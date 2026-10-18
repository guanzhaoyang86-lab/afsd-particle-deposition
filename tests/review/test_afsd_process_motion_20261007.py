"""Machine-program and finite-domain checks, not deposition validation.

Analytic volume/geometry budgets and independent quadrature check commanded
motion only. No material position, constitutive state or accepted run is made.
"""
from pathlib import Path
import copy
import json
import math
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from afsd_process_motion_20261007 import (
    ProcessProgram, domain_and_inventory_audit)


OMEGA = 10*math.pi
SPEED = .00212


def phase(name, duration, omega=0., feed=0., traverse=0., transition='smoothstep'):
    return dict(name=name, duration_s=duration, omega_rad_s=omega,
                feed_m_s=feed, traverse_m_s=traverse, transition=transition,
                status='DECLARED_ENGINEERING_ASSUMPTION',
                source='UNIT_TEST_MANUFACTURED_MACHINE_PROGRAM_NOT_VIDEO_CALIBRATION')


def specification():
    return [phase('spinup', .020, OMEGA),
            phase('feed_conditioning', .020, OMEGA, SPEED),
            phase('traverse_ramp', .020, OMEGA, SPEED, SPEED),
            phase('deposit', .100, OMEGA, SPEED, SPEED, 'hold')]


def geometry():
    source = ROOT/'outputs/contact_material_search_20261004/root/resolved_config.json'
    p = json.loads(source.read_text(encoding='utf-8'))['parameters']
    p.update(substrate_width=.010, bore_width=.002504)
    return p


class ProcessMotionTests(unittest.TestCase):
    def setUp(self):
        self.program = ProcessProgram(specification(), .0015)

    def close(self, actual, expected, atol=2e-15):
        np.testing.assert_allclose(actual, expected, atol=atol, rtol=0.)

    def test_spinup_moves_only_rotation_axis(self):
        for t in (0., .005, .010, .015, .020):
            with self.subTest(time=t):
                self.close(self.program.integral(t)[1:], [0., 0.], atol=0.)
                self.close(self.program.velocity(t)[1:], [0., 0.], atol=0.)
                self.close(self.program.global_pose(t)['origin'], [0., 0., .0015], atol=0.)
        self.assertEqual(self.program.step(0., .020).feed_displacement_m, 0.)
        self.assertGreater(self.program.integral(.020)[0], 0.)

    def test_feed_conditioning_does_not_start_traverse(self):
        self.close(self.program.velocity(.030), [OMEGA, SPEED/2, 0.])
        self.assertGreater(self.program.integral(.030)[1], 0.)
        self.assertEqual(self.program.integral(.040)[2], 0.)
        self.close(self.program.velocity(.050), [OMEGA, SPEED, SPEED/2])

    def test_cumulative_commands_have_analytic_phase_totals(self):
        # Every smooth ramp contributes half its change times its duration.
        expected = [OMEGA*.150, SPEED*.130, SPEED*.110]
        self.close(self.program.integral(self.program.end_s), expected)
        self.close(self.program.integral(.020), [OMEGA*.010, 0., 0.])
        self.close(self.program.integral(.040), [OMEGA*.030, SPEED*.010, 0.])
        self.close(self.program.integral(.060), [OMEGA*.050, SPEED*.030, SPEED*.010])

    def test_integral_matches_independent_gaussian_quadrature(self):
        # Four-node Gauss integration exactly integrates the cubic speed law;
        # it does not use the implementation's integrated quartic expression.
        nodes, weights = np.polynomial.legendre.leggauss(4)
        for p in self.program.phases:
            for left, right in ((.07, .81), (.0, .5), (.5, 1.)):
                a = p['start']+left*p['duration']
                b = p['end'] if right == 1. else p['start']+right*p['duration']
                velocities = np.array([self.program.velocity((a+b)/2+(b-a)*z/2)
                                       for z in nodes])
                numerical = (b-a)/2*np.einsum('i,ij->j', weights, velocities)
                with self.subTest(phase=p['name'], left=left, right=right):
                    self.close(self.program.integral(b)-self.program.integral(a), numerical)

    def test_position_and_speed_are_continuous_at_phase_boundaries(self):
        for p in self.program.phases[:-1]:
            t = p['end']
            before = np.nextafter(t, -math.inf)
            after = np.nextafter(t, math.inf)
            with self.subTest(time=t):
                self.close(self.program.integral(before), self.program.integral(after))
                self.close(self.program.velocity(before), self.program.velocity(after), atol=2e-13)
                self.close(self.program.global_pose(before)['origin'],
                           self.program.global_pose(after)['origin'])
                self.close(self.program.global_pose(before)['rotation'],
                           self.program.global_pose(after)['rotation'])

    def test_cross_phase_step_must_be_split_and_lands_on_boundary(self):
        for t, boundary in ((.015, .020), (.035, .040), (.055, .060)):
            with self.subTest(time=t):
                with self.assertRaises(ValueError):
                    self.program.step(t, .010)
                dt = self.program.bounded_dt(t, .010)
                self.assertEqual(dt, boundary-t)
                interval = self.program.step(t, dt)
                self.assertEqual(interval.t1, boundary)
                self.assertNotEqual(interval.phase_name, self.program.phase(boundary)['name'])

    def test_partitioned_commands_preserve_endpoints(self):
        feed = []; angles = []; travel = []
        for p in self.program.phases:
            edges = np.linspace(p['start'], p['end'], 9)
            for a, b in zip(edges[:-1], edges[1:]):
                step = self.program.step(float(a), float(b-a))
                record = step.to_dict()
                feed.append(-step.feed_displacement_m)
                angles.append(record['unwrapped_angle1_rad']-record['unwrapped_angle0_rad'])
                travel.append(record['origin1_m'][0]-record['origin0_m'][0])
        self.close([math.fsum(angles), math.fsum(feed), math.fsum(travel)],
                   self.program.integral(self.program.end_s))

    def test_step_path_has_consistent_chord_rotation_and_axial_command(self):
        interval = self.program.step(.043, .009)
        first = self.program.global_pose(interval.t0)
        last = self.program.global_pose(interval.t1)
        for t, expected in ((interval.t0, first), (interval.t1, last)):
            self.close(interval.pose(t)['origin'], expected['origin'], atol=0.)
            self.close(interval.pose(t)['rotation'], expected['rotation'], atol=0.)
        for fraction in (.2, .5, .8):
            t = interval.t0+fraction*interval.dt_s
            pose = interval.pose(t)
            self.close(pose['origin'], (1-fraction)*first['origin']+fraction*last['origin'])
            theta = (1-fraction)*interval.start[0]+fraction*interval.end[0]
            self.close(pose['rotation']@np.array([1., 0., 0.]),
                       [math.cos(theta), math.sin(theta), 0.])
            self.close(pose['rotation'].T@pose['rotation'], np.eye(3))
            self.assertAlmostEqual(np.linalg.det(pose['rotation']), 1., places=14)
            self.assertAlmostEqual(interval.prescribed_feed_at(t),
                                   fraction*interval.feed_displacement_m, places=17)

    def test_state_metadata_does_not_prescribe_free_material_trajectories(self):
        source = specification(); before = copy.deepcopy(source)
        program = ProcessProgram(source, .0015)
        snapshot = program.to_dict()
        self.assertFalse(snapshot['free_material_trajectory_prescribed'])
        self.assertEqual(source, before)
        source[0]['omega_rad_s'] = 999.
        snapshot['phases'][0]['omega_rad_s'] = 777.
        self.assertEqual(program.specification[0]['omega_rad_s'], OMEGA)
        self.close(program.integral(0.), [0., 0., 0.], atol=0.)
        self.close(program.velocity(0.), [0., 0., 0.], atol=0.)
        self.assertNotIn('material_positions', program.to_dict())

    def test_duplicate_names_missing_provenance_and_velocity_jumps_rejected(self):
        modifications = [lambda s: s[1].update(name=s[0]['name']),
                         lambda s: s[0].update(name=''),
                         lambda s: s[0].update(status='MEASURED_FROM_VIDEO'),
                         lambda s: s[0].update(source=''),
                         lambda s: s[0].update(transition='hold'),
                         lambda s: s[0].update(transition='unregistered')]
        for index, change in enumerate(modifications):
            bad = specification(); change(bad)
            with self.subTest(case=index), self.assertRaises(ValueError):
                ProcessProgram(bad, .0015)

    def test_nonfinite_negative_and_zero_inputs_rejected(self):
        for key in ('duration_s', 'omega_rad_s', 'feed_m_s', 'traverse_m_s'):
            for value in (math.nan, math.inf, -math.inf, -1.):
                bad = specification(); bad[0][key] = value
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    ProcessProgram(bad, .0015)
        bad = specification(); bad[0]['duration_s'] = 0.
        with self.assertRaises(ValueError):
            ProcessProgram(bad, .0015)
        for height in (0., -.001, math.nan, math.inf):
            with self.subTest(height=height), self.assertRaises(ValueError):
                ProcessProgram(specification(), height)
        with self.assertRaises(ValueError):
            ProcessProgram([], .0015)

    def test_time_interval_bounds_and_ambiguous_rotation_rejected(self):
        for t in (-.001, self.program.end_s+.001, math.nan, math.inf):
            for function in (self.program.integral, self.program.velocity, self.program.global_pose):
                with self.subTest(t=t, function=function.__name__), self.assertRaises(ValueError):
                    function(t)
        for dt in (0., -.001, math.nan, math.inf):
            with self.subTest(dt=dt), self.assertRaises(ValueError):
                self.program.step(0., dt)
        with self.assertRaises(ValueError):
            self.program.bounded_dt(self.program.end_s, .001)
        interval = self.program.step(0., .010)
        for t in (-.001, .011, math.nan):
            with self.subTest(t=t), self.assertRaises(ValueError):
                interval.pose(t)
        fast = ProcessProgram([phase('fast', 1., 10.)], .0015)
        with self.assertRaises(ValueError):
            fast.step(0., 1.)


class FiniteDomainTests(unittest.TestCase):
    def setUp(self):
        self.parameters = geometry()

    def program_for_distance(self, distance=0., feed=0.):
        # A one-second rest-to-speed ramp has displacement target_speed/2.
        return ProcessProgram([phase('budget', 1., feed=2*feed, traverse=2*distance)],
                              self.parameters['layer'])

    def test_ten_mm_domain_has_one_mm_command_margin(self):
        audit = domain_and_inventory_audit(self.parameters, self.program_for_distance(.001))
        self.assertTrue(audit['passed'])
        self.assertAlmostEqual(audit['required_edge_buffer_m'], .0013, places=15)
        self.assertAlmostEqual(audit['final_forward_edge_margin_m'], .0013, places=15)
        self.assertEqual(audit['particle_count'], 1350)
        over = domain_and_inventory_audit(self.parameters, self.program_for_distance(.00100001))
        self.assertFalse(over['passed'])
        self.assertEqual(over['failures'], ['SUBSTRATE_DOES_NOT_COVER_FULL_PATH_PLUS_KERNEL_BUFFER'])

    def test_eight_mm_domain_has_no_traverse_margin(self):
        self.parameters['substrate_width'] = .008
        initial = domain_and_inventory_audit(self.parameters, self.program_for_distance())
        moved = domain_and_inventory_audit(self.parameters, self.program_for_distance(.000001))
        self.assertTrue(initial['passed'])
        self.assertFalse(moved['passed'])
        self.assertIn('SUBSTRATE_DOES_NOT_COVER_FULL_PATH_PLUS_KERNEL_BUFFER', moved['failures'])

    def test_actuator_cap_reserve_limits_finite_feed(self):
        # 3 mm bar minus 1.5 mm shoulder, .5 mm actuator cap and .25 mm
        # half-particle reserve leaves .75 mm maximum commanded feed.
        exact = domain_and_inventory_audit(self.parameters, self.program_for_distance(feed=.00075))
        self.assertTrue(exact['passed'])
        self.assertAlmostEqual(exact['usable_axial_command_m'], .00075, places=15)
        depleted = domain_and_inventory_audit(self.parameters, self.program_for_distance(feed=.00075001))
        self.assertFalse(depleted['passed'])
        self.assertEqual(depleted['failures'], ['FINITE_FEEDSTOCK_ACTUATOR_RESERVE_EXHAUSTED'])
        self.assertEqual(depleted['replenishment'], 'FINITE_BAR_NO_AUTOMATIC_INSERTION')

    def test_declared_example_fits_without_predicting_material_motion(self):
        original = copy.deepcopy(self.parameters)
        audit = domain_and_inventory_audit(self.parameters, ProcessProgram(specification(), .0015))
        self.assertTrue(audit['passed'])
        self.assertEqual(self.parameters, original)
        self.assertEqual(audit['free_material_transport'], 'NOT_PREDICTED_BY_COMMAND_INTEGRALS')
        self.assertEqual(audit['domain_convergence'], 'NOT_TESTED')
        self.assertAlmostEqual(audit['programmed_feed_m'], .0002756, places=15)
        self.assertAlmostEqual(audit['programmed_traverse_m'], .0002332, places=15)

    def test_nonphysical_geometry_rejected(self):
        for key, value in (('dx', math.nan), ('substrate_width', -.01),
                           ('substrate_depth', .0017), ('bore_width', .0024),
                           ('tool_radius', .001), ('tool_top', .0015),
                           ('bar_length', .002)):
            bad = copy.deepcopy(self.parameters); bad[key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                domain_and_inventory_audit(bad, self.program_for_distance())


if __name__ == '__main__':
    unittest.main()

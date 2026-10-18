"""Independent literals and <= two-step tests; no full impact trajectory here."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import probe_normal_contact as impact
from afsd_backend.pysph.normal_penalty import normal_state, tool_state, moving_tool_rows
from afsd_backend.pysph.contact import WallContact
from afsd_scaffold.vtk_legacy import read_point_cloud


class NormalPenaltyChecks(unittest.TestCase):
    def setUp(self):
        self.p = dict(dx=.2, traverse=0., omega=0., layer=1., tool_top=2., tool_radius=2.,
                      bore_width=1., substrate_depth=1., base_temperature=293., tool_temperature=293.)

    def test_01_independent_piecewise_force_potential_and_signed_gap(self):
        s = normal_state(np.array([.2, .1, .075, -.025]), .1, 100.)
        np.testing.assert_allclose(s['signed_penetration'], [-.1, 0., .025, .125], atol=1e-16)
        np.testing.assert_array_equal(s['active'], [False, False, True, True])
        np.testing.assert_allclose(s['normal_force_N'], [0., 0., 2.5, 12.5], atol=1e-14)
        np.testing.assert_allclose(s['potential_J'], [0., 0., .03125, .78125], atol=1e-14)

    def test_02_actual_pinned_wall_normal_force_matches_literals(self):
        wall = WallContact(100., .2)
        for gap, force, U, delta in [(.2, 0., 0., 0.), (.1, 0., 0., 0.),
                                     (.075, 2.5, .03125, .025)]:
            got = wall.force(gap, [0., 0., 1.], [0., 0., -.3], .1, 1., .01)
            np.testing.assert_allclose(got[0], [0., 0., force], atol=1e-12)
            self.assertEqual(got[1], 0.)
            self.assertAlmostEqual(got[2], U, delta=1e-12)
            self.assertAlmostEqual(got[3], delta, delta=1e-12)
            self.assertAlmostEqual(float(normal_state(gap, .1, 100.)['normal_force_N']), force, delta=1e-12)

    def test_03_force_is_negative_potential_gap_gradient(self):
        # Independent centered derivative of literal U=50*(.1-gap)^2.
        gap, epsilon = .075, 1e-7
        dUdg = (50*(.1-gap-epsilon)**2-50*(.1-gap+epsilon)**2)/(2*epsilon)
        self.assertAlmostEqual(float(normal_state(gap, .1, 100.)['normal_force_N']), -dUdg, delta=1e-9)

    def test_04_nonfinite_and_nonpositive_parameters_fail(self):
        for args in [(np.nan, 1., 4.), (1., 0., 4.), (1., 1., -4.),
                     (1., np.inf, 4.), (1., 1., np.inf)]:
            with self.assertRaises(ValueError): normal_state(*args)

    def test_05_open_wall_no_normal_attraction(self):
        s = normal_state([1., 1.1, 2.], 1., 4.)
        np.testing.assert_array_equal(s['normal_force_N'], [0., 0., 0.])
        np.testing.assert_array_equal(s['potential_J'], [0., 0., 0.])
        self.assertLess(s['signed_penetration'][-1], 0.)

    def test_06_lift_position_and_velocity_piecewise_literal(self):
        before = tool_state(.2, self.p, .25, .4)
        at = tool_state(.25, self.p, .25, .4)
        after = tool_state(.5, self.p, .25, .4)
        self.assertEqual(before['layer'], 1.); self.assertEqual(at['layer'], 1.)
        self.assertAlmostEqual(after['layer'], 1.1)
        self.assertAlmostEqual(after['tool_top'], 2.1)
        np.testing.assert_array_equal(before['translation_velocity'], [0., 0., 0.])
        np.testing.assert_array_equal(at['translation_velocity'], [0., 0., .4])
        epsilon = 1e-6
        dz = (tool_state(.5+epsilon, self.p, .25, .4)['layer']-
              tool_state(.5-epsilon, self.p, .25, .4)['layer'])/(2*epsilon)
        self.assertAlmostEqual(dz, .4, delta=1e-9)

    def test_07_tool_rows_move_both_velocities_but_base_remains_fixed(self):
        points = np.array([[1., 0., 1.05], [0., 0., -.95], [.45, .2, 1.5]])
        ids = np.array([11, 12, 13])
        rows = moving_tool_rows(points, ids, .5, self.p, .25, .4)
        self.assertEqual(set(rows['feature_id']), {0, 1, 3})
        for a, feature in enumerate(rows['feature_id']):
            expected = [0., 0., 0.] if feature == 0 else [0., 0., .4]
            np.testing.assert_array_equal(rows['wall_velocity'][a], expected)
            np.testing.assert_array_equal(rows['wall_velocity_at_face'][a], expected)
            np.testing.assert_array_equal(rows['translation_velocity'][a], expected)
        self.assertEqual(self.p['layer'], 1.); self.assertEqual(self.p['tool_top'], 2.)

    def test_08_rotating_lifting_velocity_same_application_location(self):
        p = dict(self.p, omega=2., traverse=.3)
        rows = moving_tool_rows(np.array([[1., 0., 1.05]]), np.array([11]), .5, p, .25, .4)
        np.testing.assert_allclose(rows['wall_velocity'], [[.3, 1.7, .4]], atol=1e-14)
        np.testing.assert_allclose(rows['wall_velocity_at_face'], [[.3, 1.7, .4]], atol=1e-14)
        np.testing.assert_allclose(rows['body_origin'], [[.15, 0., 1.1]], atol=1e-14)

    def test_09_rotating_square_bore_stays_open(self):
        rows = moving_tool_rows(np.array([[0., 0., 1.5]]), np.array([7]), .5, self.p, .25, .4)
        self.assertEqual(len(rows['particle_id']), 0)
        self.assertEqual(rows['wall_velocity'].shape, (0, 3))

    def test_10_tool_lift_guards(self):
        for args in [(-1., .2, .4), (.1, -.1, .4), (.1, .2, -.4), (.1, .2, np.nan)]:
            with self.assertRaises(ValueError): tool_state(args[0], self.p, args[1], args[2])

    def test_11_literal_first_static_and_moving_kicks(self):
        for speed, z, v, work, initial_K in [(0., .995000125, -.99995, 0., .5),
                                            (.2, .996000125, -.79995, .00001, .32)]:
            b = impact.NormalImpact(.005, speed); record = b.advance()
            np.testing.assert_allclose(record['position_m'], [0., 0., z], atol=1e-14)
            np.testing.assert_allclose(record['velocity_m_s'], [0., 0., v], atol=1e-14)
            self.assertAlmostEqual(record['cumulative']['wall_work_J'], work, delta=1e-16)
            self.assertAlmostEqual(record['cumulative']['normal_impulse_N_s'], .00005, delta=1e-16)
            self.assertEqual(record['cumulative']['friction_raw_J'], 0.)
            self.assertAlmostEqual(record['kinetic_J']-initial_K,
                                   record['cumulative']['normal_work_J'], delta=1e-14)
            self.assertAlmostEqual(record['global_energy_residual_J'], -1.24996875e-9, delta=1e-14)

    def test_12_independent_discrete_normal_spring_defect_not_heat(self):
        pa = impact.get_particle_array(name='solid', x=[0.], y=[0.], z=[.96], w=[-.3], m=[1.])
        for name in ('ax', 'ay', 'az', 'au', 'av', 'aw', 'arho', 'ae', 'e'):
            if name not in pa.properties: pa.add_property(name)
        trace = {}
        def callback(t, dt):
            state = normal_state(pa.z[0], 1., 100.)
            result = impact.kick([0., 0., pa.w[0]], 1., np.zeros(3), [[0., 0., 1.]],
                                 [float(state['normal_force_N'])], [[0., 0., 0.]], 0., dt)
            pa.aw[:] = (result['velocity'][2]-pa.w[0])/dt
            trace.update(result)
        impact.HostLeapfrog(pa, callback).advance(0., .01)
        self.assertAlmostEqual(pa.w[0], -.2585, delta=1e-14)
        self.assertAlmostEqual(1.-pa.z[0], .0427925, delta=1e-14)
        dU = 50*(1.-pa.z[0])**2-.08
        normal_work = sum(a['normal_work_J'].sum() for a in trace['intervals'])
        defect = normal_work+dU
        self.assertAlmostEqual(defect, -2.89721875e-5, delta=1e-14)
        self.assertAlmostEqual(defect, .0025*(.5*pa.w[0]**2-.045), delta=1e-14)
        self.assertEqual(float(trace['friction_heat_J'].sum()), 0.)

    def test_13_observation_labels_do_not_change_two_step_motion(self):
        a, b = impact.NormalImpact(.005), impact.NormalImpact(.005)
        b.pa.tracer_type[:] = 999; b.pa.as_rs_tag[:] = 0
        for _ in range(2): a.advance(); b.advance()
        np.testing.assert_array_equal(impact.positions(a.pa), impact.positions(b.pa))
        np.testing.assert_array_equal(impact.velocities(a.pa), impact.velocities(b.pa))
        self.assertEqual(a.cumulative, b.cumulative)

    def test_14_nonzero_checkpoint_restart_identity(self):
        a, b = impact.NormalImpact(.0025, .2), impact.NormalImpact(.0025, .2)
        a.advance(); a.advance(); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; b.save_checkpoint(path)
            c = impact.NormalImpact.from_checkpoint(path); c.advance()
        self.assertEqual(a.cumulative, c.cumulative)
        self.assertEqual(a.time, c.time)
        for name in ('x', 'y', 'z', 'u', 'v', 'w', 'particle_id', 'tracer_type', 'as_rs_tag'):
            np.testing.assert_array_equal(getattr(a.pa, name), getattr(c.pa, name))

    def test_15_failure_rolls_back_and_marks_candidate_rejected(self):
        b = impact.NormalImpact(.005); b.advance(); before = b.record()
        with patch.object(impact, 'kick', side_effect=ValueError('INJECTED_FAILURE')):
            with self.assertRaisesRegex(ValueError, 'INJECTED_FAILURE'): b.advance()
        np.testing.assert_array_equal(b.record()['position_m'], before['position_m'])
        self.assertEqual(b.time, before['time_s']); self.assertEqual(b.cumulative, before['cumulative'])
        self.assertEqual(b.rejected_candidate['classification'], 'REJECTED')
        self.assertFalse(b.rejected_candidate['restart_allowed'])

    def test_16_output_has_no_fabricated_solid_or_thermal_fields(self):
        b = impact.NormalImpact(.005); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'frame.vtk'; b.vtk(path); result = read_point_cloud(path)
        self.assertEqual(set(result['scalars']), {'particle_id', 'tracer_type', 'as_rs_tag', 'mass_kg'})
        self.assertEqual(set(result['vectors']), {'velocity_m_per_s'})
        self.assertIn('NOT_AFSD', result['title'])

    def test_17_reject_nonregistered_case_and_explicit_configuration_provenance(self):
        with self.assertRaises(ValueError): impact.NormalImpact(.001)
        b = impact.NormalImpact(.005, .2)
        self.assertEqual(b.config['mass_kg'], 1.)
        self.assertEqual(b.config['initial_velocity_m_s'], [0., 0., -.8])
        self.assertEqual(b.config['source_configuration']['parameters']['friction']['value'], 0)
        self.assertIn('docs/NORMAL_CONTACT_CONTRACT.md', b.source_hashes)
        self.assertEqual(b.config['field_status']['temperature'], 'NOT_SOLVED')

    def test_19_rejected_or_missing_acceptance_cannot_restart(self):
        b = impact.NormalImpact(.005); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; b.save_checkpoint(path)
            pa, state = impact.restore(path)
            for key in ('accepted', 'restart_allowed'):
                for value in (False, None):
                    altered = dict(state); altered[key] = value
                    with self.assertRaisesRegex(ValueError, 'only accepted'):
                        impact.NormalImpact(.005, pa=pa, state=altered)

    def test_18_analytic_reference_literals_are_independent(self):
        # Explicit continuous reference: no normal_state/kick integration to
        # generate the expected solution. This is reference algebra, not 400 steps.
        with patch.object(impact, 'normal_state', side_effect=AssertionError('reference called law')):
            with patch.object(impact, 'kick', side_effect=AssertionError('reference called kick')):
                stop = np.pi/2
                at_max_delta, at_max_U = .5*np.sin(2*np.pi/4), .5*np.sin(2*np.pi/4)**2
                z_static, z_moving = 3.-stop, 3.-stop+.4
                initial_K, final_K, wall_work = .32, .72, .4
        self.assertAlmostEqual(at_max_delta, .5); self.assertAlmostEqual(at_max_U, .5)
        self.assertAlmostEqual(z_static, 1.4292036732051034)
        self.assertAlmostEqual(z_moving, 1.8292036732051034)
        self.assertAlmostEqual(final_K-initial_K, wall_work)


if __name__ == '__main__': unittest.main(verbosity=2)

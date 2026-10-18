"""Manufactured observer identities; no deformation or bonding validation."""
import copy
import json
import math
from types import SimpleNamespace
import unittest
import numpy as np
from afsd_process_observer_20261007 import observe_process
from afsd_backend.pysph.interval_acceptance import initial_bore_crossing_state


class Motion:
    t0, t1, dt_s = 0., 1., 1.
    def __init__(self, translation=None, omega=2.):
        self.translation = np.zeros(3) if translation is None else np.asarray(translation, float)
        self.omega = omega
    def pose(self, t):
        return dict(origin=np.array([0., 0., 1.])+t*self.translation, rotation=np.eye(3))
    def to_dict(self):
        return dict(average_omega_rad_s=self.omega)


class ProcessObserverTests(unittest.TestCase):
    def fixture(self, x, mass=None, crossed=(), ineligible=(), actuator=(), motion=None):
        x = np.asarray(x, float); n = len(x)
        mass = np.ones(n) if mass is None else np.asarray(mass, float)
        ids = np.arange(1, n+1)
        pa = SimpleNamespace(particle_id=ids, component_id=np.ones(n, int), m=mass,
            actuator=np.isin(ids, actuator).astype(int))
        initial = x.copy(); initial[:, 2] = 2.
        initial[np.isin(ids, ineligible), 2] = .5
        motion = Motion() if motion is None else motion
        h = initial_bore_crossing_state(initial, ids, pa.component_id, mass, 0.,
                                        motion.pose(0.), 2.)
        h['last_time_s'] = 1.; h['crossing_ids'] = sorted(crossed)
        h['crossing_mass_kg'] = math.fsum(mass[pid-1] for pid in crossed)
        h['events'] = [dict(particle_id=pid, time_s=.5, step_fraction=.5,
            position_m=[0., 0., 1.], tool_local_position_m=[0., 0., 0.], mass_kg=float(mass[pid-1]))
            for pid in sorted(crossed)]
        state = dict(bore_observation=h)
        result = dict(position=x, displacement=np.zeros_like(x),
                      contact=dict(force_by_feature_N=np.zeros((n, 8, 3))))
        p = dict(tool_radius=3., bore_width=2., dx=.2)
        return pa, state, result, motion, p

    def test_mass_weighted_spin_subtracts_translation_and_excludes_actuator_and_axis(self):
        args = self.fixture([[2., 1., 2.], [3., 1., 2.], [1., 1., 2.], [2., 1., 2.]],
            mass=[2., 1., 7., 50.], actuator=[4], motion=Motion([1., 1., 0.]))
        args[2]['displacement'][:] = [[1., 2., 0.], [1., 7., 0.], [1., 1., 0.], [1., 101., 0.]]
        r = observe_process(*args)['rotation']
        self.assertAlmostEqual(r['mass_weighted_material_omega_rad_s'], 5./3.)
        self.assertAlmostEqual(r['inertia_weighted_material_omega_rad_s'], 14./6.)
        self.assertEqual(r['axis_excluded_population']['particle_ids'], [3])
        self.assertEqual(r['resolved_population']['particle_ids'], [1, 2])

    def test_feature_torque_reaction_and_no_shoulder_force_leak(self):
        args = self.fixture([[2., 0., 2.], [0., 1., 2.]])
        args[2]['contact']['force_by_feature_N'][0, 2, 1] = 3.
        args[2]['contact']['force_by_feature_N'][1, 5, 0] = 2.
        args[2]['contact']['force_by_feature_N'][0, 0, 1] = 1000.
        r = observe_process(*args)['bore_wall_torque']
        self.assertEqual(r['on_all_deformable_material_z_Nm'], 4.)
        self.assertEqual(r['reaction_on_tool_z_Nm'], -4.)

    def test_only_eligible_crossed_count_for_spreading_and_retention(self):
        args = self.fixture([[.5, 0., .5], [1.5, 0., .5], [-3.2, 0., .5],
                             [1.6, 0., .5], [-3.2, 0., -.1]], mass=[1., 2., 3., 4., 5.],
                            crossed=[1, 2, 3, 5], ineligible=[4])
        r = observe_process(*args)['supplied_material']
        self.assertEqual(r['eligible_crossed']['mass_kg'], 11.)
        self.assertEqual(r['crossed_under_annular_shoulder']['particle_ids'], [2])
        self.assertEqual(r['geometric_retention_candidate']['particle_ids'], [3])
        self.assertFalse(r['geometric_retention_is_metallurgical_bonding'])

    def test_missing_feature_forces_are_not_reported_as_zero(self):
        args = self.fixture([[.5, 0., 2.]])
        args[2]['contact'] = {}
        self.assertEqual(observe_process(*args)['bore_wall_torque']['status'],
                         'NOT_AVAILABLE_NO_SAVED_FEATURE_FORCES')

    def test_rotated_square_uses_tool_coordinates(self):
        class RotatedMotion(Motion):
            def pose(self, t):
                p = super().pose(t)
                c = np.sqrt(.5)
                p['rotation'] = np.array([[c, -c, 0.], [c, c, 0.], [0., 0., 1.]])
                return p
        args = self.fixture([[1.1, 0., .5]], crossed=[1], motion=RotatedMotion())
        r = observe_process(*args)['supplied_material']
        self.assertEqual(r['crossed_inside_square_projection']['particle_ids'], [1])
        self.assertEqual(r['crossed_under_annular_shoulder']['particle_ids'], [])

    def test_stopped_tool_ratio_is_null_and_result_is_json_friendly_and_pure(self):
        args = self.fixture([[.5, 0., 2.]], motion=Motion(omega=0.))
        old = copy.deepcopy(args[1]); x = args[2]['position'].copy()
        r = observe_process(*args)
        self.assertIsNone(r['rotation']['mass_weighted_material_to_tool_ratio'])
        json.dumps(r, allow_nan=False)
        self.assertEqual(old, args[1]); np.testing.assert_array_equal(x, args[2]['position'])

    def test_inconsistent_crossing_history_is_rejected(self):
        args = self.fixture([[.5, 0., .5]], crossed=[1])
        args[1]['bore_observation']['crossing_mass_kg'] = 9.
        with self.assertRaises(ValueError):
            observe_process(*args)


if __name__ == '__main__':
    unittest.main()

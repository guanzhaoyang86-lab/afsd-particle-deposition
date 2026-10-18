"""Synthetic fixtures test passive bookkeeping, never claim physical validation."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from afsd_deposition_progress_20261010 import frame_from_saved, summarize_progress

P = dict(tool_radius=3., bore_width=2., dx=.5, initial_substrate_top_z_m=0.)


def frame(t=0., crossed=(), positions=None):
    return dict(time_s=t, ids=np.array([10, 11, 12]), component_ids=np.array([1, 1, 0]),
        position=np.array(positions if positions is not None else [[0., 0., 2.2], [0., 0., 1.], [0., 0., -.2]]),
        mass=np.array([2., 3., 5.]), ep=np.zeros(3), eligible_ids=[10], crossing_ids=list(crossed),
        tool_origin=np.array([0., 0., 2.]), tool_rotation=np.eye(3))


class ProgressTests(unittest.TestCase):
    def test_gap_initial_lower_bar_excluded(self):
        result = summarize_progress([frame()], P)
        row = result['final']
        self.assertAlmostEqual(row['nearest_eligible_uncrossed_gap_m'], .2)
        self.assertEqual(row['nearest_eligible_uncrossed_particle_id'], 10)
        self.assertEqual(row['supplied_under_shoulder']['particle_count'], 0)
        self.assertIsNone(row['new_crossing_ids_since_previous_sample'])
        self.assertFalse(result['metallurgical_bonding_proven'])
        self.assertEqual(result['physical_solver_evaluations'], 0)

    def test_gap_zero_and_negative_do_not_invent_crossing(self):
        for z in (2., 1.9):
            f = frame(); f['position'][0, 2] = z
            row = summarize_progress([f], P)['final']
            self.assertAlmostEqual(row['nearest_eligible_uncrossed_gap_m'], z-2.)
            self.assertEqual(row['supplied_crossed']['particle_count'], 0)
            self.assertEqual(row['pending_at_or_below_outlet_ids'], [10])

    def test_crossing_then_radial_spread_and_residence(self):
        f0 = frame()
        f1 = frame(1., [10]); f1['position'][0] = [0., 0., 1.9]
        f2 = frame(2., [10]); f2['position'][0] = [1.2, 0., 1.8]
        f3 = frame(3., [10]); f3['position'][0] = [1.4, 0., 1.7]; f3['ep'][0] = .2
        result = summarize_progress([f0, f1, f2, f3], P)
        self.assertEqual(result['history'][1]['new_crossing_ids_since_previous_sample'], [10])
        self.assertEqual(result['history'][1]['new_crossing_mass_since_previous_sample_kg'], 2.)
        self.assertEqual(result['final']['supplied_under_shoulder']['mass_kg'], 2.)
        residence = result['final']['under_shoulder_residence'][0]
        self.assertEqual((residence['consecutive_saved_frames'], residence['sampled_residence_s']), (2, 1.))
        self.assertEqual(result['first_observed_supplied_under_shoulder_time_s'], 2.)
        self.assertAlmostEqual(result['final']['radial_spread']['maximum_radius_change_since_first_saved_crossing_m'], 1.4)
        self.assertEqual(result['final']['radial_spread']['maximum_ep_increment_since_first_saved_crossing'], .2)
        self.assertIsNone(result['final']['nearest_eligible_uncrossed_gap_m'])

    def test_leave_reenter_resets_consecutive_residence(self):
        fs = [frame(t, [10]) for t in (1., 2., 3.)]
        fs[0]['position'][0] = [1.2, 0., 1.8]
        fs[1]['position'][0] = [0., 0., 1.8]
        fs[2]['position'][0] = [1.2, 0., 1.8]
        result = summarize_progress(fs, P)
        self.assertEqual(result['history'][1]['left_under_shoulder_ids_since_previous_sample'], [10])
        self.assertEqual(result['history'][1]['ever_observed_under_shoulder_mass_kg'], 2.)
        self.assertEqual(result['final']['under_shoulder_residence'][0]['consecutive_saved_frames'], 1)
        self.assertEqual(result['final']['under_shoulder_residence'][0]['sampled_residence_s'], 0.)

    def test_particle_order_is_irrelevant(self):
        a, b = frame(), frame(1.)
        for key in ('ids', 'component_ids', 'position', 'mass', 'ep'):
            b[key] = b[key][[2, 0, 1]]
        result = summarize_progress([a, b], P)
        self.assertEqual(result['final']['nearest_eligible_uncrossed_particle_id'], 10)

    def test_rotated_square_mask_uses_tool_pose(self):
        f = frame(1., [10]); f['position'][0] = [1.1, 0., 1.8]
        self.assertEqual(summarize_progress([f], P)['final']['supplied_under_shoulder']['particle_count'], 1)
        angle = np.pi/4
        f['tool_rotation'] = np.array([[np.cos(angle), -np.sin(angle), 0.],
            [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
        self.assertEqual(summarize_progress([f], P)['final']['supplied_under_shoulder']['particle_count'], 0)

    def test_no_mutation(self):
        f = frame(1., [10]); saved = deepcopy(f)
        summarize_progress([f], P)
        for key, value in saved.items():
            np.testing.assert_array_equal(f[key], value)

    def test_history_and_identity_rejection(self):
        mutations = [lambda f: f.update(time_s=0.), lambda f: f['mass'].__setitem__(0, 4.),
            lambda f: f.update(eligible_ids=[10, 11]), lambda f: f['component_ids'].__setitem__(0, 0),
            lambda f: f['ids'].__setitem__(0, 11)]
        for mutate in mutations:
            with self.subTest(mutation=mutate):
                f = frame(1.); mutate(f)
                with self.assertRaises(ValueError):
                    summarize_progress([frame(), f], P)
        with self.assertRaises(ValueError):
            summarize_progress([frame(0., [10]), frame(1.)], P)
        with self.assertRaises(ValueError):
            summarize_progress([frame(0., [11])], P)

    def test_finite_and_pose_validation(self):
        for key in ('position', 'mass', 'ep', 'tool_origin', 'tool_rotation'):
            with self.subTest(key=key):
                f = frame(); f[key].flat[0] = np.nan
                with self.assertRaises(ValueError): summarize_progress([f], P)
        f = frame(); f['tool_rotation'][0, 0] = -1.
        with self.assertRaises(ValueError): summarize_progress([f], P)
        with self.assertRaises(ValueError): summarize_progress([], P)

    def test_full_clearance_and_no_false_shoulder(self):
        f = frame(1., [10]); f['position'][0] = [3.3, 0., 1.]
        row = summarize_progress([f], P)['final']
        self.assertEqual(row['supplied_under_shoulder']['mass_kg'], 0.)
        self.assertEqual(row['supplied_clear_of_shoulder_projection']['mass_kg'], 2.)
        self.assertFalse(row['geometric_retention_is_metallurgical_bonding'])

    def test_saved_observer_agreement(self):
        f = frame()
        row = summarize_progress([f], P)['final']
        keys = [('eligible', 'eligible'), ('eligible_crossed', 'supplied_crossed'),
            ('eligible_not_yet_crossed', 'supplied_not_yet_crossed'),
            ('crossed_under_annular_shoulder', 'supplied_under_shoulder'),
            ('geometric_retention_candidate', 'supplied_clear_of_shoulder_projection')]
        f['saved_observer'] = dict(time1_s=0., supplied_material={a: deepcopy(row[b]) for a, b in keys})
        self.assertTrue(summarize_progress([f], P)['final']['saved_observer_matches'])
        f['saved_observer']['supplied_material']['eligible']['mass_kg'] += 1.
        with self.assertRaises(ValueError): summarize_progress([f], P)

    def test_saved_adapter_checks_identity_and_copies(self):
        f = frame()
        pa = SimpleNamespace(particle_id=f['ids'], component_id=f['component_ids'], m=f['mass'],
            x=f['position'][:, 0], y=f['position'][:, 1], z=f['position'][:, 2], ep=f['ep'])
        state = dict(time_s=0., step=0, parameters=P, bore_observation=dict(
            particle_ids=[10, 11, 12], component_ids=[1, 1, 0], masses_kg=[2., 3., 5.],
            last_time_s=0., bore_width_m=2., eligible_ids=[10], crossing_ids=[]))
        result = frame_from_saved(pa, state, dict(origin=f['tool_origin'], rotation=f['tool_rotation']))
        result['position'][0, 2] = -100.
        self.assertEqual(pa.z[0], 2.2)
        state['bore_observation']['masses_kg'][0] = 9.
        with self.assertRaises(ValueError):
            frame_from_saved(pa, state, dict(origin=f['tool_origin'], rotation=f['tool_rotation']))

    def test_saved_machine_pose_fallback(self):
        f = frame(1.)
        pa = SimpleNamespace(particle_id=f['ids'], component_id=f['component_ids'], m=f['mass'],
            x=f['position'][:, 0], y=f['position'][:, 1], z=f['position'][:, 2], ep=f['ep'])
        state = dict(time_s=1., step=2, parameters=P, bore_observation=dict(
            particle_ids=[10, 11, 12], component_ids=[1, 1, 0], masses_kg=[2., 3., 5.],
            last_time_s=1., bore_width_m=2., eligible_ids=[10], crossing_ids=[]),
            machine_history=[dict(t1_s=1., origin1_m=f['tool_origin'].tolist(), rotation1=f['tool_rotation'].tolist())])
        copied = frame_from_saved(pa, state)
        np.testing.assert_array_equal(copied['tool_origin'], f['tool_origin'])
        np.testing.assert_array_equal(copied['tool_rotation'], f['tool_rotation'])
        state['machine_history'][0]['t1_s'] = .5
        with self.assertRaises(ValueError): frame_from_saved(pa, state)
        state['machine_history'] = []
        with self.assertRaises(ValueError): frame_from_saved(pa, state)


if __name__ == '__main__':
    unittest.main()

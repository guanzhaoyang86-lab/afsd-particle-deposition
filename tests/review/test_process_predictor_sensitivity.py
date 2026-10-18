"""Independent synthetic postprocessing references; NOT physical solver tests.

Literal vector norms, event identities and signed work are specified here.
Only the NumPy postprocessor is imported; no solver/material reference calls.
"""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import compare_process_predictor_sensitivity as subject


class MetricReferences(unittest.TestCase):
    def test_vector_euclidean_and_final_not_component_maximum(self):
        result = subject.metric([[3., 4., 0.], [0., 0., 2.]],
                                [[0., 0., 0.], [0., 0., 4.]], 10., 1)
        self.assertEqual(result['denominator'], 5.)
        self.assertEqual(result['per_time_absolute_difference'], [5., 2.])
        self.assertEqual(result['maximum_normalized_difference'], 1.)
        self.assertEqual(result['final_normalized_difference'], .4)

    def test_tensor_frobenius_includes_off_diagonal_entries(self):
        a = np.zeros((2, 1, 3, 3))
        b = np.zeros_like(a)
        b[1, 0] = [[2., 3., 0.], [0., 0., 6.], [0., 0., 0.]]
        result = subject.metric(a, b, 100., 2)
        self.assertEqual(result['maximum_absolute_difference'], 7.)
        self.assertEqual(result['denominator'], 7.)
        self.assertEqual(result['per_time_absolute_difference'], [0., 7.])

    def test_signed_physical_final_values_not_replaced_by_magnitude(self):
        result = subject.metric([0., -2.], [0., 2.], 1.)
        self.assertEqual(result['coarse_final'], -2.)
        self.assertEqual(result['fine_final'], 2.)
        self.assertEqual(result['maximum_absolute_difference'], 4.)
        self.assertEqual(result['maximum_normalized_difference'], 2.)

    def test_near_zero_floor_is_exposed(self):
        result = subject.metric([0., 0.], [0., .00025], 1.)
        self.assertEqual(result['denominator'], .01)
        self.assertEqual(result['maximum_normalized_difference'], .025)
        self.assertEqual(result['status'], 'WITHIN_SCREEN')

    def test_five_percent_boundary_and_warning(self):
        self.assertEqual(subject.metric([0., 20.], [0., 19.], 1.)['status'], 'WITHIN_SCREEN')
        self.assertEqual(subject.metric([0., 20.], [0., 18.9], 1.)['status'], 'SENSITIVE_WARNING')

    def test_nonfinite_shape_and_nominal_rejected(self):
        for a, b, nominal in [([0.], [np.nan], 1.), ([0.], [0., 1.], 1.),
                              ([0.], [0.], 0.), ([0.], [0.], np.inf)]:
            with self.subTest(a=a, b=b, nominal=nominal), self.assertRaises(ValueError):
                subject.metric(a, b, nominal)


class TimeReferences(unittest.TestCase):
    def test_actual_frames_selected_without_interpolation(self):
        np.testing.assert_array_equal(subject.match_times([0., .004, .0125, .02, .025],
                                                         [0., .0125, .025]), [0, 2, 4])

    def test_missing_time_cannot_be_interpolated(self):
        with self.assertRaisesRegex(ValueError, 'missing or ambiguous'):
            subject.match_times([0., .01, .02], [.015])

    def test_two_close_accepted_frames_are_ambiguous(self):
        with self.assertRaisesRegex(ValueError, 'missing or ambiguous'):
            subject.match_times([0., .01, .01 + 1e-13, .02], [.01])

    def test_duplicate_nonmonotonic_and_nonfinite_rejected(self):
        for times in ([0., .01, .01], [0., .02, .01], [0., np.nan, .02]):
            with self.subTest(times=times), self.assertRaisesRegex(ValueError, 'invalid accepted'):
                subject.match_times(times, [0.])

    def test_frame_cannot_satisfy_two_requested_times(self):
        with self.assertRaisesRegex(ValueError, 'reused accepted'):
            subject.match_times([0., .01, .02], [.01, .01 + 1e-13])


def event(t, selected, released=()):
    return dict(physical_time_s=t,
                cone_margin_N=[0. if i in selected else 1. for i in range(3)],
                released_rows=[i in released for i in range(3)])


CONTACTS = dict(contact_particle_id=np.array([19, 24, 26]),
                contact_feature_id=np.array([1, 3, 2]))


class EventReferences(unittest.TestCase):
    def inventory(self, events):
        return subject.event_inventory(dict(events=events), CONTACTS)

    def test_simultaneous_group_may_split_and_reorder_without_index_pairing(self):
        coarse = self.inventory([event(.06, [0, 1])])
        fine = self.inventory([event(.0595, [1]), event(.0604, [0])])
        result = subject.compare_events(coarse, fine)
        self.assertEqual(result['status'], 'WITHIN_SCREEN')
        self.assertEqual(result['coarse_only'], [])
        self.assertEqual(result['fine_only'], [])
        self.assertFalse(result['sequence_equal'])
        self.assertEqual(len(result['coarse_groups']), 1)
        self.assertEqual(len(result['fine_groups']), 2)
        paired = {r['key']: r for r in result['paired']}
        self.assertAlmostEqual(paired['2:19:1:1']['signed_difference_s'], .0004)
        self.assertAlmostEqual(paired['2:24:3:1']['signed_difference_s'], -.0005)

    def test_missing_counterparts_are_reported_not_discarded(self):
        result = subject.compare_events(self.inventory([event(.06, [0, 1])]),
                                        self.inventory([event(.06, [0, 2])]))
        self.assertEqual(result['status'], 'EVENT_INVENTORY_MISMATCH')
        self.assertEqual(result['coarse_only'], ['2:24:3:1'])
        self.assertEqual(result['fine_only'], ['2:26:2:1'])
        self.assertEqual([p['key'] for p in result['paired']], ['2:19:1:1'])

    def test_phase_corner_uses_ending_phase_and_occurrence(self):
        result = self.inventory([event(.05, [0]), event(.05 + .5e-12, [0]),
                                 event(.05 + 2e-12, [0]), event(.1, [1]),
                                 event(.15, [2]), event(.2, [0])])
        self.assertEqual(set(result['times_by_key']),
                         {'1:19:1:1', '1:19:1:2', '2:19:1:1',
                          '2:24:3:1', '3:26:2:1', '4:19:1:1'})

    def test_already_released_contact_is_not_new_event(self):
        result = self.inventory([event(.06, [0, 1], released=[0])])
        self.assertEqual(list(result['times_by_key']), ['2:24:3:1'])

    def test_shift_above_frozen_screen_warns(self):
        result = subject.compare_events(self.inventory([event(.06, [0])]),
                                        self.inventory([event(.062, [0])]))
        self.assertEqual(result['status'], 'SENSITIVE_WARNING')
        self.assertAlmostEqual(result['maximum_paired_absolute_shift_s'], .002)

    def test_empty_nonfinite_and_invalid_time_events_rejected(self):
        invalid = [event(.06, []), event(0., [0]), event(.21, [0]), event(float('nan'), [0])]
        broken = event(.06, [0]); broken['cone_margin_N'][1] = float('nan')
        invalid.append(broken)
        for item in invalid:
            with self.subTest(item=item), self.assertRaises(ValueError):
                self.inventory([item])


def work_fixture():
    """Independent signed ledger: Q=1.5+1.75, Wwall=2+2, Wmat=.5+.25 J."""
    h = dict(particle_id=np.array([41, 17]), contact_particle_id=np.array([17]),
             displacement_m=np.zeros((3, 2, 3)), contact_force_N=np.zeros((3, 2, 3)),
             tangent_force_global_yz_N=np.array([[[2., 0.]], [[2., 0.]], [[-6., 0.]]]),
             interval_friction_J=np.array([[0.], [1.5], [1.75]]),
             friction_dissipation_J=np.array([[0.], [1.5], [3.25]]))
    h['displacement_m'][:, 1, 1] = [0., .25, .125]
    h['contact_force_N'][:, 1, 1] = [2., 2., -6.]
    records = [dict(wall_displacement_m=w, interval_wall_work_J=ww,
                    interval_contact_material_work_J=wm)
               for w, ww, wm in [([0., 0., 0.], 0., 0.),
                                 ([0., 1., 0.], 2., .5), ([0., 0., 0.], 2., .25)]]
    return h, records


def negative_work_fixture(amount):
    h = dict(particle_id=np.array([17]), contact_particle_id=np.array([17]),
             displacement_m=np.array([[[0., 0., 0.]], [[0., amount, 0.]]]),
             contact_force_N=np.array([[[0., 1., 0.]], [[0., 1., 0.]]]),
             tangent_force_global_yz_N=np.array([[[1., 0.]], [[1., 0.]]]),
             interval_friction_J=np.array([[0.], [-amount]]),
             friction_dissipation_J=np.array([[0.], [-amount]]))
    records = [dict(wall_displacement_m=[0., 0., 0.], interval_wall_work_J=0.,
                    interval_contact_material_work_J=0.),
               dict(wall_displacement_m=[0., 0., 0.], interval_wall_work_J=0.,
                    interval_contact_material_work_J=amount)]
    return h, records


class WorkReferences(unittest.TestCase):
    def test_signed_midpoint_work_and_persistent_id_lookup(self):
        result = subject.work_audit(*work_fixture())
        self.assertEqual(result['reconstructed_Q_J'], 3.25)
        self.assertEqual(result['reconstructed_Wwall_J'], 4.)
        self.assertEqual(result['reconstructed_Wmaterial_J'], .75)
        self.assertEqual(result['maximum_interval_contact_identity_J'], 0.)

    def test_raw_negative_roundoff_is_retained_not_absolute_or_clipped(self):
        result = subject.work_audit(*negative_work_fixture(1e-18))
        self.assertEqual(result['raw_minimum_Q_J'], -1e-18)
        self.assertEqual(result['raw_negative_components'], 1)
        self.assertEqual(result['reconstructed_Q_J'], -1e-18)

    def test_meaningful_negative_work_rejected(self):
        with self.assertRaisesRegex(ValueError, 'negative signed Q'):
            subject.work_audit(*negative_work_fixture(1e-4))

    def test_saved_interval_and_cumulative_tampering_rejected(self):
        for key in ('interval_friction_J', 'friction_dissipation_J'):
            h, records = work_fixture()
            h[key][-1, 0] += .01
            with self.subTest(key=key), self.assertRaises(ValueError):
                subject.work_audit(h, records)

    def test_contact_identity_and_saved_work_tampering_rejected(self):
        h, records = work_fixture()
        h['contact_force_N'][:, 1, 0] = 3.
        h['displacement_m'][1:, 1, 0] = .1
        with self.assertRaisesRegex(ValueError, 'contact work identity'):
            subject.work_audit(h, records)
        for key in ('interval_wall_work_J', 'interval_contact_material_work_J'):
            h, records = work_fixture()
            records[1][key] += .01
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'saved work'):
                subject.work_audit(h, records)


def minimal_comparison_pair():
    a = dict(cfg=dict(increments_per_phase=4, E_Pa=70e9), sources={'frozen.py': 'digest'},
             h={name: np.array([1, 2]) for name in subject.STATIC})
    b = copy.deepcopy(a)
    b['cfg']['increments_per_phase'] = 8
    return a, b


class GuardReferences(unittest.TestCase):
    def test_physical_parameter_mutation_rejected_before_comparison(self):
        a, b = minimal_comparison_pair()
        b['cfg']['E_Pa'] = 69e9
        with self.assertRaisesRegex(ValueError, 'physical/configuration mismatch'):
            subject.compare_loaded(a, b)

    def test_source_mutation_rejected_before_comparison(self):
        a, b = minimal_comparison_pair()
        b['sources']['frozen.py'] = 'changed'
        with self.assertRaisesRegex(ValueError, 'source inventory mismatch'):
            subject.compare_loaded(a, b)

    def test_each_static_geometry_id_or_boundary_mutation_rejected(self):
        for name in subject.STATIC:
            a, b = minimal_comparison_pair()
            b['h'][name][1] = 9
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'static geometry/ID/boundary mismatch'):
                subject.compare_loaded(a, b)

    def test_static_dtype_mutation_is_not_silently_equal(self):
        a, b = minimal_comparison_pair()
        b['h']['particle_id'] = b['h']['particle_id'].astype(float)
        with self.assertRaisesRegex(ValueError, 'static geometry/ID/boundary mismatch'):
            subject.compare_loaded(a, b)

    def test_partial_output_rejected_without_nearest_endpoint_comparison(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            fixtures = {'report.json': dict(status='BLOCKED', error='negative N'),
                        'resolved_config.json': {}, 'readback.json': {}}
            for name, data in fixtures.items():
                (path / name).write_text(json.dumps(data), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'BLOCKED_NOT_COMPARABLE: incomplete path'):
                subject.load_case(path, path / 'readback.json', 8)


if __name__ == '__main__':
    unittest.main(verbosity=2)

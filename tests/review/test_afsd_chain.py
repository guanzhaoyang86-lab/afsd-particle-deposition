"""Read-only chain evidence and independent metadata-corruption regressions.

Tiny arrays below are explicit I/O fixtures, not simulated or reference physics.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from inspect_afsd_chain import (exact_fields, id_keyed, inspect_chain,
                                join_segments, ordered_frames, validate_metadata,
                                validate_parent_receipt, accepted_target_complete)


def fixture_data():
    return dict(title='GEOMETRY_FIXTURE; NOT_PHYSICS',
                points=[[1., 2., 3.], [4., 5., 6.]],
                scalars={'particle_id': [9, 2], 'temperature_K': [701., 702.]},
                vectors={'velocity_m_per_s': [[.1, .2, .3], [.4, .5, .6]]})


def frame(step, time):
    return dict(step=step, time_s=time, file=f'fixture_{step}.vtk',
                tool_translation_m=[0., 0., 0.], tool_rotation_rad=0.,
                tool_axis=[0, 0, 1], rate_evaluation_time_s=time-.01,
                unavailable_fields_reason=None)


def fixture_segments():
    a, b = Path('chain_fixture_parent').resolve(), Path('chain_fixture_child').resolve()
    values = id_keyed(fixture_data())
    common = dict(sources={'fixture': 'not-a-real-source-digest'}, config={'fixture': True},
                  state={'contact_variant': 'fixture'}, first=values, last=values)
    return [{**deepcopy(common), 'directory': a, 'frames': [frame(8, .8), frame(9, .9)],
             'manifest': {'provenance': {'checkpoint_parent': None}}},
            {**deepcopy(common), 'directory': b, 'frames': [frame(9, .9), frame(10, 1.)],
             'manifest': {'provenance': {'checkpoint_parent': str(a/'checkpoint.npz')}}}]


class ChainBoundaryChecks(unittest.TestCase):
    def test_accepted_time_cannot_reclassify_failed_or_incomplete_run(self):
        self.assertTrue(accepted_target_complete(2e-6, 2e-6, 'COMPLETE', 'COMPLETE'))
        self.assertFalse(accepted_target_complete(1e-6, 2e-6, 'COMPLETE', 'COMPLETE'))
        for runner in ('COMPLETE', 'FAILED', 'INCOMPLETE'):
            for manifest in ('COMPLETE', 'FAILED', 'INCOMPLETE'):
                if runner == manifest == 'COMPLETE':
                    continue
                with self.subTest(runner=runner, manifest=manifest):
                    self.assertFalse(accepted_target_complete(2e-6, 2e-6, runner, manifest))

    def test_reordered_ids_compare_all_fields_exactly(self):
        original = fixture_data()
        sorted_data = deepcopy(original)
        sorted_data['points'] = original['points'][::-1]
        for group in ('scalars', 'vectors'):
            for name, values in original[group].items():
                sorted_data[group][name] = values[::-1]
        exact_fields(id_keyed(original), id_keyed(sorted_data))
        sorted_data['scalars']['temperature_K'][0] += .001
        with self.assertRaisesRegex(ValueError, 'temperature_K'):
            exact_fields(id_keyed(original), id_keyed(sorted_data))

    def test_missing_field_not_silently_ignored(self):
        a = id_keyed(fixture_data()); b = deepcopy(a)
        del b['velocity_m_per_s']
        with self.assertRaisesRegex(ValueError, 'inventory'):
            exact_fields(a, b)

    def test_nan_and_duplicate_or_fractional_ids_fail(self):
        for ids in ([1, 1], [1, 1.5], [1, float('nan')]):
            data = fixture_data(); data['scalars']['particle_id'] = ids
            with self.assertRaises(ValueError):
                id_keyed(data)
        data = fixture_data(); data['points'][0][2] = float('nan')
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            id_keyed(data)

    def test_nonzero_restart_and_single_deduplicated_boundary(self):
        result = join_segments(fixture_segments())
        self.assertEqual([f['step'] for f in result], [8, 9, 10])
        self.assertEqual([f['time_s'] for f in result], [.8, .9, 1.])

    def test_duplicate_internal_frame_is_error_not_deduplication(self):
        frames = [frame(8, .8), frame(8, .9)]
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            ordered_frames(frames)

    def test_rejected_is_never_accepted(self):
        frames = [frame(1, .1)]; frames[0]['file'] = 'REJECTED_candidate.vtk'
        with self.assertRaisesRegex(ValueError, 'REJECTED'):
            ordered_frames(frames)
        data = fixture_data(); data['title'] = 'REJECTED'
        with self.assertRaisesRegex(ValueError, 'REJECTED'):
            id_keyed(data)

    def test_broken_parent_path_fails(self):
        segments = fixture_segments()
        segments[1]['manifest']['provenance']['checkpoint_parent'] = 'wrong/checkpoint.npz'
        with self.assertRaisesRegex(ValueError, 'lineage'):
            join_segments(segments)

    def test_changed_source_config_and_tool_pose_fail(self):
        for key in ('sources', 'config'):
            segments = fixture_segments(); segments[1][key] = {'changed': True}
            with self.assertRaisesRegex(ValueError, 'source/config'):
                join_segments(segments)
        segments = fixture_segments(); segments[1]['frames'][0]['tool_rotation_rad'] = .5
        with self.assertRaisesRegex(ValueError, 'boundary metadata'):
            join_segments(segments)

    def test_boundary_time_and_unsampled_field_change_fail(self):
        segments = fixture_segments(); segments[1]['frames'][0]['time_s'] = .900001
        with self.assertRaisesRegex(ValueError, 'step/time'):
            join_segments(segments)
        segments = fixture_segments(); segments[1]['first']['velocity_m_per_s'][0, 2] += .01
        with self.assertRaisesRegex(ValueError, 'velocity_m_per_s'):
            join_segments(segments)

    def test_source_and_configuration_checkpoint_digests_independent(self):
        # File bytes and JSON fingerprint are computed in the test, not by the
        # checked metadata function.  This is evidence validation, not mechanics.
        with tempfile.TemporaryDirectory(prefix='afsd_chain_fixture_') as directory:
            path = Path(directory)/'fixture.txt'
            payload = b'explicit metadata fixture, not physics\n'
            path.write_bytes(payload)
            sources = {'fixture.txt': hashlib.sha256(payload).hexdigest()}
            config = {'parameters': {'fixture_only': True}}
            expected = hashlib.sha256(json.dumps(config['parameters'], sort_keys=True).encode()).hexdigest()
            state = {'source_hashes': sources, 'configuration_fingerprint': expected}
            self.assertEqual(validate_metadata(sources, config, state, Path(directory)), expected)
            wrong = deepcopy(state); wrong['configuration_fingerprint'] = 'bad'
            with self.assertRaisesRegex(ValueError, 'configuration fingerprint'):
                validate_metadata(sources, config, wrong, Path(directory))
            wrong = deepcopy(state); wrong['source_hashes'] = {}
            with self.assertRaisesRegex(ValueError, 'source metadata'):
                validate_metadata(sources, config, wrong, Path(directory))
            sources = {'fixture.txt': 'bad'}
            with self.assertRaisesRegex(ValueError, 'source digest'):
                validate_metadata(sources, config, state, Path(directory))


class FrozenAcceptedEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = inspect_chain(
            [ROOT/'outputs/short_window_e3_20261001/dx05_cfl01_05us'], 2e-6)

    def test_partial_target_stays_two_microseconds(self):
        r = self.report
        self.assertEqual(r['readback_status'], 'PASS')
        self.assertEqual(r['status'], 'VALID_PARTIAL_CHAIN')
        self.assertEqual(r['target_end_time_s'], 2e-6)
        self.assertFalse(r['reached_target'])
        self.assertEqual(r['accepted_end_step'], 48)
        self.assertEqual(r['accepted_end_time_s'], 5e-7)
        self.assertEqual(r['accepted_unique_frames'], 49)

    def test_origin_metrics_not_confused_with_future_segment_displacement(self):
        r = self.report
        self.assertTrue(r['origin_frame_available'])
        self.assertEqual(r['observed_origin_to_end'], r['segments'][0]['observed_segment'])
        self.assertGreater(r['observed_origin_to_end']['free_displacement_max_m'], 0.)
        self.assertEqual(r['segments'][0]['accepted_duration_s'], 5e-7)

    def test_saved_success_is_not_deposition_or_validation(self):
        r = self.report
        for key in ('deposition_validated', 'mixing_validated', 'paper_reproduction', 'production_authorized'):
            self.assertIs(r[key], False)
        self.assertEqual(r['observed_origin_to_end']['maximum_equivalent_plastic_strain'], 0.)
        self.assertEqual(r['segments'][0]['rejected_candidates_excluded'], [])

    def test_prerun_receipt_preserves_frozen_parent_and_target(self):
        receipt_path = ROOT/'logs/continuation_e3_20261001/preflight_receipt.json'
        parent = dict(directory=(ROOT/'outputs/short_window_e3_20261001/dx05_cfl01_05us').resolve(),
                      frames=[dict(step=48, time_s=5e-7)])
        report = validate_parent_receipt(receipt_path, parent, 2e-6)
        self.assertEqual(report['status'], 'PASS')
        self.assertEqual(report['parent_checkpoint_sha256'],
                         '18818a0831d7c2ad005b921362048121f00a29c934480038033687f1dd6a4b69')
        with self.assertRaisesRegex(ValueError, 'target changed'):
            validate_parent_receipt(receipt_path, parent, 5e-7)
        with tempfile.TemporaryDirectory(prefix='afsd_receipt_fixture_') as directory:
            path = Path(directory)/'tampered_receipt.json'
            content = json.loads(receipt_path.read_text())
            content['checkpoint_sha256'] = 'corrupt'
            path.write_text(json.dumps(content), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'receipt digest'):
                validate_parent_receipt(path, parent, 2e-6)


if __name__ == '__main__':
    unittest.main(verbosity=2)

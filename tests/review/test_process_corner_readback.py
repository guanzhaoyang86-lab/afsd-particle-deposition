"""Q4 saved-data regressions; corruption occurs only in temporary test fixtures."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from inspect_process_corner import inspect
from afsd_scaffold.vtk_legacy import read_point_cloud, write_point_cloud

BASE = ROOT/'outputs/process_corner_q4_20261001'


class CornerReadbackQ4(unittest.TestCase):
    def fixture(self, folder):
        target = Path(folder)/'Q4_CORRUPTION_FIXTURE_NOT_PHYSICAL'
        shutil.copytree(BASE/'j2_staggered_i8', target)
        return target

    def update_hash(self, target, name):
        path = target/'report.json'; report = json.loads(path.read_text())
        report['artifact_sha256'][name] = hashlib.sha256((target/name).read_bytes()).hexdigest()
        path.write_text(json.dumps(report), encoding='utf-8')

    def rewrite_history(self, target, change):
        path = target/'history.npz'
        with np.load(path, allow_pickle=False) as f: data = {key:f[key].copy() for key in f.files}
        change(data)
        np.savez_compressed(path, **data)
        self.update_hash(target, 'history.npz')

    def test_01_actual_four_cases_readback_and_independent_continuum_references(self):
        cases = [('elastic', 'staggered', 8), ('j2', 'symmetric', 8),
                 ('j2', 'staggered', 8), ('j2', 'staggered', 16)]
        for law, schedule, n in cases:
            with self.subTest(law=law, schedule=schedule, increments=n):
                result = inspect(BASE/f'{law}_{schedule}_i{n}')
                self.assertEqual(result['status'], 'PASS_COMPLETE_READBACK')
                self.assertFalse(result['production_authorized'])
                self.assertEqual(result['accepted_time_s'], .2)
                self.assertLess(max(result['maximum_interval_gapwork_by_wall_J']), 1e-6)
                self.assertLess(result['maximum_particle_plastic_defect_J'], 1e-10)
                if law == 'elastic':
                    self.assertEqual(result['final']['volume_mean_ep'], 0.)
                    np.testing.assert_allclose(result['final_surface_displacement_m'], 0., atol=1e-12, rtol=0.)
                    np.testing.assert_allclose(result['independent_wall_work_J'], 0., atol=1e-7, rtol=0.)
                else:
                    # Literal independent biaxial J2 references, not computed by the solver.
                    self.assertAlmostEqual(result['final']['volume_mean_ep'], .005428571428571429, delta=1e-10)
                    self.assertAlmostEqual(result['final']['plastic_work_J'], .4885714285714286, delta=1e-7)
                    np.testing.assert_allclose(result['final_surface_displacement_m'],
                        [-27.142857142857144e-6, -27.142857142857144e-6, 54.28571428571429e-6],
                        atol=2e-10, rtol=0.)
                    reference_work = .24214285714285715 if n == 8 else .24401785714285715
                    np.testing.assert_allclose(result['independent_wall_work_J'], [reference_work]*2,
                                               atol=1e-7, rtol=0.)
                    expected_times = [.13214285714285715] if schedule == 'symmetric' else [
                        .13214285714285715, .16607142857142858]
                    self.assertEqual(len(result['events']), len(expected_times))
                    for event, reference in zip(result['events'], expected_times):
                        self.assertAlmostEqual(event['physical_time_s'], reference, delta=1e-8)

    def test_02_artifact_hash_and_rehashed_id_geometry_corruption_rejected(self):
        for kind, message in [('hash', 'artifact changed'), ('id', 'identity mismatch'),
                              ('geometry', 'geometry/support mismatch')]:
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(prefix='q4_readback_') as folder:
                target = self.fixture(folder)
                if kind == 'hash':
                    path = target/'final.vtk'
                    path.write_bytes(path.read_bytes()+b'\nCORRUPTED_TEST_ONLY\n')
                else:
                    def change(data):
                        if kind == 'id': data['particle_id'][1] = data['particle_id'][0]
                        else: data['position_m'][-1, 13, 2] += 1e-6
                    self.rewrite_history(target, change)
                with self.assertRaisesRegex(ValueError, message): inspect(target)

    def test_03_opposite_per_wall_errors_cannot_cancel_in_report(self):
        for field in ('interval_wall_work_J', 'wall_work_J', 'contact_material_work_J'):
            with self.subTest(field=field), tempfile.TemporaryDirectory(prefix='q4_readback_') as folder:
                target = self.fixture(folder); path = target/'report.json'
                report = json.loads(path.read_text()); row = report['records'][1]
                previous_sum = sum(row[field])
                row[field][0] += .01; row[field][1] -= .01
                self.assertAlmostEqual(sum(row[field]), previous_sum, delta=1e-15)
                path.write_text(json.dumps(report), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'per-wall'): inspect(target)

    def test_04_release_metadata_must_be_correct_and_complete(self):
        for kind in ('wrong_feature', 'missing_events', 'duplicate_event', 'reversed_events'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(prefix='q4_readback_') as folder:
                target = self.fixture(folder); path = target/'report.json'
                report = json.loads(path.read_text())
                if kind == 'wrong_feature': report['events'][0]['near_zero_feature_ids'] = [2]
                elif kind == 'missing_events': report['events'] = []
                elif kind == 'duplicate_event': report['events'].append(report['events'][0])
                else: report['events'] = list(reversed(report['events']))
                path.write_text(json.dumps(report), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'event|release'): inspect(target)

    def test_05_rehashed_unsolved_temperature_and_velocity_cannot_be_published(self):
        for name in ('temperature_K', 'velocity_m_s'):
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix='q4_readback_') as folder:
                target = self.fixture(folder); path = target/'final.vtk'
                vtk = read_point_cloud(path)
                if name.startswith('temperature'):
                    vtk['scalars'][name] = [293.15]*len(vtk['points'])
                else:
                    vtk['vectors'][name] = [[0., 0., 0.]]*len(vtk['points'])
                write_point_cloud(path, vtk['points'], vtk['scalars'], vtk['vectors'],
                                  title='CORRUPTION_FIXTURE_NOT_PHYSICAL')
                self.update_hash(target, 'final.vtk')
                with self.assertRaisesRegex(ValueError, 'unsolved field published'): inspect(target)

    def test_06_blocked_status_cannot_be_promoted_or_authorize_production(self):
        with tempfile.TemporaryDirectory(prefix='q4_readback_') as folder:
            target = self.fixture(folder); path = target/'report.json'
            report = json.loads(path.read_text()); report['status'] = 'BLOCKED'
            path.write_text(json.dumps(report), encoding='utf-8')
            result = inspect(target)
            self.assertEqual(result['status'], 'PASS_PARTIAL_READBACK')
            self.assertEqual(result['completion_status'], 'BLOCKED')
            self.assertFalse(result['production_authorized'])
            report['production_authorized'] = True
            path.write_text(json.dumps(report), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'cannot authorize production'): inspect(target)


if __name__ == '__main__': unittest.main(verbosity=2)

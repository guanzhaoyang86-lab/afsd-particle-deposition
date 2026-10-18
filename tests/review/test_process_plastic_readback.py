"""Q2 evidence readback: deliberately damaged temporary copies, never solver input."""
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
from inspect_process_plastic import inspect
from afsd_scaffold.vtk_legacy import read_point_cloud, write_point_cloud

BASE = ROOT/'outputs/process_plastic_q2_20261001'


def scalar_trapezoid_reference(n):
    """Independent piecewise uniaxial continuum law; no tested functions called."""
    strains = [-.004*i/n for i in range(n+1)]
    strains += [-.004+.004*i/n for i in range(1, n+1)]
    stresses = []
    ey = 90e6/70e9
    for i, e in enumerate(strains):
        if i <= n:
            sigma = 70e9*e if e >= -ey else -90e6
        else:
            sigma = -90e6+70e9*(e+.004) if e+.004 <= 2*ey else 90e6
        stresses.append(sigma)
    return sum(.5*(stresses[i-1]+stresses[i])*(strains[i]-strains[i-1])*1e-6
               for i in range(1, len(strains)))


class PlasticReadbackQ2(unittest.TestCase):
    def fixture(self, folder):
        target = Path(folder)/'Q2_CORRUPTION_FIXTURE_NOT_PHYSICAL'
        shutil.copytree(BASE/'j2_n3_i8', target)
        return target

    def refresh_artifact(self, target, name, change_report=None):
        path = target/'report.json'; report = json.loads(path.read_text())
        report['artifact_sha256'][name] = hashlib.sha256((target/name).read_bytes()).hexdigest()
        if change_report is not None:
            change_report(report)
        path.write_text(json.dumps(report), encoding='utf-8')

    def test_actual_accepted_cases_read_back_and_independent_trapezoid_oracle(self):
        for n, expected_work in ((8, .4275), (16, .430)):
            result = inspect(BASE/f'j2_n3_i{n}')
            self.assertEqual(result['status'], 'PASS_COMPLETE_READBACK')
            self.assertEqual(result['accepted_end_s'], .2)
            self.assertEqual(result['accepted_steps'], 2*n)
            self.assertFalse(result['production_authorized'])
            self.assertGreater(result['final']['plastic_work_J'], 0.)
            self.assertAlmostEqual(scalar_trapezoid_reference(n), expected_work, delta=2e-15)
            # The continuum oracle itself is tested here. Coarse-space deviation
            # is reported by the review, not converted into a new precision gate.

    def test_artifact_tamper_is_rejected_before_any_physical_readback(self):
        with tempfile.TemporaryDirectory(prefix='q2_readback_') as folder:
            target = self.fixture(folder); path = target/'final.vtk'
            path.write_bytes(path.read_bytes()+b'\nCORRUPTED_TEST_COPY\n')
            with self.assertRaisesRegex(ValueError, 'artifact changed'):
                inspect(target)

    def test_rehashed_and_ledger_consistent_plastic_work_still_needs_independent_reference(self):
        with tempfile.TemporaryDirectory(prefix='q2_readback_') as folder:
            target = self.fixture(folder); path = target/'history.npz'
            with np.load(path, allow_pickle=False) as f:
                data = {key:f[key].copy() for key in f.files}
            data['plastic_work_J'] *= 2.
            np.savez_compressed(path, **data)
            def change(report):
                for row, values in zip(report['records'], data['plastic_work_J']):
                    row['plastic_work_J'] = float(values.sum())
            self.refresh_artifact(target, 'history.npz', change)
            with self.assertRaisesRegex(ValueError, 'independent material reference'):
                inspect(target)

    def test_blocked_status_at_target_is_never_promoted_to_complete(self):
        with tempfile.TemporaryDirectory(prefix='q2_readback_') as folder:
            target = self.fixture(folder); path = target/'report.json'
            report = json.loads(path.read_text()); report['status'] = 'BLOCKED'
            path.write_text(json.dumps(report), encoding='utf-8')
            result = inspect(target)
            self.assertEqual(result['status'], 'PASS_PARTIAL_READBACK')
            self.assertEqual(result['completion_status'], 'BLOCKED')
            self.assertFalse(result['production_authorized'])

    def test_truncated_ledger_cannot_skip_material_and_work_checks(self):
        with tempfile.TemporaryDirectory(prefix='q2_readback_') as folder:
            target = self.fixture(folder); path = target/'report.json'
            report = json.loads(path.read_text()); report['records'] = report['records'][:1]
            path.write_text(json.dumps(report), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'ledger/history'):
                inspect(target)

    def test_ledger_step_and_time_must_correspond_to_saved_frame(self):
        for field, value in (('step', 900), ('physical_time_s', .199)):
            with self.subTest(field=field):
                with tempfile.TemporaryDirectory(prefix='q2_readback_') as folder:
                    target = self.fixture(folder); path = target/'report.json'
                    report = json.loads(path.read_text()); report['records'][1][field] = value
                    path.write_text(json.dumps(report), encoding='utf-8')
                    with self.assertRaisesRegex(ValueError, 'ledger/history'):
                        inspect(target)

    def test_rehashed_duplicate_particle_id_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix='q2_readback_') as folder:
            target = self.fixture(folder); path = target/'history.npz'
            with np.load(path, allow_pickle=False) as f:
                data = {key:f[key].copy() for key in f.files}
            data['particle_id'][1] = data['particle_id'][0]
            np.savez_compressed(path, **data); self.refresh_artifact(target, 'history.npz')
            with self.assertRaisesRegex(ValueError, 'persistent ID'):
                inspect(target)

    def test_rehashed_unsolved_temperature_field_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix='q2_readback_') as folder:
            target = self.fixture(folder); path = target/'final.vtk'
            data = read_point_cloud(path)
            data['scalars']['temperature_K'] = [293.15]*len(data['points'])
            write_point_cloud(path, data['points'], data['scalars'], data['vectors'],
                              title='CORRUPTED_TEST_COPY_NOT_A_PHYSICS_RESULT')
            self.refresh_artifact(target, 'final.vtk')
            with self.assertRaisesRegex(ValueError, 'unsolved thermal/dynamic field'):
                inspect(target)


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Q1 archived data integrity; temporary corrupt copies are test fixtures only."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from inspect_process_time import inspect

BASE = ROOT/'outputs/process_time_q1_20261001'


class Q1Readback(unittest.TestCase):
    def temporary_copy(self, folder):
        target = Path(folder)/'fixture'
        shutil.copytree(BASE/'n3_i4', target)
        return target

    def test_four_actual_cases_reach_same_end_without_production_claim(self):
        for name, particles, frames in (('n3_i4',27,5),('n5_i4',125,5),('n5_i8',125,9),('n7_i4',343,5)):
            report, data, summary = inspect(BASE/name)
            self.assertEqual(summary['status'], 'PASS_READBACK_ONLY')
            self.assertEqual(summary['common_end_s'], .1)
            self.assertFalse(summary['production_authorized'])
            self.assertEqual(len(data['particle_id']), particles)
            self.assertEqual(len(data['physical_time_s']), frames)
            self.assertEqual(report['continuum_accuracy']['policy'], 'REPORTED_NOT_A_PRECISION_GATE')

    def test_failed_or_partial_run_cannot_be_called_common_end(self):
        with tempfile.TemporaryDirectory(prefix='q1_readback_fixture_') as folder:
            target = self.temporary_copy(folder)
            path = target/'report.json'; original = json.loads(path.read_text())
            for key, value in (('status', 'BLOCKED'), ('accepted_time_s', .025)):
                changed = dict(original); changed[key] = value
                path.write_text(json.dumps(changed), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'incomplete'):
                    inspect(target)

    def test_modified_artifact_is_detected(self):
        with tempfile.TemporaryDirectory(prefix='q1_readback_fixture_') as folder:
            target = self.temporary_copy(folder)
            path = target/'final.vtk'
            path.write_bytes(path.read_bytes()+b'\nCORRUPTED_TEST_COPY\n')
            with self.assertRaisesRegex(ValueError, 'changed artifact'):
                inspect(target)

    def test_reported_work_is_checked_against_force_displacement_not_itself(self):
        with tempfile.TemporaryDirectory(prefix='q1_readback_fixture_') as folder:
            target = self.temporary_copy(folder)
            path = target/'report.json'; data = json.loads(path.read_text())
            data['records'][-1]['external_work_J'] *= 2
            path.write_text(json.dumps(data), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'actual force-displacement'):
                inspect(target)


if __name__ == '__main__':
    unittest.main(verbosity=2)

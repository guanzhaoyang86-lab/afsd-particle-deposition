"""Q5 saved-data adversarial tests; only temporary copies are corrupted."""
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
from inspect_process_friction import inspect
from afsd_scaffold.vtk_legacy import read_point_cloud, write_point_cloud

BASE = ROOT/'outputs/process_friction_q5_20261001'


class FrictionReadbackQ5(unittest.TestCase):
    def fixture(self, directory):
        target = Path(directory)/'Q5_CORRUPTION_FIXTURE_NOT_PHYSICAL'
        shutil.copytree(BASE/'n3_i4', target)
        return target

    def report(self, target):
        return json.loads((target/'report.json').read_text(encoding='utf-8'))

    def save_report(self, target, report):
        (target/'report.json').write_text(json.dumps(report), encoding='utf-8')

    def rehash(self, target, name):
        report = self.report(target)
        report['artifact_sha256'][name] = hashlib.sha256((target/name).read_bytes()).hexdigest()
        self.save_report(target, report)

    def change_history(self, target, change):
        with np.load(target/'history.npz', allow_pickle=False) as saved:
            data = {key: saved[key].copy() for key in saved.files}
        change(data)
        np.savez_compressed(target/'history.npz', **data)
        self.rehash(target, 'history.npz')

    def test_01_actual_both_paths_complete_and_signed_slip_occurs(self):
        for increments in (4, 8):
            with self.subTest(increments=increments):
                path = BASE/f'n3_i{increments}'
                result = inspect(path)
                self.assertEqual(result['status'], 'PASS_COMPLETE_READBACK')
                self.assertEqual(result['accepted_time_s'], .2)
                self.assertFalse(result['production_authorized'])
                self.assertLess(result['maximum_cumulative_contact_work_defect_J'], 1e-6)
                self.assertGreaterEqual(result['minimum_raw_contact_increment_J'], -1e-10)
                # This is a path coverage check, not a manufactured cube continuum reference.
                with np.load(path/'history.npz', allow_pickle=False) as saved:
                    modes = saved['interval_modes']; slip = saved['interval_relative_slip_m']
                    times = saved['physical_time_s']
                    self.assertTrue(np.any((modes==-1) & (slip < -1e-11)))
                    self.assertTrue(np.any((modes==1) & (slip > 1e-11)))
                    self.assertTrue(np.any(modes[(times>.1)&(times<.15)]==0))

    def test_02_manifest_hashes_and_rehashed_identity_modes_are_checked(self):
        for kind, message in (('artifact','artifact changed'), ('source','source changed'),
                              ('id','identity'), ('mode','shape/mode')):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(prefix='q5_readback_') as directory:
                target = self.fixture(directory)
                if kind == 'artifact':
                    path = target/'final.vtk'; path.write_bytes(path.read_bytes()+b'\nTEST_CORRUPTION\n')
                elif kind == 'source':
                    report = self.report(target)
                    report['source_sha256']['scripts/probe_process_friction.py'] = '0'*64
                    self.save_report(target, report)
                else:
                    def change(data):
                        if kind == 'id': data['particle_id'][1] = data['particle_id'][0]
                        else: data['interval_modes'][1,0] = 2
                    self.change_history(target, change)
                with self.assertRaisesRegex(ValueError, message): inspect(target)

    def test_03_opposite_raw_contact_work_errors_cannot_cancel(self):
        for field, message in (('interval_friction_J','per-contact friction work'),
                               ('friction_dissipation_J','per-contact accumulated dissipation')):
            with self.subTest(field=field), tempfile.TemporaryDirectory(prefix='q5_readback_') as directory:
                target = self.fixture(directory)
                def change(data):
                    previous = data[field][-1].sum()
                    data[field][-1,0] -= .01
                    data[field][-1,1] += .01
                    self.assertAlmostEqual(float(data[field][-1].sum()), float(previous), delta=1e-15)
                    self.assertLess(float(data[field][-1,0]), 0.)
                self.change_history(target, change)
                with self.assertRaisesRegex(ValueError, message): inspect(target)

    def test_04_event_inventory_order_and_contact_modes_are_verified(self):
        for kind in ('missing', 'duplicate', 'reversed', 'wrong_mode', 'wrong_margin'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(prefix='q5_readback_') as directory:
                target = self.fixture(directory); report = self.report(target)
                self.assertGreater(len(report['events']), 1)
                if kind == 'missing': report['events'] = []
                elif kind == 'duplicate': report['events'].append(report['events'][0])
                elif kind == 'reversed': report['events'].reverse()
                elif kind == 'wrong_mode': report['events'][0]['ending_interval_modes'][0] = 99
                else: report['events'][0]['cone_margin_N'][0] += 1.
                self.save_report(target, report)
                with self.assertRaisesRegex(ValueError, 'event'): inspect(target)

    def test_05_unsolved_temperature_and_velocity_never_become_output(self):
        for name in ('temperature_K', 'velocity_m_s'):
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix='q5_readback_') as directory:
                target = self.fixture(directory); path = target/'final.vtk'
                vtk = read_point_cloud(path)
                if name.startswith('temperature'): vtk['scalars'][name] = [293.15]*27
                else: vtk['vectors'][name] = [[0.,0.,0.]]*27
                write_point_cloud(path, vtk['points'], vtk['scalars'], vtk['vectors'],
                                  title='Q5_CORRUPTION_FIXTURE_NOT_PHYSICAL')
                self.rehash(target, 'final.vtk')
                with self.assertRaisesRegex(ValueError, 'unsolved field'): inspect(target)
        with tempfile.TemporaryDirectory(prefix='q5_readback_') as directory:
            target = self.fixture(directory)
            self.change_history(target, lambda data:data.update(temperature_K=np.full((len(data['physical_time_s']),27),293.15)))
            with self.assertRaisesRegex(ValueError, 'unsolved field'): inspect(target)

    def test_06_blocked_is_partial_and_cannot_authorize_production(self):
        with tempfile.TemporaryDirectory(prefix='q5_readback_') as directory:
            target = self.fixture(directory); report = self.report(target)
            report.update(status='BLOCKED', error='TEST_FORCED_BLOCK_NOT_PHYSICAL')
            self.save_report(target, report)
            rejected = dict(classification='REJECTED', restart_allowed=False, error=report['error'],
                last_accepted_step=report['accepted_steps'], candidate_available=False,
                candidate_time_s=None, candidate_base_step=None, attempted_time_s=.2)
            (target/'REJECTED_manifest.json').write_text(json.dumps(rejected), encoding='utf-8')
            self.rehash(target, 'REJECTED_manifest.json')
            result = inspect(target)
            self.assertEqual(result['status'], 'PASS_PARTIAL_READBACK')
            self.assertEqual(result['completion_status'], 'BLOCKED')
            report = self.report(target); report['production_authorized'] = True
            self.save_report(target, report)
            with self.assertRaisesRegex(ValueError, 'cannot authorize production'): inspect(target)

    def test_07_actual_relative_motion_and_wall_work_cannot_be_fabricated(self):
        for kind, message in (('relative','relative slip'), ('wall_work','wall interval'), ('wall_motion','wall program')):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(prefix='q5_readback_') as directory:
                target = self.fixture(directory)
                if kind == 'relative':
                    self.change_history(target, lambda data:data['interval_relative_slip_m'].__setitem__((-1,0), .001))
                else:
                    report = self.report(target)
                    if kind == 'wall_work': report['records'][-1]['interval_wall_work_J'] += .01
                    else: report['records'][-1]['wall_displacement_m'][1] += .001
                    self.save_report(target, report)
                with self.assertRaisesRegex(ValueError, message): inspect(target)


if __name__ == '__main__': unittest.main(verbosity=2)

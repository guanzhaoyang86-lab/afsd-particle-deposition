"""Administrative launcher tests: Popen and supervisor are always mocked."""
from pathlib import Path
import contextlib
import io
import json
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import start_gpu_four_hours_detached_20261009 as launcher


class DetachedLauncherTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='detached_admin_test_', dir=launcher.BASE)
        self.folder = Path(self.temp.name); self.parent = self.folder/'parent'; self.parent.mkdir()
        self.output = self.folder/'campaign'
        self.qualification = self.folder/'qualification.json'; self.qualification.write_text('{}')
        self.args = ['--output', str(self.output), '--parent', str(self.parent),
                     '--runner', str(ROOT/'scripts/run_gpu_cached_preflight_20261008.py'),
                     '--runner-qualification', str(self.qualification), '--initial-dt', '.0002',
                     '--minimum-dt', '1e-7', '--max-halvings', '6']
        self.policy = dict(pid=123, in_job=False, limit_flags=0,
                           breakaway_allowed=False, kill_on_job_close=False)

    def tearDown(self):
        self.temp.cleanup()

    def test_parent_uses_exact_hidden_breakaway_flags_and_forwards_options(self):
        with patch.object(launcher, 'job_policy', return_value=self.policy), patch.object(
                launcher.subprocess, 'Popen', return_value=Mock(pid=456)) as popen, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(launcher.main(self.args), 0)
        self.assertEqual(popen.call_count, 1)
        command = popen.call_args.args[0]; options = popen.call_args.kwargs
        self.assertIn('--detached-child', command)
        for name, value in [('--wall-seconds', '14400.0'), ('--initial-dt', '0.0002'),
                            ('--minimum-dt', '1e-07'), ('--max-halvings', '6')]:
            self.assertEqual(command[command.index(name)+1], value)
        self.assertIn('--runner-qualification', command)
        self.assertEqual(options['creationflags'], 16777736)
        self.assertEqual(options['startupinfo'].wShowWindow, 0)
        self.assertTrue(options['startupinfo'].dwFlags & launcher.subprocess.STARTF_USESHOWWINDOW)
        self.assertEqual(options['stdin'], launcher.subprocess.DEVNULL)
        self.assertEqual(options['stderr'], launcher.subprocess.STDOUT)
        self.assertTrue(options['close_fds'])
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            self.assertEqual(options['env'][name], '1')
        self.assertEqual(options['env']['PYTHONPATH'], '')
        receipt = json.loads(self.output.with_name('campaign.launch.json').read_text())
        self.assertEqual(receipt['creationflags'], 16777736)
        self.assertIn('scripts/start_gpu_resume_detached_20261008.py', receipt['source_sha256'])
        self.assertEqual(receipt['wall_limit_s'], 14400.)
        self.assertFalse(self.output.exists())

    def test_child_records_policy_before_invoking_fixed_supervisor(self):
        def inspect_before_run(path, run_name):
            receipt = json.loads(self.output.with_name('campaign.bootstrap.json').read_text())
            self.assertEqual(receipt['status'], 'PASS_DETACHED_JOB_POLICY')
            self.assertEqual(path, str(launcher.SUPERVISOR))
            self.assertEqual(run_name, '__main__')
            self.assertEqual(sys.argv[0], str(launcher.SUPERVISOR))
            self.assertEqual(sys.argv[sys.argv.index('--wall-seconds')+1], '1800.0')
        with patch.object(launcher, 'job_policy', return_value=self.policy), patch.object(
                launcher.runpy, 'run_path', side_effect=inspect_before_run) as run, patch.object(
                launcher.subprocess, 'Popen') as popen, patch.object(sys, 'argv', ['launcher_test']), patch.object(sys, 'path', list(sys.path)):
            self.assertEqual(launcher.main(self.args+['--detached-child', '--wall-seconds', '1800']), 0)
        run.assert_called_once(); popen.assert_not_called()

    def test_remaining_kill_on_close_never_runs_supervisor(self):
        policy = dict(self.policy, kill_on_job_close=True, in_job=True, limit_flags=0x2000)
        with patch.object(launcher, 'job_policy', return_value=policy), patch.object(
                launcher.runpy, 'run_path') as run, patch.object(launcher.subprocess, 'Popen') as popen:
            with self.assertRaisesRegex(RuntimeError, 'kill-on-job-close'):
                launcher.main(self.args+['--detached-child'])
        run.assert_not_called(); popen.assert_not_called()
        receipt = json.loads(self.output.with_name('campaign.bootstrap.json').read_text())
        self.assertEqual(receipt['status'], 'FAIL_KILL_ON_JOB_CLOSE_STILL_PRESENT')

    def test_breakaway_launch_failure_has_no_fallback_or_success_receipt(self):
        with patch.object(launcher, 'job_policy', return_value=self.policy), patch.object(
                launcher.subprocess, 'Popen', side_effect=PermissionError('Job breakaway denied')) as popen:
            with self.assertRaises(PermissionError):
                launcher.main(self.args)
        self.assertEqual(popen.call_count, 1)
        self.assertFalse(self.output.with_name('campaign.launch.json').exists())

    def test_excess_time_invalid_steps_and_reused_output_fail_before_launch(self):
        with patch.object(launcher.subprocess, 'Popen') as popen:
            for extra in (['--wall-seconds', '14401'], ['--wall-seconds', 'nan'],
                          ['--initial-dt', '0'], ['--max-halvings', '9']):
                with self.subTest(extra=extra), self.assertRaises(ValueError):
                    launcher.main(self.args+extra)
            self.output.mkdir()
            with self.assertRaises(ValueError):
                launcher.main(self.args)
        popen.assert_not_called()

    def test_external_path_is_rejected_before_launch(self):
        with patch.object(launcher.subprocess, 'Popen') as popen:
            with self.assertRaisesRegex(ValueError, 'project-local'):
                launcher.main(self.args+['--output', str(ROOT.parent/'external_campaign')])
        popen.assert_not_called()

    def test_existing_launch_record_is_not_overwritten(self):
        log = self.output.with_name('campaign.stdout.log'); log.write_bytes(b'old evidence')
        with patch.object(launcher.subprocess, 'Popen') as popen:
            with self.assertRaisesRegex(ValueError, 'new detached launch records'):
                launcher.main(self.args)
        popen.assert_not_called(); self.assertEqual(log.read_bytes(), b'old evidence')


if __name__ == '__main__':
    unittest.main(verbosity=2)

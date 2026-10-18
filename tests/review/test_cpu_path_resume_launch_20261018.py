"""Launcher-only tests. Every Popen and controller execution is intercepted."""
from contextlib import ExitStack, redirect_stderr, redirect_stdout
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'outputs/short_bar_event_continue_20261018_v5/launch.py'
spec = importlib.util.spec_from_file_location('test_local_interval_launcher', SOURCE)
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)


class LocalIntervalLaunch(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.directory = Path(self.stack.enter_context(tempfile.TemporaryDirectory(
            prefix='launch_unit_only_', dir=ROOT / 'outputs')))
        self.campaign = self.directory / 'controller.py'
        self.campaign.write_text('# Not a physical controller; mock-only test fixture.\n', encoding='utf-8')
        for field, path in {
                'BASE': self.directory, 'CAMPAIGN': self.campaign,
                'CONFIG': self.directory / 'run.configuration.json',
                'BOOTSTRAP': self.directory / 'run.bootstrap.json',
                'RECEIPT': self.directory / 'run.launch.json',
                'LOG': self.directory / 'run.stdout.log', 'EXIT': self.directory / 'run.exit.json'}.items():
            self.stack.enter_context(patch.object(launcher, field, path))
        self.started = datetime(2026, 10, 20, 10, 0, 0, tzinfo=timezone.utc)
        self.clock = self.stack.enter_context(patch.object(launcher, 'now', return_value=self.started))
        self.policy = self.stack.enter_context(patch.object(launcher, 'job_policy',
            return_value={'in_job': False, 'kill_on_job_close': False}))
        # This guard is installed before every test, including rejection cases.
        self.spawn = self.stack.enter_context(patch.object(launcher.subprocess, 'Popen',
            side_effect=AssertionError('test must explicitly configure a mocked spawn')))
        self.controller = self.stack.enter_context(patch.object(launcher.runpy, 'run_path',
            side_effect=AssertionError('test must explicitly configure mocked controller execution')))
        self.output = self.stack.enter_context(redirect_stdout(io.StringIO()))
        self.errors = self.stack.enter_context(redirect_stderr(io.StringIO()))

    def empty_records(self):
        self.assertEqual({path.name for path in self.directory.iterdir()}, {'controller.py'})
        self.spawn.assert_not_called()
        self.controller.assert_not_called()

    def configure_mock_spawn(self):
        self.spawn.side_effect = None
        self.spawn.return_value = Mock(pid=123456)

    def child_configuration(self, seconds=18000.):
        config = dict(schema=launcher.CONFIGURATION_SCHEMA,
            started_utc=self.started.isoformat(),
            deadline_utc=(self.started+timedelta(seconds=seconds)).isoformat(), wall_seconds=seconds)
        launcher.write_new(launcher.CONFIG, config)
        return config, launcher.digest(launcher.CONFIG)

    def child_args(self, configuration_sha):
        return ['--detached-child', '--expected-campaign-sha', launcher.digest(self.campaign),
                '--expected-configuration-sha', configuration_sha]

    def test_missing_duration_rejects_without_files_or_processes(self):
        with self.assertRaises(SystemExit) as result:
            launcher.main([])
        self.assertEqual(result.exception.code, 2)
        self.empty_records()

    def test_inspect_is_read_only_before_duration_validation(self):
        self.assertEqual(launcher.main(['--inspect-job', '--wall-seconds', 'not-a-number']), 0)
        self.empty_records()
        self.assertIn('kill_on_job_close', self.output.getvalue())
        self.clock.assert_not_called()

    def test_inspect_does_not_require_duration_or_create_configuration(self):
        with patch.object(launcher, 'write_new', side_effect=AssertionError('inspection cannot write')):
            self.assertEqual(launcher.main(['--inspect-job']), 0)
        self.empty_records()

    def test_invalid_durations_reject_without_side_effects(self):
        for value in ('nan', 'inf', '-inf', '-1', '0', '649.999', 'true', 'bad'):
            with self.subTest(value=value), self.assertRaises(SystemExit):
                launcher.main(['--wall-seconds=' + value])
            self.empty_records()
        self.policy.assert_not_called()

    def test_bool_duration_rejected_by_validation(self):
        with self.assertRaises(ValueError):
            launcher.duration(True)

    def test_parent_freezes_start_plus_duration_and_hidden_child_inherits_hash(self):
        self.configure_mock_spawn()
        self.assertEqual(launcher.main(['--wall-seconds', '18000']), 0)
        config = json.loads(launcher.CONFIG.read_text())
        self.assertEqual(config, dict(schema=launcher.CONFIGURATION_SCHEMA,
            started_utc=self.started.isoformat(),
            deadline_utc=(self.started+timedelta(hours=5)).isoformat(), wall_seconds=18000.))
        args, kwargs = self.spawn.call_args
        self.assertIn('--expected-configuration-sha', args[0])
        self.assertIn(launcher.digest(launcher.CONFIG), args[0])
        self.assertNotIn('--wall-seconds', args[0])
        self.assertEqual(kwargs['startupinfo'].wShowWindow, 0)
        self.assertTrue(kwargs['creationflags'] & launcher.subprocess.CREATE_BREAKAWAY_FROM_JOB)
        self.assertTrue(kwargs['creationflags'] & launcher.subprocess.DETACHED_PROCESS)
        self.assertEqual(kwargs['env']['OMP_NUM_THREADS'], '1')
        receipt = json.loads(launcher.RECEIPT.read_text())
        self.assertEqual(receipt['authorized_deadline_utc'], config['deadline_utc'])
        self.assertEqual(receipt['source_sha256'][launcher.CONFIG.relative_to(ROOT).as_posix()],
                         launcher.digest(launcher.CONFIG))
        self.controller.assert_not_called()

    def test_no_arbitrary_twelve_hour_limit(self):
        self.configure_mock_spawn()
        self.assertEqual(launcher.main(['--wall-seconds', '172800']), 0)
        config = json.loads(launcher.CONFIG.read_text())
        self.assertEqual(config['wall_seconds'], 172800.)
        self.assertEqual(config['deadline_utc'], (self.started+timedelta(days=2)).isoformat())

    def test_unrepresentable_duration_rejects_without_configuration_or_spawn(self):
        with self.assertRaises(ValueError):
            launcher.main(['--wall-seconds', '1e308'])
        self.empty_records()

    def test_child_uses_original_deadline_without_restarting_clock(self):
        config, expected = self.child_configuration()
        original_bytes = launcher.CONFIG.read_bytes()
        self.clock.return_value = self.started+timedelta(minutes=20)
        original_argv = sys.argv
        called_argv = []
        def mock_controller(path, run_name):
            called_argv.extend(sys.argv)
            self.assertEqual(path, str(self.campaign))
            self.assertEqual(run_name, '__main__')
            raise SystemExit(0)
        self.controller.side_effect = mock_controller
        self.assertEqual(launcher.main(self.child_args(expected)), 0)
        self.assertEqual(called_argv, [str(self.campaign), '--run', '--configuration', str(launcher.CONFIG),
                                      '--expected-configuration-sha', expected])
        self.assertIs(sys.argv, original_argv)
        self.assertEqual(launcher.CONFIG.read_bytes(), original_bytes)
        for path in (launcher.BOOTSTRAP, launcher.EXIT):
            receipt = json.loads(path.read_text())
            self.assertEqual(receipt['authorized_deadline_utc'], config['deadline_utc'])
            self.assertEqual(receipt['authorized_started_utc'], config['started_utc'])
            self.assertEqual(receipt['configuration_sha256'], expected)
        self.spawn.assert_not_called()

    def test_child_rejects_new_duration_instead_of_resetting_budget(self):
        _, expected = self.child_configuration()
        with self.assertRaises(SystemExit):
            launcher.main(self.child_args(expected)+['--wall-seconds', '18000'])
        self.spawn.assert_not_called()
        self.controller.assert_not_called()
        self.assertFalse(launcher.BOOTSTRAP.exists())

    def test_child_configuration_hash_mismatch_rejects_execution(self):
        self.child_configuration()
        with self.assertRaises(RuntimeError):
            launcher.main(self.child_args('0'*64))
        self.controller.assert_not_called()
        self.spawn.assert_not_called()
        self.assertFalse(launcher.BOOTSTRAP.exists())

    def test_child_expired_configuration_is_not_extended(self):
        _, expected = self.child_configuration()
        self.clock.return_value = self.started+timedelta(hours=6)
        with self.assertRaises(RuntimeError):
            launcher.main(self.child_args(expected))
        self.controller.assert_not_called()
        self.spawn.assert_not_called()
        self.assertFalse(launcher.BOOTSTRAP.exists())

    def test_child_job_kill_on_close_policy_rejects_controller(self):
        _, expected = self.child_configuration()
        self.policy.return_value = {'in_job': True, 'kill_on_job_close': True}
        with self.assertRaises(RuntimeError):
            launcher.main(self.child_args(expected))
        self.controller.assert_not_called()
        receipt = json.loads(launcher.BOOTSTRAP.read_text())
        self.assertEqual(receipt['status'], 'FAIL_DETACHED_JOB_POLICY_OR_SOURCE')

    def test_spawn_failure_keeps_configuration_and_records_failure(self):
        self.spawn.side_effect = OSError('synthetic spawn failure')
        with self.assertRaises(OSError):
            launcher.main(['--wall-seconds', '18000'])
        self.assertTrue(launcher.CONFIG.is_file())
        self.assertTrue(launcher.LOG.is_file())
        receipt = json.loads(launcher.RECEIPT.read_text())
        self.assertEqual(receipt['status'], 'FAIL_DETACHED_CAMPAIGN_LAUNCH')
        self.assertTrue(receipt['configuration_retained'])
        self.assertFalse(receipt['physics_started'])
        self.assertEqual(receipt['configuration_sha256'], launcher.digest(launcher.CONFIG))
        self.controller.assert_not_called()

    def test_existing_configuration_prevents_duplicate_launch_without_overwrite(self):
        self.child_configuration()
        original = launcher.CONFIG.read_bytes()
        with self.assertRaises(RuntimeError):
            launcher.main(['--wall-seconds', '18000'])
        self.assertEqual(launcher.CONFIG.read_bytes(), original)
        self.spawn.assert_not_called()

    def test_inconsistent_configuration_deadline_rejected(self):
        config, _ = self.child_configuration()
        config['deadline_utc'] = (self.started+timedelta(hours=6)).isoformat()
        with self.assertRaises(ValueError):
            launcher.validate_configuration(config)

    def test_configuration_boolean_duration_and_non_utc_timestamp_rejected(self):
        config, _ = self.child_configuration()
        damaged = dict(config, wall_seconds=True)
        with self.assertRaises(ValueError):
            launcher.validate_configuration(damaged)
        damaged = dict(config, started_utc='2026-10-20T18:00:00+08:00')
        with self.assertRaises(ValueError):
            launcher.validate_configuration(damaged)


if __name__ == '__main__':
    unittest.main()

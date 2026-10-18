"""Fixed-deadline handoff and reused controller tests; no physical execution."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import ast
import contextlib
import copy
import importlib.util
import io
import sys
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import supervise_gpu_campaign_fast_20261008 as fast
import supervise_gpu_campaign_20261008 as frozen

spec = importlib.util.spec_from_file_location('frozen_campaign_components',
    ROOT/'tests/review/test_gpu_campaign_supervisor_20261008.py')
old_tests = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old_tests)
STOP = datetime.fromisoformat(fast.AUTHORIZED_STOP_AT)


class DeadlineTests(unittest.TestCase):
    def test_original_absolute_deadline_gives_only_remaining_budget(self):
        self.assertEqual(fast.remaining_seconds(fast.AUTHORIZED_STOP_AT, STOP-timedelta(seconds=321.)), 321.)
        self.assertEqual(fast.remaining_seconds(fast.AUTHORIZED_STOP_AT, STOP+timedelta(seconds=1.)), 0.)

    def test_deadline_extension_and_naive_time_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'absolute deadline'):
            fast.remaining_seconds((STOP+timedelta(seconds=1.)).isoformat(), STOP)
        with self.assertRaisesRegex(ValueError, 'timezone-aware'):
            fast.remaining_seconds('2026-10-09T02:19:45', STOP)
        with self.assertRaisesRegex(ValueError, 'aware current'):
            fast.remaining_seconds(fast.AUTHORIZED_STOP_AT, datetime(2026, 10, 8))

    def test_wall_clock_jumps_cannot_extend_either_deadline(self):
        p = dict(latest_stop_at=fast.AUTHORIZED_STOP_AT)
        with mock.patch.object(fast.time, 'perf_counter', return_value=100.), \
             mock.patch.object(fast, 'utc_now', return_value=STOP-timedelta(seconds=20.)):
            self.assertEqual(fast.budget_remaining(p, 300.), 20.)
        with mock.patch.object(fast.time, 'perf_counter', return_value=100.), \
             mock.patch.object(fast, 'utc_now', return_value=STOP-timedelta(seconds=30000.)):
            self.assertEqual(fast.budget_remaining(p, 300.), 200.)

    def test_frozen_physics_and_controller_helpers_are_reused_by_identity(self):
        for name in ('decision', 'verify_parent', 'snapshot', 'classify_failure',
                     'singleton_and_awake', 'phase_table', 'wall_now'):
            self.assertIs(getattr(fast, name), getattr(frozen, name))

    def test_supervisor_imports_only_stdlib_and_old_stdlib_supervisor(self):
        tree = ast.parse(Path(fast.__file__).read_text(encoding='utf-8'))
        imported = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(a.name.split('.')[0] for a in node.names)
            if isinstance(node, ast.ImportFrom):
                imported.append(node.module.split('.')[0])
        self.assertLessEqual(set(imported), {'pathlib', 'datetime', 'argparse', 'ctypes', 'json',
            'math', 'os', 'subprocess', 'sys', 'time', 'traceback', 'supervise_gpu_campaign_20261008'})
        self.assertEqual(fast.RUNNER.name, 'run_gpu_cached_preflight_20261008.py')

    def test_no_wall_seconds_cli_can_reset_authorization(self):
        text = Path(fast.__file__).read_text(encoding='utf-8')
        self.assertNotIn("add_argument('--wall-seconds'", text)
        self.assertIn("add_argument('--origin-campaign', required=True)", text)


class OriginTests(old_tests.TemporaryWorkspace):
    def origin(self, finished=True):
        origin = self.folder/'original'; origin.mkdir()
        parent = origin/'attempt_0001'; parent.mkdir()
        protocol = dict(schema='GPU_TEN_HOUR_SUPERVISOR_20261008', output=fast.rel(origin),
            initial_parent=fast.rel(parent), source_sha256={}, child_source_sha256={},
            wall_limit_s=36000., started_at=(STOP-timedelta(hours=10)).isoformat(),
            latest_stop_at=fast.AUTHORIZED_STOP_AT, supervisor_pid=12345)
        progress = dict(latest_stop_at=fast.AUTHORIZED_STOP_AT, latest_pass_stage=fast.rel(parent),
            status='STOPPED_USER_REQUEST' if finished else 'RUNNING_GPU_CHILD',
            finished=finished, exit_code=0 if finished else None, child_pid=None if finished else 456,
            checkpoint='synthetic.npz', checkpoint_sha256='synthetic-hash', physical_time_s=.02)
        fast.atomic_json(origin/'protocol.json', protocol)
        fast.atomic_json(origin/'progress.json', progress)
        if finished:
            fast.atomic_json(origin/'final_summary.json', progress)
        return origin, parent, protocol, progress

    def test_clean_user_stop_is_a_valid_handoff_and_binds_final_receipts(self):
        origin, parent, _, _ = self.origin()
        with mock.patch.object(fast, 'ORIGIN', origin), \
             mock.patch.object(fast, 'process_is_running', return_value=False):
            _, _, hashes = fast.read_origin(origin, parent)
        self.assertEqual(len(hashes), 3)
        self.assertEqual(hashes[fast.rel(origin/'protocol.json')], fast.sha(origin/'protocol.json'))

    def test_running_origin_is_preflight_only_and_never_launchable(self):
        origin, parent, _, _ = self.origin(finished=False)
        with mock.patch.object(fast, 'ORIGIN', origin):
            _, progress, hashes = fast.read_origin(origin, parent, allow_running=True)
            self.assertFalse(progress['finished']); self.assertEqual(len(hashes), 1)
            with self.assertRaisesRegex(ValueError, 'finish cleanly'):
                fast.read_origin(origin, parent)

    def test_stale_parent_and_live_origin_process_are_rejected(self):
        origin, parent, _, _ = self.origin()
        stale = origin/'older_PASS'; stale.mkdir()
        with mock.patch.object(fast, 'ORIGIN', origin):
            with self.assertRaisesRegex(ValueError, 'exact latest'):
                fast.read_origin(origin, stale)
            with mock.patch.object(fast, 'process_is_running', return_value=True):
                with self.assertRaisesRegex(ValueError, 'has not exited'):
                    fast.read_origin(origin, parent)

    def test_origin_budget_and_absolute_deadline_are_not_reinterpreted(self):
        origin, parent, protocol, _ = self.origin()
        with mock.patch.object(fast, 'ORIGIN', origin):
            for limit in (36000.1, float('inf'), float('nan'), -1.):
                altered = dict(protocol, wall_limit_s=limit)
                # Deliberately malformed fixture, not an artifact writer change.
                import json
                (origin/'protocol.json').write_text(json.dumps(altered), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'at most ten hours'):
                    fast.read_origin(origin, parent, allow_running=True)
            altered = dict(protocol, latest_stop_at=(STOP+timedelta(seconds=1.)).isoformat())
            fast.atomic_json(origin/'protocol.json', altered)
            with self.assertRaisesRegex(ValueError, 'deadline or duration'):
                fast.read_origin(origin, parent, allow_running=True)

    def test_failed_stop_and_mismatched_final_cannot_be_handoff(self):
        origin, parent, _, progress = self.origin()
        with mock.patch.object(fast, 'ORIGIN', origin):
            altered = dict(progress, status='STOPPED_CHILD_FAILURE', exit_code=2)
            fast.atomic_json(origin/'progress.json', altered)
            with self.assertRaisesRegex(ValueError, 'finish cleanly'):
                fast.read_origin(origin, parent)
            fast.atomic_json(origin/'progress.json', progress)
            fast.atomic_json(origin/'final_summary.json', dict(progress, physical_time_s=1.))
            with self.assertRaisesRegex(ValueError, 'final receipt'):
                fast.read_origin(origin, parent)


class ReusedChildControllerTests(old_tests.ChildControllerTests):
    """Run the frozen five controller scenarios against the new copied loop."""
    def setUp(self):
        super().setUp()
        actual_supervise = fast.supervise
        def with_bound_deadline(protocol, meta):
            protocol['latest_stop_at'] = fast.AUTHORIZED_STOP_AT
            now = STOP-timedelta(seconds=protocol['wall_limit_s'])
            with mock.patch.object(fast, 'validate_handoff'), mock.patch.object(fast, 'utc_now', return_value=now):
                return actual_supervise(protocol, meta)
        self.enterContext(mock.patch.object(old_tests, 'campaign', fast))
        self.enterContext(mock.patch.object(fast, 'supervise', side_effect=with_bound_deadline))


class NewLoopTests(old_tests.TemporaryWorkspace):
    def protocol(self):
        parent = self.folder/'parent'; parent.mkdir()
        meta = old_tests.metadata()
        meta['observer'] = dict(supplied_material=dict(
            eligible_crossed=dict(particle_count=0, mass_kg=0.),
            crossed_under_annular_shoulder=dict(particle_count=0, mass_kg=0.)))
        p = dict(output=fast.rel(self.folder/'campaign'), initial_parent=fast.rel(parent),
            initial_snapshot=dict(latest_pass_stage=fast.rel(parent), checkpoint='unchanged.npz'),
            requested_macro_dt_s=.0004, wall_limit_s=36000., source_sha256={}, child_source_sha256={},
            minimum_dt_s=1e-7, max_halvings=8, phases=copy.deepcopy(old_tests.PHASES),
            python=sys.executable, latest_stop_at=fast.AUTHORIZED_STOP_AT)
        return p, meta

    def test_supervise_recomputes_remaining_after_planning(self):
        protocol, meta = self.protocol()
        with mock.patch.object(fast, 'validate_handoff'), \
             mock.patch.object(fast, 'utc_now', return_value=STOP+timedelta(seconds=1.)), \
             mock.patch.object(fast, 'singleton_and_awake', return_value=(mock.Mock(), 123)), \
             mock.patch.object(fast.subprocess, 'Popen') as launch, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(fast.supervise(protocol, meta), 0)
        launch.assert_not_called()
        saved = fast.read(ROOT/protocol['output']/'protocol.json')
        self.assertEqual(saved['planned_remaining_wall_s'], 36000.)
        self.assertEqual(saved['wall_limit_s'], 0.)
        self.assertEqual(saved['latest_stop_at'], fast.AUTHORIZED_STOP_AT)

    def test_stop_after_stage_is_graceful_and_only_checked_between_children(self):
        protocol, meta = self.protocol()
        actual_write = fast.atomic_json
        def write_and_request(path, value):
            actual_write(path, value)
            if Path(path).name == 'progress.json' and value['status'] == 'STARTING':
                (Path(path).parent/'STOP_AFTER_STAGE').write_text('finish current stage', encoding='utf-8')
        with mock.patch.object(fast, 'validate_handoff'), \
             mock.patch.object(fast, 'utc_now', return_value=STOP-timedelta(seconds=1000.)), \
             mock.patch.object(fast, 'singleton_and_awake', return_value=(mock.Mock(), 123)), \
             mock.patch.object(fast, 'atomic_json', side_effect=write_and_request), \
             mock.patch.object(fast.subprocess, 'Popen') as launch, \
             mock.patch.object(fast, 'terminate_owned_tree') as terminate, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(fast.supervise(protocol, meta), 0)
        launch.assert_not_called(); terminate.assert_not_called()
        self.assertEqual(fast.read(ROOT/protocol['output']/'final_summary.json')['status'], 'STOPPED_AFTER_STAGE_REQUEST')
        tree = ast.parse(Path(fast.__file__).read_text(encoding='utf-8'))
        child_loop = next(node for node in ast.walk(tree) if isinstance(node, ast.While)
            and isinstance(node.test, ast.Compare))
        self.assertNotIn('STOP_AFTER_STAGE', ast.dump(child_loop))


def load_tests(loader, standard_tests, pattern):
    # These run the exact shared helper implementation without changing frozen
    # tests. The only omitted decision case targets the retired --wall-seconds.
    for name in loader.getTestCaseNames(old_tests.DecisionTests):
        if name != 'test_ten_hour_authorization_cannot_be_extended':
            standard_tests.addTest(old_tests.DecisionTests(name))
    standard_tests.addTests(loader.loadTestsFromTestCase(old_tests.FailureClassificationTests))
    standard_tests.addTests(loader.loadTestsFromTestCase(old_tests.ParentReceiptTests))
    return standard_tests


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Administrative four-hour authorization and fake-child checks; no physics."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
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
import supervise_gpu_four_hours_20261009 as four
import supervise_gpu_campaign_20261008 as frozen

spec = importlib.util.spec_from_file_location('four_hour_frozen_components',
    ROOT/'tests/review/test_gpu_campaign_supervisor_20261008.py')
old_tests = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old_tests)
START = datetime(2026, 10, 9, 10, 12, 0, 123456, tzinfo=timezone.utc)


class AuthorizationTests(unittest.TestCase):
    def test_new_authorization_has_exact_four_hour_deadline(self):
        timing = four.fresh_deadline(14400., START, 100.)
        self.assertEqual(timing['monotonic_deadline'], 14500.)
        self.assertEqual(four.aware_time(timing['started_at']), START)
        self.assertEqual(four.aware_time(timing['latest_stop_at']), START+timedelta(hours=4))
        self.assertNotEqual(timing['latest_stop_at'], '2026-10-09T02:19:45+08:00')

    def test_nonfinite_negative_boolean_and_over_four_hour_budgets_reject(self):
        for value in (14400.000001, 36000., 0., -1., True, float('nan'), float('inf')):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'four hours'):
                four.bounded_wall_seconds(value)
        self.assertEqual(four.bounded_wall_seconds(300.), 300.)

    def test_cli_plan_rejects_excess_budget_before_reading_any_parent(self):
        args = argparse.Namespace(output=str(four.BASE/'uncreated_four_hour_test'),
            parent=str(four.BASE/'unread_parent'), wall_seconds=14400.01)
        with mock.patch.object(four, 'verify_parent') as parent:
            with self.assertRaisesRegex(ValueError, 'four hours'):
                four.plan(args)
        parent.assert_not_called()

    def test_naive_start_and_nonfinite_monotonic_values_reject(self):
        with self.assertRaisesRegex(ValueError, 'aware new'):
            four.fresh_deadline(14400., START.replace(tzinfo=None), 0.)
        with self.assertRaisesRegex(ValueError, 'finite monotonic'):
            four.fresh_deadline(14400., START, float('nan'))

    def test_wall_clock_jumps_cannot_extend_either_budget(self):
        timing = four.fresh_deadline(14400., START, 100.)
        with mock.patch.object(four.time, 'perf_counter', return_value=200.), \
             mock.patch.object(four, 'utc_now', return_value=START-timedelta(days=1)):
            self.assertEqual(four.budget_remaining(timing, timing['monotonic_deadline']), 14300.)
        with mock.patch.object(four.time, 'perf_counter', return_value=200.), \
             mock.patch.object(four, 'utc_now', return_value=START+timedelta(hours=4)):
            self.assertEqual(four.budget_remaining(timing, timing['monotonic_deadline']), 0.)

    def test_initial_timestep_may_only_stay_or_shrink(self):
        parent_dt = .000474609375
        self.assertEqual(four.initial_dt(parent_dt), parent_dt)
        self.assertEqual(four.initial_dt(parent_dt, parent_dt/2), parent_dt/2)
        for value in (parent_dt*2, 0., -1., float('nan'), float('inf'), True):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'cannot grow'):
                four.initial_dt(parent_dt, value)

    def test_frozen_decisions_parent_verification_and_process_helpers_are_reused(self):
        for name in ('decision', 'verify_parent', 'snapshot', 'classify_failure',
                'singleton_and_awake', 'phase_table', 'wall_now', 'atomic_json',
                'launch_options', 'terminate_owned_tree'):
            self.assertIs(getattr(four, name), getattr(frozen, name))

    def test_module_imports_are_only_stdlib_and_frozen_administrative_helpers(self):
        tree = ast.parse(Path(four.__file__).read_text(encoding='utf-8'))
        imports = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                imports.update(x.name.split('.')[0] for x in n.names)
            elif isinstance(n, ast.ImportFrom):
                imports.add(n.module.split('.')[0])
        self.assertLessEqual(imports, {'pathlib', 'datetime', 'argparse', 'json', 'math',
            'os', 'subprocess', 'sys', 'time', 'traceback', 'supervise_gpu_campaign_20261008'})


class BridgeBindingTests(old_tests.TemporaryWorkspace):
    def fixture(self):
        runner = self.folder/'runner.py'; runner.write_text('# synthetic bridge', encoding='utf-8')
        source = self.folder/'frozen.txt'; source.write_text('immutable parent', encoding='utf-8')
        bridge = self.folder/'bridge.txt'; bridge.write_text('qualified extension', encoding='utf-8')
        parent = {four.rel(source): four.sha(source)}
        child = dict(parent, **{four.rel(bridge): four.sha(bridge)})
        record = dict(status=four.BRIDGE_STATUS, exit_code=0, runner=four.rel(runner),
            runner_sha256=four.sha(runner), parent_source_sha256=parent,
            child_source_sha256=child, source_sha256={four.rel(runner): four.sha(runner)},
            evidence_sha256={four.rel(bridge): four.sha(bridge)})
        path = self.folder/'qualification.json'; four.atomic_json(path, record)
        return runner, path, record, dict(source_sha256=child)

    def test_qualified_bridge_binds_exact_runner_parent_and_child_sources(self):
        runner, path, record, meta = self.fixture()
        with mock.patch.object(four, 'ALLOWED_RUNNERS', (runner,)):
            child, evidence, qualification = four.runner_sources(runner, path, meta)
        self.assertEqual(child, record['child_source_sha256'])
        self.assertEqual(evidence[four.rel(path)], four.sha(path))
        self.assertEqual(qualification, dict(path=four.rel(path), sha256=four.sha(path)))

    def test_unqualified_unlisted_or_mismatched_parent_bridge_rejects(self):
        runner, path, record, meta = self.fixture()
        with self.assertRaisesRegex(ValueError, 'registered cached'):
            four.runner_sources(runner, path, meta)
        with mock.patch.object(four, 'ALLOWED_RUNNERS', (runner,)):
            with self.assertRaisesRegex(ValueError, 'independent qualification'):
                four.runner_sources(runner, None, meta)
            with self.assertRaisesRegex(ValueError, 'exact qualified child'):
                four.runner_sources(runner, path, dict(source_sha256={}))

    def test_failed_qualification_or_changed_runner_rejects(self):
        runner, path, record, meta = self.fixture()
        with mock.patch.object(four, 'ALLOWED_RUNNERS', (runner,)):
            four.atomic_json(path, dict(record, exit_code=2))
            with self.assertRaisesRegex(ValueError, 'did not pass'):
                four.runner_sources(runner, path, meta)
            four.atomic_json(path, record)
            runner.write_text('# modified after qualification', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'runner binding'):
                four.runner_sources(runner, path, meta)

    def test_old_source_path_cannot_be_rebound_silently(self):
        with self.assertRaisesRegex(ValueError, 'conflicting frozen source'):
            four.merge_sources({'frozen.py': 'a'*64}, {'frozen.py': 'b'*64})
        self.assertEqual(four.merge_sources({'same.py': 'a'*64}, {'same.py': 'a'*64}),
            {'same.py': 'a'*64})

    def test_cached_runner_needs_no_bridge_and_keeps_parent_source_map(self):
        runner, path, _, meta = self.fixture()
        with mock.patch.object(four, 'ALLOWED_RUNNERS', (runner,)), \
             mock.patch.object(four, 'CACHED_RUNNER', runner):
            child, _, qualification = four.runner_sources(runner, None, meta)
            self.assertEqual(child, meta['source_sha256']); self.assertIsNone(qualification)
            with self.assertRaisesRegex(ValueError, 'does not take'):
                four.runner_sources(runner, path, meta)


class CopiedChildLoopTests(old_tests.ChildControllerTests):
    """Exercise the five frozen fake-child scenarios against the new loop."""
    def setUp(self):
        super().setUp()
        actual = four.supervise
        def run(protocol, meta):
            protocol['wall_limit_s'] = min(14400., protocol['wall_limit_s'])
            protocol['runner'] = four.rel(four.CACHED_RUNNER)
            with mock.patch.object(four, 'validate_handoff'), \
                 mock.patch.object(four, 'utc_now', return_value=START):
                return actual(protocol, meta)
        self.enterContext(mock.patch.object(old_tests, 'campaign', four))
        self.enterContext(mock.patch.object(four, 'supervise', side_effect=run))

    def test_child_command_uses_the_protocol_runner(self):
        self.scenario(kind='failure')
        self.assertEqual(Path(self.commands[0][2]), four.CACHED_RUNNER)


class StopAndLaunchTests(old_tests.TemporaryWorkspace):
    def setup_protocol(self):
        parent = self.folder/'parent'; parent.mkdir()
        meta = old_tests.metadata()
        meta['observer'] = dict(supplied_material=dict(
            eligible_crossed=dict(particle_count=0, mass_kg=0.),
            crossed_under_annular_shoulder=dict(particle_count=0, mass_kg=0.)))
        return dict(output=four.rel(self.folder/'campaign'), initial_parent=four.rel(parent),
            runner=four.rel(four.BRIDGE_RUNNER), initial_snapshot=dict(checkpoint='parent.npz'),
            requested_macro_dt_s=.0004, wall_limit_s=14400., source_sha256={},
            child_source_sha256={}, minimum_dt_s=1e-7, max_halvings=8,
            phases=copy.deepcopy(old_tests.PHASES), python=sys.executable), meta

    def test_each_new_run_starts_a_new_bounded_deadline_and_honors_boundary_stop(self):
        protocol, meta = self.setup_protocol()
        actual_write = four.atomic_json
        def write(path, value):
            actual_write(path, value)
            if Path(path).name == 'progress.json' and value['status'] == 'STARTING':
                (Path(path).parent/'STOP_AFTER_STAGE').write_text('requested boundary stop', encoding='utf-8')
        with mock.patch.object(four, 'validate_handoff'), \
             mock.patch.object(four, 'utc_now', return_value=START), \
             mock.patch.object(four, 'singleton_and_awake', return_value=(mock.Mock(), 123)), \
             mock.patch.object(four, 'atomic_json', side_effect=write), \
             mock.patch.object(four.subprocess, 'Popen') as launch, \
             contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(four.supervise(protocol, meta), 0)
        launch.assert_not_called()
        final = four.read(ROOT/protocol['output']/'final_summary.json')
        self.assertEqual(final['status'], 'STOPPED_AFTER_STAGE_REQUEST')
        self.assertEqual(four.aware_time(final['latest_stop_at'])-four.aware_time(final['started_at']),
            timedelta(hours=4))

    def test_unqualified_parent_is_rejected_before_any_child_or_output(self):
        protocol, meta = self.setup_protocol()
        with mock.patch.object(four, 'validate_handoff', side_effect=ValueError('invalid complete GPU receipt')), \
             mock.patch.object(four.subprocess, 'Popen') as launch:
            with self.assertRaisesRegex(ValueError, 'complete GPU receipt'):
                four.supervise(protocol, meta)
        launch.assert_not_called(); self.assertFalse((ROOT/protocol['output']).exists())


def load_tests(loader, standard_tests, pattern):
    for name in loader.getTestCaseNames(old_tests.DecisionTests):
        if name != 'test_ten_hour_authorization_cannot_be_extended':
            standard_tests.addTest(old_tests.DecisionTests(name))
    for cls in (old_tests.FailureClassificationTests, old_tests.ParentReceiptTests):
        standard_tests.addTests(loader.loadTestsFromTestCase(cls))
    return standard_tests


if __name__ == '__main__':
    unittest.main(verbosity=2)

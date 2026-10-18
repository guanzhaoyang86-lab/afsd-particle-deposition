"""Read-only/mocked checks for the new CPU continuation administration.

No nonlinear solve or real process is started by these tests. The immutable
physical worker is compared against the previous qualified repair controller.
"""
from pathlib import Path
import ast
import importlib.util
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / 'outputs/short_bar_event_continue_20261015_v2/controller.py'
OLD = ROOT / 'outputs/short_bar_event_repair_20261015_v1/repair.py'


def definition(path, name):
    return next(node for node in ast.parse(path.read_text(encoding='utf-8')).body
                if isinstance(node, ast.FunctionDef) and node.name == name)


class MockAdministrativeScheduling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpu_continue_admin_review', CURRENT)
        cls.controller = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.controller)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='afsd_admin_review_')
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.pair = self.directory / 'pair'
        self.pair.mkdir()
        self.ctrl = self.controller
        for target, name, value in ((self.ctrl, 'RUN', self.directory),
                                    (self.ctrl, 'rel', lambda p: str(p)),
                                    (self.ctrl.old, 'launch_options', lambda: {})):
            context = patch.object(target, name, value)
            context.start()
            self.addCleanup(context.stop)

    def invoke_branch(self, maximum=25):
        return self.ctrl.run_branch(self.pair, {'source_sha256': {'fixed.py': 'a' * 64}},
            'coarse', 'parent.npz', 'b' * 64, .000474609375, maximum,
            {'maximum_single_solve_s': 34.5509886})

    def test_pause_prevents_cost_calculation_and_process_launch(self):
        (self.directory / 'PAUSE_REQUEST.json').write_text('{}', encoding='utf-8')
        with patch.object(self.ctrl.subprocess, 'Popen') as process, \
                patch.object(self.ctrl.costs, 'forecast') as forecast:
            with self.assertRaisesRegex(ValueError, 'USER_PAUSE_REQUESTED'):
                self.invoke_branch()
            process.assert_not_called()
            forecast.assert_not_called()

    def test_deadline_reserve_prevents_process_launch(self):
        with patch.object(self.ctrl, 'remaining', return_value=649.99), \
                patch.object(self.ctrl.subprocess, 'Popen') as process, \
                patch.object(self.ctrl.costs, 'forecast') as forecast:
            with self.assertRaisesRegex(ValueError, 'ABSOLUTE_DEADLINE_RESERVE'):
                self.invoke_branch()
            process.assert_not_called()
            forecast.assert_not_called()

    def test_over_budget_forecast_prevents_process_launch(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl.subprocess, 'Popen') as process, \
                patch.object(self.ctrl.costs, 'forecast',
                             return_value={'forecast_s': 450.01, 'eligible': False}):
            with self.assertRaisesRegex(ValueError, 'BRANCH_COST_REASSESSMENT_REQUIRED'):
                self.invoke_branch()
            process.assert_not_called()

    def test_zero_remaining_pair_solves_prevents_process_launch(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl.subprocess, 'Popen') as process, \
                patch.object(self.ctrl.costs, 'forecast',
                             return_value={'forecast_s': 320., 'eligible': True}):
            with self.assertRaisesRegex(ValueError, 'BRANCH_COST_REASSESSMENT_REQUIRED'):
                self.invoke_branch(maximum=0)
            process.assert_not_called()

    def test_worker_receives_actual_single_solve_maximum(self):
        class StopBeforeAnyChild(RuntimeError):
            pass
        written = []
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl.costs, 'forecast',
                             return_value={'forecast_s': 320., 'eligible': True}), \
                patch.object(self.ctrl.old, 'verify'), \
                patch.object(self.ctrl, 'write', side_effect=lambda path, data: written.append((path, data))), \
                patch.object(self.ctrl.subprocess, 'Popen', side_effect=StopBeforeAnyChild) as process:
            with self.assertRaises(StopBeforeAnyChild):
                self.invoke_branch()
        process.assert_called_once()
        self.assertEqual(len(written), 1)
        request = written[0][1]
        self.assertEqual(request['candidate_forecast_s'], 34.5509886)
        self.assertEqual(request['forecast_s'], 320.)
        self.assertEqual(request['maximum_solves'], 12)
        self.assertEqual(request['maximum_events_per_branch'], 1)
        self.assertEqual(request['internal_limit_s'], 540.)
        self.assertEqual(request['external_limit_s'], 600.)

    def test_cost_stop_is_not_reported_as_physical_failure(self):
        expected = {
            'BRANCH_COST_REASSESSMENT_REQUIRED': 'STOPPED_COST_REASSESSMENT',
            'USER_PAUSE_REQUESTED': 'STOPPED_USER_PAUSE',
            'ABSOLUTE_DEADLINE_RESERVE': 'STOPPED_WALL_RESERVE',
            'PAIR_SCREEN_FAILED: test': 'STOPPED_PAIR_SCREEN_FAILURE',
            'BRANCH_FAILED: test': 'STOPPED_VALIDATION_OR_CHILD_FAILURE',
        }
        for message, status in expected.items():
            with self.subTest(message=message):
                self.assertEqual(self.ctrl.stopped_status(ValueError(message)), status)

    def test_wall_clock_deadline_cannot_be_extended_by_monotonic_time(self):
        with patch.object(self.ctrl, 'DEADLINE', 100.), \
                patch.object(self.ctrl, 'MONOTONIC_DEADLINE', 200.), \
                patch.object(self.ctrl.time, 'time', return_value=101.), \
                patch.object(self.ctrl.time, 'monotonic', return_value=190.):
            self.assertEqual(self.ctrl.remaining(), -1.)

    def test_wall_clock_rollback_cannot_extend_monotonic_deadline(self):
        with patch.object(self.ctrl, 'DEADLINE', 100.), \
                patch.object(self.ctrl, 'MONOTONIC_DEADLINE', 200.), \
                patch.object(self.ctrl.time, 'time', return_value=90.), \
                patch.object(self.ctrl.time, 'monotonic', return_value=201.):
            self.assertEqual(self.ctrl.remaining(), -1.)


class FrozenPhysicalAdministration(unittest.TestCase):
    def test_physical_worker_is_exactly_the_qualified_original_worker(self):
        self.assertEqual(ast.dump(definition(CURRENT, 'worker'), include_attributes=False),
                         ast.dump(definition(OLD, 'worker'), include_attributes=False))

    def test_pair_comparison_and_physical_calls_remain_identical(self):
        def selected(path):
            calls = [node for node in ast.walk(definition(path, 'run_pair'))
                     if isinstance(node, ast.Call)]
            return sorted(ast.dump(node, include_attributes=False) for node in calls
                          if ast.unparse(node.func).startswith(('physical.', 'verification.')))
        self.assertEqual(selected(CURRENT), selected(OLD))

    def test_launcher_and_controller_have_same_fixed_absolute_deadline(self):
        def assigned_deadline(path):
            tree = ast.parse(path.read_text(encoding='utf-8'))
            assignment = next(node for node in tree.body if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == 'DEADLINE' for target in node.targets))
            value = assignment.value
            if isinstance(value, ast.Call) and isinstance(value.func, ast.Attribute):
                self.assertEqual(value.func.attr, 'timestamp')
                value = value.func.value
            return ast.dump(value, include_attributes=False)
        self.assertEqual(assigned_deadline(CURRENT), assigned_deadline(CURRENT.parent / 'launch.py'))


if __name__ == '__main__':
    unittest.main()

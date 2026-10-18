"""Mocked whole-pair rollback and read-only invariants for v3 administration."""
from pathlib import Path
import ast
import copy
import importlib.util
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / 'outputs/short_bar_event_continue_20261015_v3/controller.py'
PREVIOUS = ROOT / 'outputs/short_bar_event_continue_20261015_v2/controller.py'


def definition(path, name):
    return next(node for node in ast.parse(path.read_text(encoding='utf-8')).body
                if isinstance(node, ast.FunctionDef) and node.name == name)


class AdministrativeRollback(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpu_multi_event_admin_review', CURRENT)
        cls.controller = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.controller)

    def setUp(self):
        self.ctrl = self.controller
        self.tmp = tempfile.TemporaryDirectory(prefix='afsd_multi_event_admin_')
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.written, self.progressed = [], []
        self.model = self.ctrl.costs.update(self.ctrl.costs.empty_model(), 160., 10, [10.] * 10)
        self.meta = dict(screened_continuation_checkpoint='verified_parent.npz',
            screened_continuation_sha256=self.ctrl.PARENT_SHA,
            source_sha256={}, dt_s=.000474609375,
            screened_end_time_s=.07613671875, screened_end_step=404,
            observer={'supplied_material': {'crossed': 7}})
        self.evidence = dict(branches=[{'measured_single_solve_s': 20.}], accepted=False,
            failed_candidate_used_as_parent=False)

        def fake_read(path):
            if Path(path).name == 'component_checks.json':
                return {'exit_code': 0, 'source_sha256': {}}
            if Path(path).name == 'manifest.json':
                return self.meta
            if Path(path).name == 'cost_model.json':
                return self.model
            raise AssertionError('Unexpected mocked read: ' + str(path))

        changes = (
            (self.ctrl, 'RUN', self.directory / 'run'),
            (self.ctrl, 'BASE', self.directory),
            (self.ctrl, 'rel', lambda p: str(p)),
            (self.ctrl, 'read', fake_read),
            (self.ctrl, 'sha', lambda p: self.ctrl.PARENT_SHA),
            (self.ctrl, 'write', lambda p, d: self.written.append((Path(p), copy.deepcopy(d)))),
            (self.ctrl, 'progress', lambda **d: self.progressed.append(copy.deepcopy(d))),
            (self.ctrl, 'calibrate', lambda: (copy.deepcopy(self.meta), dict(self.model))),
            (self.ctrl, 'resume_recent', lambda meta, model: (meta, model)),
            (self.ctrl, 'source_closure', lambda source: dict(source)),
            (self.ctrl.old, 'verify', lambda source: None),
            (self.ctrl.physical, 'adaptation_audit', lambda: {'mocked': True}),
        )
        for target, name, value in changes:
            context = patch.object(target, name, value)
            context.start()
            self.addCleanup(context.stop)

    def test_retry_recomputes_whole_pair_from_same_complete_parent(self):
        accepted = dict(self.meta, screened_continuation_checkpoint='new_verified.npz',
                        screened_continuation_sha256='a' * 64, screened_end_time_s=.07625537109375,
                        screened_end_step=406)
        with patch.object(self.ctrl, 'remaining', side_effect=[1000., 0.]), \
                patch.object(self.ctrl, 'retry_evidence', return_value=copy.deepcopy(self.evidence)) as retry, \
                patch.object(self.ctrl, 'run_pair', side_effect=[ValueError('BRANCH_FAILED: fine_2'),
                                                               (accepted, self.model)]) as pairs:
            self.assertEqual(self.ctrl.run(), 0)
        calls = [call.args for call in pairs.call_args_list]
        self.assertEqual([call[0] for call in calls], ['continuation_0001_h1', 'continuation_0001_h2'])
        self.assertEqual([call[1] for call in calls], ['verified_parent.npz'] * 2)
        self.assertEqual([call[2] for call in calls], [self.ctrl.PARENT_SHA] * 2)
        self.assertEqual([call[3] for call in calls], [self.meta['dt_s'] / 2., self.meta['dt_s'] / 4.])
        self.assertEqual(retry.call_args.args[1:], ('verified_parent.npz', self.ctrl.PARENT_SHA))
        promoted = [row for row in self.progressed if row['status'] == 'CONTINUATION_VERIFIED']
        self.assertEqual(len(promoted), 1)
        self.assertEqual(promoted[0]['last_complete_checkpoint'], 'new_verified.npz')

    def test_initial_halving_counts_toward_four_total_limit(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence', return_value=copy.deepcopy(self.evidence)) as retry, \
                patch.object(self.ctrl, 'run_pair', side_effect=ValueError('BRANCH_FAILED: coarse')) as pairs:
            self.assertEqual(self.ctrl.run(), 2)
        self.assertEqual(pairs.call_count, 4)
        self.assertEqual(retry.call_count, 3)
        calls = [call.args for call in pairs.call_args_list]
        self.assertEqual([call[0] for call in calls], ['continuation_0001_h%d' % level for level in range(1, 5)])
        self.assertEqual([call[3] for call in calls], [self.meta['dt_s'] / 2. ** level for level in range(1, 5)])
        self.assertFalse(any(row['status'] == 'CONTINUATION_VERIFIED' for row in self.progressed))

    def test_other_physical_failure_is_not_retried(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence') as retry, \
                patch.object(self.ctrl, 'run_pair', side_effect=ValueError('PAIR_SCREEN_FAILED: test')) as pairs:
            self.assertEqual(self.ctrl.run(), 2)
        pairs.assert_called_once()
        retry.assert_not_called()
        self.assertEqual(self.progressed[-1]['status'], 'STOPPED_PAIR_SCREEN_FAILURE')

    def test_ineligible_branch_failure_does_not_start_second_pair(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence', side_effect=ValueError('ineligible unknown event')), \
                patch.object(self.ctrl, 'run_pair', side_effect=ValueError('BRANCH_FAILED: coarse')) as pairs:
            self.assertEqual(self.ctrl.run(), 2)
        pairs.assert_called_once()
        self.assertFalse(any(row['status'] == 'CONTINUATION_VERIFIED' for row in self.progressed))

    def test_success_keeps_reduced_dt_for_next_complete_pair(self):
        first = dict(self.meta, screened_continuation_checkpoint='first_verified.npz',
                     screened_continuation_sha256='a' * 64, screened_end_time_s=.0763740234375,
                     screened_end_step=406)
        second = dict(first, screened_continuation_checkpoint='second_verified.npz',
                      screened_continuation_sha256='b' * 64, screened_end_time_s=.076611328125,
                      screened_end_step=408)
        with patch.object(self.ctrl, 'remaining', side_effect=[1000., 1000., 0.]), \
                patch.object(self.ctrl, 'retry_evidence') as retry, \
                patch.object(self.ctrl, 'run_pair', side_effect=[(first, self.model), (second, self.model)]) as pairs:
            self.assertEqual(self.ctrl.run(), 0)
        calls = [call.args for call in pairs.call_args_list]
        self.assertEqual([call[3] for call in calls], [self.meta['dt_s'] / 2.] * 2)
        self.assertEqual(calls[1][1:3], ('first_verified.npz', 'a' * 64))
        retry.assert_not_called()

    def test_failed_solve_only_raises_floor_without_complete_sample(self):
        before = copy.deepcopy(self.model)
        raised = self.ctrl.failed_solve_floor(self.model,
            {'branches': [{'measured_single_solve_s': 35.}]} )
        self.assertEqual(self.model, before)
        self.assertEqual(raised['sample_count'], before['sample_count'])
        self.assertEqual(raised['maximum_single_solve_s'], 35.)
        self.assertEqual(raised['maximum_projected_branch_s'], before['maximum_projected_branch_s'])

    def test_invalid_failed_solve_duration_is_rejected(self):
        for value in (float('nan'), float('inf'), -1., 0., True):
            with self.subTest(value=value):
                with self.assertRaises((ValueError, TypeError)):
                    self.ctrl.failed_solve_floor(self.model,
                        {'branches': [{'measured_single_solve_s': value}]})


class FrozenV3Physics(unittest.TestCase):
    def test_worker_branch_and_pair_physics_unchanged_from_v2(self):
        for name in ('worker', 'run_branch', 'run_pair'):
            with self.subTest(name=name):
                self.assertEqual(ast.dump(definition(CURRENT, name), include_attributes=False),
                                 ast.dump(definition(PREVIOUS, name), include_attributes=False))

    def test_deadline_is_not_reset_by_retry_campaign(self):
        def assigned(path):
            tree = ast.parse(path.read_text(encoding='utf-8'))
            return next(node for node in tree.body if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == 'DEADLINE' for target in node.targets))
        self.assertEqual(ast.dump(assigned(CURRENT), include_attributes=False),
                         ast.dump(assigned(PREVIOUS), include_attributes=False))


if __name__ == '__main__':
    unittest.main()

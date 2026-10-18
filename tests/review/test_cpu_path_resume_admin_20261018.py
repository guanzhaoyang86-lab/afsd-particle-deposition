"""Mock v5 continuation control without physical solves or process launches."""
from pathlib import Path
import ast
import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest
from datetime import datetime
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / 'outputs/short_bar_event_continue_20261018_v5/controller.py'
PREVIOUS = ROOT / 'outputs/short_bar_event_continue_20261015_v4/controller.py'


def definition(path, name):
    return next(node for node in ast.parse(path.read_text(encoding='utf-8')).body
                if isinstance(node, ast.FunctionDef) and node.name == name)


class LocalIntervalAdministration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpu_local_interval_admin_review', CURRENT)
        cls.controller = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.controller)

    def setUp(self):
        self.ctrl = self.controller
        self.tmp = tempfile.TemporaryDirectory(prefix='afsd_local_interval_admin_')
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.saved, self.progressed, self.completed = {}, [], []
        self.pause_at_start = False
        self.model = self.ctrl.costs.update(self.ctrl.costs.empty_model(), 160., 10, [10.] * 10)
        self.meta = dict(screened_continuation_checkpoint='verified_parent.npz',
            screened_continuation_sha256=self.ctrl.PARENT_SHA, source_sha256={},
            dt_s=.00011865234375, screened_end_time_s=.0796962890625,
            screened_end_step=454, observer={'supplied_material': {'crossed': 7}})
        self.parents = {'verified_parent.npz': self.meta}
        self.evidence = dict(branches=[{'measured_single_solve_s': 20.}], accepted=False,
                             failed_candidate_used_as_parent=False)
        self.configuration = dict(schema='CPU_LOCAL_INTERVAL_RUN_CONFIGURATION_V1',
            started_utc='2026-10-16T00:00:00+00:00', deadline_utc='2026-10-16T05:00:00+00:00',
            wall_seconds=18000.)

        def fake_read(path):
            path = Path(path)
            if path in self.saved:
                return copy.deepcopy(self.saved[path])
            if path.name == 'component_checks.json':
                return {'exit_code': 0, 'source_sha256': {}}
            if path.name == 'manifest.json':
                return copy.deepcopy(self.meta)
            if path.name == 'cost_model.json':
                return dict(self.model)
            if path.name == 'run.configuration.json':
                return dict(self.configuration)
            raise AssertionError('Unexpected mocked read: ' + str(path))

        def fake_write(path, value):
            self.saved[Path(path)] = copy.deepcopy(value)
            if Path(path).name == 'instance.lock.json' and self.pause_at_start:
                (self.ctrl.RUN / 'PAUSE_REQUEST.json').write_text('{}', encoding='utf-8')

        changes = (
            (self.ctrl, 'RUN', self.directory / 'run'),
            (self.ctrl, 'BASE', self.directory),
            (self.ctrl, 'CONFIGURATION', self.directory / 'run.configuration.json'),
            (self.ctrl, 'CONFIGURATION_SHA', 'c' * 64),
            (self.ctrl, 'DEADLINE', 1792126800.),
            (self.ctrl, 'MONOTONIC_DEADLINE', 18000.),
            (self.ctrl, 'rel', lambda p: str(p)),
            (self.ctrl, 'read', fake_read),
            (self.ctrl, 'sha', lambda p: 'c' * 64 if Path(p).name == 'run.configuration.json' else self.ctrl.PARENT_SHA),
            (self.ctrl, 'write', fake_write),
            (self.ctrl, 'progress', lambda **d: self.progressed.append(copy.deepcopy(d))),
            (self.ctrl, 'restore_parent', lambda: (copy.deepcopy(self.meta), dict(self.model))),
            (self.ctrl, 'seed_pending', lambda meta, model:
                (self.ctrl.intervals.start(meta['screened_end_time_s'], meta['dt_s']), model)),
            (self.ctrl, 'source_closure', lambda source: dict(source)),
            (self.ctrl.old, 'verify', lambda source: None),
            (self.ctrl.physical, 'adaptation_audit', lambda: {'mocked': True}),
        )
        for target, name, value in changes:
            context = patch.object(target, name, value, create=True)
            context.start()
            self.addCleanup(context.stop)
        guard = patch.object(self.ctrl.subprocess, 'Popen',
                             side_effect=AssertionError('real process forbidden in administrative tests'))
        self.process_guard = guard.start()
        self.addCleanup(guard.stop)
        physics_guard = patch.object(self.ctrl.physical, 'complete_branch',
                             side_effect=AssertionError('physical solve forbidden in administrative tests'))
        self.physical_guard = physics_guard.start()
        self.addCleanup(physics_guard.stop)

    def tearDown(self):
        self.process_guard.assert_not_called()
        self.physical_guard.assert_not_called()

    def successful_pair(self, name, parent, parent_sha, dt, scale, sources, cost, purpose):
        before = self.parents[parent]
        self.assertEqual(parent_sha, before['screened_continuation_sha256'])
        checkpoint = 'verified_%s.npz' % name
        after = dict(before, screened_continuation_checkpoint=checkpoint,
            screened_continuation_sha256=('%064x' % (len(self.completed) + 1)),
            screened_end_time_s=before['screened_end_time_s'] + dt,
            screened_end_step=before['screened_end_step'] + 2, dt_s=dt)
        self.parents[checkpoint] = after
        self.completed.append(after)
        return after, cost

    def test_six_eligible_failures_then_success_have_no_retry_counter_cap(self):
        attempts = []
        def pair(*args):
            attempts.append(args)
            if len(attempts) <= 6:
                raise ValueError('BRANCH_FAILED: fine_2')
            return self.successful_pair(*args)
        with patch.object(self.ctrl, 'remaining', side_effect=lambda: 0. if self.completed else 1000.), \
                patch.object(self.ctrl, 'retry_evidence', return_value=copy.deepcopy(self.evidence)) as retry, \
                patch.object(self.ctrl.intervals, 'accept', wraps=self.ctrl.intervals.accept) as accept, \
                patch.object(self.ctrl, 'run_pair', side_effect=pair):
            self.assertEqual(self.ctrl.run(), 0)
        self.assertEqual(len(attempts), 7)
        self.assertEqual(retry.call_count, 6)
        accept.assert_called_once()
        self.assertEqual([args[0] for args in attempts],
                         ['piece_000001_attempt_%06d' % i for i in range(1, 8)])
        self.assertEqual([args[1] for args in attempts], ['verified_parent.npz'] * 7)
        self.assertEqual([args[2] for args in attempts], [self.ctrl.PARENT_SHA] * 7)
        for earlier, later in zip(attempts, attempts[1:]):
            self.assertAlmostEqual(later[3], earlier[3] / 2., delta=2e-17)

    def test_two_completed_local_pieces_restore_nominal_macro_interval(self):
        def seed(meta, model):
            state = self.ctrl.intervals.start(meta['screened_end_time_s'], meta['dt_s'])
            return self.ctrl.intervals.split(state, meta['screened_end_time_s']), model
        with patch.object(self.ctrl, 'seed_pending', side_effect=seed), \
                patch.object(self.ctrl, 'remaining', side_effect=lambda: 0. if len(self.completed) >= 3 else 1000.), \
                patch.object(self.ctrl, 'retry_evidence') as retry, \
                patch.object(self.ctrl, 'run_pair', side_effect=self.successful_pair) as pairs:
            self.assertEqual(self.ctrl.run(), 0)
        calls = [call.args for call in pairs.call_args_list]
        self.assertEqual(len(calls), 3)
        self.assertAlmostEqual(calls[0][3], self.meta['dt_s'] / 2., delta=2e-17)
        self.assertAlmostEqual(calls[1][3], self.meta['dt_s'] / 2., delta=2e-17)
        self.assertAlmostEqual(calls[2][3], self.meta['dt_s'], delta=2e-17)
        self.assertEqual(calls[1][1], self.completed[0]['screened_continuation_checkpoint'])
        self.assertEqual(calls[2][1], self.completed[1]['screened_continuation_checkpoint'])
        retry.assert_not_called()

    def test_unknown_branch_failure_cannot_trigger_split_or_second_attempt(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence', side_effect=ValueError('unrecognized failure')), \
                patch.object(self.ctrl.intervals, 'split', wraps=self.ctrl.intervals.split) as split, \
                patch.object(self.ctrl.intervals, 'accept', wraps=self.ctrl.intervals.accept) as accept, \
                patch.object(self.ctrl, 'run_pair', side_effect=ValueError('BRANCH_FAILED: coarse')) as pairs:
            self.assertEqual(self.ctrl.run(), 2)
        pairs.assert_called_once()
        split.assert_not_called()
        accept.assert_not_called()

    def test_pair_screen_failure_never_enters_event_retry(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence') as retry, \
                patch.object(self.ctrl.intervals, 'accept', wraps=self.ctrl.intervals.accept) as accept, \
                patch.object(self.ctrl, 'run_pair', side_effect=ValueError('PAIR_SCREEN_FAILED: test')) as pairs:
            self.assertEqual(self.ctrl.run(), 2)
        pairs.assert_called_once()
        retry.assert_not_called()
        accept.assert_not_called()
        self.assertEqual(self.progressed[-1]['status'], 'STOPPED_PAIR_SCREEN_FAILURE')

    def test_cost_reassessment_does_not_trigger_geometric_retry(self):
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence') as retry, \
                patch.object(self.ctrl, 'run_pair', side_effect=ValueError('BRANCH_COST_REASSESSMENT_REQUIRED')):
            self.assertEqual(self.ctrl.run(), 2)
        retry.assert_not_called()
        self.assertEqual(self.progressed[-1]['status'], 'STOPPED_COST_REASSESSMENT')

    def test_invalid_success_endpoint_is_not_promoted_or_popped(self):
        def invalid_pair(*args):
            result, model = self.successful_pair(*args)
            result = dict(result, screened_end_time_s=self.meta['screened_end_time_s'] + args[3] / 2.)
            return result, model
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence') as retry, \
                patch.object(self.ctrl, 'run_pair', side_effect=invalid_pair) as pairs:
            self.assertEqual(self.ctrl.run(), 2)
        pairs.assert_called_once()
        retry.assert_not_called()
        promoted = [row for row in self.progressed if row.get('last_complete_checkpoint')
                    and row['last_complete_checkpoint'] != 'verified_parent.npz']
        self.assertEqual(promoted, [])

    def test_pause_prevents_first_pair(self):
        self.pause_at_start = True
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'run_pair') as pairs:
            self.assertEqual(self.ctrl.run(), 0)
        pairs.assert_not_called()
        self.assertEqual(self.progressed[-1]['status'], 'STOPPED_USER_PAUSE')

    def test_deadline_prevents_first_pair(self):
        with patch.object(self.ctrl, 'remaining', side_effect=[1000., 649.]), \
                patch.object(self.ctrl, 'run_pair') as pairs:
            self.assertEqual(self.ctrl.run(), 0)
        pairs.assert_not_called()
        self.assertEqual(self.progressed[-1]['status'], 'STOPPED_WALL_RESERVE')

    def test_pause_between_eligible_attempts_prevents_second_pair(self):
        def failed(*args):
            (self.ctrl.RUN / 'PAUSE_REQUEST.json').write_text('{}', encoding='utf-8')
            raise ValueError('BRANCH_FAILED: fine_2')
        with patch.object(self.ctrl, 'remaining', return_value=1000.), \
                patch.object(self.ctrl, 'retry_evidence', return_value=copy.deepcopy(self.evidence)), \
                patch.object(self.ctrl, 'run_pair', side_effect=failed) as pairs:
            self.assertEqual(self.ctrl.run(), 0)
        pairs.assert_called_once()
        self.assertEqual(self.progressed[-1]['status'], 'STOPPED_USER_PAUSE')

    def test_deadline_between_attempts_prevents_second_pair(self):
        attempted = []
        def failed(*args):
            attempted.append(args)
            raise ValueError('BRANCH_FAILED: fine_2')
        with patch.object(self.ctrl, 'remaining', side_effect=lambda: 0. if attempted else 1000.), \
                patch.object(self.ctrl, 'retry_evidence', return_value=copy.deepcopy(self.evidence)), \
                patch.object(self.ctrl, 'run_pair', side_effect=failed) as pairs:
            self.assertEqual(self.ctrl.run(), 0)
        pairs.assert_called_once()
        self.assertEqual(self.progressed[-1]['status'], 'STOPPED_WALL_RESERVE')


class ResumeSavedInterval(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpu_path_resume_saved_review', CURRENT)
        cls.ctrl = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.ctrl)

    def setUp(self):
        self.begin = .07969628906249966
        self.nominal = .00011865234375000002
        initial = self.ctrl.intervals.start(self.begin, self.nominal)
        for _ in range(3):
            initial = self.ctrl.intervals.split(initial, self.begin)
        self.dt = self.ctrl.intervals.current_dt(initial, self.begin)
        self.end = self.begin + self.dt
        self.pending = self.ctrl.intervals.accept(initial, self.end)
        self.meta = dict(parent_sha256=self.ctrl.ANCHOR_SHA, parent_time_s=self.begin,
            screened_end_time_s=self.end, dt_s=self.dt)
        self.records = {
            'initial_local_intervals.json': dict(parent_sha256=self.ctrl.ANCHOR_SHA,
                nominal_dt_s=self.nominal, scheduler=initial),
            'scheduler.json': dict(self.pending, pending_targets=list(self.pending['pending_targets'])),
            'source_bridge.json': dict(initial_local_intervals_sha256=
                self.ctrl.RESUME_PINS['initial_local_intervals.json']),
        }
        self.saved = {}
        patches = (
            patch.object(self.ctrl, 'pinned_resume', side_effect=lambda name: copy.deepcopy(self.records[name])),
            patch.object(self.ctrl, 'write', side_effect=lambda path, value: self.saved.update({Path(path).name: value})),
            patch.object(self.ctrl.subprocess, 'Popen', side_effect=AssertionError('no process')),
            patch.object(self.ctrl.physical, 'complete_branch', side_effect=AssertionError('no solve')),
        )
        for item in patches:
            item.start(); self.addCleanup(item.stop)

    def test_restore_exact_three_remaining_boundaries_not_last_piece_dt(self):
        pending, model = self.ctrl.seed_pending(self.meta, {'unchanged': True})
        self.assertEqual(pending, self.pending)
        self.assertEqual(len(pending['pending_targets']), 3)
        self.assertEqual(pending['nominal_dt'], self.nominal)
        self.assertGreater(pending['nominal_dt'], self.meta['dt_s'])
        self.assertEqual(model, {'unchanged': True})
        self.assertFalse(self.saved['initial_local_intervals.json']['failed_candidate_used_as_parent'])

    def test_remaining_boundaries_complete_then_original_nominal_returns(self):
        pending, _ = self.ctrl.seed_pending(self.meta, {})
        for end in tuple(pending['pending_targets']):
            pending = self.ctrl.intervals.accept(pending, end)
        self.assertEqual(pending['pending_targets'], ())
        following = self.ctrl.intervals.start(pending['accepted_time_s'], pending['nominal_dt'])
        self.assertAlmostEqual(self.ctrl.intervals.current_dt(following, following['accepted_time_s']),
                               self.nominal, delta=1e-17)

    def test_changed_saved_boundary_is_rejected(self):
        self.records['scheduler.json']['pending_targets'][0] += 1e-7
        with self.assertRaisesRegex(ValueError, 'reconstructed'):
            self.ctrl.seed_pending(self.meta, {})
        self.assertEqual(self.saved, {})

    def test_changed_saved_nominal_is_rejected(self):
        self.records['scheduler.json']['nominal_dt'] = self.dt
        with self.assertRaisesRegex(ValueError, 'reconstructed'):
            self.ctrl.seed_pending(self.meta, {})

    def test_changed_source_bridge_pin_is_rejected(self):
        self.records['source_bridge.json']['initial_local_intervals_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'anchor'):
            self.ctrl.seed_pending(self.meta, {})

    def test_partial_or_wrong_endpoint_cannot_be_used(self):
        self.meta['screened_end_time_s'] -= self.dt / 2.
        with self.assertRaisesRegex(ValueError, 'endpoint'):
            self.ctrl.seed_pending(self.meta, {})

    def test_parent_interval_does_not_match_pending_first_piece_is_rejected(self):
        self.meta['dt_s'] *= 2.
        with self.assertRaisesRegex(ValueError, 'cover'):
            self.ctrl.seed_pending(self.meta, {})

    def test_actual_parent_time_and_initial_scheduler_must_match(self):
        self.meta['parent_time_s'] += 1e-6
        with self.assertRaisesRegex(ValueError, 'identity'):
            self.ctrl.seed_pending(self.meta, {})


class PinnedResumeEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpu_path_resume_pinned_review', CURRENT)
        cls.ctrl = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.ctrl)

    def test_changed_cached_resume_hash_rejected_before_read(self):
        with patch.object(self.ctrl, 'sha', return_value='0' * 64), \
             patch.object(self.ctrl, 'read') as reader:
            with self.assertRaisesRegex(ValueError, 'evidence changed'):
                self.ctrl.pinned_resume('scheduler.json')
        reader.assert_not_called()

    def test_unregistered_resume_file_rejected_before_read(self):
        with patch.object(self.ctrl, 'sha') as digest, patch.object(self.ctrl, 'read') as reader:
            with self.assertRaisesRegex(ValueError, 'evidence changed'):
                self.ctrl.pinned_resume('failed_candidate.npz')
        digest.assert_not_called()
        reader.assert_not_called()

    def test_parent_verification_uses_original_reader_once_and_no_ancestor_loop(self):
        function = definition(CURRENT, 'restore_parent')
        calls = [node for node in ast.walk(function) if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute) and node.func.attr == 'verify_pair']
        self.assertEqual(len(calls), 1)
        self.assertEqual(ast.unparse(calls[0]), 'legacy_verification.verify_pair(PARENT)')
        self.assertNotIn('calibrate()', ast.unparse(function))
        self.assertNotIn('resume_recent(', ast.unparse(function))


class ImmediateParentCostRestore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpu_path_resume_cost_review', CURRENT)
        cls.ctrl = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.ctrl)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='afsd_v5_cost_mock_')
        self.addCleanup(self.tmp.cleanup)
        self.parent = Path(self.tmp.name) / 'parent'
        self.meta = dict(screened_continuation_sha256=self.ctrl.PARENT_SHA,
            parent_sha256=self.ctrl.ANCHOR_SHA, parent_checkpoint='complete_old_parent.npz',
            source_sha256={'original.py': 'a' * 64})
        self.model = self.ctrl.costs.empty_model()
        self.model.update(sample_count=99, maximum_projected_branch_s=224.41016451428732,
                          maximum_single_solve_s=34.55098859999998)
        self.records = {
            'source_bridge.json': dict(to_source_sha256=self.meta['source_sha256'],
                parent_sha256=self.ctrl.ANCHOR_SHA, parent_checkpoint=self.meta['parent_checkpoint'],
                recent_parent_readback_sha256=self.ctrl.RESUME_PINS['recent_parent_readback.json']),
            'recent_parent_readback.json': dict(status='PASS_COMPLETE_RECENT_PARENT',
                parent_sha256=self.ctrl.ANCHOR_SHA, model=self.model),
        }
        for role in self.ctrl.ROLES:
            branch = self.parent / role
            branch.mkdir(parents=True)
            (branch / 'manifest.json').write_text(json.dumps({'solves': 1}), encoding='utf-8')
            (branch / 'supervision.json').write_text(json.dumps({'wall_seconds': 10.}), encoding='utf-8')
            (branch / 'physical_solve_cost.json').write_text(json.dumps({'wall_seconds': 9.}), encoding='utf-8')
        self.saved = {}
        changes = (
            (self.ctrl, 'PARENT', self.parent),
            (self.ctrl, 'sha', lambda path: self.ctrl.PARENT_MANIFEST_SHA),
            (self.ctrl, 'rel', lambda path: str(path)),
            (self.ctrl, 'read', lambda path: json.loads(Path(path).read_text(encoding='utf-8'))),
            (self.ctrl, 'pinned_resume', lambda name: copy.deepcopy(self.records[name])),
            (self.ctrl, 'write', lambda path, value: self.saved.update({Path(path).name: value})),
        )
        for target, name, value in changes:
            context = patch.object(target, name, value)
            context.start(); self.addCleanup(context.stop)
        verifier = patch.object(self.ctrl.legacy_verification, 'verify_pair',
            return_value=dict(meta=self.meta, status='PASS_CPU_EVENT_PAIR_READBACK'))
        self.verifier = verifier.start(); self.addCleanup(verifier.stop)
        for target, name in ((self.ctrl.subprocess, 'Popen'), (self.ctrl.physical, 'complete_branch')):
            guard = patch.object(target, name, side_effect=AssertionError('physical execution forbidden'))
            guard.start(); self.addCleanup(guard.stop)

    def test_one_complete_parent_readback_and_three_cost_samples_restore_102(self):
        meta, model = self.ctrl.restore_parent()
        self.verifier.assert_called_once_with(self.parent)
        self.assertEqual(meta, self.meta)
        self.assertEqual(model['sample_count'], 102)
        self.assertEqual(model['maximum_projected_branch_s'], self.model['maximum_projected_branch_s'])
        self.assertEqual(model['maximum_single_solve_s'], self.model['maximum_single_solve_s'])
        self.assertAlmostEqual(self.ctrl.costs.forecast(model)['forecast_s'], 366.615246771431)
        self.assertTrue(self.saved['cost_calibration.json']['incomplete_pair_samples_excluded'])

    def test_parent_receipt_failure_prevents_cost_reuse(self):
        self.verifier.side_effect = ValueError('complete receipt invalid')
        with self.assertRaisesRegex(ValueError, 'receipt invalid'):
            self.ctrl.restore_parent()
        self.assertEqual(self.saved, {})

    def test_source_bridge_cannot_change_original_map(self):
        self.records['source_bridge.json']['to_source_sha256'] = {'original.py': 'b' * 64}
        with self.assertRaisesRegex(ValueError, 'bridge mismatch'):
            self.ctrl.restore_parent()
        self.assertEqual(self.saved, {})

    def test_cached_cost_anchor_must_match_actual_complete_parent_chain(self):
        self.records['recent_parent_readback.json']['parent_sha256'] = 'f' * 64
        with self.assertRaisesRegex(ValueError, 'cost evidence anchor'):
            self.ctrl.restore_parent()
        self.assertEqual(self.saved, {})


class FrozenLocalIntervalPhysics(unittest.TestCase):
    def test_worker_and_complete_pair_are_unchanged_from_v4(self):
        for name in ('worker', 'run_pair'):
            with self.subTest(name=name):
                self.assertEqual(ast.dump(definition(CURRENT, name), include_attributes=False),
                                 ast.dump(definition(PREVIOUS, name), include_attributes=False))

    def test_worker_branch_and_pair_have_same_execution_as_v4(self):
        for name in ('worker', 'run_branch', 'run_pair'):
            with self.subTest(name=name):
                self.assertEqual(ast.dump(definition(CURRENT, name), include_attributes=False),
                                 ast.dump(definition(PREVIOUS, name), include_attributes=False))


class ImmutableLaunchDeadline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpu_local_config_review', CURRENT)
        cls.controller = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.controller)

    def setUp(self):
        self.ctrl = self.controller
        self.tmp = tempfile.TemporaryDirectory(prefix='afsd_local_config_review_')
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name).resolve()
        self.path = self.directory / 'run.configuration.json'
        self.valid = dict(schema='CPU_LOCAL_INTERVAL_RUN_CONFIGURATION_V1',
            started_utc='2026-10-16T00:00:00+00:00', deadline_utc='2026-10-16T05:00:00+00:00',
            wall_seconds=18000.)
        changes = (
            (self.ctrl, 'BASE', self.directory),
            (self.ctrl, 'read', lambda path: json.loads(Path(path).read_text(encoding='utf-8'))),
            (self.ctrl, 'sha', lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()),
            (self.ctrl, 'CONFIGURATION', None),
            (self.ctrl, 'CONFIGURATION_SHA', None),
            (self.ctrl, 'DEADLINE', None),
            (self.ctrl, 'MONOTONIC_DEADLINE', None),
        )
        for target, name, value in changes:
            context = patch.object(target, name, value)
            context.start()
            self.addCleanup(context.stop)

    def fixture(self, **changes):
        self.path.write_text(json.dumps(dict(self.valid, **changes)), encoding='utf-8')
        return hashlib.sha256(self.path.read_bytes()).hexdigest()

    def test_reloading_configuration_never_grants_a_fresh_absolute_budget(self):
        digest = self.fixture()
        original_bytes = self.path.read_bytes()
        begin = datetime.fromisoformat(self.valid['started_utc']).timestamp()
        end = datetime.fromisoformat(self.valid['deadline_utc']).timestamp()
        with patch.object(self.ctrl.time, 'time', return_value=begin), \
                patch.object(self.ctrl.time, 'monotonic', return_value=100.):
            self.ctrl.configure(self.path, digest)
            self.assertEqual(self.ctrl.DEADLINE, end)
            first_monotonic = self.ctrl.MONOTONIC_DEADLINE
            self.assertEqual(self.ctrl.remaining(), 18000.)
        with patch.object(self.ctrl.time, 'time', return_value=begin + 600.), \
                patch.object(self.ctrl.time, 'monotonic', return_value=700.):
            self.ctrl.configure(self.path, digest)
            self.assertEqual(self.ctrl.DEADLINE, end)
            self.assertEqual(self.ctrl.MONOTONIC_DEADLINE, first_monotonic)
            self.assertEqual(self.ctrl.remaining(), 17400.)
        self.assertEqual(self.path.read_bytes(), original_bytes)

    def test_changed_configuration_sha_is_rejected(self):
        digest = self.fixture()
        self.fixture(wall_seconds=19000.)
        with self.assertRaisesRegex(ValueError, 'configuration changed'):
            self.ctrl.configure(self.path, digest)
        self.assertIsNone(self.ctrl.CONFIGURATION)
        self.assertIsNone(self.ctrl.DEADLINE)

    def test_start_duration_deadline_contradiction_is_rejected(self):
        digest = self.fixture(wall_seconds=17999.)
        with self.assertRaisesRegex(ValueError, 'disagree'):
            self.ctrl.configure(self.path, digest)
        self.assertIsNone(self.ctrl.DEADLINE)

    def test_naive_start_or_deadline_is_rejected(self):
        for changes in ({'started_utc': '2026-10-16T00:00:00'},
                        {'deadline_utc': '2026-10-16T05:00:00'},
                        {'started_utc': '2026-10-16T00:00:00', 'deadline_utc': '2026-10-16T05:00:00'}):
            with self.subTest(changes=changes):
                digest = self.fixture(**changes)
                with self.assertRaisesRegex(ValueError, 'timezone-aware'):
                    self.ctrl.configure(self.path, digest)
                self.assertIsNone(self.ctrl.DEADLINE)

    def test_invalid_wall_duration_is_rejected_before_configuration_changes(self):
        for duration in (0., -1., 649., float('nan'), float('inf'), True, '18000'):
            with self.subTest(duration=duration):
                digest = self.fixture(wall_seconds=duration)
                with self.assertRaises(ValueError):
                    self.ctrl.configure(self.path, digest)
                self.assertIsNone(self.ctrl.DEADLINE)
                self.assertIsNone(self.ctrl.CONFIGURATION)

    def test_unregistered_configuration_path_is_rejected(self):
        digest = self.fixture()
        alternate = self.directory / 'another.configuration.json'
        alternate.write_bytes(self.path.read_bytes())
        with self.assertRaisesRegex(ValueError, 'registered launch configuration'):
            self.ctrl.configure(alternate, digest)


if __name__ == '__main__':
    unittest.main()

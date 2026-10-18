"""Independent incremental-driver administration tests; no GPU solving."""
import argparse
import ast
import copy
import inspect
import json
import multiprocessing
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_incremental_20261010 as driver
import afsd_incremental_worker_20261010 as worker
import afsd_adaptive_dt_20261010 as adaptive
from afsd_locked_history_20261010 import LockedHistory, HistoryLockError

SHA = 'a'*64
PARENT_SHA = 'b'*64


def administrative_child(channel, protocol_sha):
    """Real spawn/Pipe transport; imports call no GPU environment or solver."""
    try:
        grant = worker.request_grant(channel, protocol_sha, 'start')
        channel.send(dict(administrative_child_finished=True, controller_pid=grant['controller_pid']))
    finally:
        channel.close()


class ProtocolTests(unittest.TestCase):
    def fixture(self):
        sources = {'scripts/frozen_fixture.py': SHA}
        template = dict(source_sha256=sources, reference_sha256={'old_fixture':SHA},
            search_tensor_reference_stage='old-replay', administrative_preflight_statistics={'old':True},
            search_tensor_protocol=driver.PROTOCOL)
        meta = dict(source_sha256=copy.deepcopy(sources), screened_continuation_checkpoint='outputs/parent_fixture.npz',
            screened_continuation_sha256=PARENT_SHA, maximum_candidate_wall_s=20.)
        return template, meta

    def make(self, template=None, meta=None, **changes):
        default_template, default_meta = self.fixture()
        return driver.stage_protocol(default_template if template is None else template,
            default_meta if meta is None else meta, driver.BASE/'parent_fixture',
            changes.get('dt', .00011865234375), changes.get('administrative', {'admin_fixture':SHA}),
            changes.get('session_record', {'session_fixture':SHA}), changes.get('refs', {'old_fixture':SHA}))

    def test_static_audit_reuses_exact_candidate_and_pair_gates(self):
        result = driver.static_audit()
        self.assertEqual(result['status'], 'PASS_UNCHANGED_INITIAL_RECEIPT_ASSERTIONS_AND_PHYSICAL_PAIR')
        self.assertFalse(result['runtime_global_rebinding'])
        self.assertTrue(result['inherited_physical_audit']['pair_gates_ast_identical'])
        self.assertTrue(result['inherited_physical_audit']['complete_candidate_ast_identical'])
        self.assertEqual(driver.CHECKS, driver.bridge.CHECKS)
        self.assertNotEqual(driver.INCREMENTAL_CHECKS, driver.bridge.CHECKS)

    def test_worker_calls_original_full_pair_and_no_replacement_solver(self):
        tree = ast.parse(inspect.getsource(worker.stage_worker))
        calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute) and node.func.attr == '_run_pair_original']
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].func.value.id, 'physical')
        self.assertIs(worker.physical, driver.bridge)

    def test_stage_inherits_exact_physical_source_and_complete_parent(self):
        template, meta = self.fixture(); before = copy.deepcopy(template)
        result = self.make(template, meta)
        self.assertEqual(result['source_sha256'], template['source_sha256'])
        self.assertEqual(result['parent_implementation_source_sha256'], template['source_sha256'])
        self.assertEqual(result['parent_sha256'], PARENT_SHA)
        self.assertEqual(result['parent_checkpoint'], meta['screened_continuation_checkpoint'])
        self.assertEqual(result['reference_sha256'], {'old_fixture':SHA,'admin_fixture':SHA,'session_fixture':SHA})
        self.assertNotIn('search_tensor_reference_stage', result)
        self.assertNotIn('administrative_preflight_statistics', result)
        self.assertEqual(template, before)

    def test_inconsistent_physical_source_rejects(self):
        template, meta = self.fixture(); meta['source_sha256']['scripts/frozen_fixture.py'] = 'c'*64
        with self.assertRaises(ValueError):
            self.make(template, meta)

    def test_reference_admin_and_session_digest_conflicts_reject(self):
        for changes in (dict(administrative={'old_fixture':'c'*64}),
                dict(session_record={'admin_fixture':'c'*64}),
                dict(refs={'scripts/frozen_fixture.py':'c'*64})):
            with self.subTest(changes=changes), self.assertRaisesRegex(ValueError, 'conflicting evidence'):
                self.make(**changes)

    def test_original_cost_limit_rejects(self):
        template, meta = self.fixture(); meta['maximum_candidate_wall_s'] = 94.
        with self.assertRaisesRegex(ValueError, 'forecast limit'):
            self.make(template, meta)

    def test_driver_rejects_unbounded_requested_inputs_before_creating_output(self):
        defaults = dict(output='outputs/gpu_acceleration_20261007/uncreated_driver_test_fixture',
            parent='outputs/gpu_acceleration_20261007/unread_parent', wall_seconds=600., max_stages=3,
            maximum_dt=adaptive.DEFAULT_MAXIMUM_DT_S, shoulder_observe_ms=5., target_time=.22, replay_stage=None)
        for bad in (dict(wall_seconds=14400.01), dict(wall_seconds=float('nan')), dict(max_stages=2001),
                dict(max_stages=0), dict(maximum_dt=.0051), dict(shoulder_observe_ms=20.1),
                dict(target_time=.22001), dict(target_time=0.)):
            with self.subTest(bad=bad), mock.patch.object(driver, 'LockedHistory') as pool:
                with self.assertRaises(ValueError):
                    driver.run(argparse.Namespace(**dict(defaults, **bad)))
                pool.assert_not_called()

    def test_actual_latest_parent_starts_small_and_all_ids_allow_later_growth(self):
        path = driver.BASE/'campaign_4h_20261009_v1/attempt_0052/manifest.json'
        meta = json.loads(path.read_text(encoding='utf-8'))
        state = adaptive.initial_state(meta['screened_continuation_sha256'], meta['dt_s'])
        result = adaptive.normalized_screen(meta['pair_screen'], .0005, list(range(1,1382)))
        self.assertLess(result['ratio'], .05)
        self.assertEqual(state['nominal_dt_s'], meta['dt_s'])
        self.assertEqual(state['small_error_streak'], 0)


class GrantTests(unittest.TestCase):
    def protocol(self):
        return dict(source_sha256={'source':SHA},reference_sha256={'history':SHA},
            parent_checkpoint='outputs/parent_fixture.npz',parent_sha256=PARENT_SHA)

    def pool(self):
        result = mock.Mock(); result.sha.return_value = PARENT_SHA
        return result

    def test_private_pipe_grant_binds_exact_hash_and_phase(self):
        pool = self.pool(); protocol = self.protocol()
        response = driver.authorize_grant(dict(protocol_sha256=SHA,phase='start'),SHA,'start',pool,protocol)
        self.assertTrue(response['history_locked']); self.assertEqual(response['phase'],'start')
        self.assertEqual(response['protocol_sha256'],SHA)
        self.assertEqual(pool.verify.call_args_list,
            [mock.call(protocol['source_sha256']),mock.call(protocol['reference_sha256'])])

    def test_wrong_identity_order_and_extra_requests_reject_before_history_access(self):
        for request, phase in ((dict(protocol_sha256='c'*64,phase='start'),'start'),
                (dict(protocol_sha256=SHA,phase='end'),'start'),
                (dict(protocol_sha256=SHA,phase='start',extra=True),'start'),
                (dict(protocol_sha256=SHA,phase=None),None)):
            pool = self.pool()
            with self.subTest(request=request), self.assertRaises(ValueError):
                driver.authorize_grant(request,SHA,phase,pool,self.protocol())
            pool.verify.assert_not_called()

    def test_inactive_pool_or_changed_parent_cannot_grant(self):
        pool = self.pool(); pool.verify.side_effect = HistoryLockError('closed live session')
        with self.assertRaises(HistoryLockError):
            driver.authorize_grant(dict(protocol_sha256=SHA,phase='end'),SHA,'end',pool,self.protocol())
        pool = self.pool(); pool.sha.return_value = 'c'*64
        with self.assertRaisesRegex(ValueError, 'parent changed'):
            driver.authorize_grant(dict(protocol_sha256=SHA,phase='end'),SHA,'end',pool,self.protocol())

    def test_worker_rejects_missing_wrong_or_expired_live_authorization(self):
        good = dict(status='PASS_LIVE_HISTORY_GRANT',phase='start',protocol_sha256=SHA,history_locked=True)
        for changes in (dict(status='FAIL'),dict(phase='end'),dict(protocol_sha256='c'*64),
                dict(history_locked=False),dict(history_locked=1)):
            channel = mock.Mock();channel.poll.return_value=True;channel.recv.return_value=dict(good,**changes)
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                worker.request_grant(channel,SHA,'start')
        channel=mock.Mock();channel.poll.return_value=False
        with self.assertRaisesRegex(ValueError,'did not authorize'):
            worker.request_grant(channel,SHA,'start')
        channel.recv.assert_not_called()

    def test_worker_rejects_invalid_phase_or_hash_before_sending_request(self):
        for protocol_sha,phase in ((SHA,None),(SHA,'middle'),('bad','start'),('G'*64,'start')):
            channel=mock.Mock()
            with self.subTest(protocol_sha=protocol_sha,phase=phase), self.assertRaises(ValueError):
                worker.request_grant(channel,protocol_sha,phase)
            channel.send.assert_not_called()

    def test_real_spawn_administrative_handshake_with_live_locked_history(self):
        with tempfile.TemporaryDirectory(prefix='incremental_driver_test_',dir=driver.BASE) as temporary:
            folder=Path(temporary); parent=folder/'parent.txt'; source=folder/'source.txt'; reference=folder/'reference.txt'
            for item in (parent,source,reference): item.write_text(item.name,encoding='utf-8')
            pool=LockedHistory(ROOT)
            context=multiprocessing.get_context('spawn');server,child_channel=context.Pipe(duplex=True)
            child=context.Process(target=administrative_child,args=(child_channel,SHA))
            try:
                protocol=dict(parent_checkpoint=driver.rel(parent),parent_sha256=pool.sha(parent),
                    source_sha256={driver.rel(source):pool.sha(source)},reference_sha256={driver.rel(reference):pool.sha(reference)})
                child.start();child_channel.close()
                self.assertTrue(server.poll(30.),'administrative child failed to request live grant')
                request=server.recv()
                response=driver.authorize_grant(request,SHA,'start',pool,protocol)
                server.send(response)
                self.assertTrue(server.poll(30.))
                done=server.recv(); self.assertTrue(done['administrative_child_finished'])
                self.assertEqual(done['controller_pid'],response['controller_pid'])
                child.join(timeout=30.);self.assertEqual(child.exitcode,0)
                pool.close()
                with self.assertRaises(HistoryLockError):
                    driver.authorize_grant(dict(protocol_sha256=SHA,phase='end'),SHA,'end',pool,protocol)
            finally:
                if child.pid is not None and child.is_alive(): child.terminate();child.join(timeout=5.)
                server.close();child_channel.close();pool.close()


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Bounded live-history authorization regression; no GPU solve or long waits."""
import ast
import copy
import inspect
import multiprocessing
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest import mock

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import afsd_incremental_worker_20261010 as old_worker
import afsd_incremental_worker_resume_20261010 as worker
import run_gpu_incremental_resume_20261010 as driver
from afsd_locked_history_20261010 import LockedHistory

SHA='a'*64


class Clock:
    def __init__(self,now=100.): self.now=now
    def __call__(self): return self.now


class DelayedChannel:
    def __init__(self,clock,delay=0.,response=None,recv_delay=0.,send_delay=0.):
        self.clock=clock;self.delay=delay;self.recv_delay=recv_delay;self.send_delay=send_delay
        self.response=response if response is not None else dict(status='PASS_LIVE_HISTORY_GRANT',
            phase='start',protocol_sha256=SHA,history_locked=True)
        self.sent=[];self.waits=[];self.received=0
    def send(self,value): self.sent.append(value);self.clock.now+=self.send_delay
    def poll(self,timeout):
        self.waits.append(timeout);self.clock.now+=min(timeout,self.delay)
        return self.delay<=timeout
    def recv(self):
        self.received+=1;self.clock.now+=self.recv_delay
        return self.response


def two_grant_child(channel,protocol_sha):
    try:
        deadline=time.perf_counter()+525.
        grants=[worker.request_grant(channel,protocol_sha,phase,deadline) for phase in ('start','end')]
        channel.send(dict(complete=True,phases=[g['phase'] for g in grants]))
    finally: channel.close()


class BoundedGrantTests(unittest.TestCase):
    def request(self,channel,clock,**kwargs):
        with mock.patch.object(worker.time,'perf_counter',clock):
            return worker.request_grant(channel,SHA,'start',**kwargs)

    def test_60_second_controller_delay_succeeds_and_would_exceed_old_45_seconds(self):
        clock=Clock();channel=DelayedChannel(clock,delay=60.)
        result=self.request(channel,clock,deadline=625.)
        self.assertEqual(result['status'],'PASS_LIVE_HISTORY_GRANT')
        self.assertEqual(channel.waits,[180.]);self.assertEqual(clock.now,160.)

    def test_180_second_cap_rejects_missing_grant_without_extended_stage_budget(self):
        for deadline in (None,625.,10000.):
            clock=Clock();channel=DelayedChannel(clock,delay=181.)
            with self.subTest(deadline=deadline),self.assertRaisesRegex(ValueError,'did not authorize'):
                self.request(channel,clock,deadline=deadline)
            self.assertEqual(channel.waits,[180.]);self.assertEqual(channel.received,0)

    def test_remaining_525_second_stage_budget_truncates_wait(self):
        # Stage began at 100; end grant starts after 500 seconds of stage work.
        clock=Clock(600.);channel=DelayedChannel(clock,delay=26.)
        with self.assertRaisesRegex(ValueError,'did not authorize'):
            self.request(channel,clock,deadline=100.+525.)
        self.assertEqual(channel.waits,[25.]);self.assertEqual(channel.received,0)

    def test_grant_can_pass_within_remaining_stage_budget(self):
        clock=Clock(600.);channel=DelayedChannel(clock,delay=24.)
        self.request(channel,clock,deadline=625.)
        self.assertEqual(channel.waits,[25.]);self.assertEqual(clock.now,624.)

    def test_expired_and_invalid_deadlines_reject_before_send(self):
        for deadline in (99.,100.,float('nan'),float('inf'),float('-inf'),'625'):
            clock=Clock();channel=DelayedChannel(clock)
            with self.subTest(deadline=deadline),self.assertRaises(ValueError):
                self.request(channel,clock,deadline=deadline)
            self.assertEqual(channel.sent,[]);self.assertEqual(channel.waits,[])

    def test_send_delay_uses_same_absolute_grant_and_stage_deadline(self):
        clock=Clock();channel=DelayedChannel(clock,send_delay=5.,delay=14.)
        self.request(channel,clock,deadline=120.)
        self.assertEqual(channel.waits,[15.])
        clock=Clock();channel=DelayedChannel(clock,send_delay=21.)
        with self.assertRaisesRegex(ValueError,'did not authorize'):
            self.request(channel,clock,deadline=120.)
        self.assertEqual(channel.waits,[])

    def test_response_after_deadline_is_rejected_even_if_poll_was_ready(self):
        for deadline,delay,recv_delay in ((120.,19.,2.),(625.,179.,2.),(120.,20.,0.)):
            clock=Clock();channel=DelayedChannel(clock,delay=delay,recv_delay=recv_delay)
            with self.subTest(deadline=deadline,delay=delay),self.assertRaisesRegex(ValueError,'deadline exhausted'):
                self.request(channel,clock,deadline=deadline)

    def test_wrong_response_status_identity_phase_and_lock_state_reject(self):
        good=dict(status='PASS_LIVE_HISTORY_GRANT',phase='start',protocol_sha256=SHA,history_locked=True)
        for changes in (dict(status='FAIL'),dict(phase='end'),dict(protocol_sha256='b'*64),
                dict(history_locked=False),dict(history_locked=1)):
            clock=Clock();channel=DelayedChannel(clock,response=dict(good,**changes))
            with self.subTest(changes=changes),self.assertRaisesRegex(ValueError,'invalid live history'):
                self.request(channel,clock)
        clock=Clock();channel=DelayedChannel(clock,response='PASS')
        with self.assertRaisesRegex(ValueError,'invalid live history'): self.request(channel,clock)

    def test_invalid_request_rejects_before_send(self):
        for sha,phase in ((SHA,'middle'),('bad','start'),('G'*64,'end'),(SHA,None)):
            channel=mock.Mock()
            with self.subTest(sha=sha,phase=phase),self.assertRaises(ValueError):
                worker.request_grant(channel,sha,phase,625.)
            channel.send.assert_not_called()

    def test_closed_or_broken_channel_cannot_authorize(self):
        for method,error in (('send',BrokenPipeError()),('poll',OSError('closed')),('recv',EOFError())):
            channel=mock.Mock();channel.poll.return_value=True
            getattr(channel,method).side_effect=error
            with self.subTest(method=method),self.assertRaises(type(error)):
                worker.request_grant(channel,SHA,'start')

    def test_old_physical_stage_ast_changes_only_two_grant_deadline_arguments(self):
        old=ast.parse(inspect.getsource(old_worker.stage_worker))
        new=ast.parse(inspect.getsource(worker.stage_worker))
        calls=[n for n in ast.walk(new) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)
            and n.func.id=='request_grant']
        self.assertEqual(len(calls),2)
        for call in calls:
            self.assertEqual(len(call.args),4)
            self.assertEqual(ast.dump(call.args[-1]),ast.dump(ast.parse('began+525.',mode='eval').body))
            call.args.pop()
        self.assertEqual(ast.dump(new,include_attributes=False),ast.dump(old,include_attributes=False))
        # Outside grant and stage bodies the whole module remains exactly frozen.
        old=ast.parse(inspect.getsource(old_worker));new=ast.parse(inspect.getsource(worker))
        for tree in (old,new):
            tree.body=[node for node in tree.body if not isinstance(node,ast.FunctionDef)
                or node.name not in ('request_grant','stage_worker')]
        self.assertEqual(ast.dump(new),ast.dump(old))

    def test_real_spawn_pipe_start_and_end_require_live_locked_history(self):
        with tempfile.TemporaryDirectory(prefix='grant_resume_test_',dir=driver.BASE) as temporary:
            folder=Path(temporary);parent=folder/'parent.txt';source=folder/'source.txt';reference=folder/'reference.txt'
            for item in (parent,source,reference): item.write_text(item.name,encoding='utf-8')
            pool=LockedHistory(ROOT);context=multiprocessing.get_context('spawn')
            server,child_channel=context.Pipe(duplex=True)
            child=context.Process(target=two_grant_child,args=(child_channel,SHA))
            try:
                protocol=dict(parent_checkpoint=driver.rel(parent),parent_sha256=pool.sha(parent),
                    source_sha256={driver.rel(source):pool.sha(source)},
                    reference_sha256={driver.rel(reference):pool.sha(reference)})
                child.start();child_channel.close()
                for phase in ('start','end'):
                    self.assertTrue(server.poll(30.),'child did not request grant')
                    request=server.recv()
                    response=driver.authorize_grant(request,SHA,phase,pool,protocol)
                    self.assertEqual(response['status'],'PASS_LIVE_HISTORY_GRANT')
                    server.send(response)
                self.assertTrue(server.poll(30.))
                self.assertEqual(server.recv(),dict(complete=True,phases=['start','end']))
                child.join(timeout=30.);self.assertEqual(child.exitcode,0)
            finally:
                if child.pid is not None and child.is_alive(): child.terminate();child.join(timeout=5.)
                server.close();child_channel.close();pool.close()


if __name__=='__main__': unittest.main(verbosity=2)

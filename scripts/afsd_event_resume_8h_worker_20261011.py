"""Isolated explicit event pair; original live-history authorization and limits."""
from pathlib import Path
import contextlib
import json
import os
import threading
import time
import traceback
import numpy as np
import run_gpu_visibility_event_delay_20261010 as physical


def request_grant(channel, protocol_sha, phase, deadline=None):
    physical.require(phase in ('start','end') and isinstance(protocol_sha,str) and len(protocol_sha)==64
        and all(c in '0123456789abcdef' for c in protocol_sha),'valid grant phase and protocol SHA required')
    began=time.perf_counter()
    physical.require(deadline is None or (isinstance(deadline,(int,float)) and np.isfinite(deadline)),
        'finite authorization deadline required')
    grant_deadline=began+180. if deadline is None else min(began+180.,float(deadline))
    physical.require(grant_deadline>began,'live history authorization deadline exhausted')
    channel.send(dict(protocol_sha256=protocol_sha,phase=phase))
    remaining=grant_deadline-time.perf_counter()
    physical.require(remaining>0. and channel.poll(remaining), 'live history controller did not authorize worker')
    response=channel.recv()
    physical.require(time.perf_counter()<grant_deadline,'live history authorization deadline exhausted')
    physical.require(isinstance(response,dict) and response.get('status')=='PASS_LIVE_HISTORY_GRANT' and response.get('phase')==phase
        and response.get('protocol_sha256')==protocol_sha and response.get('history_locked') is True,
        'invalid live history authorization')
    return response


def worker_limits(protocol,began,epoch_now,perf_now):
    physical.require(protocol.get('administrative_protocol')=='FROZEN_V5_EIGHT_HOUR_DYNAMIC_ADMINISTRATIVE_RESUME_20261011_V1',
        'registered eight-hour administrative protocol required')
    soft=protocol['physical_soft_deadline_s'];hard=protocol['internal_limit_s'];outer=protocol['external_limit_s']
    epoch=protocol['absolute_deadline_epoch_s']
    physical.require(all(type(v) in (int,float) and np.isfinite(v) and v>0 for v in (soft,hard,outer,epoch))
        and soft<=hard<=outer and soft<=3600. and hard<=3630. and outer<=3720.,'bounded dynamic worker limits required')
    physical.require(epoch<=1791733307.0 and epoch>epoch_now,'eight-hour absolute authorization expired or extended')
    physical.require(protocol['worker_deadline_s']==soft and protocol['dynamic_wall_budget']['physical_soft_deadline_s']==soft
        and protocol['dynamic_wall_budget']['internal_limit_s']==hard and protocol['dynamic_wall_budget']['external_limit_s']==outer,
        'worker metadata differs from registered budget')
    remaining=epoch-epoch_now
    return min(began+soft,perf_now+remaining),min(hard,remaining)


def stage_worker(directory, channel, expected_protocol_sha):
    output=Path(directory);protocol=physical.read(output/'protocol.json')
    report=dict(protocol,status='RUNNING_GPU_DIAGNOSTIC',screened_continuation_authorized=False)
    began=time.perf_counter();code=2;grants=[]
    def timeout():
        physical.core.write(output/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NOT_ACCEPTED'));os._exit(124)
    soft_deadline,hard_seconds=worker_limits(protocol,began,time.time(),time.perf_counter())
    timer=threading.Timer(hard_seconds,timeout);timer.daemon=True;timer.start()
    with (output/'run.log').open('x',encoding='utf-8',buffering=1) as stream, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
        try:
            physical.require(physical.core.sha(output/'protocol.json')==expected_protocol_sha,'worker protocol changed')
            physical.require(all(os.environ.get(k)=='1' for k in physical.shared.THREADS),'single CPU thread required')
            grants.append(request_grant(channel,expected_protocol_sha,'start',soft_deadline))
            physical.verify(protocol['source_sha256']);physical.verify(protocol['administrative_source_sha256'])
            physical.require(physical.core.sha(physical.ROOT/protocol['parent_checkpoint'])==protocol['parent_sha256'],'actual parent changed')
            report['adaptation_audit']=physical.adaptation_audit();report['gpu_environment']=physical.gpu.environment()
            physical.run_pair(output,protocol,report,soft_deadline)
            if protocol.get('incremental_replay_stage'):
                reference=physical.ROOT/protocol['incremental_replay_stage'];exact={};comparisons={}
                for role in ('coarse','fine_1','fine_2'):
                    pa,state=physical.load(output/role/'NOT_RESTARTABLE_candidate.npz')
                    prior,old_state=physical.load(reference/role/'NOT_RESTARTABLE_candidate.npz')
                    exact[role]=set(pa.properties)==set(prior.properties) and all(np.array_equal(prop.get_npy_array(),prior.properties[name].get_npy_array()) for name,prop in pa.properties.items())
                    comparisons[role]=physical.compare_cpu_endpoint(dict(pa=pa,state=state),reference/role/'NOT_RESTARTABLE_candidate.npz')
                report.update(incremental_replay_particle_arrays_bitwise_equal=exact,incremental_replay_comparisons=comparisons)
                physical.require(all(exact.values()),'administrative change altered replay particle arrays')
            code=0
        except Exception as error:
            physical.core.write(output/'WORKER_PRE_END_FAILURE.json',dict(status='RECORDED_WORKER_FAILURE_BEFORE_END_GRANT',error=repr(error)))
            report.update(status='FAIL_GPU_STAGE_NOT_ACCEPTED',error=repr(error),traceback=traceback.format_exc(),screened_continuation_authorized=False)
        finally:
            try:
                physical.verify(protocol['source_sha256']);physical.verify(protocol['administrative_source_sha256'])
                grants.append(request_grant(channel,expected_protocol_sha,'end',soft_deadline))
            except Exception as error:
                code=2;report.update(status='FAIL_GPU_STAGE_NOT_ACCEPTED',error=repr(error),traceback=traceback.format_exc(),screened_continuation_authorized=False)
            timer.cancel();channel.close()
            report.update(wall_seconds=time.perf_counter()-began,peak_working_set_bytes=physical.core.memory_bytes(),
                gpu_statistics=physical.gpu.statistics(),live_history_grants=grants,accepted=False,restart_allowed=False,accepted_time_advanced_s=0.)
            report['artifacts_sha256']={physical.rel(p):physical.core.sha(p) for p in output.rglob('*') if p.is_file() and p.name not in ('manifest.json','run.log','supervision.json')}
            physical.core.write(output/'manifest.json',report)
    raise SystemExit(code)

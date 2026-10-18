"""Administrative continuation of one existing GPU campaign and its deadline.

Only the launcher, origin handoff and remaining-wall-clock bookkeeping differ
from the frozen supervisor.  Physics helpers and GPU modules are not imported.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ctypes
import json
import math
import os
import subprocess
import sys
import time
import traceback
from supervise_gpu_campaign_20261008 import (
    require, read, sha, inside, rel, atomic_json, verify, launch_options,
    terminate_owned_tree, decision, verify_parent, snapshot, classify_failure,
    singleton_and_awake, phase_table, wall_now, ROOT, BASE, PASS, THREADS, LOCAL)

RUNNER = ROOT/'scripts/run_gpu_cached_preflight_20261008.py'
OLD_SUPERVISOR = ROOT/'scripts/supervise_gpu_campaign_20261008.py'
CONTRACT = ROOT/'docs/GPU_FAST_CONTINUATION_20261008.md'
TEST = ROOT/'tests/review/test_gpu_fast_supervisor_20261008.py'
CHECKS = BASE/'fast_supervisor_components_20261008_v1.json'
LAUNCHER_TEST = ROOT/'tests/review/test_gpu_preflight_cache_20261008.py'
LAUNCHER_CHECKS = BASE/'preflight_components_20261008_v1.json'
ORIGIN = BASE/'campaign_10h_20261008_v1'
AUTHORIZED_STOP_AT = '2026-10-09T02:19:45+08:00'


def utc_now():
    return datetime.now(timezone.utc)


def aware_time(value):
    result = datetime.fromisoformat(value)
    require(result.tzinfo is not None and result.utcoffset() is not None,
        'timezone-aware original campaign time required')
    return result.astimezone(timezone.utc)


def remaining_seconds(latest_stop_at, now=None):
    stop = aware_time(latest_stop_at)
    require(stop == aware_time(AUTHORIZED_STOP_AT), 'original absolute deadline must not change')
    now = utc_now() if now is None else now
    require(now.tzinfo is not None and now.utcoffset() is not None, 'aware current time required')
    return max(0., (stop-now.astimezone(timezone.utc)).total_seconds())


def budget_remaining(protocol, deadline):
    # The fixed monotonic limit never grows after a backward wall-clock change;
    # the original absolute UTC limit also stops a forward clock change.
    return max(0., min(deadline-time.perf_counter(), remaining_seconds(protocol['latest_stop_at'])))


def process_is_running(pid):
    require(isinstance(pid, int) and pid > 0, 'recorded positive origin supervisor PID required')
    require(os.name == 'nt', 'registered Windows GPU campaign required')
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.OpenProcess.argtypes = [ctypes.c_uint32, ctypes.c_bool, ctypes.c_uint32]
    kernel.OpenProcess.restype = ctypes.c_void_p
    kernel.GetExitCodeProcess.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint32)]
    kernel.CloseHandle.argtypes = [ctypes.c_void_p]
    handle = kernel.OpenProcess(0x1000, False, pid)
    if not handle:
        require(ctypes.get_last_error() == 87, 'cannot establish that the origin supervisor exited')
        return False
    try:
        exit_code = ctypes.c_uint32()
        require(kernel.GetExitCodeProcess(handle, ctypes.byref(exit_code)) != 0,
            'cannot read origin supervisor process state')
        return exit_code.value == 259
    finally:
        kernel.CloseHandle(handle)


def read_origin(origin, parent, allow_running=False):
    origin, parent = inside(Path(origin)), inside(Path(parent))
    require(origin == ORIGIN.resolve(), 'only the existing authorized origin campaign may be continued')
    protocol_path = origin/'protocol.json'
    origin_hash = sha(protocol_path)
    protocol, progress = read(protocol_path), read(origin/'progress.json')
    require(protocol['schema'] == 'GPU_TEN_HOUR_SUPERVISOR_20261008'
        and protocol['output'] == rel(origin), 'original campaign identity differs')
    verify(protocol['source_sha256'])
    limit = float(protocol['wall_limit_s'])
    require(math.isfinite(limit) and 0 < limit <= 36000., 'original finite budget must be at most ten hours')
    start, stop = aware_time(protocol['started_at']), aware_time(protocol['latest_stop_at'])
    require(0 < (stop-start).total_seconds() <= limit
        and stop == aware_time(AUTHORIZED_STOP_AT), 'original absolute deadline or duration differs')
    require(progress['latest_stop_at'] == protocol['latest_stop_at']
        and parent == inside(Path(progress['latest_pass_stage'])), 'parent must be the exact latest origin PASS')
    require(parent == inside(Path(protocol['initial_parent'])) or parent.is_relative_to(origin),
        'latest PASS must belong to the original campaign')
    finished = progress.get('finished') is True
    hashes = {rel(protocol_path): origin_hash}
    if not allow_running:
        require(finished and progress.get('exit_code') == 0 and progress.get('child_pid') is None,
            'origin must finish cleanly before launching a continuation')
        require(progress['status'] == 'STOPPED_USER_REQUEST',
            'this handoff requires the original clean STOPPED_USER_REQUEST boundary')
        final = origin/'final_summary.json'
        require(read(final) == progress, 'original final receipt differs from progress')
        require(not process_is_running(protocol['supervisor_pid']), 'origin supervisor has not exited')
        hashes.update({rel(origin/'progress.json'): sha(origin/'progress.json'), rel(final): sha(final)})
    require(sha(protocol_path) == origin_hash, 'original protocol changed during handoff audit')
    return protocol, progress, hashes


def validate_handoff(protocol, meta):
    origin, progress, hashes = read_origin(protocol['origin_campaign'], protocol['initial_parent'])
    require(hashes == protocol['origin_sha256'], 'origin handoff receipt changed after planning')
    require(origin['child_source_sha256'] == protocol['child_source_sha256'] == meta['source_sha256'],
        'frozen physical implementation changed at handoff')
    parent = ROOT/protocol['initial_parent']
    require(sha(parent/'manifest.json') == protocol['initial_parent_manifest_sha256']
        and progress['checkpoint_sha256'] == meta['screened_continuation_sha256']
        and progress['checkpoint'] == meta['screened_continuation_checkpoint']
        and progress['physical_time_s'] == meta['screened_end_time_s'],
        'latest complete origin parent changed after planning')


def plan(args):
    output, parent = inside(Path(args.output)), inside(Path(args.parent))
    require(output.is_relative_to(BASE) and output != BASE and not output.exists(), 'new GPU campaign folder required')
    require(0 < args.minimum_dt <= 1e-7 and 0 <= args.max_halvings <= 8,
        'registered temporal refinement bounds required')
    origin, progress, origin_hashes = read_origin(args.origin_campaign, parent,
        allow_running=bool(args.preflight_only))
    for path, expected_status in ((CHECKS, 'PASS'),
            (LAUNCHER_CHECKS, 'PASS_GPU_PREFLIGHT_CACHE_COMPONENTS')):
        checks = read(path)
        require(checks['status'] == expected_status and checks['exit_code'] == 0, 'executed administrative components required')
        verify(checks['source_sha256']); verify(checks['artifacts_sha256'])
    meta = verify_parent(parent, origin['child_source_sha256'])
    require(progress['checkpoint_sha256'] == meta['screened_continuation_sha256']
        and progress['checkpoint'] == meta['screened_continuation_checkpoint']
        and progress['physical_time_s'] == meta['screened_end_time_s'], 'origin snapshot does not identify exact PASS parent')
    phases = phase_table(read(ROOT/'configs/deposition_fill_then_track_20261007.json'))
    require(phases == origin['phases'], 'registered phase program changed')
    source = {rel(p): sha(p) for p in (Path(__file__).resolve(), CONTRACT, TEST, CHECKS,
        RUNNER, LAUNCHER_TEST, LAUNCHER_CHECKS, OLD_SUPERVISOR,
        ROOT/'scripts/supervise_afsd_deposition_20261007.py')}
    source.update(origin['source_sha256']); source.update(meta['source_sha256']); source.update(origin_hashes)
    # Planning and all source/parent checks consume the old authorization.
    remaining = remaining_seconds(origin['latest_stop_at'])
    return dict(schema='GPU_FIXED_DEADLINE_FAST_SUPERVISOR_20261008', output=rel(output),
        origin_campaign=rel(inside(Path(args.origin_campaign))), origin_sha256=origin_hashes,
        latest_stop_at=origin['latest_stop_at'], original_started_at=origin['started_at'],
        original_wall_limit_s=origin['wall_limit_s'], initial_parent=rel(parent),
        source_sha256=source, child_source_sha256=meta['source_sha256'],
        initial_parent_manifest_sha256=sha(parent/'manifest.json'),
        wall_limit_s=remaining, minimum_dt_s=args.minimum_dt, max_halvings=args.max_halvings,
        requested_macro_dt_s=meta['dt_s'], phases=phases,
        initial_snapshot=snapshot(meta,parent), first_decision=decision(meta,meta['dt_s'],remaining,phases,args.minimum_dt),
        python=str(ROOT/'.envs/pysph/python.exe'), stop_on_first_shoulder=True,
        accepted=False, formal_accepted_time_advanced_s=0., no_physical_module_imported=True), meta


def supervise(protocol,meta):
    verify(protocol['source_sha256']); validate_handoff(protocol,meta)
    began=time.perf_counter(); started=utc_now()
    remaining=remaining_seconds(protocol['latest_stop_at'],started); deadline=began+remaining
    protocol['planned_remaining_wall_s']=protocol['wall_limit_s']; protocol['wall_limit_s']=remaining
    output=ROOT/protocol['output']; output.mkdir(parents=True,exist_ok=False)
    protocol.update(started_at=started.astimezone(LOCAL).isoformat(timespec='seconds'),
        supervisor_pid=os.getpid(),sleep_prevention='PROCESS_LOCAL_SYSTEM_REQUIRED_DISPLAY_MAY_SLEEP')
    atomic_json(output/'protocol.json',protocol)
    parent=ROOT/protocol['initial_parent']; desired=protocol['requested_macro_dt_s']; completed=0; attempt=0; halvings=0
    progress=dict(status='STARTING',finished=False,completed_stage_count=0,attempt_count=0,
        started_at=protocol['started_at'],latest_stop_at=protocol['latest_stop_at'],supervisor_pid=os.getpid(),
        accepted=False,formal_accepted_time_advanced_s=0.,**protocol['initial_snapshot'])
    child_process=None; code=0;kernel=None;mutex=None
    def update(**fields):
        progress.update(fields); progress.update(updated_at=wall_now(),wall_seconds=time.perf_counter()-began,
            remaining_wall_seconds=budget_remaining(protocol,deadline))
        atomic_json(output/'progress.json',progress)
    update()
    try:
        kernel,mutex=singleton_and_awake()
        while True:
            verify(protocol['source_sha256'])
            if (output/'STOP_REQUESTED').exists():
                update(status='STOPPED_USER_REQUEST');break
            if (output/'STOP_AFTER_STAGE').exists():
                update(status='STOPPED_AFTER_STAGE_REQUEST');break
            if meta['observer']['supplied_material']['crossed_under_annular_shoulder']['mass_kg']>0:
                update(status='STOPPED_FIRST_SCREENED_NEW_MATERIAL_UNDER_SHOULDER');break
            action=decision(meta,desired,budget_remaining(protocol,deadline),protocol['phases'],protocol['minimum_dt_s'])
            if action['action']!='RUN':
                update(status=action['action'],next_decision=action);break
            attempt+=1; child=output/f'attempt_{attempt:04d}'; logfile=output/f'attempt_{attempt:04d}_launcher.log'
            command=[protocol['python'],'-B',str(RUNNER),'--mode','continue','--output',rel(child),
                '--cost',rel(BASE/'benchmark_v1'),'--validation',rel(BASE/'replay_v1'),
                '--parent',rel(parent),'--dt',repr(action['dt_s'])]
            update(status='RUNNING_GPU_CHILD',current_attempt=rel(child),attempt_count=attempt,
                current_dt_s=action['dt_s'],next_decision=action,current_parent=rel(parent),halvings_from_same_parent=halvings)
            launched=time.perf_counter(); termination=None; timed_out=False; requested_stop=False
            env=dict(os.environ,**{k:'1' for k in THREADS});env['PYTHONPATH']=''
            with logfile.open('x',encoding='utf-8') as stream:
                child_process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,**launch_options())
                update(child_pid=child_process.pid)
                while child_process.poll() is None:
                    requested_stop=(output/'STOP_REQUESTED').exists()
                    if budget_remaining(protocol,deadline)<=0. or time.perf_counter()>=launched+650. or requested_stop:
                        timed_out=not requested_stop;termination=terminate_owned_tree(child_process,stream)
                        child_process.wait(timeout=20.);break
                    wait_seconds=max(0.,min(5.,budget_remaining(protocol,deadline),launched+650.-time.perf_counter()))
                    try:child_process.wait(timeout=wait_seconds)
                    except subprocess.TimeoutExpired:update()
                exit_code=child_process.returncode
            child_process=None
            row=dict(attempt=attempt,stage=rel(child),parent=rel(parent),dt_s=action['dt_s'],command=command,
                exit_code=exit_code,wall_seconds=time.perf_counter()-launched,timeout=timed_out,termination=termination,
                log_sha256=sha(logfile),accepted=False)
            atomic_json(output/f'attempt_{attempt:04d}_receipt.json',row)
            if timed_out or requested_stop:
                update(status='STOPPED_USER_REQUEST' if requested_stop else 'STOPPED_WALL_DEADLINE_OR_CHILD_TIMEOUT');code=0 if requested_stop else 2;break
            if exit_code!=0:
                failure=read(child/'manifest.json') if (child/'manifest.json').exists() else None
                screen=read(child/'pair_screen.json') if (child/'pair_screen.json').exists() else None
                refine=classify_failure(failure,screen)=='REFINE'
                if refine:
                    # A failed screen is eligible only with intact same-parent,
                    # source, checkpoint and three successful candidate receipts.
                    verify(failure['source_sha256']);verify(failure['reference_sha256']);verify(failure['artifacts_sha256'])
                    require(failure['source_sha256']==protocol['child_source_sha256']
                        and failure['parent_sha256']==meta['screened_continuation_sha256'], 'failed attempt parent or source differs')
                    for role in ('coarse','fine_1','fine_2'):
                        require(read(child/role/'candidate_qualification.json')['status']=='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN',
                            'all candidate gates must pass before temporal refinement')
                if refine and halvings<protocol['max_halvings'] and action['dt_s']/2>=protocol['minimum_dt_s']:
                    desired=action['dt_s']/2;halvings+=1
                    update(status='RETRY_SMALLER_DT_FROM_UNCHANGED_PASS_PARENT',failed_attempt=rel(child),next_dt_s=desired);continue
                update(status='STOPPED_CHILD_FAILURE',failed_attempt=rel(child),error=failure.get('error') if failure else 'launcher/preflight failure',
                    child_exit_code=exit_code);code=2;break
            next_meta=verify_parent(child,protocol['child_source_sha256'],parent)
            require(next_meta['screened_end_time_s']>meta['screened_end_time_s'],'no screened physical progress')
            parent,meta=child,next_meta;completed+=1;halvings=0
            update(status='VERIFIED_PASS_READY_FOR_NEXT_STAGE',completed_stage_count=completed,child_pid=None,**snapshot(meta,parent))
            row.update(verified=True,manifest_sha256=sha(child/'manifest.json'),**snapshot(meta,parent))
            atomic_json(output/f'attempt_{attempt:04d}_receipt.json',row)
    except Exception as error:
        code=2;update(status='STOPPED_SUPERVISOR_VERIFICATION_FAILURE',error=repr(error),traceback=traceback.format_exc())
    finally:
        if child_process is not None and child_process.poll() is None:
            with (output/'termination.log').open('a',encoding='utf-8') as stream:terminate_owned_tree(child_process,stream)
        if kernel is not None:
            kernel.SetThreadExecutionState(0x80000000)
            kernel.CloseHandle(mutex)
        update(finished=True,exit_code=code,child_pid=None)
        atomic_json(output/'final_summary.json',progress)
    print(json.dumps(progress),flush=True)
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--parent', required=True)
    parser.add_argument('--origin-campaign', required=True)
    parser.add_argument('--minimum-dt', type=float, default=1e-7)
    parser.add_argument('--max-halvings', type=int, default=8)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    protocol, meta = plan(args)
    if args.preflight_only:
        print(json.dumps({k: protocol[k] for k in ('initial_parent', 'origin_campaign',
            'latest_stop_at', 'wall_limit_s', 'first_decision', 'initial_snapshot')}))
        return 0
    return supervise(protocol, meta)


if __name__ == '__main__':
    raise SystemExit(main())


"""Ten-hour, stdlib-only supervision of the frozen GPU deposition runner.

No physics functions or GPU modules are imported. Child stages retain every
original gate and source binding. Only a complete screened fine pair advances
the parent. All failures and partial attempts stay on disk.
"""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
import ctypes
import json
import math
import os
import subprocess
import sys
import time
import traceback
from supervise_afsd_deposition_20261007 import (
    require, read, sha, inside, rel, atomic_json, verify, launch_options, terminate_owned_tree)

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/gpu_acceleration_20261007'
RUNNER=ROOT/'scripts/run_afsd_gpu_20261007.py'
CONTRACT=ROOT/'docs/GPU_TEN_HOUR_CAMPAIGN_20261008.md'
TEST=ROOT/'tests/review/test_gpu_campaign_supervisor_20261008.py'
CHECKS=BASE/'supervisor_components_20261008_v1.json'
INITIAL=BASE/'continue_20261008_001414_714'
PASS='PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION'
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')
LOCAL=timezone(timedelta(hours=8))


def wall_now():
    return datetime.now(timezone.utc).astimezone(LOCAL).isoformat(timespec='seconds')


def phase_table(config):
    rows=[]; start=0.; previous_traverse=0.
    for p in config['phases']:
        end=start+float(p['duration_s'])
        require(math.isfinite(end) and end>start, 'invalid registered phase duration')
        rows.append(dict(start=start,end=end,name=p['name'],
            stationary=previous_traverse==0. and float(p['traverse_m_s'])==0.))
        start=end; previous_traverse=float(p['traverse_m_s'])
    return rows


def bounded_dt(current, desired, boundary):
    dt=min(desired,boundary-current)
    if dt == boundary-current:
        # Both saved half endpoints must land on the exact registered boundary;
        # adjust only representational rounding, never a saved particle/time.
        candidates=[dt]; lower=upper=dt
        for unused in range(8):
            lower=math.nextafter(lower,0.); upper=math.nextafter(upper,math.inf)
            candidates.extend((lower,upper))
        for candidate in candidates:
            if current+candidate==boundary and (current+candidate/2)+candidate/2==boundary:
                return candidate
        raise ValueError('cannot represent exact two-half phase endpoint')
    return dt


def decision(meta, proposed_dt, remaining_s, phases, minimum_dt=1e-7):
    t=float(meta['screened_end_time_s']); measured=float(meta['maximum_candidate_wall_s'])
    require(all(math.isfinite(v) for v in (t,measured,proposed_dt,remaining_s,minimum_dt))
        and t>=0 and measured>0 and 0<minimum_dt<=proposed_dt<=.005, 'finite bounded cost/time inputs required')
    forecast=30.+4.5*measured
    result=dict(forecast_s=forecast,child_timeout_s=max(0.,min(650.,remaining_s)))
    if remaining_s<forecast+20.:
        return dict(result,action='STOP_WALL_RESERVE')
    if forecast>450.:
        return dict(result,action='STOP_COST_REASSESSMENT')
    phase=next((p for p in phases if p['start']<=t<p['end']),None)
    if phase is None or not phase['stationary']:
        return dict(result,action='STOP_REGISTERED_STATIONARY_FILL_END')
    try:
        dt=bounded_dt(t,proposed_dt,phase['end'])
    except ValueError as error:
        return dict(result,action='STOP_UNREPRESENTABLE_PHASE_BOUNDARY',reason=str(error))
    if dt<minimum_dt:
        return dict(result,action='STOP_MINIMUM_DT_OR_BOUNDARY')
    return dict(result,action='RUN',dt_s=dt,phase=phase['name'],phase_end_s=phase['end'])


def classify_failure(meta, screen):
    # Only an observed temporal disagreement may halve the interval. Root,
    # energy, geometry, CUDA and source failures need separate diagnosis.
    if not isinstance(meta,dict) or meta.get('status')!='FAIL_GPU_STAGE_NOT_ACCEPTED':
        return 'STOP'
    if meta.get('error') != "ValueError('original coarse/fine screen failed; no GPU continuation')":
        return 'STOP'
    if not isinstance(screen,dict) or screen.get('passed') is not False or screen.get('crossing_ids_match') is not True:
        return 'STOP'
    if not isinstance(screen.get('global_screen'),dict) or not isinstance(screen.get('local_screen'),dict):
        return 'STOP'
    a=screen['global_screen'].get('pass'); b=screen['local_screen'].get('passed')
    if type(a) is not bool or type(b) is not bool or (a and b):
        return 'STOP'
    return 'REFINE'


def verify_parent(folder, expected_source=None, expected_parent=None):
    folder=inside(folder)
    meta=read(folder/'manifest.json'); supervision=read(folder/'supervision.json')
    require(meta['status']==PASS and supervision['exit_code']==0 and supervision['source_inputs_still_match']
        and supervision['manifest_sha256']==sha(folder/'manifest.json')
        and supervision['log_sha256']==sha(folder/'run.log'), 'invalid complete GPU receipt')
    require(meta['accepted'] is False and meta['restart_allowed'] is False and meta['diagnostic_only'] is True
        and meta['accepted_time_advanced_s']==0. and meta['screened_continuation_authorized'] is True,
        'only screened diagnostic GPU continuation is authorized')
    if expected_source is not None:
        require(meta['source_sha256']==expected_source, 'GPU implementation changed')
    verify(meta['source_sha256']); verify(meta['reference_sha256']); verify(meta['artifacts_sha256'])
    history=meta['screened_history']; screen=read(folder/'pair_screen.json')
    require(screen==meta['pair_screen'] and len(history)==2 and screen['passed']
        and screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'],
        'complete passing two-half screen required')
    require(screen['fine_checkpoint_sha256']==[r['sha256'] for r in history]
        and meta['screened_continuation_checkpoint']==history[-1]['checkpoint']
        and meta['screened_continuation_sha256']==history[-1]['sha256'], 'must continue exact last screened fine state')
    for half,row in enumerate(history,1):
        checkpoint=inside(Path(row['checkpoint']))
        require(checkpoint==(folder/f'fine_{half}'/'NOT_RESTARTABLE_candidate.npz').resolve()
            and sha(checkpoint)==row['sha256']
            and abs(row['time_s']-meta['parent_time_s']-half*meta['dt_s']/2)<=1e-13, 'half-step binding or time differs')
    require(meta['screened_end_time_s']==history[-1]['time_s']
        and sha(folder/'coarse/NOT_RESTARTABLE_candidate.npz')==screen['coarse_checkpoint_sha256'], 'coarse/end binding differs')
    require(sha(inside(Path(meta['parent_checkpoint'])))==meta['parent_sha256'], 'parent checkpoint changed')
    if expected_parent is not None:
        parent=read(expected_parent/'manifest.json')
        require(meta['parent_gpu_directory']==rel(expected_parent)
            and meta['parent_sha256']==parent['screened_continuation_sha256']
            and meta['screened_end_step']==parent['screened_end_step']+2
            and meta['parent_time_s']==parent['screened_end_time_s'], 'not the exact prior screened parent')
    return meta


def snapshot(meta,folder):
    supplied=meta['observer']['supplied_material']
    return dict(latest_pass_stage=rel(folder),physical_time_s=meta['screened_end_time_s'],
        screened_steps=meta['screened_end_step'],checkpoint=meta['screened_continuation_checkpoint'],
        checkpoint_sha256=meta['screened_continuation_sha256'],
        new_crossed_count=supplied['eligible_crossed']['particle_count'],
        new_crossed_mass_kg=supplied['eligible_crossed']['mass_kg'],
        new_under_shoulder_count=supplied['crossed_under_annular_shoulder']['particle_count'],
        new_under_shoulder_mass_kg=supplied['crossed_under_annular_shoulder']['mass_kg'])


def singleton_and_awake():
    if os.name!='nt':
        raise RuntimeError('registered Windows GPU campaign required')
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.CreateMutexW.argtypes=[ctypes.c_void_p,ctypes.c_bool,ctypes.c_wchar_p]
    kernel.CreateMutexW.restype=ctypes.c_void_p
    kernel.CloseHandle.argtypes=[ctypes.c_void_p]
    mutex=kernel.CreateMutexW(None,False,'Local\\AFSD_GPU_DEPOSITION_CAMPAIGN_20261008')
    if not mutex:
        raise RuntimeError('cannot acquire GPU campaign mutex')
    if ctypes.get_last_error()==183:
        kernel.CloseHandle(mutex)
        raise RuntimeError('another GPU campaign supervisor is already running')
    # Process-local idle-sleep prevention, display may still turn off.
    kernel.SetThreadExecutionState.argtypes=[ctypes.c_uint32]
    if kernel.SetThreadExecutionState(0x80000001)==0:
        kernel.CloseHandle(mutex)
        raise RuntimeError('cannot establish process-local awake request')
    return kernel,mutex


def plan(args):
    output=inside(Path(args.output)); parent=inside(Path(args.parent))
    require(output.is_relative_to(BASE) and output!=BASE and not output.exists(),'new GPU campaign folder required')
    require(math.isfinite(args.wall_seconds) and 650.<=args.wall_seconds<=36000.,'authorized maximum is ten hours')
    require(0<args.minimum_dt<=1e-7 and 0<=args.max_halvings<=8,'registered temporal refinement bounds required')
    checks=read(CHECKS)
    require(checks['status']=='PASS' and checks['exit_code']==0,'executed supervisor components required')
    verify(checks['source_sha256']); verify(checks['artifacts_sha256'])
    meta=verify_parent(parent)
    config=ROOT/'configs/deposition_fill_then_track_20261007.json'; phases=phase_table(read(config))
    source={rel(p):sha(p) for p in (Path(__file__).resolve(),CONTRACT,TEST,CHECKS,RUNNER,
        ROOT/'scripts/supervise_afsd_deposition_20261007.py')}
    source.update(meta['source_sha256'])
    return dict(schema='GPU_TEN_HOUR_SUPERVISOR_20261008',output=rel(output),initial_parent=rel(parent),
        source_sha256=source,child_source_sha256=meta['source_sha256'],
        initial_parent_manifest_sha256=sha(parent/'manifest.json'),
        wall_limit_s=args.wall_seconds,minimum_dt_s=args.minimum_dt,max_halvings=args.max_halvings,
        requested_macro_dt_s=meta['dt_s'],phases=phases,
        initial_snapshot=snapshot(meta,parent),first_decision=decision(meta,meta['dt_s'],args.wall_seconds,phases,args.minimum_dt),
        python=str(ROOT/'.envs/pysph/python.exe'),stop_on_first_shoulder=True,
        accepted=False,formal_accepted_time_advanced_s=0.,no_physical_module_imported=True),meta


def supervise(protocol,meta):
    began=time.perf_counter(); deadline=began+protocol['wall_limit_s']
    started=datetime.now(timezone.utc); output=ROOT/protocol['output']; output.mkdir(parents=True,exist_ok=False)
    protocol.update(started_at=started.astimezone(LOCAL).isoformat(timespec='seconds'),
        latest_stop_at=(started+timedelta(seconds=protocol['wall_limit_s'])).astimezone(LOCAL).isoformat(timespec='seconds'),
        supervisor_pid=os.getpid(),sleep_prevention='PROCESS_LOCAL_SYSTEM_REQUIRED_DISPLAY_MAY_SLEEP')
    atomic_json(output/'protocol.json',protocol)
    parent=ROOT/protocol['initial_parent']; desired=protocol['requested_macro_dt_s']; completed=0; attempt=0; halvings=0
    progress=dict(status='STARTING',finished=False,completed_stage_count=0,attempt_count=0,
        started_at=protocol['started_at'],latest_stop_at=protocol['latest_stop_at'],supervisor_pid=os.getpid(),
        accepted=False,formal_accepted_time_advanced_s=0.,**protocol['initial_snapshot'])
    child_process=None; code=0;kernel=None;mutex=None
    def update(**fields):
        progress.update(fields); progress.update(updated_at=wall_now(),wall_seconds=time.perf_counter()-began,
            remaining_wall_seconds=max(0.,deadline-time.perf_counter()))
        atomic_json(output/'progress.json',progress)
    update()
    try:
        kernel,mutex=singleton_and_awake()
        while True:
            verify(protocol['source_sha256'])
            if (output/'STOP_REQUESTED').exists():
                update(status='STOPPED_USER_REQUEST');break
            if meta['observer']['supplied_material']['crossed_under_annular_shoulder']['mass_kg']>0:
                update(status='STOPPED_FIRST_SCREENED_NEW_MATERIAL_UNDER_SHOULDER');break
            action=decision(meta,desired,deadline-time.perf_counter(),protocol['phases'],protocol['minimum_dt_s'])
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
                    if time.perf_counter()>=min(deadline,launched+650.) or requested_stop:
                        timed_out=not requested_stop;termination=terminate_owned_tree(child_process,stream)
                        child_process.wait(timeout=20.);break
                    wait_seconds=max(0.,min(5.,deadline-time.perf_counter(),launched+650.-time.perf_counter()))
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
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True);parser.add_argument('--parent',default=rel(INITIAL))
    parser.add_argument('--wall-seconds',type=float,default=36000.)
    parser.add_argument('--minimum-dt',type=float,default=1e-7);parser.add_argument('--max-halvings',type=int,default=8)
    parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    protocol,meta=plan(args)
    if args.preflight_only:
        print(json.dumps({k:protocol[k] for k in ('initial_parent','wall_limit_s','first_decision','initial_snapshot')}));return 0
    return supervise(protocol,meta)


if __name__=='__main__':
    raise SystemExit(main())

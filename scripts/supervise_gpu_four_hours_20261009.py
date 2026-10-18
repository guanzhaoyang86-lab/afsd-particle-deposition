"""New, explicitly authorized four-hour administrative GPU continuation.

The copied child loop retains the frozen supervision decisions and receipt
checks. Only a complete screened parent can advance. This module imports no
physical solver, NumPy, CuPy, PySPH, or proposed numerical bridge.
"""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
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

CONTRACT = ROOT/'docs/GPU_FOUR_HOUR_CONTINUATION_20261009.md'
TEST = ROOT/'tests/review/test_gpu_four_hour_supervisor_20261009.py'
CHECKS = BASE/'four_hour_supervisor_components_20261009_v1.json'
CACHED_RUNNER = ROOT/'scripts/run_gpu_cached_preflight_20261008.py'
BRIDGE_RUNNER = ROOT/'scripts/run_gpu_search_tensor_20261009.py'
ALLOWED_RUNNERS = (CACHED_RUNNER, BRIDGE_RUNNER)
MAXIMUM_WALL_SECONDS = 14400.
BRIDGE_STATUS = 'PASS_SEARCH_TENSOR_GPU_REPLAY_AND_RECOVERY'


def utc_now():
    return datetime.now(timezone.utc)


def aware_time(value):
    result = datetime.fromisoformat(value)
    require(result.tzinfo is not None and result.utcoffset() is not None,
        'timezone-aware campaign time required')
    return result.astimezone(timezone.utc)


def bounded_wall_seconds(value):
    require(isinstance(value, (int, float)) and not isinstance(value, bool)
        and math.isfinite(value) and 0. < value <= MAXIMUM_WALL_SECONDS,
        'new authorization permits at most four hours')
    return float(value)


def initial_dt(parent_dt, requested=None):
    chosen = parent_dt if requested is None else requested
    require(all(isinstance(v, (int, float)) and not isinstance(v, bool)
        and math.isfinite(v) for v in (parent_dt, chosen))
        and 0. < chosen <= parent_dt <= .005,
        'initial dt must be positive and cannot grow beyond the complete parent dt')
    return float(chosen)


def fresh_deadline(wall_limit_s, started=None, began=None):
    limit = bounded_wall_seconds(wall_limit_s)
    started = utc_now() if started is None else started
    require(started.tzinfo is not None and started.utcoffset() is not None,
        'aware new campaign start required')
    began = time.perf_counter() if began is None else began
    require(math.isfinite(began), 'finite monotonic start required')
    stop = started.astimezone(timezone.utc)+timedelta(seconds=limit)
    return dict(started_at=started.astimezone(LOCAL).isoformat(timespec='microseconds'),
        latest_stop_at=stop.astimezone(LOCAL).isoformat(timespec='microseconds'),
        monotonic_deadline=began+limit)


def remaining_seconds(latest_stop_at, now=None):
    stop = aware_time(latest_stop_at)
    now = utc_now() if now is None else now
    require(now.tzinfo is not None and now.utcoffset() is not None,
        'aware current time required')
    return max(0., (stop-now.astimezone(timezone.utc)).total_seconds())


def budget_remaining(protocol, deadline):
    # Neither a backward wall-clock adjustment nor a forward one can extend
    # the original monotonic or absolute deadline established for this run.
    return max(0., min(deadline-time.perf_counter(),
        remaining_seconds(protocol['latest_stop_at'])))


def merge_sources(*mappings):
    result = {}
    for mapping in mappings:
        require(isinstance(mapping, dict), 'source mapping required')
        for name, value in mapping.items():
            require(isinstance(name, str) and isinstance(value, str)
                and len(value) == 64 and all(c in '0123456789abcdef' for c in value),
                'canonical SHA256 source mapping required')
            require(name not in result or result[name] == value,
                'conflicting frozen source hashes: '+name)
            result[name] = value
    return result


def runner_sources(runner, qualification, meta):
    runner = inside(Path(runner))
    require(runner in tuple(p.resolve() for p in ALLOWED_RUNNERS),
        'runner is outside the registered cached/qualified bridge set')
    require(runner.is_file(), 'registered runner must already exist')
    parent_sources = meta['source_sha256']
    if runner == CACHED_RUNNER.resolve():
        require(qualification is None, 'cached runner does not take bridge qualification')
        return parent_sources, {rel(runner): sha(runner)}, None
    require(qualification is not None, 'new bridge requires a complete independent qualification')
    path = inside(Path(qualification)); receipt_hash = sha(path); receipt = read(path)
    require(receipt['status'] == BRIDGE_STATUS and receipt['exit_code'] == 0,
        'bridge qualification did not pass')
    require(inside(Path(receipt['runner'])) == runner and receipt['runner_sha256'] == sha(runner),
        'qualified runner binding differs')
    child = receipt['child_source_sha256']
    require(isinstance(child, dict) and child, 'qualified nonempty child source mapping required')
    require(parent_sources == child,
        'four-hour parent must already use the exact qualified child implementation')
    verify(receipt['parent_source_sha256']); verify(receipt['source_sha256'])
    verify(receipt['evidence_sha256']); verify(child)
    require(sha(path) == receipt_hash, 'bridge qualification changed during inspection')
    evidence = merge_sources(receipt['parent_source_sha256'], receipt['source_sha256'], receipt['evidence_sha256'],
        {rel(path): receipt_hash, rel(runner): sha(runner)})
    # A new file may add an explicit alternative, but old frozen files are not
    # silently rewritten to give the bridge different bytes under the same name.
    merge_sources(parent_sources, child, evidence)
    return child, evidence, dict(path=rel(path), sha256=receipt_hash)


def validate_handoff(protocol, meta):
    parent = inside(Path(protocol['initial_parent']))
    fresh = verify_parent(parent, protocol['parent_source_sha256'])
    require(sha(parent/'manifest.json') == protocol['initial_parent_manifest_sha256']
        and fresh['screened_continuation_sha256'] == meta['screened_continuation_sha256']
        and snapshot(fresh, parent) == protocol['initial_snapshot'],
        'complete parent changed after planning')
    qualification = protocol['runner_qualification']
    child, evidence, checked = runner_sources(protocol['runner'],
        qualification['path'] if qualification is not None else None, fresh)
    require(child == protocol['child_source_sha256'] and checked == qualification
        and evidence == protocol['runner_evidence_sha256'],
        'runner qualification or implementation changed after planning')
    require(initial_dt(fresh['dt_s'], protocol['requested_macro_dt_s'])
        == protocol['requested_macro_dt_s'], 'planned initial dt changed')
    bounded_wall_seconds(protocol['wall_limit_s'])


def plan(args):
    output, parent = inside(Path(args.output)), inside(Path(args.parent))
    require(output.is_relative_to(BASE) and output != BASE and not output.exists(),
        'new GPU four-hour campaign folder required')
    limit = bounded_wall_seconds(args.wall_seconds)
    require(math.isfinite(args.minimum_dt) and 0. < args.minimum_dt <= 1e-7
        and isinstance(args.max_halvings, int) and not isinstance(args.max_halvings, bool)
        and 0 <= args.max_halvings <= 8, 'registered temporal refinement bounds required')
    checks = read(CHECKS)
    require(checks['status'] == 'PASS' and checks['exit_code'] == 0,
        'executed four-hour administrative components required')
    verify(checks['source_sha256']); verify(checks['artifacts_sha256'])
    meta = verify_parent(parent)
    desired = initial_dt(meta['dt_s'], args.initial_dt)
    require(desired >= args.minimum_dt, 'initial dt is below the registered minimum')
    runner = inside(Path(args.runner))
    child, evidence, qualification = runner_sources(runner, args.runner_qualification, meta)
    config = ROOT/'configs/deposition_fill_then_track_20261007.json'
    phases = phase_table(read(config))
    source = merge_sources(meta['source_sha256'], child, evidence,
        {rel(p): sha(p) for p in (Path(__file__).resolve(), CONTRACT, TEST, CHECKS, config,
            ROOT/'scripts/supervise_gpu_campaign_20261008.py',
            ROOT/'scripts/supervise_afsd_deposition_20261007.py')})
    return dict(schema='GPU_FOUR_HOUR_CONTINUATION_SUPERVISOR_20261009',
        output=rel(output), initial_parent=rel(parent), runner=rel(runner),
        runner_qualification=qualification, runner_evidence_sha256=evidence,
        source_sha256=source, parent_source_sha256=meta['source_sha256'], child_source_sha256=child,
        initial_parent_manifest_sha256=sha(parent/'manifest.json'), wall_limit_s=limit,
        minimum_dt_s=args.minimum_dt, max_halvings=args.max_halvings,
        requested_macro_dt_s=desired, phases=phases, initial_snapshot=snapshot(meta, parent),
        first_decision=decision(meta, desired, limit, phases, args.minimum_dt),
        python=str(ROOT/'.envs/pysph/python.exe'), stop_on_first_shoulder=True,
        failed_candidate_used=False, accepted=False, formal_accepted_time_advanced_s=0.,
        no_physical_module_imported=True), meta


def supervise(protocol,meta):
    verify(protocol['source_sha256']); validate_handoff(protocol,meta)
    began=time.perf_counter(); timing=fresh_deadline(protocol['wall_limit_s'], began=began)
    deadline=timing.pop('monotonic_deadline')
    output=ROOT/protocol['output']; output.mkdir(parents=True,exist_ok=False)
    protocol.update(**timing, supervisor_pid=os.getpid(),
        sleep_prevention='PROCESS_LOCAL_SYSTEM_REQUIRED_DISPLAY_MAY_SLEEP')
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
            command=[protocol['python'],'-B',str(ROOT/protocol['runner']),'--mode','continue','--output',rel(child),
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
    parser.add_argument('--runner', required=True)
    parser.add_argument('--runner-qualification')
    parser.add_argument('--initial-dt', type=float)
    parser.add_argument('--wall-seconds', type=float, default=14400.)
    parser.add_argument('--minimum-dt', type=float, default=1e-7)
    parser.add_argument('--max-halvings', type=int, default=8)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    protocol, meta = plan(args)
    if args.preflight_only:
        print(json.dumps({k: protocol[k] for k in ('initial_parent', 'runner',
            'wall_limit_s', 'requested_macro_dt_s', 'first_decision', 'initial_snapshot',
            'runner_qualification')}))
        return 0
    return supervise(protocol, meta)


if __name__ == '__main__':
    raise SystemExit(main())

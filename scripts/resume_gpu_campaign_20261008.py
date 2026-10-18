"""One explicitly authorized interruption resume, retaining the old deadline.

The supervisor loop is a verbatim copy of the frozen fast supervisor.  Only
the origin check is different: no incomplete Newton candidate becomes a parent.
"""
import ast
import inspect
import textwrap
import supervise_gpu_campaign_fast_20261008 as frozen
from supervise_gpu_campaign_fast_20261008 import (
    Path, argparse, json, math, os, subprocess, sys, time, traceback,
    require, read, sha, inside, rel, atomic_json, verify, launch_options,
    terminate_owned_tree, decision, verify_parent, snapshot, classify_failure,
    singleton_and_awake, phase_table, wall_now, utc_now, remaining_seconds,
    budget_remaining, aware_time, process_is_running, ROOT, BASE, PASS, THREADS,
    LOCAL, RUNNER, AUTHORIZED_STOP_AT, ORIGIN as ORIGINAL_ORIGIN,
)

INTERRUPTED = BASE/'campaign_10h_20261008_fast_v1'
PARENT = INTERRUPTED/'attempt_0003'
UNFINISHED = INTERRUPTED/'attempt_0004'
PARENT_SHA256 = '5265b7ed9245a97198a3f482951d02ddc8b718e1bd6dc5f65d5f23541d705831'
CONTRACT = ROOT/'docs/GPU_INTERRUPTION_RESUME_20261008.md'
TEST = ROOT/'tests/review/test_gpu_resume_20261008.py'
CHECKS = BASE/'resume_supervisor_components_20261008_v1.json'


def static_loop_audit():
    current = ast.parse(textwrap.dedent(inspect.getsource(supervise)))
    prior = ast.parse(textwrap.dedent(inspect.getsource(frozen.supervise)))
    require(ast.dump(current) == ast.dump(prior), 'supervision loop must remain AST-identical')
    return dict(passed=True, loop_ast_identical=True, physical_solver_imported=False)


def read_interruption(origin, parent):
    origin, parent = inside(Path(origin)), inside(Path(parent))
    require(origin == INTERRUPTED.resolve() and parent == PARENT.resolve(),
        'this authorization names one interrupted campaign and its exact stage-3 PASS')
    original_path = ORIGINAL_ORIGIN/'protocol.json'
    original, prior = read(original_path), read(origin/'protocol.json')
    progress = read(origin/'progress.json')
    verify(original['source_sha256']); verify(prior['source_sha256'])
    limit = float(original['wall_limit_s'])
    require(original['schema'] == 'GPU_TEN_HOUR_SUPERVISOR_20261008'
        and original['output'] == rel(ORIGINAL_ORIGIN)
        and math.isfinite(limit) and 0 < limit <= 36000., 'original authorization identity or budget differs')
    start, stop = aware_time(original['started_at']), aware_time(original['latest_stop_at'])
    require(stop == aware_time(AUTHORIZED_STOP_AT) and 0 < (stop-start).total_seconds() <= limit,
        'original absolute deadline cannot be extended')
    require(prior['schema'] == 'GPU_FIXED_DEADLINE_FAST_SUPERVISOR_20261008'
        and prior['output'] == rel(origin) and prior['origin_campaign'] == rel(ORIGINAL_ORIGIN)
        and prior['latest_stop_at'] == progress['latest_stop_at'] == original['latest_stop_at']
        and prior['child_source_sha256'] == original['child_source_sha256'], 'interrupted campaign lineage differs')
    require(progress.get('finished') is False and progress['status'] == 'RUNNING_GPU_CHILD'
        and progress.get('exit_code') is None and progress['completed_stage_count'] == 3
        and progress['attempt_count'] == 4 and progress['latest_pass_stage'] == rel(parent)
        and progress['current_parent'] == rel(parent) and progress['current_attempt'] == rel(UNFINISHED),
        'only this missing-completion interruption is eligible')
    require(progress['supervisor_pid'] == prior['supervisor_pid'] == 11888 and progress['child_pid'] == 51152,
        'recorded interrupted process identity changed')
    for pid in (progress['supervisor_pid'], progress['child_pid']):
        require(not process_is_running(pid), 'interrupted supervisor or child is still running')
    for path in (origin/'STOP_REQUESTED', origin/'STOP_AFTER_STAGE', origin/'HARD_TIMEOUT', origin/'HARD_TIMEOUT.json',
            origin/'final_summary.json', origin/'attempt_0004_receipt.json',
            UNFINISHED/'manifest.json', UNFINISHED/'supervision.json', UNFINISHED/'pair_screen.json',
            UNFINISHED/'HARD_TIMEOUT', UNFINISHED/'HARD_TIMEOUT.json', UNFINISHED/'STOP_REQUESTED'):
        require(not path.exists(), 'completed, failed, timed-out or explicitly stopped evidence requires separate review: '+rel(path))
    require(UNFINISHED.is_dir() and not (UNFINISHED/'fine_1').exists() and not (UNFINISHED/'fine_2').exists(),
        'registered unfinished coarse attempt changed')
    interrupted_protocol = read(UNFINISHED/'protocol.json')
    require(interrupted_protocol['source_sha256'] == prior['child_source_sha256']
        and interrupted_protocol['parent_gpu_directory'] == rel(parent)
        and interrupted_protocol['parent_sha256'] == PARENT_SHA256,
        'unfinished attempt source or original PASS parent changed')
    verify(interrupted_protocol['source_sha256'])
    if 'administrative_source_sha256' in interrupted_protocol:
        verify(interrupted_protocol['administrative_source_sha256'])
    evidence_paths = [original_path, origin/'protocol.json', origin/'progress.json',
        origin/'attempt_0003_receipt.json', origin/'attempt_0003_launcher.log', origin/'attempt_0004_launcher.log']
    evidence_paths.extend(sorted(p for p in UNFINISHED.rglob('*') if p.is_file()))
    for log in (origin/'attempt_0004_launcher.log', UNFINISHED/'run.log'):
        text = log.read_text(encoding='utf-8')
        require(not any(token in text for token in ('HARD_TIMEOUT', 'EXTERNAL_TIMEOUT_NOT_ACCEPTED', 'FAIL_GPU_STAGE_NOT_ACCEPTED',
            'STOPPED_CHILD_FAILURE', 'STOPPED_WALL_DEADLINE_OR_CHILD_TIMEOUT')),
            'recorded terminal failure is not a missing-completion interruption')
    evidence = {rel(p): sha(p) for p in evidence_paths}
    meta = verify_parent(parent, prior['child_source_sha256'], origin/'attempt_0002')
    receipt = read(origin/'attempt_0003_receipt.json')
    require(receipt['verified'] is True and receipt['exit_code'] == 0 and receipt['timeout'] is False
        and receipt['termination'] is None and receipt['stage'] == rel(parent)
        and receipt['manifest_sha256'] == sha(parent/'manifest.json')
        and receipt['log_sha256'] == sha(origin/'attempt_0003_launcher.log'), 'complete verified outer PASS receipt required')
    expected = snapshot(meta, parent)
    require(all(progress.get(k) == v and receipt.get(k) == v for k, v in expected.items()),
        'interrupted progress and outer receipt must identify the exact latest complete PASS')
    require(meta['screened_continuation_sha256'] == PARENT_SHA256
        and meta['screened_end_step'] == 156, 'registered stage-3 checkpoint identity differs')
    verify(evidence)
    return original, prior, meta, receipt, evidence


def validate_handoff(protocol, meta):
    _, prior, fresh, _, evidence = read_interruption(protocol['origin_campaign'], protocol['initial_parent'])
    require(evidence == protocol['interruption_sha256'], 'interruption evidence changed after planning')
    require(prior['child_source_sha256'] == protocol['child_source_sha256'] == meta['source_sha256']
        and fresh['screened_continuation_sha256'] == meta['screened_continuation_sha256'] == PARENT_SHA256
        and sha(PARENT/'manifest.json') == protocol['initial_parent_manifest_sha256'],
        'frozen parent or implementation changed after planning')


def plan(args):
    output = inside(Path(args.output))
    require(output.is_relative_to(BASE) and output != BASE and not output.exists(), 'new GPU resume folder required')
    require(0 < args.minimum_dt <= 1e-7 and 0 <= args.max_halvings <= 8, 'registered temporal refinement bounds required')
    checks = read(CHECKS)
    require(checks['status'] == 'PASS' and checks['exit_code'] == 0, 'executed resume components required')
    verify(checks['source_sha256']); verify(checks['artifacts_sha256']); static_loop_audit()
    original, prior, meta, receipt, evidence = read_interruption(args.origin_campaign, args.parent)
    phases = phase_table(read(ROOT/'configs/deposition_fill_then_track_20261007.json'))
    require(phases == prior['phases'] == original['phases'], 'registered phase program changed')
    source = dict(prior['source_sha256'])
    source.update({rel(p): sha(p) for p in (Path(__file__).resolve(), CONTRACT, TEST, CHECKS,
        Path(frozen.__file__).resolve(), ROOT/'tests/review/test_gpu_fast_supervisor_20261008.py',
        ROOT/'scripts/start_gpu_resume_detached_20261008.py')})
    source.update(evidence)
    remaining = remaining_seconds(original['latest_stop_at'])
    action = decision(meta, meta['dt_s'], remaining, phases, args.minimum_dt)
    return dict(schema='GPU_INTERRUPTION_RESUME_SUPERVISOR_20261008', output=rel(output),
        origin_campaign=rel(INTERRUPTED), initial_parent=rel(PARENT), interruption_sha256=evidence,
        source_sha256=source, child_source_sha256=meta['source_sha256'],
        initial_parent_manifest_sha256=sha(PARENT/'manifest.json'),
        latest_stop_at=original['latest_stop_at'], original_started_at=original['started_at'],
        original_wall_limit_s=original['wall_limit_s'], wall_limit_s=remaining,
        minimum_dt_s=args.minimum_dt, max_halvings=args.max_halvings,
        requested_macro_dt_s=meta['dt_s'], phases=phases, initial_snapshot=snapshot(meta,PARENT), first_decision=action,
        measured_cost=dict(maximum_candidate_wall_s=meta['maximum_candidate_wall_s'],
            previous_complete_stage_wall_s=receipt['wall_seconds'], forecast_s=action['forecast_s']),
        python=str(ROOT/'.envs/pysph/python.exe'), stop_on_first_shoulder=True,
        interrupted_candidate_used=False, interruption_cause='UNKNOWN_MISSING_COMPLETION_NOT_PHYSICAL_PASS',
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
    parser.add_argument('--parent', default=rel(PARENT))
    parser.add_argument('--origin-campaign', default=rel(INTERRUPTED))
    parser.add_argument('--minimum-dt', type=float, default=1e-7)
    parser.add_argument('--max-halvings', type=int, default=8)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    protocol, meta = plan(args)
    if args.preflight_only:
        print(json.dumps({k: protocol[k] for k in ('initial_parent', 'latest_stop_at', 'wall_limit_s',
            'first_decision', 'measured_cost', 'initial_snapshot')}))
        return 0
    return supervise(protocol, meta)


if __name__ == '__main__':
    raise SystemExit(main())


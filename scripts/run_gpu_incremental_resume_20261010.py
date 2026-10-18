"""Locked-history session, verified step recovery, unchanged GPU physical pairs."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
import ast
import copy
import hashlib
import inspect
import json
import multiprocessing
import os
import sys
import time
import traceback
import numpy as np
import run_gpu_search_tensor_20261009 as bridge
import run_gpu_cached_preflight_20261008 as cache
import supervise_gpu_campaign_20261008 as admin
from afsd_locked_history_20261010 import LockedHistory
import afsd_adaptive_dt_20261010 as adaptive
import afsd_deposition_progress_20261010 as deposition

ROOT, BASE, CPU_BASE = bridge.ROOT, bridge.BASE, bridge.CPU_BASE
require, read, load, rel = bridge.require, bridge.read, bridge.load, bridge.rel
CONTINUE_PASS, REPLAY_PASS = bridge.CONTINUE_PASS, bridge.REPLAY_PASS
frozen_sources, adaptation_audit = bridge.frozen_sources, bridge.adaptation_audit
PROTOCOL = bridge.PROTOCOL
CONTRACT = ROOT/'docs/GPU_INCREMENTAL_RESUME_20261010.md'
CHECKS = bridge.CHECKS
INCREMENTAL_CHECKS = BASE/'incremental_resume_components_20261010_v1.json'
SESSION_PROTOCOL = 'LIVE_WINDOWS_LOCKED_HISTORY_RESUME_BOUNDED_GRANT_20261010'


def merge(*mappings):
    result = {}
    for mapping in mappings:
        for name, digest in mapping.items():
            require(name not in result or result[name] == digest, 'conflicting evidence: '+name)
            result[name] = digest
    return result


class LockedPreflight(cache.CachedPreflight):
    """All old receipt/history assertions, protected by live immutable handles."""
    def __init__(self, pool):
        super().__init__(); self.pool = pool

    def fresh_hash(self, path):
        digest = self.pool.sha(path)
        return digest, self.signature(path)

    def finish(self):
        self.ensure_open()
        self.pool.verify({rel(row['path']): row['digest'] for row in self.records.values()})
        self.finished = True
        return dict(status='PASS_LIVE_LOCKED_HISTORY_TRANSACTION', wall_seconds=time.perf_counter()-self.started,
            unique_files=len(self.records), hash_requests=self.requests,
            parsed_json_files=len(self.json_cache), decoded_checkpoints=len(self.load_cache),
            persistent_cache_requires_live_handles=True, stat_only_trust=False)


def checks(pool):
    receipt = read(INCREMENTAL_CHECKS)
    require(receipt['status'] == 'PASS_INCREMENTAL_COMPONENTS' and receipt['exit_code'] == 0,
        'executed incremental components required')
    pool.seal(receipt['source_sha256']); pool.seal(receipt['artifacts_sha256'])
    required = ['scripts/run_gpu_incremental_resume_20261010.py', 'scripts/afsd_incremental_worker_resume_20261010.py',
        'scripts/afsd_locked_history_20261010.py', 'scripts/afsd_adaptive_dt_20261010.py',
        'scripts/afsd_deposition_progress_20261010.py', rel(CONTRACT)]
    require(set(required).issubset(receipt['source_sha256']), 'incomplete administrative component binding')
    return merge(receipt['source_sha256'], receipt['artifacts_sha256'], {rel(INCREMENTAL_CHECKS):pool.sha(INCREMENTAL_CHECKS)})


def static_audit():
    left = ast.parse(inspect.getsource(bridge.preflight)); left.body[0].name = '_initial_preflight'
    left.body[0].args.args.append(ast.arg(arg='pool'))
    replaced = 0
    for node in ast.walk(left):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'CachedPreflight':
            require(isinstance(node.func.value, ast.Name) and node.func.value.id == 'cache', 'unexpected transaction provider')
            node.func = ast.Name(id='LockedPreflight', ctx=ast.Load()); node.args = [ast.Name(id='pool', ctx=ast.Load())]; replaced += 1
    require(replaced == 1 and ast.dump(left) == ast.dump(ast.parse(inspect.getsource(_initial_preflight))),
        'initial original history assertions changed')
    return dict(status='PASS_UNCHANGED_INITIAL_RECEIPT_ASSERTIONS_AND_PHYSICAL_PAIR',
        transaction_provider_replacements=replaced, physical_pair='run_gpu_search_tensor_20261009._run_pair_original',
        inherited_physical_audit=bridge.adaptation_audit(), runtime_global_rebinding=False)


def stage_protocol(template, meta, parent, dt, administrative, session_record, refs):
    require(meta['source_sha256'] == template['source_sha256'], 'physical source differs from qualified session')
    merge(template['source_sha256'], refs, administrative, session_record)
    p = copy.deepcopy(template)
    for name in ('search_tensor_reference_stage', 'administrative_preflight_statistics'):
        p.pop(name, None)
    p.update(mode='continue', parent_checkpoint=meta['screened_continuation_checkpoint'],
        parent_sha256=meta['screened_continuation_sha256'], parent_gpu_directory=rel(parent),
        parent_implementation_source_sha256=meta['source_sha256'], dt_s=dt,
        candidate_forecast_s=meta['maximum_candidate_wall_s'], forecast_s=30.+4.5*meta['maximum_candidate_wall_s'],
        reference_sha256=merge(refs, administrative, session_record), administrative_source_sha256=administrative,
        administrative_protocol=SESSION_PROTOCOL)
    require(p['forecast_s'] <= 450., 'unchanged candidate forecast limit exceeded')
    return p


def authorize_grant(request, protocol_sha, phase, pool, protocol):
    require(phase in ('start','end'), 'only two ordered live grants are permitted')
    require(request == dict(protocol_sha256=protocol_sha, phase=phase), 'worker grant identity or order mismatch')
    pool.verify(protocol['source_sha256']); pool.verify(protocol['reference_sha256'])
    require(pool.sha(ROOT/protocol['parent_checkpoint']) == protocol['parent_sha256'], 'locked parent changed')
    return dict(status='PASS_LIVE_HISTORY_GRANT', phase=phase, protocol_sha256=protocol_sha,
        controller_pid=os.getpid(), observed_at=admin.wall_now(), history_locked=True)


def run_stage(pool, output, protocol, controller_deadline, update, stop_file):
    from afsd_incremental_worker_resume_20261010 import stage_worker
    output.mkdir(); bridge.shared.write_new(output/'protocol.json', protocol)
    protocol_sha = pool.sha(output/'protocol.json')
    pool.seal(protocol['source_sha256']); pool.seal(protocol['reference_sha256'])
    context = multiprocessing.get_context('spawn'); server, child_channel = context.Pipe(duplex=True)
    child = context.Process(target=stage_worker, args=(str(output), child_channel, protocol_sha), daemon=False)
    started = time.perf_counter(); history_before = pool.stats()
    child.start(); child_channel.close(); update(child_pid=child.pid)
    grants = []; interrupted = False; error = None; last_update = started
    try:
        while child.is_alive():
            if time.perf_counter() >= min(started+600., controller_deadline) or stop_file.exists():
                interrupted = True
                with (output/'termination.log').open('a', encoding='utf-8') as stream:
                    import subprocess
                    subprocess.run(['taskkill','/PID',str(child.pid),'/T','/F'], stdout=stream, stderr=subprocess.STDOUT, timeout=20.)
                child.join(timeout=20.); break
            if len(grants) < 2 and server.poll(.25):
                request = server.recv()
                phase = ('start','end')[len(grants)] if len(grants) < 2 else None
                response = authorize_grant(request, protocol_sha, phase, pool, protocol)
                grants.append(response); server.send(response)
            elif len(grants) == 2:
                child.join(timeout=.25)
            if time.perf_counter()-last_update >= 5.:
                update(); last_update = time.perf_counter()
        child.join(timeout=5.)
        code = child.exitcode if child.exitcode is not None else 124
        require(not child.is_alive(), 'worker could not be terminated')
        pool.verify(protocol['source_sha256']); pool.verify(protocol['reference_sha256'])
        require(code != 0 or len(grants) == 2, 'successful worker lacks both live grants')
    except Exception as exc:
        error = repr(exc); code = 2
        if child.is_alive():
            child.terminate(); child.join(timeout=15.)
    finally:
        server.close()
    manifest = output/'manifest.json'; log = output/'run.log'
    bridge.shared.write_new(output/'supervision.json', dict(exit_code=code, wall_seconds=time.perf_counter()-started,
        source_inputs_still_match=error is None and len(grants)==2, source_error=error, live_history_grants=grants,
        manifest_sha256=bridge.core.sha(manifest) if manifest.exists() else None,
        log_sha256=bridge.core.sha(log) if log.exists() else None, interrupted=interrupted,
        locked_history_before=history_before, locked_history_after=pool.stats(),
        accepted=False, restart_allowed=False, administrative_protocol=SESSION_PROTOCOL))
    return code, time.perf_counter()-started


def validated_stage(pool, folder, source, parent_meta):
    transaction = LockedPreflight(pool)
    m = transaction.gpu_receipt(folder, {CONTINUE_PASS})
    require(m['source_sha256'] == source and m['parent_sha256'] == parent_meta['screened_continuation_sha256'],
        'new stage source or exact parent mismatch')
    require(m['screened_end_step'] == parent_meta['screened_end_step']+2, 'fine state step continuity failed')
    transaction.finish()
    return m


def bind_stage(pool, folder, meta, refs):
    extra = dict(meta['artifacts_sha256'])
    for name in ('manifest.json','supervision.json','run.log'):
        path = folder/name; extra[rel(path)] = pool.sha(path)
    pool.seal(extra)
    return merge(refs, extra)


def observe(meta, frames):
    if frames:
        targets = [(row['checkpoint'], read(ROOT/row['observer_path'])) for row in meta['screened_history']]
    else:
        targets = [(meta['screened_continuation_checkpoint'], meta['observer'])]
    for checkpoint, observer in targets:
        pa, state = load(ROOT/checkpoint)
        frame = deposition.frame_from_saved(pa, state, None, observer)
        frames.append(frame)
    return deposition.summarize_progress(frames, state['parameters']), state


def run(args):
    output = (ROOT/args.output).resolve(); parent = (ROOT/args.parent).resolve()
    require(output.is_relative_to(BASE) and not output.exists(), 'new project-local session output required')
    require(0 < args.wall_seconds <= 14400. and 1 <= args.max_stages <= 2000, 'bounded authorized validation campaign required')
    require(0 < args.maximum_dt <= .005 and 0 < args.shoulder_observe_ms <= 20., 'registered recovery/observation bounds required')
    require(0 < args.target_time <= .22, 'stationary filling only')
    output.mkdir(parents=True); began = time.perf_counter(); started = datetime.now(timezone.utc)
    deadline = began+args.wall_seconds; code = 0; kernel = mutex = None; frames = []
    progress = dict(status='INITIAL_FULL_HISTORY_VERIFICATION', finished=False, started_at=started.astimezone(admin.LOCAL).isoformat(),
        latest_stop_at=(started+timedelta(seconds=args.wall_seconds)).astimezone(admin.LOCAL).isoformat(),
        supervisor_pid=os.getpid(), child_pid=None, completed_stage_count=0, attempt_count=0,
        accepted=False, formal_accepted_time_advanced_s=0.)
    def update(**fields):
        progress.update(fields, updated_at=admin.wall_now(), wall_seconds=time.perf_counter()-began,
            remaining_wall_seconds=max(0.,deadline-time.perf_counter()))
        admin.atomic_json(output/'progress.json',progress)
    update()
    pool = LockedHistory(ROOT)
    try:
        kernel, mutex = admin.singleton_and_awake()
        administrative = checks(pool); audit = static_audit()
        startup_args = argparse.Namespace(mode='continue', parent=rel(parent), dt=None,
            baseline=rel(CPU_BASE/'stage_013'), cost=rel(BASE/'benchmark_v1'), validation=rel(BASE/'replay_v1'), reference_stage=None)
        template = _initial_preflight(startup_args, pool)
        transaction = LockedPreflight(pool); meta = transaction.gpu_receipt(parent,{CONTINUE_PASS}); transaction.finish()
        refs = bind_stage(pool,parent,meta,dict(template['reference_sha256']))
        reference = None
        if args.replay_stage:
            reference = (ROOT/args.replay_stage).resolve(); transaction = LockedPreflight(pool)
            rm = transaction.gpu_receipt(reference,{CONTINUE_PASS}); transaction.finish()
            require(rm['parent_sha256'] == meta['screened_continuation_sha256'] and rm['dt_s'] == meta['dt_s']
                and rm['source_sha256'] == meta['source_sha256'], 'administrative replay must match exact physical implementation and parent')
            refs = bind_stage(pool,reference,rm,refs)
        observation, state = observe(meta,frames)
        policy = adaptive.initial_state(meta['screened_continuation_sha256'],meta['dt_s'],maximum_dt_s=args.maximum_dt)
        phases = admin.phase_table(read(ROOT/'configs/deposition_fill_then_track_20261007.json'))
        session = dict(schema=SESSION_PROTOCOL,initial_parent=rel(parent),source_sha256=meta['source_sha256'],
            administrative_source_sha256=administrative,initial_preflight_statistics=template['administrative_preflight_statistics'],
            adapter_audit=audit,wall_limit_s=args.wall_seconds,maximum_dt_s=args.maximum_dt,
            target_time_s=args.target_time,shoulder_observe_ms=args.shoulder_observe_ms,replay_stage=rel(reference) if reference else None,
            started_at=progress['started_at'],latest_stop_at=progress['latest_stop_at'],accepted=False)
        bridge.shared.write_new(output/'protocol.json',session)
        session_record={rel(output/'protocol.json'):pool.sha(output/'protocol.json')}
        shoulder_started = None; shoulder_start_completed = None; completed = 0; attempt = 0
        update(status='INITIAL_VERIFICATION_COMPLETE',initial_preflight_statistics=template['administrative_preflight_statistics'],
            adaptive_state=policy,**admin.snapshot(meta,parent))
        while completed < args.max_stages:
            if (output/'STOP_REQUESTED').exists() or (output/'STOP_AFTER_STAGE').exists():
                update(status='STOPPED_USER_REQUEST'); break
            if meta['screened_end_time_s'] >= args.target_time-1e-14:
                update(status='STOPPED_STATIONARY_TARGET_REACHED'); break
            first_shoulder=observation['first_observed_supplied_under_shoulder_time_s']
            if first_shoulder is not None and shoulder_started is None:
                shoulder_started=first_shoulder; shoulder_start_completed=completed
            if shoulder_started is not None:
                if meta['screened_end_time_s']-shoulder_started >= args.shoulder_observe_ms/1000 and completed-shoulder_start_completed >= 5:
                    update(status='STOPPED_SHOULDER_OBSERVATION_WINDOW_COMPLETE'); break
            desired = policy['nominal_dt_s']
            action = admin.decision(meta,desired,max(0.,deadline-time.perf_counter()),phases,1e-7)
            if action['action'] != 'RUN': update(status=action['action'],next_decision=action); break
            dt = min(action['dt_s'],args.target_time-meta['screened_end_time_s'])
            clipped = dt < desired*(1.-1e-12)
            policy = adaptive.begin_attempt(policy,dt,phase_clipped=clipped)
            attempt += 1; child=output/f'attempt_{attempt:04d}'
            protocol=stage_protocol(template,meta,parent,dt,administrative,session_record,refs)
            if reference is not None and completed == 0: protocol['incremental_replay_stage']=rel(reference)
            update(status='RUNNING_GPU_CHILD',attempt_count=attempt,current_attempt=rel(child),current_dt_s=dt,
                current_parent=rel(parent),adaptive_state=policy)
            exit_code,wall=run_stage(pool,child,protocol,deadline,update,output/'STOP_REQUESTED')
            row=dict(attempt=attempt,stage=rel(child),parent=rel(parent),dt_s=dt,exit_code=exit_code,wall_seconds=wall)
            if exit_code != 0:
                failure=read(child/'manifest.json') if (child/'manifest.json').exists() else {}
                screen=read(child/'pair_screen.json') if (child/'pair_screen.json').exists() else None
                gates=False
                if admin.classify_failure(failure,screen)=='REFINE':
                    pool.seal(failure['source_sha256']); pool.seal(failure['reference_sha256']);pool.seal(failure['artifacts_sha256'])
                    require(failure['parent_sha256']==meta['screened_continuation_sha256'] and failure['source_sha256']==meta['source_sha256'], 'failed trial source/parent mismatch')
                    gates=all(read(child/role/'candidate_qualification.json')['status']=='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN' for role in ('coarse','fine_1','fine_2'))
                policy,event=adaptive.resolve_failure(policy,failure,screen,state['parameters']['dx'],
                    state['bore_observation']['particle_ids'],original_candidate_gates_passed=gates)
                row['adaptive_event']=event; bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json',row)
                if event['action']=='REFINE':
                    update(status='RETRY_SMALLER_DT_FROM_SAME_COMPLETE_PARENT',adaptive_state=policy,last_adaptive_event=event);continue
                update(status='STOPPED_CHILD_FAILURE',error=failure.get('error','worker failed'),failed_attempt=rel(child));code=2;break
            next_meta=validated_stage(pool,child,meta['source_sha256'],meta)
            refs=bind_stage(pool,child,next_meta,refs)
            policy,event=adaptive.resolve_pass(policy,next_meta['pair_screen'],state['parameters']['dx'],
                state['bore_observation']['particle_ids'],next_meta['screened_continuation_sha256'],original_candidate_gates_passed=True)
            parent,meta=child,next_meta;completed+=1
            observation,state=observe(meta,frames); admin.atomic_json(output/'deposition_progress.json',observation)
            row.update(verified=True,adaptive_event=event,locked_history_statistics=pool.stats(),**admin.snapshot(meta,parent))
            bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json',row)
            update(status='VERIFIED_PASS_READY_FOR_NEXT_STAGE',completed_stage_count=completed,child_pid=None,
                last_adaptive_event=event,adaptive_state=policy,**admin.snapshot(meta,parent))
        else: update(status='STOPPED_STAGE_LIMIT')
    except Exception as error:
        code=2;update(status='STOPPED_SESSION_ERROR',error=repr(error),traceback=traceback.format_exc())
    finally:
        try:
            audit=pool.final_audit(); bridge.shared.write_new(output/'final_history_audit.json',audit)
        except Exception as error:
            code=2;update(status='STOPPED_FINAL_HISTORY_VERIFICATION_FAILURE',final_audit_error=repr(error))
        pool.close()
        if kernel is not None:
            kernel.SetThreadExecutionState(0x80000000);kernel.CloseHandle(mutex)
        update(finished=True,exit_code=code,child_pid=None)
        bridge.shared.write_new(output/'final_summary.json',progress)
    print(json.dumps(progress),flush=True);return code


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True);p.add_argument('--parent',required=True)
    p.add_argument('--replay-stage');p.add_argument('--wall-seconds',type=float,default=14400.)
    p.add_argument('--max-stages',type=int,default=1000)
    p.add_argument('--maximum-dt',type=float,default=.000474609375)
    p.add_argument('--target-time',type=float,default=.22)
    p.add_argument('--shoulder-observe-ms',type=float,default=5.)
    p.add_argument('--detach',action='store_true');p.add_argument('--detached-child',action='store_true')
    args=p.parse_args()
    if args.detach or args.detached_child:
        import subprocess
        from start_gpu_resume_detached_20261008 import job_policy, write_new
        output=(ROOT/args.output).resolve()
        require(output.is_relative_to(BASE) and not output.exists(),'new project-local detached output required')
        if args.detached_child:
            policy=job_policy()
            write_new(output.with_name(output.name+'.bootstrap.json'),dict(**policy,supervisor=rel(Path(__file__)),
                status='PASS_DETACHED_JOB_POLICY' if not policy['kill_on_job_close'] else 'FAIL_JOB_POLICY'))
            require(not policy['kill_on_job_close'],'detached session retains kill-on-close')
        else:
            require(os.name=='nt','Windows detached launch required')
            forwarded=[v for v in sys.argv[1:] if v!='--detach']
            command=[sys.executable,'-B',str(Path(__file__).resolve()),*forwarded,'--detached-child']
            startup=subprocess.STARTUPINFO();startup.dwFlags|=subprocess.STARTF_USESHOWWINDOW;startup.wShowWindow=0
            flags=subprocess.CREATE_BREAKAWAY_FROM_JOB|subprocess.DETACHED_PROCESS|subprocess.CREATE_NEW_PROCESS_GROUP
            env=dict(os.environ,**{k:'1' for k in bridge.shared.THREADS});env['PYTHONPATH']=''
            with output.with_name(output.name+'.stdout.log').open('xb',buffering=0) as stream:
                process=subprocess.Popen(command,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=stream,
                    stderr=subprocess.STDOUT,startupinfo=startup,creationflags=flags,close_fds=True)
            receipt=dict(status='STARTED_INCREMENTAL_SESSION_PENDING_BOOTSTRAP',pid=process.pid,command=command,
                creationflags=flags,source_sha256={rel(Path(__file__)):bridge.core.sha(Path(__file__))})
            write_new(output.with_name(output.name+'.launch.json'),receipt);print(json.dumps(receipt));return 0
    return run(args)


# STATIC_INITIAL_PREFLIGHT_FOLLOWS


def _initial_preflight(args, pool):
    require(args.mode == 'continue' and args.parent, 'explicit complete parent required')
    transaction = LockedPreflight(pool)
    # Reuse every original CPU, GPU benchmark and GPU replay verification.
    base_args = copy.copy(args); base_args.parent = None; base_args.dt = None
    protocol = transaction.preflight(base_args)
    legacy = dict(protocol['source_sha256']); current = frozen_sources(transaction, legacy)
    valid_dir = (ROOT/args.validation).resolve()
    valid = transaction.gpu_receipt(valid_dir, {REPLAY_PASS})
    refs = protocol['reference_sha256']; folder = (ROOT/args.parent).resolve()
    chosen = None; visited = set()
    while True:
        require(folder not in visited, 'cycle in explicit search-domain continuation')
        visited.add(folder)
        meta = transaction.gpu_receipt(folder, {CONTINUE_PASS})
        require(meta['source_sha256'] in (legacy, current)
            and meta['validation_directory'] == rel(valid_dir), 'unregistered implementation or replay in parent chain')
        if meta['source_sha256'] == current:
            require(meta.get('search_tensor_protocol') == PROTOCOL, 'new implementation lacks explicit bridge identity')
        if chosen is None:
            chosen = meta
        transaction.references(folder, meta, refs)
        previous = meta.get('parent_gpu_directory')
        earlier = transaction.gpu_receipt(ROOT/previous, {CONTINUE_PASS}) if previous else valid
        require(meta['parent_sha256'] == earlier['screened_continuation_sha256'], 'explicit bridge parent chain mismatch')
        if not previous:
            break
        folder = (ROOT/previous).resolve()
    dt = args.dt if args.dt is not None else chosen['dt_s']
    require(np.isfinite(dt) and 0 < dt <= .005, 'bounded finite positive macro interval required')
    protocol.update(source_sha256=current, legacy_source_sha256=legacy,
        parent_checkpoint=chosen['screened_continuation_checkpoint'], parent_sha256=chosen['screened_continuation_sha256'],
        parent_gpu_directory=rel((ROOT/args.parent).resolve()), parent_implementation_source_sha256=chosen['source_sha256'],
        dt_s=dt, candidate_forecast_s=chosen['maximum_candidate_wall_s'],
        forecast_s=30.+4.5*chosen['maximum_candidate_wall_s'], search_tensor_protocol=PROTOCOL,
        search_tensor_components=rel(CHECKS), search_tensor_adaptation_audit=adaptation_audit())
    require(protocol['forecast_s'] <= 450., 'original measured cost forecast exceeds450s')
    if args.reference_stage:
        reference_dir = (ROOT/args.reference_stage).resolve()
        reference = transaction.gpu_receipt(reference_dir, {CONTINUE_PASS})
        require(reference['source_sha256'] == legacy and reference['parent_sha256'] == protocol['parent_sha256']
            and reference['dt_s'] == dt, 'replay must match the exact original complete parent and interval')
        transaction.references(reference_dir, reference, refs)
        protocol['search_tensor_reference_stage'] = rel(reference_dir)
    cache.bind_administration(transaction, protocol, cache.DEFAULT_CHECKS)
    transaction.verify(refs)
    protocol['administrative_preflight_statistics'] = transaction.finish()
    return protocol

if __name__ == "__main__":
    raise SystemExit(main())

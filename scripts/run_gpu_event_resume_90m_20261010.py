"""Administrative-only 90 minute resume of the fully qualified frozen v5 history."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
import ast
import copy
import inspect
import json
import math
import multiprocessing
import os
import sys
import time
import traceback
import run_gpu_event_delay_campaign_20261010 as previous

ROOT, BASE, CPU_BASE = previous.ROOT, previous.BASE, previous.CPU_BASE
frozen, bridge, admin, adaptive = previous.frozen, previous.bridge, previous.admin, previous.adaptive
require, read, load, rel = previous.require, previous.read, previous.load, previous.rel
source_generation, source_transition = previous.source_generation, previous.source_transition
cost_evidence, keep_failure = previous.cost_evidence, previous.keep_failure
measured_envelope = previous.measured_envelope
authorize_grant = previous.authorize_grant
LOCATOR_PROTOCOL = previous.LOCATOR_PROTOCOL
BLOCKED_PREFIX_PROTOCOL = previous.BLOCKED_PREFIX_PROTOCOL
VISIBLE_PATH_PROTOCOL = previous.VISIBLE_PATH_PROTOCOL
ROLES = previous.ROLES
PROTOCOL = 'FROZEN_V5_COMPLETE_PARENT_ADMINISTRATIVE_RESUME_90M_20261010_V1'
SESSION_PROTOCOL = PROTOCOL
AUTHORIZED_AT = datetime.fromisoformat('2026-10-10T23:44:34+08:00')
ABSOLUTE_STOP = datetime.fromisoformat('2026-10-11T01:14:34+08:00')
FORECAST_LIMIT_S = 1200.
PHYSICAL_SOFT_LIMIT_S = 1260.
WORKER_HARD_LIMIT_S = 1290.
CONTROLLER_LIMIT_S = 1380.
CONTRACT = ROOT/'docs/GPU_EVENT_RESUME_90M_20261010.md'
COMPONENTS = BASE/'event_resume_90m_components_20261010_v1.json'
OLD_CAMPAIGN = BASE/'event_delay_campaign_20261010_v1'
PARENT = OLD_CAMPAIGN/'attempt_0007'
PARENT_SHA256 = 'b1f950f4f5b947ed8d13f2b256f33b7ee7fa9f50b8dbbc5ba2f9b8541d8e9c67'
PARENT_TIME_S = .07517333984374985
PARENT_STEP = 499
ADMINISTRATIVE_FILES = ('scripts/run_gpu_event_resume_90m_20261010.py',
    'scripts/afsd_event_resume_90m_worker_20261010.py', 'scripts/check_gpu_event_resume_90m_20261010.py',
    'tests/review/test_gpu_event_resume_90m_20261010.py', 'docs/GPU_EVENT_RESUME_90M_20261010.md')


def authorized_wall_seconds(requested, now):
    require(type(requested) in (int, float) and math.isfinite(requested) and 0 < requested <= 5400.,
            'requested authorization must be finite and at most 5400 seconds')
    require(now.tzinfo is not None and now >= AUTHORIZED_AT, 'timezone-aware authorization start required')
    remaining = (ABSOLUTE_STOP-now).total_seconds()
    require(remaining > 0, '90 minute absolute authorization expired')
    return min(float(requested), remaining)


def checks(pool):
    inherited = previous.checks(pool)
    receipt = read(COMPONENTS)
    require(receipt['status'] == 'PASS_EVENT_RESUME_90M_COMPONENTS' and receipt['exit_code'] == 0
        and receipt['physical_source_changed'] is False and receipt['gpu_campaign_runs'] == 0,
            'executed administrative-only resume checks required')
    require(receipt['inherited_receipt'] == rel(previous.COMPONENTS)
        and receipt['inherited_receipt_sha256'] == pool.sha(previous.COMPONENTS), 'frozen v5 component receipt required')
    require(set(ADMINISTRATIVE_FILES).issubset(receipt['source_sha256'])
        and all(receipt['source_sha256'].get(k) == v for k, v in inherited.items()),
            'administrative receipt must preserve frozen v5 evidence')
    require(receipt['physical_source_sha256'] == read(PARENT/'manifest.json')['source_sha256']
        and len(receipt['physical_source_sha256']) == 636, 'exact unchanged 636-entry physical map required')
    require(receipt['parent_directory'] == rel(PARENT) and receipt['parent_checkpoint_sha256'] == PARENT_SHA256
        and receipt['old_session_directory'] == rel(OLD_CAMPAIGN) and receipt['absolute_stop'] == ABSOLUTE_STOP.isoformat()
        and receipt['registered_wall_limits_s'] == dict(forecast=FORECAST_LIMIT_S, physical_soft=PHYSICAL_SOFT_LIMIT_S,
            worker_hard=WORKER_HARD_LIMIT_S, controller=CONTROLLER_LIMIT_S), 'registered parent, deadline or budget differs')
    pool.seal(receipt['source_sha256']); pool.seal(receipt['artifacts_sha256'])
    return frozen.merge(inherited, receipt['source_sha256'], receipt['artifacts_sha256'],
                        {rel(COMPONENTS): pool.sha(COMPONENTS)})


def static_audit():
    expected = inspect.getsource(previous.run_stage).replace('afsd_event_delay_worker_20261010',
        'afsd_event_resume_90m_worker_20261010').replace('started+1080.', 'started+1380.')
    require(ast.dump(ast.parse(expected)) == ast.dump(ast.parse(inspect.getsource(run_stage))),
            'supervision differs beyond registered worker and 1380-second bound')
    old = (ROOT/'scripts/afsd_event_delay_worker_20261010.py').read_text(encoding='utf-8')
    new = (ROOT/'scripts/afsd_event_resume_90m_worker_20261010.py').read_text(encoding='utf-8')
    expected = old.replace('threading.Timer(990.,timeout)', 'threading.Timer(1290.,timeout)').replace('began+960.', 'began+1260.')
    require(ast.dump(ast.parse(expected)) == ast.dump(ast.parse(new)), 'worker differs beyond 1260/1290-second limits')
    require(ABSOLUTE_STOP < frozen.ABSOLUTE_STOP and (ABSOLUTE_STOP-AUTHORIZED_AT).total_seconds() == 5400.,
            'new absolute stop must fit the explicit 90 minute request')
    return dict(status='PASS_ADMINISTRATIVE_ONLY_90M_RESUME', physical_source_changed=False,
        frozen_v5_audit=previous.static_audit(), registered_wall_limits_s=dict(forecast=FORECAST_LIMIT_S,
            physical_soft=PHYSICAL_SOFT_LIMIT_S, worker_hard=WORKER_HARD_LIMIT_S, controller=CONTROLLER_LIMIT_S),
        grant_timeout_unchanged_s=180., absolute_stop=ABSOLUTE_STOP.isoformat())


def cost_decision(meta, dt, remaining, phases, rows):
    envelope = measured_envelope(rows)
    t = float(meta['screened_end_time_s']); minimum_dt = 1e-7; forecast = envelope['forecast_s']
    require(all(math.isfinite(v) for v in (t, dt, remaining)) and t >= 0 and minimum_dt <= dt <= .005,
            'finite bounded cost/time inputs required')
    result = dict(forecast_s=forecast, child_timeout_s=max(0., min(CONTROLLER_LIMIT_S, remaining)), event_cost_envelope=envelope)
    if remaining < forecast+20.: return dict(result, action='STOP_WALL_RESERVE')
    if forecast > FORECAST_LIMIT_S: return dict(result, action='STOP_COST_REASSESSMENT')
    phase = next((p for p in phases if p['start'] <= t < p['end']), None)
    if phase is None or not phase['stationary']: return dict(result, action='STOP_REGISTERED_STATIONARY_FILL_END')
    try:
        bounded = admin.bounded_dt(t, dt, phase['end'])
    except ValueError as error:
        return dict(result, action='STOP_UNREPRESENTABLE_PHASE_BOUNDARY', reason=str(error))
    if bounded < minimum_dt: return dict(result, action='STOP_MINIMUM_DT_OR_BOUNDARY')
    return dict(result, action='RUN', dt_s=bounded, phase=phase['name'], phase_end_s=phase['end'])


def stage_protocol(template, meta, parent, dt, administrative, session_record, refs, scale, costs):
    require(scale == .25 and meta['source_sha256'] == template['source_sha256']
        and source_generation(meta, template) == 6 and source_transition(template, meta, parent) is None,
            'same exact qualified v5 implementation and quarter tolerance required')
    envelope = measured_envelope(costs)
    require(envelope['forecast_s'] <= FORECAST_LIMIT_S, 'registered 1200-second forecast limit exceeded')
    frozen.merge(template['source_sha256'], refs, administrative, session_record)
    result = copy.deepcopy(template)
    for name in ('search_tensor_reference_stage', 'administrative_preflight_statistics'):
        result.pop(name, None)
    result.update(mode='continue', parent_checkpoint=meta['screened_continuation_checkpoint'],
        parent_sha256=meta['screened_continuation_sha256'], parent_gpu_directory=rel(parent),
        parent_implementation_source_sha256=meta['source_sha256'], dt_s=dt,
        reference_sha256=frozen.merge(refs, administrative, session_record), administrative_source_sha256=administrative,
        implementation_transition=None, event_protocol=frozen.PROTOCOL, event_locator_protocol=LOCATOR_PROTOCOL,
        blocked_prefix_protocol=BLOCKED_PREFIX_PROTOCOL, visible_path_protocol=VISIBLE_PATH_PROTOCOL,
        administrative_protocol=PROTOCOL, event_tolerance_scale=.25,
        candidate_forecast_s=envelope['maximum_single_solve_wall_s'], forecast_s=envelope['forecast_s'],
        event_cost_envelope=envelope, sensitivity_contract=rel(previous.CONTRACT), resume_contract=rel(CONTRACT),
        administrative_checks=rel(COMPONENTS), resume_authorized_at=AUTHORIZED_AT.isoformat(),
        resume_absolute_stop=ABSOLUTE_STOP.isoformat(),
        forecast_limit_s=FORECAST_LIMIT_S, internal_limit_s=WORKER_HARD_LIMIT_S,
        external_limit_s=CONTROLLER_LIMIT_S, physical_soft_deadline_s=PHYSICAL_SOFT_LIMIT_S,
        worker_deadline_s=PHYSICAL_SOFT_LIMIT_S)
    result.pop('qualification_wall_reserve_s', None)
    return result


def canonical(document):
    return json.loads(json.dumps(document, allow_nan=False))


def verify_parent(meta, policy):
    require(meta['screened_continuation_sha256'] == PARENT_SHA256
        and meta['screened_end_time_s'] == PARENT_TIME_S and meta['screened_end_step'] == PARENT_STEP,
            'exact registered complete parent checkpoint, time and step required')
    require(policy['parent_sha256'] == PARENT_SHA256 and policy['status'] == 'READY' and policy['pending'] is None,
            'restored adaptive policy must be ready at the exact parent')


def verify_restoration(policy, observation, costs, summary, saved_observation):
    require(policy == summary['adaptive_state'], 'reconstructed policy differs from final registered policy')
    require(canonical(observation) == saved_observation, 'reconstructed actual fine-leaf observation differs')
    require(canonical(measured_envelope(costs)) == summary['next_decision']['event_cost_envelope'],
            'full historical cost evidence or prediction differs')


def verify_old_session(summary, audit):
    require(summary['finished'] is True and summary['exit_code'] == 0 and summary['status'] == 'STOP_COST_REASSESSMENT'
        and summary['completed_stage_count'] == summary['attempt_count'] == 7
        and summary['latest_pass_stage'] == rel(PARENT) and summary['checkpoint_sha256'] == PARENT_SHA256
        and summary['physical_time_s'] == PARENT_TIME_S and summary['screened_steps'] == PARENT_STEP,
            'complete seven-stage v5 cost stop at the exact checkpoint required')
    require(audit['failed'] is False and audit['final_audits'] == 1
        and audit['unique_files'] == audit['expected_files'] == audit['locked_files'] > 0,
            'complete original final history audit required')


def restore_history(pool, args):
    # Reverify all pre-event history with the original verifier and one live lock pool.
    common_parent = previous.COMMON_PARENT
    startup = argparse.Namespace(mode='continue', parent=rel(common_parent), dt=None,
        baseline=rel(CPU_BASE/'stage_013'), cost=rel(BASE/'benchmark_v1'),
        validation=rel(BASE/'replay_v1'), reference_stage=None)
    template = frozen._initial_preflight(startup, pool)
    v1 = dict(template['source_sha256'])
    v2 = previous.secant_admin.secant_sources(pool, v1)
    v3 = previous.box_admin.box_sources(pool, v2)
    v4 = previous.previous_admin.prism_sources(pool, v3)
    template.update(v1_source_sha256=v1, v2_source_sha256=v2, v3_source_sha256=v3,
        v4_source_sha256=v4, source_sha256=previous.delay_sources(pool, v4),
        delay_components_sha256=pool.sha(previous.COMPONENTS), event_locator_protocol=LOCATOR_PROTOCOL,
        blocked_prefix_protocol=BLOCKED_PREFIX_PROTOCOL, visible_path_protocol=VISIBLE_PATH_PROTOCOL)
    old_template = read(PARENT/'protocol.json'); session = read(OLD_CAMPAIGN/'protocol.json')
    for key in ('source_sha256', 'legacy_source_sha256', 'search_source_sha256', 'v1_source_sha256',
                'v2_source_sha256', 'v3_source_sha256', 'v4_source_sha256'):
        require(template[key] == old_template[key], 'frozen source generation map changed: '+key)
    require(len(template['source_sha256']) == 636 and session['source_sha256'] == template['source_sha256']
        and session['schema'] == previous.PROTOCOL, 'exact frozen v5 campaign registration required')
    require(args.maximum_dt == session['maximum_dt_s'] and args.target_time == session['target_time_s']
        and args.shoulder_observe_ms == session['shoulder_observe_ms'], 'resume must preserve old policy and observation bounds')
    summary = read(OLD_CAMPAIGN/'final_summary.json')
    verify_old_session(summary, read(OLD_CAMPAIGN/'final_history_audit.json'))
    evidence = ('protocol.json', 'sensitivity_comparison.json', 'final_summary.json', 'final_history_audit.json',
        'deposition_progress.json', 'one_receipt.json', 'half_receipt.json', 'quarter_receipt.json')
    refs = dict(template['reference_sha256'])
    for name in evidence:
        path = OLD_CAMPAIGN/name; refs = frozen.merge(refs, {rel(path):pool.sha(path)})
    tx = frozen.LockedPreflight(pool); common = tx.gpu_receipt(common_parent, {frozen.CONTINUE_PASS}); tx.finish()
    require(common['source_sha256'] == v1 and common['incremental_replay_particle_arrays_bitwise_equal'] ==
        dict(coarse=True, fine_1=True, fine_2=True) and not any(common['event_records'][r] for r in ROLES),
            'original complete bitwise replay parent required')
    refs = frozen.bind_stage(pool, common_parent, common, refs)
    costs = [cost_evidence(pool, common_parent, common)]; stages = {}
    for scale, label in ((1., 'one'), (.5, 'half'), (.25, 'quarter')):
        folder = OLD_CAMPAIGN/('sensitivity_'+label)
        meta = frozen.validated_stage(pool, folder, template['source_sha256'], common)
        require(meta.get('implementation_transition') == source_transition(template, common, common_parent),
                'original sensitivity bridge changed')
        previous.sibling_identity(meta, common, common_parent, scale, common['dt_s'], template['source_sha256'])
        refs = frozen.bind_stage(pool, folder, meta, refs); costs.append(cost_evidence(pool, folder, meta)); stages[scale] = (folder, meta)
    comparison = previous.compare_group(common, common_parent, stages)
    require(canonical(comparison) == read(OLD_CAMPAIGN/'sensitivity_comparison.json'),
            'recomputed original three-scale comparison differs')
    frames = []; observation, state = frozen.observe(common, frames)
    policy = adaptive.initial_state(common['screened_continuation_sha256'], common['dt_s'], maximum_dt_s=args.maximum_dt)
    parent, meta = stages[.25]; policy, event = previous.selected_pass(policy, meta, state)
    observation, state = frozen.observe(meta, frames)
    for number in range(1, 8):
        child = OLD_CAMPAIGN/f'attempt_{number:04d}'
        next_meta = frozen.validated_stage(pool, child, template['source_sha256'], meta)
        require(source_generation(next_meta, template) == 6 and next_meta['event_tolerance_scale'] == .25
            and next_meta['parent_gpu_directory'] == rel(parent)
            and next_meta.get('implementation_transition') is None, 'ordinary v5 selected chain differs')
        policy = adaptive.begin_attempt(policy, next_meta['dt_s'],
            phase_clipped=next_meta['dt_s'] < policy['nominal_dt_s']*(1.-1e-12))
        policy, event = adaptive.resolve_pass(policy, next_meta['pair_screen'], state['parameters']['dx'],
            state['bore_observation']['particle_ids'], next_meta['screened_continuation_sha256'], original_candidate_gates_passed=True)
        receipt_path = OLD_CAMPAIGN/f'attempt_{number:04d}_receipt.json'; receipt = read(receipt_path)
        require(receipt['verified'] is True and receipt['exit_code'] == 0 and receipt['stage'] == rel(child)
            and receipt['parent'] == rel(parent) and canonical(receipt['adaptive_event']) == canonical(event),
                'original adaptive receipt differs from replayed policy')
        refs = frozen.merge(refs, {rel(receipt_path):pool.sha(receipt_path)})
        refs = frozen.bind_stage(pool, child, next_meta, refs); costs.append(cost_evidence(pool, child, next_meta))
        parent, meta = child, next_meta; observation, state = frozen.observe(meta, frames)
    verify_parent(meta, policy)
    require(parent.resolve() == PARENT.resolve(), 'restored history ended at a different parent directory')
    verify_restoration(policy, observation, costs, summary, read(OLD_CAMPAIGN/'deposition_progress.json'))
    require(state['implementation_source_sha256'] == template['source_sha256'], 'loaded parent physical implementation differs')
    pool.verify(refs)
    return dict(template=template, refs=refs, parent=parent, meta=meta, costs=costs, policy=policy,
        observation=observation, state=state, frames=frames, comparison=comparison)


def run(args):
    output = (ROOT/args.output).resolve()
    require(output.is_relative_to(BASE) and not output.exists(), 'new project-local output required')
    require(0 <= args.max_stages <= 2000, 'bounded continuation stage count required')
    started = datetime.now(timezone.utc); wall = authorized_wall_seconds(args.wall_seconds, started)
    output.mkdir(parents=True); began = time.perf_counter(); deadline = began+wall
    pool = frozen.LockedHistory(ROOT); kernel = mutex = None; code = 0
    progress = dict(status='INITIAL_FULL_HISTORY_VERIFICATION', finished=False, supervisor_pid=os.getpid(),
        child_pid=None, completed_stage_count=0, attempt_count=0, started_at=started.astimezone(admin.LOCAL).isoformat(),
        latest_stop_at=(started+timedelta(seconds=wall)).astimezone(admin.LOCAL).isoformat(),
        authorization_received_at=AUTHORIZED_AT.isoformat(), accepted=False, formal_accepted_time_advanced_s=0.)
    progress.update(admin.snapshot(read(PARENT/'manifest.json'), PARENT), initial_parent_verified_in_this_session=False)
    def update(**fields):
        progress.update(fields, updated_at=admin.wall_now(), wall_seconds=time.perf_counter()-began,
            remaining_wall_seconds=max(0., deadline-time.perf_counter()))
        admin.atomic_json(output/'progress.json', progress)
    update()
    try:
        kernel, mutex = admin.singleton_and_awake()
        administrative = checks(pool); audit = static_audit(); restored = restore_history(pool, args)
        template = restored['template']; parent = restored['parent']; meta = restored['meta']
        phases = admin.phase_table(read(ROOT/'configs/deposition_fill_then_track_20261007.json'))
        session = dict(schema=PROTOCOL, complete_parent=rel(parent), parent_sha256=PARENT_SHA256,
            physical_source_sha256=template['source_sha256'], administrative_source_sha256=administrative,
            physical_source_changed=False, inherited_qualified_campaign=rel(OLD_CAMPAIGN),
            sensitivity_reverified_without_solves=True, initial_preflight_statistics=template['administrative_preflight_statistics'],
            adapter_audit=audit, requested_wall_limit_s=args.wall_seconds, wall_limit_s=wall,
            authorization_received_at=AUTHORIZED_AT.isoformat(), absolute_stop=ABSOLUTE_STOP.isoformat(),
            latest_stop_at=progress['latest_stop_at'], restored_adaptive_state=restored['policy'],
            inherited_cost_evidence=restored['costs'], accepted=False)
        bridge.shared.write_new(output/'protocol.json', session)
        session_record = {rel(output/'protocol.json'):pool.sha(output/'protocol.json')}
        admin.atomic_json(output/'deposition_progress.json', restored['observation'])
        update(status='FULL_HISTORY_VERIFIED_READY_TO_RESUME', adaptive_state=restored['policy'],
            initial_parent_verified_in_this_session=True,
            initial_preflight_statistics=template['administrative_preflight_statistics'], **admin.snapshot(meta, parent))
        code = continue_selected(args, output, pool, deadline, update, template, administrative, session_record,
            restored['refs'], parent, meta, restored['costs'], restored['policy'], restored['observation'],
            restored['state'], restored['frames'], phases)
    except Exception as error:
        code = 2; update(status='STOPPED_SESSION_ERROR', error=repr(error), traceback=traceback.format_exc())
    finally:
        try:
            bridge.shared.write_new(output/'final_history_audit.json', pool.final_audit())
        except Exception as error:
            code = 2; update(status='STOPPED_FINAL_HISTORY_VERIFICATION_FAILURE', final_audit_error=repr(error))
        pool.close()
        if kernel is not None:
            kernel.SetThreadExecutionState(0x80000000); kernel.CloseHandle(mutex)
        update(finished=True, exit_code=code, child_pid=None)
        bridge.shared.write_new(output/'final_summary.json', progress)
    print(json.dumps(progress), flush=True)
    return code


def continue_selected(args, output, pool, deadline, update, template, administrative, session_record,
                      refs, parent, meta, costs, policy, observation, state, frames, phases):
    """Continue the selected physical history with the frozen original retry rules."""
    completed = attempt = 0; shoulder_started = None; shoulder_start_completed = None
    while completed < args.max_stages:
        if (output/'STOP_REQUESTED').exists() or (output/'STOP_AFTER_STAGE').exists():
            update(status='STOPPED_USER_REQUEST'); return 0
        if meta['screened_end_time_s'] >= args.target_time-1e-14:
            update(status='STOPPED_STATIONARY_TARGET_REACHED'); return 0
        first_shoulder = observation['first_observed_supplied_under_shoulder_time_s']
        if first_shoulder is not None and shoulder_started is None:
            shoulder_started = first_shoulder; shoulder_start_completed = completed
        if shoulder_started is not None and meta['screened_end_time_s']-shoulder_started >= args.shoulder_observe_ms/1000 and completed-shoulder_start_completed >= 5:
            update(status='STOPPED_SHOULDER_OBSERVATION_WINDOW_COMPLETE'); return 0
        desired = policy['nominal_dt_s']
        action = cost_decision(meta, desired, max(0., deadline-time.perf_counter()), phases, costs)
        if action['action'] != 'RUN':
            update(status=action['action'], next_decision=action); return 0
        dt = min(action['dt_s'], args.target_time-meta['screened_end_time_s'])
        policy = adaptive.begin_attempt(policy, dt, phase_clipped=dt < desired*(1.-1e-12))
        attempt += 1; child = output/f'attempt_{attempt:04d}'
        protocol = stage_protocol(template, meta, parent, dt, administrative, session_record, refs, .25, costs)
        update(status='RUNNING_GPU_CHILD', attempt_count=attempt, current_attempt=rel(child), current_dt_s=dt,
               current_parent=rel(parent), adaptive_state=policy, event_tolerance_scale=.25)
        exit_code, elapsed = run_stage(pool, child, protocol, deadline, update, output/'STOP_REQUESTED')
        row = dict(attempt=attempt, stage=rel(child), parent=rel(parent), dt_s=dt,
                   exit_code=exit_code, wall_seconds=elapsed)
        if exit_code:
            refs, failure = keep_failure(pool, child, template['source_sha256'], meta, refs)
            screen = read(child/'pair_screen.json') if (child/'pair_screen.json').exists() else None
            gates = False
            if admin.classify_failure(failure, screen) == 'REFINE':
                gates = all(read(child/role/'candidate_qualification.json')['status'] ==
                    'PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN' for role in ROLES)
            recovery = frozen.geometry_budget_recovery(pool, child, failure, meta, policy,
                template['source_sha256'], protocol.get('implementation_transition'))
            if recovery is None:
                policy, event = adaptive.resolve_failure(policy, failure, screen, state['parameters']['dx'],
                    state['bore_observation']['particle_ids'], original_candidate_gates_passed=gates)
            else:
                policy, event = recovery
            row['adaptive_event'] = event
            bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json', row)
            if event['action'] == 'REFINE':
                update(status='RETRY_SMALLER_DT_FROM_SAME_COMPLETE_PARENT', adaptive_state=policy, last_adaptive_event=event)
                continue
            update(status='STOPPED_CHILD_FAILURE', failed_attempt=rel(child), error=failure.get('error', 'worker failed'))
            return 2
        next_meta = frozen.validated_stage(pool, child, template['source_sha256'], meta)
        require(source_generation(next_meta, template) == 6 and next_meta['event_tolerance_scale'] == .25
            and next_meta.get('event_locator_protocol') == LOCATOR_PROTOCOL
            and next_meta.get('implementation_transition') == source_transition(template, meta, parent),
                'selected tolerance and registered secant implementation must persist during continuation')
        refs = frozen.bind_stage(pool, child, next_meta, refs); costs.append(cost_evidence(pool, child, next_meta))
        policy, event = adaptive.resolve_pass(policy, next_meta['pair_screen'], state['parameters']['dx'],
            state['bore_observation']['particle_ids'], next_meta['screened_continuation_sha256'], original_candidate_gates_passed=True)
        parent, meta = child, next_meta; completed += 1
        observation, state = frozen.observe(meta, frames)
        admin.atomic_json(output/'deposition_progress.json', observation)
        row.update(verified=True, adaptive_event=event, locked_history_statistics=pool.stats(), **admin.snapshot(meta, parent))
        bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json', row)
        update(status='VERIFIED_PASS_READY_FOR_NEXT_STAGE', completed_stage_count=completed, child_pid=None,
               adaptive_state=policy, last_adaptive_event=event, **admin.snapshot(meta, parent))
    update(status='STOPPED_STAGE_LIMIT')
    return 0

def run_stage(pool, output, protocol, controller_deadline, update, stop_file):
    from afsd_event_resume_90m_worker_20261010 import stage_worker
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
            if time.perf_counter() >= min(started+1380., controller_deadline) or stop_file.exists():
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

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='outputs/gpu_acceleration_20261007/event_resume_90m_20261010_v1')
    parser.add_argument('--wall-seconds', type=float, default=5400.)
    parser.add_argument('--max-stages', type=int, default=1000, help='continuation stages after full existing-history verification')
    parser.add_argument('--maximum-dt', type=float, default=.000474609375)
    parser.add_argument('--target-time', type=float, default=.22)
    parser.add_argument('--shoulder-observe-ms', type=float, default=5.)
    parser.add_argument('--detach', action='store_true'); parser.add_argument('--detached-child', action='store_true')
    args = parser.parse_args()
    if args.detach or args.detached_child:
        import subprocess
        from start_gpu_resume_detached_20261008 import job_policy, write_new
        output = (ROOT/args.output).resolve()
        require(output.is_relative_to(BASE) and not output.exists(), 'new detached output required')
        if args.detached_child:
            policy = job_policy()
            write_new(output.with_name(output.name+'.bootstrap.json'), dict(**policy, supervisor=rel(Path(__file__)),
                status='PASS_DETACHED_JOB_POLICY' if not policy['kill_on_job_close'] else 'FAIL_JOB_POLICY'))
            require(not policy['kill_on_job_close'], 'detached session retains kill-on-close')
        else:
            require(os.name == 'nt', 'Windows detached launch required')
            forwarded = [v for v in sys.argv[1:] if v != '--detach']
            command = [sys.executable, '-B', str(Path(__file__).resolve()), *forwarded, '--detached-child']
            startup = subprocess.STARTUPINFO(); startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW; startup.wShowWindow = 0
            flags = subprocess.CREATE_BREAKAWAY_FROM_JOB|subprocess.DETACHED_PROCESS|subprocess.CREATE_NEW_PROCESS_GROUP
            env = dict(os.environ, **{k:'1' for k in bridge.shared.THREADS}); env['PYTHONPATH'] = ''
            with output.with_name(output.name+'.stdout.log').open('xb', buffering=0) as stream:
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL, stdout=stream,
                    stderr=subprocess.STDOUT, startupinfo=startup, creationflags=flags, close_fds=True)
            receipt = dict(status='STARTED_90M_RESUME_PENDING_BOOTSTRAP', pid=process.pid, command=command,
                creationflags=flags, source_sha256={rel(Path(__file__)): bridge.core.sha(Path(__file__))})
            write_new(output.with_name(output.name+'.launch.json'), receipt); print(json.dumps(receipt)); return 0
    return run(args)

if __name__ == '__main__':
    raise SystemExit(main())

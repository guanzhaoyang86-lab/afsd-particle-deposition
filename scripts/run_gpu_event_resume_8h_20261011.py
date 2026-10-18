"""Eight-hour administrative continuation of unchanged qualified v5 physics."""
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
import run_gpu_event_resume_90m_20261010 as previous

ROOT, BASE, CPU_BASE = previous.ROOT, previous.BASE, previous.CPU_BASE
frozen, bridge, admin, adaptive = previous.frozen, previous.bridge, previous.admin, previous.adaptive
v5 = previous.previous
require, read, load, rel = previous.require, previous.read, previous.load, previous.rel
source_generation, source_transition = previous.source_generation, previous.source_transition
cost_evidence, measured_envelope = previous.cost_evidence, previous.measured_envelope
authorize_grant = previous.authorize_grant
LOCATOR_PROTOCOL = previous.LOCATOR_PROTOCOL
BLOCKED_PREFIX_PROTOCOL = previous.BLOCKED_PREFIX_PROTOCOL
VISIBLE_PATH_PROTOCOL = previous.VISIBLE_PATH_PROTOCOL
ROLES = previous.ROLES
canonical = previous.canonical
PROTOCOL = 'FROZEN_V5_EIGHT_HOUR_DYNAMIC_ADMINISTRATIVE_RESUME_20261011_V1'
SESSION_PROTOCOL = PROTOCOL
AUTHORIZED_AT = datetime.fromisoformat('2026-10-11T15:41:47+08:00')
ABSOLUTE_STOP = datetime.fromisoformat('2026-10-11T23:41:47+08:00')
CONTRACT = ROOT/'docs/GPU_EVENT_RESUME_8H_20261011.md'
COMPONENTS = BASE/'event_resume_8h_components_20261011_v1.json'
OLD_CAMPAIGN = BASE/'event_resume_90m_20261010_v1'
PARENT = OLD_CAMPAIGN/'attempt_0013'
FAILED = OLD_CAMPAIGN/'attempt_0014'
PARENT_SHA256 = '0079328d9ea46263d62cf975299497ea6774136d22e5bdf0ffba5d6cdcb91c97'
PARENT_TIME_S = .08134326171874973
PARENT_STEP = 527
ADMINISTRATIVE_FILES = ('scripts/run_gpu_event_resume_8h_20261011.py',
    'scripts/afsd_event_resume_8h_worker_20261011.py', 'scripts/check_gpu_event_resume_8h_20261011.py',
    'tests/review/test_gpu_event_resume_8h_20261011.py', 'docs/GPU_EVENT_RESUME_8H_20261011.md',
    'scripts/afsd_multi_event_refinement_20261011.py', 'tests/review/test_afsd_multi_event_refinement_20261011.py')


def authorized_wall_seconds(requested, now):
    require(type(requested) in (int, float) and math.isfinite(requested) and 0 < requested <= 28800.,
            'finite requested authorization at most 28800 seconds required')
    require(now.tzinfo is not None and now >= AUTHORIZED_AT, 'timezone-aware request time required')
    remaining = (ABSOLUTE_STOP-now).total_seconds()
    require(remaining > 0, 'eight-hour absolute authorization expired')
    return min(float(requested), remaining)


def checks(pool):
    inherited = previous.checks(pool); receipt = read(COMPONENTS)
    require(receipt['status'] == 'PASS_EVENT_RESUME_8H_COMPONENTS' and receipt['exit_code'] == 0
        and receipt['physical_source_changed'] is False and receipt['gpu_campaign_runs'] == 0,
            'executed administrative-only eight-hour checks required')
    require(receipt['inherited_receipt'] == rel(previous.COMPONENTS)
        and receipt['inherited_receipt_sha256'] == pool.sha(previous.COMPONENTS), 'frozen 90-minute receipt required')
    require(set(ADMINISTRATIVE_FILES).issubset(receipt['source_sha256'])
        and all(receipt['source_sha256'].get(k) == v for k, v in inherited.items()), 'old closure must remain exact')
    require(receipt['physical_source_sha256'] == read(PARENT/'manifest.json')['source_sha256']
        and len(receipt['physical_source_sha256']) == 636, 'unchanged 636-entry physical source required')
    require(receipt['parent_directory'] == rel(PARENT) and receipt['parent_checkpoint_sha256'] == PARENT_SHA256
        and receipt['failed_stage'] == rel(FAILED) and receipt['absolute_stop'] == ABSOLUTE_STOP.isoformat(),
            'registered full parent, failed evidence and absolute deadline required')
    pool.seal(receipt['source_sha256']); pool.seal(receipt['artifacts_sha256'])
    return frozen.merge(inherited, receipt['source_sha256'], receipt['artifacts_sha256'], {rel(COMPONENTS):pool.sha(COMPONENTS)})



def static_audit():
    import afsd_event_resume_8h_worker_20261011 as worker
    old=inspect.getsource(__import__('afsd_event_resume_90m_worker_20261010').stage_worker)
    new=inspect.getsource(worker.stage_worker)
    new=new.replace("    soft_deadline,hard_seconds=worker_limits(protocol,began,time.time(),time.perf_counter())\n",'')
    new=new.replace('threading.Timer(hard_seconds,timeout)','threading.Timer(1290.,timeout)').replace('soft_deadline','began+1260.')
    new=new.replace("            physical.core.write(output/'WORKER_PRE_END_FAILURE.json',dict(status='RECORDED_WORKER_FAILURE_BEFORE_END_GRANT',error=repr(error)))\n",'')
    require(ast.dump(ast.parse(new))==ast.dump(ast.parse(old)),'worker physical logic differs beyond dynamic clocks and failure journal')
    require(inspect.getsource(worker.request_grant)==inspect.getsource(__import__('afsd_event_resume_90m_worker_20261010').request_grant),
        'original live-history grants changed')
    require((ABSOLUTE_STOP-AUTHORIZED_AT).total_seconds()==28800.,'eight-hour fixed authorization required')
    return dict(status='PASS_DYNAMIC_ADMINISTRATIVE_8H_RESUME',physical_source_changed=False,
        frozen_audit=previous.static_audit(),dynamic_soft_formula='min(3600,max(1260,forecast+60))',hard_extra_s=30.,
        controller_extra_s=120.,all_limits_capped_by_authorization=True,forecast_is_admission_limit=False,
        shoulder_observation_is_stop_condition=False,absolute_stop=ABSOLUTE_STOP.isoformat())


def dynamic_budget(costs, remaining):
    require(type(remaining) in (int,float) and math.isfinite(remaining) and remaining >= 0., 'finite remaining wall required')
    envelope = measured_envelope(costs)
    forecast = envelope['forecast_s']; single = envelope['maximum_single_solve_wall_s']
    reserve = 1.5*single+20.
    nominal_soft = min(3600., max(1260., forecast+60.))
    return dict(event_cost_envelope=envelope, forecast_s=forecast,
        physical_soft_deadline_s=min(nominal_soft,remaining), internal_limit_s=min(nominal_soft+30.,remaining),
        external_limit_s=min(nominal_soft+120.,remaining), remaining_wall_seconds=remaining,
        minimum_start_reserve_s=reserve, best_effort=forecast+20. > remaining,
        forecast_is_admission_limit=False, action='RUN' if remaining >= reserve else 'STOP_WALL_RESERVE')


def cost_decision(meta, dt, remaining, phases, costs):
    budget = dynamic_budget(costs, remaining)
    if budget['action'] != 'RUN': return budget
    t = float(meta['screened_end_time_s'])
    require(math.isfinite(t) and t >= 0 and type(dt) in (int,float) and math.isfinite(dt) and 1e-7 <= dt <= .005,
            'finite bounded physical interval required')
    phase = next((p for p in phases if p['start'] <= t < p['end']),None)
    if phase is None or not phase['stationary']: return dict(budget,action='STOP_REGISTERED_STATIONARY_FILL_END')
    try: bounded = admin.bounded_dt(t,dt,phase['end'])
    except ValueError as error: return dict(budget,action='STOP_UNREPRESENTABLE_PHASE_BOUNDARY',reason=str(error))
    if bounded < 1e-7: return dict(budget,action='STOP_MINIMUM_DT_OR_BOUNDARY')
    return dict(budget,dt_s=bounded,phase=phase['name'],phase_end_s=phase['end'])


def stage_protocol(template, meta, parent, dt, administrative, session_record, refs, costs, remaining):
    require(meta['source_sha256'] == template['source_sha256'] and meta.get('event_tolerance_scale')==.25 and source_generation(meta,template) == 6
        and source_transition(template,meta,parent) is None, 'same qualified v5 implementation required')
    budget = dynamic_budget(costs,remaining); require(budget['action']=='RUN','insufficient measured single-solve wall reserve')
    frozen.merge(template['source_sha256'],refs,administrative,session_record)
    result=copy.deepcopy(template)
    for key in ('search_tensor_reference_stage','administrative_preflight_statistics','qualification_wall_reserve_s'):
        result.pop(key,None)
    result.update(mode='continue',parent_checkpoint=meta['screened_continuation_checkpoint'],
        parent_sha256=meta['screened_continuation_sha256'],parent_gpu_directory=rel(parent),
        parent_implementation_source_sha256=meta['source_sha256'],dt_s=dt,
        reference_sha256=frozen.merge(refs,administrative,session_record),administrative_source_sha256=administrative,
        implementation_transition=None,event_protocol=frozen.PROTOCOL,event_locator_protocol=LOCATOR_PROTOCOL,
        blocked_prefix_protocol=BLOCKED_PREFIX_PROTOCOL,visible_path_protocol=VISIBLE_PATH_PROTOCOL,
        event_tolerance_scale=.25,administrative_protocol=PROTOCOL,administrative_checks=rel(COMPONENTS),
        resume_contract=rel(CONTRACT),sensitivity_contract=rel(v5.CONTRACT),resume_authorized_at=AUTHORIZED_AT.isoformat(),
        resume_absolute_stop=ABSOLUTE_STOP.isoformat(),absolute_deadline_epoch_s=min(ABSOLUTE_STOP.timestamp(),template['session_deadline_epoch_s']),
        candidate_forecast_s=budget['event_cost_envelope']['maximum_single_solve_wall_s'],forecast_s=budget['forecast_s'],
        event_cost_envelope=budget['event_cost_envelope'],dynamic_wall_budget=budget,
        forecast_limit_s=None,forecast_is_admission_limit=False,internal_limit_s=budget['internal_limit_s'],
        external_limit_s=budget['external_limit_s'],physical_soft_deadline_s=budget['physical_soft_deadline_s'],
        worker_deadline_s=budget['physical_soft_deadline_s'])
    return result


def verify_parent(meta):
    require(meta['screened_continuation_sha256']==PARENT_SHA256 and meta['screened_end_time_s']==PARENT_TIME_S
        and meta['screened_end_step']==PARENT_STEP, 'exact complete attempt13 state required')


def retain_failure(pool, folder, protocol, parent_meta, refs):
    failure={}; supervision=read(folder/'supervision.json')
    require(supervision['exit_code'] != 0 and supervision.get('source_error') is None, 'failed stage source verification required')
    if (folder/'manifest.json').exists():
        refs,failure=v5.keep_failure(pool,folder,protocol['source_sha256'],parent_meta,refs)
    evidence={rel(path):pool.sha(path) for path in folder.rglob('*') if path.is_file()}
    pool.seal(evidence); refs=frozen.merge(refs,evidence)
    cost=dict(stage=rel(folder),supervision_sha256=evidence[rel(folder/'supervision.json')],
        supervised_wall_s=supervision['wall_seconds'],worker_wall_s=failure.get('wall_seconds'),
        complete_pair=False,usable_parent=False,retained_artifacts_sha256=evidence)
    return refs,failure,cost


def multi_event_recovery(pool, folder, failure, parent_meta, policy, source):
    import afsd_multi_event_refinement_20261011 as multi
    if failure.get('error') != "GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')": return None
    supervision=read(folder/'supervision.json')
    require(supervision['exit_code']!=0 and pool.sha(folder/'manifest.json')==supervision['manifest_sha256']
        and pool.sha(folder/'run.log')==supervision['log_sha256'] and read(folder/'manifest.json')==failure,
        'multi-event failed receipt identity mismatch')
    require(failure['source_sha256']==source==parent_meta['source_sha256']
        and failure['parent_implementation_source_sha256']==source and failure.get('implementation_transition') is None
        and failure['parent_sha256']==parent_meta['screened_continuation_sha256']
        and failure['parent_checkpoint']==parent_meta['screened_continuation_checkpoint']
        and failure['dt_s']==policy['pending']['dt_s'], 'multi-event source, parent or interval mismatch')
    pool.seal(failure['source_sha256']);pool.seal(failure['reference_sha256']);pool.seal(failure['artifacts_sha256'])
    if (folder/'pair_screen.json').exists(): return None
    def evidence(path):
        require(rel(path) in failure['artifacts_sha256'] and pool.sha(path)==failure['artifacts_sha256'][rel(path)],
                'multi-event evidence missing or changed')
        return read(path)
    candidate=None
    for role in ROLES:
        directory=folder/role
        if candidate is not None:
            require(not directory.exists(),'later candidate after failed multi-event role');continue
        require(directory.is_dir(),'missing preceding candidate role')
        qualification=directory/'candidate_qualification.json'
        if qualification.exists():
            require(evidence(qualification)['status']=='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN','preceding role not qualified')
        else: candidate=directory
    require(candidate is not None,'unqualified multi-event candidate missing')
    report_names=(
        ('visible_path','visible_path_audit.json'),('visible_shell','continuous_shell_audit.json'),
        ('geometry','geometry_audit.json'),('event_router','visibility_event_routing.json'),
        ('union_event_router','proposal/union_visibility_events.json'))
    reports={};unsupported=False
    for key,name in report_names:
        path=candidate/name
        if rel(path) not in failure['artifacts_sha256']:
            unsupported=True
        else:
            reports[key]=evidence(path)
    if unsupported: return None
    diagnosis=multi.recognize(failure,supervision,reports)
    if not diagnosis['eligible']: return None
    diagnosis['candidate_role']=candidate.name
    policy,event=multi.refine(policy,diagnosis);event['candidate_role']=candidate.name
    return policy,event


def timeout_recovery(pool, folder, failure, parent_meta, policy, protocol):
    supervision=read(folder/'supervision.json')
    if supervision.get('termination_reason') not in ('STAGE_WALL_LIMIT','WORKER_HARD_TIMEOUT'): return None
    # A captured physical/gate error always outranks an eventual wall timeout.
    if failure or (folder/'WORKER_PRE_END_FAILURE.json').exists(): return None
    require(supervision.get('source_error') is None and supervision['exit_code']!=0,
        'timeout source verification failed')
    if supervision['termination_reason']=='WORKER_HARD_TIMEOUT':
        require(supervision['exit_code']==124 and (folder/'HARD_TIMEOUT.json').exists()
            and read(folder/'HARD_TIMEOUT.json').get('status')=='TIMEOUT_NOT_ACCEPTED', 'explicit worker hard timeout marker required')
    else:
        require(supervision['interrupted'] is True and supervision['wall_seconds']>=protocol['external_limit_s'],
            'actual external stage wall limit evidence required')
    require(protocol['source_sha256']==parent_meta['source_sha256']
        and protocol['parent_sha256']==parent_meta['screened_continuation_sha256']
        and protocol['parent_checkpoint']==parent_meta['screened_continuation_checkpoint']
        and protocol['dt_s']==policy['pending']['dt_s'], 'timeout complete-parent identity mismatch')
    pool.verify(protocol['source_sha256']);pool.verify(protocol['reference_sha256'])
    adaptive.validate_state(policy); pending=policy['pending']; config=policy['config']; smaller=pending['dt_s']/2.
    result=copy.deepcopy(policy); action='STOP';reason='TIMEOUT_REFINEMENT_LIMIT'
    result.update(status='STOPPED',pending=None,small_error_streak=0,growth_pending=False)
    if policy['consecutive_halvings'] < config['max_halvings'] and smaller >= max(1e-7,config['minimum_dt_s']):
        result.update(status='READY',nominal_dt_s=smaller,consecutive_halvings=policy['consecutive_halvings']+1,
            cooldown_remaining=max(policy['cooldown_remaining'],config['cooldown_passes']),
            growth_failures=policy['growth_failures']+int(pending['growth_trial']))
        action='REFINE';reason='ADMINISTRATIVE_TIMEOUT_RETRY_UNCHANGED_PARENT'
    adaptive.validate_state(result)
    return result,dict(action=action,reason=reason,attempt=policy['attempts'],parent_before_sha256=policy['parent_sha256'],
        parent_after_sha256=result['parent_sha256'],attempted_dt_s=pending['dt_s'],next_nominal_dt_s=result['nominal_dt_s'],
        termination_reason=supervision['termination_reason'],whole_pair_qualified=False,accepted=False,restart_allowed=False,
        accepted_time_advanced_s=0.)


def restore_history(pool,args):
    restored=previous.restore_history(pool,args)
    template=restored['template'];refs=restored['refs'];parent=restored['parent'];meta=restored['meta']
    costs=restored['costs'];policy=restored['policy'];state=restored['state'];frames=restored['frames']
    session=read(OLD_CAMPAIGN/'protocol.json');summary=read(OLD_CAMPAIGN/'final_summary.json');audit=read(OLD_CAMPAIGN/'final_history_audit.json')
    require(session['schema']==previous.PROTOCOL and session['physical_source_sha256']==template['source_sha256']
        and session['restored_adaptive_state']==policy, 'original 90-minute source and initial policy differ')
    require(summary['finished'] is True and summary['exit_code']==2 and summary['status']=='STOPPED_CHILD_FAILURE'
        and summary['completed_stage_count']==13 and summary['attempt_count']==14
        and summary['latest_pass_stage']==rel(PARENT) and summary['failed_attempt']==rel(FAILED)
        and summary['checkpoint_sha256']==PARENT_SHA256 and summary['physical_time_s']==PARENT_TIME_S
        and summary['screened_steps']==PARENT_STEP,'exact failed attempt14 session termination required')
    require(audit['failed'] is False and audit['final_audits']==1
        and audit['unique_files']==audit['expected_files']==audit['locked_files']>0,'complete successful old final hash audit required')
    for name in ('protocol.json','final_summary.json','final_history_audit.json','deposition_progress.json'):
        path=OLD_CAMPAIGN/name;refs=frozen.merge(refs,{rel(path):pool.sha(path)})
    for number in range(1,14):
        child=OLD_CAMPAIGN/f'attempt_{number:04d}';new=frozen.validated_stage(pool,child,template['source_sha256'],meta)
        require(source_generation(new,template)==6 and new['event_tolerance_scale']==.25
            and new['parent_gpu_directory']==rel(parent) and new.get('implementation_transition') is None,
            'unchanged quarter-scale 90-minute chain required')
        policy=adaptive.begin_attempt(policy,new['dt_s'],phase_clipped=new['dt_s']<policy['nominal_dt_s']*(1.-1e-12))
        policy,event=adaptive.resolve_pass(policy,new['pair_screen'],state['parameters']['dx'],state['bore_observation']['particle_ids'],
            new['screened_continuation_sha256'],original_candidate_gates_passed=True)
        path=OLD_CAMPAIGN/f'attempt_{number:04d}_receipt.json';receipt=read(path)
        require(receipt['verified'] is True and receipt['exit_code']==0 and receipt['stage']==rel(child)
            and receipt['parent']==rel(parent) and canonical(receipt['adaptive_event'])==canonical(event),
            '90-minute adaptive replay receipt mismatch')
        refs=frozen.merge(refs,{rel(path):pool.sha(path)});refs=frozen.bind_stage(pool,child,new,refs)
        costs.append(cost_evidence(pool,child,new));parent,meta=child,new
        observation,state=frozen.observe(meta,frames)
    verify_parent(meta);require(parent.resolve()==PARENT.resolve(),'restored wrong complete parent directory')
    require(canonical(observation)==read(OLD_CAMPAIGN/'deposition_progress.json'),'selected actual fine-leaf history differs')
    require(canonical(event)==summary['last_adaptive_event'],'last PASS adaptive event differs')
    policy_before_failure=copy.deepcopy(policy);failure=read(FAILED/'manifest.json');old_protocol=read(FAILED/'protocol.json')
    require(old_protocol['event_cost_envelope']==canonical(measured_envelope(costs)),'all 24 completed pair costs must remain exact')
    require(old_protocol['source_sha256']==template['source_sha256']==state['implementation_source_sha256'],'loaded complete parent source differs')
    policy=adaptive.begin_attempt(policy,failure['dt_s'],phase_clipped=failure['dt_s']<policy['nominal_dt_s']*(1.-1e-12))
    require(policy==summary['adaptive_state'],'replayed failed pending attempt differs from old summary')
    refs,failure,failed_cost=retain_failure(pool,FAILED,old_protocol,meta,refs)
    old_receipt_path=OLD_CAMPAIGN/'attempt_0014_receipt.json';old_receipt=read(old_receipt_path)
    old_stopped,old_event=adaptive.resolve_failure(policy,failure,None,state['parameters']['dx'],state['bore_observation']['particle_ids'],
        original_candidate_gates_passed=False)
    require(canonical(old_event)==old_receipt['adaptive_event'] and old_receipt['exit_code']!=0,'original unsupported failure receipt differs')
    refs=frozen.merge(refs,{rel(old_receipt_path):pool.sha(old_receipt_path)})
    recovery=multi_event_recovery(pool,FAILED,failure,meta,policy,template['source_sha256'])
    require(recovery is not None and recovery[1]['action']=='REFINE','old multi-event failure must qualify exact bounded same-parent recovery')
    policy,event=recovery
    require(policy['parent_sha256']==PARENT_SHA256 and policy['status']=='READY' and policy['pending'] is None,'recovery must retain exact complete parent')
    pool.verify(refs)
    restored.update(refs=refs,parent=parent,meta=meta,costs=costs,policy=policy,observation=observation,state=state,frames=frames,
        previous_failed_costs=[failed_cost],initial_recovery=event,policy_before_failed_attempt=policy_before_failure)
    return restored


def continue_selected(args,output,pool,deadline,update,template,administrative,session_record,
                      refs,parent,meta,costs,policy,observation,state,frames,phases,failed_costs):
    completed=attempt=0;shoulder_milestone=False
    while completed < args.max_stages:
        if (output/'STOP_REQUESTED').exists() or (output/'STOP_AFTER_STAGE').exists():
            update(status='STOPPED_USER_REQUEST');return 0
        if meta['screened_end_time_s'] >= args.target_time-1e-14:
            update(status='STOPPED_STATIONARY_TARGET_REACHED');return 0
        if observation['first_observed_supplied_under_shoulder_time_s'] is not None and not shoulder_milestone:
            shoulder_milestone=True
            update(shoulder_observed=True,shoulder_first_time_s=observation['first_observed_supplied_under_shoulder_time_s'])
        desired=policy['nominal_dt_s'];remaining=max(0.,deadline-time.perf_counter())
        action=cost_decision(meta,desired,remaining,phases,costs)
        if action['action']!='RUN':
            update(status=action['action'],next_decision=action);return 0
        dt=min(action['dt_s'],args.target_time-meta['screened_end_time_s'])
        policy=adaptive.begin_attempt(policy,dt,phase_clipped=dt<desired*(1.-1e-12))
        attempt+=1;child=output/f'attempt_{attempt:04d}'
        protocol=stage_protocol(template,meta,parent,dt,administrative,session_record,refs,costs,max(0.,deadline-time.perf_counter()))
        update(status='RUNNING_GPU_CHILD',attempt_count=attempt,current_attempt=rel(child),current_dt_s=dt,
            current_parent=rel(parent),adaptive_state=policy,event_tolerance_scale=.25,dynamic_wall_budget=protocol['dynamic_wall_budget'])
        exit_code,elapsed=run_stage(pool,child,protocol,deadline,update,output/'STOP_REQUESTED')
        row=dict(attempt=attempt,stage=rel(child),parent=rel(parent),dt_s=dt,exit_code=exit_code,wall_seconds=elapsed)
        if exit_code:
            refs,failure,resource_cost=retain_failure(pool,child,protocol,meta,refs);failed_costs.append(resource_cost)
            bridge.shared.write_new(output/f'attempt_{attempt:04d}_failed_cost.json',resource_cost)
            failed_cost_path=output/f'attempt_{attempt:04d}_failed_cost.json'
            refs=frozen.merge(refs,{rel(failed_cost_path):pool.sha(failed_cost_path)})
            supervision=read(child/'supervision.json')
            if supervision.get('termination_reason') in ('USER_STOP','SESSION_DEADLINE'):
                row.update(accepted=False,retained_complete_parent_sha256=meta['screened_continuation_sha256'],termination_reason=supervision['termination_reason'])
                bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json',row)
                update(status='STOPPED_USER_REQUEST' if supervision['termination_reason']=='USER_STOP' else 'STOPPED_AUTHORIZATION_DEADLINE',
                    child_pid=None,retained_complete_parent_sha256=meta['screened_continuation_sha256']);return 0
            recovery=multi_event_recovery(pool,child,failure,meta,policy,template['source_sha256'])
            if recovery is None:
                recovery=frozen.geometry_budget_recovery(pool,child,failure,meta,policy,template['source_sha256'],None)
            if recovery is None:
                recovery=timeout_recovery(pool,child,failure,meta,policy,protocol)
            if recovery is None:
                screen=read(child/'pair_screen.json') if (child/'pair_screen.json').exists() else None
                gates=False
                if admin.classify_failure(failure,screen)=='REFINE':
                    gates=all(read(child/role/'candidate_qualification.json')['status']=='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN' for role in ROLES)
                policy,event=adaptive.resolve_failure(policy,failure,screen,state['parameters']['dx'],
                    state['bore_observation']['particle_ids'],original_candidate_gates_passed=gates)
            else: policy,event=recovery
            row['adaptive_event']=event;bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json',row)
            if event['action']=='REFINE':
                update(status='RETRY_SMALLER_DT_FROM_SAME_COMPLETE_PARENT',adaptive_state=policy,last_adaptive_event=event,child_pid=None);continue
            update(status='STOPPED_CHILD_FAILURE',failed_attempt=rel(child),adaptive_state=policy,last_adaptive_event=event,
                error=failure.get('error','worker failed without a recoverable administrative timeout'),child_pid=None);return 2
        new=frozen.validated_stage(pool,child,template['source_sha256'],meta)
        require(source_generation(new,template)==6 and new['event_tolerance_scale']==.25
            and new['parent_gpu_directory']==rel(parent) and new.get('implementation_transition') is None,
            'unchanged qualified v5 generation and quarter tolerance must persist')
        refs=frozen.bind_stage(pool,child,new,refs);costs.append(cost_evidence(pool,child,new))
        policy,event=adaptive.resolve_pass(policy,new['pair_screen'],state['parameters']['dx'],state['bore_observation']['particle_ids'],
            new['screened_continuation_sha256'],original_candidate_gates_passed=True)
        parent,meta=child,new;completed+=1;observation,state=frozen.observe(meta,frames)
        admin.atomic_json(output/'deposition_progress.json',observation)
        row.update(verified=True,adaptive_event=event,locked_history_statistics=pool.stats(),**admin.snapshot(meta,parent))
        bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json',row)
        update(status='VERIFIED_PASS_READY_FOR_NEXT_STAGE',completed_stage_count=completed,child_pid=None,
            adaptive_state=policy,last_adaptive_event=event,**admin.snapshot(meta,parent))
    update(status='STOPPED_STAGE_LIMIT');return 0


def run(args):
    output=(ROOT/args.output).resolve();require(output.is_relative_to(BASE) and not output.exists(),'new project-local output required')
    require(type(args.max_stages) is int and 0 <= args.max_stages <= 100000,'bounded continuation count required')
    started=datetime.now(timezone.utc);wall=authorized_wall_seconds(args.wall_seconds,started)
    output.mkdir(parents=True);began=time.perf_counter();deadline=began+wall
    pool=frozen.LockedHistory(ROOT);kernel=mutex=None;code=0
    progress=dict(status='INITIAL_FULL_HISTORY_VERIFICATION',finished=False,supervisor_pid=os.getpid(),child_pid=None,
        completed_stage_count=0,attempt_count=0,started_at=started.astimezone(admin.LOCAL).isoformat(),
        latest_stop_at=(started+timedelta(seconds=wall)).astimezone(admin.LOCAL).isoformat(),authorization_received_at=AUTHORIZED_AT.isoformat(),
        initial_parent_verified_in_this_session=False,accepted=False,formal_accepted_time_advanced_s=0.,**admin.snapshot(read(PARENT/'manifest.json'),PARENT))
    def update(**fields):
        progress.update(fields,updated_at=admin.wall_now(),wall_seconds=time.perf_counter()-began,remaining_wall_seconds=max(0.,deadline-time.perf_counter()))
        admin.atomic_json(output/'progress.json',progress)
    update()
    try:
        kernel,mutex=admin.singleton_and_awake();administrative=checks(pool);audit=static_audit();restored=restore_history(pool,args)
        template=restored['template'];template['session_deadline_epoch_s']=(started+timedelta(seconds=wall)).timestamp()
        parent,meta=restored['parent'],restored['meta'];phases=admin.phase_table(read(ROOT/'configs/deposition_fill_then_track_20261007.json'))
        session=dict(schema=PROTOCOL,complete_parent=rel(parent),parent_sha256=PARENT_SHA256,
            physical_source_sha256=template['source_sha256'],administrative_source_sha256=administrative,physical_source_changed=False,
            inherited_campaign=rel(OLD_CAMPAIGN),sensitivity_reverified_without_solves=True,
            initial_preflight_statistics=template['administrative_preflight_statistics'],adapter_audit=audit,
            requested_wall_limit_s=args.wall_seconds,wall_limit_s=wall,authorization_received_at=AUTHORIZED_AT.isoformat(),
            absolute_stop=ABSOLUTE_STOP.isoformat(),latest_stop_at=progress['latest_stop_at'],restored_adaptive_state=restored['policy'],
            inherited_cost_evidence=restored['costs'],inherited_failed_costs=restored['previous_failed_costs'],
            old_failure=rel(FAILED),initial_multi_event_recovery=restored['initial_recovery'],
            forecast_is_admission_limit=False,shoulder_observation_is_stop_condition=False,accepted=False)
        bridge.shared.write_new(output/'protocol.json',session);session_record={rel(output/'protocol.json'):pool.sha(output/'protocol.json')}
        admin.atomic_json(output/'deposition_progress.json',restored['observation'])
        update(status='FULL_HISTORY_VERIFIED_READY_TO_RESUME',initial_parent_verified_in_this_session=True,
            adaptive_state=restored['policy'],last_adaptive_event=restored['initial_recovery'],
            initial_preflight_statistics=template['administrative_preflight_statistics'],**admin.snapshot(meta,parent))
        code=continue_selected(args,output,pool,deadline,update,template,administrative,session_record,restored['refs'],parent,meta,
            restored['costs'],restored['policy'],restored['observation'],restored['state'],restored['frames'],phases,restored['previous_failed_costs'])
    except Exception as error:
        code=2;update(status='STOPPED_SESSION_ERROR',error=repr(error),traceback=traceback.format_exc())
    finally:
        try: bridge.shared.write_new(output/'final_history_audit.json',pool.final_audit())
        except Exception as error:
            code=2;update(status='STOPPED_FINAL_HISTORY_VERIFICATION_FAILURE',final_audit_error=repr(error))
        pool.close()
        if kernel is not None: kernel.SetThreadExecutionState(0x80000000);kernel.CloseHandle(mutex)
        update(finished=True,exit_code=code,child_pid=None);bridge.shared.write_new(output/'final_summary.json',progress)
    print(json.dumps(progress),flush=True);return code


def run_stage(pool, output, protocol, controller_deadline, update, stop_file):
    from afsd_event_resume_8h_worker_20261011 import stage_worker
    output.mkdir();bridge.shared.write_new(output/'protocol.json',protocol);protocol_sha=pool.sha(output/'protocol.json')
    pool.seal(protocol['source_sha256']);pool.seal(protocol['reference_sha256'])
    context=multiprocessing.get_context('spawn');server,child_channel=context.Pipe(duplex=True)
    child=context.Process(target=stage_worker,args=(str(output),child_channel,protocol_sha),daemon=False)
    started=time.perf_counter();history_before=pool.stats();child.start();child_channel.close();update(child_pid=child.pid)
    grants=[];interrupted=False;error=None;last_update=started;termination_reason=None
    try:
        while child.is_alive():
            now=time.perf_counter()
            if now >= min(started+protocol['external_limit_s'],controller_deadline) or stop_file.exists():
                interrupted=True
                termination_reason='USER_STOP' if stop_file.exists() else ('SESSION_DEADLINE' if now>=controller_deadline else 'STAGE_WALL_LIMIT')
                with (output/'termination.log').open('a',encoding='utf-8') as stream:
                    import subprocess
                    subprocess.run(['taskkill','/PID',str(child.pid),'/T','/F'],stdout=stream,stderr=subprocess.STDOUT,timeout=20.)
                child.join(timeout=20.);break
            if len(grants)<2 and server.poll(.25):
                try: request=server.recv()
                except EOFError:
                    child.join(timeout=.25);continue
                phase=('start','end')[len(grants)] if len(grants)<2 else None
                response=authorize_grant(request,protocol_sha,phase,pool,protocol);grants.append(response);server.send(response)
            elif len(grants)==2: child.join(timeout=.25)
            if time.perf_counter()-last_update>=5.: update();last_update=time.perf_counter()
        child.join(timeout=5.);code=child.exitcode if child.exitcode is not None else 124
        require(not child.is_alive(),'worker could not be terminated')
        pool.verify(protocol['source_sha256']);pool.verify(protocol['reference_sha256'])
        require(code!=0 or len(grants)==2,'successful worker lacks both live grants')
        if termination_reason is None and code==124 and (output/'HARD_TIMEOUT.json').exists():
            marker=read(output/'HARD_TIMEOUT.json')
            require(marker.get('status')=='TIMEOUT_NOT_ACCEPTED','unknown timeout marker')
            termination_reason='WORKER_HARD_TIMEOUT'
    except Exception as exc:
        error=repr(exc);code=2
        if child.is_alive(): child.terminate();child.join(timeout=15.)
    finally: server.close()
    manifest=output/'manifest.json';log=output/'run.log'
    bridge.shared.write_new(output/'supervision.json',dict(exit_code=code,wall_seconds=time.perf_counter()-started,
        source_inputs_still_match=error is None and len(grants)==2,source_error=error,live_history_grants=grants,
        manifest_sha256=bridge.core.sha(manifest) if manifest.exists() else None,
        log_sha256=bridge.core.sha(log) if log.exists() else None,interrupted=interrupted,termination_reason=termination_reason,
        locked_history_before=history_before,locked_history_after=pool.stats(),accepted=False,restart_allowed=False,
        administrative_protocol=SESSION_PROTOCOL))
    return code,time.perf_counter()-started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='outputs/gpu_acceleration_20261007/event_resume_8h_20261011_v1')
    parser.add_argument('--wall-seconds', type=float, default=28800.)
    parser.add_argument('--max-stages', type=int, default=100000, help='continuation stages after full existing-history verification')
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
            receipt = dict(status='STARTED_8H_RESUME_PENDING_BOOTSTRAP', pid=process.pid, command=command,
                creationflags=flags, source_sha256={rel(Path(__file__)): bridge.core.sha(Path(__file__))})
            write_new(output.with_name(output.name+'.launch.json'), receipt); print(json.dumps(receipt)); return 0
    return run(args)

if __name__ == '__main__':
    raise SystemExit(main())

"""Strict nonlinear-failure continuation inside the existing eight-hour authorization."""
from pathlib import Path
from datetime import datetime,timedelta,timezone
import argparse
import ast
import copy
import inspect
import json
import multiprocessing
import math
import numpy as np
import os
import sys
import time
import traceback
import run_gpu_event_resume_8h_20261011 as previous

ROOT,BASE,CPU_BASE=previous.ROOT,previous.BASE,previous.CPU_BASE
frozen,bridge,admin,adaptive=previous.frozen,previous.bridge,previous.admin,previous.adaptive
v5=previous.v5
require,read,load,rel=previous.require,previous.read,previous.load,previous.rel
source_generation,source_transition=previous.source_generation,previous.source_transition
cost_evidence,measured_envelope=previous.cost_evidence,previous.measured_envelope
cost_decision,dynamic_budget=previous.cost_decision,previous.dynamic_budget
authorize_grant,retain_failure=previous.authorize_grant,previous.retain_failure
multi_event_recovery,timeout_recovery=previous.multi_event_recovery,previous.timeout_recovery
canonical=previous.canonical
LOCATOR_PROTOCOL=previous.LOCATOR_PROTOCOL
BLOCKED_PREFIX_PROTOCOL=previous.BLOCKED_PREFIX_PROTOCOL
VISIBLE_PATH_PROTOCOL=previous.VISIBLE_PATH_PROTOCOL
ROLES=previous.ROLES
AUTHORIZED_AT,ABSOLUTE_STOP=previous.AUTHORIZED_AT,previous.ABSOLUTE_STOP
authorized_wall_seconds=previous.authorized_wall_seconds
PROTOCOL='FROZEN_V5_STRICT_NONLINEAR_RECOVERY_20261011_V1'
SESSION_PROTOCOL=PROTOCOL
CONTRACT=ROOT/'docs/GPU_NONLINEAR_RESUME_20261011.md'
COMPONENTS=BASE/'nonlinear_resume_components_20261011_v1.json'
OLD_CAMPAIGN=BASE/'event_resume_8h_20261011_v1'
PARENT=OLD_CAMPAIGN/'attempt_0039'
FAILED=OLD_CAMPAIGN/'attempt_0040'
PARENT_SHA256='aedf8bfed9c67f38ec12bf38e1e9ec77e37cbba39914f63064dfedbfb537c474'
PARENT_TIME_S=.08958959960937454
PARENT_STEP=599
ADMINISTRATIVE_FILES=('scripts/run_gpu_nonlinear_resume_20261011.py','scripts/afsd_nonlinear_resume_worker_20261011.py',
    'scripts/check_gpu_nonlinear_resume_20261011.py','tests/review/test_gpu_nonlinear_resume_20261011.py',
    'docs/GPU_NONLINEAR_RESUME_20261011.md','scripts/afsd_nonlinear_refinement_20261011.py',
    'tests/review/test_afsd_nonlinear_refinement_20261011.py')


def checks(pool):
    inherited=previous.checks(pool);receipt=read(COMPONENTS)
    require(receipt['status']=='PASS_NONLINEAR_RESUME_COMPONENTS' and receipt['exit_code']==0
        and receipt['physical_source_changed'] is False and receipt['gpu_campaign_runs']==0,'executed administrative-only nonlinear recovery checks required')
    require(receipt['inherited_receipt']==rel(previous.COMPONENTS)
        and receipt['inherited_receipt_sha256']==pool.sha(previous.COMPONENTS),'frozen eight-hour receipt required')
    require(set(ADMINISTRATIVE_FILES).issubset(receipt['source_sha256'])
        and all(receipt['source_sha256'].get(k)==v for k,v in inherited.items()),'frozen eight-hour closure changed')
    require(receipt['physical_source_sha256']==read(PARENT/'manifest.json')['source_sha256']
        and len(receipt['physical_source_sha256'])==636,'exact unchanged physical implementation required')
    require(receipt['parent_directory']==rel(PARENT) and receipt['parent_checkpoint_sha256']==PARENT_SHA256
        and receipt['failed_stage']==rel(FAILED) and receipt['absolute_stop']==ABSOLUTE_STOP.isoformat(),'registered parent/failure/absolute deadline differs')
    pool.seal(receipt['source_sha256']);pool.seal(receipt['artifacts_sha256'])
    return frozen.merge(inherited,receipt['source_sha256'],receipt['artifacts_sha256'],{rel(COMPONENTS):pool.sha(COMPONENTS)})


def static_audit():
    old=(ROOT/'scripts/afsd_event_resume_8h_worker_20261011.py').read_text(encoding='utf-8')
    new=(ROOT/'scripts/afsd_nonlinear_resume_worker_20261011.py').read_text(encoding='utf-8')
    require(ast.dump(ast.parse(old.replace(previous.PROTOCOL,PROTOCOL)))==ast.dump(ast.parse(new)),
        'worker changed beyond explicit administrative protocol identity')
    expected=inspect.getsource(previous.run_stage).replace('afsd_event_resume_8h_worker_20261011','afsd_nonlinear_resume_worker_20261011')
    require(ast.dump(ast.parse(expected))==ast.dump(ast.parse(inspect.getsource(run_stage))),'stage supervision changed beyond worker name')
    expected=inspect.getsource(previous.continue_selected).replace(
        '            if recovery is None:\n                screen=',
        "            if recovery is None:\n                recovery=nonlinear_recovery(pool,child,failure,meta,policy,template['source_sha256'])\n            if recovery is None:\n                screen=")
    require(ast.dump(ast.parse(expected))==ast.dump(ast.parse(inspect.getsource(continue_selected))),
        'ordinary loop changed beyond the explicit nonlinear recovery branch')
    return dict(status='PASS_STRICT_NONLINEAR_ADMINISTRATIVE_BRIDGE',physical_source_changed=False,
        frozen_eight_hour_audit=previous.static_audit(),absolute_stop=ABSOLUTE_STOP.isoformat(),
        new_time_authorization=False,dynamic_cost_and_deadline_policy_unchanged=True)


def stage_protocol(template,meta,parent,dt,administrative,session_record,refs,costs,remaining):
    result=previous.stage_protocol(template,meta,parent,dt,administrative,session_record,refs,costs,remaining)
    result.update(administrative_protocol=PROTOCOL,administrative_checks=rel(COMPONENTS),resume_contract=rel(CONTRACT))
    return result


def verify_parent(meta):
    require(meta['screened_continuation_sha256']==PARENT_SHA256 and meta['screened_end_time_s']==PARENT_TIME_S
        and meta['screened_end_step']==PARENT_STEP,'exact complete old39 checkpoint/time/step required')



def verify_old_audit():
    # Prove this specific final count gap from the preceding sealed snapshot,
    # immutable failed-stage inventory and the exact final byte-read totals.
    baseline=read(OLD_CAMPAIGN/'attempt_0039_receipt.json')['locked_history_statistics']
    supervision=read(FAILED/'supervision.json');before=supervision['locked_history_before'];after=supervision['locked_history_after']
    final=read(OLD_CAMPAIGN/'final_history_audit.json');cost_path=OLD_CAMPAIGN/'attempt_0040_failed_cost.json';cost=read(cost_path)
    names={'manifest.json','protocol.json','run.log','supervision.json','WORKER_PRE_END_FAILURE.json',
        'coarse/proposal/broyden/iterations.json','coarse/proposal/broyden/last_UNCOMMITTED_candidate.json',
        'coarse/proposal/broyden/last_UNCOMMITTED_candidate.npz','coarse/proposal/broyden/last_UNCOMMITTED_candidate_graph.json',
        'coarse/proposal/broyden/last_UNCOMMITTED_candidate_graph.npz','coarse/proposal/broyden/secant_failure.json',
        'coarse/proposal/newton_retry/iterations.json','coarse/proposal/newton_retry/last_UNCOMMITTED_candidate.json',
        'coarse/proposal/newton_retry/last_UNCOMMITTED_candidate.npz'}
    paths={path.relative_to(FAILED).as_posix():path for path in FAILED.rglob('*') if path.is_file()}
    require(set(paths)==names,'exact fourteen-file final failed-stage inventory required')
    hashes={rel(path):bridge.core.sha(path) for path in paths.values()}
    failure=read(FAILED/'manifest.json')
    expected=dict(stage=rel(FAILED),supervision_sha256=hashes[rel(FAILED/'supervision.json')],
        supervised_wall_s=supervision['wall_seconds'],worker_wall_s=failure['wall_seconds'],complete_pair=False,usable_parent=False,
        retained_artifacts_sha256=hashes)
    require(cost==expected,'final failed-cost ledger does not exactly bind all fourteen files')
    require(hashes[rel(FAILED/'manifest.json')]==supervision['manifest_sha256']
        and hashes[rel(FAILED/'run.log')]==supervision['log_sha256'],'old failed supervision identity mismatch')
    base=baseline['locked_files']
    require(type(base) is int and base>0 and baseline['unique_files']==baseline['expected_files']==baseline['hash_reads']==base
        and baseline['failed'] is False and baseline['final_audits']==0,'complete preceding locked snapshot required')
    proto_bytes=(FAILED/'protocol.json').stat().st_size
    for row in (before,after):
        require(row['locked_files']==row['unique_files']==row['hash_reads']==base+1 and row['expected_files']==base
            and row['hash_bytes']==baseline['hash_bytes']+proto_bytes and row['failed'] is False and row['final_audits']==0,
            'only failed protocol may be newly hashed before retaining failure')
    retained_bytes=sum(path.stat().st_size for name,path in paths.items() if name!='protocol.json')
    cost_bytes=cost_path.stat().st_size;held=base+len(names)+1;expected_count=base+len(names)
    expected_bytes=2*(after['hash_bytes']+retained_bytes+cost_bytes)
    require(final['failed'] is False and final['final_audits']==1 and final['unique_files']==final['locked_files']==held
        and final['expected_files']==expected_count and final['hash_reads']==2*held and final['hash_bytes']==expected_bytes
        and final['authority']=='LIVE_WINDOWS_DENY_WRITE_DELETE_HANDLES_V1','specific final failed-cost audit arithmetic differs')
    digest=bridge.core.sha(cost_path)
    require(digest=='6a0e7bcf4f97d41cb0563d66340a64c2cd3988eeb159f6888abe32965e053a3a','exact old final failed-cost SHA required')
    return dict(status='PASS_EXACT_FINAL_ADDITIONAL_HASHED_COST_FILE',additional_hashed_file=rel(cost_path),sha256=digest,
        preceding_locked_files=base,failed_stage_files=len(names),final_locked_files=held,final_expected_files=expected_count,
        expected_final_hash_reads=2*held,expected_final_hash_bytes=expected_bytes,
        final_audit_covers_all_held_files=True,new_session_explicit_seal_required=True)


def validate_candidate_arrays(arrays,report,parent_pa,parent_state,dt):
    ids=np.asarray(parent_pa.particle_id);n=len(ids);X=np.column_stack((parent_pa.x,parent_pa.y,parent_pa.z))
    shapes=dict(particle_id=(n,),position_m=(n,3),displacement_m=(n,3),velocity_m_s=(n,3),density_kg_m3=(n,),
        stress_Pa=(n,3,3),temperature_K=(n,),plastic_strain=(n,),force_N=(n,3),internal_force_N=(n,3),
        contact_force_N=(n,3),contact_history=(n,8,3))
    require(set(arrays)==set(shapes)|{'accepted','restart_allowed'},'exact uncommitted diagnostic array inventory required')
    for name,shape in shapes.items():
        require(arrays[name].shape==shape and np.all(np.isfinite(arrays[name])),'nonfinite or malformed uncommitted arrays')
    require(np.array_equal(arrays['particle_id'],ids),'candidate persistent IDs differ from complete parent')
    require(all(arrays[name].shape==() and arrays[name].dtype==np.dtype(bool) and arrays[name].item() is False
        for name in ('accepted','restart_allowed')),'uncommitted state must remain unaccepted')
    require(np.all(arrays['density_kg_m3']>0) and np.all(arrays['temperature_K']>0) and np.all(arrays['plastic_strain']>=0),
        'candidate material fields outside finite positive domain')
    require(np.array_equal(arrays['position_m'],X+arrays['displacement_m'])
        and np.array_equal(arrays['velocity_m_s'],arrays['displacement_m']/dt),'candidate position/velocity do not match exact parent and dt')
    require(np.array_equal(arrays['force_N'],arrays['internal_force_N']+arrays['contact_force_N'])
        and np.array_equal(arrays['contact_force_N'],np.asarray(report['contact']['force_N'])),'candidate force decomposition differs')
    free=np.asarray(parent_state['free']);require(free.dtype==np.dtype(bool) and free.shape==(n,3) and np.any(free),'complete parent free mask required')
    limit=1e-5+1e-7*max(1.,float(np.max(np.abs(arrays['internal_force_N']))),float(np.max(np.abs(arrays['contact_force_N']))))
    residual=float(np.max(np.abs(arrays['force_N'][free])));metrics=report['metrics']
    require(metrics['force_limit_N']==limit and metrics['force_residual_N']==residual and metrics['scaled_residual']==residual/limit,
        'original force threshold/residual not independently reproduced')
    require(residual>limit and report['contact']['dt_s']==dt,'only unconverged original same-dt force residual eligible')
    return dict(force_limit_N=limit,force_residual_N=residual,scaled_residual=residual/limit,particle_count=n)


def nonlinear_recovery(pool,folder,failure,parent_meta,policy,source):
    import afsd_nonlinear_refinement_20261011 as nonlinear
    if failure.get('error')!="RefinableCandidateError('complete residual line search/evaluation budget failed')":return None
    if (folder/'pair_screen.json').exists() or (folder/'fine_1').exists() or (folder/'fine_2').exists():return None
    candidate=folder/'coarse/proposal'
    required={'secant_failure':candidate/'broyden/secant_failure.json','broyden_iterations':candidate/'broyden/iterations.json',
        'broyden_candidate':candidate/'broyden/last_UNCOMMITTED_candidate.json',
        'broyden_graph':candidate/'broyden/last_UNCOMMITTED_candidate_graph.json',
        'newton_iterations':candidate/'newton_retry/iterations.json','newton_candidate':candidate/'newton_retry/last_UNCOMMITTED_candidate.json'}
    if not all(rel(path) in failure.get('artifacts_sha256',{}) for path in required.values()):return None
    require(not any(folder.rglob('NOT_RESTARTABLE_candidate.npz')) and not any(folder.rglob('candidate_qualification.json'))
        and not any(folder.rglob('SOLVED_UNCOMMITTED_candidate.npz')),'unsupported solved/qualified or nested numerical failure layout')
    supervision=read(folder/'supervision.json')
    require(supervision['exit_code']!=0 and supervision.get('source_error') is None
        and pool.sha(folder/'manifest.json')==supervision['manifest_sha256'] and pool.sha(folder/'run.log')==supervision['log_sha256']
        and read(folder/'manifest.json')==failure,'nonlinear failed receipt identity mismatch')
    require(failure['source_sha256']==source==parent_meta['source_sha256'] and failure['parent_implementation_source_sha256']==source
        and failure.get('implementation_transition') is None and failure['parent_sha256']==parent_meta['screened_continuation_sha256']
        and failure['parent_checkpoint']==parent_meta['screened_continuation_checkpoint']
        and failure['dt_s']==policy['pending']['dt_s'] and policy['parent_sha256']==parent_meta['screened_continuation_sha256'],
        'nonlinear source/full parent/dt differs')
    pool.seal(failure['source_sha256']);pool.seal(failure['reference_sha256']);pool.seal(failure['artifacts_sha256'])
    def evidence(path):
        require(rel(path) in failure['artifacts_sha256'] and pool.sha(path)==failure['artifacts_sha256'][rel(path)],'nonlinear evidence missing or SHA changed')
        return path
    reports={key:read(evidence(path)) for key,path in required.items()}
    diagnosis=nonlinear.recognize(failure,supervision,reports)
    if not diagnosis['eligible']:return None
    parent_path=ROOT/parent_meta['screened_continuation_checkpoint']
    require(pool.sha(parent_path)==parent_meta['screened_continuation_sha256'],'actual complete parent changed')
    pa,state=load(parent_path)
    require(state['implementation_source_sha256']==source,'loaded complete parent implementation differs')
    array_evidence={};saved={}
    for key,directory in (('broyden_candidate',candidate/'broyden'),('newton_candidate',candidate/'newton_retry')):
        path=evidence(directory/'last_UNCOMMITTED_candidate.npz')
        with np.load(path,allow_pickle=False) as npz: arrays={name:npz[name].copy() for name in npz.files}
        array_evidence[key]=validate_candidate_arrays(arrays,reports[key],pa,state,failure['dt_s']);saved[key]=arrays
    path=evidence(candidate/'broyden/last_UNCOMMITTED_candidate_graph.npz')
    with np.load(path,allow_pickle=False) as npz: graph={name:npz[name].copy() for name in npz.files}
    X=np.column_stack((pa.x,pa.y,pa.z))
    require(np.array_equal(graph['particle_id'],pa.particle_id) and np.array_equal(graph['position0'],X)
        and np.array_equal(graph['position1'],saved['broyden_candidate']['position_m'])
        and all(np.all(np.isfinite(value)) for value in graph.values()),'graph diagnostic does not match same-parent finite candidate')
    require(np.array_equal(graph['i'],graph['raw_i']) and np.array_equal(graph['j'],graph['raw_j'])
        and np.array_equal(graph['g'],graph['raw_g']) and graph['blocked_raw_edge_indices'].size==0,
        'graph evidence contains unrecognized visibility filtering')
    from visible_quasistatic import edge_fingerprint
    require(edge_fingerprint(graph['i'],graph['j'],pa.particle_id)==reports['broyden_graph']['graph_fingerprint'],'actual candidate graph fingerprint differs')
    diagnosis['independent_array_evidence']=array_evidence
    result,event=nonlinear.refine(policy,diagnosis)
    event['candidate_role']='coarse';event['proposal_directory']=rel(candidate)
    return result,event


def restore_history(pool,args):
    restored=previous.restore_history(pool,args)
    template=restored['template'];refs=restored['refs'];parent=restored['parent'];meta=restored['meta']
    costs=restored['costs'];policy=restored['policy'];state=restored['state'];frames=restored['frames']
    failed_costs=restored['previous_failed_costs']
    session=read(OLD_CAMPAIGN/'protocol.json');summary=read(OLD_CAMPAIGN/'final_summary.json')
    require(session['schema']==previous.PROTOCOL and session['physical_source_sha256']==template['source_sha256']
        and session['restored_adaptive_state']==policy,'old eight-hour source and initial policy differ')
    require(summary['finished'] is True and summary['exit_code']==2 and summary['status']=='STOPPED_CHILD_FAILURE'
        and summary['completed_stage_count']==35 and summary['attempt_count']==40
        and summary['latest_pass_stage']==rel(PARENT) and summary['failed_attempt']==rel(FAILED)
        and summary['checkpoint_sha256']==PARENT_SHA256 and summary['physical_time_s']==PARENT_TIME_S
        and summary['screened_steps']==PARENT_STEP,'exact old eight-hour failed40 termination required')
    audit_evidence=verify_old_audit()
    for name in ('protocol.json','final_summary.json','final_history_audit.json','deposition_progress.json'):
        path=OLD_CAMPAIGN/name;refs=frozen.merge(refs,{rel(path):pool.sha(path)})
    last_pending=None;policy_before_failure=None;last_failure=None
    for number in range(1,41):
        child=OLD_CAMPAIGN/f'attempt_{number:04d}';receipt_path=OLD_CAMPAIGN/f'attempt_{number:04d}_receipt.json'
        receipt=read(receipt_path);protocol=read(child/'protocol.json')
        require(protocol['source_sha256']==template['source_sha256'] and protocol['parent_sha256']==meta['screened_continuation_sha256']
            and protocol['parent_gpu_directory']==rel(parent) and receipt['parent']==rel(parent)
            and receipt['dt_s']==protocol['dt_s'] and receipt['stage']==rel(child),'old stage sequence identity differs')
        ready=copy.deepcopy(policy)
        policy=adaptive.begin_attempt(policy,protocol['dt_s'],phase_clipped=protocol['dt_s']<policy['nominal_dt_s']*(1.-1e-12))
        if receipt['exit_code']==0:
            require(number not in (1,2,4,6,40),'registered failed attempt changed to PASS')
            new=frozen.validated_stage(pool,child,template['source_sha256'],meta)
            require(source_generation(new,template)==6 and new['event_tolerance_scale']==.25
                and new['parent_gpu_directory']==rel(parent) and new.get('implementation_transition') is None,'unchanged old quarter-scale history required')
            policy,event=adaptive.resolve_pass(policy,new['pair_screen'],state['parameters']['dx'],state['bore_observation']['particle_ids'],
                new['screened_continuation_sha256'],original_candidate_gates_passed=True)
            require(receipt.get('verified') is True,'original PASS receipt not verified')
            refs=frozen.bind_stage(pool,child,new,refs);costs.append(cost_evidence(pool,child,new));parent,meta=child,new
            observation,state=frozen.observe(meta,frames)
        else:
            require(number in (1,2,4,6,40),'unexpected old failed attempt')
            refs,failure,resource_cost=retain_failure(pool,child,protocol,meta,refs)
            cost_path=OLD_CAMPAIGN/f'attempt_{number:04d}_failed_cost.json'
            require(canonical(resource_cost)==read(cost_path),'retained failed cost or complete artifact map differs')
            refs=frozen.merge(refs,{rel(cost_path):pool.sha(cost_path)});pool.seal({rel(cost_path):pool.sha(cost_path)})
            failed_costs.append(resource_cost)
            if number==40:
                verify_parent(meta);require(parent.resolve()==PARENT.resolve(),'last complete parent directory differs')
                require(protocol['event_cost_envelope']==canonical(measured_envelope(costs)),'all 59 complete pair costs differ')
                last_pending=copy.deepcopy(policy);policy_before_failure=ready;last_failure=failure
                policy,event=adaptive.resolve_failure(policy,failure,None,state['parameters']['dx'],state['bore_observation']['particle_ids'],
                    original_candidate_gates_passed=False)
                require(policy==summary['adaptive_state'] and canonical(event)==summary['last_adaptive_event'],
                    'replayed final STOPPED policy or event differs')
            else:
                recovery=multi_event_recovery(pool,child,failure,meta,policy,template['source_sha256'])
                require(recovery is not None and recovery[1]['action']=='REFINE','historical multi-event retry no longer qualifies')
                policy,event=recovery
        require(canonical(event)==receipt['adaptive_event'],'historical adaptive receipt differs')
        refs=frozen.merge(refs,{rel(receipt_path):pool.sha(receipt_path)})
    require(canonical(observation)==read(OLD_CAMPAIGN/'deposition_progress.json'),'all 120 selected fine-leaf observations differ')
    require(len(costs)==59 and len(frames)==120 and len(failed_costs)==6,'complete historical cost/frame inventory differs')
    require(state['implementation_source_sha256']==template['source_sha256']==meta['source_sha256'],'complete parent physical source differs')
    recovery=nonlinear_recovery(pool,FAILED,last_failure,meta,last_pending,template['source_sha256'])
    require(recovery is not None and recovery[1]['action']=='REFINE','old numerical rejection does not qualify strict same-parent refinement')
    policy,event=recovery
    require(policy['parent_sha256']==PARENT_SHA256 and policy['status']=='READY' and policy['pending'] is None,
        'new numerical recovery must retain the full old39 parent')
    pool.seal(refs)
    restored.update(refs=refs,parent=parent,meta=meta,costs=costs,policy=policy,state=state,frames=frames,observation=observation,
        previous_failed_costs=failed_costs,initial_recovery=event,policy_before_failed_attempt=policy_before_failure,
        original_final_audit_evidence=audit_evidence)
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
                recovery=nonlinear_recovery(pool,child,failure,meta,policy,template['source_sha256'])
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
            old_failure=rel(FAILED),initial_nonlinear_recovery=restored['initial_recovery'],
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
    from afsd_nonlinear_resume_worker_20261011 import stage_worker
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
    parser.add_argument('--output', default='outputs/gpu_acceleration_20261007/nonlinear_resume_20261011_v1')
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
            receipt = dict(status='STARTED_NONLINEAR_RESUME_PENDING_BOOTSTRAP', pid=process.pid, command=command,
                creationflags=flags, source_sha256={rel(Path(__file__)): bridge.core.sha(Path(__file__))})
            write_new(output.with_name(output.name+'.launch.json'), receipt); print(json.dumps(receipt)); return 0
    return run(args)

if __name__=='__main__':
    raise SystemExit(main())

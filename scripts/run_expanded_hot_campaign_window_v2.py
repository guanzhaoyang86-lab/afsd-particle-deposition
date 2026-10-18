"""Explicit single-window expanded process with immutable budget registrations.

--plan-only registers a coarse protocol without physics. Each later invocation
runs one stage only. All physical work stays in the frozen candidate helper.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import run_visible_support_process as old
import expanded_hot_continuation_candidate as qualified
import expanded_hot_domain as model
from qualify_expanded_hot_domain import immutable_hashes

core,post,event=old.core,old.post,old.event
ROOT=core.ROOT
CONTRACT=ROOT/'docs/EXPANDED_HOT_CAMPAIGN_WINDOW_V2_CONTRACT_20261005.md'
QUALIFICATION=ROOT/'outputs/expanded_hot_qualification_20261005/v1/summary.json'
BRANCH_ORIGIN_SHA='31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3'
COMPONENTS=ROOT/'outputs/visible_top_components_20261004/v1/summary.json'
COMPONENTS_SHA='7fbb5fdca371335ff528f4132fb840948dc3239de8eb50e420b3533014505cab'
PATH_COMPONENTS=ROOT/'outputs/visible_path_tight_20261004/qualification_v1/summary.json'
PATH_COMPONENTS_SHA='3d461593a2fab30900a83bcde5b2835b353bec3d66783393b3500ce493b750bc'
PATH_QUALIFICATION='FIXED_RETAINED_GRAPH_RODRIGUES_CANCELLATION_ENCLOSURE_V1'
DT_CHOICES={'1.953125':1.953125e-6,'3.90625':3.90625e-6,'7.8125':7.8125e-6,'15.625':15.625e-6,'31.25':31.25e-6,'62.5':62.5e-6,'125':125e-6}
PASS_COARSE='QUALIFIED_EXPANDED_CAMPAIGN_V2_COARSE_NOT_ACCEPTED'
PASS_FINE='PASS_EXPANDED_CAMPAIGN_V2_SELECTED_FINE_SCOPE_NOT_PHYSICAL_VALIDATION'
PARENT_AUDITS={
 'PASS_EXPANDED_HOT_BUDGET_PREFIX_READBACK_NOT_PHYSICAL_VALIDATION':'scripts/readback_expanded_hot_budget_prefix.py',
 'PASS_EXPANDED_HOT_CAMPAIGN_READBACK_NOT_PHYSICAL_VALIDATION':'scripts/readback_expanded_hot_campaign_window.py',
 'PASS_EXPANDED_HOT_CAMPAIGN_V2_READBACK_NOT_PHYSICAL_VALIDATION':'scripts/readback_expanded_hot_campaign_window_v2.py'}
PROTOCOL='EXPANDED_HOT_SINGLE_WINDOW_PROTOCOL_V2'
SELECTION='EXPANDED_HOT_FINE_PREFIX_SELECTION_V2'
TIME_STUDY=ROOT/'outputs/contact_onset_time_study_20261005/v1/readback_eight/summary.json'
TIME_STUDY_SHA='a0a9c0b99c86c6f2b02182bbc5320f87231296464d1afe1539515edde137a67f'
TIME_READER=ROOT/'scripts/readback_contact_onset_time_study.py'
TIME_READER_SHA='759fad998a9e51f59690864d2558ae63fe3b09ce69ff451358af4a3610b5a92b'
require=old.require


def relative(path):return Path(path).resolve().relative_to(ROOT).as_posix()

def project(path):
    path=Path(path);path=(ROOT/path).resolve() if not path.is_absolute() else path.resolve()
    require(path.is_relative_to(ROOT),'path escaped project');return path

def write_new(path,value):
    with Path(path).open('x',encoding='utf-8') as stream:
        json.dump(core.clean(value),stream,indent=2,allow_nan=False);stream.write('\n')

def bind_map(frozen,mapping):
    for name,digest in mapping.items():
        path=project(name);require(core.sha(path)==digest,'bound source/input drift: '+name)
        require(name not in frozen or frozen[name]==digest,'conflicting source binding '+name);frozen[name]=digest

def time_resolution_qualification():
    """Fixed saved study receipt, never a diagnostic restart or physics call."""
    require(core.sha(TIME_STUDY)==TIME_STUDY_SHA and core.sha(TIME_READER)==TIME_READER_SHA,
        'bounded temporal-study reader/receipt changed')
    q=post.read(TIME_STUDY)
    require(q['status']=='PASS_CONTACT_ONSET_TIME_STUDY_READBACK_NOT_PHYSICAL_VALIDATION' and q['exit_code']==0 and
        q['physical_solver_evaluations']==q['nnps_evaluations']==q['material_evaluations']==q['contact_evaluations']==0 and
        q['accepted'] is False and q['restart_allowed'] is False and q['new_accepted_steps']==0 and q['accepted_time_advanced_s']==0. and
        q['actual_accepted_parent_sha256']=='12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b' and
        q['actual_accepted_global_step']==14 and abs(q['actual_accepted_time_s']-62.5e-6)<1e-15 and
        q['evidence_sha256'].get(relative(TIME_READER))==TIME_READER_SHA,
        'independent bounded temporal-study scope or true parent differs')
    inputs=dict(q['evidence_sha256']);inputs[relative(TIME_STUDY)]=TIME_STUDY_SHA
    bind_map({},inputs);measurements=[]
    require([s['stage'] for s in q['stages']]==['four','eight'],'study resolution sequence differs')
    for stage,count in zip(q['stages'],(4,8)):
        require(stage['physical_trajectory_completed'] and stage['exit_code']==0 and stage['candidate_count']==count and
            stage['accepted_steps']==0 and stage['accepted_time_advanced_s']==0. and stage['full_prefix_preserved'] and
            stage['diagnostic_global_step']==14+count and abs(stage['diagnostic_candidate_end_time_s']-70.3125e-6)<1e-15 and
            len(stage['screens'])==count//2,'study trajectory incomplete or retrospectively accepted')
        folder=TIME_STUDY.parent.parent/stage['stage']
        meta=post.read(folder/'manifest.json');rows=post.read(folder/'candidate_history.json')
        require(len(rows)==count and core.sha(folder/'manifest.json')==stage['manifest_sha256'] and
            meta['wall_seconds']==stage['wall_seconds'] and
            [r['wall_seconds'] for r in rows]==[r['wall_seconds'] for r in stage['candidate_summaries']],
            'study complete candidate costs differ from independent readback')
        costs=[r['wall_seconds'] for r in rows]
        require(all(np.isfinite(v) and v>0 for v in costs),'invalid study cost')
        measurements.append(dict(stage='contact_onset_'+stage['stage'],directory=relative(folder),
            manifest_sha256=core.sha(folder/'manifest.json'),candidate_history_sha256=core.sha(folder/'candidate_history.json'),
            complete_candidate_count=count,complete_candidate_wall_seconds=costs,maximum_complete_candidate_s=max(costs),
            stage_wall_seconds=meta['wall_seconds'],stage_wall_per_complete_candidate_s=meta['wall_seconds']/count,
            terminal_status=meta['status']))
    require(q['stages'][0]['temporal_screens_all_pass'] is False and
        q['stages'][1]['temporal_screens_all_pass'] is True and all(s['passed'] for s in q['stages'][1]['screens']),
        'original four-step failure or passing eight-versus-four screen differs')
    metadata=dict(qualification='BOUNDED_CONTACT_ONSET_TIME_SCREEN_NOT_GLOBAL_DT_VALIDATION_V1',
        receipt=relative(TIME_STUDY),receipt_sha256=TIME_STUDY_SHA,reader_sha256=TIME_READER_SHA,
        study_parent_sha256=q['actual_accepted_parent_sha256'],study_start_s=q['actual_accepted_time_s'],
        study_end_s=q['stages'][1]['diagnostic_candidate_end_time_s'],
        passing_coarse_dt_s=1.953125e-6,passing_fine_dt_s=.9765625e-6,
        diagnostic_states_remain_unaccepted=True,whole_history_convergence_established=False,physics_changed=False)
    return metadata,measurements,inputs


def parent_context(parent_path,audit_path):
    """Saved evidence, state arrays and hashes only; no NNPS or residual."""
    parent_path,audit_path=project(parent_path),project(audit_path)
    require(parent_path.name=='accepted.npz' and parent_path.parent.parent.name=='fine','only an actual fine publication can be a parent')
    pa,state=core.restore(parent_path);audit=post.read(audit_path);inputs={}
    def bind(path):
        path=project(path);inputs[relative(path)]=core.sha(path);return path
    bind(parent_path);bind(audit_path)
    require(audit['status'] in PARENT_AUDITS and audit['exit_code']==0 and
        audit['physical_solver_evaluations']==0 and audit['nnps_evaluations']==0,'eligible independent saved-evidence readback required')
    reader=ROOT/PARENT_AUDITS[audit['status']]
    require(audit['evidence_sha256'].get(relative(reader))==core.sha(reader),'parent reader implementation binding missing')
    bind_map({},audit['evidence_sha256'])
    p,provenance=model.parameters(.008)
    require(state['classification']==model.CLASSIFICATION and state['model_scope']==model.SCOPE and
        state['operator_variant']=='visible' and state['independent_initialization'] is True and
        state['accepted'] is True and state['restart_allowed'] is True and len(pa.x)==918 and
        int(state['free'].sum())==1961 and state['parameters']==p and state['parameter_provenance']==provenance and
        state['source_sha256']==model.source_hashes() and 'spatial_model_transition' not in state and
        state['branch_origin_sha256']==BRANCH_ORIGIN_SHA,'wrong independent accepted physical model/origin')
    for name,prop in pa.properties.items():require(np.isfinite(prop.get_npy_array()).all(),'nonfinite parent field: '+name)
    require(audit['accepted_history'] and audit['accepted_history'][-1]['checkpoint']==relative(parent_path) and
        audit['accepted_history'][-1]['sha256']==core.sha(parent_path) and audit['accepted_time_s']==state['time_s'] and
        audit['accepted_global_steps']==state['step'] and audit['branch_origin_sha256']==BRANCH_ORIGIN_SHA and
        audit['window_origin_sha256']==state['window_origin_sha256'],'readback is not for this exact final accepted endpoint')
    window=project(audit['window_directory']);fine=window/'fine'
    require(parent_path.parent.parent==fine,'parent and readback window differ')
    published=post.read(bind(fine/'accepted_history.json'));local_index=int(parent_path.parent.name.removeprefix('step_'))
    require(local_index==len(published)==len(audit['accepted_history']) and local_index%2==0 and
        published[-1]['sha256']==core.sha(parent_path) and published[-1]['global_step']==state['step'] and
        published[-1]['time_s']==state['time_s'],'parent accepted history is incomplete or not a complete pair')
    previous=project(audit['parent_checkpoint']);previous_hash=audit['parent_sha256']
    for index,row in enumerate(published,1):
        path=fine/row['checkpoint'];bind(path)
        require(path==fine/f'step_{index:03d}/accepted.npz' and row['sha256']==core.sha(path) and
            row['parent_sha256']==previous_hash and audit['accepted_history'][index-1]['sha256']==row['sha256'],
            'parent published direct ancestry is inconsistent')
        previous_hash=row['sha256']
    dt=float(published[-1]['dt_s']);require(np.isfinite(dt) and dt>0,'invalid last accepted fine step')
    screens=[row for row in audit['pair_screens'] if row['sha256']==state['window_temporal_screen_sha256']]
    require(len(screens)==1 and screens[0]['passed'] is True and screens[0]['time_s']==state['time_s'] and
        screens[0]['pair']==local_index//2,'last accepted endpoint lacks its actual passing pair screen')
    screen=post.read(bind(fine/f'pair_{local_index//2:03d}_temporal_screen.json'))
    require(core.sha(fine/f'pair_{local_index//2:03d}_temporal_screen.json')==screens[0]['sha256'] and
        screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'] and
        screen['time_s']==state['time_s'] and screen['cumulative_event_and_blocked_ids']==screens[0]['cumulative_event_and_blocked_ids'],
        'saved parent passing screen differs')
    # Later failed screens are retained as evidence, never inherited as seeds.
    seeds=list(screens[0]['cumulative_event_and_blocked_ids'])
    model.discrete_process.validate_interval_energy_state(state['interval_energy'])
    model.discrete_process.validate_bore_crossing_state(state['bore_observation'])
    require(state['interval_energy']['steps']==state['step'] and len(state['interval_energy']['records'])==state['step'],
        'parent lost part of the global energy prefix')
    ledger=bind(audit['trajectory']['complete_interval_energy_path'])
    require(core.sha(ledger)==audit['trajectory']['complete_interval_energy_sha256'] and
        post.read(ledger)==core.clean(state['interval_energy']),'full readback energy ledger differs from parent')
    metadata=state['expanded_initial_qualification'];q=post.read(bind(QUALIFICATION))
    require(metadata==audit['independent_initial_qualification'] and metadata['initialization_sha256']==BRANCH_ORIGIN_SHA and
        project(metadata['receipt'])==QUALIFICATION and metadata['receipt_sha256']==core.sha(QUALIFICATION) and
        q['qualification_passed'] is True and q['initialization_sha256']==BRANCH_ORIGIN_SHA and
        q['parameters']==p and q['parameter_provenance']==provenance and immutable_hashes(pa,state)==q['immutable_array_sha256'],
        'initial rank/map/identity qualification changed')
    bind_map({},q['source_and_input_sha256']);bind_map({},q['output_sha256'])
    require(state['controller_source_sha256']==state['implementation_source_sha256'],'parent source maps differ')
    frozen={};bind_map(frozen,state['implementation_source_sha256'])
    measurements=[]
    coarse=project(audit['budget_selection']['original_coarse_directory']) if 'budget_selection' in audit else window/'coarse'
    for name,folder in (('coarse',coarse),('fine',fine)):
        manifest=post.read(bind(folder/'manifest.json'));rows=post.read(bind(folder/'candidate_history.json'))
        receipt=post.read(bind(manifest['supervision_receipt']))
        require(manifest['status']!='RUNNING' and receipt['external_wall_s']<=605 and
            manifest['candidate_steps']==len(rows) and manifest['wall_seconds']>0 and len(rows)>0,
            'complete terminal parent candidate measurements required')
        costs=[row['wall_seconds'] for row in rows]
        require(all(np.isfinite(v) and v>0. for v in costs),'invalid measured complete candidate cost')
        measurements.append(dict(stage=name,directory=relative(folder),manifest_sha256=core.sha(folder/'manifest.json'),
            candidate_history_sha256=core.sha(folder/'candidate_history.json'),complete_candidate_count=len(rows),
            complete_candidate_wall_seconds=costs,maximum_complete_candidate_s=max(costs),
            stage_wall_seconds=manifest['wall_seconds'],stage_wall_per_complete_candidate_s=manifest['wall_seconds']/len(rows),
            terminal_status=manifest['status']))
    temporal,study_costs,study_inputs=time_resolution_qualification()
    measurements.extend(study_costs);bind_map(frozen,study_inputs)
    extras=[Path(__file__).resolve(),CONTRACT,Path(qualified.__file__),ROOT/'scripts/expanded_hot_continuation_acceptance.py']
    bind_map(frozen,{relative(path):core.sha(path) for path in extras});bind_map(frozen,inputs)
    outcome={key:audit.get(key) for key in ('status','window_completed','budget_prefix_completed',
        'registered_window_completed','selected_fine_scope_completed','window_outcome','accepted_time_s','accepted_global_steps')}
    outcome['stage_statuses']={key:value['status'] for key,value in audit['stages'].items()}
    return dict(pa=pa,state=state,p=p,provenance=provenance,parent=parent_path,parent_hash=core.sha(parent_path),
        audit=audit_path,audit_hash=core.sha(audit_path),parent_audit=audit,parent_outcome=outcome,
        parent_evidence=audit['evidence_sha256'],frozen=frozen,measurements=measurements,
        last_fine_dt=dt,seeds=seeds,prefix=state['interval_energy']['records'],time_resolution_qualification=temporal)


def measured_cost(measurements):
    return max(max(row['maximum_complete_candidate_s'],row['stage_wall_per_complete_candidate_s']) for row in measurements)


def affordable_pairs(maximum,per_step):
    return max((n for n in range(1,maximum+1) if 30+1.5*(2*n)*per_step<=450.),default=0)


def protocol_value(parent,audit,dt_key,requested,window):
    require(dt_key in DT_CHOICES and 1<=requested<=8,'registered coarse dt and 1..8 pair upper bound required')
    c=parent_context(parent,audit);dt=DT_CHOICES[dt_key]
    require(dt<=4*c['last_fine_dt']*(1+16*np.finfo(float).eps),'coarse dt growth exceeds four previous accepted fine steps')
    cost=measured_cost(c['measurements']);per_step=cost+1.;count=affordable_pairs(requested,per_step)
    state=c['state'];start=state['time_s']
    budget=dict(formula='30 + 1.5 * planned_candidate_count * per_candidate_cost_s',forecast_limit_s=450.,
        measurements=c['measurements'],measured_cost_s=cost,extra_allowance_s=1.,per_candidate_cost_s=per_step,
        requested_pair_count=requested,selected_pair_count=count,coarse_forecast_s=30+1.5*count*per_step,
        fine_forecast_s=30+1.5*(2*count)*per_step)
    return dict(schema=PROTOCOL,registered_before_coarse_worker=True,
        status='REGISTERED_SINGLE_WINDOW' if count else 'BUDGET_BLOCKED_NO_COARSE_RUN',
        window_directory=relative(window),contract=relative(CONTRACT),contract_sha256=core.sha(CONTRACT),
        source_sha256=c['frozen'],parent_evidence_sha256=c['parent_evidence'],parent_outcome=c['parent_outcome'],
        time_resolution_qualification=c['time_resolution_qualification'],
        parent_checkpoint=relative(c['parent']),parent_sha256=c['parent_hash'],parent_readback=relative(c['audit']),
        parent_readback_sha256=c['audit_hash'],branch_origin_sha256=BRANCH_ORIGIN_SHA,window_origin_sha256=c['parent_hash'],
        parent_time_s=start,parent_global_step=state['step'],parent_last_accepted_fine_dt_s=c['last_fine_dt'],
        independent_qualification=relative(QUALIFICATION),independent_qualification_sha256=core.sha(QUALIFICATION),
        parameters=c['p'],parameter_provenance=c['provenance'],inherited_seed_ids=c['seeds'],
        coarse_dt_us=dt_key,coarse_dt_s=dt,fine_dt_s=.5*dt,requested_pair_count=requested,registered_pair_count=count,
        requested_end_s=start+requested*dt,registered_end_s=start+count*dt,budget=budget,
        internal_limit_s=540.,external_limit_s=600.,maximum_coarse_to_previous_fine_ratio=4.,
        accepted=False,restart_allowed=False,physical_solver_evaluations=0,nnps_evaluations=0)


def load_context(window):
    protocol_path=window/'window_protocol.json';protocol=post.read(protocol_path)
    require(protocol['schema']==PROTOCOL and protocol['status']=='REGISTERED_SINGLE_WINDOW' and
        protocol['window_directory']==relative(window),'invalid/unaffordable immutable window protocol')
    expected=protocol_value(protocol['parent_checkpoint'],protocol['parent_readback'],protocol['coarse_dt_us'],
        protocol['requested_pair_count'],window)
    require(expected==protocol,'protocol differs from actual parent/source/budget evidence')
    c=parent_context(protocol['parent_checkpoint'],protocol['parent_readback'])
    c.update(protocol=protocol,protocol_path=protocol_path,protocol_sha256=core.sha(protocol_path),window=window,
        qualification=QUALIFICATION,qualification_hash=core.sha(QUALIFICATION),refs={},
        start=protocol['parent_time_s'],start_step=protocol['parent_global_step'],branch_origin=BRANCH_ORIGIN_SHA,
        dt=protocol['coarse_dt_s'],half=protocol['fine_dt_s'],registered_pairs=protocol['registered_pair_count'])
    return c


def coarse_references(c):
    folder=c['window']/'coarse';meta=post.read(folder/'manifest.json');rows=post.read(folder/'candidate_history.json')
    supervision=post.read(project(meta['supervision_receipt']));refs={}
    require(meta['status'] in (PASS_COARSE,'FAILED_NO_NEW_ACCEPTED_STATE') and meta['source_sha256']==c['frozen'] and
        meta['parent_sha256']==c['parent_hash'] and meta['window_protocol_sha256']==c['protocol_sha256'] and
        meta['accepted_steps']==0 and meta['candidate_steps']==len(rows) and
        supervision['exit_code']==(0 if meta['status']==PASS_COARSE else 2), 'terminal matching coarse stage required')
    require(not list(folder.glob('step_*/accepted.npz')) and not post.read(folder/'accepted_history.json'),'coarse became accepted')
    require(0<=len(rows)<=c['registered_pairs'],'too many coarse candidates')
    if meta['status']==PASS_COARSE:require(len(rows)==c['registered_pairs'],'successful coarse stage incomplete')
    previous=c['parent_hash'];prefix=c['prefix']
    for index,row in enumerate(rows,1):
        require(row['step']==index and row['parent_sha256']==previous and row['dt_s']==c['dt'] and
            row['global_step']==c['start_step']+index and abs(row['time_s']-(c['start']+index*c['dt']))<1e-13,
            'coarse reference ancestry/time differs')
        for name,digest in row['artifacts'].items():
            path=(folder/name).resolve();require(path.is_relative_to((folder/f'step_{index:03d}').resolve()) and
                core.sha(path)==digest,'coarse artifact hash/path differs')
        path=(folder/row['checkpoint']).resolve()
        require(path==folder/f'step_{index:03d}/NOT_RESTARTABLE_candidate.npz' and core.sha(path)==row['sha256'],
            'coarse checkpoint path/hash differs')
        saved=core.load(str(path));state=saved['solver_data']
        model.discrete_process.validate_interval_energy_state(state['interval_energy'])
        require(state['accepted'] is False and state['restart_allowed'] is False and
            state['implementation_source_sha256']==c['frozen'] and state['classification']==model.CLASSIFICATION and
            state['interval_energy']['records'][:-1]==prefix and state['time_s']==row['time_s'],
            'coarse reference source/field/history gate differs')
        previous=row['sha256'];prefix=state['interval_energy']['records']
    # Include failed/unregistered artifacts as immutable evidence, not references.
    for path in sorted(folder.rglob('*')):
        if path.is_file():refs[relative(path)]=core.sha(path)
    for path in (project(meta['supervision_receipt']),c['window']/'execution/coarse.log'):
        refs[relative(path)]=core.sha(path)
    return meta,rows,refs


def selection_value(c):
    meta,rows,refs=coarse_references(c);original_cost=c['protocol']['budget']['measured_cost_s']
    costs=[r['wall_seconds'] for r in rows]
    if rows:
        require(all(np.isfinite(v) and v>0 for v in costs),'nonpositive coarse cost')
        actual=max(max(costs),meta['wall_seconds']/len(rows));cost=max(original_cost,actual)
    else:actual=None;cost=original_cost
    per_step=cost+1.;count=affordable_pairs(len(rows),per_step)
    if not count:reason='NO_AFFORDABLE_QUALIFIED_FINE_PREFIX'
    elif count<len(rows):reason='MEASURED_FINE_BUDGET_PREFIX'
    elif len(rows)<c['registered_pairs']:reason='FAILED_COARSE_QUALIFIED_PREFIX'
    else:reason='FULL_REGISTERED_COARSE_SCOPE'
    return dict(schema=SELECTION,registered_before_fine_worker=True,
        status='REGISTERED_FINE_PREFIX' if count else 'BLOCKED_NO_FINE_RUN',reason=reason,
        window_protocol=relative(c['protocol_path']),window_protocol_sha256=c['protocol_sha256'],
        source_sha256=c['frozen'],parent_checkpoint=relative(c['parent']),parent_sha256=c['parent_hash'],
        coarse_directory=relative(c['window']/'coarse'),coarse_status=meta['status'],coarse_reference_input_sha256=refs,
        registered_pair_count=c['registered_pairs'],qualified_coarse_count=len(rows),selected_pair_count=count,
        selected_pair_indices=list(range(1,count+1)),selected_fine_step_count=2*count,
        registered_end_s=c['protocol']['registered_end_s'],selected_end_s=c['start']+count*c['dt'],
        cost=dict(prior_measured_cost_s=original_cost,coarse_complete_candidate_wall_seconds=costs,
            coarse_stage_wall_seconds=meta['wall_seconds'],coarse_stage_wall_per_complete_candidate_s=
                meta['wall_seconds']/len(rows) if rows else None,measured_cost_s=cost,extra_allowance_s=1.,
            per_candidate_cost_s=per_step,all_qualified_fine_forecast_s=30+1.5*2*len(rows)*per_step,
            selected_fine_forecast_s=30+1.5*2*count*per_step,forecast_limit_s=450.),
        unused_coarse_references=[dict(step=row['step'],global_step=row['global_step'],time_s=row['time_s'],
            checkpoint=relative(c['window']/'coarse'/row['checkpoint']),sha256=row['sha256'],status='UNUSED_NOT_ACCEPTED')
            for row in rows[count:]],unqualified_registered_pair_indices=list(range(len(rows)+1,c['registered_pairs']+1)),
        accepted=False,restart_allowed=False,accepted_time_advanced_s=0.,physical_solver_evaluations=0,nnps_evaluations=0)


def preflight(window,stage):
    c=load_context(window);protocol=c['protocol'];c['coarse']=[];c['coarse_status']=None;c['selection']=None
    if stage=='coarse':
        c.update(pairs=c['registered_pairs'],per_root=protocol['budget']['measured_cost_s'],
            per_step=protocol['budget']['per_candidate_cost_s'],forecast=protocol['budget']['coarse_forecast_s'])
    else:
        path=window/'fine_selection.json';selection=post.read(path)
        require(selection==selection_value(c) and selection['status']=='REGISTERED_FINE_PREFIX','fine selection absent, unaffordable or changed')
        meta,rows,refs=coarse_references(c);c['refs']=refs
        c.update(pairs=selection['selected_pair_count'],coarse=rows[:selection['selected_pair_count']],coarse_status=meta['status'],
            selection=selection,selection_path=path,selection_sha256=core.sha(path),
            per_root=selection['cost']['measured_cost_s'],per_step=selection['cost']['per_candidate_cost_s'],
            forecast=selection['cost']['selected_fine_forecast_s'])
    require(c['forecast']<=450.,'selected stage exceeds forecast cap');return c


def check_context(c):
    bind_map({},c['frozen']);bind_map({},c['parent_evidence']);bind_map({},c['refs'])
    require(model.parameters(.008)==(c['p'],c['provenance']) and core.sha(c['parent'])==c['parent_hash'] and
        core.sha(c['audit'])==c['audit_hash'] and core.sha(c['qualification'])==c['qualification_hash'] and
        core.sha(c['protocol_path'])==c['protocol_sha256'],'source/parent/protocol drift')
    if c['selection'] is not None:require(core.sha(c['selection_path'])==c['selection_sha256'],'fine selection drift')

def publish(item,c,direct_hash,screen_path,coarse_hash):
    pa,state,trial,r=item['pa'],item['state'],item['trial'],item['result'];sd=item['directory']
    require(state['classification']==model.CLASSIFICATION and 'spatial_model_transition' not in state and
        core.clean(state['independent_domain_initialization'])==core.clean(c['state']['independent_domain_initialization']),
        'independent model origin changed')
    screen=post.read(screen_path)
    require(screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'] and
        post.read(sd/'visible_path_audit.json')['status']==old.CERTIFIED,'unqualified publication')
    state.update(accepted=True,restart_allowed=True,event_candidate_accepted=True,
        path_certificate_qualification=dict(qualification=PATH_QUALIFICATION,contract_sha256=core.sha(CONTRACT),
            component_receipt_sha256=PATH_COMPONENTS_SHA,physics_changed=False),
        contact_event_qualification=dict(qualification=qualified.QUALIFICATION,allowed_nonbore_features=sorted(qualified.ALLOWED_FEATURES),
            contract_sha256=core.sha(CONTRACT),component_receipt_sha256=COMPONENTS_SHA,physics_changed=False),
        implementation='INDEPENDENT_EXPANDED_HOT_VISIBLE_SUPPORT_V1',source_sha256=c['state']['source_sha256'],
        implementation_source_sha256=c['frozen'],controller_source_sha256=c['frozen'],
        continuation_parent_sha256=direct_hash,branch_origin_sha256=c['branch_origin'],window_origin_sha256=c['parent_hash'],
        window_temporal_screen_sha256=core.sha(screen_path),window_coarse_reference_sha256=coarse_hash,
        campaign_window_protocol_sha256=c['protocol_sha256'],campaign_fine_selection_sha256=c['selection_sha256'],
        parameters=c['p'],parameter_provenance=c['provenance'])
    core.checkpoint(sd/'UNPUBLISHED_candidate.npz',pa,dict(state,restart_allowed=False,accepted=False));core.verify_snapshot(sd/'UNPUBLISHED_candidate.npz',pa,dict(state,restart_allowed=False,accepted=False))
    saved=core.load(str(sd/'UNPUBLISHED_candidate.npz'))['solver_data']
    for key in ('interval_energy','bore_observation','independent_domain_initialization','expanded_initial_qualification',
                'implementation_source_sha256','continuation_parent_sha256'):
        require(core.clean(saved[key])==core.clean(state[key]),'staged metadata differs '+key)
    check_context(c)
    core.vtk(sd/'fields.vtk',pa,trial.p,trial.bulk,state['time_s'],model.CLASSIFICATION)
    core.write(sd/'tool_pose.json',core.pose(trial.p,state['time_s']));core.write(sd/'transport.json',state['bore_observation'])
    core.checkpoint(sd/'accepted.npz',pa,state);core.verify_snapshot(sd/'accepted.npz',pa,state)


def worker(window,stage):
    folder=window/stage;folder.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    rows=[];history=[];report=dict(status='RUNNING',stage=stage,window_directory=relative(window),accepted_steps=0,candidate_steps=0,
        classification=model.CLASSIFICATION,physical_accuracy='ENGINEERING_APPROXIMATE_NOT_DEPOSITION_OR_EXPERIMENTAL_VALIDATION')
    core.write(folder/'manifest.json',report)
    def stop():
        try:core.write(folder/'HARD_TIMEOUT.json',dict(accepted_history=history,candidate_history=rows,
            published_checkpoints_requiring_readback=[x.relative_to(folder).as_posix() for x in folder.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        c=preflight(window,stage);pa,state=c['pa'],c['state'];per_root=c['per_root'];per_step=c['per_step'];seeds=set(c['seeds'])
        report.update(path_certificate_qualification=PATH_QUALIFICATION,
            time_resolution_qualification=c['time_resolution_qualification'],
            path_component_receipt=relative(PATH_COMPONENTS),path_component_receipt_sha256=PATH_COMPONENTS_SHA,
            contact_event_qualification=qualified.QUALIFICATION,allowed_nonbore_features=sorted(qualified.ALLOWED_FEATURES),
            contact_event_contract=relative(CONTRACT),contact_event_contract_sha256=core.sha(CONTRACT),
            component_receipt=relative(COMPONENTS),component_receipt_sha256=COMPONENTS_SHA,
            source_sha256=c['frozen'],provenance=c['provenance'],parent_checkpoint=relative(c['parent']),parent_sha256=c['parent_hash'],
            branch_origin_sha256=c['branch_origin'],window_origin_sha256=c['parent_hash'],
            parent_readback=relative(c['audit']),parent_readback_sha256=c['audit_hash'],
            independent_qualification=relative(c['qualification']),independent_qualification_sha256=c['qualification_hash'],
            start_time_s=c['start'],start_global_step=c['start_step'],accepted_time_s=c['start'],candidate_time_s=c['start'],
            registered_end_s=c['protocol']['registered_end_s'],reference_end_s=c['start']+c['pairs']*c['dt'],
            registered_pair_count=c['registered_pairs'],reference_pair_count=c['pairs'],
            requested_pair_count=c['protocol']['requested_pair_count'],requested_end_s=c['protocol']['requested_end_s'],
            coarse_dt_s=c['dt'],fine_dt_s=c['half'],inherited_seed_ids=sorted(seeds),
            coarse_original_status=c['coarse_status'],forecast_wall_s=c['forecast'],complete_stage_cost_per_step_s=per_step,
            window_protocol=relative(c['protocol_path']),window_protocol_sha256=c['protocol_sha256'],
            fine_selection=None if c['selection'] is None else relative(c['selection_path']),
            fine_selection_sha256=None if c['selection'] is None else c['selection_sha256'],
            supervision_receipt=relative(window/'execution'/(stage+'_supervision.json')),
            candidate_history_path=str(folder/'candidate_history.json'))
        core.write(folder/'manifest.json',report);core.write(folder/'accepted_history.json',history);core.write(folder/'candidate_history.json',rows)
        core.checkpoint(folder/'CONTINUATION_PARENT_NOT_RESTARTABLE.npz',pa,dict(state,restart_allowed=False))
        direct_hash=c['parent_hash']
        for pair in range(1,c['pairs']+1):
            pending=[]
            for half in range(1,3 if stage=='fine' else 2):
                require(time.perf_counter()+1.5*per_step+20<=begin+525,'remaining complete-stage budget exhausted')
                number_in_stage=2*(pair-1)+half if stage=='fine' else pair;dt=c['half'] if stage=='fine' else c['dt']
                sd=folder/f'step_{number_in_stage:03d}';sd.mkdir()
                item=qualified.candidate(pa,state,c,dt,sd,begin+525,per_root);pending.append(item)
                row=old.candidate_row(item,number_in_stage,direct_hash if len(pending)==1 else core.sha(pending[-2]['checkpoint']),folder)
                rows.append(row);core.write(folder/'candidate_history.json',rows)
                pa,state=item['pa'],item['state'];per_root=max(per_root,item['wall_seconds'])
                seeds.update(e['particle_id'] for e in row['events']);seeds.update(row['blocked_particle_ids']);seeds.update(state['bore_observation']['crossing_ids'])
                report.update(candidate_steps=len(rows),candidate_time_s=state['time_s']);core.write(folder/'manifest.json',report)
                print(json.dumps(dict(step=number_in_stage,time_s=state['time_s'],blocked=row['blocked_edge_count'],wall_s=item['wall_seconds'],accepted=False)),flush=True)
            if stage=='fine':
                cr=c['coarse'][pair-1];coarsefolder=folder.parent/'coarse';saved=core.load(str(coarsefolder/cr['checkpoint']))
                coarse,cs=saved['arrays']['solid'],saved['solver_data'];require(abs(cs['time_s']-state['time_s'])<1e-15,'pair time mismatch')
                seeds.update(e['particle_id'] for e in cr['events']);seeds.update(cr['blocked_particle_ids']);seeds.update(cs['bore_observation']['crossing_ids'])
                graphs=[]
                with np.load(coarsefolder/f'step_{pair:03d}'/'visible_graph.npz') as data:
                    for a,b in (('i','j'),('raw_i','raw_j')):graphs.append((data['particle_id'].copy(),data[a].copy(),data[b].copy()))
                for item in pending:
                    for a,b in (('i','j'),('raw_i','raw_j')):graphs.append((item['pa'].particle_id,item['result'][a],item['result'][b]))
                selected=event.one_hop_ids(seeds,graphs);glob=post.compare(pa,coarse,c['p']);local=event.local_compare(pa,coarse,selected)
                match=state['bore_observation']['crossing_ids']==cs['bore_observation']['crossing_ids'];screen=folder/f'pair_{pair:03d}_temporal_screen.json'
                core.write(screen,dict(pair=pair,time_s=state['time_s'],parent_sha256=c['parent_hash'],coarse_checkpoint_sha256=cr['sha256'],
                    fine_candidate_sha256=[core.sha(x['checkpoint']) for x in pending],cumulative_event_and_blocked_ids=sorted(seeds),
                    global_screen=glob,local_screen=local,crossing_ids_match=match,fine_crossing_ids=state['bore_observation']['crossing_ids'],
                    coarse_crossing_ids=cs['bore_observation']['crossing_ids'],endpoint_graphs_match=pending[-1]['result']['graph_fingerprint']==cr['graph_fingerprint']))
                require(glob['pass'] and local['passed'] and match,'same-parent global/local/strict-cross comparison failed')
                for item in pending:
                    publish(item,c,direct_hash,screen,cr['sha256']);path=item['directory']/'accepted.npz'
                    entry=dict(step=len(history)+1,global_step=item['state']['step'],time_s=item['state']['time_s'],dt_s=c['half'],
                        parent_sha256=direct_hash,checkpoint=path.relative_to(folder).as_posix(),sha256=core.sha(path),
                        temporal_screen_sha256=core.sha(screen),coarse_reference_sha256=cr['sha256'])
                    history.append(entry);direct_hash=entry['sha256'];core.write(folder/'accepted_history.json',history)
                    report.update(accepted_steps=len(history),accepted_time_s=entry['time_s'],accepted_time_advanced_s=entry['time_s']-c['start'])
                    core.write(folder/'manifest.json',report)
            else:direct_hash=core.sha(pending[-1]['checkpoint'])
            check_context(c);per_step=max(per_step,(time.perf_counter()-begin)/len(rows),per_root)
        top_ids=sorted({e['particle_id'] for row in rows for e in row['events'] if e['feature']=='TOP'})
        report.update(top_event_observed_in_candidates=bool(top_ids),top_particle_ids=top_ids)
        report.update(status=PASS_COARSE if stage=='coarse' else PASS_FINE,
            selected_fine_scope_completed=stage=='fine',
            registered_window_completed=stage=='fine' and c['pairs']==c['registered_pairs'] and c['coarse_status']==PASS_COARSE,
            requested_target_completed=stage=='fine' and c['pairs']==c['protocol']['requested_pair_count'] and c['coarse_status']==PASS_COARSE,
            interval_energy=state['interval_energy'],crossing_ids=state['bore_observation']['crossing_ids'],crossing_mass_kg=state['crossing_mass_kg'])
        code=0
    except Exception as error:
        traceback.print_exc();published=[x.relative_to(folder).as_posix() for x in folder.glob('step_*/accepted.npz')]
        status='FAILED_PUBLISHED_CHECKPOINTS_REQUIRE_READBACK' if len(published)!=len(history) else ('FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE')
        report.update(status=status,error=repr(error),published_checkpoints_requiring_readback=published);code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes());core.write(folder/'manifest.json',report)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('coarse','fine'),required=True);parser.add_argument('--output',required=True)
    parser.add_argument('--parent');parser.add_argument('--parent-audit');parser.add_argument('--coarse-dt-us',choices=tuple(DT_CHOICES))
    parser.add_argument('--pairs',type=int);parser.add_argument('--plan-only',action='store_true');parser.add_argument('--worker',action='store_true')
    args=parser.parse_args();window=project(args.output)
    require(window.is_relative_to(ROOT/'outputs') and window!=ROOT/'outputs','new child directory under outputs required')
    if args.worker:
        require(not args.plan_only and not any((args.parent,args.parent_audit,args.coarse_dt_us,args.pairs)),
            'worker reads only immutable protocol inputs')
        return worker(window,args.stage)
    protocol_path=window/'window_protocol.json'
    if args.stage=='coarse':
        if not protocol_path.exists():
            require(all((args.parent,args.parent_audit,args.coarse_dt_us,args.pairs)) and not window.exists(),
                'new coarse registration requires explicit parent/audit/dt/pair bound and new output')
            require(not project(args.parent).is_relative_to(window) and not project(args.parent_audit).is_relative_to(window),
                'new output may not contain/overwrite its own inputs')
            protocol=protocol_value(args.parent,args.parent_audit,args.coarse_dt_us,args.pairs,window)
            window.mkdir(parents=True,exist_ok=False);write_new(protocol_path,protocol)
        else:
            require(not any((args.parent,args.parent_audit,args.coarse_dt_us,args.pairs)),
                'existing protocol must be consumed without replacement input arguments')
            protocol=post.read(protocol_path)
        if protocol['status']!='REGISTERED_SINGLE_WINDOW':
            print(json.dumps(dict(status=protocol['status'],output=relative(window))));return 2
        load_context(window)
        if args.plan_only:
            print(json.dumps(dict(status='PLAN_REGISTERED_NO_PHYSICS',output=relative(window),budget=protocol['budget']),indent=2));return 0
    else:
        require(not any((args.parent,args.parent_audit,args.coarse_dt_us,args.pairs,args.plan_only)),
                'fine consumes the existing immutable window protocol')
        c=load_context(window);path=window/'fine_selection.json'
        require(not path.exists(),'refuse immutable fine selection overwrite')
        selection=selection_value(c);write_new(path,selection)
        if selection['status']!='REGISTERED_FINE_PREFIX':
            print(json.dumps(dict(status=selection['status'],output=relative(window))));return 2
    c=preflight(window,args.stage)
    execution=window/'execution';execution.mkdir(exist_ok=True);log=execution/(args.stage+'.log')
    require(not log.exists() and not (window/args.stage).exists(),'refuse stage evidence overwrite or re-run')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--output',relative(window),'--worker']
    start=time.perf_counter()
    with log.open('x',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    error=None
    try:check_context(c)
    except Exception as problem:error=repr(problem);code=2
    write_new(execution/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,
        external_wall_s=time.perf_counter()-start,limit_s=600,log_sha256=core.sha(log),
        source_inputs_still_match=error is None,source_error=error,window_protocol_sha256=c['protocol_sha256'],
        fine_selection_sha256=None if c['selection'] is None else c['selection_sha256']))
    print(json.dumps(dict(stage=args.stage,exit_code=code,output=relative(window))));return code


if __name__=='__main__':raise SystemExit(main())

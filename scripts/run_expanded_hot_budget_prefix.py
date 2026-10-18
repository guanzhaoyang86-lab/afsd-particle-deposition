"""Budget-selected seven-pair fine prefix using all-audited existing coarse references."""
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

core,post,event=old.core,old.post,old.event
ROOT=core.ROOT
BASE=ROOT/'outputs/expanded_hot_budget_prefix_20261005'
LOGS=ROOT/'logs/expanded_hot_budget_prefix_20261005'
CONTRACT=ROOT/'docs/EXPANDED_HOT_BUDGET_PREFIX_CONTRACT_20261005.md'
ORIGINAL_CONTROLLER=ROOT/'scripts/run_expanded_hot_continuation.py'
ORIGINAL_CONTRACT=ROOT/'docs/EXPANDED_HOT_CONTINUATION_CONTRACT_20261005.md'
COARSE_DIRECTORY=ROOT/'outputs/expanded_hot_continuation_20261005/window_001/coarse'
ORIGINAL_FINE=COARSE_DIRECTORY.parent/'fine'
ORIGINAL_PAIRS=8
FIRST_WINDOW=ROOT/'outputs/expanded_hot_startup_20261005/window_001'
FIRST=FIRST_WINDOW/'fine/step_002/accepted.npz'
FIRST_SHA='82d3ed6b1c0e64786526c6cd9e1a6f1994f8bb97a10017591c63d67ac0120420'
FIRST_AUDIT=FIRST_WINDOW/'readback/summary.json'
FIRST_AUDIT_SHA='52aa6f27a293d00b0d4e934ea2c49d198ebd70a80392a5f1eed7e412bb74e360'
BRANCH_ORIGIN_SHA='31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3'
START,START_STEP=.000015625,2
QUALIFICATION=ROOT/'outputs/expanded_hot_qualification_20261005/v1/summary.json'
COMPONENTS=ROOT/'outputs/visible_top_components_20261004/v1/summary.json'
COMPONENTS_SHA='7fbb5fdca371335ff528f4132fb840948dc3239de8eb50e420b3533014505cab'
PATH_COMPONENTS=ROOT/'outputs/visible_path_tight_20261004/qualification_v1/summary.json'
PATH_COMPONENTS_SHA='3d461593a2fab30900a83bcde5b2835b353bec3d66783393b3500ce493b750bc'
PATH_QUALIFICATION='FIXED_RETAINED_GRAPH_RODRIGUES_CANCELLATION_ENCLOSURE_V1'
DT,HALF,PAIRS=.0000078125,.00000390625,7
require=old.require
PASS_COARSE='QUALIFIED_EXPANDED_CONTINUATION_COARSE_NOT_ACCEPTED'
PASS_FINE='PASS_EXPANDED_BUDGET_SELECTED_FINE_PREFIX_NOT_FULL_WINDOW'


def relative(path):return Path(path).relative_to(ROOT).as_posix()


def bind_map(frozen,mapping):
    post.check_sources(mapping)
    for name,digest in mapping.items():
        require(name not in frozen or frozen[name]==digest,'conflicting source binding '+name)
        frozen[name]=digest


def preflight(number,stage):
    require(number==1 and stage=='fine','only the registered seven-pair fine budget prefix is authorized')
    require(not ORIGINAL_FINE.exists(),'original sixteen-step fine must remain NOT_RUN')
    require(core.sha(FIRST)==FIRST_SHA and core.sha(FIRST_AUDIT)==FIRST_AUDIT_SHA,
        'real accepted parent/readback drift')
    pa,state=core.restore(FIRST);p,provenance=model.parameters(.008)
    require(state['classification']==model.CLASSIFICATION and state['model_scope']==model.SCOPE and
        state['operator_variant']=='visible' and state['independent_initialization'] is True and
        state['restart_allowed'] is True and state['accepted'] is True and state['time_s']==START and
        state['step']==START_STEP and state['parameters']==p and state['parameter_provenance']==provenance and
        state['source_sha256']==model.source_hashes() and 'spatial_model_transition' not in state and
        state['branch_origin_sha256']==BRANCH_ORIGIN_SHA and state['window_origin_sha256']==BRANCH_ORIGIN_SHA,
        'actual accepted independent startup prefix required; no reset/old-state transition')
    model.discrete_process.validate_interval_energy_state(state['interval_energy'])
    model.discrete_process.validate_bore_crossing_state(state['bore_observation'])
    require(state['interval_energy']['steps']==START_STEP and len(state['interval_energy']['records'])==START_STEP,
        'complete genuine initial energy prefix required')
    model.make_trial(pa,state)  # constructor only, no NNPS/evaluation
    audit=post.read(FIRST_AUDIT)
    require(audit['status']=='PASS_EXPANDED_HOT_STARTUP_READBACK_NOT_PHYSICAL_VALIDATION' and
        audit['exit_code']==0 and audit['window_completed'] is False and audit['accepted_time_s']==START and
        audit['accepted_global_steps']==START_STEP and audit['added_accepted_steps']==START_STEP and
        audit['accepted_history'][-1]['sha256']==FIRST_SHA and audit['accepted_history'][-1]['checkpoint']==relative(FIRST) and
        audit['stages']['coarse']['status']=='FAILED_NO_NEW_ACCEPTED_STATE' and
        audit['stages']['fine']['status']=='PASS_QUALIFIED_FINE_PREFIX_OF_FAILED_COARSE' and
        len(audit['pair_screens'])==1 and audit['pair_screens'][0]['passed'],
        'independent readback must preserve the startup failure and its exact accepted prefix')
    post.check_sources(audit['evidence_sha256'])
    ledger=post.read(ROOT/audit['trajectory']['complete_interval_energy_path'])
    require(core.sha(ROOT/audit['trajectory']['complete_interval_energy_path'])==
        audit['trajectory']['complete_interval_energy_sha256'] and core.clean(state['interval_energy'])==ledger,
        'parent saved complete interval energy differs from its readback')
    q=post.read(QUALIFICATION)
    require(q['status']=='PASS_EXPANDED_INITIAL_ELASTIC_RANK_AND_SAME_MAP_NOT_PHYSICAL_VALIDATION' and
        q['qualification_passed'] is True and q['accepted'] is False and q['restart_allowed'] is False and
        q['accepted_time_advanced_s']==0. and q['initialization_sha256']==BRANCH_ORIGIN_SHA and
        q['parameters']==p and q['parameter_provenance']==provenance and
        state['expanded_initial_qualification']['receipt_sha256']==core.sha(QUALIFICATION),
        'preserved initial qualification differs')
    frozen={};bind_map(frozen,state['implementation_source_sha256'])
    require(state['controller_source_sha256']==frozen,'parent controller source map differs')
    extras=[ORIGINAL_CONTROLLER,Path(qualified.__file__),ROOT/'scripts/expanded_hot_continuation_acceptance.py',
        ORIGINAL_CONTRACT,FIRST,FIRST_AUDIT]
    measured=[]
    for prior_stage in ('coarse','fine'):
        parent_meta=post.read(FIRST_WINDOW/prior_stage/'manifest.json')
        parent_rows=post.read(FIRST_WINDOW/prior_stage/'candidate_history.json')
        require(parent_meta['source_sha256']==state['implementation_source_sha256'] and parent_rows and
            parent_meta['candidate_steps']==len(parent_rows),'parent complete-candidate cost history differs')
        for row in parent_rows:
            require(np.isfinite(row['wall_seconds']) and row['wall_seconds']>0.,'invalid measured candidate cost')
            measured.append(row['wall_seconds'])
        extras.extend(FIRST_WINDOW/prior_stage/name for name in ('manifest.json','candidate_history.json'))
    bind_map(frozen,{relative(path):core.sha(path) for path in extras})
    original_frozen=dict(frozen)
    folder=COARSE_DIRECTORY;meta=post.read(folder/'manifest.json');all_coarse=post.read(folder/'candidate_history.json')
    supervision=post.read(ROOT/meta['supervision_receipt'])
    require(meta['status']==PASS_COARSE and meta['source_sha256']==original_frozen and
        meta['parent_sha256']==FIRST_SHA and meta['accepted_steps']==0 and meta['accepted_time_s']==START and
        meta['candidate_steps']==ORIGINAL_PAIRS and len(all_coarse)==ORIGINAL_PAIRS and
        meta['registered_pair_count']==ORIGINAL_PAIRS and meta['reference_pair_count']==ORIGINAL_PAIRS and
        abs(meta['candidate_time_s']-(START+ORIGINAL_PAIRS*DT))<1e-15 and
        supervision['exit_code']==0 and supervision['limit_s']==600 and supervision['external_wall_s']<=600.,
        'all eight original coarse references must be complete, supervised and unaccepted')
    require(not list(folder.glob('step_*/accepted.npz')) and not post.read(folder/'accepted_history.json'),
        'original coarse improperly published')
    refs={};previous=FIRST_SHA;previous_energy=state['interval_energy']['records']
    for index,row in enumerate(all_coarse,1):
        require(row['step']==index and row['parent_sha256']==previous and row['dt_s']==DT and
            row['global_step']==START_STEP+index and abs(row['time_s']-(START+index*DT))<1e-15,
            'original coarse chain/time differs')
        for name,digest in row['artifacts'].items():
            path=(folder/name).resolve()
            require(path.is_relative_to((folder/f'step_{index:03d}').resolve()) and core.sha(path)==digest,
                'original coarse artifact drift/path escape');refs[path]=digest
        path=(folder/row['checkpoint']).resolve()
        require(path==folder/f'step_{index:03d}/NOT_RESTARTABLE_candidate.npz' and
            core.sha(path)==row['sha256'],'original coarse checkpoint/hash differs')
        candidate=core.load(str(path))['solver_data']
        model.discrete_process.validate_interval_energy_state(candidate['interval_energy'])
        require(candidate['restart_allowed'] is False and candidate['accepted'] is False and
            candidate['classification']==model.CLASSIFICATION and candidate['time_s']==row['time_s'] and
            candidate['step']==START_STEP+index and candidate['implementation_source_sha256']==original_frozen and
            candidate['interval_energy']['records'][:-1]==previous_energy and
            candidate['interval_energy']['records'][:START_STEP]==state['interval_energy']['records'],
            'original coarse accepted, changed sources or reset a complete energy prefix')
        previous=row['sha256'];previous_energy=candidate['interval_energy']['records']
    for name in ('manifest.json','candidate_history.json','accepted_history.json'):
        refs[folder/name]=core.sha(folder/name)
    refs[ROOT/meta['supervision_receipt']]=core.sha(ROOT/meta['supervision_receipt'])
    # The old references predate this new controller. Their exact old source map
    # is checked first; only the new fine candidates acquire these two files.
    bind_map(frozen,{relative(path):core.sha(path) for path in (Path(__file__),CONTRACT)})
    count=PAIRS;coarse=all_coarse[:count];coarse_status=meta['status']
    per_root=max(max(measured),max(row['wall_seconds'] for row in all_coarse));per_step=per_root+1.
    original_forecast=30.+1.5*(2*ORIGINAL_PAIRS)*per_step
    forecast=30.+1.5*(2*count)*per_step
    require(original_forecast>450. and forecast<=450.,'registered first-seven budget selection no longer applies')
    selection=dict(reason='ORIGINAL_16_FINE_STEP_FORECAST_EXCEEDS_450_SECONDS',
        original_coarse_directory=relative(folder),original_coarse_status=coarse_status,
        original_coarse_manifest_sha256=core.sha(folder/'manifest.json'),original_reference_pair_count=ORIGINAL_PAIRS,
        selected_pair_indices=list(range(1,PAIRS+1)),selected_pair_count=PAIRS,selected_fine_step_count=2*PAIRS,
        original_16_step_forecast_s=original_forecast,selected_14_step_forecast_s=forecast,forecast_limit_s=450.,
        original_fine_directory=relative(ORIGINAL_FINE),original_fine_status='NOT_RUN_BUDGET_PREFLIGHT_REJECTED',
        original_registered_end_s=START+ORIGINAL_PAIRS*DT,selected_reference_end_s=START+PAIRS*DT,
        original_full_window_completed=False,
        unused_coarse_references=[dict(step=row['step'],global_step=row['global_step'],time_s=row['time_s'],
            checkpoint=relative(folder/row['checkpoint']),sha256=row['sha256'],status='UNUSED_NOT_ACCEPTED')
            for row in all_coarse[PAIRS:]])
    seeds=set(audit['pair_screens'][-1]['cumulative_event_and_blocked_ids'])
    return dict(pa=pa,state=state,p=p,provenance=provenance,frozen=frozen,refs=refs,
        parent=FIRST,parent_hash=FIRST_SHA,audit=FIRST_AUDIT,audit_hash=FIRST_AUDIT_SHA,
        qualification=QUALIFICATION,qualification_hash=core.sha(QUALIFICATION),
        seeds=seeds,prefix=state['interval_energy']['records'],
        coarse=coarse,coarse_status=coarse_status,pairs=count,per_root=per_root,per_step=per_step,
        forecast=forecast,start=START,start_step=START_STEP,branch_origin=BRANCH_ORIGIN_SHA,
        original_frozen=original_frozen,budget_selection=selection,coarse_directory=COARSE_DIRECTORY,
        cost_evidence=dict(measured_complete_candidate_seconds=measured,extra_allowance_s=1.,
            inherited_complete_energy_records=START_STEP))


def check_context(c):
    require(not ORIGINAL_FINE.exists(),'original sixteen-step fine must remain NOT_RUN')
    if 'budget_registration' in c:
        require(core.sha(c['budget_registration'])==c['budget_registration_sha256'],
                'before-worker budget registration changed')
    post.check_sources(c['frozen'])
    require(model.parameters(.008)==(c['p'],c['provenance']),'expanded registered inputs changed')
    for path,digest in c['refs'].items():require(core.sha(path)==digest,'bound coarse reference changed')
    require(core.sha(c['parent'])==c['parent_hash'] and core.sha(c['audit'])==c['audit_hash'] and
        core.sha(c['qualification'])==c['qualification_hash'],'independent parent/qualification drift')


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


def worker(number,stage):
    folder=BASE/f'window_{number:03d}'/stage;folder.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    rows=[];history=[];report=dict(status='RUNNING',stage=stage,window=number,accepted_steps=0,candidate_steps=0,
        classification=model.CLASSIFICATION,physical_accuracy='ENGINEERING_APPROXIMATE_NOT_DEPOSITION_OR_EXPERIMENTAL_VALIDATION')
    core.write(folder/'manifest.json',report)
    def stop():
        try:core.write(folder/'HARD_TIMEOUT.json',dict(accepted_history=history,candidate_history=rows,
            published_checkpoints_requiring_readback=[x.relative_to(folder).as_posix() for x in folder.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        c=preflight(number,stage);pa,state=c['pa'],c['state'];per_root=c['per_root'];per_step=c['per_step'];seeds=set(c['seeds'])
        registration=LOGS/f'window_{number:03d}'/'fine_budget_registration.json'
        registered=post.read(registration)
        require(registered['registered_before_worker'] is True and registered['source_sha256']==c['frozen'] and
            registered['budget_selection']==c['budget_selection'] and registered['contract_sha256']==core.sha(CONTRACT) and
            registered['coarse_reference_input_sha256']=={relative(path):digest for path,digest in c['refs'].items()},
            'before-worker budget/source/reference registration differs')
        c.update(budget_registration=registration,budget_registration_sha256=core.sha(registration))
        report.update(path_certificate_qualification=PATH_QUALIFICATION,
            path_component_receipt=relative(PATH_COMPONENTS),path_component_receipt_sha256=PATH_COMPONENTS_SHA,
            contact_event_qualification=qualified.QUALIFICATION,allowed_nonbore_features=sorted(qualified.ALLOWED_FEATURES),
            contact_event_contract=relative(CONTRACT),contact_event_contract_sha256=core.sha(CONTRACT),
            component_receipt=relative(COMPONENTS),component_receipt_sha256=COMPONENTS_SHA,
            source_sha256=c['frozen'],provenance=c['provenance'],parent_checkpoint=relative(c['parent']),parent_sha256=c['parent_hash'],
            branch_origin_sha256=c['branch_origin'],window_origin_sha256=c['parent_hash'],
            parent_readback=relative(c['audit']),parent_readback_sha256=c['audit_hash'],
            independent_qualification=relative(c['qualification']),independent_qualification_sha256=c['qualification_hash'],
            start_time_s=c['start'],start_global_step=c['start_step'],accepted_time_s=c['start'],candidate_time_s=c['start'],
            registered_end_s=c['start']+PAIRS*DT,reference_end_s=c['start']+c['pairs']*DT,
            registered_pair_count=PAIRS,reference_pair_count=c['pairs'],inherited_seed_ids=sorted(seeds),
            coarse_original_status=c['coarse_status'],forecast_wall_s=c['forecast'],complete_stage_cost_per_step_s=per_step,
            cost_evidence=c['cost_evidence'],budget_selection=c['budget_selection'],
            original_coarse_source_sha256=c['original_frozen'],
            coarse_reference_directory=relative(c['coarse_directory']),
            coarse_reference_input_sha256={relative(path):digest for path,digest in c['refs'].items()},
            original_registered_pair_count=ORIGINAL_PAIRS,original_registered_end_s=c['start']+ORIGINAL_PAIRS*DT,
            original_full_window_completed=False,
            budget_registration=relative(registration),budget_registration_sha256=core.sha(registration),
            supervision_receipt=relative(LOGS/f'window_{number:03d}'/(stage+'_supervision.json')),
            candidate_history_path=str(folder/'candidate_history.json'))
        core.write(folder/'manifest.json',report);core.write(folder/'accepted_history.json',history);core.write(folder/'candidate_history.json',rows)
        core.checkpoint(folder/'CONTINUATION_PARENT_NOT_RESTARTABLE.npz',pa,dict(state,restart_allowed=False))
        direct_hash=c['parent_hash']
        for pair in range(1,c['pairs']+1):
            pending=[]
            for half in range(1,3 if stage=='fine' else 2):
                require(time.perf_counter()+1.5*per_step+20<=begin+525,'remaining complete-stage budget exhausted')
                number_in_stage=2*(pair-1)+half if stage=='fine' else pair;dt=HALF if stage=='fine' else DT
                sd=folder/f'step_{number_in_stage:03d}';sd.mkdir()
                item=qualified.candidate(pa,state,c,dt,sd,begin+525,per_root);pending.append(item)
                row=old.candidate_row(item,number_in_stage,direct_hash if len(pending)==1 else core.sha(pending[-2]['checkpoint']),folder)
                rows.append(row);core.write(folder/'candidate_history.json',rows)
                pa,state=item['pa'],item['state'];per_root=max(per_root,item['wall_seconds'])
                seeds.update(e['particle_id'] for e in row['events']);seeds.update(row['blocked_particle_ids']);seeds.update(state['bore_observation']['crossing_ids'])
                report.update(candidate_steps=len(rows),candidate_time_s=state['time_s']);core.write(folder/'manifest.json',report)
                print(json.dumps(dict(step=number_in_stage,time_s=state['time_s'],blocked=row['blocked_edge_count'],wall_s=item['wall_seconds'],accepted=False)),flush=True)
            if stage=='fine':
                cr=c['coarse'][pair-1];coarsefolder=c['coarse_directory'];saved=core.load(str(coarsefolder/cr['checkpoint']))
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
                    entry=dict(step=len(history)+1,global_step=item['state']['step'],time_s=item['state']['time_s'],dt_s=HALF,
                        parent_sha256=direct_hash,checkpoint=path.relative_to(folder).as_posix(),sha256=core.sha(path),
                        temporal_screen_sha256=core.sha(screen),coarse_reference_sha256=cr['sha256'])
                    history.append(entry);direct_hash=entry['sha256'];core.write(folder/'accepted_history.json',history)
                    report.update(accepted_steps=len(history),accepted_time_s=entry['time_s'],accepted_time_advanced_s=entry['time_s']-c['start'])
                    core.write(folder/'manifest.json',report)
            else:direct_hash=core.sha(pending[-1]['checkpoint'])
            check_context(c);per_step=max(per_step,(time.perf_counter()-begin)/len(rows),per_root)
        top_ids=sorted({e['particle_id'] for row in rows for e in row['events'] if e['feature']=='TOP'})
        report.update(top_event_observed_in_candidates=bool(top_ids),top_particle_ids=top_ids)
        report.update(status=PASS_FINE,budget_selected_prefix_completed=True,original_full_window_completed=False,
            interval_energy=state['interval_energy'],crossing_ids=state['bore_observation']['crossing_ids'],crossing_mass_kg=state['crossing_mass_kg'])
        code=0
    except Exception as error:
        traceback.print_exc();published=[x.relative_to(folder).as_posix() for x in folder.glob('step_*/accepted.npz')]
        status='FAILED_PUBLISHED_CHECKPOINTS_REQUIRE_READBACK' if len(published)!=len(history) else ('FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE')
        report.update(status=status,error=repr(error),published_checkpoints_requiring_readback=published);code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes());core.write(folder/'manifest.json',report)
    return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--window',type=int,required=True);parser.add_argument('--stage',choices=['fine'],required=True)
    parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.window,args.stage))
    require(args.window==1 and args.stage=='fine','only the registered seven-pair fine prefix is permitted')
    logs=LOGS/f'window_{args.window:03d}';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/f'window_{args.window:03d}'/args.stage).exists(),'refuse evidence overwrite')
    registration=logs/'fine_budget_registration.json'
    require(not registration.exists(),'refuse before-worker budget registration overwrite')
    context=preflight(args.window,args.stage)
    core.write(registration,dict(registered_before_worker=True,contract=relative(CONTRACT),contract_sha256=core.sha(CONTRACT),
        source_sha256=context['frozen'],budget_selection=context['budget_selection'],
        coarse_reference_input_sha256={relative(path):digest for path,digest in context['refs'].items()},
        physical_solver_evaluations=0,nnps_evaluations=0,accepted=False,restart_allowed=False,
        accepted_time_advanced_s=0.,internal_limit_s=540.,external_limit_s=600.))
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--window',str(args.window),'--stage',args.stage,'--worker'];start=time.perf_counter()
    with log.open('w',encoding='utf-8') as out:
        try:code=subprocess.run(command,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,external_wall_s=time.perf_counter()-start,limit_s=600))
    print(json.dumps(dict(window=args.window,stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

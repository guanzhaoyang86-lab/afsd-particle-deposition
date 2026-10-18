"""One bounded continuation of the real independent expanded accepted prefix."""
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
BASE=ROOT/'outputs/expanded_hot_continuation_20261005'
LOGS=ROOT/'logs/expanded_hot_continuation_20261005'
CONTRACT=ROOT/'docs/EXPANDED_HOT_CONTINUATION_CONTRACT_20261005.md'
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
DT,HALF,PAIRS=.0000078125,.00000390625,8
require=old.require
PASS_COARSE='QUALIFIED_EXPANDED_CONTINUATION_COARSE_NOT_ACCEPTED'
PASS_FINE='PASS_EXPANDED_CONTINUATION_WINDOW_NOT_DEPOSITION_VALIDATION'


def relative(path):return Path(path).relative_to(ROOT).as_posix()


def bind_map(frozen,mapping):
    post.check_sources(mapping)
    for name,digest in mapping.items():
        require(name not in frozen or frozen[name]==digest,'conflicting source binding '+name)
        frozen[name]=digest


def preflight(number,stage):
    require(number==1,'only the registered 15.625--78.125us continuation window is authorized')
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
    extras=[Path(__file__),Path(qualified.__file__),ROOT/'scripts/expanded_hot_continuation_acceptance.py',
        CONTRACT,FIRST,FIRST_AUDIT]
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
    count=PAIRS;refs={};coarse=[];coarse_status=None
    if stage=='fine':
        folder=BASE/'window_001/coarse';meta=post.read(folder/'manifest.json');coarse=post.read(folder/'candidate_history.json')
        require(meta['status'] in (PASS_COARSE,'FAILED_NO_NEW_ACCEPTED_STATE') and meta['source_sha256']==frozen and
            meta['parent_sha256']==FIRST_SHA and meta['accepted_steps']==0 and 0<len(coarse)<=PAIRS,
            'no qualified same-real-parent coarse prefix')
        require(not list(folder.glob('step_*/accepted.npz')),'coarse improperly published')
        if meta['status']==PASS_COARSE:require(len(coarse)==PAIRS,'complete coarse count differs')
        count=len(coarse);coarse_status=meta['status'];previous=FIRST_SHA
        for index,row in enumerate(coarse,1):
            require(row['step']==index and row['parent_sha256']==previous and row['dt_s']==DT and
                row['global_step']==START_STEP+index,'coarse numerical chain differs')
            for name,digest in row['artifacts'].items():
                path=(folder/name).resolve()
                require(path.is_relative_to((folder/f'step_{index:03d}').resolve()) and
                    core.sha(path)==digest,'coarse artifact drift/path escape');refs[path]=digest
            path=folder/row['checkpoint'];require(core.sha(path)==row['sha256'],'coarse checkpoint hash differs')
            candidate=core.load(str(path))['solver_data']
            require(candidate['restart_allowed'] is False and candidate['accepted'] is False and
                candidate['classification']==model.CLASSIFICATION and
                candidate['interval_energy']['records'][:START_STEP]==state['interval_energy']['records'],
                'coarse incorrectly accepted or reset the inherited energy')
            previous=row['sha256']
        for name in ('manifest.json','candidate_history.json'):refs[folder/name]=core.sha(folder/name)
    per_root=max(measured);per_step=per_root+1.
    if coarse:per_root=max(per_root,max(row['wall_seconds'] for row in coarse));per_step=max(per_step,per_root+1.)
    forecast=30.+1.5*count*(2 if stage=='fine' else 1)*per_step
    require(forecast<=450.,'complete stage forecast exceeds 450s')
    seeds=set(audit['pair_screens'][-1]['cumulative_event_and_blocked_ids'])
    return dict(pa=pa,state=state,p=p,provenance=provenance,frozen=frozen,refs=refs,
        parent=FIRST,parent_hash=FIRST_SHA,audit=FIRST_AUDIT,audit_hash=FIRST_AUDIT_SHA,
        qualification=QUALIFICATION,qualification_hash=core.sha(QUALIFICATION),
        seeds=seeds,prefix=state['interval_energy']['records'],
        coarse=coarse,coarse_status=coarse_status,pairs=count,per_root=per_root,per_step=per_step,
        forecast=forecast,start=START,start_step=START_STEP,branch_origin=BRANCH_ORIGIN_SHA,
        cost_evidence=dict(measured_complete_candidate_seconds=measured,extra_allowance_s=1.,
            inherited_complete_energy_records=START_STEP))


def check_context(c):
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
            cost_evidence=c['cost_evidence'],
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
        report.update(status=PASS_COARSE if stage=='coarse' else (PASS_FINE if c['pairs']==PAIRS else 'PASS_QUALIFIED_FINE_PREFIX_OF_FAILED_COARSE'),
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
    parser=argparse.ArgumentParser();parser.add_argument('--window',type=int,required=True);parser.add_argument('--stage',choices=['coarse','fine'],required=True)
    parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.window,args.stage))
    require(args.window==1,'only the registered 15.625--78.125us continuation is permitted')
    logs=LOGS/f'window_{args.window:03d}';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/f'window_{args.window:03d}'/args.stage).exists(),'refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--window',str(args.window),'--stage',args.stage,'--worker'];start=time.perf_counter()
    with log.open('w',encoding='utf-8') as out:
        try:code=subprocess.run(command,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,external_wall_s=time.perf_counter()-start,limit_s=600))
    print(json.dumps(dict(window=args.window,stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

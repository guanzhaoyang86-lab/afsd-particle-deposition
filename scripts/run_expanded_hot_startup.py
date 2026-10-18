"""Bounded independent 918-particle startup; paired acceptance only."""
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
import expanded_hot_candidate as qualified
import expanded_hot_domain as model

core,post,event=old.core,old.post,old.event
ROOT=core.ROOT
BASE=ROOT/'outputs/expanded_hot_startup_20261005'
LOGS=ROOT/'logs/expanded_hot_startup_20261005'
CONTRACT=ROOT/'docs/EXPANDED_HOT_STARTUP_CONTRACT_20261005.md'
FIRST=ROOT/'outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
FIRST_SHA='31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3'
FIRST_AUDIT=FIRST.parent/'summary.json'
QUALIFICATION=ROOT/'outputs/expanded_hot_qualification_20261005/v1/summary.json'
COMPONENTS=ROOT/'outputs/visible_top_components_20261004/v1/summary.json'
COMPONENTS_SHA='7fbb5fdca371335ff528f4132fb840948dc3239de8eb50e420b3533014505cab'
PATH_COMPONENTS=ROOT/'outputs/visible_path_tight_20261004/qualification_v1/summary.json'
PATH_COMPONENTS_SHA='3d461593a2fab30900a83bcde5b2835b353bec3d66783393b3500ce493b750bc'
PATH_QUALIFICATION='FIXED_RETAINED_GRAPH_RODRIGUES_CANCELLATION_ENCLOSURE_V1'
ROOTS=[ROOT/'outputs/expanded_hot_domain_cost_20261004/v3'/name/'root/manifest.json'
    for name in ('dt_15_625us','dt_7_8125us')]
DT,HALF,PAIRS=.000015625,.0000078125,8
require=old.require
PASS_COARSE='QUALIFIED_EXPANDED_STARTUP_COARSE_NOT_ACCEPTED'
PASS_FINE='PASS_EXPANDED_STARTUP_WINDOW_NOT_DEPOSITION_VALIDATION'


def relative(path):return Path(path).relative_to(ROOT).as_posix()


def bind_map(frozen,mapping):
    post.check_sources(mapping)
    for name,digest in mapping.items():
        require(name not in frozen or frozen[name]==digest,'conflicting source binding '+name)
        frozen[name]=digest


def preflight(number,stage):
    require(number==1,'only the registered independent startup window is authorized by this controller')
    require(core.sha(FIRST)==FIRST_SHA,'independent initial checkpoint drift')
    saved=core.load(str(FIRST));pa,state=saved['arrays']['solid'],saved['solver_data']
    p,provenance=model.parameters(.008)
    require(state['classification']==model.CLASSIFICATION and state['restart_allowed'] is False and
        state['accepted'] is False and state['time_s']==0. and state['step']==0 and
        state['parameters']==p and state['parameter_provenance']==provenance and
        state['source_sha256']==model.source_hashes() and 'spatial_model_transition' not in state and
        'branch_origin_sha256' not in state,'explicit independent zero-time origin required')
    model.initial_geometry_and_field_audit(pa,state);model.make_trial(pa,state)
    q=post.read(QUALIFICATION)
    require(q['status']=='PASS_EXPANDED_INITIAL_ELASTIC_RANK_AND_SAME_MAP_NOT_PHYSICAL_VALIDATION' and
        q['qualification_passed'] is True and q['accepted'] is False and q['restart_allowed'] is False and
        q['accepted_time_advanced_s']==0. and q['initialization_sha256']==FIRST_SHA and
        q['parameters']==p and q['parameter_provenance']==provenance,'new-domain initial qualification absent')
    frozen={}
    bind_map(frozen,state['source_sha256'])
    for key in ('source_and_input_sha256','output_sha256'):bind_map(frozen,q[key])
    require(core.sha(COMPONENTS)==COMPONENTS_SHA and core.sha(PATH_COMPONENTS)==PATH_COMPONENTS_SHA,
        'contact/path component receipt drift')
    components=post.read(COMPONENTS)
    require(components['status']=='PASS_TOP_CONTACT_COMPONENTS_NOT_PROCESS_VALIDATION' and
        components['tests_run']==15 and components['exit_code']==0,'contact component qualification absent')
    for key in ('source_sha256','output_sha256'):bind_map(frozen,components[key])
    path=post.read(PATH_COMPONENTS)
    require(path['qualification_passed'] and path['components']['tests_run']==26 and
        all(path['components'][k]==0 for k in ('failures','errors','skipped')) and
        path['replay_status']==old.CERTIFIED and path['accepted_time_added_s']==0.,'path qualification absent')
    bind_map(frozen,path['source_and_input_sha256'])
    bind_map(frozen,{path['replay_file']:path['replay_sha256']})
    for row in path['workers']:
        require(row['exit_code']==0 and row['external_wall_seconds']<=row['external_limit_seconds'],'component worker failed')
        bind_map(frozen,{row['log']:row['log_sha256']})
    operators=ROOT/'outputs/visible_graph_20261004/operators_v1/summary.json'
    op=post.read(operators)
    require(op['status']=='PASS_REGISTERED_COMPONENT_BENCHMARKS_ONLY' and op['exit_code']==0,'visible operator components absent')
    bind_map(frozen,op['source_sha256'])
    if 'input_sha256' in op:bind_map(frozen,op['input_sha256'])
    roots=[]
    for root_path,dt in zip(ROOTS,(DT,HALF)):
        root=post.read(root_path);receipt=post.read(root_path.parent/'supervision.json')
        require(root['status']=='PASS_EXPANDED_SINGLE_ROOT_DIAGNOSTIC_NOT_PROCESS_ACCEPTANCE' and
            root['exit_code']==0 and root['dt_s']==dt and root['accepted_steps']==0 and
            receipt['exit_code']==0 and receipt['manifest_sha256']==core.sha(root_path) and
            receipt['external_wall_s']<=600 and receipt['before_after_hashes_match'],'new-domain root qualification absent')
        bind_map(frozen,root['source_and_input_sha256']);bind_map(frozen,root['artifact_sha256'])
        bind_map(frozen,{relative(root_path.parent/'run.log'):receipt['log_sha256']})
        for artifact in (root_path,root_path.parent/'supervision.json'):frozen[relative(artifact)]=core.sha(artifact)
        roots.append(root)
    extras=[Path(__file__),Path(qualified.__file__),ROOT/'scripts/expanded_hot_acceptance.py',
        ROOT/'scripts/visible_support_path_tight.py',ROOT/'scripts/run_visible_support_process.py',
        CONTRACT,FIRST,FIRST_AUDIT,QUALIFICATION,COMPONENTS,PATH_COMPONENTS,operators]
    bind_map(frozen,{relative(path):core.sha(path) for path in extras})
    count=PAIRS;refs={};coarse=[];coarse_status=None
    if stage=='fine':
        folder=BASE/'window_001/coarse';meta=post.read(folder/'manifest.json');coarse=post.read(folder/'candidate_history.json')
        require(meta['status'] in (PASS_COARSE,'FAILED_NO_NEW_ACCEPTED_STATE') and meta['source_sha256']==frozen and
            meta['parent_sha256']==FIRST_SHA and meta['accepted_steps']==0 and 0<len(coarse)<=PAIRS,
            'no qualified independent coarse prefix')
        require(not list(folder.glob('step_*/accepted.npz')),'coarse improperly published')
        if meta['status']==PASS_COARSE:require(len(coarse)==PAIRS,'complete coarse count differs')
        count=len(coarse);coarse_status=meta['status'];previous=FIRST_SHA
        for index,row in enumerate(coarse,1):
            require(row['step']==index and row['parent_sha256']==previous and row['dt_s']==DT,'coarse numerical chain differs')
            for name,digest in row['artifacts'].items():
                path=folder/name;require(core.sha(path)==digest,'coarse artifact drift');refs[path]=digest
            candidate=core.load(str(folder/row['checkpoint']))['solver_data']
            require(candidate['restart_allowed'] is False and candidate['accepted'] is False and
                candidate['classification']==model.CLASSIFICATION,'coarse incorrectly accepted')
            previous=row['sha256']
        for name in ('manifest.json','candidate_history.json'):refs[folder/name]=core.sha(folder/name)
    root=roots[stage=='fine'];maps=[r['wall_seconds'] for r in root['calls'] if r['kind']=='complete_map' and r['completed']]
    per_root=float(root['wall_seconds']);per_step=per_root+2*max(maps)+1.
    if coarse:per_root=max(per_root,max(row['wall_seconds'] for row in coarse));per_step=max(per_step,per_root)
    forecast=30.+1.5*count*(2 if stage=='fine' else 1)*per_step
    require(forecast<=450.,'complete stage forecast exceeds 450s')
    return dict(pa=pa,state=state,p=p,provenance=provenance,frozen=frozen,refs=refs,
        parent=FIRST,parent_hash=FIRST_SHA,audit=FIRST_AUDIT,audit_hash=core.sha(FIRST_AUDIT),
        qualification=QUALIFICATION,qualification_hash=core.sha(QUALIFICATION),
        seeds=set(pa.particle_id[pa.actuator.astype(bool)].astype(int).tolist()),prefix=[],
        coarse=coarse,coarse_status=coarse_status,pairs=count,per_root=per_root,per_step=per_step,
        forecast=forecast,start=0.,start_step=0)


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
        continuation_parent_sha256=direct_hash,branch_origin_sha256=c['parent_hash'],window_origin_sha256=c['parent_hash'],
        window_temporal_screen_sha256=core.sha(screen_path),window_coarse_reference_sha256=coarse_hash,
        parameters=c['p'],parameter_provenance=c['provenance'])
    core.checkpoint(sd/'UNPUBLISHED_candidate.npz',pa,dict(state,restart_allowed=False,accepted=False));core.verify_snapshot(sd/'UNPUBLISHED_candidate.npz',pa,state)
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
            branch_origin_sha256=c['parent_hash'],window_origin_sha256=c['parent_hash'],
            parent_initialization_summary=relative(c['audit']),parent_initialization_sha256=c['audit_hash'],
            independent_qualification=relative(c['qualification']),independent_qualification_sha256=c['qualification_hash'],
            start_time_s=c['start'],start_global_step=c['start_step'],accepted_time_s=c['start'],candidate_time_s=c['start'],
            registered_end_s=c['start']+PAIRS*DT,reference_end_s=c['start']+c['pairs']*DT,
            registered_pair_count=PAIRS,reference_pair_count=c['pairs'],inherited_seed_ids=sorted(seeds),
            coarse_original_status=c['coarse_status'],forecast_wall_s=c['forecast'],complete_stage_cost_per_step_s=per_step,
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
    require(args.window==1,'only the registered independent startup window is permitted')
    logs=LOGS/f'window_{args.window:03d}';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/f'window_{args.window:03d}'/args.stage).exists(),'refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--window',str(args.window),'--stage',args.stage,'--worker'];start=time.perf_counter()
    with log.open('w',encoding='utf-8') as out:
        try:code=subprocess.run(command,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,external_wall_s=time.perf_counter()-start,limit_s=600))
    print(json.dumps(dict(window=args.window,stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

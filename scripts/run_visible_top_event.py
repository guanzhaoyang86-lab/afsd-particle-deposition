"""Explicit TOP qualification and bounded continuation of unchanged physics.

Static adaptation of run_visible_deposition; original failures remain immutable.
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
import visible_top_candidate as qualified

core,post,event,model=old.core,old.post,old.event,old.model
ROOT=core.ROOT;BASE=ROOT/'outputs/visible_top_event_20261004'
LOGS=ROOT/'logs/visible_top_event_20261004'
CONTRACT=ROOT/'docs/VISIBLE_TOP_EVENT_PROCESS_CONTRACT_20261004.md'
PREVIOUS=ROOT/'outputs/visible_deposition_20261004/window_001'
FIRST=PREVIOUS/'fine/step_040/accepted.npz'
FIRST_SHA='97eb0b8ed76892a8dd6c85100033c083a650ca1f0c3c8c854d32c877cdadcf1d'
FIRST_AUDIT=PREVIOUS/'readback/summary.json'
FIRST_AUDIT_SHA='9e0866f0433c4321470a0a7dcb25a26c55e9ad6e70aec46ec123398892f2f3c0'
COMPONENTS=ROOT/'outputs/visible_top_components_20261004/v1/summary.json'
COMPONENTS_SHA='7fbb5fdca371335ff528f4132fb840948dc3239de8eb50e420b3533014505cab'
DT,HALF,PAIRS=.000125,.0000625,16
require=old.require
PASS_COARSE='QUALIFIED_CONTINUATION_COARSE_NOT_ACCEPTED'
PASS_FINE='PASS_VISIBLE_CONTINUATION_WINDOW_NOT_DEPOSITION_VALIDATION'


def relative(path):return Path(path).relative_to(ROOT).as_posix()


def parent_for(number):
    require(isinstance(number,int) and 1<=number<=100,'bounded numbered window required')
    if number==1:
        parent,audit=FIRST,FIRST_AUDIT
        require(core.sha(parent)==FIRST_SHA and core.sha(audit)==FIRST_AUDIT_SHA,'first audited parent changed')
        screen=post.read(PREVIOUS/'fine/pair_020_temporal_screen.json')
        costs={s:post.read(PREVIOUS/s/'manifest.json') for s in ('coarse','fine')}
    else:
        previous=BASE/f'window_{number-1:03d}';audit=previous/'readback/summary.json'
        meta=post.read(previous/'fine/manifest.json');history=post.read(previous/'fine/accepted_history.json')
        require(meta['status']==PASS_FINE and len(history)==2*PAIRS,'previous full window incomplete')
        require(post.read(audit)['top_event_process_screen_passed'],
            'previous TOP window has no audited real event')
        parent=previous/'fine'/history[-1]['checkpoint'];require(core.sha(parent)==history[-1]['sha256'],'previous parent hash mismatch')
        screen=post.read(previous/'fine'/f'pair_{PAIRS:03d}_temporal_screen.json')
        costs={s:post.read(previous/s/'manifest.json') for s in ('coarse','fine')}
    reviewed=post.read(audit)
    require(reviewed['status'] in ('PASS_SAVED_EVIDENCE_READBACK_NOT_FRESH_PHYSICS','PASS_VISIBLE_CONTINUATION_READBACK_NOT_PHYSICAL_VALIDATION','PASS_VISIBLE_TOP_EVENT_READBACK_NOT_PHYSICAL_VALIDATION'),
        'independent parent readback did not pass')
    require(reviewed['accepted_history'][-1]['sha256']==core.sha(parent),'parent is not audited accepted endpoint')
    # Verify the prior audit's evidence once at entry; its receipt and exact
    # parent/source are then immutable links, avoiding reset of physics history.
    post.check_sources(reviewed['evidence_sha256'])
    return parent,audit,screen,costs


def preflight(number,stage):
    parent,audit,screen,costs=parent_for(number);pa,state=core.restore(parent);p,provenance=core.parameters()
    require(state['classification']==model.CLASSIFICATION and state['restart_allowed'] is True and
        state['parameters']==p and state['parameter_provenance']==provenance and state['source_sha256']==core.frozen_sources(),
        'continuation requires unchanged accepted visible model and inputs')
    model.VisibleSupportTrial(pa,state,p)
    frozen=dict(state['implementation_source_sha256']);post.check_sources(frozen)
    require(core.sha(COMPONENTS)==COMPONENTS_SHA,'TOP component receipt drift')
    components=post.read(COMPONENTS)
    require(components['status']=='PASS_TOP_CONTACT_COMPONENTS_NOT_PROCESS_VALIDATION' and
        components['tests_run']==15 and components['exit_code']==0,'TOP component qualification absent')
    for key in ('source_sha256','output_sha256'):
        post.check_sources(components[key]);frozen.update(components[key])
    extra=[Path(__file__),Path(qualified.__file__),CONTRACT,COMPONENTS,parent,audit,ROOT/'scripts/run_visible_support_process.py']
    frozen.update({relative(path):core.sha(path) for path in extra})
    seeds=set(screen['cumulative_event_and_blocked_ids'])|set(state['bore_observation']['crossing_ids'])
    count=PAIRS;refs={};coarse=[];coarse_status=None
    if stage=='fine':
        folder=BASE/f'window_{number:03d}'/'coarse';meta=post.read(folder/'manifest.json');coarse=post.read(folder/'candidate_history.json')
        require(meta['status'] in (PASS_COARSE,'FAILED_NO_NEW_ACCEPTED_STATE') and meta['source_sha256']==frozen and
            meta['parent_sha256']==core.sha(parent) and meta['accepted_steps']==0 and 0<len(coarse)<=PAIRS,
            'no fully qualified same-parent coarse reference prefix')
        require(not list(folder.glob('step_*/accepted.npz')),'coarse improperly published')
        if meta['status']==PASS_COARSE:require(len(coarse)==PAIRS,'complete coarse status has wrong count')
        count=len(coarse);coarse_status=meta['status'];previous=core.sha(parent)
        for index,row in enumerate(coarse,1):
            require(row['step']==index and row['parent_sha256']==previous and row['dt_s']==DT,'coarse numerical chain differs')
            for name,digest in row['artifacts'].items():
                path=folder/name;require(core.sha(path)==digest,'coarse artifact changed');refs[path]=digest
            saved=core.load(str(folder/row['checkpoint']))['solver_data']
            require(saved['restart_allowed'] is False and saved['classification']==model.CLASSIFICATION,'coarse source is restartable')
            require(saved['interval_energy']['records'][:state['step']]==state['interval_energy']['records'],'coarse lost parent prefix')
            previous=row['sha256']
        for name in ('manifest.json','candidate_history.json'):refs[folder/name]=core.sha(folder/name)
    per_root=max(row['wall_seconds'] for row in post.read(Path(costs[stage].get('candidate_history_path',
        (PREVIOUS/stage/'candidate_history.json') if number==1 else BASE/f'window_{number-1:03d}'/stage/'candidate_history.json'))))
    measured=float(costs[stage]['wall_seconds'])/int(costs[stage]['candidate_steps'])
    if coarse:per_root=max(per_root,max(row['wall_seconds'] for row in coarse))
    per_step=max(measured,per_root);steps=count*(2 if stage=='fine' else 1)
    forecast=30.+1.5*steps*per_step
    require(forecast<=450.,'complete stage cost forecast exceeds 450s')
    return dict(pa=pa,state=state,p=p,provenance=provenance,frozen=frozen,refs=refs,
        parent=parent,parent_hash=core.sha(parent),audit=audit,audit_hash=core.sha(audit),
        seeds=seeds,prefix=list(state['interval_energy']['records']),coarse=coarse,coarse_status=coarse_status,
        pairs=count,per_root=per_root,per_step=per_step,forecast=forecast,start=state['time_s'],start_step=state['step'])


def check_context(c):
    old.check_context(c)
    require(core.sha(c['parent'])==c['parent_hash'] and core.sha(c['audit'])==c['audit_hash'],'qualified parent/audit drift')


def publish(item,c,direct_hash,screen_path,coarse_hash):
    pa,state,trial,r=item['pa'],item['state'],item['trial'],item['result'];sd=item['directory']
    require(state['classification']==model.CLASSIFICATION and not state['spatial_model_transition']['diagnostic_only'] and
        core.clean(state['spatial_model_transition'])==core.clean(c['state']['spatial_model_transition']),
        'continuation altered model transition or diagnostic qualification')
    screen=post.read(screen_path)
    require(screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'] and
        post.read(sd/'visible_path_audit.json')['status']==old.CERTIFIED,'unqualified candidate publication')
    state.update(contact_event_qualification=dict(qualification=qualified.QUALIFICATION,
        allowed_nonbore_features=sorted(qualified.ALLOWED_FEATURES),
        contract_sha256=core.sha(CONTRACT),component_receipt_sha256=COMPONENTS_SHA,physics_changed=False),
        implementation='UNIFIED_VISIBLE_SUPPORT_REBUILT_ALL_OPERATORS_V1',restart_allowed=True,
        source_sha256=core.frozen_sources(),implementation_source_sha256=c['frozen'],controller_source_sha256=c['frozen'],
        continuation_parent_sha256=direct_hash,branch_origin_sha256=c['state']['branch_origin_sha256'],window_origin_sha256=c['parent_hash'],
        window_temporal_screen_sha256=core.sha(screen_path),window_coarse_reference_sha256=coarse_hash,
        parameters=c['p'],parameter_provenance=c['provenance'])
    core.checkpoint(sd/'UNPUBLISHED_candidate.npz',pa,dict(state,restart_allowed=False));core.verify_snapshot(sd/'UNPUBLISHED_candidate.npz',pa,state)
    saved=core.load(str(sd/'UNPUBLISHED_candidate.npz'))['solver_data']
    for key in ('interval_energy','bore_observation','spatial_model_transition','implementation_source_sha256','continuation_parent_sha256'):
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
        report.update(contact_event_qualification=qualified.QUALIFICATION,allowed_nonbore_features=sorted(qualified.ALLOWED_FEATURES),
            contact_event_contract=relative(CONTRACT),contact_event_contract_sha256=core.sha(CONTRACT),
            component_receipt=relative(COMPONENTS),component_receipt_sha256=COMPONENTS_SHA,
            source_sha256=c['frozen'],provenance=c['provenance'],parent_checkpoint=relative(c['parent']),parent_sha256=c['parent_hash'],
            branch_origin_sha256=state['branch_origin_sha256'],window_origin_sha256=c['parent_hash'],
            parent_readback_summary=relative(c['audit']),parent_readback_sha256=c['audit_hash'],
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
        require(number!=1 or top_ids,'registered TOP event not observed; no first-window qualification')
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
    require(1<=args.window<=100,'bounded registered window index required')
    logs=LOGS/f'window_{args.window:03d}';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/f'window_{args.window:03d}'/args.stage).exists(),'refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--window',str(args.window),'--stage',args.stage,'--worker'];start=time.perf_counter()
    with log.open('w',encoding='utf-8') as out:
        try:code=subprocess.run(command,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,external_wall_s=time.perf_counter()-start,limit_s=600))
    print(json.dumps(dict(window=args.window,stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

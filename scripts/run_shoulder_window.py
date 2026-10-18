"""22->30 ms: nonrestartable coarse reference, pair-screened fine commits."""
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
import run_first_shoulder_enclosed as enclosed
from contact_shell_enclosure import audit_contact_shell_enclosure

event,batch,core,post=enclosed.original,enclosed.batch,enclosed.core,enclosed.post
BASE=batch.BASE/'shoulder_window_22_30ms'
PARENT_DIR=enclosed.BASE/'halves'
PARENT_SHA='5368f67e717b4ce42d3927b8eca4d03a502a7ab59c5899d3d9b0b776c71d2d95'
CONTRACT=core.ROOT/'docs/SHOULDER_WINDOW_CONTRACT_20261004.md'
CHAIN_GATE=core.ROOT/'outputs/coupled_cost_20261004/chain_shell_21_875ms/manifest.json'
CHAIN_RECEIPT=core.ROOT/'outputs/coupled_cost_20261004/chain_shell_21_875ms_supervision.json'
require,checked_file=event.require,event.checked_file
COARSE_DT,FINE_DT,PAIRS=.00025,.000125,32


def preflight(stage):
    prior=enclosed.preflight('halves')
    p,provenance,baseline=prior['p'],prior['provenance'],core.frozen_sources()
    meta_file=PARENT_DIR/'manifest.json';history_file=PARENT_DIR/'accepted_history.json'
    meta,rows=post.read(meta_file),post.read(history_file)
    require(meta['status']=='PASS_FIRST_SHOULDER_LOCAL_TIME_SCREEN_NOT_TRANSPORT' and
            meta['global_time_screen']['pass'] and meta['local_time_screen']['passed'] and
            meta['source_sha256']==prior['frozen'] and meta['provenance']==provenance and
            meta['accepted_steps']==len(rows)==2 and abs(meta['accepted_time_s']-.022)<1e-15,
            'qualified 22ms parent time screen absent')
    post.check_sources(meta['source_sha256'])
    references={meta_file:core.sha(meta_file),history_file:core.sha(history_file)}
    prefix=list(prior['state']['interval_energy']['records']);previous=prior['parent_hash'];parent=None
    for number,row in enumerate(rows,1):
        parent=PARENT_DIR/row['checkpoint'];require(core.sha(parent)==row['sha256'] and row['parent_sha256']==previous,'22ms parent chain differs')
        pa,state=core.restore(parent)
        require(state['step']==151+number and state['time_s']==row['time_s'] and
                state['source_sha256']==baseline and state['implementation_source_sha256']==prior['frozen'] and
                state['parameters']==p and state['parameter_provenance']==provenance and
                state['continuation_parent_sha256']==previous and state['branch_origin_sha256']==prior['parent_hash'] and
                state['interval_energy']['records'][:-1]==prefix,'22ms parent checkpoint history/model differs')
        batch.make_trial(pa,state,p)
        shell_file=parent.parent/'continuous_shell_audit.json';shell=post.read(shell_file)
        require(shell['status']==event.CERTIFIED and shell['additional_finite_patch_queries']==0,
                'one of the last two actual segments lacks enclosure certification')
        references.update({parent:row['sha256'],shell_file:core.sha(shell_file)})
        prefix=list(state['interval_energy']['records']);previous=row['sha256']
    require(previous==PARENT_SHA and state['step']==153 and abs(state['time_s']-.022)<1e-15,'registered final parent differs')
    chain,receipt=post.read(CHAIN_GATE),post.read(CHAIN_RECEIPT)
    require(receipt['exit_code']==0 and receipt['physical_evaluations']==0 and
            chain['status']==event.CERTIFIED and chain['total_segment_count']==chain['audited_segment_count']==chain['certified_segment_count']==151 and
            chain['unexamined_segment_count']==0 and chain['physical_evaluations']==0 and
            abs(chain['last_audited_time_s']-.021875)<1e-15 and
            chain['parameters']==p and chain['parameter_provenance']==provenance,
            'actual 0->21.875ms chain certificate absent')
    post.check_sources(chain['source_sha256'])
    chain_inputs={checked_file(name,digest).resolve():digest for name,digest in chain['input_sha256'].items()}
    require(chain_inputs.get(prior['parent'].resolve())==prior['parent_hash'],'151-segment certificate ends at a different parent')
    references.update(chain_inputs);references.update({CHAIN_GATE:core.sha(CHAIN_GATE),CHAIN_RECEIPT:core.sha(CHAIN_RECEIPT)})
    local_file=PARENT_DIR/'temporal_comparison.json';local=post.read(local_file)
    require(local['global_screen']['pass'] and local['local_screen']['passed'],'saved parent local screen differs')
    seeds=set(int(v) for v in local['event_particle_ids']);require(seeds,'initial shoulder IDs absent')
    references[local_file]=core.sha(local_file)
    extra=[Path(__file__).resolve(),CONTRACT,CHAIN_GATE,CHAIN_RECEIPT,meta_file,history_file,local_file]
    frozen=dict(prior['frozen']);frozen.update(chain['source_sha256'])
    frozen.update({path.relative_to(core.ROOT).as_posix():core.sha(path) for path in extra})
    macro=post.read(enclosed.BASE/'macro/manifest.json')
    per_root=max(macro['candidate_wall_seconds'],max(row['wall_seconds'] for row in rows))
    coarse_rows=None
    if stage=='fine':
        coarse_meta_file=BASE/'coarse/manifest.json';coarse_history_file=BASE/'coarse/candidate_history.json'
        coarse_meta,coarse_rows=post.read(coarse_meta_file),post.read(coarse_history_file)
        require(coarse_meta['status']=='NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED' and
                coarse_meta['source_sha256']==frozen and coarse_meta['provenance']==provenance and
                coarse_meta['parent_sha256']==PARENT_SHA and coarse_meta['candidate_steps']==len(coarse_rows)==PAIRS and
                coarse_meta['accepted_steps']==0 and abs(coarse_meta['accepted_time_s']-.022)<1e-15 and
                abs(coarse_meta['candidate_time_s']-.030)<1e-15,'complete same-parent nonaccepted coarse reference absent')
        require(not list((BASE/'coarse').glob('step_*/accepted.npz')),'coarse reference was wrongly published')
        references.update({coarse_meta_file:core.sha(coarse_meta_file),coarse_history_file:core.sha(coarse_history_file)})
        previous=PARENT_SHA
        for number,row in enumerate(coarse_rows,1):
            require(row['step']==number and row['dt_s']==COARSE_DT and row['parent_sha256']==previous and
                    abs(row['time_s']-(.022+number*COARSE_DT))<1e-15,'coarse reference numerical chain differs')
            for path_key,hash_key in (('checkpoint','sha256'),('evidence','evidence_sha256'),('events','events_sha256')):
                path=BASE/'coarse'/row[path_key];require(core.sha(path)==row[hash_key],'coarse reference artifact drift')
                references[path]=row[hash_key]
            audits={}
            for name,artifact in row['audit_artifacts'].items():
                path=BASE/'coarse'/artifact['path'];require(core.sha(path)==artifact['sha256'],'coarse audit artifact drift')
                references[path]=artifact['sha256'];audits[name]=post.read(path)
            require(set(audits)=={'continuous_shell_audit.json','scalar_batch_equivalence.json','geometry_audit.json'} and
                    audits['continuous_shell_audit.json']['status']==event.CERTIFIED and
                    audits['scalar_batch_equivalence.json']['status']=='PASS_SAME_CANDIDATE_COMPLETE_MAP' and
                    not audits['geometry_audit.json']['transient_contacts_missed_by_endpoints'] and
                    set(audits['geometry_audit.json']['new_non_bore_features'])<= {'SHOULDER'},
                    'coarse candidate lacks matching shell/equivalence/feature evidence')
            previous=row['sha256']
        per_root=max(per_root,max(row['wall_seconds'] for row in coarse_rows))
    count=PAIRS*(2 if stage=='fine' else 1);forecast=30+1.5*count*per_root
    require(forecast<=450,'registered shoulder-window forecast exceeds 450s')
    return dict(p=p,provenance=provenance,pa=pa,state=state,parent=parent,parent_hash=PARENT_SHA,
        origin=parent,origin_hash=PARENT_SHA,frozen=frozen,references=references,per_root=per_root,
        forecast=forecast,coarse_rows=coarse_rows,seeds=seeds,prefix=prefix)


def candidate_step(pa,state,p,dt,directory,frozen,deadline,per_root):
    started=time.perf_counter();trial=batch.make_trial(pa,state,p)
    newpa,newstate,result,metrics=event.solve_candidate(trial,dt,directory,deadline,per_root)
    candidate,evidence=event.save_candidate(directory,trial,newpa,newstate,result,frozen)
    geometry,shell,events=event.audit_candidate(trial,result,directory)
    features=set(e['feature'] for e in events)
    require(features<= {'SHOULDER'} and set(geometry['new_non_bore_features'])<= {'SHOULDER'},
            'new non-shoulder finite feature outside qualified scope')
    event.scalar_equivalence(trial,result,directory)
    return dict(pa=newpa,state=newstate,result=result,trial=trial,metrics=metrics,directory=directory,
                candidate=candidate,evidence=evidence,events=events,wall_seconds=time.perf_counter()-started)


def crossing_summary(state):
    observer=state['bore_observation']
    return dict(crossing_ids=list(observer['crossing_ids']),crossing_mass_kg=observer['crossing_mass_kg'],
                events=list(observer['events']),last_time_s=observer['last_time_s'])


def worker(stage):
    # Explicit helper injection only; inherited physics and acceptance remain
    # unchanged, and every resulting certificate retains its original report.
    event.audit_contact_shell=audit_contact_shell_enclosure
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    history=[];candidates=[];screens=[];first_passage=False
    report=dict(stage=stage,status='RUNNING',accepted_steps=0,accepted_time_s=None,candidate_steps=0,
        physical_accuracy='ENGINEERING_PAIRWISE_TIME_SCREEN_NOT_EXPERIMENTAL_VALIDATION',
        controller='FROZEN_SHOULDER_HELPERS_WITH_EXPLICIT_ENCLOSURE_AUDIT',
        coarse_reference_restart_allowed=False)
    def stop():
        try:
            core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,candidate_history=candidates,
                published_checkpoints=[str(x.relative_to(directory)) for x in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        c=preflight(stage);pa,state,p=c['pa'],c['state'],c['p'];per_root=c['per_root'];seeds=set(c['seeds'])
        report.update(classification=core.CLASSIFICATION,source_sha256=c['frozen'],provenance=c['provenance'],
            parent_checkpoint=c['parent'].relative_to(core.ROOT).as_posix(),parent_sha256=PARENT_SHA,
            starting_time_s=state['time_s'],accepted_time_s=state['time_s'],candidate_time_s=state['time_s'],
            requested_time_s=.030,forecast_wall_s=c['forecast'],dt_s=FINE_DT if stage=='fine' else COARSE_DT,
            parent_certificate_scope='151 actual chain segments plus the two certified first-shoulder half segments',
            thread_environment={key:os.environ.get(key) for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',history)
        core.write(directory/'candidate_history.json',candidates)
        initial,initial_state=core.restore(batch.OLD/'quarters/independent_initial.npz')
        core.checkpoint(directory/'independent_initial.npz',initial,initial_state)
        core.checkpoint(directory/'continuation_initial.npz',pa,state)
        direct_file,direct_hash=c['parent'],PARENT_SHA
        for pair in range(1,PAIRS+1):
            pending=[]
            for half in range(1,3 if stage=='fine' else 2):
                require(time.perf_counter()-begin+1.5*per_root+15<=525,'remaining shoulder-window budget exhausted')
                number=2*(pair-1)+half if stage=='fine' else pair
                dt=FINE_DT if stage=='fine' else COARSE_DT
                stepdir=directory/f'step_{number:03d}';stepdir.mkdir()
                item=candidate_step(pa,state,p,dt,stepdir,c['frozen'],begin+525,per_root)
                require(item['state']['interval_energy']['records'][:153]==c['prefix'],'candidate reset 153-record parent prefix')
                pending.append(item);pa,state=item['pa'],item['state'];seeds.update(e['particle_id'] for e in item['events'])
                event_file=stepdir/'sampled_nonbore_events.json'
                candidate_parent=(core.sha(pending[-2]['candidate']) if stage=='fine' and len(pending)>1 else direct_hash)
                row=dict(step=number,global_step=state['step'],time_s=state['time_s'],dt_s=dt,
                    parent_sha256=candidate_parent,checkpoint=f'step_{number:03d}/NOT_RESTARTABLE_candidate.npz',
                    sha256=core.sha(item['candidate']),evidence=f'step_{number:03d}/candidate_evidence.npz',
                    evidence_sha256=core.sha(item['evidence']),events=f'step_{number:03d}/sampled_nonbore_events.json',
                    events_sha256=core.sha(event_file),cumulative_shoulder_ids=sorted(seeds),
                    bore_observation=crossing_summary(state),
                    wall_seconds=item['wall_seconds'],restart_allowed=False,**item['metrics'])
                row['audit_artifacts']={name:dict(path=f'step_{number:03d}/'+name,sha256=core.sha(stepdir/name))
                    for name in ('continuous_shell_audit.json','scalar_batch_equivalence.json','geometry_audit.json')}
                candidates.append(row);core.write(directory/'candidate_history.json',candidates)
                report.update(candidate_steps=len(candidates),candidate_time_s=state['time_s'])
                core.write(directory/'manifest.json',report);per_root=max(per_root,item['wall_seconds'])
            if stage=='coarse':
                event.check_references(c);require(core.sha(direct_file)==direct_hash,'coarse numerical parent drift')
                direct_file=pending[-1]['candidate'];direct_hash=core.sha(direct_file)
                print(json.dumps(dict(candidate_step=pair,time_s=state['time_s'],accepted=False)),flush=True)
                continue
            coarse_row=c['coarse_rows'][pair-1];coarse_path=BASE/'coarse'/coarse_row['checkpoint']
            require(core.sha(coarse_path)==coarse_row['sha256'],'coarse same-time reference drift')
            raw=core.load(str(coarse_path));coarse,coarse_state=raw['arrays']['solid'],raw['solver_data']
            require(coarse_state['restart_allowed'] is False and coarse_state['implementation_source_sha256']==c['frozen'] and
                    coarse_state['parameters']==p and coarse_state['parameter_provenance']==c['provenance'] and
                    coarse_state['interval_energy']['records'][:153]==c['prefix'] and
                    abs(coarse_state['time_s']-state['time_s'])<1e-15,'coarse reference endpoint/model differs')
            seeds.update(coarse_row['cumulative_shoulder_ids'])
            with np.load(BASE/'coarse'/coarse_row['evidence']) as data:
                graphs=[(data['particle_id'].copy(),data['i'].copy(),data['j'].copy())]
            graphs += [(item['pa'].particle_id,item['result']['i'],item['result']['j']) for item in pending]
            selected=event.one_hop_ids(seeds,graphs)
            global_screen=post.compare(pa,coarse,p);local_screen=event.local_compare(pa,coarse,selected)
            screen=dict(pair=pair,time_s=state['time_s'],parent_sha256=PARENT_SHA,
                coarse_checkpoint_sha256=coarse_row['sha256'],cumulative_shoulder_ids=sorted(seeds),
                global_screen=global_screen,local_screen=local_screen,
                coarse_bore_observation=crossing_summary(coarse_state),fine_bore_observation=crossing_summary(state))
            screen['crossing_ids_match']=screen['coarse_bore_observation']['crossing_ids']==screen['fine_bore_observation']['crossing_ids']
            screen_file=directory/f'pair_{pair:03d}_temporal_screen.json';core.write(screen_file,screen)
            report['last_screen']=screen
            require(screen['crossing_ids_match'],'strict bore crossing IDs differ: crossing time unresolved; pair uncommitted')
            require(global_screen['pass'] and local_screen['passed'],'paired shoulder-window global/local time screen failed')
            event.check_references(c)
            for half,item in enumerate(pending,1):
                require(core.sha(direct_file)==direct_hash,'accepted fine parent drift')
                item['state'].update(window_temporal_screen_sha256=core.sha(screen_file),
                    window_coarse_reference_sha256=coarse_row['sha256'],window_pair_index=pair)
                started=time.perf_counter()
                batch.publish(item['directory'],item['pa'],item['state'],item['trial'],item['result'],
                    c['frozen'],c['provenance'],direct_hash,PARENT_SHA)
                path=item['directory']/'accepted.npz';number=2*(pair-1)+half
                row=dict(step=number,global_step=item['state']['step'],time_s=item['state']['time_s'],dt_s=FINE_DT,
                    parent_sha256=direct_hash,checkpoint=f'step_{number:03d}/accepted.npz',sha256=core.sha(path),
                    temporal_screen_sha256=core.sha(screen_file),coarse_reference_sha256=coarse_row['sha256'],
                    wall_seconds=item['wall_seconds']+time.perf_counter()-started,**item['metrics'])
                history.append(row);core.write(directory/'accepted_history.json',history)
                direct_file=path;direct_hash=row['sha256'];per_root=max(per_root,row['wall_seconds'])
                report.update(accepted_steps=len(history),accepted_time_s=item['state']['time_s'],
                    interval_energy=item['state']['interval_energy'],crossing_mass_kg=item['state']['crossing_mass_kg'])
                core.write(directory/'manifest.json',report)
            screens.append(dict(pair=pair,screen=screen_file.name,sha256=core.sha(screen_file)))
            core.write(directory/'screen_history.json',screens);print(json.dumps(dict(pair=pair,time_s=state['time_s'],accepted_steps=len(history))),flush=True)
            if set(state['bore_observation']['crossing_ids'])-set(c['state']['bore_observation']['crossing_ids']):
                first_passage=True
                report.update(status='PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION',
                    first_passage_pair=pair,bore_observation=crossing_summary(state))
                break
        if not first_passage:
            require(abs(state['time_s']-.030)<1e-15,'window terminal time differs')
        if stage=='coarse':
            require(not history and not list(directory.glob('step_*/accepted.npz')),'nonaccepted coarse branch was published')
            report['status']='NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED'
        elif not first_passage:
            require(len(history)==64 and len(screens)==32,'incomplete pairwise acceptance history')
            report.update(status='PASS_WINDOW_NOT_TRANSPORT',max_plastic_strain=float(pa.ep.max()))
        event.check_references(c);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['coarse','fine'],required=True)
    parser.add_argument('--worker',action='store_true');parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    if args.preflight_only:
        c=preflight(args.stage);print(json.dumps(dict(stage=args.stage,parent_sha256=PARENT_SHA,parent_time_s=c['state']['time_s'],
            parent_global_step=c['state']['step'],forecast_wall_s=c['forecast'],physical_evaluations=0)));raise SystemExit(0)
    if args.worker:raise SystemExit(worker(args.stage))
    logs=core.ROOT/'logs/coupled_cost_20261004/batched/shoulder_window_22_30ms';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/args.stage).exists(),'refuse window evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-started,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

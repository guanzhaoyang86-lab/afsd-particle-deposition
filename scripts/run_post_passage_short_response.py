"""One four-millisecond response window after a qualified first centre passage.

All existing physics and qualification helpers remain frozen. Crossed-particle
observations add a stricter local screen, never a force or thermal source.
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
import run_refined_shoulder_continuation as series

finite=series.finite
event,batch,core,post=series.event,series.batch,series.core,series.post
require=event.require
candidate_step,crossing_summary=finite.candidate_step,finite.crossing_summary
audit_contact_shell_finite_enclosure=finite.audit_contact_shell_finite_enclosure
BASE=batch.BASE/'post_passage_short_response'
CONTRACT=core.ROOT/'docs/POST_PASSAGE_SHORT_RESPONSE_CONTRACT_20261004.md'
COARSE_DT,FINE_DT,PAIRS=.000125,.0000625,32
DURATION=.004
FIRST_PASSAGE='PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION'


def passage_parent(passage_window):
    require(passage_window in series.WINDOWS,'first passage must belong to a registered continuation window')
    directory=series.BASE/f'window_{passage_window}'/'fine'
    meta_file=directory/'manifest.json';meta=post.read(meta_file)
    require(meta['status']==FIRST_PASSAGE,'a completed first-passage special PASS is required; RUNNING/FAIL/coarse forbidden')
    before=series.preflight(passage_window,'fine');references=dict(before['references'])
    rows_file=directory/'accepted_history.json';candidates_file=directory/'candidate_history.json';screens_file=directory/'screen_history.json'
    rows,candidates,screens=post.read(rows_file),post.read(candidates_file),post.read(screens_file)
    pairs=meta['first_passage_pair'];start=before['state']['time_s'];end=start+pairs*COARSE_DT
    require(isinstance(pairs,int) and 1<=pairs<=series.PAIRS and
            meta['accepted_steps']==len(rows)==meta['candidate_steps']==len(candidates)==2*pairs and len(screens)==pairs and
            meta['source_sha256']==before['frozen'] and meta['provenance']==before['provenance'] and
            meta['parent_sha256']==before['parent_hash'] and abs(meta['starting_time_s']-start)<1e-15 and
            abs(meta['accepted_time_s']-end)<1e-15 and meta['crossing_mass_kg']>0 and
            meta['bore_observation']['crossing_ids'] and len(list(directory.glob('step_*/accepted.npz')))==2*pairs,
            'first-passage manifest, chain length or qualification differs')
    references.update({path:core.sha(path) for path in (meta_file,rows_file,candidates_file,screens_file)})
    baseline=core.frozen_sources();prefix=list(before['state']['interval_energy']['records'])
    previous=before['parent_hash'];seeds=set(before['seeds']);coarse=before['coarse_rows']
    for index,(row,candidate) in enumerate(zip(rows,candidates),1):
        parent=directory/row['checkpoint']
        require(row['step']==candidate['step']==index and row['dt_s']==candidate['dt_s']==FINE_DT and
                core.sha(parent)==row['sha256'] and row['parent_sha256']==previous and
                candidate['parent_sha256']==(previous if index%2 else candidates[index-2]['sha256']),
                'first-passage accepted/candidate parent chain differs')
        pa,state=core.restore(parent)
        require(state['step']==row['global_step']==candidate['global_step']==before['state']['step']+index and
                state['time_s']==row['time_s']==candidate['time_s'] and abs(state['time_s']-(start+index*FINE_DT))<1e-15 and
                state['source_sha256']==baseline and state['implementation_source_sha256']==before['frozen'] and
                state['parameters']==before['p'] and state['parameter_provenance']==before['provenance'] and
                state['continuation_parent_sha256']==previous and state['branch_origin_sha256']==before['parent_hash'] and
                state['interval_energy']['records'][:-1]==prefix,'first-passage source/parameters/time/ledger differs')
        batch.make_trial(pa,state,before['p'])
        finite.check_artifacts(directory,candidate,references);core.verify_snapshot(directory/candidate['checkpoint'],pa,state)
        saved=core.load(str(directory/candidate['checkpoint']))['solver_data']
        require(saved['restart_allowed'] is False and saved['implementation_source_sha256']==before['frozen'] and
                saved['source_sha256']==baseline and saved['parameters']==before['p'] and
                saved['parameter_provenance']==before['provenance'] and saved['interval_energy']==state['interval_energy'] and
                saved['bore_observation']==state['bore_observation'],'first-passage uncommitted candidate metadata differs')
        pair=(index-1)//2+1;entry=screens[pair-1];screen_file=directory/entry['screen'];screen=post.read(screen_file)
        require(entry['pair']==screen['pair']==pair and core.sha(screen_file)==entry['sha256']==row['temporal_screen_sha256']==state['window_temporal_screen_sha256'] and
                abs(screen['time_s']-(start+pair*COARSE_DT))<1e-15 and screen['parent_sha256']==before['parent_hash'] and
                screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'] and
                screen['coarse_bore_observation']['crossing_ids']==screen['fine_bore_observation']['crossing_ids'] and
                screen['coarse_checkpoint_sha256']==coarse[pair-1]['sha256']==row['coarse_reference_sha256']==state['window_coarse_reference_sha256'],
                'first-passage same-time pair screen differs')
        ids=screen['fine_bore_observation']['crossing_ids']
        require(bool(ids)==(pair==pairs),'first centre crossing was not stopped at its first screened pair')
        if pair<pairs:require(not state['bore_observation']['crossing_ids'],'crossing preceded registered first pair')
        seeds.update(screen['cumulative_shoulder_ids']);references[screen_file]=entry['sha256'];references[parent]=row['sha256']
        prefix=list(state['interval_energy']['records']);previous=row['sha256']
    require(crossing_summary(state)==meta['bore_observation'] and state['bore_observation']['crossing_ids']==ids and
            state['crossing_mass_kg']==state['bore_observation']['crossing_mass_kg']==meta['crossing_mass_kg'] and
            len(prefix)==state['step'],'first-passage final observer/complete energy ledger differs')
    extra=[Path(__file__).resolve(),CONTRACT,meta_file,rows_file,candidates_file,screens_file,parent]
    frozen=dict(before['frozen']);frozen.update({path.relative_to(core.ROOT).as_posix():core.sha(path) for path in extra})
    references.update({path:core.sha(path) for path in extra})
    return dict(pa=pa,state=state,parent=parent,parent_hash=previous,origin=parent,origin_hash=previous,
        p=before['p'],provenance=before['provenance'],seeds=seeds,prefix=prefix,frozen=frozen,references=references,
        guard_cost=before['guard_cost'],per_root=max(before['per_root'],max(row['wall_seconds'] for row in rows),
        max(row['wall_seconds'] for row in coarse)),parent_crossed_ids=list(ids),passage_window=passage_window,
        start_time_s=state['time_s'],end_time_s=state['time_s']+DURATION)


def qualified_coarse(context):
    directory=BASE/'coarse';meta_file=directory/'manifest.json';rows_file=directory/'candidate_history.json'
    meta,rows=post.read(meta_file),post.read(rows_file)
    require(meta['status']=='NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED' and
            meta['source_sha256']==context['frozen'] and meta['provenance']==context['provenance'] and
            meta['parent_sha256']==context['parent_hash'] and meta['candidate_steps']==len(rows)==PAIRS and
            meta['accepted_steps']==0 and abs(meta['accepted_time_s']-context['start_time_s'])<1e-15 and
            abs(meta['candidate_time_s']-context['end_time_s'])<1e-15 and not list(directory.glob('step_*/accepted.npz')),
            'complete nonaccepted same-parent post-passage coarse reference absent')
    references=context['references'];references.update({meta_file:core.sha(meta_file),rows_file:core.sha(rows_file)})
    previous=context['parent_hash']
    for index,row in enumerate(rows,1):
        require(row['step']==index and row['dt_s']==COARSE_DT and row['parent_sha256']==previous and
                abs(row['time_s']-(context['start_time_s']+index*COARSE_DT))<1e-15,'post-passage coarse numerical chain differs')
        finite.check_artifacts(directory,row,references)
        metric_file=directory/row['parent_crossed_response'];require(core.sha(metric_file)==row['parent_crossed_response_sha256'],'post-passage observation drift')
        references[metric_file]=row['parent_crossed_response_sha256'];previous=row['sha256']
    return rows


def parent_crossed_response(pa,state,context):
    initial=context['pa'];initial_index={int(pid):i for i,pid in enumerate(initial.particle_id)}
    index={int(pid):i for i,pid in enumerate(pa.particle_id)}
    pose=core.pose(context['p'],state['time_s']);initial_pose=core.pose(context['p'],context['state']['time_s'])
    body=(np.column_stack((pa.x,pa.y,pa.z))-np.asarray(pose['origin']))@np.asarray(pose['rotation'])
    body0=(np.column_stack((initial.x,initial.y,initial.z))-np.asarray(initial_pose['origin']))@np.asarray(initial_pose['rotation'])
    rows=[]
    for pid in context['parent_crossed_ids']:
        i,j=index[int(pid)],initial_index[int(pid)]
        r,r0=float(np.linalg.norm(body[i,:2])),float(np.linalg.norm(body0[j,:2]))
        rows.append(dict(particle_id=int(pid),plastic_strain=float(pa.ep[i]),plastic_strain_increment_from_parent=float(pa.ep[i]-initial.ep[j]),
            temperature_K=float(pa.temperature[i]),temperature_increment_from_parent_K=float(pa.temperature[i]-initial.temperature[j]),
            tool_body_position_m=body[i].tolist(),tool_radial_m=r,tool_radial_increment_from_parent_m=r-r0,
            tool_z_m=float(body[i,2]),tool_z_increment_from_parent_m=float(body[i,2]-body0[j,2])))
    return dict(classification='OBSERVATION_ONLY_NOT_DEPOSITION_OR_FORCE_INPUT',time_s=state['time_s'],
        parent_time_s=context['state']['time_s'],parent_sha256=context['parent_hash'],tool_pose=pose,
        parent_crossed_ids=context['parent_crossed_ids'],current_bore_observation=crossing_summary(state),rows=rows)


def preflight(passage_window,stage):
    context=passage_parent(passage_window);coarse=None
    if stage=='fine':
        coarse=qualified_coarse(context);context['per_root']=max(context['per_root'],max(row['wall_seconds'] for row in coarse))
    context['per_root']+=context['guard_cost']
    context['forecast']=30+1.5*PAIRS*(2 if stage=='fine' else 1)*context['per_root']
    require(context['forecast']<=450,'post-passage forecast exceeds 450 s')
    context['coarse_rows']=coarse
    return context


def worker(passage_window,stage):
    # Explicit helper injection only; inherited physics and acceptance remain
    # unchanged, and every resulting certificate retains its original report.
    event.audit_contact_shell=audit_contact_shell_finite_enclosure
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    history=[];candidates=[];screens=[]
    report=dict(stage=stage,status='RUNNING',accepted_steps=0,accepted_time_s=None,candidate_steps=0,
        physical_accuracy='ENGINEERING_PAIRWISE_TIME_SCREEN_NOT_EXPERIMENTAL_VALIDATION',
        controller='ONE_FOUR_MS_POST_PASSAGE_RESPONSE_WITH_STRICTER_OBSERVATION_SCOPE',
        coarse_reference_restart_allowed=False)
    def stop():
        try:
            core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,candidate_history=candidates,
                published_checkpoints=[str(x.relative_to(directory)) for x in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        c=preflight(passage_window,stage);pa,state,p=c['pa'],c['state'],c['p'];per_root=c['per_root'];seeds=set(c['seeds'])
        report.update(classification=core.CLASSIFICATION,source_sha256=c['frozen'],provenance=c['provenance'],
            parent_checkpoint=c['parent'].relative_to(core.ROOT).as_posix(),parent_sha256=c['parent_hash'],
            starting_time_s=state['time_s'],accepted_time_s=state['time_s'],candidate_time_s=state['time_s'],
            requested_time_s=c['end_time_s'],registered_refinement_level=1,registered_start_s=c['start_time_s'],forecast_wall_s=c['forecast'],dt_s=FINE_DT if stage=='fine' else COARSE_DT,
            parent_certificate_scope='Full qualified first-passage chain with no state or ledger reset',
            parent_crossed_ids=c['parent_crossed_ids'],passage_window=passage_window,
            endpoint_visibility_scope='ACTUAL_ENDPOINT_NNPS_SEGMENTS_ONLY',endpoint_visibility_minimum_length_m=1e-9,
            measured_guard_cost_s=c['guard_cost'],
            thread_environment={key:os.environ.get(key) for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',history)
        core.write(directory/'candidate_history.json',candidates)
        initial,initial_state=core.restore(batch.OLD/'quarters/independent_initial.npz')
        core.checkpoint(directory/'independent_initial.npz',initial,initial_state)
        core.checkpoint(directory/'continuation_initial.npz',pa,state)
        direct_file,direct_hash=c['parent'],c['parent_hash']
        for pair in range(1,PAIRS+1):
            pending=[]
            for half in range(1,3 if stage=='fine' else 2):
                require(time.perf_counter()-begin+1.5*per_root+15<=525,'remaining shoulder-window budget exhausted')
                number=2*(pair-1)+half if stage=='fine' else pair
                dt=FINE_DT if stage=='fine' else COARSE_DT
                stepdir=directory/f'step_{number:03d}';stepdir.mkdir()
                item=candidate_step(pa,state,p,dt,stepdir,c['frozen'],begin+525,per_root)
                require(item['state']['interval_energy']['records'][:len(c['prefix'])]==c['prefix'],'candidate reset the complete registered parent prefix')
                pending.append(item);pa,state=item['pa'],item['state'];seeds.update(e['particle_id'] for e in item['events'])
                observed=time.perf_counter();response_file=stepdir/'parent_crossed_response.json'
                core.write(response_file,parent_crossed_response(pa,state,c))
                item['wall_seconds']+=time.perf_counter()-observed
                event_file=stepdir/'sampled_nonbore_events.json'
                candidate_parent=(core.sha(pending[-2]['candidate']) if stage=='fine' and len(pending)>1 else direct_hash)
                row=dict(step=number,global_step=state['step'],time_s=state['time_s'],dt_s=dt,
                    parent_sha256=candidate_parent,checkpoint=f'step_{number:03d}/NOT_RESTARTABLE_candidate.npz',
                    sha256=core.sha(item['candidate']),evidence=f'step_{number:03d}/candidate_evidence.npz',
                    evidence_sha256=core.sha(item['evidence']),events=f'step_{number:03d}/sampled_nonbore_events.json',
                    events_sha256=core.sha(event_file),cumulative_shoulder_ids=sorted(seeds),
                    bore_observation=crossing_summary(state),
                    parent_crossed_response=f'step_{number:03d}/parent_crossed_response.json',
                    parent_crossed_response_sha256=core.sha(response_file),
                    wall_seconds=item['wall_seconds'],restart_allowed=False,**item['metrics'])
                row['audit_artifacts']={name:dict(path=f'step_{number:03d}/'+name,sha256=core.sha(stepdir/name))
                    for name in ('continuous_shell_audit.json','scalar_batch_equivalence.json','geometry_audit.json','solid_edge_visibility.json')}
                candidates.append(row);core.write(directory/'candidate_history.json',candidates)
                report.update(candidate_steps=len(candidates),candidate_time_s=state['time_s'])
                core.write(directory/'manifest.json',report);per_root=max(per_root,item['wall_seconds'])
            if stage=='coarse':
                event.check_references(c);require(core.sha(direct_file)==direct_hash,'coarse numerical parent drift')
                direct_file=pending[-1]['candidate'];direct_hash=core.sha(direct_file)
                print(json.dumps(dict(candidate_step=pair,time_s=state['time_s'],accepted=False)),flush=True)
                continue
            coarse_row=c['coarse_rows'][pair-1];coarse_path=directory.parent/'coarse'/coarse_row['checkpoint']
            require(core.sha(coarse_path)==coarse_row['sha256'],'coarse same-time reference drift')
            raw=core.load(str(coarse_path));coarse,coarse_state=raw['arrays']['solid'],raw['solver_data']
            require(coarse_state['restart_allowed'] is False and coarse_state['implementation_source_sha256']==c['frozen'] and
                    coarse_state['parameters']==p and coarse_state['parameter_provenance']==c['provenance'] and
                    coarse_state['interval_energy']['records'][:len(c['prefix'])]==c['prefix'] and
                    abs(coarse_state['time_s']-state['time_s'])<1e-15,'coarse reference endpoint/model differs')
            seeds.update(coarse_row['cumulative_shoulder_ids'])
            with np.load(directory.parent/'coarse'/coarse_row['evidence']) as data:
                graphs=[(data['particle_id'].copy(),data['i'].copy(),data['j'].copy())]
            graphs += [(item['pa'].particle_id,item['result']['i'],item['result']['j']) for item in pending]
            crossed=set(state['bore_observation']['crossing_ids'])|set(coarse_state['bore_observation']['crossing_ids'])
            for item in pending:crossed.update(item['state']['bore_observation']['crossing_ids'])
            selected=event.one_hop_ids(seeds|crossed,graphs)
            global_screen=post.compare(pa,coarse,p);local_screen=event.local_compare(pa,coarse,selected)
            screen=dict(pair=pair,time_s=state['time_s'],parent_sha256=c['parent_hash'],
                coarse_checkpoint_sha256=coarse_row['sha256'],cumulative_shoulder_ids=sorted(seeds),
                global_screen=global_screen,local_screen=local_screen,cumulative_crossing_ids=sorted(crossed),
                local_scope='CUMULATIVE_SHOULDER_AND_CROSSING_IDS_ONE_HOP_UNION',
                coarse_parent_crossed_response=parent_crossed_response(coarse,coarse_state,c),
                fine_parent_crossed_response=parent_crossed_response(pa,state,c),
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
                    c['frozen'],c['provenance'],direct_hash,c['parent_hash'])
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
        require(abs(state['time_s']-c['end_time_s'])<1e-15,'post-passage terminal time differs')
        if stage=='coarse':
            require(not history and not list(directory.glob('step_*/accepted.npz')),'nonaccepted coarse branch was published')
            report['status']='NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED'
        else:
            require(len(history)==2*PAIRS and len(screens)==PAIRS,'incomplete 64-step 32-pair post-passage history')
            report.update(status='POST_PASSAGE_SHORT_RESPONSE_NOT_DEPOSITION',max_plastic_strain=float(pa.ep.max()),
                final_parent_crossed_response=parent_crossed_response(pa,state,c),bore_observation=crossing_summary(state))
        event.check_references(c);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--passage-window',type=int,choices=[1,2,3],default=1)
    parser.add_argument('--stage',choices=['coarse','fine'],required=True)
    parser.add_argument('--worker',action='store_true');parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    if args.preflight_only:
        c=preflight(args.passage_window,args.stage)
        observed=parent_crossed_response(c['pa'],c['state'],c)
        require(all(row['plastic_strain_increment_from_parent']==0 and row['temperature_increment_from_parent_K']==0 and
            row['tool_radial_increment_from_parent_m']==0 and row['tool_z_increment_from_parent_m']==0 for row in observed['rows']),
            'parent observation must have exact zero increments relative to itself')
        print(json.dumps(dict(stage=args.stage,parent_sha256=c['parent_hash'],parent_time_s=c['state']['time_s'],
            parent_global_step=c['state']['step'],prefix_records=len(c['prefix']),forecast_wall_s=c['forecast'],
            requested_time_s=c['end_time_s'],registered_pairs=PAIRS,parent_crossed_ids=c['parent_crossed_ids'],physical_evaluations=0)))
        raise SystemExit(0)
    if args.worker:raise SystemExit(worker(args.passage_window,args.stage))
    logs=core.ROOT/'logs/coupled_cost_20261004/batched/post_passage_short_response'
    logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/args.stage).exists(),'refuse post-passage evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--passage-window',str(args.passage_window),'--stage',args.stage,'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-started,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

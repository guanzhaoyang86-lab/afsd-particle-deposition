"""Three bounded windows, gated by the complete same-parent 30 ms refinement.

Only registered 30->38->46->54 ms windows at the existing refinement level.
The frozen candidate helpers and explicit pair transaction are retained;
physical runs are separate CLI invocations after every parent qualifies.
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
import run_finite_shoulder_refinement as refinement

finite=refinement.finite
event,batch,core,post=finite.event,finite.batch,finite.core,finite.post
require=event.require
candidate_step,crossing_summary=finite.candidate_step,finite.crossing_summary
audit_contact_shell_finite_enclosure=finite.audit_contact_shell_finite_enclosure
BASE=batch.BASE/'refined_shoulder_continuation'
CONTRACT=core.ROOT/'docs/REFINED_SHOULDER_CONTINUATION_CONTRACT_20261004.md'
COARSE_DT,FINE_DT,PAIRS=.000125,.0000625,64
WINDOWS={1:(.030,.038),2:(.038,.046),3:(.046,.054)}


def qualified_coarse(number,parent_hash,p,provenance,frozen,references):
    directory=BASE/f'window_{number}'/'coarse';meta_file=directory/'manifest.json';rows_file=directory/'candidate_history.json'
    meta,rows=post.read(meta_file),post.read(rows_file);start,end=WINDOWS[number]
    require(meta['status']=='NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED' and
            meta['source_sha256']==frozen and meta['provenance']==provenance and
            meta['parent_sha256']==parent_hash and meta['candidate_steps']==len(rows)==PAIRS and
            meta['accepted_steps']==0 and abs(meta['accepted_time_s']-start)<1e-15 and
            abs(meta['candidate_time_s']-end)<1e-15 and not list(directory.glob('step_*/accepted.npz')),
            'complete uncommitted same-origin 64-step coarse reference absent')
    references.update({meta_file:core.sha(meta_file),rows_file:core.sha(rows_file)})
    previous=parent_hash
    for index,row in enumerate(rows,1):
        require(row['step']==index and row['dt_s']==COARSE_DT and row['parent_sha256']==previous and
                abs(row['time_s']-(start+index*COARSE_DT))<1e-15,'continuation coarse numerical chain differs')
        finite.check_artifacts(directory,row,references);previous=row['sha256']
    return rows


def complete_fine(directory,before,start,end,expected_frozen,coarse,references):
    """Read and qualify the entire 128-step parent window, never a failed prefix."""
    meta_file=directory/'manifest.json';rows_file=directory/'accepted_history.json'
    candidates_file=directory/'candidate_history.json';screens_file=directory/'screen_history.json'
    meta,rows=post.read(meta_file),post.read(rows_file)
    candidates,screens=post.read(candidates_file),post.read(screens_file)
    require(meta['status']=='PASS_WINDOW_NOT_TRANSPORT' and meta['accepted_steps']==len(rows)==2*PAIRS and
            meta['candidate_steps']==len(candidates)==2*PAIRS and len(screens)==PAIRS and
            meta['source_sha256']==expected_frozen and meta['provenance']==before['provenance'] and
            meta['parent_sha256']==before['parent_hash'] and abs(meta['starting_time_s']-start)<1e-15 and
            abs(meta['accepted_time_s']-end)<1e-15 and meta['crossing_mass_kg']==0 and
            len(list(directory.glob('step_*/accepted.npz')))==2*PAIRS,
            'parent window must completely PASS without crossing; failed/incomplete/first-passage parent forbidden')
    references.update({path:core.sha(path) for path in (meta_file,rows_file,candidates_file,screens_file)})
    baseline=core.frozen_sources();prefix=list(before['state']['interval_energy']['records'])
    previous=before['parent_hash'];seeds=set(before['seeds'])
    for index,(row,candidate) in enumerate(zip(rows,candidates),1):
        parent=directory/row['checkpoint']
        require(row['step']==candidate['step']==index and row['dt_s']==candidate['dt_s']==FINE_DT and
                core.sha(parent)==row['sha256'] and row['parent_sha256']==previous and
                candidate['parent_sha256']==(previous if index%2 else candidates[index-2]['sha256']),
                'parent fine direct/candidate chain differs')
        pa,state=core.restore(parent)
        require(state['step']==row['global_step']==candidate['global_step']==before['state']['step']+index and
                state['time_s']==row['time_s']==candidate['time_s'] and abs(state['time_s']-(start+index*FINE_DT))<1e-15 and
                state['source_sha256']==baseline and state['implementation_source_sha256']==expected_frozen and
                state['parameters']==before['p'] and state['parameter_provenance']==before['provenance'] and
                state['continuation_parent_sha256']==previous and state['branch_origin_sha256']==before['parent_hash'] and
                state['interval_energy']['records'][:-1]==prefix and not state['bore_observation']['crossing_ids'] and
                state['crossing_mass_kg']==0 and state['bore_observation']['crossing_mass_kg']==0,
                'parent fine source/parameters/prefix/time/crossing differs')
        batch.make_trial(pa,state,before['p'])  # Constructor checks serialized history only; no evaluate.
        finite.check_artifacts(directory,candidate,references)
        core.verify_snapshot(directory/candidate['checkpoint'],pa,state)
        saved=core.load(str(directory/candidate['checkpoint']))['solver_data']
        require(saved['restart_allowed'] is False and saved['implementation_source_sha256']==expected_frozen and
                saved['source_sha256']==baseline and saved['parameters']==before['p'] and
                saved['parameter_provenance']==before['provenance'] and
                saved['interval_energy']==state['interval_energy'] and saved['bore_observation']==state['bore_observation'],
                'parent uncommitted candidate source/parameters/ledger/observer differs')
        pair=(index-1)//2+1;entry=screens[pair-1];screen_file=directory/entry['screen']
        require(entry['pair']==pair and core.sha(screen_file)==entry['sha256']==row['temporal_screen_sha256']==state['window_temporal_screen_sha256'],
                'parent pair screen hash differs')
        screen=post.read(screen_file)
        require(screen['pair']==pair and abs(screen['time_s']-(start+pair*COARSE_DT))<1e-15 and
                screen['parent_sha256']==before['parent_hash'] and screen['global_screen']['pass'] and screen['local_screen']['passed'] and
                screen['crossing_ids_match'] and not screen['fine_bore_observation']['crossing_ids'] and
                not screen['coarse_bore_observation']['crossing_ids'] and
                screen['fine_bore_observation']['crossing_mass_kg']==screen['coarse_bore_observation']['crossing_mass_kg']==0 and
                screen['coarse_checkpoint_sha256']==coarse[pair-1]['sha256']==row['coarse_reference_sha256']==state['window_coarse_reference_sha256'],
                'parent pair lacks complete equal-time qualification without crossing')
        seeds.update(screen['cumulative_shoulder_ids']);references[screen_file]=entry['sha256'];references[parent]=row['sha256']
        prefix=list(state['interval_energy']['records']);previous=row['sha256']
    return dict(pa=pa,state=state,parent=parent,parent_hash=previous,seeds=seeds,p=before['p'],provenance=before['provenance'],
        per_root=max(before['per_root'],max(row['wall_seconds'] for row in rows),max(row['wall_seconds'] for row in coarse)))


def refinement_parent():
    # This existing gate preserves every earlier failure and all finite/visibility
    # sources, while returning the original 22 ms parent and completed coarse.
    before=refinement.preflight('fine');references=dict(before['references'])
    directory=refinement.BASE/'fine'
    parent=complete_fine(directory,before,.022,.030,before['frozen'],before['coarse_rows'],references)
    require(parent['state']['step']==281 and len(parent['state']['interval_energy']['records'])==281,
            '30 ms origin must retain all 153+128 accepted energy records')
    extra=[Path(__file__).resolve(),CONTRACT,directory/'manifest.json',directory/'accepted_history.json',
           directory/'candidate_history.json',directory/'screen_history.json',parent['parent']]
    frozen=dict(before['frozen']);frozen.update({path.relative_to(core.ROOT).as_posix():core.sha(path) for path in extra})
    references.update({path:core.sha(path) for path in extra})
    return parent,frozen,references,before['guard_cost']


def registered_parent(number,origin,frozen,references):
    if number==1:return origin
    before=registered_parent(number-1,origin,frozen,references)
    coarse=qualified_coarse(number-1,before['parent_hash'],before['p'],before['provenance'],frozen,references)
    return complete_fine(BASE/f'window_{number-1}'/'fine',before,*WINDOWS[number-1],frozen,coarse,references)


def preflight(window_number,stage):
    require(window_number in WINDOWS,'only registered continuation windows 1,2,3 exist')
    origin,frozen,references,guard_cost=refinement_parent()
    parent=registered_parent(window_number,origin,frozen,references);start,end=WINDOWS[window_number]
    state=parent['state'];require(abs(state['time_s']-start)<1e-15 and not state['bore_observation']['crossing_ids'],
                               'continuation requires its complete registered parent and no earlier crossing')
    per_root=parent['per_root'];coarse_rows=None
    if stage=='fine':
        coarse_rows=qualified_coarse(window_number,parent['parent_hash'],parent['p'],parent['provenance'],frozen,references)
        per_root=max(per_root,max(row['wall_seconds'] for row in coarse_rows))
    per_root+=guard_cost
    forecast=30+1.5*PAIRS*(2 if stage=='fine' else 1)*per_root
    require(forecast<=450,'registered continuation forecast exceeds 450 s')
    return dict(parent,frozen=frozen,references=references,origin=parent['parent'],origin_hash=parent['parent_hash'],
        prefix=list(state['interval_energy']['records']),per_root=per_root,guard_cost=guard_cost,
        coarse_rows=coarse_rows,forecast=forecast,start_time_s=start,end_time_s=end,window_number=window_number)


def worker(window_number,stage):
    # Explicit helper injection only; inherited physics and acceptance remain
    # unchanged, and every resulting certificate retains its original report.
    event.audit_contact_shell=audit_contact_shell_finite_enclosure
    directory=BASE/f'window_{window_number}'/stage;directory.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    history=[];candidates=[];screens=[];first_passage=False
    report=dict(stage=stage,status='RUNNING',accepted_steps=0,accepted_time_s=None,candidate_steps=0,
        physical_accuracy='ENGINEERING_PAIRWISE_TIME_SCREEN_NOT_EXPERIMENTAL_VALIDATION',
        controller='SAME_REFINEMENT_LEVEL_BOUNDED_CONTINUATION_FROM_30MS',
        coarse_reference_restart_allowed=False)
    def stop():
        try:
            core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,candidate_history=candidates,
                published_checkpoints=[str(x.relative_to(directory)) for x in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        c=preflight(window_number,stage);pa,state,p=c['pa'],c['state'],c['p'];per_root=c['per_root'];seeds=set(c['seeds'])
        report.update(classification=core.CLASSIFICATION,source_sha256=c['frozen'],provenance=c['provenance'],
            parent_checkpoint=c['parent'].relative_to(core.ROOT).as_posix(),parent_sha256=c['parent_hash'],
            starting_time_s=state['time_s'],accepted_time_s=state['time_s'],candidate_time_s=state['time_s'],
            requested_time_s=c['end_time_s'],registered_window=window_number,registered_start_s=c['start_time_s'],forecast_wall_s=c['forecast'],dt_s=FINE_DT if stage=='fine' else COARSE_DT,
            parent_certificate_scope='Complete 281-step 30ms refinement chain plus fully screened registered continuation windows',
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
            selected=event.one_hop_ids(seeds,graphs)
            global_screen=post.compare(pa,coarse,p);local_screen=event.local_compare(pa,coarse,selected)
            screen=dict(pair=pair,time_s=state['time_s'],parent_sha256=c['parent_hash'],
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
            if set(state['bore_observation']['crossing_ids'])-set(c['state']['bore_observation']['crossing_ids']):
                first_passage=True
                report.update(status='PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION',
                    first_passage_pair=pair,bore_observation=crossing_summary(state))
                break
        if not first_passage:
            require(abs(state['time_s']-c['end_time_s'])<1e-15,'window terminal time differs')
        if stage=='coarse':
            require(not history and not list(directory.glob('step_*/accepted.npz')),'nonaccepted coarse branch was published')
            report['status']='NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED'
        elif not first_passage:
            require(len(history)==2*PAIRS and len(screens)==PAIRS,'incomplete 128-step 64-pair continuation history')
            report.update(status='PASS_WINDOW_NOT_TRANSPORT',max_plastic_strain=float(pa.ep.max()))
        event.check_references(c);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--window',type=int,choices=[1,2,3],default=1)
    parser.add_argument('--stage',choices=['coarse','fine'],required=True)
    parser.add_argument('--worker',action='store_true');parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    if args.preflight_only:
        c=preflight(args.window,args.stage)
        print(json.dumps(dict(window=args.window,stage=args.stage,parent_sha256=c['parent_hash'],parent_time_s=c['state']['time_s'],
            parent_global_step=c['state']['step'],prefix_records=len(c['prefix']),forecast_wall_s=c['forecast'],registered_pairs=PAIRS,physical_evaluations=0)))
        raise SystemExit(0)
    if args.worker:raise SystemExit(worker(args.window,args.stage))
    logs=core.ROOT/f'logs/coupled_cost_20261004/batched/refined_shoulder_continuation/window_{args.window}'
    logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/f'window_{args.window}'/args.stage).exists(),'refuse refined continuation evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--window',str(args.window),'--stage',args.stage,'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-started,limit_s=600))
    print(json.dumps(dict(window=args.window,stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

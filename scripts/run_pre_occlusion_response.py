"""One 2 ms fine comparison restricted to the qualified pre-occlusion range.

The old four-ms reference remains FAILED. Its first 16 uncommitted candidates
are read-only comparisons; all new states start from the true 37.375 ms parent.
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
import run_post_passage_short_response as postpass

finite=postpass.finite
event,batch,core,post=postpass.event,postpass.batch,postpass.core,postpass.post
require=event.require
candidate_step,crossing_summary=postpass.candidate_step,postpass.crossing_summary
parent_crossed_response=postpass.parent_crossed_response
audit_contact_shell_finite_enclosure=postpass.audit_contact_shell_finite_enclosure
BASE=batch.BASE/'pre_occlusion_response'
OLD_COARSE=postpass.BASE/'coarse'
CONTRACT=core.ROOT/'docs/PRE_OCCLUSION_RESPONSE_CONTRACT_20261004.md'
COARSE_DT,FINE_DT,PAIRS=.000125,.0000625,16
PARENT_SHA='9f78a757bb9897e97598b09bee5842f29f6d6d7a2d217dfb39add0fb4a096b67'


def preflight():
    context=postpass.passage_parent(1);coarse_frozen=dict(context['frozen']);references=context['references']
    require(context['parent_hash']==PARENT_SHA and context['state']['step']==len(context['prefix'])==399 and
            abs(context['state']['time_s']-.037375)<1e-15,'registered 37.375 ms accepted parent differs')
    meta_file=OLD_COARSE/'manifest.json';rows_file=OLD_COARSE/'candidate_history.json'
    meta,rows=post.read(meta_file),post.read(rows_file)
    require(meta['status']=='FAILED_NO_NEW_ACCEPTED_STATE' and meta['accepted_steps']==0 and
            meta['candidate_steps']==len(rows)==19 and abs(meta['accepted_time_s']-.037375)<1e-15 and
            abs(meta['candidate_time_s']-.03975)<1e-15 and meta['source_sha256']==coarse_frozen and
            meta['parent_sha256']==PARENT_SHA and meta['provenance']==context['provenance'] and
            not list(OLD_COARSE.glob('step_*/accepted.npz')) and
            len(list(OLD_COARSE.glob('step_*/NOT_RESTARTABLE_candidate.npz')))==20,
            'old four-ms failure or its 19 qualified plus one rejected candidates changed')
    extra=[Path(__file__).resolve(),CONTRACT,meta_file,rows_file]
    previous=PARENT_SHA;prefix=list(context['prefix'])
    for index,row in enumerate(rows,1):
        finite.check_artifacts(OLD_COARSE,row,references)
        path=OLD_COARSE/row['checkpoint'];raw=core.load(str(path));state=raw['solver_data']
        require(row['step']==index and row['global_step']==state['step']==399+index and row['dt_s']==COARSE_DT and
                row['parent_sha256']==previous and row['time_s']==state['time_s'] and
                abs(row['time_s']-(.037375+index*COARSE_DT))<1e-15 and
                state['restart_allowed'] is False and state['source_sha256']==core.frozen_sources() and
                state['implementation_source_sha256']==coarse_frozen and state['parameters']==context['p'] and
                state['parameter_provenance']==context['provenance'] and state['interval_energy']['records'][:-1]==prefix,
                'old numerical reference source, parameters, prefix or chain differs')
        observation_file=OLD_COARSE/row['parent_crossed_response']
        require(core.sha(observation_file)==row['parent_crossed_response_sha256'],'old observation drift')
        references[observation_file]=row['parent_crossed_response_sha256'];extra.append(path)
        previous=row['sha256'];prefix=list(state['interval_energy']['records'])
    rejected=OLD_COARSE/'step_020';blocked=post.read(rejected/'solid_edge_visibility.json')
    require(blocked['status']=='FAIL_SOLID_OCCLUDED_ENDPOINT_EDGE' and blocked['occluded_edge_count']==1 and
            blocked['max_occluded_interval_length_m']>1e-9 and blocked['max_inset_occluded_interval_length_m']>1e-9 and
            abs(blocked['endpoint_time_s']-.039875)<1e-15 and
            {blocked['occluded_edges'][0]['particle_id_i'],blocked['occluded_edges'][0]['particle_id_j']}=={219,251} and
            core.sha(rejected/'NOT_RESTARTABLE_candidate.npz')==blocked['candidate_sha256'] and
            core.sha(rejected/'candidate_evidence.npz')==blocked['evidence_sha256'],
            'actual 39.875 ms solid-occluded candidate evidence differs')
    failed=core.load(str(rejected/'NOT_RESTARTABLE_candidate.npz'))['solver_data']
    require(failed['restart_allowed'] is False and failed['implementation_source_sha256']==coarse_frozen and
            failed['parameters']==context['p'] and failed['parameter_provenance']==context['provenance'] and
            failed['interval_energy']['records'][:-1]==prefix and failed['step']==419 and
            abs(failed['time_s']-.039875)<1e-15,'rejected candidate was changed or made restartable')
    extra.extend(rejected/name for name in ('NOT_RESTARTABLE_candidate.npz','candidate_evidence.npz',
        'FULL_UNCOMMITTED_candidate.npz','FULL_UNCOMMITTED_candidate.json','solid_edge_visibility.json',
        'continuous_shell_audit.json','geometry_audit.json','scalar_batch_equivalence.json','energy.json','sampled_nonbore_events.json'))
    frozen=dict(coarse_frozen);frozen.update({path.relative_to(core.ROOT).as_posix():core.sha(path) for path in extra})
    references.update({path:core.sha(path) for path in extra})
    per_root=max(context['per_root'],max(row['wall_seconds'] for row in rows))+context['guard_cost']
    forecast=30+1.5*(2*PAIRS)*per_root;require(forecast<=450,'pre-occlusion two-ms forecast exceeds 450 s')
    require(abs(rows[PAIRS-1]['time_s']-.039375)<1e-15,'reference prefix does not end at registered 39.375 ms')
    context.update(frozen=frozen,coarse_frozen=coarse_frozen,per_root=per_root,forecast=forecast,
        coarse_rows=rows[:PAIRS],start_time_s=.037375,end_time_s=.039375)
    return context


def worker():
    stage='fine';passage_window=1
    # Explicit helper injection only; inherited physics and acceptance remain
    # unchanged, and every resulting certificate retains its original report.
    event.audit_contact_shell=audit_contact_shell_finite_enclosure
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    history=[];candidates=[];screens=[]
    report=dict(stage=stage,status='RUNNING',accepted_steps=0,accepted_time_s=None,candidate_steps=0,
        physical_accuracy='ENGINEERING_PAIRWISE_TIME_SCREEN_NOT_EXPERIMENTAL_VALIDATION',
        controller='ONE_PRE_OCCLUSION_TWO_MS_COMPARISON_OLD_FAILURE_PRESERVED',
        coarse_reference_restart_allowed=False)
    def stop():
        try:
            core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,candidate_history=candidates,
                published_checkpoints=[str(x.relative_to(directory)) for x in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        c=preflight();pa,state,p=c['pa'],c['state'],c['p'];per_root=c['per_root'];seeds=set(c['seeds'])
        report.update(classification=core.CLASSIFICATION,source_sha256=c['frozen'],provenance=c['provenance'],
            parent_checkpoint=c['parent'].relative_to(core.ROOT).as_posix(),parent_sha256=c['parent_hash'],
            starting_time_s=state['time_s'],accepted_time_s=state['time_s'],candidate_time_s=state['time_s'],
            requested_time_s=c['end_time_s'],registered_refinement_level=1,registered_start_s=c['start_time_s'],forecast_wall_s=c['forecast'],dt_s=FINE_DT if stage=='fine' else COARSE_DT,
            parent_certificate_scope='Full qualified first-passage chain with no state or ledger reset',
            parent_crossed_ids=c['parent_crossed_ids'],passage_window=passage_window,
            coarse_reference_directory=OLD_COARSE.relative_to(core.ROOT).as_posix(),
            coarse_reference_source_sha256=c['coarse_frozen'],original_four_ms_status='FAILED_NO_NEW_ACCEPTED_STATE',
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
            coarse_row=c['coarse_rows'][pair-1];coarse_path=OLD_COARSE/coarse_row['checkpoint']
            require(core.sha(coarse_path)==coarse_row['sha256'],'coarse same-time reference drift')
            raw=core.load(str(coarse_path));coarse,coarse_state=raw['arrays']['solid'],raw['solver_data']
            require(coarse_state['restart_allowed'] is False and coarse_state['implementation_source_sha256']==c['coarse_frozen'] and
                    coarse_state['parameters']==p and coarse_state['parameter_provenance']==c['provenance'] and
                    coarse_state['interval_energy']['records'][:len(c['prefix'])]==c['prefix'] and
                    abs(coarse_state['time_s']-state['time_s'])<1e-15,'coarse reference endpoint/model differs')
            seeds.update(coarse_row['cumulative_shoulder_ids'])
            with np.load(OLD_COARSE/coarse_row['evidence']) as data:
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
            require(len(history)==2*PAIRS and len(screens)==PAIRS,'incomplete 32-step 16-pair pre-occlusion history')
            report.update(status='PASS_PRE_OCCLUSION_2MS_RESPONSE_NOT_DEPOSITION',max_plastic_strain=float(pa.ep.max()),
                final_parent_crossed_response=parent_crossed_response(pa,state,c),bore_observation=crossing_summary(state))
        event.check_references(c);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--worker',action='store_true')
    parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    if args.preflight_only:
        c=preflight()
        print(json.dumps(dict(stage='fine',parent_sha256=c['parent_hash'],parent_time_s=c['state']['time_s'],
            parent_global_step=c['state']['step'],prefix_records=len(c['prefix']),forecast_wall_s=c['forecast'],
            requested_time_s=c['end_time_s'],registered_pairs=PAIRS,parent_crossed_ids=c['parent_crossed_ids'],
            original_four_ms_status='FAILED_NO_NEW_ACCEPTED_STATE',physical_evaluations=0)))
        raise SystemExit(0)
    if args.worker:raise SystemExit(worker())
    logs=core.ROOT/'logs/coupled_cost_20261004/batched/pre_occlusion_response'
    logs.mkdir(parents=True,exist_ok=True);log=logs/'fine.log'
    require(not log.exists() and not (BASE/'fine').exists(),'refuse pre-occlusion evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/'fine_supervision.json',dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-started,limit_s=600))
    print(json.dumps(dict(stage='fine',exit_code=code,log=str(log))));raise SystemExit(code)

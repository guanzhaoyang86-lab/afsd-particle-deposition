"""One registered whole-window refinement from the unchanged 22 ms parent.

The frozen finite-window candidate helpers are reused without changing them.
The pair transaction is copied explicitly with 64 pairs, 0.125/0.0625 ms steps.
There is no recursive refinement, continuation from the failed prefix, or
change of a physical equation, acceptance threshold, or event scope.
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
import run_finite_shoulder_windows as finite

event,batch,core,post=finite.event,finite.batch,finite.core,finite.post
require=event.require
candidate_step,crossing_summary=finite.candidate_step,finite.crossing_summary
audit_contact_shell_finite_enclosure=finite.audit_contact_shell_finite_enclosure
BASE=batch.BASE/'finite_shoulder_refinement'
CONTRACT=core.ROOT/'docs/FINITE_SHOULDER_REFINEMENT_CONTRACT_20261004.md'
OLD=finite.BASE/'window_1'
COARSE_DT,FINE_DT,PAIRS=.000125,.0000625,64
START,END=.022,.030


def old_failure_gate(context):
    """Keep the actual failed window and accepted 26-step prefix immutable."""
    references=dict(context['references']);directory=OLD/'fine'
    meta_file=directory/'manifest.json';rows_file=directory/'accepted_history.json'
    candidates_file=directory/'candidate_history.json';screens_file=directory/'screen_history.json'
    failed_screen_file=directory/'pair_014_temporal_screen.json'
    meta,rows=post.read(meta_file),post.read(rows_file)
    candidates,screens=post.read(candidates_file),post.read(screens_file)
    failed=post.read(failed_screen_file)
    require(meta['status']=='FAILED_WITH_ACCEPTED_PREFIX' and meta['accepted_steps']==len(rows)==26 and
            meta['candidate_steps']==len(candidates)==28 and len(screens)==13 and
            abs(meta['accepted_time_s']-.02525)<1e-15 and abs(meta['candidate_time_s']-.0255)<1e-15 and
            meta['source_sha256']==context['frozen'] and meta['provenance']==context['provenance'] and
            meta['parent_sha256']==context['parent_hash'] and meta['last_screen']==failed and
            failed['pair']==14 and abs(failed['time_s']-.0255)<1e-15 and
            failed['global_screen']['pass'] and not failed['local_screen']['passed'] and failed['crossing_ids_match'] and
            not failed['coarse_bore_observation']['crossing_ids'] and not failed['fine_bore_observation']['crossing_ids'] and
            any(row['particle_id']==275 and row['stress_relative']>.1 for row in failed['local_screen']['rows']),
            'original window failure and its 26-step accepted prefix differ')
    old_coarse=finite.qualified_coarse(1,context['parent_hash'],context['p'],context['provenance'],context['frozen'],references)
    require(failed['coarse_checkpoint_sha256']==old_coarse[13]['sha256'] and
            len(list(directory.glob('step_*/accepted.npz')))==26 and
            not (directory/'step_027/accepted.npz').exists() and not (directory/'step_028/accepted.npz').exists(),
            'original failed pair was committed or its coarse reference changed')
    previous=context['parent_hash']
    for index,row in enumerate(rows,1):
        path=directory/row['checkpoint']
        require(row['step']==index and row['global_step']==153+index and row['dt_s']==.000125 and
                abs(row['time_s']-(START+index*.000125))<1e-15 and row['parent_sha256']==previous and
                core.sha(path)==row['sha256'],'original accepted prefix chain changed')
        references[path]=row['sha256'];previous=row['sha256']
        entry=screens[(index-1)//2];screen_file=directory/entry['screen'];screen=post.read(screen_file)
        require(core.sha(screen_file)==entry['sha256']==row['temporal_screen_sha256'] and
                screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'],
                'original accepted pair qualification changed')
        references[screen_file]=entry['sha256']
    for index,row in enumerate(candidates,1):
        finite.check_artifacts(directory,row,references)
        expected_parent=(context['parent_hash'] if index==1 else
            candidates[index-2]['sha256'] if index%2==0 else rows[index-2]['sha256'])
        require(row['step']==index and row['global_step']==153+index and row['dt_s']==.000125 and
                row['parent_sha256']==expected_parent and abs(row['time_s']-(START+index*.000125))<1e-15,
                'original candidate chain differs')
    for row in candidates[-2:]:
        raw=core.load(str(directory/row['checkpoint']));state=raw['solver_data']
        require(state['restart_allowed'] is False and state['implementation_source_sha256']==context['frozen'] and
                state['parameters']==context['p'] and state['parameter_provenance']==context['provenance'] and
                state['interval_energy']['records'][:153]==context['prefix'],
                'old failed candidate was changed or made restartable')
    extra=[Path(__file__).resolve(),CONTRACT,meta_file,rows_file,candidates_file,screens_file,failed_screen_file,
           OLD/'coarse/manifest.json',OLD/'coarse/candidate_history.json',directory/rows[-1]['checkpoint']]
    for row in candidates[-2:]:
        extra.extend(directory/row[key] for key in ('checkpoint','evidence','events'))
        extra.extend(directory/artifact['path'] for artifact in row['audit_artifacts'].values())
    frozen=dict(context['frozen']);frozen.update({path.relative_to(core.ROOT).as_posix():core.sha(path) for path in extra})
    references.update({path:core.sha(path) for path in extra})
    per_root=max(context['per_root'],max(row['wall_seconds'] for row in old_coarse),
                 max(row['wall_seconds'] for row in candidates),max(row['wall_seconds'] for row in rows))
    return frozen,references,per_root


def qualified_coarse(context,frozen,references):
    directory=BASE/'coarse';meta_file=directory/'manifest.json';rows_file=directory/'candidate_history.json'
    meta,rows=post.read(meta_file),post.read(rows_file)
    require(meta['status']=='NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED' and
            meta['source_sha256']==frozen and meta['provenance']==context['provenance'] and
            meta['parent_sha256']==context['parent_hash'] and meta['candidate_steps']==len(rows)==PAIRS and
            meta['accepted_steps']==0 and abs(meta['accepted_time_s']-START)<1e-15 and
            abs(meta['candidate_time_s']-END)<1e-15 and not list(directory.glob('step_*/accepted.npz')),
            'complete uncommitted 64-step same-parent refinement coarse absent')
    references.update({meta_file:core.sha(meta_file),rows_file:core.sha(rows_file)})
    previous=context['parent_hash']
    for index,row in enumerate(rows,1):
        require(row['step']==index and row['global_step']==153+index and row['dt_s']==COARSE_DT and
                row['parent_sha256']==previous and abs(row['time_s']-(START+index*COARSE_DT))<1e-15,
                'refinement coarse chain/time differs')
        finite.check_artifacts(directory,row,references);previous=row['sha256']
    return rows


def preflight(stage):
    context=finite.preflight(1,'coarse')
    require(context['state']['step']==153 and len(context['prefix'])==153 and
            abs(context['state']['time_s']-START)<1e-15 and context['parent_hash']==finite.PARENT_SHA,
            'refinement requires original whole-window 22 ms parent')
    frozen,references,per_root=old_failure_gate(context);coarse_rows=None
    if stage=='fine':
        coarse_rows=qualified_coarse(context,frozen,references)
        per_root=max(per_root,max(row['wall_seconds'] for row in coarse_rows))
    per_root+=context['guard_cost']
    forecast=30+1.5*PAIRS*(2 if stage=='fine' else 1)*per_root
    require(forecast<=450,'single-level refinement forecast exceeds 450 s')
    context.update(frozen=frozen,references=references,per_root=per_root,forecast=forecast,
        coarse_rows=coarse_rows,start_time_s=START,end_time_s=END)
    return context


def worker(stage):
    # Explicit helper injection only; inherited physics and acceptance remain
    # unchanged, and every resulting certificate retains its original report.
    event.audit_contact_shell=audit_contact_shell_finite_enclosure
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False);begin=time.perf_counter()
    history=[];candidates=[];screens=[];first_passage=False
    report=dict(stage=stage,status='RUNNING',accepted_steps=0,accepted_time_s=None,candidate_steps=0,
        physical_accuracy='ENGINEERING_PAIRWISE_TIME_SCREEN_NOT_EXPERIMENTAL_VALIDATION',
        controller='ONE_LEVEL_WHOLE_WINDOW_REFINEMENT_WITH_FROZEN_FINITE_HELPERS',
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
            parent_checkpoint=c['parent'].relative_to(core.ROOT).as_posix(),parent_sha256=c['parent_hash'],
            starting_time_s=state['time_s'],accepted_time_s=state['time_s'],candidate_time_s=state['time_s'],
            requested_time_s=c['end_time_s'],registered_refinement_level=1,registered_start_s=c['start_time_s'],forecast_wall_s=c['forecast'],dt_s=FINE_DT if stage=='fine' else COARSE_DT,
            parent_certificate_scope='Original 153 segments; no records taken from old post-22ms failed window',
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
            require(len(history)==2*PAIRS and len(screens)==PAIRS,'incomplete 128-step 64-pair refinement history')
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
        c=preflight(args.stage)
        print(json.dumps(dict(stage=args.stage,parent_sha256=c['parent_hash'],parent_time_s=c['state']['time_s'],
            parent_global_step=c['state']['step'],prefix_records=len(c['prefix']),forecast_wall_s=c['forecast'],
            registered_pairs=PAIRS,coarse_dt_s=COARSE_DT,fine_dt_s=FINE_DT,physical_evaluations=0)))
        raise SystemExit(0)
    if args.worker:raise SystemExit(worker(args.stage))
    logs=core.ROOT/'logs/coupled_cost_20261004/batched/finite_shoulder_refinement'
    logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/args.stage).exists(),'refuse refinement evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-started,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

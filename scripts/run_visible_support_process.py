"""One registered visible-support window; immutable evidence and bounded roots.

Cost/coarse are never restartable. Fine commits only after same-parent pair
screens, scalar-map identity and independent whole-path geometry certificates.
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
import run_first_shoulder_event as event
import visible_quasistatic as model
import visibility_equilibrium as solver
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from contact_shell_finite_enclosure import audit_contact_shell_finite_enclosure
from visible_support_path import audit_visible_support_path, CERTIFIED

core,post,batch=event.core,event.post,event.batch
ROOT=core.ROOT
BASE=ROOT/'outputs/visible_support_20261004/process_v1'
PARENT=ROOT/'outputs/coupled_cost_20261004/batched/pre_occlusion_response/fine/step_032/accepted.npz'
PARENT_SHA='553f57b1dd5f0555d5a3726f7942705bd485ffd688ef3bf34f16ff342fd9a804'
SUMMARY=ROOT/'outputs/coupled_cost_20261004/batched/readback_pre_occlusion_39_375ms/summary.json'
SUMMARY_SHA='4632c3cfc7430ad00942aeea58edb227607333c5b078610ca66af1838595b600'
QUAL=ROOT/'outputs/visible_support_20261004/qualification_v3/manifest.json'
OPERATORS=ROOT/'outputs/visible_graph_20261004/operators_v1/summary.json'
GEOMETRY=ROOT/'outputs/coupled_cost_20261004/visible_support_geometry_saved.json'
PATH_GATE=ROOT/'outputs/coupled_cost_20261004/visible_support_last_accepted_path_v2.json'
PATH_COMPONENTS=ROOT/'outputs/coupled_cost_20261004/visible_support_geometry_components_v2.json'
COMPONENT_GATE=ROOT/'outputs/visible_support_20261004/component_receipt.json'
START,END,COARSE_DT,FINE_DT,PAIRS=.039375,.040375,.000125,.0000625,8
require=event.require


class RecordedTrial(model.VisibleSupportTrial):
    """Retain last actual complete map for failure evidence, with no reevaluation."""
    def evaluate(self,displacement,dt):
        result=super().evaluate(displacement,dt)
        self.last_complete_result=result
        return result


def preflight(stage):
    require(core.sha(PARENT)==PARENT_SHA and core.sha(SUMMARY)==SUMMARY_SHA,'registered true parent changed')
    pa,state=core.restore(PARENT);p,provenance=core.parameters();baseline=core.frozen_sources()
    require(state['step']==431 and abs(state['time_s']-START)<1e-15 and state['source_sha256']==baseline and
        state['parameters']==p and state['parameter_provenance']==provenance,'parent fields or provenance differ')
    prior=post.read(SUMMARY)
    require(prior['accepted_steps']==431 and prior['inherited_interval_never_reset'],'qualified old prefix absent')
    bound=dict(prior['pre_occlusion_source_and_evidence_sha256'])
    post.check_sources(bound);post.check_sources(state['implementation_source_sha256'])
    bound.update(state['implementation_source_sha256']);bound.update(baseline)
    qual,operators,geometry,pathgate,components,pathcomponents=[post.read(v) for v in (QUAL,OPERATORS,GEOMETRY,PATH_GATE,COMPONENT_GATE,PATH_COMPONENTS)]
    require(qual['status']=='PASS_DIAGNOSTICS_NOT_PROCESS_ACCEPTANCE' and qual['accepted_time_advanced_s']==0.,'same-map gate absent')
    require(operators['status']=='PASS_REGISTERED_COMPONENT_BENCHMARKS_ONLY' and operators['exit_code']==0,'operator gate absent')
    require(geometry['status']=='PASS_SAVED_ENDPOINT_GRAPH_IDENTITY_AND_OLD_BLOCKED_EDGE_DIAGNOSTIC','endpoint geometry gate absent')
    require(pathgate['status']==CERTIFIED and abs(pathgate['end_time_s']-START)<1e-15 and pathgate['physical_solver_evaluations']==0,'saved path gate absent')
    require(pathcomponents['status']=='PASS_GEOMETRY_COMPONENTS' and pathcomponents['component_test_count']==28 and
        all(row['exit_code']==0 and row['pass_receipt'] for row in pathcomponents['rows']),'whole-path component gate absent')
    require(components['status']=='PASS_VISIBLE_GEOMETRY_AND_SEARCH_COMPONENTS_ONLY','geometry/search receipt absent')
    for value in (qual,operators,geometry,pathgate,components,pathcomponents):
        post.check_sources(value['source_sha256']);bound.update(value['source_sha256'])
        if 'input_sha256' in value:post.check_sources(value['input_sha256']);bound.update(value['input_sha256'])
    extra=[Path(__file__),PARENT,SUMMARY,QUAL,OPERATORS,GEOMETRY,PATH_GATE,COMPONENT_GATE,PATH_COMPONENTS,
        ROOT/'scripts/visible_quasistatic.py',ROOT/'scripts/visible_support.py',ROOT/'scripts/visible_support_path.py',
        ROOT/'scripts/visible_support_numerics.py',ROOT/'scripts/run_batched_process.py',ROOT/'scripts/run_discrete_post_start.py',
        ROOT/'scripts/visibility_equilibrium.py',ROOT/'scripts/run_first_shoulder_event.py',
        ROOT/'scripts/contact_shell_finite_enclosure.py',ROOT/'docs/VISIBLE_SUPPORT_PROCESS_CONTRACT_20261004.md']
    bound.update({path.relative_to(ROOT).as_posix():core.sha(path) for path in extra});post.check_sources(bound)
    newstate=model.explicit_transition(state,PARENT_SHA)
    for key in state:
        if key!='classification':require(core.clean(newstate[key])==core.clean(state[key]),'transition reset history '+key)
    newstate['spatial_model_transition']['qualification']='PASS_PARENT_MAP_TANGENT_AND_COMPONENT_GATES'
    parent_screen=post.read(PARENT.parent.parent/'pair_016_temporal_screen.json')
    require(core.sha(PARENT.parent.parent/'pair_016_temporal_screen.json')==
        '38b653f34ea8b2d89d98990e847d09307ac7495046c91b7e6619ef3a74578563','inherited event seeds changed')
    seeds=set(parent_screen['cumulative_shoulder_ids'])|set(state['bore_observation']['crossing_ids'])
    refs={};coarse=[]
    if stage=='cost':
        q=qual['parent_equivalence']
        per_root=1.5*(206*max(q['new_candidate_wall_s'],q['old_candidate_wall_s'])+
                      28*max(q['new_tangent_wall_s'],q['old_tangent_wall_s']))+10.
    else:
        cost=post.read(BASE/'cost/manifest.json')
        require(cost['status']=='PASS_COST_ROOT_NOT_ACCEPTED' and cost['source_sha256']==bound and
            cost['parent_sha256']==PARENT_SHA,'matching measured root cost absent')
        per_root=cost['candidate_wall_s'];refs[BASE/'cost/manifest.json']=core.sha(BASE/'cost/manifest.json')
    if stage=='fine':
        meta=post.read(BASE/'coarse/manifest.json');coarse=post.read(BASE/'coarse/candidate_history.json')
        require(meta['status']=='QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED' and meta['source_sha256']==bound and
            meta['parent_sha256']==PARENT_SHA and meta['accepted_steps']==0 and len(coarse)==PAIRS and
            abs(meta['candidate_time_s']-END)<1e-15,'complete same-origin coarse window absent')
        require(not list((BASE/'coarse').glob('step_*/accepted.npz')),'coarse became restartable')
        previous=PARENT_SHA
        for index,row in enumerate(coarse,1):
            require(row['step']==index and row['parent_sha256']==previous and row['dt_s']==COARSE_DT,'coarse chain differs')
            for path,digest in row['artifacts'].items():
                actual=BASE/'coarse'/path;require(core.sha(actual)==digest,'coarse evidence drift');refs[actual]=digest
            saved=core.load(str(BASE/'coarse'/row['checkpoint']))['solver_data']
            require(saved['restart_allowed'] is False and saved['classification']==model.CLASSIFICATION,'coarse metadata differs')
            previous=row['sha256']
        for file in ('manifest.json','candidate_history.json'):
            refs[BASE/'coarse'/file]=core.sha(BASE/'coarse'/file)
        per_root=max(per_root,max(r['wall_seconds'] for r in coarse))
    forecast=30+1.5*(1 if stage=='cost' else PAIRS*(2 if stage=='fine' else 1))*per_root
    require(forecast<=450.,'registered stage cost forecast exceeds 450s')
    return dict(pa=pa,state=newstate,p=p,provenance=provenance,frozen=bound,refs=refs,per_root=per_root,
        forecast=forecast,seeds=seeds,prefix=list(state['interval_energy']['records']),coarse=coarse)


def check_context(c):
    post.check_sources(c['frozen'])
    require(core.parameters()==(c['p'],c['provenance']),'registered input changed')
    for path,digest in c['refs'].items():require(core.sha(path)==digest,'bound reference changed')


def save_graph(directory,trial,r):
    graph=r['visibility_graph']
    core.write(directory/'visible_graph.json',dict(diagnostics=graph['diagnostics'],fingerprints=graph['fingerprints'],
        parent_time_s=trial.state['time_s'],endpoint_time_s=trial.state['time_s']+r['dt_s'],
        graph_fingerprint=r['graph_fingerprint'],numerical_domain=r['visibility_numerical_domain'],accepted=False))
    np.savez_compressed(directory/'visible_graph.npz',particle_id=trial.pa.particle_id,
        raw_i=r['raw_i'],raw_j=r['raw_j'],raw_g=r['raw_g'],i=r['i'],j=r['j'],g=r['g'],
        kept_raw_edge_indices=graph['kept_raw_edge_indices'],blocked_raw_edge_indices=graph['blocked_raw_edge_indices'],
        position0=trial.x0,position1=r['position'],correction=r['correction'],volume=r['volume'])


def candidate(pa,state,c,dt,directory,deadline,per_root):
    started=time.perf_counter();trial=RecordedTrial(pa,state,c['p']);guess=batch.guess_for(trial,dt)
    attempt=directory/'broyden';attempt.mkdir();fallback=None
    try:r,metrics=solver.solve(trial,dt,guess,attempt,deadline,24)
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
        require(time.perf_counter()+per_root+15<=deadline,'insufficient budget for one same-parent Newton retry')
        fallback=repr(error);attempt=directory/'newton_retry';attempt.mkdir()
        try:r,metrics=core.solve(trial,dt,guess,attempt,deadline,24)
        except Exception:
            if hasattr(trial,'last_complete_result'):
                evidence=attempt/'LAST_FULL_EVALUATION_NOT_NECESSARILY_BEST_ITERATE';evidence.mkdir()
                solver.journal(evidence,trial,trial.last_complete_result,'last_FULL_UNCOMMITTED_candidate')
                save_graph(evidence,trial,trial.last_complete_result)
            raise
        metrics['controller']='ORIGINAL_MODIFIED_NEWTON_SAME_PARENT_VISIBLE_FULL_MAP'
    metrics['fallback_reason']=fallback
    core.journal(directory,trial,r,'SOLVED_UNCOMMITTED_candidate');save_graph(directory,trial,r)
    # Fresh complete residual and all original mechanical/thermal/prefix gates.
    nextpa,nextstate,fresh=model.accepted_copy(trial,r)
    require(np.array_equal(r['i'],fresh['i']) and np.array_equal(r['j'],fresh['j']) and
        r['graph_fingerprint']==fresh['graph_fingerprint'],'fresh graph changed at identical displacement')
    require(nextstate['interval_energy']['records'][:len(c['prefix'])]==c['prefix'],'complete old prefix was reset')
    target,evidence=event.save_candidate(directory,trial,nextpa,nextstate,fresh,c['frozen'])
    save_graph(directory,trial,fresh)
    scalar=model.VisibleSupportTrial(pa,state,c['p']);p=c['p']
    scalar.contact=ChainRuleDiscreteContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'],particle_radius=p['dx']/2,normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'],friction=p['friction'],base_z=-p['substrate_depth'])
    sr=scalar.evaluate(fresh['displacement'],dt)
    differences=batch.candidate_equivalence(dict(sr,energy_comparison=core.energy_diagnostics(scalar,sr)),
        dict(fresh,energy_comparison=core.energy_diagnostics(trial,fresh)))
    require(sr['graph_fingerprint']==fresh['graph_fingerprint'] and np.array_equal(sr['g'],fresh['g']) and
        np.array_equal(sr['thermal']['pair_conductance_W_K'],fresh['thermal']['pair_conductance_W_K']),
        'scalar map did not use same graph/heat')
    core.write(directory/'scalar_batch_equivalence.json',dict(status='PASS_SAME_FILTERED_COMPLETE_MAP',differences=differences))
    geometry=post.geometry_audit(trial,fresh);core.write(directory/'geometry_audit.json',geometry)
    require(not geometry['transient_contacts_missed_by_endpoints'] and
        set(geometry['new_non_bore_features'])<= {'SHOULDER'},'unqualified contact event')
    t=state['time_s'];pose0,pose1=core.pose(p,t),core.pose(p,t+dt)
    shell=audit_contact_shell_finite_enclosure(trial.contact,trial.x0,fresh['position'],pose0,pose1)
    core.write(directory/'continuous_shell_audit.json',shell)
    require(shell['status']==event.CERTIFIED,'continuous contact shell unresolved')
    before=time.perf_counter()
    path=audit_visible_support_path(trial.x0,fresh['position'],fresh['i'],fresh['j'],pose0,pose1,p,pa.particle_id)
    path.update(wall_seconds=time.perf_counter()-before,candidate_sha256=core.sha(target),graph_sha256=core.sha(directory/'visible_graph.npz'))
    core.write(directory/'visible_path_audit.json',path)
    require(path['status']==CERTIFIED,'retained support graph not certified over complete step')
    events=[]
    for fraction in (0.,.25,.5,.75,1.):
        endpoint=trial.contact.geometry(trial.x0+fraction*fresh['displacement'],core.pose(p,t+fraction*dt))
        for i,j in np.argwhere(endpoint['active']):
            if j in (0,1,6,7):events.append(dict(fraction=fraction,particle_id=int(pa.particle_id[i]),feature=trial.contact.features[j]))
    require(set(e['feature'] for e in events)<= {'SHOULDER'},'sampled nonshoulder event')
    core.write(directory/'sampled_nonbore_events.json',events)
    return dict(pa=nextpa,state=nextstate,trial=trial,result=fresh,metrics=metrics,directory=directory,
        checkpoint=target,evidence=evidence,events=events,wall_seconds=time.perf_counter()-started)


def candidate_row(item,number,direct_hash,directory):
    sd=item['directory'];state=item['state'];r=item['result']
    artifacts={path.relative_to(directory).as_posix():core.sha(path) for path in sd.iterdir() if path.is_file()}
    return dict(step=number,global_step=state['step'],time_s=state['time_s'],dt_s=r['dt_s'],
        parent_sha256=direct_hash,checkpoint=item['checkpoint'].relative_to(directory).as_posix(),sha256=core.sha(item['checkpoint']),
        artifacts=artifacts,wall_seconds=item['wall_seconds'],events=item['events'],
        blocked_particle_ids=sorted(set(r['visibility_graph']['blocked_edge_pid_pairs'].ravel().tolist())),
        graph_fingerprint=r['graph_fingerprint'],blocked_edge_count=r['visibility_graph']['diagnostics']['blocked_edge_count'],
        force_residual_N=r['force_residual_N'],**item['metrics'])


def publish(item,c,direct_hash,screen_file,coarse_hash):
    pa,state,trial,r=item['pa'],item['state'],item['trial'],item['result'];directory=item['directory']
    require(state['classification']==model.CLASSIFICATION and not state['spatial_model_transition']['diagnostic_only'] and
        state['spatial_model_transition']['qualification']=='PASS_PARENT_MAP_TANGENT_AND_COMPONENT_GATES',
        'model or diagnostic state cannot be published')
    require(post.read(directory/'visible_path_audit.json')['status']==CERTIFIED,'whole-path proof absent')
    screen=post.read(screen_file)
    require(screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'],'pair qualification absent')
    state.update(implementation='UNIFIED_VISIBLE_SUPPORT_REBUILT_ALL_OPERATORS_V1',restart_allowed=True,
        source_sha256=core.frozen_sources(),implementation_source_sha256=c['frozen'],controller_source_sha256=c['frozen'],
        parameters=c['p'],parameter_provenance=c['provenance'],continuation_parent_sha256=direct_hash,
        branch_origin_sha256=PARENT_SHA,window_temporal_screen_sha256=core.sha(screen_file),window_coarse_reference_sha256=coarse_hash)
    core.checkpoint(directory/'UNPUBLISHED_candidate.npz',pa,dict(state,restart_allowed=False))
    core.verify_snapshot(directory/'UNPUBLISHED_candidate.npz',pa,state)
    staged=core.load(str(directory/'UNPUBLISHED_candidate.npz'))['solver_data']
    for key in ('interval_energy','bore_observation','spatial_model_transition','implementation_source_sha256','continuation_parent_sha256'):
        require(core.clean(staged[key])==core.clean(state[key]),'staged metadata differs '+key)
    check_context(c)
    core.vtk(directory/'fields.vtk',pa,trial.p,trial.bulk,state['time_s'],model.CLASSIFICATION)
    core.write(directory/'tool_pose.json',core.pose(trial.p,state['time_s']))
    core.write(directory/'transport.json',state['bore_observation'])
    core.checkpoint(directory/'accepted.npz',pa,state);core.verify_snapshot(directory/'accepted.npz',pa,state)


def worker(stage):
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False);start_wall=time.perf_counter()
    history=[];rows=[];report=dict(status='RUNNING',stage=stage,accepted_steps=0,accepted_time_s=START,
        accepted_time_advanced_s=0.,classification=model.CLASSIFICATION,physical_accuracy='ENGINEERING_APPROXIMATE_NOT_VALIDATED')
    core.write(directory/'manifest.json',report)
    def stop():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,candidate_history=rows,
            published_checkpoints_requiring_readback=[p.relative_to(directory).as_posix() for p in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        c=preflight(stage);pa,state=c['pa'],c['state'];per_root=c['per_root'];seeds=set(c['seeds'])
        report.update(source_sha256=c['frozen'],provenance=c['provenance'],parent_sha256=PARENT_SHA,
            parent_checkpoint=PARENT.relative_to(ROOT).as_posix(),forecast_wall_s=c['forecast'],candidate_time_s=START,
            registered_end_s=END,thread_environment={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')})
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',history)
        core.write(directory/'candidate_history.json',rows)
        # This transfer snapshot is explicitly not a new accepted physics step.
        core.checkpoint(directory/'MODEL_TRANSITION_NOT_RESTARTABLE.npz',pa,dict(state,restart_allowed=False))
        direct_hash=PARENT_SHA
        for pair in range(1,2 if stage=='cost' else PAIRS+1):
            pending=[]
            for half in range(1,3 if stage=='fine' else 2):
                require(time.perf_counter()+1.5*per_root+15<=start_wall+525,'remaining bounded window budget exhausted')
                number=2*(pair-1)+half if stage=='fine' else pair;dt=FINE_DT if stage=='fine' else COARSE_DT
                sd=directory/f'step_{number:03d}';sd.mkdir()
                item=candidate(pa,state,c,dt,sd,start_wall+525,per_root);pending.append(item)
                row=candidate_row(item,number,direct_hash if len(pending)==1 else core.sha(pending[-2]['checkpoint']),directory)
                rows.append(row);core.write(directory/'candidate_history.json',rows)
                pa,state=item['pa'],item['state'];per_root=max(per_root,item['wall_seconds'])
                seeds.update(e['particle_id'] for e in item['events']);seeds.update(row['blocked_particle_ids'])
                seeds.update(state['bore_observation']['crossing_ids'])
                report.update(candidate_steps=len(rows),candidate_time_s=state['time_s']);core.write(directory/'manifest.json',report)
                print(json.dumps(dict(step=number,time_s=state['time_s'],wall_s=item['wall_seconds'],blocked=row['blocked_edge_count'],accepted=False)),flush=True)
            if stage=='fine':
                cr=c['coarse'][pair-1];raw=core.load(str(BASE/'coarse'/cr['checkpoint']));coarse,cstate=raw['arrays']['solid'],raw['solver_data']
                require(abs(cstate['time_s']-state['time_s'])<1e-15,'coarse/fine endpoint times differ')
                seeds.update(e['particle_id'] for e in cr['events']);seeds.update(cr['blocked_particle_ids'])
                seeds.update(cstate['bore_observation']['crossing_ids']);graphs=[]
                with np.load(BASE/'coarse'/f'step_{pair:03d}'/'visible_graph.npz') as saved:
                    for a,b in (('i','j'),('raw_i','raw_j')):graphs.append((saved['particle_id'].copy(),saved[a].copy(),saved[b].copy()))
                for item in pending:
                    for a,b in (('i','j'),('raw_i','raw_j')):graphs.append((item['pa'].particle_id,item['result'][a],item['result'][b]))
                selected=event.one_hop_ids(seeds,graphs);global_screen=post.compare(pa,coarse,c['p']);local_screen=event.local_compare(pa,coarse,selected)
                match=state['bore_observation']['crossing_ids']==cstate['bore_observation']['crossing_ids']
                screen_file=directory/f'pair_{pair:03d}_temporal_screen.json'
                core.write(screen_file,dict(pair=pair,time_s=state['time_s'],parent_sha256=PARENT_SHA,
                    coarse_checkpoint_sha256=cr['sha256'],fine_candidate_sha256=[core.sha(x['checkpoint']) for x in pending],
                    cumulative_event_and_blocked_ids=sorted(seeds),global_screen=global_screen,local_screen=local_screen,
                    crossing_ids_match=match,fine_crossing_ids=state['bore_observation']['crossing_ids'],
                    coarse_crossing_ids=cstate['bore_observation']['crossing_ids'],
                    endpoint_graphs_match=pending[-1]['result']['graph_fingerprint']==cr['graph_fingerprint']))
                require(global_screen['pass'] and local_screen['passed'] and match,'same-parent global/local pair screen failed')
                for item in pending:
                    publish(item,c,direct_hash,screen_file,cr['sha256']);path=item['directory']/'accepted.npz'
                    entry=dict(step=len(history)+1,global_step=item['state']['step'],time_s=item['state']['time_s'],dt_s=FINE_DT,
                        parent_sha256=direct_hash,checkpoint=path.relative_to(directory).as_posix(),sha256=core.sha(path),
                        temporal_screen_sha256=core.sha(screen_file),coarse_reference_sha256=cr['sha256'])
                    history.append(entry);direct_hash=entry['sha256'];core.write(directory/'accepted_history.json',history)
                    report.update(accepted_steps=len(history),accepted_time_s=entry['time_s'],accepted_time_advanced_s=entry['time_s']-START)
                    core.write(directory/'manifest.json',report)
            else:direct_hash=core.sha(pending[-1]['checkpoint'])
            check_context(c)
        report['status']={'cost':'PASS_COST_ROOT_NOT_ACCEPTED','coarse':'QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED',
            'fine':'PASS_VISIBLE_SUPPORT_WINDOW_TIME_SCREEN_NOT_PHYSICAL_VALIDATION'}[stage]
        report.update(candidate_wall_s=max(r['wall_seconds'] for r in rows),crossing_ids=state['bore_observation']['crossing_ids'],
            crossing_mass_kg=state['crossing_mass_kg'],interval_energy=state['interval_energy'])
        code=0
    except Exception as error:
        traceback.print_exc()
        published=[p.relative_to(directory).as_posix() for p in directory.glob('step_*/accepted.npz')]
        status='FAILED_PUBLISHED_CHECKPOINTS_REQUIRE_READBACK' if len(published)!=len(history) else (
            'FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE')
        report.update(status=status,error=repr(error),published_checkpoints_requiring_readback=published);code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-start_wall,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',required=True,choices=['cost','coarse','fine']);parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.stage))
    logs=ROOT/'logs/visible_support_20261004';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/args.stage).exists(),'refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker'];begin=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,external_wall_s=time.perf_counter()-begin,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

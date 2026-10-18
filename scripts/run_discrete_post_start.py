"""Bounded post-start accuracy and cost probe; preserves v3 physical equations."""
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
import run_discrete_contact_process as core

BASE=core.ROOT/'outputs/discrete_contact_process_20261004/v3'
CONTRACT=core.ROOT/'docs/DISCRETE_POST_START_CONTRACT_20261004.md'


def read(path):return json.loads(path.read_text())


def check_sources(mapping):
    for path,digest in mapping.items():
        if core.sha(core.ROOT/path)!=digest:raise ValueError('source drift: '+path)


def compare(pa,coarse,p):
    if not np.array_equal(pa.particle_id,coarse.particle_id):raise ValueError('comparison IDs differ')
    norm=np.linalg.norm
    result=dict(position_max_m=float(norm(np.column_stack((pa.x-coarse.x,pa.y-coarse.y,pa.z-coarse.z)),axis=1).max()),
        temperature_max_K=float(np.abs(pa.temperature-coarse.temperature).max()),
        stress_relative=float(norm(pa.S-coarse.S)/max(norm(pa.S),1e6*np.sqrt(len(pa.x)))),
        plastic_relative=float(norm(pa.ep-coarse.ep)/max(norm(pa.ep),1e-5*np.sqrt(len(pa.x)))))
    result['pass']=bool(result['position_max_m']<=.01*p['dx'] and result['temperature_max_K']<=.5 and
                        result['stress_relative']<=.1 and result['plastic_relative']<=.2)
    return result


def geometry_audit(trial,r):
    """Samples supplement the existing continuous solid sweep, not certify it."""
    t=trial.state['time_s'];dt=r['dt_s'];contact=trial.contact
    first=contact.geometry(trial.x0,core.pose(trial.p,t))
    last=contact.geometry(r['position'],core.pose(trial.p,t+dt))
    transient=[];counts=[]
    active=first['active']|last['active']
    for fraction in (.25,.5,.75):
        middle=contact.geometry(trial.x0+fraction*r['displacement'],core.pose(trial.p,t+fraction*dt))
        missed=middle['active']&~first['active']&~last['active']
        transient.extend(dict(fraction=fraction,particle_id=int(trial.pa.particle_id[i]),feature=contact.features[j]) for i,j in np.argwhere(missed))
        counts.append(middle['active'].sum(axis=0).tolist());active|=middle['active']
    neighbors=np.bincount(np.r_[r['i'],r['j']],minlength=len(trial.x0))
    result=dict(endpoint_counts=last['active'].sum(axis=0).tolist(),interior_counts=counts,
        features=list(contact.features),transient_contacts_missed_by_endpoints=transient,
        maximum_penetration_m=float(np.max(last['overlap_m'])),min_neighbor_count=int(neighbors.min()),
        maximum_correction_condition=float(np.linalg.cond(r['correction']).max()),
        new_non_bore_features=[contact.features[j] for j in (0,1,6,7) if active[:,j].any()],
        limitation='Three samples do not certify absence of all penalty-shell contact events.')
    return result


def worker(stage):
    directory=BASE/'post_start'/stage;directory.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter();history=[]
    def stop():
        core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history));os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    report=dict(classification=core.CLASSIFICATION,stage=stage,accepted_steps=0,accepted_time_s=None,
        production_authorized=False,status='RUNNING_POST_START_DIAGNOSTIC',physical_accuracy='NOT_EXPERIMENTALLY_VALIDATED')
    try:
        baseline=core.frozen_sources();p,provenance=core.parameters()
        cost=read(BASE/'cost/manifest.json');tests=read(BASE/'tests/manifest.json')
        if cost['status']!='PASS_MEASURED_COST_NOT_PHYSICAL_VALIDATION' or tests['status']!='PASS_SELECTED_COMPONENT_TESTS_NOT_FULL_REGRESSION':
            raise ValueError('component/cost gate absent')
        for meta in (cost,tests):
            if meta['source_sha256']!=baseline or meta['provenance']!=provenance:raise ValueError('component/cost source or input mismatch')
        quarters=read(BASE/'quarters/manifest.json')
        check_sources(quarters['source_sha256'])
        if quarters['status']!='PASS_SHORT_TIME_REFINEMENT_NOT_TRANSPORT' or not quarters['temporal_comparison']['pass'] or quarters['accepted_steps']!=4:
            raise ValueError('startup refinement has not passed')
        if quarters['provenance']!=provenance or abs(quarters['accepted_time_s']-.0005)>1e-15:raise ValueError('startup input/time mismatch')
        parent=BASE/'quarters';parent_history=read(parent/'accepted_history.json')
        per_root=max(r['wall_seconds'] for r in parent_history)*3.
        dts=[.0005] if stage=='macro' else [.00025,.00025]
        if stage in ('halves','pilot'):
            macro=read(BASE/'post_start/macro/manifest.json')
            if macro['status']!='PASS_POST_START_ROOT' or macro['baseline_source_sha256']!=baseline or macro['provenance']!=provenance:
                raise ValueError('matching post-start macro absent')
            check_sources(macro['source_sha256']);per_root=macro['wall_seconds']
        if stage=='pilot':
            parent=BASE/'post_start/halves';prior=read(parent/'manifest.json')
            if prior['status']!='PASS_POST_START_TIME_SCREEN' or not prior['temporal_comparison']['pass'] or prior['accepted_steps']!=2:
                raise ValueError('matching post-start temporal screen absent')
            if prior['baseline_source_sha256']!=baseline or prior['provenance']!=provenance or abs(prior['accepted_time_s']-.001)>1e-15:
                raise ValueError('post-start parent version/time mismatch')
            check_sources(prior['source_sha256']);parent_history=read(parent/'accepted_history.json')
            per_root=max(per_root,max(r['wall_seconds'] for r in parent_history));dts=[.0005]*8
        forecast=30+1.5*len(dts)*per_root
        if forecast>450:raise ValueError('measured forecast exceeds 450s')
        parent_row=parent_history[-1];parent_file=parent/parent_row['checkpoint'];parent_hash=core.sha(parent_file)
        if parent_hash!=parent_row['sha256']:raise ValueError('parent checkpoint hash mismatch')
        pa,state=core.restore(parent_file)
        expected_time=.001 if stage=='pilot' else .0005
        expected_step=6 if stage=='pilot' else 4
        expected_history_count=2 if stage=='pilot' else 4
        if len(parent_history)!=expected_history_count or state['step']!=expected_step or abs(state['time_s']-expected_time)>1e-15 or abs(parent_row['time_s']-state['time_s'])>1e-15:
            raise ValueError('parent step/time/history mismatch')
        if stage=='halves' and macro['parent_sha256']!=parent_hash:raise ValueError('refinement branches have different parents')
        if state['source_sha256']!=baseline or state['parameters']!=p or state['parameter_provenance']!=provenance:
            raise ValueError('checkpoint source/parameters mismatch')
        core.DiscreteContactTrial(pa,state,p) # validate all inherited histories
        extra={str(path.relative_to(core.ROOT)).replace('\\','/'):core.sha(path) for path in (Path(__file__).resolve(),CONTRACT)}
        frozen=dict(baseline,**extra)
        report.update(source_sha256=frozen,baseline_source_sha256=baseline,provenance=provenance,
            parent_checkpoint=str(parent_file.relative_to(core.ROOT)),parent_sha256=parent_hash,
            starting_time_s=state['time_s'],accepted_time_s=state['time_s'],forecast_wall_s=forecast,
            particle_count=len(pa.x),requested_dts_s=dts,
            numerical_thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',[])
        core.checkpoint(directory/'continuation_initial.npz',pa,state)
        initial_pa,initial_state=core.restore(BASE/'quarters/independent_initial.npz')
        core.checkpoint(directory/'independent_initial.npz',initial_pa,initial_state)
        for number,dt in enumerate(dts,1):
            if time.perf_counter()-begin+1.5*per_root+15>525:raise ValueError('remaining budget exhausted')
            stepdir=directory/f'step_{number:03d}';stepdir.mkdir()
            trial=core.DiscreteContactTrial(pa,state,p)
            guess=np.column_stack((pa.u,pa.v,pa.w))*dt;fixed=~state['free'];guess[fixed]=trial.prescribed(dt)[fixed]
            started=time.perf_counter();r,metrics=core.solve(trial,dt,guess,stepdir,begin+525,cost['maximum_iterations'])
            audit=geometry_audit(trial,r);core.write(stepdir/'geometry_audit.json',audit)
            if audit['transient_contacts_missed_by_endpoints']:raise ValueError('unresolved interior penalty-shell contact; needs event refinement')
            if audit['new_non_bore_features']:raise ValueError('new finite contact feature requires a separate local time screen')
            newpa,newstate,fresh=core.accepted_copy(trial,r)
            newstate.update(controller_source_sha256=extra,continuation_parent_sha256=parent_hash)
            check_sources(frozen)
            if core.sha(parent_file)!=parent_hash:raise ValueError('parent drift')
            core.publish(stepdir,newpa,newstate,trial,fresh,baseline,provenance)
            check_sources(frozen)
            _,saved=core.restore(stepdir/'accepted.npz')
            if saved['controller_source_sha256']!=extra or saved['continuation_parent_sha256']!=parent_hash:raise ValueError('controller provenance lost')
            pa,state=newpa,newstate
            row=dict(step=number,time_s=state['time_s'],dt_s=dt,checkpoint=f'step_{number:03d}/accepted.npz',**metrics,
                sha256=core.sha(stepdir/'accepted.npz'),wall_seconds=time.perf_counter()-started,
                force_residual_N=fresh['force_residual_N'],energy_raw_J=state['energy']['global_raw_residual_J'])
            history.append(row);core.write(directory/'accepted_history.json',history);print(json.dumps(core.clean(row)),flush=True)
            report.update(accepted_steps=len(history),accepted_time_s=state['time_s'],interval_energy=state['interval_energy'],
                crossing_ids=state['crossing_ids'],crossing_mass_kg=state['crossing_mass_kg'],max_plastic_strain=float(pa.ep.max()))
            core.write(directory/'manifest.json',report)
        report['status']='PASS_POST_START_ROOT' if stage=='macro' else 'PASS_BOUNDED_5MS_PILOT_NOT_TRANSPORT_VALIDATION'
        if stage=='halves':
            coarse,coarse_state=core.restore(BASE/'post_start/macro/step_001/accepted.npz')
            if abs(coarse_state['time_s']-state['time_s'])>1e-15:raise ValueError('comparison endpoint time differs')
            comparison=compare(pa,coarse,p);core.write(directory/'temporal_comparison.json',comparison)
            report['temporal_comparison']=comparison
            if not comparison['pass']:raise ValueError('post-start temporal screen failed')
            report['status']='PASS_POST_START_TIME_SCREEN'
        check_sources(frozen);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_FACTS_RETAINED' if history else 'REJECTED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['macro','halves','pilot'],required=True);parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.stage))
    logs=core.ROOT/'logs/discrete_contact_process_20261004/v3/post_start';logs.mkdir(parents=True,exist_ok=True)
    log=logs/(args.stage+'.log')
    if log.exists() or (BASE/'post_start'/args.stage).exists():raise ValueError('refuse evidence overwrite')
    cmd=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker']
    child_env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
    begin=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(cmd,cwd=core.ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,env=child_env).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=cmd,exit_code=code,external_wall_s=time.perf_counter()-begin,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

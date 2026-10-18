"""Equivalent batched contact: tests, same-root replay, bounded 5->6 ms screen."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import unittest
import numpy as np
import run_discrete_contact_process as core
import run_discrete_post_start as post
import secant_equilibrium as solver

BASE=core.ROOT/'outputs/coupled_cost_20261004/batched'
OLD=post.BASE
CONTRACT=core.ROOT/'docs/BATCHED_PROCESS_CONTRACT_20261004.md'


def sources():
    paths=[Path(__file__).resolve(),CONTRACT,core.ROOT/'scripts/batched_contact.py',
           core.ROOT/'scripts/secant_equilibrium.py',core.ROOT/'scripts/run_discrete_post_start.py',
           post.CONTRACT,core.ROOT/'tests/review/test_batched_contact.py',
           core.ROOT/'tests/review/test_secant_equilibrium.py']
    return dict(core.frozen_sources(),**{p.relative_to(core.ROOT).as_posix():core.sha(p) for p in paths})


def make_trial(pa,state,p):
    from batched_contact import BatchedChainRuleContact
    trial=core.DiscreteContactTrial(pa,state,p)
    trial.contact=BatchedChainRuleContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'],particle_radius=p['dx']/2,
        normal_stiffness=p['normal_stiffness'],tangential_stiffness=p['normal_stiffness'],
        friction=p['friction'],base_z=-p['substrate_depth'])
    return trial


def check_inputs():
    p,provenance=core.parameters();baseline=core.frozen_sources()
    meta=post.read(OLD/'cost/manifest.json')
    if meta['source_sha256']!=baseline or meta['provenance']!=provenance:raise ValueError('legacy source/input drift')
    return p,provenance,baseline


def parent_for(stage):
    if stage=='replay':return OLD/'post_start/pilot',6
    return OLD/'post_start/pilot_refinement',15


def guess_for(trial,dt):
    pa=trial.pa;u=np.column_stack((pa.u,pa.v,pa.w))*dt
    fixed=~trial.state['free'];u[fixed]=trial.prescribed(dt)[fixed];return u


def candidate_equivalence(old,new):
    rows={}
    for key,atol in (('force',1e-9),('stress',1e-5),('ep',1e-12),('temperature',1e-9),
                     ('density',1e-9),('gradient',1e-10),('position',1e-15)):
        a,b=np.asarray(old[key]),np.asarray(new[key])
        rows[key]=float(np.max(np.abs(a-b)))
        if not np.allclose(a,b,rtol=1e-10,atol=atol):raise ValueError('candidate mismatch '+key)
    for key,a in old['contact'].items():
        if key not in new['contact']:raise ValueError('missing contact field '+key)
        b=new['contact'][key]
        if isinstance(a,np.ndarray):
            if a.dtype.kind in 'biu':
                if not np.array_equal(a,b):raise ValueError('contact branch mismatch '+key)
            elif not np.allclose(a,b,rtol=1e-10,atol=1e-11,equal_nan=False):raise ValueError('contact field mismatch '+key)
    for key in ('plastic_work_J','global_raw_residual_J','absolute_work_scale_J'):
        # caller attaches same-function raw ledgers to bind equivalence to power
        a,b=old['energy_comparison'][key],new['energy_comparison'][key]
        if not np.isclose(a,b,rtol=1e-10,atol=1e-12):raise ValueError('raw energy mismatch '+key)
    return rows


def publish(directory,pa,state,trial,r,frozen,provenance,parent_hash,branch_hash=None):
    state.update(implementation='NUMPY_BATCHED_FINITE_PATCH_EQUIVALENT',
        implementation_source_sha256=frozen,controller_source_sha256=frozen,
        continuation_parent_sha256=parent_hash,branch_origin_sha256=branch_hash or parent_hash)
    state.update(source_sha256=core.frozen_sources(),parameters=trial.p,parameter_provenance=provenance)
    staged=dict(state,restart_allowed=False)
    core.checkpoint(directory/'UNPUBLISHED_candidate.npz',pa,staged)
    core.verify_snapshot(directory/'UNPUBLISHED_candidate.npz',pa,staged)
    saved=core.load(str(directory/'UNPUBLISHED_candidate.npz'))['solver_data']
    for key in ('implementation_source_sha256','continuation_parent_sha256','branch_origin_sha256','interval_energy','bore_observation'):
        if core.clean(saved[key])!=core.clean(staged[key]):raise ValueError('staged metadata differs: '+key)
    initial,_=core.restore(directory.parent/'independent_initial.npz')
    if not np.array_equal(initial.particle_id,pa.particle_id):raise ValueError('ID drift')
    core.vtk(directory/'fields.vtk',pa,trial.p,trial.bulk,state['time_s'],core.CLASSIFICATION)
    core.write(directory/'energy.json',dict(mechanical=state['energy'],interval=state['interval_energy'],thermal=r['thermal']['ledger']))
    core.write(directory/'tool_pose.json',core.pose(trial.p,state['time_s']))
    core.write(directory/'transport.json',state['bore_observation'])
    post.check_sources(frozen)
    if core.parameters()!=(trial.p,provenance):raise ValueError('input drift before commit')
    core.checkpoint(directory/'accepted.npz',pa,state)


def worker(stage):
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter();history=[]
    report=dict(stage=stage,status='RUNNING',classification=core.CLASSIFICATION,accepted_time_advanced_s=0.)
    def stop():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,published_checkpoints=[str(x.relative_to(directory)) for x in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        frozen=sources();p,provenance,baseline=check_inputs()
        report.update(source_sha256=frozen,baseline_source_sha256=baseline,provenance=provenance,
            thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        core.write(directory/'manifest.json',report)
        if stage=='tests':
            suite=unittest.TestSuite()
            for pattern in ('test_quasistatic*.py','test_discrete*.py','test_interval_acceptance.py','test_secant_equilibrium.py','test_batched_contact.py'):
                suite.addTests(unittest.defaultTestLoader.discover(str(core.ROOT/'tests/review'),pattern))
            result=unittest.TextTestRunner(verbosity=2,stream=sys.stdout).run(suite)
            report['tests_run']=result.testsRun
            if not result.wasSuccessful():raise ValueError('component gate failed')
            report['status']='PASS_SELECTED_COMPONENTS_NOT_FULL_REGRESSION'
        else:
            gate=post.read(BASE/'tests/manifest.json')
            if gate['status']!='PASS_SELECTED_COMPONENTS_NOT_FULL_REGRESSION' or gate['source_sha256']!=frozen:raise ValueError('matching tests absent')
            parent_dir,index=parent_for(stage);rows=post.read(parent_dir/'accepted_history.json')
            parent_meta=post.read(parent_dir/'manifest.json')
            expected_status='PASS_BOUNDED_5MS_PILOT_NOT_TRANSPORT_VALIDATION' if stage=='replay' else 'PASS_1_TO_5MS_TIME_SCREEN_NOT_TRANSPORT'
            expected_count=8 if stage=='replay' else 16
            if parent_meta['status']!=expected_status or len(rows)!=expected_count or parent_meta['accepted_steps']!=expected_count or abs(parent_meta['accepted_time_s']-.005)>1e-15:
                raise ValueError('parent branch incomplete or unqualified')
            if stage!='replay' and not parent_meta['temporal_comparison']['pass']:raise ValueError('parent time screen absent')
            if parent_meta['baseline_source_sha256']!=baseline or parent_meta['provenance']!=provenance:raise ValueError('parent manifest provenance differs')
            post.check_sources(parent_meta['source_sha256'])
            parent_file=parent_dir/rows[index]['checkpoint'];parent_hash=core.sha(parent_file)
            if parent_hash!=rows[index]['sha256']:raise ValueError('parent hash drift')
            pa,state=core.restore(parent_file)
            if state['source_sha256']!=baseline or state['parameters']!=p or state['parameter_provenance']!=provenance:raise ValueError('parent input mismatch')
            report.update(parent_checkpoint=parent_file.relative_to(core.ROOT).as_posix(),parent_sha256=parent_hash,
                starting_time_s=state['time_s'],accepted_time_s=state['time_s'],particle_count=len(pa.x))
            if stage=='replay':
                trial=make_trial(pa,state,p);old_trial=core.DiscreteContactTrial(pa,state,p);dt=.0005;guess=guess_for(trial,dt)
                t=time.perf_counter();old=old_trial.evaluate(guess,dt);old_wall=time.perf_counter()-t
                t=time.perf_counter();new=trial.evaluate(guess,dt);new_wall=time.perf_counter()-t
                old['energy_comparison']=core.energy_diagnostics(old_trial,old);new['energy_comparison']=core.energy_diagnostics(trial,new)
                differences=candidate_equivalence(old,new)
                t=time.perf_counter();core.newton_direction(trial,new);tangent_wall=time.perf_counter()-t
                forecast=30+1.5*(100*new_wall+3*tangent_wall)
                if forecast>180:raise ValueError('replay forecast exceeds 180s')
                started=time.perf_counter();r,metrics=solver.solve(trial,dt,guess,directory,begin+525,24)
                next_pa,next_state,fresh=core.accepted_copy(trial,r)
                audit=post.geometry_audit(trial,fresh)
                if audit['transient_contacts_missed_by_endpoints'] or audit['new_non_bore_features']:raise ValueError('replay geometry event')
                elapsed=time.perf_counter()-started
                reference_file=parent_dir/rows[index+1]['checkpoint']
                if core.sha(reference_file)!=rows[index+1]['sha256']:raise ValueError('reference drift')
                reference,_=core.restore(reference_file);comparison=post.compare(next_pa,reference,p)
                comparison['pass']=bool(comparison['position_max_m']<=1e-8 and comparison['temperature_max_K']<=1e-4 and comparison['stress_relative']<=1e-3 and comparison['plastic_relative']<=1e-3)
                next_state.update(restart_allowed=False,implementation_source_sha256=frozen)
                core.checkpoint(directory/'NOT_RESTARTABLE_replay.npz',next_pa,next_state)
                core.write(directory/'energy.json',dict(mechanical=next_state['energy'],interval=next_state['interval_energy'],thermal=fresh['thermal']['ledger']))
                report.update(candidate_differences=differences,original_candidate_wall_s=old_wall,candidate_wall_s=new_wall,
                    tangent_linear_wall_s=tangent_wall,root_wall_s=elapsed,solver_metrics=metrics,comparison=comparison,
                    reference_sha256=core.sha(reference_file),forecast_wall_s=forecast)
                if not comparison['pass']:raise ValueError('same-root screen failed')
                report['status']='PASS_EQUIVALENT_MAP_AND_SAME_ROOT_MEASURED_COST'
            else:
                if state['step']!=22 or abs(state['time_s']-.005)>1e-15:raise ValueError('expected fine 5ms parent')
                cost=post.read(BASE/'replay/manifest.json')
                if cost['status']!='PASS_EQUIVALENT_MAP_AND_SAME_ROOT_MEASURED_COST' or cost['source_sha256']!=frozen:raise ValueError('equivalence/cost gate absent')
                dts=[.001] if stage=='macro' else [.0005,.0005]
                per_root=3*cost['root_wall_s']
                if stage=='halves':
                    macro=post.read(BASE/'macro/manifest.json')
                    if macro['status']!='PASS_REGISTERED_1MS_ROOT' or macro['source_sha256']!=frozen or macro['parent_sha256']!=parent_hash:raise ValueError('matching macro parent absent')
                    per_root=macro['wall_seconds']
                forecast=30+1.5*len(dts)*per_root
                if forecast>450:raise ValueError('root sequence forecast exceeds 450s')
                report['forecast_wall_s']=forecast;core.write(directory/'manifest.json',report)
                initial,initial_state=core.restore(OLD/'quarters/independent_initial.npz')
                core.checkpoint(directory/'independent_initial.npz',initial,initial_state)
                core.checkpoint(directory/'continuation_initial.npz',pa,state)
                direct_parent_file=parent_file;direct_parent_hash=parent_hash
                for number,dt in enumerate(dts,1):
                    if time.perf_counter()-begin+1.5*per_root+15>525:raise ValueError('remaining budget exhausted')
                    stepdir=directory/f'step_{number:03d}';stepdir.mkdir();t=time.perf_counter()
                    trial=make_trial(pa,state,p);r,metrics=solver.solve(trial,dt,guess_for(trial,dt),stepdir,begin+525,24)
                    audit=post.geometry_audit(trial,r);core.write(stepdir/'geometry_audit.json',audit)
                    if audit['transient_contacts_missed_by_endpoints'] or audit['new_non_bore_features']:raise ValueError('new geometry event requires local validation')
                    newpa,newstate,fresh=core.accepted_copy(trial,r)
                    if core.sha(parent_file)!=parent_hash or core.sha(direct_parent_file)!=direct_parent_hash:raise ValueError('parent drift')
                    publish(stepdir,newpa,newstate,trial,fresh,frozen,provenance,direct_parent_hash,parent_hash)
                    pa,state=newpa,newstate
                    row=dict(step=number,time_s=state['time_s'],dt_s=dt,parent_sha256=direct_parent_hash,checkpoint=f'step_{number:03d}/accepted.npz',sha256=core.sha(stepdir/'accepted.npz'),wall_seconds=time.perf_counter()-t,**metrics)
                    history.append(row);core.write(directory/'accepted_history.json',history);print(json.dumps(row),flush=True)
                    direct_parent_file=stepdir/'accepted.npz';direct_parent_hash=row['sha256']
                    per_root=max(per_root,row['wall_seconds'])
                    report.update(accepted_steps=len(history),accepted_time_s=state['time_s'],accepted_time_advanced_s=state['time_s']-.005,
                        interval_energy=state['interval_energy'],crossing_mass_kg=state['crossing_mass_kg'],max_plastic_strain=float(pa.ep.max()))
                    core.write(directory/'manifest.json',report)
                report['status']='PASS_REGISTERED_1MS_ROOT'
                if stage=='halves':
                    macro_rows=post.read(BASE/'macro/accepted_history.json');macro_file=BASE/'macro/step_001/accepted.npz'
                    if len(macro_rows)!=1 or core.sha(macro_file)!=macro_rows[0]['sha256']:raise ValueError('macro endpoint hash drift')
                    coarse,coarse_state=core.restore(macro_file)
                    if abs(coarse_state['time_s']-.006)>1e-15 or abs(coarse_state['time_s']-state['time_s'])>1e-15 or coarse_state['branch_origin_sha256']!=parent_hash:
                        raise ValueError('comparison origin/time mismatch')
                    comparison=post.compare(pa,coarse,p)
                    core.write(directory/'temporal_comparison.json',comparison);report['temporal_comparison']=comparison
                    if not comparison['pass']:raise ValueError('5 to 6ms temporal screen failed')
                    report['status']='PASS_5_TO_6MS_TIME_SCREEN_NOT_TRANSPORT'
        post.check_sources(frozen);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['tests','replay','macro','halves'],required=True);parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.stage))
    logs=core.ROOT/'logs/coupled_cost_20261004/batched';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    if log.exists() or (BASE/args.stage).exists():raise ValueError('refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker'];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');t=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-t,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

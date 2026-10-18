"""Tests -> cost -> DG force root -> equal-time refinement, all supervised.

Each tag is an immutable evidence version; no prior result is overwritten.
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
import unittest
import numpy as np
from run_hot_quasistatic_probe import (ROOT, parameters, sources, write, sha, clean,
    summary, journal, verify_snapshot, memory_bytes, initial_audit, pose,
    checkpoint, restore, vtk, energy_diagnostics, newton_direction)
from afsd_backend.pysph.discrete_process import (
    CLASSIFICATION, DiscreteContactTrial, independent_initial_state, accepted_copy)
from pysph.solver.utils import load


def frozen_sources():
    extra = [Path(__file__).resolve(), ROOT/'docs/DISCRETE_CONTACT_PROCESS_CONTRACT_20261004.md',
             ROOT/'docs/CHAIN_CONTACT_ADDENDUM_20261004.md']
    extra += list((ROOT/'tests/review').glob('test_discrete*.py'))
    extra += list((ROOT/'tests/review').glob('test_interval_acceptance.py'))
    return dict(sources(), **{str(p.relative_to(ROOT)).replace('\\', '/'):sha(p) for p in extra})


def solve(trial, dt, guess, directory, deadline, maximum=24):
    r = trial.evaluate(guess, dt); count = 1; rows = []
    for iteration in range(maximum+1):
        row = dict(iteration=iteration, evaluations=count, **summary(r)); rows.append(row)
        write(directory/'iterations.json', rows)
        if r['scaled_residual'] <= 1:
            journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            return r, dict(iterations=iteration, evaluations=count)
        if iteration == maximum or time.perf_counter() > deadline:
            journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            raise ValueError('registered iteration/time budget exhausted')
        try:
            delta, linear = newton_direction(trial, r)
        except Exception:
            journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            raise
        merit = np.linalg.norm(r['force'][trial.state['free']]); found = False
        for backtrack in range(8):
            alpha = .5**backtrack
            try:
                next_r = trial.evaluate(r['displacement']+alpha*delta, dt); count += 1
                next_merit = float(np.linalg.norm(next_r['force'][trial.state['free']]))
                row.setdefault('line_search', []).append(dict(alpha=alpha, merit_N=next_merit))
                if next_merit < (1-1e-4*alpha)*merit:
                    r = next_r; found = True; break
            except (ValueError, ArithmeticError) as error:
                count += 1
                row.setdefault('line_search', []).append(dict(alpha=alpha, rejected=str(error)))
            if count >= 100 or time.perf_counter() > deadline: break
        row['linear'] = linear; write(directory/'iterations.json', rows)
        if not found or count >= 100:
            journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            raise ValueError('complete residual line search/evaluation budget failed')


def publish(directory, pa, state, trial, r, frozen, provenance):
    state.update(source_sha256=frozen, parameters=trial.p, parameter_provenance=provenance)
    staged = dict(state, restart_allowed=False)
    checkpoint(directory/'UNPUBLISHED_candidate.npz', pa, staged)
    verify_snapshot(directory/'UNPUBLISHED_candidate.npz', pa, staged)
    saved = load(str(directory/'UNPUBLISHED_candidate.npz'))['solver_data']
    for key in ('interval_energy', 'bore_observation'):
        if clean(saved[key]) != clean(staged[key]):
            raise ValueError('checkpoint new-model history mismatch: '+key)
    read_pa, _ = restore(directory.parent/'independent_initial.npz')
    if not np.array_equal(pa.particle_id, read_pa.particle_id): raise ValueError('ID drift')
    vtk(directory/'fields.vtk', pa, trial.p, trial.bulk, state['time_s'], CLASSIFICATION)
    write(directory/'energy.json', dict(mechanical=state['energy'], interval=state['interval_energy'],
                                      thermal=r['thermal']['ledger']))
    write(directory/'tool_pose.json', pose(trial.p,state['time_s']))
    write(directory/'transport.json', state['bore_observation'])
    if frozen_sources() != frozen or parameters() != (trial.p, provenance):
        raise RuntimeError('source/input drift before commit')
    checkpoint(directory/'accepted.npz', pa, state)


def worker(stage, base):
    directory = base/stage; directory.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter(); history = []
    def timeout():
        write(directory/'HARD_TIMEOUT.json', dict(accepted_history=history, restart_only_published=True))
        os._exit(124)
    timer = threading.Timer(540.,timeout); timer.daemon=True; timer.start()
    report = dict(classification=CLASSIFICATION, stage=stage, accepted_steps=0, accepted_time_s=0.,
        production_authorized=False, physical_accuracy='NOT_VALIDATED_BY_NUMERICAL_GATES',
        original_dynamic_time_s=1.3503936329109039e-6, unrelated_old_QS_time_s=8.28399658203125e-5,
        numerical_thread_environment={key:os.environ.get(key) for key in
            ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
    try:
        frozen = frozen_sources(); p, provenance = parameters()
        report.update(source_sha256=frozen, provenance=provenance, particle_count=342)
        write(directory/'manifest.json', dict(report,status='RUNNING'))
        if stage == 'tests':
            suite = unittest.TestSuite()
            for pattern in ('test_quasistatic*.py','test_discrete*.py','test_interval_acceptance.py'):
                suite.addTests(unittest.defaultTestLoader.discover(str(ROOT/'tests/review'), pattern))
            count = suite.countTestCases()
            result = unittest.TextTestRunner(verbosity=2,stream=sys.stdout).run(suite)
            if not result.wasSuccessful(): raise ValueError('selected component tests failed')
            report.update(status='PASS_SELECTED_COMPONENT_TESTS_NOT_FULL_REGRESSION', tests_run=count)
        else:
            gate = json.loads((base/'tests/manifest.json').read_text())
            if gate['status'] != 'PASS_SELECTED_COMPONENT_TESTS_NOT_FULL_REGRESSION' or gate['source_sha256'] != frozen:
                raise ValueError('matching tests have not passed')
            pa, state = independent_initial_state(p)
            checkpoint(directory/'independent_initial.npz',pa,state)
            write(directory/'accepted_history.json',[])
            write(directory/'resolved_config.json',dict(parameters=p,provenance=provenance))
            if stage == 'cost':
                trial = DiscreteContactTrial(pa,state,p)
                audit = initial_audit(trial); write(directory/'initial_audit.json',audit)
                t0=time.perf_counter(); r=trial.evaluate(trial.prescribed(.0005),.0005)
                trial_time=time.perf_counter()-t0
                journal(directory,trial,r,'UNCOMMITTED_cost_candidate')
                t0=time.perf_counter(); _, linear=newton_direction(trial,r)
                tangent_time=time.perf_counter()-t0
                maximum=min(24,int((450-30)/(1.5*(tangent_time+8*trial_time))))
                if maximum < 2: raise ValueError('measured cost cannot fit a useful root')
                report.update(status='PASS_MEASURED_COST_NOT_PHYSICAL_VALIDATION',
                    trial_wall_s=trial_time,tangent_linear_wall_s=tangent_time,linear=linear,
                    maximum_iterations=maximum,root_forecast_s=30+1.5*maximum*(tangent_time+8*trial_time),
                    stable_timestep_s='NOT_ESTABLISHED',requested_step_s=.0005,initial_audit=audit)
            else:
                cost=json.loads((base/'cost/manifest.json').read_text())
                if cost['source_sha256'] != frozen or cost['status'] != 'PASS_MEASURED_COST_NOT_PHYSICAL_VALIDATION':
                    raise ValueError('matching measured cost gate absent')
                if cost['provenance'] != provenance: raise ValueError('cost inputs differ')
                if stage == 'macro':
                    dts=[.0005]; forecast=cost['root_forecast_s']
                else:
                    macro=json.loads((base/'macro/manifest.json').read_text())
                    if macro['status']!='PASS_REGISTERED_ROOT' or macro['source_sha256'] != frozen:
                        raise ValueError('macro root must pass before halves')
                    dts=[.00025,.00025]; forecast=30+1.5*2*macro['wall_seconds']
                if forecast>450: raise ValueError('root sequence forecast exceeds 450 seconds')
                report['forecast_s']=forecast
                for number,dt in enumerate(dts,1):
                    stepdir=directory/f'step_{number:03d}'; stepdir.mkdir()
                    trial=DiscreteContactTrial(pa,state,p)
                    guess=np.column_stack((pa.u,pa.v,pa.w))*dt
                    fixed=~state['free']; guess[fixed]=trial.prescribed(dt)[fixed]
                    t0=time.perf_counter()
                    r,metrics=solve(trial,dt,guess,stepdir,started+525,cost['maximum_iterations'])
                    next_pa,next_state,fresh=accepted_copy(trial,r)
                    publish(stepdir,next_pa,next_state,trial,fresh,frozen,provenance)
                    pa,state=next_pa,next_state
                    entry=dict(step=number,time_s=state['time_s'],dt_s=dt,
                        checkpoint=f'step_{number:03d}/accepted.npz',**metrics)
                    history.append(entry); report.update(accepted_steps=len(history),accepted_time_s=state['time_s'])
                    entry.update(sha256=sha(stepdir/'accepted.npz'),wall_seconds=time.perf_counter()-t0,
                                 force_residual_N=fresh['force_residual_N'],
                                 energy_raw_J=state['energy']['global_raw_residual_J'])
                    write(directory/'accepted_history.json',history)
                    print(json.dumps(clean(entry)),flush=True)
                report.update(status='PASS_REGISTERED_ROOT',interval_energy=state['interval_energy'],
                    crossing_ids=state['crossing_ids'],crossing_mass_kg=state['crossing_mass_kg'],
                    max_plastic_strain=float(pa.ep.max()))
                if stage=='halves':
                    coarse,_=restore(base/'macro/step_001/accepted.npz')
                    dx=np.column_stack((pa.x-coarse.x,pa.y-coarse.y,pa.z-coarse.z))
                    norm=np.linalg.norm
                    comparison=dict(position_max_m=float(norm(dx,axis=1).max()),
                        temperature_max_K=float(np.max(np.abs(pa.temperature-coarse.temperature))),
                        stress_relative=float(norm(pa.S-coarse.S)/max(norm(pa.S),1e6*np.sqrt(len(pa.x)))),
                        plastic_relative=float(norm(pa.ep-coarse.ep)/max(norm(pa.ep),1e-5*np.sqrt(len(pa.x)))))
                    comparison['pass']=bool(comparison['position_max_m']<=.01*p['dx'] and
                        comparison['temperature_max_K']<=.5 and comparison['stress_relative']<=.1 and
                        comparison['plastic_relative']<=.2)
                    write(directory/'temporal_comparison.json',comparison)
                    report['temporal_comparison']=comparison
                    if not comparison['pass']: raise ValueError('registered equal-time refinement screen failed')
                    report['status']='PASS_SHORT_INTERVAL_REFINEMENT_NOT_PROCESS_VALIDATION'
        if frozen_sources()!=frozen: raise RuntimeError('source changed during stage')
        code=0
    except Exception as error:
        traceback.print_exc(); report.update(status='FAILED_WITH_ACCEPTED_FACTS_RETAINED' if history else 'REJECTED_NO_ACCEPTED_STATE',error=repr(error))
        code=2
    finally: timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
                  accepted_steps=len(history),accepted_time_s=history[-1]['time_s'] if history else 0.)
    write(directory/'manifest.json',report)
    print(json.dumps(clean({k:v for k,v in report.items() if k not in ('source_sha256','interval_energy')})),flush=True)
    return code


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['tests','cost','macro','halves'],required=True)
    parser.add_argument('--tag',default='v1');parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if not args.tag.isalnum(): raise ValueError('alphanumeric evidence tag required')
    base=ROOT/'outputs/discrete_contact_process_20261004'/args.tag
    if args.worker: raise SystemExit(worker(args.stage,base))
    logs=ROOT/'logs/discrete_contact_process_20261004'/args.tag;logs.mkdir(parents=True,exist_ok=True)
    log=logs/(args.stage+'.log')
    if log.exists() or (base/args.stage).exists(): raise ValueError('refuse evidence overwrite')
    cmd=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--tag',args.tag,'--worker']
    child_env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
    begin=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try: code=subprocess.run(cmd,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,env=child_env).returncode
        except subprocess.TimeoutExpired: code=124
    write(logs/(args.stage+'_supervision.json'),dict(command=cmd,exit_code=code,external_wall_s=time.perf_counter()-begin,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))))
    raise SystemExit(code)

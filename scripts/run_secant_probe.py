"""Frozen same-root solver replay; never advances the accepted process time."""
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

BASE=post.BASE/'secant_replay'
CONTRACT=core.ROOT/'docs/SECANT_REPLAY_CONTRACT_20261004.md'


def sources():
    extra=[Path(__file__).resolve(),CONTRACT,core.ROOT/'scripts/secant_equilibrium.py',
           core.ROOT/'tests/review/test_secant_equilibrium.py',core.ROOT/'scripts/run_discrete_post_start.py']
    return dict(core.frozen_sources(),**{str(p.relative_to(core.ROOT)).replace('\\','/'):core.sha(p) for p in extra})


def worker(stage):
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    def stop():
        core.write(directory/'HARD_TIMEOUT.json',dict(accepted_time_advanced_s=0));os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    report=dict(stage=stage,accepted_time_advanced_s=0.,status='RUNNING',classification=core.CLASSIFICATION)
    try:
        frozen=sources();report['source_sha256']=frozen
        core.write(directory/'manifest.json',report)
        if stage=='tests':
            suite=unittest.defaultTestLoader.discover(str(core.ROOT/'tests/review'),'test_secant_equilibrium.py')
            result=unittest.TextTestRunner(verbosity=2,stream=sys.stdout).run(suite)
            report['tests_run']=result.testsRun
            if not result.wasSuccessful():raise ValueError('solver fixture tests failed')
            report['status']='PASS_SOLVER_FIXTURES_NOT_AFSD_VALIDATION'
        else:
            import secant_equilibrium
            test=post.read(BASE/'tests/manifest.json')
            if test['status']!='PASS_SOLVER_FIXTURES_NOT_AFSD_VALIDATION' or test['source_sha256']!=frozen:raise ValueError('matching solver test gate absent')
            p,provenance=core.parameters();cost=post.read(post.BASE/'cost/manifest.json')
            if cost['source_sha256']!=core.frozen_sources() or cost['provenance']!=provenance:raise ValueError('physical cost source/input mismatch')
            forecast=30+1.5*(100*cost['trial_wall_s']+3*cost['tangent_linear_wall_s'])
            if forecast>180:raise ValueError('same-root replay forecast exceeds 180s')
            pilot=post.BASE/'post_start/pilot';manifest=post.read(pilot/'manifest.json');history=post.read(pilot/'accepted_history.json')
            if manifest['status']!='PASS_BOUNDED_5MS_PILOT_NOT_TRANSPORT_VALIDATION' or len(history)!=8:raise ValueError('pilot not complete')
            parent_file=pilot/history[-2]['checkpoint'];reference_file=pilot/history[-1]['checkpoint']
            for path,row in ((parent_file,history[-2]),(reference_file,history[-1])):
                if core.sha(path)!=row['sha256']:raise ValueError('replay input changed')
            parent,parent_state=core.restore(parent_file);reference,reference_state=core.restore(reference_file)
            if parent_state['source_sha256']!=core.frozen_sources() or parent_state['parameters']!=p:raise ValueError('parent model version mismatch')
            dt=history[-1]['dt_s']
            if abs(parent_state['time_s']+dt-reference_state['time_s'])>1e-15:raise ValueError('replay times differ')
            trial=core.DiscreteContactTrial(parent,parent_state,p)
            guess=np.column_stack((parent.u,parent.v,parent.w))*dt;fixed=~parent_state['free'];guess[fixed]=trial.prescribed(dt)[fixed]
            solve_started=time.perf_counter()
            result,metrics=secant_equilibrium.solve(trial,dt,guess,directory,started+525,24)
            next_pa,next_state,fresh=core.accepted_copy(trial,result)
            audit=post.geometry_audit(trial,fresh);core.write(directory/'geometry_audit.json',audit)
            if audit['transient_contacts_missed_by_endpoints'] or audit['new_non_bore_features']:raise ValueError('replay geometry event')
            solve_wall=time.perf_counter()-solve_started
            comparison=post.compare(next_pa,reference,p)
            comparison['pass']=bool(comparison['position_max_m']<=1e-8 and comparison['temperature_max_K']<=1e-4 and
                                    comparison['stress_relative']<=1e-3 and comparison['plastic_relative']<=1e-3)
            speedup=history[-1]['wall_seconds']/solve_wall
            next_state.update(restart_allowed=False,solver_controller_source_sha256=frozen,
                replay_parent_sha256=core.sha(parent_file),replay_reference_sha256=core.sha(reference_file),
                accepted_time_advanced_s=0.)
            core.checkpoint(directory/'NOT_RESTARTABLE_replay_candidate.npz',next_pa,next_state)
            core.verify_snapshot(directory/'NOT_RESTARTABLE_replay_candidate.npz',next_pa,next_state)
            core.write(directory/'energy.json',dict(mechanical=next_state['energy'],interval=next_state['interval_energy'],thermal=fresh['thermal']['ledger']))
            report.update(comparison=comparison,speedup=speedup,solve_wall_s=solve_wall,reference_root_wall_s=history[-1]['wall_seconds'],
                forecast_wall_s=forecast,solver_metrics=metrics,force_residual_N=fresh['force_residual_N'],
                parent_checkpoint=str(parent_file.relative_to(core.ROOT)),parent_sha256=core.sha(parent_file),
                reference_checkpoint=str(reference_file.relative_to(core.ROOT)),reference_sha256=core.sha(reference_file),
                numerical_thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')},
                status='PASS_SAME_ROOT_AND_COST_NOT_PROCESS_VALIDATION' if comparison['pass'] and speedup>=2 else 'FAIL_SAME_ROOT_OR_COST_SCREEN')
        if sources()!=frozen:raise ValueError('source drift')
        code=0 if not report['status'].startswith('FAIL') else 3
    except Exception as error:
        traceback.print_exc();report.update(status='REPLAY_OR_FIXTURE_FAILED_NO_TIME_ADVANCED',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['tests','replay'],required=True);parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.stage))
    logs=core.ROOT/'logs/discrete_contact_process_20261004/v3/secant_replay';logs.mkdir(parents=True,exist_ok=True)
    log=logs/(args.stage+'.log')
    if log.exists() or (BASE/args.stage).exists():raise ValueError('refuse evidence overwrite')
    cmd=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker']
    child_env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
    begin=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(cmd,cwd=core.ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,env=child_env).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=cmd,exit_code=code,external_wall_s=time.perf_counter()-begin,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

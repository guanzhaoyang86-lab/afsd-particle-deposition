"""Additional registered h/4 or h/8 diagnostic; no long-process authority."""
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


def worker(count):
    base=core.ROOT/'outputs/discrete_contact_process_20261004/v3'
    name='quarters' if count==4 else 'eighths';directory=base/name
    directory.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter();history=[]
    def stop():
        core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history));os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    report=dict(classification=core.CLASSIFICATION,stage=name,accepted_steps=0,accepted_time_s=0.,
        status='RUNNING_TIME_REFINEMENT_NOT_TRANSPORT',requested_time_s=.0005,substeps=count)
    try:
        baseline=core.frozen_sources();p,provenance=core.parameters()
        cost=json.loads((base/'cost/manifest.json').read_text())
        if baseline!=cost['source_sha256']:raise ValueError('cost source mismatch')
        previous_name='halves' if count==4 else 'quarters'
        previous=json.loads((base/previous_name/'manifest.json').read_text())
        if previous['accepted_time_s']!=.0005:raise ValueError('reference did not complete equal-time accepted sequence')
        per_root=previous['wall_seconds']/(count/2)
        forecast=1.5*per_root*count+30.
        if forecast>450.:raise ValueError('refinement sequence exceeds cost limit')
        extra=[Path(__file__).resolve(),core.ROOT/'docs/DISCRETE_REFINEMENT_EXTENSION_20261004.md']
        frozen=dict(baseline,**{str(p.relative_to(core.ROOT)).replace('\\','/'):core.sha(p) for p in extra})
        report.update(source_sha256=frozen,provenance=provenance,forecast_wall_s=forecast,
            numerical_thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',[])
        pa,state=core.independent_initial_state(p)
        core.checkpoint(directory/'independent_initial.npz',pa,state)
        for number in range(1,count+1):
            if time.perf_counter()-begin+1.5*per_root+15>525:raise ValueError('remaining cost budget exhausted')
            stepdir=directory/f'step_{number:03d}';stepdir.mkdir()
            dt=.0005/count;trial=core.DiscreteContactTrial(pa,state,p)
            guess=np.column_stack((pa.u,pa.v,pa.w))*dt
            fixed=~state['free'];guess[fixed]=trial.prescribed(dt)[fixed]
            started=time.perf_counter()
            r,metrics=core.solve(trial,dt,guess,stepdir,begin+525,cost['maximum_iterations'])
            newpa,newstate,fresh=core.accepted_copy(trial,r)
            # core publication checks baseline source version, supplemented by
            # the wrapper/contract checks before and after publication here.
            for path in extra:
                if core.sha(path)!=frozen[str(path.relative_to(core.ROOT)).replace('\\','/')]:raise RuntimeError('wrapper source drift')
            core.publish(stepdir,newpa,newstate,trial,fresh,baseline,provenance)
            pa,state=newpa,newstate
            row=dict(step=number,time_s=state['time_s'],dt_s=dt,checkpoint=f'step_{number:03d}/accepted.npz',**metrics)
            history.append(row)
            row.update(sha256=core.sha(stepdir/'accepted.npz'),wall_seconds=time.perf_counter()-started,
                       force_residual_N=fresh['force_residual_N'],energy_raw_J=state['energy']['global_raw_residual_J'])
            core.write(directory/'accepted_history.json',history);print(json.dumps(core.clean(row)),flush=True)
        coarse,_=core.restore(base/previous_name/f'step_{count//2:03d}/accepted.npz')
        norm=np.linalg.norm;delta=np.column_stack((pa.x-coarse.x,pa.y-coarse.y,pa.z-coarse.z))
        comparison=dict(position_max_m=float(norm(delta,axis=1).max()),
            temperature_max_K=float(np.abs(pa.temperature-coarse.temperature).max()),
            stress_relative=float(norm(pa.S-coarse.S)/max(norm(pa.S),1e6*np.sqrt(len(pa.x)))),
            plastic_relative=float(norm(pa.ep-coarse.ep)/max(norm(pa.ep),1e-5*np.sqrt(len(pa.x)))))
        comparison['pass']=bool(comparison['position_max_m']<=.01*p['dx'] and comparison['temperature_max_K']<=.5 and
                                 comparison['stress_relative']<=.1 and comparison['plastic_relative']<=.2)
        core.write(directory/'temporal_comparison.json',comparison)
        report.update(temporal_comparison=comparison,interval_energy=state['interval_energy'],
            status='PASS_SHORT_TIME_REFINEMENT_NOT_TRANSPORT' if comparison['pass'] else 'FAIL_SHORT_TIME_REFINEMENT',
            naive_0_13s_process_forecast_s=1.5*(time.perf_counter()-begin)/.0005*.13+30.)
        code=0 if comparison['pass'] else 3
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_FACTS_RETAINED' if history else 'REJECTED_NO_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(accepted_steps=len(history),accepted_time_s=history[-1]['time_s'] if history else 0.,
        wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes(),
        production_authorized=False)
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--substeps',choices=[4,8],type=int,required=True);parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.substeps))
    name='quarters' if args.substeps==4 else 'eighths'
    logs=core.ROOT/'logs/discrete_contact_process_20261004/v3';log=logs/(name+'.log')
    if log.exists():raise ValueError('refuse evidence overwrite')
    cmd=[sys.executable,'-B',str(Path(__file__).resolve()),'--substeps',str(args.substeps),'--worker']
    child_env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
    begin=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(cmd,cwd=core.ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,env=child_env).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(name+'_supervision.json'),dict(command=cmd,exit_code=code,external_wall_s=time.perf_counter()-begin,limit_s=600))
    print(json.dumps(dict(stage=name,exit_code=code,log=str(log))));raise SystemExit(code)

"""At most two additional equal-time screens, never an unbounded halving loop."""
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
import run_batched_window as window
batch,core,post=window.batch,window.core,window.post
BASE=window.BASE/'refinement'
CONTRACT=core.ROOT/'docs/WINDOW_REFINEMENT_CONTRACT_20261004.md'


def worker(level):
    name='quarter' if level==4 else 'eighth';directory=BASE/name;directory.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter();history=[];report=dict(status='RUNNING',accepted_steps=0,accepted_time_s=None,level=level)
    def stop():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,published_checkpoints=[str(p.relative_to(directory)) for p in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        context=window.preflight('coarse');pa,state=context['pa'],context['state'];p=context['p'];parent=context['parent'];parent_hash=context['parent_hash']
        previous=window.BASE/'fine' if level==4 else BASE/'quarter'
        meta=post.read(previous/'manifest.json');previous_rows=post.read(previous/'accepted_history.json')
        if abs(meta['accepted_time_s']-.020)>1e-15 or len(previous_rows)!=14*level//2 or meta['parent_sha256']!=parent_hash or meta['temporal_comparison']['pass']:
            raise ValueError('complete same-parent failed preceding time screen required')
        post.check_sources(meta['source_sha256'])
        previous_file=previous/previous_rows[-1]['checkpoint'];previous_hash=core.sha(previous_file)
        if previous_hash!=previous_rows[-1]['sha256']:raise ValueError('previous endpoint hash mismatch')
        per_root=max(row['wall_seconds'] for row in previous_rows);count=14*level;dt=.001/level
        forecast=30+1.5*count*per_root
        if forecast>450:raise ValueError('refinement forecast exceeds 450s')
        extra={path.relative_to(core.ROOT).as_posix():core.sha(path) for path in (Path(__file__).resolve(),CONTRACT)}
        frozen=dict(context['frozen'],**extra)
        report.update(source_sha256=frozen,provenance=context['provenance'],classification=core.CLASSIFICATION,
            parent_sha256=parent_hash,parent_checkpoint=parent.relative_to(core.ROOT).as_posix(),starting_time_s=state['time_s'],
            previous_endpoint=previous_file.relative_to(core.ROOT).as_posix(),previous_sha256=previous_hash,
            dt_s=dt,requested_time_s=.020,forecast_wall_s=forecast,particle_count=len(pa.x),
            thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',history)
        initial,initial_state=core.restore(batch.OLD/'quarters/independent_initial.npz')
        core.checkpoint(directory/'independent_initial.npz',initial,initial_state)
        core.checkpoint(directory/'continuation_initial.npz',pa,state)
        direct_file=parent;direct_hash=parent_hash
        for number in range(1,count+1):
            if time.perf_counter()-begin+1.5*per_root+15>525:raise ValueError('remaining budget exhausted')
            stepdir=directory/f'step_{number:03d}';stepdir.mkdir();started=time.perf_counter()
            trial=batch.make_trial(pa,state,p);guess=batch.guess_for(trial,dt);attempt=stepdir/'broyden';attempt.mkdir();fallback=None
            try:r,metrics=batch.solver.solve(trial,dt,guess,attempt,begin+525,24)
            except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
                if time.perf_counter()-begin+per_root+15>525:raise
                fallback=repr(error);attempt=stepdir/'newton_retry';attempt.mkdir()
                r,metrics=core.solve(trial,dt,guess,attempt,begin+525,24);metrics['controller']='ORIGINAL_MODIFIED_NEWTON_RETRY_SAME_PARENT'
            audit=post.geometry_audit(trial,r);core.write(stepdir/'geometry_audit.json',audit)
            if audit['transient_contacts_missed_by_endpoints'] or audit['new_non_bore_features']:
                core.journal(stepdir,trial,r,'geometry_event_UNCOMMITTED');raise ValueError('finite feature event requires separate validation')
            newpa,newstate,fresh=core.accepted_copy(trial,r)
            if core.sha(direct_file)!=direct_hash or core.sha(parent)!=parent_hash or core.sha(previous_file)!=previous_hash:raise ValueError('reference/parent drift')
            batch.publish(stepdir,newpa,newstate,trial,fresh,frozen,context['provenance'],direct_hash,parent_hash)
            pa,state=newpa,newstate
            row=dict(step=number,global_step=state['step'],time_s=state['time_s'],dt_s=dt,parent_sha256=direct_hash,
                checkpoint=f'step_{number:03d}/accepted.npz',sha256=core.sha(stepdir/'accepted.npz'),
                wall_seconds=time.perf_counter()-started,fallback_reason=fallback,force_residual_N=fresh['force_residual_N'],**metrics)
            history.append(row);core.write(directory/'accepted_history.json',history)
            direct_file=stepdir/'accepted.npz';direct_hash=row['sha256'];per_root=max(per_root,row['wall_seconds'])
            report.update(accepted_steps=len(history),accepted_time_s=state['time_s'],interval_energy=state['interval_energy'],
                crossing_mass_kg=state['crossing_mass_kg'],max_plastic_strain=float(pa.ep.max()))
            core.write(directory/'manifest.json',report);print(json.dumps(row),flush=True)
        other,other_state=core.restore(previous_file)
        if abs(other_state['time_s']-state['time_s'])>1e-15:raise ValueError('comparison end time mismatch')
        comparison=post.compare(pa,other,p);core.write(directory/'temporal_comparison.json',comparison)
        report.update(temporal_comparison=comparison,status='PASS_WINDOW_REFINEMENT_NOT_TRANSPORT' if comparison['pass'] else 'FAIL_WINDOW_REFINEMENT')
        post.check_sources(frozen);code=0 if comparison['pass'] else 3
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--level',choices=[4,8],type=int,required=True);parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.level))
    name='quarter' if args.level==4 else 'eighth';logs=core.ROOT/'logs/coupled_cost_20261004/batched/window_6_20ms/refinement';logs.mkdir(parents=True,exist_ok=True);log=logs/(name+'.log')
    if log.exists() or (BASE/name).exists():raise ValueError('refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--level',str(args.level),'--worker'];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');t=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(name+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-t,limit_s=600))
    print(json.dumps(dict(stage=name,exit_code=code,log=str(log))));raise SystemExit(code)

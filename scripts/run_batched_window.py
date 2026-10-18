"""One bounded, same-parent 6->20ms process window and temporal comparison."""
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
import run_batched_process as batch
core,post=batch.core,batch.post
BASE=batch.BASE/'window_6_20ms'
CONTRACT=core.ROOT/'docs/BATCHED_WINDOW_CONTRACT_20261004.md'


def preflight(stage):
    p,provenance,baseline=batch.check_inputs();old_frozen=batch.sources()
    gate=post.read(batch.BASE/'halves/manifest.json')
    if gate['status']!='PASS_5_TO_6MS_TIME_SCREEN_NOT_TRANSPORT' or gate['source_sha256']!=old_frozen or not gate['temporal_comparison']['pass']:
        raise ValueError('6ms matching time screen absent')
    history=post.read(batch.BASE/'halves/accepted_history.json')
    if len(history)!=2:raise ValueError('parent history size mismatch')
    parent=batch.BASE/'halves'/history[-1]['checkpoint'];parent_hash=core.sha(parent)
    if parent_hash!=history[-1]['sha256']:raise ValueError('parent hash mismatch')
    pa,state=core.restore(parent)
    if state['step']!=24 or abs(state['time_s']-.006)>1e-15 or state['source_sha256']!=baseline or state['implementation_source_sha256']!=old_frozen or state['parameters']!=p:
        raise ValueError('parent history/model mismatch')
    post.check_sources(old_frozen);batch.make_trial(pa,state,p)
    macro=post.read(batch.BASE/'macro/accepted_history.json')
    per_root=2*max([r['wall_seconds'] for r in history+macro]);count,dt=14,.001
    coarse_file=None;coarse_hash=None
    if stage=='fine':
        meta=post.read(BASE/'coarse/manifest.json');coarse_rows=post.read(BASE/'coarse/accepted_history.json')
        if meta['status']!='PASS_BOUNDED_WINDOW_NOT_TRANSPORT_VALIDATION' or meta['parent_sha256']!=parent_hash or len(coarse_rows)!=14:
            raise ValueError('complete same-parent coarse window absent')
        post.check_sources(meta['source_sha256'])
        per_root=max(per_root,max(r['wall_seconds'] for r in coarse_rows));count,dt=28,.0005
        coarse_file=BASE/'coarse'/coarse_rows[-1]['checkpoint'];coarse_hash=core.sha(coarse_file)
        if coarse_hash!=coarse_rows[-1]['sha256']:raise ValueError('coarse endpoint hash mismatch')
    forecast=30+1.5*count*per_root
    if forecast>450:raise ValueError('registered window forecast exceeds 450s')
    extra={path.relative_to(core.ROOT).as_posix():core.sha(path) for path in (Path(__file__).resolve(),CONTRACT)}
    frozen=dict(old_frozen,**extra)
    return dict(p=p,provenance=provenance,pa=pa,state=state,parent=parent,parent_hash=parent_hash,
        per_root=per_root,count=count,dt=dt,forecast=forecast,frozen=frozen,
        coarse_file=coarse_file,coarse_hash=coarse_hash)


def worker(stage):
    directory=BASE/stage;directory.mkdir(parents=True,exist_ok=False);begin=time.perf_counter();history=[]
    report=dict(stage=stage,status='RUNNING',accepted_steps=0,accepted_time_s=None)
    def stop():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(accepted_history=history,published_checkpoints=[str(x.relative_to(directory)) for x in directory.glob('step_*/accepted.npz')]))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        context=preflight(stage);p=context['p'];pa,state=context['pa'],context['state'];dt=context['dt'];frozen=context['frozen']
        per_root=context['per_root'];parent_file=context['parent'];parent_hash=context['parent_hash']
        report.update(classification=core.CLASSIFICATION,source_sha256=frozen,provenance=context['provenance'],
            parent_checkpoint=parent_file.relative_to(core.ROOT).as_posix(),parent_sha256=parent_hash,
            starting_time_s=state['time_s'],accepted_time_s=state['time_s'],requested_time_s=.020,
            particle_count=len(pa.x),forecast_wall_s=context['forecast'],dt_s=dt,
            thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',history)
        initial,initial_state=core.restore(batch.OLD/'quarters/independent_initial.npz')
        core.checkpoint(directory/'independent_initial.npz',initial,initial_state)
        core.checkpoint(directory/'continuation_initial.npz',pa,state)
        direct_file=parent_file;direct_hash=parent_hash
        for number in range(1,context['count']+1):
            if time.perf_counter()-begin+1.5*per_root+15>525:raise ValueError('remaining budget exhausted')
            stepdir=directory/f'step_{number:03d}';stepdir.mkdir();started=time.perf_counter()
            trial=batch.make_trial(pa,state,p);guess=batch.guess_for(trial,dt)
            first=stepdir/'broyden';first.mkdir();fallback=None
            try:r,metrics=batch.solver.solve(trial,dt,guess,first,begin+525,24)
            except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
                if time.perf_counter()-begin+per_root+15>525:raise
                fallback=repr(error);second=stepdir/'newton_retry';second.mkdir()
                r,metrics=core.solve(trial,dt,guess,second,begin+525,24)
                metrics['controller']='ORIGINAL_MODIFIED_NEWTON_RETRY_SAME_PARENT'
            audit=post.geometry_audit(trial,r);core.write(stepdir/'geometry_audit.json',audit)
            if audit['transient_contacts_missed_by_endpoints'] or audit['new_non_bore_features']:
                core.journal(stepdir,trial,r,'geometry_event_UNCOMMITTED');raise ValueError('first new finite contact event requires local validation')
            newpa,newstate,fresh=core.accepted_copy(trial,r)
            if core.sha(direct_file)!=direct_hash or core.sha(parent_file)!=parent_hash:raise ValueError('parent drift')
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
        if abs(state['time_s']-.020)>1e-15:raise ValueError('end time differs')
        report['status']='PASS_BOUNDED_WINDOW_NOT_TRANSPORT_VALIDATION'
        if stage=='fine':
            if core.sha(context['coarse_file'])!=context['coarse_hash']:raise ValueError('coarse endpoint drift')
            coarse,coarse_state=core.restore(context['coarse_file'])
            if abs(coarse_state['time_s']-state['time_s'])>1e-15 or coarse_state['branch_origin_sha256']!=parent_hash:raise ValueError('comparison time/origin mismatch')
            comparison=post.compare(pa,coarse,p);core.write(directory/'temporal_comparison.json',comparison);report['temporal_comparison']=comparison
            if not comparison['pass']:raise ValueError('6 to 20ms window time screen failed')
            report['status']='PASS_WINDOW_TIME_SCREEN_NOT_TRANSPORT_VALIDATION'
        post.check_sources(frozen);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['coarse','fine'],required=True);parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.stage))
    logs=core.ROOT/'logs/coupled_cost_20261004/batched/window_6_20ms';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    if log.exists() or (BASE/args.stage).exists():raise ValueError('refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker'];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-started,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

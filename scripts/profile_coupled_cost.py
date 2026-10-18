"""Bounded profiling of an already accepted root; advances no physical time."""
from pathlib import Path
import cProfile
import io
import json
import os
import pstats
import subprocess
import sys
import threading
import time
import numpy as np
import run_discrete_contact_process as core

OUT=core.ROOT/'outputs/coupled_cost_20261004/profile'


def worker():
    OUT.mkdir(parents=True,exist_ok=False)
    timer=threading.Timer(110.,lambda:os._exit(124));timer.daemon=True;timer.start()
    base=core.ROOT/'outputs/discrete_contact_process_20261004/v3/post_start/pilot_refinement'
    history=json.loads((base/'accepted_history.json').read_text())
    frozen=core.frozen_sources();p,provenance=core.parameters()
    files=[base/row['checkpoint'] for row in history[-2:]]
    for path,row in zip(files,history[-2:]):
        if core.sha(path)!=row['sha256']:raise ValueError('input hash drift')
    pa,state=core.restore(files[0]);end,end_state=core.restore(files[1])
    if state['source_sha256']!=frozen or state['parameters']!=p:raise ValueError('source mismatch')
    trial=core.DiscreteContactTrial(pa,state,p)
    u=np.column_stack((end.x-pa.x,end.y-pa.y,end.z-pa.z));dt=history[-1]['dt_s']
    u[~state['free']]=trial.prescribed(dt)[~state['free']]
    start=time.perf_counter();r=trial.evaluate(u,dt);eval_wall=time.perf_counter()-start
    if r['scaled_residual']>1:raise ValueError('saved accepted displacement did not replay')
    start=time.perf_counter();core.newton_direction(trial,r);tangent_wall=time.perf_counter()-start
    profiles={}
    for name,func in (('residual',lambda:trial.evaluate(u,dt)),('tangent',lambda:core.newton_direction(trial,r))):
        profiler=cProfile.Profile();profiler.enable();func();profiler.disable()
        profiler.dump_stats(str(OUT/(name+'.prof')))
        stream=io.StringIO();stats=pstats.Stats(profiler,stream=stream).strip_dirs().sort_stats('cumulative')
        stats.print_stats(45);(OUT/(name+'.txt')).write_text(stream.getvalue(),encoding='utf-8')
        profiles[name]=[dict(file=file,line=line,function=function,calls=value[1],self_s=value[2],cumulative_s=value[3])
            for (file,line,function),value in sorted(profiler.getstats() and pstats.Stats(profiler).stats.items(),key=lambda kv:kv[1][3],reverse=True)[:45]]
    if frozen!=core.frozen_sources():raise ValueError('source changed during profiling')
    report=dict(status='PROFILE_ONLY_NO_TIME_ADVANCE',unprofiled_residual_wall_s=eval_wall,
        unprofiled_tangent_linear_wall_s=tangent_wall,source_sha256=frozen,provenance=provenance,
        controller_sha256=core.sha(Path(__file__)),inputs={str(path.relative_to(core.ROOT)):core.sha(path) for path in files},
        force_residual_N=r['force_residual_N'],profiles=profiles,
        thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
    core.write(OUT/'manifest.json',report);timer.cancel()
    print(json.dumps(dict(residual_s=eval_wall,tangent_s=tangent_wall)));return 0


if __name__=='__main__':
    if '--worker' in sys.argv:raise SystemExit(worker())
    logs=core.ROOT/'logs/coupled_cost_20261004';logs.mkdir(parents=True,exist_ok=True)
    log=logs/'profile.log'
    if log.exists() or OUT.exists():raise ValueError('refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
    start=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=120).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/'profile_receipt.json',dict(command=command,exit_code=code,wall_s=time.perf_counter()-start,limit_s=120))
    print(log.read_text(encoding='utf-8'));raise SystemExit(code)

"""Bounded nonrestartable dt replay against actual smaller-step TOP trajectory."""
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
import run_visible_top_event as top

core,post,event,old=top.core,top.post,top.event,top.old
ROOT=core.ROOT
BASE=ROOT/'outputs/visible_top_timestep_20261004'
LOGS=ROOT/'logs/visible_top_timestep_20261004'
CONTRACT=ROOT/'docs/VISIBLE_TOP_TIMESTEP_SCREEN_CONTRACT_20261004.md'
REFERENCE=top.BASE/'window_001'
require=core.require if hasattr(core,'require') else top.require


def relative(path):return Path(path).relative_to(ROOT).as_posix()


def preflight(multiplier):
    require(multiplier in (4,8,16),'unregistered timestep')
    c=top.preflight(1,'coarse')
    auditpath=REFERENCE/'readback/summary.json';audit=post.read(auditpath)
    require(audit['status']=='PASS_VISIBLE_TOP_EVENT_READBACK_NOT_PHYSICAL_VALIDATION' and
        audit['window_completed'] and audit['top_event_process_screen_passed'] and
        audit['added_accepted_steps']==32 and audit['parent_sha256']==c['parent_hash'],
        'complete independent smaller-step reference absent')
    post.check_sources(audit['evidence_sha256'])
    manifest=post.read(REFERENCE/'fine/manifest.json')
    history=post.read(REFERENCE/'fine/accepted_history.json')
    rows=post.read(REFERENCE/'fine/candidate_history.json')
    require(manifest['status']==top.PASS_FINE and len(rows)==len(history)==32,'reference is incomplete')
    previous=c['parent_hash']
    for number,(entry,row) in enumerate(zip(history,rows),1):
        path=REFERENCE/'fine'/entry['checkpoint']
        require(entry['step']==number and entry['parent_sha256']==previous and
            entry['dt_s']==top.HALF and core.sha(path)==entry['sha256'],'reference accepted chain drift')
        saved=core.load(str(path))['solver_data']
        require(saved['restart_allowed'] is True and saved['interval_energy']['records'][:c['start_step']]==c['prefix'],
            'reference accepted history differs')
        previous=entry['sha256']
        c['refs'][path]=entry['sha256']
        for name,value in row['artifacts'].items():c['refs'][REFERENCE/'fine'/name]=value
    for path in (auditpath,REFERENCE/'fine/manifest.json',REFERENCE/'fine/accepted_history.json',
                 REFERENCE/'fine/candidate_history.json'):
        c['refs'][path]=core.sha(path)
    c['frozen'].update({relative(path):core.sha(path) for path in (Path(__file__),CONTRACT,auditpath)})
    count=32//multiplier
    fullcost=float(manifest['wall_seconds'])/32
    maxroot=max(r['wall_seconds'] for r in rows)
    predicted=max(fullcost,maxroot)*multiplier
    forecast=30+1.5*count*predicted
    require(forecast<=450,'larger-dt complete stage forecast exceeds 450s')
    c.update(reference_history=history,reference_rows=rows,audit_reference=auditpath,
        count=count,dt=top.HALF*multiplier,forecast=forecast,per_root=predicted,per_step=predicted)
    top.check_context(c)
    return c


def worker(multiplier):
    folder=BASE/f'mult_{multiplier:02d}';folder.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter();rows=[];screens=[]
    report=dict(status='RUNNING',classification=top.model.CLASSIFICATION,diagnostic_only=True,
        restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.,multiplier=multiplier)
    core.write(folder/'manifest.json',report)
    def stop():
        try:core.write(folder/'HARD_TIMEOUT.json',dict(candidate_history=rows,screens=screens))
        finally:os._exit(124)
    timer=threading.Timer(540,stop);timer.daemon=True;timer.start()
    try:
        c=preflight(multiplier);pa,state=c['pa'],c['state'];seeds=set(c['seeds']);previous=c['parent_hash']
        report.update(source_sha256=c['frozen'],parent_checkpoint=relative(c['parent']),parent_sha256=previous,
            start_time_s=c['start'],start_global_step=c['start_step'],dt_s=c['dt'],candidate_count=c['count'],
            reference_readback=relative(c['audit_reference']),reference_readback_sha256=core.sha(c['audit_reference']),
            forecast_wall_s=c['forecast'],reference_sha256={relative(p):v for p,v in c['refs'].items()},
            supervision_receipt=relative(LOGS/f'mult_{multiplier:02d}_supervision.json'))
        core.write(folder/'manifest.json',report)
        for number in range(1,c['count']+1):
            require(time.perf_counter()+1.5*c['per_step']+20<begin+525,'remaining stage budget exhausted')
            directory=folder/f'step_{number:03d}';directory.mkdir()
            item=top.qualified.candidate(pa,state,c,c['dt'],directory,begin+525,c['per_root'])
            row=old.candidate_row(item,number,previous,folder);rows.append(row)
            core.write(folder/'candidate_history.json',rows)
            pa,state=item['pa'],item['state'];previous=row['sha256']
            seeds.update(e['particle_id'] for e in row['events']);seeds.update(row['blocked_particle_ids'])
            seeds.update(state['bore_observation']['crossing_ids'])
            graphs=[]
            for a,b in (('i','j'),('raw_i','raw_j')):
                graphs.append((pa.particle_id,item['result'][a],item['result'][b]))
            first=(number-1)*multiplier;last=number*multiplier
            for refrow in c['reference_rows'][first:last]:
                seeds.update(e['particle_id'] for e in refrow['events']);seeds.update(refrow['blocked_particle_ids'])
                graphpath=REFERENCE/'fine'/f"step_{refrow['step']:03d}"/'visible_graph.npz'
                with np.load(graphpath,allow_pickle=False) as graph:
                    for a,b in (('i','j'),('raw_i','raw_j')):graphs.append((graph['particle_id'],graph[a],graph[b]))
            reference=c['reference_history'][last-1];refpath=REFERENCE/'fine'/reference['checkpoint']
            saved=core.load(str(refpath));refpa,rs=saved['arrays']['solid'],saved['solver_data']
            require(abs(rs['time_s']-state['time_s'])<1e-14,'comparison time mismatch')
            seeds.update(rs['bore_observation']['crossing_ids'])
            selected=event.one_hop_ids(seeds,graphs)
            # The established fine argument supplies stress/ep denominators.
            glob=post.compare(refpa,pa,c['p']);local=event.local_compare(refpa,pa,selected)
            match=state['bore_observation']['crossing_ids']==rs['bore_observation']['crossing_ids']
            screen=dict(step=number,time_s=state['time_s'],candidate_sha256=row['sha256'],
                reference_checkpoint=relative(refpath),reference_sha256=reference['sha256'],
                reference_fine_steps=list(range(first+1,last+1)),cumulative_event_and_blocked_ids=sorted(seeds),
                global_screen=glob,local_screen=local,crossing_ids_match=match)
            screens.append(screen);core.write(directory/'temporal_screen.json',screen)
            require(glob['pass'] and local['passed'] and match,'larger physical timestep disagrees with reference')
            top.check_context(c)
            report.update(completed_candidates=len(rows),candidate_end_s=state['time_s'])
            core.write(folder/'manifest.json',report)
            print(json.dumps(dict(multiplier=multiplier,step=number,time_s=state['time_s'],global_screen=glob,
                local_pass=local['passed'],accepted=False)),flush=True)
        topids=sorted({e['particle_id'] for r in rows for e in r['events'] if e['feature']=='TOP'})
        require(topids,'TOP event was not exercised')
        report.update(status='PASS_LARGER_TIMESTEP_REFERENCE_SCREEN_NOT_RESTARTABLE',top_particle_ids=topids,
            screens=screens,interval_energy=state['interval_energy'],source_unchanged=True,
            comparisons_are_against_62_5us=True,contact_history_convergence='NOT_ESTABLISHED')
        code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_NONRESTARTABLE_TIMESTEP_SCREEN',error=repr(error));code=2
    finally:timer.cancel()
    report.update(exit_code=code,wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(folder/'manifest.json',report)
    return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--multiplier',type=int,choices=(4,8,16),required=True)
    parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    if args.worker:raise SystemExit(worker(args.multiplier))
    LOGS.mkdir(parents=True,exist_ok=True);log=LOGS/f'mult_{args.multiplier:02d}.log'
    require(not log.exists() and not (BASE/f'mult_{args.multiplier:02d}').exists(),'refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--multiplier',str(args.multiplier),'--worker']
    begin=time.perf_counter()
    with log.open('w',encoding='utf-8') as out:
        try:code=subprocess.run(command,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(LOGS/f'mult_{args.multiplier:02d}_supervision.json',dict(command=command,exit_code=code,
        external_wall_s=time.perf_counter()-begin,limit_s=600))
    print(json.dumps(dict(multiplier=args.multiplier,exit_code=code,log=str(log))));raise SystemExit(code)

"""Two predeclared refinements from one real parent; diagnostic, never restartable."""
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
import run_expanded_hot_campaign_window as campaign

ROOT=campaign.ROOT
core,post,event,old=campaign.core,campaign.post,campaign.event,campaign.old
qualified=campaign.qualified
require=campaign.require
rel=campaign.relative
BASELINE=ROOT/'outputs/expanded_hot_campaign_20261005/window_001'
AUDIT=BASELINE/'readback/summary.json'
AUDIT_SHA='89f999f398a9f3de70535500121bbba3a359ab3b39cbd351a0e8b00960e88658'
PARENT_SHA='12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b'
CONTRACT=ROOT/'docs/CONTACT_ONSET_TIME_STUDY_CONTRACT_20261005.md'
COUNTS={'four':4,'eight':8}
SPAN=7.8125e-6
COMPLETE='COMPLETED_CONTACT_ONSET_RESOLUTION_DIAGNOSTIC_NOT_ACCEPTED'
THREADS=('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')


def bind(mapping):
    campaign.bind_map({},mapping)


def references(folder,count,dt,c):
    rows=post.read(folder/'candidate_history.json')
    require(len(rows)==count,'reference trajectory is not the registered full length')
    prefix=c['prefix'];previous=c['parent_hash']
    for index,row in enumerate(rows,1):
        require(row['step']==index and row['dt_s']==dt and row['parent_sha256']==previous and
            row['global_step']==c['start_step']+index and
            abs(row['time_s']-(c['start']+index*dt))<1e-15,'reference time/ancestry differs')
        for name,digest in row['artifacts'].items():
            path=(folder/name).resolve()
            require(path.is_relative_to((folder/f'step_{index:03d}').resolve()) and core.sha(path)==digest,
                'reference artifact changed')
        path=folder/row['checkpoint']
        require(path==folder/f'step_{index:03d}/NOT_RESTARTABLE_candidate.npz' and
            core.sha(path)==row['sha256'],'reference checkpoint identity differs')
        saved=core.load(str(path));state=saved['solver_data']
        require(state['accepted'] is False and state['restart_allowed'] is False and
            state['parameters']==c['p'] and state['classification']==campaign.model.CLASSIFICATION and
            state['time_s']==row['time_s'] and state['step']==row['global_step'] and
            state['interval_energy']['records'][:-1]==prefix,'reference physical history changed')
        campaign.model.discrete_process.validate_interval_energy_state(state['interval_energy'])
        prefix=state['interval_energy']['records'];previous=row['sha256']
    return rows


def preflight(stage,base):
    require(stage in COUNTS,'only the two registered resolutions exist')
    c=campaign.load_context(BASELINE)
    require(c['parent_hash']==PARENT_SHA and c['start_step']==14 and abs(c['start']-62.5e-6)<1e-15,
        'exact accepted parent required; no failed candidate restart')
    audit=post.read(AUDIT)
    require(core.sha(AUDIT)==AUDIT_SHA and audit['status']=='PASS_EXPANDED_HOT_CAMPAIGN_READBACK_NOT_PHYSICAL_VALIDATION'
        and audit['added_accepted_steps']==0 and audit['accepted_global_steps']==14 and
        audit['window_completed'] is False and audit['stages']['fine']['status']=='FAILED_NO_NEW_ACCEPTED_STATE',
        'original failed same-parent study must remain failed')
    bind(audit['evidence_sha256'])
    refs=dict(audit['evidence_sha256']);refs[rel(AUDIT)]=AUDIT_SHA
    measurements=list(c['protocol']['budget']['measurements'])
    for name in ('coarse','fine'):
        folder=BASELINE/name;manifest=post.read(folder/'manifest.json');rows=post.read(folder/'candidate_history.json')
        require(manifest['candidate_steps']==len(rows)==2,'baseline complete candidate count changed')
        measurements.append(dict(stage='baseline_'+name,maximum_complete_candidate_s=max(r['wall_seconds'] for r in rows),
            stage_wall_per_complete_candidate_s=manifest['wall_seconds']/len(rows)))
    baseline_rows=references(BASELINE/'fine',2,SPAN/2,c)
    folder=BASELINE/'fine'
    if stage=='eight':
        folder=base/'four';manifest=post.read(folder/'manifest.json');receipt=post.read(folder/'supervision.json')
        require(manifest['status']==COMPLETE and manifest['physical_trajectory_completed'] and
            manifest['parent_sha256']==PARENT_SHA and manifest['accepted_steps']==0 and
            receipt['exit_code']==0 and receipt['source_inputs_still_match'] and
            receipt['manifest_sha256']==core.sha(folder/'manifest.json') and
            receipt['log_sha256']==core.sha(folder/'run.log'),'four-step diagnostic must physically complete first')
        require(manifest['source_sha256']==post.read(folder/'protocol.json')['source_sha256'],
            'previous diagnostic source map differs')
        bind(manifest['source_sha256']);bind(manifest['reference_sha256'])
        refs.update({rel(path):core.sha(path) for path in folder.rglob('*') if path.is_file()})
        prior=references(folder,4,SPAN/4,c)
        measurements.append(dict(stage='four',maximum_complete_candidate_s=max(r['wall_seconds'] for r in prior),
            stage_wall_per_complete_candidate_s=manifest['wall_seconds']/len(prior)))
    else:prior=baseline_rows
    frozen=dict(c['frozen'])
    campaign.bind_map(frozen,{rel(path):core.sha(path) for path in (Path(__file__).resolve(),CONTRACT,AUDIT)})
    measured=max(max(r['maximum_complete_candidate_s'],r['stage_wall_per_complete_candidate_s']) for r in measurements)
    count=COUNTS[stage];cost=measured+1.;forecast=30+1.5*count*cost
    protocol=dict(schema='BOUNDED_CONTACT_ONSET_TIME_STUDY_V1',stage=stage,
        source_sha256=frozen,reference_sha256=refs,parent_checkpoint=rel(c['parent']),parent_sha256=PARENT_SHA,
        independent_origin_sha256=campaign.BRANCH_ORIGIN_SHA,parent_time_s=c['start'],parent_global_step=c['start_step'],
        baseline_readback=rel(AUDIT),baseline_readback_sha256=AUDIT_SHA,
        reference_directory=rel(folder),reference_count=len(prior),step_count=count,dt_s=SPAN/count,
        endpoint_s=c['start']+SPAN,measurements=measurements,measured_cost_s=measured,per_candidate_cost_s=cost,
        forecast_s=forecast,forecast_limit_s=450.,internal_limit_s=540.,external_limit_s=600.,
        status='REGISTERED_DIAGNOSTIC' if forecast<=450. else 'BLOCKED_COST_NO_RUN',
        accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.,
        physical_evaluations_in_preflight=0,nnps_evaluations_in_preflight=0)
    c.update(frozen=frozen,refs=refs,selection=None,study_protocol=protocol,reference_rows=prior,reference_directory=folder,
        study_count=count,study_dt=SPAN/count,study_per_step=cost)
    return c


def compare_pair(c,number,pending,seeds,directory):
    row=c['reference_rows'][number-1];folder=c['reference_directory']
    data=core.load(str(folder/row['checkpoint']));coarse,cs=data['arrays']['solid'],data['solver_data']
    fine,state=pending[-1]['pa'],pending[-1]['state']
    require(abs(cs['time_s']-state['time_s'])<1e-15,'comparison endpoints differ')
    seeds.update(e['particle_id'] for e in row['events']);seeds.update(row['blocked_particle_ids'])
    seeds.update(cs['bore_observation']['crossing_ids']);graphs=[]
    with np.load(folder/f'step_{number:03d}/visible_graph.npz') as z:
        for a,b in (('i','j'),('raw_i','raw_j')):graphs.append((z['particle_id'].copy(),z[a].copy(),z[b].copy()))
    for item in pending:
        for a,b in (('i','j'),('raw_i','raw_j')):graphs.append((item['pa'].particle_id,item['result'][a],item['result'][b]))
    selected=event.one_hop_ids(seeds,graphs);glob=post.compare(fine,coarse,c['p']);local=event.local_compare(fine,coarse,selected)
    index=np.flatnonzero(np.isin(fine.particle_id,selected))
    absolute=[]
    for k in index:
        absolute.append(dict(particle_id=int(fine.particle_id[k]),stress_norm_difference_Pa=float(np.linalg.norm(
            fine.S.reshape(-1,9)[k]-coarse.S.reshape(-1,9)[k])),fine_stress_norm_Pa=float(np.linalg.norm(fine.S.reshape(-1,9)[k])),
            plastic_strain_difference=float(abs(fine.ep[k]-coarse.ep[k])),fine_plastic_strain=float(fine.ep[k])))
    same=state['bore_observation']['crossing_ids']==cs['bore_observation']['crossing_ids']
    result=dict(pair=number,time_s=state['time_s'],parent_sha256=PARENT_SHA,
        coarse_reference_sha256=row['sha256'],fine_candidate_sha256=[core.sha(v['checkpoint']) for v in pending],
        cumulative_event_and_blocked_ids=sorted(seeds),global_screen=glob,local_screen=local,
        local_absolute_differences=absolute,crossing_ids_match=same,
        passed=glob['pass'] and local['passed'] and same,accepted=False,restart_allowed=False)
    core.write(directory/f'pair_{number:03d}_screen.json',result)
    return result


def reduction_observation(stage,base,current):
    previous_path=(BASELINE/'fine/pair_001_temporal_screen.json' if stage=='four' else base/'four/pair_002_screen.json')
    previous=post.read(previous_path)
    def metrics(screen):
        out={key:screen['global_screen'][key] for key in ('position_max_m','temperature_max_K','stress_relative','plastic_relative')}
        out.update({'local_'+key:max(row[key] for row in screen['local_screen']['rows'])
            for key in ('stress_relative','plastic_relative')})
        return out
    before,after=metrics(previous),metrics(current)
    return dict(previous_screen=rel(previous_path),previous_screen_sha256=core.sha(previous_path),
        previous_errors=before,current_errors=after,
        current_over_previous={key:after[key]/before[key] if before[key]>0 else None for key in before},
        identical_local_id_set=previous['local_screen']['particle_ids']==current['local_screen']['particle_ids'],
        denominator_floors=dict(stress_Pa=1e6,plastic_strain=1e-5),
        scope='EMPIRICAL_ADJACENT_RESOLUTION_DIFFERENCES_NOT_ASSUMED_ORDER_OR_EXTRAPOLATED_SOLUTION')


def worker(stage,base):
    directory=base/stage;protocol=post.read(directory/'protocol.json');started=time.perf_counter()
    rows=[];screens=[];code=2;report=dict(status='RUNNING_DIAGNOSTIC',stage=stage,accepted=False,restart_allowed=False,
        accepted_steps=0,accepted_time_advanced_s=0.,physical_trajectory_completed=False)
    def stop():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_DIAGNOSTIC_NOT_ACCEPTED',completed_candidates=len(rows)))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        require(all(os.environ.get(k)=='1' for k in THREADS),'single thread environment required')
        c=preflight(stage,base);require(c['study_protocol']==protocol and protocol['status']=='REGISTERED_DIAGNOSTIC',
            'immutable study protocol differs or is unaffordable')
        report.update(protocol,physical_trajectory_completed=False,status='RUNNING_DIAGNOSTIC')
        pa,state=c['pa'],c['state'];parent_hash=PARENT_SHA;seeds=set(c['seeds']);pending=[]
        per_step=c['study_per_step'];per_root=protocol['measured_cost_s']
        core.write(directory/'manifest.json',report);core.write(directory/'candidate_history.json',rows)
        for number in range(1,c['study_count']+1):
            require(time.perf_counter()+1.5*per_step+20<=started+525.,'remaining measured stage budget exhausted')
            sd=directory/f'step_{number:03d}';sd.mkdir()
            item=qualified.candidate(pa,state,c,c['study_dt'],sd,started+525.,per_root)
            require(item['state']['accepted'] is False and item['state']['restart_allowed'] is False,
                'diagnostic attempted publication')
            row=old.candidate_row(item,number,parent_hash,directory);rows.append(row);pending.append(item)
            core.write(directory/'candidate_history.json',rows)
            pa,state=item['pa'],item['state'];parent_hash=row['sha256']
            seeds.update(e['particle_id'] for e in row['events']);seeds.update(row['blocked_particle_ids'])
            seeds.update(state['bore_observation']['crossing_ids'])
            if number%2==0:
                screen=compare_pair(c,number//2,pending,seeds,directory);screens.append(screen);pending=[]
            campaign.check_context(c);per_root=max(per_root,item['wall_seconds'])
            per_step=max(per_step,per_root,(time.perf_counter()-started)/len(rows))
            report.update(candidate_steps=len(rows),candidate_time_s=state['time_s'])
            core.write(directory/'manifest.json',report)
            print(json.dumps(dict(candidate=number,time_s=state['time_s'],wall_s=item['wall_seconds'],accepted=False)),flush=True)
        require(abs(state['time_s']-protocol['endpoint_s'])<1e-15,'registered endpoint not reached')
        report.update(status=COMPLETE,physical_trajectory_completed=True,
            temporal_screens_all_pass=all(v['passed'] for v in screens),pair_screens=screens,
            terminal_error_reduction=reduction_observation(stage,base,screens[-1]),
            interval_energy=state['interval_energy'],crossing_ids=state['bore_observation']['crossing_ids'],
            candidate_end_checkpoint=rel(item['checkpoint']),candidate_end_sha256=core.sha(item['checkpoint']))
        code=0
    except Exception as error:
        report.update(status='FAIL_CONTACT_ONSET_DIAGNOSTIC_NOT_ACCEPTED',error=repr(error),traceback=traceback.format_exc())
    finally:
        timer.cancel();report.update(wall_seconds=time.perf_counter()-started,candidate_steps=len(rows),
            completed_screen_count=len(screens),peak_working_set_bytes=core.memory_bytes())
        report['artifacts_sha256']={rel(p):core.sha(p) for p in directory.rglob('*') if p.is_file() and
            p.name not in ('manifest.json','run.log','supervision.json')}
        core.write(directory/'manifest.json',report)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',choices=tuple(COUNTS),required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS);args=parser.parse_args()
    base=campaign.project(args.output);require(base.is_relative_to(ROOT/'outputs') and base!=ROOT/'outputs','output outside project')
    if args.worker:return worker(args.stage,base)
    c=preflight(args.stage,base);protocol=c['study_protocol']
    if args.preflight_only:
        print(json.dumps({k:protocol[k] for k in ('status','step_count','dt_s','forecast_s','parent_sha256')},indent=2));return 0
    directory=base/args.stage;directory.mkdir(parents=True,exist_ok=False);campaign.write_new(directory/'protocol.json',protocol)
    if protocol['status']!='REGISTERED_DIAGNOSTIC':print(protocol['status']);return 2
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--output',rel(base),'--worker']
    started=time.perf_counter();code=124
    with (directory/'run.log').open('x',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,env=dict(os.environ,**{k:'1' for k in THREADS}),
            stdout=stream,stderr=subprocess.STDOUT,timeout=600.).returncode
        except subprocess.TimeoutExpired:stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    error=None
    try:bind(protocol['source_sha256']);bind(protocol['reference_sha256'])
    except Exception as problem:error=repr(problem);code=2
    manifest=directory/'manifest.json'
    campaign.write_new(directory/'supervision.json',dict(command=command,exit_code=code,
        wall_seconds=time.perf_counter()-started,source_inputs_still_match=error is None,source_error=error,
        log_sha256=core.sha(directory/'run.log'),manifest_sha256=core.sha(manifest) if manifest.exists() else None,
        worker_limit_s=540.,external_limit_s=600.,accepted=False,restart_allowed=False))
    print(json.dumps(dict(stage=args.stage,exit_code=code,output=rel(directory),
        status=post.read(manifest)['status'] if manifest.exists() else 'NO_COMPLETE_MANIFEST')));return code


if __name__=='__main__':raise SystemExit(main())

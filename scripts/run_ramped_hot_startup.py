"""One bounded independent startup experiment; diagnostic states never publish."""
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
import ramped_hot_candidate as complete
import ramped_hot_trial as model
import qualify_expanded_hot_domain as qualification

core,post,event,old = complete.core,complete.post,complete.event,complete.old
ROOT=core.ROOT
require=complete.require
CONTRACT=ROOT/'docs/RAMPED_HOT_STARTUP_DIAGNOSTIC_CONTRACT_20261005.md'
CHECK=ROOT/'outputs/ramped_hot_motion_20261005/v1/receipt.json'
OLD_COST=ROOT/'outputs/expanded_hot_domain_cost_20261004/v3/dt_15_625us/cost/manifest.json'
OLD_RANK=ROOT/'outputs/expanded_hot_qualification_20261005/v1/summary.json'
OLD_ROOT=ROOT/'outputs/expanded_hot_campaign_v2_20261005/window_001/fine/candidate_history.json'
COUNTS={'coarse':2,'fine':4}
DTS={'coarse':.0005,'fine':.00025}
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')
PASS='PASS_RAMP_PHYSICAL_STAGE_NOT_ACCEPTED'


def rel(path): return Path(path).resolve().relative_to(ROOT).as_posix()
def read(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def load(path):
    data=core.load(str(path));return data['arrays']['solid'],data['solver_data']
def new(path,value):
    with Path(path).open('x',encoding='utf-8') as stream: json.dump(core.clean(value),stream,indent=2,allow_nan=False);stream.write('\n')
def verify(mapping): post.check_sources(mapping)


def source_inputs():
    # Import closure plus every backend module; lazy imports live in that backend.
    paths={Path(__file__).resolve(),CONTRACT,CHECK,OLD_COST,OLD_RANK,OLD_ROOT}
    paths.add(ROOT/'scripts/check_ramped_hot_motion.py')
    for module in list(sys.modules.values()):
        file=getattr(module,'__file__',None)
        if file:
            path=Path(file).resolve()
            if path.suffix=='.py' and path.is_relative_to(ROOT/'scripts'):paths.add(path)
    mapping={}
    def merge(values):
        verify(values)
        for name,digest in values.items():
            require(name not in mapping or mapping[name]==digest,'conflicting frozen source: '+name)
            mapping[name]=digest
    merge(core.frozen_sources());merge(model.source_hashes())
    for receipt_path in (CHECK,OLD_COST,OLD_RANK):
        receipt=read(receipt_path)
        for key in ('source_sha256','source_and_input_sha256','output_sha256'):
            merge(receipt.get(key,{}))
    merge({rel(path):core.sha(path) for path in paths})
    verify(mapping)
    check=read(CHECK)
    require(check.get('status')=='PASS_RAMPED_MOTION_COMPONENT_QUALIFICATION_NOT_PROCESS_ACCEPTANCE' and
        all(v['status']=='PASS' for v in check['checks']),'ramp motion qualification required')
    return mapping


def stage_dir(base,stage,chunk):
    return base/'cost' if stage=='cost' else base/f'chunk_{chunk:03d}'/stage


def receipt(folder):
    manifest=read(folder/'manifest.json');supervision=read(folder/'supervision.json')
    require(supervision['exit_code']==0 and supervision['source_inputs_still_match'] and
        supervision['manifest_sha256']==core.sha(folder/'manifest.json') and
        supervision['log_sha256']==core.sha(folder/'run.log'),'failed or changed stage receipt')
    verify(manifest['source_sha256']);verify(manifest['artifacts_sha256'])
    return manifest,supervision


def preflight(base,stage,chunk):
    require(1<=chunk<=4,'only four registered chunks')
    frozen=source_inputs(); refs={}; measured=max(r['wall_seconds'] for r in read(OLD_ROOT))
    elapsed=0.
    for path in base.glob('**/supervision.json') if base.exists() else ():
        d=read(path);elapsed+=d['wall_seconds']
        refs[rel(path)]=core.sha(path)
        manifest_path=path.parent/'manifest.json'
        if manifest_path.exists():
            m=read(manifest_path)
            refs[rel(manifest_path)]=core.sha(manifest_path)
            rows=read(path.parent/'candidate_history.json') if (path.parent/'candidate_history.json').exists() else []
            if (path.parent/'candidate_history.json').exists():
                refs[rel(path.parent/'candidate_history.json')]=core.sha(path.parent/'candidate_history.json')
            if rows:measured=max(measured,max(r['wall_seconds'] for r in rows),m['wall_seconds']/len(rows))
    if stage=='cost':
        require(chunk==1 and not base.exists(),'independent experiment requires a new output')
        c=read(OLD_COST);r=read(OLD_RANK)
        forecast=30+1.5*(2*c['full_map_wall_s']+c['tangent_linear_wall_s']+r['wall_seconds'])+120
        parent=None;count=0;dt=DTS['coarse']
    else:
        cost,_=receipt(base/'cost');require(cost['status']=='PASS_RAMP_COST_INITIAL_RANK_SAME_MAP_NOT_ACCEPTED','cost qualification absent')
        require(cost['source_sha256']==frozen,'experiment implementation changed')
        parent=base/'initial/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
        if chunk>1:
            audit=base/f'chunk_{chunk-1:03d}/readback/summary.json';a=read(audit)
            require(a['status']=='PASS_RAMP_CHUNK_SAVED_READBACK_NOT_ACCEPTED' and a['temporal_screens_all_pass'], 'previous chunk lacks independent readback')
            verify(a['evidence_sha256'])
            for name,digest in a['evidence_sha256'].items():
                require(name not in refs or refs[name]==digest,'conflicting prior audit evidence')
                refs[name]=digest
            refs[rel(audit)]=core.sha(audit)
            prior=stage_dir(base,stage,chunk-1);m,_=receipt(prior)
            require(m['status']==PASS,'previous path did not physically complete')
            parent=ROOT/m['candidate_end_checkpoint']
        if stage=='fine':
            coarse=stage_dir(base,'coarse',chunk);m,_=receipt(coarse)
            require(m['status']==PASS,'coarse path failed; no automatic alternate attempt')
            refs.update({rel(p):core.sha(p) for p in coarse.rglob('*') if p.is_file()})
        refs.update({rel(p):core.sha(p) for p in (base/'cost').rglob('*') if p.is_file()})
        refs[rel(parent)]=core.sha(parent)
        pa,state=load(parent)
        require(state['classification']==model.CLASSIFICATION and state['accepted'] is False and state['restart_allowed'] is False,
            'only this diagnostic path may be chained inside the experiment')
        require(abs(state['time_s']-(chunk-1)*.001)<1e-15 and state['step']==(chunk-1)*COUNTS[stage], 'diagnostic path ancestry differs')
        count=COUNTS[stage];dt=DTS[stage];forecast=30+1.5*count*(measured+1.)
    require(forecast<=450. and elapsed+forecast<=1200.,'bounded stage or cumulative experiment budget exceeded')
    external_limit=min(600.,1200.-elapsed)
    internal_limit=min(540.,external_limit-5.)
    require(forecast<=internal_limit-15.,'insufficient hard cumulative budget margin')
    protocol=dict(schema='RAMPED_HOT_STARTUP_DIAGNOSTIC_V1',stage=stage,chunk=chunk,source_sha256=frozen,
        reference_sha256=refs,parent_checkpoint=rel(parent) if parent else None,parent_sha256=core.sha(parent) if parent else None,
        step_count=count,dt_s=dt,endpoint_s=chunk*.001 if stage!='cost' else 0.,
        comparison_scope='COMMON_INITIAL_CONTINUOUS_TWO_RESOLUTION',measured_candidate_cost_s=measured,
        forecast_s=forecast,previous_supervised_wall_s=elapsed,total_wall_limit_s=1200.,
        internal_limit_s=internal_limit,external_limit_s=external_limit,accepted=False,restart_allowed=False,accepted_time_advanced_s=0.)
    return protocol


def cost_stage(base,directory,report):
    pa,state=model.independent_initial_state()
    init=base/'initial';init.mkdir()
    path=init/'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
    core.checkpoint(path,pa,state);core.verify_snapshot(path,pa,state)
    core.write(init/'state.json',state)
    trial=model.make_trial(pa,state,DTS['coarse'])
    # Old initial-rank function only queries t=0 pose, checked by motion qualification.
    rank=qualification.initial_rank(trial,directory)
    begin=time.perf_counter();r=trial.evaluate(trial.prescribed(DTS['coarse']),DTS['coarse']);map_wall=time.perf_counter()-begin
    core.journal(directory,trial,r,'UNSOLVED_COST_ONLY');old.save_graph(directory,trial,r)
    begin=time.perf_counter();J=trial.tangent(r);tangent_wall=time.perf_counter()-begin
    require(np.isfinite(J).all(),'nonfinite initial tangent')
    scalar=model.make_trial(pa,state,DTS['coarse']);scalar.contact=qualification.contact_factory(old.ChainRuleDiscreteContact,trial.p)
    sr=scalar.evaluate(r['displacement'],DTS['coarse'])
    differences=complete.batch.candidate_equivalence(dict(sr,energy_comparison=core.energy_diagnostics(scalar,sr)),
        dict(r,energy_comparison=core.energy_diagnostics(trial,r)))
    require(sr['graph_fingerprint']==r['graph_fingerprint'],'initial scalar graph differs')
    report.update(status='PASS_RAMP_COST_INITIAL_RANK_SAME_MAP_NOT_ACCEPTED',classification=model.CLASSIFICATION,
        initial_elastic_rank=rank,full_map_wall_s=map_wall,tangent_wall_s=tangent_wall,scalar_batch_differences=differences,
        initial_checkpoint=rel(path),initial_sha256=core.sha(path),particle_count=len(pa.x),free_dof=int(state['free'].sum()),
        dt_s=DTS['coarse'],stable_timestep_s='NOT_ESTABLISHED',physical_full_map_evaluations=2,
        initial_state_unchanged=core.clean(load(path)[1])==core.clean(state),static_adaptation=model.static_adaptation_audit())


def pair_screen(base,chunk,number,fine_item,pending):
    folder=stage_dir(base,'coarse',chunk);row=read(folder/'candidate_history.json')[number-1]
    coarse,cs=load(folder/row['checkpoint']);fine,state=fine_item['pa'],fine_item['state']
    require(abs(cs['time_s']-state['time_s'])<1e-15,'shared endpoints differ')
    global_screen=post.compare(fine,coarse,state['parameters'])
    local_screen=event.local_compare(fine,coarse,sorted(int(v) for v in fine.particle_id))
    same=state['bore_observation']['crossing_ids']==cs['bore_observation']['crossing_ids']
    return dict(pair=number,time_s=state['time_s'],comparison_scope='COMMON_INITIAL_CONTINUOUS_TWO_RESOLUTION',
        coarse_checkpoint=rel(folder/row['checkpoint']),coarse_sha256=row['sha256'],
        fine_checkpoint_sha256=[core.sha(v['checkpoint']) for v in pending],global_screen=global_screen,
        local_screen=local_screen,crossing_ids_match=same,passed=bool(global_screen['pass'] and local_screen['passed'] and same),
        accepted=False,restart_allowed=False)


def worker(base,stage,chunk):
    directory=stage_dir(base,stage,chunk);started=time.perf_counter();rows=[];screens=[];code=2
    protocol=read(directory/'protocol.json')
    report=dict(protocol,status='RUNNING_RAMP_DIAGNOSTIC',accepted_steps=0,accepted_time_advanced_s=0.)
    def timeout():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NOT_ACCEPTED',completed_candidates=len(rows)))
        finally:os._exit(124)
    deadline=started+protocol['internal_limit_s']-15.
    timer=threading.Timer(protocol['internal_limit_s'],timeout);timer.daemon=True;timer.start()
    try:
        require(all(os.environ.get(k)=='1' for k in THREADS),'single thread environment required')
        verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
        require(source_inputs()==protocol['source_sha256'],'frozen inputs changed')
        core.write(directory/'manifest.json',report);core.write(directory/'candidate_history.json',rows)
        if stage=='cost':cost_stage(base,directory,report)
        else:
            pa,state=load(ROOT/protocol['parent_checkpoint']);direct_hash=protocol['parent_sha256'];pending=[]
            per_root=protocol['measured_candidate_cost_s'];dt=protocol['dt_s']
            for number in range(1,protocol['step_count']+1):
                require(time.perf_counter()+1.5*(per_root+1)+20<=deadline,'remaining measured stage budget exhausted')
                sd=directory/f'step_{number:03d}';sd.mkdir()
                item=complete.candidate(pa,state,dt,sd,deadline,per_root,protocol['source_sha256'])
                require(item['state']['accepted'] is False and item['state']['restart_allowed'] is False,'diagnostic publication attempted')
                row=old.candidate_row(item,number,direct_hash,directory);row['response']=item['response']
                rows.append(row);pending.append(item);core.write(directory/'candidate_history.json',rows)
                pa,state=item['pa'],item['state'];direct_hash=row['sha256']
                if stage=='fine' and number%2==0:
                    screen=pair_screen(base,chunk,number//2,item,pending);screens.append(screen);pending=[]
                    core.write(directory/f'pair_{number//2:03d}_screen.json',screen)
                    require(screen['passed'],'registered common-initial temporal screen failed')
                per_root=max(per_root,item['wall_seconds'],(time.perf_counter()-started)/len(rows))
                verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
                report.update(candidate_steps=len(rows),candidate_time_s=state['time_s']);core.write(directory/'manifest.json',report)
                print(json.dumps(dict(candidate=number,time_s=state['time_s'],wall_s=item['wall_seconds'],response=item['response'],accepted=False)),flush=True)
            require(abs(state['time_s']-protocol['endpoint_s'])<1e-15,'endpoint not reached')
            report.update(status=PASS,candidate_end_checkpoint=rel(item['checkpoint']),candidate_end_sha256=direct_hash,
                interval_energy=state['interval_energy'],crossing_ids=state['bore_observation']['crossing_ids'],
                temporal_screens_all_pass=all(v['passed'] for v in screens),pair_screens=screens)
        verify(protocol['source_sha256']);verify(protocol['reference_sha256']);code=0
    except Exception as error:
        report.update(status='FAIL_RAMP_STARTUP_DIAGNOSTIC_NOT_ACCEPTED',error=repr(error),traceback=traceback.format_exc())
    finally:
        timer.cancel();report.update(wall_seconds=time.perf_counter()-started,candidate_steps=len(rows),pair_screens=screens,
            peak_working_set_bytes=core.memory_bytes(),accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.)
        report['artifacts_sha256']={rel(p):core.sha(p) for p in directory.rglob('*') if p.is_file() and
            p.name not in ('manifest.json','run.log','supervision.json')}
        if stage=='cost' and (base/'initial').exists():report['artifacts_sha256'].update({rel(p):core.sha(p) for p in (base/'initial').iterdir() if p.is_file()})
        core.write(directory/'manifest.json',report)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('cost','coarse','fine'),required=True)
    parser.add_argument('--chunk',type=int,default=1);parser.add_argument('--output',required=True)
    parser.add_argument('--preflight-only',action='store_true');parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();base=(ROOT/args.output).resolve()
    require(base.is_relative_to(ROOT/'outputs') and base!=ROOT/'outputs','new project output required')
    if args.worker:return worker(base,args.stage,args.chunk)
    protocol=preflight(base,args.stage,args.chunk)
    if args.preflight_only:print(json.dumps({k:protocol[k] for k in ('stage','chunk','step_count','dt_s','forecast_s','previous_supervised_wall_s')}));return 0
    directory=stage_dir(base,args.stage,args.chunk);directory.mkdir(parents=True,exist_ok=False);new(directory/'protocol.json',protocol)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--chunk',str(args.chunk),'--output',rel(base),'--worker']
    begin=time.perf_counter();code=124
    with (directory/'run.log').open('x',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,env=dict(os.environ,**{k:'1' for k in THREADS}),stdout=stream,stderr=subprocess.STDOUT,timeout=protocol['external_limit_s']).returncode
        except subprocess.TimeoutExpired:stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    error=None
    try:verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
    except Exception as problem:error=repr(problem);code=2
    manifest=directory/'manifest.json'
    new(directory/'supervision.json',dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-begin,
        source_inputs_still_match=error is None,source_error=error,manifest_sha256=core.sha(manifest) if manifest.exists() else None,
        log_sha256=core.sha(directory/'run.log'),internal_limit_s=protocol['internal_limit_s'],external_limit_s=protocol['external_limit_s'],accepted=False,restart_allowed=False))
    print(json.dumps(dict(stage=args.stage,chunk=args.chunk,exit_code=code,status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    return code


if __name__=='__main__':raise SystemExit(main())

"""Bounded same-parent timestep study; frozen ramp physics, no publication."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import ramped_hot_candidate as complete
import ramped_hot_trial as model
import bounded_time_control as policy

core,post,event,old=complete.core,complete.post,complete.event,complete.old
ROOT=core.ROOT
require=complete.require
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')
ORIGIN=ROOT/'outputs/ramped_hot_startup_20261005/v1/chunk_003/fine/step_004/NOT_RESTARTABLE_candidate.npz'
ORIGIN_SHA='43b2957db88cb3db23f41dadd4a0b6a9c17f3b73aa92ebb4741f51b4ff03775f'
ORIGIN_AUDIT=ROOT/'outputs/ramped_hot_startup_20261005/v1/chunk_003/readback/summary.json'
ORIGIN_AUDIT_SHA='ef2beb0c45691ae755c342191a7077c7f3c4c9328eaa5d173fe8694979dad859'
OLD_FAILURE=ROOT/'outputs/ramped_hot_startup_20261005/v1/chunk_004/readback/summary.json'
OLD_FAILURE_SHA='25746d994a5b2d98c68c3e6929dee005b507e3f3527f52ae23ed154b815a50be'
CHECK=ROOT/'outputs/bounded_time_control_20261005/v2/receipt.json'
CONTRACT=ROOT/'docs/LOCAL_TIME_CONTROL_CONTRACT_20261005.md'
READBACK=ROOT/'scripts/readback_local_time_control.py'
OLD_ROOT=ROOT/'outputs/expanded_hot_campaign_v2_20261005/window_001/fine/candidate_history.json'
DIAGNOSTIC=dict(accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.)


def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def load(path):
    data=core.load(str(path));return data['arrays']['solid'],data['solver_data']
def verify(mapping):post.check_sources(mapping)
def new(path,value):
    with Path(path).open('x',encoding='utf-8') as stream:
        json.dump(core.clean(value),stream,indent=2,allow_nan=False);stream.write('\n')
def merge(mapping,values):
    verify(values)
    for name,digest in values.items():
        require(name not in mapping or mapping[name]==digest,'conflicting bound evidence: '+name)
        mapping[name]=digest


def inputs():
    require(core.sha(ORIGIN)==ORIGIN_SHA and core.sha(ORIGIN_AUDIT)==ORIGIN_AUDIT_SHA and
        core.sha(OLD_FAILURE)==OLD_FAILURE_SHA,'registered prior evidence changed')
    pa,state=load(ORIGIN);model.validate_ramped_state(state)
    require(len(pa.x)==918 and state['step']==12 and abs(state['time_s']-.003)<1e-15 and
        state['accepted'] is False and state['restart_allowed'] is False,'wrong diagnostic origin')
    audit=read(ORIGIN_AUDIT);failure=read(OLD_FAILURE)
    require(audit['temporal_screens_all_pass'] and not failure['temporal_screens_all_pass'],
        'prior success/failure distinction changed')
    sources={};refs={}
    merge(sources,core.frozen_sources());merge(sources,state['implementation_source_sha256'])
    for path in (ORIGIN_AUDIT,OLD_FAILURE):
        merge(refs,read(path)['evidence_sha256']);merge(refs,{rel(path):core.sha(path)})
    merge(refs,{rel(ORIGIN):ORIGIN_SHA,rel(OLD_ROOT):core.sha(OLD_ROOT)})
    check=read(CHECK)
    require(check['status']=='PASS_BOUNDED_TIME_CONTROLLER_LOGIC_NOT_PROCESS_ACCEPTANCE' and
        check['exit_code']==0,'pure controller checks required')
    merge(sources,check['source_sha256'])
    paths={Path(__file__).resolve(),Path(policy.__file__).resolve(),CHECK,CONTRACT,READBACK,
        ROOT/'scripts/readback_ramped_hot_startup.py',ROOT/'scripts/readback_expanded_hot_campaign_window_v2.py',
        ROOT/'scripts/readback_visible_support_process.py'}
    for module in list(sys.modules.values()):
        file=getattr(module,'__file__',None)
        if file:
            path=Path(file).resolve()
            if path.suffix=='.py' and path.is_relative_to(ROOT/'scripts'):paths.add(path)
    merge(sources,{rel(path):core.sha(path) for path in paths})
    return sources,refs,state


def audited_stage(folder,refs):
    manifest=read(folder/'manifest.json');supervision=read(folder/'supervision.json')
    audit=read(folder/'readback/summary.json')
    require(supervision['exit_code']==0 and supervision['source_inputs_still_match'] and
        supervision['manifest_sha256']==core.sha(folder/'manifest.json') and
        supervision['log_sha256']==core.sha(folder/'run.log'),'stage supervision failed')
    require(audit['status']=='PASS_SAVED_LOCAL_TIME_READBACK_NOT_ACCEPTED','independent stage readback absent')
    merge(refs,audit['evidence_sha256']);merge(refs,{rel(folder/'readback/summary.json'):core.sha(folder/'readback/summary.json')})
    verify(manifest['source_sha256']);verify(manifest['artifacts_sha256'])
    return manifest,supervision


def preflight(base,phase,index):
    require(1<=index<=10,'bounded stage index required')
    frozen,refs,origin=inputs()
    controller=policy.LocalTimeController(ORIGIN_SHA,origin['time_s'])
    parent=ORIGIN;completed=[];primary=[];elapsed=0.;count=0
    measured=max(row['wall_seconds'] for row in read(OLD_ROOT))
    folders=sorted(base.glob('attempt_*')) if base.exists() else []
    require(all(p.is_dir() for p in folders),'unexpected attempt file')
    if phase=='primary':require(len(folders)==index-1,'primary attempts must be consecutive and fresh')
    for number,folder in enumerate(folders,1):
        require(folder.name==f'attempt_{number:03d}','missing primary attempt')
        m,s=audited_stage(folder,refs)
        require(m['source_sha256']==frozen and m['controller_before']==controller.snapshot(),'changed sources or controller ancestry')
        controller=policy.LocalTimeController.from_snapshot(m['controller_after'])
        rows=read(folder/'candidate_history.json')
        completed.append((m,s,rows))
        if m['status']=='PASS_LOCAL_PAIR_DIAGNOSTIC_NOT_ACCEPTED':
            parent=folder/'fine_002'/rows[-1]['checkpoint'];primary.append((folder,m,rows))
        else:require(m['status']=='TEMPORAL_REJECT_DIAGNOSTIC_NOT_ACCEPTED','physical failure forbids continuation')
    replay_plan=[];offset=0
    if phase=='primary':
        require(not list(base.glob('replay_*')),'primary cannot resume after independent replay')
        before=controller.snapshot();proposal=controller.propose();n=3
    else:
        require(controller.status=='PRIMARY_COMPLETE_INDEPENDENT_REPLAY_REQUIRED','primary interval not complete')
        before=controller.snapshot();proposal=None
        for grid in controller.replay_grid():
            folder,m,rows=primary[grid['primary_interval']-1]
            row=rows[grid['primary_fine_step']]
            reference=folder/row['role']/row['checkpoint']
            replay_plan.append(dict(grid,reference_checkpoint=rel(reference),reference_sha256=row['sha256']))
            merge(refs,{rel(reference):row['sha256']})
        parent=ORIGIN
        replays=sorted(base.glob('replay_*'))
        require(len(replays)==index-1,'replay stages must be consecutive and fresh')
        for number,folder in enumerate(replays,1):
            require(folder.name==f'replay_{number:03d}','missing replay stage')
            m,s=audited_stage(folder,refs);rows=read(folder/'candidate_history.json')
            require(m['source_sha256']==frozen and m['status']=='PASS_REPLAY_STAGE_DIAGNOSTIC_NOT_ACCEPTED' and
                m['replay_plan']==replay_plan and m['replay_offset']==offset,'replay failed or grid changed')
            offset+=len(rows);parent=folder/rows[-1]['role']/rows[-1]['checkpoint'];completed.append((m,s,rows))
        require(offset<len(replay_plan),'independent replay already complete')
        n=min(4,len(replay_plan)-offset)
    for m,s,rows in completed:
        elapsed+=s['wall_seconds'];count+=len(rows)
        if rows:measured=max(measured,max(r['wall_seconds'] for r in rows),s['wall_seconds']/len(rows))
    budget=policy.stage_budget(n,measured,elapsed,count)
    require(core.sha(parent)==(proposal['parent_sha256'] if proposal else core.sha(parent)),'selected parent mismatch')
    merge(refs,{rel(parent):core.sha(parent)})
    return dict(schema='LOCAL_TIME_CONTROL_DIAGNOSTIC_V1',phase=phase,index=index,source_sha256=frozen,
        reference_sha256=refs,parent_checkpoint=rel(parent),parent_sha256=core.sha(parent),
        controller_before=before,proposal=proposal,replay_plan=replay_plan,replay_offset=offset,
        step_count=n,**budget)


def compare(fine,fs,coarse,cs,fine_path,coarse_path,scope):
    require(abs(fs['time_s']-cs['time_s'])<1e-15,'shared endpoints differ')
    glob=post.compare(fine,coarse,fs['parameters'])
    local=event.local_compare(fine,coarse,sorted(int(v) for v in fine.particle_id))
    same=fs['crossing_ids']==cs['crossing_ids']
    return dict(time_s=fs['time_s'],comparison_scope=scope,coarse_checkpoint=rel(coarse_path),
        coarse_sha256=core.sha(coarse_path),fine_checkpoint=rel(fine_path),fine_sha256=core.sha(fine_path),
        global_screen=glob,local_screen=local,crossing_ids_match=same,
        passed=bool(glob['pass'] and local['passed'] and same),
        normalized_error=policy.normalized_screen_error(glob,local,fs['parameters']['dx']),**DIAGNOSTIC)


def worker(base,phase,index):
    directory=base/f"{'attempt' if phase=='primary' else 'replay'}_{index:03d}"
    protocol=read(directory/'protocol.json');started=time.perf_counter();rows=[];screens=[];code=2
    controller=policy.LocalTimeController.from_snapshot(protocol['controller_before'])
    proposal=controller.propose() if phase=='primary' else None
    report=dict(protocol,status='RUNNING_LOCAL_TIME_DIAGNOSTIC',controller_after=protocol['controller_before'],screens=[])
    def timeout():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NOT_ACCEPTED',candidate_steps=len(rows)))
        finally:os._exit(124)
    timer=threading.Timer(protocol['internal_limit_s'],timeout);timer.daemon=True;timer.start()
    deadline=started+protocol['internal_limit_s']-15.;per_root=protocol['measured_candidate_cost_s']
    def solve(role,path,dt):
        nonlocal per_root
        require(time.perf_counter()+1.5*(per_root+1)+20<=deadline,'insufficient candidate deadline margin')
        pa,state=load(path);role_dir=directory/role;sd=role_dir/'step_001';sd.mkdir(parents=True)
        item=complete.candidate(pa,state,dt,sd,deadline,per_root,protocol['source_sha256'])
        row=old.candidate_row(item,1,core.sha(path),role_dir);row.update(role=role,response=item['response'])
        rows.append(row);core.write(directory/'candidate_history.json',rows)
        per_root=max(per_root,row['wall_seconds'],(time.perf_counter()-started)/len(rows))
        verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
        print(json.dumps(dict(role=role,time_s=state['time_s']+dt,wall_seconds=row['wall_seconds'],
            new_plastic_count=len(item['response']['new_plastic_ids']),tool_loaded_count=len(item['response']['nonbase_loaded_ids']),accepted=False)),flush=True)
        return item
    try:
        require(all(os.environ.get(k)=='1' for k in THREADS),'single thread environment required')
        verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
        require(inputs()[0]==protocol['source_sha256'],'implementation changed after preflight')
        core.write(directory/'manifest.json',report);core.write(directory/'candidate_history.json',rows)
        parent=ROOT/protocol['parent_checkpoint']
        if phase=='primary':
            require(proposal==protocol['proposal'],'controller proposal changed')
            coarse=solve('coarse',parent,proposal['coarse_dt_s'])
            fine1=solve('fine_001',parent,proposal['fine_dt_s'])
            fine2=solve('fine_002',fine1['checkpoint'],proposal['fine_dt_s'])
            screen=compare(fine2['pa'],fine2['state'],coarse['pa'],coarse['state'],fine2['checkpoint'],coarse['checkpoint'],policy.COMPARISON_SCOPE)
            screens.append(screen);core.write(directory/'screen.json',screen)
            controller.resolve(proposal,physical_passed=True,temporal_passed=screen['passed'],
                normalized_error=screen['normalized_error'],fine2_sha256=core.sha(fine2['checkpoint']) if screen['passed'] else None)
            report.update(status='PASS_LOCAL_PAIR_DIAGNOSTIC_NOT_ACCEPTED' if screen['passed'] else 'TEMPORAL_REJECT_DIAGNOSTIC_NOT_ACCEPTED',
                screen=screen,controller_after=controller.snapshot())
        else:
            report['status']='PASS_REPLAY_STAGE_DIAGNOSTIC_NOT_ACCEPTED'
            for k,grid in enumerate(protocol['replay_plan'][protocol['replay_offset']:protocol['replay_offset']+protocol['step_count']],1):
                item=solve(f'replay_{k:03d}',parent,grid['dt_s']);parent=item['checkpoint']
                require(abs(item['state']['time_s']-grid['t1_s'])<1e-15,'replay grid endpoint changed')
                if grid['shared_primary_fine_endpoint']:
                    reference=ROOT/grid['reference_checkpoint'];cp,cs=load(reference)
                    screen=compare(item['pa'],item['state'],cp,cs,parent,reference,'COMMON_3MS_ORIGIN_INDEPENDENT_HALF_PRIMARY_FINE_GRID')
                    screens.append(screen);core.write(directory/f'screen_{k:03d}.json',screen)
                    if not screen['passed']:
                        report['status']='TEMPORAL_REJECT_REPLAY_NOT_QUALIFIED';break
        verify(protocol['source_sha256']);verify(protocol['reference_sha256']);code=0
    except Exception as error:
        if phase=='primary' and controller.pending is not None:
            controller.resolve(proposal,physical_passed=False,temporal_passed=False)
        report.update(status='FAIL_LOCAL_TIME_PHYSICAL_NOT_ACCEPTED',error=repr(error),traceback=traceback.format_exc(),
            controller_after=controller.snapshot())
    finally:
        timer.cancel();report.update(wall_seconds=time.perf_counter()-started,candidate_steps=len(rows),screens=screens,
            peak_working_set_bytes=core.memory_bytes(),**DIAGNOSTIC)
        report['artifacts_sha256']={rel(p):core.sha(p) for p in directory.rglob('*') if p.is_file() and
            p.name not in ('manifest.json','run.log','supervision.json')}
        core.write(directory/'manifest.json',report)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True);parser.add_argument('--phase',choices=('primary','replay'),required=True)
    parser.add_argument('--index',type=int,required=True);parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS);args=parser.parse_args()
    base=(ROOT/args.output).resolve()
    require(base.is_relative_to(ROOT/'outputs') and base!=ROOT/'outputs','project output directory required')
    if args.worker:return worker(base,args.phase,args.index)
    protocol=preflight(base,args.phase,args.index)
    if args.preflight_only:
        print(json.dumps({k:protocol[k] for k in ('phase','index','step_count','forecast_s','previous_supervised_wall_s','proposal')}));return 0
    directory=base/f"{'attempt' if args.phase=='primary' else 'replay'}_{args.index:03d}"
    directory.mkdir(parents=True,exist_ok=False);new(directory/'protocol.json',protocol)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--output',rel(base),'--phase',args.phase,'--index',str(args.index),'--worker']
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
        log_sha256=core.sha(directory/'run.log'),internal_limit_s=protocol['internal_limit_s'],external_limit_s=protocol['external_limit_s'],**DIAGNOSTIC))
    print(json.dumps(dict(phase=args.phase,index=args.index,exit_code=code,status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    return code


if __name__=='__main__':raise SystemExit(main())

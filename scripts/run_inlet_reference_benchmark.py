"""Exactly one maximum-inventory cost probe and two four-root diagnostic branches."""
from pathlib import Path
import argparse
import json
import math
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import inlet_reference_benchmark as model
import run_two_insertion_diagnostic as prior
import visibility_equilibrium as solver

ROOT=model.ROOT
THREADS=prior.THREADS
OLD_BASE=ROOT/'outputs/two_insertion_diagnostic_20261005/v1'
READBACK=ROOT/'outputs/two_insertion_diagnostic_20261005/readback_v1/summary.json'
READBACK_SHA='29ff03b3f4197fcbd2df2804acfbc7dcda9ea7c35fe347852ac84f4392694a1a'
COMPONENT=ROOT/'outputs/fixed_inventory_source_work_20261005/v1/receipt.json'
COST_PASS='PASS_MEASURED_968_REFERENCE_COST_NOT_PROCESS_ACCEPTANCE'
BRANCH_PASS='PASS_FOUR_STEP_REFERENCE_BRANCH_NOT_PROCESS_ACCEPTANCE'
require,write,read,sha,rel,verify=model.require,model.write,model.read,model.sha,model.rel,model.verify


def bindings():
    frozen=prior.source_bindings();references={}
    component=read(COMPONENT)
    require(component['exit_code']==0 and component['test_count']==6 and component['status']=='PASS_FIXED_INVENTORY_WORK_COMPONENT_ONLY',
        'fixed-inventory source component qualification absent')
    verify(component['source_sha256']);verify(component['artifact_sha256'])
    frozen.update(component['source_sha256']);frozen[rel(COMPONENT)]=sha(COMPONENT)
    require(sha(READBACK)==READBACK_SHA,'qualified two-birth readback changed')
    old=read(OLD_BASE/'diagnostic/manifest.json');sup=read(OLD_BASE/'diagnostic/supervision.json');back=read(READBACK)
    require(old['status']==prior.PASS and sup['exit_code']==0 and sup['before_after_hashes_match'] and
        sup['manifest_sha256']==sha(OLD_BASE/'diagnostic/manifest.json') and sup['log_sha256']==sha(OLD_BASE/'diagnostic/run.log'),
        'original two-birth diagnostic receipt absent')
    require(back['status']=='PASS_SAVED_TWO_INSERTION_DIAGNOSTIC_READBACK_NOT_JOINT_OR_PROCESS_ACCEPTANCE', 'independent original readback failed')
    verify(old['source_sha256']);verify(old['artifacts_sha256']);verify(back['source_and_input_sha256'])
    references.update(old['artifacts_sha256']);references.update(back['source_and_input_sha256'])
    for path in (OLD_BASE/'diagnostic/manifest.json',OLD_BASE/'diagnostic/supervision.json',OLD_BASE/'diagnostic/run.log',READBACK):
        references[rel(path)]=sha(path)
    for name in ('scripts/inlet_reference_benchmark.py','scripts/run_inlet_reference_benchmark.py',
        'scripts/compare_inlet_reference_benchmark.py','scripts/fixed_inventory_source_work.py',
        'docs/INLET_REFERENCE_BENCHMARK_CONTRACT_20261005.md','docs/TWO_INSERTION_REFERENCE_PROPOSAL_20261005.md'):
        frozen[name]=sha(ROOT/name)
    verify(frozen);verify(references);return frozen,references


def root_budget(e,t,rank):
    best=None
    for count in range(2,25):
        maps=min(100,1+8*count+16)+2
        forecast=30+1.5*4*(maps*e+max(count,3)*t+rank+10)
        if forecast<=450.:best=dict(maximum_iterations=count,root_map_allowance=maps,forecast_s=forecast)
    return dict(best or dict(maximum_iterations=0,root_map_allowance=None,forecast_s=None),useful_four_roots_fit=best is not None,
        full_map_wall_s=e,tangent_linear_wall_s=t,reference_rank_wall_s=rank,forecast_limit_s=450.,
        formula='30+1.5*4*((min(100,1+8*M+16)+2)*e+max(M,3)*t+rank+10)',root_count=4,fallback_roots=0)


def preflight(stage,base):
    frozen,refs=bindings()
    if stage=='cost_968':
        cost=read(OLD_BASE/'cost/manifest.json');sup=read(OLD_BASE/'cost/supervision.json')
        require(cost['status']==prior.COST_PASS and sup['exit_code']==0 and sup['manifest_sha256']==sha(OLD_BASE/'cost/manifest.json'),
            'maximum-inventory probe requires measured predecessor cost')
        forecast=30+1.5*(2061/2011)**3*(cost['full_map_wall_s']+cost['tangent_linear_wall_s']+cost['reference_rank_wall_s'])
        require(forecast<=450.,'probe forecast exceeds budget')
        budget=dict(forecast_s=forecast,forecast_limit_s=450.,basis='MEASURED_943_MAP_TANGENT_RANK_WITH_CUBIC_DOF_ESTIMATE')
        refs.update({rel(path):sha(path) for path in (OLD_BASE/'cost/manifest.json',OLD_BASE/'cost/supervision.json')})
    else:
        directory=base/'cost_968';cost=read(directory/'manifest.json');sup=read(directory/'supervision.json')
        require(cost['status']==COST_PASS and cost['source_sha256']==frozen and cost['dt_s']==model.H and
            cost['particle_count']==968 and sup['exit_code']==0 and sup['before_after_hashes_match'] and
            sup['manifest_sha256']==sha(directory/'manifest.json') and sup['log_sha256']==sha(directory/'run.log'),
            'matching maximum-inventory cost qualification absent')
        verify(cost['artifacts_sha256']);budget=cost['four_root_budget']
        require(budget==root_budget(cost['full_map_wall_s'],cost['tangent_linear_wall_s'],cost['reference_rank_wall_s']) and
            budget['useful_four_roots_fit'] and budget['forecast_s']<=450.,'four roots do not fit measured budget')
        refs.update({rel(path):sha(path) for path in directory.rglob('*') if path.is_file()})
        if stage=='preloaded_fine':
            prerequisite=base/'inlet_fine';meta=read(prerequisite/'manifest.json');external=read(prerequisite/'supervision.json')
            require(meta['status']==BRANCH_PASS and meta['temporal_screen']['passed'] and external['exit_code']==0 and
                external['manifest_sha256']==sha(prerequisite/'manifest.json') and external['before_after_hashes_match'],
                'preloaded branch requires successful inlet fine and temporal screen')
            verify(meta['artifacts_sha256']);refs.update({rel(path):sha(path) for path in prerequisite.rglob('*') if path.is_file()})
    return dict(stage=stage,source_sha256=frozen,reference_sha256=refs,budget=budget,study_scope=model.STUDY,
        accepted=False,restart_allowed=False,dt_s=model.H,physical_evaluations_in_preflight=0)


def initial_snapshot(pa,state,directory):
    path=directory/'INDEPENDENT_STUDY_INITIAL_NOT_RESTARTABLE.npz'
    require(not path.exists(),'study initial output exists');model.old.core.checkpoint(path,pa,state)
    saved=model.old.domain.load(str(path));actual=saved['arrays']['solid'];actualstate=saved['solver_data']
    require(model.old.clean(actualstate)==model.old.clean(state) and
        model.old.inlet.inventory._array_digest(actual)==model.old.inlet.inventory._array_digest(pa),'study initial readback differs')
    require(model.old.inlet.restore_acceptance_state(json.dumps(actualstate['inlet_acceptance']),actual)==state['inlet_acceptance'],
        'initial V2 registry not bound to actual PA')
    state['parent_checkpoint_sha256']=sha(path);state['parent_checkpoint_classification']=state['classification']
    return path


def worker(stage,base):
    directory=base/stage;protocol=read(directory/'protocol.json');started=time.perf_counter();trial=None;rows=[];calculations=[];code=2
    report=dict(status='RUNNING_INDEPENDENT_REFERENCE_DIAGNOSTIC',stage=stage,study_scope=model.STUDY,dt_s=model.H,
        source_sha256=protocol['source_sha256'],reference_sha256=protocol['reference_sha256'],budget=protocol['budget'],
        accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.,diagnostic_steps=0,
        joint_process_inlet_qualification='NOT_GRANTED',control_release='NOT_EXERCISED',continued_windows=0)
    def timeout():
        try:write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NOT_ACCEPTED',completed_diagnostic_steps=len(rows),wall_seconds=time.perf_counter()-started))
        finally:os._exit(124)
    timer=threading.Timer(540.,timeout);timer.daemon=True;timer.start()
    try:
        require(all(os.environ.get(k)=='1' for k in THREADS),'single-thread cost environment required')
        verify(protocol['source_sha256']);verify(protocol['reference_sha256']);require(protocol==preflight(stage,base),'worker preflight changed')
        branch='preloaded_fine' if stage=='cost_968' else stage
        parent,state=model.initialize(branch,protocol['source_sha256']);initial=initial_snapshot(parent,state,directory)
        report.update(initial_checkpoint=rel(initial),initial_checkpoint_sha256=sha(initial),initial_mass_kg=state['study_initial_mass_kg'],
            initial_sensible_heat_J=state['study_initial_sensible_heat_J'],initial_particle_count=len(parent.x))
        if stage=='cost_968':
            prepared,inner,event,storage,geometry=model.prepare(parent,state,1);trial=model.old.InletTrial(prepared,inner,state['parameters'])
            begin=time.perf_counter();r=trial.evaluate(model.old.batch.guess_for(trial,model.H),model.H);e=time.perf_counter()-begin
            model.old.save_map(directory,trial,r)
            begin=time.perf_counter();direction,linear=model.old.core.newton_direction(trial,r);t=time.perf_counter()-begin
            np.savez_compressed(directory/'actual_search_direction.npz',particle_id=prepared.particle_id,direction=direction)
            write(directory/'actual_search_linear.json',linear);rank=model.old.reference_rank(trial,r,directory)
            report.update(status=COST_PASS,particle_count=968,free_dof=int(inner['free'].sum()),full_map_wall_s=e,
                tangent_linear_wall_s=t,reference_rank_wall_s=rank['wall_seconds'],four_root_budget=root_budget(e,t,rank['wall_seconds']),
                actual_scaled_force_residual=r['scaled_residual'],root_solved=False)
        else:
            deadline=started+525.;maximum=protocol['budget']['maximum_iterations']
            for step in range(1,5):
                verify(protocol['source_sha256']);verify(protocol['reference_sha256']);sd=directory/('step_'+str(step).zfill(3));sd.mkdir()
                prepared,inner,event,storage,geometry=model.prepare(parent,state,step)
                write(sd/'source_event.json',event);write(sd/'prepared_geometry.json',geometry)
                trial=model.old.InletTrial(prepared,inner,state['parameters'])
                result,metrics=solver.solve(trial,model.H,model.old.batch.guess_for(trial,model.H),sd/'root',deadline,maximum)
                reserve=1.5*(2*protocol['budget']['full_map_wall_s']+protocol['budget']['reference_rank_wall_s']+10)
                require(time.perf_counter()+reserve<deadline,'fresh/scalar/rank/source transaction lacks measured budget')
                actual,nextstate,fresh,calculated=model.qualified_step(parent,trial,result,event,storage,geometry,sd,protocol['source_sha256'])
                calculations.append(calculated);write(sd/'cumulative_open_energy.json',prior.cumulative_source_ledger(calculations,nextstate))
                expected_births=(1 if step<=2 else 2) if stage=='inlet_fine' else 0
                require(len(nextstate['inlet_acceptance']['source_events'])==expected_births and
                    nextstate['inlet_acceptance']['mechanical_step_count']==step,'birth/mechanical history count changed')
                require(nextstate['inlet_acceptance']['scope_lock']==model.old.SCOPE,'diagnostic scope lock changed')
                previous=state['parent_checkpoint_sha256'];path=sd/'DIAGNOSTIC_NOT_RESTARTABLE.npz'
                parent,state=model.old.write_and_reload(actual,nextstate,path)
                row=dict(step=step,time_s=state['time_s'],dt_s=model.H,particle_count=len(parent.x),parent_checkpoint_sha256=previous,
                    checkpoint=rel(path),checkpoint_sha256=sha(path),source_kind='FIXED_INVENTORY' if event is None else 'BIRTH',
                    graph_fingerprint=fresh['graph_fingerprint'],metrics=metrics,accepted=False,restart_allowed=False)
                write(sd/'transaction.json',row);rows.append(row);report['diagnostic_steps']=step
            require(state['time_s']==4*model.H and len(parent.x)==968,'unexpected branch endpoint')
            report.update(status=BRANCH_PASS,diagnostic_time_s=state['time_s'],particle_count=len(parent.x),
                final_checkpoint_sha256=state['parent_checkpoint_sha256'],transactions=rows,
                cumulative_source_ledger=prior.cumulative_source_ledger(calculations,state),final_bore_observation=state['bore_observation'])
            if stage=='inlet_fine':
                from compare_inlet_reference_benchmark import temporal_screen
                screen=temporal_screen(base);write(directory/'temporal_screen.json',screen);report['temporal_screen']=screen
                require(screen['passed'],'inlet fine versus original coarse temporal screen failed')
        verify(protocol['source_sha256']);verify(protocol['reference_sha256']);code=0
    except Exception as error:
        report.update(status='FAIL_REFERENCE_'+stage.upper(),error=repr(error),traceback=traceback.format_exc())
        if trial is not None and trial.last_complete_result is not None:
            try:
                last=directory/'LAST_FULL_UNCOMMITTED';last.mkdir();model.old.save_map(last,trial,trial.last_complete_result)
            except Exception as saving:report['last_candidate_save_error']=repr(saving)
    finally:
        timer.cancel();report.update(wall_seconds=time.perf_counter()-started,completed_transactions=rows,
            last_trial_calls=[] if trial is None else trial.calls,peak_working_set_bytes=model.old.core.memory_bytes())
        report['artifacts_sha256']={rel(path):sha(path) for path in directory.rglob('*') if path.is_file() and
            path.name not in ('manifest.json','supervision.json','run.log')};write(directory/'manifest.json',report)
    return code


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=('cost_968','inlet_fine','preloaded_fine','compare'),required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--preflight-only',action='store_true');parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();base=Path(args.output)
    if not base.is_absolute():base=ROOT/base
    base=base.resolve();require(base.is_relative_to(ROOT/'outputs'),'new output must remain under project outputs')
    if args.stage=='compare':
        from compare_inlet_reference_benchmark import compare
        require(not args.worker and not args.preflight_only,'compare is saved-only and cannot spawn physics')
        compare(base);return
    if args.worker:raise SystemExit(worker(args.stage,base))
    protocol=preflight(args.stage,base)
    if args.preflight_only:print(json.dumps(dict(status='PASS_READ_ONLY_PREFLIGHT',**protocol)));return
    directory=base/args.stage;require(not directory.exists(),'immutable stage exists');directory.mkdir(parents=True)
    write(directory/'protocol.json',protocol)
    command=[str(ROOT/'.envs/pysph/python.exe'),str(Path(__file__).resolve()),'--stage',args.stage,'--output',str(base),'--worker']
    env=dict(os.environ,**{key:'1' for key in THREADS});started=time.perf_counter();code=124
    with (directory/'run.log').open('x',encoding='utf-8') as log:
        try:code=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600.).returncode
        except subprocess.TimeoutExpired:log.write('\nEXTERNAL_600_SECOND_TIMEOUT_NOT_ACCEPTED\n')
    unchanged=True
    try:verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
    except Exception:unchanged=False
    manifest=directory/'manifest.json';write(directory/'supervision.json',dict(command=command,cwd=str(ROOT),exit_code=code,
        wall_seconds=time.perf_counter()-started,before_after_hashes_match=unchanged,manifest_sha256=sha(manifest) if manifest.exists() else None,
        log_sha256=sha(directory/'run.log'),worker_limit_s=540.,external_limit_s=600.,accepted=False,restart_allowed=False))
    print(json.dumps(dict(output=rel(directory),exit_code=code,status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    if code or not unchanged:raise SystemExit(code or 2)


if __name__=='__main__':main()

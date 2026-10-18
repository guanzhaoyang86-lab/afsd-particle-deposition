"""Supervised cost then two true coupled inlet steps; NEVER process acceptance."""
from pathlib import Path
import argparse
import math
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import two_insertion_diagnostic as model
import visibility_equilibrium as solver

ROOT=model.ROOT
THREADS=('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')
require,sha,read,write,verify,rel=model.require,model.sha,model.read,model.write,model.verify,model.rel
CORE=model.core
COST_PASS='PASS_MEASURED_INLET_MAP_TANGENT_RANK_COST_NOT_PHYSICAL_VALIDATION'
PASS='PASS_TWO_INSERTION_DIAGNOSTIC_NOT_JOINT_OR_PROCESS_ACCEPTANCE'


def source_bindings():
    mapping={}
    def merge(values):
        verify(values)
        for name,digest in values.items():
            require(name not in mapping or mapping[name]==digest,'conflicting source binding '+name)
            mapping[name]=digest
    merge(CORE.frozen_sources());merge(model.domain.source_hashes())
    receipts=[
        ('outputs/inlet_acceptance_v2_20261004/component_receipt.json','PASS_INLET_ACCEPTANCE_V2_COMPONENT_ONLY',
         '6f730d4f2d165db6b8c488696fe3379bcf5f995931f3aed242a43d91895b66a5'),
        ('outputs/inlet_source_work_20261005/component_receipt.json','PASS_INLET_SOURCE_WORK_COMPONENT_ONLY',
         '7a91f63d421a7b3893d8c5a8164862c545618e70e691d20a613c9df6efa483fb'),
        ('outputs/visible_path_tight_20261004/qualification_v1/summary.json',
         'PASS_VISIBLE_PATH_TIGHT_COMPONENTS_AND_SAVED_REPLAY_NOT_PHYSICAL_VALIDATION',
         '3d461593a2fab30900a83bcde5b2835b353bec3d66783393b3500ce493b750bc'),
        ('outputs/visible_top_components_20261004/v1/summary.json','PASS_TOP_CONTACT_COMPONENTS_NOT_PROCESS_VALIDATION',
         '7fbb5fdca371335ff528f4132fb840948dc3239de8eb50e420b3533014505cab'),
        ('outputs/visible_graph_20261004/operators_v1/summary.json','PASS_REGISTERED_COMPONENT_BENCHMARKS_ONLY',None),
        ('outputs/expanded_hot_qualification_20261005/v1/summary.json',
         'PASS_EXPANDED_INITIAL_ELASTIC_RANK_AND_SAME_MAP_NOT_PHYSICAL_VALIDATION',None)]
    for name,status,digest in receipts:
        path=ROOT/name;receipt=read(path)
        require(receipt['status']==status and (digest is None or sha(path)==digest),'component qualification changed: '+name)
        if name.endswith('component_receipt.json'):
            require(receipt['exit_code']==0 and receipt['test_count']==11,'component receipt failed')
        for key in ('source_sha256','source_sha256_after','source_and_input_sha256','artifact_sha256','output_sha256'):
            if isinstance(receipt.get(key),dict):merge(receipt[key])
        mapping[name]=sha(path)
    init=read(model.INITIAL.parent/'summary.json')
    require(init['status']=='PASS_INITIAL_ARRAYS_AND_GEOMETRY_ONLY_NOT_PHYSICAL_QUALIFICATION' and
        init['initialization_sha256']==model.INITIAL_SHA and sha(model.INITIAL)==model.INITIAL_SHA,'initial_v2 receipt differs')
    merge(init['source_sha256'])
    names=['scripts/two_insertion_diagnostic.py','scripts/run_two_insertion_diagnostic.py',
        'scripts/inlet_acceptance_v2.py','scripts/inlet_source_work.py','scripts/visibility_equilibrium.py',
        'scripts/secant_equilibrium.py','scripts/run_batched_process.py','scripts/run_discrete_post_start.py',
        'scripts/run_visible_support_process.py','scripts/run_first_shoulder_event.py',
        'scripts/contact_shell_audit.py','scripts/contact_shell_enclosure.py','scripts/contact_shell_finite_enclosure.py',
        'scripts/visible_support_path.py','scripts/visible_support_path_tight.py','scripts/visible_support_numerics.py']
    files=[ROOT/name for name in names]+[model.CONTRACT,model.INITIAL,model.INITIAL.parent/'summary.json']
    merge({rel(path):sha(path) for path in files});return mapping


def budget_for_two_roots(e,t,rank):
    # Second N=968 is not yet measured: cubic free-DOF growth is an explicit estimate.
    growth=(2061/2011)**3;best=None
    for updates in range(2,25):
        # Includes up to two refresh-related extra line searches, capped at solver's 100 maps.
        maps=min(100,1+8*updates+16)+2
        forecast=30+1.5*growth*2*(maps*e+max(updates,3)*t+rank+10)
        if forecast<=450.:
            best=dict(maximum_iterations=updates,forecast_s=forecast,root_map_allowance=maps)
    return dict(best or dict(maximum_iterations=0,forecast_s=None,root_map_allowance=None),
        useful_two_roots_fit=best is not None,forecast_limit_s=450.,second_inventory_cost_measured=False,
        particle_count_first=943,particle_count_second=968,free_dof_first=2011,free_dof_second=2061,
        second_inventory_cubic_dof_growth=growth,measured_full_map_s=e,measured_tangent_linear_s=t,
        measured_reference_rank_s=rank,
        formula='30+1.5*(2061/2011)^3*2*((min(100,1+8*M+16)+2)*e+max(M,3)*t+rank+10)',
        fallback_roots_allowed=0,acceptance=False)


def preflight(stage,base):
    frozen=source_bindings();references={}
    if stage=='cost':
        old=read(ROOT/'outputs/expanded_hot_qualification_20261005/v1/summary.json')
        # Qualification actually measured three maps, two tangents and the same-size reference eigensolve.
        measured=float(old['initial_elastic_rank']['wall_seconds'])+float(old['initial_same_map']['wall_seconds'])
        forecast=30+1.5*(2011/1961)**3*measured
        require(forecast<=450.,'cost-stage forecast exceeds 450 s')
        budget=dict(forecast_s=forecast,forecast_limit_s=450.,source='OLD_918_QUALIFICATION_MEASURED_TOTAL_TIMES_CUBIC_DOF_ESTIMATE',
            prior_qualification_wall_s=measured,not_new_inventory_measurement=True)
    else:
        directory=base/'cost';manifest=read(directory/'manifest.json');supervision=read(directory/'supervision.json')
        require(manifest['status']==COST_PASS and manifest['source_sha256']==frozen and
            manifest['initialization_sha256']==model.INITIAL_SHA and manifest['dt_s']==model.DT,
            'matching measured inlet cost missing')
        require(supervision['exit_code']==0 and supervision['before_after_hashes_match'] and
            supervision['manifest_sha256']==sha(directory/'manifest.json') and supervision['log_sha256']==sha(directory/'run.log'),
            'cost external receipt mismatch')
        verify(manifest['artifacts_sha256']);budget=manifest['two_root_budget']
        require(budget==budget_for_two_roots(manifest['full_map_wall_s'],manifest['tangent_linear_wall_s'],
            manifest['reference_rank_wall_s']) and budget['useful_two_roots_fit'] and budget['forecast_s']<=450.,
            'two-root forecast is absent or exceeds registered budget')
        references={rel(path):sha(path) for path in directory.rglob('*') if path.is_file()}
    return dict(source_sha256=frozen,reference_sha256=references,budget=budget,stage=stage,
        classification=model.CLASSIFICATION,model_scope=model.SCOPE,dt_s=model.DT,initialization_sha256=model.INITIAL_SHA,
        accepted=False,restart_allowed=False,physical_evaluations_in_preflight=0)


def cumulative_source_ledger(rows,state):
    require(len(rows)==state['step']==len(state['inlet_acceptance']['source_work_events']),'source work prefix missing')
    # Each term below is one actual step; inventory changes telescope across births.
    sums={key:math.fsum(row['open_energy'][key] for row in rows) for key in
        ('total_storage_increment_J','tool_work_J','plastic_unthermalized_J','wall_heat_out_J',
         'mechanical_raw_residual_J','thermal_external_raw_residual_J','open_raw_residual_J')}
    sums.update({key:math.fsum(row['computed_terms'][key]['value'] for row in rows) for key in model.inlet.WORK_TERMS})
    heat=math.fsum(row['source_increment']['sensible_heat_J'] for row in rows)
    recomputed=(sums['mechanical_reaction_work']+sums['tool_work_J']+sums['elastic_convected_energy']+heat-
        sums['total_storage_increment_J']-sums['plastic_unthermalized_J']-sums['wall_heat_out_J'])
    bounds=math.fsum(row['open_energy']['algebra_roundoff_bound_J'] for row in rows)
    require(abs(recomputed-sums['open_raw_residual_J'])<=bounds,'cumulative open-energy algebra mismatch')
    for previous,current in zip(rows,rows[1:]):
        for name in ('material','sensible','contact'):
            a=previous['open_energy'][name+'_storage_after_J'];b=current['open_energy'][name+'_storage_before_J']
            require(abs(a-b)<=512*np.finfo(float).eps*max(abs(a)+abs(b),np.finfo(float).tiny),
                'source energy inventory not continuous between births: '+name)
    return dict(status='COMPUTED_CUMULATIVE_OPEN_INVENTORY_DIAGNOSTIC_NOT_PROCESS_ACCEPTANCE',
        steps=len(rows),sum=sums,source_sensible_heat_J=heat,
        source_registry=state['inlet_acceptance']['protocol_registry']['cumulative'],
        raw_absolute_sum_J=math.fsum(abs(row['open_energy']['open_raw_residual_J']) for row in rows),
        recomputed_open_raw_J=recomputed,identity_roundoff_bound_J=bounds,
        kinetic_and_momentum_diagnostics_not_in_qs_energy=True,reference_temperature_K=state['parameters']['room'],
        retained_mechanical_interval_energy=state['interval_energy'])


def worker(stage,base):
    directory=base/stage;protocol=read(directory/'protocol.json');started=time.perf_counter()
    trial=None;state=None;rows=[];source_rows=[];code=2
    report=dict(status='RUNNING_UNACCEPTED_INLET_DIAGNOSTIC',stage=stage,classification=model.CLASSIFICATION,
        model_scope=model.SCOPE,source_sha256=protocol['source_sha256'],reference_sha256=protocol['reference_sha256'],
        initialization_sha256=model.INITIAL_SHA,dt_s=model.DT,accepted=False,restart_allowed=False,
        accepted_steps=0,accepted_time_advanced_s=0.,diagnostic_steps=0,diagnostic_time_s=0.,
        joint_inlet_benchmark='NOT_RUN',preloaded_long_bar_comparison='NOT_RUN',temporal_convergence='NOT_RUN',
        continuous_mass_flow_rate='NOT_MODELED',control_release='NOT_EXERCISED',budget=protocol['budget'])
    def timeout():
        try:write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NOT_ACCEPTED',wall_seconds=time.perf_counter()-started,
            accepted=False,restart_allowed=False,completed_diagnostic_steps=len(rows)))
        finally:os._exit(124)
    timer=threading.Timer(540.,timeout);timer.daemon=True;timer.start()
    try:
        require(all(os.environ.get(key)=='1' for key in THREADS),'single-thread pinned budget environment required')
        verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
        require(protocol==preflight(stage,base),'worker preflight changed')
        parent,state=model.initial_state(protocol['source_sha256'])
        initial_digest=model.inlet.inventory._array_digest(parent)
        if stage=='cost':
            prepared,inner,event,storage,birth=model.prepare_birth(parent,state,1)
            trial=model.InletTrial(prepared,inner,state['parameters'])
            write(directory/'birth_geometry.json',birth);write(directory/'source_event.json',event)
            guess=model.batch.guess_for(trial,model.DT)
            begin=time.perf_counter();r=trial.evaluate(guess,model.DT);elapsed=time.perf_counter()-begin
            model.save_map(directory,trial,r)
            begin=time.perf_counter();direction,linear=CORE.newton_direction(trial,r);tangent=time.perf_counter()-begin
            np.savez_compressed(directory/'actual_tangent_direction.npz',particle_id=prepared.particle_id,
                direction_m=direction,free=inner['free'],prescribed_guess_m=guess)
            write(directory/'actual_tangent_linear.json',linear)
            rank=model.reference_rank(trial,r,directory)
            budget=budget_for_two_roots(elapsed,tangent,rank['wall_seconds'])
            require(model.inlet.inventory._array_digest(parent)==initial_digest,'cost map mutated independent parent')
            report.update(status=COST_PASS,full_map_wall_s=elapsed,tangent_linear_wall_s=tangent,
                reference_rank_wall_s=rank['wall_seconds'],two_root_budget=budget,actual_free_dof=int(inner['free'].sum()),
                particle_count=len(prepared.x),cost_candidate_scaled_force_residual=r['scaled_residual'],
                root_solved=False,source_work='NOT_COMPUTED_COST_ONLY')
        else:
            deadline=started+525.;maximum=protocol['budget']['maximum_iterations']
            for index in (1,2):
                verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
                stepdir=directory/('step_'+str(index).zfill(3));stepdir.mkdir()
                prepared,inner,event,storage,birth=model.prepare_birth(parent,state,index)
                write(stepdir/'source_event.json',event);write(stepdir/'birth_geometry.json',birth)
                trial=model.InletTrial(prepared,inner,state['parameters'])
                rootdir=stepdir/'root';rootdir.mkdir()
                result,metrics=solver.solve(trial,model.DT,model.batch.guess_for(trial,model.DT),rootdir,deadline,maximum)
                reserve=2*protocol['budget']['measured_full_map_s']+protocol['budget']['measured_reference_rank_s']+10.
                require(time.perf_counter()+1.5*protocol['budget']['second_inventory_cubic_dof_growth']*reserve<deadline,
                    'insufficient measured budget for fresh/scalar/rank/source transaction')
                actual,nextstate,fresh,calculated=model.qualify_candidate(parent,trial,result,event,storage,birth,
                    stepdir,protocol['source_sha256'])
                require(nextstate['inlet_acceptance']['mechanical_step_count']==index and
                    len(nextstate['inlet_acceptance']['source_events'])==index,'source/physical step sequence mismatch')
                require(not nextstate['inlet_acceptance']['control_events'][-1]['released_ids'],
                    'control release outside this two-step diagnostic scope')
                source_rows.append(calculated)
                write(stepdir/'cumulative_open_energy.json',cumulative_source_ledger(source_rows,nextstate))
                checkpoint=stepdir/'DIAGNOSTIC_NOT_RESTARTABLE.npz'
                previous_sha=state['parent_checkpoint_sha256']
                parent,state=model.write_and_reload(actual,nextstate,checkpoint)
                row=dict(step=index,time_s=state['time_s'],dt_s=model.DT,parent_checkpoint_sha256=previous_sha,
                    checkpoint=rel(checkpoint),checkpoint_sha256=sha(checkpoint),particle_count=len(parent.x),
                    receipt_sha256=sha(stepdir/'mechanical_receipt.json'),source_work_sha256=sha(stepdir/'source_work.json'),
                    source_event_sha256=sha(stepdir/'source_event.json'),graph_fingerprint=fresh['graph_fingerprint'],metrics=metrics,
                    accepted=False,restart_allowed=False,diagnostic_only=True)
                rows.append(row);write(stepdir/'transaction.json',row)
                report.update(diagnostic_steps=index,diagnostic_time_s=state['time_s'])
            require(state['time_s']==2*model.DT and len(parent.x)==968,'registered two-step endpoint differs')
            require(state['inlet_acceptance']['protocol_registry']['issued_particle_count']==50 and
                state['inlet_acceptance']['scope_lock']==model.SCOPE,'birth source count/scope changed')
            report.update(status=PASS,particle_count=len(parent.x),diagnostic_end_checkpoint_sha256=state['parent_checkpoint_sha256'],
                final_bore_observation=state['inlet_acceptance']['bore_observation'],
                cumulative_source_ledger=cumulative_source_ledger(source_rows,state),diagnostic_transactions=rows)
        verify(protocol['source_sha256']);verify(protocol['reference_sha256']);code=0
    except Exception as error:
        report.update(status='FAIL_UNACCEPTED_TWO_INSERTION_'+stage.upper(),error=repr(error),traceback=traceback.format_exc())
        if trial is not None and trial.last_complete_result is not None:
            try:
                evidence=directory/'LAST_FULL_UNCOMMITTED_EVALUATION';evidence.mkdir()
                model.save_map(evidence,trial,trial.last_complete_result)
            except Exception as saving:report['last_candidate_save_error']=repr(saving)
    finally:
        timer.cancel();report.update(wall_seconds=time.perf_counter()-started,completed_diagnostic_transactions=rows,
            last_trial_calls=[] if trial is None else trial.calls,peak_working_set_bytes=CORE.memory_bytes())
        report['artifacts_sha256']={rel(path):sha(path) for path in directory.rglob('*')
            if path.is_file() and path.name not in ('run.log','manifest.json','supervision.json')}
        write(directory/'manifest.json',report)
    return code


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=('cost','diagnostic'),required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS);args=parser.parse_args()
    base=Path(args.output)
    if not base.is_absolute():base=ROOT/base
    base=base.resolve();require(base.is_relative_to(ROOT/'outputs'),'new diagnostic output must remain under project outputs')
    if args.worker:raise SystemExit(worker(args.stage,base))
    protocol=preflight(args.stage,base)
    if args.preflight_only:
        print(__import__('json').dumps(dict(status='PASS_READ_ONLY_PREFLIGHT',**protocol)));return
    directory=base/args.stage;require(not directory.exists(),'immutable stage exists; choose a new output version')
    directory.mkdir(parents=True,exist_ok=False);write(directory/'protocol.json',protocol)
    command=[str(ROOT/'.envs/pysph/python.exe'),str(Path(__file__).resolve()),'--stage',args.stage,'--output',str(base),'--worker']
    env=dict(os.environ,**{key:'1' for key in THREADS});started=time.perf_counter();code=124
    with (directory/'run.log').open('x',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600.).returncode
        except subprocess.TimeoutExpired:stream.write('\nEXTERNAL_600_SECOND_TIMEOUT_NOT_ACCEPTED\n')
    unchanged=True
    try:verify(protocol['source_sha256']);verify(protocol['reference_sha256'])
    except Exception:unchanged=False
    manifest=directory/'manifest.json'
    write(directory/'supervision.json',dict(command=command,cwd=str(ROOT),exit_code=code,wall_seconds=time.perf_counter()-started,
        before_after_hashes_match=unchanged,log_sha256=sha(directory/'run.log'),manifest_sha256=sha(manifest) if manifest.exists() else None,
        worker_limit_s=540.,external_limit_s=600.,accepted=False,restart_allowed=False))
    print(__import__('json').dumps(dict(stage=args.stage,output=rel(directory),exit_code=code,
        status=read(manifest)['status'] if manifest.exists() else 'NO_COMPLETE_MANIFEST',accepted=False,restart_allowed=False)))
    if code or not unchanged:raise SystemExit(code or 2)


if __name__=='__main__':main()

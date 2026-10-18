"""Scheduled independent 918-particle cost/root diagnostics, never acceptance."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import math
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import expanded_hot_domain as domain
import run_batched_process as batch
import run_discrete_contact_process as core
import run_discrete_post_start as post
from afsd_backend.pysph.interval_acceptance import assess_interval_energy,observe_square_bore_crossings
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from contact_shell_finite_enclosure import audit_contact_shell_finite_enclosure
from contact_shell_audit import CERTIFIED as SHELL_CERTIFIED
from visible_support_path_tight import audit_visible_support_path,CERTIFIED

ROOT=domain.ROOT
INITIAL=ROOT/'outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2'
CONTRACT=ROOT/'docs/EXPANDED_HOT_DOMAIN_COST_CONTRACT_20261004.md'
TIGHT_QUAL=ROOT/'outputs/visible_path_tight_20261004/qualification_v1/summary.json'
BASE=ROOT/'outputs/expanded_hot_domain_cost_20261004'
DT_CHOICES={'15.625':15.625e-6,'7.8125':7.8125e-6}
THREAD_KEYS=('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')


def require(value,message):
    if not bool(value):raise ValueError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def rel(path):return Path(path).relative_to(ROOT).as_posix()
def write(path,value):
    with Path(path).open('x',encoding='utf-8') as stream:json.dump(core.clean(value),stream,indent=2,allow_nan=False);stream.write('\n')
def verify(mapping):
    for name,digest in mapping.items():require(sha(ROOT/name)==digest,'source/input drift: '+name)


class RecordedTrial:
    """Instrumentation proxy around the actual factory, without altered maps."""
    def __init__(self,delegate):
        self.delegate=delegate;self.calls=[];self.last_complete_result=None
    def __getattr__(self,key):return getattr(self.delegate,key)
    def evaluate(self,displacement,dt):
        begin=time.perf_counter();row=dict(kind='complete_map',dt_s=dt,completed=False)
        self.calls.append(row)
        try:
            r=self.delegate.evaluate(displacement,dt);self.last_complete_result=r
            row.update(completed=True,force_residual_N=r['force_residual_N'],edge_count=len(r['i']))
            return r
        except Exception as error:row['error']=repr(error);raise
        finally:row['wall_seconds']=time.perf_counter()-begin
    def tangent(self,result):
        begin=time.perf_counter();row=dict(kind='frozen_tangent',completed=False);self.calls.append(row)
        try:
            tangent=self.delegate.tangent(result);row.update(completed=True,shape=list(tangent.shape))
            return tangent
        except Exception as error:row['error']=repr(error);raise
        finally:row['wall_seconds']=time.perf_counter()-begin


def baseline_bindings():
    initial=read(INITIAL/'summary.json')
    require(initial['status']=='PASS_INITIAL_ARRAYS_AND_GEOMETRY_ONLY_NOT_PHYSICAL_QUALIFICATION',
        'independent initialization receipt absent')
    verify(initial['source_sha256'])
    require(initial['particle_count']==918 and initial['free_displacement_dof']==1961,'wrong independent population')
    snapshot=ROOT/initial['initialization_checkpoint']
    require(sha(snapshot)==initial['initialization_sha256'],'independent initial snapshot changed')
    qualification=read(TIGHT_QUAL)
    require(qualification['qualification_passed'] and qualification['status']==
        'PASS_VISIBLE_PATH_TIGHT_COMPONENTS_AND_SAVED_REPLAY_NOT_PHYSICAL_VALIDATION','tight certificate qualification absent')
    verify(qualification['source_and_input_sha256'])
    names=['scripts/probe_expanded_hot_domain_cost.py','scripts/probe_expanded_hot_domain_cost_v2.py','scripts/probe_expanded_hot_domain_cost_v3.py',
        'docs/EXPANDED_HOT_DOMAIN_STARTUP_V3_ADDENDUM_20261005.md',
        'docs/EXPANDED_HOT_DOMAIN_COST_V2_ADDENDUM_20261005.md','scripts/run_batched_process.py',
        'scripts/run_discrete_post_start.py','scripts/contact_shell_finite_enclosure.py',
        'scripts/contact_shell_audit.py','scripts/visible_support_path_tight.py']
    files=[ROOT/name for name in names]+[CONTRACT,TIGHT_QUAL,snapshot,INITIAL/'summary.json',
        INITIAL/'resolved_parameters.json',INITIAL/'initialization.json',INITIAL/'initial_particle_fields.npz']
    # Bind old cost and failed energy evidence; no old result is replaced.
    for dt_key in ('62_5','31_25'):
        prior=BASE/'v2'/('dt_'+dt_key+'us')
        for name in ('cost/manifest.json','cost/supervision.json','root/manifest.json','root/supervision.json',
                     'root/fresh_basic_gates.json','root/fresh/ledger_and_geometry.json'):
            files.append(prior/name)
    # All transitive backend sources and the old solve/journal helpers are bound.
    mapping={}
    for sources in (core.frozen_sources(),domain.source_hashes(),qualification['source_and_input_sha256']):
        for name,digest in sources.items():
            require(name not in mapping or mapping[name]==digest,'conflicting dependency SHA: '+name)
            mapping[name]=digest
    mapping.update({rel(p):sha(p) for p in files});verify(mapping)
    return mapping


def load_independent_initial():
    initial=read(INITIAL/'summary.json');snapshot=ROOT/initial['initialization_checkpoint']
    require(sha(snapshot)==initial['initialization_sha256'],'independent snapshot hash differs')
    restored=domain.load(str(snapshot));pa=restored['arrays']['solid'];state=restored['solver_data']
    p,provenance=domain.parameters(.008)
    require(state['parameters']==p and state['parameter_provenance']==provenance,'independent parameters/provenance differ')
    require(state['operator_variant']=='visible' and state['time_s']==0. and state['step']==0 and
        state['restart_allowed'] is False and state['accepted'] is False,'diagnostic must start at independent time zero')
    require('spatial_model_transition' not in state and 'branch_origin_sha256' not in state,'fabricated old parent is forbidden')
    audit=domain.initial_geometry_and_field_audit(pa,state)
    require(audit['particle_count']==918 and audit['free_displacement_dof']==1961 and
        abs(audit['initial_energy']['sensible_heat_J']-119.902275)<1e-10,'independent initial audit differs')
    return pa,state,p,provenance,audit,snapshot


def root_budget(e,t):
    maximum=0;forecast=None
    for count in range(1,25):
        candidate=30+1.5*((1+8*count+2)*e+count*t)
        if candidate<=450:maximum=count;forecast=candidate
    return dict(maximum_iterations=maximum,forecast_s=forecast,useful_root_fits=maximum>=2,
        formula='30 + 1.5*((1+8*M+2)*full_map_wall_s + M*tangent_linear_wall_s)',
        forecast_limit_s=450.,fresh_and_scalar_maps_reserved=2)


def save_full_map(directory,trial,r,label):
    directory.mkdir(parents=True,exist_ok=False)
    core.journal(directory,trial,r,label)
    np.savez_compressed(directory/'material_and_operators.npz',accepted=False,restart_allowed=False,
        particle_id=trial.pa.particle_id,component_id=trial.pa.component_id,mass_kg=trial.pa.m,
        free_axes=trial.state['free'],velocity_m_s=r['displacement']/r['dt_s'],
        gradient_per_s=r['gradient'],correction=r['correction'],volume_m3=r['volume'],
        stress_Pa=r['material'][0],plastic_strain=r['material'][1],plastic_work_density_J_m3=r['material'][2],
        return_dissipation_density_J_m3=r['material'][3],effective_rate_per_s=r['material'][4],
        material_temperature_K=r['material_temperature'],plastic_strain_rate_per_s=(r['ep']-trial.pa.ep)/r['dt_s'])
    np.savez_compressed(directory/'visible_graph.npz',accepted=False,restart_allowed=False,
        particle_id=trial.pa.particle_id,position0=trial.x0,position1=r['position'],
        i=r['i'],j=r['j'],g=r['g'],raw_i=r['raw_i'],raw_j=r['raw_j'],raw_g=r['raw_g'],
        thermal_pair_conductance_W_K=r['thermal']['pair_conductance_W_K'])
    energy=core.energy_diagnostics(trial,r)
    write(directory/'ledger_and_geometry.json',dict(accepted=False,restart_allowed=False,
        mechanical=energy,thermal=r['thermal'],interval_energy_assessment=assess_interval_energy(trial.state['interval_energy'],energy),
        visibility=r['visibility_graph']['diagnostics'],visibility_numerical_domain=r['visibility_numerical_domain'],
        graph_fingerprint=r['graph_fingerprint'],max_correction_condition=float(np.linalg.cond(r['correction']).max()),
        tool_pose0=core.pose(trial.p,0.),tool_pose1=core.pose(trial.p,r['dt_s']),
        initial_energy=trial.state['initial_array_geometry_audit']['initial_energy']))


def check_fresh_root(trial,result,directory,deadline,map_wall):
    require(time.perf_counter()+3*map_wall+15<deadline,'insufficient measured budget for fresh/scalar gates')
    fresh=trial.evaluate(result['displacement'],result['dt_s']);p=trial.p;pa=trial.pa;dt=fresh['dt_s']
    save_full_map(directory/'fresh',trial,fresh,'FRESH_ROOT_NOT_ACCEPTED')
    energy=core.energy_diagnostics(trial,fresh);interval=assess_interval_energy(trial.state['interval_energy'],energy)
    gates=dict(fresh_force=fresh['scaled_residual']<=1 and float(np.max(np.abs(fresh['force'][trial.state['free']])))<=fresh['force_limit_N'],
        contact_sweep=bool(fresh['contact']['sweep'].get('certified',False)),
        raw_energy=abs(energy['global_raw_residual_J'])<=energy['registered_energy_limit_J'],
        interval_energy=bool(interval['accepted']),penetration=float(np.max(fresh['contact']['penetration_m']))<=.1*p['dx'],
        identical_fresh_graph=result['graph_fingerprint']==fresh['graph_fingerprint'] and
            np.array_equal(result['i'],fresh['i']) and np.array_equal(result['j'],fresh['j']),
        prescribed_axes=np.allclose(fresh['displacement'][~trial.state['free']],trial.prescribed(dt)[~trial.state['free']],atol=1e-16,rtol=0.),
        finite_fields=all(np.isfinite(fresh[key]).all() for key in ('force','stress','ep','temperature','density','gradient','position')))
    write(directory/'fresh_basic_gates.json',dict(gates=gates,accepted=False,restart_allowed=False))
    require(all(gates.values()),'fresh root force/energy/geometry gates failed')
    require(time.perf_counter()+1.5*map_wall+15<deadline,'insufficient budget for scalar full map')
    scalar=domain.make_trial(pa,trial.state)
    scalar.contact=ChainRuleDiscreteContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'],particle_radius=p['dx']/2,normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'],friction=p['friction'],base_z=-p['substrate_depth'])
    sr=scalar.evaluate(fresh['displacement'],dt)
    difference=batch.candidate_equivalence(dict(sr,energy_comparison=core.energy_diagnostics(scalar,sr)),
        dict(fresh,energy_comparison=energy))
    require(sr['graph_fingerprint']==fresh['graph_fingerprint'] and np.array_equal(sr['g'],fresh['g']) and
        np.array_equal(sr['thermal']['pair_conductance_W_K'],fresh['thermal']['pair_conductance_W_K']),
        'scalar map graph/heat differs')
    write(directory/'scalar_batch_equivalence.json',dict(status='PASS_SAME_FILTERED_FULL_MAP',differences=difference,
        scalar_full_map_evaluations=1,scalar_nnps_evaluations=1,wall_seconds=sr['wall_seconds']))
    geometry=post.geometry_audit(trial,fresh);write(directory/'geometry_audit.json',geometry)
    require(not geometry['transient_contacts_missed_by_endpoints'] and
        set(geometry['new_non_bore_features'])<=set(('SHOULDER','TOP')),'unqualified finite contact event')
    shell=audit_contact_shell_finite_enclosure(trial.contact,trial.x0,fresh['position'],core.pose(p,0.),core.pose(p,dt))
    write(directory/'continuous_shell_audit.json',shell)
    require(shell['status']==SHELL_CERTIFIED,'continuous shell unresolved')
    path=audit_visible_support_path(trial.x0,fresh['position'],fresh['i'],fresh['j'],core.pose(p,0.),core.pose(p,dt),p,pa.particle_id,
        maximum_depth=32,maximum_queries=256)
    write(directory/'visible_path_audit.json',path);require(path['status']==CERTIFIED,'whole retained graph path unresolved')
    observation=observe_square_bore_crossings(trial.state['bore_observation'],trial.x0,fresh['position'],
        pa.particle_id,pa.component_id,pa.m,0.,dt,lambda t:core.pose(p,t))
    write(directory/'UNACCEPTED_crossing_observation.json',dict(accepted=False,restart_allowed=False,observation=observation))
    return dict(status='PASS_FRESH_ROOT_RESIDUAL_AND_GEOMETRY_DIAGNOSTICS_NOT_PROCESS_ACCEPTANCE',gates=gates,
        initial_elastic_rank='NOT_RUN',temporal_and_domain_convergence='NOT_RUN',accepted_time_advanced_s=0.)


def worker(stage,dt_key,base):
    directory=base/stage;protocol=read(directory/'protocol.json');verify(protocol['source_and_input_sha256'])
    begin=time.perf_counter();trial=None;initial_fields=None;state=None;report=dict(status='RUNNING_UNACCEPTED_DIAGNOSTIC',
        stage=stage,dt_s=DT_CHOICES[dt_key],classification=domain.CLASSIFICATION,model_scope=domain.SCOPE,
        accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.,initial_elastic_rank='NOT_RUN',
        temporal_domain_and_paper_validation='NOT_RUN',source_and_input_sha256=protocol['source_and_input_sha256'],
        model_source_and_input_sha256=protocol['model_source_and_input_sha256'],
        thread_environment={key:os.environ.get(key) for key in THREAD_KEYS})
    def timeout():
        try:write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NO_ACCEPTANCE',elapsed_s=time.perf_counter()-begin,
            accepted=False,restart_allowed=False,calls=[] if trial is None else trial.calls))
        finally:os._exit(124)
    timer=threading.Timer(540.,timeout);timer.daemon=True;timer.start();code=2
    try:
        require(all(os.environ.get(key)=='1' for key in THREAD_KEYS),'single-thread environment required')
        pa,state,p,provenance,audit,snapshot=load_independent_initial();dt=DT_CHOICES[dt_key]
        initial_fields={key:prop.get_npy_array().copy() for key,prop in pa.properties.items()};initial_state=copy.deepcopy(state)
        trial=RecordedTrial(domain.make_trial(pa,state));guess=batch.guess_for(trial,dt)
        write(directory/'independent_initial_reference.json',dict(checkpoint=rel(snapshot),sha256=sha(snapshot),
            parameters=p,provenance=provenance,initial_audit=audit,parent_checkpoint=None,reused_342_state=False))
        np.savez_compressed(directory/'guess.npz',accepted=False,restart_allowed=False,particle_id=pa.particle_id,
            displacement_m=guess,dt_s=dt,free_axes=state['free'])
        report.update(particle_count=len(pa.x),free_displacement_dof=int(state['free'].sum()),initial_energy=audit['initial_energy'],
            theoretical_dense_allocations=audit['theoretical_float64_dense_allocations'],initialization_sha256=sha(snapshot))
        if stage=='cost':
            start=time.perf_counter();r=trial.evaluate(guess,dt);eval_wall=time.perf_counter()-start
            save_full_map(directory/'full_map',trial,r,'FULL_MAP_COST_NOT_ACCEPTED')
            report.update(full_map_wall_s=eval_wall,full_map_metrics=core.summary(r))
            write(directory/'completed_full_map_receipt.json',dict(full_map_wall_s=eval_wall,metrics=core.summary(r),
                peak_process_working_set_bytes=core.memory_bytes(),accepted=False,restart_allowed=False))
            require(time.perf_counter()<begin+525,'budget exhausted before tangent probe')
            start=time.perf_counter();direction,linear=core.newton_direction(trial,r);tangent_wall=time.perf_counter()-start
            np.savez_compressed(directory/'tangent_search_direction.npz',accepted=False,restart_allowed=False,
                displacement_correction_m=direction,particle_id=pa.particle_id)
            write(directory/'linear_probe.json',dict(linear=linear,tangent_plus_linear_wall_s=tangent_wall,
                initial_elastic_rank='NOT_RUN',search_tangent_is_not_full_jacobian=True))
            report.update(status='PASS_MEASURED_EXPANDED_FULL_MAP_COST_NOT_PHYSICAL_VALIDATION',
                tangent_linear_wall_s=tangent_wall,linear=linear,root_budget=root_budget(eval_wall,tangent_wall))
        else:
            cost=read(base/'cost/manifest.json')
            require(cost['status']=='PASS_MEASURED_EXPANDED_FULL_MAP_COST_NOT_PHYSICAL_VALIDATION' and
                cost['model_source_and_input_sha256']==protocol['model_source_and_input_sha256'] and cost['dt_s']==dt,
                'matching completed independent cost receipt required')
            budget=root_budget(cost['full_map_wall_s'],cost['tangent_linear_wall_s'])
            require(budget['useful_root_fits'] and budget['forecast_s']<=450,'measured root forecast exceeds budget')
            report['root_budget']=budget
            attempt=directory/'root_attempt';attempt.mkdir()
            result,metrics=core.solve(trial,dt,guess,attempt,begin+525,budget['maximum_iterations'])
            save_full_map(directory/'solved_map',trial,result,'SOLVED_ROOT_NOT_ACCEPTED')
            gates=check_fresh_root(trial,result,directory,begin+525,cost['full_map_wall_s'])
            report.update(status='PASS_EXPANDED_SINGLE_ROOT_DIAGNOSTIC_NOT_PROCESS_ACCEPTANCE',root_metrics=metrics,fresh_gates=gates)
        for key,values in initial_fields.items():require(np.array_equal(pa.properties[key].get_npy_array(),values),'initial particle fields were mutated: '+key)
        require(core.clean(state)==core.clean(initial_state),'independent initial metadata/history mutated')
        verify(protocol['source_and_input_sha256']);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_EXPANDED_DIAGNOSTIC_NO_ACCEPTANCE',error=repr(error))
        if trial is not None and trial.last_complete_result is not None:
            try:save_full_map(directory/'last_complete_map',trial,trial.last_complete_result,'LAST_COMPLETE_NOT_NECESSARILY_BEST_NOT_ACCEPTED')
            except Exception as evidence_error:report['failure_evidence_error']=repr(evidence_error)
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes(),exit_code=code,
        calls=[] if trial is None else trial.calls,
        completed_full_map_count=0 if trial is None else sum(c['kind']=='complete_map' and c['completed'] for c in trial.calls),
        nnps_count_note='Each completed full map has one real NNPS. Failed attempts may have reached NNPS; count is not inferred.',
        artifact_sha256={rel(f):sha(f) for f in sorted(directory.rglob('*')) if f.is_file() and f.name!='run.log'})
    write(directory/'manifest.json',report);print(json.dumps(core.clean(report),indent=2));return code


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('cost','root'),required=True)
    parser.add_argument('--dt-us',choices=tuple(DT_CHOICES),required=True)
    parser.add_argument('--tag',default='v3');parser.add_argument('--worker',action='store_true')
    args=parser.parse_args();require(args.tag.isalnum(),'alphanumeric new evidence tag required')
    base=BASE/args.tag/('dt_'+args.dt_us.replace('.','_')+'us');directory=base/args.stage
    if args.worker:return worker(args.stage,args.dt_us,base)
    require(not directory.exists(),'refuse diagnostic evidence overwrite')
    baseline=baseline_bindings();bindings=dict(baseline)
    if args.stage=='root':
        cost=read(base/'cost/manifest.json');receipt=read(base/'cost/supervision.json')
        require(cost['status']=='PASS_MEASURED_EXPANDED_FULL_MAP_COST_NOT_PHYSICAL_VALIDATION' and receipt['exit_code']==0 and
            cost['model_source_and_input_sha256']==baseline and cost['root_budget']['useful_root_fits'],
            'matching completed affordable cost required before scheduling root')
        for f in sorted((base/'cost').rglob('*')):
            if f.is_file():bindings[rel(f)]=sha(f)
    directory.mkdir(parents=True,exist_ok=False)
    write(directory/'protocol.json',dict(contract=rel(CONTRACT),stage=args.stage,dt_s=DT_CHOICES[args.dt_us],
        source_and_input_sha256=bindings,model_source_and_input_sha256=baseline,registered_before_worker=True,
        internal_limit_s=540.,external_limit_s=600.,root_forecast_limit_s=450.,accepted=False,restart_allowed=False,
        independent_initialization=rel(INITIAL/'summary.json'),no_old_accepted_copy=True,no_accepted_checkpoint_publication=True))
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--dt-us',args.dt_us,'--tag',args.tag,'--worker']
    begin=time.perf_counter()
    with (directory/'run.log').open('x',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,**{key:'1' for key in THREAD_KEYS})).returncode
        except subprocess.TimeoutExpired:code=124
    source_error=None
    try:verify(bindings)
    except Exception as error:source_error=repr(error);code=2
    write(directory/'supervision.json',dict(command=command,exit_code=code,external_wall_s=time.perf_counter()-begin,
        external_limit_s=600.,log_sha256=sha(directory/'run.log'),manifest_sha256=sha(directory/'manifest.json') if
        (directory/'manifest.json').exists() else None,before_after_hashes_match=source_error is None,source_error=source_error))
    print(json.dumps(dict(stage=args.stage,dt_us=args.dt_us,exit_code=code,directory=rel(directory))))
    return code


if __name__=='__main__':raise SystemExit(main())

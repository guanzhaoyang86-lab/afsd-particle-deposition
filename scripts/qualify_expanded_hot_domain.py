"""Supervised independent 8 mm initial rank/maps; no process acceptance."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import expanded_hot_domain as domain
import run_batched_process as batch
from check_visible_trial import same_map
from afsd_backend.pysph.mechanics import operators, positions
from afsd_backend.pysph.quasistatic_process import _geometry
from afsd_backend.pysph.quasistatic_tangent import (
    gradient_map, material_tangent, assemble_frozen_force_tangent)
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from batched_contact import BatchedChainRuleContact
from visible_support import filter_visible_support, _hash_array
from visible_support_numerics import require_visible_support_numerical_domain

ROOT = domain.ROOT
core = batch.core
BASE = ROOT/'outputs/expanded_hot_qualification_20261005'
CONTRACT = ROOT/'docs/EXPANDED_HOT_DOMAIN_QUALIFICATION_CONTRACT_20261005.md'
INITIAL = ROOT/'outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2'
DISCRETE_INITIAL = ROOT/'outputs/expanded_hot_domain_20261004/initial_8mm_discrete_v2'
COST = ROOT/'outputs/expanded_hot_domain_cost_20261004/v3/dt_15_625us/cost'
INITIAL_SHA = '31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3'
DISCRETE_SHA = 'bceace049bba901796576577f7937744331e03a04738bead7978dd6d5153524c'
DT = 15.625e-6
QUALIFICATION = 'EXPANDED_8MM_INITIAL_ELASTIC_RANK_AND_SAME_MAP_V1'
PASS = 'PASS_EXPANDED_INITIAL_ELASTIC_RANK_AND_SAME_MAP_NOT_PHYSICAL_VALIDATION'
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')
IMMUTABLE = ('particle_id', 'x0', 'y0', 'z0', 'm', 'h', 'material_id',
             'component_id', 'tracer_type', 'wholebar_tag', 'actuator')


def require(value, message):
    if not bool(value): raise ValueError(message)


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rel(path): return Path(path).resolve().relative_to(ROOT).as_posix()
def read(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def write(path, value): domain.write_new(Path(path), value)


def verify(mapping):
    for name, digest in mapping.items():
        path = (ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and sha(path) == digest, 'source/input drift: '+name)


def immutable_hashes(pa, state):
    return dict({name:_hash_array(getattr(pa,name)) for name in IMMUTABLE},
                free_axes=_hash_array(np.asarray(state['free'])))


def contact_factory(cls, p):
    return cls(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2,
        normal_stiffness=p['normal_stiffness'], tangential_stiffness=p['normal_stiffness'],
        friction=p['friction'], base_z=-p['substrate_depth'])


def inputs():
    mapping = {}
    def merge(values):
        for name, digest in values.items():
            require(name not in mapping or mapping[name] == digest, 'conflicting source binding: '+name)
            mapping[name] = digest
    merge(core.frozen_sources()); merge(domain.source_hashes())
    files = [Path(__file__).resolve(), CONTRACT, ROOT/'scripts/expanded_hot_acceptance.py',
        ROOT/'scripts/check_visible_trial.py', ROOT/'scripts/run_batched_process.py',
        ROOT/'scripts/run_discrete_post_start.py', ROOT/'scripts/secant_equilibrium.py']
    for directory, expected in ((INITIAL, INITIAL_SHA), (DISCRETE_INITIAL, DISCRETE_SHA)):
        receipt = read(directory/'summary.json')
        require(receipt['status'] == 'PASS_INITIAL_ARRAYS_AND_GEOMETRY_ONLY_NOT_PHYSICAL_QUALIFICATION',
                'independent initial array/geometry receipt missing')
        require(receipt['initialization_sha256'] == expected and
                sha(ROOT/receipt['initialization_checkpoint']) == expected, 'fixed independent snapshot differs')
        merge(receipt['source_sha256'])
        files.extend(directory/name for name in ('summary.json', 'initialization.json',
            'resolved_parameters.json', 'initial_particle_fields.npz', 'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'))
    cost = read(COST/'manifest.json'); supervision = read(COST/'supervision.json')
    require(cost['status'] == 'PASS_MEASURED_EXPANDED_FULL_MAP_COST_NOT_PHYSICAL_VALIDATION'
        and cost['dt_s'] == DT and cost['initialization_sha256'] == INITIAL_SHA
        and supervision['exit_code'] == 0 and supervision['before_after_hashes_match']
        and supervision['manifest_sha256'] == sha(COST/'manifest.json')
        and supervision['log_sha256'] == sha(COST/'run.log'), 'matching measured cost receipt required')
    merge(cost['source_and_input_sha256'])
    files.extend(COST/name for name in ('manifest.json', 'supervision.json', 'protocol.json', 'run.log'))
    merge({rel(path):sha(path) for path in files}); verify(mapping)
    forecast = 30+1.5*(3*cost['full_map_wall_s']+2*cost['tangent_linear_wall_s'])+120
    require(forecast <= 450., 'qualification forecast exceeds registered budget')
    return mapping, dict(forecast_s=forecast, forecast_limit_s=450.,
        measured_full_map_s=cost['full_map_wall_s'], measured_tangent_linear_s=cost['tangent_linear_wall_s'],
        elastic_eigen_and_io_allowance_s=120., eigen_time_is_not_previously_measured=True,
        formula='30 + 1.5*(3*measured_full_map_s + 2*measured_tangent_linear_s) + 120')


def load_initial(directory, digest):
    snapshot = directory/'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
    require(sha(snapshot) == digest, 'independent initial snapshot hash differs')
    saved = domain.load(str(snapshot)); pa, state = saved['arrays']['solid'], saved['solver_data']
    p, provenance = domain.parameters(.008)
    require(state['parameters'] == p and state['parameter_provenance'] == provenance, 'resolved inputs differ')
    audit = domain.initial_geometry_and_field_audit(pa,state)
    require(audit['particle_count'] == 918 and audit['free_displacement_dof'] == 1961, 'wrong independent population')
    require('spatial_model_transition' not in state and 'branch_origin_sha256' not in state, 'invented parent/transition forbidden')
    return pa, state, audit


def initial_rank(trial, directory):
    start = time.perf_counter(); pa, p = trial.pa, trial.p; x = positions(pa); zero = np.zeros_like(x)
    contact_rows = {}
    for name, cls in (('scalar', ChainRuleDiscreteContact), ('batched', BatchedChainRuleContact)):
        c = contact_factory(cls,p).evaluate(x,x,core.pose(p,0.),core.pose(p,0.),DT,trial.state['contact_history'])
        row = dict(maximum_force_N=float(np.max(np.abs(c['force_N']))),
                   normal_potential_J=float(np.sum(c['normal_potential1_J'])))
        contact_rows[name] = row
        require(row['maximum_force_N'] <= 1e-10 and row['normal_potential_J'] <= 1e-20,
                'initial unloaded contact equilibrium differs')
    i,j,g,_,_ = operators(pa)
    numerical = require_visible_support_numerical_domain(x,i,j,core.pose(p,0.),p,pa.particle_id)
    graph = filter_visible_support(x,i,j,g,core.pose(p,0.),p,pa.particle_id)
    require(graph['diagnostics']['blocked_edge_count'] == 0 and all(
        np.array_equal(graph[key], a) for key,a in (('i',i),('j',j),('g',g))), 'initial graphs are not identical')
    V = pa.m/pa.rho; H,B = _geometry(x,zero,pa.m,pa.rho,i,j,g)
    moment = np.zeros((len(x),3,3)); pair = (x[j]-x[i])[:,:,None]*g[:,None,:]
    np.add.at(moment,i,V[j,None,None]*pair); np.add.at(moment,j,V[i,None,None]*pair)
    mev = np.linalg.eigvalsh(moment); condition = mev[:,-1]/mev[:,0]
    require(np.all(mev[:,0]>0) and np.all(condition<=1e6), 'initial local correction rank/condition failed')
    require(np.allclose(moment@B,np.eye(3),rtol=1e-12,atol=1e-12), 'moment/correction identity failed')
    D = gradient_map(V,i,j,g,B)
    C = material_tangent(trial.s0,pa.ep,H,pa.temperature,DT,p['shear'],trial.law)
    require(not np.any(C['plastic_active']) and not np.any(C['rate_per_s']), 'initial branch is not zero-rate elastic')
    elastic = np.zeros((9,9))
    for a in range(3):
        for b in range(3):
            basis = np.zeros((3,3)); basis[a,b] = 1.
            dev = .5*(basis+basis.T)-np.trace(basis)*np.eye(3)/3
            elastic[:,3*a+b] = (2*p['shear']*DT*dev).ravel()
    require(np.allclose(C['tangent_Pa_s'],elastic,rtol=1e-12,atol=0.), 'initial material tangent differs from elastic derivative')
    J = assemble_frozen_force_tangent(D,V,C['tangent_Pa_s'],pa.rho,p['rho0'],trial.bulk,DT,zero)
    free = np.asarray(trial.state['free']).ravel(); K = -J[np.ix_(free,free)]; del J
    asym = float(np.max(np.abs(K-K.T))); scale = float(np.max(np.abs(K)))
    require(asym <= 1e-12*scale, 'elastic stiffness asymmetry exceeds roundoff gate')
    eigen_start = time.perf_counter(); eigen, vectors = np.linalg.eigh(K)
    eigen_wall = time.perf_counter()-eigen_start; least = vectors[:,0].copy(); del vectors
    residual = float(np.max(np.abs(K@least-eigen[0]*least)))
    relative_residual = residual/(np.linalg.norm(K,np.inf)*np.max(np.abs(least)))
    u = np.zeros(3*len(x)); u[free] = 1e-9*least
    strain = (D@u).reshape(-1,3,3); symmetric = .5*(strain+strain.swapaxes(1,2))
    trace = np.trace(strain,axis1=1,axis2=2); dev = symmetric-trace[:,None,None]*np.eye(3)/3
    direct_energy = float(np.sum(V*(p['shear']*np.sum(dev*dev,axis=(1,2))+.5*trial.bulk*trace**2)))
    matrix_energy = float(.5*u[free]@K@u[free])
    np.savez_compressed(directory/'initial_elastic_rank_arrays.npz',accepted=False,restart_allowed=False,
        particle_id=pa.particle_id,position_m=x,volume_m3=V,free_axes=trial.state['free'],i=i,j=j,g=g,
        moment=moment,moment_eigenvalues=mev,correction=B,free_stiffness_N_m=K,
        stiffness_eigenvalues_N_m=eigen,least_unit_eigenvector=least,least_mode_displacement_m=u.reshape(-1,3))
    row = dict(status='PASS_INITIAL_FREE_ELASTIC_RANK_NOT_FUTURE_FLOW_STABILITY',
        raw_visible_graph_identical=True,graph=graph['diagnostics'],numerical_domain=numerical,
        free_displacement_dof=int(free.sum()),fixed_displacement_dof=int((~free).sum()),
        dimensionless_moment_min_eigenvalue=float(mev.min()),maximum_moment_condition=float(condition.max()),
        stiffness_min_eigenvalue_N_m=float(eigen[0]),stiffness_max_eigenvalue_N_m=float(eigen[-1]),
        relative_min_eigenvalue=float(eigen[0]/eigen[-1]),rank_relative_threshold=1e-12,
        stiffness_max_asymmetry_N_m=asym,least_eigenpair_residual_N_m=residual,
        least_eigenpair_relative_residual=float(relative_residual),matrix_least_mode_energy_J=matrix_energy,
        direct_strain_least_mode_energy_J=direct_energy,least_mode_amplitude_m=1e-9,
        initial_contact_equilibrium=contact_rows,nnps_evaluations=1,full_map_evaluations=0,
        mass_shift=False,pseudoinverse=False,artificial_axis_constraints=False,matrix_symmetrized=False,
        eigen_wall_seconds=eigen_wall,wall_seconds=time.perf_counter()-start,
        D_shape=list(D.shape),D_bytes=D.nbytes,free_K_bytes=K.nbytes)
    valid = eigen[0] > 1e-12*eigen[-1] > 0 and relative_residual <= 1e-10 and np.isclose(
        matrix_energy,direct_energy,rtol=1e-9,atol=0.)
    if not valid: row['status'] = 'FAILED_INITIAL_ELASTIC_RANK_OR_NUMERICAL_CONSISTENCY'
    write(directory/'initial_elastic_rank.json',row)
    require(valid,'initial free elastic rank/eigen/energy consistency failed')
    return row


def save_map(directory, trial, result):
    directory.mkdir(exist_ok=False)
    core.journal(directory,trial,result,'FULL_MAP_NOT_ACCEPTED')
    np.savez_compressed(directory/'operators_and_fields.npz',accepted=False,restart_allowed=False,
        particle_id=trial.pa.particle_id,component_id=trial.pa.component_id,mass_kg=trial.pa.m,
        free_axes=trial.state['free'],i=result['i'],j=result['j'],g=result['g'],
        raw_i=result.get('raw_i',result['i']),raw_j=result.get('raw_j',result['j']),raw_g=result.get('raw_g',result['g']),
        correction=result['correction'],volume_m3=result['volume'],gradient_per_s=result['gradient'],
        material_temperature_K=result['material_temperature'],plastic_work_density_J_m3=result['material'][2],
        return_dissipation_density_J_m3=result['material'][3],effective_rate_per_s=result['rate'],
        pair_conductance_W_K=result['thermal']['pair_conductance_W_K'],pair_heat_into_i_J=result['thermal']['pair_heat_into_i_J'])
    write(directory/'thermal_and_graph.json',dict(thermal=result['thermal'],
        visibility=result.get('visibility_graph',{}),numerical_domain=result.get('visibility_numerical_domain',{})))


def initial_same_map(visible, discrete, directory):
    start=time.perf_counter(); guess=batch.guess_for(visible,DT)
    require(np.array_equal(guess,batch.guess_for(discrete,DT)), 'independent prescribed guesses differ')
    discrete.contact=contact_factory(BatchedChainRuleContact,discrete.p)
    a=discrete.evaluate(guess,DT); save_map(directory/'discrete_batched',discrete,a)
    b=visible.evaluate(guess,DT); save_map(directory/'visible_batched',visible,b)
    require(b['visibility_graph']['diagnostics']['blocked_edge_count']==0,'prescribed forecast graph is occluded')
    first=same_map(a,b,discrete,visible)
    started=time.perf_counter(); Ja=discrete.tangent(a); Jb=visible.tangent(b)
    tangent_wall=time.perf_counter()-started
    require(np.array_equal(Ja,Jb),'identical initial maps have different actual search tangents')
    tangent_difference=float(np.max(np.abs(Ja-Jb))); del Ja,Jb
    scalar=domain.make_trial(visible.pa,visible.state);scalar.contact=contact_factory(ChainRuleDiscreteContact,scalar.p)
    c=scalar.evaluate(guess,DT); save_map(directory/'visible_scalar',scalar,c)
    second=same_map(c,b,scalar,visible)
    row=dict(status='PASS_INDEPENDENT_INITIAL_SAME_GRAPH_COMPLETE_MAP_AND_TANGENT',
        reference_dt_s=DT,prescribed_guess_only_not_root=True,discrete_visible=first,scalar_batched_visible=second,
        blocked_edge_count=0,edge_count=len(b['i']),full_map_evaluations=3,nnps_evaluations=3,
        tangent_evaluations=2,tangent_maximum_difference_N_m=tangent_difference,
        tangent_wall_seconds=tangent_wall,wall_seconds=time.perf_counter()-start,
        raw_energy_is_compared_not_accepted=True,graph_fingerprint=b['graph_fingerprint'])
    write(directory/'initial_same_map.json',row);return row


def worker(directory):
    begin=time.perf_counter(); protocol=read(directory/'protocol.json'); code=2
    report=dict(status='RUNNING_EXPANDED_INITIAL_DIAGNOSTIC',qualification=QUALIFICATION,
        qualification_passed=False,classification=domain.CLASSIFICATION,model_scope=domain.SCOPE,
        reference_dt_s=DT,accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.,
        source_and_input_sha256=protocol['source_and_input_sha256'],budget=protocol['budget'],
        initial_elastic_rank='NOT_COMPLETED',initial_same_map='NOT_COMPLETED')
    def timeout():
        try:write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NO_ACCEPTANCE',wall_seconds=time.perf_counter()-begin,
            accepted=False,restart_allowed=False))
        finally:os._exit(124)
    timer=threading.Timer(540.,timeout);timer.daemon=True;timer.start()
    try:
        require(all(os.environ.get(key)=='1' for key in THREADS),'single-thread environment required')
        verify(protocol['source_and_input_sha256'])
        pa,state,audit=load_initial(INITIAL,INITIAL_SHA);other,otherstate,_=load_initial(DISCRETE_INITIAL,DISCRETE_SHA)
        require(state['operator_variant']=='visible' and otherstate['operator_variant']=='discrete','initial variants differ')
        require(pa.properties.keys()==other.properties.keys() and all(np.array_equal(
            pa.properties[key].get_npy_array(),other.properties[key].get_npy_array()) for key in pa.properties),
            'independent initial physical arrays differ')
        states=[copy.deepcopy(state),copy.deepcopy(otherstate)]
        arrays=[{key:prop.get_npy_array().copy() for key,prop in obj.properties.items()} for obj in (pa,other)]
        normalized=copy.deepcopy(otherstate);normalized['classification']=state['classification']
        normalized['operator_variant']=state['operator_variant']
        normalized['independent_domain_initialization']['spatial_operator_classification']=state['independent_domain_initialization']['spatial_operator_classification']
        # The factory records the variant inside its pure initial audit as well.
        if 'initial_array_geometry_audit' in normalized:
            normalized['initial_array_geometry_audit']['classification']=state['initial_array_geometry_audit']['classification']
        require(core.clean(normalized)==core.clean(state),'independent initial histories differ beyond declared operator labels')
        report.update(particle_count=len(pa.x),free_displacement_dof=int(state['free'].sum()),parameters=state['parameters'],
            parameter_provenance=state['parameter_provenance'],initialization_checkpoint=rel(INITIAL/'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'),
            initialization_sha256=INITIAL_SHA,discrete_initialization_checkpoint=rel(DISCRETE_INITIAL/'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'),
            discrete_initialization_sha256=DISCRETE_SHA,independent_domain_initialization=state['independent_domain_initialization'],
            immutable_array_sha256=immutable_hashes(pa,state),initial_energy=audit['initial_energy'],
            theoretical_float64_dense_allocations=audit['theoretical_float64_dense_allocations'])
        visible=domain.make_trial(pa,state);discrete=domain.make_trial(other,otherstate)
        report['initial_elastic_rank']=initial_rank(visible,directory)
        report['initial_same_map']=initial_same_map(visible,discrete,directory)
        for obj,current,before,fields in zip((pa,other),(state,otherstate),states,arrays):
            require(core.clean(current)==core.clean(before),'initial state mutated')
            for key,value in fields.items():require(np.array_equal(obj.properties[key].get_npy_array(),value),'initial field mutated: '+key)
        verify(protocol['source_and_input_sha256']);report.update(status=PASS,qualification_passed=True,
            full_map_evaluations=3,nnps_evaluations=4,initial_states_unchanged=True,
            root_solve='NOT_RUN_BY_THIS_DIAGNOSTIC',temporal_domain_inlet_and_paper_validation='NOT_RUN')
        code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_EXPANDED_INITIAL_QUALIFICATION_NO_ACCEPTANCE',error=repr(error))
    finally:timer.cancel()
    report.update(exit_code=code,wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes(),
        output_sha256={rel(path):sha(path) for path in sorted(directory.rglob('*')) if path.is_file() and path.name!='run.log'})
    write(directory/'summary.json',report);print(json.dumps(core.clean(report),indent=2));return code


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--tag',default='v1');parser.add_argument('--worker',action='store_true')
    args=parser.parse_args();require(args.tag.isalnum(),'new alphanumeric evidence tag required');directory=BASE/args.tag
    if args.worker:return worker(directory)
    require(not directory.exists(),'refuse qualification evidence overwrite');mapping,budget=inputs()
    directory.mkdir(parents=True,exist_ok=False)
    write(directory/'protocol.json',dict(contract=rel(CONTRACT),contract_sha256=sha(CONTRACT),qualification=QUALIFICATION,
        registered_before_worker=True,reference_dt_s=DT,source_and_input_sha256=mapping,budget=budget,
        internal_limit_s=540.,external_limit_s=600.,accepted=False,restart_allowed=False))
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--tag',args.tag,'--worker'];begin=time.perf_counter()
    with (directory/'run.log').open('x',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=600,
            env=dict(os.environ,**{key:'1' for key in THREADS})).returncode
        except subprocess.TimeoutExpired:code=124
    error=None
    try:verify(mapping)
    except Exception as problem:error=repr(problem);code=2
    write(directory/'supervision.json',dict(command=command,exit_code=code,external_wall_s=time.perf_counter()-begin,
        external_limit_s=600.,log_sha256=sha(directory/'run.log'),summary_sha256=sha(directory/'summary.json') if
        (directory/'summary.json').exists() else None,before_after_hashes_match=error is None,source_error=error))
    print(json.dumps(dict(status='QUALIFICATION_WORKER_COMPLETED',exit_code=code,directory=rel(directory))));return code


if __name__=='__main__':raise SystemExit(main())

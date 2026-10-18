"""Bounded, independently initialized inlet diagnostic helpers; never process acceptance."""
from pathlib import Path
import copy
import json
import time
import numpy as np
import expanded_hot_domain as domain
import inlet_acceptance_v2 as inlet
import inlet_source_work as source
import run_batched_process as batch
import run_discrete_contact_process as core
import run_discrete_post_start as post
import run_visible_support_process as visible_run
from visible_quasistatic import VisibleSupportTrial
from batched_contact import BatchedChainRuleContact
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from afsd_backend.pysph.discrete_process import accepted_copy as fixed_inventory_copy
from afsd_backend.pysph.mechanics import positions, bulk_modulus
from afsd_backend.pysph.material import FKS
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.quasistatic_tangent import gradient_map, assemble_frozen_force_tangent
from afsd_backend.pysph.interval_acceptance import validate_interval_energy_state
from contact_shell_finite_enclosure import audit_contact_shell_finite_enclosure
from contact_shell_audit import CERTIFIED as SHELL_CERTIFIED
from visible_support_path_tight import audit_visible_support_path, CERTIFIED as PATH_CERTIFIED

ROOT=domain.ROOT
CLASSIFICATION='HOT_QS_TWO_INVENTORY_BIRTH_DIAGNOSTIC_V1'
SCOPE=inlet.DIAGNOSTIC_SCOPE
DT=3.90625e-6
CONTRACT=ROOT/'docs/TWO_INSERTION_DIAGNOSTIC_CONTRACT_20261005.md'
INITIAL=ROOT/'outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
INITIAL_SHA='31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3'
require=domain.require
sha=domain.sha
write=domain.write_new
clean=domain.clean


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()
def verify(mapping):
    for name,digest in mapping.items():
        path=(ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and sha(path)==digest,'source/input changed: '+name)


def contact_factory(cls,p):
    return cls(tool_radius=p['tool_radius'],bore_width=p['bore_width'],tool_height=p['tool_top']-p['layer'],
        particle_radius=p['dx']/2,normal_stiffness=p['normal_stiffness'],tangential_stiffness=p['normal_stiffness'],
        friction=p['friction'],base_z=-p['substrate_depth'])


class InletTrial(VisibleSupportTrial):
    """Reuse the complete visible map; no old process classification adapter."""
    def __init__(self,pa,state,p):
        require(state.get('classification')==CLASSIFICATION and state.get('model_scope')==SCOPE and
            state.get('accepted') is False and state.get('restart_allowed') is False,'independent diagnostic state required')
        require('spatial_model_transition' not in state,'no fabricated visible process transition')
        validate_interval_energy_state(state['interval_energy'])
        require(state['interval_energy']['steps']==state['step'],'full diagnostic energy prefix missing')
        pending=state['inlet_acceptance'];aligned=inlet.process_state_for(pa,pending)
        require(pending['phase']=='PENDING' and pending['mechanical_step_count']==state['step'] and
            state['time_s']==pending['time_s'],'pending inlet/time/count differs')
        require(np.array_equal(state['free'],aligned['free']) and np.array_equal(state['contact_history'],aligned['contact_history']),
            'physical history not aligned by persistent ID')
        require(state['bore_observation']==inlet._legacy_observer(pending['bore_observation']),
            'only fixed-inventory internal legacy observer allowed')
        self.pa,self.state,self.p=pa,state,dict(p)
        self.x0=positions(pa);self.s0=pa.S.reshape(-1,3,3).copy()
        self.bulk=bulk_modulus(p['shear'],p['poisson']);self.law=FKS(**{key:p[key] for key in FKS.__dataclass_fields__})
        self.contact=contact_factory(BatchedChainRuleContact,p)
        self.calls=[];self.last_complete_result=None

    def evaluate(self,displacement,dt):
        started=time.perf_counter();row=dict(kind='full_map',completed=False);self.calls.append(row)
        try:
            result=super().evaluate(displacement,dt);self.last_complete_result=result
            row.update(completed=True,scaled_residual=result['scaled_residual'],graph_fingerprint=result['graph_fingerprint'])
            return result
        except Exception as error:row['error']=repr(error);raise
        finally:row['wall_seconds']=time.perf_counter()-started

    def tangent(self,result):
        started=time.perf_counter();row=dict(kind='actual_search_tangent',completed=False);self.calls.append(row)
        try:
            value=super().tangent(result);row['completed']=True;return value
        except Exception as error:row['error']=repr(error);raise
        finally:row['wall_seconds']=time.perf_counter()-started


def initial_state(frozen):
    require(sha(INITIAL)==INITIAL_SHA,'registered independent origin changed')
    saved=domain.load(str(INITIAL));pa=saved['arrays']['solid'];old=saved['solver_data']
    p,provenance=domain.parameters(.008)
    require(old['parameters']==p and old['parameter_provenance']==provenance,'independent physical parameters changed')
    require(old['classification']==domain.CLASSIFICATION and old['model_scope']==domain.SCOPE and
        old['time_s']==0. and old['step']==0 and old['accepted'] is False and old['restart_allowed'] is False,
        'independent time-zero origin required')
    require('spatial_model_transition' not in old and 'branch_origin_sha256' not in old,'old process transition forbidden')
    audit=domain.initial_geometry_and_field_audit(pa,old)
    registry=inlet.inventory.initialize_registry(pa,rel(CONTRACT),acknowledgement=inlet.inventory.ACKNOWLEDGEMENT)
    plane=dict(frame='WORLD_Z',plane_z_m=.0025,unit='m',source=rel(CONTRACT))
    v2=inlet.initialize_acceptance_state(pa,registry,old,plane,rel(CONTRACT))
    state=dict(classification=CLASSIFICATION,model_scope=SCOPE,accepted=False,restart_allowed=False,
        process_acceptance=False,independent_initialization=True,time_s=0.,step=0,parameters=p,parameter_provenance=provenance,
        inlet_acceptance=v2,bore_observation=v2['bore_observation'],interval_energy=copy.deepcopy(old['interval_energy']),
        crossing_ids=[],crossing_mass_kg=0.,source_sha256=frozen,
        diagnostic_origin=dict(checkpoint=rel(INITIAL),sha256=INITIAL_SHA,classification=old['classification'],
            model_scope=old['model_scope'],time_s=old['time_s'],step=old['step'],accepted=False,restart_allowed=False),
        parent_checkpoint_sha256=INITIAL_SHA,parent_checkpoint_classification=old['classification'],initial_array_audit=audit)
    return pa,state


def prepare_birth(parent,state,index):
    require(index in (1,2) and state['step']==index-1 and state['time_s']==(index-1)*DT,'registered birth order/time required')
    p=state['parameters'];xy=np.arange(-2,3)*p['dx'];z=(.00325,.00375)[index-1]
    x=np.array([[a,b,z] for a in xy for b in xy])
    incoming=make_array(x,np.ones(25,dtype=int),p);incoming.particle_id[:]=0
    incoming.u[:]=0.;incoming.v[:]=0.;incoming.w[:]=-p['feed'];incoming.temperature[:]=723.15
    fields={key:dict(unit=inlet.inventory.REQUIRED.get(key,'1'),source=rel(CONTRACT),
        status='INITIAL_CONDITION' if key in inlet.inventory.REQUIRED else 'UNRESOLVED') for key in incoming.properties}
    spec=dict(event_id='diagnostic-inventory-birth-'+str(index),time_s=state['time_s'],source=rel(CONTRACT),fields=fields,
        heat_capacity=dict(value=p['cp'],unit='J/(kg K)',source=rel(CONTRACT),status='DECLARED_ASSUMPTION'),
        reference_temperature=dict(value=p['room'],unit='K',source=rel(CONTRACT),status='DECLARED_ASSUMPTION'))
    pa,pending,event=inlet.propose_insertion(parent,state['inlet_acceptance'],incoming,spec,core.pose(p,state['time_s']))
    aligned=inlet.process_state_for(pa,pending)
    inner=copy.deepcopy(state);inner.update(aligned,inlet_acceptance=pending,
        bore_observation=inlet._legacy_observer(pending['bore_observation']),
        observer_adapter='LEGACY_VIEW_ONLY_INSIDE_ONE_FIXED_INVENTORY_STEP_NOT_PUBLISHED')
    require(len(pa.x)==918+25*index and int(aligned['free'].sum())==1961+50*index,'birth inventory/free DOF differs')
    require(len(pending['control_active_ids'])==25*(index+1),'unexpected inlet control selection')
    old=inlet.process_state_for(parent,state['inlet_acceptance'])
    rows={int(pid):k for k,pid in enumerate(pa.particle_id)}
    require(np.array_equal(aligned['contact_history'][[rows[int(pid)] for pid in parent.particle_id]],old['contact_history']),
        'birth changed retained physical history')
    born=np.isin(pa.particle_id,event['allocated_particle_ids']);pos=positions(pa)
    squared=np.sum((pos[born,None]-pos[None,:])**2,axis=2)
    for k,j in enumerate(np.flatnonzero(born)):squared[k,j]=np.inf
    separation=float(np.sqrt(squared.min()));margin=128*np.finfo(float).eps*max(p['substrate_width'],z)
    require(separation>=p['dx']-margin,'newborn nominal cells overlap existing inventory')
    contact=contact_factory(BatchedChainRuleContact,p);geometry=contact.geometry(pos,core.pose(p,state['time_s']))
    require(not np.any(geometry['potential_J'][born]) and not np.any(aligned['contact_history'][born]),
        'newborn initial contact storage must be exactly zero')
    parent_geometry=contact.geometry(positions(parent),core.pose(p,state['time_s']))
    parent_storage=np.column_stack((np.sum(parent_geometry['potential_J'],axis=1),
        .5*contact.kt*np.sum(old['contact_history']**2,axis=(1,2))))
    geometry_evidence=dict(status='PASS',minimum_newborn_pair_distance_m=separation,roundoff_margin_m=margin,
        born_particle_ids=event['allocated_particle_ids'],newborn_contact_storage_J=0.,
        initial_stress_zero=bool(not np.any(pa.S.reshape(-1,3,3)[born])),
        initial_density_reference=bool(np.all(pa.rho[born]==p['rho0'])),old_contact_history_unchanged=True,
        initial_kinematics_m_s=[0.,0.,-p['feed']],control_event=pending['control_events'][-1])
    return pa,inner,event,parent_storage,geometry_evidence


def reference_rank(trial,r,directory):
    """Current actual graph, unloaded elastic reference only; no nonlinear stability claim."""
    begin=time.perf_counter();p=trial.p;n=len(trial.pa.x);V=r['volume'];B=r['correction']
    moment=np.zeros((n,3,3));pair=(r['position'][r['j']]-r['position'][r['i']])[:,:,None]*r['g'][:,None,:]
    np.add.at(moment,r['i'],V[r['j'],None,None]*pair);np.add.at(moment,r['j'],V[r['i'],None,None]*pair)
    ev=np.linalg.eigvalsh(moment)
    require(np.all(ev[:,0]>0) and np.all(ev[:,-1]<=1e6*ev[:,0]),'actual correction rank/condition failed')
    require(np.allclose(moment@B,np.eye(3),rtol=1e-12,atol=1e-12),'actual correction does not invert actual moment')
    D=gradient_map(V,r['i'],r['j'],r['g'],B);elastic=np.zeros((9,9))
    for a in range(3):
        for b in range(3):
            basis=np.zeros((3,3));basis[a,b]=1.
            dev=.5*(basis+basis.T)-np.trace(basis)*np.eye(3)/3
            elastic[:,3*a+b]=(2*p['shear']*DT*dev).ravel()
    C=np.broadcast_to(elastic,(n,9,9));zero=np.zeros((n,3))
    J=assemble_frozen_force_tangent(D,V,C,np.full(n,p['rho0']),p['rho0'],trial.bulk,DT,zero)
    free=trial.state['free'].ravel();K=-J[np.ix_(free,free)];del J
    asym=float(np.max(np.abs(K-K.T)));require(asym<=1e-12*np.max(np.abs(K)),'elastic reference asymmetric')
    eig,vec=np.linalg.eigh(K);least=vec[:,0].copy();del vec
    residual=float(np.max(np.abs(K@least-eig[0]*least))/(np.linalg.norm(K,np.inf)*np.max(np.abs(least))))
    u=np.zeros(3*n);u[free]=1e-9*least;strain=(D@u).reshape(n,3,3)
    sym=.5*(strain+strain.swapaxes(1,2));tr=np.trace(strain,axis1=1,axis2=2);dev=sym-tr[:,None,None]*np.eye(3)/3
    direct=float(np.sum(V*(p['shear']*np.sum(dev*dev,axis=(1,2))+.5*trial.bulk*tr*tr)))
    matrix=float(.5*u[free]@K@u[free])
    valid=bool(eig[0]>1e-12*eig[-1]>0 and residual<=1e-10 and np.isclose(matrix,direct,rtol=1e-9,atol=0.))
    np.savez_compressed(directory/'reference_rank_arrays.npz',particle_id=trial.pa.particle_id,
        i=r['i'],j=r['j'],g=r['g'],volume=V,correction=B,free=trial.state['free'],
        stiffness_eigenvalues_N_m=eig,least_unit_eigenvector=least)
    row=dict(status='PASS' if valid else 'FAIL',scope='ACTUAL_GEOMETRY_UNLOADED_ELASTIC_REFERENCE_NOT_NONLINEAR_STABILITY',
        graph_fingerprint=r['graph_fingerprint'],particle_count=n,free_dof=int(free.sum()),
        moment_min_eigenvalue=float(ev.min()),moment_max_condition=float((ev[:,-1]/ev[:,0]).max()),
        stiffness_min_eigenvalue_N_m=float(eig[0]),stiffness_max_eigenvalue_N_m=float(eig[-1]),
        relative_min_eigenvalue=float(eig[0]/eig[-1]),least_eigenpair_relative_residual=residual,
        direct_strain_energy_J=direct,matrix_energy_J=matrix,maximum_asymmetry_N_m=asym,
        matrix_symmetrized=False,pseudoinverse=False,mass_or_modulus_shift=False,wall_seconds=time.perf_counter()-begin)
    write(directory/'reference_rank.json',row);require(valid,'actual graph elastic-reference rank failed');return row


def save_map(directory,trial,r):
    core.journal(directory,trial,r,'FULL_MAP_DIAGNOSTIC_NOT_ACCEPTED')
    visible_run.save_graph(directory,trial,r)
    np.savez_compressed(directory/'actual_fields.npz',particle_id=trial.pa.particle_id,position=r['position'],
        displacement=r['displacement'],density=r['density'],stress=r['stress'],ep=r['ep'],temperature=r['temperature'],
        force=r['force'],free=trial.state['free'],gradient=r['gradient'],volume=r['volume'],
        plastic_work_density=r['material'][2],return_dissipation_density=r['material'][3],rate=r['material'][4],
        contact_history=r['contact']['history'],pair_conductance=r['thermal']['pair_conductance_W_K'])
    write(directory/'actual_energy.json',dict(mechanical=core.energy_diagnostics(trial,r),thermal=r['thermal'],
        pose0=core.pose(trial.p,trial.state['time_s']),pose1=core.pose(trial.p,trial.state['time_s']+r['dt_s'])))


def qualify_candidate(parent,trial,r,event,parent_storage,birth_geometry,directory,frozen):
    """All existing gates first; new source accounting and V2 receipt second."""
    pa,state,fresh=fixed_inventory_copy(trial,r)
    require(state['interval_energy']['records'][:-1]==trial.state['interval_energy']['records'],'energy prefix reset')
    require(np.array_equal(r['i'],fresh['i']) and np.array_equal(r['j'],fresh['j']) and
        r['graph_fingerprint']==fresh['graph_fingerprint'],'fresh map graph changed')
    save_map(directory,trial,fresh);energy=core.energy_diagnostics(trial,fresh)
    scalar=InletTrial(trial.pa,trial.state,trial.p);scalar.contact=contact_factory(ChainRuleDiscreteContact,trial.p)
    sr=scalar.evaluate(fresh['displacement'],DT)
    diff=batch.candidate_equivalence(dict(sr,energy_comparison=core.energy_diagnostics(scalar,sr)),
        dict(fresh,energy_comparison=energy))
    require(sr['graph_fingerprint']==fresh['graph_fingerprint'] and np.array_equal(sr['g'],fresh['g']) and
        np.array_equal(sr['thermal']['pair_conductance_W_K'],fresh['thermal']['pair_conductance_W_K']),
        'scalar complete-map actual graph/heat changed')
    write(directory/'scalar_equivalence.json',dict(status='PASS',differences=diff,calls=scalar.calls))
    geometry=post.geometry_audit(trial,fresh);write(directory/'sampled_geometry.json',geometry)
    require(not geometry['transient_contacts_missed_by_endpoints'] and
        set(geometry['new_non_bore_features'])<=set(('SHOULDER','TOP')),'unqualified finite contact event')
    p=trial.p;t=trial.state['time_s'];pose0,pose1=core.pose(p,t),core.pose(p,t+DT)
    shell=audit_contact_shell_finite_enclosure(trial.contact,trial.x0,fresh['position'],pose0,pose1)
    write(directory/'continuous_shell.json',shell);require(shell['status']==SHELL_CERTIFIED,'continuous shell not certified')
    path=audit_visible_support_path(trial.x0,fresh['position'],fresh['i'],fresh['j'],pose0,pose1,p,trial.pa.particle_id,
        maximum_depth=32,maximum_queries=256)
    write(directory/'visible_path.json',path);require(path['status']==PATH_CERTIFIED,'kept support path not certified')
    rank=reference_rank(trial,fresh,directory)
    port=dict(contract=source.PORT_CONTRACT,surface='EXTERIOR_BIRTH_BOUNDARY_NOT_INTERNAL_CONTROL_PLANE',
        internal_interface_work='NOT_ADDED_INTERNAL_PAIR_WORK',traction_Pa=[0.,0.,0.],source=rel(CONTRACT))
    calculated=source.evaluate_source_work(parent,trial.pa,pa,event,parameters=p,free=trial.state['free'],result=fresh,
        mechanical_energy=energy,parent_contact_storage_J=parent_storage,external_port=port)
    write(directory/'source_work.json',calculated);work=source.source_work_record(calculated,directory/'source_work.json')
    write(directory/'source_work_v2.json',work)
    pending=trial.state['inlet_acceptance'];history=fresh['contact']['history'];free=trial.state['free']
    binding=inlet.mechanical_binding(pa,pending,history,free,state['time_s'])
    actual_force=float(np.max(np.abs(fresh['force'][free])))
    require(actual_force<=fresh['force_limit_N'] and fresh['scaled_residual']<=1.,'actual fresh force failed')
    details=dict(force=dict(actual_free_force_N=actual_force,force_limit_N=fresh['force_limit_N']),
        heat=dict(thermal=fresh['thermal']['ledger'],material_feedback_residual_K=fresh['thermal_residual_K']),
        energy=dict(energy=energy,interval_energy=state['interval_energy']),geometry_overlap=birth_geometry,
        NNPS=dict(graph_fingerprint=fresh['graph_fingerprint'],edge_count=len(fresh['i']),
            complete_map_rebuilt_actual_3d_nnps=True,visibility=fresh['visibility_graph']['diagnostics']),
        visible_support=dict(path_status=path['status'],path_sha256=sha(directory/'visible_path.json')),
        gradient_rank=rank,contact=dict(shell_status=shell['status'],shell_sha256=sha(directory/'continuous_shell.json'),
            scalar_sha256=sha(directory/'scalar_equivalence.json'),sweep=fresh['contact']['sweep']))
    gates={}
    for key,value in details.items():
        path=directory/('gate_'+key+'.json');write(path,dict(status='PASS',binding=binding,evidence=value,
            scope=SCOPE,source_sha256=frozen,source_work_report_sha256=sha(directory/'source_work.json')))
        gates[key]=dict(status='PASS',evidence_sha256=sha(path),artifact=rel(path))
    audit_artifacts={rel(path):sha(path) for path in directory.rglob('*') if path.is_file()}
    receipt=dict(scope=SCOPE,step_id='real-coupled-inlet-diagnostic-'+str(state['step']),binding=binding,
        gates=gates,source_work_sha256=inlet._sha(work),joint_inlet_benchmark='NOT_RUN',accepted=False,restart_allowed=False)
    receipt['audit_artifacts_sha256']=audit_artifacts
    write(directory/'mechanical_receipt.json',receipt)

    def external_verifier(given):
        require(given==read(directory/'mechanical_receipt.json'),'mechanical receipt file changed')
        verify(frozen);verify(given['audit_artifacts_sha256'])
        for key in inlet.GATES:
            gate=given['gates'][key];path=ROOT/gate['artifact'];obj=read(path)
            require(sha(path)==gate['evidence_sha256'] and obj['status']=='PASS' and obj['binding']==binding and
                obj['scope']==SCOPE and obj['source_sha256']==frozen,'qualified actual gate artifact mismatch: '+key)
            require(obj['source_work_report_sha256']==sha(directory/'source_work.json'),'source calculation drift')
        require(source.source_work_record(calculated,directory/'source_work.json')==work and
            read(directory/'source_work_v2.json')==work,'actual source-work artifact mismatch')
        require(binding==inlet.mechanical_binding(pa,pending,history,free,state['time_s']),'actual PA/history binding drift')
        return True

    sealed=inlet.seal_mechanical_acceptance(pa,pending,history,free,state['time_s'],receipt,work,
        external_verifier,lambda value:core.pose(p,value))
    state.update(inlet_acceptance=sealed,bore_observation=sealed['bore_observation'],accepted=False,restart_allowed=False,
        process_acceptance=False,observer_adapter='BIRTH_AWARE_V2_PUBLISHED; LEGACY_VIEW_WAS_INTERNAL_FIXED_N_ONLY',
        mechanical_receipt_sha256=sha(directory/'mechanical_receipt.json'),source_work_report_sha256=sha(directory/'source_work.json'))
    require(inlet._legacy_observer(sealed['bore_observation'])==core.clean(fresh_observer(trial,pa)),
        'sealed V2 observer differs from actual fixed-inventory observation')
    return pa,state,fresh,calculated


def fresh_observer(trial,pa):
    from afsd_backend.pysph.interval_acceptance import observe_square_bore_crossings
    return observe_square_bore_crossings(trial.state['bore_observation'],trial.x0,positions(pa),pa.particle_id,
        pa.component_id,pa.m,trial.state['time_s'],trial.state['time_s']+DT,lambda t:core.pose(trial.p,t))['next_state']


def write_and_reload(pa,state,path):
    require(not path.exists(),'diagnostic checkpoint already exists')
    core.checkpoint(path,pa,state)
    saved=domain.load(str(path));restored=saved['arrays']['solid'];actual=saved['solver_data']
    require(clean(actual)==clean(state),'diagnostic checkpoint metadata changed')
    require(inlet.inventory._array_digest(restored)==inlet.inventory._array_digest(pa),'checkpoint actual PA changed')
    v2=inlet.restore_acceptance_state(json.dumps(actual['inlet_acceptance']),restored)
    require(v2==state['inlet_acceptance'] and actual['bore_observation']==v2['bore_observation'],'V2 checkpoint observer changed')
    permutation=np.arange(len(restored.x)-1,-1,-1,dtype=np.uint32);sorted_pa=restored.extract_particles(permutation)
    v2=inlet.restore_acceptance_state(json.dumps(v2),sorted_pa);aligned=inlet.process_state_for(sorted_pa,v2)
    require(inlet.inventory._array_digest(sorted_pa)==inlet.inventory._array_digest(pa),'sorting changed physical fields')
    actual.update(contact_history=aligned['contact_history'],free=aligned['free'],inlet_acceptance=v2,
        parent_checkpoint_sha256=sha(path),parent_checkpoint_classification=actual['classification'])
    write(path.parent/'reload_and_sort.json',dict(status='PASS_DIAGNOSTIC_RELOAD_NOT_PROCESS_RESTART',
        checkpoint=rel(path),checkpoint_sha256=sha(path),actual_pa_sha256=inlet.inventory._array_digest(sorted_pa),
        source_events=len(v2['source_events']),mechanical_steps=v2['mechanical_step_count'],
        particle_ids_before=restored.particle_id.tolist(),particle_ids_after=sorted_pa.particle_id.tolist(),
        history_restored_by_persistent_id=True,observer_sha256=inlet._sha(v2['bore_observation']),
        process_restore_called=False,accepted=False,restart_allowed=False))
    return sorted_pa,actual

"""Computed source work for a declared stress-free inventory birth port.

Reads saved fields and an already evaluated coupled map; performs no physical
solve, NNPS or acceptance. Keeps the original QS raw energy defect unchanged.
"""
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.inlet_protocol import _array_digest
from afsd_backend.pysph.mechanics import bulk_modulus,volumetric_energy
from inlet_acceptance_v2 import WORK_CONVENTION,WORK_TERMS

CLASSIFICATION='STRESS_FREE_BIRTH_SOURCE_WORK_COMPONENT_ONLY_V1'
PORT_CONTRACT='STRESS_FREE_EXTERIOR_INVENTORY_BIRTH_PORT_V1'


def _json(value):
    def convert(item):
        if isinstance(item,(np.ndarray,np.generic)):return item.tolist()
        raise TypeError('unsupported source-work evidence type: '+type(item).__name__)
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False,default=convert)
def _sha(value):return hashlib.sha256(_json(value).encode()).hexdigest()
def _sum(value):return math.fsum(float(x) for x in np.asarray(value).ravel())
def _pos(pa):return np.column_stack((pa.x,pa.y,pa.z))
def _vel(pa):return np.column_stack((pa.u,pa.v,pa.w))


def _finite(value,name,shape=None):
    out=np.asarray(value,float)
    if not np.isfinite(out).all() or (shape is not None and out.shape!=shape):
        raise ValueError('finite matching '+name+' required')
    return out


def _near(actual,expected,name,scale=0.):
    a=np.asarray(actual,float);b=np.asarray(expected,float)
    if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('invalid '+name)
    bound=512*np.finfo(float).eps*max(float(np.max(np.abs(a),initial=0.)),
        float(np.max(np.abs(b),initial=0.)),float(scale),np.finfo(float).tiny)
    if np.any(np.abs(a-b)>bound):raise ValueError('inconsistent '+name)
    return bound


def _fields(pa):
    digest=_array_digest(pa);n=pa.get_number_of_particles()
    m=_finite(pa.m,'mass',(n,));rho=_finite(pa.rho,'density',(n,));T=_finite(pa.temperature,'temperature',(n,))
    S=_finite(pa.S,'stress',(9*n,)).reshape(n,3,3)
    if np.any(m<=0) or np.any(rho<=0) or np.any(T<=0):raise ValueError('positive mass/density/temperature required')
    return dict(digest=digest,ids=pa.particle_id.astype(int),m=m,rho=rho,T=T,S=S)


def _storage(fields,G,K,rho_reference):
    dev=fields['m']/fields['rho']*np.sum(fields['S']**2,axis=(1,2))/(4*G)
    bulk=fields['m']*volumetric_energy(fields['rho'],rho_reference,K)
    _finite(dev,'deviatoric storage');_finite(bulk,'volumetric storage')
    return dev+bulk


def evaluate_source_work(parent_pa,prepared_pa,actual_pa,source_event,*,parameters,free,result,
        mechanical_energy,parent_contact_storage_J,external_port):
    """Compute three source terms and an open-inventory QS energy identity.

    Result arrays and ``free`` align with prepared_pa; actual_pa can be sorted.
    Parent contact storage rows align with parent_pa: normal/tangential [J].
    This deliberately supports ONLY zero external traction and stress-free,
    reference-density newborn material with no initial contact storage.
    """
    B,P,A=(_fields(pa) for pa in (parent_pa,prepared_pa,actual_pa))
    n=len(P['ids']);parent_set=set(B['ids']);prepared_set=set(P['ids'])
    if not parent_set.issubset(prepared_set) or set(A['ids'])!=prepared_set:raise ValueError('birth inventory lost/replaced IDs')
    born_ids=sorted(prepared_set-parent_set)
    if not born_ids or born_ids!=source_event.get('allocated_particle_ids'):raise ValueError('source event birth IDs differ')
    if source_event.get('parent_particle_state_sha256')!=B['digest']:raise ValueError('source event parent PA binding differs')
    pi={int(pid):k for k,pid in enumerate(P['ids'])};ai={int(pid):k for k,pid in enumerate(A['ids'])}
    old_rows=np.array([pi[int(pid)] for pid in B['ids']]);born=np.array([pi[pid] for pid in born_ids])
    actual_rows=np.array([ai[int(pid)] for pid in P['ids']])
    if set(parent_pa.properties)!=set(prepared_pa.properties) or set(actual_pa.properties)!=set(prepared_pa.properties):
        raise ValueError('PA property schemas differ')
    if (set(parent_pa.constants)!=set(prepared_pa.constants) or set(actual_pa.constants)!=set(prepared_pa.constants)
            or any(not np.array_equal(parent_pa.constants[key],prepared_pa.constants[key])
                or not np.array_equal(actual_pa.constants[key],prepared_pa.constants[key]) for key in parent_pa.constants)):
        raise ValueError('birth or mechanical update changed PA constants')
    for name,prop in parent_pa.properties.items():
        stride=int(parent_pa.stride.get(name,1));before=prop.get_npy_array().reshape(len(B['ids']),stride)
        prepared=prepared_pa.properties[name].get_npy_array().reshape(n,stride)
        if name!='actuator' and not np.array_equal(before,prepared[old_rows]):
            raise ValueError('birth changed an old particle field: '+name)
    for name in ('m','x0','y0','z0','material_id','component_id','tracer_type','wholebar_tag','actuator'):
        if not np.array_equal(getattr(prepared_pa,name),getattr(actual_pa,name)[actual_rows]):
            raise ValueError('actual immutable identity/control changed: '+name)
    G=float(parameters['shear']);rho_reference=float(parameters['rho0']);cp=float(parameters['cp'])
    poisson=float(parameters['poisson'])
    if not -1<poisson<.5:raise ValueError('positive isotropic bulk modulus requires -1 < poisson < 0.5')
    K=float(bulk_modulus(G,poisson))
    reference=float(source_event['reference_temperature']['value'])
    if not np.isfinite([G,K,rho_reference,cp,reference]).all() or min(G,K,rho_reference,cp,reference)<=0:
        raise ValueError('positive finite moduli/reference density/cp/reference temperature required')
    if (source_event['heat_capacity']['value']!=cp or source_event['heat_capacity']['unit']!='J/(kg K)'
            or source_event['reference_temperature']['unit']!='K'):
        raise ValueError('incoming heat convention differs from coupled model')
    if np.any(P['S'][born]!=0.) or np.any(P['rho'][born]!=rho_reference):
        raise ValueError('source route requires exactly stress-free reference-density birth')
    for name in ('ep','plastic_work','return_dissipation'):
        if np.any(getattr(prepared_pa,name)[born]!=0.):raise ValueError('source route requires virgin birth history: '+name)
    if (external_port.get('contract')!=PORT_CONTRACT
            or external_port.get('surface')!='EXTERIOR_BIRTH_BOUNDARY_NOT_INTERNAL_CONTROL_PLANE'
            or external_port.get('internal_interface_work')!='NOT_ADDED_INTERNAL_PAIR_WORK'
            or not isinstance(external_port.get('source'),str) or not external_port['source'].strip()
            or not np.array_equal(_finite(external_port.get('traction_Pa'), 'declared exterior traction',(3,)),np.zeros(3))):
        raise ValueError('explicit zero-traction exterior birth port required; no internal-interface pressure work')
    f=np.asarray(free)
    if f.shape!=(n,3) or f.dtype.kind!='b':raise ValueError('boolean prepared free axes required')
    u=_finite(result['displacement'],'actual displacement',(n,3));force=_finite(result['force'],'actual force',(n,3))
    if not np.array_equal(_pos(prepared_pa)+u,_finite(result['position'],'result position',(n,3))):
        raise ValueError('actual result is not the declared displacement from prepared PA')
    for label,value,target in [('position',result['position'],_pos(actual_pa)[actual_rows]),
            ('density',result['density'],A['rho'][actual_rows]),('stress',result['stress'],A['S'][actual_rows]),
            ('temperature',result['temperature'],A['T'][actual_rows])]:
        if not np.array_equal(value,target):raise ValueError('actual PA/result mismatch: '+label)
    if not math.isfinite(float(result['dt_s'])) or float(result['dt_s'])<=0:raise ValueError('positive actual time interval required')
    volume=_finite(result['volume'],'result volume',(n,));_near(volume,P['m']/A['rho'][actual_rows],'actual volume')
    if not np.isin(prepared_pa.actuator,[0,1]).all():raise ValueError('binary inlet actuator mask required')
    actuator=np.asarray(prepared_pa.actuator).astype(bool)
    if np.any(f[actuator,2]) or np.any(~f[actuator,:2]):raise ValueError('inlet actuator must constrain only z')
    moving_axes=np.zeros((n,3),bool);moving_axes[actuator,2]=True
    if np.any(u[(~f)&(~moving_axes)]!=0.):raise ValueError('additional moving constrained axes unsupported by source-work split')
    reaction=-force*(~f);Wact=_sum(reaction[:,2]*u[:,2]*actuator)
    _near(_sum(reaction*u),Wact,'no double-counted or extra boundary work')
    contact=result['contact']
    def contact_rows(name):
        a=_finite(contact[name],name)
        if a.ndim!=2 or a.shape[0]!=n:raise ValueError('contact face rows must align with prepared PA')
        return a
    n0,n1,t0,t1=(contact_rows(key) for key in ('normal_potential0_J','normal_potential1_J',
        'tangential_potential0_J','tangential_potential1_J'))
    if np.any(n0[born]!=0.) or np.any(t0[born]!=0.):raise ValueError('newborn initial contact storage must be zero')
    parent_contact=_finite(parent_contact_storage_J,'parent contact storage',(len(B['ids']),2))
    _near(parent_contact,np.column_stack((n0.sum(axis=1)[old_rows],t0.sum(axis=1)[old_rows])),
        'old contact storage across birth')
    Cb=_sum(parent_contact);Ca=_sum(n1)+_sum(t1)
    normal_change=_sum(n1)-_sum(n0);tangent_change=_sum(t1)-_sum(t0)
    Wtool=_sum(contact_rows('wall_work_J'));friction=_sum(contact_rows('friction_dissipation_J'))
    plastic_density=_finite(result['material'][2],'plastic work density',(n,));Wp=_sum(volume*plastic_density)
    Ub,Up,Ua=(_storage(fields,G,K,rho_reference) for fields in (B,P,A))
    elastic_birth=_sum(Up[born])
    momentum=np.sum(P['m'][born,None]*_vel(prepared_pa)[born],axis=0)
    sensible_birth=_sum(P['m'][born]*cp*(P['T'][born]-reference))
    kinetic_birth=_sum(.5*P['m'][born]*np.sum(_vel(prepared_pa)[born]**2,axis=1))
    source=dict(mass_kg=_sum(P['m'][born]),momentum_kg_m_s=momentum.tolist(),sensible_heat_J=sensible_birth,
        kinetic_energy_J=kinetic_birth,volume_m3=_sum(P['m'][born]/P['rho'][born]))
    for key,value in source.items():_near(source_event['source_increment'][key],value,'inlet '+key)
    mechanical=dict(material_storage_increment_J=_sum(Ua)-_sum(Up),plastic_work_J=Wp,
        normal_storage_increment_J=normal_change,tangential_storage_increment_J=tangent_change,
        friction_sliding_J=friction,prescribed_boundary_work_J=Wact,tool_work_J=Wtool)
    mech_raw=Wact+Wtool-(mechanical['material_storage_increment_J']+normal_change+tangent_change+Wp+friction)
    mechanical['global_raw_residual_J']=mech_raw
    mech_scale=sum(abs(value) for key,value in mechanical.items() if key!='global_raw_residual_J')
    # These are consistency-roundoff bounds, NEVER process acceptance limits.
    # Core uses np.sum and separate deviatoric/bulk sums; this reader uses fsum.
    # Subtracting two large stored energies needs their absolute magnitudes in
    # the roundoff scale even when their actual step increment is tiny.
    storage_scales=dict(material_storage_increment_J=_sum(np.abs(Ua))+_sum(np.abs(Up)),
        normal_storage_increment_J=_sum(np.abs(n1))+_sum(np.abs(n0)),
        tangential_storage_increment_J=_sum(np.abs(t1))+_sum(np.abs(t0)))
    storage_scales['global_raw_residual_J']=sum(storage_scales.values())
    consistency=[]
    for key,value in mechanical.items():
        magnitude=mech_scale+storage_scales.get(key,0.)
        bound=_near(mechanical_energy[key],value,'actual mechanical ledger '+key,magnitude)
        consistency.append(dict(field=key,provided_J=float(mechanical_energy[key]),recomputed_J=float(value),
            absolute_roundoff_bound_J=bound,absolute_storage_terms_scale_J=storage_scales.get(key,0.),
            scale_basis='STEP_MAGNITUDES_PLUS_ABSOLUTE_SUBTRACTED_STORAGE_TERMS',acceptance_threshold=False))
    chi=float(parameters['plastic_heat_fraction']);partition=float(parameters['heat_partition'])
    if not 0<=chi<=1 or not 0<=partition<=1:raise ValueError('heat fractions outside [0,1]')
    thermal=result['thermal']['ledger'];contact_in=float(thermal['contact_in_J'])
    dT=_sum(P['m']*cp*(A['T'][actual_rows]-P['T']))
    thermo=dict(plastic_work_J=Wp,plastic_heat_J=chi*Wp,plastic_unthermalized_J=(1-chi)*Wp,
        friction_raw_J=friction,friction_material_heat_J=partition*friction,
        friction_wall_heat_J=(1-partition)*friction,measured_thermal_increment_J=dT,
        wall_heat_out_J=(1-partition)*friction-contact_in)
    thermal_scale=sum(abs(value) for value in thermo.values())+abs(contact_in)
    for key,value in thermo.items():_near(thermal[key],value,'actual thermal ledger '+key,thermal_scale)
    internal_conduction=float(thermal['internal_conduction_J']);local_thermal_residual=float(thermal['mc_delta_T_residual_J'])
    thermal_raw=dT-(chi*Wp+partition*friction+contact_in)
    _near(local_thermal_residual+internal_conduction,thermal_raw,'global versus local thermal balance',thermal_scale)
    ETb=_sum(B['m']*cp*(B['T']-reference));ETa=_sum(A['m']*cp*(A['T']-reference))
    delta_total_storage=_sum(Ua)-_sum(Ub)+ETa-ETb+Ca-Cb
    raw_open=Wact+Wtool+elastic_birth+sensible_birth-delta_total_storage-thermo['plastic_unthermalized_J']-thermo['wall_heat_out_J']
    identity_residual=raw_open-(mech_raw-thermal_raw)
    identity_scale=sum(abs(v) for v in (Wact,Wtool,elastic_birth,sensible_birth,delta_total_storage,
        thermo['plastic_unthermalized_J'],thermo['wall_heat_out_J'],ETa,ETb,_sum(Ua),_sum(Ub),Ca,Cb))
    identity_bound=_near(raw_open,mech_raw-thermal_raw,'open-inventory energy identity',identity_scale)
    evidence=dict(parent_pa_sha256=B['digest'],prepared_pa_sha256=P['digest'],actual_pa_sha256=A['digest'],
        source_event_sha256=_sha(source_event),external_port_sha256=_sha(external_port),parameters_sha256=_sha(parameters),
        force_displacement_free_sha256=_sha(dict(force_N=force.tolist(),displacement_m=u.tolist(),free=f.tolist())),
        mechanical_energy_sha256=_sha(mechanical_energy),thermal_ledger_sha256=_sha(thermal),
        parent_contact_storage_sha256=_sha(parent_contact.tolist()),
        contact_energy_rows_sha256=_sha({key:contact[key] for key in ('normal_potential0_J','normal_potential1_J',
            'tangential_potential0_J','tangential_potential1_J','wall_work_J','friction_dissipation_J')}),
        material_volume_work_dt_sha256=_sha(dict(volume_m3=volume,plastic_work_density_Pa=plastic_density,dt_s=float(result['dt_s']))))
    report=dict(classification=CLASSIFICATION,status='COMPUTED_SOURCE_WORK_NOT_PROCESS_ACCEPTANCE',
        physical_solver_evaluations=0,nnps_evaluations=0,force_evaluations=0,thermal_evaluations=0,
        born_particle_ids=born_ids,source_increment=source,external_port=external_port,
        thermal_inventory_reference=dict(reference_temperature_K=reference,heat_capacity_J_kg_K=cp,
            convention='ONE_REFERENCE_FOR_PARENT_PREPARED_ACTUAL_AND_BIRTH_INVENTORIES',
            cumulative_history_requires_same_reference=True),
        computed_terms=dict(pressure_flow_work=dict(status='COMPUTED',unit='J',value=0.,basis='DECLARED_ZERO_EXTERNAL_BIRTH_TRACTION'),
            elastic_convected_energy=dict(status='COMPUTED',unit='J',value=elastic_birth,basis='ACTUAL_BIRTH_MATERIAL_STORAGE'),
            mechanical_reaction_work=dict(status='COMPUTED',unit='J',value=Wact,basis='ACTUAL_COUPLED_REACTION_DOT_DISPLACEMENT_ONCE')),
        open_energy=dict(material_storage_before_J=_sum(Ub),material_storage_after_J=_sum(Ua),
            sensible_storage_before_J=ETb,sensible_storage_after_J=ETa,contact_storage_before_J=Cb,contact_storage_after_J=Ca,
            total_storage_increment_J=delta_total_storage,tool_work_J=Wtool,
            plastic_unthermalized_J=thermo['plastic_unthermalized_J'],wall_heat_out_J=thermo['wall_heat_out_J'],
            mechanical_raw_residual_J=mech_raw,thermal_external_raw_residual_J=thermal_raw,
            provided_mechanical_raw_residual_J=float(mechanical_energy['global_raw_residual_J']),
            mechanical_recomputed_minus_provided_J=mech_raw-float(mechanical_energy['global_raw_residual_J']),
            internal_conduction_sum_J=internal_conduction,local_thermal_residual_sum_J=local_thermal_residual,
            open_raw_residual_J=raw_open,algebra_identity_residual_J=identity_residual,algebra_roundoff_bound_J=identity_bound),
        mechanical_ledger_consistency_roundoff=consistency,
        quasistatic_kinetic_diagnostic=dict(source_kinetic_energy_J=kinetic_birth,
            inventory_change_J=_sum(.5*A['m']*np.sum(_vel(actual_pa)**2,axis=1))-_sum(.5*B['m']*np.sum(_vel(parent_pa)**2,axis=1)),
            included_in_quasistatic_energy=False,converted_to_heat=False),input_sha256=evidence,
        limitations=['Computed source accounting and algebra only; existing coupled and interval gates remain unchanged.',
            'Stress-free, reference-density, zero-traction birth port only; no arbitrary pressure-driven inlet is implemented.',
            'Original endpoint-material raw defects are retained, not corrected or turned into heat.',
            'Reported QS storage diagnostic is not a finite-strain free-energy validation.',
            'A full cumulative open-process ledger and real two-insertion benchmark remain caller responsibilities.'])
    report=json.loads(_json(report))
    report['calculation_sha256']=_sha(report)
    return report


def source_work_record(report,artifact_path):
    """Read back an immutable computed report, then bind V2 work to its file SHA."""
    path=Path(artifact_path);payload=path.read_bytes();stored=json.loads(payload)
    if stored!=report:raise ValueError('source-work artifact content differs from supplied calculation')
    unsealed=dict(stored);expected=unsealed.pop('calculation_sha256',None)
    if expected!=_sha(unsealed) or stored.get('classification')!=CLASSIFICATION:
        raise ValueError('source-work calculation seal/classification mismatch')
    evidence=hashlib.sha256(payload).hexdigest();out=dict(convention=WORK_CONVENTION)
    for key in WORK_TERMS:
        row=stored['computed_terms'][key]
        out[key]=dict(status=row['status'],unit=row['unit'],value=row['value'],evidence_sha256=evidence)
    return out

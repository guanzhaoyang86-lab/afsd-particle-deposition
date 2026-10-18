"""No-birth source work from actual saved fields, not a fabricated inlet event."""
import json
import hashlib
from pathlib import Path
import numpy as np
import inlet_source_work as birth
from inlet_acceptance_v2 import WORK_CONVENTION,WORK_TERMS
from afsd_backend.pysph.mechanics import bulk_modulus

CLASSIFICATION='FIXED_INVENTORY_ZERO_EXTERIOR_TRACTION_SOURCE_WORK_DIAGNOSTIC_V1'


def evaluate(parent,prepared,actual,*,parameters,free,result,mechanical_energy,parent_contact_storage_J,external_port):
    """Zero mass/enthalpy/advection input only after identical inventory checks.

    Actual mechanical reaction work is computed, never zero-filled. The
    explicit exterior traction is zero; no internal control-plane work is
    added a second time. No particle, force, NNPS or heat solve is performed.
    """
    B,P,A=(birth._fields(pa) for pa in (parent,prepared,actual));n=len(P['ids'])
    if set(B['ids'])!=set(P['ids']) or set(A['ids'])!=set(P['ids']):raise ValueError('fixed-inventory route forbids birth/removal/replacement')
    pi={int(pid):k for k,pid in enumerate(P['ids'])};ai={int(pid):k for k,pid in enumerate(A['ids'])}
    oldrows=np.array([pi[int(pid)] for pid in B['ids']]);actualrows=np.array([ai[int(pid)] for pid in P['ids']])
    if set(parent.properties)!=set(prepared.properties) or set(actual.properties)!=set(prepared.properties):raise ValueError('property schema changed')
    for key,prop in parent.properties.items():
        stride=int(parent.stride.get(key,1));old=prop.get_npy_array().reshape(n,stride)
        new=prepared.properties[key].get_npy_array().reshape(n,stride)
        if not np.array_equal(old,new[oldrows]):raise ValueError('fixed-inventory preparation changed old field: '+key)
    if (set(parent.constants)!=set(prepared.constants) or set(actual.constants)!=set(prepared.constants) or
        any(not np.array_equal(parent.constants[k],prepared.constants[k]) or not np.array_equal(actual.constants[k],prepared.constants[k]) for k in parent.constants)):
        raise ValueError('PA constants changed')
    for key in ('m','x0','y0','z0','material_id','component_id','tracer_type','wholebar_tag','actuator'):
        if not np.array_equal(getattr(prepared,key),getattr(actual,key)[actualrows]):raise ValueError('actual immutable field changed: '+key)
    p=parameters;G=float(p['shear']);rho0=float(p['rho0']);cp=float(p['cp']);reference=float(p['room'])
    K=bulk_modulus(G,p['poisson'])
    if not np.isfinite([G,K,rho0,cp,reference]).all() or min(G,K,rho0,cp,reference)<=0:raise ValueError('invalid physical inventory constants')
    if not np.isfinite(float(result['dt_s'])) or float(result['dt_s'])<=0:raise ValueError('positive actual interval required')
    if (external_port.get('contract')!=birth.PORT_CONTRACT or external_port.get('surface')!='EXTERIOR_BIRTH_BOUNDARY_NOT_INTERNAL_CONTROL_PLANE' or
        external_port.get('internal_interface_work')!='NOT_ADDED_INTERNAL_PAIR_WORK' or not external_port.get('source') or
        not np.array_equal(birth._finite(external_port.get('traction_Pa'),'exterior traction',(3,)),np.zeros(3))):
        raise ValueError('explicit zero exterior traction required')
    f=np.asarray(free);u=birth._finite(result['displacement'],'displacement',(n,3));force=birth._finite(result['force'],'force',(n,3))
    if f.shape!=(n,3) or f.dtype.kind!='b':raise ValueError('boolean free axes required')
    if not np.array_equal(birth._pos(prepared)+u,result['position']):raise ValueError('actual displacement endpoint differs')
    for key,value in (('position',birth._pos(actual)[actualrows]),('stress',A['S'][actualrows]),('density',A['rho'][actualrows]),('temperature',A['T'][actualrows])):
        if not np.array_equal(result[key],value):raise ValueError('actual PA/result differs: '+key)
    if not np.isin(prepared.actuator,[0,1]).all():raise ValueError('binary physical control required')
    actuator=prepared.actuator.astype(bool);moving=np.zeros((n,3),bool);moving[actuator,2]=True
    if np.any(f[actuator,2]) or np.any(~f[actuator,:2]) or np.any(u[(~f)&(~moving)]!=0):raise ValueError('unregistered moving constraint')
    Wact=birth._sum(-force*(~f)*u);contact=result['contact']
    def face(key):
        a=birth._finite(contact[key],key)
        if a.shape!=(n,8):raise ValueError('eight actual contact faces required')
        return a
    n0,n1,t0,t1=(face(key) for key in ('normal_potential0_J','normal_potential1_J','tangential_potential0_J','tangential_potential1_J'))
    parent_contact=birth._finite(parent_contact_storage_J,'parent contact',(n,2))
    birth._near(parent_contact,np.column_stack((n0.sum(axis=1),t0.sum(axis=1)))[oldrows],'retained contact inventory')
    V=birth._finite(result['volume'],'volume',(n,));birth._near(V,P['m']/A['rho'][actualrows],'actual volume')
    Ub,Ua=(birth._storage(fields,G,K,rho0) for fields in (B,A))
    dU=birth._sum(Ua)-birth._sum(Ub);dN=birth._sum(n1)-birth._sum(n0);dC=birth._sum(t1)-birth._sum(t0)
    Wp=birth._sum(V*result['material'][2]);Wf=birth._sum(face('friction_dissipation_J'));Wt=birth._sum(face('wall_work_J'))
    raw=Wact+Wt-dU-dN-dC-Wp-Wf
    recomputed=dict(material_storage_increment_J=dU,normal_storage_increment_J=dN,tangential_storage_increment_J=dC,
        plastic_work_J=Wp,friction_sliding_J=Wf,prescribed_boundary_work_J=Wact,tool_work_J=Wt,global_raw_residual_J=raw)
    storage_scale=birth._sum(np.abs(Ua))+birth._sum(np.abs(Ub))+birth._sum(np.abs(n0))+birth._sum(np.abs(n1))+birth._sum(np.abs(t0))+birth._sum(np.abs(t1))
    ledger_scale=storage_scale+sum(abs(value) for value in recomputed.values());bounds=[]
    for key,value in recomputed.items():
        bound=birth._near(mechanical_energy[key],value,'actual fixed-inventory mechanical '+key,ledger_scale)
        bounds.append(dict(field=key,roundoff_bound_J=bound,acceptance_threshold=False))
    thermal=result['thermal']['ledger'];dT=birth._sum(P['m']*cp*(A['T'][actualrows]-P['T']))
    chi=p['plastic_heat_fraction'];eta=p['heat_partition'];contact_in=float(thermal['contact_in_J'])
    if not 0<=chi<=1 or not 0<=eta<=1:raise ValueError('heat partition outside physical range')
    thermo=dict(plastic_work_J=Wp,plastic_heat_J=chi*Wp,plastic_unthermalized_J=(1-chi)*Wp,
        friction_raw_J=Wf,friction_material_heat_J=eta*Wf,friction_wall_heat_J=(1-eta)*Wf,
        measured_thermal_increment_J=dT,wall_heat_out_J=(1-eta)*Wf-contact_in)
    for key,value in thermo.items():birth._near(thermal[key],value,'actual fixed-inventory thermal '+key,sum(abs(x) for x in thermo.values()))
    thermal_raw=dT-(chi*Wp+eta*Wf+contact_in)
    birth._near(thermal_raw,thermal['mc_delta_T_residual_J']+thermal['internal_conduction_J'],'local/global heat identity',abs(dT)+abs(contact_in))
    ETb=birth._sum(B['m']*cp*(B['T']-reference));ETa=birth._sum(A['m']*cp*(A['T']-reference))
    Cb=birth._sum(parent_contact);Ca=birth._sum(n1)+birth._sum(t1)
    delta=birth._sum(Ua)-birth._sum(Ub)+ETa-ETb+Ca-Cb
    raw_open=Wact+Wt-delta-thermo['plastic_unthermalized_J']-thermo['wall_heat_out_J']
    bound=birth._near(raw_open,raw-thermal_raw,'fixed-inventory open balance identity',
        abs(ETa)+abs(ETb)+storage_scale+sum(abs(x) for x in (Wact,Wt,delta,*thermo.values())))
    zero=dict(mass_kg=0.,momentum_kg_m_s=[0.,0.,0.],sensible_heat_J=0.,kinetic_energy_J=0.,volume_m3=0.)
    terms=dict(pressure_flow_work=dict(status='COMPUTED',unit='J',value=0.,basis='EXPLICIT_ZERO_EXTERNAL_TRACTION_NO_BIRTH_INTERVAL'),
        elastic_convected_energy=dict(status='COMPUTED',unit='J',value=0.,basis='IDENTICAL_PERSISTENT_INVENTORY_NO_ADVECTION_SOURCE'),
        mechanical_reaction_work=dict(status='COMPUTED',unit='J',value=Wact,basis='ACTUAL_FORCE_DOT_PRESCRIBED_DISPLACEMENT_ONCE'))
    report=dict(classification=CLASSIFICATION,status='COMPUTED_FIXED_INVENTORY_SOURCE_WORK_NOT_PROCESS_ACCEPTANCE',
        source_increment=zero,birth_event=None,inventory_change_kind='NO_BIRTH_REMOVAL_OR_ID_REPLACEMENT',computed_terms=terms,
        thermal_inventory_reference=dict(reference_temperature_K=reference,heat_capacity_J_kg_K=cp),external_port=external_port,
        open_energy=dict(material_storage_before_J=birth._sum(Ub),material_storage_after_J=birth._sum(Ua),
            sensible_storage_before_J=ETb,sensible_storage_after_J=ETa,contact_storage_before_J=Cb,contact_storage_after_J=Ca,
            total_storage_increment_J=delta,tool_work_J=Wt,plastic_unthermalized_J=thermo['plastic_unthermalized_J'],wall_heat_out_J=thermo['wall_heat_out_J'],
            mechanical_raw_residual_J=raw,thermal_external_raw_residual_J=thermal_raw,open_raw_residual_J=raw_open,
            algebra_identity_residual_J=raw_open-(raw-thermal_raw),algebra_roundoff_bound_J=bound),
        mechanical_consistency_roundoff=bounds,physical_evaluations=0,
        quasistatic_kinetic_diagnostic=dict(source_kinetic_energy_J=0.,
            inventory_change_J=birth._sum(.5*A['m']*np.sum(birth._vel(actual)**2,axis=1))-
                birth._sum(.5*B['m']*np.sum(birth._vel(parent)**2,axis=1)),
            included_in_quasistatic_energy=False,converted_to_heat=False),
        input_sha256=dict(parent_pa_sha256=B['digest'],prepared_pa_sha256=P['digest'],actual_pa_sha256=A['digest'],
            parameters_sha256=birth._sha(p),mechanical_energy_sha256=birth._sha(mechanical_energy),thermal_ledger_sha256=birth._sha(thermal),
            actual_force_displacement_sha256=birth._sha(dict(force=force,displacement=u,free=f))))
    report=json.loads(birth._json(report));report['calculation_sha256']=birth._sha(report);return report


def record(report,path):
    data=Path(path).read_bytes();stored=json.loads(data);unsealed=dict(stored);expected=unsealed.pop('calculation_sha256',None)
    if stored!=report or expected!=birth._sha(unsealed) or stored['classification']!=CLASSIFICATION:raise ValueError('fixed-inventory source artifact differs')
    result=dict(convention=WORK_CONVENTION);digest=hashlib.sha256(data).hexdigest()
    for key in WORK_TERMS:
        row=report['computed_terms'][key];result[key]=dict(status=row['status'],unit='J',value=row['value'],evidence_sha256=digest)
    return result

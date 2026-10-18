"""Unified centre-collocated rigid contacts and original-material axial drive.

Uses the existing finite open-square feature geometry and coupled event kick.
This is NOT a complete rigid exterior/edge closest-point model. Contact point
projections are diagnostic; forces and wall velocities stay collocated at x_i.
`time` is the caller's frozen kick geometry/control evaluation time. The caller
must split integration steps at actuator programme discontinuities.

contact_power_W is ONLY signed thermal CONDUCTION from reservoirs into material.
It does not include friction heat; the thermal update alone partitions raw Q.
No tracer/AS/RS or current-position actuator selector enters this interface.
"""
import numpy as np
from .rigid_contact_kinematics import contact_rows
from .contact_kick import kick
from .axial_actuator import control


class BoundaryEvaluationError(RuntimeError):
    """Uncommitted boundary candidate, with original particle/feature context."""


def normal_snapshot(x, particle_ids, time, p):
    """Actual active finite features and raw penalty potential, not a cohort.

    Feature disappearance can change this potential discontinuously at an
    unsupported finite edge. The caller must retain before/after identities
    and report the normal temporal residual; it must never become heat.
    """
    stiffness=float(p['normal_stiffness'])
    if not np.isfinite(stiffness) or stiffness<=0:
        raise ValueError('positive finite normal stiffness required')
    rows=contact_rows(x,particle_ids,time,p)
    normal=stiffness*rows['penetration']
    potential=.5*stiffness*rows['penetration']**2
    positive=normal>0.
    return dict(rows=rows,normal_force_N=normal,potential_by_contact_J=potential,
        normal_potential_J=float(np.sum(potential)),
        max_penetration_m=float(np.max(rows['penetration'],initial=0.)),
        active_contact_ids=[str(value) for value in rows['contact_id'][positive]],
        potential_by_contact_id_J={str(cid):float(value) for cid,value in zip(rows['contact_id'],potential)},
        classification='FINITE_FEATURE_CENTRE_COLLOCATED_APPROXIMATE_NOT_COMPLETE_EDGES')


def _array(value,shape,name):
    result=np.asarray(value,float)
    if result.shape!=shape or not np.isfinite(result).all():
        raise ValueError('invalid/nonfinite '+name)
    return result


def evaluate_boundary(x,v,mass,T,time,dt,p,force,z0,particle_ids):
    """Return an uncommitted all-particle kick, raw work and conduction flux.

    Every particle calls the same constrained kick, even with zero contacts.
    All x/y DOFs remain free; only the selected initial top layer's z is driven.
    Each simultaneous face shares the old kick's projected mobility P/m.
    """
    x=np.asarray(x,float)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():raise ValueError('invalid particle positions')
    count=len(x);v=_array(v,(count,3),'velocities');mass=_array(mass,(count,),'mass')
    temperature=_array(T,(count,),'temperature');force=_array(force,(count,3),'internal force')
    z0=_array(z0,(count,),'original z')
    if np.any(mass<=0) or not np.isfinite(dt) or dt<=0:raise ValueError('positive mass and dt required')
    mu=float(p['friction']);conductance=float(p['contact_conductance'])
    if not np.isfinite([mu,conductance]).all() or mu<0 or conductance<0:
        raise ValueError('finite nonnegative friction and contact conductance required')
    snapshot=normal_snapshot(x,particle_ids,time,p);rows=snapshot['rows']
    command=control(z0,time,0.,p['bar_length'],p['actuator_depth'],p['feed'],p['release_time'])
    after=np.empty_like(v);work_velocity=np.empty_like(v);q_particle=np.zeros(count)
    conduction=np.zeros(count);thermal_row=np.zeros(count);traces=[]
    ledger=dict(friction_raw_J=0.,normal_material_work_J=0.,normal_wall_work_J=0.,
        normal_relative_work_J=0.,wall_work_J=0.,actuator_work_J=0.,internal_work_J=0.,
        delta_kinetic_J=0.,kinetic_work_residual_J=0.,contact_kinetic_residual_J=0.,
        contact_conduction_into_material_J=0.,wall_conduction_out_J=0.,
        normal_potential_J=snapshot['normal_potential_J'],max_penetration_m=snapshot['max_penetration_m'],
        contacts=int(np.count_nonzero(snapshot['normal_force_N']>0.)),contact_subintervals=0,
        contact_rejected_subintervals=0,contact_zero_events=0,minimum_raw_friction_J=None)
    for body in ('base','tool','actuator'):
        ledger[body+'_impulse_Ns']=np.zeros(3)
        ledger[body+'_angular_impulse_Nms']=np.zeros(3)
    for index,pid in enumerate(np.asarray(particle_ids)):
        selected=np.flatnonzero((rows['particle_index']==index)&(snapshot['normal_force_N']>0.))
        faces=[{key:value if isinstance(value,str) else np.asarray(value)[r].copy()
                for key,value in rows.items()} for r in selected]
        normals=rows['normal'][selected];walls=rows['wall_velocity'][selected]
        normal_forces=snapshot['normal_force_N'][selected]
        wall_temperatures=np.where(rows['feature_id'][selected]==0,p['base_temperature'],p['tool_temperature']).astype(float)
        for a,face in enumerate(faces):
            face['normal_force_N']=float(normal_forces[a]);face['normal_potential_J']=float(snapshot['potential_by_contact_J'][selected[a]])
            face['wall_temperature_K']=float(wall_temperatures[a])
        free=command['free_axes'][index];prescribed=command['prescribed_velocity'][index]
        try:
            result=kick(v[index],mass[index],force[index],normals,normal_forces,walls,mu,dt,
                free_axes=free,prescribed_velocity=prescribed)
        except Exception as error:
            failure=BoundaryEvaluationError('particle_id='+str(int(pid))+' contacts='+
                ','.join(str(face['contact_id']) for face in faces)+'; '+str(error))
            failure.contact_context=dict(classification='REJECTED_BOUNDARY_CANDIDATE',restart_allowed=False,
                particle_id=int(pid),particle_index=index,contact_ids=[str(face['contact_id']) for face in faces],
                feature_ids=rows['feature_id'][selected].copy(),faces=faces,position_m=x[index].copy(),
                velocity_before_m_s=v[index].copy(),internal_force_N=force[index].copy(),mass_kg=float(mass[index]),
                time_mid_s=float(time),dt_s=float(dt),free_axes=free.copy(),prescribed_velocity=prescribed.copy(),
                actuator_active=bool(command['active'][index]),normal_forces_N=normal_forces.copy(),
                wall_velocities_m_s=walls.copy(),completed_particle_traces_in_rejected_candidate=traces,
                failed_cone=getattr(error,'cone_context',None),failed_interval=getattr(error,'interval_context',None))
            raise failure from error
        after[index]=result['velocity']
        work_velocity[index]=sum(a['dt_s']*.5*(a['v_before']+a['v_after']) for a in result['intervals'])/dt
        q_particle[index]=float(np.sum(result['friction_heat_J']))
        face_h=conductance*p['dx']**2
        flux=face_h*(wall_temperatures-temperature[index])
        conduction[index]=float(np.sum(flux));thermal_row[index]=len(selected)*face_h
        normal_material=0.;normal_wall=0.
        for interval in result['intervals']:
            normal_material+=float(np.sum(interval['normal_work_J']))
            normal_wall+=float(np.sum(interval['normal_impulse']*walls))
            raw=np.asarray(interval['friction_heat_raw_J'])
            if raw.size:
                minimum=float(raw.min());previous=ledger['minimum_raw_friction_J']
                ledger['minimum_raw_friction_J']=minimum if previous is None else min(previous,minimum)
        for name,value in [('friction_raw_J',q_particle[index]),('normal_material_work_J',normal_material),
            ('normal_wall_work_J',normal_wall),('normal_relative_work_J',normal_material-normal_wall),
            ('wall_work_J',result['wall_work_J']),('actuator_work_J',result['constraint_work_J']),
            ('internal_work_J',result['other_work_J']),('delta_kinetic_J',result['delta_kinetic_J']),
            ('contact_kinetic_residual_J',result['kinetic_work_residual_J'])]:ledger[name]+=float(value)
        for body,mask in [('base',rows['feature_id'][selected]==0),('tool',rows['feature_id'][selected]!=0)]:
            impulse=np.sum((result['normal_impulse']+result['tangential_impulse'])[mask],axis=0)
            ledger[body+'_impulse_Ns']+=impulse
            ledger[body+'_angular_impulse_Nms']+=np.cross(x[index],impulse)
        ledger['actuator_impulse_Ns']+=result['constraint_impulse']
        ledger['actuator_angular_impulse_Nms']+=np.cross(x[index],result['constraint_impulse'])
        ledger['contact_subintervals']+=len(result['intervals'])
        ledger['contact_rejected_subintervals']+=int(result['rejected_intervals'])
        ledger['contact_zero_events']+=len(result['events'])
        traces.append(dict(particle_id=int(pid),particle_index=index,stage='MIDPOINT_FROZEN_GEOMETRY_KICK',
            time_mid_s=float(time),dt_s=float(dt),position_m=x[index].copy(),mass_kg=float(mass[index]),
            material_velocity_before_m_s=v[index].copy(),internal_force_N=force[index].copy(),
            free_axes=free.copy(),prescribed_velocity=prescribed.copy(),actuator_active=bool(command['active'][index]),
            actuator_force_on_material_N=result['constraint_impulse']/dt,contact_rows=faces,result=result,
            contact_conduction_W=float(conduction[index]),wall_row_W_K=float(thermal_row[index])))
    ledger['contact_conduction_into_material_J']=float(dt*np.sum(conduction))
    ledger['wall_conduction_out_J']=-ledger['contact_conduction_into_material_J']
    ledger['kinetic_work_residual_J']=(ledger['delta_kinetic_J']+ledger['friction_raw_J']-
        ledger['wall_work_J']-ledger['actuator_work_J']-ledger['internal_work_J']-ledger['normal_relative_work_J'])
    return dict(velocity=after,work_velocity=work_velocity,traces=traces,command=command,rows=rows,
        friction_raw_J=q_particle,contact_power_W=conduction,contact_conduction_W=conduction.copy(),
        wall_row_W_K=thermal_row,ledger=ledger,normal_snapshot=snapshot,
        contact_power_definition='SIGNED_CONTACT_CONDUCTION_ONLY; FRICTION_PARTITION_BELONGS_TO_THERMAL_UPDATE')

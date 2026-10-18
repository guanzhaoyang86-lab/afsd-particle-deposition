"""Same manufactured rigid planes with explicit one-sided event callbacks.

No AFSD geometry or free-material trajectory is prescribed by this callback.
Outside selected zero events the frozen slow callback formula is unchanged.
"""
import copy
import numpy as np
from afsd_backend.pysph.mechanics import velocities


def boundary_for(p,event_context=None):
    context=copy.deepcopy(event_context or {})
    def boundary(pa,t):
        angle=2*np.pi*p['normal_frequency_Hz']*t
        stroke=p['normal_stroke']*np.sin(angle/2)**2
        speed=p['normal_stroke']*np.pi*p['normal_frequency_Hz']*np.sin(angle)
        planes=[-.5*p['dx']+p['normal_preload']+stroke,
            p['length']+.5*p['dx']-p['normal_preload']-stroke]
        vv=velocities(pa);rows=[];seen=set()
        for face in (0,1):
            sign=1. if face==0 else -1.;normal=np.array([sign,0.,0.])
            wall=np.array([sign*speed,0.,p['wall_speed']])
            for index in np.flatnonzero(pa.x0==(0. if face==0 else p['length'])):
                delta=.5*p['dx']+sign*(planes[face]-pa.x[index])
                relative=vv[index]-wall;tangent=relative-(relative@normal)*normal;norm=np.linalg.norm(tangent)
                cid=f'{face}:{int(pa.particle_id[index])}';seen.add(cid)
                if delta<=0:raise ValueError('OPENING_BOUNDARY '+cid)
                fn=p['normal_stiffness']*delta
                if cid in context:
                    entry=context[cid];kind=entry['endpoint'];incoming=np.asarray(entry['incoming_direction'])
                    if kind not in ('left_zero','right_zero','private_continuation'):
                        raise ValueError('unknown event callback mode')
                    limit=np.asarray(entry['outgoing_direction']) if kind=='right_zero' else incoming
                    bound=512*np.finfo(float).eps*max(1.,np.linalg.norm(vv[index]),np.linalg.norm(wall))
                    if norm<=bound:
                        ft=-p['friction']*fn*limit
                    else:
                        ft=-p['friction']*fn*tangent/norm
                        if kind=='private_continuation' and tangent@incoming<0:ft=-ft
                else:
                    if norm==0:raise ValueError('EXACT_ZERO_SLIP_BOUNDARY '+cid)
                    ft=-p['friction']*fn*tangent/norm
                rows.append(dict(particle_index=int(index),contact_id=cid,normal=normal,
                    normal_force_N=fn,tangential_force_N=ft,wall_velocity_m_s=wall,
                    wall_temperature_K=p['wall_temperature_K'][face],
                    conductance_W_K=p['contact_conductance']*p['dx']**2,penetration_m=delta,
                    normal_potential_J=.5*p['normal_stiffness']*delta*delta,normal_stiffness_N_m=p['normal_stiffness']))
        if set(context)-seen:raise ValueError('selected event identity not in actual active contacts')
        return dict(force_definition='NONCONTACT_APPLIED_FORCE_ONLY',force_N=np.zeros((len(pa.x),3)),contacts=rows)
    return boundary


def actual_wall_accelerations(p,t,contacts):
    scalar=p['normal_stroke']*2*np.pi**2*p['normal_frequency_Hz']**2*np.cos(2*np.pi*p['normal_frequency_Hz']*t)
    return {row['contact_id']:row['normal']*scalar for row in contacts}

"""APPROXIMATE frozen-geometry event-aware coupled Coulomb kick.
No geometry/material/NNPS changes. Force and wall velocity are collocated at
particle centres. Normal penalty and internal forces are frozen for this kick.
All simultaneous contacts share projected mobility P/m; no independent mass cap.
"""
import numpy as np

class ContactConvergenceError(RuntimeError):
    pass

def tangent_bases(normals):
    """Orthonormal columns span each plane without subtracting n*n^T from I."""
    normals = np.asarray(normals, float).reshape(-1, 3)
    bases = np.empty((len(normals), 3, 2))
    for k, n in enumerate(normals):
        unit = n/np.linalg.norm(n)
        axis = np.eye(3)[int(np.argmin(np.abs(unit)))]
        first = np.cross(unit, axis)
        first /= np.linalg.norm(first)
        bases[k, :, 0] = first
        bases[k, :, 1] = np.cross(unit, first)
    return bases

def cone_solve_projected_legacy(base, bases, mobility, coefficient, limits, tolerance=2e-12):
    """Coupled disk solve in 2D tangent coordinates; returns 3D face forces.

    base is already tangent (N,2). The Hessian blocks are coefficient times
    T_a^T mobility T_b, including off-diagonal faces and constrained DOFs.
    No large normal base velocity is repeatedly added then projected away.
    """
    count = len(limits)
    if count == 0: return np.empty((0,3))
    f = np.zeros((count,2))
    lipschitz = coefficient*count*np.linalg.norm(mobility, 2)
    if lipschitz == 0:
        # No unconstrained response: friction still opposes prescribed slip.
        for k in range(count):
            u = base[k]
            speed = np.linalg.norm(u)
            if speed: f[k] = -limits[k]*u/speed
        return np.einsum('nki,ni->nk', bases, f)
    scale = max(float(np.max(limits)), np.finfo(float).tiny)
    last_updates = []
    for iteration in range(2000):
        resultant = np.einsum('nki,ni->k', bases, f)
        common = coefficient*mobility@resultant
        gradient = base+np.einsum('nki,k->ni', bases, common)
        trial = f-gradient/lipschitz
        norm = np.linalg.norm(trial, axis=1)
        outside = norm > limits
        trial[outside] *= (limits[outside]/norm[outside])[:,None]
        change = float(np.linalg.norm(trial-f, axis=1).max())
        if iteration >= 1996:
            last_updates.append(dict(iteration=iteration+1, change_N=change,
                                     force_before=f.copy(), force_after=trial.copy()))
        if change <= tolerance*scale:
            return np.einsum('nki,ni->nk', bases, trial)
        f = trial
    failure = ContactConvergenceError('coupled Coulomb cone iteration did not converge')
    failure.cone_context = dict(base_tangent=base,tangent_bases=bases,mobility=mobility,
        coefficient=coefficient,limits_N=limits,threshold_N=tolerance*scale,
        lipschitz=lipschitz,last_updates=last_updates)
    raise failure

def cone_solve(base,bases,mobility,coefficient,limits,tolerance=2e-12):
    """Same convex problem, acceptance independent of the chosen force split."""
    from .contact_cone import solve,ConeCertificateError
    try:
        return solve(base,bases,mobility,coefficient,limits,tolerance)
    except ConeCertificateError as error:
        failure=ContactConvergenceError(str(error));failure.cone_context=error.context
        raise failure from error


def kick(v0, mass, other_force, normals, normal_forces, wall_velocities,
         mu, dt, free_axes=None, prescribed_velocity=None,
         rtol=2e-7, atol=1e-12, max_intervals=4096):
    """Return impulses, subinterval work and endpoint; no particle mutation.
    Exact straight-slip stop event; static feasibility then reverse sliding.
    Curved slip: adaptive midpoint, step-doubling, angle guard. Not a global
    implicit method. No negative heat removal; rejected intervals are uncommitted.
    """
    v0 = np.asarray(v0, float); external = np.asarray(other_force, float)
    n = np.asarray(normals, float).reshape(-1,3)
    fn = np.asarray(normal_forces, float)
    walls = np.asarray(wall_velocities, float).reshape(-1,3)
    nc = len(fn)
    if mass <= 0 or dt <= 0 or mu < 0 or np.any(fn < 0):
        raise ValueError('invalid contact kick input')
    if n.shape != walls.shape or n.shape != (nc,3):
        raise ValueError('contact array mismatch')
    if nc and np.max(np.abs(np.linalg.norm(n,axis=1)-1)) > 1e-12:
        raise ValueError('unit contact normals required')
    free = np.ones(3) if free_axes is None else np.asarray(free_axes, float)
    if np.any((free != 0)&(free != 1)): raise ValueError('free axes must be 0/1')
    P = np.diag(free)
    prescribed = v0 if prescribed_velocity is None else np.asarray(prescribed_velocity,float)
    if np.linalg.norm((1-free)*(v0-prescribed)) > 1e-12:
        raise ValueError('inconsistent prescribed initial velocity; explicit actuator event required')
    # Translating wall frame: common lab boosts must not change slip thresholds.
    # This changes numerical coordinates, not the physical trajectory.
    frame=walls[0].copy() if nc else v0.copy()
    v0=v0-frame;walls=walls-frame
    bases = tangent_bases(n)
    R = bases@bases.swapaxes(1,2)
    limits = mu*fn
    normal_vectors = fn[:,None]*n
    total_other = external+normal_vectors.sum(axis=0)
    mobility = P/mass
    coupling = np.einsum('aij,jk,bkl->abil', R, mobility, R)
    speed_scale = max(np.linalg.norm(v0), float(np.linalg.norm(walls)) if nc else 0.,
                      dt*np.linalg.norm(mobility@total_other), np.finfo(float).tiny)
    zero_tolerance = 128*np.finfo(float).eps*speed_scale

    def slip(v):
        return np.einsum('nki,nk->ni',bases,v-walls)

    def instantaneous(v):
        u = slip(v); speed = np.linalg.norm(u,axis=1)
        fixed = np.zeros((nc,3))
        sliding = speed > zero_tolerance
        directions = np.einsum('nki,ni->nk',bases,u)
        fixed[sliding] = -limits[sliding,None]*directions[sliding]/speed[sliding,None]
        zero = ~sliding
        if np.any(zero):
            base = np.einsum('nki,k->ni',bases[zero],mobility@(total_other+fixed.sum(axis=0)))
            fixed[zero] = cone_solve(base,bases[zero],mobility,1.,limits[zero])
        return fixed

    def midpoint(v,h):
        base = slip(v)+.5*h*np.einsum('nki,k->ni',bases,mobility@total_other)
        f = cone_solve(base,bases,mobility,.5*h,limits)
        return v+h*mobility@(total_other+f.sum(axis=0)), f

    def interval(v,w,h,f,t,kind):
        midpoint_relative = .5*(v+w)
        midpoint_velocity = midpoint_relative+frame
        jn = h*normal_vectors; jt = h*f
        reaction = mass*(w-v)-h*external-jn.sum(axis=0)-jt.sum(axis=0)
        umid = slip(midpoint_relative)
        tangent_impulse = np.einsum('nki,nk->ni',bases,jt)
        heat = -np.einsum('ij,ij->i',tangent_impulse,umid)
        # Roundoff tolerance only for the rejection decision; heat remains raw.
        roundoff = 256*np.finfo(float).eps*np.linalg.norm(jt,axis=1)*(
            speed_scale+np.linalg.norm(v)+np.linalg.norm(w)+2*np.linalg.norm(walls,axis=1))
        if np.any(heat < -roundoff):
            failure=ContactConvergenceError('negative accepted-subinterval friction work')
            failure.interval_context=dict(start_s=t,dt_s=h,friction_heat_raw_J=heat,
                rounding_bound_J=roundoff,initial_relative_velocity=v0,
                before_relative_velocity=v,after_relative_velocity=w,
                tangential_midpoint=umid,tangential_impulse=jt,
                zero_tolerance=zero_tolerance,speed_scale=speed_scale)
            raise failure
        work_other = float(h*external@midpoint_velocity)
        work_normal = np.einsum('ij,j->i',jn,midpoint_velocity)
        wall_work = float(np.sum((jn+jt)*(walls+frame)))
        constraint_work = float(reaction@midpoint_velocity)
        dk = mass*float((w-v)@midpoint_velocity)
        residual = dk-work_other-constraint_work-float(np.sum((jn+jt)*midpoint_velocity))
        return dict(start_s=t, dt_s=h, kind=kind, v_before=v+frame, v_after=w+frame,
            normal_impulse=jn, tangential_impulse=jt, constraint_impulse=reaction,
            tangential_relative_midpoint=np.einsum('nki,ni->nk',bases,umid),
            tangent_coordinates_midpoint=umid, friction_heat_raw_J=heat,
            friction_rounding_bound_J=roundoff,relative_speed_scale_m_s=speed_scale,
            reference_frame_velocity=frame.copy(),
            other_work_J=work_other, normal_work_J=work_normal,
            constraint_work_J=constraint_work, wall_work_J=wall_work,
            delta_kinetic_J=dk, kinetic_work_residual_J=residual)

    v=v0.copy(); t=0.; records=[]; rejected=0; events=[]; attempts=0
    next_h=dt
    previous_speed=np.linalg.norm(slip(v),axis=1)
    while t < dt:
        attempts += 1
        if attempts > max_intervals or len(records) > max_intervals:
            raise ContactConvergenceError('contact subinterval budget exceeded')
        remaining=dt-t; h=min(next_h,remaining)
        f0=instantaneous(v); acceleration=mobility@(total_other+f0.sum(axis=0))
        u=slip(v); du=np.einsum('nki,k->ni',bases,acceleration)
        now_speed=np.linalg.norm(u,axis=1)
        reached=np.flatnonzero((previous_speed>zero_tolerance)&(now_speed<=zero_tolerance)).tolist()
        fresh=[c for c in reached if not any(c in e['contacts'] and
               abs(e['time_s']-t)<=16*np.finfo(float).eps*dt for e in events)]
        if fresh:
            events.append(dict(time_s=t,contacts=fresh,event='ZERO_SLIP_WITHIN_ROUNDOFF',
                residual_slip_speed=now_speed[fresh],zero_tolerance=zero_tolerance,
                following='STICK_FEASIBILITY; VELOCITY_NOT_OVERWRITTEN'))
        exact=True; stop=remaining+1.; stops=[]
        for c in range(nc):
            speed=np.linalg.norm(u[c])
            if speed <= zero_tolerance: continue
            transverse=du[c]-u[c]*float(u[c]@du[c])/(speed*speed)
            if np.linalg.norm(transverse) > 1e-10*max(np.linalg.norm(du[c]),speed/dt):
                exact=False
                if np.linalg.norm(du[c]):
                    h=min(h,.2*speed/np.linalg.norm(du[c]))
            elif float(u[c]@du[c]) < 0:
                tau=-speed*speed/float(u[c]@du[c])
                if tau < stop: stop=tau; stops=[c]
                elif abs(tau-stop) < 1e-12*dt: stops.append(c)
        if exact:
            h=min(remaining,stop)
            w=v+h*acceleration
            records.append(interval(v,w,h,f0,t,'EXACT_CONSTANT_DIRECTION'))
            if stop <= remaining:
                events.append(dict(time_s=t+h,contacts=stops,event='ZERO_SLIP',
                                   following='STICK_FEASIBILITY_AT_NEXT_INTERVAL'))
        else:
            # A scale-based cutoff can reject a representable stop interval:
            # the recorded counterexample has h=124 ulps of current t.
            # Check actual half/full time advancement instead; no tolerance or
            # velocity is enlarged/snapped and non-progress still rejects.
            if h <= 0 or t+h/2 == t or t+h == t:
                failure=ContactConvergenceError('contact interval underflow, no silent stick')
                failure.interval_context=dict(t_s=t,h_s=h,full_dt_s=dt,
                    time_ulp_s=float(np.spacing(t)),slip_speed=now_speed,
                    zero_tolerance=zero_tolerance,acceleration=acceleration)
                raise failure
            full,_=midpoint(v,h)
            half,f1=midpoint(v,h/2); w,f2=midpoint(half,h/2)
            tolerance=atol+rtol*max(np.linalg.norm(v),np.linalg.norm(w),speed_scale)
            error=np.linalg.norm(w-full)/3
            if error > tolerance:
                next_h=h/2; rejected+=1; continue
            records.append(interval(v,half,h/2,f1,t,'ADAPTIVE_MIDPOINT'))
            records.append(interval(half,w,h/2,f2,t+h/2,'ADAPTIVE_MIDPOINT'))
            next_h=min(remaining,2*h) if error < tolerance/8 else h
        previous_speed=now_speed
        t+=h; v=w
        if dt-t <= 16*np.finfo(float).eps*dt: t=dt
    return dict(velocity=v+frame, intervals=records, events=events, rejected_intervals=rejected,
        contact_mobility=coupling, constraint_projector=P,tangent_bases=bases,
        friction_heat_J=np.sum([a['friction_heat_raw_J'] for a in records],axis=0),
        tangential_impulse=np.sum([a['tangential_impulse'] for a in records],axis=0),
        normal_impulse=np.sum([a['normal_impulse'] for a in records],axis=0),
        constraint_impulse=np.sum([a['constraint_impulse'] for a in records],axis=0),
        other_work_J=sum(a['other_work_J'] for a in records),
        constraint_work_J=sum(a['constraint_work_J'] for a in records),
        wall_work_J=sum(a['wall_work_J'] for a in records),
        delta_kinetic_J=mass*float((v-v0)@(.5*(v+v0)+frame)),
        kinetic_work_residual_J=sum(a['kinetic_work_residual_J'] for a in records))


def domain_kick(wall, x, velocity, mass, temperature, time, dt, p,
                internal_force, actuator, particle_ids):
    """Small adapter: existing analytic faces, unchanged normal penalty/thermal BC."""
    from .geometry import surfaces, assert_outside_tool
    from .kick_diagnostics import FEATURES
    assert_outside_tool(x,time,p)
    grouped={}
    for feature,(indices,gaps,normals,wall_v,wall_t) in zip(FEATURES,surfaces(x,time,p)):
        for a,idx in enumerate(indices):
            overlap=p['dx']/2-gaps[a]
            if overlap <= 0: continue
            grouped.setdefault(int(idx),[]).append(dict(feature=feature,gap=float(gaps[a]),
                penetration=float(overlap),normal=normals[a],wall_velocity=wall_v[a],wall_temperature=wall_t))
    acceleration=internal_force/mass[:,None]
    after=velocity+dt*acceleration
    cap=np.asarray(actuator)==1
    after[cap,0]=p['traverse'];after[cap,2]=-p['feed']
    midpoint=.5*(velocity+after)
    reaction=mass[:,None]*(after-velocity)-dt*internal_force
    actuator_work=float(np.sum(reaction*midpoint))
    internal_work=float(np.sum(dt*internal_force*midpoint))
    heat=np.zeros(len(x))
    ledger=dict(wall_work_rate=0.,wall_heat_rate=0.,friction_power=0.,contact_heat_rate=0.,
        contact_potential=0.,max_penetration=0.,contacts=0,
        contact_subintervals=0,contact_rejected_subintervals=0,contact_zero_events=0,
        contact_min_raw_heat_J=None,contact_kinetic_work_residual_J=0.)
    traces=[]
    for idx,faces in grouped.items():
        normals=[a['normal'] for a in faces]
        normal_forces=wall.law.kn*np.array([a['penetration'] for a in faces])
        walls=[a['wall_velocity'] for a in faces]
        free=[0.,1.,0.] if cap[idx] else [1.,1.,1.]
        prescribed=[p['traverse'],velocity[idx,1],-p['feed']] if cap[idx] else velocity[idx]
        try:
            result=kick(velocity[idx],mass[idx],internal_force[idx],normals,
                normal_forces,walls,wall.law.mu,dt,free_axes=free,prescribed_velocity=prescribed)
        except (ValueError,ContactConvergenceError) as error:
            failure = ContactConvergenceError('particle_id='+str(int(particle_ids[idx]))+
                ' features='+','.join(a['feature'] for a in faces)+'; '+str(error))
            failure.contact_context = dict(particle_id=int(particle_ids[idx]), faces=faces,
                position=x[idx].copy(), velocity_before=velocity[idx].copy(),
                internal_force=internal_force[idx].copy(), mass_kg=float(mass[idx]),
                time_mid_s=time, dt_s=dt, free_axes=free, prescribed_velocity=prescribed,
                failed_cone=getattr(error, 'cone_context', None),
                completed_local_traces_in_rejected_candidate=traces)
            raise failure from error
        actuator_work-=float(reaction[idx]@midpoint[idx])
        actuator_work+=result['constraint_work_J']
        internal_work-=float(dt*internal_force[idx]@midpoint[idx])
        internal_work+=result['other_work_J']
        after[idx]=result['velocity']
        q=result['friction_heat_J']
        flux=sum(p['contact_conductance']*p['dx']**2*(a['wall_temperature']-temperature[idx]) for a in faces)
        heat[idx]=p['heat_partition']*float(q.sum())/dt+flux
        ledger['wall_work_rate']+=result['wall_work_J']/dt
        ledger['wall_heat_rate']+=(1-p['heat_partition'])*float(q.sum())/dt-flux
        ledger['friction_power']+=float(q.sum())/dt
        ledger['contact_heat_rate']+=float(flux)
        ledger['contact_potential']+=float(np.sum(normal_forces**2)/(2*wall.law.kn))
        ledger['max_penetration']=max(ledger['max_penetration'],max(a['penetration'] for a in faces))
        ledger['contacts']+=len(faces)
        ledger['contact_subintervals']+=len(result['intervals'])
        ledger['contact_rejected_subintervals']+=result['rejected_intervals']
        ledger['contact_zero_events']+=len(result['events'])
        minimum=min(float(a['friction_heat_raw_J'].min()) for a in result['intervals'])
        previous=ledger['contact_min_raw_heat_J']
        ledger['contact_min_raw_heat_J']=minimum if previous is None else min(previous,minimum)
        ledger['contact_kinetic_work_residual_J']+=result['kinetic_work_residual_J']
        traces.append(dict(particle_id=int(particle_ids[idx]),position=x[idx].copy(),
            faces=faces,internal_force=internal_force[idx].copy(),mass_kg=float(mass[idx]),
            time_mid_s=time,dt_s=dt,local_kick=result))
    ledger['frozen_internal_impulse_work_J']=internal_work
    return (after-velocity)/dt,heat,ledger,actuator_work,traces



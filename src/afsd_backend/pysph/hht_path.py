"""Compatible HHT path, explicit force-filter work, and whole-path contact audit.

No forces are invented here and no state is committed. These are numerical
path certificates; they cannot replace a force-driven coupled endpoint solve.
"""
import numpy as np


def parameters(weight=.9):
    a=float(weight)
    if not np.isfinite(a) or not 2/3<=a<=1:
        raise ValueError('HHT force weight outside registered linear range')
    return 1.5-a, (2-a)**2/4


def endpoint(x0,v0,acc0,acc1,dt,weight=.9):
    h=float(dt);gamma,beta=parameters(weight)
    if not np.isfinite(h) or h<=0:raise ValueError('positive finite dt')
    arrays=[np.asarray(a,float) for a in (x0,v0,acc0,acc1)]
    if any(a.shape!=arrays[0].shape or not np.isfinite(a).all() for a in arrays):
        raise ValueError('matching finite state arrays')
    x,v,a0,a1=arrays
    return x+h*v+h*h*((.5-beta)*a0+beta*a1),v+h*((1-gamma)*a0+gamma*a1)


def sample(x0,v0,acc0,acc1,dt,fraction,weight=.9):
    s=float(fraction)
    if not np.isfinite(s) or not 0<=s<=1:raise ValueError('path fraction outside interval')
    x1,v1=endpoint(x0,v0,acc0,acc1,dt,weight)
    gamma,beta=parameters(weight);d=beta-gamma/2
    v0=np.asarray(v0);da=np.asarray(acc1)-acc0;dv=v1-v0
    v=v0+s*dv+6*dt*d*da*s*(1-s)
    x=np.asarray(x0)+dt*(s*v0+.5*s*s*dv+dt*d*da*(3*s*s-2*s**3))
    return x,v


def linear_work(force0,force1,v0,v1,acc0,acc1,dt,weight=.9):
    f0,f1,v0,v1,a0,a1=map(lambda v:np.asarray(v,float),(force0,force1,v0,v1,acc0,acc1))
    if f0.ndim<1 or f0.shape[-1]!=3 or any(v.shape!=f0.shape or not np.isfinite(v).all() for v in (f0,f1,v0,v1,a0,a1)):
        raise ValueError('matching finite force/path arrays')
    if not np.isfinite(dt) or dt<=0:raise ValueError('positive finite dt')
    gamma,beta=parameters(weight);d=beta-gamma/2
    nodes=(.5-1/np.sqrt(12),.5+1/np.sqrt(12));rows=[]
    for s in nodes:
        f=f0+s*(f1-f0);v=v0+s*(v1-v0)+6*dt*d*(a1-a0)*s*(1-s)
        impulse=.5*dt*f
        rows.append(dict(fraction=s,impulse_Ns=impulse,velocity_m_s=v,
                         work_J=np.sum(impulse*v,axis=-1)))
    mean=.5*(v0+v1)+dt*d*(a1-a0)
    return dict(impulse_Ns=.5*dt*(f0+f1),work_J=sum(r['work_J'] for r in rows),
        quadrature=rows,mean_velocity_m_s=mean,interpolation='LINEAR_FORCE_QUADRATIC_VELOCITY')


def kinetic_identity(mass,v0,v1,acc0,acc1,force0,force1,dt,weight=.9):
    m=np.asarray(mass,float);v0=np.asarray(v0,float);v1=np.asarray(v1,float)
    a0=np.asarray(acc0,float);a1=np.asarray(acc1,float);f0=np.asarray(force0,float);f1=np.asarray(force1,float)
    if m.shape!=v0.shape[:-1] or not np.isfinite(m).all() or np.any(m<=0):raise ValueError('positive finite mass per vector')
    gamma,beta=parameters(weight);d=beta-gamma/2;eta=1-weight;q=eta*(.5-eta)
    r0=m[...,None]*a0-f0;eq=m[...,None]*a1-(weight*f1+(1-weight)*f0)
    alg=dt*((1-gamma)*r0+q*(f1-f0))
    residual_impulse=dt*gamma*eq
    work=linear_work(f0,f1,v0,v1,a0,a1,dt,weight)
    mean=.5*(v0+v1);df=f1-f0
    filter_work=np.sum(alg*mean,axis=-1)-dt*dt*d*np.sum(.5*(f0+f1)*(a1-a0),axis=-1)-dt/12*np.sum(df*(v1-v0),axis=-1)
    root_work=np.sum(residual_impulse*mean,axis=-1)
    deltaK=.5*m*np.sum((v1-v0)*(v1+v0),axis=-1)
    return dict(physical_impulse_Ns=work['impulse_Ns'],algorithmic_filter_impulse_Ns=alg,
        root_residual_impulse_Ns=residual_impulse,
        momentum_identity_residual_Ns=m[...,None]*(v1-v0)-work['impulse_Ns']-alg-residual_impulse,
        physical_linear_force_work_J=work['work_J'],algorithmic_filter_work_J=filter_work,
        root_residual_work_J=root_work,delta_kinetic_J=deltaK,
        kinetic_identity_residual_J=deltaK-work['work_J']-filter_work-root_work,
        algorithmic_work_is_not_heat=True,algorithmic_work_sign_not_assumed=True)


def _roots(coeff):
    c=np.trim_zeros(np.asarray(coeff,float),'b')
    if len(c)<=1:return []
    values=np.roots(c[::-1]);result=[]
    for value in values:
        if abs(value.imag)<=1e-9 and 0<float(value.real)<1:result.append(float(value.real))
    return sorted(result)


def contact_audit(normal,force0,force1,wall0,wall1,v0,v1,acc0,acc1,dt,
                  capacity0,capacity1,weight=.9):
    """Fixed-normal force/slip polynomial audit, NEVER an endpoint solver.

Checks power everywhere, not just endpoints or quadrature. Exact interior
Coulomb direction is NOT implied for turning slip. Zero-slip needs caller events.
"""
    if not np.isfinite(dt) or dt<=0:raise ValueError('positive finite dt')
    n,f0,f1,w0,w1,v0,v1,a0,a1=map(lambda x:np.asarray(x,float),
        (normal,force0,force1,wall0,wall1,v0,v1,acc0,acc1))
    if any(a.shape!=(3,) or not np.isfinite(a).all() for a in (n,f0,f1,w0,w1,v0,v1,a0,a1)):
        raise ValueError('finite Cartesian contact vectors required')
    if abs(np.linalg.norm(n)-1)>1e-12:raise ValueError('unit contact normal')
    caps=np.asarray([capacity0,capacity1],float)
    if not np.isfinite(caps).all() or np.any(caps<0):raise ValueError('nonnegative capacity')
    eps=np.finfo(float).eps;fs=max(np.linalg.norm(f0),np.linalg.norm(f1),float(caps.max()),np.finfo(float).tiny)
    for f,cap in zip((f0,f1),caps):
        if abs(f@n)>256*eps*fs or np.linalg.norm(f)>cap+256*eps*fs:
            raise ValueError('traction outside tangent cone')
    gamma,beta=parameters(weight);d=beta-gamma/2;da=a1-a0
    def tangent(v):return v-(v@n)*n
    u=np.array([tangent(v0-w0),tangent(v1-v0+6*dt*d*da-(w1-w0)),tangent(-6*dt*d*da)])
    df=f1-f0
    power=np.array([f0@u[0],f0@u[1]+df@u[0],f0@u[2]+df@u[1],df@u[2]])
    candidates=[0.,1.]+_roots([power[1],2*power[2],3*power[3]])
    vals=np.polynomial.polynomial.polyval(candidates,power)
    squared=np.array([u[0]@u[0],2*u[0]@u[1],u[1]@u[1]+2*u[0]@u[2],2*u[1]@u[2],u[2]@u[2]])
    minima=[0.,1.]+_roots([squared[1],2*squared[2],3*squared[3],4*squared[4]])
    # A tangential touch makes the quartic derivative have a multiple root.
    # Its generic cubic roots are ill-conditioned. Component quadratics and
    # their vertices provide independent candidates; EVERY candidate is still
    # checked against the full vector slip, never snapped or accepted by name.
    for axis in range(3):
        minima+=_roots(u[:,axis])
        if u[2,axis]!=0:
            vertex=-u[1,axis]/(2*u[2,axis])
            if 0<vertex<1:minima.append(float(vertex))
    speed_scale=max(float(np.linalg.norm(u,axis=1).sum()),np.finfo(float).tiny)
    zero_bound=256*eps*speed_scale
    speeds=[float(np.linalg.norm(u[0]+s*u[1]+s*s*u[2])) for s in minima]
    zero_events=[]
    for s,speed in sorted(zip(minima,speeds)):
        if 0<s<1 and speed<=zero_bound and (not zero_events or abs(s-zero_events[-1])>1e-8):
            zero_events.append(s)
    stick=bool(np.linalg.norm(u,axis=1).sum()<=zero_bound)
    work=linear_work(f0,f1,v0-w0,v1-w1,a0,a1,dt,weight)
    q=-float(work['work_J']);bound=256*eps*fs*speed_scale
    for row in work['quadrature']:
        row['relative_velocity_m_s']=row.pop('velocity_m_s')
        row['friction_heat_raw_J']=-float(row.pop('work_J'))
    return dict(whole_path_passive=bool(float(max(vals))<=bound),
        requires_event_split=bool(zero_events),interior_zero_fractions=zero_events,
        all_interval_stick=stick,minimum_checked_slip_m_s=min(speeds),
        zero_speed_bound_m_s=zero_bound,maximum_friction_power_W=float(max(vals)),
        power_rounding_bound_W=bound,friction_heat_raw_J=q,
        friction_rounding_bound_J=dt*bound,quadrature=work['quadrature'],
        power_stationary_fractions=candidates,slip_coefficients_m_s=u,power_coefficients_W=power,
        interior_coulomb_direction='NOT_CERTIFIED_BY_PASSIVITY',accepted=False)

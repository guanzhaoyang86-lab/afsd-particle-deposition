"""Independent RK4 located-stop/reverse extension, never a tested solver oracle.

Only the old independent all-pair spatial formula and literal inputs are
imported. Frozen parent reference states/ledgers are verified and preserved.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
from hht_slow_inputs import inputs,REFERENCE_STEPS_S,SNAPSHOT_STEP_S
from reference_hht_block import pair_fields

ROOT=Path(__file__).resolve().parents[1]
N=27; LEDGER_SIZE=72; EPS=np.finfo(float).eps
SLOTS=dict(wall_work=slice(0,2),normal_material_work=slice(2,4),
    tangent_material_work=slice(4,6),reservoir_heat=slice(6,8),
    friction_raw=slice(8,10),impulse=slice(10,16),stress_work=16,
    plastic_work=17,particle_friction=slice(18,45),particle_wall_heat=slice(45,72))
CONTRACT=ROOT/'docs/PROCESS_STOP_REFERENCE_DESIGN.md'
ADDENDUM=ROOT/'docs/PROCESS_STOP_SEQUENTIAL_ADDENDUM.md'
ALIAS_ADDENDUM=ROOT/'docs/PROCESS_STOP_ENDPOINT_ALIAS_ADDENDUM.md'


class ReferenceStop(RuntimeError):
    def __init__(self,kind,details=None,state=None):
        super().__init__(kind);self.kind=kind;self.details=details or {};self.state=state


def jsonable(obj):
    if isinstance(obj,np.ndarray):return obj.tolist()
    if isinstance(obj,np.generic):return obj.item()
    if isinstance(obj,dict):return {k:jsonable(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):return [jsonable(v) for v in obj]
    return obj


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unpack(y):
    return y[:81].reshape(27,3),y[81:162].reshape(27,3),y[162:405].reshape(27,3,3),y[405:432],y[432:459],y[459:]


def zero_bound(slip):return 512*EPS*max(1.,float(np.max(np.abs(slip))))


def static_disk(driving_tangent,capacity):
    force=np.asarray(driving_tangent,float);cap=np.asarray(capacity,float)
    if force.ndim!=2 or cap.shape!=force.shape[:-1] or not np.isfinite(force).all() or not np.isfinite(cap).all() or np.any(cap<0):
        raise ValueError('finite tangent vectors and nonnegative capacities')
    required=np.linalg.norm(force,axis=-1);feasible=required<=cap
    direction=np.zeros_like(force);moving=required>0
    direction[moving]=force[moving]/required[moving,None]
    traction=-cap[:,None]*direction;traction[feasible]=-force[feasible]
    return dict(feasible=feasible,required_norm_N=required,capacity_N=cap.copy(),
        force_margin_N=cap-required,right_limit_direction=direction,
        friction_force_N=traction,resultant_tangent_force_N=force+traction)


def rk4_step(fun,t,y,h):
    if not np.isfinite(h) or h<=0:raise ValueError('positive RK4 interval')
    k1=fun(t,y);k2=fun(t+h/2,y+h*k1/2);k3=fun(t+h/2,y+h*k2/2);k4=fun(t+h,y+h*k3)
    result=y+h*(k1+2*k2+2*k3+k4)/6
    coeff=np.stack((y,h*k1,h*(-1.5*k1+k2+k3-.5*k4),h*(2/3)*(k1-k2-k3+k4)))
    return result,coeff


def locate_zero(stepper,t0,y0,width,slips,directions,maximum_width=None,deadline=None,allow_near_incoming=False):
    """Locate mean scalar root, then certify EVERY full vector; never snap y."""
    before=np.array(y0,copy=True);directions=np.asarray(directions,float)
    u0=np.asarray(slips(y0),float)
    if (u0.shape!=directions.shape or u0.ndim!=2 or not np.isfinite(u0).all()
            or not np.allclose(np.linalg.norm(directions,axis=1),1.,atol=1e-14,rtol=0)):
        raise ValueError('finite matching full-slip vectors/unit incoming directions')
    if not np.isfinite(width) or width<=0:raise ValueError('positive event bracket duration')
    cap_width=width if maximum_width is None else float(maximum_width)
    duration_tolerance=max(5e-19,128*EPS*max(abs(t0),width))
    bound=zero_bound(u0);g0=float(np.mean(np.sum(u0*directions,axis=1)))
    if g0<=(0. if allow_near_incoming else bound):
        raise ReferenceStop('INITIAL_EVENT_NOT_STRICTLY_INCOMING',dict(g0=g0,bound=bound))
    trials=[]
    def evaluate(tau):
        if len(trials)>=32:raise ReferenceStop('EVENT_ROOT_BUDGET',dict(candidates=trials))
        if deadline is not None and time.perf_counter()>deadline:raise TimeoutError('event reference budget')
        state=stepper(tau);u=np.asarray(slips(state),float)
        if not np.isfinite(state).all() or u.shape!=u0.shape:raise ReferenceStop('INVALID_EVENT_CANDIDATE')
        value=float(np.mean(np.sum(u*directions,axis=1)))
        row=dict(status='REJECTED_PRIVATE_INCOMING_BRANCH_CANDIDATE',accepted=False,restart_allowed=False,
            time_s=float(t0+tau),duration_s=float(tau),mean_incoming_velocity_m_s=value,
            tangential_velocity_m_s=u.copy(),state=state.copy())
        trials.append(row);return tau,value,state,u
    low=(0.,g0,before,u0);high=evaluate(width)
    for _ in range(4):
        if high[1]<0:break
        width=min(cap_width,2*high[0])
        if width<=high[0]:break
        high=evaluate(width)
    if high[1]>=0:raise ReferenceStop('EVENT_NOT_BRACKETED',dict(candidates=trials))
    while len(trials)<32:
        proposed=low[0]-low[1]*(high[0]-low[0])/(high[1]-low[1])
        if not low[0]<proposed<high[0]:proposed=.5*(low[0]+high[0])
        point=evaluate(proposed)
        if abs(point[1])<=bound/16:
            if np.max(np.linalg.norm(point[3],axis=1))>bound:
                raise ReferenceStop('SCALAR_ROOT_NOT_FULL_VECTOR_ZERO',dict(candidates=trials,
                    tangential_velocity_m_s=point[3],bound_m_s=bound))
            radius=.49*duration_tolerance
            left=evaluate(max(0.,point[0]-radius));right=evaluate(point[0]+radius)
            if left[1]>=0 and right[1]<=0 and right[0]-left[0]<=duration_tolerance:
                if not np.array_equal(y0,before):raise RuntimeError('event locator mutated input')
                return dict(time_s=float(t0+point[0]),state=point[2].copy(),
                    interval_s=[float(t0+left[0]),float(t0+right[0])],
                    mean_residual_m_s=point[1],full_slip_m_s=point[3].copy(),
                    full_zero_bound_m_s=bound,time_width_bound_s=duration_tolerance,
                    evaluations=len(trials),candidates=trials,velocity_projection_used=False)
        if point[1]>0:low=point
        else:high=point
    raise ReferenceStop('EVENT_ROOT_BUDGET',dict(candidates=trials))


def wall_program(t,p):
    omega=2*np.pi*p['normal_frequency_Hz'];d=p['normal_stroke']*np.sin(omega*t/2)**2
    rate=.5*p['normal_stroke']*omega*np.sin(omega*t)
    return np.array([-p['dx']/2+p['normal_preload']+d,p['length']+p['dx']/2-p['normal_preload']-d]),np.array([[rate,0,.2],[-rate,0,.2]])


def fields(t,y,data,mode=None):
    x,v,S,rho,T,ell=unpack(y);p=data['parameters']
    if not np.isfinite(y).all() or np.any(rho<=0) or np.any(T<p['room']) or np.any(T>=p['melt']):
        raise ReferenceStop('INVALID_PHYSICAL_DOMAIN',dict(time_s=t),y.copy())
    f=pair_fields(x,v,S,rho,data);L=f['gradient'];D=.5*(L+L.transpose(0,2,1));W=.5*(L-L.transpose(0,2,1))
    tr=np.trace(D,axis1=1,axis2=2);dev=D-tr[:,None,None]*np.eye(3)/3
    rate=np.sqrt(2/3*np.sum(dev*dev,axis=(1,2)));theta=(T-p['room'])/(p['melt']-p['room'])
    thermal=-np.expm1(-np.logaddexp(0,-p['c1']*theta)/p['c2'])
    ys=p['a1']*(1+p['b1']*theta**p['b2']*np.log(np.maximum(rate,p['rate_floor'])/p['reference_rate']))*thermal
    margin=ys-np.sqrt(1.5*np.sum(S*S,axis=(1,2)))
    if margin.min()<=0:raise ReferenceStop('YIELD',dict(time_s=t,minimum_margin_Pa=float(margin.min())),y.copy())
    planes,w=wall_program(t,p);fn=np.zeros((2,N,3));ft=fn.copy();heat=np.zeros((2,N));Q=heat.copy();gap=heat.copy();limits=[]
    for face,ids in enumerate(data['face_nodes']):
        n=np.array([1. if face==0 else -1.,0,0]);overlap=(planes[face]-x[ids,0])*n[0]+p['dx']/2
        if overlap.min()<=0:raise ReferenceStop('OPENING',dict(time_s=t,face=face,particle_id=int(data['particle_id'][ids[np.argmin(overlap)]])),y.copy())
        gap[face,ids]=overlap
        for z,j in enumerate(ids):
            load=p['normal_stiffness']*overlap[z];cap=p['friction']*load
            u=v[j]-w[face];u-=np.dot(u,n)*n;speed=float(np.linalg.norm(u));bound=zero_bound(u)
            override=False
            incoming_override=bool(mode and j in mode['directions'] and
                (mode['kind']=='INCOMING_PRIVATE' or j in mode.get('incoming_indices',())))
            epoch=mode.get('event_epochs',{}).get(j,mode.get('time_s',float('-inf'))) if mode else float('-inf')
            if incoming_override:
                direction=mode['directions'][j];override=True
            elif mode and j in mode['directions'] and speed<=bound and abs(t-epoch)<=mode['time_bound_s']:
                direction=mode['directions'][j];override=True
                limits.append(dict(time_s=float(t),particle_id=int(data['particle_id'][j]),
                    slip_m_s=u.copy(),speed_bound_m_s=bound,kind='EVENT_RIGHT_LIMIT_WITHOUT_VELOCITY_PROJECTION'))
            elif (mode and mode.get('allow_small_outgoing',False) and j in mode.get('event_epochs',{})
                    and t>=epoch and speed>0 and float(u@mode['directions'][j])>0):
                direction=u/speed
            elif speed<=bound:
                raise ReferenceStop('NEXT_FULL_VECTOR_ZERO',dict(time_s=t,particle_id=int(data['particle_id'][j]),face=face,slip_m_s=u),y.copy())
            else:direction=u/speed
            fn[face,j]=load*n;ft[face,j]=-cap*direction;Q[face,j]=-float(ft[face,j]@u)
            if Q[face,j]<0 and not override:
                raise ReferenceStop('NEGATIVE_RAW_FRICTION_POWER',dict(time_s=t,particle_id=int(data['particle_id'][j]),power_W=float(Q[face,j])),y.copy())
            if limits and limits[-1]['particle_id']==int(data['particle_id'][j]):
                limits[-1]['raw_friction_power_W']=float(Q[face,j])
                limits[-1]['root_roundoff_power_bound_W']=float(cap*bound)
                if Q[face,j]<-cap*bound:raise ReferenceStop('RIGHT_LIMIT_POWER_OUTSIDE_ROOT_ROUNDOFF')
            heat[face,j]=p['contact_conductance']*p['dx']**2*(data['wall_temperature_K'][face]-T[j])
    force=fn+ft
    b=dict(force=force.sum(axis=0),face_force=force.sum(axis=1),
        wall_work=np.einsum('fna,fa->f',force,w),normal_material_work=np.einsum('fna,na->f',fn,v),
        tangent_material_work=np.einsum('fna,na->f',ft,v),wall_heat=heat,friction=Q,penetration=gap,
        normal_potential=.5*p['normal_stiffness']*float(np.sum(gap*gap)),right_limit_records=limits)
    return f,b,dict(L=L,dev=dev,W=W,trace=tr,yield_margin=margin)


def rhs(t,y,data,mode=None,limit_records=None):
    x,v,S,rho,T,ell=unpack(y);p=data['parameters'];m=data['mass_kg'];f,b,k=fields(t,y,data,mode)
    if limit_records is not None:limit_records.extend(b['right_limit_records'])
    pairheat=f['conductance']*(T[f['b']]-T[f['a']]);heat=np.zeros(N)
    np.add.at(heat,f['a'],pairheat);np.add.at(heat,f['b'],-pairheat)
    q=b['friction'].sum(axis=0);wall=b['wall_heat'].sum(axis=0)
    dT=(heat+wall+p['heat_partition']*q)/(m*p['cp']);dS=2*p['shear']*k['dev']+k['W']@S-S@k['W']
    dl=np.zeros(LEDGER_SIZE)
    for key in ('wall_work','normal_material_work','tangent_material_work'):dl[SLOTS[key]]=b[key]
    dl[SLOTS['reservoir_heat']]=b['wall_heat'].sum(axis=1)
    dl[SLOTS['friction_raw']]=b['friction'].sum(axis=1);dl[SLOTS['impulse']]=b['face_force'].ravel()
    dl[SLOTS['stress_work']]=np.sum((m/rho)[:,None,None]*f['stress']*k['L'])
    dl[SLOTS['particle_friction']]=q;dl[SLOTS['particle_wall_heat']]=wall
    return np.r_[v.ravel(),((f['force']+b['force'])/m[:,None]).ravel(),dS.ravel(),-rho*k['trace'],dT,dl]


def slip_rows(y,indices,data):
    # Walls have fixed x normal and constant tangent velocity .2e_z.
    u=unpack(y)[1][indices].copy();u[:,0]=0;u[:,2]-=data['parameters']['wall_speed'];return u


def event_group(y,data,parent_terminal):
    detail=parent_terminal['details'].get('dense_event',parent_terminal['details'])
    pid=int(detail['particle_id']);index=int(np.flatnonzero(data['particle_id']==pid)[0])
    x=data['reference_position_m'];mirror=x[index].copy();mirror[0]=data['parameters']['length']-mirror[0]
    partners=np.flatnonzero(np.max(np.abs(x-mirror),axis=1)<1e-14)
    if len(partners)!=1 or partners[0]==index:raise ReferenceStop('EVENT_MIRROR_NOT_UNIQUE')
    indices=np.array(sorted([index,int(partners[0])]),int);u=slip_rows(y,indices,data)
    if np.max(np.abs(u[:,1]))>zero_bound(u):raise ReferenceStop('EVENT_GROUP_TRANSVERSE_SYMMETRY_NOT_PRESENT')
    if np.any(u[:,2]>=0):raise ReferenceStop('PARENT_EVENT_GROUP_NOT_INCOMING')
    # Direction is the measured incoming limit in the mechanical symmetry
    # plane. No physical velocity component is set or projected.
    directions=np.tile([0.,0.,-1.],(len(indices),1))
    return indices,directions


def real_roots(coeff):
    coeff=np.asarray(coeff,float);nz=np.flatnonzero(coeff!=0)
    if len(nz)==0 or nz[-1]==0:return []
    return [float(z.real) for z in np.polynomial.polynomial.polyroots(coeff[:nz[-1]+1]) if abs(z.imag)<=1e-8 and 0<z.real<1]


def endpoint_alias_certificate(poly,h,start_time,event_time,direction,s,bound):
    """Certify a positive outgoing component-root alias, without changing data."""
    poly=np.asarray(poly,float);direction=np.asarray(direction,float)
    if poly.shape!=(4,3) or direction.shape!=(3,) or not 0<s<=1 or h<=0:return None
    value=poly[0]+s*(poly[1]+s*(poly[2]+s*poly[3]))
    if np.linalg.norm(value)>bound:return None
    elapsed=start_time+s*h-event_time
    if elapsed<0:return None
    signed=poly@direction;derivative=np.array([signed[1],2*signed[2],3*signed[3]])/h
    # Polynomial extrema, not samples, bound the whole one-sided interval.
    velocity_points=[0.,s]+[r for r in real_roots(derivative) if r<s]
    acceleration_points=[0.,s]
    if derivative[2]!=0:
        vertex=-derivative[1]/(2*derivative[2])
        if 0<vertex<s:acceleration_points.append(float(vertex))
    minimum_velocity=min(signed[0]+q*(signed[1]+q*(signed[2]+q*signed[3])) for q in velocity_points)
    final_velocity=signed[0]+s*(signed[1]+s*(signed[2]+s*signed[3]))
    minimum_acceleration=min(derivative[0]+q*(derivative[1]+q*derivative[2]) for q in acceleration_points)
    if minimum_velocity<0 or final_velocity<=0 or minimum_acceleration<=0:return None
    width=bound/minimum_acceleration
    if elapsed>width:return None
    return dict(kind='CERTIFIED_POSITIVE_OUTGOING_ENDPOINT_ROUNDOFF_ALIAS',
        full_slip_coefficients_m_s=poly.copy(),candidate_slip_m_s=value,
        event_time_s=event_time,elapsed_since_event_s=elapsed,original_full_zero_bound_m_s=bound,
        derivative_derived_roundoff_window_s=width,projected_velocity_min_m_s=float(minimum_velocity),
        projected_velocity_at_candidate_m_s=float(final_velocity),projected_acceleration_min_m_s2=float(minimum_acceleration),
        velocity_extrema_fractions=velocity_points,acceleration_extrema_fractions=acceleration_points)


def dense_check(t,h,coeff,data,mode):
    vc=coeff[:,81:162].reshape(4,N,3).copy();vc[0,:,2]-=.2;vc[:,:,0]=0.
    aliases=[]
    for face,ids in enumerate(data['face_nodes']):
        for j in ids:
            poly=vc[:,j];norm2=sum((np.convolve(poly[:,k],poly[:,k]) for k in (1,2)),start=np.zeros(7))
            candidates=[1.]+real_roots(np.arange(1,7)*norm2[1:])+real_roots(poly[:,1])+real_roots(poly[:,2])
            for s in candidates:
                u=poly[0]+s*(poly[1]+s*(poly[2]+s*poly[3]))
                if np.linalg.norm(u)<=zero_bound(poly):
                    epoch=mode.get('event_epochs',{}).get(j,mode.get('time_s',float('-inf')))
                    if j in mode['directions'] and abs(t+s*h-epoch)<=mode['time_bound_s']:
                        aliases.append(dict(time_s=t+s*h,particle_id=int(data['particle_id'][j]),kind='INITIAL_EVENT_ROUNDOFF_ALIAS'))
                    else:
                        certificate=endpoint_alias_certificate(poly,h,t,epoch,mode['directions'][j],s,zero_bound(poly)) if j in mode.get('event_epochs',{}) else None
                        if certificate:
                            certificate.update(time_s=t+s*h,particle_id=int(data['particle_id'][j]));aliases.append(certificate)
                        else:raise ReferenceStop('NEXT_DENSE_FULL_VECTOR_ZERO',dict(time_s=t+s*h,particle_id=int(data['particle_id'][j]),face=face),coeff[0]+s*(coeff[1]+s*(coeff[2]+s*coeff[3])))
    for s in (.25,.5,.75):fields(t+s*h,coeff[0]+s*(coeff[1]+s*(coeff[2]+s*coeff[3])),data,mode)
    return aliases


def readback(t,y,data):
    x,v,S,rho,T,ell=unpack(y);f,b,k=fields(t,y,data)
    p=data['parameters'];m=data['mass_kg']
    result=dict(time_s=float(t),face_force_N=b['face_force'],minimum_yield_margin_Pa=float(k['yield_margin'].min()),
        minimum_overlap_m=float(b['penetration'][b['penetration']>0].min()),condition_max=f['condition'],
        thermal_balance_residual_J=float(np.sum(m*p['cp']*(T-data['temperature_K']))-ell[6:8].sum()-p['heat_partition']*ell[8:10].sum()),
        kinetic_J=float(.5*np.sum(m[:,None]*v*v)),normal_potential_J=b['normal_potential'])
    for key,slot in SLOTS.items():result[key]=ell[slot].copy() if isinstance(slot,slice) else float(ell[slot])
    return result


def select_first_event(probes):
    """Measured times, not labels, determine order; overlapping brackets block."""
    ordered=sorted(probes,key=lambda item:item[1]['time_s'])
    if len(ordered)>1 and ordered[0][1]['interval_s'][1]>=ordered[1][1]['interval_s'][0]:
        raise ReferenceStop('INDIVIDUAL_EVENT_ORDER_NOT_CERTIFIED',dict(
            measured=[dict(identity=key,time_s=row['time_s'],interval_s=row['interval_s']) for key,row in ordered]))
    return ordered[0]


def extend(parent,name,data,step,deadline,end=8e-6):
    archive=np.load(parent/(name+'_history.npz'));report=json.loads((parent/'report.json').read_text())
    times=list(archive['time_s']);states=list(archive['states']);y=states[-1].copy();t=float(times[-1]);initial=y.copy()
    old_terminal=report[name+'_terminal'];indices,direction=event_group(y,data,old_terminal)
    def fun(mode,records=None):
        def f(at,state):
            if time.perf_counter()>deadline:raise TimeoutError('combined stop-reference55s budget')
            return rhs(at,state,data,mode,records)
        return f
    remaining=list(map(int,indices));directions={int(j):d.copy() for j,d in zip(indices,direction)}
    epochs={};events=[];all_candidates=[];limits=[];mode=None
    while remaining:
        # Only the not-yet-switched incoming branches are formally continued.
        incoming=dict(kind='SEQUENTIAL_PRIVATE',directions={j:d.copy() for j,d in directions.items()},
            incoming_indices=remaining.copy(),event_epochs=epochs.copy(),time_bound_s=5e-19,
            allow_small_outgoing=True)
        start_state=y.copy();start_time=t;width=old_terminal['interval_s'][1]-t
        probes=[]
        for j in remaining:
            located=locate_zero(lambda tau:rk4_step(fun(incoming),start_time,start_state,tau)[0],
                start_time,start_state,width,lambda z,j=j:slip_rows(z,[j],data),directions[j][None,:],
                maximum_width=step,deadline=deadline,allow_near_incoming=bool(events))
            for row in located['candidates']:
                row['root_particle_id']=int(data['particle_id'][j]);row['branch_event_count']=len(events)
                all_candidates.append(row)
            probes.append((j,located))
        j,located=select_first_event(probes)
        event_time=located['time_s'];event_state=located['state'].copy()
        f,b,k=fields(event_time,event_state,data,incoming)
        drive=f['force'][[j]].copy();drive[:,0]=0
        face=0 if j in data['face_nodes'][0] else 1
        capacity=np.array([data['parameters']['friction']*data['parameters']['normal_stiffness']*b['penetration'][face,j]])
        disk=static_disk(drive,capacity)
        if np.any(disk['feasible']):
            raise ReferenceStop('FEASIBLE_3D_STICK_NOT_IMPLEMENTED',dict(located=located,static_disk=disk),event_state)
        left_derivative=rhs(event_time,event_state,data,incoming)
        old_direction=directions[j].copy();directions[j]=disk['right_limit_direction'][0].copy()
        remaining.remove(j);epochs[j]=event_time
        mode=dict(kind='SEQUENTIAL_RIGHT',directions={key:d.copy() for key,d in directions.items()},
            incoming_indices=remaining.copy(),event_epochs=epochs.copy(),time_bound_s=located['time_width_bound_s'],
            allow_small_outgoing=True)
        right_derivative=rhs(event_time,event_state,data,mode,limits)
        event=dict(time_s=event_time,interval_s=located['interval_s'],particle_id=int(data['particle_id'][j]),
            particle_index=j,full_slip_m_s=located['full_slip_m_s'],full_zero_bound_m_s=located['full_zero_bound_m_s'],
            time_width_bound_s=located['time_width_bound_s'],mean_residual_m_s=located['mean_residual_m_s'],
            evaluations=located['evaluations'],static_disk=disk,
            incoming_direction=old_direction,right_limit_direction=disk['right_limit_direction'][0],
            acceleration_left_m_s2=left_derivative[81:162].reshape(N,3),
            acceleration_right_m_s2=right_derivative[81:162].reshape(N,3),
            physical_state_before=event_state.copy(),physical_state_after=event_state.copy(),
            physical_impulse_jump_Ns=0.,physical_work_jump_J=0.,physical_heat_jump_J=0.,velocity_projection_used=False)
        events.append(event);y=event_state;t=event_time;times.append(t);states.append(y.copy())
    shared=np.load(parent/'reference_history.npz');sample_t=list(shared['time_s'])
    sample_y=list(shared['coarse_states' if name=='coarse' else 'states'])
    sample_index=len(sample_t);next_sample=sample_index*SNAPSHOT_STEP_S;new_rejected=[];aliases=[];post_steps=0;terminal=None
    while t<end-1e-20:
        target=min(t+step,next_sample,end);h=target-t
        try:
            candidate,coeff=rk4_step(fun(mode,limits),t,y,h)
            fields(target,candidate,data,mode);aliases.extend(dense_check(t,h,coeff,data,mode))
        except ReferenceStop as exc:
            new_rejected.append(dict(status='REJECTED',accepted=False,restart_allowed=False,
                interval_s=[t,target],kind=exc.kind,details=exc.details,state=exc.state))
            terminal=dict(status='NEXT_UNSUPPORTED_BRANCH',kind=exc.kind,interval_s=[t,target],details=exc.details);break
        y=candidate;t=target;times.append(t);states.append(y.copy());post_steps+=1
        if abs(t-next_sample)<1e-20:
            sample_t.append(next_sample);sample_y.append(y.copy());sample_index+=1;next_sample=sample_index*SNAPSHOT_STEP_S
    return dict(time_s=np.asarray(times),states=np.asarray(states),sample_time_s=np.asarray(sample_t),sample_states=np.asarray(sample_y),
        events=events,event_candidates=all_candidates,right_limit_records=limits,start_aliases=aliases,
        rejected=new_rejected,terminal=terminal or dict(status='COMPLETE_8US',kind=None),
        post_event_steps=post_steps,last_accepted_time_s=t,parent_accepted_state_unchanged=bool(np.array_equal(initial,archive['states'][-1])))


def compare(coarse,fine,data):
    count=min(len(coarse['sample_time_s']),len(fine['sample_time_s']));ta=coarse['sample_time_s'][:count];tb=fine['sample_time_s'][:count]
    if not np.allclose(ta,tb,rtol=0,atol=1e-19):raise ValueError('different shared times')
    a=coarse['sample_states'][:count];b=fine['sample_states'][:count];delta=a-b;V=data['reference_volume_m3'];p=data['parameters']
    diff=dict(position_m=float(np.max(np.sqrt(np.sum(V[None,:,None]*delta[:,:81].reshape(-1,N,3)**2,axis=(1,2))/V.sum()))),
        velocity_m_s=float(np.max(np.sqrt(np.sum(V[None,:,None]*delta[:,81:162].reshape(-1,N,3)**2,axis=(1,2))/V.sum()))),
        temperature_K=float(np.max(abs(delta[:,432:459]))))
    records=[];stress=[]
    for t,ya,yb in zip(ta,a,b):
        records.append((readback(t,ya,data),readback(t,yb,data)))
        Sa=unpack(ya)[2]-p['bulk']*(unpack(ya)[3]/p['rho0']-1)[:,None,None]*np.eye(3)
        Sb=unpack(yb)[2]-p['bulk']*(unpack(yb)[3]/p['rho0']-1)[:,None,None]*np.eye(3)
        stress.append(float(np.sqrt(np.sum(V[:,None,None]*(Sa-Sb)**2)/V.sum())))
    diff['stress_Pa']=max(stress)
    for key in ('face_force_N','impulse','wall_work','friction_raw','reservoir_heat'):
        diff[key]=float(max(np.max(abs(np.asarray(ra[key])-rb[key])) for ra,rb in records))
    bounds=dict(position_m=1e-11,velocity_m_s=2e-6,stress_Pa=70.,temperature_K=1e-9,
        face_force_N=.007,impulse=1.4e-9,wall_work=2e-9,friction_raw=2e-9,reservoir_heat=2e-10)
    return dict(passed=all(diff[k]<=bounds[k] for k in bounds),differences=diff,bounds=bounds,
        common_sample_count=count,common_end_time_s=float(ta[-1]),new_post_event_samples=count-96)


def scalar_reference(other=-2.,mass=1.,capacity=1.,u0=1.,end=1.):
    """Independent force-driven constant-force scalar test; no velocity snap."""
    y=np.array([0.,u0,0.,0.]);incoming=np.sign(u0);ft=-capacity*incoming
    def scalar_rhs(force):
        return lambda t,z:np.array([z[1],(other+force)/mass,-force*z[1],other*z[1]])
    located=locate_zero(lambda h:rk4_step(scalar_rhs(ft),0.,y,h)[0],0.,y,end,
        lambda z:np.array([[z[1]]]),np.array([[incoming]]))
    disk=static_disk(np.array([[other]]),np.array([capacity]));right_force=float(disk['friction_force_N'][0,0])
    state=located['state'].copy();remaining=end-located['time_s']
    final=rk4_step(scalar_rhs(right_force),located['time_s'],state,remaining)[0]
    return dict(event_time_s=located['time_s'],event_state=state,final=final,static=disk,
        root=located,physical_jump=np.zeros_like(state),delta_kinetic_J=.5*mass*(final[1]**2-u0**2))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.perf_counter();data=inputs()
    source={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),CONTRACT,ADDENDUM,ALIAS_ADDENDUM,
        ROOT/'scripts/hht_slow_inputs.py',ROOT/'scripts/reference_hht_block.py')}
    report=dict(status='RUNNING',classification='INDEPENDENT_LOCATED_STOP_REVERSE_REFERENCE_NOT_AFSD',
        production_authorized=False,new_afsd_steps=0,source_sha256=source,
        configuration_sha256=hashlib.sha256(json.dumps(jsonable(data),sort_keys=True).encode()).hexdigest())
    try:
        parent=args.parent.resolve();old=json.loads((parent/'report.json').read_text())
        if old['status']!='PASS_PREFIX_REFERENCE_EVENT_BLOCKED':raise ValueError('verified parent prefix/event required')
        for name,digest in old['source_sha256'].items():
            if sha(ROOT/name)!=digest:raise ValueError('old source changed '+name)
        for name,digest in old['artifact_sha256'].items():
            if sha(parent/name)!=digest:raise ValueError('old artifact changed '+name)
        if report['configuration_sha256']!=old['configuration_sha256']:raise ValueError('physical inputs changed')
        report['parent_report_sha256']=sha(parent/'report.json');report['parent_artifact_sha256']=old['artifact_sha256']
        runs=[]
        for name,step in zip(('coarse','fine'),REFERENCE_STEPS_S):
            run=extend(parent,name,data,step,started+55);runs.append(run)
            np.savez_compressed(args.output/(name+'_history.npz'),time_s=run['time_s'],states=run['states'],
                sample_time_s=run['sample_time_s'],sample_states=run['sample_states'],particle_id=data['particle_id'])
            np.savez_compressed(args.output/(name+'_event_candidates.npz'),states=np.asarray([r['state'] for r in run['event_candidates']]))
            np.savez_compressed(args.output/(name+'_event_state.npz'),
                time_s=np.asarray([event['time_s'] for event in run['events']]),
                states_before=np.asarray([event['physical_state_before'] for event in run['events']]),
                states_after=np.asarray([event['physical_state_after'] for event in run['events']]),particle_id=data['particle_id'])
            for key in ('event_candidates','rejected'):
                metadata=[]
                for row in run[key]:metadata.append({k:v for k,v in row.items() if k!='state'})
                (args.output/(name+'_'+key+'.json')).write_text(json.dumps(jsonable(metadata),indent=2,allow_nan=False)+'\n',encoding='utf-8')
            if run['rejected']:
                np.savez_compressed(args.output/(name+'_rejected.npz'),states=np.asarray([r['state'] for r in run['rejected']]))
        cmp=compare(*runs,data)
        maps=[{event['particle_id']:event for event in run['events']} for run in runs]
        event_differences={str(pid):abs(maps[0][pid]['time_s']-maps[1][pid]['time_s']) for pid in maps[0] if pid in maps[1]}
        agreement=bool(maps[0].keys()==maps[1].keys() and len(event_differences)==2 and
            all(value<=1e-12 for value in event_differences.values()) and
            all(not np.any(event['static_disk']['feasible']) for mapping in maps for event in mapping.values()))
        full=all(r['terminal']['status']=='COMPLETE_8US' for r in runs)
        status=('PASS_EVENT_REFERENCE_8US' if full else 'PASS_EVENT_REFERENCE_PREFIX_NEXT_BRANCH_BLOCKED') if cmp['passed'] and agreement and cmp['new_post_event_samples']>0 else 'FAIL_EVENT_REFERENCE_COMPARISON'
        count=cmp['common_sample_count']
        np.savez_compressed(args.output/'reference_history.npz',time_s=runs[1]['sample_time_s'][:count],states=runs[1]['sample_states'][:count],
            coarse_states=runs[0]['sample_states'][:count],particle_id=data['particle_id'],mass_kg=data['mass_kg'],reference_position_m=data['reference_position_m'])
        public=[]
        for run in runs:
            item={k:v for k,v in run.items() if k not in ('time_s','states','sample_time_s','sample_states','event_candidates','rejected')}
            item['events']=[{k:v for k,v in event.items() if k not in ('physical_state_before','physical_state_after')} for event in item['events']]
            public.append(item)
        report.update(status=status,prefix_comparison=cmp,certified_common_end_time_s=cmp['common_end_time_s'] if cmp['passed'] and agreement else 0.,
            event_agreement_passed=agreement,event_time_difference_s=event_differences,event_time_difference_bound_s=1e-12,
            reference_complete_8us=full,event_status='TWO_SEQUENTIAL_FULL_ZEROS_LOCATED_AND_REVERSE_SLIDING',reference_runs=public,
            desired_1ms_status='NOT_RUN',general_multiface_actuator_status='NOT_IMPLEMENTED')
        if source!={name:sha(ROOT/name) for name in source}:raise RuntimeError('source changed during reference')
        for name,digest in old['artifact_sha256'].items():
            if sha(parent/name)!=digest:raise RuntimeError('parent artifact changed during reference')
    except Exception as exc:
        report.update(status='FAIL_EVENT_REFERENCE',error=type(exc).__name__+': '+str(exc),diagnostic=getattr(exc,'details',None),
            prefix_comparison=dict(passed=False),certified_common_end_time_s=0.,event_status='UNVERIFIED')
    report['wall_seconds']=time.perf_counter()-started
    report['artifact_sha256']={p.name:sha(p) for p in args.output.iterdir() if p.is_file()}
    (args.output/'report.json').write_text(json.dumps(jsonable(report),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(jsonable({k:v for k,v in report.items() if k not in ('reference_runs','diagnostic')}),indent=2))
    return 0 if report['status'].startswith('PASS') else 1


if __name__=='__main__':raise SystemExit(main())

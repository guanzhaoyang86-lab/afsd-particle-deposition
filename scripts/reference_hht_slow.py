"""Independent current-geometry RK4 slow-loading startup/event reference.

No tested backend/path/material/thermal routines are imported. The old oracle's
independent all-pair algebra is reused, not its fast wall programme or guards.
An event-approach state is NOT a solution of the branch after that event.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np

from hht_slow_inputs import inputs, CLASSIFICATION, REFERENCE_STEPS_S, SNAPSHOT_STEP_S
from reference_hht_block import (N, ENERGY_SLOTS, LEDGER_SIZE, pair_fields,
    unpack, packed_initial, jsonable)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT/'docs/PROCESS_SLOW_CONTRACT.md'
EPS = np.finfo(float).eps
SOURCES = (Path(__file__).resolve(), ROOT/'scripts/hht_slow_inputs.py',
    ROOT/'scripts/reference_hht_block.py', CONTRACT)


class DomainBoundary(RuntimeError):
    def __init__(self, kind, t, state, **details):
        super().__init__(kind)
        self.kind=kind; self.t=float(t); self.state=state.copy(); self.details=details


def wall_program(t, p):
    omega=2*np.pi*p['normal_frequency_Hz']
    # sin squared avoids cancellation for the deliberately slow normal stroke.
    stroke=p['normal_stroke']*np.sin(.5*omega*t)**2
    rate=.5*p['normal_stroke']*omega*np.sin(omega*t)
    plane=np.array([-p['dx']/2+p['normal_preload']+stroke,
        p['length']+p['dx']/2-p['normal_preload']-stroke])
    return plane,np.array([[rate,0.,p['wall_speed']],[-rate,0.,p['wall_speed']]])


def boundary_fields(t,x,v,T,data,check=True):
    p=data['parameters']; planes,w=wall_program(t,p)
    normal=np.array([[1.,0.,0.],[-1.,0.,0.]])
    fn=np.zeros((2,N,3)); ft=fn.copy(); heat=np.zeros((2,N)); q=heat.copy()
    gaps=np.zeros((2,N)); tangent=np.zeros((2,N,3)); minimum=np.inf
    details=[]
    for face,ids in enumerate(data['face_nodes']):
        gap=(planes[face]-x[ids,0])*normal[face,0]+p['dx']/2
        rel=v[ids]-w[face]; ut=rel-(rel@normal[face])[:,None]*normal[face]
        speed=np.linalg.norm(ut,axis=1)
        for j,pid in enumerate(ids):
            details.append(dict(face=face,particle_id=int(data['particle_id'][pid]),
                particle_index=int(pid),overlap_m=float(gap[j]),
                tangential_relative_velocity_m_s=ut[j].copy(),speed_m_s=float(speed[j])))
        if check and np.any(gap<=0):
            j=int(np.argmin(gap)); raise ValueError(('OPENING',face,int(ids[j]),float(gap[j])))
        if check and np.any(speed<=0):
            j=int(np.argmin(speed)); raise ValueError(('ZERO_SLIP_SINGULAR',face,int(ids[j]),float(speed[j])))
        if np.any(speed==0):
            # No division regularization, static branch or invented direction.
            raise ValueError(('ZERO_SLIP_SINGULAR',face,int(ids[np.argmin(speed)]),0.))
        gaps[face,ids]=gap; tangent[face,ids]=ut
        load=p['normal_stiffness']*gap
        fn[face,ids]=load[:,None]*normal[face]
        ft[face,ids]=-p['friction']*load[:,None]*ut/speed[:,None]
        q[face,ids]=-np.sum(ft[face,ids]*ut,axis=1)
        heat[face,ids]=p['contact_conductance']*p['dx']**2*(data['wall_temperature_K'][face]-T[ids])
        minimum=min(minimum,float(speed.min()))
    force=fn+ft
    return dict(force=force.sum(axis=0),face_force=force.sum(axis=1),
        wall_work=np.einsum('fna,fa->f',force,w),
        normal_material_work=np.einsum('fna,na->f',fn,v),
        tangent_material_work=np.einsum('fna,na->f',ft,v),
        friction=q,wall_heat=heat,penetration=gaps,tangent=tangent,
        minimum_speed_m_s=minimum,contacts=details,
        normal_potential=.5*p['normal_stiffness']*float(np.sum(gaps*gaps)))


def fields(t,y,data):
    x,v,S,rho,T,ell=unpack(y); p=data['parameters']
    if not np.all(np.isfinite(y)) or np.any(rho<=0) or np.any(T<p['room']) or np.any(T>=p['melt']):
        raise DomainBoundary('INVALID_PHYSICAL_DOMAIN',t,y)
    try: f=pair_fields(x,v,S,rho,data)
    except ValueError as exc: raise DomainBoundary('INVALID_GEOMETRY',t,y,error=str(exc)) from exc
    L=f['gradient']; D=.5*(L+L.transpose(0,2,1)); W=.5*(L-L.transpose(0,2,1))
    trace=np.trace(D,axis1=1,axis2=2); dev=D-trace[:,None,None]*np.eye(3)/3
    rate=np.sqrt(2/3*np.sum(dev*dev,axis=(1,2))); theta=(T-p['room'])/(p['melt']-p['room'])
    thermal=-np.expm1(-np.logaddexp(0,-p['c1']*theta)/p['c2'])
    ys=p['a1']*(1+p['b1']*theta**p['b2']*np.log(np.maximum(rate,p['rate_floor'])/p['reference_rate']))*thermal
    margin=ys-np.sqrt(1.5*np.sum(S*S,axis=(1,2)))
    if np.any(margin<=0):
        j=int(np.argmin(margin)); raise DomainBoundary('YIELD',t,y,particle_index=j,particle_id=int(data['particle_id'][j]),yield_margin_Pa=float(margin[j]))
    try: b=boundary_fields(t,x,v,T,data)
    except ValueError as exc:
        kind,face,j,value=exc.args[0]
        raise DomainBoundary(kind,t,y,face=face,particle_index=j,particle_id=int(data['particle_id'][j]),margin=value) from exc
    return f,b,dict(L=L,D=D,W=W,trace=trace,dev=dev,yield_margin=margin)


def rhs(t,y,data):
    x,v,S,rho,T,ell=unpack(y); p=data['parameters']; m=data['mass_kg']
    f,b,k=fields(t,y,data)
    pairheat=f['conductance']*(T[f['b']]-T[f['a']]); conduction=np.zeros(N)
    np.add.at(conduction,f['a'],pairheat); np.add.at(conduction,f['b'],-pairheat)
    q=b['friction'].sum(axis=0); wall=b['wall_heat'].sum(axis=0)
    dT=(conduction+wall+p['heat_partition']*q)/(m*p['cp'])
    dS=2*p['shear']*k['dev']+k['W']@S-S@k['W']
    dl=np.zeros(LEDGER_SIZE)
    for key in ('wall_work','normal_material_work','tangent_material_work'): dl[ENERGY_SLOTS[key]]=b[key]
    dl[ENERGY_SLOTS['reservoir_heat']]=b['wall_heat'].sum(axis=1)
    dl[ENERGY_SLOTS['friction_raw']]=b['friction'].sum(axis=1)
    dl[ENERGY_SLOTS['impulse']]=b['face_force'].ravel()
    dl[ENERGY_SLOTS['stress_work']]=np.sum((m/rho)[:,None,None]*f['stress']*k['L'])
    dl[ENERGY_SLOTS['particle_friction']]=q; dl[ENERGY_SLOTS['particle_wall_heat']]=wall
    return np.r_[v.ravel(),((f['force']+b['force'])/m[:,None]).ravel(),
        dS.ravel(),-rho*k['trace'],dT,dl]


def dense_coefficients(y,h,k):
    """Classical RK4 cubic extension, coefficients in ascending powers."""
    k1,k2,k3,k4=k
    return np.stack((y,h*k1,h*(-1.5*k1+k2+k3-.5*k4),h*(2/3)*(k1-k2-k3+k4)))


def dense_value(coeff,s):
    return coeff[0]+s*(coeff[1]+s*(coeff[2]+s*coeff[3]))


def real_roots(coeff):
    coeff=np.asarray(coeff,float)
    nz=np.flatnonzero(coeff!=0)
    if len(nz)==0 or nz[-1]==0:return []
    roots=np.polynomial.polynomial.polyroots(coeff[:nz[-1]+1])
    return sorted(float(z.real) for z in roots if abs(z.imag)<=1e-8 and 0<z.real<1)


def tangent_events(coeff,initial_positions,data):
    """No z-only stop: full vector checked at roots and full-norm extrema."""
    vc=coeff[:,3*N:6*N].reshape(4,N,3).copy()
    vc[0,:,2]-=data['parameters']['wall_speed']; vc[:,:,0]=0.
    events=[]; turns=[]; minimum=np.inf
    for face,ids in enumerate(data['face_nodes']):
        for j in ids:
            poly=vc[:,j]; norm2=np.zeros(7)
            for axis in (1,2): norm2+=np.convolve(poly[:,axis],poly[:,axis])
            derivative=np.arange(1,7)*norm2[1:]
            zroots=real_roots(poly[:,2]); yroots=real_roots(poly[:,1])
            sample=sorted(set([0.,1.]+zroots+yroots+real_roots(derivative)))
            scale=max(1.,float(np.max(np.abs(poly))),data['parameters']['wall_speed'])
            tol=512*EPS*scale
            for s in sample:
                u=poly[0]+s*(poly[1]+s*(poly[2]+s*poly[3])); speed=float(np.linalg.norm(u))
                minimum=min(minimum,speed)
                if 0<s<1 and speed<=tol:
                    events.append(dict(kind='ZERO_SLIP_DENSE_ROOT',face=face,particle_index=int(j),
                        particle_id=int(data['particle_id'][j]),fraction=s,speed_m_s=speed,
                        tangential_relative_velocity_m_s=u,roundoff_bound_m_s=tol,
                        symmetry_centerline=bool(initial_positions[j,1]==data['parameters']['length']/2)))
            for s in zroots:
                u=poly[0]+s*(poly[1]+s*(poly[2]+s*poly[3])); speed=float(np.linalg.norm(u))
                turns.append(dict(face=face,particle_id=int(data['particle_id'][j]),fraction=s,
                    full_speed_m_s=speed,transverse_velocity_m_s=float(u[1]),
                    stop_certified=speed<=tol))
    return dict(events=events,component_turns=turns,minimum_dense_speed_m_s=minimum)


def dense_audit(t,h,coeff,data):
    p=data['parameters']; result=tangent_events(coeff,data['reference_position_m'],data)
    omega=2*np.pi*p['normal_frequency_Hz']; amp=.5*p['normal_stroke']
    # Cubic Taylor wall plus an explicit fourth derivative remainder bound.
    wall=np.array([2*amp*np.sin(.5*omega*t)**2,
        amp*omega*np.sin(omega*t)*h,.5*amp*omega**2*np.cos(omega*t)*h*h,
        -amp*omega**3*np.sin(omega*t)*h**3/6])
    remainder=amp*omega**4*h**4/24
    xc=coeff[:,:3*N].reshape(4,N,3)
    for face,ids in enumerate(data['face_nodes']):
        nx=1 if face==0 else -1
        plane0=-p['dx']/2+p['normal_preload'] if face==0 else p['length']+p['dx']/2-p['normal_preload']
        for j in ids:
            gap=-nx*xc[:,j,0]+wall
            gap[0]+=nx*plane0+p['dx']/2
            points=[0.,1.]+real_roots(np.arange(1,4)*gap[1:])
            values=[float(np.polynomial.polynomial.polyval(s,gap)) for s in points]
            if min(values)<=remainder:
                roots=real_roots(gap)
                s=roots[0] if roots else points[int(np.argmin(values))]
                result['events'].append(dict(kind='OPENING_DENSE_ROOT',face=face,particle_index=int(j),
                    particle_id=int(data['particle_id'][j]),fraction=s,minimum_overlap_m=min(values),
                    wall_taylor_remainder_bound_m=remainder))
    # Dense checks cover geometric/thermal/yield domain at fixed additional
    # nodes. They do not claim a whole-polynomial nonlinear FKS certificate.
    for s in (.25,.5,.75): fields(t+s*h,dense_value(coeff,s),data)
    return result


def rk4_candidate(t,y,h,data,deadline):
    k=[]
    for frac,z in ((0.,y),):
        if time.perf_counter()>deadline: raise TimeoutError('combined independent reference 55 s budget')
        k.append(rhs(t,z,data))
    for frac,z in ((.5,y+.5*h*k[0]),):
        k.append(rhs(t+frac*h,z,data))
    k.append(rhs(t+.5*h,y+.5*h*k[1],data)); k.append(rhs(t+h,y+h*k[2],data))
    candidate=y+h/6*(k[0]+2*k[1]+2*k[2]+k[3])
    fields(t+h,candidate,data)
    coeff=dense_coefficients(y,h,k)
    audit=dense_audit(t,h,coeff,data)
    if audit['events']:
        event=min(audit['events'],key=lambda row:row['fraction'])
        raise DomainBoundary(event['kind'],t+h,candidate,dense_event=event,
            candidate_interval_s=[t,t+h],dense_coefficient=coeff)
    return candidate,audit


def record(t,y,data):
    x,v,S,rho,T,ell=unpack(y); f,b,k=fields(t,y,data); p=data['parameters']; m=data['mass_kg']
    result=dict(time_s=float(t),face_force_N=b['face_force'],condition_max=f['condition'],
        minimum_speed_m_s=b['minimum_speed_m_s'],minimum_overlap_m=min(row['overlap_m'] for row in b['contacts']),
        minimum_yield_margin_Pa=float(k['yield_margin'].min()),contacts=b['contacts'],
        kinetic_J=float(.5*np.sum(m[:,None]*v*v)),normal_potential_J=b['normal_potential'],
        temperature_change_max_K=float(np.max(np.abs(T-data['temperature_K']))),
        thermal_balance_residual_J=float(np.sum(m*p['cp']*(T-data['temperature_K']))-
            ell[ENERGY_SLOTS['reservoir_heat']].sum()-p['heat_partition']*ell[ENERGY_SLOTS['friction_raw']].sum()),
        internal_net_force_N=f['force'].sum(axis=0),internal_net_torque_Nm=np.cross(x,f['force']).sum(axis=0),
        spatial_power_residual_W=float(np.sum(f['force']*v)+np.sum((m/rho)[:,None,None]*f['stress']*k['L'])))
    for key,slot in ENERGY_SLOTS.items(): result[key]=ell[slot].copy() if isinstance(slot,slice) else float(ell[slot])
    return result


def run_reference(data,step_s,deadline=None):
    if step_s not in REFERENCE_STEPS_S: raise ValueError('unregistered reference dt')
    started=time.perf_counter(); deadline=started+55 if deadline is None else deadline
    end=data['parameters']['end_time']; total=int(round(end/step_s)); y=packed_initial(data)
    states=[y.copy()]; times=[0.]; sample_states=[y.copy()]; sample_times=[0.]
    stride=int(round(SNAPSHOT_STEP_S/step_s)); turns=[]; rejected=[]; t=0.; terminal=None
    for index in range(total):
        h=step_s; nominal_end=(index+1)*step_s; event_mode=False
        actual_halvings=0; rejected_right=None
        # At most 12 actual attempted halvings. Accepted prefixes approach an
        # event but are not allowed to turn into a post-event integration path.
        for subdivision in range(13):
            try: candidate,audit=rk4_candidate(t,y,h,data,deadline)
            except DomainBoundary as exc:
                event_mode=True
                rejected.append(dict(status='REJECTED',accepted=False,restart_allowed=False,
                    start_time_s=t,candidate_end_s=t+h,failed_evaluation_time_s=exc.t,
                    kind=exc.kind,details=exc.details,state=exc.state.copy(),
                    candidate_attempt_index=subdivision,actual_halvings_before=actual_halvings))
                rejected_right=t+h if rejected_right is None else min(rejected_right,t+h)
                if subdivision==12:
                    terminal=dict(status='EVENT_APPROACH_BLOCKED',kind=exc.kind,
                        interval_s=[t,rejected_right],details={k:v for k,v in exc.details.items() if k!='dense_coefficient'},
                        candidate_attempts=subdivision+1,actual_halvings=actual_halvings)
                    break
                h*=.5; actual_halvings+=1
                continue
            if event_mode and t+h>=rejected_right:
                # A formerly rejected interval was traversed by positive
                # candidates. It may have been RK-stage overshoot, not an
                # event. Do not fabricate an inverted or nonexistent bracket.
                rejected.append(dict(status='REJECTED',accepted=False,restart_allowed=False,
                    start_time_s=t,candidate_end_s=t+h,failed_evaluation_time_s=t+h,
                    kind='EVENT_LOCALIZATION_INCONSISTENT',
                    details=dict(previous_rejected_right_s=rejected_right),state=candidate.copy(),
                    candidate_attempt_index=subdivision,actual_halvings_before=actual_halvings))
                terminal=dict(status='EVENT_LOCALIZATION_INCONSISTENT',kind='EVENT_LOCALIZATION_INCONSISTENT',
                    interval_s=None,details=dict(previous_rejected_right_s=rejected_right,
                        last_accepted_s=t,positive_candidate_end_s=t+h),
                    candidate_attempts=subdivision+1,actual_halvings=actual_halvings)
                break
            for turn in audit['component_turns']:
                turns.append(dict(turn,time_s=t+h*turn['fraction']))
            y=candidate; t=t+h; states.append(y.copy()); times.append(t)
            if event_mode:
                # This positive state narrows the right bracket retained by a
                # previous rejected candidate; resume with the SAME half width.
                if subdivision==12:
                    terminal=dict(status='EVENT_APPROACH_BLOCKED',kind=rejected[-1]['kind'],
                        interval_s=[t,rejected_right],
                        details={k:v for k,v in rejected[-1]['details'].items() if k!='dense_coefficient'},
                        candidate_attempts=subdivision+1,actual_halvings=actual_halvings)
                    break
                continue
            if (index+1)%stride==0: sample_times.append(nominal_end); sample_states.append(y.copy())
            break
        if event_mode:
            if terminal is None:
                terminal=dict(status='EVENT_APPROACH_BLOCKED',kind=rejected[-1]['kind'],
                    interval_s=[t,rejected_right],
                    details={k:v for k,v in rejected[-1]['details'].items() if k!='dense_coefficient'},
                    candidate_attempts=13,actual_halvings=actual_halvings)
            break
    terminal=terminal or dict(status='COMPLETE_REGISTERED_SCOUT',kind=None,interval_s=None)
    return dict(times_s=np.asarray(times),states=np.asarray(states),
        sample_times_s=np.asarray(sample_times),sample_states=np.asarray(sample_states),
        rejected=rejected,component_turns=turns,terminal=terminal,
        last_accepted_record=record(t,y,data),step_s=step_s,
        accepted_steps=len(times)-1,wall_seconds=time.perf_counter()-started,
        method='CLASSICAL_RK4_CURRENT_GEOMETRY_INDEPENDENT',restart_after_terminal_allowed=False)


def comparison(coarse,fine,data):
    count=min(len(coarse['sample_times_s']),len(fine['sample_times_s']))
    times=coarse['sample_times_s'][:count]
    np.testing.assert_allclose(times,fine['sample_times_s'][:count],rtol=0,atol=1e-20)
    a=coarse['sample_states'][:count]; b=fine['sample_states'][:count]; delta=a-b
    V=data['reference_volume_m3']; rows=[]; stress=[]; p=data['parameters']
    for ya,yb,t in zip(a,b,times):
        rows.append((record(t,ya,data),record(t,yb,data)))
        sa=unpack(ya)[2]-p['bulk']*(unpack(ya)[3]/p['rho0']-1)[:,None,None]*np.eye(3)
        sb=unpack(yb)[2]-p['bulk']*(unpack(yb)[3]/p['rho0']-1)[:,None,None]*np.eye(3)
        stress.append(float(np.sqrt(np.sum(V[:,None,None]*(sa-sb)**2)/V.sum())))
    diff=dict(position_m=float(np.max(np.sqrt(np.sum(V[None,:,None]*delta[:,:3*N].reshape(-1,N,3)**2,axis=(1,2))/V.sum()))),
        velocity_m_s=float(np.max(np.sqrt(np.sum(V[None,:,None]*delta[:,3*N:6*N].reshape(-1,N,3)**2,axis=(1,2))/V.sum()))),
        stress_Pa=max(stress),temperature_K=float(np.max(np.abs(delta[:,16*N:17*N]))))
    for key in ('face_force_N','impulse','wall_work','friction_raw','reservoir_heat'):
        diff[key]=float(max(np.max(np.abs(np.asarray(ra[key])-rb[key])) for ra,rb in rows))
    bound=dict(position_m=1e-11,velocity_m_s=2e-6,stress_Pa=70.,temperature_K=1e-9,
        face_force_N=.007,impulse=1.4e-9,wall_work=2e-9,friction_raw=2e-9,reservoir_heat=2e-10)
    return dict(common_sample_count=count,common_end_s=float(times[-1]),differences=diff,bounds=bound,
        passed=count>1 and all(diff[key]<=bound[key] for key in bound))


def event_agreement(coarse,fine):
    a=coarse['terminal']; b=fine['terminal']
    if any(e['status']=='EVENT_LOCALIZATION_INCONSISTENT' for e in (a,b)):
        return dict(passed=False,status='BLOCKED_INCONSISTENT_EVENT_LOCALIZATION')
    if a['kind'] is None or b['kind'] is None:
        return dict(passed=a['kind'] is None and b['kind'] is None,status='NO_EVENT_BOTH_COMPLETE' if a['kind'] is None and b['kind'] is None else 'BLOCKED_EVENT_MISMATCH')
    def identity(event):
        details=event['details'].get('dense_event',event['details'])
        return ('OPENING' if event['kind'].startswith('OPENING') else 'ZERO_SLIP' if event['kind'].startswith('ZERO_SLIP') else event['kind'],details.get('face'),details.get('particle_id'))
    ia,ib=identity(a),identity(b); aa=a['interval_s']; bb=b['interval_s']
    allowance=max(coarse['step_s'],aa[1]-aa[0],bb[1]-bb[0]); separation=max(0.,aa[0]-bb[1],bb[0]-aa[1])
    return dict(passed=ia==ib and separation<=allowance,status='APPROACH_ONLY_NO_POST_EVENT_BRANCH',
        identities=[ia,ib],intervals_s=[aa,bb],separation_s=separation,allowance_s=allowance)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter(); data=inputs(); source={str(p.relative_to(ROOT)):sha(p) for p in SOURCES}
    report=dict(status='RUNNING',classification=CLASSIFICATION,production_authorized=False,
        source_sha256=source,configuration_sha256=hashlib.sha256(json.dumps(jsonable(data),sort_keys=True).encode()).hexdigest(),
        desired_1ms_status='NOT_RUN_COST_EXCEEDS_REFERENCE_CAP',positive_plastic_accuracy='NOT_RUN_ELASTIC_ORACLE')
    try:
        runs=[]
        for name,h in zip(('coarse','fine'),REFERENCE_STEPS_S):
            result=run_reference(data,h,started+55.); runs.append(result)
            np.savez_compressed(args.output/(name+'_history.npz'),time_s=result['times_s'],states=result['states'],
                sample_time_s=result['sample_times_s'],sample_states=result['sample_states'],
                particle_id=data['particle_id'],mass_kg=data['mass_kg'],reference_position_m=data['reference_position_m'])
            rejected=result['rejected']
            np.savez_compressed(args.output/(name+'_rejected.npz'),states=np.asarray([r['state'] for r in rejected]))
            metadata=[]
            for r in rejected:
                item={k:v for k,v in r.items() if k!='state'}
                item['details']={k:v for k,v in item['details'].items() if k!='dense_coefficient'}
                metadata.append(item)
            (args.output/(name+'_rejected.json')).write_text(json.dumps(jsonable(metadata),indent=2,allow_nan=False)+'\n',encoding='utf-8')
        cmp=comparison(*runs,data); agreement=event_agreement(*runs)
        complete=all(r['terminal']['kind'] is None for r in runs)
        status=('PASS_COMPLETE_REFERENCE' if complete else 'PASS_PREFIX_REFERENCE_EVENT_BLOCKED') if cmp['passed'] and agreement['passed'] else 'FAIL_REFERENCE'
        common=cmp['common_sample_count']
        np.savez_compressed(args.output/'reference_history.npz',time_s=runs[1]['sample_times_s'][:common],
            states=runs[1]['sample_states'][:common],coarse_states=runs[0]['sample_states'][:common],
            particle_id=data['particle_id'],mass_kg=data['mass_kg'],reference_position_m=data['reference_position_m'])
        report.update(status=status,comparison=cmp,prefix_comparison=cmp,event_agreement=agreement,
            reference_runs=[{k:v for k,v in r.items() if k not in ('states','times_s','sample_states','sample_times_s','rejected')} for r in runs],
            reference_complete_20us=complete,certified_common_prefix_end_s=cmp['common_end_s'],
            certified_common_end_time_s=cmp['common_end_s'] if cmp['passed'] and agreement['passed'] else 0.,
            event_status='NO_EVENT' if complete else 'EVENT_APPROACH_BLOCKED',
            coarse_terminal=runs[0]['terminal'],fine_terminal=runs[1]['terminal'],
            full_slip_event_scope='DENSE_RK4_ROUNDOFF_ROOT_OR_APPROACH_ONLY; NO_POST_EVENT_SOLUTION',
            nonlinear_dense_domain_scope='RK_STAGES_ENDPOINT_AND_THREE_DENSE_SAMPLES_NOT_WHOLE_PATH_FKS_PROOF')
        if source!={str(p.relative_to(ROOT)):sha(p) for p in SOURCES}: raise RuntimeError('source changed during reference')
    except Exception as exc: report.update(status='FAIL_REFERENCE',prefix_comparison=dict(passed=False),
        certified_common_end_time_s=0.,event_status='UNVERIFIED',error=type(exc).__name__+': '+str(exc))
    report['wall_seconds']=time.perf_counter()-started
    report['artifact_sha256']={p.name:sha(p) for p in args.output.iterdir() if p.is_file()}
    (args.output/'report.json').write_text(json.dumps(jsonable(report),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(jsonable({k:v for k,v in report.items() if k!='reference_runs'}),indent=2))
    return 0 if report['status'] in ('PASS_COMPLETE_REFERENCE','PASS_PREFIX_REFERENCE_EVENT_BLOCKED') else 1


if __name__=='__main__': raise SystemExit(main())

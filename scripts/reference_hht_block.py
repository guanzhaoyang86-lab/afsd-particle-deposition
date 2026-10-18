"""Independent all-pair elastic/contact/thermal LOCAL_BENCHMARK reference.

Never imports the tested backend, HHT, material return, force or heat solver.
This is a bounded reference ODE, not a physical AFSD driver.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CLASSIFICATION = 'LOCAL_BENCHMARK_HHT_CURRENT_3D_TWO_FACE_THERMAL'
CONTRACT = ROOT/'docs/PROCESS_HHT_FIXTURE_CONTRACT.md'
N = 27
PHYSICAL_SIZE = 17*N
LEDGER_SIZE = 18+2*N
ENERGY_SLOTS = dict(wall_work=slice(0,2), normal_material_work=slice(2,4),
    tangent_material_work=slice(4,6), reservoir_heat=slice(6,8),
    friction_raw=slice(8,10), impulse=slice(10,16), stress_work=16,
    plastic_work=17, particle_friction=slice(18,18+N),
    particle_wall_heat=slice(18+N,18+2*N))


def fixture():
    """Literal inputs only; the tested driver may consume these arrays."""
    p = dict(length=.01, dx=.005, h_ratio=1.3, rho0=2700., young=70e9,
        shear=35e9, poisson=0., bulk=70e9/3, end_time=2e-6,
        normal_preload=2e-6, normal_stroke=2e-7, normal_stiffness=7e6,
        friction=.2, wall_speed=.2, conductivity=167., cp=900.,
        contact_conductance=1000., heat_partition=.5, plastic_heat_fraction=1.,
        a1=317.6e6, a2=1e6, a3=70., b1=12.44, b2=5.65,
        c1=13., c2=.01, room=293.15, melt=855.15,
        reference_rate=1., rate_floor=1.)
    ijk = np.indices((3,3,3)).reshape(3,-1).T
    x = p['dx']*ijk.astype(float)
    V = p['dx']**3*np.prod(np.where((ijk==0)|(ijk==2),.5,1.),axis=1)
    v = np.zeros((N,3)); v[:,2] = .02
    return dict(classification=CLASSIFICATION, production_authorized=False,
        parameters=p, particle_id=np.arange(1,N+1), reference_position_m=x.copy(),
        position_m=x, velocity_m_s=v, mass_kg=p['rho0']*V,
        reference_volume_m3=V, density_kg_m3=np.full(N,p['rho0']),
        smoothing_length_m=np.full(N,p['h_ratio']*p['dx']),
        deviatoric_stress_Pa=np.tile(np.diag([2e4,-1e4,-1e4]),(N,1,1)),
        ep=np.zeros(N), temperature_K=400+2*x[:,0]/p['length'],
        face_nodes=(np.flatnonzero(ijk[:,0]==0),np.flatnonzero(ijk[:,0]==2)),
        wall_temperature_K=np.array([410.,420.]))


def wall_program(t, p, side=None):
    phase=2*np.pi*t/p['end_time']
    stroke=.5*p['normal_stroke']*(1-np.cos(phase))
    speed=p['normal_stroke']*np.pi/p['end_time']*np.sin(phase)
    planes=np.array([-p['dx']/2+p['normal_preload']+stroke,
        p['length']+p['dx']/2-p['normal_preload']-stroke])
    right=(t>=p['end_time']/2) if side is None else side
    w=np.array([[speed,0.,p['wall_speed']*(-1 if right else 1)],
                [-speed,0.,p['wall_speed']*(-1 if right else 1)]])
    return planes,w


def fixture_data():
    """Shared physical data only; no shared force/boundary implementation."""
    data=fixture(); p=dict(data['parameters'])
    p.update(wall_temperature_K=[410.,420.], contact_radius=.0025,
        minimum_slip_m_s=.05)
    return dict(X=data['reference_position_m'].copy(),mass=data['mass_kg'].copy(),
        h=data['smoothing_length_m'].copy(),S0=data['deviatoric_stress_Pa'].copy(),
        rho0=data['density_kg_m3'].copy(),T0=data['temperature_K'].copy(),
        v0=data['velocity_m_s'].copy(),parameters=p,
        particle_id=data['particle_id'].copy(),face_nodes=data['face_nodes'])


def pair_fields(x,v,S,rho,data):
    """Independently enumerated scalar polynomial and current moments."""
    p=data['parameters']; m=data['mass_kg']; h=data['smoothing_length_m']
    a,b=np.triu_indices(N,1); rvec=x[a]-x[b]; r=np.linalg.norm(rvec,axis=1)
    if np.any(r<=0): raise ValueError('coincident independent-reference pair')
    hp=.5*(h[a]+h[b]); q=r/hp; active=q<2
    a,b,rvec,r,hp,q=(z[active] for z in (a,b,rvec,r,hp,q))
    radial=15/(62*np.pi*hp**4)*np.where(q<1,3*q*q-6,-3*(2-q)**2)
    g=radial[:,None]*rvec/r[:,None]; vol=m/rho
    moment=np.zeros((N,3,3)); term=-rvec[:,:,None]*g[:,None,:]
    np.add.at(moment,a,vol[b,None,None]*term)
    np.add.at(moment,b,vol[a,None,None]*term)
    eig=np.linalg.eigvalsh(moment)
    if np.any(eig[:,0]<=0) or np.any(eig[:,-1]>1e6*eig[:,0]):
        raise ValueError('invalid current reference correction moment')
    B=np.linalg.inv(moment)
    ga=np.einsum('nji,nj->ni',B[a],g)
    gb=np.einsum('nji,nj->ni',B[b],g)
    L=np.zeros((N,3,3)); dv=v[b]-v[a]
    np.add.at(L,a,vol[b,None,None]*dv[:,:,None]*ga[:,None,:])
    np.add.at(L,b,vol[a,None,None]*dv[:,:,None]*gb[:,None,:])
    pressure=p['bulk']*(rho/p['rho0']-1)
    sigma=S-pressure[:,None,None]*np.eye(3)
    pair=vol[a,None]*vol[b,None]*(np.einsum('nij,nj->ni',sigma[a],ga)+
        np.einsum('nij,nj->ni',sigma[b],gb))
    force=np.zeros((N,3)); np.add.at(force,a,pair); np.add.at(force,b,-pair)
    conductance=vol[a]*vol[b]*2*p['conductivity']*(-radial/r)
    return dict(a=a,b=b,g=g,B=B,gradient=L,stress=sigma,force=force,
        conductance=conductance,condition=float(np.max(eig[:,-1]/eig[:,0])))


def boundary_fields(t,x,v,T,data,side):
    p=data['parameters']; planes,w=wall_program(t,p,side)
    normals=np.array([[1.,0.,0.],[-1.,0.,0.]])
    fn=np.zeros((2,N,3)); ft=fn.copy(); heat=np.zeros((2,N)); q=heat.copy()
    penetration=np.zeros((2,N)); slip_margin=np.inf
    for face,ids in enumerate(data['face_nodes']):
        d=(planes[face]-x[ids,0])*normals[face,0]+p['dx']/2
        if np.any(d<=0): raise ValueError('reference contact opened; new event branch required')
        penetration[face,ids]=d; normal=p['normal_stiffness']*d
        rel=v[ids]-w[face]; tangent=rel-(rel@normals[face])[:,None]*normals[face]
        speed=np.linalg.norm(tangent,axis=1)
        # Both registered halves must maintain their signed z sliding branch.
        expected=1 if side else -1
        margin=expected*tangent[:,2]
        if np.any(margin<=.05) or np.any(speed<=.05):
            raise ValueError('reference zero/reversed slip inside smooth half; event branch required')
        slip_margin=min(slip_margin,float(margin.min()))
        fn[face,ids]=normal[:,None]*normals[face]
        ft[face,ids]=-p['friction']*normal[:,None]*tangent/speed[:,None]
        q[face,ids]=p['friction']*normal*speed
        heat[face,ids]=p['contact_conductance']*p['dx']**2*(data['wall_temperature_K'][face]-T[ids])
    force=fn+ft
    return dict(force=force.sum(axis=0),face_force=force.sum(axis=1),
        wall_work=np.einsum('fna,fa->f',force,w),
        normal_material_work=np.einsum('fna,na->f',fn,v),
        tangent_material_work=np.einsum('fna,na->f',ft,v),
        friction=q,wall_heat=heat,penetration=penetration,slip_margin=slip_margin,
        normal_potential=.5*p['normal_stiffness']*float(np.sum(penetration**2)))


def unpack(y):
    return (y[:3*N].reshape(N,3),y[3*N:6*N].reshape(N,3),
        y[6*N:15*N].reshape(N,3,3),y[15*N:16*N],y[16*N:17*N],y[17*N:])


def packed_initial(data):
    return np.r_[data['position_m'].ravel(),data['velocity_m_s'].ravel(),
        data['deviatoric_stress_Pa'].ravel(),data['density_kg_m3'],
        data['temperature_K'],np.zeros(LEDGER_SIZE)]


def state_scaling():
    scale=np.r_[np.full(3*N,.01),np.full(3*N,.2),np.full(9*N,1e6),
        np.full(N,2700.),np.full(N,10.),np.full(LEDGER_SIZE,1e-4)]
    offset=np.zeros_like(scale); offset[16*N:17*N]=400.
    return offset,scale


def rhs(t,y,data,side):
    x,v,S,rho,T,ledger=unpack(y); p=data['parameters']; m=data['mass_kg']
    if np.any(rho<=0) or np.any(T<p['room']) or np.any(T>=p['melt']):
        raise ValueError('reference density/temperature physical-domain violation')
    f=pair_fields(x,v,S,rho,data); L=f['gradient']; D=.5*(L+L.transpose(0,2,1))
    W=.5*(L-L.transpose(0,2,1)); trace=np.trace(D,axis1=1,axis2=2)
    dev=D-trace[:,None,None]*np.eye(3)/3
    rate=np.sqrt(2/3*np.sum(dev*dev,axis=(1,2)))
    theta=(T-p['room'])/(p['melt']-p['room'])
    thermal=-np.expm1(-np.logaddexp(0,-p['c1']*theta)/p['c2'])
    yield_stress=p['a1']*(1+p['b1']*theta**p['b2']*
        np.log(np.maximum(rate,p['rate_floor'])/p['reference_rate']))*thermal
    equivalent=np.sqrt(1.5*np.sum(S*S,axis=(1,2)))
    if np.any(equivalent>=yield_stress): raise ValueError('reference yielded; elastic oracle invalid')
    b=boundary_fields(t,x,v,T,data,side)
    pairheat=f['conductance']*(T[f['b']]-T[f['a']]); conduction=np.zeros(N)
    np.add.at(conduction,f['a'],pairheat); np.add.at(conduction,f['b'],-pairheat)
    q=b['friction'].sum(axis=0); wall=b['wall_heat'].sum(axis=0)
    dT=(conduction+wall+p['heat_partition']*q)/(m*p['cp'])
    dS=2*p['shear']*dev+W@S-S@W
    dl=np.zeros(LEDGER_SIZE)
    for key in ('wall_work','normal_material_work','tangent_material_work'):
        dl[ENERGY_SLOTS[key]]=b[key]
    dl[ENERGY_SLOTS['reservoir_heat']]=b['wall_heat'].sum(axis=1)
    dl[ENERGY_SLOTS['friction_raw']]=b['friction'].sum(axis=1)
    dl[ENERGY_SLOTS['impulse']]=b['face_force'].ravel()
    dl[ENERGY_SLOTS['stress_work']]=np.sum((m/rho)[:,None,None]*f['stress']*L)
    dl[ENERGY_SLOTS['particle_friction']]=q
    dl[ENERGY_SLOTS['particle_wall_heat']]=wall
    return np.r_[v.ravel(),((f['force']+b['force'])/m[:,None]).ravel(),
        dS.ravel(),-rho*trace,dT,dl]


def run_reference(data,steps_per_half=128,deadline=None):
    """Classical RK4, independently refined; no tested time-integration map."""
    if steps_per_half not in (128,256): raise ValueError('unregistered reference grid')
    y=packed_initial(data); start=time.perf_counter()
    if deadline is None: deadline=start+55.
    times=np.linspace(0,data['parameters']['end_time'],33)
    states=[y.copy()]; calls=0
    for side in (False,True):
        begin=data['parameters']['end_time']*.5*int(side)
        step=data['parameters']['end_time']*.5/steps_per_half
        def fun(t,z):
            nonlocal calls
            if time.perf_counter()>deadline: raise TimeoutError('independent reference 55 s total budget')
            calls+=1
            return rhs(t,z,data,side)
        for n in range(steps_per_half):
            t=begin+n*step
            k1=fun(t,y); k2=fun(t+.5*step,y+.5*step*k1)
            k3=fun(t+.5*step,y+.5*step*k2); k4=fun(t+step,y+step*k3)
            y=y+step/6*(k1+2*k2+2*k3+k4)
            if (n+1)%(steps_per_half//16)==0: states.append(y.copy())
    return dict(times_s=times,states=np.asarray(states),nfev=calls,
        wall_seconds=time.perf_counter()-start,method='CLASSICAL_RK4',
        steps_per_half=steps_per_half,step_s=step)


def readback_record(t,y,data,side=None):
    x,v,S,rho,T,ell=unpack(y)
    if side is None: side=t>=data['parameters']['end_time']/2
    f=pair_fields(x,v,S,rho,data); b=boundary_fields(t,x,v,T,data,side)
    p=data['parameters']; m=data['mass_kg']; L=f['gradient']; D=.5*(L+L.transpose(0,2,1))
    dev=D-np.trace(D,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
    rate=np.sqrt(2/3*np.sum(dev*dev,axis=(1,2))); th=(T-p['room'])/(p['melt']-p['room'])
    yield_stress=p['a1']*(1+p['b1']*th**p['b2']*np.log(np.maximum(rate,p['rate_floor']))) * (-np.expm1(-np.logaddexp(0,-p['c1']*th)/p['c2']))
    measured=float(np.sum(m*p['cp']*(T-data['temperature_K'])))
    supplied=float(ell[ENERGY_SLOTS['reservoir_heat']].sum()+p['heat_partition']*ell[ENERGY_SLOTS['friction_raw']].sum())
    r=dict(time_s=float(t),side='RIGHT' if side else 'LEFT',
        face_force_N=b['face_force'],condition_max=f['condition'],
        minimum_slip_margin_m_s=b['slip_margin'],
        minimum_overlap_m=float(b['penetration'][b['penetration']>0].min()),
        minimum_yield_margin_Pa=float(np.min(yield_stress-np.sqrt(1.5*np.sum(S*S,axis=(1,2))))),
        kinetic_J=float(.5*np.sum(m[:,None]*v*v)),
        normal_potential_J=b['normal_potential'],temperature_change_max_K=float(np.max(np.abs(T-data['temperature_K']))),
        thermal_increment_J=measured,thermal_balance_residual_J=measured-supplied,
        internal_net_force_N=f['force'].sum(axis=0),
        internal_net_torque_Nm=np.cross(x,f['force']).sum(axis=0),
        spatial_power_residual_W=float(np.sum(f['force']*v)+np.sum((m/rho)[:,None,None]*f['stress']*L)))
    for key,slot in ENERGY_SLOTS.items(): r[key]=ell[slot].copy() if isinstance(slot,slice) else float(ell[slot])
    return r


def jsonable(obj):
    if isinstance(obj,np.ndarray): return obj.tolist()
    if isinstance(obj,np.generic): return obj.item()
    if isinstance(obj,dict): return {k:jsonable(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)): return [jsonable(v) for v in obj]
    return obj


def convergence(coarse,fine,data):
    delta=coarse['states']-fine['states']; vol=data['reference_volume_m3']; rows=[]
    for a,b,t in zip(coarse['states'],fine['states'],fine['times_s']):
        ra=readback_record(t,a,data); rb=readback_record(t,b,data)
        rows.append((ra,rb))
    difference=dict(position_m=float(np.max(np.sqrt(np.sum(vol[None,:,None]*delta[:,:3*N].reshape(-1,N,3)**2,axis=(1,2))/vol.sum()))),
        velocity_m_s=float(np.max(np.sqrt(np.sum(vol[None,:,None]*delta[:,3*N:6*N].reshape(-1,N,3)**2,axis=(1,2))/vol.sum()))),
        temperature_K=float(np.max(np.abs(delta[:,16*N:17*N]))))
    stress=[]
    for a,b in zip(coarse['states'],fine['states']):
        xa,va,Sa,ra,_,_=unpack(a); xb,vb,Sb,rb,_,_=unpack(b)
        sa=Sa-data['parameters']['bulk']*(ra/data['parameters']['rho0']-1)[:,None,None]*np.eye(3)
        sb=Sb-data['parameters']['bulk']*(rb/data['parameters']['rho0']-1)[:,None,None]*np.eye(3)
        stress.append(np.sqrt(np.sum(vol[:,None,None]*(sa-sb)**2)/vol.sum()))
    difference['stress_Pa']=float(max(stress))
    for key in ('face_force_N','impulse','wall_work','friction_raw','reservoir_heat'):
        difference[key]=float(max(np.max(np.abs(np.asarray(a[key])-b[key])) for a,b in rows))
    bounds=dict(position_m=1e-11,velocity_m_s=2e-6,temperature_K=1e-9,stress_Pa=70.,
        face_force_N=.007,impulse=1.4e-9,wall_work=2e-9,friction_raw=2e-9,reservoir_heat=2e-10)
    return dict(differences=difference,bounds=bounds,passed=all(difference[k]<=bounds[k] for k in bounds))


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter(); data=fixture(); report=dict(classification=CLASSIFICATION,
        production_authorized=False,joint_hht_status='NOT_RUN',status='RUNNING')
    report['configuration_sha256']=hashlib.sha256(json.dumps(jsonable(data),sort_keys=True).encode()).hexdigest()
    report['source_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__).resolve(),CONTRACT)}
    try:
        coarse=run_reference(data,deadline=started+55.)
        fine=run_reference(data,steps_per_half=256,deadline=started+55.)
        cmp=convergence(coarse,fine,data); records=[readback_record(t,y,data) for t,y in zip(fine['times_s'],fine['states'])]
        switch_left=readback_record(fine['times_s'][16],fine['states'][16],data,False)
        np.savez_compressed(args.output/'reference_history.npz',time_s=fine['times_s'],
            states=fine['states'],coarse_states=coarse['states'],particle_id=data['particle_id'],
            reference_position_m=data['reference_position_m'],mass_kg=data['mass_kg'])
        report.update(status='PASS_INDEPENDENT_ELASTIC_REFERENCE' if cmp['passed'] else 'FAIL_REFERENCE_REFINEMENT',
            comparison=cmp,records=records,switch_left_record=switch_left,
            reference_runs=[{k:v for k,v in r.items() if k not in ('times_s','states')} for r in (coarse,fine)],
            particle_count=N,minimum_yield_margin_Pa=min(r['minimum_yield_margin_Pa'] for r in records),
            minimum_slip_margin_m_s=min(r['minimum_slip_margin_m_s'] for r in records),
            minimum_overlap_m=min(r['minimum_overlap_m'] for r in records),
            maximum_thermal_balance_residual_J=max(abs(r['thermal_balance_residual_J']) for r in records),
            material_energy_status='CURRENT_VOLUME_HYPOELASTIC_DIAGNOSTIC_NOT_CERTIFIED',
            positive_plastic_work_accuracy='NOT_TESTED_ELASTIC_PATH',
            reference_scope='SEMIDISCRETE_TEMPORAL_REFERENCE_NOT_CONTINUUM_VALIDATION')
    except Exception as exc:
        report.update(status='FAIL_REFERENCE',error=type(exc).__name__+': '+str(exc))
    report['wall_seconds']=time.perf_counter()-started
    (args.output/'report.json').write_text(json.dumps(jsonable(report),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(jsonable({k:v for k,v in report.items() if k not in ('records','switch_left_record')}),indent=2))
    return 0 if report['status'].startswith('PASS') else 1


if __name__=='__main__': raise SystemExit(main())

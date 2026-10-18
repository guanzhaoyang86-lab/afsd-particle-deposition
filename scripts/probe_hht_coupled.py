"""Bounded real 3D local HHT coupling; never advances the AFSD checkpoint."""
from pathlib import Path
import argparse
import copy
import gzip
import hashlib
import json
import os
import sys
import threading
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from reference_hht_block import fixture_data
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.mechanics import positions,velocities,pressure
from afsd_backend.pysph.material import equivalent,FKS
from afsd_backend.pysph.hht_coupled import initialize,advance,restart_at_control_event,restore_checkpoint
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.output import checkpoint,restore
from afsd_backend.pysph.application import memory_bytes
from afsd_backend.pysph import hht_path
from afsd_scaffold.vtk_legacy import write_point_cloud

CLASS='LOCAL_BENCHMARK_HHT_CURRENT_3D_TWO_FACE_THERMAL_NOT_AFSD'
SOURCES=['scripts/probe_hht_coupled.py','scripts/reference_hht_block.py',
    'src/afsd_backend/pysph/hht_coupled.py','src/afsd_backend/pysph/hht_path.py',
    'src/afsd_backend/pysph/hht_linear.py','src/afsd_backend/pysph/particles.py',
    'src/afsd_backend/pysph/mechanics.py','src/afsd_backend/pysph/kernel.py',
    'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/implicit_thermal.py',
    'src/afsd_backend/pysph/structural_mobility.py','src/afsd_backend/pysph/output.py',
    'docs/PROCESS_HHT_FIXTURE_CONTRACT.md','docs/PROCESS_HHT_COUPLED_CONTRACT.md',
    'docs/PROCESS_HHT_PATH_CONTRACT.md','docs/PROCESS_HHT_RUN_CONTRACT.md',
    'scripts/restart_hht_probe.py','tests/review/test_hht_restart.py','docs/PROCESS_HHT_RESTART_CONTRACT.md']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def clean(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:clean(a) for k,a in v.items()}
    if isinstance(v,(tuple,list)):return [clean(a) for a in v]
    return v
def write(path,value):
    path=Path(path);pending=path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    pending.replace(path)
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def require(b,message):
    if not b:raise ValueError(message)


def setup():
    f=fixture_data();p=copy.deepcopy(f['parameters'])
    p.update(hht_weight=.9,contact_coulomb_relative_tolerance=.01,
        width=.01,marker_low=.003,marker_high=.004,bar_length=.01,actuator_depth=.001,
        temperature=400.,traverse=0.,feed=0.)
    pa=make_array(f['X'],np.zeros(len(f['X']),dtype=int),p)
    pa.m[:]=f['mass'];pa.h[:]=f['h'];pa.rho[:]=f['rho0'];pa.S[:]=f['S0'].ravel()
    pa.temperature[:]=f['T0'];pa.u[:],pa.v[:],pa.w[:]=f['v0'].T
    # Manufactured passive markers; neither wall law nor mechanics reads them.
    pa.tracer_type[:]=np.arange(len(pa.x))%3;pa.wholebar_tag[:]=pa.tracer_type
    pa.set_output_arrays(list(pa.properties))
    return pa,p


def boundary_for(p,right):
    """Independent implementation of the registered moving penalty planes."""
    def boundary(pa,t):
        cycle=2*np.pi*t/p['end_time'];stroke=.5*p['normal_stroke']*(1-np.cos(cycle))
        speed=p['normal_stroke']*np.pi/p['end_time']*np.sin(cycle)
        planes=[-.5*p['dx']+p['normal_preload']+stroke,
                p['length']+.5*p['dx']-p['normal_preload']-stroke]
        vv=velocities(pa);xx=positions(pa);rows=[]
        for face in (0,1):
            sign=1. if face==0 else -1.;normal=np.array([sign,0.,0.])
            wall=np.array([sign*speed,0.,p['wall_speed']*(-1 if right else 1)])
            indices=np.flatnonzero(pa.x0==(0. if face==0 else p['length']))
            for index in indices:
                penetration=.5*p['dx']+sign*(planes[face]-xx[index,0])
                require(penetration>0,'contact opening needs a new event branch')
                slip=vv[index]-wall;slip=slip-(slip@normal)*normal
                require((1 if right else -1)*slip[2]>.05,'left declared smooth sliding branch')
                fn=p['normal_stiffness']*penetration
                tangent=-p['friction']*fn*slip/np.linalg.norm(slip)
                rows.append(dict(particle_index=int(index),contact_id=f'{face}:{int(pa.particle_id[index])}',
                    normal=normal,normal_force_N=fn,tangential_force_N=tangent,
                    wall_velocity_m_s=wall,wall_temperature_K=p['wall_temperature_K'][face],
                    conductance_W_K=p['contact_conductance']*p['dx']**2,
                    penetration_m=penetration,normal_potential_J=.5*p['normal_stiffness']*penetration**2,
                    normal_stiffness_N_m=p['normal_stiffness']))
        return dict(force_definition='NONCONTACT_APPLIED_FORCE_ONLY',force_N=np.zeros((len(pa.x),3)),contacts=rows)
    return boundary


def inventory():
    watched=set(SOURCES)|{p.relative_to(ROOT).as_posix() for p in (ROOT/'src/afsd_backend/pysph').glob('*.py')}
    return {s:sha(ROOT/s) for s in sorted(watched)}


def verify_receipts(preflight,reference,sources):
    a=read(preflight);b=read(reference)
    require(a['status']=='PASS_HHT_COMPONENT_PREFLIGHT_NOT_AFSD','component preflight required')
    for s,h in a['source_sha256'].items():require(sha(ROOT/s)==h,'stale preflight '+s)
    for s,h in sources.items():require(a['source_sha256'].get(s)==h,'source missing preflight '+s)
    for s,h in a['artifact_sha256'].items():require(sha(Path(preflight).parent/s)==h,'preflight artifact drift')
    require(b['status']=='PASS_INDEPENDENT_ELASTIC_REFERENCE','independent reference required')
    for s,h in b['source_sha256'].items():require(sha(ROOT/s)==h,'reference source drift')
    return dict(preflight_sha256=sha(preflight),reference_sha256=sha(reference),
        reference_history_sha256=sha(Path(reference).parent/'reference_history.npz'))


def frame(pa,state,integrals):
    face=np.zeros((2,3))
    for c in state['boundary_snapshot']['contacts']:
        face[int(c['contact_id'].split(':')[0])]+=c['normal_force_N']*c['normal']+c['tangential_force_N']
    return dict(time_s=state['time_s'],position_m=positions(pa),velocity_m_s=velocities(pa),
        stress_Pa=pa.S.reshape(-1,3,3)-pressure(pa.rho,2700.,70e9/3)[:,None,None]*np.eye(3),
        density_kg_m3=pa.rho.copy(),temperature_K=pa.temperature.copy(),ep=pa.ep.copy(),
        acceleration_m_s2=state['acceleration_m_s2'].copy(),force_N=state['force_N'].copy(),
        face_force_N=face,integrals=copy.deepcopy(integrals))


def add_integrals(total,old,new,diagnostic,p):
    h=diagnostic['dt_s'];a0=old['acceleration_m_s2'];a1=new['acceleration_m_s2']
    contacts0={c['contact_id']:c for c in old['boundary_snapshot']['contacts']}
    contacts1={c['contact_id']:c for c in new['boundary_snapshot']['contacts']}
    qnode=diagnostic['thermal']['particle_energy_J']['contact_in_J']
    for audit in diagnostic['contact']['contacts']:
        cid=audit['contact_id'];c0=contacts0[cid];c1=contacts1[cid];side=int(cid.split(':')[0]);idx=c0['particle_index']
        total['friction_raw'][side]+=audit['friction_heat_raw_J']
        total['reservoir_heat'][side]+=qnode[idx]
        for node in audit['quadrature']:
            s=node['fraction'];wall=(1-s)*c0['wall_velocity_m_s']+s*c1['wall_velocity_m_s']
            normal=(1-s)*c0['normal_force_N']*c0['normal']+s*c1['normal_force_N']*c1['normal']
            impulse=node['impulse_Ns']+.5*h*normal
            total['impulse'][side]+=impulse
            total['wall_work'][side]+=float(impulse@wall)
    total['algorithmic_filter_work_J']+=float(np.sum(diagnostic['kinetic']['algorithmic_filter_work_J']))
    total['root_residual_work_J']+=float(np.sum(diagnostic['kinetic']['root_residual_work_J']))


def physical_checks(pa,old,initial,state,diag,p):
    require(np.array_equal(pa.particle_id,initial.particle_id) and np.array_equal(pa.m,initial.m),'ID/mass changed')
    require(np.max(np.linalg.norm(positions(pa)-positions(initial),axis=1))<=5e-6,'benchmark displacement validity')
    require(np.max(np.abs(pa.rho/p['rho0']-1))<=.005,'benchmark density validity')
    require(max(c['penetration_m'] for c in state['boundary_snapshot']['contacts'])<=5e-6,'benchmark penetration validity')
    require(np.all(pa.ep==0),'elastic benchmark yielded')
    kin=diag['kinetic'];ws=np.abs(kin['physical_linear_force_work_J']).sum()+np.abs(kin['algorithmic_filter_work_J']).sum()
    require(np.max(np.abs(kin['kinetic_identity_residual_J']))<=1e-12+1e-9*ws,'kinetic identity')
    js=np.abs(kin['physical_impulse_Ns']).sum()
    require(np.max(np.abs(kin['momentum_identity_residual_Ns']))<=1e-12+1e-9*js,'momentum identity')
    thermal=diag['thermal']['particle_energy_J']
    actual_balance=pa.m*p['cp']*(pa.temperature-old.temperature)-thermal['total_heat_input_J']
    require(np.max(np.abs(actual_balance))<=1e-12+1e-9*np.abs(thermal['total_heat_input_J']).sum(),'saved temperature thermal identity')
    # Tangential wall speed is constant within each registered half. Certify
    # its signed quadratic slip margin, including the interior extremum.
    right=state['boundary_snapshot']['contacts'][0]['wall_velocity_m_s'][2]<0;minimum_slip=float('inf')
    for row in diag['contact']['contacts']:
        u=np.asarray(row['slip_coefficients_m_s'])[:,2]*(1 if right else -1)
        samples=[0.,1.]
        if u[2]!=0 and 0 < -u[1]/(2*u[2]) < 1:samples.append(-u[1]/(2*u[2]))
        minimum_slip=min(minimum_slip,min(u[0]+s*u[1]+s*s*u[2] for s in samples))
    require(minimum_slip>.05,'interior registered sliding branch margin')
    # Bound the actual cubic particle path against the whole cosine wall
    # stroke [0, normal_stroke]; no endpoint-only opening certificate.
    h=diag['dt_s']
    # Endpoint displacement and velocities determine this cubic exactly.
    x0=positions(old);v0=velocities(old);v1=velocities(pa)
    c0=x0[:,0];c1=h*v0[:,0];dx=pa.x-old.x
    c2=3*dx-h*(2*v0[:,0]+v1[:,0]);c3=-2*dx+h*(v0[:,0]+v1[:,0])
    lower=float('inf');upper=0.
    for row in state['boundary_snapshot']['contacts']:
        i=row['particle_index'];side=int(row['contact_id'].split(':')[0]);sgn=1 if side==0 else -1
        fractions=[0.,1.]
        for root in np.polynomial.polynomial.polyroots([c1[i],2*c2[i],3*c3[i]]):
            if abs(root.imag)<1e-12 and 0<root.real<1:fractions.append(float(root.real))
        signed_x=[sgn*(c0[i]+s*c1[i]+s*s*c2[i]+s**3*c3[i]) for s in fractions]
        base=-.5*p['dx']+p['normal_preload'] if side==0 else p['length']+.5*p['dx']-p['normal_preload']
        lower=min(lower,.5*p['dx']+sgn*base-max(signed_x))
        upper=max(upper,.5*p['dx']+sgn*base-min(signed_x)+p['normal_stroke'])
    require(lower>0 and upper<=5e-6,'whole-path conservative overlap bounds')
    law=FKS(**{k:p[k] for k in FKS.__dataclass_fields__ if k in p})
    margin=law.yield_stress(pa.ep,np.zeros(len(pa.x)),pa.temperature)-equivalent(pa.S.reshape(-1,3,3))
    diag['publication_checks']=dict(saved_temperature_balance_J=actual_balance,
        minimum_signed_slip_m_s=minimum_slip,whole_path_overlap_lower_m=lower,whole_path_overlap_upper_m=upper,
        minimum_low_rate_yield_margin_Pa=float(np.min(margin)),
        maximum_temperature_change_K=float(np.max(np.abs(pa.temperature-initial.temperature))))


def save_stage(directory,pa,state,p,diag,integrals):
    directory.mkdir()
    checkpoint(directory/'accepted.npz',pa,state)
    loaded,stored=restore(directory/'accepted.npz');checked,again=restore_checkpoint(loaded,stored,p)
    require(particle_hash(checked)==particle_hash(pa) and again['state_hash']==state['state_hash'],'full checkpoint readback')
    scalars=dict(particle_id=pa.particle_id,temperature_K=pa.temperature,density_kg_m3=pa.rho,
        pressure_Pa=pressure(pa.rho,p['rho0'],p['bulk']),equivalent_stress_Pa=equivalent(pa.S.reshape(-1,3,3)),
        plastic_strain=pa.ep,tracer_type=pa.tracer_type,mass_kg=pa.m)
    sigma=pa.S.reshape(-1,3,3)-pressure(pa.rho,p['rho0'],p['bulk'])[:,None,None]*np.eye(3)
    for i,j,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):scalars['stress_'+name+'_Pa']=sigma[:,i,j]
    write_point_cloud(directory/'accepted.vtk',positions(pa).tolist(),clean(scalars),
        clean(dict(velocity_m_s=velocities(pa),acceleration_m_s2=state['acceleration_m_s2'])),title=CLASS)
    with gzip.open(directory/'diagnostic.json.gz','wt',encoding='utf-8') as stream:json.dump(clean(diag),stream,allow_nan=False)
    write(directory/'integrals.json',integrals)
    return dict(time_s=state['time_s'],checkpoint=directory.name+'/accepted.npz',
        checkpoint_sha256=sha(directory/'accepted.npz'),state_hash=state['state_hash'])


def trajectory(directory,steps,maximum_steps=None):
    directory.mkdir();pa,p=setup();initial=pa.extract_particles(np.arange(len(pa.x),dtype=np.uint32))
    state=initialize(pa,p,boundary_for(p,False));h=p['end_time']/steps
    integrals=dict(friction_raw=np.zeros(2),reservoir_heat=np.zeros(2),wall_work=np.zeros(2),
        impulse=np.zeros((2,3)),algorithmic_filter_work_J=0.,root_residual_work_J=0.)
    frames=[frame(pa,state,integrals)];history=[];metrics=[];started=time.perf_counter();sources=inventory()
    limit=steps if maximum_steps is None else min(steps,maximum_steps)
    try:
        for index in range(limit):
            right=index>=steps//2;old=state
            candidate,candidate_state,diag=advance(pa,state,p,boundary_for(p,right),h)
            physical_checks(candidate,pa,initial,candidate_state,diag,p)
            candidate_integrals=copy.deepcopy(integrals);add_integrals(candidate_integrals,old,candidate_state,diag,p)
            left_frame=None
            if index+1==steps//2:
                left=copy.deepcopy(candidate_state);left_frame=frame(candidate,candidate_state,candidate_integrals)
                candidate_state,event=restart_at_control_event(candidate,candidate_state,p,boundary_for(p,True))
                diag['control_switch']=event;diag['left_state_hash_at_switch']=left['state_hash']
            require(inventory()==sources,'source drift before candidate save')
            entry=save_stage(directory/f'stage_{index+1:03d}',candidate,candidate_state,p,diag,candidate_integrals)
            candidate_frame=frame(candidate,candidate_state,candidate_integrals)
            if left_frame is not None:candidate_frame['switch_left_face_force_N']=left_frame['face_force_N']
            require(inventory()==sources,'source drift before accepted publication')
            next_history=history+[entry];write(directory/'accepted_history.json',next_history)
            # Publication is the sole acceptance transaction. A saved but
            # unpublished NPZ is not a restart point, even if locally valid.
            pa,state,integrals,history=candidate,candidate_state,candidate_integrals,next_history
            frames.append(candidate_frame);metrics.append(dict(step=index+1,evaluations=diag['evaluations'],
                step_seconds=diag['elapsed_s'],maximum_contact_defect=diag['maximum_sampled_coulomb_relative_defect'],
                normal_temporal_residual_J=diag['normal_temporal_residual_J'],residuals=diag['residuals'],
                publication_checks=diag['publication_checks']))
            print(json.dumps(clean(metrics[-1])),flush=True)
        write(directory/'frames.json',frames);write(directory/'step_metrics.json',metrics)
        require(abs(state['time_s']-limit*h)<1e-19,'trajectory common end time')
        return dict(status='COMPLETE_LOCAL_TRAJECTORY_NOT_AFSD',steps=limit,planned_steps=steps,
            time_s=state['time_s'],wall_seconds=time.perf_counter()-started,metrics=metrics,frames=frames)
    except Exception as exc:
        write(directory/'REJECTED.json',dict(accepted=False,restart_allowed=False,error=repr(exc),
            diagnostic=getattr(exc,'diagnostics',None),last_accepted_time_s=state['time_s'],published_steps=len(history)))
        raise


def compare(runs,reference):
    ref=read(reference);history=np.load(Path(reference).parent/'reference_history.npz',allow_pickle=False)
    states=history['states'];times=history['time_s'];vol=fixture_data()['mass']/2700.;rows=[]
    n=len(vol)
    def rms(a):return float(np.sqrt(np.sum(vol.reshape((n,)+(1,)*(a.ndim-1))*a*a)/vol.sum()))
    for run in runs:
        maximum=dict(position_m=0.,velocity_m_s=0.,stress_Pa=0.,temperature_K=0.,face_force_N=0.,
            impulse=0.,wall_work=0.,friction_raw=0.,reservoir_heat=0.)
        for f in run['frames']:
            idx=int(np.argmin(np.abs(times-f['time_s'])));require(abs(times[idx]-f['time_s'])<1e-19,'reference common time')
            y=states[idx];x=y[:3*n].reshape(n,3);v=y[3*n:6*n].reshape(n,3);S=y[6*n:15*n].reshape(n,3,3)
            rho=y[15*n:16*n];T=y[16*n:17*n];sigma=S-(70e9/3)*(rho/2700-1)[:,None,None]*np.eye(3)
            errors=dict(position_m=rms(f['position_m']-x),velocity_m_s=rms(f['velocity_m_s']-v),
                stress_Pa=rms(f['stress_Pa']-sigma),temperature_K=float(np.max(np.abs(f['temperature_K']-T))),
                face_force_N=float(np.max(np.abs(f['face_force_N']-ref['records'][idx]['face_force_N']))))
            if 'switch_left_face_force_N' in f:
                errors['face_force_N']=max(errors['face_force_N'],float(np.max(np.abs(
                    f['switch_left_face_force_N']-ref['switch_left_record']['face_force_N']))))
            for k in ('impulse','wall_work','friction_raw','reservoir_heat'):
                a=np.asarray(f['integrals'][k]).ravel();b=np.asarray(ref['records'][idx][k]).ravel()
                errors[k]=float(np.max(np.abs(a-b)))
            for k,value in errors.items():maximum[k]=max(maximum[k],value)
        rows.append(dict(steps=run['steps'],errors=maximum))
    bounds=dict(position_m=1e-9,velocity_m_s=2e-4,stress_Pa=7000.,temperature_K=1e-7,
        face_force_N=.7,impulse=1.4e-7,wall_work=2e-7,friction_raw=2e-7,reservoir_heat=2e-8)
    uncertainty=ref['comparison']['differences'];decrease={}
    for k in bounds:
        e=[row['errors'][k] for row in rows]
        decrease[k]=all(b<=.7*a if a>10*uncertainty[k] else b<=10*uncertainty[k] for a,b in zip(e[:-1],e[1:]))
    return dict(rows=rows,finest_bounds=bounds,refinement_pass=decrease,
        passed=all(rows[-1]['errors'][k]<=bounds[k] and decrease[k] for k in bounds),
        oracle_scope='INDEPENDENT_CURRENT_GEOMETRY_SEMIDISCRETE_NOT_CONTINUUM_VALIDATION')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--preflight',type=Path,required=True);ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--mode',choices=['pilot','sweep'],required=True);ap.add_argument('--pilot',type=Path)
    a=ap.parse_args();out=a.output.resolve();require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);sources=inventory();started=time.perf_counter()
    def timeout():
        write(out/'HARD_TIMEOUT.json',dict(status='WALL_BUDGET_STOP',accepted=False));os._exit(124)
    timer=threading.Timer(60. if a.mode=='pilot' else 55.,timeout);timer.daemon=True;timer.start()
    report=dict(classification=CLASS,status='RUNNING',production_authorized=False,afsd_progress=False,
        mode=a.mode,source_sha256=sources,new_afsd_steps=0)
    try:
        receipts=verify_receipts(a.preflight,a.reference,sources);report['receipts']=receipts
        _,p=setup();write(out/'resolved_config.json',p);report['config_sha256']=sha(out/'resolved_config.json')
        if a.mode=='pilot':
            run=trajectory(out/'pilot',8,1);report.update(status='PASS_LOCAL_COST_PILOT_NOT_COMMON_END_TIME',
                pilot_wall_seconds=run['wall_seconds'],pilot_evaluations=sum(v['evaluations'] for v in run['metrics']))
        else:
            require(a.pilot is not None,'same-source pilot required');pilot=read(a.pilot)
            require(pilot['status']=='PASS_LOCAL_COST_PILOT_NOT_COMMON_END_TIME' and pilot['source_sha256']==sources
                and pilot['receipts']==receipts,'stale/failed pilot')
            estimate=2*56*max(pilot['pilot_wall_seconds'],.01)+10
            report['forecast_seconds']=estimate;require(estimate<=50.,'NOT_RUN_SWEEP_COST_EXCEEDS_REGISTERED_LIMIT')
            runs=[trajectory(out/f'grid_{steps:03d}',steps) for steps in (8,16,32)]
            comparison=compare(runs,a.reference);write(out/'comparison.json',comparison)
            report.update(status='PASS_LOCAL_HHT_TEMPORAL_BENCHMARK_NOT_AFSD' if comparison['passed'] else 'FAIL_LOCAL_TEMPORAL_ACCURACY',
                comparison=comparison,completed_grids=[dict(steps=r['steps'],time_s=r['time_s'],wall_seconds=r['wall_seconds']) for r in runs])
        require(inventory()==sources,'source drift during run')
        code=0 if report['status'].startswith('PASS') else 2
    except Exception as exc:report.update(status='FAILED_LOCAL_BENCHMARK',error=repr(exc));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    write(out/'report.json',report);print(json.dumps(clean(report),allow_nan=False),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

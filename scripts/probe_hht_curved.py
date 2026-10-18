"""Bounded curved-contact slow scout; never advances the AFSD checkpoint."""
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
from hht_slow_inputs import inputs,schedule,SCOUT_END_S
CLASSIFICATION='LOCAL_CURVED_CONTACT_IMPULSE_HHT_SLOW_STARTUP_NOT_AFSD'
from probe_hht_coupled import clean,sha,write,read,require,frame
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.mechanics import positions,velocities,pressure
from afsd_backend.pysph.material import equivalent
from afsd_backend.pysph.hht_curved_coupled import initialize,advance,restore_checkpoint
from afsd_backend.pysph.output import checkpoint,restore
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.application import memory_bytes
from afsd_scaffold.vtk_legacy import write_point_cloud

SOURCES=['scripts/probe_hht_curved.py','scripts/probe_hht_curved_regression.py','docs/PROCESS_CURVED_CONTRACT.md','docs/PROCESS_CURVED_DERIVATION.md','scripts/probe_hht_slow.py','scripts/hht_slow_inputs.py',
    'scripts/reference_hht_slow.py','scripts/reference_hht_block.py',
    'scripts/probe_hht_coupled.py','docs/PROCESS_SLOW_CONTRACT.md']


def inventory():
    names=set(SOURCES)|{p.relative_to(ROOT).as_posix() for p in (ROOT/'src/afsd_backend/pysph').glob('*.py')}
    return {name:sha(ROOT/name) for name in sorted(names)}


from probe_hht_slow import setup, boundary_for
from afsd_backend.pysph import hht_path


def slow_gap_certificate(pa,acc0,acc1,dt,p,left,right):
    """Conservative cubic particle/analytic wall enclosure, independent of nodes."""
    x0=positions(pa);v0=velocities(pa)
    x1,v1=hht_path.endpoint(x0,v0,acc0,acc1,dt,p['hht_weight'])
    delta=x1[:,0]-x0[:,0]
    c0=x0[:,0];c1=dt*v0[:,0]
    c2=3*delta-dt*(2*v0[:,0]+v1[:,0]);c3=-2*delta+dt*(v0[:,0]+v1[:,0])
    lo=float('inf');hi=0.
    for row in left['contacts']:
        idx=row['particle_index'];face=int(row['contact_id'].split(':')[0]);sgn=1 if face==0 else -1
        samples=[0.,1.]
        for r in np.polynomial.polynomial.polyroots([c1[idx],2*c2[idx],3*c3[idx]]):
            if abs(r.imag)<1e-12 and 0<r.real<1:samples.append(float(r.real))
        xx=[sgn*(c0[idx]+s*c1[idx]+s*s*c2[idx]+s**3*c3[idx]) for s in samples]
        base=-.5*p['dx']+p['normal_preload'] if face==0 else p['length']+.5*p['dx']-p['normal_preload']
        lo=min(lo,.5*p['dx']+sgn*base-max(xx))
        hi=max(hi,.5*p['dx']+sgn*base-min(xx)+p['normal_stroke'])
    require(lo>0 and hi<=5e-6,'conservative whole-path contact gap domain')
    return dict(minimum_overlap_m=lo,maximum_overlap_m=hi,fixed_normals=True,
        affine_tangential_wall_velocity=True,
        source='CUBIC_PARTICLE_X_EXTREMA_PLUS_ANALYTIC_5HZ_STROKE_RANGE',
        definition='FIXED_X_PLANES_CONSTANT_PLUS_Z_TANGENTIAL_WALL_SPEED')


def add_integrals(total,old,new,diagnostic,p):
    """Actual integrated path forces, not endpoint interpolation or two-node weights."""
    qwall=diagnostic['thermal']['particle_energy_J']['contact_in_J']
    seen=set()
    for row in diagnostic['contact']['contacts']:
        side=int(row['contact_id'].split(':')[0]);idx=row['particle_index']
        require(idx not in seen,'disjoint-face reservoir attribution required')
        seen.add(idx)
        total['friction_raw'][side]+=row['friction_heat_raw_J']
        total['reservoir_heat'][side]+=qwall[idx]
        total['impulse'][side]+=row['total_actual_impulse_Ns']
        total['wall_work'][side]+=row['actual_wall_work_J']
    total['algorithmic_filter_work_J']+=float(np.sum(diagnostic['kinetic']['algorithmic_filter_work_J']))
    total['root_residual_work_J']+=float(np.sum(diagnostic['kinetic']['root_residual_work_J']))


def checks(candidate,old,initial,state,diag,p):
    require(np.array_equal(candidate.particle_id,initial.particle_id) and np.array_equal(candidate.m,initial.m),'ID/mass changed')
    x=positions(candidate);reference=positions(initial);disp=np.linalg.norm(x-reference,axis=1)
    require(float(disp.max())<=p['maximum_translation_m'],'manufactured translation domain')
    i,j=np.triu_indices(27,1);r0=np.linalg.norm(reference[i]-reference[j],axis=1)
    shape=float(np.max(np.abs(np.linalg.norm(x[i]-x[j],axis=1)/r0-1)))
    require(shape<=p['maximum_pair_relative_change'],'manufactured small-strain domain')
    require(np.max(np.abs(candidate.rho/p['rho0']-1))<=.005,'density domain')
    require(np.all(candidate.ep==0),'positive plastic response outside registered scope')
    kin=diag['kinetic'];scale=np.abs(kin['physical_force_work_J']).sum()+np.abs(kin['algorithmic_filter_work_J']).sum()
    require(np.max(np.abs(kin['kinetic_identity_residual_J']))<=1e-12+1e-9*scale,'kinetic identity')
    require(np.max(np.abs(kin['momentum_identity_residual_Ns']))<=1e-12+1e-9*np.abs(kin['physical_impulse_Ns']).sum(),'momentum identity')
    heat=diag['thermal']['particle_energy_J']['total_heat_input_J']
    balance=candidate.m*p['cp']*(candidate.temperature-old.temperature)-heat
    require(np.max(np.abs(balance))<=1e-12+1e-9*np.abs(heat).sum(),'saved temperature heat identity')
    # Conservative whole-cubic-path overlap enclosure, not endpoint-only.
    h=diag['dt_s'];v0=velocities(old);v1=velocities(candidate);dx=candidate.x-old.x
    c0=old.x;c1=h*v0[:,0];c2=3*dx-h*(2*v0[:,0]+v1[:,0]);c3=-2*dx+h*(v0[:,0]+v1[:,0])
    low=float('inf');high=0.;minimum_speed=float('inf')
    for contact in state['boundary_snapshot']['contacts']:
        index=contact['particle_index'];face=int(contact['contact_id'].split(':')[0]);sgn=1 if face==0 else -1
        samples=[0.,1.]
        for root in np.polynomial.polynomial.polyroots([c1[index],2*c2[index],3*c3[index]]):
            if abs(root.imag)<1e-12 and 0<root.real<1:samples.append(float(root.real))
        sx=[sgn*(c0[index]+s*c1[index]+s*s*c2[index]+s**3*c3[index]) for s in samples]
        base=-.5*p['dx']+p['normal_preload'] if face==0 else p['length']+.5*p['dx']-p['normal_preload']
        low=min(low,.5*p['dx']+sgn*base-max(sx));high=max(high,.5*p['dx']+sgn*base-min(sx)+p['normal_stroke'])
        v=velocities(candidate)[index]-contact['wall_velocity_m_s'];v-=np.dot(v,contact['normal'])*contact['normal']
        minimum_speed=min(minimum_speed,float(np.linalg.norm(v)))
    require(low>0 and high<=5e-6,'conservative full-path contact overlap domain')
    diag['publication_checks']=dict(maximum_translation_m=float(disp.max()),maximum_pair_relative_change=shape,
        whole_path_overlap_lower_m=low,whole_path_overlap_upper_m=high,
        minimum_endpoint_vector_slip_m_s=minimum_speed,saved_temperature_balance_J=balance,
        signed_z_reversal_is_failure=False)


def save_stage(directory,pa,state,p,diag,integrals):
    directory.mkdir();checkpoint(directory/'accepted.npz',pa,state)
    loaded,stored=restore(directory/'accepted.npz');again,verified=restore_checkpoint(loaded,stored,p)
    require(particle_hash(again)==particle_hash(pa) and verified['state_hash']==state['state_hash'],'full checkpoint readback')
    scalar=dict(particle_id=pa.particle_id,mass_kg=pa.m,tracer_type=pa.tracer_type,
        temperature_K=pa.temperature,density_kg_m3=pa.rho,pressure_Pa=pressure(pa.rho,p['rho0'],p['bulk']),
        equivalent_stress_Pa=equivalent(pa.S.reshape(-1,3,3)),plastic_strain=pa.ep)
    sigma=pa.S.reshape(-1,3,3)-scalar['pressure_Pa'][:,None,None]*np.eye(3)
    for i,j,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):scalar['stress_'+name+'_Pa']=sigma[:,i,j]
    write_point_cloud(directory/'accepted.vtk',positions(pa).tolist(),clean(scalar),
        clean(dict(velocity_m_s=velocities(pa),acceleration_m_s2=state['acceleration_m_s2'])),
        title=CLASSIFICATION+' time_s='+str(state['time_s']))
    with gzip.open(directory/'diagnostic.json.gz','wt',encoding='utf-8') as stream:json.dump(clean(diag),stream,allow_nan=False)
    write(directory/'integrals.json',integrals)
    return dict(time_s=state['time_s'],checkpoint=directory.name+'/accepted.npz',
        checkpoint_sha256=sha(directory/'accepted.npz'),state_hash=state['state_hash'])


def receipts(preflight,reference,sources):
    check=read(preflight);ref=read(reference)
    require(check['status']=='PASS_CURVED_COMPONENT_PREFLIGHT_NOT_PROCESS','new component preflight required')
    for name,h in check['source_sha256'].items():require(sha(ROOT/name)==h,'stale preflight '+name)
    for name,h in sources.items():require(check['source_sha256'].get(name)==h,'missing preflight source '+name)
    for name,h in check['artifact_sha256'].items():require(sha(preflight.parent/name)==h,'changed preflight artifact')
    require(ref['status'] in ('PASS_COMPLETE_REFERENCE','PASS_PREFIX_REFERENCE_EVENT_BLOCKED'),'validated independent reference prefix required')
    require(ref['prefix_comparison']['passed'],'reference uncertainty failed')
    for name,h in ref['source_sha256'].items():require(sha(ROOT/name)==h,'stale independent source '+name)
    expected_config=hashlib.sha256(json.dumps(clean(inputs()),sort_keys=True).encode()).hexdigest()
    require(ref['configuration_sha256']==expected_config,'reference physical inputs differ')
    require('reference_history.npz' in ref['artifact_sha256'],'reference history missing artifact binding')
    for name,h in ref['artifact_sha256'].items():
        target=(reference.parent/name).resolve()
        require(target.is_relative_to(reference.parent.resolve()) and sha(target)==h,'independent artifact drift '+name)
    require(ref['certified_common_end_time_s']>=125e-9,'insufficient independently resolved startup')
    return ref,dict(preflight_sha256=sha(preflight),reference_sha256=sha(reference),
        reference_history_sha256=sha(reference.parent/'reference_history.npz'))


def trajectory(directory,refinement,reference_end,sources,maximum_steps=None):
    directory.mkdir();pa,p=setup();initial=pa.extract_particles(np.arange(27,dtype=np.uint32));boundary_journal=[]
    boundary=boundary_for(p,boundary_journal);state=initialize(pa,p,boundary)
    integrals=dict(friction_raw=np.zeros(2),reservoir_heat=np.zeros(2),wall_work=np.zeros(2),
        impulse=np.zeros((2,3)),algorithmic_filter_work_J=0.,root_residual_work_J=0.)
    history=[];metrics=[];frames=[frame(pa,state,integrals)];plan=schedule(refinement);started=time.perf_counter()
    if maximum_steps is not None:plan=plan[:maximum_steps]
    try:
        for index,(begin,end) in enumerate(plan):
            require(abs(state['time_s']-begin)<1e-19,'schedule/published time mismatch')
            boundary_journal.clear();candidate,candidate_state,diag=advance(pa,state,p,boundary,end-begin,slow_gap_certificate)
            checks(candidate,pa,initial,candidate_state,diag,p)
            if end>reference_end+1e-19:
                checkpoint(directory/'REJECTED_UNREFERENCED.npz',candidate,dict(classification='REJECTED',
                    accepted=False,restart_allowed=False,candidate_state=candidate_state,
                    reason='UNREFERENCED_INTERVAL',last_accepted_state_hash=state['state_hash']))
                raise ValueError('UNREFERENCED_INTERVAL: candidate not published beyond certified reference prefix')
            cumulative=copy.deepcopy(integrals);add_integrals(cumulative,state,candidate_state,diag,p)
            require(inventory()==sources,'source drift before candidate save')
            item=save_stage(directory/f'stage_{index+1:03d}',candidate,candidate_state,p,diag,cumulative)
            require(inventory()==sources,'source drift before publication')
            updated=history+[item];write(directory/'accepted_history.json',updated)
            pa,state,integrals,history=candidate,candidate_state,cumulative,updated
            frames.append(frame(pa,state,integrals));metrics.append(dict(step=index+1,time_s=state['time_s'],dt_s=end-begin,
                evaluations=diag['evaluations'],step_seconds=diag['elapsed_s'],residuals=diag['residuals'],
                maximum_contact_defect=diag['maximum_sampled_coulomb_relative_defect'],
                normal_temporal_residual_J=diag['normal_temporal_residual_J'],checks=diag['publication_checks']))
            write(directory/'frames.json',frames);write(directory/'step_metrics.json',metrics)
            print(json.dumps(dict(step=index+1,time_s=state['time_s'],evaluations=diag['evaluations'])),flush=True)
        return dict(status='COMPLETE_LOCAL_SCOUT_NOT_AFSD' if maximum_steps is None else 'PASS_LOCAL_COST_PILOT',
            time_s=state['time_s'],published_steps=len(history),planned_steps=len(schedule(refinement)),
            wall_seconds=time.perf_counter()-started,metrics=metrics)
    except Exception as error:
        write(directory/'REJECTED.json',dict(classification='REJECTED',accepted=False,restart_allowed=False,
            error=repr(error),last_accepted_time_s=state['time_s'],last_accepted_state_hash=state['state_hash'],
            published_steps=len(history),candidate_interval_s=[begin,end],reference_certified_end_time_s=reference_end,
            diagnostic=getattr(error,'diagnostics',None)))
        with gzip.open(directory/'rejected_boundary_evaluations.json.gz','wt',encoding='utf-8') as stream:
            json.dump(clean(boundary_journal),stream,allow_nan=False)
        return dict(status='STOPPED_LOCAL_SCOUT_NOT_PROCESS_FEASIBLE',error=repr(error),time_s=state['time_s'],
            published_steps=len(history),planned_steps=len(schedule(refinement)),
            wall_seconds=time.perf_counter()-started,metrics=metrics,failed_interval_s=[begin,end])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--preflight',type=Path,required=True);ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--mode',choices=['pilot','S0'],required=True);ap.add_argument('--pilot',type=Path)
    ap.add_argument('--regression',type=Path,required=True)
    a=ap.parse_args();out=a.output.resolve();require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);started=time.perf_counter();sources=inventory();attempt_started=False
    def timeout():
        write(out/'HARD_TIMEOUT.json',dict(status='WALL_BUDGET_STOP',accepted=False,restart_allowed=False));os._exit(124)
    timer=threading.Timer(55.,timeout);timer.daemon=True;timer.start()
    report=dict(classification=CLASSIFICATION,status='NOT_RUN_GATES_PENDING',mode=a.mode,source_sha256=sources,
        production_authorized=False,new_afsd_steps=0,process_transport=False,desired_1ms_status='NOT_RUN_REFERENCE_COST')
    try:
        ref,receipt=receipts(a.preflight,a.reference,sources)
        regression=read(a.regression)
        require(regression['status']=='PASS_PRIVATE_SAVED_RED_REGRESSION_NOT_TRAJECTORY' and
            regression['accepted'] is False and regression['restart_allowed'] is False and
            regression['new_afsd_steps']==0,'passed private saved RED regression required')
        require(regression['source_sha256']==sources and regression['receipts']==receipt,'stale saved RED regression')
        require(abs(regression['candidate_time_s']-2.5e-6)<1e-19,'wrong regression interval')
        for name,digest in regression['artifact_sha256'].items():
            target=(a.regression.parent/name).resolve()
            require(target.is_relative_to(a.regression.parent.resolve()) and sha(target)==digest,'changed private regression artifact')
        receipt['regression_sha256']=sha(a.regression);report['receipts']=receipt
        _,p=setup();write(out/'resolved_config.json',p);report['configuration_sha256']=sha(out/'resolved_config.json')
        refinement=0 if a.mode=='pilot' else int(a.mode[1]);steps=len(schedule(refinement))
        if a.mode!='pilot':
            require(a.pilot is not None,'measured same-source pilot required');pilot=read(a.pilot)
            require(pilot['status']=='PASS_LOCAL_COST_PILOT' and pilot['source_sha256']==sources and pilot['receipts']==receipt,'stale/failed pilot')
            estimate=2*steps*max(pilot['run']['wall_seconds'],.01)+10
            report['forecast_seconds']=estimate;require(estimate<=50,'NOT_RUN_COST_EXCEEDS_REGISTERED_LIMIT')
        attempt_started=True
        run=trajectory(out/'trajectory',refinement,ref['certified_common_end_time_s'],sources,1 if a.mode=='pilot' else None)
        report.update(status=run['status'],run=run,reference_event_status=ref['event_status'],
            comparison_status='NOT_RUN_NO_THREE_COMPLETED_COMMON_END_TIMES')
        require(inventory()==sources,'source drift during scout')
        code=0 if report['status'] in ('COMPLETE_LOCAL_SCOUT_NOT_AFSD','PASS_LOCAL_COST_PILOT') else 2
    except Exception as error:
        report.update(status='FAIL_SCOUT_EXECUTION_OR_INTEGRITY' if attempt_started else 'NOT_RUN_GATE_REJECTED',error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={path.relative_to(out).as_posix():sha(path) for path in out.rglob('*') if path.is_file()})
    write(out/'report.json',report);print(json.dumps(clean({k:v for k,v in report.items() if k not in ('source_sha256','artifact_sha256','run')})),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

"""One bounded local continuation from the PUBLISHED5us parent, not AFSD."""
from pathlib import Path
import argparse
import copy
import gzip
import hashlib
import json
import os
import threading
import time
import numpy as np
import probe_hht_curved as base
from afsd_backend.pysph import hht_stop_stage as stage
from afsd_backend.pysph.hht_stop_flow import certify_interval,discover_events,locate_duration,restart_sliding
from afsd_backend.pysph.output import checkpoint,restore
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.mechanics import positions,velocities,pressure
from afsd_backend.pysph.material import equivalent
from afsd_backend.pysph.application import memory_bytes
from afsd_scaffold.vtk_legacy import write_point_cloud
from hht_stop_boundary import boundary_for,actual_wall_accelerations
from hht_slow_inputs import inputs

ROOT=base.ROOT


def save_stage(directory,pa,state,p,diagnostic,integrals):
    directory.mkdir();checkpoint(directory/'accepted.npz',pa,state)
    loaded,stored=restore(directory/'accepted.npz');again,verified=stage.restore_checkpoint(loaded,stored,p)
    base.require(particle_hash(again)==particle_hash(pa) and verified['state_hash']==state['state_hash'],'checkpoint readback')
    scalar=dict(particle_id=pa.particle_id,mass_kg=pa.m,tracer_type=pa.tracer_type,
        temperature_K=pa.temperature,density_kg_m3=pa.rho,pressure_Pa=pressure(pa.rho,p['rho0'],p['bulk']),
        equivalent_stress_Pa=equivalent(pa.S.reshape(-1,3,3)),plastic_strain=pa.ep)
    sigma=pa.S.reshape(-1,3,3)-scalar['pressure_Pa'][:,None,None]*np.eye(3)
    for i,j,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):scalar['stress_'+name+'_Pa']=sigma[:,i,j]
    write_point_cloud(directory/'accepted.vtk',positions(pa).tolist(),base.clean(scalar),
        base.clean(dict(velocity_m_s=velocities(pa),acceleration_m_s2=state['acceleration_m_s2'])),
        title='LOCAL_LOCATED_STOP_NOT_AFSD time_s='+str(state['time_s']))
    with gzip.open(directory/'diagnostic.json.gz','wt',encoding='utf-8') as stream:json.dump(base.clean(diagnostic),stream,allow_nan=False)
    base.write(directory/'integrals.json',integrals)
    return dict(time_s=state['time_s'],checkpoint=directory.name+'/accepted.npz',
        checkpoint_sha256=base.sha(directory/'accepted.npz'),state_hash=state['state_hash'])


def validate_receipt(path,status=None):
    value=base.read(path)
    if status is not None:base.require(value['status'] in status,'required receipt failed '+str(path))
    for name,digest in value['source_sha256'].items():base.require(base.sha(ROOT/name)==digest,'source drift '+name)
    for name,digest in value['artifact_sha256'].items():
        target=(path.parent/name).resolve();base.require(target.is_relative_to(path.parent.resolve()) and base.sha(target)==digest,'artifact drift '+name)
    return value


def main():
    from check_process_stop import inventory,verify_previous
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('output','preflight','reference'):ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args();out=a.output.resolve();base.require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);started=time.perf_counter();sources=inventory();code=2
    report=dict(status='NOT_RUN_GATES_PENDING',source_sha256=sources,accepted_afsd_steps=0,new_afsd_steps=0,
        production_authorized=False,classification='LOCAL_STOP_REVERSE_SCOUT_NOT_AFSD',physical_geometry_unchanged=True)
    history=[];pa=None;state=None;current_context={};started_actual=False
    def timeout():
        base.write(out/'HARD_TIMEOUT.json',dict(status='WALL_BUDGET_STOP',accepted=False,restart_allowed=False));os._exit(124)
    timer=threading.Timer(55.,timeout);timer.daemon=True;timer.start()
    try:
        report['prior_before']=verify_previous()
        gate=validate_receipt(a.preflight,{'PASS_STOP_COMPONENT_PREFLIGHT_NOT_PROCESS'})
        base.require(gate['source_sha256']==sources,'stale preflight inventory')
        reference=validate_receipt(a.reference,{'PASS_EVENT_REFERENCE_8US','PASS_EVENT_REFERENCE_PREFIX_NEXT_BRANCH_BLOCKED'})
        base.require(reference['prefix_comparison']['passed'] and reference['event_agreement_passed'],'reference uncertainties failed')
        expected_config=hashlib.sha256(json.dumps(base.clean(inputs()),sort_keys=True).encode()).hexdigest()
        base.require(reference['configuration_sha256']==expected_config,'reference physical configuration mismatch')
        base.require('reference_history.npz' in reference['artifact_sha256'],'missing bound reference trajectory')
        base.require(reference['certified_common_end_time_s']>=6.5e-6,'insufficient postevent reference horizon')
        parent_root=ROOT/'outputs/process_curved_20261003/S0_initial'
        parent=validate_receipt(parent_root/'report.json')
        published=base.read(parent_root/'trajectory/accepted_history.json');last=published[-1]
        base.require(len(published)==13 and abs(last['time_s']-5e-6)<1e-19,'published5us parent only')
        path=parent_root/'trajectory'/last['checkpoint'];base.require(base.sha(path)==last['checkpoint_sha256'],'history checkpoint hash')
        raw,oldstate=restore(path);initial,p=base.setup()
        raw,oldstate=base.restore_checkpoint(raw,oldstate,p)
        base.require(oldstate['state_hash']==last['state_hash'],'published semantic hash')
        pa=raw;state=stage.migrate_published(pa,oldstate,p)
        integrals=base.read(path.parent/'integrals.json')
        for key in ('friction_raw','reservoir_heat','wall_work','impulse'):integrals[key]=np.asarray(integrals[key],float)
        # Actual previous same-problem15-evaluation step, inflated to24 and2x.
        measured=parent['run']['metrics'][-1]['step_seconds']
        forecast=2.*(12+3+2)*measured*(24./15.)+4.
        report.update(forecast_seconds=forecast,forecast_full_coupled_solves_max=17,
            measured_parent_step_seconds=measured,registered_inner_evaluations_max=24,
            parent_report_sha256=base.sha(parent_root/'report.json'),parent_checkpoint_sha256=base.sha(path),
            parent_state_hash=oldstate['state_hash'],preflight_sha256=base.sha(a.preflight),reference_sha256=base.sha(a.reference),
            configuration_sha256=reference['configuration_sha256'])
        base.require(forecast<=50.,'NOT_RUN_COST_EXCEEDS50S_FORECAST')
        base.write(out/'report.json',report)
        print(json.dumps(dict(forecast_seconds=forecast,maximum_local_end_s=8e-6,maximum_root_trials=12)),flush=True)
        trajectory=out/'trajectory';trajectory.mkdir();private=out/'private';private.mkdir()
        frames=[base.frame(pa,state,integrals)];calls=0;event_receipt=None;failures=[]
        base.require(inventory()==sources,'source drift before migration publication')
        migration=save_stage(trajectory/'stage_000_transition',pa,state,p,state['last_event'],integrals)
        migration['record_kind']='SAME_PHYSICAL_FIELDS_NUMERICAL_VARIANT_TRANSITION'
        history=[migration];base.write(trajectory/'accepted_history.json',history)
        base.write(trajectory/'frames.json',frames)
        def persist(candidate,next_state,diag,tag):
            nonlocal pa,state,integrals,frames,history
            base.checks(candidate,pa,initial,next_state,diag,p)
            base.require(next_state['time_s']<=reference['certified_common_end_time_s']+1e-19,'unreferenced candidate')
            total=copy.deepcopy(integrals);base.add_integrals(total,state,next_state,diag,p)
            base.require(inventory()==sources,'source drift before publication')
            item=save_stage(trajectory/tag,candidate,next_state,p,diag,total)
            updated=history+[item];base.write(trajectory/'accepted_history.json',updated)
            pa,state,integrals,history=candidate,next_state,total,updated
            frames.append(base.frame(pa,state,integrals));base.write(trajectory/'frames.json',frames)
        def trial(h,context):
            nonlocal calls
            base.require(calls<17,'registered full coupled solve count exhausted')
            calls+=1;name=f'trial_{calls:03d}';bc=boundary_for(p,context)
            try:
                result,rs,diag=stage.advance_private(pa,state,p,bc,h,base.slow_gap_certificate,context)
                checkpoint(private/(name+'.npz'),result,rs)
                with gzip.open(private/(name+'.json.gz'),'wt',encoding='utf-8') as stream:json.dump(base.clean(diag),stream,allow_nan=False)
                return result,rs,diag
            except Exception as error:
                receipt=dict(status='REJECTED',accepted=False,restart_allowed=False,dt_s=h,
                    parent_state_hash=state['state_hash'],error=str(error),diagnostic=getattr(error,'diagnostics',None))
                base.write(private/(name+'_REJECTED.json'),receipt);failures.append(receipt);raise
        started_actual=True
        horizon=min(8e-6,reference['certified_common_end_time_s'])
        targets=sorted(set([t for t in (6e-6,7e-6,8e-6) if t<=horizon]+[horizon]))
        report['targets_s']=targets
        for target in targets:
            duration=target-state['time_s'];bc=boundary_for(p,current_context)
            try:
                candidate,cs,diag=trial(duration,current_context)
                accepted,detail=certify_interval(pa,state,candidate,cs,diag,p,bc,base.slow_gap_certificate,current_context)
                persist(candidate,accepted,detail,f'stage_{len(history)+1:03d}')
                current_context={}
            except stage.CoupledStageError as original:
                if event_receipt is not None:raise
                pending=original.diagnostics.get('pending_unequilibrated_candidate')
                base.require(pending is not None,'failed stage has no dense-path predictor')
                context=discover_events(pa,state,p,boundary_for(p),duration,pending['acceleration_m_s2'])
                base.write(out/'event_predictor.json',dict(context=context,source='ACTUAL_FAILED_CANDIDATE_NOT_ORACLE'))
                # First bracket estimate uses this solver's actual failed path.
                # It is not an imposed event time and each duration fully resolves.
                upper=min(duration,1.2*min(row['predictor_dt_s'] for row in context.values()))
                located,ls,ld,root=locate_duration(lambda h:trial(h,context),pa,state,context,upper)
                strict={cid:dict(entry,endpoint='left_zero') for cid,entry in context.items()}
                accepted,detail=certify_interval(pa,state,located,ls,ld,p,boundary_for(p,strict),base.slow_gap_certificate,strict)
                detail['event_location']=root
                persist(located,accepted,detail,f'stage_{len(history)+1:03d}_event_left')
                right,outgoing,event=restart_sliding(pa,state,p,strict,lambda ctx:boundary_for(p,ctx),
                    np.repeat((pa.actuator==0)[:,None],3,axis=1).astype(int),
                    actual_wall_accelerations(p,state['time_s'],state['boundary_snapshot']['contacts']))
                base.require(inventory()==sources,'source drift before zero-duration restart publication')
                base.require(right['time_s']<=reference['certified_common_end_time_s']+1e-19,'unreferenced restart')
                restarted_integrals=copy.deepcopy(integrals)
                restarted_integrals['acceleration_restart_storage_change_J']=integrals.get('acceleration_restart_storage_change_J',0.)+event['signed_acceleration_storage']['signed_jump_J']
                item=save_stage(trajectory/f'stage_{len(history)+1:03d}_event_right',pa,right,p,event,restarted_integrals)
                item['record_kind']='ZERO_DURATION_FORCE_LIMIT_RESTART_NO_PHYSICAL_IMPULSE_OR_HEAT'
                history=history+[item];base.write(trajectory/'accepted_history.json',history);state=right;integrals=restarted_integrals
                frames.append(base.frame(pa,state,integrals));base.write(trajectory/'frames.json',frames)
                event_receipt=dict(location=root,restart=event);base.write(out/'event.json',event_receipt)
                current_context=outgoing
                remainder=target-state['time_s']
                base.require(remainder>0,'event outside attempted target interval')
                candidate,cs,diag=trial(remainder,current_context)
                accepted,detail=certify_interval(pa,state,candidate,cs,diag,p,boundary_for(p,current_context),base.slow_gap_certificate,current_context)
                persist(candidate,accepted,detail,f'stage_{len(history)+1:03d}');current_context={}
            print(json.dumps(dict(time_s=state['time_s'],published_records=len(history),private_calls=calls)),flush=True)
        base.require(event_receipt is not None,'no first fullzero/restart was resolved')
        report.update(status='PASS_LOCAL_LOCATED_STOP_REVERSE_NOT_AFSD',event=event_receipt,
            actual_private_coupled_calls=calls,published_records=len(history),time_s=state['time_s'],
            initial_failures_retained=len(failures),process_feasible=False);code=0
    except Exception as error:
        report.update(status='STOPPED_LOCAL_EVENT_SCOUT_NOT_PROCESS' if started_actual else 'NOT_RUN_GATE_REJECTED',
            error=repr(error),diagnostic=getattr(error,'diagnostics',None),published_records=len(history),
            time_s=None if state is None else state['time_s'],last_accepted_state_hash=None if state is None else state['state_hash'])
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={q.relative_to(out).as_posix():base.sha(q) for q in out.rglob('*') if q.is_file() and q.name!='report.json'})
    base.write(out/'report.json',report)
    print(json.dumps(base.clean({k:v for k,v in report.items() if k not in ('source_sha256','artifact_sha256','diagnostic','event')})),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

"""One real342 same-history response measurement, never a physical commit."""
from pathlib import Path
import argparse,gzip,json,os,sys,threading,time
import numpy as np
import probe_updated_stage as old
from check_process_bridge import verify_previous,inventory
from afsd_process_bridge import (direction,run_four,perturb,forecast,EPSILON,FIELDS,
    candidate_copy,candidate_hash,BridgeFailure)
from afsd_bridge_audit import audit
from afsd_backend.pysph.process_thermal_trial import ProcessThermalTrial
from afsd_backend.pysph.axial_actuator import control
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.output import restore
from afsd_backend.pysph.mechanics import positions
from afsd_backend.pysph.stage_iteration import residual_norms
from afsd_backend.pysph.application import memory_bytes

ROOT=old.ROOT;DT=5e-8;CLASSIFICATION='REAL_AFSD342_SAME_HISTORY_RESPONSE_MEASUREMENT_NOT_PROCESS'


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path,value):
    path=Path(path);pending=path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(old.clean(value),indent=2,allow_nan=False)+'\n',encoding='utf-8');pending.replace(path)


def compressed(path,value):
    path=Path(path);pending=path.with_suffix(path.suffix+'.pending')
    with gzip.open(pending,'wt',encoding='utf-8') as stream:json.dump(old.clean(value),stream,allow_nan=False,separators=(',',':'))
    pending.replace(path)


def hashes_match(values):
    for name,digest in values.items():
        path=(ROOT/name).resolve();old.require(path.is_relative_to(ROOT) and old.sha(path)==digest,'source/input drift '+name)


def receipt(path,status):
    path=Path(path).resolve();data=read(path);old.require(data['status']==status,'required receipt failed '+str(path))
    hashes_match(data.get('source_sha256',{}))
    for name,digest in data.get('artifact_sha256',{}).items():
        child=(path.parent/name).resolve();old.require(child.is_relative_to(path.parent) and old.sha(child)==digest,'artifact drift '+name)
    return data


def signature(e,pa,p):
    traces=[]
    for trace in e['boundary']['traces']:
        result=trace['result']
        # Event times are allowed to move; changes in the type/order/partition
        # are recorded as changed numerical branches, not smoothness proofs.
        event_sequence=[]
        for event in result['events']:
            contacts=[]
            for index in event['contacts']:
                old.require(type(index) in (int,np.int32,np.int64) and 0<=index<len(trace['contact_rows']),
                    'event contact index outside current face inventory')
                face=trace['contact_rows'][index]
                contacts.append((str(face['contact_id']),int(face['feature_id'])))
            event_sequence.append(dict(event=event['event'],contacts=contacts,keys=sorted(event)))
        traces.append(dict(particle_id=int(trace['particle_id']),
            faces=[(str(f['contact_id']),int(f['feature_id'])) for f in trace['contact_rows']],
            free_axes=np.asarray(trace['free_axes'],int),
            interval_kinds=[q['kind'] for q in result['intervals']],
            events=event_sequence))
    mat=e['material']
    return dict(traces=traces,plastic_active=np.asarray(mat['ep'])>pa.ep,
        rate_floor_active=np.asarray(mat['strain_rate_per_s'])<=p['rate_floor'],
        normal_before=e['normal_before']['active_contact_ids'],normal_after=e['normal_after']['active_contact_ids'])


def mechanism_coverage(summary):
    return (summary['positive_plastic_particle_count']>0 and summary['positive_plastic_work_J']>0
        and summary['contact_count']>0 and summary['multicontact_particle_count']>0
        and summary['actuator_active_particle_count']>0 and summary['nonzero_constraint_impulse_particle_count']>0
        and summary['constraint_work_J']!=0)


def response_arrays(e):
    b=e['boundary'];mat=e['material'];th=e['thermal'];traces=b['traces']
    face={};constraint=np.zeros_like(e['internal_force_N'])
    for trace in traces:
        result=trace['result'];idx=trace['particle_index'];constraint[idx]=result['constraint_impulse']
        for j,f in enumerate(trace['contact_rows']):
            face[str(f['contact_id'])]=dict(normal_impulse_Ns=np.asarray(result['normal_impulse'])[j],
                tangential_impulse_Ns=np.asarray(result['tangential_impulse'])[j],raw_Q_J=result['friction_heat_J'][j])
    return dict(faces=face,constraint_impulse_Ns=constraint,
        internal_force_N=e['internal_force_N'],plastic_work_J=mat['plastic_work_J'],
        ep=mat['ep'],return_dissipation_J=mat['return_dissipation_J'],
        friction_raw_J=b['friction_raw_J'],thermal_total_input_J=th['particle_energy_J']['total_heat_input_J'],
        actuator_work_J=b['ledger']['actuator_work_J'],wall_work_J=b['ledger']['wall_work_J'])


def numerical_digest(e):
    # No clocks occur in current evaluation output. A later model adding them
    # must explicitly amend scope rather than exclude arbitrary differing data.
    return old.fingerprint(e)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--preflight',type=Path,required=True)
    args=ap.parse_args();out=args.output.resolve();old.require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output')
    out.mkdir(parents=True);started=time.perf_counter()
    report=dict(status='NOT_RUN_GATES_PENDING',classification=CLASSIFICATION,accepted=False,restart_allowed=False,
        can_publish=False,production_authorized=False,new_afsd_steps=0,new_transport_events=0,
        actual_evaluations=0,actual_evaluations_attempted=0,measurements=[],
        physical_parameters_or_laws_changed=False,
        temporal_variant='EXISTING_NAMED_BE_THERMAL_MAP_VERSUS_HISTORICAL_EXPLICIT_MAP',
        internal_limit_seconds=110.,external_limit_seconds=120.)
    sources=None;inputs=None;pa=None;state=None;trial=None;before=None;state_before=None
    code=2;evidence={};first_digest=None;reentry=False;input_candidate=None
    def timeout():
        write(out/'HARD_TIMEOUT.json',dict(status='REJECTED_WALL_TIMEOUT',accepted=False,restart_allowed=False,
            completed_evaluations=report['actual_evaluations']));os._exit(124)
    timer=threading.Timer(110.,timeout);timer.daemon=True;timer.start()
    try:
        report['prior_before']=verify_previous();sources=inventory();report['source_sha256']=sources
        gate=receipt(args.preflight,'PASS_BRIDGE_PREFLIGHT_NOT_PROCESS')
        old.require(gate['source_sha256']==sources,'stale preflight');report['preflight_sha256']=old.sha(args.preflight)
        pa,state,p,inputs,parent_sources=old.inputs();hashes_match(parent_sources)
        before=particle_hash(pa);state_before=old.fingerprint(state)
        folder=ROOT/'outputs/nonlinear_continuation_20261003/continuation'
        archive=receipt(folder/'manifest.json','COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION')
        hashes_match(archive['input_sha256'])
        history=read(folder/'accepted_history.json');old.require(len(history)==4,'four recorded stages required')
        old.require(history[-1]['time_s']==archive['accepted_time_s']==850.393632910904e-9,'wrong latest physical time')
        old.require(history[-1]['snapshot']==archive['last_accepted_snapshot'],'wrong latest physical snapshot')
        for row in (history[0],history[-1]):
            path=folder/row['snapshot'];old.require(old.sha(path)==row['snapshot_sha256'],'published checkpoint drift')
            check_pa,check_state=restore(path)
            old.require(particle_hash(check_pa)==row['particle_state_sha256'],'published particle hash mismatch')
            old.require(check_state['parameters']==p and check_state['t']==row['time_s'],'published parameter/time mismatch')
        old.require(history[0]['dt_s']==DT and history[0]['time_s']==state['t']+DT,'wrong historical slice')
        input_candidate=folder/'stage_001/accepted_trial.json.gz'
        with gzip.open(input_candidate,'rt',encoding='utf-8') as stream:prior=json.load(stream)
        c=candidate_copy(prior['candidate']);old_eval=prior['evaluation']
        old.require(np.array_equal(old_eval['particle_id'],pa.particle_id) and
            old_eval['accepted_state_hash']==before and old_eval['time_before_s']==state['t'], 'candidate source history mismatch')
        stage_pa,stage_state=restore(folder/history[0]['snapshot'])
        for key,actual in [('position_m',positions(stage_pa)),('velocity_m_s',np.column_stack((stage_pa.u,stage_pa.v,stage_pa.w))),
                           ('density_kg_m3',stage_pa.rho),('temperature_K',stage_pa.temperature)]:
            old.require(np.array_equal(c[key],actual),'accepted candidate/checkpoint field mismatch '+key)
        inputs=dict(inputs)
        for name,digest in archive['artifact_sha256'].items():inputs[(folder/name).relative_to(ROOT).as_posix()]=digest
        inputs[(folder/'manifest.json').relative_to(ROOT).as_posix()]=old.sha(folder/'manifest.json')
        inputs.update(archive['source_sha256']);hashes_match(inputs)
        trial=ProcessThermalTrial(pa,p,state['t'])
        cmd=control(pa.z0,state['t']+.5*DT,0.,p['bar_length'],p['actuator_depth'],p['feed'],p['release_time'])
        D=direction(positions(pa),pa.m,pa.particle_id,cmd['free_axes'])
        config=dict(parameters=p,dt_s=DT,epsilon_m_s=EPSILON,classification=CLASSIFICATION,
            parent_configuration_hash=state['configuration_hash'],direction_definition='FROZEN_GEOMETRIC_FREE_AXIS_AFFINE_A',
            thermal_variant='UNCHANGED_EXISTING_PROCESS_THERMAL_TRIAL_BE',candidate_source=str(input_candidate.relative_to(ROOT)))
        report.update(input_sha256=inputs,configuration_sha256=old.fingerprint(config),parent_particle_hash=before,
            parent_state_hash=state_before,parent_time_s=state['t'],diagnostic_endpoint_s=state['t']+DT,
            latest_accepted_time_s=history[-1]['time_s'],particle_count=len(pa.x),
            historical_candidate_sha256=old.sha(input_candidate),baseline_candidate_hash=candidate_hash(c),
            initial_forecast_seconds=forecast(0.,4,[]),physical_dt_s=DT,epsilon_m_s=EPSILON,
            old_maximum_scaled_residual=history[0]['maximum_scaled_residual'],
            basis_not_a_zero_residual_reference=True)
        write(out/'configuration.json',config)
        compressed(out/'parent_state.json.gz',dict(state=state,properties={k:v.get_npy_array().copy() for k,v in pa.properties.items()}))
        write(out/'direction.json',dict(particle_id=pa.particle_id,D=D,free_axes=cmd['free_axes'],
            actual_plus_minus_delta={k:perturb(c,D,DT,1)[k]-perturb(c,D,DT,-1)[k] for k in FIELDS}))
        write(out/'last_accepted_reference.json',dict(latest_checkpoint=str((folder/history[-1]['snapshot']).relative_to(ROOT)),
            latest_time_s=history[-1]['time_s'],historical_parent=str((old.CASE/'plastic_accepted.npz').relative_to(ROOT)),
            historical_parent_time_s=state['t'],history_rewound=False,new_physical_state=False))
        write(out/'report.json',report)
        print(json.dumps(dict(forecast_seconds=report['initial_forecast_seconds'],remaining_forecast_seconds=forecast(time.perf_counter()-started,4,[]),
            maximum_evaluations=4,particle_count=len(pa.x),physical_steps=0,internal_seconds=110,external_seconds=120)),flush=True)
        def completed(index,label,candidate,e,seconds):
            nonlocal first_digest,reentry
            report['actual_evaluations']=index
            compressed(out/(label+'_candidate.json.gz'),dict(accepted=False,restart_allowed=False,can_publish=False,
                candidate=candidate,evaluation=e,evaluation_seconds=seconds))
            check=audit(pa,p,state['t'],DT,e);write(out/(label+'_audit.json'),check)
            old.require(particle_hash(pa)==before and trial.accepted_hash()==before and old.fingerprint(state)==state_before,'parent history mutation')
            sig=signature(e,pa,p);write(out/(label+'_branches.json'),sig)
            digest=numerical_digest(e)
            if label=='base':first_digest=digest
            if label=='base_repeat':reentry=digest==first_digest
            norms,scaled=residual_norms(candidate,e['map'],pa.rho)
            item=dict(label=label,evaluation_seconds=seconds,audit_status=check['status'],audit_passed=check['passed'],
                audit_summary=check['summary'],residual_norms=norms,maximum_scaled_residual=scaled,
                numerical_output_hash=digest,branch_hash=old.fingerprint(sig))
            report['measurements'].append(item);write(out/'report.json',report);print(json.dumps(old.clean(item)),flush=True)
            old.require(check['passed'],'independent mapped-account audit failed '+label)
            if label in ('base','base_repeat'):
                old.require(mechanism_coverage(check['summary']),'registered base does not cover plastic/multicontact/actuator mechanisms')
                report[label+'_mechanism_coverage']=True
            if label=='base_repeat':old.require(reentry,'same-history base reentry differs')
            if label in ('base','plus','minus'):
                evidence[label]=dict(residual={k:np.asarray(e['residual'][k]).copy() for k in FIELDS},
                    responses=response_arrays(e),branches=old.fingerprint(sig))
            return dict(audit_passed=True,numerical_output_hash=digest,branch_hash=item['branch_hash'])
        def evaluate(candidate,dt):
            report['actual_evaluations_attempted']+=1
            return trial.evaluate(candidate,dt)
        result=run_four(evaluate,c,D,DT,time.perf_counter,completed,started=started,limit_seconds=110.)
        report['orchestration']=result
        if result['completed']==4:
            a,b=evidence['plus'],evidence['minus'];same=a['branches']==b['branches']==evidence['base']['branches']
            quotient=dict(classification='ONE_EPSILON_DIRECTIONAL_QUOTIENT_NOT_CERTIFIED_JACOBIAN',
                epsilon_m_s=EPSILON,same_branch_signature=same,accepted=False,can_publish=False,
                residual_derivative={k:(a['residual'][k]-b['residual'][k])/(2*EPSILON) for k in FIELDS},
                residual_signed_difference={k:a['residual'][k]-b['residual'][k] for k in FIELDS},
                units=dict(position_m='s',velocity_m_s='dimensionless',mean_velocity_m_s='dimensionless',
                    density_kg_m3='kg*s/m^4',temperature_K='K*s/m'))
            for key in a['responses']:
                if key=='faces':continue
                quotient[key+'_derivative_per_m_s']=(np.asarray(a['responses'][key])-np.asarray(b['responses'][key]))/(2*EPSILON)
            common=set(a['responses']['faces']) & set(b['responses']['faces'])
            quotient['face_id_sets_equal']=set(a['responses']['faces'])==set(b['responses']['faces'])
            quotient['face_derivative_per_m_s']={cid:{key:(np.asarray(a['responses']['faces'][cid][key])-
                np.asarray(b['responses']['faces'][cid][key]))/(2*EPSILON) for key in a['responses']['faces'][cid]} for cid in sorted(common)}
            write(out/'directional_response.json',quotient)
            report.update(status='PASS_RESPONSE_MEASUREMENT_NOT_PROCESS' if same else 'NONSMOOTH_NOT_CERTIFIED',
                same_branch_signature=same,base_reentry_exact=reentry,actual_map_response_measured=True);code=0
        else:report.update(status='PARTIAL_BUDGET_NOT_RUN_REMAINDER',not_run=result['not_run'])
    except Exception as error:
        context=getattr(error,'original',error)
        compressed(out/'REJECTED_DIAGNOSTIC.json.gz',dict(status='REJECTED',accepted=False,restart_allowed=False,
            error=repr(error),candidate=getattr(error,'candidate',None),completed=getattr(error,'records',None),
            failed_call=getattr(error,'label',None),failure_context=getattr(context,'contact_context',getattr(context,'feature_context',None))))
        report.update(status='FAIL_BRIDGE_DIAGNOSTIC' if report['actual_evaluations_attempted'] else 'NOT_RUN_GATE_REJECTED',error=repr(error))
    try:
        if sources is not None:old.require(inventory()==sources,'source drift after measurement')
        if inputs is not None:hashes_match(inputs)
        if pa is not None:
            report['parent_unchanged']=particle_hash(pa)==before and old.fingerprint(state)==state_before
            old.require(report['parent_unchanged'],'parent changed')
        report['prior_after']=verify_previous()
    except Exception as error:report.update(status='FAIL_BRIDGE_INTEGRITY',integrity_error=repr(error));code=2
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={q.relative_to(out).as_posix():old.sha(q) for q in out.rglob('*') if q.is_file() and q.name!='report.json'})
    write(out/'report.json',report)
    print(json.dumps(old.clean({k:v for k,v in report.items() if k not in ('source_sha256','input_sha256','artifact_sha256','measurements','orchestration')})),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

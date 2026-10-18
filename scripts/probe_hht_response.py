"""One saved PRIVATE numerical search experiment; no physical publication."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import os
import threading
import time
import numpy as np
import probe_hht_curved as base
from probe_hht_stop import validate_receipt
from hht_stop_boundary import boundary_for
from hht_slow_inputs import inputs
from afsd_backend.pysph import hht_stop_stage as original
from afsd_backend.pysph.hht_response_stage import advance_private_response
from afsd_backend.pysph.contact_response_preconditioner import STRATEGY
from afsd_backend.pysph.output import checkpoint,restore
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.application import memory_bytes

ROOT=base.ROOT


def object_hash(value):
    return hashlib.sha256(json.dumps(base.clean(value),sort_keys=True,allow_nan=False).encode()).hexdigest()


def main():
    from check_process_response import inventory,verify_previous
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('output','preflight','limits'):ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args();out=args.output.resolve()
    base.require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);started=time.perf_counter();sources=inventory()
    report=dict(status='NOT_RUN_GATES_PENDING',source_sha256=sources,accepted=False,
        restart_allowed=False,can_publish=False,production_authorized=False,new_afsd_steps=0,
        new_published_local_steps=0,actual_private_coupled_calls=0,
        classification=STRATEGY,continuous_cross_zero_jacobian_certified=False)
    code=2;pa=None;state=None;parent_particle=None;parent_object=None
    def timeout():
        base.write(out/'HARD_TIMEOUT.json',dict(status='WALL_BUDGET_STOP',accepted=False,restart_allowed=False))
        os._exit(124)
    timer=threading.Timer(55.,timeout);timer.daemon=True;timer.start()
    try:
        report['prior_before']=verify_previous()
        preflight=validate_receipt(args.preflight,{'PASS_RESPONSE_COMPONENT_PREFLIGHT_NOT_PROCESS'})
        base.require(preflight['source_sha256']==sources,'stale response preflight')
        limits=base.read(args.limits)
        base.require(limits['status']=='PASS_MATHEMATICAL_NONDIFFERENTIABILITY_EVIDENCE_NOT_SOLVER_VALIDATION',
            'missing independent continuous cusp evidence')
        for name,digest in limits['source_sha256'].items():
            base.require(base.sha(ROOT/name)==digest,'limits source drift '+name)
        base.require(limits['checks_run']==130 and all(row['passed'] for row in limits['checks']),
            'independent mathematical checks incomplete')
        parent_root=ROOT/'outputs/process_stop_20261003/scout_initial'
        prior=validate_receipt(parent_root/'report.json',{'STOPPED_LOCAL_EVENT_SCOUT_NOT_PROCESS'})
        history=base.read(parent_root/'trajectory/accepted_history.json');last=history[-1]
        base.require(len(history)==2 and last['time_s']==6e-6,'only published6us parent permitted')
        parent_path=parent_root/'trajectory'/last['checkpoint']
        base.require(base.sha(parent_path)==last['checkpoint_sha256'],'parent file hash mismatch')
        pa,stored=restore(parent_path);initial,p=base.setup()
        pa,state=original.restore_checkpoint(pa,stored,p)
        base.require(state['state_hash']==last['state_hash'],'parent semantic hash mismatch')
        parent_particle=particle_hash(pa);parent_object=object_hash(state)
        old_trial_path=parent_root/'private/trial_005_REJECTED.json';old_trial=base.read(old_trial_path)
        base.require(old_trial['parent_state_hash']==state['state_hash'],'saved failure parent differs')
        dt=old_trial['dt_s'];base.require(dt==2.7254762326048973e-7,'saved interval changed')
        context_path=parent_root/'event_predictor.json';context=base.read(context_path)['context']
        base.require(context and all(row['endpoint']=='private_continuation' for row in context.values()),
            'saved private context missing')
        configuration_hash=object_hash(inputs())
        base.require(prior['configuration_sha256']==configuration_hash,'physical configuration changed')
        elapsed=old_trial['diagnostic']['elapsed_s'];evaluations=old_trial['diagnostic']['evaluations']
        forecast=2.*elapsed*(24./evaluations)+5.
        report.update(preflight_sha256=base.sha(args.preflight),limits_receipt_sha256=base.sha(args.limits),
            configuration_sha256=configuration_hash,parent_report_sha256=base.sha(parent_root/'report.json'),
            parent_checkpoint=str(parent_path.relative_to(ROOT)),parent_checkpoint_sha256=base.sha(parent_path),
            parent_state_hash=state['state_hash'],parent_particle_hash=parent_particle,
            parent_state_object_hash=parent_object,parent_time_s=state['time_s'],
            saved_failure_sha256=base.sha(old_trial_path),event_context_sha256=base.sha(context_path),
            dt_s=dt,attempted_time_s=state['time_s']+dt,
            forecast_seconds=forecast,old_measured_seconds=elapsed,old_evaluations=evaluations,
            maximum_evaluations=24,old_raw_residuals=old_trial['diagnostic']['last_metrics'])
        base.write(out/'last_accepted_state.json',dict(checkpoint=report['parent_checkpoint'],
            checkpoint_sha256=report['parent_checkpoint_sha256'],state_hash=state['state_hash'],
            time_s=state['time_s'],new_physical_state=False))
        base.write(out/'configuration.json',dict(parameters=p,canonical_inputs=inputs(),
            event_context=context,dt_s=dt,classification=STRATEGY))
        base.require(forecast<=50.,'NOT_RUN_COST_EXCEEDS50S_FORECAST')
        base.require(inventory()==sources,'source drift before actual call')
        base.write(out/'report.json',report)
        print(json.dumps(dict(forecast_seconds=forecast,maximum_coupled_calls=1,
            maximum_evaluations=24,physical_steps=0,private_dt_s=dt)),flush=True)
        report['actual_private_coupled_calls']=1
        candidate,cstate,diag=advance_private_response(pa,state,p,boundary_for(p,context),dt,
            base.slow_gap_certificate,context)
        base.require(not cstate['accepted'] and not cstate['restart_allowed'] and not diag['can_publish'],
            'private output has publication authority')
        target=out/'PRIVATE_NONRESTARTABLE.npz';checkpoint(target,candidate,cstate)
        from pysph.solver.utils import load
        raw=load(str(target))
        base.require(particle_hash(raw['arrays']['solid'])==particle_hash(candidate), 'debug checkpoint field readback')
        base.require(object_hash(raw['solver_data'])==object_hash(cstate),'debug checkpoint history readback')
        try:restore(target)
        except ValueError:report['restart_rejection_verified']=True
        else:raise RuntimeError('private output incorrectly restartable')
        with gzip.open(out/'PRIVATE_DIAGNOSTIC.json.gz','wt',encoding='utf-8') as stream:
            json.dump(base.clean(diag),stream,allow_nan=False)
        report.update(status='PASS_PRIVATE_NUMERICAL_ROOT_NOT_EVENT_OR_PROCESS',
            evaluations=diag['evaluations'],actual_solve_seconds=diag['elapsed_s'],
            raw_residuals=diag['residuals'],response_updates=len(diag['contact_response_history']),
            raw_total_friction_J=float(np.sum(diag['contact']['friction_raw_J'])),
            candidate_state_hash=cstate['state_hash'],candidate_particle_hash=particle_hash(candidate),
            old_velocity_residual_over_new=(report['old_raw_residuals']['velocity_residual_m_s']/
                diag['residuals']['velocity_residual_m_s']))
        code=0
    except Exception as error:
        detail=getattr(error,'diagnostics',None)
        if report['actual_private_coupled_calls']:
            base.write(out/'REJECTED_CANDIDATE.json',dict(status='REJECTED',accepted=False,
                restart_allowed=False,error=repr(error),diagnostic=detail))
            report.update(status='FAIL_PRIVATE_NUMERICAL_DIAGNOSTIC_NOT_PROCESS',error=repr(error))
            if detail:
                report.update(evaluations=detail.get('evaluations'),actual_solve_seconds=detail.get('elapsed_s'),
                    raw_residuals=detail.get('last_metrics'),
                    response_updates=len(detail.get('contact_response_history',[])))
        else:report.update(status='NOT_RUN_GATE_REJECTED',error=repr(error))
    finally:timer.cancel()
    try:
        if pa is not None and parent_particle is not None:
            report['parent_particle_unchanged']=particle_hash(pa)==parent_particle
            report['parent_state_unchanged']=object_hash(state)==parent_object
            base.require(report['parent_particle_unchanged'] and report['parent_state_unchanged'],
                'accepted parent was mutated')
        base.require(inventory()==sources,'source drift during diagnostic')
        report['prior_after']=verify_previous()
    except Exception as error:
        report.update(status='FAIL_INTEGRITY',integrity_error=repr(error));code=2
    report.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={q.relative_to(out).as_posix():base.sha(q) for q in out.rglob('*')
            if q.is_file() and q.name!='report.json'})
    base.write(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','artifact_sha256')}),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

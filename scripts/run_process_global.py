"""One preregistered 500 ns joint root; never deposition/production validation."""
from pathlib import Path
import argparse
import json
import os
import sys
import threading
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from probe_updated_stage import sha,require,compressed
from run_coupled_afsd import write_json,clean,fingerprint
from run_nonlinear_continuation import accumulate
from afsd_backend.pysph.continuation_state import restore as restore_parent
from afsd_backend.pysph.process_thermal_trial import ProcessThermalTrial,acceptance
from afsd_backend.pysph.stage_continuity_preconditioner import ContinuityStagePreconditioner
from afsd_backend.pysph.stage_structural_root import solve_stage
from afsd_backend.pysph.stage_iteration import initial_candidate,StageFailure,residual_norms
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.process_global_state import CLASSIFICATION,commit,state,save,restore
from afsd_backend.pysph.transport_observer import observe
from afsd_backend.pysph.application import memory_bytes

CASE=ROOT/'outputs/nonlinear_continuation_20261003/continuation'
SETTINGS=dict(dt_s=5e-7,maximum_steps=1,maximum_evaluations=32,beta=.75,history_depth=6,
    qr_condition_limit=1e10,growth_limit=2.,forecast_multiplier=1.5,overhead_s=30.,
    launch_limit_s=450.,inner_limit_s=540.,thermal_variant='BACKWARD_EULER',
    preconditioner='FROZEN_CONTINUITY_SCHUR_ELASTIC_NOT_FULL_TANGENT')
MODULES=['structural_mobility','implicit_thermal','stage_preconditioner','stage_continuity_preconditioner',
         'process_thermal_trial','stage_structural_root','process_global_state']
SOURCES=['scripts/run_process_global.py','scripts/check_process_global.py',
    'docs/PROCESS_GLOBAL_ROOT_CONTRACT.md','docs/PROCESS_GLOBAL_THERMAL_VARIANT.md',
    'docs/PROCESS_GLOBAL_PRECONDITIONER.md']+['src/afsd_backend/pysph/'+m+'.py' for m in MODULES]


def setup():
    manifest=json.loads((CASE/'manifest.json').read_text(encoding='utf-8'))
    require(manifest['status']=='COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION','published parent required')
    watched={(CASE/'manifest.json').relative_to(ROOT).as_posix():sha(CASE/'manifest.json')}
    for rel,h in manifest['artifact_sha256'].items():
        path=(CASE/rel).resolve();require(path.is_relative_to(CASE) and sha(path)==h,'parent artifact drift '+rel)
        watched[path.relative_to(ROOT).as_posix()]=h
    for rel,h in manifest['source_sha256'].items():require(sha(ROOT/rel)==h,'parent source drift '+rel)
    history=json.loads((CASE/'accepted_history.json').read_text(encoding='utf-8'));entry=history[-1]
    require(len(history)==4 and entry['snapshot']==manifest['last_accepted_snapshot'],'parent history mismatch')
    p=json.loads((CASE/'resolved_config.json').read_text(encoding='utf-8'))['parameters']
    path=CASE/entry['snapshot']
    pa,old=restore_parent(path,p,manifest['source_sha256'],manifest['configuration_sha256'],entry['snapshot_sha256'])
    require(old['t']==entry['time_s']==manifest['accepted_time_s'] and len(pa.x)==342,'parent time/size mismatch')
    sources=dict(manifest['source_sha256'],**{s:sha(ROOT/s) for s in SOURCES})
    config=fingerprint(dict(classification=CLASSIFICATION,parameters=p,settings=SETTINGS))
    parent=dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),source_time_s=old['t'],source_stage=4)
    return pa,old,p,watched,sources,config,parent


def preflight_receipt(path,sources):
    path=Path(path).resolve();require(path.is_relative_to(ROOT),'local preflight required')
    r=json.loads(path.read_text(encoding='utf-8'))
    require(r['status']=='PASS_GLOBAL_PREFLIGHT_NOT_PROCESS_VALIDATION','global preflight must pass')
    for rel,h in r['source_sha256'].items():require(sha(ROOT/rel)==h,'stale preflight '+rel)
    for rel in SOURCES:require(r['source_sha256'].get(rel)==sources[rel],'source not tested '+rel)
    for rel,h in r['artifact_sha256'].items():
        item=(path.parent/rel).resolve();require(item.is_relative_to(path.parent) and sha(item)==h,'preflight artifact drift')
    return dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),tests_run=r['tests_run'])


def budget(path,sources,config,preflight,watched,parent):
    path=Path(path).resolve();require(path.is_relative_to(ROOT),'local cost probe required')
    r=json.loads(path.read_text(encoding='utf-8'))
    require(r['status']=='COMPLETE_UNCOMMITTED_COST_PROBE' and r['classification']==CLASSIFICATION
        and r['settings']==SETTINGS and r['accepted_steps']==0 and r['production_authorized'] is False,'wrong probe')
    require(r['source_sha256']==sources and r['configuration_sha256']==config and r['preflight']==preflight,'stale probe')
    require(r['input_sha256']==watched and r['parent_reference']==parent,'cost probe initial history changed')
    for rel,h in r['artifact_sha256'].items():
        item=(path.parent/rel).resolve();require(item.is_relative_to(path.parent) and sha(item)==h,'probe artifact drift')
    build=r['preconditioner_build_wall_seconds'];cost=r['trial_and_proposal_wall_seconds']
    require(np.isfinite([build,cost]).all() and build>0 and cost>0,'invalid measured costs')
    predicted=SETTINGS['forecast_multiplier']*(build+SETTINGS['maximum_evaluations']*cost)+SETTINGS['overhead_s']
    return dict(probe=path.relative_to(ROOT).as_posix(),probe_sha256=sha(path),forecast_wall_seconds=predicted,
        allowed=predicted<=SETTINGS['launch_limit_s'],maximum_steps=1,maximum_evaluations=32,
        requested_physical_increment_s=SETTINGS['dt_s'])


def run(directory,preflight_path,probe=False,probe_manifest=None):
    directory=Path(directory).resolve()
    require(directory.is_relative_to(ROOT) and not directory.exists(),'fresh project-local output required')
    directory.mkdir(parents=True);started=time.perf_counter();candidate=evaluation=None;trace=[];published=False
    manifest=dict(classification=CLASSIFICATION,status='RUNNING_UNCOMMITTED',settings=SETTINGS,
        mode='COST_PROBE' if probe else 'ONE_LARGER_STEP_ROOT',accepted_steps=0,accepted_progress=False,
        production_authorized=False,process_accuracy='NOT_VALIDATED',backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5')
    def deadline():
        try:write_json(directory/'HARD_TIMEOUT.json',dict(status='WALL_BUDGET_STOP',accepted=False,restart_allowed=False))
        finally:os._exit(124)
    timer=threading.Timer(80. if probe else SETTINGS['inner_limit_s'],deadline);timer.daemon=True;timer.start()
    try:
        pa,old,p,watched,sources,config,parent=setup();before=particle_hash(pa)
        preflight=preflight_receipt(preflight_path,sources)
        manifest.update(source_sha256=sources,input_sha256=watched,configuration_sha256=config,preflight=preflight,
            parent_reference=parent,source_time_s=old['t'],accepted_time_s=old['t'],particle_count=len(pa.x),
            requested_end_s=old['t']+SETTINGS['dt_s'])
        write_json(directory/'resolved_config.json',dict(classification=CLASSIFICATION,parameters=p,settings=SETTINGS,
            physical_sources='UNCHANGED_PARENT_PARAMETERS',configuration_sha256=config))
        write_json(directory/'last_accepted_reference.json',parent);write_json(directory/'accepted_history.json',[])
        write_json(directory/'manifest.json',manifest)
        if not probe:
            require(probe_manifest is not None,'measured cost required')
            receipt=budget(probe_manifest,sources,config,preflight,watched,parent);write_json(directory/'budget.json',receipt)
            manifest['budget']=receipt;require(receipt['allowed'],'NOT_RUN_COST_EXCEEDS_BUDGET')
        tick=time.perf_counter();pre=ContinuityStagePreconditioner(pa,p,SETTINGS['dt_s'],old['t'])
        manifest['preconditioner_build_wall_seconds']=time.perf_counter()-tick
        trial=ProcessThermalTrial(pa,p,old['t'])
        if probe:
            candidate=initial_candidate(trial,SETTINGS['dt_s']);tick=time.perf_counter()
            evaluation=trial.evaluate(candidate,SETTINGS['dt_s'])
            proposal=pre.apply(candidate,evaluation['map'])
            compressed(directory/'probe_trial.json.gz',dict(candidate=candidate,evaluation=evaluation,numerical_proposal_only=proposal))
            norms,scaled=residual_norms(candidate,evaluation['map'],pa.rho)
            row=dict(evaluation=1,residual_norms=norms,maximum_scaled_residual=scaled,
                candidate_disposition='UNCOMMITTED_COST_PROBE_NOT_ACCEPTED')
            compressed(directory/'probe_candidate.json.gz',dict(candidate=candidate,row=row))
            write_json(directory/'probe_iterations.json',[row])
            # Include same-size complete journal, candidate and trace IO as a real evaluation.
            manifest.update(status='COMPLETE_UNCOMMITTED_COST_PROBE',trial_and_proposal_wall_seconds=time.perf_counter()-tick)
            require(particle_hash(pa)==before,'parent state mutated')
            for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'source/input drift '+rel)
        else:
            def journal(row,c,e):
                nonlocal candidate,evaluation
                candidate=c;evaluation=e;trace.append(row)
                compressed(directory/'latest_trial.json.gz',dict(candidate=c,evaluation=e))
                compressed(directory/f'candidate_{len(trace):03d}.json.gz',dict(classification='UNACCEPTED_NUMERICAL_CANDIDATE',candidate=c,row=row))
                write_json(directory/'iterations.json',trace)
                print(json.dumps(clean(row),allow_nan=False),flush=True)
            solved=solve_stage(trial,pre,SETTINGS['dt_s'],SETTINGS['maximum_evaluations'],
                max(1.,SETTINGS['inner_limit_s']-(time.perf_counter()-started)),journal)
            candidate,evaluation=solved['candidate'],solved['evaluation']
            checks=acceptance(pa,p,trial,candidate,evaluation)
            after=commit(pa,p,candidate,evaluation,checks)
            transport=observe(pa,after,old['t'],evaluation['time_after_s'],p,old['transport_cumulative'])
            ledger=dict(boundary={},thermal={},energy_diagnostic_residual_J=0.)
            accumulate(ledger,evaluation,checks)
            s=state(after,p,evaluation['time_after_s'],1,sources,config,parent,transport['cumulative'],ledger)
            for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'pre-publication drift '+rel)
            save(directory,after,p,s)
            loaded,_=restore(directory/'accepted.npz',p,sources,config,sha(directory/'accepted.npz'))
            require(particle_hash(loaded)==particle_hash(after),'accepted roundtrip failed')
            write_json(directory/'acceptance.json',checks);write_json(directory/'transport.json',transport)
            write_json(directory/'tool_pose.json',s['tool_pose'])
            compressed(directory/'accepted_trial.json.gz',dict(candidate=candidate,evaluation=evaluation))
            entry=dict(stage_index=1,time_s=s['t'],dt_s=SETTINGS['dt_s'],evaluations=len(trace),snapshot='accepted.npz',
                snapshot_sha256=sha(directory/'accepted.npz'),particle_state_sha256=particle_hash(after),
                maximum_scaled_residual=checks['maximum_scaled_residual'],transport=transport['crossing'],
                checkpoint_readback='IDENTICAL_FULL_HISTORY; RECOMPUTED_RESTART_EQUIVALENCE_NOT_RUN')
            require(particle_hash(pa)==before,'parent state mutated')
            for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'source/input drift before publication '+rel)
            write_json(directory/'accepted_history.json',[entry])
            published=True
            manifest.update(status='COMPLETE_ONE_NUMERICAL_ROOT_NOT_PROCESS_VALIDATION',accepted_steps=1,
                accepted_progress=True,accepted_time_s=s['t'],last_accepted_snapshot='accepted.npz',
                total_evaluations=len(trace),transport_status=transport['crossing']['status'])
        code=0
    except Exception as error:
        if isinstance(error,StageFailure):candidate,evaluation,trace=error.candidate,error.evaluation,error.trace
        compressed(directory/'REJECTED_candidate.json.gz',dict(classification='REJECTED',accepted=False,restart_allowed=False,
            error=repr(error),candidate=candidate,last_completed_evaluation=evaluation,iterations=trace,
            failure_context=getattr(error,'original_context',getattr(error,'contact_context',getattr(error,'feature_context',None)))))
        manifest.update(status='FAILED_AFTER_PUBLISHED_ACCEPTED_STATE' if published else 'BLOCKED_NO_NEW_ACCEPTED_STATE',
            error=repr(error),total_evaluations=len(trace));code=2
    finally:timer.cancel()
    manifest.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={f.relative_to(directory).as_posix():sha(f) for f in directory.rglob('*')
            if f.is_file() and f!=directory/'manifest.json'})
    write_json(directory/'manifest.json',manifest,rejected=code!=0)
    print(json.dumps(clean({k:v for k,v in manifest.items() if k not in ('source_sha256','input_sha256','artifact_sha256')},rejected=code!=0)),flush=True)
    return code


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--preflight-manifest',type=Path,required=True)
    mode=ap.add_mutually_exclusive_group(required=True);mode.add_argument('--probe',action='store_true')
    mode.add_argument('--probe-manifest',type=Path)
    a=ap.parse_args();return run(a.output,a.preflight_manifest,a.probe,a.probe_manifest)


if __name__=='__main__':raise SystemExit(main())

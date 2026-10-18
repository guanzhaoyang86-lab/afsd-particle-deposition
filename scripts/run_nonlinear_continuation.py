"""Bounded same-physics global roots, transactional continuation and transport."""
from pathlib import Path
import argparse
import json
import os
import sys
import threading
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from probe_updated_stage import inputs, acceptance, compressed, sha, require
from run_coupled_afsd import clean, write_json, fingerprint
from afsd_backend.pysph.updated_trial import UpdatedTrial,particle_hash
from afsd_backend.pysph.stage_iteration import initial_candidate,StageFailure
from afsd_backend.pysph.stage_anderson import solve_stage
from afsd_backend.pysph.continuation_state import CLASSIFICATION,commit,state,save,restore
from afsd_backend.pysph.transport_observer import observe
from afsd_backend.pysph.application import memory_bytes

SETTINGS=dict(dt_s=5e-8,maximum_steps=4,maximum_evaluations_per_step=32,total_evaluation_budget=72,
    beta=.75,history_depth=6,qr_condition_limit=1e10,growth_limit=2.,
    forecast_multiplier=1.5,overhead_s=30.,launch_limit_s=450.,inner_limit_s=540.)
SOURCES=['scripts/run_nonlinear_continuation.py','src/afsd_backend/pysph/stage_anderson.py',
    'src/afsd_backend/pysph/continuation_state.py','src/afsd_backend/pysph/transport_observer.py',
    'docs/NONLINEAR_CONTINUATION_CONTRACT.md']


def setup():
    pa,old,p,watched,sources=inputs()
    frozen=json.loads((ROOT/'docs/UPDATED_TRIAL_FREEZE.json').read_text(encoding='utf-8'))['files_sha256']
    for rel,h in sources.items():
        if not rel.startswith('third_party/'):
            require(frozen.get(rel)==h,'prior updated-source freeze mismatch '+rel)
    sources.update({s:sha(ROOT/s) for s in SOURCES})
    watched['docs/UPDATED_TRIAL_FREEZE.json']=sha(ROOT/'docs/UPDATED_TRIAL_FREEZE.json')
    config=fingerprint(dict(classification=CLASSIFICATION,parameters=p,settings=SETTINGS))
    parent=dict(path='outputs/coupled_afsd_20261002/short_07us/plastic_accepted.npz',
                sha256=sha(ROOT/'outputs/coupled_afsd_20261002/short_07us/plastic_accepted.npz'),
                source_time_s=old['t'],source_step=old['count'])
    return pa,old,p,watched,sources,config,parent


def budget(path,sources,config):
    path=Path(path).resolve();require(path.is_relative_to(ROOT),'project-local cost probe required')
    report=json.loads(path.read_text(encoding='utf-8'))
    require(report['status']=='COMPLETE_UNCOMMITTED_COST_PROBE' and report['mode']=='COST_PROBE' and
        report['classification']==CLASSIFICATION and report['production_authorized'] is False and
        report['accepted_steps']==0 and report['settings']==SETTINGS,'wrong cost-probe scope')
    require(report['source_sha256']==sources and report['configuration_sha256']==config,'cost source/config drift')
    require({'resolved_config.json','last_accepted_reference.json','probe_trial.json.gz'} <= set(report['artifact_sha256']),
        'incomplete probe evidence')
    for name,h in report['artifact_sha256'].items():
        item=(path.parent/name).resolve()
        require(item.is_relative_to(path.parent) and sha(item)==h,'probe artifact drift')
    elapsed=float(report['complete_trial_wall_seconds']);require(np.isfinite(elapsed) and elapsed>0,'invalid measured cost')
    affordable=int(np.floor((SETTINGS['launch_limit_s']-SETTINGS['overhead_s'])/(SETTINGS['forecast_multiplier']*elapsed)))
    evaluation_budget=min(SETTINGS['total_evaluation_budget'],affordable)
    predicted=SETTINGS['forecast_multiplier']*elapsed*evaluation_budget+SETTINGS['overhead_s']
    return dict(probe=path.relative_to(ROOT).as_posix(),probe_sha256=sha(path),forecast_wall_seconds=predicted,
        allowed=evaluation_budget>=SETTINGS['maximum_evaluations_per_step'] and predicted<=SETTINGS['launch_limit_s'],
        evaluation_budget=evaluation_budget,maximum_physical_steps=SETTINGS['maximum_steps'],
        maximum_physical_increment_s=SETTINGS['dt_s']*SETTINGS['maximum_steps'])


def accumulate(ledger,evaluation,checks):
    boundary=evaluation['boundary']['ledger']
    nonincrements={'contacts','max_penetration_m','normal_potential_J','minimum_raw_friction_J'}
    counts={'contact_subintervals','contact_rejected_subintervals','contact_zero_events'}
    for key,value in boundary.items():
        if key in nonincrements:continue
        require(key.endswith(('_J','_Ns','_Nms')) or key in counts,'unclassified boundary ledger field '+key)
        a=np.asarray(value);ledger['boundary'][key]=clean(np.asarray(ledger['boundary'].get(key,np.zeros_like(a)))+a)
    for key,value in evaluation['thermal']['ledger'].items():
        require(key.endswith('_J'),'nonincremental thermal ledger field '+key)
        ledger['thermal'][key]=float(ledger['thermal'].get(key,0.))+float(value)
    ledger['boundary_state_latest']={key:boundary[key] for key in nonincrements}
    ledger['maximum_penetration_m']=max(ledger.get('maximum_penetration_m',0.),boundary['max_penetration_m'])
    q=boundary['minimum_raw_friction_J']
    if q is not None:ledger['minimum_raw_friction_J']=q if ledger.get('minimum_raw_friction_J') is None else min(ledger['minimum_raw_friction_J'],q)
    ledger['energy_diagnostic_residual_J']+=checks['total_energy_diagnostic_residual_J']


def run(directory,probe=False,probe_manifest=None):
    directory=Path(directory).resolve()
    require(directory.is_relative_to(ROOT) and not directory.exists(),'fresh project-local output directory required')
    directory.mkdir(parents=True);started=time.perf_counter()
    def deadline():
        try:
            write_json(directory/'HARD_TIMEOUT.json',dict(status='WALL_BUDGET_STOP_UNCOMMITTED_TRIAL_REJECTED',
                restart='Only a published accepted_history entry is eligible for a future same-source restart.',production_authorized=False))
        finally:os._exit(124)
    timer=threading.Timer(80. if probe else SETTINGS['inner_limit_s'],deadline);timer.daemon=True;timer.start()
    manifest=dict(classification=CLASSIFICATION,status='RUNNING_UNCOMMITTED',mode='COST_PROBE' if probe else 'BOUNDED_CONTINUATION',
        production_authorized=False,settings=SETTINGS,accepted_steps=0,accepted_progress=False,
        backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5')
    history=[];iterations=[];evaluations=0;candidate=None;evaluation=None;stage_dir=directory
    try:
        pa,old,p,watched,sources,config,parent=setup();time_s=old['t'];initial_hash=particle_hash(pa)
        manifest.update(source_sha256=sources,input_sha256=watched,configuration_sha256=config,
            parent_reference=parent,source_time_s=time_s,accepted_time_s=time_s,particle_count=len(pa.x),
            requested_end_s=time_s+SETTINGS['dt_s']*SETTINGS['maximum_steps'])
        write_json(directory/'resolved_config.json',dict(classification=CLASSIFICATION,parameters=p,settings=SETTINGS,
            physical_sources='UNCHANGED_FROZEN_PARENT_RESOLVED_CONFIG',configuration_sha256=config))
        write_json(directory/'last_accepted_reference.json',parent);write_json(directory/'accepted_history.json',history)
        write_json(directory/'manifest.json',manifest)
        if probe:
            trial=UpdatedTrial(pa,p,time_s);candidate=initial_candidate(trial,SETTINGS['dt_s']);tick=time.perf_counter()
            evaluation=trial.evaluate(candidate,SETTINGS['dt_s'])
            compressed(directory/'probe_trial.json.gz',dict(candidate=candidate,evaluation=evaluation))
            manifest.update(status='COMPLETE_UNCOMMITTED_COST_PROBE',complete_trial_wall_seconds=time.perf_counter()-tick)
            require(particle_hash(pa)==initial_hash,'probe mutated accepted history')
        else:
            require(probe_manifest is not None,'same-source measured cost required')
            receipt=budget(probe_manifest,sources,config);write_json(directory/'budget.json',receipt)
            manifest['budget']=receipt;require(receipt['allowed'],'NOT_RUN_COST_EXCEEDS_BUDGET')
            cumulative=None;ledger=dict(boundary={},thermal={},energy_diagnostic_residual_J=0.)
            for step in range(1,SETTINGS['maximum_steps']+1):
                remaining=receipt['evaluation_budget']-evaluations
                require(remaining>0,'total nonlinear evaluation budget exhausted')
                stage_dir=directory/f'stage_{step:03d}';stage_dir.mkdir()
                trial=UpdatedTrial(pa,p,time_s);iterations=[];candidate=None;evaluation=None
                def journal(row,c,e):
                    nonlocal candidate,evaluation,evaluations
                    candidate=c;evaluation=e;evaluations+=1;iterations.append(row)
                    if e is not None:
                        index={int(pid):i for i,pid in enumerate(e['particle_id'])}
                        row['material_rate_branch_diagnostic']={str(pid):dict(
                            strain_rate_per_s=float(e['material']['strain_rate_per_s'][index[pid]]),
                            rate_minus_fks_floor_per_s=float(e['material']['strain_rate_per_s'][index[pid]]-p['rate_floor']),
                            plastic_increment=float(e['material']['ep'][index[pid]]-pa.ep[index[pid]]),
                            signed_velocity_residual_m_s=(c['velocity_m_s'][index[pid]]-e['map']['velocity_m_s'][index[pid]]).tolist(),
                            signed_temperature_residual_K=float(c['temperature_K'][index[pid]]-e['map']['temperature_K'][index[pid]]))
                            for pid in (268,269) if pid in index}
                    compressed(stage_dir/'latest_trial.json.gz',dict(candidate=c,evaluation=e))
                    # Every attempted numerical candidate is retained, including safeguarded rejection.
                    compressed(stage_dir/f'candidate_{len(iterations):03d}.json.gz',dict(candidate=c,row=row))
                    write_json(stage_dir/'iterations.json',iterations)
                    print(json.dumps(clean(dict(stage_index=step,total_evaluations=evaluations,**row)),allow_nan=False),flush=True)
                solved=solve_stage(trial,SETTINGS['dt_s'],min(remaining,SETTINGS['maximum_evaluations_per_step']),
                    max(1.,SETTINGS['inner_limit_s']-(time.perf_counter()-started)),journal)
                candidate,evaluation=solved['candidate'],solved['evaluation']
                checks=acceptance(pa,p,trial,candidate,evaluation)
                after=commit(pa,p,candidate,evaluation,checks)
                transport=observe(pa,after,time_s,evaluation['time_after_s'],p,cumulative)
                for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'source/input drift before commit '+rel)
                accumulate(ledger,evaluation,checks)
                cumulative=transport['cumulative']
                newstate=state(after,p,evaluation['time_after_s'],step,sources,config,parent,cumulative,ledger)
                save(stage_dir,after,p,newstate)
                loaded,loaded_state=restore(stage_dir/'accepted.npz',p,sources,config,sha(stage_dir/'accepted.npz'))
                require(particle_hash(loaded)==particle_hash(after),'roundtrip physical history mismatch')
                write_json(stage_dir/'acceptance.json',checks);write_json(stage_dir/'transport.json',transport)
                write_json(stage_dir/'tool_pose.json',newstate['tool_pose'])
                compressed(stage_dir/'accepted_trial.json.gz',dict(candidate=candidate,evaluation=evaluation))
                for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'source/input drift before publication '+rel)
                entry=dict(stage_index=step,time_s=newstate['t'],dt_s=SETTINGS['dt_s'],evaluations=len(iterations),
                    snapshot=(stage_dir/'accepted.npz').relative_to(directory).as_posix(),
                    snapshot_sha256=sha(stage_dir/'accepted.npz'),particle_state_sha256=particle_hash(after),
                    checkpoint_readback='PASS_IDENTICAL_FULL_HISTORY; RECOMPUTED_RESTART_EQUIVALENCE_NOT_RUN',
                    maximum_scaled_residual=checks['maximum_scaled_residual'],transport=transport['crossing'])
                new_history=history+[entry];write_json(directory/'accepted_history.json',new_history);history=new_history
                # Only a fully validated/published/roundtripped state supplies the next stage.
                pa=loaded;time_s=loaded_state['t']
                manifest.update(accepted_steps=len(history),accepted_progress=True,accepted_time_s=time_s,
                    last_accepted_snapshot=entry['snapshot'],total_evaluations=evaluations)
                write_json(directory/'manifest.json',manifest)
            manifest.update(status='COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION',transport_status=transport['crossing']['status'])
        for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'source/input drift '+rel)
        code=0
    except Exception as error:
        if isinstance(error,StageFailure):candidate,evaluation,iterations=error.candidate,error.evaluation,error.trace
        compressed(stage_dir/'REJECTED_candidate.json.gz',dict(classification='REJECTED',accepted=False,restart_allowed=False,
            error=repr(error),candidate=candidate,last_completed_evaluation=evaluation,iterations=iterations,
            failure_context=getattr(error,'original_context',getattr(error,'contact_context',getattr(error,'feature_context',None)))))
        manifest.update(status='BLOCKED_WITH_ACCEPTED_PREFIX' if history else 'BLOCKED_NO_ACCEPTED_PROGRESS',
            error=repr(error),accepted_steps=len(history),accepted_progress=bool(history),total_evaluations=evaluations)
        code=2
    finally:timer.cancel()
    manifest.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={f.relative_to(directory).as_posix():sha(f) for f in directory.rglob('*')
            if f.is_file() and f!=directory/'manifest.json'})
    write_json(directory/'manifest.json',manifest,rejected=code!=0)
    print(json.dumps(clean({k:v for k,v in manifest.items() if k not in ('source_sha256','input_sha256','artifact_sha256')},rejected=code!=0)),flush=True)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--probe',action='store_true');mode.add_argument('--probe-manifest',type=Path)
    args=parser.parse_args();return run(args.output,args.probe,args.probe_manifest)


if __name__=='__main__':raise SystemExit(main())

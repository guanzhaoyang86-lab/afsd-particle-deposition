"""Registered event-response candidate solve; unchanged physical stage map."""
from pathlib import Path
import argparse
import gzip
import json
import os
import sys
import threading
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from run_process_global import setup as parent_setup
from probe_updated_stage import sha,require,compressed
from run_coupled_afsd import write_json,clean,fingerprint
from run_nonlinear_continuation import accumulate
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.process_thermal_trial import ProcessThermalTrial,acceptance
from afsd_backend.pysph.stage_structural_root import solve_stage
from afsd_backend.pysph.stage_iteration import initial_candidate,StageFailure,residual_norms
from afsd_backend.pysph.event_response import response_from_trace
from afsd_backend.pysph.stage_event_preconditioner import StageEventPreconditioner
from afsd_backend.pysph.process_global_state import CLASSIFICATION,commit,state,save,restore
from afsd_backend.pysph.transport_observer import observe
from afsd_backend.pysph.application import memory_bytes
from afsd_backend.pysph.axial_actuator import control

ANCHOR=ROOT/'outputs/process_global_20261003/root_attempt'
SETTINGS=dict(dt_s=5e-7,maximum_steps=1,maximum_evaluations=32,beta=.75,history_depth=6,
    qr_condition_limit=1e10,growth_limit=2.,forecast_multiplier=1.5,overhead_s=30.,
    launch_limit_s=450.,inner_limit_s=540.,thermal_variant='UNCHANGED_BACKWARD_EULER',
    preconditioner='FROZEN_EVENT_MEAN_AND_ENDPOINT_CONTINUITY_SCHUR',
    numerical_anchor='REJECTED_CANDIDATE_ONLY_NOT_PHYSICAL_INITIAL_STATE')
SOURCES=['scripts/run_process_event.py','scripts/check_process_event.py',
    'src/afsd_backend/pysph/event_response.py','src/afsd_backend/pysph/stage_event_preconditioner.py',
    'docs/PROCESS_EVENT_ROOT_CONTRACT.md','docs/PROCESS_EVENT_RESPONSE_CONTRACT.md',
    'docs/PROCESS_EVENT_SCHUR_CONTRACT.md']


def setup():
    pa,old,p,watched,sources,_,parent=parent_setup()
    m=json.loads((ANCHOR/'manifest.json').read_text(encoding='utf-8'))
    require(m['status']=='BLOCKED_NO_NEW_ACCEPTED_STATE' and m['accepted_steps']==0,'registered rejected anchor required')
    require(m['parent_reference']==parent and m['settings']['dt_s']==SETTINGS['dt_s'],'anchor physical parent/dt mismatch')
    require(json.loads((ANCHOR/'accepted_history.json').read_text(encoding='utf-8'))==[],'rejected anchor has published state')
    watched[(ANCHOR/'manifest.json').relative_to(ROOT).as_posix()]=sha(ANCHOR/'manifest.json')
    for rel,h in m['source_sha256'].items():require(sha(ROOT/rel)==h,'anchor source drift '+rel)
    for rel,h in m['artifact_sha256'].items():
        item=(ANCHOR/rel).resolve();require(item.is_relative_to(ANCHOR) and sha(item)==h,'anchor artifact drift '+rel)
        watched[item.relative_to(ROOT).as_posix()]=h
    with gzip.open(ANCHOR/'latest_trial.json.gz','rt',encoding='utf-8') as stream:anchor=json.load(stream)
    e=anchor['evaluation']
    require(e['accepted_state_hash']==particle_hash(pa) and e['dt_s']==SETTINGS['dt_s'],'anchor initial history mismatch')
    require(np.array_equal(e['particle_id'],pa.particle_id) and e['accepted'] is False and e['restart_allowed'] is False,
        'anchor ID/disposition mismatch')
    require(e['time_before_s']==old['t'] and e['time_after_s']==old['t']+SETTINGS['dt_s'],'anchor time mismatch')
    controls=control(pa.z0,e['time_mid_s'],0.,p['bar_length'],p['actuator_depth'],p['feed'],p['release_time'])
    initial_velocity=np.column_stack((pa.u,pa.v,pa.w))
    require(len(e['boundary']['traces'])==len(pa.x),'anchor trace inventory mismatch')
    for index,record in enumerate(e['boundary']['traces']):
        require(record['particle_index']==index and record['particle_id']==pa.particle_id[index], 'anchor trace ID mismatch')
        require(np.array_equal(record['material_velocity_before_m_s'],initial_velocity[index])
            and record['mass_kg']==pa.m[index] and np.array_equal(record['free_axes'],controls['free_axes'][index])
            and np.array_equal(record['prescribed_velocity'],controls['prescribed_velocity'][index]),'anchor physical kick history/control mismatch')
    require(json.loads((ANCHOR/'resolved_config.json').read_text(encoding='utf-8'))['parameters']==p,'anchor physical parameters changed')
    sources.update({s:sha(ROOT/s) for s in SOURCES})
    config=fingerprint(dict(classification=CLASSIFICATION,parameters=p,settings=SETTINGS))
    anchor_ref=dict(path=(ANCHOR/'latest_trial.json.gz').relative_to(ROOT).as_posix(),
        sha256=sha(ANCHOR/'latest_trial.json.gz'),scope='REJECTED_NUMERICAL_LINEARIZATION_ONLY')
    return pa,old,p,watched,sources,config,parent,anchor,anchor_ref


def receipt(path,sources):
    path=Path(path).resolve();require(path.is_relative_to(ROOT),'local preflight required')
    r=json.loads(path.read_text(encoding='utf-8'))
    require(r['status']=='NUMERICAL_EVENT_PREFLIGHT_PASS_WITH_CONTINUOUS_REFERENCE_FAIL',
        'explicit numerical-only preflight with retained accuracy failure required')
    for rel,h in r['source_sha256'].items():require(sha(ROOT/rel)==h,'stale preflight '+rel)
    for rel in SOURCES:require(r['source_sha256'].get(rel)==sources[rel],'source not covered by preflight '+rel)
    for rel,h in r['artifact_sha256'].items():
        item=(path.parent/rel).resolve();require(item.is_relative_to(path.parent) and sha(item)==h,'preflight artifact drift')
    return dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),tests_run=r['tests_run'],
        scope='NUMERICAL_DERIVATIVE_ONLY; REGISTERED_CONTINUOUS_ACCURACY_FAIL_RETAINED')


def checked_report(path,sources,config,watched,parent,anchor_ref,preflight,status):
    path=Path(path).resolve();require(path.is_relative_to(ROOT),'local component report required')
    r=json.loads(path.read_text(encoding='utf-8'))
    require(r['status']==status and r['settings']==SETTINGS and r['classification']==CLASSIFICATION
        and r['accepted_steps']==0 and r['production_authorized'] is False,'wrong component scope')
    require(r['source_sha256']==sources and r['configuration_sha256']==config and r['input_sha256']==watched
        and r['parent_reference']==parent and r['numerical_anchor']==anchor_ref and r['preflight']==preflight,'stale component/parent/anchor')
    for rel,h in r['artifact_sha256'].items():
        item=(path.parent/rel).resolve();require(item.is_relative_to(path.parent) and sha(item)==h,'component artifact drift')
    return r


def response_arrays(path,pa,sources,config,watched,parent,anchor_ref,preflight):
    path=Path(path).resolve()
    r=checked_report(path,sources,config,watched,parent,anchor_ref,preflight,'COMPLETE_FROZEN_EVENT_RESPONSE_NOT_PHYSICAL_STATE')
    require(r['completed_particles']==len(pa.x),'incomplete response inventory')
    require('response_blocks.npz' in r['artifact_sha256'],'response artifact missing')
    with np.load(path.parent/'response_blocks.npz',allow_pickle=False) as d:
        require(np.array_equal(d['particle_id'],pa.particle_id),'response IDs changed')
        b=d['B'].copy();c=d['C'].copy()
    return b,c,dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),
        blocks_sha256=sha(path.parent/'response_blocks.npz'),particle_id=pa.particle_id.copy(),
        anchor=anchor_ref,scope='FROZEN_NUMERICAL_DERIVATIVE_NOT_PHYSICAL_CONTACT_REPLACEMENT')


def cost_budget(r):
    build=float(r['preconditioner_build_wall_seconds']);cost=float(r['trial_proposal_journal_wall_seconds'])
    require(np.isfinite([build,cost]).all() and build>0 and cost>0,'invalid measured costs')
    forecast=SETTINGS['forecast_multiplier']*(build+SETTINGS['maximum_evaluations']*cost)+SETTINGS['overhead_s']
    return dict(forecast_wall_seconds=forecast,allowed=forecast<=SETTINGS['launch_limit_s'],
        maximum_steps=1,maximum_evaluations=32,physical_increment_s=SETTINGS['dt_s'])


def run(directory,preflight_path,mode,response_manifest=None,cost_manifest=None):
    directory=Path(directory).resolve()
    require(directory.is_relative_to(ROOT) and not directory.exists(),'fresh project-local output required')
    require(mode in ('response','cost','solve'),'invalid mode')
    directory.mkdir(parents=True);start=time.perf_counter();candidate=evaluation=None;trace=[];published=False
    m=dict(classification=CLASSIFICATION,experiment='EVENT_RESPONSE_CANDIDATE_SCHUR',status='RUNNING_UNCOMMITTED',
        settings=SETTINGS,mode=mode,accepted_steps=0,accepted_progress=False,production_authorized=False,
        process_accuracy='NOT_VALIDATED',backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5')
    def deadline():
        try:write_json(directory/'HARD_TIMEOUT.json',dict(status='WALL_BUDGET_STOP',accepted=False,restart_allowed=False))
        finally:os._exit(124)
    timer=threading.Timer({'response':180.,'cost':80.,'solve':540.}[mode],deadline);timer.daemon=True;timer.start()
    try:
        pa,old,p,watched,sources,config,parent,anchor,anchor_ref=setup();before=particle_hash(pa)
        preflight=receipt(preflight_path,sources)
        m.update(source_sha256=sources,input_sha256=watched,configuration_sha256=config,preflight=preflight,
            parent_reference=parent,numerical_anchor=anchor_ref,source_time_s=old['t'],accepted_time_s=old['t'],
            particle_count=len(pa.x),requested_end_s=old['t']+SETTINGS['dt_s'])
        write_json(directory/'resolved_config.json',dict(classification=CLASSIFICATION,parameters=p,settings=SETTINGS,
            physical_sources='UNCHANGED_ACCEPTED_PARENT',configuration_sha256=config))
        write_json(directory/'last_accepted_reference.json',parent);write_json(directory/'accepted_history.json',[])
        write_json(directory/'manifest.json',m)
        if mode=='response':
            blocks_b=[];blocks_c=[];response_rows=[];tick=time.perf_counter()
            for index,record in enumerate(anchor['evaluation']['boundary']['traces']):
                require(record['particle_index']==index and record['particle_id']==pa.particle_id[index],'anchor trace IDs/order')
                local=time.perf_counter();r=response_from_trace(record,dt=SETTINGS['dt_s'],mu=p['friction'])
                compressed(directory/f'particle_{int(pa.particle_id[index]):06d}.json.gz',r)
                row=dict(particle_id=int(pa.particle_id[index]),available=bool(r['available']),
                    kick_calls=r['kick_calls'],wall_seconds=time.perf_counter()-local,
                    reason=r.get('reason'),contact_count=len(record['contact_rows']))
                response_rows.append(row);write_json(directory/'response_progress.json',response_rows)
                if not r['available']:
                    m['failed_response_particle_id']=row['particle_id'];raise ValueError('UNAVAILABLE_EVENT_RESPONSE '+str(row))
                blocks_b.append(r['B']);blocks_c.append(r['C'])
                m['completed_particles']=len(response_rows)
                if record['contact_rows']:print(json.dumps(row),flush=True)
            np.savez_compressed(directory/'response_blocks.npz',B=np.asarray(blocks_b),C=np.asarray(blocks_c),particle_id=pa.particle_id)
            m.update(status='COMPLETE_FROZEN_EVENT_RESPONSE_NOT_PHYSICAL_STATE',response_wall_seconds=time.perf_counter()-tick,
                total_kick_calls=sum(r['kick_calls'] for r in response_rows))
        else:
            require(response_manifest is not None,'complete local event response required')
            b,c,provenance=response_arrays(response_manifest,pa,sources,config,watched,parent,anchor_ref,preflight)
            m['response_receipt']=clean(provenance)
            if mode=='solve':
                require(cost_manifest is not None,'same-source cost probe required')
                cost=checked_report(cost_manifest,sources,config,watched,parent,anchor_ref,preflight,'COMPLETE_EVENT_COST_PROBE_UNCOMMITTED')
                require(cost['response_receipt']==clean(provenance),'different event response cost anchor')
                budget=cost_budget(cost);write_json(directory/'budget.json',budget);m['budget']=budget
                require(budget['allowed'],'NOT_RUN_COST_EXCEEDS_BUDGET')
            tick=time.perf_counter();pre=StageEventPreconditioner(pa,p,SETTINGS['dt_s'],old['t'],b,c,provenance)
            np.savez_compressed(directory/'numerical_preconditioner.npz',matrix=pre.matrix,B=b,C=c,particle_id=pa.particle_id)
            m['preconditioner_build_wall_seconds']=time.perf_counter()-tick
            trial=ProcessThermalTrial(pa,p,old['t'])
            if mode=='cost':
                candidate=initial_candidate(trial,SETTINGS['dt_s']);tick=time.perf_counter()
                evaluation=trial.evaluate(candidate,SETTINGS['dt_s']);proposal=pre.apply(candidate,evaluation['map'])
                compressed(directory/'probe_trial.json.gz',dict(candidate=candidate,evaluation=evaluation,numerical_proposal_only=proposal))
                norms,scaled=residual_norms(candidate,evaluation['map'],pa.rho)
                row=dict(evaluation=1,residual_norms=norms,maximum_scaled_residual=scaled)
                compressed(directory/'probe_candidate.json.gz',dict(candidate=candidate,row=row));write_json(directory/'iterations.json',[row])
                m.update(status='COMPLETE_EVENT_COST_PROBE_UNCOMMITTED',trial_proposal_journal_wall_seconds=time.perf_counter()-tick)
            else:
                def journal(row,c,e):
                    nonlocal candidate,evaluation
                    candidate=c;evaluation=e;trace.append(row)
                    compressed(directory/'latest_trial.json.gz',dict(candidate=c,evaluation=e))
                    compressed(directory/f'candidate_{len(trace):03d}.json.gz',dict(classification='UNACCEPTED',candidate=c,row=row))
                    write_json(directory/'iterations.json',trace);print(json.dumps(clean(row),allow_nan=False),flush=True)
                solved=solve_stage(trial,pre,SETTINGS['dt_s'],SETTINGS['maximum_evaluations'],
                    max(1.,SETTINGS['inner_limit_s']-(time.perf_counter()-start)),journal)
                candidate,evaluation=solved['candidate'],solved['evaluation']
                checks=acceptance(pa,p,trial,candidate,evaluation);after=commit(pa,p,candidate,evaluation,checks)
                transport=observe(pa,after,old['t'],evaluation['time_after_s'],p,old['transport_cumulative'])
                ledger=dict(boundary={},thermal={},energy_diagnostic_residual_J=0.);accumulate(ledger,evaluation,checks)
                s=state(after,p,evaluation['time_after_s'],1,sources,config,parent,transport['cumulative'],ledger)
                require(particle_hash(pa)==before,'parent history changed before save')
                for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'source/input drift before save '+rel)
                save(directory,after,p,s);loaded,_=restore(directory/'accepted.npz',p,sources,config,sha(directory/'accepted.npz'))
                require(particle_hash(loaded)==particle_hash(after),'accepted roundtrip failed')
                write_json(directory/'acceptance.json',checks);write_json(directory/'transport.json',transport)
                write_json(directory/'tool_pose.json',s['tool_pose']);compressed(directory/'accepted_trial.json.gz',dict(candidate=candidate,evaluation=evaluation))
                require(particle_hash(pa)==before,'parent history changed')
                for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'pre-publication source/input drift '+rel)
                entry=dict(stage_index=1,time_s=s['t'],dt_s=SETTINGS['dt_s'],evaluations=len(trace),snapshot='accepted.npz',
                    snapshot_sha256=sha(directory/'accepted.npz'),particle_state_sha256=particle_hash(after),
                    maximum_scaled_residual=checks['maximum_scaled_residual'],transport=transport['crossing'],
                    checkpoint_readback='FULL_HISTORY_IDENTICAL; RECOMPUTED_EQUIVALENCE_NOT_RUN')
                write_json(directory/'accepted_history.json',[entry]);published=True
                m.update(status='COMPLETE_EVENT_ROOT_NOT_PROCESS_VALIDATION',accepted_steps=1,accepted_progress=True,
                    accepted_time_s=s['t'],last_accepted_snapshot='accepted.npz',total_evaluations=len(trace))
        if not published:
            require(particle_hash(pa)==before,'parent state mutated')
            for rel,h in dict(watched,**sources).items():require(sha(ROOT/rel)==h,'source/input drift '+rel)
        code=0
    except Exception as error:
        if isinstance(error,StageFailure):candidate,evaluation,trace=error.candidate,error.evaluation,error.trace
        compressed(directory/'REJECTED_candidate.json.gz',dict(classification='REJECTED',accepted=False,restart_allowed=False,
            error=repr(error),candidate=candidate,last_completed_evaluation=evaluation,iterations=trace,
            failure_context=getattr(error,'original_context',getattr(error,'contact_context',None))))
        m.update(status='FAILED_AFTER_PUBLISHED_STATE' if published else 'BLOCKED_NO_NEW_ACCEPTED_STATE',
            error=repr(error),total_evaluations=len(trace));code=2
    finally:timer.cancel()
    m.update(wall_seconds=time.perf_counter()-start,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={f.relative_to(directory).as_posix():sha(f) for f in directory.rglob('*') if f.is_file() and f!=directory/'manifest.json'})
    write_json(directory/'manifest.json',m,rejected=code!=0)
    print(json.dumps(clean({k:v for k,v in m.items() if k not in ('source_sha256','input_sha256','artifact_sha256','response_receipt')},rejected=code!=0)),flush=True)
    return code


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--preflight-manifest',type=Path,required=True);ap.add_argument('--mode',choices=['response','cost','solve'],required=True)
    ap.add_argument('--response-manifest',type=Path);ap.add_argument('--cost-manifest',type=Path)
    a=ap.parse_args();return run(a.output,a.preflight_manifest,a.mode,a.response_manifest,a.cost_manifest)


if __name__=='__main__':raise SystemExit(main())

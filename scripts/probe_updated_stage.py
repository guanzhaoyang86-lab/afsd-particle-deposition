"""One registered 50 ns event-average global trial; never a production run."""
from pathlib import Path
import argparse
import gzip
import json
import os
import sys
import threading
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from assess_process_route import checked_inputs, CASE, sha, require
from run_coupled_afsd import clean, write_json, fingerprint
from afsd_backend.pysph.updated_trial import UpdatedTrial, particle_hash, CLASSIFICATION
from afsd_backend.pysph.stage_iteration import initial_candidate, solve_stage, StageFailure, residual_norms
from afsd_backend.pysph.output import restore, checkpoint
from afsd_backend.pysph.mechanics import positions, velocities, pressure, volumetric_energy
from afsd_backend.pysph.material import equivalent
from afsd_backend.pysph.application import memory_bytes
from afsd_scaffold.vtk_legacy import write_point_cloud

NEW_SOURCES = ['scripts/probe_updated_stage.py', 'scripts/assess_process_route.py',
    'src/afsd_backend/pysph/updated_trial.py', 'src/afsd_backend/pysph/stage_iteration.py',
    'docs/UPDATED_TRIAL_CONTRACT.md']
SETTINGS = dict(dt_s=5e-8, maximum_evaluations=32, relaxation=.75,
    forecast_multiplier=1.5, overhead_s=30., launch_limit_s=450., inner_limit_s=540.)


def inputs():
    manifest, hashes = checked_inputs(CASE)
    path = CASE/'plastic_accepted.npz'
    pa, state = restore(path)
    p = json.loads((CASE/'resolved_config.json').read_text(encoding='utf-8'))['parameters']
    require(state['accepted'] is True and state['count'] == 62 and len(pa.x) == 342,
            'registered nonzero plastic snapshot required')
    require(state['source_hashes'] == manifest['source_sha256'] and state['parameters'] == p and
            state['configuration_hash'] == manifest['configuration_sha256'], 'snapshot provenance mismatch')
    require(np.any(pa.ep > 0) and np.any(pa.S != 0), 'plastic/stress history missing')
    sources = dict(manifest['source_sha256'], **{s: sha(ROOT/s) for s in NEW_SOURCES})
    return pa, state, p, hashes, sources


def compressed(path, data):
    path = Path(path); pending = path.with_suffix(path.suffix+'.pending')
    with gzip.open(pending, 'wt', encoding='utf-8') as stream:
        json.dump(clean(data, rejected=True), stream, allow_nan=False, separators=(',', ':'))
    pending.replace(path)


def budget(probe, source, config_hash):
    probe = Path(probe).resolve()
    require(probe.is_relative_to(ROOT), 'project-local probe required')
    m = json.loads(probe.read_text(encoding='utf-8'))
    require(m['status'] == 'COMPLETE_COST_PROBE_UNCOMMITTED', 'complete successful cost probe required')
    require(m['classification'] == CLASSIFICATION and m['mode'] == 'COST_PROBE' and
            m['accepted_progress'] is False and m['production_authorized'] is False and
            m['settings'] == SETTINGS, 'wrong cost-probe scope/settings')
    require(m['source_sha256'] == source and m['configuration_sha256'] == config_hash, 'stale cost probe')
    require({'resolved_config.json','last_accepted_reference.json','uncommitted_trial.json.gz'} <=
            set(m['artifact_sha256']), 'incomplete cost-probe evidence')
    for name, digest in m['artifact_sha256'].items():
        require(Path(name).name == name and sha(probe.parent/name) == digest, 'probe artifact drift')
    elapsed = float(m['complete_trial_wall_seconds'])
    require(np.isfinite(elapsed) and elapsed > 0, 'invalid probe elapsed time')
    predicted = SETTINGS['forecast_multiplier']*elapsed*SETTINGS['maximum_evaluations']+SETTINGS['overhead_s']
    return dict(probe=str(probe.relative_to(ROOT)), probe_sha256=sha(probe),
                forecast_wall_seconds=predicted, allowed=predicted <= SETTINGS['launch_limit_s'])


def acceptance(pa, p, trial, c, e):
    """Checks refer to the evaluated candidate, never snap it to an unevaluated map."""
    require(e['classification']==CLASSIFICATION and e['accepted'] is False and e['restart_allowed'] is False,
            'uncommitted registered trial metadata required')
    norms, scaled = residual_norms(c, e['map'], pa.rho)
    require(scaled <= 1, 'global residual not converged')
    require(trial.accepted_hash() == particle_hash(pa) == e['accepted_state_hash'], 'accepted history changed')
    require(np.array_equal(e['particle_id'], pa.particle_id), 'ID correspondence lost')
    require(e['material_scope'] == 'UNCHANGED_FKS_APPROXIMATE_PARAMETERS', 'manufactured material cannot be accepted')
    s = e['spatial']; b = e['boundary']; ledger = b['ledger']; th = e['thermal']
    require(abs(s['power_conjugacy_residual_W']) <= 1e-12+1e-10*s['power_absolute_scale_W'], 'spatial power')
    require(np.linalg.norm(s['internal_net_force_N']) <= 1e-12+1e-10*s['force_absolute_scale_N'], 'internal net force')
    require(np.linalg.norm(s['internal_net_torque_Nm']) <= 1e-14+1e-10*s['torque_absolute_scale_Nm'], 'internal net torque')
    source_heat = th['particle_energy_J']['total_heat_input_J']
    thermal_residual = float(np.sum(pa.m*p['cp']*(c['temperature_K']-pa.temperature)-source_heat))
    require(abs(thermal_residual) <= 1e-12+1e-9*float(np.abs(source_heat).sum()), 'saved candidate thermal balance')
    require(abs(th['ledger']['internal_conduction_J']) <= 1e-15+1e-10*float(np.abs(th['pair_power_W']).sum()*e['dt_s']), 'pair heat balance')
    work_scale = sum(abs(ledger[k]) for k in ('delta_kinetic_J','wall_work_J','actuator_work_J','internal_work_J','normal_relative_work_J','friction_raw_J'))
    for key in ('kinetic_work_residual_J', 'contact_kinetic_residual_J'):
        require(abs(ledger[key]) <= 1e-15+1e-9*work_scale, key)
    seen = set(); max_momentum = 0.; intervals = 0
    for trace in b['traces']:
        idx = trace['particle_index']; require(idx not in seen and trace['particle_id'] == pa.particle_id[idx], 'trace IDs')
        seen.add(idx); result = trace['result']; free = np.asarray(trace['free_axes'])
        reaction = np.asarray(result['constraint_impulse'])
        require(np.max(np.abs(free*reaction)) <= 1e-17, 'reaction on free axis')
        require(np.max(np.abs((1-free)*(c['velocity_m_s'][idx]-trace['prescribed_velocity']))) <= 1e-12, 'axial endpoint constraint')
        impulse = np.asarray(result['normal_impulse']).reshape(-1,3).sum(axis=0)+np.asarray(result['tangential_impulse']).reshape(-1,3).sum(axis=0)+reaction
        defect = pa.m[idx]*(c['velocity_m_s'][idx]-velocities(pa)[idx])-e['dt_s']*e['internal_force_N'][idx]-impulse
        max_momentum = max(max_momentum, float(np.max(np.abs(defect))))
        require(np.max(np.abs(defect)) <= pa.m[idx]*1e-9+1e-17, 'candidate momentum identity')
        mean = np.zeros(3); covered = 0.
        for interval in result['intervals']:
            intervals += 1; h = interval['dt_s']; covered += h
            require(h > 0, 'nonpositive contact subinterval')
            mean += h*.5*(np.asarray(interval['v_before'])+interval['v_after'])/e['dt_s']
            require(np.all(np.asarray(interval['friction_heat_raw_J']) >= -np.asarray(interval['friction_rounding_bound_J'])), 'negative raw friction beyond roundoff')
        require(abs(covered-e['dt_s']) <= 128*np.finfo(float).eps*e['dt_s'], 'incomplete contact interval')
        require(np.max(np.abs(mean-c['mean_velocity_m_s'][idx])) <= 1e-9, 'event-average velocity mismatch')
    require(len(seen) == len(pa.x), 'missing particle traces')
    require(np.all(e['material']['ep'] >= pa.ep), 'plastic history decreased')
    require(e['dt_s'] <= th['thermal_dt_bound_s'], 'thermal dt bound')
    v0 = velocities(pa); v1 = c['velocity_m_s']; rho1 = c['density_kg_m3']; S0 = pa.S.reshape(-1,3,3); S1 = e['material']['deviatoric_stress_Pa']
    dK = float(.5*np.sum(pa.m[:,None]*(v1*v1-v0*v0)))
    dEs = float(np.sum(pa.m/rho1*np.sum(S1*S1,axis=(1,2))/(4*p['shear'])-pa.m/pa.rho*np.sum(S0*S0,axis=(1,2))/(4*p['shear'])))
    dEv = float(np.sum(pa.m*(volumetric_energy(rho1,p['rho0'],trial.bulk)-volumetric_energy(pa.rho,p['rho0'],trial.bulk))))
    dU = e['normal_after']['normal_potential_J']-e['normal_before']['normal_potential_J']
    dEth = float(np.sum(pa.m*p['cp']*(c['temperature_K']-pa.temperature)))
    diagnostic = dK+dEs+dEv+dU+dEth+th['ledger']['plastic_unthermalized_J']-ledger['wall_work_J']-ledger['actuator_work_J']+th['ledger']['wall_heat_out_J']
    x0=positions(pa); x1=c['position_m']; xmid=e['midpoint']['position_m']
    angular_delta=np.sum(pa.m[:,None]*(np.cross(x1,v1)-np.cross(x0,v0)),axis=0)
    external_angular=sum((np.asarray(ledger[k+'_angular_impulse_Nms']) for k in ('base','tool','actuator')),np.zeros(3))
    extra_kinematic=np.sum(pa.m[:,None]*np.cross(x1-x0,.5*(v0+v1)),axis=0)
    midpoint_internal_torque=e['dt_s']*np.cross(xmid,e['internal_force_N']).sum(axis=0)
    nonlinear_work=e['dt_s']*float(np.sum(e['internal_force_N']*(c['mean_velocity_m_s']-b['work_velocity'])))
    return dict(status='ACCEPTED_NUMERICAL_STAGE_APPROXIMATE_NOT_ENERGY_CERTIFIED', residual_norms=norms,
        maximum_scaled_residual=scaled, momentum_residual_max_Ns=max_momentum, thermal_residual_J=thermal_residual,
        contact_subintervals=intervals, spatial=s, contact_ledger=ledger, thermal_ledger=th['ledger'],
        delta_K_J=dK, delta_Es_diagnostic_J=dEs, delta_volumetric_diagnostic_J=dEv,
        candidate_minus_contact_map_delta_K_J=dK-ledger['delta_kinetic_J'],
        delta_normal_potential_J=dU, delta_thermal_J=dEth, total_energy_diagnostic_residual_J=diagnostic,
        normal_temporal_residual_J=dU+ledger['normal_relative_work_J'],
        nonlinear_internal_work_defect_J=nonlinear_work,
        angular_momentum_change_kg_m2_s=angular_delta,external_midpoint_angular_impulse_Nms=external_angular,
        angular_residual_Nms=angular_delta-external_angular,
        event_average_kinematic_angular_term_Nms=extra_kinematic,
        angular_identity_defect_Nms=angular_delta-external_angular-midpoint_internal_torque-extra_kinematic,
        angular_status='FROZEN_MIDPOINT_MOMENT_APPROXIMATE_EXTRA_KINEMATIC_TERM_REPORTED_NOT_CORRECTED',
        total_energy_status='DIAGNOSTIC_HYPOELASTIC_NOT_FREE_ENERGY_CERTIFIED_NO_RESIDUAL_HEAT',
        plastic_particle_ids=pa.particle_id[e['material']['ep']>pa.ep], production_authorized=False)


def save_stage(directory, pa, state, p, c, e, checks, sources, config_hash, bulk):
    accepted = pa.extract_particles(np.arange(len(pa.x), dtype=np.uint32))
    accepted.x[:],accepted.y[:],accepted.z[:] = c['position_m'].T
    accepted.u[:],accepted.v[:],accepted.w[:] = c['velocity_m_s'].T
    accepted.rho[:] = c['density_kg_m3']; accepted.temperature[:] = c['temperature_K']
    accepted.S[:] = e['material']['deviatoric_stress_Pa'].ravel(); accepted.ep[:] = e['material']['ep']
    accepted.strain_rate[:] = e['material']['strain_rate_per_s']; accepted.ep_rate[:] = (accepted.ep-pa.ep)/e['dt_s']
    gradient=e['midpoint']['velocity_gradient_per_s']; accepted.D[:] = (.5*(gradient+gradient.swapaxes(1,2))).ravel()
    # Store actual history increments; do not reuse old per-volume work arrays with new semantics.
    record = dict(classification=CLASSIFICATION, accepted=True, restart_allowed=False,
        t=e['time_after_s'], count=state['count']+1, parameters=p, source_sha256=sources,
        configuration_sha256=config_hash, checks=clean(checks), production_authorized=False,
        restart_status='SNAPSHOT_ONLY_NEW_INTEGRATOR_RESTART_NOT_VALIDATED',
        valid_properties=['x','y','z','u','v','w','rho','temperature','S','D','ep','strain_rate','ep_rate',
            'particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','as_rs_tag','m','h','x0','y0','z0'],
        inherited_scratch_fields='p/e/plastic_work/return_dissipation and derivative scratch are NOT solved stage outputs')
    checkpoint(directory/'accepted_snapshot.npz', accepted, record)
    scalar = {key: getattr(accepted,key).tolist() for key in ('particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','as_rs_tag')}
    pr = pressure(accepted.rho,p['rho0'],bulk); S = accepted.S.reshape(-1,3,3); sigma = S-pr[:,None,None]*np.eye(3)
    scalar.update(mass_kg=accepted.m.tolist(),density_kg_per_m3=accepted.rho.tolist(),temperature_K=accepted.temperature.tolist(),
        pressure_Pa=pr.tolist(),equivalent_stress_Pa=equivalent(S).tolist(),equivalent_plastic_strain=accepted.ep.tolist())
    for a,b,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')): scalar['stress_'+name+'_Pa']=sigma[:,a,b].tolist()
    write_point_cloud(directory/'accepted_stage.vtk',positions(accepted).tolist(),scalar,
        dict(velocity_m_per_s=velocities(accepted).tolist(),interval_mean_velocity_m_per_s=c['mean_velocity_m_s'].tolist()),
        title=CLASSIFICATION+'; NOT_PAPER_VALIDATION; time_s='+str(e['time_after_s']))


def run(directory, probe=False, probe_manifest=None):
    directory = Path(directory).resolve()
    require(directory.is_relative_to(ROOT) and not directory.exists(), 'fresh project-local output directory required')
    directory.mkdir(parents=True); started=time.perf_counter()
    def timeout():
        write_json(directory/'HARD_TIMEOUT.json',dict(status='REJECTED_WALL_TIMEOUT',accepted=False,production_authorized=False))
        os._exit(124)
    timer=threading.Timer(80. if probe else SETTINGS['inner_limit_s'],timeout); timer.daemon=True; timer.start()
    m=dict(classification=CLASSIFICATION,status='RUNNING_UNCOMMITTED',accepted_progress=False,production_authorized=False,
           settings=SETTINGS, backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5',
           mode='COST_PROBE' if probe else 'SINGLE_GLOBAL_ATTEMPT')
    candidate=None; evaluated=None; trace=[]; pa=None; state=None
    try:
        pa,state,p,hashes,sources=inputs(); trial=UpdatedTrial(pa,p,state['t'])
        config_hash=fingerprint(dict(parameters=p,stage_settings=SETTINGS,classification=CLASSIFICATION))
        m.update(source_sha256=sources,input_sha256=hashes,configuration_sha256=config_hash,
            parent_configuration_sha256=state['configuration_hash'],accepted_initial_hash=particle_hash(pa),
            particle_count=len(pa.x),source_time_s=state['t'],requested_time_s=state['t']+SETTINGS['dt_s'])
        write_json(directory/'resolved_config.json',dict(parameters=p,settings=SETTINGS,source='UNCHANGED_REGISTERED_COUPLED_AFSD'))
        write_json(directory/'last_accepted_reference.json',dict(path=str((CASE/'plastic_accepted.npz').relative_to(ROOT)),
            sha256=sha(CASE/'plastic_accepted.npz'),time_s=state['t'],step=state['count'],classification=state['classification']))
        write_json(directory/'manifest.json',m)
        candidate=initial_candidate(trial,SETTINGS['dt_s']); tick=time.perf_counter()
        if probe:
            evaluated=trial.evaluate(candidate,SETTINGS['dt_s'])
            compressed(directory/'uncommitted_trial.json.gz',dict(candidate=candidate,evaluation=evaluated))
            m.update(status='COMPLETE_COST_PROBE_UNCOMMITTED',complete_trial_wall_seconds=time.perf_counter()-tick,
                candidate_metrics=dict(spatial=evaluated['spatial'],contact=evaluated['boundary']['ledger'],thermal=evaluated['thermal']['ledger']),
                root_status='NOT_RUN_COST_PROBE_ONLY')
        else:
            require(probe_manifest is not None,'measured same-source cost probe required')
            receipt=budget(probe_manifest,sources,config_hash); m['budget']=receipt
            write_json(directory/'budget.json',receipt)
            require(receipt['allowed'],'NOT_RUN_COST_EXCEEDS_REGISTERED_BUDGET')
            def journal(row,c,e):
                nonlocal candidate,evaluated,trace
                candidate=c; evaluated=e; trace.append(row)
                compressed(directory/'latest_uncommitted_trial.json.gz',dict(candidate=c,evaluation=e))
                write_json(directory/'iterations.json',trace)
                print(json.dumps(clean(row)),flush=True)
            solved=solve_stage(trial,SETTINGS['dt_s'],SETTINGS['maximum_evaluations'],SETTINGS['relaxation'],
                max(1.,SETTINGS['inner_limit_s']-(time.perf_counter()-started)),journal)
            candidate=solved['candidate']; evaluated=solved['evaluation']
            checks=acceptance(pa,p,trial,candidate,evaluated)
            # Verify provenance again BEFORE any accepted output is written.
            for rel,h in dict(hashes,**sources).items(): require(sha(ROOT/rel)==h,'input/source drift '+rel)
            save_stage(directory,pa,state,p,candidate,evaluated,checks,sources,config_hash,trial.bulk)
            write_json(directory/'acceptance.json',checks)
            m.update(status='ACCEPTED_ONE_APPROXIMATE_GLOBAL_STAGE_NOT_SUSTAINED_DEPOSITION',accepted_progress=True,
                     stage_acceptance=checks,accepted_time_s=evaluated['time_after_s'],evaluations=len(trace))
        require(particle_hash(pa)==m['accepted_initial_hash'],'caller accepted state mutated')
        for rel,h in dict(hashes,**sources).items(): require(sha(ROOT/rel)==h,'input/source drift '+rel)
        code=0
    except Exception as error:
        if isinstance(error,StageFailure): candidate,error_trace,evaluated=error.candidate,error.trace,error.evaluation; trace=error_trace
        context=getattr(error,'original_context',getattr(error,'contact_context',getattr(error,'feature_context',None)))
        compressed(directory/'REJECTED_candidate.json.gz',dict(classification='REJECTED',accepted=False,restart_allowed=False,
            error=repr(error),candidate=candidate,last_completed_evaluation=evaluated,iterations=trace,failure_context=context,
            last_completed_evaluation_may_refer_to_previous_candidate=True))
        m.update(status='BLOCKED_REJECTED_GLOBAL_STAGE' if not probe else 'BLOCKED_REJECTED_COST_PROBE',
            error=repr(error),accepted_progress=False,evaluations=len(trace),root_status='NOT_RUN' if probe else 'REJECTED')
        code=2
    finally:
        timer.cancel()
    m.update(wall_seconds=time.perf_counter()-started,peak_process_working_set_bytes=memory_bytes(),
        artifact_sha256={f.name:sha(f) for f in directory.iterdir() if f.is_file() and f.name!='manifest.json'})
    write_json(directory/'manifest.json',m,rejected=code!=0)
    print(json.dumps(clean({k:v for k,v in m.items() if k not in ('source_sha256','input_sha256','artifact_sha256','candidate_metrics','stage_acceptance')},rejected=code!=0)),flush=True)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--probe',action='store_true');mode.add_argument('--probe-manifest',type=Path)
    args=parser.parse_args();return run(args.output,args.probe,args.probe_manifest)


if __name__=='__main__':raise SystemExit(main())

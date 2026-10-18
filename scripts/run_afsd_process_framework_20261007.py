"""Supervised independent cost/probe stages; never publish a physical restart.

The 0.16 s machine program is only a plan. This command executes a separately
bounded startup diagnostic with common-initial coarse/fine comparisons.
"""
from pathlib import Path
import argparse
import json
import math
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import afsd_process_framework_20261007 as model
import ramped_hot_candidate as inherited
from afsd_process_observer_20261007 import observe_process
from visible_support_events_20261007 import audit_visibility_events
from visible_support_path_tight import audit_visible_support_path

core, post, event, batch, solver, old = (inherited.core, inherited.post, inherited.event,
                                      inherited.batch, inherited.solver, inherited.old)
ROOT = model.ROOT
CONTRACT = ROOT/'docs/PROCESS_FRAMEWORK_CONTRACT_20261007.md'
DEFAULT_CHECKS = ROOT/'outputs/process_framework_20261007/components_v1.json'
REFERENCE_COST = ROOT/'outputs/expanded_hot_domain_cost_20261004/v3/dt_15_625us/cost/manifest.json'
DTS = dict(coarse=.0005, fine=.00025)
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')
require = model.require
COST_PASS = 'PASS_PROCESS_FRAMEWORK_MEASURED_MAP_COST_NOT_PHYSICAL_VALIDATION'
PROBE_PASS = 'PASS_PROCESS_FRAMEWORK_SHORT_TIME_SCREEN_NOT_ACCEPTED'


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def load(path):
    saved = core.load(str(path))
    return saved['arrays']['solid'], saved['solver_data']


def write_new(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(core.clean(value), stream, ensure_ascii=False, indent=2, allow_nan=False)


def verify(mapping):
    for name, digest in mapping.items():
        path = (ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and core.sha(path) == digest, 'source/evidence drift: '+name)


def source_inputs(config, checks):
    require(checks.is_file(), 'executed component/operator receipt required before a physical stage')
    receipt = read(checks)
    require(str(receipt.get('status', '')).startswith('PASS'), 'component receipt did not pass')
    rows = receipt.get('checks', [])
    require(isinstance(rows, list) and rows and all(row.get('status') == 'PASS'
            and row.get('exit_code', 0) == 0 for row in rows), 'executed component checks required')
    require(receipt.get('operator_reuse_audit', {}).get('status') ==
            'PASS_UNCHANGED_INHERITED_THERMOMECHANICAL_OPERATORS', 'operator reuse receipt required')
    require(receipt.get('source_sha256'), 'component source binding absent')
    mapping = dict(model.source_hashes(config))
    for key in ('source_sha256', 'output_sha256', 'evidence_sha256'):
        values = receipt.get(key, {})
        verify(values)
        for name, digest in values.items():
            require(name not in mapping or mapping[name] == digest, 'conflicting source binding: '+name)
            mapping[name] = digest
    files = {Path(__file__).resolve(), CONTRACT, config, checks,
        ROOT/'scripts/afsd_process_observer_20261007.py', ROOT/'scripts/visible_support_events_20261007.py',
        ROOT/'scripts/test_afsd_process_observer_20261007.py', ROOT/'scripts/test_visible_support_events_20261007.py',
        ROOT/'docs/VISIBLE_SUPPORT_EVENT_ROUTING_20261007.md'}
    # Capture the actual imported script closure and lazy backend module set.
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT/'scripts'):
                files.add(path)
    files.update((ROOT/'src/afsd_backend/pysph').glob('*.py'))
    for path in files:
        name, digest = rel(path), core.sha(path)
        require(name not in mapping or mapping[name] == digest, 'source receipt changed: '+name)
        mapping[name] = digest
    verify(mapping)
    model.operator_reuse_audit()
    return mapping


def cost_receipt(base):
    directory = base/'cost'
    meta, supervision = read(directory/'manifest.json'), read(directory/'supervision.json')
    require(meta['status'] == COST_PASS and supervision['exit_code'] == 0
        and supervision['source_inputs_still_match']
        and supervision['manifest_sha256'] == core.sha(directory/'manifest.json')
        and supervision['log_sha256'] == core.sha(directory/'run.log'), 'successful matching cost receipt required')
    verify(meta['source_sha256']); verify(meta['artifacts_sha256'])
    return meta


def preflight(base, stage, config, checks, end_time):
    frozen = source_inputs(config, checks)
    _, p, program, domain = model.load_configuration(config)
    require(math.isfinite(end_time) and 0 < end_time <= program.end_s,
            'finite positive short endpoint within the machine program required')
    counts = {name: round(end_time/dt) for name, dt in DTS.items()}
    require(all(counts[name] >= 1 and abs(counts[name]*dt-end_time) <= 1e-14
        for name, dt in DTS.items()) and sum(counts.values()) <= 64,
        'endpoint must be a coarse-step multiple with <=64 total coarse+fine steps')
    # Both diagnostic branches must respect actual phase boundaries.
    for name, dt in DTS.items():
        t = 0.
        for _ in range(counts[name]):
            program.step(t, dt); t += dt
    references = {}
    if stage == 'cost':
        require(not base.exists(), 'cost stage requires a new output base')
        prior = read(REFERENCE_COST)
        require(str(prior['status']).startswith('PASS') and prior['particle_count'] > 0,
                'executed prior CPU cost evidence required for initial forecast')
        ratio = domain['particle_count']/prior['particle_count']
        estimated_map = prior['full_map_wall_s']*max(1., ratio**2)
        estimated_tangent = prior['tangent_linear_wall_s']*max(1., ratio**3)
        forecast = 60.+1.5*(3*estimated_map+2*estimated_tangent)
        measured = None
        references[rel(REFERENCE_COST)] = core.sha(REFERENCE_COST)
        estimate = dict(kind='PRIOR_MEASURED_COST_WITH_DECLARED_N_SQUARED_MAP_N_CUBED_DENSE_TANGENT_SCALING',
                        reference=rel(REFERENCE_COST), reference_particle_count=prior['particle_count'],
                        particle_count_ratio=ratio, estimated_map_s=estimated_map,
                        estimated_tangent_linear_s=estimated_tangent, fixed_io_margin_s=60.)
    else:
        meta = cost_receipt(base)
        require(meta['source_sha256'] == frozen, 'implementation changed after cost stage')
        require(not (base/'probe').exists(), 'refuse to overwrite an existing probe')
        measured = float(meta['candidate_forecast_s'])
        forecast = 30.+1.5*sum(counts.values())*measured
        estimate = dict(kind='MEASURED_CURRENT_MAP_TANGENT_LINEAR_COST_WITH_DECLARED_ROOT_ALLOWANCE',
            formula='30 + 1.5 * total_candidate_count * (12*map_s + 2*tangent_linear_s + 3)',
            candidate_forecast_s=measured, predicted_candidates_per_second=1./measured,
            full_map_wall_s=meta['full_map_wall_s'], tangent_linear_wall_s=meta['tangent_linear_wall_s'],
            stable_timestep='NOT_ESTABLISHED_REQUIRES_THIS_PROBE_AND_LOCAL_REFINEMENT')
        for path in (base/'cost').rglob('*'):
            if path.is_file():
                references[rel(path)] = core.sha(path)
        for path in (base/'initial').rglob('*'):
            if path.is_file():
                references[rel(path)] = core.sha(path)
    require(forecast <= 450., 'predicted stage exceeds 450 s; choose a shorter explicit probe')
    n = domain['particle_count']
    return dict(schema='INDEPENDENT_PROCESS_FRAMEWORK_DIAGNOSTIC_20261007', stage=stage,
        configuration=rel(config), component_receipt=rel(checks), source_sha256=frozen,
        reference_sha256=references, domain_inventory_plan=domain, end_time_s=end_time,
        counts=counts, dts_s=DTS, total_candidate_count=sum(counts.values()) if stage == 'probe' else 0,
        forecast_s=forecast, forecast_limit_s=450., cost_estimate=estimate,
        measured_candidate_forecast_s=measured, internal_limit_s=540., external_limit_s=600.,
        particle_count=n, dense_full_tangent_bytes=8*(3*n)**2,
        memory_plan='At least full dense tangent plus free matrix/inverse; actual process peak recorded.',
        classification=model.CLASSIFICATION, accepted=False, restart_allowed=False,
        accepted_steps=0, accepted_time_advanced_s=0., full_machine_program_executed=False)


def scalar_equivalence(pa, state, dt, trial, result, directory):
    scalar = model.ProcessTrial(pa, state, dt)
    p = trial.p
    scalar.contact = old.ChainRuleDiscreteContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2,
        normal_stiffness=p['normal_stiffness'], tangential_stiffness=p['normal_stiffness'],
        friction=p['friction'], base_z=-p['substrate_depth'])
    sr = scalar.evaluate(result['displacement'], dt)
    differences = batch.candidate_equivalence(
        dict(sr, energy_comparison=core.energy_diagnostics(scalar, sr)),
        dict(result, energy_comparison=core.energy_diagnostics(trial, result)))
    require(sr['graph_fingerprint'] == result['graph_fingerprint']
        and np.array_equal(sr['g'], result['g'])
        and np.array_equal(sr['thermal']['pair_conductance_W_K'], result['thermal']['pair_conductance_W_K']),
        'scalar and batch must share exact graph and heat conductances')
    core.write(directory/'scalar_batch_equivalence.json', dict(status='PASS_SAME_COMPLETE_PROCESS_MAP', differences=differences))
    return differences


def complete_candidate(pa, state, dt, directory, deadline, per_root, frozen):
    """Solve and screen one detached diagnostic candidate; never accept it."""
    started = time.perf_counter()
    trial = model.ProcessTrial(pa, state, dt)
    guess = batch.guess_for(trial, dt)
    attempt = directory/'broyden'; attempt.mkdir()
    fallback = None
    try:
        result, metrics = solver.solve(trial, dt, guess, attempt, deadline, 24)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        require(time.perf_counter()+1.5*per_root+20 < deadline,
                'insufficient budget for one same-parent unchanged Newton retry')
        fallback = repr(error); attempt = directory/'newton_retry'; attempt.mkdir()
        result, metrics = core.solve(trial, dt, guess, attempt, deadline, 24)
        metrics['controller'] = 'ORIGINAL_NEWTON_SAME_PARENT_PROCESS_MAP'
    metrics['fallback_reason'] = fallback
    core.journal(directory, trial, result, 'SOLVED_UNCOMMITTED_candidate')
    old.save_graph(directory, trial, result)
    nextpa, nextstate, fresh = model.advance_candidate(trial, result)
    require(np.array_equal(pa.particle_id, nextpa.particle_id)
        and np.array_equal(pa.m, nextpa.m)
        and abs(nextpa.m.sum()/state['initial_mass_kg']-1.) <= 1e-12, 'candidate mass/identity changed')
    target, evidence = event.save_candidate(directory, trial, nextpa, nextstate, fresh, frozen)
    old.save_graph(directory, trial, fresh)
    core.write(directory/'observer.json', observe_process(pa, nextstate, fresh, trial.motion, trial.p))
    core.write(directory/'motion.json', trial.motion.to_dict())
    core.write(directory/'transport.json', nextstate['bore_observation'])
    scalar_equivalence(pa, state, dt, trial, fresh, directory)
    geometry, events = inherited.geometry_audit(trial, fresh)
    core.write(directory/'geometry_audit.json', geometry)
    core.write(directory/'sampled_events.json', events)
    require(not geometry['transient_contacts_missed_by_endpoints']
        and set(geometry['new_non_bore_features']) <= {'SHOULDER', 'TOP'},
        'unqualified or endpoint-hidden contact feature')
    pose0, pose1 = trial.motion.pose(trial.motion.t0), trial.motion.pose(trial.motion.t1)
    shell = old.audit_contact_shell_finite_enclosure(trial.contact, trial.x0, fresh['position'], pose0, pose1)
    core.write(directory/'continuous_shell_audit.json', shell)
    require(shell['status'] == event.CERTIFIED, 'continuous contact shell unresolved')
    path = audit_visible_support_path(trial.x0, fresh['position'], fresh['i'], fresh['j'],
                                     pose0, pose1, trial.p, pa.particle_id)
    path.update(candidate_sha256=core.sha(target), graph_sha256=core.sha(directory/'visible_graph.npz'))
    core.write(directory/'visible_path_audit.json', path)
    if path['status'] != old.CERTIFIED:
        try:
            routing = audit_visibility_events(trial.x0, fresh['position'], fresh['i'], fresh['j'],
                                               pose0, pose1, trial.p, pa.particle_id)
        except Exception as error:
            routing = dict(status='EVENT_DIAGNOSTIC_FAILED_REQUIRES_RESOLVE', error=repr(error),
                           requires_resolve=True, acceptance_allowed=False)
        core.write(directory/'visibility_event_routing.json', routing)
        raise ValueError('whole-path visibility rejected; see event routing, coupled resolve still required')
    core.write(directory/'candidate_qualification.json', dict(status='PASS_LOCAL_CANDIDATE_GATES_REQUIRES_TIME_SCREEN',
        accepted=False, restart_allowed=False, accepted_time_advanced_s=0., metrics=metrics))
    core.vtk(directory/'fields.vtk', nextpa, trial.p, trial.bulk, nextstate['time_s'], model.CLASSIFICATION)
    return dict(pa=nextpa, state=nextstate, checkpoint=target, evidence=evidence,
                metrics=metrics, wall_seconds=time.perf_counter()-started)


def cost_stage(base, directory, protocol, report):
    pa, state = model.independent_initial_state(ROOT/protocol['configuration'])
    initial = base/'initial'; initial.mkdir()
    path = initial/'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
    core.checkpoint(path, pa, state); core.verify_snapshot(path, pa, state)
    core.write(initial/'state.json', state)
    trial = model.ProcessTrial(pa, state, DTS['coarse'])
    guess = batch.guess_for(trial, DTS['coarse'])
    started = time.perf_counter(); result = trial.evaluate(guess, DTS['coarse']); map_wall = time.perf_counter()-started
    core.journal(directory, trial, result, 'UNSOLVED_COST_ONLY')
    old.save_graph(directory, trial, result)
    started = time.perf_counter(); tangent = trial.tangent(result); tangent_wall = time.perf_counter()-started
    require(np.isfinite(tangent).all(), 'nonfinite initial tangent')
    free = state['free'].ravel()
    matrix = tangent[np.ix_(free, free)]
    started = time.perf_counter()
    inverse = np.linalg.inv(matrix)
    direction = np.linalg.solve(matrix, -result['force'][state['free']])
    linear_wall = time.perf_counter()-started
    condition = float(np.linalg.norm(matrix, np.inf)*np.linalg.norm(inverse, np.inf))
    require(np.isfinite(direction).all() and np.isfinite(condition) and condition <= 1e12,
            'initial free tangent near singular; no regularization authorized')
    rhs = -result['force'][state['free']]
    residual = float(np.linalg.norm(matrix@direction-rhs, np.inf)/max(1., np.linalg.norm(rhs, np.inf)))
    require(residual <= 1e-8, 'initial linear direction residual failed')
    differences = scalar_equivalence(pa, state, DTS['coarse'], trial, result, directory)
    candidate_forecast = 12*map_wall+2*(tangent_wall+linear_wall)+3.
    report.update(status=COST_PASS, initial_checkpoint=rel(path), initial_sha256=core.sha(path),
        full_map_wall_s=map_wall, tangent_wall_s=tangent_wall, linear_wall_s=linear_wall,
        tangent_linear_wall_s=tangent_wall+linear_wall, full_map_evaluations_per_second=1./map_wall,
        candidate_forecast_s=candidate_forecast,
        cost_forecast_formula='12*measured_map_s + 2*measured_tangent_linear_s + 3 s IO margin',
        force_root_cost_measured=False, physical_full_map_evaluations=2,
        initial_free_tangent_condition_infinity_estimate=condition, initial_linear_relative_residual=residual,
        scalar_batch_differences=differences, dt_s=DTS['coarse'],
        stable_timestep_s='NOT_ESTABLISHED_BY_UNSOLVED_COST_MAP',
        registered_maximum_increment_m=.1*trial.p['dx'],
        evaluated_maximum_increment_m=float(np.linalg.norm(result['displacement'], axis=1).max()),
        particle_count=len(pa.x), free_dof=int(free.sum()),
        operator_reuse_audit=model.operator_reuse_audit())


def probe_stage(base, directory, protocol, report, rows, screens, deadline):
    cost = cost_receipt(base)
    initial = ROOT/cost['initial_checkpoint']
    require(core.sha(initial) == cost['initial_sha256'], 'independent parent changed')
    coarse_rows = []
    per_root = float(protocol['measured_candidate_forecast_s'])
    for role in ('coarse', 'fine'):
        folder = directory/role; folder.mkdir()
        pa, state = load(initial)
        require(state['time_s'] == 0. and state['step'] == 0, 'both paths require the identical initial state')
        parent = core.sha(initial)
        for number in range(1, protocol['counts'][role]+1):
            require(time.perf_counter()+1.5*per_root+20 < deadline, 'remaining measured stage budget exhausted')
            step = folder/f'step_{number:03d}'; step.mkdir()
            item = complete_candidate(pa, state, DTS[role], step, deadline, per_root, protocol['source_sha256'])
            pa, state = item['pa'], item['state']
            require(not state['accepted'] and not state['restart_allowed'], 'diagnostic publication forbidden')
            row = dict(role=role, step=number, time_s=state['time_s'], dt_s=DTS[role],
                checkpoint=rel(item['checkpoint']), sha256=core.sha(item['checkpoint']), parent_sha256=parent,
                wall_seconds=item['wall_seconds'], metrics=item['metrics'], accepted=False, restart_allowed=False)
            rows.append(row); parent = row['sha256']
            if role == 'coarse':
                coarse_rows.append(row)
            elif number % 2 == 0:
                coarse_row = coarse_rows[number//2-1]
                cp, cs = load(ROOT/coarse_row['checkpoint'])
                require(abs(cs['time_s']-state['time_s']) <= 1e-14, 'coarse/fine endpoint time mismatch')
                global_screen = post.compare(pa, cp, state['parameters'])
                local_screen = event.local_compare(pa, cp, sorted(int(v) for v in pa.particle_id))
                same_crossings = state['bore_observation']['crossing_ids'] == cs['bore_observation']['crossing_ids']
                screen = dict(pair=number//2, time_s=state['time_s'],
                    scope='COMMON_INITIAL_CONTINUOUS_TWO_RESOLUTION_ALL_PARTICLE_IDS',
                    coarse_checkpoint_sha256=coarse_row['sha256'], fine_checkpoint_sha256=parent,
                    global_screen=global_screen, local_screen=local_screen, crossing_ids_match=same_crossings,
                    passed=bool(global_screen['pass'] and local_screen['passed'] and same_crossings), accepted=False)
                screens.append(screen); core.write(directory/f'pair_{number//2:03d}_screen.json', screen)
                require(screen['passed'], 'common-initial global/local/passage time screen failed')
            core.write(directory/'candidate_history.json', rows)
            per_root = max(per_root, item['wall_seconds'])
            report.update(candidate_steps=len(rows), current_role=role, candidate_time_s=state['time_s'])
            core.write(directory/'manifest.json', report)
            verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
            print(json.dumps(dict(role=role, step=number, time_s=state['time_s'],
                wall_seconds=item['wall_seconds'], accepted=False)), flush=True)
        require(abs(state['time_s']-protocol['end_time_s']) <= 1e-14, 'requested diagnostic endpoint not reached')
    report.update(status=PROBE_PASS, candidate_end_checkpoint=row['checkpoint'],
        candidate_end_sha256=row['sha256'], temporal_screens_all_pass=bool(screens) and all(v['passed'] for v in screens),
        diagnostic_endpoint_s=state['time_s'], diagnostic_steps_per_second=len(rows)/sum(r['wall_seconds'] for r in rows),
        interval_energy=state['interval_energy'], crossing_ids=state['bore_observation']['crossing_ids'],
        stable_timestep_s='ONLY_THIS_SHORT_WINDOW_SCREENED_NOT_GENERAL_STABILITY',
        full_machine_program_executed=False)


def worker(base, stage):
    directory = base/stage; protocol = read(directory/'protocol.json')
    started = time.perf_counter(); rows, screens = [], []; code = 2
    report = dict(protocol, status='RUNNING_INDEPENDENT_PROCESS_FRAMEWORK_DIAGNOSTIC')
    def timeout():
        try:
            core.write(directory/'HARD_TIMEOUT.json', dict(status='TIMEOUT_NOT_ACCEPTED', completed_candidates=len(rows)))
        finally:
            os._exit(124)
    timer = threading.Timer(protocol['internal_limit_s'], timeout); timer.daemon = True; timer.start()
    deadline = started+protocol['internal_limit_s']-15.
    try:
        require(all(os.environ.get(k) == '1' for k in THREADS), 'single-thread CPU environment required')
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
        require(source_inputs(ROOT/protocol['configuration'], ROOT/protocol['component_receipt']) ==
                protocol['source_sha256'], 'complete source binding differs')
        core.write(directory/'manifest.json', report); core.write(directory/'candidate_history.json', rows)
        if stage == 'cost':
            cost_stage(base, directory, protocol, report)
        else:
            probe_stage(base, directory, protocol, report, rows, screens, deadline)
        verify(protocol['source_sha256']); verify(protocol['reference_sha256']); code = 0
    except Exception as error:
        report.update(status='FAIL_PROCESS_FRAMEWORK_DIAGNOSTIC_NOT_ACCEPTED', error=repr(error),
                      traceback=traceback.format_exc(), full_machine_program_executed=False)
    finally:
        timer.cancel()
        core.write(directory/'candidate_history.json', rows)
        report.update(wall_seconds=time.perf_counter()-started, candidate_steps=len(rows), pair_screens=screens,
            peak_working_set_bytes=core.memory_bytes(), accepted=False, restart_allowed=False,
            accepted_steps=0, accepted_time_advanced_s=0.)
        report['artifacts_sha256'] = {rel(p): core.sha(p) for p in directory.rglob('*') if p.is_file()
            and p.name not in ('manifest.json', 'run.log', 'supervision.json')}
        if stage == 'cost' and (base/'initial').exists():
            report['artifacts_sha256'].update({rel(p): core.sha(p) for p in (base/'initial').iterdir() if p.is_file()})
        core.write(directory/'manifest.json', report)
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('cost', 'probe'), required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--checks', type=Path, default=DEFAULT_CHECKS)
    parser.add_argument('--config', type=Path, default=model.DEFAULT_CONFIG)
    parser.add_argument('--end-time', type=float, default=.003)
    parser.add_argument('--preflight-only', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    base = (ROOT/args.output).resolve(); config = (ROOT/args.config).resolve(); checks = (ROOT/args.checks).resolve()
    require(base.is_relative_to(ROOT/'outputs') and base != ROOT/'outputs', 'new project output directory required')
    if args.worker:
        return worker(base, args.stage)
    protocol = preflight(base, args.stage, config, checks, args.end_time)
    if args.preflight_only:
        print(json.dumps({k: protocol[k] for k in ('stage', 'particle_count', 'counts', 'forecast_s', 'cost_estimate')}, indent=2))
        return 0
    directory = base/args.stage; directory.mkdir(parents=True, exist_ok=False)
    write_new(directory/'protocol.json', protocol)
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--stage', args.stage,
               '--output', rel(base), '--worker']
    began = time.perf_counter(); code = 124
    environment = dict(os.environ, **{k: '1' for k in THREADS})
    environment['PYTHONPATH'] = ''
    with (directory/'run.log').open('x', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=ROOT, env=environment, stdout=stream,
                                  stderr=subprocess.STDOUT, timeout=protocol['external_limit_s']).returncode
        except subprocess.TimeoutExpired:
            stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    source_error = None
    try:
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
    except Exception as error:
        source_error = repr(error); code = 2
    manifest = directory/'manifest.json'
    write_new(directory/'supervision.json', dict(command=command, exit_code=code,
        wall_seconds=time.perf_counter()-began, source_inputs_still_match=source_error is None,
        source_error=source_error, manifest_sha256=core.sha(manifest) if manifest.exists() else None,
        log_sha256=core.sha(directory/'run.log'), internal_limit_s=540., external_limit_s=600.,
        accepted=False, restart_allowed=False, accepted_time_advanced_s=0.))
    print(json.dumps(dict(stage=args.stage, exit_code=code,
        status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST', output=rel(directory))))
    return code


if __name__ == '__main__':
    raise SystemExit(main())

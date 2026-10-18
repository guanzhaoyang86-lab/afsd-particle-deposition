"""Bounded, adaptive, screened diagnostic continuation of the deposition setup.

Every macro interval compares one coarse solve with two half solves from the
same unchanged parent. Only the screened fine branch may seed a later stage.
Formal acceptance/restart flags remain false. No old global module is patched.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import afsd_deposition_remote_setup_20261007 as model
import run_afsd_process_framework_20261007 as shared
from visible_support_remote_domain_20261007 import audit_visible_support_path as remote_audit_visible_support_path
from visible_support_remote_events_20261007 import audit_visibility_events as remote_audit_visibility_events

core, post, event, batch, solver, old = (shared.core, shared.post, shared.event,
                                      shared.batch, shared.solver, shared.old)
ROOT, require = model.ROOT, model.require
CONTRACT = ROOT/'docs/DEPOSITION_STAGE_RUNNER_CONTRACT_20261007.md'
DEFAULT_CHECKS = ROOT/'outputs/deposition_campaign_20261007/components_remote_v1.json'
REFERENCE_COST = ROOT/'outputs/process_framework_20261007/run_v1/cost/manifest.json'
COST_PASS = 'PASS_DEPOSITION_COST_MAP_NOT_PHYSICAL_VALIDATION'
STAGE_PASS = 'PASS_SCREENED_DEPOSITION_STAGE_DIAGNOSTIC_CONTINUATION'
SHOULDER_PASS = 'PASS_FIRST_NEW_MATERIAL_UNDER_SHOULDER_DIAGNOSTIC_CONTINUATION'
PREFIX_PASS = 'PASS_SCREENED_PREFIX_BUDGET_BOUND_NOT_TARGET'
read, load, rel, verify = shared.read, shared.load, shared.rel, shared.verify


class GeometryPathStop(RuntimeError):
    """Observed or unresolved whole-path geometry is not auto-refined here."""


class RefinableCandidateError(RuntimeError):
    """Failed full physical/numerical candidate; may retry from unchanged parent."""


class RefinementExhausted(RuntimeError):
    pass


class StageBudgetStop(RuntimeError):
    def __init__(self, reason, next_dt):
        self.reason, self.next_dt = reason, next_dt
        super().__init__(reason)


def digest(value):
    return hashlib.sha256(json.dumps(core.clean(value), sort_keys=True,
                                    allow_nan=False).encode()).hexdigest()


def parent_fingerprint(pa, state):
    values = {}
    for name, prop in pa.properties.items():
        a = np.ascontiguousarray(prop.get_npy_array())
        values[name] = hashlib.sha256(str((a.shape, str(a.dtype))).encode()+a.tobytes()).hexdigest()
    return digest(dict(properties=values, state=core.clean(state)))


def bounded_interval(time_s, desired_dt, end_time_s, phase_end_s, minimum_dt):
    values = (time_s, desired_dt, end_time_s, phase_end_s, minimum_dt)
    if not all(math.isfinite(float(v)) for v in values) or min(desired_dt, minimum_dt) <= 0:
        raise ValueError('finite positive time increments required')
    dt = min(desired_dt, end_time_s-time_s, phase_end_s-time_s)
    if dt < minimum_dt:
        raise ValueError('remaining commanded interval below minimum dt or outside phase')
    return float(dt)


def adaptive_pair(dt, minimum_dt, maximum_halvings, attempt, record):
    """Controller only. ``attempt`` must close over the SAME immutable parent.

    Geometry exceptions propagate immediately. Failed temporal screens and
    explicitly refinable candidate errors halve dt without changing thresholds.
    """
    if not isinstance(maximum_halvings, int) or not 0 <= maximum_halvings <= 12:
        raise ValueError('maximum_halvings must be0..12')
    if not (math.isfinite(dt) and math.isfinite(minimum_dt) and dt >= minimum_dt > 0):
        raise ValueError('finite dt at or above a positive minimum required')
    trial_dt = float(dt)
    for level in range(maximum_halvings+1):
        try:
            outcome = attempt(trial_dt, level)
            passed = bool(outcome['screen']['passed'])
            row = dict(level=level, dt_s=trial_dt, passed=passed,
                       status='PASS_PAIR' if passed else 'FAIL_ORIGINAL_TEMPORAL_THRESHOLDS')
        except RefinableCandidateError as error:
            outcome = None
            row = dict(level=level, dt_s=trial_dt, passed=False,
                       status='FAIL_CANDIDATE_NO_STATE_COMMIT', error=str(error))
        record(row)
        if row['passed']:
            return outcome, trial_dt
        if level == maximum_halvings or trial_dt/2 < minimum_dt:
            raise RefinementExhausted('same-parent refinement exhausted; failed candidates remain uncommitted')
        trial_dt /= 2
    raise AssertionError('unreachable adaptive controller state')


def screen_utilization(global_result, local_result, parameters):
    """Largest fraction of an UNCHANGED comparison limit actually consumed."""
    ratios = [global_result['position_max_m']/(.01*parameters['dx']),
        global_result['temperature_max_K']/.5, global_result['stress_relative']/.1,
        global_result['plastic_relative']/.2]
    for row in local_result['rows']:
        ratios.extend((row['position_error_m']/local_result['position_limit_m'],
            row['temperature_error_K']/local_result['temperature_limit_K'],
            row['stress_relative']/local_result['stress_relative_limit'],
            row['plastic_relative']/local_result['plastic_relative_limit']))
    require(all(math.isfinite(v) and v >= 0 for v in ratios), 'nonfinite temporal utilization')
    return max(ratios)


def next_step_size(passed_dt, utilization, maximum_dt, growth_factor):
    if not (0 < passed_dt <= maximum_dt and 1 <= growth_factor <= 2
            and math.isfinite(utilization) and 0 <= utilization <= 1):
        raise ValueError('valid passed interval and bounded growth required')
    return min(maximum_dt, passed_dt*growth_factor) if utilization <= .25 else passed_dt


def source_inputs(config, checks):
    receipt = read(checks)
    require(str(receipt.get('status', '')).startswith('PASS') and receipt.get('checks')
        and all(row.get('status') == 'PASS' and row.get('exit_code', 0) == 0 for row in receipt['checks']),
        'executed deposition component tests required')
    require(receipt.get('operator_reuse_audit', {}).get('status') ==
            'PASS_UNCHANGED_INHERITED_THERMOMECHANICAL_OPERATORS', 'new setup operator reuse receipt required')
    mapping = dict(model.source_hashes(config))
    for key in ('source_sha256', 'output_sha256', 'evidence_sha256'):
        verify(receipt.get(key, {}))
        for name, sha in receipt.get(key, {}).items():
            require(name not in mapping or mapping[name] == sha, 'conflicting component source binding')
            mapping[name] = sha
    required_tests = ('tests/review/test_afsd_deposition_setup_20261007.py',
                      'scripts/test_afsd_deposition_stage_20261007.py',
                      'scripts/test_afsd_deposition_evidence_20261007.py',
                      'scripts/test_visible_support_remote_domain_20261007.py',
                      'scripts/test_afsd_deposition_remote_adaptation_20261007.py')
    for key in required_tests:
        require(receipt.get('source_sha256', {}).get(key) == core.sha(ROOT/key),
                'component receipt must bind '+key)
    files = {Path(__file__).resolve(), CONTRACT, Path(config).resolve(), Path(checks).resolve()}
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT/'scripts'):
                files.add(path)
    for path in files:
        mapping[rel(path)] = core.sha(path)
    verify(mapping); model.operator_reuse_audit()
    return mapping


def checked_receipt(folder, statuses, frozen=None):
    meta, supervision = read(folder/'manifest.json'), read(folder/'supervision.json')
    require(meta['status'] in statuses and supervision['exit_code'] == 0
        and supervision['source_inputs_still_match']
        and supervision['manifest_sha256'] == core.sha(folder/'manifest.json')
        and supervision['log_sha256'] == core.sha(folder/'run.log'), 'stage receipt absent, failed or changed')
    verify(meta['source_sha256']); verify(meta['artifacts_sha256'])
    if frozen is not None:
        require(meta['source_sha256'] == frozen, 'new implementation requires a new independent cost stage')
    require(meta['classification'] == model.CLASSIFICATION and not meta['accepted']
        and not meta['restart_allowed'] and meta['accepted_time_advanced_s'] == 0.,
        'only unchanged deposition diagnostic states may continue')
    return meta


def bind_receipt(folder, meta, references):
    references.update(meta['artifacts_sha256'])
    for name in ('manifest.json', 'supervision.json', 'run.log'):
        path = folder/name; references[rel(path)] = core.sha(path)


def continuation_parent(cost_folder, parent_folder, frozen):
    cost = checked_receipt(cost_folder, {COST_PASS}, frozen)
    references = {}; bind_receipt(cost_folder, cost, references)
    initial = ROOT/cost['initial_checkpoint']
    require(core.sha(initial) == cost['initial_sha256'], 'independent initial checkpoint changed')
    target = initial; parent_manifest = None; visited = set()
    current = parent_folder
    while current is not None:
        current = current.resolve()
        require(current not in visited, 'cycle in diagnostic continuation history')
        visited.add(current)
        meta = checked_receipt(current, {STAGE_PASS, SHOULDER_PASS, PREFIX_PASS}, frozen)
        require(meta['screened_continuation_authorized'] and meta['screened_history']
            and meta['initial_sha256'] == cost['initial_sha256']
            and meta['initial_checkpoint'] == cost['initial_checkpoint'], 'stage origin or screening authorization differs')
        require(meta['screened_continuation_checkpoint'] == meta['screened_history'][-1]['checkpoint']
            and meta['screened_continuation_sha256'] == meta['screened_history'][-1]['sha256'],
            'only the last screened fine state may continue')
        require(len(meta['screened_history']) == 2*len(meta['pair_screens'])
            and all(row['passed'] for row in meta['pair_screens']), 'all selected fine pairs must pass')
        for row in meta['screened_history']:
            require(core.sha(ROOT/row['checkpoint']) == row['sha256']
                and core.sha(ROOT/row['observer_path']) == row['observer_sha256']
                and core.sha(ROOT/row['screen_path']) == row['screen_sha256'], 'screened history evidence changed')
        for offset in range(0, len(meta['screened_history']), 2):
            first, last = meta['screened_history'][offset:offset+2]
            require(first['half_index'] == 1 and last['half_index'] == 2
                and first['macro_index'] == last['macro_index']
                and first['screen_path'] == last['screen_path'], 'screened pair row identity mismatch')
            actual_screen = read(ROOT/first['screen_path'])
            require(actual_screen == meta['pair_screens'][offset//2]
                and actual_screen['passed'] and actual_screen['global_screen']['pass']
                and actual_screen['local_screen']['passed'] and actual_screen['crossing_ids_match']
                and actual_screen['fine_checkpoint_sha256'] == [first['sha256'], last['sha256']]
                and actual_screen['time_s'] == last['time_s'], 'PASS screen does not bind these exact fine candidates')
        if parent_manifest is None:
            parent_manifest = meta; target = ROOT/meta['screened_continuation_checkpoint']
        bind_receipt(current, meta, references)
        prior = meta.get('parent_stage_manifest')
        if prior is None:
            require(meta['parent_checkpoint_sha256'] == cost['initial_sha256'], 'first stage must start at independent initial state')
            current = None
        else:
            prior_path = ROOT/prior
            require(core.sha(prior_path) == meta['parent_stage_manifest_sha256'], 'parent stage manifest changed')
            prior_meta = read(prior_path)
            require(meta['parent_checkpoint_sha256'] == prior_meta['screened_continuation_sha256'], 'stage chain skipped or mixed parents')
            current = prior_path.parent
    pa, state = load(target); model.validate_state(state)
    require(state['classification'] == model.CLASSIFICATION and state['diagnostic_only']
            and not state['accepted'] and not state['restart_allowed'], 'forbidden formal/other-model parent')
    if parent_manifest is not None:
        require(state['time_s'] == parent_manifest['screened_end_time_s']
            and state['step'] == parent_manifest['screened_end_step'], 'parent checkpoint endpoint differs')
    return cost, target, pa, state, references, parent_manifest


def preflight(args):
    output = (ROOT/args.output).resolve()
    require(output.is_relative_to(ROOT/'outputs') and output != ROOT/'outputs' and not output.exists(), 'new project output required')
    config, checks = (ROOT/args.config).resolve(), (ROOT/args.checks).resolve()
    frozen = source_inputs(config, checks)
    _, _, program, domain = model.load_configuration(config)
    references = {}; parent_info = {}
    if args.stage == 'cost':
        prior = read(REFERENCE_COST); references[rel(REFERENCE_COST)] = core.sha(REFERENCE_COST)
        require(str(prior['status']).startswith('PASS'), 'prior measured cost required')
        ratio = domain['particle_count']/prior['particle_count']
        forecast = 60.+1.5*(3*prior['full_map_wall_s']*max(1., ratio**2)
            +2*prior['tangent_linear_wall_s']*max(1., ratio**3))
        cost_info = dict(kind='PRIOR_MEASURED_COST_N2_MAP_N3_DENSE_LINEAR_SCALING', particle_count_ratio=ratio)
    else:
        require(args.cost is not None and args.end_time is not None, 'stage requires --cost and explicit global --end-time')
        require(1 <= args.max_macrosteps <= 32 and 3 <= args.max_candidate_attempts <= 192 and 0 <= args.max_halvings <= 12
            and math.isfinite(args.max_dt) and math.isfinite(args.min_dt)
            and 0 < args.min_dt <= args.max_dt <= .02 and 1 <= args.growth_factor <= 2,
            'bounded steps/refinements and positive dt<=20ms with growth<=2 required')
        cost_folder = (ROOT/args.cost).resolve(); parent_folder = (ROOT/args.parent).resolve() if args.parent else None
        cost, target, pa, state, references, parent_meta = continuation_parent(cost_folder, parent_folder, frozen)
        desired_dt = args.dt if args.dt is not None else (parent_meta.get('next_proposed_macro_dt_s', .0005) if parent_meta else .0005)
        require(math.isfinite(desired_dt) and args.min_dt <= desired_dt <= args.max_dt,
                'initial dt must lie within explicit minimum/maximum')
        program = model.validate_state(state); phase = program.phase(state['time_s'])
        require(math.isfinite(args.end_time) and state['time_s'] < args.end_time <= program.end_s,
                'global end-time must advance within the finite program')
        endpoint = min(args.end_time, phase['end'])
        bounded_interval(state['time_s'], desired_dt, endpoint, phase['end'], args.min_dt)
        measured = max(float(cost['candidate_forecast_s']),
                       float(parent_meta.get('maximum_candidate_wall_s', 0.)) if parent_meta else 0.)
        forecast = 30.+1.5*args.max_candidate_attempts*measured
        cost_info = dict(kind='CURRENT_MEASURED_COST_AND_PRIOR_SCREENED_CANDIDATE_MAX', candidate_forecast_s=measured,
            formula='30 + 1.5 * max_candidate_attempts * per_candidate_forecast; every actual candidate including retries consumes quota')
        parent_info = dict(cost_directory=rel(cost_folder), initial_checkpoint=cost['initial_checkpoint'], initial_sha256=cost['initial_sha256'],
            parent_checkpoint=rel(target), parent_checkpoint_sha256=core.sha(target), parent_time_s=state['time_s'], parent_step=state['step'],
            parent_stage_manifest=rel(parent_folder/'manifest.json') if parent_folder else None,
            parent_stage_manifest_sha256=core.sha(parent_folder/'manifest.json') if parent_folder else None,
            parent_state_fingerprint=parent_fingerprint(pa, state), parent_machine_history_sha256=digest(state['machine_history']),
            parent_energy_history_sha256=digest(state['interval_energy']['records']),
            requested_end_time_s=args.end_time, stage_end_time_s=endpoint, phase_name=phase['name'], phase_end_time_s=phase['end'],
            dt_s=desired_dt, min_dt_s=args.min_dt, max_dt_s=args.max_dt, growth_factor=args.growth_factor,
            max_halvings=args.max_halvings, max_macrosteps=args.max_macrosteps,
            max_candidate_attempts=args.max_candidate_attempts,
            stop_on_first_shoulder=bool(args.stop_on_first_shoulder), measured_candidate_forecast_s=measured,
            traverse_spreading_gate_required=bool(phase['target'][2] > 0 or phase['initial'][2] > 0))
        if parent_info['traverse_spreading_gate_required']:
            require(parent_meta is not None, 'traverse requires actual screened material spreading observations')
            observer = read(ROOT/parent_meta['screened_history'][-1]['observer_path'])
            parent_info['observed_spreading_before_traverse_kg'] = observer['supplied_material']['crossed_under_annular_shoulder']['mass_kg']
            supplied_ids = observer['supplied_material']['crossed_under_annular_shoulder']['particle_ids']
            initial_pa, _ = load(ROOT/cost['initial_checkpoint'])
            require(np.array_equal(initial_pa.particle_id, pa.particle_id), 'initial/current material ID order differs')
            plastic = np.isin(pa.particle_id, supplied_ids) & (pa.ep > initial_pa.ep)
            parent_info['plastic_supplied_shoulder_ids'] = pa.particle_id[plastic].astype(int).tolist()
            parent_info['plastic_supplied_shoulder_mass_kg'] = math.fsum(pa.m[plastic].tolist())
            parent_info['plastic_spreading_gate_scope'] = 'CURRENT_ELIGIBLE_CROSSED_ANNULAR_IDS_WITH_POSITIVE_ACCUMULATED_EP_RELATIVE_TO_INDEPENDENT_INITIAL'
    require(forecast <= 450., 'predicted bounded stage exceeds450s; choose fewer macrosteps')
    return dict(schema='SCREENED_DEPOSITION_STAGE_20261007', stage=args.stage, classification=model.CLASSIFICATION,
        configuration=rel(config), component_receipt=rel(checks), source_sha256=frozen, reference_sha256=references,
        domain_inventory_plan=domain, particle_count=domain['particle_count'], cost_estimate=cost_info,
        forecast_s=forecast, forecast_limit_s=450., internal_limit_s=540., external_limit_s=600.,
        accepted=False, restart_allowed=False, diagnostic_only=True, accepted_steps=0, accepted_time_advanced_s=0.,
        **parent_info)


def scalar_equivalence(trial, result, directory):
    scalar = model.ProcessTrial(trial.pa, trial.state, result['dt_s']); p = trial.p
    scalar.contact = old.ChainRuleDiscreteContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2, normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'], friction=p['friction'], base_z=-p['substrate_depth'])
    sr = scalar.evaluate(result['displacement'], result['dt_s'])
    differences = batch.candidate_equivalence(dict(sr, energy_comparison=core.energy_diagnostics(scalar, sr)),
        dict(result, energy_comparison=core.energy_diagnostics(trial, result)))
    require(sr['graph_fingerprint'] == result['graph_fingerprint'] and np.array_equal(sr['g'], result['g'])
        and np.array_equal(sr['thermal']['pair_conductance_W_K'], result['thermal']['pair_conductance_W_K']),
        'scalar/batch exact common graph and heat differ')
    core.write(directory/'scalar_batch_equivalence.json', dict(status='PASS_COMPLETE_MAP_EQUIVALENCE', differences=differences))
    return differences


def complete_candidate(pa, state, dt, directory, deadline, per_candidate, frozen):
    """Explicit new-model candidate; no rebinding of any inherited global."""
    began = time.perf_counter(); trial = model.ProcessTrial(pa, state, dt)
    attempt = directory/'broyden'; attempt.mkdir(); guess = batch.guess_for(trial, dt)
    try:
        try:
            result, metrics = solver.solve(trial, dt, guess, attempt, deadline, 24)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            require(time.perf_counter()+1.5*per_candidate+20 < deadline, 'no budget for one unchanged Newton retry')
            attempt = directory/'newton_retry'; attempt.mkdir()
            result, metrics = core.solve(trial, dt, guess, attempt, deadline, 24)
            metrics['fallback_reason'] = repr(error)
        core.journal(directory, trial, result, 'SOLVED_UNCOMMITTED_candidate'); old.save_graph(directory, trial, result)
        nextpa, nextstate, fresh = model.advance_candidate(trial, result)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        raise RefinableCandidateError(str(error)) from error
    require(np.array_equal(pa.particle_id, nextpa.particle_id) and np.array_equal(pa.m, nextpa.m)
        and model.identity_digest(pa) == model.identity_digest(nextpa), 'immutable material inventory changed')
    target, _ = event.save_candidate(directory, trial, nextpa, nextstate, fresh, frozen)
    old.save_graph(directory, trial, fresh)
    observer = shared.observe_process(pa, nextstate, fresh, trial.motion, trial.p)
    core.write(directory/'observer.json', observer); core.write(directory/'motion.json', trial.motion.to_dict())
    core.write(directory/'transport.json', nextstate['bore_observation'])
    scalar_equivalence(trial, fresh, directory)
    geometry, events = shared.inherited.geometry_audit(trial, fresh)
    core.write(directory/'geometry_audit.json', geometry); core.write(directory/'sampled_events.json', events)
    pose0, pose1 = trial.motion.pose(trial.motion.t0), trial.motion.pose(trial.motion.t1)
    shell = old.audit_contact_shell_finite_enclosure(trial.contact, trial.x0, fresh['position'], pose0, pose1)
    core.write(directory/'continuous_shell_audit.json', shell)
    path = remote_audit_visible_support_path(trial.x0, fresh['position'], fresh['i'], fresh['j'],
                                             pose0, pose1, trial.p, pa.particle_id)
    core.write(directory/'visible_path_audit.json', path)
    if (geometry['transient_contacts_missed_by_endpoints'] or set(geometry['new_non_bore_features'])-{'SHOULDER', 'TOP'}
            or shell['status'] != event.CERTIFIED or path['status'] != old.CERTIFIED):
        try:
            routing = remote_audit_visibility_events(trial.x0, fresh['position'], fresh['i'], fresh['j'],
                                                      pose0, pose1, trial.p, pa.particle_id)
        except Exception as error:
            routing = dict(status='UNRESOLVED_EVENT_DIAGNOSTIC', error=repr(error), requires_resolve=True, acceptance_allowed=False)
        core.write(directory/'visibility_event_routing.json', routing)
        raise GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')
    core.write(directory/'candidate_qualification.json', dict(status='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN',
        accepted=False, restart_allowed=False, metrics=metrics))
    core.vtk(directory/'fields.vtk', nextpa, trial.p, trial.bulk, nextstate['time_s'], model.CLASSIFICATION)
    return dict(pa=nextpa, state=nextstate, checkpoint=target, observer_path=directory/'observer.json', observer=observer,
                metrics=metrics, wall_seconds=time.perf_counter()-began)


def cost_stage(output, protocol, report):
    pa, state = model.independent_initial_state(ROOT/protocol['configuration'])
    initial = output/'INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
    core.checkpoint(initial, pa, state); core.verify_snapshot(initial, pa, state)
    core.write(output/'initial_state.json', state)
    trial = model.ProcessTrial(pa, state, .0005)
    started = time.perf_counter(); result = trial.evaluate(trial.prescribed(.0005), .0005); map_wall = time.perf_counter()-started
    core.journal(output, trial, result, 'UNSOLVED_COST_ONLY'); old.save_graph(output, trial, result)
    started = time.perf_counter(); tangent = trial.tangent(result); tangent_wall = time.perf_counter()-started
    require(np.isfinite(tangent).all(), 'nonfinite initial tangent')
    free = state['free'].ravel(); matrix = tangent[np.ix_(free, free)]
    started = time.perf_counter(); inverse = np.linalg.inv(matrix)
    direction = np.linalg.solve(matrix, -result['force'][state['free']]); linear_wall = time.perf_counter()-started
    condition = float(np.linalg.norm(matrix, np.inf)*np.linalg.norm(inverse, np.inf))
    require(np.isfinite(direction).all() and condition <= 1e12, 'initial tangent near singular; no regularization')
    rhs = -result['force'][state['free']]
    linear_residual = float(np.linalg.norm(matrix@direction-rhs, np.inf)/max(1., np.linalg.norm(rhs, np.inf)))
    require(linear_residual <= 1e-8, 'initial linear direction residual failed')
    differences = scalar_equivalence(trial, result, output)
    report.update(status=COST_PASS, initial_checkpoint=rel(initial), initial_sha256=core.sha(initial),
        full_map_wall_s=map_wall, tangent_wall_s=tangent_wall, linear_wall_s=linear_wall,
        tangent_linear_wall_s=tangent_wall+linear_wall, free_dof=int(free.sum()),
        candidate_forecast_s=12*map_wall+2*(tangent_wall+linear_wall)+3.,
        initial_free_tangent_condition_infinity_estimate=condition, initial_linear_relative_residual=linear_residual,
        dt_s=.0005, force_root_cost_measured=False, scalar_batch_differences=differences,
        stable_timestep_s='NOT_ESTABLISHED_BY_UNSOLVED_MAP', registered_maximum_increment_m=.1*trial.p['dx'],
        evaluated_maximum_increment_m=float(np.linalg.norm(result['displacement'], axis=1).max()),
        physical_full_map_evaluations=2, operator_reuse_audit=model.operator_reuse_audit())


def stage_work(output, protocol, report, screened, attempts, screens, deadline):
    pa, state = load(ROOT/protocol['parent_checkpoint'])
    require(parent_fingerprint(pa, state) == protocol['parent_state_fingerprint'], 'full parent state changed')
    program = model.validate_state(state)
    if protocol['traverse_spreading_gate_required'] and (protocol.get('observed_spreading_before_traverse_kg', 0.) <= 0.
            or protocol.get('plastic_supplied_shoulder_mass_kg', 0.) <= 0.):
        raise GeometryPathStop('traverse blocked: no screened new supplied material with accumulated plastic response under annular shoulder')
    prefix_machine = list(state['machine_history']); prefix_energy = list(state['interval_energy']['records'])
    per_candidate = protocol['measured_candidate_forecast_s']; current_dt = protocol['dt_s']; reason = None
    maximum_wall = 0.; candidate_attempt_count = 0
    def call_candidate(candidate_pa, candidate_state, dt, directory):
        nonlocal maximum_wall, per_candidate, candidate_attempt_count
        require(candidate_attempt_count < protocol['max_candidate_attempts'], 'candidate quota unexpectedly exhausted')
        candidate_attempt_count += 1
        report['candidate_attempt_count'] = candidate_attempt_count
        began = time.perf_counter()
        try:
            return complete_candidate(candidate_pa, candidate_state, dt, directory,
                                      deadline, per_candidate, protocol['source_sha256'])
        finally:
            maximum_wall = max(maximum_wall, time.perf_counter()-began)
            per_candidate = max(per_candidate, maximum_wall)
    for macro in range(1, protocol['max_macrosteps']+1):
        if state['time_s'] >= protocol['stage_end_time_s']-64*np.finfo(float).eps*max(1., protocol['stage_end_time_s']):
            reason = 'REQUESTED_ENDPOINT_OR_PHASE_BOUNDARY'; break
        if time.perf_counter()+4.5*per_candidate+20 >= deadline:
            require(bool(screened), 'no budget for one complete comparison')
            reason = 'REMAINING_MEASURED_WALL_BUDGET'; break
        if candidate_attempt_count+3 > protocol['max_candidate_attempts']:
            require(bool(screened), 'candidate quota cannot cover a complete first comparison')
            reason = 'DECLARED_TOTAL_CANDIDATE_QUOTA'; break
        dt = bounded_interval(state['time_s'], current_dt, protocol['stage_end_time_s'],
                              program.phase(state['time_s'])['end'], protocol['min_dt_s'])
        parent_pa, parent_state = pa, state
        fingerprint = parent_fingerprint(parent_pa, parent_state)
        macro_folder = output/f'macro_{macro:04d}'; macro_folder.mkdir()
        def attempt(trial_dt, level):
            if candidate_attempt_count+3 > protocol['max_candidate_attempts']:
                raise StageBudgetStop('DECLARED_TOTAL_CANDIDATE_QUOTA', trial_dt)
            if time.perf_counter()+4.5*per_candidate+20 >= deadline:
                raise StageBudgetStop('REMAINING_MEASURED_WALL_BUDGET', trial_dt)
            require(parent_fingerprint(parent_pa, parent_state) == fingerprint, 'same-parent retry state mutated')
            folder = macro_folder/f'attempt_{level:02d}'; folder.mkdir()
            manifest = dict(macro_index=macro, refinement_level=level, dt_s=trial_dt,
                parent_fingerprint=fingerprint, parent_time_s=parent_state['time_s'], status='RUNNING_PAIR', accepted=False)
            try:
                items = []
                for role, candidate_pa, candidate_state, step_dt in (
                        ('coarse', parent_pa, parent_state, trial_dt), ('fine_1', parent_pa, parent_state, trial_dt/2)):
                    directory = folder/role; directory.mkdir()
                    items.append(call_candidate(candidate_pa, candidate_state, step_dt, directory))
                directory = folder/'fine_2'; directory.mkdir()
                items.append(call_candidate(items[1]['pa'], items[1]['state'], trial_dt/2, directory))
                coarse, first, last = items
                require(parent_fingerprint(parent_pa, parent_state) == fingerprint, 'candidate modified its original parent')
                require(abs(last['state']['time_s']-coarse['state']['time_s']) <= 1e-13, 'pair endpoint time mismatch')
                global_result = post.compare(last['pa'], coarse['pa'], state['parameters'])
                local_result = event.local_compare(last['pa'], coarse['pa'], sorted(int(v) for v in last['pa'].particle_id))
                matching = last['state']['bore_observation']['crossing_ids'] == coarse['state']['bore_observation']['crossing_ids']
                screen = dict(macro_index=macro, refinement_level=level, parent_fingerprint=fingerprint,
                    time_s=last['state']['time_s'], dt_s=trial_dt, scope='ONE_COARSE_VERSUS_TWO_HALVES_IDENTICAL_FULL_PARENT_ALL_IDS',
                    global_screen=global_result, local_screen=local_result, crossing_ids_match=matching,
                    coarse_checkpoint_sha256=core.sha(coarse['checkpoint']),
                    fine_checkpoint_sha256=[core.sha(first['checkpoint']), core.sha(last['checkpoint'])],
                    passed=bool(global_result['pass'] and local_result['passed'] and matching), accepted=False)
                screen['maximum_original_tolerance_utilization'] = screen_utilization(global_result, local_result, state['parameters'])
                screen_path = folder/'pair_screen.json'; core.write(screen_path, screen)
                manifest.update(status='PASS_PAIR' if screen['passed'] else 'FAIL_ORIGINAL_TIME_SCREEN', screen_path=rel(screen_path))
                return dict(fine=[first, last], screen=screen, screen_path=screen_path)
            except Exception as error:
                manifest.update(status='FAIL_PAIR_NO_STATE_COMMIT', error=repr(error), traceback=traceback.format_exc())
                raise
            finally:
                require(parent_fingerprint(parent_pa, parent_state) == fingerprint, 'rejected trial changed parent')
                core.write(folder/'attempt.json', manifest)
        def record(row):
            attempts.append(dict(macro_index=macro, **row)); core.write(output/'attempt_history.json', attempts)
        try:
            outcome, used_dt = adaptive_pair(dt, protocol['min_dt_s'], protocol['max_halvings'], attempt, record)
        except StageBudgetStop as bound:
            require(bool(screened), 'budget exhausted before any screened prefix: '+bound.reason)
            reason, current_dt = bound.reason, bound.next_dt
            break
        current_dt = next_step_size(used_dt, outcome['screen']['maximum_original_tolerance_utilization'],
                                    protocol['max_dt_s'], protocol['growth_factor'])
        screen_path = outcome['screen_path']; screens.append(outcome['screen'])
        for half, item in enumerate(outcome['fine'], 1):
            s = item['state']; n = len(prefix_machine)
            require(s['machine_history'][:n] == prefix_machine and s['interval_energy']['records'][:n] == prefix_energy,
                    'full stage parent machine/energy prefix changed')
            screened.append(dict(checkpoint=rel(item['checkpoint']), sha256=core.sha(item['checkpoint']), time_s=s['time_s'],
                observer_path=rel(item['observer_path']), observer_sha256=core.sha(item['observer_path']),
                macro_index=macro, half_index=half, screen_path=rel(screen_path), screen_sha256=core.sha(screen_path),
                crossed_mass_kg=item['observer']['supplied_material']['eligible_crossed']['mass_kg'],
                under_shoulder_mass_kg=item['observer']['supplied_material']['crossed_under_annular_shoulder']['mass_kg'],
                wall_seconds=item['wall_seconds']))
        pa, state = outcome['fine'][-1]['pa'], outcome['fine'][-1]['state']
        core.write(output/'screened_history.json', screened)
        report.update(screened_history=screened, pair_screens=screens, screened_end_time_s=state['time_s'],
            screened_end_step=state['step'], maximum_candidate_wall_s=maximum_wall,
            last_passed_macro_dt_s=used_dt, next_proposed_macro_dt_s=current_dt)
        core.write(output/'manifest.json', report); verify(protocol['source_sha256'])
        print(json.dumps(dict(macro=macro, time_s=state['time_s'], dt_s=used_dt, screened_steps=len(screened),
            crossing_mass_kg=screened[-1]['crossed_mass_kg'], under_shoulder_mass_kg=screened[-1]['under_shoulder_mass_kg'])), flush=True)
        if protocol['stop_on_first_shoulder'] and any(r['under_shoulder_mass_kg'] > 0 for r in screened[-2:]):
            reason = 'FIRST_SCREENED_NEW_MATERIAL_UNDER_ANNULAR_SHOULDER'; break
    require(bool(screened), 'stage produced no screened diagnostic continuation')
    status = (SHOULDER_PASS if reason == 'FIRST_SCREENED_NEW_MATERIAL_UNDER_ANNULAR_SHOULDER' else
              PREFIX_PASS if reason in ('DECLARED_TOTAL_CANDIDATE_QUOTA', 'REMAINING_MEASURED_WALL_BUDGET') else STAGE_PASS)
    report.update(status=status,
        stop_reason=reason or 'DECLARED_MACROSTEP_BOUND', screened_continuation_authorized=True,
        continuation_scope='THIS_NEW_MODEL_AND_UNCHANGED_SOURCE_SCREENED_DIAGNOSTIC_ONLY',
        screened_continuation_checkpoint=screened[-1]['checkpoint'], screened_continuation_sha256=screened[-1]['sha256'],
        screened_end_time_s=state['time_s'], screened_end_step=state['step'], maximum_candidate_wall_s=maximum_wall,
        candidate_attempt_count=candidate_attempt_count, next_proposed_macro_dt_s=current_dt,
        stage_target_reached=abs(state['time_s']-protocol['stage_end_time_s']) <= 1e-13,
        failed_attempts_are_not_screened=True,
        parent_prefix_unchanged=True, formal_restart_allowed=False, full_program_completed=False)


def worker(output):
    protocol = read(output/'protocol.json'); began = time.perf_counter(); code = 2
    screened, attempts, screens = [], [], []
    report = dict(protocol, status='RUNNING_DEPOSITION_DIAGNOSTIC', screened_history=screened,
                  pair_screens=screens, screened_continuation_authorized=False)
    def timeout():
        try:
            core.write(output/'HARD_TIMEOUT.json', dict(status='TIMEOUT_NOT_ACCEPTED', screened_steps=len(screened)))
        finally:
            os._exit(124)
    timer = threading.Timer(540., timeout); timer.daemon = True; timer.start()
    try:
        require(all(os.environ.get(k) == '1' for k in shared.THREADS), 'single-thread CPU environment required')
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
        require(source_inputs(ROOT/protocol['configuration'], ROOT/protocol['component_receipt']) == protocol['source_sha256'],
                'full source closure changed')
        core.write(output/'manifest.json', report)
        if protocol['stage'] == 'cost':
            cost_stage(output, protocol, report)
        else:
            stage_work(output, protocol, report, screened, attempts, screens, began+525.)
        verify(protocol['source_sha256']); verify(protocol['reference_sha256']); code = 0
    except Exception as error:
        report.update(status='FAIL_DEPOSITION_STAGE_NOT_ACCEPTED', error=repr(error), traceback=traceback.format_exc(),
                      screened_continuation_authorized=False)
    finally:
        timer.cancel(); core.write(output/'screened_history.json', screened); core.write(output/'attempt_history.json', attempts)
        report.update(wall_seconds=time.perf_counter()-began, peak_working_set_bytes=core.memory_bytes(),
            screened_history=screened, pair_screens=screens, attempts=attempts,
            accepted=False, restart_allowed=False, accepted_steps=0, accepted_time_advanced_s=0.)
        report['artifacts_sha256'] = {rel(p): core.sha(p) for p in output.rglob('*') if p.is_file()
            and p.name not in ('manifest.json', 'run.log', 'supervision.json')}
        core.write(output/'manifest.json', report)
    return code


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=('cost', 'stage'), required=True); p.add_argument('--output', required=True)
    p.add_argument('--config', type=Path, default=model.DEFAULT_CONFIG); p.add_argument('--checks', type=Path, default=DEFAULT_CHECKS)
    p.add_argument('--cost', type=Path); p.add_argument('--parent', type=Path)
    p.add_argument('--end-time', type=float); p.add_argument('--dt', '--initial-dt', dest='dt', type=float)
    p.add_argument('--max-dt', type=float, default=.005); p.add_argument('--growth-factor', type=float, default=1.5)
    p.add_argument('--min-dt', type=float, default=1e-7); p.add_argument('--max-halvings', type=int, default=8)
    p.add_argument('--max-macrosteps', type=int, default=4); p.add_argument('--stop-on-first-shoulder', action='store_true')
    p.add_argument('--max-candidate-attempts', type=int, default=18)
    p.add_argument('--preflight-only', action='store_true'); p.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    return p


def main():
    args = parser().parse_args(); output = (ROOT/args.output).resolve()
    require(output.is_relative_to(ROOT/'outputs') and output != ROOT/'outputs', 'project output required')
    if args.worker:
        return worker(output)
    protocol = preflight(args)
    if args.preflight_only:
        print(json.dumps({k: protocol.get(k) for k in ('stage', 'particle_count', 'forecast_s', 'parent_time_s',
            'requested_end_time_s', 'stage_end_time_s', 'max_macrosteps', 'cost_estimate')}, indent=2)); return 0
    output.mkdir(parents=True, exist_ok=False); shared.write_new(output/'protocol.json', protocol)
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--stage', args.stage, '--output', rel(output), '--worker']
    env = dict(os.environ, **{k: '1' for k in shared.THREADS}); env['PYTHONPATH'] = ''
    began = time.perf_counter(); code = 124
    with (output/'run.log').open('x', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=600.).returncode
        except subprocess.TimeoutExpired:
            stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    error = None
    try:
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
    except Exception as problem:
        error = repr(problem); code = 2
    manifest = output/'manifest.json'
    shared.write_new(output/'supervision.json', dict(command=command, exit_code=code, wall_seconds=time.perf_counter()-began,
        source_inputs_still_match=error is None, source_error=error, manifest_sha256=core.sha(manifest) if manifest.exists() else None,
        log_sha256=core.sha(output/'run.log'), internal_limit_s=540., external_limit_s=600.,
        accepted=False, restart_allowed=False, accepted_time_advanced_s=0.))
    print(json.dumps(dict(exit_code=code, output=rel(output), status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    return code


if __name__ == '__main__':
    raise SystemExit(main())

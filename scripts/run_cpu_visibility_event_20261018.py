"""Explicit bounded activation-event integration of the original CPU physics.

No interpolated point becomes a state. Each event bracket endpoint is solved
again from exactly the same parent with one-sided common thermomechanical
support. All used support retains the original continuous-path certificate.
"""
from pathlib import Path
import ast
import inspect
import json
import shutil
import time
import numpy as np
import run_afsd_deposition_stage_remote_20261007 as original
import afsd_cpu_event_trial_20261015 as model
import visible_support_anchored_20261018 as path_proof
from afsd_event_locator_delay_20261010 import locate_activation as _legacy_locate_activation
from afsd_event_blocked_box_20261010 import audit_blocked_prefix
from visible_support import filter_visible_support

cpu = original
core, post, event, batch, old, shared = (original.core, original.post,
    original.event, original.batch, original.old, original.shared)
solver = original.solver
ROOT = original.ROOT
require, read, load, rel, verify = original.require, original.read, original.load, original.rel, original.verify
CONTINUE_PASS = 'PASS_CPU_SCREENED_EVENT_DIAGNOSTIC_CONTINUATION'
GeometryPathStop, RefinableCandidateError = original.GeometryPathStop, original.RefinableCandidateError
PROTOCOL = model.EVENT_PROTOCOL
LOCATOR_PROTOCOL = 'REAL_ENDPOINT_AND_PATH_DELAY_SECANT_20261010_V1'
BLOCKED_PREFIX_PROTOCOL = 'FIXED_WITNESS_BOX_CONTINUOUS_INCOMING_BLOCKED_PREFIX_20261010_V1'
VISIBLE_PATH_PROTOCOL = 'ENDPOINT_ANCHORED_SWEPT_SEGMENT_CERTIFICATE_20261018_V1'
INHERITED_CONTRACT = ROOT/'docs/GPU_VISIBILITY_EVENT_DELAY_20261010.md'
# Three proposals plus at most two independent ten-probe event brackets and
# their two post-event solves. The registered administrative deadline binds.
MAXIMUM_PHYSICAL_SOLVES = 25
MAXIMUM_EVENTS_PER_BRANCH = 2
EVENT_TIME_TOLERANCE_S = 1e-7
EVENT_RELATIVE_TIME_TOLERANCE = 1./1024.
EVENT_POSITION_TOLERANCE_M = 5e-8


def resolved_routing(route):
    if route.get('unresolved_samples'):
        return False
    changes = route.get('transition_brackets', [])
    if not changes:
        return route.get('status') == 'NO_OBSERVED_EVENT_REQUIRES_WHOLE_PATH_CERTIFICATE'
    return (len(changes) == 1 and route.get('status') == 'VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE'
        and changes[0].get('status') == 'BRACKET_LOCALIZED_NOT_UNIQUE_OR_PHYSICALLY_RESOLVED'
        and 0 <= changes[0]['bracket_width_fraction'] <= route['fraction_tolerance'])


def probe_prefix_is_ready(report):
    """Only a successful original proof can finish a probe; pure budget may retry."""
    if report['status'] == 'CERTIFIED_BLOCKED_PREFIX':
        require(report['failure_kind'] is None and not report['visible_samples']
            and not report['degenerate_edges'] and not report['numerical_domain_failures']
            and not report['unresolved_intervals']
            and not report.get('error') and report['maximum_queries'] == 256
            and 0 < report['analytic_batch_queries'] <= 256,
            'inconsistent successful event probe prefix proof')
        return True
    require(report['status'] == 'UNRESOLVED'
        and report['failure_kind'] == 'MAXIMUM_ANALYTIC_BATCH_QUERIES'
        and report['analytic_batch_queries'] == report['maximum_queries'] == 256
        and not report['visible_samples'] and not report['degenerate_edges']
        and report['unresolved_intervals']
        and all(row['reason'] == 'MAXIMUM_ANALYTIC_BATCH_QUERIES'
                for row in report['unresolved_intervals'])
        and not report['numerical_domain_failures'] and not report.get('error'),
        'event probe blocked-prefix proof has a non-budget failure')
    return False


def adaptation_audit():
    """Audit the inherited v5 geometry/event gates without importing GPU code."""
    inherited_path = ROOT/'scripts/run_gpu_visibility_event_delay_20261010.py'
    inherited = ast.parse(inherited_path.read_text(encoding='utf-8'))
    functions = {node.name: node for node in inherited.body if isinstance(node, ast.FunctionDef)}
    unchanged = ('resolved_routing', 'probe_prefix_is_ready', 'scalar_equivalence',
        'event_probe_guess', 'routing', 'qualify', 'leaf_row', 'copy_macro_endpoint',
        'classify_activation')
    for name in unchanged:
        actual = ast.parse(inspect.getsource(globals()[name])).body[0]
        if name == 'qualify':
            bindings = [node for node in actual.body if isinstance(node, ast.Assign)
                and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)
                and node.value.func.id == 'bind_parent_metadata']
            require(len(bindings) == 1, 'exactly one explicit path parent metadata binding required')
            actual.body.remove(bindings[0])
        require(ast.dump(actual) == ast.dump(functions[name]),
            'inherited event or physical qualification changed: '+name)
    require(solver is original.solver and core is original.core,
        'original CPU nonlinear/Newton providers required')
    return dict(status='PASS_CPU_EXPLICIT_COUPLED_EVENT_INTEGRATION', protocol=PROTOCOL,
        trial=model.operator_reuse_audit(), unchanged_event_functions=list(unchanged),
        inherited_event_driver=rel(inherited_path), inherited_event_driver_sha256=core.sha(inherited_path),
        cpu_nonlinear_provider=solver.__name__, cpu_newton_provider=core.__name__,
        event_locator_protocol=LOCATOR_PROTOCOL,
        probe_initial_guess='SAME_PARENT_MASKED_DISPLACEMENT_SCALED_WITH_PRESCRIBED_AXES_RESET',
        endpoint_policy_unchanged=True, activation_delay_explicit=True,
        event_time_tolerance_s=EVENT_TIME_TOLERANCE_S,
        event_relative_time_tolerance=EVENT_RELATIVE_TIME_TOLERANCE,
        event_position_tolerance_m=EVENT_POSITION_TOLERANCE_M,
        original_module_globals_rebound=False, gpu_provider_used=False)


def bounded_budget_limit(budget, key, default):
    """Administrative callers may tighten inherited budgets, never enlarge."""
    value = budget.get(key, default)
    require(type(value) is int and 1 <= value <= default, 'invalid bounded '+key)
    return value



def bind_parent_metadata(report, parent_checkpoint, parent_sha, parent_fingerprint):
    """Relabel stale fixture provenance; never alter a geometry result or gate."""
    require(core.sha(parent_checkpoint) == parent_sha, 'geometry report parent SHA changed')
    result = dict(report)
    legacy = result.pop('frozen_parent_sha256')
    result.update(parent_checkpoint=rel(parent_checkpoint), parent_sha256=parent_sha,
        actual_parent_fingerprint=parent_fingerprint,
        legacy_fixture_parent_sha256=legacy, legacy_fixture_parent_is_current=False,
        parent_binding='CALLER_VERIFIED_ACTUAL_CHECKPOINT_SHA256',
        inherited_geometry_algorithm_unchanged=True)
    return result


def locate_activation(solve_probe, full_dt, initial_fraction, *, parent_sha256, **kwargs):
    """Keep the legacy locator algorithm; explicitly bind its actual parent."""
    require(isinstance(parent_sha256, str) and len(parent_sha256) == 64
        and all(c in '0123456789abcdef' for c in parent_sha256),
        'actual parent SHA256 required for CPU event locator')
    location = _legacy_locate_activation(solve_probe, full_dt, initial_fraction, **kwargs)
    legacy = location.pop('frozen_parent_sha256')
    location.update(parent_sha256=parent_sha256,
        legacy_fixture_parent_sha256=legacy,
        legacy_fixture_parent_is_current=False,
        parent_binding='CALLER_VERIFIED_ACTUAL_CHECKPOINT_SHA256',
        inherited_locator_algorithm_unchanged=True)
    return location


def scalar_equivalence(trial, result, directory):
    scalar = model.ProcessTrial(trial.pa, trial.state, result['dt_s'], trial.incoming_pid_pairs)
    p = trial.p
    scalar.contact = old.ChainRuleDiscreteContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2, normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'], friction=p['friction'], base_z=-p['substrate_depth'])
    sr = scalar.evaluate(result['displacement'], result['dt_s'])
    differences = batch.candidate_equivalence(dict(sr, energy_comparison=core.energy_diagnostics(scalar, sr)),
        dict(result, energy_comparison=core.energy_diagnostics(trial, result)))
    require(sr['graph_fingerprint'] == result['graph_fingerprint'] and np.array_equal(sr['g'], result['g'])
        and np.array_equal(sr['thermal']['pair_conductance_W_K'], result['thermal']['pair_conductance_W_K']),
        'scalar/batch exact common event graph and heat differ')
    core.write(directory/'scalar_batch_equivalence.json', dict(status='PASS_COMPLETE_MAP_EQUIVALENCE', differences=differences))
    return differences



def event_probe_guess(trial, dt, prior):
    """A numerical initial guess only; every later state is independently solved."""
    previous = prior['trial']
    require(cpu.parent_fingerprint(previous.pa, previous.state) == cpu.parent_fingerprint(trial.pa, trial.state),
        'warm-start probe must have the same immutable parent')
    require(tuple(previous.incoming_pid_pairs) == tuple(trial.incoming_pid_pairs)
        and bool(trial.incoming_pid_pairs), 'warm-start event identity differs')
    previous_dt = previous.motion.dt_s
    require(np.isfinite(previous_dt) and previous_dt > 0. and np.isfinite(dt) and dt > 0.,
        'positive finite warm-start interval required')
    guess = np.array(prior['fresh']['displacement'], dtype=float, copy=True)*(dt/previous_dt)
    require(guess.shape == trial.x0.shape and np.isfinite(guess).all(), 'invalid prior probe displacement')
    fixed = ~trial.state['free']
    guess[fixed] = trial.prescribed(dt)[fixed]
    require(np.max(np.linalg.norm(guess,axis=1)) <= .1*trial.p['dx'],
        'scaled probe guess exceeds original displacement domain')
    return guess


def solve_physical(pa, state, dt, directory, deadline, per_candidate, frozen, budget, incoming=(), initial_guess=None):
    require(time.perf_counter()+5. < deadline, 'no time for another physical solve')
    maximum_solves = bounded_budget_limit(budget, 'maximum_solves', MAXIMUM_PHYSICAL_SOLVES)
    require(type(budget.get('solves')) is int and 0 <= budget['solves'] < maximum_solves,
        'event physical solve budget exhausted')
    budget['solves'] += 1
    directory.mkdir(parents=True)
    began = time.perf_counter()
    fingerprint = cpu.parent_fingerprint(pa, state)
    trial = model.ProcessTrial(pa, state, dt, incoming)
    attempt = directory/'broyden'; attempt.mkdir()
    guess = batch.guess_for(trial, dt)
    if initial_guess is not None:
        require(bool(incoming), 'warm start is restricted to explicitly masked event probes')
        guess = event_probe_guess(trial, dt, initial_guess)
        core.write(directory/'probe_initial_guess.json', dict(
            status='UNACCEPTED_NONLINEAR_INITIAL_GUESS_ONLY',
            source_checkpoint=rel(initial_guess['checkpoint']),
            source_sha256=core.sha(initial_guess['checkpoint']),
            source_dt_s=initial_guess['trial'].motion.dt_s, target_dt_s=dt,
            scale=dt/initial_guess['trial'].motion.dt_s,
            parent_fingerprint=fingerprint, prescribed_axes_reset=True,
            state_interpolated_or_accepted=False))
    try:
        try:
            result, metrics = solver.solve(trial, dt, guess, attempt, deadline, 24)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            require(time.perf_counter()+1.5*per_candidate+20 < deadline, 'no budget for unchanged Newton retry')
            attempt = directory/'newton_retry'; attempt.mkdir()
            result, metrics = core.solve(trial, dt, guess, attempt, deadline, 24)
            metrics['fallback_reason'] = repr(error)
        core.journal(directory, trial, result, 'SOLVED_UNCOMMITTED_candidate')
        old.save_graph(directory, trial, result)
        nextpa, nextstate, fresh = model.advance_candidate(trial, result)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        raise RefinableCandidateError(str(error)) from error
    require(cpu.parent_fingerprint(pa, state) == fingerprint, 'event probe mutated parent')
    require(np.array_equal(pa.particle_id, nextpa.particle_id) and np.array_equal(pa.m, nextpa.m)
        and model.identity_digest(pa) == model.identity_digest(nextpa), 'immutable material inventory changed')
    target, _ = event.save_candidate(directory, trial, nextpa, nextstate, fresh, frozen)
    # Continue from the exact saved state, including its registered source map
    # and diagnostic flags, rather than the pre-save in-memory metadata.
    nextpa, nextstate = load(target)
    require(nextstate['implementation_source_sha256'] == frozen,
        'saved physical substep implementation identity differs')
    old.save_graph(directory, trial, fresh)
    observer = shared.observe_process(pa, nextstate, fresh, trial.motion, trial.p)
    core.write(directory/'observer.json', observer); core.write(directory/'motion.json', trial.motion.to_dict())
    core.write(directory/'transport.json', nextstate['bore_observation'])
    core.write(directory/'physical_solve_cost.json', dict(
        status='MEASURED_COMPLETED_PHYSICAL_SOLVE', wall_seconds=time.perf_counter()-began,
        checkpoint=rel(target), checkpoint_sha256=core.sha(target), parent_fingerprint=fingerprint,
        dt_s=dt, incoming_pid_pairs=[list(pair) for pair in incoming],
        scope='NONLINEAR_SOLVE_FRESH_ADVANCE_AND_SAVED_STATE_EXCLUDES_LATER_GEOMETRY_QUALIFICATION',
        accepted=False))
    return dict(pa=nextpa, state=nextstate, checkpoint=target, directory=directory, trial=trial, fresh=fresh,
        observer_path=directory/'observer.json', observer=observer, metrics=metrics,
        physical_solve_wall_seconds=time.perf_counter()-began, began=began)


def routing(item, union=False):
    trial, fresh = item['trial'], item['fresh']
    if union:
        si, sj, _, _, _ = model.operators(trial.pa)
        pairs = np.unique(np.concatenate((np.column_stack((si, sj)),
            np.column_stack((fresh['raw_i'], fresh['raw_j'])))), axis=0)
        i, j = pairs[:, 0], pairs[:, 1]
    else:
        i, j = fresh['i'], fresh['j']
    result = original.remote_audit_visibility_events(trial.x0, fresh['position'], i, j,
        trial.motion.pose(trial.motion.t0), trial.motion.pose(trial.motion.t1), trial.p, trial.pa.particle_id)
    name = 'union_visibility_events.json' if union else 'visibility_event_routing.json'
    core.write(item['directory']/name, result)
    return result


def qualify(item):
    trial, fresh, directory = item['trial'], item['fresh'], item['directory']
    union = routing(item, union=True)
    incoming = set(trial.incoming_pid_pairs)
    union_valid = (resolved_routing(union) and all(
        r['direction'] == 'BLOCKED_TO_VISIBLE'
        and tuple(sorted((r['particle_id_i'],r['particle_id_j']))) in incoming
        for r in union['transition_brackets']))
    if not incoming:
        union_valid = union_valid and union['status'] == 'NO_OBSERVED_EVENT_REQUIRES_WHOLE_PATH_CERTIFICATE'
    item['union_routing'] = union
    scalar_equivalence(trial, fresh, directory)
    geometry, events = shared.inherited.geometry_audit(trial, fresh)
    core.write(directory/'geometry_audit.json', geometry); core.write(directory/'sampled_events.json', events)
    pose0, pose1 = trial.motion.pose(trial.motion.t0), trial.motion.pose(trial.motion.t1)
    shell = old.audit_contact_shell_finite_enclosure(trial.contact, trial.x0, fresh['position'], pose0, pose1)
    core.write(directory/'continuous_shell_audit.json', shell)
    path = path_proof.audit_visible_support_path(trial.x0, fresh['position'], fresh['i'], fresh['j'],
        pose0, pose1, trial.p, trial.pa.particle_id)
    path = bind_parent_metadata(path, item['parent_checkpoint'], item['parent_sha256'],
        cpu.parent_fingerprint(trial.pa, trial.state))
    core.write(directory/'visible_path_audit.json', path)
    item.update(geometry=geometry, shell=shell, path=path)
    valid = (union_valid and not geometry['transient_contacts_missed_by_endpoints']
        and not set(geometry['new_non_bore_features'])-{'SHOULDER', 'TOP'}
        and shell['status'] == event.CERTIFIED and path['status'] == old.CERTIFIED)
    if not valid:
        item['routing'] = routing(item)
        return False
    qualification = dict(status='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN', accepted=False,
        restart_allowed=False, metrics=item['metrics'], event_protocol=PROTOCOL,
        delayed_activation_pairs=[list(p) for p in trial.incoming_pid_pairs])
    core.write(directory/'candidate_qualification.json', qualification)
    core.vtk(directory/'fields.vtk', item['pa'], trial.p, trial.bulk, item['state']['time_s'], model.CLASSIFICATION)
    item['wall_seconds'] = time.perf_counter()-item['began']
    return True


def leaf_row(item, parent_checkpoint, parent_sha, parent_state):
    return dict(checkpoint=rel(item['checkpoint']), sha256=core.sha(item['checkpoint']),
        parent_checkpoint=rel(parent_checkpoint), parent_sha256=parent_sha,
        parent_time_s=parent_state['time_s'], parent_step=parent_state['step'],
        time_s=item['state']['time_s'], step=item['state']['step'], dt_s=item['trial'].motion.dt_s,
        observer_path=rel(item['observer_path']), qualification_path=rel(item['directory']/'candidate_qualification.json'))


def copy_macro_endpoint(item, directory):
    checkpoint = directory/'NOT_RESTARTABLE_candidate.npz'
    shutil.copyfile(item['checkpoint'], checkpoint)
    shutil.copyfile(item['observer_path'], directory/'observer.json')
    core.write(directory/'candidate_qualification.json', dict(
        status='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN', accepted=False, restart_allowed=False,
        event_protocol=PROTOCOL, physical_endpoint_sha256=core.sha(checkpoint),
        original_candidate_gates_passed=True))
    return dict(item, checkpoint=checkpoint, observer_path=directory/'observer.json')


def classify_activation(item):
    route = item.get('routing') or routing(item)
    require(item['shell']['status'] == event.CERTIFIED
        and not item['geometry']['transient_contacts_missed_by_endpoints']
        and not set(item['geometry']['new_non_bore_features'])-{'SHOULDER', 'TOP'},
        'event resolver cannot repair shell/contact failure')
    require(item['path']['failure_kind'] == 'STRICT_SOLID_INTERSECTION_ON_FIXED_RETAINED_PATH'
        and not item['path']['unresolved_intervals'] and not item['path']['degenerate_edges']
        and not item['path']['numerical_domain_failures'], 'not a resolved activation candidate')
    require(resolved_routing(route) and route['status'] == 'VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE'
        and len(route['transition_brackets']) == 1 and len(route['changed_edges']) == 1
        and not route['unresolved_samples'], 'only one observed activation per event solve')
    bracket = route['transition_brackets'][0]
    require(bracket['direction'] == 'BLOCKED_TO_VISIBLE', 'deactivation requires separate integration')
    pair = tuple(sorted((bracket['particle_id_i'], bracket['particle_id_j'])))
    require(all(tuple(sorted((r['particle_id_i'], r['particle_id_j']))) == pair and r['fraction'] == 0.
        for r in item['path']['found_edges']), 'additional strict intersection cannot be repaired')
    return (pair,), .5*(bracket['start_fraction']+bracket['end_fraction'])


def complete_branch(pa, state, dt, directory, deadline, per_candidate, frozen, budget,
                    parent_checkpoint, parent_sha, tolerance_scale=1., depth=0):
    maximum_events = bounded_budget_limit(budget, 'maximum_events_per_branch', MAXIMUM_EVENTS_PER_BRANCH)
    require(core.sha(parent_checkpoint) == parent_sha, 'actual branch parent checkpoint SHA changed')
    began = time.perf_counter(); directory.mkdir()
    candidate = solve_physical(pa, state, dt, directory/'proposal', deadline, per_candidate, frozen, budget)
    candidate.update(parent_checkpoint=parent_checkpoint, parent_sha256=parent_sha)
    if qualify(candidate):
        row = leaf_row(candidate, parent_checkpoint, parent_sha, state)
        result = copy_macro_endpoint(candidate, directory)
        result.update(leaves=[row], event_records=[], wall_seconds=time.perf_counter()-began)
        return result
    try:
        incoming, initial = classify_activation(candidate)
    except Exception:
        # Preserve ordinary root-level evidence for the original pure-budget
        # recovery recognizer; no rejected state is used as a parent.
        for name in ('geometry_audit.json', 'continuous_shell_audit.json', 'visible_path_audit.json',
                     'visibility_event_routing.json'):
            shutil.copyfile(candidate['directory']/name, directory/name)
        raise GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')
    require(depth < maximum_events, 'event count limit reached')
    union = routing(candidate, union=True)
    target = set(incoming)
    require(resolved_routing(union) and union['status'] == 'VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE'
        and not union['unresolved_samples']
        and all(r['direction'] == 'BLOCKED_TO_VISIBLE' and tuple(sorted((r['particle_id_i'],r['particle_id_j']))) in target
                for r in union['transition_brackets'])
        and len(union['transition_brackets']) == 1, 'unknown union event or deactivation')
    probes = []
    def solve_probe(probe_dt):
        initial_guess = min(probes, key=lambda q:abs(q['trial'].motion.dt_s-probe_dt)) if probes else None
        probe = solve_physical(pa, state, probe_dt, directory/f'event_probe_{len(probes)+1:02d}',
            deadline, per_candidate, frozen, budget, incoming, initial_guess)
        probe.update(parent_checkpoint=parent_checkpoint, parent_sha256=parent_sha)
        probes.append(probe)
        diagnostic = routing(probe, union=True)
        require(resolved_routing(diagnostic), 'unresolved event probe geometry')
        changes = diagnostic['transition_brackets']
        require(len(changes) <= 1 and all(r['direction'] == 'BLOCKED_TO_VISIBLE'
            and tuple(sorted((r['particle_id_i'],r['particle_id_j']))) in target for r in changes),
            'event probe has unknown event or reversal')
        trial, fresh = probe['trial'], probe['fresh']
        ids = {int(pid):n for n,pid in enumerate(pa.particle_id)}
        i = np.array([min(ids[p[0]],ids[p[1]]) for p in incoming])
        j = np.array([max(ids[p[0]],ids[p[1]]) for p in incoming])
        endpoint = filter_visible_support(fresh['position'], i, j, np.zeros((len(i),3)),
            trial.motion.pose(trial.motion.t1), trial.p, pa.particle_id)
        blocked = bool(len(endpoint['blocked_raw_edge_indices']))
        require(blocked or len(changes) == 1, 'visible endpoint lacks localized activation')
        fraction = None if blocked else .5*(changes[0]['start_fraction']+changes[0]['end_fraction'])
        measure = max((float(row['strict_maximum_contiguous_length_m'])
            for row in endpoint['diagnostics']['blocked_edges']), default=0.)
        activation_delay = None if blocked else probe_dt*(1.-changes[0]['start_fraction'])
        prefix_ready = False
        if activation_delay is not None and activation_delay <= tolerance:
            prefix_end = 1.-tolerance/probe_dt
            require(0. < prefix_end < 1., 'bounded activation probe requires a positive blocked prefix')
            prefix = audit_blocked_prefix(trial.x0, fresh['position'], i, j,
                trial.motion.pose(trial.motion.t0), trial.motion.pose(trial.motion.t1),
                trial.p, pa.particle_id, prefix_end)
            prefix = bind_parent_metadata(prefix, parent_checkpoint, parent_sha,
                cpu.parent_fingerprint(pa, state))
            core.write(probe['directory']/'probe_incoming_blocked_prefix.json', prefix)
            prefix_ready = probe_prefix_is_ready(prefix)
        return dict(blocked=blocked, blocked_measure=measure, event_fraction=fraction, payload=probe,
            activation_delay_s=activation_delay, activation_prefix_certified=prefix_ready,
            positions=fresh['position'], event_count=len(changes), event_kind='BLOCKED_TO_VISIBLE', outgoing=False)
    tolerance = min(EVENT_TIME_TOLERANCE_S, dt*EVENT_RELATIVE_TIME_TOLERANCE)*tolerance_scale
    location = locate_activation(solve_probe, dt, initial, parent_sha256=parent_sha, time_tolerance=tolerance,
        maximum_solves=10, minimum_dt=1e-10, maximum_position_difference=EVENT_POSITION_TOLERANCE_M)
    audit = {k:v for k,v in location.items() if k not in ('lo','hi','attempts')}
    audit['attempts'] = [dict(**{k:v for k,v in q.items() if k not in ('payload','positions')},
        checkpoint=rel(q['payload']['checkpoint']), sha256=core.sha(q['payload']['checkpoint']))
        for q in location['attempts'] if q.get('payload') is not None]
    for side in ('lo','hi'):
        if location.get(side):
            q=location[side];audit[side]=dict(dt_s=q['dt_s'], blocked=q['blocked'],
                checkpoint=rel(q['payload']['checkpoint']), sha256=core.sha(q['payload']['checkpoint']))
    core.write(directory/'event_location.json', audit)
    require(location['status'] == 'BRACKETED_ACTIVATION', 'real physical event bracket unresolved')
    low, high = location['lo']['payload'], location['hi']['payload']
    high_route = routing(high, union=True)
    require(resolved_routing(high_route) and len(high_route['transition_brackets']) == 1,
        'selected high endpoint lacks resolved activation on its own path')
    high_delay = high['trial'].motion.dt_s*(1.-high_route['transition_brackets'][0]['start_fraction'])
    require(high_delay <= tolerance, 'selected high path activation delay exceeds registered time bound')
    prefix_end = 1.-tolerance/high['trial'].motion.dt_s
    require(0. < prefix_end < 1., 'bounded activation requires a positive blocked prefix')
    index = {int(pid):n for n,pid in enumerate(pa.particle_id)}
    prefix_i = np.array([min(index[a],index[b]) for a,b in incoming])
    prefix_j = np.array([max(index[a],index[b]) for a,b in incoming])
    prefix = audit_blocked_prefix(high['trial'].x0, high['fresh']['position'], prefix_i, prefix_j,
        high['trial'].motion.pose(high['trial'].motion.t0), high['trial'].motion.pose(high['trial'].motion.t1),
        high['trial'].p, pa.particle_id, prefix_end)
    prefix = bind_parent_metadata(prefix, parent_checkpoint, parent_sha,
        cpu.parent_fingerprint(pa, state))
    prefix_path = high['directory']/'incoming_blocked_prefix.json'
    core.write(prefix_path, prefix)
    require(prefix['status'] == 'CERTIFIED_BLOCKED_PREFIX',
        'incoming edge not continuously certified blocked before bounded activation window')
    require(qualify(low) and qualify(high), 'one-sided event substep failed original geometry/physical gates')
    oldrow = leaf_row(high, parent_checkpoint, parent_sha, state)
    remaining = state['time_s']+dt-high['state']['time_s']
    require(remaining > 1e-10, 'event leaves no positive post-event interval')
    after = complete_branch(high['pa'], high['state'], remaining, directory/'after_event', deadline,
        per_candidate, frozen, budget, high['checkpoint'], core.sha(high['checkpoint']), tolerance_scale, depth+1)
    first_graph = np.load(ROOT/after['leaves'][0]['checkpoint'].replace('NOT_RESTARTABLE_candidate.npz','visible_graph.npz'))
    try:
        ids = first_graph['particle_id']; restored = {tuple(sorted((int(ids[i]),int(ids[j]))))
            for i,j in zip(first_graph['i'],first_graph['j'])}
        require(target <= restored, 'post-event physical graph did not restore incoming support')
    finally:
        first_graph.close()
    resolution = dict(status='PASS_BOUNDED_PHYSICALLY_RESOLVED_ACTIVATION', protocol=PROTOCOL,
        incoming_pid_pairs=[list(v) for v in incoming], location=rel(directory/'event_location.json'),
        location_sha256=core.sha(directory/'event_location.json'), tolerance_s=tolerance,
        low_checkpoint=rel(low['checkpoint']), high_checkpoint=rel(high['checkpoint']),
        bracket_width_s=location['bracket_width_s'],
        maximum_position_difference_m=location['maximum_position_difference_m'],
        first_post_event_checkpoint=after['leaves'][0]['checkpoint'], incoming_restored=True,
        original_active_path_certified=True, accepted=False, restart_allowed=False,
        interpolation_accepted=False, real_solution_bracket_width_s=location['bracket_width_s'],
        observed_high_path_activation_delay_s=high_delay, activation_delay_bound_s=tolerance,
        blocked_prefix_path=rel(prefix_path), blocked_prefix_sha256=core.sha(prefix_path),
        blocked_prefix_end_fraction=prefix_end)
    core.write(high['directory']/'event_resolution.json', resolution)
    qualification = read(high['directory']/'candidate_qualification.json')
    qualification.update(event_resolution_path=rel(high['directory']/'event_resolution.json'),
        event_resolution_sha256=core.sha(high['directory']/'event_resolution.json'))
    core.write(high['directory']/'candidate_qualification.json', qualification)
    result = copy_macro_endpoint(after, directory)
    result.update(leaves=[oldrow]+after['leaves'], event_records=[resolution]+after['event_records'],
        wall_seconds=time.perf_counter()-began)
    return result


def run_pair(output, protocol, report, deadline):
    require(core.sha(ROOT/protocol['parent_checkpoint']) == protocol['parent_sha256'],
        'actual pair parent checkpoint SHA changed')
    pa, state = load(ROOT/protocol['parent_checkpoint']); dt = protocol['dt_s']
    fingerprint = cpu.parent_fingerprint(pa, state)
    require(protocol.get('event_tolerance_scale',1.) in (1.,.5,.25), 'unregistered event tolerance scale')
    program = model.validate_state(state); phase = program.phase(state['time_s'])
    require(state['time_s']+dt <= phase['end']+1e-15, 'CPU interval crosses machine phase')
    require(phase['target'][2] <= 0 and phase['initial'][2] <= 0, 'stationary fill only')
    items=[]; budget=dict(solves=0)
    for role in ('coarse','fine_1','fine_2'):
        parent = items[1] if role == 'fine_2' else None
        pp, ss = (parent['pa'],parent['state']) if parent else (pa,state)
        checkpoint = parent['checkpoint'] if parent else ROOT/protocol['parent_checkpoint']
        sha = core.sha(checkpoint) if parent else protocol['parent_sha256']
        item=complete_branch(pp,ss,dt if role=='coarse' else dt/2,output/role,deadline,
            protocol['candidate_forecast_s'],protocol['source_sha256'],budget,checkpoint,sha,
            protocol.get('event_tolerance_scale',1.))
        items.append(item)
        require(cpu.parent_fingerprint(pa,state)==fingerprint,'original parent mutated')
        print(json.dumps(dict(role=role,wall_s=item['wall_seconds'],time_s=item['state']['time_s'],
            physical_substeps=len(item['leaves']),activation_events=len(item['event_records']))),flush=True)
    coarse,first,last=items
    global_screen=post.compare(last['pa'],coarse['pa'],state['parameters'])
    local_screen=event.local_compare(last['pa'],coarse['pa'],sorted(int(v) for v in pa.particle_id))
    matching=last['state']['bore_observation']['crossing_ids']==coarse['state']['bore_observation']['crossing_ids']
    screen=dict(global_screen=global_screen,local_screen=local_screen,crossing_ids_match=matching,
        coarse_checkpoint_sha256=core.sha(coarse['checkpoint']),
        fine_checkpoint_sha256=[core.sha(first['checkpoint']),core.sha(last['checkpoint'])],
        passed=bool(global_screen['pass'] and local_screen['passed'] and matching),accepted=False)
    core.write(output/'pair_screen.json',screen)
    require(screen['passed'],'original coarse/fine screen failed; no CPU event continuation')
    prefix_n=len(state['machine_history'])
    require(last['state']['machine_history'][:prefix_n]==state['machine_history']
        and last['state']['interval_energy']['records'][:prefix_n]==state['interval_energy']['records'],
        'history prefix changed')
    report.update(status=CONTINUE_PASS,pair_screen=screen,screened_continuation_checkpoint=rel(last['checkpoint']),
        screened_continuation_sha256=core.sha(last['checkpoint']),screened_end_time_s=last['state']['time_s'],
        screened_end_step=last['state']['step'],parent_time_s=state['time_s'],parent_fingerprint=fingerprint,
        screened_continuation_authorized=True,parent_prefix_unchanged=True,
        maximum_candidate_wall_s=max(v['wall_seconds'] for v in items),candidate_wall_s=[v['wall_seconds'] for v in items],
        screened_history=[dict(checkpoint=rel(v['checkpoint']),sha256=core.sha(v['checkpoint']),
            time_s=v['state']['time_s'],observer_path=rel(v['observer_path'])) for v in items[1:]],
        event_substeps=dict(schema='EXPLICIT_COUPLED_EVENT_SUBSTEPS_V1',
            **{role:item['leaves'] for role,item in zip(('coarse','fine_1','fine_2'),items)}),
        event_records={role:item['event_records'] for role,item in zip(('coarse','fine_1','fine_2'),items)},
        event_protocol=PROTOCOL,event_locator_protocol=LOCATOR_PROTOCOL,physical_solve_count=budget['solves'],
        observer=last['observer'],full_program_completed=False)

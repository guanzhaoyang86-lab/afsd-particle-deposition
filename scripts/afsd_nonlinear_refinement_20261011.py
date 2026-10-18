"""Administrative half-step retry of a narrowly evidenced nonlinear failure.

Only the first coarse/proposal Broyden-then-Newton layout is supported. The
caller must bind these six original JSON reports and associated NPZ files to
the sealed manifest, exact sources, complete parent, pending dt and parent free
mask. It must verify finite NPZ arrays, original force threshold/residual and
unchanged inventory, and disallow later, committed or qualified candidates and
pair screens. This module performs no physical solve or file I/O. A rejected
candidate remains rejected; convergence and all original gates are still
required after a genuine smaller solve from the same complete parent.
"""
import copy
import math

import afsd_adaptive_dt_20261010 as adaptive


PROTOCOL = 'UNCHANGED_PARENT_BROYDEN_NEWTON_NONCONVERGENCE_REFINEMENT_20261011_V1'
NONLINEAR_ERROR = "RefinableCandidateError('complete residual line search/evaluation budget failed')"
SECANT_ERROR = "ValueError('complete residual Armijo search failed; refresh budget exhausted')"
CONTROLLER = 'BOUNDED_BROYDEN_DISCARD_ON_VISIBLE_GRAPH_CHANGE_V1'
REPORT_NAMES = ('secant_failure', 'broyden_iterations', 'broyden_candidate',
    'broyden_graph', 'newton_iterations', 'newton_candidate')
MINIMUM_DT_S = 1e-7
COOLDOWN_PASSES = 8
METRICS = ('force_residual_N', 'force_limit_N', 'scaled_residual', 'density_residual',
    'thermal_residual_K', 'virtual_work_residual_J', 'material_endpoint_work_J',
    'density_iterations', 'thermal_iterations', 'wall_seconds')


def _require(condition, reason):
    if not condition:
        raise ValueError(reason)


def _integer(report, key, minimum=0):
    value = report.get(key)
    _require(type(value) is int and value >= minimum, 'INVALID_INTEGER_'+key)
    return value


def _finite_tree(value):
    """Reject nonfinite saved JSON numerics without trusting string conversion."""
    if type(value) in (int, float):
        _require(math.isfinite(value), 'NONFINITE_DIAGNOSTIC_VALUE')
    elif type(value) is dict:
        for child in value.values():
            _finite_tree(child)
    elif type(value) is list:
        for child in value:
            _finite_tree(child)
    else:
        _require(value is None or type(value) in (str, bool), 'NON_JSON_DIAGNOSTIC_VALUE')


def _metrics(row):
    values = {key: adaptive.finite(row[key], key) for key in METRICS}
    force, limit = values['force_residual_N'], values['force_limit_N']
    _require(limit >= 1.01e-5 and force > limit
        and values['scaled_residual'] == force/limit, 'NOT_ORIGINAL_UNCONVERGED_FORCE_RESIDUAL')
    _require(0. <= values['density_residual'] <= 1e-10
        and 0. <= values['thermal_residual_K'] <= 1e-8
        and abs(values['virtual_work_residual_J']) <= 1e-12+1e-10*abs(values['material_endpoint_work_J'])
        and 1 <= _integer(row, 'density_iterations') <= 40
        and 1 <= _integer(row, 'thermal_iterations') <= 30
        and values['wall_seconds'] >= 0., 'ORIGINAL_INNER_PHYSICAL_DOMAIN_FAILED')


def _iterations(rows, broyden):
    _require(type(rows) is list and 2 <= len(rows) <= 24, 'INCOMPLETE_NONLINEAR_ITERATIONS')
    _finite_tree(rows)
    evaluations = 1
    builds = []
    for index, row in enumerate(rows):
        _require(type(row) is dict and _integer(row, 'iteration') == index
            and _integer(row, 'evaluations', 1) == evaluations, 'INCONSISTENT_NONLINEAR_ITERATION_ACCOUNTING')
        _metrics(row)
        line = row['line_search']
        _require(type(line) is list and 1 <= len(line) <= 8, 'INCOMPLETE_OR_EXTENDED_LINE_SEARCH')
        for backtrack, trial in enumerate(line):
            _require(type(trial) is dict and 'rejected' not in trial and 'error' not in trial
                and adaptive.finite(trial['alpha'], 'backtrack alpha') == .5**backtrack
                and adaptive.finite(trial['merit_N'], 'line search merit') >= 0.,
                'PHYSICAL_DOMAIN_REJECTION_OR_INVALID_BACKTRACK')
            if broyden:
                adaptive.digest(trial['graph_fingerprint'])
                _require(type(trial.get('graph_changed')) is bool
                    and trial['graph_changed'] == (trial['graph_fingerprint'] != row['graph_fingerprint'])
                    and 1 <= _integer(trial, 'jacobian_build_number') <= 3,
                    'INVALID_BROYDEN_SEARCH_GRAPH_OR_BUILD')
        evaluations += len(line)
        _require(evaluations < 100, 'EVALUATION_OR_TIME_BUDGET_NOT_IN_SCOPE')
        linear = row['linear']
        if broyden:
            _require(row.get('controller') == CONTROLLER
                and adaptive.finite(row['actual_free_force_residual_N'], 'actual force') == row['force_residual_N']
                and type(linear) is list and len(linear) == 1, 'INVALID_BROYDEN_IDENTITY_OR_LINEAR_ROUTE')
            adaptive.digest(row['graph_fingerprint'])
            for build in row.get('jacobian_builds', []):
                _require(build.get('reason') in ('INITIAL_FULL_TANGENT', 'ACTUAL_SUPPORT_GRAPH_CHANGED',
                    'POOR_PROGRESS', 'HEAVY_DAMPING', 'ARMIJO_FAILURE')
                    and 0. < adaptive.finite(build['condition_infinity'], 'condition') <= 1e12,
                    'UNSAFE_OR_UNKNOWN_TANGENT_REFRESH')
                builds.append(_integer(build, 'build_number', 1))
            linear = linear[0]
            _require(linear.get('fixed_axis_increment_zero') is True
                and linear.get('condition_policy') == 'INITIAL_INVERSE_NORM_THEN_SHERMAN_MORRISON_GUARDED_ESTIMATE'
                and linear.get('jacobian_build_number') == line[0]['jacobian_build_number'],
                'UNREGISTERED_BROYDEN_LINEAR_POLICY')
        else:
            _require(type(linear) is dict
                and linear.get('coercivity') == 'NOT_CERTIFIED_BY_LU; INITIAL_ELASTIC_RANK_AUDIT_REQUIRED'
                and linear.get('inertia') is False and 'terminal_error' not in row,
                'UNREGISTERED_NEWTON_LINEAR_POLICY')
        _require(0. <= adaptive.finite(linear['linear_backward_error_N'], 'linear backward error')
            <= 1e-9*max(1., row['force_residual_N']) and _integer(linear, 'nnz', 1) > 0,
            'ORIGINAL_LINEAR_ACCURACY_FAILED')
        if index < len(rows)-1:
            _require('terminal_error' not in row and 'evaluations_at_failure' not in row,
                'EARLIER_TERMINAL_FAILURE')
            if broyden:
                _require(_integer(row, 'evaluations_after_search') == evaluations
                    and adaptive.finite(row['accepted_alpha'], 'accepted alpha') == line[-1]['alpha']
                    and 0. <= adaptive.finite(row['merit_ratio'], 'merit ratio') < 1.-1e-4*line[-1]['alpha'],
                    'INCONSISTENT_BROYDEN_ACCEPTED_ITERATE')
                secant = row['secant']
                _require(secant.get('status') in ('GOOD_BROYDEN_RANK_ONE', 'DISCARDED_ACROSS_GRAPH_CHANGE'),
                    'UNSAFE_SECANT_UPDATE')
                if secant['status'] == 'DISCARDED_ACROSS_GRAPH_CHANGE':
                    _require(adaptive.digest(secant['previous_graph']) == row['graph_fingerprint']
                        and adaptive.digest(secant['current_graph']) == line[-1]['graph_fingerprint']
                        and secant['previous_graph'] != secant['current_graph'], 'INCONSISTENT_GRAPH_DISCARD')
        else:
            _require(len(line) == 8 and 'accepted_alpha' not in row
                and 'evaluations_after_search' not in row and 'secant' not in row,
                'LAST_LINE_SEARCH_NOT_EXHAUSTED')
            # The preceding accepted line-search endpoint is this iterate, so
            # its recorded L2 merit certifies that all eight final trials fail
            # the unchanged Armijo inequality; no physical state interpolation.
            merit = adaptive.finite(rows[-2]['line_search'][-1]['merit_N'], 'last iterate merit')
            _require(merit > 0. and all(trial['merit_N'] >= (1.-1e-4*trial['alpha'])*merit
                for trial in line), 'LAST_TRIAL_ACTUALLY_SATISFIES_ARMIJO')
            if broyden:
                _require(row.get('terminal_error') == SECANT_ERROR
                    and _integer(row, 'evaluations_at_failure') == evaluations,
                    'NOT_EXACT_BROYDEN_TERMINAL_FAILURE')
    if broyden:
        _require(builds == [1, 2, 3], 'ORIGINAL_TWO_REFRESH_BUDGET_NOT_EXHAUSTED')
    return evaluations


def _candidate(candidate, last, dt):
    _require(type(candidate) is dict and candidate.get('accepted') is False
        and candidate.get('restart_allowed') is False, 'CANDIDATE_NOT_UNCOMMITTED')
    _finite_tree(candidate)
    metrics = candidate['metrics']
    _require(type(metrics) is dict and set(metrics) == set(METRICS)
        and all(type(metrics[key]) is type(last[key]) and metrics[key] == last[key] for key in METRICS),
        'MIXED_LAST_ITERATE_AND_CANDIDATE')
    contact = candidate['contact']
    _require(contact.get('classification') == 'HOT_START_QUASISTATIC_CHAIN_RULE_DISCRETE_NORMAL_BE_FRICTION_CONTACT_APPROXIMATE_COMPONENT'
        and adaptive.finite(contact['dt_s'], 'contact dt') == dt
        and contact['sweep'].get('certified') is True
        and contact['sweep'].get('interpolation') == 'LINEAR_PARTICLE_AND_ORIGIN_SHORTEST_SO3_TOOL_ROTATION',
        'UNBOUND_DT_OR_UNCERTIFIED_CENTRE_SWEEP')
    _require(type(candidate['energy']) is dict and type(candidate['thermal']) is dict,
        'INCOMPLETE_PHYSICAL_DIAGNOSTICS')
    energy = candidate['energy']
    scale = adaptive.finite(energy['absolute_work_scale_J'], 'energy scale')
    limit = adaptive.finite(energy['registered_energy_limit_J'], 'energy limit')
    _require(scale >= 0. and limit == 1e-7+.05*scale
        and abs(adaptive.finite(energy['global_raw_residual_J'], 'energy residual')) <= limit,
        'ORIGINAL_RAW_ENERGY_DIAGNOSTIC_FAILED')


def recognize(manifest, supervision, reports):
    """Fail closed; eligible means retry only, never accepted or restartable."""
    diagnosis = dict(protocol=PROTOCOL, eligible=False, reason='UNRECOGNIZED_NONLINEAR_FAILURE',
        diagnostic_only=True, accepted=False, restart_allowed=False, accepted_time_advanced_s=0.)
    try:
        _require(type(manifest) is dict and type(supervision) is dict and type(reports) is dict
            and set(reports) == set(REPORT_NAMES), 'INCOMPLETE_OR_UNEXPECTED_REPORT_SET')
        _require(manifest.get('status') == 'FAIL_GPU_STAGE_NOT_ACCEPTED'
            and manifest.get('error') == NONLINEAR_ERROR and manifest.get('accepted') is False
            and manifest.get('restart_allowed') is False and manifest.get('diagnostic_only') is True
            and manifest.get('screened_continuation_authorized') is False
            and adaptive.finite(manifest['accepted_time_advanced_s'], 'advanced time') == 0.,
            'NOT_EXACT_REJECTED_NONLINEAR_CANDIDATE')
        parent = adaptive.digest(manifest['parent_sha256'])
        dt = adaptive.finite(manifest['dt_s'], 'attempted dt')
        _require(dt > 0., 'INVALID_REJECTED_DT')
        _require(type(supervision.get('exit_code')) is int and supervision['exit_code'] == 2
            and supervision.get('source_inputs_still_match') is True
            and 'source_error' in supervision and supervision['source_error'] is None
            and supervision.get('interrupted') is False and supervision.get('accepted') is False
            and supervision.get('restart_allowed') is False and supervision.get('error') is None
            and supervision.get('termination_reason') is None, 'UNTRUSTED_OR_INTERRUPTED_SUPERVISION')
        adaptive.digest(supervision['manifest_sha256'])
        grants = supervision['live_history_grants']
        _require(type(grants) is list and len(grants) == 2
            and manifest.get('live_history_grants') == grants, 'MISSING_OR_INCONSISTENT_HISTORY_GRANTS')
        for grant, phase in zip(grants, ('start', 'end')):
            _require(type(grant) is dict and grant.get('status') == 'PASS_LIVE_HISTORY_GRANT'
                and grant.get('phase') == phase and grant.get('history_locked') is True,
                'INVALID_HISTORY_GRANT')
            adaptive.digest(grant['protocol_sha256'])
            _integer(grant, 'controller_pid', 1)
        _require(grants[0]['protocol_sha256'] == grants[1]['protocol_sha256']
            and grants[0]['controller_pid'] == grants[1]['controller_pid'], 'MIXED_HISTORY_GRANTS')
        protocols = [v for k, v in manifest['artifacts_sha256'].items() if k.endswith('/protocol.json')]
        _require(protocols == [grants[0]['protocol_sha256']], 'UNBOUND_GRANT_PROTOCOL')
        failure = reports['secant_failure']
        _require(type(failure) is dict and failure.get('controller') == CONTROLLER
            and failure.get('error') == SECANT_ERROR and failure.get('accepted') is False
            and failure.get('restart_allowed') is False
            and _integer(failure, 'jacobian_builds') == 3
            and _integer(failure, 'jacobian_refreshes') == 2, 'NOT_EXACT_BROYDEN_ARMIJO_REFRESH_FAILURE')
        b, n = reports['broyden_iterations'], reports['newton_iterations']
        b_evaluations = _iterations(b, True)
        n_evaluations = _iterations(n, False)
        _require(_integer(failure, 'evaluations') == b_evaluations
            and all(b[0][key] == n[0][key] for key in METRICS if key != 'wall_seconds'),
            'MIXED_SOLVE_INITIAL_STATE_OR_FAILURE_COUNT')
        _candidate(reports['broyden_candidate'], b[-1], dt)
        _candidate(reports['newton_candidate'], n[-1], dt)
        graph = reports['broyden_graph']
        _require(type(graph) is dict and graph.get('accepted') is False and graph.get('restart_allowed') is False
            and adaptive.digest(graph['graph_fingerprint']) == b[-1]['graph_fingerprint'],
            'MIXED_OR_ACCEPTED_BROYDEN_GRAPH')
        _finite_tree(graph)
        diagnostics = graph['diagnostics']
        _require(diagnostics.get('classification') == 'ENGINEERING_APPROXIMATE_ENDPOINT_VISIBLE_SUPPORT_V1'
            and diagnostics.get('policy') == 'BLOCK_ANY_POSITIVE_STRICT_INTERSECTION_OR_ROUNDOFF_EXPANDED_SOLID_HIT'
            and diagnostics.get('scope') == 'ENDPOINT_CANONICAL_SUPPORT_GRAPH_ONLY'
            and diagnostics.get('old_1nm_guard_not_used_for_filter') is True
            and diagnostics.get('canonical_orientation_preserved') is True
            and diagnostics.get('raw_edge_order_preserved') is True
            and _integer(diagnostics, 'raw_edge_count', 1) == _integer(diagnostics, 'visible_edge_count', 1)
            and all(_integer(diagnostics, key) == 0 for key in ('blocked_edge_count', 'strict_blocked_edge_count',
                'roundoff_only_blocked_edge_count', 'nnps_evaluations', 'physical_solver_evaluations'))
            and diagnostics.get('blocked_edges') == []
            and adaptive.finite(diagnostics['physical_solid_intersection_length_tolerance_m'], 'solid tolerance') == 0.
            and adaptive.finite(diagnostics['maximum_strict_intersection_length_m'], 'strict intersection') == 0.,
            'UNKNOWN_OR_MODIFIED_ENDPOINT_GRAPH_DOMAIN')
        for key in ('particle_id_order_sha256', 'raw_edge_pid_order_sha256', 'raw_topology_sha256',
                'visible_topology_sha256', 'blocked_topology_sha256', 'raw_gradient_sha256',
                'visible_gradient_sha256', 'world_positions_sha256', 'tool_origin_sha256',
                'tool_rotation_sha256', 'geometry_policy_sha256'):
            adaptive.digest(graph['fingerprints'][key])
        diagnosis.update(eligible=True, reason='COMPLETE_FINITE_BROYDEN_AND_NEWTON_ARMIJO_FAILURE_ONLY',
            parent_sha256=parent, attempted_dt_s=dt, scope='FIRST_COARSE_PROPOSAL_ONLY',
            broyden_evaluations=b_evaluations, newton_evaluations=n_evaluations,
            broyden_scaled_residual=b[-1]['scaled_residual'], newton_scaled_residual=n[-1]['scaled_residual'],
            whole_path_certified=False, nonlinear_converged=False, failed_candidate_retained=True)
    except (ValueError, KeyError, TypeError, OverflowError, AttributeError, IndexError) as error:
        diagnosis['reason'] = str(error)
    return diagnosis


def refine(state, diagnosis):
    """Halve actual rejected dt with original bounds and no accepted progress."""
    adaptive.validate_state(state)
    _require(state['status'] == 'RUNNING', 'MATCHING_PENDING_ATTEMPT_REQUIRED')
    pending, config = state['pending'], state['config']
    result = copy.deepcopy(state)
    result.update(status='STOPPED', pending=None, small_error_streak=0, growth_pending=False)
    action, reason = 'STOP', 'UNRECOGNIZED_NONLINEAR_FAILURE_REQUIRES_DIAGNOSIS'
    eligible = (type(diagnosis) is dict and diagnosis.get('protocol') == PROTOCOL
        and diagnosis.get('eligible') is True and diagnosis.get('parent_sha256') == state['parent_sha256']
        and diagnosis.get('attempted_dt_s') == pending['dt_s']
        and diagnosis.get('accepted') is False and diagnosis.get('restart_allowed') is False
        and diagnosis.get('whole_path_certified') is False and diagnosis.get('nonlinear_converged') is False
        and diagnosis.get('failed_candidate_retained') is True
        and diagnosis.get('scope') == 'FIRST_COARSE_PROPOSAL_ONLY'
        and type(diagnosis.get('accepted_time_advanced_s')) in (int, float)
        and diagnosis['accepted_time_advanced_s'] == 0.)
    smaller = pending['dt_s']/2.
    if eligible:
        if state['consecutive_halvings'] >= config['max_halvings']:
            reason = 'NONLINEAR_REFINEMENT_LIMIT'
        elif smaller < max(MINIMUM_DT_S, config['minimum_dt_s']):
            reason = 'MINIMUM_DT_REACHED'
        else:
            result.update(status='READY', nominal_dt_s=smaller,
                consecutive_halvings=state['consecutive_halvings']+1,
                cooldown_remaining=max(state['cooldown_remaining'], COOLDOWN_PASSES, config['cooldown_passes']),
                growth_failures=state['growth_failures']+int(pending['growth_trial']))
            action, reason = 'REFINE', 'NONLINEAR_NONCONVERGENCE_RETRY_UNCHANGED_PARENT'
    adaptive.validate_state(result)
    event = dict(protocol=PROTOCOL, action=action, reason=reason,
        attempt=state['attempts'], parent_before_sha256=state['parent_sha256'],
        parent_after_sha256=result['parent_sha256'], attempted_dt_s=pending['dt_s'],
        nominal_before_s=state['nominal_dt_s'], next_nominal_dt_s=result['nominal_dt_s'],
        phase_clipped=pending['phase_clipped'], growth_trial=pending['growth_trial'],
        small_error_streak=result['small_error_streak'], cooldown_remaining=result['cooldown_remaining'],
        normalized_screen=None, diagnosis=copy.deepcopy(diagnosis),
        diagnostic_only=True, accepted=False, restart_allowed=False, accepted_time_advanced_s=0.)
    return result, event

"""Retry only exhausted geometry certification, by solving smaller physical steps.

This is an administrative policy, not a geometry or physics implementation.
The caller must verify the manifest/supervision/report hashes, select the first
unqualified candidate, require all preceding candidates qualified, and reject
later candidates or an existing pair screen. The four reports passed here are
the original candidate JSON objects, not summaries. A recognized failure stays
FAIL: no candidate, interpolated point, or failed checkpoint becomes a parent.
Sampling without observed events does not certify the rejected whole path.
"""
import copy

import afsd_adaptive_dt_20261010 as adaptive


PROTOCOL = 'UNCHANGED_PARENT_GEOMETRY_QUERY_BUDGET_REFINEMENT_20261010_V1'
GEOMETRY_ERROR = "GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')"
BUDGET_REASON = 'MAXIMUM_ANALYTIC_BATCH_QUERIES'
REPORT_NAMES = ('visible_path', 'visible_shell', 'geometry', 'event_router')
MINIMUM_DT_S = 1e-7
COOLDOWN_PASSES = 8


def _require(condition, reason):
    if not condition:
        raise ValueError(reason)


def _empty_lists(report, fields):
    return all(type(report.get(name)) is list and report[name] == [] for name in fields)


def _zero(report, field):
    return type(report.get(field)) is int and report[field] == 0


def _integer(report, field, minimum=0):
    value = report.get(field)
    _require(type(value) is int and value >= minimum, 'INVALID_INTEGER_'+field)
    return value


def recognize(manifest, supervision, reports):
    """Fail closed unless the sole observed failure is the original 256-query cap.

    reports keys: visible_path, visible_shell, geometry, event_router. Filesystem
    trust and candidate ordering belong to the caller; this function has no I/O.
    Returned parent/dt bind a recognized rejection to one pending policy attempt.
    """
    diagnosis = dict(protocol=PROTOCOL, eligible=False,
        reason='UNRECOGNIZED_GEOMETRY_FAILURE', diagnostic_only=True,
        accepted=False, restart_allowed=False, accepted_time_advanced_s=0.)
    try:
        _require(isinstance(manifest, dict) and isinstance(supervision, dict)
            and isinstance(reports, dict) and set(reports) == set(REPORT_NAMES),
            'INCOMPLETE_OR_UNEXPECTED_REPORT_SET')
        _require(manifest.get('status') == 'FAIL_GPU_STAGE_NOT_ACCEPTED'
            and manifest.get('error') == GEOMETRY_ERROR
            and manifest.get('accepted') is False
            and manifest.get('restart_allowed') is False
            and manifest.get('diagnostic_only') is True
            and type(manifest.get('accepted_time_advanced_s')) in (int, float)
            and manifest['accepted_time_advanced_s'] == 0., 'NOT_EXACT_REJECTED_GEOMETRY_CANDIDATE')
        parent = adaptive.digest(manifest['parent_sha256'])
        dt = adaptive.finite(manifest['dt_s'], 'attempted dt')
        _require(dt > 0., 'INVALID_REJECTED_DT')
        _require(type(supervision.get('exit_code')) is int and supervision['exit_code'] == 2
            and supervision.get('source_inputs_still_match') is True
            and 'source_error' in supervision and supervision['source_error'] is None
            and supervision.get('interrupted') is False
            and supervision.get('accepted') is False
            and supervision.get('restart_allowed') is False
            and supervision.get('error') is None, 'UNTRUSTED_OR_INTERRUPTED_SUPERVISION')
        adaptive.digest(supervision['manifest_sha256'])
        grants = supervision['live_history_grants']
        _require(type(grants) is list and len(grants) == 2
            and manifest.get('live_history_grants') == grants, 'MISSING_OR_INCONSISTENT_HISTORY_GRANTS')
        for grant, phase in zip(grants, ('start', 'end')):
            _require(isinstance(grant, dict) and grant.get('status') == 'PASS_LIVE_HISTORY_GRANT'
                and grant.get('phase') == phase and grant.get('history_locked') is True,
                'INVALID_HISTORY_GRANT')
            adaptive.digest(grant['protocol_sha256'])
            _integer(grant, 'controller_pid', 1)
        _require(grants[0]['protocol_sha256'] == grants[1]['protocol_sha256']
            and grants[0]['controller_pid'] == grants[1]['controller_pid'], 'MIXED_HISTORY_GRANTS')
        artifacts = manifest['artifacts_sha256']
        protocol_hashes = [value for name, value in artifacts.items() if name.endswith('/protocol.json')]
        _require(protocol_hashes == [grants[0]['protocol_sha256']], 'UNBOUND_GRANT_PROTOCOL')

        path, shell, geometry, router = (reports[name] for name in REPORT_NAMES)
        _require(all(isinstance(report, dict) for report in (path, shell, geometry, router)),
            'REPORTS_MUST_BE_ORIGINAL_OBJECTS')
        _require(path.get('classification') == 'FIXED_RETAINED_GRAPH_TIGHT_PATH_FINITE_TOOL_REMOTE_DOMAIN_V1'
            and path.get('status') == 'UNRESOLVED_VISIBLE_SUPPORT_PATH'
            and path.get('failure_kind') == BUDGET_REASON
            and _integer(path, 'maximum_queries') == 256
            and _integer(path, 'analytic_batch_queries') == 256
            and _integer(path, 'maximum_depth') == 32
            and _integer(path, 'deepest_subdivision') <= 32,
            'NOT_ORIGINAL_VISIBLE_PATH_QUERY_BUDGET')
        _require(all(_zero(path, field) for field in ('physical_solver_evaluations', 'nnps_evaluations', 'edge_deletions'))
            and _empty_lists(path, ('found_edges', 'degenerate_edges', 'numerical_domain_failures')),
            'OBSERVED_INTERSECTION_DOMAIN_OR_GRAPH_FAILURE')
        _require(type(path['query_log']) is list and len(path['query_log']) == 256
            and _integer(path, 'strict_geometry_batch_queries')
                + _integer(path, 'enclosure_batch_queries') == 256,
            'INCONSISTENT_PATH_QUERY_ACCOUNTING')
        edges = _integer(path, 'edge_count', 1)
        unresolved = path['unresolved_intervals']
        _require(type(unresolved) is list and unresolved, 'MISSING_UNRESOLVED_INTERVALS')
        for interval in unresolved:
            _require(isinstance(interval, dict) and interval.get('reason') == BUDGET_REASON,
                'NON_BUDGET_UNRESOLVED_INTERVAL')
            a = adaptive.finite(interval['start_fraction'], 'start fraction')
            b = adaptive.finite(interval['end_fraction'], 'end fraction')
            _require(0. <= a < b <= 1. and _integer(interval, 'depth') <= 32,
                'INVALID_UNRESOLVED_INTERVAL')
            indices, ids = interval['edge_indices'], interval['edge_particle_ids']
            _require(type(indices) is list and indices and len(set(indices)) == len(indices)
                and all(type(i) is int and 0 <= i < edges for i in indices)
                and type(ids) is list and len(ids) == len(indices)
                and all(type(pair) is list and len(pair) == 2
                    and all(type(pid) is int for pid in pair) and pair[0] != pair[1] for pair in ids),
                'INVALID_UNRESOLVED_EDGE_IDENTITIES')

        _require(shell.get('status') == 'CERTIFIED_NO_HIDDEN_SHELL_CONTACT'
            and shell.get('failure_kind') is None
            and _empty_lists(shell, ('hidden_contacts', 'unresolved_intervals', 'finite_patch_unresolved_intervals')),
            'SHELL_NOT_COMPLETELY_CERTIFIED')
        _require(_empty_lists(geometry, ('transient_contacts_missed_by_endpoints',))
            and type(geometry.get('new_non_bore_features')) is list
            and set(geometry['new_non_bore_features']) <= {'SHOULDER', 'TOP'},
            'UNALLOWED_CONTACT_GEOMETRY_EVENT')
        _require(router.get('classification') == 'OBSERVED_VISIBLE_SUPPORT_REMOTE_DOMAIN_EVENT_ROUTER_20261007'
            and router.get('status') == 'NO_OBSERVED_EVENT_REQUIRES_WHOLE_PATH_CERTIFICATE'
            and router.get('requires_resolve') is False and router.get('acceptance_allowed') is False
            and router.get('earliest_observed_bracket') is None
            and _empty_lists(router, ('endpoint_changed_edges', 'changed_edges', 'transition_brackets',
                'unresolved_samples', 'proposed_intervals'))
            and all(_zero(router, key) for key in ('physical_solver_evaluations', 'nnps_evaluations',
                'force_evaluations', 'material_updates', 'edge_deletions',
                'start_blocked_edge_count', 'end_blocked_edge_count')),
            'OBSERVED_OR_UNRESOLVED_TOPOLOGY_EVENT')
        _require(router.get('edge_count') == edges
            and router.get('particle_count') == manifest.get('particle_count') == shell.get('particle_count')
            and type(router.get('query_log')) is list and router['query_log']
            and all(isinstance(row, dict) and _zero(row, 'blocked_edge_count') for row in router['query_log']),
            'INCONSISTENT_EVENT_ROUTER')
        path_inputs, router_inputs = path['input_sha256'], router['input_sha256']
        for key in ('world_position0', 'world_position1', 'particle_id_order',
                'origin0', 'origin1', 'rotation0', 'rotation1'):
            _require(adaptive.digest(path_inputs[key]) == adaptive.digest(router_inputs[key]),
                'MIXED_PATH_AND_ROUTER_INPUTS')
        for path_key, router_key in (('retained_edge_pid_order', 'supplied_edge_pid_order'),
                ('retained_topology', 'supplied_topology')):
            _require(adaptive.digest(path_inputs[path_key]) == adaptive.digest(router_inputs[router_key]),
                'MIXED_PATH_AND_ROUTER_TOPOLOGY')
        diagnosis.update(eligible=True, reason='ONLY_UNRESOLVED_ORIGINAL_GEOMETRY_QUERY_BUDGET',
            parent_sha256=parent, attempted_dt_s=dt, original_query_limit=256,
            unresolved_interval_count=len(unresolved), observed_intersections=0,
            whole_path_certified=False, failed_candidate_retained=True)
    except (ValueError, KeyError, TypeError, OverflowError) as error:
        diagnosis['reason'] = str(error)
    return diagnosis


def refine(state, diagnosis):
    """Halve a recognized rejection from the exact same complete parent.

    The original adaptive state remains compatible with begin_attempt and
    resolve_pass. Shared consecutive-halving count never exceeds the original
    eight; a successful complete pair is still required for any advancement.
    """
    adaptive.validate_state(state)
    _require(state['status'] == 'RUNNING', 'MATCHING_PENDING_ATTEMPT_REQUIRED')
    pending, config = state['pending'], state['config']
    result = copy.deepcopy(state)
    result.update(status='STOPPED', pending=None, small_error_streak=0, growth_pending=False)
    action, reason = 'STOP', 'UNRECOGNIZED_GEOMETRY_FAILURE_REQUIRES_DIAGNOSIS'
    eligible = (isinstance(diagnosis, dict) and diagnosis.get('protocol') == PROTOCOL
        and diagnosis.get('eligible') is True and diagnosis.get('parent_sha256') == state['parent_sha256']
        and diagnosis.get('attempted_dt_s') == pending['dt_s']
        and diagnosis.get('accepted') is False and diagnosis.get('restart_allowed') is False
        and diagnosis.get('whole_path_certified') is False)
    smaller = pending['dt_s']/2.
    if eligible:
        if state['consecutive_halvings'] >= config['max_halvings']:
            reason = 'GEOMETRY_REFINEMENT_LIMIT'
        elif smaller < max(MINIMUM_DT_S, config['minimum_dt_s']):
            reason = 'MINIMUM_DT_REACHED'
        else:
            result.update(status='READY', nominal_dt_s=smaller,
                consecutive_halvings=state['consecutive_halvings']+1,
                cooldown_remaining=max(state['cooldown_remaining'], COOLDOWN_PASSES, config['cooldown_passes']),
                growth_failures=state['growth_failures']+int(pending['growth_trial']))
            action, reason = 'REFINE', 'GEOMETRY_BUDGET_RETRY_UNCHANGED_PARENT'
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

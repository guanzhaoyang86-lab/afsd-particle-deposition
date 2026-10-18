"""Administrative smaller-step retry for multiple observed opening events.

The five original reports are required, including the start/end union router:
an endpoint-only graph can omit a closing edge. The caller binds every file to
the original manifest, source map, complete parent, timestep and first failed
candidate, and disallows later candidates or a pair screen. This policy does
not certify event uniqueness, the failed path, or any failed physical state.
It only permits a fresh half-step solve from the unchanged complete parent.
"""
import copy

import afsd_adaptive_dt_20261010 as adaptive


PROTOCOL = 'UNCHANGED_PARENT_MULTIPLE_OPENING_EVENT_REFINEMENT_20261011_V1'
GEOMETRY_ERROR = "GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')"
REPORT_NAMES = ('visible_path', 'visible_shell', 'geometry', 'event_router', 'union_event_router')
MINIMUM_DT_S = 1e-7
COOLDOWN_PASSES = 8
FRACTION_TOLERANCE = 2.**-30


def _require(condition, reason):
    if not condition:
        raise ValueError(reason)


def _integer(report, field, minimum=0):
    value = report.get(field)
    _require(type(value) is int and value >= minimum, 'INVALID_INTEGER_'+field)
    return value


def _empty_lists(report, fields):
    return all(type(report.get(name)) is list and report[name] == [] for name in fields)


def _zero(report, field):
    return type(report.get(field)) is int and report[field] == 0


def _identities(rows, particles, edges):
    _require(type(rows) is list and rows, 'MISSING_EDGE_IDENTITIES')
    answer = {}
    used_indices = set()
    for row in rows:
        _require(isinstance(row, dict), 'INVALID_EDGE_OBJECT')
        i, j = _integer(row, 'i'), _integer(row, 'j')
        a, b = _integer(row, 'particle_id_i'), _integer(row, 'particle_id_j')
        edge = _integer(row, 'edge_index')
        _require(i < j < particles and a != b and edge < edges,
            'INVALID_EDGE_IDENTITY')
        pair = tuple(sorted((a, b)))
        _require(pair not in answer and edge not in used_indices, 'DUPLICATE_EDGE_IDENTITY')
        answer[pair] = row
        used_indices.add(edge)
    return answer


def _same_identity(a, b, *, edge_index=True):
    keys = ('i', 'j', 'particle_id_i', 'particle_id_j') + (('edge_index',) if edge_index else ())
    return all(a[k] == b[k] for k in keys)


def _router(report, found, particles, active):
    _require(report.get('classification') == 'OBSERVED_VISIBLE_SUPPORT_REMOTE_DOMAIN_EVENT_ROUTER_20261007'
        and report.get('status') == 'VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE'
        and report.get('requires_resolve') is True and report.get('acceptance_allowed') is False
        and _empty_lists(report, ('unresolved_samples',))
        and all(_zero(report, key) for key in ('physical_solver_evaluations', 'nnps_evaluations',
            'force_evaluations', 'material_updates', 'edge_deletions')),
        'UNRESOLVED_OR_UNREGISTERED_EVENT_ROUTER')
    n = len(found)
    edges = _integer(report, 'edge_count', 1)
    _require(_integer(report, 'particle_count', 1) == particles
        and _integer(report, 'start_blocked_edge_count') == n
        and _zero(report, 'end_blocked_edge_count')
        and _integer(report, 'sample_intervals') == 16
        and _integer(report, 'maximum_evaluations') == 256
        and adaptive.finite(report['fraction_tolerance'], 'fraction tolerance') == FRACTION_TOLERANCE,
        'UNREGISTERED_EVENT_COUNTS_OR_LIMITS')
    tables = [_identities(report[key], particles, edges) for key in
        ('endpoint_changed_edges', 'changed_edges', 'transition_brackets')]
    _require(all(set(table) == set(found) for table in tables), 'MIXED_OR_UNKNOWN_CHANGED_PAIRS')
    for pair, original in found.items():
        endpoint, changed, bracket = (table[pair] for table in tables)
        _require(_same_identity(original, endpoint, edge_index=active)
            and _same_identity(endpoint, changed) and _same_identity(endpoint, bracket),
            'MIXED_EVENT_IDENTITIES')
        _require(endpoint.get('start_blocked') is True and endpoint.get('end_blocked') is False
            and bracket.get('start_blocked') is True and bracket.get('end_blocked') is False
            and bracket.get('direction') == 'BLOCKED_TO_VISIBLE', 'NON_OPENING_EVENT')
        _require(adaptive.finite(endpoint['start_strict_intersection_length_m'], 'start intersection')
                == original['total_inside_length_m']
            and adaptive.finite(endpoint['end_strict_intersection_length_m'], 'end intersection') == 0.,
            'INCONSISTENT_ENDPOINT_INTERSECTION')
        a = adaptive.finite(bracket['start_fraction'], 'bracket start')
        b = adaptive.finite(bracket['end_fraction'], 'bracket end')
        width = adaptive.finite(bracket['bracket_width_fraction'], 'bracket width')
        _require(0. <= a < b <= 1. and width == b-a and width <= FRACTION_TOLERANCE
            and bracket.get('status') == 'BRACKET_LOCALIZED_NOT_UNIQUE_OR_PHYSICALLY_RESOLVED'
            and bracket.get('single_event_certified') is False
            and bracket.get('physical_substeps_solved') is False, 'UNLOCALIZED_OR_MISREPRESENTED_EVENT')
    brackets = report['transition_brackets']
    earliest = min(brackets, key=lambda row: row['start_fraction'])
    _require(report.get('earliest_observed_bracket') == earliest, 'INCONSISTENT_EARLIEST_BRACKET')
    a, b = earliest['start_fraction'], earliest['end_fraction']
    expected = [dict(start_fraction=lo, end_fraction=hi, role=role,
        status='DIAGNOSTIC_ONLY_NOT_SOLVED_NOT_CERTIFIED') for lo, hi, role in
        ((0., a, 'BEFORE_OBSERVED_EVENT'), (a, b, 'OBSERVED_EVENT_BRACKET'),
         (b, 1., 'AFTER_OBSERVED_EVENT'))]
    _require(report.get('proposed_intervals') == expected, 'UNREGISTERED_PROPOSED_INTERVALS')
    log = report['query_log']
    evaluations = _integer(report, 'geometry_evaluations', 1)
    _require(type(log) is list and len(log) == evaluations <= 256
        and _integer(report, 'analytic_partition_queries') == 2*evaluations,
        'INCONSISTENT_ROUTER_QUERY_ACCOUNTING')
    samples = {}
    for row in log:
        f = adaptive.finite(row['fraction'], 'sample fraction')
        count = _integer(row, 'blocked_edge_count')
        _require(0. <= f <= 1. and count <= n and f not in samples,
            'INVALID_EVENT_SAMPLE')
        samples[f] = count
    _require(samples.get(0.) == n and samples.get(1.) == 0
        and all(k/16 in samples for k in range(17)), 'INCOMPLETE_ORIGINAL_ROUTER_SAMPLES')
    ordered = sorted(samples)
    _require(all(samples[a] >= samples[b] for a, b in zip(ordered, ordered[1:])),
        'OBSERVED_CLOSING_OR_UNKNOWN_EVENT')
    for row in brackets:
        a, b = row['start_fraction'], row['end_fraction']
        _require(a in samples and b in samples and samples[a] > samples[b],
            'UNSUPPORTED_EVENT_BRACKET')
    return tables[2]


def recognize(manifest, supervision, reports):
    """Return eligible=False for anything outside the narrow opening-only case.

    reports contains visible_path, visible_shell, geometry, event_router and
    union_event_router original JSON objects. No filesystem I/O occurs here.
    A positive diagnosis authorizes only administrative refinement, never PASS.
    """
    diagnosis = dict(protocol=PROTOCOL, eligible=False,
        reason='UNRECOGNIZED_MULTI_EVENT_FAILURE', diagnostic_only=True,
        accepted=False, restart_allowed=False, accepted_time_advanced_s=0.)
    try:
        _require(isinstance(manifest, dict) and isinstance(supervision, dict)
            and isinstance(reports, dict) and set(reports) == set(REPORT_NAMES),
            'INCOMPLETE_OR_UNEXPECTED_REPORT_SET')
        _require(manifest.get('status') == 'FAIL_GPU_STAGE_NOT_ACCEPTED'
            and manifest.get('error') == GEOMETRY_ERROR
            and manifest.get('accepted') is False and manifest.get('restart_allowed') is False
            and manifest.get('diagnostic_only') is True
            and manifest.get('screened_continuation_authorized') is False
            and adaptive.finite(manifest['accepted_time_advanced_s'], 'advanced time') == 0.,
            'NOT_EXACT_REJECTED_GEOMETRY_CANDIDATE')
        parent = adaptive.digest(manifest['parent_sha256'])
        dt = adaptive.finite(manifest['dt_s'], 'attempted dt')
        _require(dt > 0., 'INVALID_REJECTED_DT')
        _require(type(supervision.get('exit_code')) is int and supervision['exit_code'] == 2
            and supervision.get('source_inputs_still_match') is True
            and 'source_error' in supervision and supervision['source_error'] is None
            and supervision.get('interrupted') is False
            and supervision.get('accepted') is False and supervision.get('restart_allowed') is False
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
        protocol_hashes = [value for name, value in manifest['artifacts_sha256'].items()
            if name.endswith('/protocol.json')]
        _require(protocol_hashes == [grants[0]['protocol_sha256']], 'UNBOUND_GRANT_PROTOCOL')
        path, shell, geometry, router, union = (reports[name] for name in REPORT_NAMES)
        _require(all(isinstance(report, dict) for report in reports.values()), 'INVALID_REPORT_OBJECT')
        particles = _integer(manifest, 'particle_count', 1)
        edges = _integer(path, 'edge_count', 1)
        _require(path.get('classification') == 'FIXED_RETAINED_GRAPH_CURVATURE_PRISM_SAT_20261010_V1'
            and path.get('status') == 'SOLID_OCCLUSION_FOUND_ON_RETAINED_EDGE_PATH'
            and path.get('failure_kind') == 'STRICT_SOLID_INTERSECTION_ON_FIXED_RETAINED_PATH'
            and _integer(path, 'maximum_queries') == 256 and _integer(path, 'maximum_depth') == 32
            and _integer(path, 'analytic_batch_queries') == 2
            and _integer(path, 'strict_geometry_batch_queries') == 2
            and all(_zero(path, key) for key in ('enclosure_batch_queries', 'prism_batch_queries',
                'deepest_subdivision', 'physical_solver_evaluations', 'nnps_evaluations',
                'edge_deletions', 'certified_edge_interval_count'))
            and _empty_lists(path, ('degenerate_edges', 'numerical_domain_failures',
                'unresolved_intervals', 'certificates'))
            and path.get('evaluated_fractions') == [0., 1.]
            and path.get('query_log') == [dict(kind='strict', fraction=f, edge_count=edges) for f in (0., 1.)],
            'NOT_EXACT_INITIAL_BLOCKAGE_ON_ENDPOINT_GRAPH')
        for expected, value, query in zip((0., 1.), path['evaluated_fractions'], path['query_log']):
            _require(adaptive.finite(value, 'evaluated fraction') == expected
                and adaptive.finite(query['fraction'], 'strict query fraction') == expected
                and _integer(query, 'edge_count', 1) == edges, 'INVALID_INITIAL_PATH_SAMPLE')
        found = _identities(path['found_edges'], particles, edges)
        _require(len(found) >= 2, 'MULTIPLE_OPENING_EVENTS_REQUIRED')
        for row in found.values():
            total = adaptive.finite(row['total_inside_length_m'], 'strict intersection')
            maximum = adaptive.finite(row['maximum_contiguous_inside_length_m'], 'maximum intersection')
            _require(adaptive.finite(row['fraction'], 'found fraction') == 0.
                and 0. < maximum <= total, 'NON_INITIAL_OR_INVALID_STRICT_INTERSECTION')
            intervals = row['strict_intervals']
            _require(type(intervals) is list and intervals, 'MISSING_STRICT_INTERSECTION_INTERVALS')
            for interval in intervals:
                _require(type(interval) is list and len(interval) == 2, 'INVALID_STRICT_INTERVAL')
                a, b = (adaptive.finite(v, 'strict interval') for v in interval)
                _require(0. <= a < b <= 1., 'INVALID_STRICT_INTERVAL')
        _require(shell.get('status') == 'CERTIFIED_NO_HIDDEN_SHELL_CONTACT'
            and shell.get('failure_kind') is None
            and (shell.get('solid_path_valid') is None or shell.get('solid_path_valid') is True)
            and _integer(shell, 'particle_count', 1) == particles
            and _empty_lists(shell, ('hidden_contacts', 'unresolved_intervals', 'finite_patch_unresolved_intervals')),
            'SHELL_NOT_COMPLETELY_CERTIFIED')
        _require(_empty_lists(geometry, ('transient_contacts_missed_by_endpoints',))
            and type(geometry.get('new_non_bore_features')) is list
            and set(geometry['new_non_bore_features']) <= {'SHOULDER', 'TOP'},
            'UNALLOWED_CONTACT_GEOMETRY_EVENT')
        active_brackets = _router(router, found, particles, True)
        union_brackets = _router(union, found, particles, False)
        _require(router['edge_count'] == edges and union['edge_count'] >= edges,
            'INCONSISTENT_ACTIVE_OR_UNION_GRAPH_SIZE')
        for pair in found:
            for key in ('start_fraction', 'end_fraction', 'bracket_width_fraction'):
                _require(active_brackets[pair][key] == union_brackets[pair][key],
                    'MIXED_ACTIVE_AND_UNION_EVENT_BRACKETS')
        inputs = path['input_sha256']
        adaptive.digest(inputs['geometry'])
        for other in (router, union):
            for key in ('world_position0', 'world_position1', 'particle_id_order',
                    'origin0', 'origin1', 'rotation0', 'rotation1'):
                _require(adaptive.digest(inputs[key]) == adaptive.digest(other['input_sha256'][key]),
                    'MIXED_PATH_AND_ROUTER_INPUTS')
            for key in ('supplied_edge_pid_order', 'supplied_topology'):
                adaptive.digest(other['input_sha256'][key])
        for a, b in (('retained_edge_pid_order', 'supplied_edge_pid_order'),
                ('retained_topology', 'supplied_topology')):
            _require(adaptive.digest(inputs[a]) == router['input_sha256'][b],
                'MIXED_ACTIVE_PATH_TOPOLOGY')
        diagnosis.update(eligible=True, reason='MULTIPLE_LOCALIZED_OPENINGS_REQUIRE_SMALLER_PHYSICAL_STEP',
            parent_sha256=parent, attempted_dt_s=dt, observed_opening_event_count=len(found),
            opening_particle_id_pairs=[list(pair) for pair in sorted(found)],
            transition_brackets=copy.deepcopy(router['transition_brackets']),
            original_query_limit=256, whole_path_certified=False, event_uniqueness_certified=False,
            physical_substeps_solved=False, failed_candidate_retained=True)
    except (ValueError, KeyError, TypeError, OverflowError, AttributeError) as error:
        diagnosis['reason'] = str(error)
    return diagnosis


def refine(state, diagnosis):
    """Halve actual rejected dt; never advance a parent or manufacture a PASS."""
    adaptive.validate_state(state)
    _require(state['status'] == 'RUNNING', 'MATCHING_PENDING_ATTEMPT_REQUIRED')
    pending, config = state['pending'], state['config']
    result = copy.deepcopy(state)
    result.update(status='STOPPED', pending=None, small_error_streak=0, growth_pending=False)
    action, reason = 'STOP', 'UNRECOGNIZED_MULTI_EVENT_FAILURE_REQUIRES_DIAGNOSIS'
    eligible = (isinstance(diagnosis, dict) and diagnosis.get('protocol') == PROTOCOL
        and diagnosis.get('eligible') is True and diagnosis.get('parent_sha256') == state['parent_sha256']
        and diagnosis.get('attempted_dt_s') == pending['dt_s']
        and diagnosis.get('accepted') is False and diagnosis.get('restart_allowed') is False
        and diagnosis.get('whole_path_certified') is False
        and diagnosis.get('event_uniqueness_certified') is False
        and diagnosis.get('physical_substeps_solved') is False
        and diagnosis.get('failed_candidate_retained') is True
        and type(diagnosis.get('accepted_time_advanced_s')) in (int, float)
        and diagnosis.get('accepted_time_advanced_s') == 0.)
    smaller = pending['dt_s']/2.
    if eligible:
        if state['consecutive_halvings'] >= config['max_halvings']:
            reason = 'MULTI_EVENT_REFINEMENT_LIMIT'
        elif smaller < max(MINIMUM_DT_S, config['minimum_dt_s']):
            reason = 'MINIMUM_DT_REACHED'
        else:
            result.update(status='READY', nominal_dt_s=smaller,
                consecutive_halvings=state['consecutive_halvings']+1,
                cooldown_remaining=max(state['cooldown_remaining'], COOLDOWN_PASSES, config['cooldown_passes']),
                growth_failures=state['growth_failures']+int(pending['growth_trial']))
            action, reason = 'REFINE', 'MULTIPLE_OPENING_EVENTS_RETRY_UNCHANGED_PARENT'
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

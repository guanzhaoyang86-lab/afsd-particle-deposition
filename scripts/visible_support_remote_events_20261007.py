"""Read-only routing of observed support-visibility changes on a trial chord.

This diagnostic complements, and never replaces, the fixed-graph whole-path
certificate.  It localizes observed changes in the existing conservative
endpoint policy.  It cannot certify absence or uniqueness of hidden events,
does not choose a force/thermal graph, and does not accept a physical step.
"""
import numpy as np

from visible_support import _inputs, _hash_array, _topology_hash, filter_visible_support
from visible_support_remote_domain_20261007 import (
    require_visible_support_numerical_domain, VisibleSupportNumericalDomainError, require_remote_filter_consistency,
)
from visible_support_path_tight import _rotation_path


CLASSIFICATION = 'OBSERVED_VISIBLE_SUPPORT_REMOTE_DOMAIN_EVENT_ROUTER_20261007'
EVENT = 'VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE'
UNRESOLVED = 'UNRESOLVED_VISIBILITY_EVENT_DIAGNOSTIC'
NO_OBSERVED_EVENT = 'NO_OBSERVED_EVENT_REQUIRES_WHOLE_PATH_CERTIFICATE'


def audit_visibility_events(x0, x1, i, j, pose0, pose1, p, particle_ids,
                            sample_intervals=16, fraction_tolerance=2.**-30,
                            maximum_evaluations=256):
    """Locate sampled policy switches, with no mutation or solver evaluation.

    ``i,j`` may be the endpoint-retained graph (diagnosing an old rejection)
    or a canonical union of trial support pairs (finding activation and loss).
    NNPS/support-radius changes are outside this API. Fractions always refer to
    the supplied world chord and shortest SO(3) tool path, not a future resolved
    material trajectory. Bisection keeps opposite policy states at its bracket
    ends; it proves neither a single crossing nor a valid graph on either side.

    All paths need the existing continuous certificate even when
    ``requires_resolve`` is false. ``acceptance_allowed`` is ALWAYS false.
    """
    if (isinstance(sample_intervals, bool) or not isinstance(sample_intervals, int)
            or not 1 <= sample_intervals <= 64):
        raise ValueError('sample_intervals must be an integer in 1..64')
    if (isinstance(maximum_evaluations, bool) or not isinstance(maximum_evaluations, int)
            or not sample_intervals+1 <= maximum_evaluations <= 4096):
        raise ValueError('maximum_evaluations must cover the grid and be <=4096')
    if not np.isfinite(fraction_tolerance) or not 2.**-48 <= fraction_tolerance <= .125:
        raise ValueError('fraction_tolerance must be finite in [2**-48, .125]')
    first = np.asarray(x0, float)
    last = np.asarray(x1, float)
    if first.shape != last.shape:
        raise ValueError('matching endpoint position arrays required')
    parsed0 = _inputs(first, i, j, pose0, p, particle_ids)
    parsed1 = _inputs(last, i, j, pose1, p, particle_ids)
    _, i, j, ids, o0, R0, _, _, _, _, _ = parsed0
    _, _, _, _, o1, R1, _, _, _, _, _ = parsed1
    angle, rotation_at = _rotation_path(R0, R1)
    dx, do = last-first, o1-o0
    pairs = np.column_stack((ids[i], ids[j])).astype('<i8', copy=False)
    report = dict(
        classification=CLASSIFICATION, status=NO_OBSERVED_EVENT,
        scope='SAMPLED_POLICY_EVENTS_ON_SUPPLIED_PAIRS_AND_DECLARED_CHORD_ONLY',
        requires_resolve=False, acceptance_allowed=False,
        physical_solver_evaluations=0, nnps_evaluations=0, force_evaluations=0,
        material_updates=0, edge_deletions=0, particle_count=len(first), edge_count=len(i),
        sample_intervals=sample_intervals, maximum_evaluations=maximum_evaluations,
        fraction_tolerance=fraction_tolerance, geometry_evaluations=0,
        analytic_partition_queries=0, rotation_angle_rad=float(angle),
        endpoint_changed_edges=[], changed_edges=[], transition_brackets=[],
        unresolved_samples=[], query_log=[], earliest_observed_bracket=None,
        proposed_intervals=[], input_sha256=dict(
            world_position0=_hash_array(first), world_position1=_hash_array(last),
            particle_id_order=_hash_array(ids.astype('<i8')),
            supplied_edge_pid_order=_hash_array(pairs), supplied_topology=_topology_hash(pairs),
            origin0=_hash_array(o0), origin1=_hash_array(o1),
            rotation0=_hash_array(R0), rotation1=_hash_array(R1)),
        required_resolution=[
            'Reject the unresolved full candidate; keep the last accepted parent and all old FAIL evidence.',
            'Solve physical substeps again from an accepted parent; do not accept interpolated diagnostic chord points.',
            'Use one declared graph policy in density/correction, material gradient, paired force and heat.',
            'Check geometry, common-graph identity, mass, contact, full energy history and time refinement on newly solved substeps.',
            'Recompute event locations after resolving; a saved candidate bracket is not a prescribed event time.',
        ],
        limitations=[
            'No absence or uniqueness of hidden events is certified by finite sampling/bisection.',
            'Fraction tolerance is diagnostic localization precision, not a relaxed physical intersection threshold.',
            'No edge is deleted and no mechanical, thermal, constitutive or accepted state is changed.',
            'No neighbor-radius/support-membership event is detected unless the caller supplies that pair.',
            'A sampled numerical ambiguity or degenerate segment fails closed.',
            'This router does not implement the coupled topology-event integrator.',
        ])
    cache = {}

    def edge_row(edge):
        return dict(edge_index=int(edge), i=int(i[edge]), j=int(j[edge]),
                    particle_id_i=int(ids[i[edge]]), particle_id_j=int(ids[j[edge]]))

    def sample(fraction):
        fraction = float(fraction)
        if fraction in cache:
            return cache[fraction]
        if report['geometry_evaluations'] >= maximum_evaluations:
            return None
        report['geometry_evaluations'] += 1
        x = first+fraction*dx
        current_pose = dict(origin=o0+fraction*do, rotation=rotation_at(fraction))
        try:
            numerical_domain = require_visible_support_numerical_domain(x, i, j, current_pose, p, ids)
            # The diagnostic passes dummy gradients solely because the frozen
            # filter returns a subset; these never enter mechanics or heat.
            filtered = filter_visible_support(x, i, j, np.zeros((len(i), 3)),
                                              current_pose, p, ids)
            require_remote_filter_consistency(numerical_domain, filtered)
            report['analytic_partition_queries'] += 2
            blocked = np.zeros(len(i), dtype=bool)
            blocked[filtered['blocked_raw_edge_indices']] = True
            strict_length = np.zeros(len(i))
            for row in filtered['diagnostics']['blocked_edges']:
                strict_length[row['raw_edge_index']] = row['strict_total_length_m']
            result = dict(blocked=blocked, strict_length=strict_length)
            report['query_log'].append(dict(fraction=fraction,
                                           blocked_edge_count=int(blocked.sum())))
        except (ValueError, ArithmeticError) as error:
            record = dict(fraction=fraction, reason=str(error))
            if isinstance(error, VisibleSupportNumericalDomainError):
                record['numerical_domain'] = error.report
            report['unresolved_samples'].append(record)
            result = None
        cache[fraction] = result
        return result

    grid = np.linspace(0., 1., sample_intervals+1)
    samples = [sample(f) for f in grid]
    begin, end = samples[0], samples[-1]
    if begin is not None and end is not None:
        report['start_blocked_edge_count'] = int(begin['blocked'].sum())
        report['end_blocked_edge_count'] = int(end['blocked'].sum())
        for edge in np.flatnonzero(begin['blocked'] != end['blocked']):
            row = edge_row(edge)
            row.update(start_blocked=bool(begin['blocked'][edge]),
                       end_blocked=bool(end['blocked'][edge]),
                       start_strict_intersection_length_m=float(begin['strict_length'][edge]),
                       end_strict_intersection_length_m=float(end['strict_length'][edge]))
            report['endpoint_changed_edges'].append(row)
    changed = set()
    for slot, (left, right) in enumerate(zip(samples[:-1], samples[1:])):
        if left is None or right is None:
            continue
        for edge in np.flatnonzero(left['blocked'] != right['blocked']):
            changed.add(int(edge))
            lo, hi = float(grid[slot]), float(grid[slot+1])
            low_state, high_state = bool(left['blocked'][edge]), bool(right['blocked'][edge])
            outcome = 'BRACKET_LOCALIZED_NOT_UNIQUE_OR_PHYSICALLY_RESOLVED'
            while hi-lo > fraction_tolerance:
                mid = .5*(lo+hi)
                middle = sample(mid)
                if middle is None:
                    outcome = ('NUMERICAL_OR_DEGENERATE_SAMPLE_UNRESOLVED'
                               if mid in cache else 'EVALUATION_BUDGET_EXHAUSTED')
                    break
                if bool(middle['blocked'][edge]) == low_state:
                    lo = mid
                else:
                    hi = mid
            row = edge_row(edge)
            row.update(start_fraction=lo, end_fraction=hi,
                       bracket_width_fraction=hi-lo,
                       start_blocked=low_state, end_blocked=high_state,
                       direction='BLOCKED_TO_VISIBLE' if low_state else 'VISIBLE_TO_BLOCKED',
                       status=outcome, single_event_certified=False,
                       physical_substeps_solved=False)
            report['transition_brackets'].append(row)
    report['changed_edges'] = [edge_row(edge) for edge in sorted(changed)]
    if report['transition_brackets']:
        report['status'] = EVENT
        report['requires_resolve'] = True
        earliest = min(report['transition_brackets'], key=lambda r: r['start_fraction'])
        report['earliest_observed_bracket'] = dict(earliest)
        lo, hi = earliest['start_fraction'], earliest['end_fraction']
        for a, b, role in ((0., lo, 'BEFORE_OBSERVED_EVENT'), (lo, hi, 'OBSERVED_EVENT_BRACKET'),
                           (hi, 1., 'AFTER_OBSERVED_EVENT')):
            if b > a:
                report['proposed_intervals'].append(dict(start_fraction=a, end_fraction=b,
                    role=role, status='DIAGNOSTIC_ONLY_NOT_SOLVED_NOT_CERTIFIED'))
    if (report['unresolved_samples'] or any(row['status'] !=
            'BRACKET_LOCALIZED_NOT_UNIQUE_OR_PHYSICALLY_RESOLVED' for row in report['transition_brackets'])):
        report['status'] = UNRESOLVED
        report['requires_resolve'] = True
    return report

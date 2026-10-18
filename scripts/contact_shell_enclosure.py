"""Supplement an unchanged Lipschitz audit with analytic planar enclosures.

The old audit still owns its original depth and finite-patch query limits.
Only a MAXIMUM_GEOMETRY_QUERIES / MAXIMUM_SUBDIVISION_DEPTH remainder may be
upgraded. For tool patches 0..5, distance to the finite patch is at least
distance to its infinite plane. For q(s)=R(s).T*(x(s)-o(s)),
 ||q''|| <= 2*angle*||dx-do|| + angle**2*max(||r0||,||r1||).
Each component lies within M*h*h/8 of its endpoint chord on an interval of
width h. Same-side plane endpoints therefore give a conservative lower bound.
Every accepted subinterval is proved by that bound, never by sampling alone.

No additional closest-point/finite-patch geometry queries, contact forces,
material updates, or history changes occur. Outer/base remainders are left to
the original certificate; unknowns cannot become accepted by this supplement.
"""
from copy import deepcopy

import numpy as np

from contact_shell_audit import audit_contact_shell, CERTIFIED, FOUND, UNRESOLVED
from afsd_backend.pysph.quasistatic_contact import _pose, _rotation_path


CLASSIFICATION = 'LIPSCHITZ_PLUS_FINITE_PATCH_PLANE_CURVATURE_ENCLOSURE_V1'


def audit_contact_shell_enclosure(contact, x0, x1, pose0, pose1,
                                  maximum_depth=32, maximum_queries=256,
                                  maximum_plane_evaluations=256):
    """Run the old audit once, then prove only its recorded unresolved part.

    The old report is preserved verbatim under ``original_audit`` and is
    generated here on the exact same inputs; an unrelated external report
    cannot accidentally qualify a different trajectory. Analytic plane path
    evaluations have a separate explicit cap and use no finite-patch queries.
    """
    if (isinstance(maximum_plane_evaluations, bool) or
            not isinstance(maximum_plane_evaluations, int) or
            not 2 <= maximum_plane_evaluations <= 256):
        raise ValueError('analytic plane path evaluation budget 2..256 required')
    original = audit_contact_shell(contact, x0, x1, pose0, pose1,
        maximum_depth=maximum_depth, maximum_queries=maximum_queries)
    result = deepcopy(original)
    result.update(classification=CLASSIFICATION, original_audit=original,
        additional_finite_patch_queries=0, maximum_plane_evaluations=maximum_plane_evaluations,
        plane_path_evaluations=0, plane_evaluated_fractions=[], planar_certificates=[],
        planar_unresolved_intervals=[], original_failure_kind=original['failure_kind'])
    if original['status'] != UNRESOLVED:
        return result
    if original['failure_kind'] not in ('MAXIMUM_GEOMETRY_QUERIES', 'MAXIMUM_SUBDIVISION_DEPTH'):
        return result
    if not original['unresolved_intervals']:
        result['failure_kind'] = 'MISSING_ORIGINAL_UNRESOLVED_INTERVALS'
        return result
    first, last = np.asarray(x0, float), np.asarray(x1, float)
    origin0, rotation0 = _pose(pose0)
    origin1, rotation1 = _pose(pose1)
    angle, rotation_at = _rotation_path(rotation0, rotation1)
    dx, do = last-first, origin1-origin0
    r0, r1 = first-origin0, last-origin1
    max_radius = np.maximum(np.linalg.norm(r0, axis=1), np.linalg.norm(r1, axis=1))
    scale = np.maximum(max_radius, max(contact.radius, contact.height, contact.particle_radius))
    rotation_defect = max(np.linalg.norm(rotation0.T@rotation0-np.eye(3), ord=np.inf),
                          np.linalg.norm(rotation1.T@rotation1-np.eye(3), ord=np.inf),
                          np.linalg.norm(rotation_at(1.)-rotation1, ord=np.inf))
    arithmetic_factor = 2048*np.finfo(float).eps+8*rotation_defect
    coordinate_scale = np.maximum(scale, np.maximum(np.linalg.norm(first, axis=1), np.linalg.norm(last, axis=1)))
    coordinate_scale = np.maximum(coordinate_scale, max(np.linalg.norm(origin0), np.linalg.norm(origin1)))
    margin = 2048*np.finfo(float).eps*coordinate_scale+8*rotation_defect*scale
    curvature = (2*angle*np.linalg.norm(dx-do, axis=1)+angle**2*max_radius)*(1.+arithmetic_factor)
    result.update(maximum_path_second_derivative_m=float(np.max(curvature, initial=0.)),
                  maximum_planar_roundoff_margin_m=float(np.max(margin, initial=0.)),
                  planar_bound='MIN_SAME_SIDE_ENDPOINT_PLANE_DISTANCE_MINUS_M_H_SQUARED_OVER_8_MINUS_ROUNDOFF')
    planes = ((2, 0.), (2, contact.height), (0, -contact.half_bore),
              (0, contact.half_bore), (1, -contact.half_bore), (1, contact.half_bore))
    cache = {}

    class AnalyticBudget(Exception):
        pass

    def local_at(fraction):
        if fraction in cache:
            return cache[fraction]
        if len(cache) >= maximum_plane_evaluations:
            raise AnalyticBudget
        rotation = rotation0 if fraction == 0. else rotation1 if fraction == 1. else rotation_at(fraction)
        # Same world linear interpolation as the frozen certificate. A linear
        # interpolation of endpoint BODY coordinates would be a different path.
        local = (first+fraction*dx-origin0-fraction*do)@rotation
        cache[fraction] = local
        return local

    def record(i, j, lo, hi, depth, **extra):
        return dict(particle_index=int(i), feature_index=int(j), feature=contact.features[j],
                    start_fraction=float(lo), end_fraction=float(hi), depth=int(depth), **extra)

    stack = []
    for interval in reversed(original['unresolved_intervals']):
        for i, j in reversed(interval['pairs']):
            stack.append((int(i), int(j), float(interval['start_fraction']),
                          float(interval['end_fraction']), int(interval['depth'])))
    while stack:
        i, j, lo, hi, depth = stack.pop()
        if j >= 6:
            result['planar_unresolved_intervals'].append(record(i, j, lo, hi, depth,
                reason='OUTER_OR_BASE_REMAINDER_REQUIRES_ORIGINAL_BOUND'))
            continue
        try:
            q0, q1 = local_at(lo)[i], local_at(hi)[i]
        except AnalyticBudget:
            result['planar_unresolved_intervals'].append(record(i, j, lo, hi, depth, reason='MAXIMUM_ANALYTIC_PLANE_EVALUATIONS'))
            result['planar_unresolved_intervals'].extend(record(*item, reason='MAXIMUM_ANALYTIC_PLANE_EVALUATIONS') for item in stack)
            break
        axis, offset = planes[j]
        signed0, signed1 = float(q0[axis]-offset), float(q1[axis]-offset)
        same_side = (signed0 > 0. and signed1 > 0.) or (signed0 < 0. and signed1 < 0.)
        deviation = float(curvature[i]*(hi-lo)**2/8.)
        lower = min(abs(signed0), abs(signed1))-deviation-float(margin[i]) if same_side else -deviation-float(margin[i])
        evidence = record(i, j, lo, hi, depth, plane_axis=axis, plane_offset_m=float(offset),
            endpoint_signed_plane_distances_m=[signed0, signed1], same_side=bool(same_side),
            second_derivative_bound_m=float(curvature[i]), chord_deviation_bound_m=deviation,
            roundoff_margin_m=float(margin[i]), finite_patch_distance_lower_bound_m=float(lower),
            particle_radius_m=float(contact.particle_radius))
        if same_side and lower > contact.particle_radius:
            result['planar_certificates'].append(evidence)
            continue
        if depth >= maximum_depth:
            result['planar_unresolved_intervals'].append(dict(evidence, reason='MAXIMUM_SUBDIVISION_DEPTH'))
            continue
        middle = .5*(lo+hi)
        stack.append((i, j, middle, hi, depth+1))
        stack.append((i, j, lo, middle, depth+1))
    result['plane_path_evaluations'] = len(cache)
    result['plane_evaluated_fractions'] = sorted(cache)
    if result['planar_unresolved_intervals']:
        result['status'] = UNRESOLVED
        result['failure_kind'] = 'PLANAR_ENCLOSURE_UNRESOLVED'
    else:
        result['status'] = CERTIFIED
        result['failure_kind'] = None
        result['unresolved_intervals'] = []
        result['certified_interval_count'] += len(result['planar_certificates'])
    return result

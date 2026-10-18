"""Generic finite-patch distance lower bounds for declared curved tool paths.

The frozen Lipschitz audit runs first with its unchanged query/depth budgets.
Only its unresolved budget/depth intervals can receive a second certificate.
For q(s)=R(s).T*(x(s)-o(s)), ||q''|| <= M, where
 M=2*angle*||dx-do||+angle**2*max(||x0-o0||,||x1-o1||).
On an interval of width h, every component lies between its endpoint extrema
expanded by M*h*h/8. Expanding further for arithmetic and rotation-endpoint
alias error and rounding outwards gives a box containing the entire path.

Distance from that box to a patch is a lower bound for the path distance:
* Bore walls are exact finite axis-aligned rectangles. Their box distance is
  the norm of the three coordinate-interval gaps (one rectangle axis is fixed).
* Each cap is (closed disk) intersect (outside the open square), on a fixed z
  plane. Dropping the circle gives a union of four planar half-planes. Its box
  distance uses the smaller gap to either outside-square coordinate set and
  the z gap. Dropping the square instead gives a disk-distance bound using the
  minimum box radius. The maximum of these two lower bounds is still safe.
* The outer cylinder is a circle times [0,H]. The box's exact radial range is
  [hypot(gap(x,0),gap(y,0)), hypot(maxabs(x),maxabs(y))]; radial and axial gaps
  give a conservative cylinder-distance bound, with outward radial rounding.
* The base retains the original WORLD-coordinate affine certificate. A tool
  body box must never be applied to that stationary laboratory surface.

All final distance bounds lose an additional arithmetic margin. No force,
temperature, history, or physical dissipation is clipped or changed. A bound
that fails to prove separation is UNRESOLVED, never a sampled absence claim.
There are no additional finite-patch queries and no particle-specific rules.
"""
from copy import deepcopy
import math

import numpy as np

from contact_shell_audit import audit_contact_shell, CERTIFIED, UNRESOLVED
from afsd_backend.pysph.quasistatic_contact import _pose, _rotation_path


CLASSIFICATION = 'LIPSCHITZ_PLUS_FINITE_PATCH_CURVATURE_BOX_ENCLOSURE_V1'


def finite_patch_box_lower_bound(contact, feature, lower, upper):
    """Conservative finite-patch distance for one closed BODY-coordinate box.

    This helper deliberately handles tool features 0..6 only. The base uses a
    different coordinate frame and remains with the original affine check.
    """
    lo, hi = np.asarray(lower, float), np.asarray(upper, float)
    if (lo.shape != (3,) or hi.shape != (3,) or not np.isfinite(lo).all() or
            not np.isfinite(hi).all() or np.any(lo > hi) or
            isinstance(feature, bool) or not isinstance(feature, (int, np.integer)) or not 0 <= feature <= 6):
        raise ValueError('finite ordered (3,) box and tool feature 0..6 required')
    a, R, H = contact.half_bore, contact.radius, contact.height
    scale = max(float(np.max(np.abs(lo))), float(np.max(np.abs(hi))), a, R, H, contact.particle_radius)
    arithmetic_margin = 128*np.finfo(float).eps*scale

    def gap(axis, low, high):
        return max(low-hi[axis], lo[axis]-high, 0.)

    # Radial extrema of the entire box, not merely endpoint point norms.
    radial_min = max(0., float(np.nextafter(math.hypot(gap(0, 0., 0.), gap(1, 0., 0.)), -np.inf)))
    radial_max = float(np.nextafter(math.hypot(max(abs(lo[0]), abs(hi[0])),
                                             max(abs(lo[1]), abs(hi[1]))), np.inf))
    if feature in (0, 1):
        z = 0. if feature == 0 else H
        zgap = gap(2, z, z)
        xoutside = max(a-max(abs(lo[0]), abs(hi[0])), 0.)
        youtside = max(a-max(abs(lo[1]), abs(hi[1])), 0.)
        outside_square_bound = math.hypot(zgap, min(xoutside, youtside))
        disk_bound = math.hypot(zgap, max(radial_min-R, 0.))
        raw = max(outside_square_bound, disk_bound)
    elif feature <= 5:
        axis = 0 if feature <= 3 else 1
        fixed = -a if feature in (2, 4) else a
        raw = math.hypot(gap(axis, fixed, fixed), gap(1-axis, -a, a), gap(2, 0., H))
    else:
        radial_gap = max(radial_min-R, R-radial_max, 0.)
        raw = math.hypot(radial_gap, gap(2, 0., H))
    if not math.isfinite(raw) or not math.isfinite(arithmetic_margin):
        raise ValueError('nonfinite distance enclosure arithmetic')
    # Zero is a mathematically valid lower bound for geometric distance; this
    # does not alter a physical force, energy, temperature, or dissipation.
    return max(0., float(np.nextafter(raw-arithmetic_margin, -np.inf)))


def audit_contact_shell_finite_enclosure(contact, x0, x1, pose0, pose1,
                                        maximum_depth=32, maximum_queries=256,
                                        maximum_path_evaluations=256):
    """Preserve the original audit, then bound every remaining tool feature.

    The original audit is recomputed on these exact inputs, not supplied by a
    caller with possibly unrelated geometry. Its FOUND/invalid-solid outcomes
    cannot be upgraded. The base remains governed by its original affine rule.
    """
    if (isinstance(maximum_path_evaluations, bool) or not isinstance(maximum_path_evaluations, int) or
            not 2 <= maximum_path_evaluations <= 256):
        raise ValueError('analytic path budget 2..256 required')
    original = audit_contact_shell(contact, x0, x1, pose0, pose1,
        maximum_depth=maximum_depth, maximum_queries=maximum_queries)
    result = deepcopy(original)
    result.update(classification=CLASSIFICATION, original_audit=original,
        original_failure_kind=original['failure_kind'], additional_finite_patch_queries=0,
        maximum_path_evaluations=maximum_path_evaluations, analytic_path_evaluations=0,
        analytic_evaluated_fractions=[], finite_patch_certificates=[], finite_patch_unresolved_intervals=[])
    if original['status'] != UNRESOLVED or original['failure_kind'] not in (
            'MAXIMUM_GEOMETRY_QUERIES', 'MAXIMUM_SUBDIVISION_DEPTH'):
        return result
    if not original['unresolved_intervals']:
        result['failure_kind'] = 'MISSING_ORIGINAL_UNRESOLVED_INTERVALS'
        return result
    first, last = np.asarray(x0, float), np.asarray(x1, float)
    origin0, rotation0 = _pose(pose0)
    origin1, rotation1 = _pose(pose1)
    angle, rotation_at = _rotation_path(rotation0, rotation1)
    dx, do = last-first, origin1-origin0
    max_radius = np.maximum(np.linalg.norm(first-origin0, axis=1), np.linalg.norm(last-origin1, axis=1))
    scale = np.maximum(max_radius, max(contact.radius, contact.height, contact.particle_radius))
    rotation_defect = max(np.linalg.norm(rotation0.T@rotation0-np.eye(3), ord=np.inf),
                          np.linalg.norm(rotation1.T@rotation1-np.eye(3), ord=np.inf),
                          np.linalg.norm(rotation_at(1.)-rotation1, ord=np.inf))
    arithmetic_factor = 2048*np.finfo(float).eps+8*rotation_defect
    coordinate_scale = np.maximum(scale, np.maximum(np.linalg.norm(first, axis=1), np.linalg.norm(last, axis=1)))
    coordinate_scale = np.maximum(coordinate_scale, max(np.linalg.norm(origin0), np.linalg.norm(origin1)))
    margin = 2048*np.finfo(float).eps*coordinate_scale+8*rotation_defect*scale
    curvature = (2*angle*np.linalg.norm(dx-do, axis=1)+angle**2*max_radius)*(1.+arithmetic_factor)
    if not np.isfinite(curvature).all() or not np.isfinite(margin).all():
        result['failure_kind'] = 'NONFINITE_CURVATURE_ENCLOSURE'
        return result
    result.update(maximum_path_second_derivative_m=float(np.max(curvature, initial=0.)),
                  maximum_coordinate_roundoff_margin_m=float(np.max(margin, initial=0.)),
                  box_bound='ENDPOINT_COORDINATE_EXTREMA_EXPANDED_BY_M_H_SQUARED_OVER_8_AND_ARITHMETIC')
    cache = {}

    class AnalyticBudget(Exception):
        pass

    def local_at(fraction):
        if fraction in cache:
            return cache[fraction]
        if len(cache) >= maximum_path_evaluations:
            raise AnalyticBudget
        rotation = rotation0 if fraction == 0. else rotation1 if fraction == 1. else rotation_at(fraction)
        # Keep the same actual world interpolation and endpoint aliases as the
        # old geometry. The rotation defect explicitly bounds R(1)-R1 mismatch.
        x = first+fraction*dx
        origin = origin0+fraction*do
        cache[fraction] = (x-origin)@rotation
        return cache[fraction]

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
        if j == 7:
            result['finite_patch_unresolved_intervals'].append(record(i, j, lo, hi, depth,
                reason='BASE_REMAINS_WITH_ORIGINAL_WORLD_AFFINE_BOUND'))
            continue
        try:
            q0, q1 = local_at(lo)[i], local_at(hi)[i]
        except AnalyticBudget:
            result['finite_patch_unresolved_intervals'].append(record(i, j, lo, hi, depth,
                reason='MAXIMUM_ANALYTIC_PATH_EVALUATIONS'))
            result['finite_patch_unresolved_intervals'].extend(record(*item,
                reason='MAXIMUM_ANALYTIC_PATH_EVALUATIONS') for item in stack)
            break
        deviation = float(curvature[i]*(hi-lo)**2/8.)
        lower = np.nextafter(np.minimum(q0, q1)-deviation-margin[i], -np.inf)
        upper = np.nextafter(np.maximum(q0, q1)+deviation+margin[i], np.inf)
        distance_bound = finite_patch_box_lower_bound(contact, j, lower, upper)
        evidence = record(i, j, lo, hi, depth, endpoint_body_positions_m=[q0.tolist(), q1.tolist()],
            second_derivative_bound_m=float(curvature[i]), chord_deviation_bound_m=deviation,
            coordinate_roundoff_margin_m=float(margin[i]), body_box_lower_m=lower.tolist(), body_box_upper_m=upper.tolist(),
            finite_patch_distance_lower_bound_m=distance_bound, particle_radius_m=float(contact.particle_radius))
        if distance_bound > contact.particle_radius:
            result['finite_patch_certificates'].append(evidence)
            continue
        if depth >= maximum_depth:
            result['finite_patch_unresolved_intervals'].append(dict(evidence, reason='MAXIMUM_SUBDIVISION_DEPTH'))
            continue
        middle = .5*(lo+hi)
        stack.append((i, j, middle, hi, depth+1))
        stack.append((i, j, lo, middle, depth+1))
    result['analytic_path_evaluations'] = len(cache)
    result['analytic_evaluated_fractions'] = sorted(cache)
    if result['finite_patch_unresolved_intervals']:
        result['status'] = UNRESOLVED
        result['failure_kind'] = 'FINITE_PATCH_ENCLOSURE_UNRESOLVED'
    else:
        result['status'] = CERTIFIED
        result['failure_kind'] = None
        result['unresolved_intervals'] = []
        result['certified_interval_count'] += len(result['finite_patch_certificates'])
    return result

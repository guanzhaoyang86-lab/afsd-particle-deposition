"""Read-only, bounded certificate against endpoint-hidden penalty contact.

This does not integrate contact, solve equilibrium, or validate fields. For
each finite patch inactive at BOTH endpoints it bounds the actual linear
world-particle/origin and shortest-SO(3) tool path. Distance to a closed fixed
patch is 1-Lipschitz. Unproved intervals are subdivided; exhaustion is never
reported as absence of contact. Active-endpoint pairs are outside this audit.

The frozen solver's independent solid-sweep requirement remains necessary.
The result certifies the declared interpolation only, not the unknown true
motion between time samples. Conservative floating margins are included; this
is floating-point geometry with safety bounds, not interval-arithmetic proof.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))

from afsd_backend.pysph.quasistatic_contact import (
    ContactGeometryError, _pose, _rotation_path)


CERTIFIED = 'CERTIFIED_NO_HIDDEN_SHELL_CONTACT'
FOUND = 'HIDDEN_CONTACT_FOUND'
UNRESOLVED = 'UNRESOLVED'


def audit_contact_shell(contact, x0, x1, pose0, pose1,
                        maximum_depth=32, maximum_queries=256):
    """Return a pure diagnostic; particle indices refer to input row order.

    ``maximum_queries`` includes endpoints and counts full-particle geometry
    evaluations, cached at each dyadic fraction. The three statuses are not
    interchangeable. INVALID_SOLID_PATH is explicitly an UNRESOLVED failure,
    never a clean shell certificate. No particle/history array is mutated.
    """
    first = np.asarray(x0, float)
    last = np.asarray(x1, float)
    if (first.ndim != 2 or first.shape[1] != 3 or last.shape != first.shape or
            not np.isfinite(first).all() or not np.isfinite(last).all()):
        raise ValueError('finite equal (N,3) positions required')
    if (isinstance(maximum_depth, bool) or not isinstance(maximum_depth, int) or
            maximum_depth < 0 or maximum_depth > 32 or
            isinstance(maximum_queries, bool) or not isinstance(maximum_queries, int) or
            maximum_queries < 2 or maximum_queries > 256):
        raise ValueError('depth 0..32 and full-geometry query budget 2..256 required')
    origin0, rotation0 = _pose(pose0)
    origin1, rotation1 = _pose(pose1)
    angle, rotation_at = _rotation_path(rotation0, rotation1)
    n = len(first)
    dx = last-first
    do = origin1-origin0
    radius = float(contact.particle_radius)
    scale = np.maximum(np.maximum(np.linalg.norm(first-origin0, axis=1),
                                   np.linalg.norm(last-origin1, axis=1)),
                       max(contact.radius, contact.height, radius))
    # Account for the admitted rotation-matrix orthonormality tolerance and
    # finite endpoint mismatch of the same inherited Rodrigues interpolation.
    rotation_defect = max(np.linalg.norm(rotation0.T@rotation0-np.eye(3), ord=np.inf),
                          np.linalg.norm(rotation1.T@rotation1-np.eye(3), ord=np.inf),
                          np.linalg.norm(rotation_at(1.)-rotation1, ord=np.inf))
    arithmetic_factor = 2048*np.finfo(float).eps+8*rotation_defect
    coordinate_scale = np.maximum(scale, np.maximum(np.linalg.norm(first, axis=1),
                                                     np.linalg.norm(last, axis=1)))
    coordinate_scale = np.maximum(coordinate_scale,
                                  max(np.linalg.norm(origin0), np.linalg.norm(origin1)))
    margin = 2048*np.finfo(float).eps*coordinate_scale+8*rotation_defect*scale
    raw_speed = (np.linalg.norm(dx-do, axis=1)+angle*np.maximum(
        np.linalg.norm(first-origin0, axis=1), np.linalg.norm(last-origin1, axis=1)))
    speed = raw_speed*(1.+arithmetic_factor)+margin
    velocities = np.repeat(speed[:, None], 8, axis=1)
    velocities[:, 7] = np.abs(dx[:, 2])*(1.+arithmetic_factor)+margin
    stationary = np.zeros((n, 8), bool)
    stationary[:, :7] = ((np.all(dx-do == 0., axis=1)) & (angle == 0.))[:, None]
    stationary[:, 7] = dx[:, 2] == 0.
    cache = {}
    queried_fractions = []
    report = dict(status=UNRESOLVED,
        interpolation='LINEAR_PARTICLE_AND_ORIGIN_SHORTEST_SO3_TOOL_ROTATION',
        features=list(contact.features), particle_count=n,
        scope='PATCH_PAIRS_INACTIVE_AT_BOTH_ENDPOINTS_ONLY',
        maximum_depth=maximum_depth, maximum_queries=maximum_queries,
        geometry_queries=0, deepest_subdivision=0, evaluated_fractions=[],
        endpoint_inactive_pair_count=0, checked_pairs=[], hidden_contacts=[],
        unresolved_intervals=[], certified_interval_count=0,
        solid_path_valid=None, failure_kind=None,
        rotation_angle_rad=float(angle), rotation_defect=float(rotation_defect),
        maximum_clearance_roundoff_margin_m=float(np.max(margin, initial=0.)),
        limitation='Declared interpolation only; active-endpoint contact accuracy and the independent solid sweep are not certified.')

    class QueryBudget(Exception):
        pass

    class InvalidPath(Exception):
        def __init__(self, fraction, particle, reason):
            self.fraction, self.particle, self.reason = fraction, particle, reason

    def distances(fraction):
        if fraction in cache:
            return cache[fraction]
        if len(queried_fractions) >= maximum_queries:
            raise QueryBudget
        queried_fractions.append(fraction)
        x = first+fraction*dx
        origin = origin0+fraction*do
        # Endpoints exactly match the frozen force evaluations. Any small
        # interpolation endpoint mismatch is covered by rotation_defect above.
        rotation = rotation0 if fraction == 0. else rotation1 if fraction == 1. else rotation_at(fraction)
        local = (x-origin)@rotation
        gaps = np.full((n, 8), np.inf)
        try:
            if hasattr(contact, '_patches_batch'):
                gaps[:, :7] = contact._patches_batch(local, reject_active_ties=False)[0]
            else:
                for i, point in enumerate(local):
                    gaps[i, :7] = contact._patches(point, reject_active_ties=False)[0]
        except ContactGeometryError as error:
            # A failed batch has no row metadata. Resolve only this diagnostic;
            # never relax the invalid solid/boundary geometry into a distance.
            particle = None
            for i, point in enumerate(local):
                try:
                    contact._patches(point, reject_active_ties=False)
                except ContactGeometryError:
                    particle = i
                    break
            raise InvalidPath(fraction, particle, str(error))
        if contact.base_z is not None:
            gaps[:, 7] = x[:, 2]-contact.base_z
            below = np.flatnonzero(gaps[:, 7] < 0.)
            if len(below):
                raise InvalidPath(fraction, int(below[0]), 'particle centre below rigid base')
        cache[fraction] = gaps
        return gaps

    def interval_record(lo, hi, pairs, depth):
        return dict(start_fraction=float(lo), end_fraction=float(hi), depth=depth,
                    pairs=[[int(i), int(j)] for i, j in pairs])

    try:
        gaps0, gaps1 = distances(0.), distances(1.)
        needed = (gaps0 >= radius) & (gaps1 >= radius) & np.isfinite(gaps0) & np.isfinite(gaps1)
        pairs = np.argwhere(needed)
        report['endpoint_inactive_pair_count'] = len(pairs)
        report['checked_pairs'] = pairs.tolist()
        stack = [(0., 1., pairs, 0)] if len(pairs) else []
        while stack:
            lo, hi, pairs, depth = stack.pop()
            report['deepest_subdivision'] = max(report['deepest_subdivision'], depth)
            ii, jj = pairs[:, 0], pairs[:, 1]
            clear0 = cache[lo][ii, jj]-radius
            clear1 = cache[hi][ii, jj]-radius
            # A known stationary pair has exactly the same geometry throughout
            # and needs no strictly-positive margin (e.g. a base tangent point).
            certified = stationary[ii, jj] | (np.maximum(clear0, clear1) >
                velocities[ii, jj]*(hi-lo)+margin[ii])
            # The base distance is affine, unlike finite tool patch distances.
            certified |= (jj == 7) & (np.minimum(clear0, clear1) > margin[ii])
            report['certified_interval_count'] += int(np.count_nonzero(certified))
            unresolved = pairs[~certified]
            if not len(unresolved):
                continue
            if depth >= maximum_depth:
                report['unresolved_intervals'].append(interval_record(lo, hi, unresolved, depth))
                report['failure_kind'] = 'MAXIMUM_SUBDIVISION_DEPTH'
                continue
            middle = .5*(lo+hi)
            try:
                middle_gaps = distances(middle)
            except QueryBudget:
                report['unresolved_intervals'].append(interval_record(lo, hi, unresolved, depth))
                report['unresolved_intervals'].extend(interval_record(*entry) for entry in stack)
                report['failure_kind'] = 'MAXIMUM_GEOMETRY_QUERIES'
                break
            mi, mj = unresolved[:, 0], unresolved[:, 1]
            hidden = middle_gaps[mi, mj] < radius
            if np.any(hidden):
                for i, j in unresolved[hidden]:
                    gap = float(middle_gaps[i, j])
                    report['hidden_contacts'].append(dict(particle_index=int(i), feature_index=int(j),
                        feature=contact.features[j], fraction=float(middle), gap_m=gap,
                        clearance_m=gap-radius, within_roundoff_margin=bool(radius-gap <= margin[i])))
                report['status'] = FOUND
                report['failure_kind'] = 'ENDPOINT_INACTIVE_INTERIOR_ACTIVE'
                break
            stack.append((middle, hi, unresolved, depth+1))
            stack.append((lo, middle, unresolved, depth+1))
        else:
            if not report['unresolved_intervals']:
                report['status'] = CERTIFIED
        # This routine only detects invalid solids at queried fractions. It
        # intentionally does not claim a second complete solid-sweep proof.
    except InvalidPath as error:
        report.update(status=UNRESOLVED, failure_kind='INVALID_SOLID_PATH', solid_path_valid=False,
                      invalid_geometry=dict(fraction=float(error.fraction), particle_index=error.particle,
                                            reason=error.reason))
    report['geometry_queries'] = len(queried_fractions)
    report['evaluated_fractions'] = sorted(queried_fractions)
    return report

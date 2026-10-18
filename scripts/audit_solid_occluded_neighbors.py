"""Pure endpoint segment visibility diagnostic for the finite cylindrical tool.

No NNPS, constitutive, force, thermal, or time-step evaluation is performed.
The guard uses strict solid interiors; a separately inset solid is diagnostic.
"""
import numpy as np

MINIMUM_OCCLUDED_LENGTH_M = 1e-9
INTERIOR_CLEARANCE_DIAGNOSTIC_M = 1e-10


def _partitions(start, direction, radius, half_bore, bottom, top):
    n = len(start)
    cuts = np.ones((n, 10), dtype=float)
    cuts[:, 0] = 0.
    # All possible state changes of slab, cylinder, and open square bore.
    for column, (axis, boundary) in enumerate(((0, -half_bore), (0, half_bore),
            (1, -half_bore), (1, half_bore), (2, bottom), (2, top)), start=2):
        root = np.full(n, np.inf)
        np.divide(boundary-start[:, axis], direction[:, axis], out=root,
                  where=direction[:, axis] != 0.)
        cuts[:, column] = np.where((root > 0.) & (root < 1.), root, 1.)
    aa = np.sum(direction[:, :2]**2, axis=1)
    bb = 2.*np.sum(start[:, :2]*direction[:, :2], axis=1)
    cc = np.sum(start[:, :2]**2, axis=1)-radius**2
    disc = bb*bb-4.*aa*cc
    real = (aa > 0.) & (disc >= 0.)
    sqrt_disc = np.zeros(n)
    np.sqrt(disc, out=sqrt_disc, where=real)
    stable_q = -.5*(bb+np.copysign(sqrt_disc, bb))
    root1 = np.full(n, np.inf)
    root2 = np.full(n, np.inf)
    np.divide(stable_q, aa, out=root1, where=real)
    np.divide(cc, stable_q, out=root2, where=real & (stable_q != 0.))
    # At q=0 the only double root is zero; it is already a segment endpoint.
    for column, root in ((8, root1), (9, root2)):
        cuts[:, column] = np.where((root > 0.) & (root < 1.), root, 1.)
    cuts.sort(axis=1)
    spans = np.diff(cuts, axis=1)
    mid = start[:, None, :]+.5*(cuts[:, :-1]+cuts[:, 1:])[:, :, None]*direction[:, None, :]
    inside = ((mid[:, :, 2] > bottom) & (mid[:, :, 2] < top)
              & (np.sum(mid[:, :, :2]**2, axis=2) < radius**2)
              & (np.max(np.abs(mid[:, :, :2]), axis=2) > half_bore))
    lengths = spans*np.linalg.norm(direction, axis=1)[:, None]
    maximum = np.zeros(n)
    running = np.zeros(n)
    for slot in range(spans.shape[1]):
        positive = spans[:, slot] > 0.
        running = np.where(positive, np.where(inside[:, slot], running+lengths[:, slot], 0.), running)
        maximum = np.maximum(maximum, running)
    total = np.sum(np.where(inside, lengths, 0.), axis=1)
    return cuts, inside, maximum, total


def audit_solid_edge_visibility(x, i, j, pose, p):
    """Audit actual endpoint NNPS segments; caller supplies matching raw rows.

    ``pose`` is a dict or object with origin and rotation; local rows are
    (x-origin) @ rotation. ``p`` supplies bore_width, tool_radius, layer,
    tool_top in metres. A strict solid interval longer than 1 nm fails.
    Invalid geometry/indices raise ValueError, so callers must fail closed.
    Output indices refer to the supplied raw particle ordering, not PID order.
    """
    x = np.asarray(x, dtype=float)
    i, j = np.asarray(i), np.asarray(j)
    if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
        raise ValueError('x must be finite N by 3 endpoint positions')
    if i.ndim != 1 or j.shape != i.shape or i.dtype.kind not in 'iu' or j.dtype.kind not in 'iu':
        raise ValueError('i,j must be same-length integer vectors')
    if np.any(i < 0) or np.any(j >= len(x)) or np.any(i >= j):
        raise ValueError('edges require 0 <= i < j < particle count')
    if len(i) and len(np.unique(np.column_stack((i, j)), axis=0)) != len(i):
        raise ValueError('duplicate NNPS edges')
    origin = np.asarray(pose['origin'] if isinstance(pose, dict) else pose.origin, dtype=float)
    rotation = np.asarray(pose['rotation'] if isinstance(pose, dict) else pose.rotation, dtype=float)
    if (origin.shape != (3,) or rotation.shape != (3, 3) or not np.isfinite(origin).all()
            or not np.isfinite(rotation).all() or not np.allclose(rotation.T@rotation, np.eye(3), rtol=0., atol=1e-12)
            or abs(np.linalg.det(rotation)-1.) > 1e-12):
        raise ValueError('pose must be finite proper orthogonal rotation and origin')
    radius = float(p['tool_radius'])
    half_bore = .5*float(p['bore_width'])
    height = float(p['tool_top'])-float(p['layer'])
    delta = INTERIOR_CLEARANCE_DIAGNOSTIC_M
    if not np.isfinite([radius, half_bore, height]).all() or min(radius, half_bore, .5*height) <= delta:
        raise ValueError('invalid finite tool dimensions')
    local = (x-origin)@rotation
    start, direction = local[i], local[j]-local[i]
    cuts, inside, maximum, total = _partitions(start, direction, radius, half_bore, 0., height)
    _, _, inset_max, inset_total = _partitions(start, direction, radius-delta, half_bore+delta, delta, height-delta)
    occluded = maximum > MINIMUM_OCCLUDED_LENGTH_M
    inset_occluded = inset_max > MINIMUM_OCCLUDED_LENGTH_M
    rows = []
    for edge in np.flatnonzero(occluded):
        intervals = []
        for slot in range(inside.shape[1]):
            lo, hi = float(cuts[edge, slot]), float(cuts[edge, slot+1])
            if hi <= lo or not inside[edge, slot]:
                continue
            if intervals and intervals[-1][1] == lo:
                intervals[-1][1] = hi
            else:
                intervals.append([lo, hi])
        rows.append(dict(edge_index=int(edge), i=int(i[edge]), j=int(j[edge]),
            strict_interior_parameter_intervals=intervals,
            maximum_contiguous_inside_length_m=float(maximum[edge]),
            total_inside_length_m=float(total[edge]),
            inset_maximum_contiguous_inside_length_m=float(inset_max[edge]),
            inset_total_inside_length_m=float(inset_total[edge]),
            local_start_m=start[edge].tolist(), local_end_m=(start[edge]+direction[edge]).tolist()))
    return dict(status='FAIL_SOLID_OCCLUDED_ENDPOINT_EDGE' if rows else 'PASS_NO_SOLID_OCCLUDED_ENDPOINT_EDGE',
        scope='ACTUAL_ENDPOINT_NNPS_SEGMENTS_ONLY', edge_count=int(len(i)),
        occluded_edge_count=int(occluded.sum()),
        max_occluded_interval_length_m=float(maximum.max(initial=0.)),
        total_strict_inside_length_m=float(total.sum()),
        minimum_occluded_length_m=MINIMUM_OCCLUDED_LENGTH_M,
        strict_guard_geometric_inset_m=0., inner_clearance_diagnostic_m=delta,
        inset_occluded_edge_count=int(inset_occluded.sum()),
        boundary_sensitive_occluded_edge_count=int(np.count_nonzero(occluded & ~inset_occluded)),
        max_inset_occluded_interval_length_m=float(inset_max.max(initial=0.)),
        occluded_edges=rows,
        limitations=['Endpoint straight segments only; no swept graph or future states certified.',
            'Floating-point analytic intersections; strict inside interval must exceed 1 nm.',
            '0.1 nm inset is a separate sensitivity diagnostic and never suppresses the guard.',
            'Clear visibility does not validate free-surface kernel consistency or physical transport.'])

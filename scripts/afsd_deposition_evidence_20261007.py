"""Read-only, persistent-ID evidence for supplied material and trailing residue.

The caller must separately qualify the physical states. These observations do
not qualify a solver, prove adhesion, or prescribe any material trajectory.
"""
import math
import numpy as np


def summarize_frames(frames, parameters):
    """Summarize chronological saved states; never infer a path between them.

    Frame keys: time_s, ids, position, mass, ep, eligible_ids, crossing_ids,
    tool_origin, tool_rotation. Positive x is the declared straight traverse.
    Trailing candidates must lie entirely behind the shoulder projection and
    have centres between the original substrate top and the shoulder plane.
    Consecutive residence is a sampled geometric observation, not bonding.
    """
    radius = float(parameters['tool_radius'])
    dx = float(parameters['dx'])
    top = float(parameters.get('initial_substrate_top_z_m', 0.))
    if not math.isfinite(radius + dx + top) or min(radius, dx) <= 0:
        raise ValueError('finite positive tool and particle dimensions required')
    if not frames:
        raise ValueError('at least one actual saved frame required')
    identity = None
    previous_time = -math.inf
    cross_ep = {}
    residence = {}
    seen_crossed = set()
    rows = []
    for frame in frames:
        t = float(frame['time_s'])
        ids = np.asarray(frame['ids'])
        x = np.asarray(frame['position'], float)
        mass = np.asarray(frame['mass'], float)
        ep = np.asarray(frame['ep'], float)
        origin = np.asarray(frame['tool_origin'], float)
        rotation = np.asarray(frame['tool_rotation'], float)
        n = len(ids)
        if (not math.isfinite(t) or t < 0 or t <= previous_time or
                ids.dtype.kind not in 'iu' or len(np.unique(ids)) != n or
                x.shape != (n, 3) or mass.shape != (n,) or ep.shape != (n,) or
                origin.shape != (3,) or rotation.shape != (3, 3) or
                not all(np.isfinite(a).all() for a in (x, mass, ep, origin, rotation)) or
                np.any(mass <= 0) or np.any(ep < 0) or
                not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12, rtol=0) or
                abs(float(np.linalg.det(rotation))-1) > 1e-12 or
                not np.allclose(rotation[:, 2], [0., 0., 1.], atol=1e-12, rtol=0)):
            raise ValueError('finite chronological physical arrays and vertical tool required')
        eligible = set(map(int, frame['eligible_ids']))
        crossed = set(map(int, frame['crossing_ids']))
        material_ids = set(map(int, ids))
        if not crossed <= eligible <= material_ids or not seen_crossed <= crossed:
            raise ValueError('persistent eligible passage history required')
        order = np.argsort(ids)
        current_identity = (ids[order].tolist(), mass[order].tolist(), sorted(eligible))
        if identity is None:
            identity = current_identity
        elif current_identity != identity:
            raise ValueError('persistent IDs, per-ID masses or initial eligibility changed')
        supplied = np.isin(ids, sorted(crossed))
        r = x - origin
        local = r @ rotation
        height = (x[:, 2] > top) & (local[:, 2] < 0.)
        outside_square = np.max(np.abs(local[:, :2]), axis=1) >= parameters['bore_width']/2
        under = supplied & height & outside_square & (np.sum(r[:, :2]**2, axis=1) <= radius**2)
        trailing = supplied & height & (r[:, 0] < -radius-dx/2)
        # Lateral extrusion is deliberately excluded from trailing observations.
        records = []
        new_residence = {}
        for k, pid in enumerate(ids):
            pid = int(pid)
            if pid in crossed:
                cross_ep.setdefault(pid, float(ep[k]))
            if trailing[k]:
                start, start_x, count = residence.get(pid, (t, x[k].copy(), 0))
                new_residence[pid] = (start, start_x, count+1)
                records.append(dict(particle_id=pid, mass_kg=float(mass[k]),
                    consecutive_saved_frames=count+1, sampled_residence_s=t-start,
                    displacement_since_first_trailing_sample_m=float(np.linalg.norm(x[k]-start_x)),
                    ep_increment_since_first_saved_crossing=float(ep[k])-cross_ep[pid]))
        residence = new_residence
        def population(mask):
            return dict(particle_ids=ids[mask].astype(int).tolist(), particle_count=int(mask.sum()),
                        mass_kg=math.fsum(mass[mask].tolist()))
        rows.append(dict(time_s=t, tool_origin_m=origin.tolist(),
            supplied_crossed=population(supplied), supplied_under_shoulder=population(under),
            supplied_trailing=population(trailing), trailing_residence=records,
            maximum_ep=float(ep.max(initial=0.)),
            maximum_saved_post_crossing_ep_increment=max(
                (float(ep[k])-cross_ep[int(pid)] for k, pid in enumerate(ids) if int(pid) in crossed),
                default=0.)))
        previous_time = t
        seen_crossed = crossed
    return dict(classification='SAMPLED_SUPPLIED_MATERIAL_DEPOSITION_EVIDENCE_V1',
        physical_solver_evaluations=0, frame_count=len(rows), history=rows, final=rows[-1],
        first_supplied_shoulder_time_s=next((r['time_s'] for r in rows
            if r['supplied_under_shoulder']['particle_count']), None),
        first_supplied_trailing_time_s=next((r['time_s'] for r in rows
            if r['supplied_trailing']['particle_count']), None),
        solver_qualification='CALLER_MUST_VERIFY_SAVED_PHYSICAL_EVIDENCE',
        metallurgical_bonding_proven=False, complete_deposition_proven=False,
        limitations=[
            'Saved endpoints do not certify residence between samples.',
            'Initial material below the outlet is excluded by immutable eligibility.',
            'Behind the shoulder means the negative world-x half-plane for a positive-x straight pass.',
            'Geometric residence and plastic strain do not establish adhesion or metallurgical bonding.',
            'Post-crossing strain uses the first saved crossed state, not an interpolated event value.',
        ])

"""Passive supplied-material progress from actual, caller-qualified saved states.

This module does not call a solver, modify particles, prescribe a trajectory,
infer unsaved crossings, or certify physical acceptance. Geometry is explicitly
separate from metallurgical bonding. Histories may start at a continuation
checkpoint: their residence observations apply only to the supplied samples.
"""
from copy import deepcopy
import math
import numpy as np


CLASSIFICATION = 'PASSIVE_SAMPLED_SUPPLIED_MATERIAL_PROGRESS_20261010'


def _ids(values, name):
    values = np.asarray(values)
    if values.ndim != 1 or (values.size and values.dtype.kind not in 'iu'):
        raise ValueError(name + ' must be an integer ID vector')
    result = list(map(int, values))
    if len(result) != len(set(result)):
        raise ValueError(name + ' contains duplicate IDs')
    return result


def frame_from_saved(pa, state, pose=None, observer=None):
    """Copy an actual saved particle/state/pose triplet without changing it.

    ``pose`` must be the caller's actual machine pose at ``state['time_s']``.
    If omitted, the exact last saved machine endpoint is used; an initial
    state without machine history requires an explicit pose.
    The optional saved observer is checked against independently recomputed
    endpoint populations in ``summarize_progress``. No solver is imported.
    """
    if pose is None:
        machine = state.get('machine_history', [])
        if not machine or float(machine[-1]['t1_s']) != float(state['time_s']):
            raise ValueError('matching saved machine endpoint or explicit pose required')
        pose = dict(origin=machine[-1]['origin1_m'], rotation=machine[-1]['rotation1'])
    history = state['bore_observation']
    ids = np.asarray(pa.particle_id).copy()
    order = np.argsort(ids)
    components = np.asarray(pa.component_id).copy()
    mass = np.asarray(pa.m, float).copy()
    if (ids[order].astype(int).tolist() != history['particle_ids']
            or components[order].astype(int).tolist() != history['component_ids']
            or mass[order].tolist() != history['masses_kg']
            or float(history['last_time_s']) != float(state['time_s'])
            or float(history['bore_width_m']) != float(state['parameters']['bore_width'])):
        raise ValueError('saved particle identity/time/bore history mismatch')
    frame = dict(time_s=float(state['time_s']), step=int(state['step']), ids=ids,
        component_ids=components, position=np.column_stack((pa.x, pa.y, pa.z)),
        mass=mass, ep=np.asarray(pa.ep, float).copy(),
        eligible_ids=deepcopy(history['eligible_ids']),
        crossing_ids=deepcopy(history['crossing_ids']),
        tool_origin=np.asarray(pose['origin'], float).copy(),
        tool_rotation=np.asarray(pose['rotation'], float).copy())
    if observer is not None:
        frame['saved_observer'] = deepcopy(observer)
    return frame


def _validated(frame, previous_time, identity, previous_crossed):
    t = float(frame['time_s'])
    ids = _ids(frame['ids'], 'particle IDs')
    eligible = set(_ids(frame['eligible_ids'], 'eligible IDs'))
    crossed = set(_ids(frame['crossing_ids'], 'crossing IDs'))
    component = np.asarray(frame['component_ids'])
    x, mass, ep = (np.asarray(frame[k], float) for k in ('position', 'mass', 'ep'))
    origin, rotation = (np.asarray(frame[k], float) for k in ('tool_origin', 'tool_rotation'))
    n = len(ids)
    if (not n or not math.isfinite(t) or t < 0 or t <= previous_time
            or component.shape != (n,) or component.dtype.kind not in 'iu'
            or x.shape != (n, 3) or mass.shape != (n,) or ep.shape != (n,)
            or origin.shape != (3,) or rotation.shape != (3, 3)
            or not all(np.isfinite(v).all() for v in (x, mass, ep, origin, rotation))
            or np.any(mass <= 0) or np.any(ep < 0)
            or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12, rtol=0)
            or abs(float(np.linalg.det(rotation))-1) > 1e-12
            or not np.allclose(rotation[:, 2], [0., 0., 1.], atol=1e-12, rtol=0)):
        raise ValueError('finite chronological physical arrays and vertical-axis pose required')
    if not previous_crossed <= crossed <= eligible <= {pid for pid, c in zip(ids, component) if c == 1}:
        raise ValueError('persistent supplied-material eligibility and crossing history required')
    by_id = sorted(zip(ids, component.astype(int).tolist(), mass.tolist()))
    this_identity = (by_id, sorted(eligible))
    if identity is not None and this_identity != identity:
        raise ValueError('persistent IDs, components, masses or initial eligibility changed')
    return t, np.asarray(ids, dtype=np.int64), x, mass, ep, origin, rotation, eligible, crossed, this_identity


def summarize_progress(frames, parameters):
    """Measure progress over chronological saved endpoints; return JSON data.

    Gap is signed centre height above the outlet for eligible IDs not yet in
    the saved crossing history. It is not clamped to zero. Crossing counts
    come only from saved persistent-ID events, never from gap or geometry.
    Shoulder residence means consecutive *sampled* membership; positions
    between samples and metallurgical bonding remain unproven.
    """
    radius, bore, dx = (float(parameters[k]) for k in ('tool_radius', 'bore_width', 'dx'))
    top = float(parameters.get('initial_substrate_top_z_m', 0.))
    if not np.isfinite([radius, bore, dx, top]).all() or min(radius, bore, dx) <= 0:
        raise ValueError('finite positive dimensions and finite initial substrate top required')
    if not frames:
        raise ValueError('at least one actual saved frame required')
    previous_time, identity = -math.inf, None
    previous_crossed, ever_under, previous_under = set(), set(), set()
    residence, cross_radius, cross_ep = {}, {}, {}
    rows = []
    for frame in frames:
        t, ids, x, mass, ep, origin, rotation, eligible, crossed, identity = _validated(
            frame, previous_time, identity, previous_crossed)
        relative = x-origin
        local = relative @ rotation
        radial = np.linalg.norm(relative[:, :2], axis=1)
        eligible_mask = np.isin(ids, sorted(eligible))
        supplied = np.isin(ids, sorted(crossed))
        pending = eligible_mask & ~supplied
        height = (x[:, 2] > top) & (local[:, 2] < 0.)
        outside_square = np.max(np.abs(local[:, :2]), axis=1) >= bore/2
        under = supplied & height & outside_square & (radial <= radius)
        clearance = supplied & height & (radial > radius+dx/2)

        def population(mask):
            return dict(particle_ids=sorted(map(int, ids[mask])), particle_count=int(mask.sum()),
                        mass_kg=math.fsum(mass[mask].tolist()))

        under_ids = set(map(int, ids[under]))
        ever_under |= under_ids
        pending_indices = np.flatnonzero(pending)
        closest = min(pending_indices, key=lambda k: (local[k, 2], ids[k])) if pending_indices.size else None
        new_residence, residence_rows = {}, []
        for k, pid_value in enumerate(ids):
            pid = int(pid_value)
            if pid in crossed:
                cross_radius.setdefault(pid, float(radial[k]))
                cross_ep.setdefault(pid, float(ep[k]))
            if pid not in under_ids:
                continue
            start, start_x, count = residence.get(pid, (t, x[k].copy(), 0))
            new_residence[pid] = (start, start_x, count+1)
            residence_rows.append(dict(particle_id=pid, mass_kg=float(mass[k]),
                first_consecutive_under_shoulder_sample_s=start,
                consecutive_saved_frames=count+1, sampled_residence_s=t-start,
                displacement_since_first_consecutive_sample_m=float(np.linalg.norm(x[k]-start_x))))
        residence = new_residence
        spread = dict(status='OBSERVED_SAVED_CENTRES' if crossed else 'NO_SUPPLIED_CROSSED_PARTICLES',
            minimum_radius_m=float(radial[supplied].min()) if crossed else None,
            maximum_radius_m=float(radial[supplied].max()) if crossed else None,
            mass_weighted_radius_m=float(np.sum(mass[supplied]*radial[supplied])/math.fsum(mass[supplied].tolist())) if crossed else None,
            world_xy_extent_m=np.ptp(x[supplied, :2], axis=0).tolist() if crossed else None,
            maximum_radius_change_since_first_saved_crossing_m=max(
                (float(radial[k])-cross_radius[int(pid)] for k, pid in enumerate(ids) if int(pid) in crossed), default=None),
            maximum_ep_increment_since_first_saved_crossing=max(
                (float(ep[k])-cross_ep[int(pid)] for k, pid in enumerate(ids) if int(pid) in crossed), default=None))
        row = dict(time_s=t, step=frame.get('step'), tool_origin_m=origin.tolist(),
            eligible=population(eligible_mask), supplied_crossed=population(supplied),
            supplied_not_yet_crossed=population(pending), supplied_under_shoulder=population(under),
            supplied_clear_of_shoulder_projection=population(clearance),
            nearest_eligible_uncrossed_gap_m=float(local[closest, 2]) if closest is not None else None,
            nearest_eligible_uncrossed_particle_id=int(ids[closest]) if closest is not None else None,
            pending_at_or_below_outlet_ids=sorted(map(int, ids[pending & (local[:, 2] <= 0.)])),
            new_crossing_ids_since_previous_sample=sorted(crossed-previous_crossed) if rows else None,
            new_crossing_mass_since_previous_sample_kg=math.fsum(mass[np.isin(ids, sorted(crossed-previous_crossed))].tolist()) if rows else None,
            entered_under_shoulder_ids_since_previous_sample=sorted(under_ids-previous_under) if rows else None,
            left_under_shoulder_ids_since_previous_sample=sorted(previous_under-under_ids) if rows else None,
            ever_observed_under_shoulder_ids=sorted(ever_under),
            ever_observed_under_shoulder_mass_kg=math.fsum(mass[np.isin(ids, sorted(ever_under))].tolist()),
            under_shoulder_residence=sorted(residence_rows, key=lambda r: r['particle_id']),
            radial_spread=spread, geometric_retention_is_metallurgical_bonding=False)
        saved = frame.get('saved_observer')
        if saved is not None:
            if not math.isclose(float(saved['time1_s']), t, rel_tol=0., abs_tol=1e-12):
                raise ValueError('saved observer endpoint time differs')
            for observer_key, row_key in (('eligible', 'eligible'), ('eligible_crossed', 'supplied_crossed'),
                    ('eligible_not_yet_crossed', 'supplied_not_yet_crossed'),
                    ('crossed_under_annular_shoulder', 'supplied_under_shoulder'),
                    ('geometric_retention_candidate', 'supplied_clear_of_shoulder_projection')):
                actual = saved['supplied_material'][observer_key]
                expected = row[row_key]
                if (sorted(actual['particle_ids']) != expected['particle_ids']
                        or actual['particle_count'] != expected['particle_count']
                        or actual['mass_kg'] != expected['mass_kg']):
                    raise ValueError('saved observer population differs: '+observer_key)
            row['saved_observer_matches'] = True
        else:
            row['saved_observer_matches'] = None
        rows.append(row)
        previous_time, previous_crossed, previous_under = t, crossed, under_ids
    return dict(classification=CLASSIFICATION, status='PASSIVE_SAVED_STATE_OBSERVATION_ONLY',
        physical_solver_evaluations=0, material_trajectory_prescribed=False,
        acceptance_changed=False, formal_accepted_time_advanced_s=0.,
        frame_count=len(rows), first_sample_time_s=rows[0]['time_s'], final=rows[-1], history=rows,
        first_observed_supplied_crossing_time_s=next((r['time_s'] for r in rows if r['supplied_crossed']['particle_count']), None),
        first_observed_supplied_under_shoulder_time_s=next((r['time_s'] for r in rows if r['supplied_under_shoulder']['particle_count']), None),
        metallurgical_bonding_proven=False, complete_deposition_proven=False,
        solver_qualification='CALLER_MUST_VERIFY_SAVED_PHYSICAL_EVIDENCE',
        limitations=[
            'Gap zero marks a centre at the outlet plane, not a certified passage event or deposited layer.',
            'Counts use immutable eligible IDs and saved crossing history; initial lower bar is excluded.',
            'Residence is observed only at supplied saved samples and does not certify motion between them.',
            'The first sample is a baseline: historical crossings are not claimed as new events within this window.',
            'Geometry, radial spread and plastic strain do not establish adhesion or metallurgical bonding.',
            'Hot-start, finite-inventory, constitutive and spatial-resolution limitations remain with the physical model.',
        ])

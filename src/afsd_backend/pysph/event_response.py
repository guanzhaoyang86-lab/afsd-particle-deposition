"""Detached candidate-only derivatives of the unchanged physical contact kick.

All forces are in N; both returned Jacobians have units s/kg. This module never
commits a state and never substitutes its response for the physical mobility.
"""
import hashlib
import json
import time

import numpy as np

from .contact_kick import kick


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    return value


def _digest(value):
    payload = json.dumps(_plain(value), sort_keys=True, separators=(',', ':'),
                         allow_nan=False).encode('utf-8')
    return hashlib.sha256(payload).hexdigest()


def _mean(result, dt):
    return sum(float(a['dt_s']) * .5 * (np.asarray(a['v_before']) +
               np.asarray(a['v_after'])) for a in result['intervals']) / dt


def _summary(result, dt, normals, capacities, walls, labels):
    """Retain physical branch and signed work evidence, not step-count branches."""
    endpoint = np.asarray(result['velocity'], float)
    intervals = result['intervals']
    relative_scale = max(float(a['relative_speed_scale_m_s']) for a in intervals)
    zero = 512 * np.finfo(float).eps * relative_scale
    slip = endpoint[None, :] - walls
    slip -= np.sum(slip * normals, axis=1)[:, None] * normals
    speed = np.linalg.norm(slip, axis=1)
    end_force = np.asarray(intervals[-1]['tangential_impulse']) / intervals[-1]['dt_s']
    cone = []
    for label, s, f, limit in zip(labels, speed, end_force, capacities):
        if limit == 0:
            status = 'ZERO_CAPACITY'
        elif s > zero:
            status = 'SLIDING'
        elif np.linalg.norm(f) < limit * (1 - 2e-8):
            status = 'ZERO_SLIP_CONE_INTERIOR'
        else:
            status = 'ZERO_SLIP_CONE_BOUNDARY'
        cone.append([label, status])
    events = []
    for event in result['events']:
        ids = sorted(labels[int(c)] for c in event['contacts'])
        t = float(event['time_s'])
        if events and abs(t - events[-1]['time_s']) <= 1e-9 * dt:
            events[-1]['contacts'] = sorted(set(events[-1]['contacts'] + ids))
        else:
            events.append(dict(time_s=t, contacts=ids, event='ZERO_SLIP'))
    raw = np.concatenate([np.asarray(a['friction_heat_raw_J']) for a in intervals])
    bounds = np.concatenate([np.asarray(a['friction_rounding_bound_J']) for a in intervals])
    rebuilt = np.concatenate([-np.sum(np.asarray(a['tangential_impulse']) *
        np.asarray(a['tangential_relative_midpoint']), axis=1) for a in intervals])
    rebuild_scale = np.concatenate([np.linalg.norm(a['tangential_impulse'], axis=1)*
        np.linalg.norm(a['tangential_relative_midpoint'], axis=1) for a in intervals])
    q = np.asarray(result['friction_heat_J'])
    discrepancy = float(np.max(np.abs(rebuilt - raw), initial=0.))
    rebuild_bounds = bounds + 1024*np.finfo(float).eps*(np.abs(raw)+rebuild_scale)
    qcheck = bool(np.all(raw >= -bounds) and np.all(np.abs(rebuilt-raw) <= rebuild_bounds))
    return dict(result_sha256=_digest(result), velocity_m_s=endpoint.copy(),
        mean_velocity_m_s=_mean(result, dt), friction_raw_by_face_J=q.copy(),
        friction_raw_J=float(np.sum(q)), minimum_raw_interval_friction_J=(
            float(np.min(raw)) if raw.size else None),
        minimum_raw_plus_rounding_bound_J=(float(np.min(raw + bounds)) if raw.size else None),
        maximum_friction_rounding_bound_J=float(np.max(bounds, initial=0.)),
        reconstructed_heat_max_difference_J=discrepancy,
        reconstructed_heat_max_bound_J=float(np.max(rebuild_bounds, initial=0.)),
        reconstructed_heat_max_excess_J=float(np.max(np.abs(rebuilt-raw)-rebuild_bounds, initial=0.)),
        raw_heat_certificate_pass=qcheck,
        kinetic_work_residual_J=float(result['kinetic_work_residual_J']),
        events=events, endpoint_cone=sorted(cone), endpoint_slip_vectors_m_s=slip.copy(),
        endpoint_zero_tolerance_m_s=zero, relative_speed_scale_m_s=relative_scale,
        interval_count=len(intervals), rejected_intervals=int(result['rejected_intervals']))


def _branch_change(anchor, probe, labels):
    if [e['contacts'] for e in anchor['events']] != [e['contacts'] for e in probe['events']]:
        return 'PHYSICAL_ZERO_SLIP_SEQUENCE_CHANGED'
    if anchor['endpoint_cone'] != probe['endpoint_cone']:
        return 'ENDPOINT_ACTIVE_CONE_CHANGED'
    status = dict(anchor['endpoint_cone'])
    for i, label in enumerate(labels):
        if status[label] == 'SLIDING' and float(
                anchor['endpoint_slip_vectors_m_s'][i] @ probe['endpoint_slip_vectors_m_s'][i]) <= 0:
            return 'ENDPOINT_SLIDING_DIRECTION_REVERSED'
    return None


def response(v0, mass, other_force, normals, normal_forces, wall_velocities, mu, dt,
             free_axes=None, prescribed_velocity=None, *, baseline_result=None,
             contact_ids=None, force_step_relative=1e-3, derivative_rtol=2e-3,
             derivative_atol_factor=2e-5, keep_raw_traces=False):
    """Return available B,C or explicit unavailable evidence, always accepted=False.

    h and h/2 central differences use precisely the legacy kick defaults. A
    changed physical branch or inconsistent refinement never yields a fallback.
    """
    started = time.perf_counter()
    v0 = np.array(v0, float, copy=True)
    force = np.array(other_force, float, copy=True)
    n = np.array(normals, float, copy=True).reshape(-1, 3)
    fn = np.array(normal_forces, float, copy=True)
    walls = np.array(wall_velocities, float, copy=True).reshape(-1, 3)
    free = np.ones(3) if free_axes is None else np.array(free_axes, float, copy=True)
    prescribed = v0.copy() if prescribed_velocity is None else np.array(prescribed_velocity, float, copy=True)
    if (v0.shape != (3,) or force.shape != (3,) or free.shape != (3,) or
            prescribed.shape != (3,) or fn.shape != (len(n),) or walls.shape != n.shape):
        raise ValueError('event response input shape mismatch')
    arrays = (v0, force, n, fn, walls, free, prescribed)
    scalars = (mass, mu, dt, force_step_relative, derivative_rtol, derivative_atol_factor)
    if (not all(np.isfinite(a).all() for a in arrays) or not np.isfinite(scalars).all()
            or mass <= 0 or dt <= 0 or mu < 0 or force_step_relative <= 0
            or derivative_rtol < 0 or derivative_atol_factor < 0 or np.any(fn < 0)
            or np.any((free != 0) & (free != 1))):
        raise ValueError('invalid event response physical input or derivative tolerance')
    if len(n) and np.max(np.abs(np.linalg.norm(n, axis=1)-1)) > 1e-12:
        raise ValueError('unit contact normals required')
    if np.linalg.norm((1-free)*(v0-prescribed)) > 1e-12:
        raise ValueError('inconsistent prescribed initial velocity')
    labels = [str(i) for i in range(len(n))] if contact_ids is None else list(map(str, contact_ids))
    if len(labels) != len(n) or len(set(labels)) != len(labels):
        raise ValueError('unique contact identities required')
    result = dict(available=False, accepted=False, B=None, C=None, units='s/kg',
        scope='NUMERICAL_CANDIDATE_RESPONSE_ONLY', probes=[], kick_calls=0,
        kick_seconds=0., contact_ids=labels, free_axes=free.copy(),
        constrained_columns_exact_zero=np.flatnonzero(free == 0).tolist(),
        derivative_rtol=derivative_rtol, derivative_atol_s_kg=derivative_atol_factor*dt/mass)
    result['fixed_inputs_sha256'] = _digest(dict(v0=v0, mass=mass, normals=n,
        normal_forces=fn, walls=walls, mu=mu, dt=dt, free=free, prescribed=prescribed))

    def finish(reason, available=False):
        result.update(reason=reason, available=available, response_seconds=time.perf_counter()-started)
        return result

    def replay(probe_force, label, axis=None, step=None):
        entry = dict(label=label, force_N=probe_force.copy(), axis=axis, step_N=step)
        result['probes'].append(entry)
        tick = time.perf_counter()
        result['kick_calls'] += 1
        try:
            physical = kick(v0, mass, probe_force, n, fn, walls, mu, dt,
                            free_axes=free, prescribed_velocity=prescribed)
            entry.update(_summary(physical, dt, n, mu*fn, walls, labels))
            if keep_raw_traces:
                entry['raw_kick'] = physical
            entry['status'] = 'PROBE_COMPLETED'
        except Exception as error:
            entry.update(status='PROBE_FAILED', error_type=type(error).__name__, error=str(error),
                cone_context=getattr(error, 'cone_context', None),
                interval_context=getattr(error, 'interval_context', None))
            return None
        finally:
            entry['kick_seconds'] = time.perf_counter()-tick
            result['kick_seconds'] += entry['kick_seconds']
        return entry

    if len(n) == 0:
        endpoint = v0 + dt * free * force / mass
        mean = .5 * (v0 + endpoint)
        result.update(method='ANALYTIC_NO_CONTACT', baseline=dict(velocity_m_s=endpoint,
            mean_velocity_m_s=mean, friction_raw_J=0., friction_raw_by_face_J=np.empty(0)),
            B=np.diag(free)*dt/(2*mass), C=np.diag(free)*dt/mass,
            force_step_N=np.zeros(3))
    else:
        anchor = replay(force, 'baseline')
        if anchor is None:
            return finish('BASELINE_KICK_FAILED')
        result.update(method='ACTUAL_EVENT_CENTRAL_H_AND_HALF_H', baseline=anchor)
        if not anchor['raw_heat_certificate_pass']:
            return finish('BASELINE_RAW_HEAT_CERTIFICATE_FAILED')
    anchor = result['baseline']
    if baseline_result is not None:
        saved_mean = _mean(baseline_result, dt)
        endpoint_difference = float(np.max(np.abs(anchor['velocity_m_s'] - baseline_result['velocity'])))
        mean_difference = float(np.max(np.abs(anchor['mean_velocity_m_s'] - saved_mean)))
        saved_q = np.asarray(baseline_result['friction_heat_J'])
        q_difference = float(np.max(np.abs(anchor['friction_raw_by_face_J'] - saved_q), initial=0.))
        speed_scale = anchor.get('relative_speed_scale_m_s', float(np.linalg.norm(dt*free*force/mass)))
        vbound = 2e-12 + 1e-9*speed_scale
        qbound = 1e-24 + 1e-9*float(np.sum(np.abs(saved_q)))
        match = endpoint_difference <= vbound and mean_difference <= vbound and q_difference <= qbound
        result['saved_replay_comparison'] = dict(passed=bool(match),
            saved_result_sha256=_digest(baseline_result), endpoint_max_difference_m_s=endpoint_difference,
            mean_max_difference_m_s=mean_difference, raw_Q_max_difference_J=q_difference,
            velocity_tolerance_m_s=vbound, raw_Q_tolerance_J=qbound)
        if not match:
            result.update(B=None, C=None)
            return finish('SAVED_JOURNAL_REPLAY_MISMATCH')
    else:
        result['saved_replay_comparison'] = dict(passed=None, reason='NO_SAVED_JOURNAL_PROVIDED')
    if len(n) == 0:
        return finish('ANALYTIC_NO_CONTACT_AVAILABLE', True)
    relative_speed = max(float(np.linalg.norm(v0-walls[0])),
                         float(np.max(np.linalg.norm(walls-walls[0], axis=1), initial=0.)))
    fscale = max(float(np.linalg.norm(force)), float(np.sum(mu*fn)),
                 mass*relative_speed/dt, np.finfo(float).tiny)
    step = force_step_relative*fscale
    result.update(force_scale_N=fscale, force_step_relative=force_step_relative,
                  force_step_N=step*free, column_evidence=[])
    B = np.zeros((3, 3)); C = np.zeros((3, 3))
    for axis in np.flatnonzero(free):
        estimates = []
        for factor in (1., .5):
            pair = []
            for sign in (1., -1.):
                shifted = force.copy(); shifted[axis] += sign*factor*step
                if shifted[axis] == force[axis]:
                    return finish('FORCE_STEP_NOT_REPRESENTABLE')
                entry = replay(shifted, 'central', int(axis), sign*factor*step)
                if entry is None:
                    return finish('DERIVATIVE_KICK_FAILED')
                if not entry['raw_heat_certificate_pass']:
                    return finish('DERIVATIVE_RAW_HEAT_CERTIFICATE_FAILED')
                changed = _branch_change(anchor, entry, labels)
                entry['physical_branch_change'] = changed
                if changed:
                    return finish(changed)
                pair.append(entry)
            estimates.append(tuple((pair[0][key]-pair[1][key])/(2*factor*step)
                for key in ('mean_velocity_m_s', 'velocity_m_s')))
        coarse = np.stack(estimates[0]); fine = np.stack(estimates[1])
        error = np.abs(coarse-fine)
        bound = derivative_rtol*np.maximum(np.abs(coarse), np.abs(fine)) + derivative_atol_factor*dt/mass
        passed = bool(np.isfinite(fine).all() and np.isfinite(coarse).all() and np.all(error <= bound))
        forward = np.stack([(pair[0][key]-anchor[key])/(.5*step)
            for key in ('mean_velocity_m_s', 'velocity_m_s')])
        backward = np.stack([(anchor[key]-pair[1][key])/(.5*step)
            for key in ('mean_velocity_m_s', 'velocity_m_s')])
        sided_error = np.abs(forward-backward)
        sided_bound = derivative_rtol*np.maximum(np.abs(forward), np.abs(backward)) + derivative_atol_factor*dt/mass
        sided_pass = bool(np.isfinite(forward).all() and np.isfinite(backward).all()
                          and np.all(sided_error <= sided_bound))
        result['column_evidence'].append(dict(axis=int(axis), coarse_B_s_kg=coarse[0].copy(),
            fine_B_s_kg=fine[0].copy(), coarse_C_s_kg=coarse[1].copy(), fine_C_s_kg=fine[1].copy(),
            absolute_refinement_difference_s_kg=error, allowed_difference_s_kg=bound,
            refinement_pass=passed, forward_half_h_s_kg=forward, backward_half_h_s_kg=backward,
            one_sided_difference_s_kg=sided_error, one_sided_allowed_s_kg=sided_bound,
            one_sided_pass=sided_pass))
        if not passed:
            return finish('CENTRAL_DIFFERENCE_REFINEMENT_INCONSISTENT')
        if not sided_pass:
            return finish('ONE_SIDED_DIFFERENCES_INCONSISTENT')
        B[:, axis], C[:, axis] = fine
    prescribed_rows = np.flatnonzero(free == 0)
    raw_constrained = np.stack((B[prescribed_rows], C[prescribed_rows]))
    constraint_error = float(np.max(np.abs(raw_constrained), initial=0.))
    result['prescribed_response_certificate'] = dict(raw_derivative_rows_s_kg=raw_constrained,
        maximum_magnitude_s_kg=constraint_error, tolerance_s_kg=128*np.finfo(float).eps*dt/mass)
    if constraint_error > 128*np.finfo(float).eps*dt/mass:
        return finish('PRESCRIBED_OUTPUT_DEPENDS_ON_FORCE')
    B[prescribed_rows] = 0.; C[prescribed_rows] = 0.
    result.update(B=B.copy(), C=C.copy())
    return finish('SAME_PHYSICAL_BRANCH_REFINED_RESPONSE_AVAILABLE', True)


def response_from_trace(trace, dt=None, mu=None, **options):
    """Consume the existing coupled-boundary journal without changing it."""
    if mu is None:
        raise ValueError('frozen friction coefficient mu is required')
    interval = float(trace['dt_s']) if dt is None else float(dt)
    if interval != float(trace['dt_s']):
        raise ValueError('response dt must equal the saved physical kick dt')
    faces = trace['contact_rows']
    answer = response(trace['material_velocity_before_m_s'], trace['mass_kg'],
        trace['internal_force_N'], [f['normal'] for f in faces],
        [f['normal_force_N'] for f in faces], [f['wall_velocity'] for f in faces], mu,
        interval, free_axes=trace['free_axes'], prescribed_velocity=trace['prescribed_velocity'],
        baseline_result=trace['result'], contact_ids=[f['contact_id'] for f in faces], **options)
    answer['particle_id'] = int(trace['particle_id'])
    answer['particle_index'] = int(trace['particle_index'])
    return answer

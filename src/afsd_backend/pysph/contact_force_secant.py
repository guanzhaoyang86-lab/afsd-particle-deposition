"""Fixed-amplitude numerical force secants of the unchanged contact kick.

Branches may differ. This is explicitly NOT a certified classical Jacobian,
nor a replacement for the physical contact map or its work/heat.
"""
import time
import numpy as np
from .contact_kick import kick
from .event_response import _summary, _branch_change, _digest

CLASSIFICATION = 'CONTACT_FORCE_SECANT_APPROXIMATE_NOT_JACOBIAN'
RELATIVE_STEP = 1e-3


def response(v0, mass, force, normals, normal_forces, walls, mu, dt,
             free_axes=None, prescribed_velocity=None, baseline_result=None,
             contact_ids=None):
    start = time.perf_counter()
    v0 = np.array(v0, float, copy=True)
    force = np.array(force, float, copy=True)
    normals = np.array(normals, float, copy=True).reshape(-1, 3)
    fn = np.array(normal_forces, float, copy=True)
    walls = np.array(walls, float, copy=True).reshape(-1, 3)
    free = np.ones(3) if free_axes is None else np.array(free_axes, float, copy=True)
    prescribed = v0.copy() if prescribed_velocity is None else np.array(prescribed_velocity, float, copy=True)
    nc = len(normals)
    if (v0.shape != (3,) or force.shape != (3,) or free.shape != (3,) or
            prescribed.shape != (3,) or fn.shape != (nc,) or walls.shape != (nc, 3) or
            not all(np.isfinite(a).all() for a in (v0, force, normals, fn, walls, free, prescribed)) or
            not np.isfinite([mass, mu, dt]).all() or mass <= 0 or dt <= 0 or mu < 0 or
            np.any(fn < 0) or np.any((free != 0) & (free != 1)) or
            np.max(abs(np.linalg.norm(normals, axis=1)-1), initial=0.) > 1e-12 or
            np.linalg.norm((1-free)*(v0-prescribed)) > 1e-12):
        raise ValueError('invalid contact secant input')
    labels = list(map(str, range(nc) if contact_ids is None else contact_ids))
    if len(labels) != nc or len(set(labels)) != nc:
        raise ValueError('unique contact IDs required')
    records = []

    def replay(f, axis=None, step=0.):
        result = kick(v0, mass, f, normals, fn, walls, mu, dt,
                      free_axes=free, prescribed_velocity=prescribed)
        summary = _summary(result, dt, normals, mu*fn, walls, labels)
        if not summary['raw_heat_certificate_pass']:
            raise ValueError('secant physical kick raw heat certificate failed')
        summary.update(axis=axis, step_N=step, force_N=f.copy())
        records.append(summary)
        return summary

    if nc:
        anchor = replay(force)
    else:
        end = v0+dt*free*force/mass
        anchor = dict(velocity_m_s=end, mean_velocity_m_s=.5*(v0+end),
                      friction_raw_by_face_J=np.empty(0), friction_raw_J=0.)
    if baseline_result is not None:
        mean = sum(a['dt_s']*.5*(np.asarray(a['v_before'])+np.asarray(a['v_after']))
                   for a in baseline_result['intervals'])/dt
        scale = anchor.get('relative_speed_scale_m_s', np.linalg.norm(dt*free*force/mass))
        if (np.max(abs(anchor['velocity_m_s']-baseline_result['velocity'])) > 2e-12+1e-9*scale or
                np.max(abs(anchor['mean_velocity_m_s']-mean)) > 2e-12+1e-9*scale or
                np.max(abs(anchor['friction_raw_by_face_J']-baseline_result['friction_heat_J']), initial=0.) >
                1e-24+1e-9*np.sum(abs(np.asarray(baseline_result['friction_heat_J'])))):
            raise ValueError('saved contact baseline replay mismatch')
    B = np.zeros((3, 3)); C = B.copy(); changes = []
    scale = 0.; step = 0.
    if not nc:
        B = np.diag(free)*dt/(2*mass); C = 2*B
    else:
        speed = max(np.linalg.norm(v0-walls[0]),
                    np.max(np.linalg.norm(walls-walls[0], axis=1), initial=0.))
        scale = max(np.linalg.norm(force), np.sum(mu*fn), mass*speed/dt, np.finfo(float).tiny)
        step = RELATIVE_STEP*scale
        for axis in np.flatnonzero(free):
            pair = []
            for sign in (1., -1.):
                shifted = force.copy(); shifted[axis] += sign*step
                if shifted[axis] == force[axis]:
                    raise ValueError('unrepresentable force secant')
                probe = replay(shifted, int(axis), sign*step)
                change = _branch_change(anchor, probe, labels)
                probe['branch_change'] = change
                if change:
                    changes.append(dict(axis=int(axis), sign=sign, reason=change))
                pair.append(probe)
            B[:, axis] = (pair[0]['mean_velocity_m_s']-pair[1]['mean_velocity_m_s'])/(2*step)
            C[:, axis] = (pair[0]['velocity_m_s']-pair[1]['velocity_m_s'])/(2*step)
    if not np.isfinite(B).all() or not np.isfinite(C).all() or np.any(B[free == 0] != 0.) or np.any(C[free == 0] != 0.):
        raise ValueError('invalid contact secant or prescribed response')
    return dict(classification=CLASSIFICATION, accepted=False, restart_allowed=False,
        full_tangent=False, B=B, C=C, units='s/kg', force_scale_N=float(scale),
        force_step_N=float(step), force_step_relative=RELATIVE_STEP, branch_changes=changes,
        baseline=anchor, probes=records, kick_calls=len(records),
        wall_seconds=time.perf_counter()-start,
        input_sha256=_digest(dict(v0=v0, mass=mass, force=force, normals=normals,
            fn=fn, walls=walls, mu=mu, dt=dt, free=free, prescribed=prescribed)))


def response_from_trace(trace, mu):
    faces = trace['contact_rows']
    return response(trace['material_velocity_before_m_s'], trace['mass_kg'],
        trace['internal_force_N'], [f['normal'] for f in faces],
        [f['normal_force_N'] for f in faces], [f['wall_velocity'] for f in faces],
        mu, trace['dt_s'], trace['free_axes'], trace['prescribed_velocity'],
        trace['result'], [f['contact_id'] for f in faces])

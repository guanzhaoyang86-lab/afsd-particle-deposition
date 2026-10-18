"""Explicit undamped unilateral normal spring and translating-tool companion.

Same normal law as pinned PySPH RigidBodyWallCollision with en=1 (see the
unchanged WallContact adapter); not a new restitution/adhesion/material model.
Only geometric overlap is positive-part selected. Work and heat are NEVER clipped.
"""
import numpy as np

from .rigid_contact_kinematics import contact_rows, rigid_velocity


def normal_state(gap, radius, stiffness):
    """SI: gap/radius [m], stiffness [N/m], force [N], energy [J]."""
    gap = np.asarray(gap, dtype=float)
    if (not np.isfinite(gap).all() or not np.isfinite(radius) or
            not np.isfinite(stiffness) or radius <= 0. or stiffness <= 0.):
        raise ValueError('finite gap and positive finite radius/stiffness required')
    signed = radius-gap
    active = signed > 0.
    penetration = np.where(active, signed, 0.)
    return dict(signed_penetration=signed, penetration=penetration, active=active,
                normal_force_N=stiffness*penetration,
                potential_J=.5*stiffness*penetration**2)


def tool_state(t, p, lift_start, lift_speed):
    """Piecewise prescribed rigid lift; base stays fixed, no material motion rule.

    At the lift breakpoint the right-hand wall velocity is reported. Geometry
    remains continuous. The moving-tool stepper must resolve this breakpoint.
    """
    if (not np.isfinite([t, lift_start, lift_speed]).all() or
            t < 0. or lift_start < 0. or lift_speed < 0.):
        raise ValueError('finite nonnegative time/lift start/lift speed required')
    if not np.isfinite([p['layer'], p['tool_top'], p['traverse'], p['omega']]).all():
        raise ValueError('nonfinite rigid tool parameter')
    offset = max(t-lift_start, 0.)*lift_speed
    vz = lift_speed if t >= lift_start else 0.
    return dict(layer=p['layer']+offset, tool_top=p['tool_top']+offset,
                translation_velocity=np.array([p['traverse'], 0., vz]),
                origin=np.array([p['traverse']*t, 0., p['layer']+offset]),
                omega=np.array([0., 0., p['omega']]), offset=offset)


def moving_tool_rows(x, ids, t, p, lift_start, lift_speed):
    """Existing centre-collocated rows at lifted geometry, with true rigid speed.

    The face velocity is still diagnostic: no new face-force reconstruction.
    BASE feature 0 remains fixed. Existing square opening/features are preserved.
    """
    state = tool_state(t, p, lift_start, lift_speed)
    moved = dict(p, layer=state['layer'], tool_top=state['tool_top'])
    rows = contact_rows(x, ids, t, moved)
    for a, feature in enumerate(rows['feature_id']):
        if feature == 0:
            continue
        rows['translation_velocity'][a] = state['translation_velocity']
        rows['body_origin'][a] = state['origin']
        for key, point in (('wall_velocity', 'application_point'),
                           ('wall_velocity_at_face', 'contact_point')):
            rows[key][a] = rigid_velocity(rows[point][a:a+1], state['origin'],
                                         state['translation_velocity'], state['omega'])[0]
    return rows

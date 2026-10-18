"""SI midpoint equation assembly only; NOT an implicit or quasistatic solver.

The caller supplies forces from the same candidate and accepted history.
Contact events must be resolved separately; endpoint work is not friction heat.
No candidate is advanced, accepted, committed or written by this module.
"""
import numpy as np


def _real(value, name, shape=None):
    array = np.asarray(value)
    if array.dtype.kind not in 'iuf':
        raise ValueError(name + ' must contain real numeric values')
    if shape is not None and array.shape != shape:
        raise ValueError(name + ' must have shape ' + str(shape))
    array = np.asarray(array, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError(name + ' must contain finite values')
    return array


def midpoint_residual(mass, x0, x1, v0, v1, dt, internal_force,
                      contact_force, external_force, reaction_force,
                      free_axes, prescribed_velocity):
    """Return seven N-by-3 residual/force arrays, without mutating inputs.

    mass is positive [kg], dt is positive [s]; positions [m], velocities
    [m/s] and forces [N] are strictly N-by-3. free_axes is N-by-3 boolean
    or exact numeric 0/1, NOT a dense projection matrix. Nonzero reaction
    on a free axis remains visible. Rq is diagnostic only; Rm retains mass.

    An interval containing unresolved slip/constraint/feature events is not
    certified by these equations. This interface supplies no heat or energy.
    """
    m = _real(mass, 'mass')
    if m.ndim != 1 or len(m) == 0 or np.any(m <= 0):
        raise ValueError('mass must be a nonempty positive N vector')
    shape = (len(m), 3)
    a = _real(x0, 'x0', shape)
    b = _real(x1, 'x1', shape)
    u = _real(v0, 'v0', shape)
    w = _real(v1, 'v1', shape)
    fi = _real(internal_force, 'internal_force', shape)
    fc = _real(contact_force, 'contact_force', shape)
    fe = _real(external_force, 'external_force', shape)
    reaction = _real(reaction_force, 'reaction_force', shape)
    vp = _real(prescribed_velocity, 'prescribed_velocity', shape)
    h = _real(dt, 'dt', ())
    if h <= 0:
        raise ValueError('dt must be positive')
    p = np.asarray(free_axes)
    if (p.shape != shape or p.dtype.kind not in 'biuf' or
            not np.isfinite(p).all() or not np.all((p == 0) | (p == 1))):
        raise ValueError('free_axes must be N-by-3 boolean or exact numeric 0/1')
    p = p.astype(float, copy=True)
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            physical = fi + fc + fe
            inertial = m[:, None] * ((w - u) / h)
            result = {
                'kinematic_residual_m': b - a - h * (u + w) / 2.,
                'dynamic_residual_N': inertial - (physical + reaction),
                'quasistatic_free_residual_N': -p * physical,
                'constraint_velocity_residual_m_s': (1. - p) * (w - vp),
                'free_reaction_force_N': p * reaction,
                'inertial_force_N': inertial,
                'physical_force_N': physical,
            }
    except FloatingPointError as error:
        raise ValueError('residual arithmetic became nonfinite; no clipping applied') from error
    if not all(np.isfinite(value).all() for value in result.values()):
        raise ValueError('residual arithmetic became nonfinite; no clipping applied')
    return result

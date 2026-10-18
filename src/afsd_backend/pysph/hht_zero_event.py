"""Local zero-slip limits, not an event integrator or a velocity projector.

This component covers one active fixed-normal contact on a free material
particle. Its use of the real mass is the instantaneous Newton balance,
NEVER an independent finite-step contact impulse cap. Shared contacts and
actuator constraints require their coupled mobility/reactions and reject here.
Every returned direction uses the complete Cartesian tangential vector.
"""
import numpy as np

from .hht_path import parameters


ROUNDING_MULTIPLIER = 512.


class ZeroEventError(ValueError):
    """Unsupported scope, invalid input, or an uncertified zero endpoint."""


def _vector(value, name):
    result = np.array(value, dtype=float, copy=True)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ZeroEventError('finite Cartesian '+name+' required')
    return result


def _scalar(value, name, positive=False):
    if np.ndim(value) != 0:
        raise ZeroEventError('scalar '+name+' required')
    result = float(value)
    if not np.isfinite(result) or (result <= 0 if positive else result < 0):
        raise ZeroEventError(('positive' if positive else 'nonnegative')+' finite '+name+' required')
    return result


def zero_event_limits(mass_kg, normal, normal_force_N, friction,
                      other_force_N, wall_acceleration_m_s2,
                      event_slip_m_s, velocity_scale_m_s, incoming_direction,
                      *, other_contact_count=0, free_axes=(1, 1, 1)):
    """Classify static feasibility and compute unchanged-state force limits.

``other_force_N`` excludes THIS contact's normal and tangential forces.
It includes actual internal and externally applied force at the event state.
No other contact or actuator is supported. ``event_slip_m_s`` is already
the complete relative TANGENTIAL velocity, not one selected component.
``velocity_scale_m_s`` must bound the caller's velocities and h*accelerations;
the zero tolerance is512eps*max(1m/s, supplied scale), not a slip speed floor.
The raw slip is returned unchanged; nothing is advanced or snapped to a wall.

With P=I-n*n.T, the required static traction is P*(m*a_wall-F_other).
If its norm exceeds mu*Fn, the force-selected outgoing direction is the
opposite of that required traction. This need not reverse the incoming
direction in every mathematical input; that observation is returned separately.
Feasibility at the Coulomb-disk boundary is non-strict and separately flagged
as roundoff-marginal. A later driver must not infer a robust branch margin
from that flag or claim this algebra integrates a sticking interval.
"""
    mass = _scalar(mass_kg, 'mass_kg', positive=True)
    fn = _scalar(normal_force_N, 'normal_force_N', positive=True)
    mu = _scalar(friction, 'friction')
    scale = _scalar(velocity_scale_m_s, 'velocity_scale_m_s')
    if (isinstance(other_contact_count, (bool, np.bool_)) or
            np.ndim(other_contact_count) != 0 or other_contact_count != 0):
        raise ZeroEventError('disjoint single-contact scope: shared contacts require coupled mobility')
    axes = _vector(free_axes, 'free_axes')
    if not np.array_equal(axes, np.ones(3)):
        raise ZeroEventError('all physical axes must be free; actuator constraints require coupled reactions')
    n = _vector(normal, 'normal')
    incoming = _vector(incoming_direction, 'incoming_direction')
    other = _vector(other_force_N, 'other_force_N')
    wall_a = _vector(wall_acceleration_m_s2, 'wall_acceleration_m_s2')
    slip = _vector(event_slip_m_s, 'event_slip_m_s')
    eps = np.finfo(float).eps
    unit_bound = ROUNDING_MULTIPLIER*eps
    if abs(np.linalg.norm(n)-1.) > unit_bound:
        raise ZeroEventError('unit contact normal required')
    if abs(np.linalg.norm(incoming)-1.) > unit_bound or abs(incoming@n) > unit_bound:
        raise ZeroEventError('unit tangential incoming direction required')
    zero_bound = ROUNDING_MULTIPLIER*eps*max(1., scale)
    if abs(slip@n) > ROUNDING_MULTIPLIER*eps*max(np.linalg.norm(slip), np.finfo(float).tiny):
        raise ZeroEventError('event_slip_m_s must be tangential, not unprojected relative velocity')
    if np.linalg.norm(slip) > zero_bound:
        raise ZeroEventError('full-vector zero-slip residual exceeds kinematic roundoff bound')
    tangent = lambda value: value-(value@n)*n
    required = tangent(mass*wall_a-other)
    required_norm = float(np.linalg.norm(required))
    capacity = mu*fn
    if (not np.isfinite(capacity) or not np.isfinite(required).all() or
            not np.isfinite(required_norm)):
        raise ZeroEventError('nonfinite static-force balance')
    margin = float(capacity-required_norm)
    force_scale = max(capacity, required_norm, np.linalg.norm(other),
                      mass*np.linalg.norm(wall_a), np.finfo(float).tiny)
    if not np.isfinite(force_scale):
        raise ZeroEventError('nonfinite force-balance norm')
    force_bound = ROUNDING_MULTIPLIER*eps*force_scale
    stick = bool(required_norm <= capacity)
    outgoing = None if stick else -required/required_norm
    left_ft = -capacity*incoming
    right_ft = required.copy() if stick else -capacity*outgoing
    left_a = (other+fn*n+left_ft)/mass
    right_a = (other+fn*n+right_ft)/mass
    relative_a = tangent(right_a-wall_a)
    if not all(np.isfinite(value).all() for value in (left_ft, right_ft, left_a, right_a, relative_a)):
        raise ZeroEventError('nonfinite one-sided acceleration')
    return dict(classification='STICK_FEASIBLE' if stick else 'RESTART_SLIDING',
        scope='INSTANTANEOUS_SINGLE_FIXED_NORMAL_FREE_PARTICLE_NOT_COUPLED_EVENT_INTEGRATOR',
        accepted=False, restart_allowed=False, stick_feasible=stick,
        reverse_relative_to_incoming=bool(not stick and outgoing@incoming < 0.),
        incoming_direction=incoming, outgoing_direction=outgoing,
        required_static_traction_N=required, capacity_N=float(capacity),
        static_margin_N=margin, static_force_roundoff_bound_N=float(force_bound),
        feasibility_margin_is_roundoff_marginal=bool(abs(margin) <= force_bound),
        left_tangential_force_N=left_ft, right_tangential_force_N=right_ft,
        left_acceleration_m_s2=left_a, right_acceleration_m_s2=right_a,
        right_relative_tangential_acceleration_m_s2=relative_a,
        event_slip_m_s=slip, zero_slip_bound_m_s=float(zero_bound),
        mass_kg=mass, normal=n, normal_force_N=fn, friction=mu,
        other_force_N=other, wall_acceleration_m_s2=wall_a,
        event_physical_impulse_Ns=np.zeros(3), event_physical_work_J=0.,
        event_friction_heat_J=0., velocity_was_modified=False,
        finite_step_effective_mass_certified=False, numerical_terms_are_material_heat=False)


def one_sided_slip_direction(coefficients, endpoint, tolerance=None):
    """Limit from INSIDE a quadratic-slip interval, with its raw endpoint.

For u(s)=c0+s*c1+s*s*c2, the left-end inside expansion is
c1*epsilon+c2*epsilon**2; the right-end inside expansion is
-(c1+2*c2)*epsilon+c2*epsilon**2. If the computed linear coefficient is
exactly zero, the quadratic coefficient is used. A nonzero but unresolved
linear cancellation rejects, rather than pretending a higher-order touch.
An all-stick path has no sliding direction and is rejected.

The caller must separately certify no earlier full-vector zero and the
normal/tangential-wall contract. A scalar component root is never enough.
``tolerance`` may tighten, never exceed,512eps*max(1,sum(norm(ci)))m/s.
This function does not overwrite c0, c1, c2 or the endpoint velocity.
"""
    c = np.array(coefficients, dtype=float, copy=True)
    if c.shape != (3, 3) or not np.isfinite(c).all():
        raise ZeroEventError('finite (3,3) ascending-power Cartesian slip coefficients required')
    if isinstance(endpoint, (bool, np.bool_)) or np.ndim(endpoint) != 0 or endpoint not in (0, 1):
        raise ZeroEventError('endpoint must be0(left) or1(right)')
    scale = float(np.linalg.norm(c, axis=1).sum())
    maximum_zero_bound = ROUNDING_MULTIPLIER*np.finfo(float).eps*max(1., scale)
    bound = maximum_zero_bound if tolerance is None else _scalar(tolerance, 'tolerance')
    if bound > maximum_zero_bound:
        raise ZeroEventError('zero-event tolerance cannot exceed declared kinematic roundoff scale')
    raw = c[0].copy() if endpoint == 0 else c[0]+c[1]+c[2]
    if np.linalg.norm(raw) > bound:
        raise ZeroEventError('full-vector endpoint is not a certified zero-slip limit')
    first = c[1].copy() if endpoint == 0 else -(c[1]+2*c[2])
    second = c[2].copy()
    # Derivative cancellation is judged at derivative scale, not the1m/s
    # endpoint scale: genuinely tiny but resolved sliding is not discarded.
    derivative_scale = max(float(np.linalg.norm(c[1])+2*np.linalg.norm(c[2])), np.finfo(float).tiny)
    derivative_bound = ROUNDING_MULTIPLIER*np.finfo(float).eps*derivative_scale
    first_norm = float(np.linalg.norm(first))
    if first_norm > derivative_bound:
        leading, order = first, 1
    elif first_norm > 0.:
        raise ZeroEventError('one-sided derivative unresolved at floating-point cancellation scale')
    elif np.linalg.norm(second) > 0.:
        leading, order = second, 2
    else:
        raise ZeroEventError('all-stick path has no one-sided sliding direction')
    norm = float(np.linalg.norm(leading))
    if not np.isfinite(norm) or norm == 0.:
        raise ZeroEventError('unresolved one-sided sliding direction')
    return dict(direction=leading/norm, leading_coefficient=leading,
        leading_order=order, endpoint=int(endpoint), endpoint_slip_m_s=raw,
        zero_slip_bound_m_s=float(bound), derivative_roundoff_bound_m_s=float(derivative_bound),
        coefficients_m_s=c, accepted=False, velocity_was_modified=False,
        whole_interval_no_other_zero_certified=False)


def acceleration_storage_jump(old_acceleration, new_acceleration, mass,
                              last_dt, weight=.9):
    """Signed acceleration-history diagnostic at a zero-duration restart.

Z=.5*h_last**2*(beta-gamma/2)*sum(m*A*A). This is NOT a proof of the
complete nonlinear algorithmic energy or stability, and NEVER material heat.
The reference step remains the accepted LEFT interval's dt on both sides.
No force, position, velocity or thermodynamic history is changed here.
"""
    old = np.array(old_acceleration, dtype=float, copy=True)
    new = np.array(new_acceleration, dtype=float, copy=True)
    masses = np.array(mass, dtype=float, copy=True)
    h = _scalar(last_dt, 'last_dt', positive=True)
    if (old.ndim not in (1, 2) or old.shape[-1] != 3 or new.shape != old.shape or
            masses.shape != old.shape[:-1] or not np.isfinite(old).all() or
            not np.isfinite(new).all() or not np.isfinite(masses).all() or np.any(masses <= 0)):
        raise ZeroEventError('matching finite Cartesian accelerations and positive physical mass required')
    gamma, beta = parameters(weight)
    d = beta-gamma/2
    value = .5*h*h*d*float(np.sum(masses[..., None]*(new-old)*(new+old)))
    if not np.isfinite(value):
        raise ZeroEventError('nonfinite acceleration storage jump')
    return dict(signed_jump_J=value, d=float(d), last_dt_s=h, weight=float(weight),
        definition='0.5*h_last^2*(beta-gamma/2)*sum(m*(Aplus^2-Aminus^2))',
        full_algorithmic_energy_jump_certified=False, numerical_terms_are_material_heat=False,
        accepted=False)

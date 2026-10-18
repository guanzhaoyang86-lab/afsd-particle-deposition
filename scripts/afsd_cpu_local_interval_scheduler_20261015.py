"""Pure administrative subdivision of already requested physical intervals.

This module never solves physics, accepts a candidate, or classifies failures.
The caller must authorize a split from complete failure evidence and may call
``accept`` only after the entire coarse/fine pair has independently passed.
Returned mappings are new objects; their pending boundaries are immutable tuples.
There is no refinement counter and no fixed minimum time step.
"""
import math
from collections.abc import Mapping
from numbers import Real


SCHEMA = 'CPU_LOCAL_PENDING_INTERVALS_20261015_V1'
_FIELDS = {'schema', 'accepted_time_s', 'nominal_dt', 'target_time_s', 'pending_targets'}


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(name + ' must be a finite real number')
    try:
        result = float(value)
    except (ValueError, OverflowError, TypeError) as error:
        raise ValueError(name + ' must be finite') from error
    if not math.isfinite(result) or result < 0. or (positive and result <= 0.):
        raise ValueError(name + ' must be finite and ' + ('positive' if positive else 'nonnegative'))
    return result


def _interval(begin, end):
    """Require that the original two fine steps can still advance time."""
    dt = end - begin
    if not math.isfinite(dt) or dt <= 0.:
        raise ValueError('interval must have a finite positive duration')
    half = dt / 2.
    middle = begin + half
    if half <= 0. or not begin < middle < end or not middle < middle + half:
        raise ValueError('FLOATING_POINT_FINE_INTERVAL_NOT_REPRESENTABLE')
    return dt, middle


def _validated(state):
    if not isinstance(state, Mapping) or set(state) != _FIELDS or state.get('schema') != SCHEMA:
        raise ValueError('invalid pending-interval state schema')
    accepted = _number(state['accepted_time_s'], 'accepted_time_s')
    nominal = _number(state['nominal_dt'], 'nominal_dt', positive=True)
    target = _number(state['target_time_s'], 'target_time_s', positive=True)
    if accepted > target:
        raise ValueError('accepted time exceeds phase target')
    values = state['pending_targets']
    if not isinstance(values, (list, tuple)):
        raise ValueError('pending_targets must be a list or tuple')
    pending = tuple(_number(value, 'pending boundary', positive=True) for value in values)
    prior = accepted
    for end in pending:
        if not prior < end <= target:
            raise ValueError('pending boundaries must be strictly ordered within the phase')
        _interval(prior, end)
        prior = end
    return accepted, nominal, target, pending


def _state(accepted, nominal, target, pending):
    result = dict(schema=SCHEMA, accepted_time_s=accepted, nominal_dt=nominal,
                  target_time_s=target, pending_targets=tuple(pending))
    _validated(result)
    return result


def start(time_s, nominal_dt, target=.22):
    """Create one nominal interval, clipped to the registered phase target.

    Once its pending boundaries are empty, the caller can call ``start`` again
    from the actual verified endpoint using the same nominal_dt. This restores
    the nominal interval after local refinement without growing the registered
    maximum interval or retaining an unnecessarily small global time step.
    """
    begin = _number(time_s, 'time_s')
    nominal = _number(nominal_dt, 'nominal_dt', positive=True)
    phase_end = _number(target, 'target', positive=True)
    if begin > phase_end:
        raise ValueError('start time exceeds phase target')
    if begin == phase_end:
        return _state(begin, nominal, phase_end, ())
    end = min(phase_end, begin + min(nominal, phase_end - begin))
    _interval(begin, end)
    return _state(begin, nominal, phase_end, (end,))


def current_dt(state, accepted_time):
    """Return the next requested dt, only for this exact accepted parent time."""
    accepted, _, _, pending = _validated(state)
    actual = _number(accepted_time, 'accepted_time')
    if actual != accepted:
        raise ValueError('accepted time does not match the pending state parent')
    if not pending:
        raise ValueError('no pending interval; start a new nominal interval if needed')
    return _interval(actual, pending[0])[0]


def split(state, accepted_time):
    """Split only the current interval; retain its original end and siblings.

    A split that cannot support the original fine half steps is rejected by
    floating-point representability, rather than by a counter or fixed dt.
    """
    current_dt(state, accepted_time)
    accepted, nominal, target, pending = _validated(state)
    _, middle = _interval(accepted, pending[0])
    _interval(accepted, middle)
    _interval(middle, pending[0])
    return _state(accepted, nominal, target, (middle,) + pending)


def accept(state, end_time):
    """Consume one boundary after independent verification of its whole pair.

    Saved fine/event endpoints can differ from the scheduled end by roundoff.
    Permit at most eight ULPs, additionally bounded by one eighth of the actual
    interval. This tolerance cannot consume the parent, a subsequent boundary,
    or time beyond the phase. The actual saved time is retained unchanged.
    """
    accepted, nominal, target, pending = _validated(state)
    actual = _number(end_time, 'end_time')
    if not pending:
        raise ValueError('no pending interval to accept')
    expected = pending[0]
    dt, _ = _interval(accepted, expected)
    tolerance = min(8. * max(math.ulp(expected), math.ulp(actual)), dt / 8.)
    if (not accepted < actual <= target or abs(actual - expected) > tolerance
            or (len(pending) > 1 and actual >= pending[1])):
        raise ValueError('verified endpoint does not match exactly one pending interval')
    return _state(actual, nominal, target, pending[1:])

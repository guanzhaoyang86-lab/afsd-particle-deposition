"""Bound an incoming visibility event using actual same-parent physical solves.

The callback owns the fixed parent, target edge identity, physical equations,
full candidate gates and deadline.  A returned router fraction only proposes
a time: it cannot establish either side of a bracket.  This module neither
interpolates particle states nor promotes a diagnostic state to accepted.
"""
import math
import numbers
import numpy as np


BRACKETED = 'BRACKETED_ACTIVATION'
STOP = 'STOP_EVENT_LOCATION'
MAXIMUM_POSITION_DIFFERENCE_M = 5e-8


def _finite_real(value):
    return (not isinstance(value, (bool, np.bool_))
        and isinstance(value, numbers.Real) and math.isfinite(float(value)))


def locate_activation(solve_probe, full_dt, initial_fraction, *, time_tolerance,
        maximum_solves=10, minimum_dt=1e-10,
        maximum_position_difference=MAXIMUM_POSITION_DIFFERENCE_M):
    """Return real blocked/visible endpoints or a bounded, explicit stop.

    ``solve_probe(dt)`` returns ``blocked`` (native bool), ``event_fraction``
    (None or a finite fraction in [0,1]), ``positions`` (finite N by 3), and an
    opaque ``payload``.  Optional event_count, event_kind and outgoing fields
    are checked when present.  Every call must independently restart from the
    SAME parent and classify the SAME incoming event.  Callback exceptions
    propagate unchanged; the callback also owns the wall-clock deadline.

    Payload objects retain their identity and are never copied or inspected.
    Position snapshots are private, read-only evidence.  The full endpoint is
    physically probed first (a callback may reuse a proven same-parent solve).
    No time-zero assumed endpoint can serve as the returned lower endpoint.
    """
    if not callable(solve_probe):
        raise TypeError('solve_probe must be callable')
    for name, value in [('full_dt', full_dt), ('time_tolerance', time_tolerance),
            ('minimum_dt', minimum_dt),
            ('maximum_position_difference', maximum_position_difference)]:
        if not _finite_real(value) or float(value) <= 0.:
            raise ValueError(name+' must be positive finite')
    full_dt, time_tolerance, minimum_dt = map(float,
        (full_dt, time_tolerance, minimum_dt))
    maximum_position_difference = float(maximum_position_difference)
    if full_dt <= minimum_dt:
        raise ValueError('full_dt must exceed minimum_dt')
    if maximum_position_difference > MAXIMUM_POSITION_DIFFERENCE_M:
        raise ValueError('the original position difference bound cannot be enlarged')
    if not _finite_real(initial_fraction) or not 0. <= float(initial_fraction) <= 1.:
        raise ValueError('initial_fraction must be finite in [0,1]')
    if type(maximum_solves) is not int or not 1 <= maximum_solves <= 10:
        raise ValueError('maximum_solves must be an integer in 1..10')

    attempts, cache = [], {}
    lo, hi = None, None
    shape = None
    stop_reason = None

    def finish(reason=None):
        width = None if lo is None or hi is None else hi['dt_s']-lo['dt_s']
        delta = None
        if lo is not None and hi is not None:
            with np.errstate(over='ignore', invalid='ignore'):
                delta = float(np.max(np.linalg.norm(hi['positions']-lo['positions'], axis=1)))
            if not math.isfinite(delta):
                delta = None
        success = (reason is None and lo is not None and hi is not None
            and lo['dt_s'] > 0. and lo['blocked'] is True and hi['blocked'] is False
            and 0. < width <= time_tolerance and delta is not None
            and delta <= maximum_position_difference)
        return dict(status=BRACKETED if success else STOP,
            reason='REAL_ENDPOINT_BRACKET_WITHIN_TIME_AND_POSITION_BOUNDS' if success
                else reason or 'BRACKET_NOT_WITHIN_BOUNDS',
            lo=lo, hi=hi, attempts=attempts, solve_count=len(attempts),
            bracket_width_s=width, maximum_position_difference_m=delta,
            full_dt_s=full_dt, time_tolerance_s=time_tolerance,
            minimum_dt_s=minimum_dt, maximum_solves=maximum_solves,
            position_difference_limit_m=maximum_position_difference,
            accepted=False, restart_allowed=False, diagnostic_only=True,
            interpolated_states=0, assumed_zero_endpoint_used_as_result=False)

    def probe(dt):
        nonlocal shape, stop_reason
        dt = float(dt)
        if not minimum_dt <= dt <= full_dt:
            stop_reason = 'PROPOSED_TIME_OUTSIDE_DECLARED_BOUNDS'
            return None
        if dt in cache:
            return cache[dt]
        if len(attempts) >= maximum_solves:
            stop_reason = 'MAXIMUM_PHYSICAL_SOLVES'
            return None
        # Record the physical call before invoking it; do not catch or convert
        # callback failures into a usable event location.
        entry = dict(dt_s=dt, valid=False)
        attempts.append(entry)
        raw = solve_probe(dt)
        if not isinstance(raw, dict):
            stop_reason = 'PROBE_IS_NOT_A_DICTIONARY'
            return None
        entry['payload'] = raw.get('payload')
        if type(raw.get('blocked')) is not bool:
            stop_reason = 'ENDPOINT_BLOCKED_IS_NOT_A_NATIVE_BOOL'
            return None
        if 'event_fraction' not in raw or 'positions' not in raw or 'payload' not in raw:
            stop_reason = 'PROBE_REQUIRED_FIELD_MISSING'
            return None
        if 'event_count' in raw and (type(raw['event_count']) is not int
                or not 0 <= raw['event_count'] <= 1):
            stop_reason = 'MULTIPLE_OR_INVALID_EVENT_COUNT'
            return None
        if 'outgoing' in raw and (type(raw['outgoing']) is not bool or raw['outgoing']):
            stop_reason = 'OUTGOING_OR_INVALID_EVENT_DIRECTION'
            return None
        if 'event_kind' in raw and raw['event_kind'] not in (None, 'BLOCKED_TO_VISIBLE'):
            stop_reason = 'OUTGOING_OR_UNSUPPORTED_EVENT_KIND'
            return None
        if 'dt_s' in raw and (not _finite_real(raw['dt_s']) or float(raw['dt_s']) != dt):
            stop_reason = 'PROBE_TIME_DOES_NOT_MATCH_REQUEST'
            return None
        fraction = raw['event_fraction']
        if fraction is not None and (not _finite_real(fraction) or not 0. <= float(fraction) <= 1.):
            stop_reason = 'EVENT_FRACTION_OUTSIDE_FINITE_UNIT_INTERVAL'
            return None
        if raw['blocked'] and fraction is not None:
            stop_reason = 'ACTIVATION_REPORTED_AT_BLOCKED_ENDPOINT'
            return None
        if raw.get('event_count') == 0 and fraction is not None:
            stop_reason = 'EVENT_COUNT_AND_FRACTION_DISAGREE'
            return None
        if 'event_count' in raw and raw['event_count'] != (0 if raw['blocked'] else 1):
            stop_reason = 'EVENT_COUNT_AND_ENDPOINT_DISAGREE'
            return None
        try:
            data = np.asarray(raw['positions'])
            if data.dtype.kind not in 'fiu' or data.ndim != 2 or data.shape[1] != 3 or not len(data):
                raise ValueError('finite N by 3 positions required')
            position = np.array(data, dtype=float, copy=True)
            if not np.isfinite(position).all() or (shape is not None and position.shape != shape):
                raise ValueError('finite same-shape positions required')
        except (TypeError, ValueError, OverflowError):
            stop_reason = 'INVALID_OR_INCONSISTENT_POSITIONS'
            return None
        for previous in cache.values():
            if ((raw['blocked'] and not previous['blocked'] and previous['dt_s'] <= dt)
                    or (not raw['blocked'] and previous['blocked'] and previous['dt_s'] >= dt)):
                stop_reason = 'NONMONOTONE_ENDPOINT_CLASSIFICATION'
                return None
        position.setflags(write=False)
        shape = position.shape
        entry.update(valid=True, blocked=raw['blocked'],
            event_fraction=None if fraction is None else float(fraction), positions=position)
        cache[dt] = entry
        return entry

    hi = probe(full_dt)
    if hi is None:
        return finish(stop_reason)
    if hi['blocked']:
        hi = None
        return finish('FULL_ENDPOINT_IS_NOT_VISIBLE')

    # Use a freshly solved visible endpoint's router fraction as a guess only.
    fraction = hi['event_fraction']
    center = full_dt*(float(initial_fraction) if fraction is None else fraction)
    half_width = .45*time_tolerance
    lower_dt = max(minimum_dt, min(full_dt, center-half_width))
    lower = probe(lower_dt)
    if lower is None:
        return finish(stop_reason)
    if lower['blocked']:
        lo = lower
        # Seek the real visible side near the guessed event before bisecting.
        offset = half_width
        while True:
            upper_dt = min(hi['dt_s'], max(minimum_dt, center+offset))
            if upper_dt <= lo['dt_s']:
                upper_dt = hi['dt_s']
            upper = probe(upper_dt)
            if upper is None:
                return finish(stop_reason)
            if not upper['blocked']:
                hi = upper
                break
            lo = upper
            offset *= 2.
            if not math.isfinite(offset):
                return finish('NONFINITE_SEARCH_TIME')
    else:
        hi = lower
        # The first lower guess was already visible. Expand toward the known
        # blocked parent, while letting each NEW visible solve guide the guess.
        offset = half_width
        while lo is None:
            hint = None if hi['event_fraction'] is None else hi['dt_s']*hi['event_fraction']
            if hint is not None and minimum_dt < hint < hi['dt_s']-half_width:
                center, offset = hint, half_width
            else:
                offset *= 2.
            lower_dt = max(minimum_dt, min(hi['dt_s'], center-offset))
            if lower_dt >= hi['dt_s']:
                lower_dt = max(minimum_dt, .5*hi['dt_s'])
            if lower_dt == hi['dt_s']:
                return finish('NO_POSITIVE_BLOCKED_ENDPOINT_ABOVE_MINIMUM_DT')
            lower = probe(lower_dt)
            if lower is None:
                return finish(stop_reason)
            if lower['blocked']:
                lo = lower
            else:
                hi = lower
                if lower_dt == minimum_dt:
                    return finish('NO_POSITIVE_BLOCKED_ENDPOINT_ABOVE_MINIMUM_DT')

        # A later visible solve can correct a poor initial guess by a large
        # amount.  Its newly proposed upper side still needs a REAL solve;
        # bisecting against the much older, distant hi would waste the budget.
        upper_dt = min(hi['dt_s'], max(minimum_dt, center+offset))
        if lo['dt_s'] < upper_dt < hi['dt_s']:
            upper = probe(upper_dt)
            if upper is None:
                return finish(stop_reason)
            if upper['blocked']:
                lo = upper
            else:
                hi = upper

    while True:
        current = finish()
        if current['status'] == BRACKETED:
            return current
        if current['maximum_position_difference_m'] is None:
            return finish('NONFINITE_ENDPOINT_POSITION_DIFFERENCE')
        midpoint = lo['dt_s']+.5*(hi['dt_s']-lo['dt_s'])
        if not lo['dt_s'] < midpoint < hi['dt_s'] or midpoint < minimum_dt:
            return finish('TIME_BRACKET_CANNOT_BE_REFINED')
        middle = probe(midpoint)
        if middle is None:
            return finish(stop_reason)
        if middle['blocked']:
            lo = middle
        else:
            hi = middle

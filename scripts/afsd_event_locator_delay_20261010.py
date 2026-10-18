"""Delay-qualified incoming event location using actual same-parent solves.

A nonnegative original strict intersection length can guide a proposal from
two blocked endpoints. The measured decrease is NOT an acceptance condition.
Every returned endpoint is physically solved; original v1 evidence is intact.

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
FROZEN_PARENT_SHA256 = '07341bc16976e59e2e7a516a44371f62530bc6fe11eb56a7e5d651599970a56c'


def _finite_real(value):
    return (not isinstance(value, (bool, np.bool_))
        and isinstance(value, numbers.Real) and math.isfinite(float(value)))


def locate_activation(solve_probe, full_dt, initial_fraction, *, time_tolerance,
        maximum_solves=10, minimum_dt=1e-10,
        maximum_position_difference=MAXIMUM_POSITION_DIFFERENCE_M):
    """Return real blocked/visible endpoints or a bounded, explicit stop.

    ``solve_probe(dt)`` returns ``blocked`` (native bool), ``event_fraction``
    (None or a finite fraction in [0,1]), ``positions`` (finite N by 3), and an
    opaque ``payload``. Required ``activation_delay_s`` is None for blocked
    endpoints; for visible endpoints it is the finite nonnegative conservative
    delay from that candidate's own router start fraction. The callback owns
    that geometric evidence; delay is never used to interpolate or accept a
    state. The returned high must satisfy its own delay <= time_tolerance.
    Required ``activation_prefix_certified`` is a native bool, False at blocked
    endpoints. A visible high can finish only when it is True: the callback
    must certify that candidate's original blocked prefix under the same
    registered budget. False can only request further real bracket refinement;
    it cannot classify an event or substitute a state. Geometry/domain failures
    propagate from the callback and must not be disguised as budget exhaustion.
    Optional event_count, event_kind and outgoing fields
    are checked when present.  Every call must independently restart from the
    SAME parent and classify the SAME incoming event.  Callback exceptions
    propagate unchanged; the callback also owns the wall-clock deadline.

    Payload objects retain their identity and are never copied or inspected.
    Optional ``blocked_measure`` is the original filter strict positive-length
    intersection measure in metres (visible requires zero; a roundoff-only
    blocked zero cannot guide a secant). Negative/nonfinite values are invalid.
    A trustworthy decreasing secant predicts times strictly inside the real
    bracket; it never predicts or interpolates a physical state.

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

    attempts, cache, proposals = [], {}, []
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
            and delta <= maximum_position_difference
            and hi['activation_delay_s'] <= time_tolerance
            and hi['activation_prefix_certified'] is True)
        return dict(status=BRACKETED if success else STOP,
            reason='REAL_ENDPOINT_BRACKET_WITHIN_TIME_POSITION_DELAY_AND_PREFIX_BOUNDS' if success
                else reason or 'BRACKET_NOT_WITHIN_BOUNDS',
            lo=lo, hi=hi, attempts=attempts, solve_count=len(attempts),
            proposal_history=proposals,
            location_method="REAL_ENDPOINTS_WITH_DECREASING_BLOCKED_MEASURE_SECANT_AND_HIGH_PATH_DELAY_PREFIX",
            frozen_parent_sha256=FROZEN_PARENT_SHA256,
            high_activation_delay_s=None if hi is None else hi['activation_delay_s'],
            activation_delay_limit_s=time_tolerance,
            activation_delay_within_bound=False if hi is None else hi['activation_delay_s']<=time_tolerance,
            high_activation_prefix_certified=False if hi is None else hi['activation_prefix_certified'],
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
        if ('event_fraction' not in raw or 'positions' not in raw or 'payload' not in raw
                or 'activation_delay_s' not in raw or 'activation_prefix_certified' not in raw):
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
        measure = raw.get('blocked_measure')
        if 'blocked_measure' in raw and (not _finite_real(measure) or float(measure) < 0.):
            stop_reason = 'INVALID_BLOCKED_INTERSECTION_MEASURE'
            return None
        if measure is not None and not raw['blocked'] and float(measure) != 0.:
            stop_reason = 'VISIBLE_ENDPOINT_HAS_POSITIVE_BLOCKED_MEASURE'
            return None
        delay = raw['activation_delay_s']
        if raw['blocked'] and delay is not None:
            stop_reason = 'BLOCKED_ENDPOINT_HAS_ACTIVATION_DELAY'
            return None
        if not raw['blocked'] and (not _finite_real(delay) or float(delay) < 0.):
            stop_reason = 'INVALID_ACTIVATION_DELAY'
            return None
        prefix = raw['activation_prefix_certified']
        if type(prefix) is not bool:
            stop_reason = 'ACTIVATION_PREFIX_CERTIFIED_IS_NOT_A_NATIVE_BOOL'
            return None
        if raw['blocked'] and prefix:
            stop_reason = 'BLOCKED_ENDPOINT_HAS_CERTIFIED_ACTIVATION_PREFIX'
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
            event_fraction=None if fraction is None else float(fraction), positions=position,
            activation_delay_s=None if delay is None else float(delay),
            activation_prefix_certified=prefix)
        if measure is not None:
            entry['blocked_measure'] = float(measure)
        cache[dt] = entry
        return entry

    def solved(dt, proposal):
        nonlocal lo, hi
        value = probe(dt)
        if value is None:
            return None
        proposals.append(dict(dt_s=float(dt), **proposal))
        if value['blocked']:
            if lo is None or value['dt_s'] > lo['dt_s']:
                lo = value
        elif hi is None or value['dt_s'] < hi['dt_s']:
            hi = value
        return value

    full = solved(full_dt, dict(kind='FULL_REAL_ENDPOINT'))
    if full is None:
        return finish(stop_reason)
    if full['blocked']:
        return finish('FULL_ENDPOINT_IS_NOT_VISIBLE')

    fraction = hi['event_fraction']
    center = full_dt*(float(initial_fraction) if fraction is None else fraction)
    half_width = .45*time_tolerance
    expansion_stride = .9*time_tolerance

    def real_pair_at_guess(center, kind, extra=None):
        """Try real lower/upper endpoints around a proposal, respecting hi/lo."""
        extra = {} if extra is None else extra
        lower_bound = 0. if lo is None else lo['dt_s']
        lower_dt = max(minimum_dt, center-half_width)
        if not lower_bound < lower_dt < hi['dt_s']:
            return False
        value = solved(lower_dt, dict(kind=kind+'_LOWER', center_s=center, **extra))
        if value is None:
            return True
        if finish()['status'] == BRACKETED:
            return True
        if value['blocked']:
            upper_dt = center+half_width
            if lo['dt_s'] < upper_dt < hi['dt_s']:
                solved(upper_dt, dict(kind=kind+'_UPPER', center_s=center, **extra))
        return True

    # The original first narrow pair is retained. Two genuinely blocked
    # measurements can already provide a much better event-time proposal.
    if not real_pair_at_guess(center, 'INITIAL_ROUTER'):
        first_dt = max(minimum_dt, .5*hi['dt_s'])
        if first_dt >= hi['dt_s']:
            return finish('NO_POSITIVE_BLOCKED_ENDPOINT_ABOVE_MINIMUM_DT')
        solved(first_dt, dict(kind='INITIAL_BOUNDED_MIDPOINT'))
    if stop_reason:
        return finish(stop_reason)

    while True:
        current = finish()
        if current['status'] == BRACKETED:
            return current
        if len(attempts) >= maximum_solves:
            return finish('MAXIMUM_PHYSICAL_SOLVES')
        if lo is not None and current['maximum_position_difference_m'] is None:
            return finish('NONFINITE_ENDPOINT_POSITION_DIFFERENCE')
        if hi['dt_s'] <= minimum_dt and lo is None:
            return finish('NO_POSITIVE_BLOCKED_ENDPOINT_ABOVE_MINIMUM_DT')

        lower_bound = 0. if lo is None else lo['dt_s']
        width = hi['dt_s']-lower_bound
        proposed = False
        # Inside a narrow time bracket, ordinary bisection can still be needed
        # to meet the independent original 50 nm position difference bound
        # or the visible endpoint own-path delay and blocked-prefix bounds.
        if width > time_tolerance:
            blocked = sorted((q for q in cache.values()
                if q['blocked'] and q.get('blocked_measure', 0.) > 0.),
                key=lambda q:q['dt_s'])
            if len(blocked) >= 2:
                first, second = blocked[-2:]
                change = first['blocked_measure']-second['blocked_measure']
                interval = second['dt_s']-first['dt_s']
                arithmetic = 64*np.finfo(float).eps*max(first['blocked_measure'],second['blocked_measure'])
                if interval > 0. and change > arithmetic:
                    decrease = change/interval
                    if math.isfinite(decrease) and decrease > 0.:
                        estimate = second['dt_s']+second['blocked_measure']/decrease
                        if math.isfinite(estimate) and lower_bound < estimate < hi['dt_s']:
                            proposed = real_pair_at_guess(estimate, 'BLOCKED_MEASURE_SECANT',
                                dict(reference_dt_s=[first['dt_s'],second['dt_s']],
                                     reference_measure_m=[first['blocked_measure'],second['blocked_measure']],
                                     positive_measure_decrease_m_s=decrease))
            if not proposed and hi['event_fraction'] is not None:
                estimate = hi['dt_s']*hi['event_fraction']
                if lower_bound < estimate < hi['dt_s']:
                    proposed = real_pair_at_guess(estimate, 'LATEST_REAL_VISIBLE_ROUTER')
        if stop_reason:
            return finish(stop_reason)
        if proposed:
            continue
        # Without trustworthy geometry slope or an interior router hint,
        # first expand the latest real endpoint locally. This preserves the
        # useful narrow-guess behavior when the optional measure is absent.
        # Expansion is only a proposed real solve and stays inside hi/lo.
        if width > time_tolerance:
            expanded = (hi['dt_s']-expansion_stride if lo is None
                        else lo['dt_s']+expansion_stride)
            expanded = max(minimum_dt, expanded)
            expansion_stride *= 2.
            if lower_bound < expanded < hi['dt_s']:
                solved(expanded, dict(kind='BOUNDED_REAL_ENDPOINT_EXPANSION'))
                if stop_reason:
                    return finish(stop_reason)
                continue
        middle_dt = lower_bound+.5*width
        middle_dt = max(minimum_dt,middle_dt)
        if not lower_bound < middle_dt < hi['dt_s']:
            return finish('TIME_BRACKET_CANNOT_BE_REFINED')
        solved(middle_dt, dict(kind='BOUNDED_REAL_BRACKET_BISECTION'))
        if stop_reason:
            return finish(stop_reason)

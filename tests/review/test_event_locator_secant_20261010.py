"""Secant locator: independent physics, negative probes and measured-data diagnostic fit.

The saved-v1 data oracle is a fit outside the observed blocked sample range;
it is an algorithm diagnostic, never evidence that new physical solves passed.
"""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from afsd_event_locator_secant_20261010 import locate_activation, BRACKETED, STOP


class OpaqueState:
    def __deepcopy__(self, memo):
        raise AssertionError('physical payload must not be deep-copied')


class ForcedParticle:
    """Exact m*x''=F from one fixed parent; activation at x=x(event_time)."""
    def __init__(self, event_time=.6e-3, fraction_bias_s=0., initial_hint=None,
            velocity=.002, acceleration=.5):
        self.parent = dict(x=np.array([[.001, 0., 0.]]), v=velocity, a=acceleration, mass=2.)
        self.original = copy.deepcopy(self.parent)
        self.event_time = event_time
        self.fraction_bias_s = fraction_bias_s
        self.initial_hint = initial_hint
        self.calls = []
        self.payloads = []

    def __call__(self, dt):
        np.testing.assert_array_equal(self.parent['x'], self.original['x'])
        assert self.parent['v'] == self.original['v']
        v, a = self.parent['v'], self.parent['a']
        displacement = v*dt+.5*a*dt*dt
        positions = self.parent['x']+np.array([[displacement, 0., 0.]])
        payload = OpaqueState()
        payload.dt = dt
        payload.position = positions.copy()
        payload.velocity = v+a*dt
        payload.work = self.parent['mass']*a*displacement
        payload.kinetic_change = .5*self.parent['mass']*(payload.velocity**2-v**2)
        self.calls.append(dt)
        self.payloads.append(payload)
        # This manufactured probe has an independently known physical event.
        # Deliberately biased fractions model a diagnostic locator's inexact
        # proposed event time, without perturbing any physical endpoint.
        blocked = bool(dt < self.event_time)
        proposed = self.event_time+self.fraction_bias_s
        if len(self.calls) == 1 and self.initial_hint is not None:
            proposed = self.initial_hint
        fraction = None if blocked else min(1., max(0., proposed/dt))
        return dict(blocked=blocked, event_fraction=fraction, positions=positions,
            payload=payload, event_count=0 if blocked else 1,
            event_kind=None if blocked else 'BLOCKED_TO_VISIBLE', outgoing=False)


def locate(probe, **kwargs):
    return locate_activation(probe, .001, .6, time_tolerance=1e-7, **kwargs)


class PhysicalBracketTests(unittest.TestCase):
    def assert_physical_bracket(self, result, model):
        self.assertEqual(result['status'], BRACKETED, result['reason'])
        lo, hi = result['lo'], result['hi']
        self.assertGreater(lo['dt_s'], 0.)
        self.assertLess(lo['dt_s'], model.event_time)
        self.assertGreaterEqual(hi['dt_s'], model.event_time)
        self.assertLessEqual(hi['dt_s']-lo['dt_s'], result['time_tolerance_s'])
        self.assertLessEqual(result['maximum_position_difference_m'], 5e-8)
        self.assertFalse(result['accepted'])
        self.assertFalse(result['restart_allowed'])
        self.assertEqual(result['interpolated_states'], 0)
        for bound in (lo, hi):
            self.assertTrue(any(bound['payload'] is item for item in model.payloads))
            np.testing.assert_array_equal(bound['positions'], bound['payload'].position)
            roundoff = 32*np.finfo(float).eps*model.parent['mass']*(
                model.parent['v']**2+bound['payload'].velocity**2)
            self.assertAlmostEqual(bound['payload'].work, bound['payload'].kinetic_change,
                delta=roundoff)
        self.assertEqual(result['solve_count'], len(model.calls))
        self.assertLessEqual(result['solve_count'], 10)

    def test_exact_one_sided_solution_is_bracketed_with_two_real_endpoints(self):
        model = ForcedParticle()
        result = locate(model)
        self.assert_physical_bracket(result, model)
        self.assertEqual(result['solve_count'], 3)

    def test_bad_initial_router_fraction_is_only_a_guess(self):
        model = ForcedParticle(initial_hint=.0008)
        result = locate(model)
        self.assert_physical_bracket(result, model)
        self.assertGreater(result['solve_count'], 3)

    def test_consistently_early_diagnostic_guess_requires_real_expansion(self):
        model = ForcedParticle(fraction_bias_s=-1.1e-7)
        result = locate(model)
        self.assert_physical_bracket(result, model)
        self.assertGreater(result['solve_count'], 3)

    def test_time_tolerance_alone_cannot_accept_large_position_difference(self):
        model = ForcedParticle(velocity=4., acceleration=0.)
        result = locate(model)
        self.assert_physical_bracket(result, model)
        self.assertLess(result['bracket_width_s'], 1e-7)
        self.assertGreater(result['solve_count'], 3)

    def test_stricter_tolerances_converge_and_preserve_work(self):
        widths = []
        for tolerance in (1e-6, 1e-7, 1e-8):
            model = ForcedParticle()
            result = locate_activation(model, .001, .6, time_tolerance=tolerance)
            self.assert_physical_bracket(result, model)
            widths.append(result['bracket_width_s'])
        self.assertTrue(all(a > b for a, b in zip(widths, widths[1:])))

    def test_parent_is_not_modified_and_payload_is_not_copied(self):
        model = ForcedParticle()
        result = locate(model)
        self.assert_physical_bracket(result, model)
        np.testing.assert_array_equal(model.parent['x'], model.original['x'])
        self.assertFalse(result['lo']['positions'].flags.writeable)

    def test_snapshots_do_not_alias_callback_buffers(self):
        model = ForcedParticle()
        buffer = np.zeros((1, 3))
        def reused(dt):
            value = model(dt)
            buffer[:] = value['positions']
            value['positions'] = buffer
            return value
        result = locate(reused)
        before = result['lo']['positions'].copy()
        buffer[:] = 99.
        np.testing.assert_array_equal(result['lo']['positions'], before)
        self.assert_physical_bracket(result, model)

    def test_maximum_solves_cannot_be_bypassed_by_exact_diagnostic_fraction(self):
        model = ForcedParticle()
        result = locate(model, maximum_solves=2)
        self.assertEqual(result['status'], STOP)
        self.assertEqual(result['reason'], 'MAXIMUM_PHYSICAL_SOLVES')
        self.assertEqual(len(model.calls), 2)

    def test_no_positive_blocked_endpoint_cannot_use_assumed_zero(self):
        model = ForcedParticle(event_time=5e-11)
        result = locate(model)
        self.assertEqual(result['status'], STOP)
        self.assertIsNone(result['lo'])
        self.assertEqual(result['reason'], 'NO_POSITIVE_BLOCKED_ENDPOINT_ABOVE_MINIMUM_DT')
        self.assertTrue(all(dt >= 1e-10 for dt in model.calls))

    def test_full_endpoint_must_actually_be_visible(self):
        result = locate(ForcedParticle(event_time=.002))
        self.assertEqual(result['status'], STOP)
        self.assertEqual(result['reason'], 'FULL_ENDPOINT_IS_NOT_VISIBLE')

    def test_continuous_physical_family_never_returns_a_false_bracket(self):
        successes = 0
        for event_time in (1.2e-4, 5.5e-4, 9.4e-4):
            for bias in (-1.3e-7, 0., 1.3e-7):
                for velocity in (.002, .08, 2.):
                    model = ForcedParticle(event_time=event_time,
                        fraction_bias_s=bias, velocity=velocity)
                    result = locate(model)
                    self.assertLessEqual(len(model.calls), 10)
                    self.assertTrue(all(1e-10 <= dt <= .001 for dt in model.calls))
                    if result['status'] == BRACKETED:
                        successes += 1
                        self.assert_physical_bracket(result, model)
                    else:
                        self.assertEqual(result['reason'], 'MAXIMUM_PHYSICAL_SOLVES')
        self.assertGreaterEqual(successes, 20)

    def test_position_jump_cannot_be_hidden_by_an_arbitrarily_small_time_bracket(self):
        calls = []
        def discontinuous(dt):
            calls.append(dt)
            blocked = dt < .0006
            return dict(blocked=blocked, event_fraction=None if blocked else .0006/dt,
                positions=np.array([[0. if blocked else 1e-6, 0., 0.]]), payload=object())
        result = locate(discontinuous)
        self.assertEqual(result['status'], STOP)
        self.assertEqual(result['reason'], 'MAXIMUM_PHYSICAL_SOLVES')
        self.assertEqual(len(calls), 10)
        self.assertGreater(result['maximum_position_difference_m'], 5e-8)


class RejectionTests(unittest.TestCase):
    def test_invalid_probe_directions_counts_fields_and_coordinates_stop(self):
        cases = [
            ({'blocked': 0}, 'ENDPOINT_BLOCKED_IS_NOT_A_NATIVE_BOOL'),
            ({'blocked': np.bool_(False)}, 'ENDPOINT_BLOCKED_IS_NOT_A_NATIVE_BOOL'),
            ({'event_count': 2}, 'MULTIPLE_OR_INVALID_EVENT_COUNT'),
            ({'event_count': 0, 'event_fraction': None}, 'EVENT_COUNT_AND_ENDPOINT_DISAGREE'),
            ({'outgoing': True}, 'OUTGOING_OR_INVALID_EVENT_DIRECTION'),
            ({'event_kind': 'VISIBLE_TO_BLOCKED'}, 'OUTGOING_OR_UNSUPPORTED_EVENT_KIND'),
            ({'event_fraction': float('nan')}, 'EVENT_FRACTION_OUTSIDE_FINITE_UNIT_INTERVAL'),
            ({'event_fraction': 1.0001}, 'EVENT_FRACTION_OUTSIDE_FINITE_UNIT_INTERVAL'),
            ({'event_fraction': -1e-4}, 'EVENT_FRACTION_OUTSIDE_FINITE_UNIT_INTERVAL'),
            ({'positions': [[np.inf, 0., 0.]]}, 'INVALID_OR_INCONSISTENT_POSITIONS'),
            ({'positions': [0., 0., 0.]}, 'INVALID_OR_INCONSISTENT_POSITIONS'),
            ({'positions': np.zeros((1, 3), complex)}, 'INVALID_OR_INCONSISTENT_POSITIONS'),
            ({'dt_s': .002}, 'PROBE_TIME_DOES_NOT_MATCH_REQUEST'),
        ]
        for changes, reason in cases:
            with self.subTest(changes=changes):
                model = ForcedParticle()
                def invalid(dt):
                    result = model(dt)
                    result.update(changes)
                    return result
                result = locate(invalid)
                self.assertEqual(result['status'], STOP)
                self.assertEqual(result['reason'], reason)
                self.assertEqual(result['solve_count'], 1)

    def test_blocked_endpoint_cannot_also_report_an_activation(self):
        model = ForcedParticle()
        def inconsistent(dt):
            result = model(dt)
            if result['blocked']:
                result['event_fraction'] = .5
            return result
        result = locate(inconsistent)
        self.assertEqual(result['reason'], 'ACTIVATION_REPORTED_AT_BLOCKED_ENDPOINT')

    def test_changed_position_shape_is_rejected(self):
        model = ForcedParticle()
        def inconsistent(dt):
            result = model(dt)
            if len(model.calls) > 1:
                result['positions'] = np.zeros((2, 3))
            return result
        self.assertEqual(locate(inconsistent)['reason'], 'INVALID_OR_INCONSISTENT_POSITIONS')

    def test_callback_physical_failure_and_deadline_propagate(self):
        failure = RuntimeError('real physical candidate failed')
        def failed(dt):
            raise failure
        with self.assertRaises(RuntimeError) as raised:
            locate(failed)
        self.assertIs(raised.exception, failure)

    def test_invalid_configuration_is_rejected_before_any_solve(self):
        for kwargs in [dict(maximum_solves=11), dict(maximum_solves=True),
                dict(minimum_dt=.001), dict(maximum_position_difference=5.1e-8)]:
            with self.subTest(kwargs=kwargs):
                model = ForcedParticle()
                with self.assertRaises(ValueError):
                    locate(model, **kwargs)
                self.assertEqual(model.calls, [])


class AcceleratingVisibility(ForcedParticle):
    """Exact constant-force motion; geometry uses its own linear endpoint chord."""
    def __call__(self, dt):
        result = super().__call__(dt)
        v, a = self.parent['v'], self.parent['a']
        event_position = v*self.event_time+.5*a*self.event_time**2
        actual_position = v*dt+.5*a*dt**2
        result['blocked_measure'] = max(0., event_position-actual_position)
        result['event_fraction'] = None if result['blocked'] else event_position/actual_position
        return result


class SecantPhysicsTests(unittest.TestCase):
    assert_physical_bracket = PhysicalBracketTests.assert_physical_bracket

    def test_nonlinear_physical_motion_uses_geometry_only_as_a_guess(self):
        for event_time, velocity, acceleration in [(.000249, .00055, 2.11),
                (.000020, .002, 100.), (.000340, .0004, 5.)]:
            with self.subTest(event_time=event_time):
                model = AcceleratingVisibility(event_time=event_time,
                    velocity=velocity, acceleration=acceleration)
                result = locate_activation(model, .000474609375, .35550876846536994,
                    time_tolerance=1e-7)
                self.assert_physical_bracket(result, model)
                self.assertTrue(any(p['kind'].startswith('BLOCKED_MEASURE_SECANT')
                    for p in result['proposal_history']))
                self.assertTrue(all(0. < t <= .000474609375 for t in model.calls))

    def test_noisy_increasing_measure_is_not_trusted_as_a_decreasing_secant(self):
        model = ForcedParticle(initial_hint=.0005998)
        def wrong_slope(dt):
            result = model(dt)
            result['blocked_measure'] = dt if result['blocked'] else 0.
            return result
        result = locate(wrong_slope)
        self.assert_physical_bracket(result, model)
        self.assertFalse(any(p['kind'].startswith('BLOCKED_MEASURE_SECANT')
            for p in result['proposal_history']))

    def test_outside_bracket_secant_does_not_extrapolate_a_physical_solve(self):
        model = ForcedParticle(initial_hint=.0005998)
        def distant_root(dt):
            result = model(dt)
            result['blocked_measure'] = 1.-1e-6*dt if result['blocked'] else 0.
            return result
        result = locate(distant_root)
        self.assert_physical_bracket(result, model)
        self.assertTrue(all(1e-10 <= t <= .001 for t in model.calls))
        self.assertFalse(any(p['kind'].startswith('BLOCKED_MEASURE_SECANT')
            for p in result['proposal_history']))

    def test_roundoff_only_blocked_zero_measure_uses_bounded_fallback(self):
        model = ForcedParticle(initial_hint=.0005998)
        def zero_measure(dt):
            result = model(dt)
            result['blocked_measure'] = 0.
            return result
        result = locate(zero_measure)
        self.assert_physical_bracket(result, model)
        self.assertFalse(any(p['kind'].startswith('BLOCKED_MEASURE_SECANT')
            for p in result['proposal_history']))


class MeasuredSequenceDiagnosticTests(unittest.TestCase):
    """Use actual v1 blocked measurements, with explicitly unvalidated extrapolation.

    Source: event_validation_20261010_v1/attempt_0002/coarse/event_probe_02..10.
    Values were independently recalculated with the frozen original endpoint
    filter on PID 1334--1371 from each saved graph/motion/checkpoint.  The full
    visible probe's own router fraction is also measured.  A quadratic fit
    outside these observed blocked times is ONLY a diagnostic oracle; success
    here neither accepts an old probe nor establishes a real new event time.
    """
    times = np.array([168.68279440836892, 168.77279440836892, 168.81779440836893,
        168.9077944083689, 169.08779440836892, 169.44779440836893,
        170.16779440836893, 171.60779440836892, 174.48779440836892])*1e-6
    measures = np.array([8.007990104919844e-8, 7.999833689432452e-8,
        7.995754840748444e-8, 7.987595861556766e-8, 7.971267349107168e-8,
        7.938600365268643e-8, 7.873187158099601e-8,
        7.742029584364029e-8, 7.478392821376174e-8])

    def test_measured_initial_failure_sequence_guides_a_bounded_diagnostic_fit(self):
        x, y = self.times, self.measures
        coefficients = np.polyfit((x-x[0])/1e-6, y, 2)
        fit_event = x[0]+max(np.roots(coefficients))*1e-6
        initial_clearance = np.polyval(coefficients, -x[0]/1e-6)
        actual_calls = []
        def diagnostic_fit(dt):
            signed = float(np.polyval(coefficients, (dt-x[0])/1e-6))
            # Retain the original measurement at every stored sample time.
            # Intermediate/outside values are labelled interpolation/fit only.
            exact = np.flatnonzero(np.abs(x-dt) < 1e-16)
            if len(exact):
                signed = float(y[int(exact[0])])
            blocked = bool(signed > 0.)
            fraction = None if blocked else float(initial_clearance/(initial_clearance-signed))
            if not actual_calls:
                self.assertEqual(dt, .000474609375)
                blocked, fraction = False, .35550876846536994
            actual_calls.append(dt)
            return dict(blocked=blocked, event_fraction=fraction,
                blocked_measure=max(0., signed) if blocked else 0.,
                positions=np.array([[.001+initial_clearance-signed, 0., 0.]]),
                payload=dict(classification='DIAGNOSTIC_FIT_NOT_PHYSICAL_SOLUTION', dt_s=dt))
        result = locate_activation(diagnostic_fit, .000474609375, .35550876846536994,
            time_tolerance=1e-7)
        self.assertEqual(result['status'], BRACKETED, result['reason'])
        self.assertLess(result['lo']['dt_s'], fit_event)
        self.assertGreaterEqual(result['hi']['dt_s'], fit_event)
        self.assertLessEqual(result['solve_count'], 10)
        self.assertLessEqual(result['bracket_width_s'], 1e-7)
        self.assertLessEqual(result['maximum_position_difference_m'], 5e-8)
        self.assertFalse(result['accepted'])
        np.testing.assert_allclose(actual_calls[1:3], x[:2], rtol=0., atol=1e-16)
        self.assertTrue(any(p['kind'].startswith('BLOCKED_MEASURE_SECANT')
            for p in result['proposal_history']))


class SecantRejectionTests(unittest.TestCase):
    def test_invalid_and_contradictory_measure_stops_without_state_promotion(self):
        for value, reason in [(float('nan'), 'INVALID_BLOCKED_INTERSECTION_MEASURE'),
                (float('inf'), 'INVALID_BLOCKED_INTERSECTION_MEASURE'),
                (-1e-10, 'INVALID_BLOCKED_INTERSECTION_MEASURE'),
                (True, 'INVALID_BLOCKED_INTERSECTION_MEASURE'),
                (1e-9, 'VISIBLE_ENDPOINT_HAS_POSITIVE_BLOCKED_MEASURE')]:
            with self.subTest(measure=value):
                model = ForcedParticle()
                def invalid(dt):
                    result = model(dt)
                    result['blocked_measure'] = value
                    return result
                result = locate(invalid)
                self.assertEqual(result['status'], STOP)
                self.assertEqual(result['reason'], reason)
                self.assertEqual(len(model.calls), 1)


if __name__ == '__main__':
    unittest.main()

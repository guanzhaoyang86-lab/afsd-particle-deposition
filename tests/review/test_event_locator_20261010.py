"""Independent constant-force particle probes for bounded event location."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from afsd_event_locator_20261010 import locate_activation, BRACKETED, STOP


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


if __name__ == '__main__':
    unittest.main()

"""Static legacy fixture adaptation plus independent high-path delay checks.

The saved-v1 data oracle is a fit outside the observed blocked sample range;
it is an algorithm diagnostic, never evidence that new physical solves passed.
"""
from pathlib import Path
import copy
import ast
import hashlib
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from afsd_event_locator_delay_20261010 import locate_activation, BRACKETED, STOP
# Static fixture adaptation: only the required callback delay evidence is added.
# Original physical equations, opaque-state assertions and rejection tests remain.


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
            activation_delay_s=None if blocked else max(0.,dt-self.event_time),
            activation_prefix_certified=not blocked,
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
        self.assertLessEqual(result['hi']['activation_delay_s'],result['time_tolerance_s'])
        self.assertTrue(result['activation_delay_within_bound'])
        self.assertTrue(result['high_activation_prefix_certified'])
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
                activation_delay_s=None if blocked else 0.,activation_prefix_certified=not blocked,
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
        result['activation_delay_s'] = None if result['blocked'] else dt*(1.-result['event_fraction'])
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
                activation_delay_s=None if blocked else dt*(1.-fraction),activation_prefix_certified=not blocked,
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
        self.assertLessEqual(result['hi']['activation_delay_s'],result['time_tolerance_s'])
        self.assertTrue(result['activation_delay_within_bound'])
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


class SavedFailureShape:
    """Synthetic callbacks following v4 half-tolerance's measured proposals.

    The first eight times/measures/fractions describe the saved failure shape;
    this oracle then supplies a ninth *synthetic* independently computed state.
    It does not claim a GPU solve or promote an old diagnostic state.
    """
    rows=[
        (.000474609375,False,.35550876846536994,0.,.00030588158081263825),
        (.00016870529440836892,True,None,8.005952235947122e-8,None),
        (.00016875029440836892,True,None,8.001872202461743e-8,None),
        (.0002569830134895916,False,.9545198115520179,0.,1.1687636001103987e-5),
        (.0002452728776081547,True,None,4.4299188073777186e-9,None),
        (.00024531787760815466,True,None,4.381764206028607e-9,None),
        (.00024939009304025963,True,None,6.6813479740122535e-12,None),
        (.0002494350930402596,False,.9997591120190918,0.,6.008603208237255e-8),
    ]

    def __init__(self,prefix_ready_call=9,high_delay_override=None):
        self.calls=[];self.payloads=[]
        self.prefix_ready_call=prefix_ready_call;self.high_delay_override=high_delay_override

    def __call__(self,dt):
        if len(self.calls)<8:
            expected,blocked,fraction,measure,delay=self.rows[len(self.calls)]
            if abs(expected-dt)>1e-18:
                raise AssertionError('Original secant proposal sequence changed')
        else:
            old_lo,old_hi=self.rows[-2][0],self.rows[-1][0]
            event=old_lo+.1*(old_hi-old_lo)
            blocked=bool(dt<event)
            delay=None if blocked else max(0.,self.rows[-1][-1]+dt-old_hi)
            fraction=None if blocked else 1.-delay/dt
            measure=max(0.,event-dt)*.001 if blocked else 0.
        if len(self.calls)==7 and self.high_delay_override is not None:
            delay=self.high_delay_override
        positions=np.array([[.001+.04455*dt,0.,0.]])
        payload=OpaqueState();payload.dt=dt;payload.position=positions.copy()
        self.calls.append(dt);self.payloads.append(payload)
        return dict(blocked=blocked,event_fraction=fraction,blocked_measure=measure,
            positions=positions,payload=payload,activation_delay_s=delay,
            activation_prefix_certified=not blocked and len(self.calls)>=self.prefix_ready_call,
            event_count=0 if blocked else 1,event_kind=None if blocked else 'BLOCKED_TO_VISIBLE')


class OwnPathDelayTests(unittest.TestCase):
    def test_delay_pass_but_prefix_false_requires_ninth_solved_endpoint(self):
        capped=locate_activation(SavedFailureShape(high_delay_override=4e-8),
            .000474609375,.35550876846536994,time_tolerance=5e-8,maximum_solves=8)
        self.assertEqual(capped['status'],STOP)
        self.assertTrue(capped['activation_delay_within_bound'])
        self.assertFalse(capped['high_activation_prefix_certified'])
        self.assertLessEqual(capped['bracket_width_s'],5e-8)
        result=locate_activation(SavedFailureShape(high_delay_override=4e-8),
            .000474609375,.35550876846536994,time_tolerance=5e-8)
        self.assertEqual(result['status'],BRACKETED)
        self.assertEqual(result['solve_count'],9)
        self.assertTrue(result['high_activation_prefix_certified'])

    def test_prefix_false_at_ninth_probe_uses_tenth_or_stops_at_nine_cap(self):
        capped=locate_activation(SavedFailureShape(prefix_ready_call=10),
            .000474609375,.35550876846536994,time_tolerance=5e-8,maximum_solves=9)
        self.assertEqual(capped['status'],STOP)
        self.assertEqual(capped['reason'],'MAXIMUM_PHYSICAL_SOLVES')
        self.assertTrue(capped['activation_delay_within_bound'])
        self.assertFalse(capped['high_activation_prefix_certified'])
        model=SavedFailureShape(prefix_ready_call=10)
        result=locate_activation(model,.000474609375,.35550876846536994,time_tolerance=5e-8)
        self.assertEqual(result['status'],BRACKETED)
        self.assertEqual(result['solve_count'],10)
        self.assertIs(result['hi']['payload'],model.payloads[-1])
        self.assertTrue(result['hi']['activation_prefix_certified'])
        self.assertEqual(result['interpolated_states'],0)
        self.assertEqual(result['proposal_history'][-1]['kind'],'BOUNDED_REAL_BRACKET_BISECTION')

    def test_invalid_missing_or_blocked_true_prefix_is_rejected(self):
        model=ForcedParticle()
        def missing(dt):
            row=model(dt);del row['activation_prefix_certified'];return row
        self.assertEqual(locate(missing)['reason'],'PROBE_REQUIRED_FIELD_MISSING')
        for value in (None,0,1,np.bool_(True),'True'):
            model=ForcedParticle()
            def invalid(dt):
                row=model(dt);row['activation_prefix_certified']=value;return row
            self.assertEqual(locate(invalid)['reason'],'ACTIVATION_PREFIX_CERTIFIED_IS_NOT_A_NATIVE_BOOL')
        model=ForcedParticle()
        def blocked_claim(dt):
            row=model(dt);row['activation_prefix_certified']=True;return row
        self.assertEqual(locate(blocked_claim)['reason'],'BLOCKED_ENDPOINT_HAS_CERTIFIED_ACTIVATION_PREFIX')

    def test_prefix_failure_never_changes_endpoint_classification(self):
        model=ForcedParticle()
        def unresolved(dt):
            row=model(dt);row['activation_prefix_certified']=False;return row
        result=locate(unresolved)
        self.assertEqual(result['status'],STOP)
        self.assertEqual(result['solve_count'],10)
        self.assertFalse(result['hi']['blocked'])
        self.assertTrue(result['lo']['blocked'])
        self.assertLess(result['lo']['dt_s'],model.event_time)
        self.assertGreaterEqual(result['hi']['dt_s'],model.event_time)
        self.assertFalse(result['high_activation_prefix_certified'])

    def test_saved_45ns_bracket_60ns_delay_shape_uses_ninth_real_callback(self):
        model=SavedFailureShape()
        result=locate_activation(model,.000474609375,.35550876846536994,time_tolerance=5e-8)
        self.assertEqual(result['status'],BRACKETED,result['reason'])
        self.assertEqual(result['solve_count'],9)
        old_lo,old_hi=model.rows[-2][0],model.rows[-1][0]
        self.assertEqual(model.calls[-1],old_lo+.5*(old_hi-old_lo))
        self.assertEqual(result['proposal_history'][-1]['kind'],'BOUNDED_REAL_BRACKET_BISECTION')
        self.assertLessEqual(result['bracket_width_s'],5e-8)
        self.assertLessEqual(result['maximum_position_difference_m'],5e-8)
        self.assertLessEqual(result['high_activation_delay_s'],5e-8)
        self.assertGreater(model.rows[-1][-1],5e-8)
        self.assertIs(result['hi']['payload'],model.payloads[-1])
        np.testing.assert_array_equal(result['hi']['positions'],model.payloads[-1].position)
        self.assertEqual(result['interpolated_states'],0)
        self.assertFalse(result['accepted'] or result['restart_allowed'])

    def test_same_eight_probe_cap_stops_instead_of_accepting_short_bracket(self):
        model=SavedFailureShape()
        result=locate_activation(model,.000474609375,.35550876846536994,
            time_tolerance=5e-8,maximum_solves=8)
        self.assertEqual(result['status'],STOP)
        self.assertEqual(result['reason'],'MAXIMUM_PHYSICAL_SOLVES')
        self.assertEqual(result['solve_count'],8)
        self.assertAlmostEqual(result['bracket_width_s'],45e-9,delta=1e-18)
        self.assertLess(result['maximum_position_difference_m'],5e-8)
        self.assertGreater(result['high_activation_delay_s'],5e-8)
        self.assertFalse(result['activation_delay_within_bound'])

    def test_invalid_visible_delays_are_not_coerced_or_ignored(self):
        for value in (None,True,np.bool_(False),float('nan'),float('inf'),-1e-20,'0'):
            with self.subTest(value=value):
                model=ForcedParticle()
                def invalid(dt):
                    row=model(dt);row['activation_delay_s']=value;return row
                result=locate(invalid)
                self.assertEqual(result['status'],STOP)
                self.assertEqual(result['reason'],'INVALID_ACTIVATION_DELAY')
                self.assertEqual(result['solve_count'],1)

    def test_missing_required_delay_and_blocked_non_none_delay_stop(self):
        model=ForcedParticle()
        def missing(dt):
            row=model(dt);del row['activation_delay_s'];return row
        self.assertEqual(locate(missing)['reason'],'PROBE_REQUIRED_FIELD_MISSING')
        for bad in (0.,1e-9,False):
            model=ForcedParticle()
            def invalid(dt):
                row=model(dt)
                if row['blocked']:row['activation_delay_s']=bad
                return row
            result=locate(invalid)
            self.assertEqual(result['reason'],'BLOCKED_ENDPOINT_HAS_ACTIVATION_DELAY')
            self.assertEqual(result['solve_count'],2)

    def test_delay_gate_is_independent_of_position_and_time_gates(self):
        # Zero delay does not manufacture a real blocked endpoint or narrow
        # the independently observed time bracket.
        model=ForcedParticle()
        def zero_delay(dt):
            row=model(dt)
            if not row['blocked']:row['activation_delay_s']=0.
            return row
        result=locate(zero_delay,maximum_solves=2)
        self.assertEqual(result['status'],STOP)
        self.assertTrue(result['activation_delay_within_bound'])
        self.assertGreater(result['bracket_width_s'],result['time_tolerance_s'])
        # Nor can a zero delay hide a finite discontinuity in position.
        model=ForcedParticle()
        def jump(dt):
            row=model(dt)
            row['positions']=np.array([[0. if row['blocked'] else 1e-6,0.,0.]])
            if not row['blocked']:row['activation_delay_s']=0.
            return row
        result=locate(jump)
        self.assertEqual(result['status'],STOP)
        self.assertTrue(result['activation_delay_within_bound'])
        self.assertLessEqual(result['bracket_width_s'],result['time_tolerance_s'])
        self.assertGreater(result['maximum_position_difference_m'],5e-8)

    def test_delay_remains_too_large_even_after_refinement_cannot_pass(self):
        model=ForcedParticle()
        def persistent_delay(dt):
            row=model(dt)
            if not row['blocked']:row['activation_delay_s']=2e-7
            return row
        result=locate(persistent_delay)
        self.assertEqual(result['status'],STOP)
        self.assertEqual(result['reason'],'MAXIMUM_PHYSICAL_SOLVES')
        self.assertEqual(result['solve_count'],10)
        self.assertLessEqual(result['bracket_width_s'],1e-7)
        self.assertLessEqual(result['maximum_position_difference_m'],5e-8)
        self.assertFalse(result['activation_delay_within_bound'])

    def test_zero_delay_is_valid_for_visible_endpoint(self):
        model=ForcedParticle()
        def zero(dt):
            row=model(dt)
            if not row['blocked']:row['activation_delay_s']=0.
            return row
        result=locate(zero)
        self.assertEqual(result['status'],BRACKETED)
        self.assertEqual(result['high_activation_delay_s'],0.)

    def test_original_secant_and_bisection_algorithm_are_bytecode_ast_equivalent(self):
        import afsd_event_locator_delay_20261010 as current
        old_path=ROOT/'scripts/afsd_event_locator_secant_20261010.py'
        self.assertEqual(hashlib.sha256(old_path.read_bytes()).hexdigest(),current.FROZEN_PARENT_SHA256)
        def fn(path):
            return next(n for n in ast.parse(path.read_text()).body
                if isinstance(n,ast.FunctionDef) and n.name=='locate_activation')
        old,new=fn(old_path),fn(ROOT/'scripts/afsd_event_locator_delay_20261010.py')
        # Proposal and bracket-update routines are untouched. Only probe
        # validation/evidence and finish acceptance know about own-path delay.
        for name in ('solved','real_pair_at_guess'):
            a=next(n for n in old.body if isinstance(n,ast.FunctionDef) and n.name==name)
            b=next(n for n in new.body if isinstance(n,ast.FunctionDef) and n.name==name)
            self.assertEqual(ast.dump(a),ast.dump(b))
        a=next(n for n in old.body if isinstance(n,ast.While))
        b=next(n for n in new.body if isinstance(n,ast.While))
        self.assertEqual(ast.dump(a),ast.dump(b))


if __name__ == '__main__':
    unittest.main()

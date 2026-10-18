"""Controller invariants; these fixtures do not validate physical deposition."""
import unittest
from run_afsd_deposition_stage_20261007 import (adaptive_pair, bounded_interval,
    next_step_size, screen_utilization, GeometryPathStop, RefinableCandidateError,
    RefinementExhausted, StageBudgetStop)


class DepositionStageControllerTests(unittest.TestCase):
    def test_temporal_failure_retries_same_parent_without_selecting_failed_fine(self):
        parent = object(); parents, dts, rows = [], [], []
        failed_fine, passed_fine = object(), object()
        def attempt(dt, level):
            parents.append(parent); dts.append(dt)
            return dict(screen=dict(passed=level == 1), fine=passed_fine if level == 1 else failed_fine)
        result, dt = adaptive_pair(.001, .0001, 4, attempt, rows.append)
        self.assertEqual(dts, [.001, .0005]); self.assertTrue(all(v is parent for v in parents))
        self.assertIs(result['fine'], passed_fine); self.assertEqual(dt, .0005)
        self.assertEqual([row['passed'] for row in rows], [False, True])

    def test_geometry_failure_stops_without_halving(self):
        calls, rows = [], []
        def attempt(dt, level):
            calls.append(dt); raise GeometryPathStop('path intersects tool')
        with self.assertRaises(GeometryPathStop):
            adaptive_pair(.001, 1e-6, 8, attempt, rows.append)
        self.assertEqual(calls, [.001]); self.assertEqual(rows, [])

    def test_budget_after_failed_attempt_never_returns_an_unscreened_state(self):
        rows = []
        def attempt(dt, level):
            if level == 0:
                return dict(screen=dict(passed=False), fine='FAILED_FINE')
            raise StageBudgetStop('DECLARED_TOTAL_CANDIDATE_QUOTA', dt)
        with self.assertRaises(StageBudgetStop) as caught:
            adaptive_pair(.001, 1e-6, 8, attempt, rows.append)
        self.assertEqual(caught.exception.next_dt, .0005)
        self.assertEqual(len(rows), 1); self.assertFalse(rows[0]['passed'])

    def test_explicit_refinable_candidate_failure_can_retry(self):
        rows = []
        def attempt(dt, level):
            if level == 0:
                raise RefinableCandidateError('force root rejected')
            return dict(screen=dict(passed=True))
        _, dt = adaptive_pair(.001, .0001, 4, attempt, rows.append)
        self.assertEqual(dt, .0005); self.assertEqual(rows[0]['status'], 'FAIL_CANDIDATE_NO_STATE_COMMIT')

    def test_nonrefinable_identity_error_is_not_hidden_by_retry(self):
        def attempt(dt, level):
            raise ValueError('material identity changed')
        with self.assertRaisesRegex(ValueError, 'identity'):
            adaptive_pair(.001, .0001, 4, attempt, lambda r: None)

    def test_minimum_increment_and_retry_count_are_hard_limits(self):
        for minimum, maximum, expected in ((.0004, 8, 2), (.000001, 1, 2)):
            rows = []
            with self.assertRaises(RefinementExhausted):
                adaptive_pair(.001, minimum, maximum, lambda dt, level: dict(screen=dict(passed=False)), rows.append)
            self.assertEqual(len(rows), expected)

    def test_interval_ends_at_phase_boundary_and_never_extrapolates(self):
        self.assertAlmostEqual(bounded_interval(.019, .005, .1, .020, 1e-6), .001)
        with self.assertRaises(ValueError):
            bounded_interval(.020, .005, .1, .020, 1e-6)

    def test_growth_only_when_all_actual_errors_have_headroom_and_is_bounded(self):
        self.assertEqual(next_step_size(.001, .1, .0012, 2.), .0012)
        self.assertEqual(next_step_size(.001, .3, .002, 2.), .001)
        self.assertEqual(next_step_size(.001, .1, .002, 1.), .001)
        with self.assertRaises(ValueError):
            next_step_size(.001, 1.1, .002, 2.)

    def test_local_particle_error_controls_growth_even_with_small_global_norm(self):
        g = dict(position_max_m=0., temperature_max_K=0., stress_relative=.001, plastic_relative=0.)
        l = dict(position_limit_m=5e-6, temperature_limit_K=.5, stress_relative_limit=.1,
            plastic_relative_limit=.2, rows=[dict(position_error_m=0., temperature_error_K=0.,
                                                stress_relative=.08, plastic_relative=0.)])
        self.assertAlmostEqual(screen_utilization(g, l, dict(dx=.0005)), .8)


if __name__ == '__main__':
    unittest.main()

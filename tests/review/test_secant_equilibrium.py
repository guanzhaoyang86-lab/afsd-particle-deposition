"""Independent coupled nonlinear search fixtures, not AFSD physics validation."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import secant_equilibrium as search


class CoupledFixture:
    def __init__(self, mode='smooth', reported_zero=False):
        self.mode = mode; self.reported_zero = reported_zero
        self.state = dict(free=np.array([[True, True, False]]),
            contact_history=np.array([[[.003, -.002, 0.]]]), time_s=.2)
        self.root = np.array([.4, -.3]) if mode == 'smooth' else np.array([.2, .1])
        self.evaluations = 0; self.tangents = 0; self.attempts = []

    def function(self, x):
        if self.mode == 'log':
            if x[0] <= 0.: raise ValueError('log-domain contact-like guard')
            return np.array([np.log(x[0])+.25*x[1], x[1]+.3*x[0]])
        if self.mode == 'singular': return np.array([x[0], 1e-14*x[1]])
        if self.mode == 'constant': return np.ones(2)
        return np.array([x[0]+.2*x[0]*x[1]+.1*x[1]**2,
                         x[1]+.3*np.sin(x[0])+.1*x[0]**2])

    def evaluate(self, displacement, dt):
        u = np.asarray(displacement).copy(); self.evaluations += 1; self.attempts.append(u)
        if u[0, 2] != .125: raise ValueError('prescribed fixture axis changed')
        f = self.function(u[0, :2])
        if self.mode != 'constant': f -= self.function(self.root)
        force = np.zeros((1, 3)); force[0, :2] = f
        maximum = float(np.max(np.abs(f))); limit = 1e-11
        return dict(displacement=u, force=force, force_residual_N=maximum,
            force_limit_N=limit, scaled_residual=0. if self.reported_zero else maximum/limit,
            density_residual=0., thermal_residual_K=0., virtual_work_residual_J=0.,
            material_endpoint_work_J=0., density_iterations=1, thermal_iterations=1,
            wall_seconds=0., dt_s=dt)

    def tangent(self, r):
        self.tangents += 1; x, y = r['displacement'][0, :2]
        J = np.eye(3)
        if self.mode == 'log': J[:2, :2] = [[1/x, .25], [.3, 1.]]
        elif self.mode == 'singular': J[:2, :2] = np.diag([1., 1e-14])
        elif self.mode != 'constant':
            J[:2, :2] = [[1+.2*y, .2*x+.2*y], [.3*np.cos(x)+.2*x, 1.]]
        return J


class SecantEquilibriumTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.directory = Path(self.temp.name)
        self.journals = []
        def journal(directory, trial, r, name):
            self.journals.append(copy.deepcopy(r))
            search.write(directory/(name+'.json'), dict(accepted=False,
                restart_allowed=False, force_N=r['force'], displacement=r['displacement']))
        self.mock = patch.object(search, 'journal', side_effect=journal); self.mock.start()

    def tearDown(self):
        self.mock.stop(); self.temp.cleanup()

    def solve(self, fixture, guess=None, maximum=24):
        if guess is None: guess = np.array([[0., 0., .125]])
        return search.solve(fixture, .1, guess, self.directory, time.perf_counter()+30., maximum)

    def test_coupled_nonlinear_root_uses_full_force_and_preserves_axes_history(self):
        fixture = CoupledFixture(); before = copy.deepcopy(fixture.state)
        guess = np.array([[0., 0., .125]]); original_guess = guess.copy()
        r, metrics = self.solve(fixture, guess)
        np.testing.assert_allclose(r['displacement'][0, :2], fixture.root, atol=1e-10, rtol=0)
        self.assertLessEqual(np.max(np.abs(r['force'][fixture.state['free']])), r['force_limit_N'])
        self.assertGreater(metrics['broyden_updates'], 1)
        self.assertLess(metrics['jacobian_builds'], metrics['iterations'])
        self.assertEqual(metrics['evaluations'], fixture.evaluations)
        self.assertLessEqual(metrics['evaluations'], 100)
        self.assertLessEqual(metrics['iterations'], 24)
        self.assertLessEqual(metrics['jacobian_refreshes'], 2)
        for u in fixture.attempts: self.assertEqual(u[0, 2], .125)
        np.testing.assert_array_equal(guess, original_guess)
        np.testing.assert_array_equal(fixture.state['contact_history'], before['contact_history'])
        np.testing.assert_array_equal(fixture.state['free'], before['free'])
        self.assertEqual(fixture.state['time_s'], before['time_s'])
        self.assertTrue((self.directory/'last_UNCOMMITTED_candidate.json').exists())

    def test_invalid_domain_line_trials_are_logged_counted_and_backtracked(self):
        fixture = CoupledFixture('log')
        r, metrics = self.solve(fixture, np.array([[2., 0., .125]]))
        np.testing.assert_allclose(r['displacement'][0, :2], fixture.root, atol=1e-10, rtol=0)
        rows = json.loads((self.directory/'iterations.json').read_text())
        errors = [entry['rejected'] for row in rows for entry in row.get('line_search', []) if 'rejected' in entry]
        self.assertTrue(any('log-domain' in error for error in errors))
        self.assertEqual(metrics['evaluations'], fixture.evaluations)
        self.assertTrue(any(row.get('accepted_alpha', 1.) < 1. for row in rows))

    def test_stale_reported_convergence_cannot_replace_actual_force(self):
        fixture = CoupledFixture(reported_zero=True)
        r, metrics = self.solve(fixture)
        self.assertGreater(metrics['iterations'], 0)
        self.assertLessEqual(np.max(np.abs(r['force'])), r['force_limit_N'])
        np.testing.assert_allclose(r['displacement'][0, :2], fixture.root, atol=1e-10, rtol=0)

    def test_rejected_root_has_uncommitted_journal_and_unchanged_history(self):
        fixture = CoupledFixture('constant'); before = fixture.state['contact_history'].copy()
        with self.assertRaisesRegex(ValueError, 'Armijo search failed'):
            self.solve(fixture)
        np.testing.assert_array_equal(fixture.state['contact_history'], before)
        self.assertEqual(fixture.state['time_s'], .2)
        failure = json.loads((self.directory/'secant_failure.json').read_text())
        self.assertFalse(failure['accepted']); self.assertFalse(failure['restart_allowed'])
        self.assertEqual(failure['jacobian_refreshes'], 2)
        self.assertLessEqual(failure['evaluations'], 100)
        self.assertEqual(len(self.journals), 1)
        np.testing.assert_array_equal(self.journals[0]['force'][0, :2], [1., 1.])

    def test_near_singular_initial_matrix_is_rejected_without_shift_or_pseudoinverse(self):
        fixture = CoupledFixture('singular')
        with self.assertRaisesRegex(ValueError, 'near-singular'):
            self.solve(fixture)
        self.assertEqual(fixture.evaluations, 1)
        self.assertTrue((self.directory/'secant_failure.json').exists())
        np.testing.assert_array_equal(self.journals[0]['displacement'], [[0., 0., .125]])

    def test_rank_one_singular_update_is_rejected_and_tiny_secant_is_skipped(self):
        with self.assertRaisesRegex(ValueError, 'determinant factor'):
            search._broyden(np.eye(2), np.eye(2), np.array([1., 0.]), np.zeros(2), 1.)
        A = np.array([[2., .2], [.1, 1.]])
        same, _, report = search._broyden(A, np.linalg.inv(A), np.array([1e-20, 0.]),
                                        np.array([2e-20, 1e-21]), 1.)
        np.testing.assert_array_equal(same, A)
        self.assertEqual(report['status'], 'SKIPPED_UNRESOLVED_SECANT_STEP')

    def test_iteration_budget_keeps_last_real_candidate_without_commit(self):
        fixture = CoupledFixture()
        with self.assertRaisesRegex(search.SearchBudgetError, 'iteration budget'):
            self.solve(fixture, maximum=1)
        self.assertEqual(len(self.journals), 1)
        self.assertGreater(np.max(np.abs(self.journals[0]['force'])), 1e-11)
        self.assertEqual(fixture.state['time_s'], .2)

    def test_expired_deadline_does_not_evaluate_or_forge_candidate(self):
        fixture = CoupledFixture()
        with self.assertRaisesRegex(search.SearchBudgetError, 'time budget'):
            search.solve(fixture, .1, np.array([[0., 0., .125]]), self.directory,
                         time.perf_counter()-1.)
        self.assertEqual(fixture.evaluations, 0); self.assertEqual(len(self.journals), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)

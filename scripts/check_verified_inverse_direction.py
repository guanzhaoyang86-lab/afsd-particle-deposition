"""Bounded matrix and manufactured-force qualification, never AFSD physics."""
import argparse
import ast
import copy
import hashlib
import io
import json
import os
from pathlib import Path
import statistics
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import numpy as np
import verified_inverse_direction as linear
import visibility_inverse_equilibrium as search

ROOT = Path(__file__).resolve().parents[1]
MEASUREMENTS = {}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen_broyden_algebra():
    """Execute only the two original pure algebra definitions, not its imports."""
    path = ROOT/'scripts/secant_equilibrium.py'
    tree = ast.parse(path.read_text(encoding='utf-8'))
    selected = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('_condition', '_broyden')]
    if len(selected) != 2: raise ValueError('original algebra definitions missing')
    namespace = dict(np=np, CONDITION_LIMIT=1e12)
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['_broyden']


class DirectionTests(unittest.TestCase):
    def test_nonsymmetric_inverse_application_matches_original_dense(self):
        rng = np.random.default_rng(521)
        u, _ = np.linalg.qr(rng.normal(size=(48, 48)))
        v, _ = np.linalg.qr(rng.normal(size=(48, 48)))
        a = (u*np.geomspace(1e7, 250., 48))@v.T
        h = np.linalg.inv(a); b = rng.normal(size=48)*.05
        d, r = linear.verified_inverse_direction(a, h, b)
        ref = np.linalg.solve(a, b)
        self.assertFalse(r['dense_fallback_used'])
        self.assertLessEqual(np.max(np.abs(a@d-b)), r['linear_backward_limit_N'])
        self.assertLess(np.linalg.norm(d-ref)/np.linalg.norm(ref), 1e-10)
        MEASUREMENTS['nonsymmetric'] = r

    def test_twenty_four_original_broyden_updates(self):
        rng = np.random.default_rng(124)
        u, _ = np.linalg.qr(rng.normal(size=(64, 64)))
        v, _ = np.linalg.qr(rng.normal(size=(64, 64)))
        a = (u*np.geomspace(1e7, 250., 64))@v.T; h = np.linalg.inv(a)
        update = frozen_broyden_algebra(); results = []
        for _ in range(24):
            step = 1e-9*rng.normal(size=64)
            change = a@step + .1*np.linalg.norm(step)*rng.normal(size=64)
            a, h, secant = update(a, h, step, change, 1e-6)
            self.assertEqual(secant['status'], 'GOOD_BROYDEN_RANK_ONE')
            b = rng.normal(size=64)*.05
            d, r = linear.verified_inverse_direction(a, h, b)
            ref = np.linalg.solve(a, b)
            self.assertLessEqual(np.max(np.abs(a@d-b)), r['linear_backward_limit_N'])
            self.assertLess(np.linalg.norm(d-ref)/np.linalg.norm(ref), 1e-10)
            results.append(r)
        MEASUREMENTS['broyden_24_updates'] = dict(steps=24,
            dense_fallbacks=sum(r['dense_fallback_used'] for r in results),
            maximum_backward_error_N=max(r['linear_backward_error_N'] for r in results),
            maximum_condition_estimate=max(r['condition_infinity_estimate'] for r in results))

    def test_one_beneficial_correction(self):
        a = np.eye(5); h = (1.-1e-6)*a; b = np.full(5, .05)
        d, r = linear.verified_inverse_direction(a, h, b)
        self.assertTrue(r['correction_attempted'] and r['correction_beneficial'])
        self.assertFalse(r['dense_fallback_used'])
        self.assertLessEqual(np.max(np.abs(d-b)), 1e-9)
        MEASUREMENTS['beneficial_correction'] = r

    def test_stale_inverse_falls_back_after_only_one_correction(self):
        a = np.diag([1., 2., 3., 4.]); h = .5*np.linalg.inv(a); b = np.full(4, .05)
        with patch.object(linear.np.linalg, 'solve', wraps=np.linalg.solve) as solve:
            d, r = linear.verified_inverse_direction(a, h, b)
            self.assertEqual(solve.call_count, 1)
        self.assertTrue(r['correction_attempted'] and r['dense_fallback_used'])
        np.testing.assert_array_equal(d, np.linalg.solve(a, b))
        MEASUREMENTS['stale_inverse_fallback'] = r

    def test_unhelpful_correction_is_discarded(self):
        a = np.eye(3); b = np.array([.05, -.03, .02])
        d, r = linear.verified_inverse_direction(a, -a, b)
        self.assertTrue(r['correction_attempted'] and r['dense_fallback_used'])
        self.assertFalse(r['correction_beneficial'])
        np.testing.assert_array_equal(d, b)

    def test_two_original_probe_vectors_do_not_replace_rhs_check(self):
        a = np.eye(4); v = np.array([1., 1., -1., -1.])/2
        h = a+.1*np.outer(np.array([1., 0., 0., 0.]), v)
        probe = np.column_stack((np.ones(4), (-1.)**np.arange(4)))
        np.testing.assert_allclose(a@h@probe, probe, rtol=0, atol=1e-16)
        d, r = linear.verified_inverse_direction(a, h, .05*v)
        self.assertTrue(r['dense_fallback_used'])
        np.testing.assert_array_equal(d, .05*v)
        MEASUREMENTS['probe_blind_stale_inverse'] = r

    def test_existing_condition_limit_rejects_without_regularization(self):
        a = np.diag([1., 1e-13]); h = np.diag([1., 1e13])
        with patch.object(linear.np.linalg, 'solve', side_effect=AssertionError('should not solve')):
            with self.assertRaisesRegex(ValueError, 'near-singular'):
                linear.verified_inverse_direction(a, h, np.ones(2))

    def test_singular_matrix_cannot_be_hidden_by_stale_inverse(self):
        with self.assertRaises(np.linalg.LinAlgError):
            linear.verified_inverse_direction(np.ones((2, 2)), np.eye(2), np.array([1., -1.]))

    def test_readonly_inputs_preserved_in_success_and_fallback(self):
        a = np.array([[2., 1.], [-.2, 3.]]); b = np.array([.05, .08])
        for factor in (1., .5):
            h = factor*np.linalg.inv(a)
            before = [v.copy() for v in (a, h, b)]
            for v in (a, h, b): v.setflags(write=False)
            linear.verified_inverse_direction(a, h, b)
            for original, actual in zip(before, (a, h, b)):
                np.testing.assert_array_equal(original, actual)

    def test_expired_deadline_does_not_start_dense_fallback(self):
        with patch.object(linear.np.linalg, 'solve', side_effect=AssertionError('should not solve')):
            with self.assertRaises(linear.LinearDirectionBudgetError):
                linear.verified_inverse_direction(np.eye(2), np.eye(2), np.ones(2), deadline=time.perf_counter()-1)

    def test_nonfinite_and_wrong_shape_rejected(self):
        for h, b in ((np.full((2, 2), np.nan), np.ones(2)), (np.eye(3), np.ones(2)), (np.eye(2), np.ones((2, 1)))):
            with self.assertRaises(ValueError): linear.verified_inverse_direction(np.eye(2), h, b)

    def test_static_solver_changes_only_direction_route_and_diagnostics(self):
        old = (ROOT/'scripts/visibility_equilibrium.py').read_text(encoding='utf-8')
        new = (ROOT/'scripts/visibility_inverse_equilibrium.py').read_text(encoding='utf-8')
        old = old[old.index('from pathlib import Path'):]
        new = new[new.index('from pathlib import Path'):]
        new = new.replace('from verified_inverse_direction import verified_inverse_direction\n', '', 1)
        new = new.replace('BOUNDED_BROYDEN_VERIFIED_INVERSE_VISIBLE_GRAPH_CHANGE_V1',
                          'BOUNDED_BROYDEN_DISCARD_ON_VISIBLE_GRAPH_CHANGE_V1', 1)
        new = new.replace("                    direction, verified_linear = verified_inverse_direction(\n"
            "                        matrix, inverse, rhs, deadline=deadline)\n"
            "                    backward = verified_linear['linear_backward_error_N']",
            "                    direction = np.linalg.solve(matrix, rhs)\n"
            "                    backward = float(np.max(np.abs(matrix@direction-rhs)))\n"
            "                    if not np.isfinite(direction).all() or backward > 1e-9*max(1., float(np.max(np.abs(rhs)))):\n"
            "                        raise ValueError('secant linear solve backward error')", 1)
        new = new.replace("                    condition_policy='INITIAL_INVERSE_NORM_THEN_SHERMAN_MORRISON_GUARDED_ESTIMATE',\n"
            "                    verified_inverse_direction=verified_linear))",
            "                    condition_policy='INITIAL_INVERSE_NORM_THEN_SHERMAN_MORRISON_GUARDED_ESTIMATE'))", 1)
        self.assertEqual(new, old)


class ManufacturedForce:
    def __init__(self, oscillating=False, forged_zero=False):
        self.state = {'free':np.array([[True, False, False]]), 'time_s':0., 'history':np.array([3.])}
        self.oscillating, self.forged_zero, self.calls, self.builds = oscillating, forged_zero, 0, []
    def evaluate(self, displacement, dt):
        u = np.array(displacement, copy=True); self.calls += 1; x = float(u[0, 0])
        f = np.exp(x)-3 if self.oscillating else x+2*max(x-.5, 0)-1.4
        graph = str(self.calls) if self.oscillating else ('left' if x <= .5 else 'right')
        return dict(displacement=u, force=np.array([[f, 0., 0.]]), force_limit_N=1e-11,
            scaled_residual=0. if self.forged_zero else abs(f)/1e-11, graph_fingerprint=graph, dt_s=dt)
    def tangent(self, result):
        x = float(result['displacement'][0, 0]); self.builds.append(result['graph_fingerprint'])
        return np.diag([np.exp(x) if self.oscillating else 1. if x <= .5 else 3., 1., 1.])


class ManufacturedGraphTests(unittest.TestCase):
    def run_case(self, trial):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        path = Path(self.temp.name)
        # Replace output-only adapters on this NEW solver, never old globals.
        with patch.object(search, 'journal'), patch.object(search, 'summary', return_value={}):
            result = search.solve(trial, .1, np.array([[0., .2, .3]]), path, time.perf_counter()+20, 24)
        return result, json.loads((path/'iterations.json').read_text())
    def test_graph_change_discards_update_and_rebuilds(self):
        trial = ManufacturedForce(); before = copy.deepcopy(trial.state)
        (result, metrics), rows = self.run_case(trial)
        self.assertAlmostEqual(result['displacement'][0, 0], .8, places=11)
        self.assertEqual(trial.builds, ['left', 'right'])
        self.assertEqual(metrics['jacobian_refreshes'], 1)
        self.assertEqual(sum(r.get('secant', {}).get('status') == 'DISCARDED_ACROSS_GRAPH_CHANGE' for r in rows), 1)
        np.testing.assert_array_equal(result['displacement'][0, 1:], [.2, .3])
        np.testing.assert_array_equal(trial.state['history'], before['history'])
    def test_full_force_overrides_forged_convergence(self):
        (result, metrics), _ = self.run_case(ManufacturedForce(forged_zero=True))
        self.assertGreater(metrics['iterations'], 0)
        self.assertLessEqual(abs(result['force'][0, 0]), 1e-11)
    def test_refresh_budget_exhaustion_still_rejects(self):
        trial = ManufacturedForce(oscillating=True)
        with self.assertRaisesRegex(ValueError, 'refresh budget exhausted'): self.run_case(trial)
        self.assertEqual(len(trial.builds), 3)
        report = json.loads((Path(self.temp.name)/'secant_failure.json').read_text())
        self.assertFalse(report['accepted']); self.assertFalse(report['restart_allowed'])
        self.assertEqual(report['jacobian_refreshes'], 2)
        self.assertLessEqual(report['evaluations'], 100)


def small_timing():
    """Small algebra timing only; do not extrapolate this to 1961-DOF physics."""
    rng = np.random.default_rng(1961); n = 192
    a = rng.normal(size=(n, n)); a += 2*n*np.eye(n)
    h = np.linalg.inv(a); b = rng.normal(size=n)*.05
    def dense():
        d = np.linalg.solve(a, b)
        error = float(np.max(np.abs(a@d-b)))
        if error > 1e-9*max(1., float(np.max(np.abs(b)))): raise AssertionError('dense backward gate')
        return d
    dense(); linear.verified_inverse_direction(a, h, b)
    values = {'dense_solve_and_backward_check':[], 'verified_inverse_with_condition_and_backward_check':[]}
    for _ in range(9):
        start = time.perf_counter(); dense(); values['dense_solve_and_backward_check'].append(time.perf_counter()-start)
        start = time.perf_counter(); linear.verified_inverse_direction(a, h, b)
        values['verified_inverse_with_condition_and_backward_check'].append(time.perf_counter()-start)
    return dict(matrix_size=n, repetitions=9, includes_inverse_construction=False,
        purpose='SMALL_ALGEBRA_ONLY_NOT_PROCESS_SPEEDUP', timings=values,
        median_seconds={k:statistics.median(v) for k,v in values.items()})


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); args = parser.parse_args()
    output = (ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT): raise ValueError('output outside project')
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), ROOT/'scripts/verified_inverse_direction.py', ROOT/'scripts/visibility_inverse_equilibrium.py',
             ROOT/'scripts/visibility_equilibrium.py', ROOT/'scripts/secant_equilibrium.py']
    hashes = {p.relative_to(ROOT).as_posix():digest(p) for p in files}
    begin = time.perf_counter(); stream = io.StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (output/'tests.log').write_text(stream.getvalue(), encoding='utf-8')
    timing = small_timing() if result.wasSuccessful() else None
    unchanged = all(digest(ROOT/path) == expected for path, expected in hashes.items())
    passed = result.wasSuccessful() and unchanged
    report = dict(status='PASS_ALGEBRA_AND_MANUFACTURED_SEARCH_ONLY' if passed else 'FAIL_COMPONENT_QUALIFICATION',
        exit_code=0 if passed else 2, test_count=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        skipped=len(result.skipped), wall_seconds=time.perf_counter()-begin, source_sha256=hashes,
        source_unchanged=unchanged, measurements=MEASUREMENTS, small_matrix_timing=timing,
        physical_solver_evaluations=0, nnps_evaluations=0, process_time_advanced_s=0.,
        maximum_matrix_test_dimension=192, manufactured_scalar_force_tests=3,
        numpy_version=np.__version__, thread_environment={k:os.environ.get(k) for k in
            ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
        limitations=['No actual AFSD residual, tangent, process root or checkpoint acceptance was run.',
            'Small matrix timings do not measure 1961-DOF speedup or amortize inverse construction.',
            'Actual same-parent A/B root qualification is still required before process use.'])
    (output/'receipt.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(stream.getvalue()); print(json.dumps({k:report[k] for k in ('status','exit_code','test_count','wall_seconds')}, indent=2))
    return report['exit_code']


if __name__ == '__main__': raise SystemExit(main())

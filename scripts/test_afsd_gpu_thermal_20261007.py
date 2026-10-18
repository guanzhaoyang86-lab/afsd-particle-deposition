"""Static adapter tests and optional explicit CPU/GPU thermal components.

The default command requires no CUDA. --gpu really executes the GPU component;
it must fail if CUDA is absent and never falls back silently to CPU. These are
small thermal fixtures, not a full AFSD time-step or acceleration claim.

Comparison revision 3 keeps the v1/v2 failed receipts. Near-zero rounding residuals
need not equal one another: CPU and GPU must each meet the original physical
heat-balance limit. All temperature and per-particle heat comparisons are
unchanged. Global ledger comparisons propagate those per-particle differences
through fsum, because cancellation can leave a near-zero global sum.
"""
from pathlib import Path
import argparse
import ast
import copy
import inspect
import json
import math
import sys
import textwrap
import time
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import implicit_thermal as cpu
import afsd_gpu_thermal_20261007 as gpu
import afsd_gpu_trial_20261007 as trial_adapter


def arguments(n=2):
    x = np.column_stack((np.arange(n), np.zeros(n), np.zeros(n)))
    return dict(mass=np.ones(n), rho=np.ones(n), T=np.linspace(300., 400., n),
        conductivity=.5, x=x, i=np.arange(n-1, dtype=int), j=np.arange(1, n, dtype=int),
        gradient=np.tile([1., 0., 0.], (n-1, 1)), dt=1., cp=1.,
        plastic_work_J=np.zeros(n), friction_raw_J=np.zeros(n),
        friction_rounding_bound_J=np.zeros(n), wall_row_W_K=np.zeros(n),
        wall_weighted_temperature_W=np.zeros(n), chi=.9, partition=.25,
        room=293.15, melt=855.15)


def arguments_from_trial_result(trial, result):
    """Read-only final thermal fixed-point inputs of a completed residual map."""
    p, pa, c = trial.p, trial.pa, result['contact']
    conductance = np.asarray(c['active'], bool).astype(float)*p['contact_conductance']*p['dx']**2
    bath = np.full(8, p['tool_temperature']); bath[-1] = p['base_temperature']
    return dict(mass=pa.m.copy(), rho=result['density'].copy(), T=pa.temperature.copy(),
        conductivity=p['conductivity'], x=result['position'].copy(),
        i=result['i'].copy(), j=result['j'].copy(), gradient=result['g'].copy(),
        dt=result['dt_s'], cp=p['cp'], plastic_work_J=result['volume']*result['material'][2],
        friction_raw_J=np.asarray(c['friction_dissipation_J']).sum(axis=1),
        friction_rounding_bound_J=np.zeros(len(pa.x)), wall_row_W_K=conductance.sum(axis=1),
        wall_weighted_temperature_W=conductance@bath, chi=p['plastic_heat_fraction'],
        partition=p['heat_partition'], room=p['room'], melt=p['melt'])


def check_residuals(result):
    """Check residual magnitudes against heat amounts, never another residual."""
    particle = result['particle_energy_J']
    total = np.asarray(particle['total_heat_input_J'])
    # Exactly the existing production per-particle physical heat-balance gate.
    balance_limit = 1e-12+1e-10*float(np.abs(total).sum())
    residual = np.asarray(particle['mc_delta_T_residual_J'])
    if residual.shape != total.shape or not np.isfinite(residual).all():
        raise AssertionError('invalid particle thermal balance residual')
    if float(np.max(np.abs(residual))) > balance_limit:
        raise AssertionError('particle thermal balance residual exceeds original limit')
    ledger_residual = float(result['ledger']['mc_delta_T_residual_J'])
    if not np.isfinite(ledger_residual) or abs(ledger_residual) > balance_limit:
        raise AssertionError('ledger thermal balance residual exceeds original limit')
    # The matrix residual is also in joules. Its scale is the signed source
    # ledger's absolute heat amounts, before cancellation, with the same
    # absolute/relative coefficients. This is a test-only backward check;
    # no production equation or heat-balance threshold is changed.
    source_scale = sum(float(np.abs(particle[key]).sum()) for key in
        ('plastic_heat_J', 'friction_material_heat_J', 'contact_in_J', 'internal_conduction_J'))
    linear_limit = 1e-12+1e-10*source_scale
    linear = np.asarray(result['linear_system_residual_J'])
    if linear.shape != total.shape or not np.isfinite(linear).all():
        raise AssertionError('invalid linear thermal residual')
    if float(np.max(np.abs(linear))) > linear_limit:
        raise AssertionError('linear thermal residual exceeds source-scaled limit')
    return dict(heat_balance_limit_J=balance_limit,
        particle_residual_max_J=float(np.max(np.abs(residual))), ledger_residual_J=ledger_residual,
        linear_residual_max_J=float(np.max(np.abs(linear))), linear_residual_limit_J=linear_limit)


def compare_ledgers(reference, actual):
    """Validate each own fsum, then propagate already-checked particle errors."""
    maximum = 0.
    for field in reference['ledger']:
        terms = (('friction_wall_heat_J', 1.), ('contact_in_J', -1.)) if field == 'wall_heat_out_J' else ((field, 1.),)
        source_scales = []
        for result in (reference, actual):
            arrays = [np.asarray(result['particle_energy_J'][name]) for name, _ in terms]
            scale = math.fsum(float(value) for values in arrays for value in np.abs(values))
            own_sum = math.fsum(sign*math.fsum(float(value) for value in values)
                               for values, (_, sign) in zip(arrays, terms))
            recorded = float(result['ledger'][field])
            rounding = 64*np.finfo(float).eps*max(1., scale)
            if not np.isfinite(recorded) or abs(recorded-own_sum) > rounding:
                raise AssertionError('ledger differs from its own particle fsum: '+field)
            source_scales.append(scale)
        # |sum(a)-sum(b)| <= sum(|a-b|), plus sum/subtraction rounding.
        # No relative comparison with a possibly cancelled global sum is used.
        propagation = math.fsum(float(value) for name, _ in terms for value in
            np.abs(np.asarray(reference['particle_energy_J'][name])-
                   np.asarray(actual['particle_energy_J'][name])))
        limit = propagation+64*np.finfo(float).eps*max(1., *source_scales)
        difference = abs(float(reference['ledger'][field])-float(actual['ledger'][field]))
        if difference > limit:
            raise AssertionError('ledger difference exceeds propagated particle error: '+field)
        maximum = max(maximum, difference)
    return maximum


def compare_results(reference, actual):
    """Compare heat/temperature; independently check both rounding residuals."""
    if set(reference) != set(actual):
        raise AssertionError('thermal output keys differ')
    differences = dict(residual_checks=dict(cpu=check_residuals(reference), gpu=check_residuals(actual)))
    for key in reference:
        a, b = reference[key], actual[key]
        if isinstance(a, dict):
            if set(a) != set(b):
                raise AssertionError('thermal ledger keys differ: '+key)
            if key == 'ledger':
                differences[key] = compare_ledgers(reference, actual)
                continue
            for field in a:
                if field == 'mc_delta_T_residual_J':
                    continue  # Each implementation was checked against the original balance limit.
                np.testing.assert_allclose(a[field], b[field], rtol=1e-10, atol=1e-12,
                                           err_msg=key+'.'+field)
            differences[key] = max(float(np.max(np.abs(np.asarray(a[f])-np.asarray(b[f])))) for f in a)
        elif isinstance(a, np.ndarray):
            atol = 1e-9 if key == 'temperature' else 1e-12
            if key != 'linear_system_residual_J':
                np.testing.assert_allclose(a, b, rtol=1e-10, atol=atol, err_msg=key)
            differences[key] = float(np.max(np.abs(a-b))) if a.size else 0.
        elif isinstance(a, (str, bool)):
            if a != b:
                raise AssertionError('thermal metadata differs: '+key)
        else:
            np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-12, err_msg=key)
    return differences


def compare_fixture(params):
    before = copy.deepcopy(params)
    started = time.perf_counter(); reference = cpu.advance_implicit_thermal(**params)
    cpu_seconds = time.perf_counter()-started
    started = time.perf_counter(); actual = gpu.advance_implicit_thermal(**params)
    gpu_seconds = time.perf_counter()-started
    differences = compare_results(reference, actual)
    for key in params:
        np.testing.assert_array_equal(params[key], before[key], err_msg='input mutated: '+key)
    return dict(status='PASS_THERMAL_CPU_GPU_COMPONENT_EQUIVALENCE', n=len(params['T']),
        cpu_wall_s=cpu_seconds, gpu_wall_s=gpu_seconds, differences=differences,
        accepted=False, physical_time_advanced_s=0., speed_claim=False)


class StaticAdapterTests(unittest.TestCase):
    def test_only_one_thermal_linear_provider_differs(self):
        self.assertEqual(gpu.static_adaptation_audit()['gpu_solve_calls'], 1)

    def test_trial_equations_and_other_providers_unchanged(self):
        self.assertTrue(trial_adapter.operator_reuse_audit()['endpoint_equations_ast_identical'])

    def test_changed_heat_balance_tolerance_is_rejected(self):
        tree = ast.parse(textwrap.dedent(inspect.getsource(gpu.advance_implicit_thermal))).body[0]
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == 1e-12:
                node.value = 1e-11; break
        else:
            self.fail('thermal balance tolerance not found')
        with self.assertRaisesRegex(ValueError, 'thermal equations/gates differ'):
            gpu.static_adaptation_audit(tree)

    def test_changed_thermal_solver_arguments_are_rejected(self):
        tree = ast.parse(textwrap.dedent(inspect.getsource(gpu.advance_implicit_thermal))).body[0]
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == '_gpu_solve':
                node.args = list(reversed(node.args)); break
        with self.assertRaisesRegex(ValueError, 'thermal equations/gates differ'):
            gpu.static_adaptation_audit(tree)

    def test_changed_force_tolerance_is_rejected(self):
        tree = ast.parse(textwrap.dedent(inspect.getsource(trial_adapter.ProcessTrial.evaluate))).body[0]
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == 1e-5:
                node.value = 1e-4; break
        else:
            self.fail('force tolerance not found')
        with self.assertRaisesRegex(ValueError, 'GPU endpoint residual AST differs'):
            trial_adapter.operator_reuse_audit(tree)

    def test_trial_remains_compatible_with_original_candidate_guards(self):
        self.assertTrue(issubclass(trial_adapter.ProcessTrial, trial_adapter.parent.ProcessTrial))
        self.assertIs(trial_adapter.advance_candidate, trial_adapter.parent.advance_candidate)


class CPUFixtureTests(unittest.TestCase):
    def test_independent_two_body_literal_reference(self):
        r = cpu.advance_implicit_thermal(**arguments())
        np.testing.assert_allclose(r['temperature'], [1000/3, 1100/3], rtol=0, atol=2e-12)

    def test_comparison_detects_temperature_change(self):
        r = cpu.advance_implicit_thermal(**arguments()); changed = copy.deepcopy(r)
        changed['temperature'][0] += 1e-4
        with self.assertRaises(AssertionError):
            compare_results(r, changed)

    def test_comparison_detects_energy_change(self):
        r = cpu.advance_implicit_thermal(**arguments()); changed = copy.deepcopy(r)
        changed['particle_energy_J']['total_heat_input_J'][0] += 1e-7
        with self.assertRaises(AssertionError):
            compare_results(r, changed)

    def test_comparison_detects_each_particle_balance_residual(self):
        r = cpu.advance_implicit_thermal(**arguments()); changed = copy.deepcopy(r)
        changed['particle_energy_J']['mc_delta_T_residual_J'][0] = 2*check_residuals(r)['heat_balance_limit_J']
        for pair in ((r, changed), (changed, r)):
            with self.assertRaisesRegex(AssertionError, 'particle thermal balance residual'):
                compare_results(*pair)

    def test_comparison_detects_each_ledger_balance_residual(self):
        r = cpu.advance_implicit_thermal(**arguments()); changed = copy.deepcopy(r)
        changed['ledger']['mc_delta_T_residual_J'] = 2*check_residuals(r)['heat_balance_limit_J']
        for pair in ((r, changed), (changed, r)):
            with self.assertRaisesRegex(AssertionError, 'ledger thermal balance residual'):
                compare_results(*pair)

    def test_comparison_detects_each_linear_residual(self):
        r = cpu.advance_implicit_thermal(**arguments()); changed = copy.deepcopy(r)
        changed['linear_system_residual_J'][0] = 2*check_residuals(r)['linear_residual_limit_J']
        for pair in ((r, changed), (changed, r)):
            with self.assertRaisesRegex(AssertionError, 'linear thermal residual'):
                compare_results(*pair)

    def test_comparison_detects_each_ledger_fsum_change(self):
        r = cpu.advance_implicit_thermal(**arguments())
        for field in ('measured_thermal_increment_J', 'wall_heat_out_J'):
            changed = copy.deepcopy(r)
            changed['ledger'][field] += 1e-7
            for pair in ((r, changed), (changed, r)):
                with self.assertRaisesRegex(AssertionError, 'own particle fsum'):
                    compare_results(*pair)


GPU_REPORTS = []


class GPUFixtureTests(unittest.TestCase):
    def test_two_body_and_unequal_capacity(self):
        for cp in (1., np.array([2., 1.])):
            params = arguments(); params['cp'] = cp
            GPU_REPORTS.append(compare_fixture(params))

    def test_single_wall_with_heat_sources(self):
        params = arguments(1)
        params.update(cp=2., dt=.5, plastic_work_J=np.array([4.]), chi=1.,
            wall_row_W_K=np.array([3.]), wall_weighted_temperature_W=np.array([1500.]))
        GPU_REPORTS.append(compare_fixture(params))

    def test_large_dt_and_particle_permutation(self):
        params = arguments(3); params['T'] = np.array([400., 300., 500.]); params['dt'] = 100.
        GPU_REPORTS.append(compare_fixture(params))
        perm = np.array([2, 0, 1]); inv = np.argsort(perm); permuted = copy.deepcopy(params)
        for key in ('mass', 'rho', 'T', 'x', 'plastic_work_J', 'friction_raw_J',
                    'friction_rounding_bound_J', 'wall_row_W_K', 'wall_weighted_temperature_W'):
            permuted[key] = permuted[key][perm]
        permuted['i'] = inv[params['i']]; permuted['j'] = inv[params['j']]
        GPU_REPORTS.append(compare_fixture(permuted))

    def test_signed_roundoff_friction_and_rejection(self):
        params = arguments(1); params['friction_raw_J'][0] = -1e-26
        params['friction_rounding_bound_J'][0] = 2e-26
        GPU_REPORTS.append(compare_fixture(params))
        params['friction_rounding_bound_J'][0] = 0.
        with self.assertRaisesRegex(ValueError, 'not clipped'):
            gpu.advance_implicit_thermal(**params)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--gpu', action='store_true')
    parser.add_argument('--output', type=Path); args = parser.parse_args()
    if args.output and args.output.exists():
        raise ValueError('refuse existing receipt overwrite')
    suite = unittest.TestSuite()
    for cls in (StaticAdapterTests, CPUFixtureTests) + ((GPUFixtureTests,) if args.gpu else ()):
        suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(cls))
    started = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    receipt = dict(status='PASS' if result.wasSuccessful() else 'FAIL', tests_run=result.testsRun,
        gpu_requested=args.gpu, gpu_fixture_reports=GPU_REPORTS, wall_seconds=time.perf_counter()-started,
        accepted=False, physical_time_advanced_s=0.,
        scope='Static adapters and explicitly selected thermal components only; no AFSD root or speed claim.')
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())

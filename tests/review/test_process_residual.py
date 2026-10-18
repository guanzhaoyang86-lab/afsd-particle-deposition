"""Independent algebra/time references, NOT an AFSD trajectory or solver test."""
import inspect
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from afsd_backend.pysph.process_residual import midpoint_residual


def fixture(count=1):
    z = np.zeros((count, 3))
    return dict(mass=np.ones(count), x0=z.copy(), x1=z.copy(),
                v0=z.copy(), v1=z.copy(), dt=.2, internal_force=z.copy(),
                contact_force=z.copy(), external_force=z.copy(),
                reaction_force=z.copy(), free_axes=np.ones((count, 3)),
                prescribed_velocity=z.copy())


def spring_fixture():
    a = fixture()
    a.update(x0=np.array([[1., 0., -1.]]), v0=np.array([[0., 2., 1.]]),
             v1=np.array([[-35 / 52, 162 / 109, 4.]]),
             external_force=np.array([[.5, -1., 3.]]))
    a['x1'] = a['x0'] + .1 * (a['v0'] + a['v1'])
    a['internal_force'] = -np.array([4., 9., 16.]) * (a['x0'] + a['x1']) / 2
    return a


class ProcessResidual(unittest.TestCase):
    def close(self, actual, expected):
        np.testing.assert_allclose(actual, expected, rtol=0., atol=1e-12)

    def test_01_three_dimensional_spring_literal_solution(self):
        result = midpoint_residual(**spring_fixture())
        self.close(result['kinematic_residual_m'], 0.)
        self.close(result['dynamic_residual_N'], 0.)
        self.close(result['constraint_velocity_residual_m_s'], 0.)
        self.assertTrue(np.all(np.abs(result['inertial_force_N']) > 0))

    def test_02_independent_spring_energy_external_work(self):
        a = spring_fixture()
        k = np.array([4., 9., 16.])
        e0 = .5 * np.sum(a['v0'] ** 2 + k * a['x0'] ** 2)
        e1 = .5 * np.sum(a['v1'] ** 2 + k * a['x1'] ** 2)
        work = np.sum(a['external_force'] * (a['x1'] - a['x0']))
        self.assertAlmostEqual(e1 - e0, work, delta=1e-12)
        self.close(midpoint_residual(**a)['dynamic_residual_N'], 0.)

    def test_03_literal_three_axis_constant_force(self):
        a = fixture()
        a.update(mass=np.array([2.]), dt=.5,
                 internal_force=np.array([[1., -2., 3.]]),
                 contact_force=np.array([[3., 2., -1.]]),
                 external_force=np.array([[4., 4., -6.]]),
                 v1=np.array([[2., 1., -1.]]), x1=np.array([[.5, .25, -.25]]))
        r = midpoint_residual(**a)
        self.close(r['physical_force_N'], [[8., 4., -4.]])
        self.close(r['inertial_force_N'], [[8., 4., -4.]])
        self.close(r['dynamic_residual_N'], 0.)
        self.close(r['kinematic_residual_m'], 0.)
        self.close(r['quasistatic_free_residual_N'], [[-8., -4., 4.]])

    def test_04_constrained_axis_reaction_is_not_dropped(self):
        a = fixture()
        a.update(dt=.5, mass=np.array([2.]), free_axes=np.array([[1, 1, 0]]),
                 v0=np.array([[0., 0., -.5]]), v1=np.array([[1., -2., -.5]]),
                 x1=np.array([[.25, -.5, -.25]]),
                 external_force=np.array([[4., -8., 6.]]),
                 reaction_force=np.array([[0., 0., -6.]]),
                 prescribed_velocity=np.array([[99., 99., -.5]]))
        r = midpoint_residual(**a)
        self.close(r['dynamic_residual_N'], 0.)
        self.close(r['kinematic_residual_m'], 0.)
        self.close(r['constraint_velocity_residual_m_s'], 0.)
        self.close(r['free_reaction_force_N'], 0.)
        self.close(r['quasistatic_free_residual_N'], [[-4., 8., 0.]])

    def test_05_free_axis_reaction_is_exposed_not_projected_away(self):
        a = fixture()
        a['reaction_force'] = np.array([[3., -2., 5.]])
        a['free_axes'] = np.array([[1, 0, 1]])
        r = midpoint_residual(**a)
        self.close(r['free_reaction_force_N'], [[3., 0., 5.]])
        self.close(r['dynamic_residual_N'], [[-3., 2., -5.]])
        self.close(r['quasistatic_free_residual_N'], 0.)

    def test_06_literal_static_equilibrium(self):
        a = fixture()
        a['internal_force'] = np.array([[2., -3., 4.]])
        a['contact_force'] = np.array([[-1., 1., -1.]])
        a['external_force'] = np.array([[-1., 2., -3.]])
        r = midpoint_residual(**a)
        self.close(r['dynamic_residual_N'], 0.)
        self.close(r['quasistatic_free_residual_N'], 0.)

    def test_07_quasistatic_zero_does_not_remove_unbalanced_inertia(self):
        a = fixture()
        a['v1'] = np.array([[.2, -.4, .6]])
        r = midpoint_residual(**a)
        self.close(r['quasistatic_free_residual_N'], 0.)
        self.close(r['dynamic_residual_N'], [[1., -2., 3.]])

    def test_08_endpoint_constraint_error_keeps_sign(self):
        a = fixture()
        a['free_axes'] = np.array([[1, 0, 0]])
        a['v1'] = np.array([[9., 2., -3.]])
        a['prescribed_velocity'] = np.array([[-9., 4., -4.]])
        self.close(midpoint_residual(**a)['constraint_velocity_residual_m_s'], [[0., -2., 1.]])

    def test_09_particle_permutation_equivariance(self):
        rng = np.random.default_rng(121)
        a = fixture(5)
        for name in a:
            if name not in ('mass', 'dt', 'free_axes'):
                a[name] = rng.normal(size=(5, 3))
        a['mass'] = np.arange(1., 6.)
        a['free_axes'] = rng.integers(0, 2, size=(5, 3))
        order = np.array([3, 0, 4, 2, 1])
        b = {key: (value if key == 'dt' else value[order].copy()) for key, value in a.items()}
        ra, rb = midpoint_residual(**a), midpoint_residual(**b)
        for key in ra:
            self.close(rb[key], ra[key][order])

    def test_10_inputs_unchanged_and_readonly_arrays_supported(self):
        a = spring_fixture()
        copies = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in a.items()}
        for v in a.values():
            if isinstance(v, np.ndarray):
                v.setflags(write=False)
        midpoint_residual(**a)
        for key in copies:
            np.testing.assert_array_equal(a[key], copies[key])

    def test_11_output_arrays_do_not_alias_inputs_or_each_other(self):
        a = spring_fixture()
        r = midpoint_residual(**a)
        values = list(r.values())
        for i, out in enumerate(values):
            for inp in a.values():
                if isinstance(inp, np.ndarray):
                    self.assertFalse(np.shares_memory(out, inp))
            for other in values[i + 1:]:
                self.assertFalse(np.shares_memory(out, other))

    def test_12_bool_and_exact_numeric_projection_equivalent(self):
        a = spring_fixture()
        a['free_axes'] = np.array([[True, False, True]])
        r = midpoint_residual(**a)
        a['free_axes'] = np.array([[1., 0., 1.]])
        s = midpoint_residual(**a)
        for key in r:
            np.testing.assert_array_equal(r[key], s[key])

    def test_13_interface_has_no_tags_thermal_or_commit(self):
        expected = ('mass', 'x0', 'x1', 'v0', 'v1', 'dt', 'internal_force',
                    'contact_force', 'external_force', 'reaction_force',
                    'free_axes', 'prescribed_velocity')
        self.assertEqual(tuple(inspect.signature(midpoint_residual).parameters), expected)
        expected_output = {'kinematic_residual_m', 'dynamic_residual_N',
                           'quasistatic_free_residual_N', 'constraint_velocity_residual_m_s',
                           'free_reaction_force_N', 'inertial_force_N', 'physical_force_N'}
        self.assertEqual(set(midpoint_residual(**fixture())), expected_output)

    def test_14_array_shape_errors_rejected_without_broadcast(self):
        keys = ('x0', 'x1', 'v0', 'v1', 'internal_force', 'contact_force',
                'external_force', 'reaction_force', 'prescribed_velocity')
        for key in keys:
            for bad in (0., np.zeros(3), np.zeros((1, 2)), np.zeros((2, 3))):
                with self.subTest(key=key, shape=np.shape(bad)):
                    a = fixture(); a[key] = bad
                    with self.assertRaises(ValueError): midpoint_residual(**a)

    def test_15_mass_requires_positive_nonempty_N_vector(self):
        for bad in (1., [], [[1.]], [0.], [-1.], [1., 2.], [True], ['1']):
            with self.subTest(bad=bad):
                a = fixture(); a['mass'] = bad
                with self.assertRaises(ValueError): midpoint_residual(**a)

    def test_16_dt_requires_positive_real_scalar(self):
        for bad in (0., -1., [1.], True, '1', 1. + 0.j, np.inf, np.nan):
            with self.subTest(bad=bad):
                a = fixture(); a['dt'] = bad
                with self.assertRaises(ValueError): midpoint_residual(**a)

    def test_17_projection_requires_exact_zero_one_shape(self):
        for bad in (1., [1, 1, 1], [[1, 1, .5]], [[1, 2, 0]], [[-1, 1, 0]],
                    [[1, 1, np.nan]], [[1, 1, np.inf]], [['1', '0', '1']],
                    [[1 + 0j, 0j, 1 + 0j]], np.eye(3)):
            with self.subTest(bad=bad):
                a = fixture(); a['free_axes'] = bad
                with self.assertRaises(ValueError): midpoint_residual(**a)

    def test_18_every_particle_input_rejects_nonfinite(self):
        for key in fixture():
            if key == 'dt': continue
            for value in (np.nan, np.inf, -np.inf):
                with self.subTest(key=key, value=value):
                    a = fixture(); a[key].flat[0] = value
                    with self.assertRaises(ValueError): midpoint_residual(**a)

    def test_19_complex_string_object_arrays_rejected(self):
        for key in ('x0', 'x1', 'v0', 'v1', 'internal_force', 'contact_force',
                    'external_force', 'reaction_force', 'prescribed_velocity'):
            for dtype in (complex, str, object, bool):
                with self.subTest(key=key, dtype=dtype):
                    a = fixture(); a[key] = np.zeros((1, 3), dtype=dtype)
                    with self.assertRaises(ValueError): midpoint_residual(**a)

    def test_20_arithmetic_overflow_is_explicit_failure(self):
        a = fixture()
        a['v1'][0, 0] = 1e308
        a['mass'][0] = 1e308
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            midpoint_residual(**a)

    def test_21_same_candidate_repeated_calls_are_pure(self):
        a = spring_fixture()
        first = midpoint_residual(**a)
        midpoint_residual(**fixture())
        second = midpoint_residual(**a)
        for key in first:
            np.testing.assert_array_equal(first[key], second[key])

    def test_22_independent_oscillator_second_order_phase_and_energy(self):
        reports = []
        # Reference state advances are obtained from an independently written
        # two-equation linear system, NOT by solving a residual from this API.
        for h in (.25, .125, .0625):
            x, v = 1., 0.
            max_energy_error, max_residual = 0., 0.
            count = round(1. / h)
            for step in range(count):
                x1, v1 = np.linalg.solve(np.array([[1., -h / 2], [2 * h, 1.]]),
                                        np.array([x + h * v / 2, v - 2 * h * x]))
                a = fixture()
                a.update(dt=h, x0=np.array([[x, 0., 0.]]), x1=np.array([[x1, 0., 0.]]),
                         v0=np.array([[v, 0., 0.]]), v1=np.array([[v1, 0., 0.]]),
                         internal_force=np.array([[-2 * (x + x1), 0., 0.]]))
                r = midpoint_residual(**a)
                max_residual = max(max_residual, float(np.max(np.abs(r['dynamic_residual_N']))),
                                   float(np.max(np.abs(r['kinematic_residual_m']))))
                max_energy_error = max(max_energy_error, abs(.5 * v1 ** 2 + 2 * x1 ** 2 - 2.))
                x, v = x1, v1
            phase = count * 2 * np.arctan(h)
            self.assertAlmostEqual(x, np.cos(phase), delta=1e-12)
            self.assertAlmostEqual(v, -2 * np.sin(phase), delta=1e-12)
            self.assertLess(max_energy_error, 1e-12)
            self.assertLess(max_residual, 1e-12)
            reports.append(dict(dt_s=h, steps=count, final_x_m=x, final_v_m_s=v,
                                continuous_x_error_m=abs(x - np.cos(2.)),
                                continuous_v_error_m_s=abs(v + 2 * np.sin(2.)),
                                phase_error_rad=abs(phase - 2.),
                                max_energy_error_J=max_energy_error,
                                max_residual_SI=max_residual))
        for coarse, fine in zip(reports, reports[1:]):
            for key in ('continuous_x_error_m', 'continuous_v_error_m_s', 'phase_error_rad'):
                self.assertLess(fine[key], coarse[key])
                self.assertLess(fine[key] / coarse[key], .27)
        print('OSCILLATOR_REFERENCE ' + json.dumps(dict(status='INDEPENDENT_LINEAR_REFERENCE_ONLY',
              continuous_x_m=float(np.cos(2.)), continuous_v_m_s=float(-2 * np.sin(2.)),
              rows=reports)), flush=True)

    def test_23_unresolved_slip_event_endpoint_work_is_not_heat(self):
        # Independent piecewise Coulomb reference: a=-3 until t=1/3,
        # then a=-1 for remaining 2/3. Never call a contact routine here.
        stop = 1 / 3
        forward_dx = .5 * 1. * stop
        reverse_dx = -.5 * (2 / 3) * (1 - stop)
        true_dx = forward_dx + reverse_dx
        q_true = forward_dx - reverse_dx
        j_t = -stop + (1 - stop)
        u_final = -2 / 3
        u_endpoint_mid = .5 * (1 + u_final)
        q_endpoint_surrogate = -j_t * u_endpoint_mid
        self.assertAlmostEqual(q_true, 7 / 18, delta=1e-12)
        self.assertAlmostEqual(j_t, 1 / 3, delta=1e-12)
        self.assertAlmostEqual(true_dx, -1 / 18, delta=1e-12)
        self.assertAlmostEqual(q_endpoint_surrogate, -1 / 18, delta=1e-12)
        self.assertAlmostEqual(-2 * true_dx, 1 / 9, delta=1e-12)
        self.assertAlmostEqual(.5 * (u_final ** 2 - 1), -5 / 18, delta=1e-12)
        a = fixture()
        a.update(dt=1., v0=np.array([[1., 0., 0.]]), v1=np.array([[u_final, 0., 0.]]),
                 x1=np.array([[true_dx, 0., 0.]]), external_force=np.array([[-2., 0., 0.]]),
                 contact_force=np.array([[j_t, 0., 0.]]))
        r = midpoint_residual(**a)
        self.close(r['dynamic_residual_N'], 0.)
        self.close(r['kinematic_residual_m'], [[-2 / 9, 0., 0.]])
        self.assertLess(q_endpoint_surrogate, 0.)
        self.assertGreater(q_true, 0.)
        print('SLIP_EVENT_LIMITATION ' + json.dumps(dict(status='ENDPOINT_SURROGATE_IS_NOT_HEAT',
              Q_true_J=q_true, Q_endpoint_surrogate_J=q_endpoint_surrogate,
              displacement_true_m=true_dx, displacement_endpoint_mid_m=u_endpoint_mid,
              total_tangential_impulse_Ns=j_t)), flush=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)

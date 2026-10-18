"""Independent fixed-branch FKS/J2 numerical-response checks, no trajectory.

References use literal affine algebra, scalar FKS arithmetic and an independent
100-bisection return. Preregistration: MATERIAL_RESPONSE_TEST_CONTRACT.md.
"""
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import material_response as tested
from afsd_backend.pysph.material import FKS, update as original_update

EPSILONS = (1e-5, 5e-6)
STRESS_ATOL, DERIVATIVE_RTOL = .05, 2e-5
PLASTIC_ATOL = 2e-12
AXIAL = np.diag([-1/3, -1/3, 2/3])
SHEAR = 26.3e9


def emit(name, **data):
    print('MATERIAL_RESPONSE_EVIDENCE '+json.dumps(dict(name=name,
          classification='MATERIAL_POINT_COMPONENT_NOT_AFSD', **data),
          sort_keys=True, allow_nan=False), flush=True)


def literal_yield(law, ep, rate, temperature):
    """Scalar arithmetic only; no production yield/partial/tangent calls."""
    theta_star = (temperature-law.room)/(law.melt-law.room)
    hardening = law.a1+law.a2*math.atan(law.a3*ep)
    multiplier = 1+law.b1*theta_star**law.b2*math.log(
        max(rate, law.rate_floor)/law.reference_rate)
    # Stable scalar form of 1-(1+exp(-c1*Tstar))**(-1/c2).
    thermal = -math.expm1(-math.log1p(math.exp(-law.c1*theta_star))/law.c2)
    return hardening*multiplier*thermal


def literal_return(trial, ep0, rate, temperature, shear, law):
    qtrial = math.sqrt(1.5*float(np.sum(trial*trial)))
    lo, hi = 0., qtrial/(3*shear)
    if qtrial <= literal_yield(law, ep0, rate, temperature):
        increment = 0.
    else:
        for _ in range(100):
            mid = (lo+hi)/2
            residual = qtrial-3*shear*mid-literal_yield(
                law, ep0+mid, rate, temperature)
            if residual > 0:
                lo = mid
            else:
                hi = mid
        increment = (lo+hi)/2
    equivalent = qtrial-3*shear*increment
    return trial*(equivalent/qtrial), increment, equivalent


def affine_partials(shape=()):
    return {key: np.full(shape, value) for key, value in (
        ('yield_stress_Pa', 9.275), ('d_ep_Pa', 3.),
        ('d_rate_Pa_s', 2.), ('d_temperature_Pa_K', -.25))}


def dev_symmetric(matrix):
    symmetric = .5*(matrix+matrix.T)
    return symmetric-np.trace(symmetric)*np.eye(3)/3


def spectral_rotation(spin_dt):
    """Independent exponential, not Rodrigues or production rotation."""
    eigenvalues, eigenvectors = np.linalg.eigh(1j*spin_dt)
    result = (eigenvectors*np.exp(-1j*eigenvalues)) @ eigenvectors.conj().T
    if np.max(np.abs(result.imag)) > 2e-14:
        raise AssertionError('independent antisymmetric exponential lost reality')
    return result.real


def literal_trial(stress, gradient, dt, shear):
    deviator = dev_symmetric(gradient)
    spin = .5*(gradient-gradient.T)
    half = spectral_rotation(spin*dt/2)
    full = spectral_rotation(spin*dt)
    result = full@stress@full.T+2*shear*dt*(half@deviator@half.T)
    rate = math.sqrt(2/3*float(np.sum(deviator*deviator)))
    return result, rate


class MaterialResponseTests(unittest.TestCase):
    def assert_exact(self, observed, expected):
        np.testing.assert_allclose(observed, expected, atol=2e-12, rtol=2e-12)

    def test_affine_three_orthogonal_implicit_directions(self):
        trial = 14*AXIAL
        # Each expected scalar is a literal closed-form affine return value.
        cases = ((AXIAL, 0., 0., 1/24, 1/8),
                 (np.zeros((3, 3)), 1., 0., -1/12, 7/4),
                 (np.zeros((3, 3)), 0., 1., 1/96, -7/32))
        for dtrial, dr, dT, expected_dep, expected_q in cases:
            with self.subTest(dr=dr, dT=dT):
                result = tested.radial_direction(trial, .225, dtrial, dr, dT,
                                                 7., affine_partials(), True)
                self.assert_exact(result['delta_plastic_strain'], expected_dep)
                self.assert_exact(result['delta_equivalent_Pa'], expected_q)
                self.assert_exact(result['delta_stress_Pa'], expected_q*AXIAL)
                self.assert_exact(result['base_return_residual_Pa'], 0.)
                self.assertFalse(result['physical_update'])
        emit('affine_three_directions', plastic_increment=.225,
             returned_equivalent_Pa=9.275)

    def test_affine_nonproportional_transverse_tensor(self):
        direction = np.array([[0., .2, 0.], [.2, 0., .1], [0., .1, 0.]])
        result = tested.radial_direction(14*AXIAL, .225, direction, 0., 0.,
                                         7., affine_partials(), True)
        self.assert_exact(result['delta_stress_Pa'], (9.275/14)*direction)
        self.assert_exact(result['delta_plastic_strain'], 0.)
        self.assert_exact(result['delta_equivalent_Pa'], 0.)

    def test_batch_broadcast_elastic_and_input_purity(self):
        trial = np.stack((14*AXIAL, 6*AXIAL))[None]
        increments = np.array([[.225, 0.]])
        active = np.array([[True, False]])
        partials = affine_partials((1, 2))
        direction = np.array([[.1, .2, 0.], [.2, -.3, .1], [0., .1, .2]])
        saved = copy.deepcopy((trial, increments, active, partials, direction))
        result = tested.radial_direction(trial, increments, direction, .3, -.2,
                                         7., partials, active)
        self.assertEqual(result['delta_stress_Pa'].shape, (1, 2, 3, 3))
        for index in range(2):
            scalar = tested.radial_direction(trial[0, index], increments[0, index],
                direction, .3, -.2, 7., affine_partials(), bool(active[0, index]))
            for key in ('delta_stress_Pa', 'delta_plastic_strain', 'delta_equivalent_Pa'):
                self.assert_exact(result[key][0, index], scalar[key])
        self.assert_exact(result['delta_stress_Pa'][0, 1], direction)
        self.assert_exact(result['delta_plastic_strain'][0, 1], 0.)
        for current, before in zip((trial, increments, active, direction),
                                   (saved[0], saved[1], saved[2], saved[4])):
            np.testing.assert_array_equal(current, before)
            self.assertFalse(np.shares_memory(result['delta_stress_Pa'], current))
        for key in partials:
            np.testing.assert_array_equal(partials[key], saved[3][key])

    def test_literal_fks_values_and_three_partial_differences(self):
        law = FKS()
        cases = ((.01, 1.05, 723.15), (.2, 3., 400.), (.015, .4, 600.))
        values = np.array(cases)
        snapshot = values.copy()
        result = tested.fks_partials(law, values[:, 0], values[:, 1], values[:, 2])
        np.testing.assert_allclose(result['yield_stress_Pa'],
            [literal_yield(law, *case) for case in cases], atol=1e-6, rtol=2e-12)
        max_errors = dict(d_ep_Pa=0., d_rate_Pa_s=0., d_temperature_Pa_K=0.)
        for i, case in enumerate(cases):
            for axis, key, step in ((0, 'd_ep_Pa', 1e-7),
                    (1, 'd_rate_Pa_s', 1e-5), (2, 'd_temperature_Pa_K', 1e-3)):
                for epsilon in (step, step/2):
                    plus, minus = list(case), list(case)
                    plus[axis] += epsilon
                    minus[axis] -= epsilon
                    reference = (literal_yield(law, *plus)-literal_yield(law, *minus))/(2*epsilon)
                    np.testing.assert_allclose(result[key][i], reference, atol=.2, rtol=2e-6)
                    max_errors[key] = max(max_errors[key], abs(float(result[key][i])-reference))
        np.testing.assert_array_equal(values, snapshot)
        emit('literal_fks_partials', maximum_absolute_errors=max_errors)

    def test_floor_exact_rejection_and_domain_guards(self):
        law = FKS()
        result = tested.fks_partials(law, .01, np.array([0., .4, .999]), 723.15)
        np.testing.assert_array_equal(result['d_rate_Pa_s'], np.zeros(3))
        with self.assertRaises(tested.NonDifferentiableMaterial):
            tested.fks_partials(law, .01, 1., 723.15)
        for ep, rate, temperature in ((-.01, 1.1, 723.15), (.01, -.1, 723.15),
                (.01, 1.1, law.room-1), (.01, 1.1, law.melt),
                (.01, float('nan'), 723.15), (.01, True, 723.15)):
            with self.subTest(ep=ep, rate=rate, temperature=temperature):
                with self.assertRaises(ValueError):
                    tested.fks_partials(law, ep, rate, temperature)

    def test_radial_against_independent_fks_scalar_return(self):
        law = FKS()
        trial, ep0, rate, temperature = 5e6*AXIAL, .01, 1.05, 723.15
        _, increment, _ = literal_return(trial, ep0, rate, temperature, SHEAR, law)
        self.assertGreater(increment, 0.)
        partials = tested.fks_partials(law, ep0+increment, rate, temperature)
        mixed = 1e5*np.array([[.2, .4, -.1], [.4, -.3, .2], [-.1, .2, .1]])
        cases = ((1e5*AXIAL, 0., 0.), (np.zeros((3, 3)), .2, 0.),
                 (np.zeros((3, 3)), 0., 10.), (mixed, 0., 0.))
        errors = []
        for index, (dtrial, dr, dT) in enumerate(cases):
            actual = tested.radial_direction(trial, increment, dtrial, dr, dT,
                                              SHEAR, partials, True)
            for epsilon in EPSILONS:
                plus = literal_return(trial+epsilon*dtrial, ep0, rate+epsilon*dr,
                                      temperature+epsilon*dT, SHEAR, law)
                minus = literal_return(trial-epsilon*dtrial, ep0, rate-epsilon*dr,
                                       temperature-epsilon*dT, SHEAR, law)
                references = [(a-b)/(2*epsilon) for a, b in zip(plus, minus)]
                for key, expected, atol in zip(
                        ('delta_stress_Pa', 'delta_plastic_strain', 'delta_equivalent_Pa'),
                        references, (STRESS_ATOL, PLASTIC_ATOL, STRESS_ATOL)):
                    np.testing.assert_allclose(actual[key], expected, atol=atol,
                                               rtol=DERIVATIVE_RTOL)
                errors.append(dict(direction=index, epsilon=epsilon,
                    stress_max_Pa=float(np.max(np.abs(actual['delta_stress_Pa']-references[0]))),
                    plastic_abs=float(abs(actual['delta_plastic_strain']-references[1]))))
        emit('independent_fks_return', differences=errors)

    def test_original_48_return_same_history_composition(self):
        law, dt, ep0, temperature = FKS(), 5e-8, .01, 723.15
        stress = 5e6*AXIAL
        gradient = np.diag([-.525, -.525, 1.05])
        saved_stress, saved_gradient = stress.copy(), gradient.copy()

        def physical(g, temp):
            # Each call starts from identical base history; never accumulate.
            result = original_update(stress[None].copy(), np.array([ep0]),
                g[None].copy(), np.array([temp]), dt, SHEAR, law)
            return result[0][0], float(result[1][0]), float(result[4][0])

        base_stress, base_ep, rate = physical(gradient, temperature)
        trial = stress+2*SHEAR*dt*gradient
        reference, reference_dep, _ = literal_return(trial, ep0, rate, temperature, SHEAR, law)
        np.testing.assert_allclose(base_stress, reference, atol=1e-6, rtol=2e-12)
        np.testing.assert_allclose(base_ep, ep0+reference_dep, atol=2e-15, rtol=0.)
        partials = tested.fks_partials(law, base_ep, rate, temperature)
        dgradient = np.array([[-.15, .12, 0.], [.12, -.15, .07], [0., .07, .3]])
        errors = []
        for label, direction, dT in (('gradient', dgradient, 0.),
                                     ('temperature', np.zeros((3, 3)), 10.)):
            dr = (2/(3*rate))*float(np.sum(gradient*direction))
            actual = tested.radial_direction(trial, base_ep-ep0,
                2*SHEAR*dt*direction, dr, dT, SHEAR, partials, True)
            for epsilon in EPSILONS:
                plus = physical(gradient+epsilon*direction, temperature+epsilon*dT)
                minus = physical(gradient-epsilon*direction, temperature-epsilon*dT)
                ds, dep = (plus[0]-minus[0])/(2*epsilon), (plus[1]-minus[1])/(2*epsilon)
                np.testing.assert_allclose(actual['delta_stress_Pa'], ds,
                    atol=STRESS_ATOL, rtol=DERIVATIVE_RTOL)
                np.testing.assert_allclose(actual['delta_plastic_strain'], dep,
                    atol=PLASTIC_ATOL, rtol=DERIVATIVE_RTOL)
                errors.append(dict(direction=label, epsilon=epsilon,
                    stress_max_Pa=float(np.max(np.abs(actual['delta_stress_Pa']-ds))),
                    plastic_abs=float(abs(actual['delta_plastic_strain']-dep))))
        np.testing.assert_array_equal(stress, saved_stress)
        np.testing.assert_array_equal(gradient, saved_gradient)
        emit('unchanged_48_bisection_same_history', differences=errors)

    def test_radial_branch_and_norm_guards(self):
        baseline = [14*AXIAL, .225, AXIAL, 0., 0., 7., affine_partials(), True]
        changes = ((1, 0.), (1, -1.), (5, 0.), (7, 1),
                   (1, 1.), (0, np.full((3, 3), np.nan)),
                   (2, np.eye(3)), (6, {}))
        for index, value in changes:
            arguments = copy.deepcopy(baseline)
            arguments[index] = value
            with self.subTest(index=index, value=str(value)):
                with self.assertRaises(ValueError):
                    tested.radial_direction(*arguments)
        for trial in (np.zeros((3, 3)), np.diag([-3., -3., 6.])):
            boundary_partials = affine_partials()
            boundary_partials['yield_stress_Pa'] = 9.
            with self.assertRaises(tested.NonDifferentiableMaterial):
                tested.radial_direction(trial, 0., AXIAL, 0., 0.,
                                         7., boundary_partials, False)
        with self.assertRaises(ValueError):
            tested.radial_direction(14*AXIAL, 0., AXIAL, 0., 0.,
                                     7., affine_partials(), False)
        with self.assertRaises(ValueError):
            tested.radial_direction(6*AXIAL, .1, AXIAL, 0., 0.,
                                     7., affine_partials(), False)

    def test_trial_pure_rotation_of_nonzero_prestress(self):
        stress = np.diag([3., -1., -2.])
        generator = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 0.]])
        dt, omega = .2, .7
        angle = omega*dt
        rotation = np.array([[math.cos(angle), -math.sin(angle), 0.],
                             [math.sin(angle), math.cos(angle), 0.], [0., 0., 1.]])
        expected = rotation@stress@rotation.T
        derivative = dt*(generator@expected-expected@generator)
        actual = tested.trial_direction(stress, omega*generator, generator, dt, 7.)
        self.assert_exact(actual['trial_stress_Pa'], expected)
        self.assert_exact(actual['delta_trial_stress_Pa'], derivative)
        self.assert_exact(actual['rate_per_s'], 0.)
        self.assert_exact(actual['delta_rate_per_s'], 0.)

    def test_trial_noncommuting_spin_strain_and_batch(self):
        stress = np.array([[7., 2., -1.], [2., -4., .5], [-1., .5, -3.]])
        gradient = np.array([[.3, -.7, .2], [.4, -.1, -.3], [.1, .6, -.2]])
        direction = np.array([[-.2, .3, .4], [-.6, .35, .2], [.15, -.25, -.15]])
        saved = [a.copy() for a in (stress, gradient, direction)]
        actual = tested.trial_direction(stress, gradient, direction, .15, 7.)
        reference, rate = literal_trial(stress, gradient, .15, 7.)
        self.assert_exact(actual['trial_stress_Pa'], reference)
        self.assert_exact(actual['rate_per_s'], rate)
        errors = []
        for epsilon in EPSILONS:
            plus, rate_plus = literal_trial(stress, gradient+epsilon*direction, .15, 7.)
            minus, rate_minus = literal_trial(stress, gradient-epsilon*direction, .15, 7.)
            expected, dr = (plus-minus)/(2*epsilon), (rate_plus-rate_minus)/(2*epsilon)
            np.testing.assert_allclose(actual['delta_trial_stress_Pa'], expected,
                                       atol=2e-7, rtol=2e-6)
            np.testing.assert_allclose(actual['delta_rate_per_s'], dr, atol=2e-9, rtol=2e-6)
            errors.append(float(np.max(np.abs(actual['delta_trial_stress_Pa']-expected))))
        batch = tested.trial_direction(stress, np.stack((gradient, gradient)),
                                       direction, .15, 7.)
        for key in ('trial_stress_Pa', 'delta_trial_stress_Pa', 'rate_per_s', 'delta_rate_per_s'):
            self.assert_exact(batch[key][0], actual[key])
            self.assert_exact(batch[key][1], actual[key])
        for array, before in zip((stress, gradient, direction), saved):
            np.testing.assert_array_equal(array, before)
            self.assertFalse(np.shares_memory(actual['trial_stress_Pa'], array))
            self.assertFalse(np.shares_memory(actual['delta_trial_stress_Pa'], array))
        emit('spectral_trial_direction', maximum_stress_derivative_errors_Pa=errors)

    def test_trial_zero_rate_and_invalid_domain(self):
        stress = np.diag([3., -1., -2.])
        with self.assertRaises(tested.NonDifferentiableMaterial):
            tested.trial_direction(stress, np.zeros((3, 3)), AXIAL, .1, 7.)
        for dt, shear in ((0., 7.), (.1, 0.), (True, 7.), (.1, float('inf'))):
            with self.assertRaises(ValueError):
                tested.trial_direction(stress, AXIAL, AXIAL, dt, shear)


if __name__ == '__main__':
    emit('source_receipt', sha256={str(path.relative_to(ROOT)):
        hashlib.sha256(path.read_bytes()).hexdigest() for path in (
            Path(__file__), ROOT/'src/afsd_backend/pysph/material_response.py',
            ROOT/'src/afsd_backend/pysph/material.py',
            ROOT/'docs/MATERIAL_RESPONSE_TEST_CONTRACT.md')})
    unittest.main(verbosity=2)

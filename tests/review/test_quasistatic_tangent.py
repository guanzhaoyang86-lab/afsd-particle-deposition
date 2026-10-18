"""Frozen QS search derivatives; component evidence, not AFSD transport."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.material import FKS, update
from afsd_backend.pysph.material_response import NonDifferentiableMaterial
from afsd_backend.pysph.mechanics import corrected_operators, internal_force, pressure
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.quasistatic_tangent import (
    gradient_map, material_tangent, assemble_frozen_force_tangent)

G = 26.3e9
K = 2*G*(1+.33)/(3*(1-2*.33))
AXIAL = np.diag([-1/3, -1/3, 2/3])


class MaterialDerivativeTests(unittest.TestCase):
    def test_zero_elastic_state_has_literal_tensor_derivative(self):
        dt = .001
        response = material_tangent(np.zeros((2, 3, 3)), np.zeros(2),
            np.zeros((2, 3, 3)), 723.15, dt, G, FKS())
        direction = np.array([[2., 3., -4.], [-7., 1., 8.], [2., 4., -3.]])
        dev = .5*(direction+direction.T)-np.trace(direction)*np.eye(3)/3
        got = (response['tangent_Pa_s']@direction.ravel()).reshape(2, 3, 3)
        np.testing.assert_allclose(got, np.broadcast_to(2*G*dt*dev, got.shape),
                                   rtol=2e-14, atol=1e-7)
        np.testing.assert_array_equal(response['temperature_tangent_Pa_K'], 0.)
        self.assertFalse(response['plastic_active'].any())
        self.assertFalse(response['full_tangent'])

    def _direction_check(self, S0, L, dt, direction, temperature_direction=0.):
        S0, L, direction = (np.asarray(a)[None] for a in (S0, L, direction))
        ep = np.array([.01])
        T = np.array([723.15])
        law = FKS()
        physical = update(S0, ep, L, T, dt, G, law)
        response = material_tangent(S0, ep, L, T, dt, G, law, material=physical)
        derivative = (response['tangent_Pa_s'][0]@direction.ravel()
                      +response['temperature_tangent_Pa_K'][0]*temperature_direction).reshape(3, 3)
        for epsilon in (1e-5, 5e-6):
            plus = update(S0, ep, L+epsilon*direction,
                          T+epsilon*temperature_direction, dt, G, law)[0][0]
            minus = update(S0, ep, L-epsilon*direction,
                           T-epsilon*temperature_direction, dt, G, law)[0][0]
            fd = (plus-minus)/(2*epsilon)
            np.testing.assert_allclose(derivative, fd, rtol=2e-5, atol=.1)
        return response

    def test_zero_state_directional_difference(self):
        self._direction_check(np.zeros((3, 3)), np.zeros((3, 3)), .001,
                              .01*AXIAL)

    def test_rate_sensitive_plastic_combined_direction(self):
        direction = np.array([[2., -3., 4.], [-1., -2., 1.], [3., -2., 1.]])
        response = self._direction_check(3e6*AXIAL,
            np.array([[-20., 12., 3.], [-8., -20., 1.], [0., 2., 40.]]),
            1e-4, direction, 5.)
        self.assertTrue(response['plastic_active'][0])
        self.assertGreater(response['rate_per_s'][0], 1.)

    def test_zero_rate_plastic_composite_derivative(self):
        response = self._direction_check(10e6*AXIAL, np.zeros((3, 3)),
            .001, np.array([[.01, .02, 0.], [-.03, -.01, 0.], [0., 0., 0.]]), 3.)
        self.assertTrue(response['plastic_active'][0])
        self.assertEqual(response['rate_per_s'][0], 0.)

    def test_elastic_rate_floor_norm_switch_does_not_change_stress(self):
        L = 1.5*AXIAL  # sqrt(2/3 D:D) == 1 exactly.
        response = material_tangent(np.zeros((1, 3, 3)), np.zeros(1),
            L[None], 723.15, 1e-9, G, FKS())
        self.assertFalse(response['plastic_active'][0])
        self.assertEqual(response['rate_per_s'][0], 1.)

    def test_plastic_rate_floor_switch_rejects(self):
        with self.assertRaises(NonDifferentiableMaterial):
            material_tangent((10e6*AXIAL)[None], np.zeros(1),
                (1.5*AXIAL)[None], 723.15, 1e-4, G, FKS())

    def test_wrong_material_history_rejects(self):
        S = np.zeros((1, 3, 3)); ep = np.zeros(1)
        returned = list(update(S, ep, S, np.array([723.15]), .001, G, FKS()))
        returned[1] = np.array([.1])
        with self.assertRaisesRegex(ValueError, 'branch'):
            material_tangent(S, ep, S, 723.15, .001, G, FKS(), returned)


class SparseForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        x = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                                 [-.0005, 0., .0005], indexing='ij')).reshape(3, -1).T
        p = dict(dx=.0005, h_ratio=1.3, rho0=2700., width=.0015,
                 marker_low=-.001, marker_high=.001, bar_length=.01,
                 actuator_depth=.0005, temperature=723.15, traverse=.00212, feed=.00212)
        cls.pa = make_array(x, np.zeros(27, int), p)
        cls.x = x
        cls.i, cls.j, cls.g, _, _, cls.B = corrected_operators(cls.pa)
        cls.V = cls.pa.m/cls.pa.rho
        cls.D = gradient_map(cls.V, cls.i, cls.j, cls.g, cls.B)

    def test_gradient_and_force_match_actual_pair_operators(self):
        velocity = self.x@np.array([[2., 3., 4.], [1., 2., -1.], [-3., 1., 7.]]).T
        private = self.pa.extract_particles(np.arange(27, dtype=np.uint32))
        private.u[:], private.v[:], private.w[:] = velocity.T
        actual = corrected_operators(private)[3]
        np.testing.assert_allclose((self.D@velocity.ravel()).reshape(27, 3, 3), actual,
                                   rtol=2e-13, atol=1e-12)
        stress = np.random.default_rng(41).normal(size=(27, 3, 3))*1e6
        stress = .5*(stress+stress.swapaxes(1, 2))
        expected = internal_force(self.pa, stress, self.i, self.j, self.g, self.B)
        actual_force = (-self.D.T@(self.V[:, None, None]*stress).ravel()).reshape(27, 3)
        np.testing.assert_allclose(actual_force, expected, rtol=2e-13, atol=1e-12)

    def test_frozen_endpoint_force_directional_derivative(self):
        n = 27; dt = 1e-4
        u = self.x@np.array([[.0008, .0006, 0.], [0., -.0002, .0003], [0., 0., .0004]]).T
        direction = np.random.default_rng(99).normal(size=(n, 3))*1e-7
        S0 = np.broadcast_to(3e6*AXIAL, (n, 3, 3)).copy()
        ep0 = np.full(n, .01)
        L = (self.D@u.ravel()).reshape(n, 3, 3)/dt
        response = material_tangent(S0, ep0, L, 723.15, dt, G, FKS())
        J = assemble_frozen_force_tangent(self.D, self.V, response['tangent_Pa_s'],
                                         self.pa.rho, 2700., K, dt, u)
        def force(v):
            strain = (self.D@v.ravel()).reshape(n, 3, 3)
            rho = self.pa.rho/(1+np.trace(strain, axis1=1, axis2=2))
            stress = update(S0, ep0, strain/dt, np.full(n, 723.15), dt, G, FKS())[0]
            sigma = stress-pressure(rho, 2700., K)[:, None, None]*np.eye(3)
            return -self.D.T@(self.V[:, None, None]*sigma).ravel()
        expected = J@direction.ravel()
        for epsilon in (1e-4, 5e-5):
            actual = (force(u+epsilon*direction)-force(u-epsilon*direction))/(2*epsilon)
            np.testing.assert_allclose(expected, actual, rtol=3e-5, atol=1e-7)
        self.assertTrue(np.isfinite(J).all())

    def test_elastic_zero_force_block_has_no_inertia(self):
        n = 27
        matrices = []
        for dt in (1e-5, 1e-3):
            response = material_tangent(np.zeros((n, 3, 3)), np.zeros(n),
                np.zeros((n, 3, 3)), 723.15, dt, G, FKS())
            matrices.append(assemble_frozen_force_tangent(self.D, self.V,
                response['tangent_Pa_s'], self.pa.rho, 2700., K, dt,
                np.zeros((n, 3))))
        np.testing.assert_allclose(matrices[0], matrices[1], rtol=2e-14, atol=1e-7)
        rigid_translation = np.tile([1., 2., -3.], n)
        scale = np.linalg.norm(matrices[0])*np.linalg.norm(rigid_translation)
        self.assertLess(np.linalg.norm(matrices[0]@rigid_translation), 1e-13*scale)


if __name__ == '__main__':
    unittest.main(verbosity=2)

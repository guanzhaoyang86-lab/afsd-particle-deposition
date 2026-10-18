"""Independent manufactured references for the isolated fixed-cube prototype."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from probe_boundary_sbp import evaluate, fixture, extract_original, gradient, adjoint_force
from afsd_backend.pysph.mechanics import operators


class FixedCubeSBP(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report, cls.state = evaluate()
        cls.x = cls.state['position_m']
        cls.V = cls.state['volume_m3']
        cls.D = cls.state['prototype_D_per_m']
        cls.e = cls.state['signed_face_area_m2']

    def test_independent_tensor_difference_reference(self):
        # Written from endpoint/central stencil formulas, NOT a tested helper.
        dx = .5
        d = np.zeros((5, 5))
        d[0, :2] = [-1/dx, 1/dx]
        d[-1, -2:] = [-1/dx, 1/dx]
        for row in range(1, 4):
            d[row, row-1], d[row, row+1] = -.5/dx, .5/dx
        eye = np.eye(5)
        reference = np.array([np.kron(np.kron(d, eye), eye),
                              np.kron(np.kron(eye, d), eye),
                              np.kron(np.kron(eye, eye), d)])
        for k in range(3):
            q = self.V[:, None]*reference[k]
            np.testing.assert_allclose(q+q.T, np.diag(self.e[k]), atol=1e-12, rtol=0.)
            np.testing.assert_allclose(reference[k]@self.x,
                np.tile(np.eye(3)[k], (125, 1)), atol=1e-12, rtol=0.)
        # Reference demonstrates existence on this grid; NOT substituted for SPH.
        self.assertGreater(np.linalg.norm(self.D-reference), 1e-8)

    def test_independent_geometry_and_actual_source_bridge(self):
        self.assertAlmostEqual(self.V.sum(), 8., delta=1e-12)
        for k in range(3):
            for sign in (-1, 1):
                mask = self.x[:, k] == sign
                self.assertAlmostEqual(np.sum(self.e[k, mask]), sign*4., delta=1e-12)
            np.testing.assert_array_equal(self.e[k, np.abs(self.x[:, k]) < 1], 0.)
        for error in self.report['bridge'].values():
            self.assertLessEqual(error, 1e-10)
        self.assertGreater(self.report['operators']['original']['constant_patch_relative'], 1e-8)

    def test_projection_compatibility_sbp_and_original_support(self):
        for name in ('compatibility_max', 'linear_constraint_max'):
            self.assertLessEqual(self.report['projection'][name], 1e-10)
        support = np.eye(125, dtype=bool)
        i, j = self.state['neighbor_i'], self.state['neighbor_j']
        support[i, j] = support[j, i] = True
        np.testing.assert_array_equal(self.D[:, ~support], 0.)
        for k in range(3):
            q = self.V[:, None]*self.D[k]
            np.testing.assert_allclose(q+q.T, np.diag(self.e[k]), atol=1e-10, rtol=0.)

    def test_full_three_dimensional_affine_reproduction(self):
        A = np.array([[.3, -.7, .2], [.5, -.1, .8], [-.4, .9, .6]])
        velocity = self.x@A.T + np.array([.2, -.3, .7])
        expected = np.broadcast_to(A, (125, 3, 3))
        np.testing.assert_allclose(gradient(self.D, velocity), expected, atol=1e-10, rtol=0.)

    def test_independent_constant_stress_patch_and_virtual_power(self):
        sigma = np.array([[2., .4, -.3], [.4, 3., .2], [-.3, .2, 1.]])
        stress = np.broadcast_to(sigma, (125, 3, 3))
        # Compute geometric face loads directly, without reading saved Fext.
        ext = np.zeros((125, 3))
        for node, point in enumerate(self.x):
            for axis in range(3):
                if abs(point[axis]) == 1.:
                    area = 1.
                    for tangent in range(3):
                        if tangent != axis:
                            area *= .25 if abs(point[tangent]) == 1. else .5
                    ext[node] += area*point[axis]*sigma[:, axis]
        fi = adjoint_force(self.D, self.V, stress)
        self.assertLessEqual(np.linalg.norm(fi+ext)/np.linalg.norm(ext), 1e-10)
        np.testing.assert_allclose(ext, self.state['external_constant_stress_traction_N'], atol=1e-12)
        x, y, z = self.x.T
        bubble = np.column_stack((x*(1-x*x)*(1-y*y)*(1-z*z), np.zeros(125), np.zeros(125)))
        for v, exact in ((self.x, 48.), (bubble, 0.)):
            work = np.einsum('i,iak,iak', self.V, stress, gradient(self.D, v))
            self.assertAlmostEqual(work, exact, delta=1e-10)
            self.assertAlmostEqual(np.sum(ext*v), exact, delta=1e-10)
        # A different constant symmetric stress ensures it was not fitted to sigma.
        other = np.broadcast_to(np.array([[1., -.6, .8], [-.6, -2., .3], [.8, .3, .5]]),
                                stress.shape)
        force = adjoint_force(self.D, self.V, other)
        loads = np.einsum('ki,iak->ia', self.e, other)
        self.assertLessEqual(np.linalg.norm(force+loads)/np.linalg.norm(loads), 1e-10)

    def test_arbitrary_nodal_ibp_power_total_force_and_torque(self):
        rng = np.random.default_rng(4007)  # Independent holdout from probe seed.
        f, g = rng.normal(size=(2, 125))
        for d, e in zip(self.D, self.e):
            self.assertAlmostEqual(f@(self.V*(d@g))+(d@f)@(self.V*g), f@(e*g), delta=1e-10)
        velocity = rng.normal(size=(125, 3))
        stress = rng.normal(size=(125, 3, 3))
        stress += stress.transpose(0, 2, 1)
        force = adjoint_force(self.D, self.V, stress)
        power = np.einsum('i,iak,iak', self.V, stress, gradient(self.D, velocity))
        self.assertAlmostEqual(power, -np.sum(velocity*force), delta=1e-10)
        np.testing.assert_allclose(force.sum(axis=0), 0., atol=1e-10, rtol=0.)
        np.testing.assert_allclose(np.cross(self.x, force).sum(axis=0), 0., atol=1e-10, rtol=0.)

    def test_rigid_virtual_rotation_nonzero_stress(self):
        omega = np.array([.3, -.2, .5])
        velocity = np.cross(omega, self.x)
        reference = np.array([[0., -.5, -.2], [.5, 0., -.3], [.2, .3, 0.]])
        L = gradient(self.D, velocity)
        np.testing.assert_allclose(L, np.broadcast_to(reference, L.shape), atol=1e-10, rtol=0.)
        stress = np.broadcast_to(np.diag([1., 2., 3.]), L.shape)
        self.assertAlmostEqual(np.einsum('i,iak,iak', self.V, stress, L), 0., delta=1e-10)
        # This is a virtual-work test, not a Jaumann time-integration test.

    def test_no_additional_fixed_grid_zero_modes(self):
        scalar = self.state['prototype_scalar_gradient_singular_values']
        strain = self.state['prototype_symmetric_gradient_singular_values']
        self.assertEqual(np.sum(scalar <= 1e-10*scalar[0]), 1)
        self.assertEqual(np.sum(strain <= 1e-10*strain[0]), 6)

    def test_degenerate_moment_still_fails_explicitly(self):
        for kind in ('planar_distinct', 'thin_columns', 'thin_distinct_xy'):
            with self.subTest(kind=kind):
                pa, _, _, _ = fixture()
                if kind != 'thin_columns':
                    # Keep particles distinct so the moment guard is reached.
                    pa.x[:] += .07*pa.z
                if kind == 'planar_distinct':
                    pa.z[:] = 0.
                else:
                    pa.z[:] *= 1e-4
                # Thickness alone is NOT the condition number for a cusp kernel.
                # Independently assemble M before expecting the guard to reject.
                i, j, g, _, _ = operators(pa)
                x = np.column_stack((pa.x, pa.y, pa.z))
                volume = pa.m/pa.rho
                moment = np.zeros((125, 3, 3))
                for a, b, grad in zip(i, j, g):
                    term = np.outer(x[b]-x[a], grad)
                    moment[a] += volume[b]*term
                    moment[b] += volume[a]*term
                eig = np.linalg.eigvalsh(moment)
                bad = np.any((eig[:, 0] <= 0) | (eig[:, -1] > 1e6*eig[:, 0]))
                self.assertEqual(bool(bad), kind != 'thin_columns')
                if bad:
                    with self.assertRaisesRegex(ValueError, 'correction condition limit=1e\\+06'):
                        extract_original(pa)
                else:
                    _, _, condition = extract_original(pa)
                    self.assertLessEqual(condition, 1e6)

    def test_coincident_particles_fail_before_moment_correction(self):
        pa, _, _, _ = fixture()
        pa.z[:] = 0.  # Rejection-only input; NOT a deformed-grid validation.
        with self.assertRaisesRegex(ValueError, 'distinct coincident particles'):
            extract_original(pa)


if __name__ == '__main__':
    unittest.main(verbosity=2)


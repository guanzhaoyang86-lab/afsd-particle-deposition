"""Independent references for a frozen elastic PRECONDITIONER, not AFSD."""
from pathlib import Path
import hashlib
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.structural_mobility import (
    CLASSIFICATION, ConstrainedResponse, build_frozen_elastic)
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.mechanics import corrected_operators, internal_force, positions

G = 26.3e9
BULK = 2*G*(1+.33)/(3*(1-2*.33))


def cloud(distorted=False):
    xyz = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                               [-.0005, 0., .0005], indexing='ij')).reshape(3, -1).T
    if distorted:
        rng = np.random.default_rng(20261003)
        xyz += rng.uniform(-.00003, .00003, xyz.shape)
    p = dict(dx=.0005, h_ratio=1.3, rho0=2700., width=.0015,
        marker_low=-.001, marker_high=.001, bar_length=.01,
        actuator_depth=.0005, temperature=400., traverse=.00212, feed=.00212)
    return make_array(xyz, np.zeros(27, int), p)


def pa_hash(pa):
    h = hashlib.sha256()
    for key in sorted(pa.properties):
        a = pa.properties[key].get_npy_array()
        h.update(key.encode()); h.update(a.tobytes())
    return h.hexdigest()


def legacy_elastic_force(pa, displacement):
    """Independent existing operator path; never uses new D, K or stress."""
    q = pa.extract_particles(np.arange(len(pa.x), dtype=np.uint32))
    q.u[:], q.v[:], q.w[:] = displacement.T
    i, j, g, l, _, b = corrected_operators(q)
    strain = .5*(l+l.swapaxes(1, 2))
    volumetric = np.trace(strain, axis1=1, axis2=2)
    dev = strain-volumetric[:, None, None]*np.eye(3)/3
    sigma = 2*G*dev+BULK*volumetric[:, None, None]*np.eye(3)
    return internal_force(q, sigma, i, j, g, correction=b), l


class MobilityAlgebraTests(unittest.TestCase):
    def spring(self, free=(True, True)):
        return ConstrainedResponse([2., 3.], [[4., -4.], [-4., 4.]], 1., free)

    def test_01_two_dof_literal_cross_contact_response(self):
        op = self.spring()
        expected = np.array([[4., 1.], [1., 3.]])/11.
        np.testing.assert_allclose(op.contact_response(np.eye(2)), expected, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(op.contact_response(np.eye(2), at='midpoint'), expected/2,
                                   rtol=2e-13, atol=2e-13)
        self.assertGreater(op.contact_response(np.eye(2))[0, 1], .09)
        got = op.response([1., 2.])
        np.testing.assert_allclose(got['delta_velocity_m_s'], [6/11, 7/11], rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(got['reaction_impulse_Ns'], [0., 0.], atol=2e-13)

    def test_02_prescribed_increment_literal_reaction(self):
        op = self.spring((True, False))
        got = op.response([1., 0.], [0., .2])
        np.testing.assert_allclose(got['delta_velocity_m_s'], [.4, .2], atol=2e-13, rtol=2e-13)
        np.testing.assert_allclose(got['reaction_impulse_Ns'], [0., .4], atol=2e-13, rtol=2e-13)
        np.testing.assert_allclose(op.contact_response(np.eye(2)), [[1/3, 0.], [0., 0.]], atol=2e-13)

    def test_03_mass_only_exact_not_epsilon_stiffness(self):
        op = ConstrainedResponse([2., 4., 5.], np.zeros((3, 3)), 10., [True, False, True])
        contact = np.array([[1., 2., 0.], [1., 0., -2.]])
        expected = np.array([[.5, .5], [.5, 1.3]])
        np.testing.assert_allclose(op.contact_response(contact), expected, atol=2e-13, rtol=2e-13)
        result = op.response([2., 3., -5.])
        np.testing.assert_allclose(result['delta_velocity_m_s'], [1., 0., -1.], atol=2e-13)
        np.testing.assert_allclose(result['reaction_impulse_Ns'], [0., -3., 0.], atol=2e-13)
        self.assertEqual(op.diagnostics['stiffness_near_zero_eigenvalues'], 3)

    def test_04_all_constrained_dofs_retained(self):
        op = self.spring((False, False))
        result = op.response([1., -2.], [.2, -.3])
        np.testing.assert_allclose(result['delta_velocity_m_s'], [.2, -.3], atol=0.)
        # A=[[3,-1],[-1,4]], so A*[.2,-.3]-[1,-2]=[-.1,.6].
        np.testing.assert_allclose(result['reaction_impulse_Ns'], [-.1, .6], atol=2e-13)
        np.testing.assert_array_equal(op.contact_response(np.eye(2)), np.zeros((2, 2)))
        self.assertIsNone(op.diagnostics['free_matrix_condition'])

    def test_05_arbitrary_constraint_rows_independent_balance(self):
        op = ConstrainedResponse([1., 2., 3.], [[3., -1., -1.], [-1., 3., -1.], [-1., -1., 3.]],
                                 .4, [False, True, False])
        got = op.response([1., 2., 3.], [.1, 0., -.2])
        # A21=A23=-.04 and A22=2.12: 2.12*dv2=2-.004=1.996.
        v = np.array([.1, 1.996/2.12, -.2])
        a = np.array([[1.12, -.04, -.04], [-.04, 2.12, -.04], [-.04, -.04, 3.12]])
        np.testing.assert_allclose(got['delta_velocity_m_s'], v, atol=2e-13)
        np.testing.assert_allclose(got['reaction_impulse_Ns'], a@v-[1., 2., 3.], atol=2e-13)

    def test_06_no_input_alias_and_invalid_requests_reject(self):
        mass = np.array([2., 3.]); k = np.array([[4., -4.], [-4., 4.]])
        op = ConstrainedResponse(mass, k, 1., [True, False])
        mass[:] = 100.; k[:] = 0.
        np.testing.assert_allclose(op.matrix, [[3., -1.], [-1., 4.]], atol=0.)
        with self.assertRaisesRegex(ValueError, 'free axis'): op.response([0., 0.], [1., 0.])
        with self.assertRaises(ValueError): op.contact_response(np.eye(2), at='guess')
        with self.assertRaises(ValueError): op.contact_response(np.eye(3))
        with self.assertRaises(ValueError): op.response([np.nan, 0.])
        for bad in ([0., 3.], [-1., 3.], [np.nan, 3.]):
            with self.assertRaises(ValueError): ConstrainedResponse(bad, np.eye(2), 1., [1, 1])
        with self.assertRaisesRegex(ValueError, 'nonsymmetric'):
            ConstrainedResponse([1., 1.], [[1., 1.], [0., 1.]], 1., [1, 1])
        with self.assertRaisesRegex(ValueError, 'negative elastic'):
            ConstrainedResponse([1., 1.], np.diag([-1., 1.]), 1., [1, 1])
        with self.assertRaises(ValueError): ConstrainedResponse([1., 1.], np.eye(2), 0., [1, 1])
        with self.assertRaises(ValueError): ConstrainedResponse([1., 1.], np.eye(2), 1., [.5, 1])


class ActualSPHComponentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pa = cloud(); cls.original = pa_hash(cls.pa)
        cls.op = build_frozen_elastic(cls.pa, G, BULK)
        cls.k_scale = np.linalg.norm(cls.op.stiffness, 2)

    def test_07_three_dimensional_affine_gradient_stress_energy_literal(self):
        a = np.array([[.001, .002, -.001], [-.003, -.0005, .002], [.001, -.002, .0002]])
        u = positions(self.pa)@a.T+np.array([1e-5, -2e-5, 3e-5])
        expected = np.broadcast_to(a, (27, 3, 3))
        np.testing.assert_allclose(self.op.gradient(u), expected, rtol=2e-11, atol=2e-12)
        strain = .5*(a+a.T); trace = np.trace(a); dev = strain-trace*np.eye(3)/3
        sigma = 2*G*dev+BULK*trace*np.eye(3)
        np.testing.assert_allclose(self.op.stress(u), np.broadcast_to(sigma, (27, 3, 3)),
                                   rtol=2e-11, atol=1e-5)
        reference = (G*np.sum(dev*dev)+.5*BULK*trace*trace)*(27*.0005**3)
        self.assertLessEqual(abs(self.op.energy(u)-reference), 1e-15+2e-11*abs(reference))

    def test_08_nonzero_rigid_translation_and_rotation_null_modes(self):
        x = positions(self.pa)
        spin = np.array([[0., -2., -3.], [2., 0., -1.], [3., 1., 0.]])
        for u in (np.broadcast_to([.1, -.2, .3], x.shape).copy(), x@spin.T):
            residual = np.linalg.norm(self.op.internal_force(u))
            self.assertLessEqual(residual, 2e-11*self.k_scale*np.linalg.norm(u))
        rotation_gradient = self.op.gradient(x@spin.T)
        np.testing.assert_allclose(rotation_gradient, np.broadcast_to(spin, (27, 3, 3)),
                                   rtol=2e-11, atol=2e-12)

    def test_09_actual_adjoint_force_energy_and_directional_derivative(self):
        rng = np.random.default_rng(203)
        u = rng.normal(size=(27, 3))*1e-6
        direction = rng.normal(size=(27, 3))*1e-6
        force, gradient = legacy_elastic_force(self.pa, u)
        error = np.linalg.norm(self.op.internal_force(u)-force)
        self.assertLessEqual(error, 1e-12+2e-11*np.linalg.norm(force))
        np.testing.assert_allclose(self.op.gradient(u), gradient, rtol=2e-11, atol=2e-12)
        reference_energy = -.5*float(np.sum(u*force))
        self.assertLessEqual(abs(self.op.energy(u)-reference_energy), 1e-15+2e-11*abs(reference_energy))
        epsilon = 1e-6
        plus, _ = legacy_elastic_force(self.pa, u+epsilon*direction)
        minus, _ = legacy_elastic_force(self.pa, u-epsilon*direction)
        finite_difference = -(plus-minus)/(2*epsilon)
        action = (self.op.stiffness@direction.ravel()).reshape(27, 3)
        self.assertLessEqual(np.linalg.norm(action-finite_difference), 1e-12+2e-8*np.linalg.norm(action))

    def test_10_free_surface_distorted_neighbors_and_nonuniform_density(self):
        pa = cloud(distorted=True)
        pa.rho[:] *= np.linspace(.99, 1.01, 27)
        before = pa_hash(pa); op = build_frozen_elastic(pa, G, BULK)
        a = np.array([[.001, -.0002, .0003], [.0004, .0005, -.0006], [.0007, .0008, -.0009]])
        u = positions(pa)@a.T
        np.testing.assert_allclose(op.gradient(u), np.broadcast_to(a, (27, 3, 3)), rtol=2e-11, atol=2e-12)
        force, _ = legacy_elastic_force(pa, u)
        self.assertLessEqual(np.linalg.norm(op.internal_force(u)-force), 1e-12+2e-11*np.linalg.norm(force))
        self.assertEqual(before, pa_hash(pa))

    def test_11_true_mass_inertia_positive_does_not_hide_stiffness_nullspace(self):
        response = self.op.response_operator(1e-7, np.ones((27, 3), bool))
        np.testing.assert_array_equal(response.mass_per_dof, np.repeat(self.pa.m, 3))
        expected = np.diag(np.repeat(self.pa.m, 3))+1e-14*self.op.stiffness/4
        np.testing.assert_allclose(response.matrix, expected, atol=1e-20, rtol=2e-13)
        self.assertGreater(np.linalg.eigvalsh(response.matrix).min(), 0.)
        zero = response.diagnostics['stiffness_near_zero_eigenvalues']
        self.assertGreaterEqual(zero, 6)
        print('ACTUAL_SPH_STIFFNESS_DIAGNOSTIC '+repr(response.diagnostics), flush=True)
        self.assertTrue(response.diagnostics['inertia_does_not_certify_stiffness_coercivity'])

    def test_12_no_particle_mutation_or_label_force_coupling(self):
        self.assertEqual(pa_hash(self.pa), self.original)
        changed = self.pa.extract_particles(np.arange(27, dtype=np.uint32))
        changed.tracer_type[:] = np.arange(27) % 3
        changed.wholebar_tag[:] = np.arange(27) % 4
        altered = build_frozen_elastic(changed, G, BULK)
        np.testing.assert_array_equal(altered.stiffness, self.op.stiffness)
        np.testing.assert_array_equal(altered.gradient_map, self.op.gradient_map)
        self.assertEqual(altered.classification, CLASSIFICATION)
        self.assertFalse(self.op.response_operator(1e-7, np.ones((27, 3), bool)).diagnostics['full_tangent'])

    def test_13_planar_neighborhood_rejected_no_pseudoinverse(self):
        subset = self.pa.extract_particles(np.flatnonzero(self.pa.z == 0.).astype(np.uint32))
        with self.assertRaisesRegex(ValueError, 'degenerate/ill-conditioned'):
            build_frozen_elastic(subset, G, BULK)

    def test_14_particle_permutation_preserves_gradient_and_stiffness_id_correspondence(self):
        order = np.random.default_rng(44).permutation(27).astype(np.uint32)
        q = self.pa.extract_particles(order)
        other = build_frozen_elastic(q, G, BULK)
        dofs = (3*order[:, None]+np.arange(3)).ravel()
        expected = self.op.stiffness[np.ix_(dofs, dofs)]
        np.testing.assert_allclose(other.stiffness, expected, rtol=2e-11, atol=1e-5)
        np.testing.assert_array_equal(other.particle_ids, self.op.particle_ids[order])


if __name__ == '__main__':
    unittest.main(verbosity=2)

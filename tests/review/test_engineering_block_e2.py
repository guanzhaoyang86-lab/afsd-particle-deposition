"""E2 algebra, provenance and short-step regressions; not AFSD acceptance."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from probe_engineering_block import load_config, make_problem, reference
from afsd_backend.pysph.mechanics import corrected_operators, internal_force, positions, velocities
from afsd_backend.pysph.output import checkpoint, restore
from pysph.solver.utils import load


class EngineeringBlockE2(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config()

    def test_independent_literal_continuum_end_values(self):
        # Manufactured PDE oracle, not a call to the SPH gradient/material law.
        ell, rho, young, alpha = .01, 2700., 70e9, 1e-5
        c = np.sqrt(young/rho)
        k = np.pi/(2*ell)
        x = np.array([[0., 0., 0.], [.005, .003, .007], [.01, .01, .01]])
        ref = reference(self.cfg, x, ell/c)
        displacement = np.zeros((3, 3))
        displacement[:, 0] = alpha/k*np.array([0., np.sqrt(.5), 1.])
        velocity = np.zeros((3, 3))
        velocity[:, 0] = alpha*c*np.array([0., np.sqrt(.5), 1.])
        stress = np.zeros((3, 3, 3))
        stress[:, 0, 0] = 700000.*np.array([1., np.sqrt(.5), 0.])
        np.testing.assert_allclose(ref['displacement_m'], displacement, rtol=1e-11, atol=1e-20)
        np.testing.assert_allclose(ref['velocity_m_s'], velocity, rtol=1e-11, atol=1e-14)
        np.testing.assert_allclose(ref['stress_Pa'], stress, rtol=1e-11, atol=1e-8)
        self.assertAlmostEqual(ref['reaction_N'], -70., delta=1e-10)
        self.assertAlmostEqual(ref['work_J'], 3.5e-6, delta=1e-16)
        self.assertAlmostEqual(ref['kinetic_J'], 1.75e-6, delta=1e-16)
        self.assertAlmostEqual(ref['elastic_J'], 1.75e-6, delta=1e-16)
        initial = reference(self.cfg, x, 0.)
        for field in ('displacement_m', 'velocity_m_s', 'stress_Pa'):
            np.testing.assert_array_equal(initial[field], np.zeros_like(initial[field]))
        self.assertIs(self.cfg['production_authorized'], False)

    def test_fixed_domain_quadrature_mass_and_load(self):
        for n in (3, 5):
            b = make_problem(self.cfg, n, .1)
            self.assertEqual(b.pa.get_number_of_particles(), n**3)
            np.testing.assert_array_equal(b.X.min(axis=0), [0., 0., 0.])
            np.testing.assert_allclose(b.X.max(axis=0), [.01, .01, .01], rtol=0, atol=1e-16)
            self.assertAlmostEqual(b.V.sum()/1e-6, 1., delta=1e-12)
            self.assertAlmostEqual(b.pa.m.sum()/.0027, 1., delta=1e-12)
            np.testing.assert_allclose(b.pa.m, 2700.*b.V, rtol=1e-12, atol=0)
            np.testing.assert_array_equal(b.fixed, b.X[:, 0] == 0.)
            self.assertEqual(np.count_nonzero(b.fixed), n*n)
            expected = b.pa.m*1e-5*(70e9/2700.)*(np.pi/.02)*np.sin(np.pi*b.X[:, 0]/.02)
            np.testing.assert_allclose(b.body[:, 0], expected, rtol=1e-11, atol=1e-14)
            np.testing.assert_array_equal(b.body[:, 1:], np.zeros((n**3, 2)))
            self.assertLessEqual(b.dt*b.omega_bound, 1.+1e-12)
            self.assertAlmostEqual(b.steps*b.dt, .01/np.sqrt(70e9/2700.), delta=1e-20)

    def test_all_three_dimensional_gradient_and_actual_interface(self):
        b = make_problem(self.cfg, 3, .1)
        affine = np.array([[.3, -.7, .4], [1.2, .1, -.6], [-.2, .9, .8]])
        field = b.X@affine.T+np.array([.02, -.01, .03])
        np.testing.assert_allclose(b.gradient(field), np.broadcast_to(affine, (27, 3, 3)), rtol=0, atol=1e-11)
        # The extracted dense derivative must retain the actual NNPS/corrected
        # operator, not only reproduce an affine special case.
        field = np.column_stack((b.X[:, 0]*b.X[:, 2], b.X[:, 1]**2, np.sin(17*b.X[:, 2])))
        b.ref.u[:], b.ref.v[:], b.ref.w[:] = field.T
        actual = corrected_operators(b.ref)[3]
        np.testing.assert_allclose(b.gradient(field), actual, rtol=1e-11, atol=1e-12)
        self.assertLessEqual(b.condition, 1e6)

    def test_actual_force_adjoint_power_force_and_torque(self):
        b = make_problem(self.cfg, 3, .1)
        x = b.X/.01
        sigma = np.zeros((27, 3, 3))
        sigma[:, 0, 0] = 2.+x[:, 0]
        sigma[:, 1, 1] = -1.+x[:, 1]
        sigma[:, 2, 2] = 3.-x[:, 2]
        sigma[:, 0, 1] = sigma[:, 1, 0] = .7+x[:, 2]
        sigma[:, 0, 2] = sigma[:, 2, 0] = -.4+x[:, 1]
        sigma[:, 1, 2] = sigma[:, 2, 1] = .2+x[:, 0]
        force = internal_force(b.ref, sigma, b.i, b.j, b.g, correction=b.B)
        independent_adjoint = -np.einsum('kij,i,iak->ja', b.D, b.V, sigma)
        np.testing.assert_allclose(force, independent_adjoint, rtol=1e-11, atol=1e-16)
        velocity = np.column_stack((x[:, 0]*x[:, 1], x[:, 2]**2, x[:, 0]-.3*x[:, 2]))
        kinetic_power = float(np.sum(force*velocity))
        internal_power = float(np.einsum('i,iak,kij,ja->', b.V, sigma, b.D, velocity))
        scale = float(np.sum(np.linalg.norm(force, axis=1)))
        self.assertLess(abs(kinetic_power+internal_power)/scale, 1e-10)
        self.assertLess(np.linalg.norm(force.sum(axis=0))/scale, 1e-10)
        self.assertLess(np.linalg.norm(np.cross(b.X, force).sum(axis=0))/(.01*scale), 1e-10)

    def test_passive_labels_cannot_change_three_actual_steps(self):
        a = make_problem(self.cfg, 3, .1)
        b = make_problem(self.cfg, 3, .1)
        b.pa.tracer_type[:] = 23-b.pa.tracer_type
        initial_ids = a.pa.particle_id.copy()
        for _ in range(3):
            a.advance()
            b.advance()
        for left, right in ((positions(a.pa), positions(b.pa)),
                            (velocities(a.pa), velocities(b.pa)),
                            (a.sigma, b.sigma), (a.pa.ep, b.pa.ep), (a.body, b.body)):
            np.testing.assert_array_equal(left, right)
        np.testing.assert_array_equal(a.pa.particle_id, initial_ids)
        np.testing.assert_array_equal(b.pa.particle_id, initial_ids)
        np.testing.assert_array_equal(a.pa.ep, np.zeros(27))
        np.testing.assert_array_equal(a.pa.x[a.fixed], a.X[a.fixed, 0])
        self.assertGreater(np.linalg.norm(positions(a.pa)-a.X), 0.)
        self.assertEqual(a.constraint_work, 0.)
        self.assertEqual(a.stress_time, a.time)

    def test_roller_does_not_freeze_transverse_dofs(self):
        b = make_problem(self.cfg, 3, .1)
        b.pa.v[:] = .01
        b.pa.w[:] = -.02
        b.advance()
        displacement = positions(b.pa)-b.X
        np.testing.assert_allclose(displacement[:, 1], .01*b.dt, rtol=0, atol=1e-15)
        np.testing.assert_allclose(displacement[:, 2], -.02*b.dt, rtol=0, atol=1e-15)
        np.testing.assert_allclose(b.pa.v[b.fixed], .01, rtol=0, atol=1e-12)
        np.testing.assert_allclose(b.pa.w[b.fixed], -.02, rtol=0, atol=1e-12)
        np.testing.assert_array_equal(displacement[b.fixed, 0], np.zeros(9))
        self.assertGreater(float(b.pa.u[~b.fixed].max()), 0.)

    def test_diagnostic_checkpoint_preserves_ids_but_is_not_restart(self):
        b = make_problem(self.cfg, 3, .1)
        b.advance()
        with tempfile.TemporaryDirectory(prefix='afsd_e2_checkpoint_') as folder:
            path = Path(folder)/'diagnostic.npz'
            checkpoint(path, b.pa, dict(classification='E2_DIAGNOSTIC', restart_allowed=False,
                                      time_s=b.time, production_authorized=False))
            data = load(str(path))
            saved = data['arrays']['solid']
            for key in ('particle_id', 'tracer_type', 'x', 'y', 'z', 'u', 'v', 'w', 'S', 'ep'):
                np.testing.assert_array_equal(getattr(saved, key), getattr(b.pa, key))
            self.assertIs(data['solver_data']['restart_allowed'], False)
            with self.assertRaises(ValueError):
                restore(path)


if __name__ == '__main__':
    unittest.main(verbosity=2)

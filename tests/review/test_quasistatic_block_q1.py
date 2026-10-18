"""Independent Q1 references and trial-state contracts; not AFSD validation."""
import copy
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from probe_process_time import QuasiStaticBlock, continuum_reference, load_config
from afsd_backend.pysph.mechanics import internal_force


def accepted_snapshot(block):
    """Capture full accepted particle state, not just selected physical arrays."""
    particles = {}
    for prefix, pa in (('actual', block.pa), ('reference', block.mesh.ref)):
        for key, prop in pa.properties.items():
            array = prop.get_npy_array()
            particles[prefix+'/'+key] = (array.dtype.str, array.shape, array.tobytes())
    arrays = {key: getattr(block, key).copy() for key in (
        'displacement', 'stress', 'external', 'reaction', 'rate')}
    scalars = {key: copy.deepcopy(getattr(block, key)) for key in (
        'time', 'count', 'external_work', 'internal_work', 'constraint_work',
        'last_iterations', 'last_residual', 'last_force_defect',
        'last_torque_defect', 'last_trace')}
    return particles, arrays, scalars


class QuasiStaticBlockQ1(unittest.TestCase):
    def assertAcceptedEqual(self, snapshot, block):
        after = accepted_snapshot(block)
        self.assertEqual(snapshot[0], after[0])
        self.assertEqual(snapshot[2], after[2])
        for key in snapshot[1]:
            np.testing.assert_array_equal(snapshot[1][key], after[1][key])

    def test_literal_continuum_oracle_and_not_a_force_input(self):
        x = np.array([[0., 0., 0.], [.005, .003, .007], [.01, .01, .01]])
        for factor in (0., .5, 1.):
            ref = continuum_reference(x, factor)
            u = np.zeros((3, 3))
            u[:, 0] = factor*np.array([0., np.sqrt(2.)*1e-7/np.pi, 2e-7/np.pi])
            sigma = np.zeros((3, 3, 3))
            sigma[:, 0, 0] = factor*np.array([700000., 700000./np.sqrt(2.), 0.])
            np.testing.assert_allclose(ref['displacement_m'], u, rtol=1e-13, atol=1e-22)
            np.testing.assert_allclose(ref['stress_Pa'], sigma, rtol=1e-13, atol=1e-9)
            np.testing.assert_allclose(ref['reaction_N'], [-70.*factor, 0., 0.], atol=1e-13)
            self.assertAlmostEqual(ref['elastic_J'], factor*factor*1.75e-6, delta=1e-20)
        self.assertIs(load_config()['production_authorized'], False)

    def test_physical_domain_mass_rollers_and_independent_load_quadrature(self):
        for n in (3, 5):
            b = QuasiStaticBlock(n)
            self.assertAlmostEqual(float(b.V.sum()), 1e-6, delta=1e-20)
            self.assertAlmostEqual(float(b.pa.m.sum()), .0027, delta=1e-16)
            np.testing.assert_array_equal(b.X.min(axis=0), [0., 0., 0.])
            np.testing.assert_allclose(b.X.max(axis=0), [.01, .01, .01], atol=1e-17)
            np.testing.assert_array_equal(b.fixed, b.X == 0.)
            self.assertEqual(int(b.fixed.sum()), 3*n*n)
            self.assertEqual(int(b.free.sum()), 3*n*n*(n-1))
            self.assertTrue(np.all(b.fixed[b.X[:, 0] == 0., 0]))
            self.assertTrue(np.any(~b.fixed[b.X[:, 0] == 0., 1]))
            self.assertTrue(np.any(~b.fixed[b.X[:, 0] == 0., 2]))
            q = np.pi/(4*(n-1))
            expected_force = 70.*q/np.tan(q)
            self.assertAlmostEqual(float(b.mesh.body[:, 0].sum()), expected_force, delta=1e-11)
            np.testing.assert_array_equal(b.mesh.body[:, 1:], np.zeros((n**3, 2)))
            self.assertGreater(b.eigen_min, 0.)

    def test_three_dimensional_tangent_against_literal_affine_stress_and_actual_pair(self):
        b = QuasiStaticBlock(3)
        # The stress oracle is continuum sigma=E*sym(A), not the tested D or K.
        A = 1e-6*np.array([[.3, -.7, .4], [1.2, .1, -.6], [-.2, .9, .8]])
        u = b.X@A.T
        sigma = np.broadcast_to(35e9*(A+A.T), (27, 3, 3))
        pair = internal_force(b.mesh.ref, sigma, b.mesh.i, b.mesh.j,
                              b.mesh.g, correction=b.mesh.B)
        matrix = -(b.K@u.ravel()).reshape(27, 3)
        np.testing.assert_allclose(matrix, pair, rtol=1e-10, atol=1e-12)
        expected_energy = .5*70e9*1e-6*np.sum((.5*(A+A.T))**2)
        self.assertAlmostEqual(float(.5*u.ravel()@b.K@u.ravel()), expected_energy, delta=1e-19)
        np.testing.assert_allclose(pair.sum(axis=0), [0., 0., 0.], atol=1e-12)
        np.testing.assert_allclose(np.cross(b.X, pair).sum(axis=0), [0., 0., 0.], atol=1e-14)
        # All three translations and rotations must remain null before support elimination.
        for axis in np.eye(3):
            for rigid in (np.broadcast_to(axis, b.X.shape), np.cross(axis, b.X)):
                defect = np.linalg.norm(b.K@rigid.ravel())
                scale = np.linalg.norm(b.K)*np.linalg.norm(rigid)
                self.assertLess(defect/scale, 1e-13)

    def test_smooth_physical_schedule_and_time_domain(self):
        b = QuasiStaticBlock(3)
        self.assertEqual(b.load_factor(0.), 0.)
        self.assertEqual(b.load_factor(.05), .5)
        self.assertEqual(b.load_factor(.1), 1.)
        self.assertAlmostEqual(b.load_factor(.025), .15625, delta=1e-16)
        for bad in (-.01, .1001, np.nan, np.inf):
            with self.assertRaises(ValueError):
                b.load_factor(bad)

    def test_trial_is_repeatable_and_all_accepted_state_unchanged(self):
        b = QuasiStaticBlock(3)
        before = accepted_snapshot(b)
        u = 2e-6*b.X
        u[b.fixed] = 0.
        first = b.trial(u, .05)
        b.trial(.7*u, .075)
        again = b.trial(u, .05)
        for key in first:
            np.testing.assert_array_equal(first[key], again[key])
        self.assertAcceptedEqual(before, b)
        second_time = b.trial(u, .1)
        np.testing.assert_allclose(first['S'], second_time['S'], rtol=1e-13, atol=1e-10)
        np.testing.assert_allclose(first['strain_rate'], 2.*second_time['strain_rate'], rtol=1e-13, atol=1e-20)
        self.assertEqual(first['dt_s'], .05)
        self.assertEqual(second_time['dt_s'], .1)

    def test_failure_and_unconverged_candidate_roll_back_everything(self):
        b = QuasiStaticBlock(3)
        before = accepted_snapshot(b)
        with self.assertRaisesRegex(RuntimeError, 'did not converge'):
            b.solve_candidate(.1, maximum_iterations=0)
        self.assertAcceptedEqual(before, b)
        with self.assertRaisesRegex(ValueError, 'unconverged'):
            b.commit(b.trial(np.zeros_like(b.X), .1))
        self.assertAcceptedEqual(before, b)
        u = np.zeros_like(b.X); u[b.fixed] = 1e-9
        with self.assertRaisesRegex(ValueError, 'support'):
            b.trial(u, .1)
        self.assertAcceptedEqual(before, b)

    def test_commit_rechecks_mutable_candidate_and_rejects_stale_state(self):
        b = QuasiStaticBlock(3)
        before = accepted_snapshot(b)
        candidate = b.solve_candidate(.1)
        self.assertAcceptedEqual(before, b)
        expected = b.trial(candidate['displacement'], .1)
        candidate['stress'][:] = 1e90
        candidate['S'][:] = -1e90
        candidate['residual'][:] = 0.
        b.commit(candidate)
        np.testing.assert_array_equal(b.stress, expected['stress'])
        np.testing.assert_array_equal(b.pa.S.reshape(-1, 3, 3), expected['S'])
        after = accepted_snapshot(b)
        with self.assertRaisesRegex(ValueError, 'stale'):
            b.commit(candidate)
        self.assertAcceptedEqual(after, b)

    def test_equilibrium_reaction_work_and_no_hidden_dynamic_velocity(self):
        b = QuasiStaticBlock(3)
        ids = b.pa.particle_id.copy()
        for t in (.025, .05, .075, .1):
            b.advance(t)
        record = b.record()
        expected_reaction = -70.*(np.pi/8)/np.tan(np.pi/8)
        self.assertAlmostEqual(float(b.reaction[:, 0].sum()), expected_reaction, delta=1e-7)
        np.testing.assert_allclose((b.external+b.reaction).sum(axis=0), np.zeros(3), atol=1e-7)
        self.assertLess(abs(record['work_identity_residual_J']), 1.75e-13)
        self.assertEqual(record['constraint_work_J'], 0.)
        self.assertGreater(record['external_work_J'], 0.)
        self.assertGreater(np.linalg.norm(b.displacement), 0.)
        np.testing.assert_array_equal(b.displacement[b.fixed], np.zeros(int(b.fixed.sum())))
        np.testing.assert_array_equal(b.pa.ep, np.zeros(27))
        np.testing.assert_array_equal(b.pa.particle_id, ids)
        for name in ('u', 'v', 'w'):
            np.testing.assert_array_equal(getattr(b.pa, name), np.zeros(27))
        self.assertGreater(np.linalg.norm(b.rate), 0.)
        self.assertEqual(b.time, .1)
        self.assertEqual(b.count, 4)

    def test_passive_label_changes_do_not_change_equilibrium(self):
        a, b = QuasiStaticBlock(3), QuasiStaticBlock(3)
        b.pa.tracer_type[:] = 43-7*b.pa.tracer_type
        a.advance(.1); b.advance(.1)
        for name in ('displacement', 'stress', 'external', 'reaction', 'rate'):
            np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
        np.testing.assert_array_equal(a.pa.S, b.pa.S)
        np.testing.assert_array_equal(a.pa.ep, b.pa.ep)
        self.assertEqual(a.record(), b.record())


if __name__ == '__main__':
    unittest.main(verbosity=2)

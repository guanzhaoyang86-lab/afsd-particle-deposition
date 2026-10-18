"""Independent Q2 material references and transactions; NOT coupled AFSD evidence."""
import copy
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from probe_process_plastic import PlasticBlock, CandidateFailure
from probe_engineering_material import PerfectPlasticBenchmark
from afsd_backend.pysph.material import update, FKS
from afsd_backend.pysph.mechanics import internal_force

E, G, K, Y, VOLUME = 70e9, 35e9, 70e9/3, 90e6, 1e-6


def independent_uniaxial(e, loading):
    """Scalar continuum J2 solution; no solver/contact/gradient calls."""
    ey = 90e6/70e9
    if loading:
        sigma = 70e9*e if e >= -ey else -90e6
        ep = 0. if e >= -ey else -e-ey
    else:
        travel = e+.004
        sigma = -90e6+70e9*travel if travel <= 2*ey else 90e6
        ep = .004-ey if travel <= 2*ey else .004-ey+travel-2*ey
    signed_plastic = e-sigma/70e9
    strain = np.diag([e, -.5*signed_plastic, -.5*signed_plastic])
    return strain, np.diag([sigma, 0., 0.]), ep


def material_path(points):
    S = np.zeros((1, 3, 3)); ep = np.zeros(1); strain = np.zeros((3, 3))
    P = returned = 0.; history = []
    for e, loading in points:
        target, reference, reference_ep = independent_uniaxial(e, loading)
        S, ep, plastic, ret, rate = update(S, ep, ((target-strain)/.01)[None],
            np.array([293.15]), .01, G, PerfectPlasticBenchmark(Y, 293.15))
        sigma = S[0]+K*np.trace(target)*np.eye(3)
        P += float(plastic[0])*VOLUME; returned += float(ret[0])*VOLUME
        history.append((sigma.copy(), float(ep[0]), reference, reference_ep,
                        target.copy(), P, returned))
        strain = target
    return history


def snapshot(b):
    particles = {}
    for prefix, pa in (('actual', b.pa), ('reference', b.mesh.ref)):
        for name, prop in pa.properties.items():
            a = prop.get_npy_array()
            particles[prefix+'/'+name] = (a.dtype.str, a.shape, a.tobytes())
    arrays = {name: getattr(b, name).copy() for name in
              ('displacement', 'stress', 'external', 'reaction', 'rate')}
    scalar = {name: copy.deepcopy(getattr(b, name)) for name in
              ('time', 'count', 'external_work', 'internal_work', 'constraint_work',
               'last_iterations', 'last_trace', 'last_residual', 'last_force_defect',
               'last_torque_defect', 'accepted_floor_count')}
    return particles, arrays, scalar


class IndependentMaterialQ2(unittest.TestCase):
    def test_literal_compression_and_reverse_yield_values(self):
        for subdivisions in (8, 16):
            points = [(e, True) for e in np.linspace(0., -.004, subdivisions+1)[1:]]
            points += [(e, False) for e in np.linspace(-.004, 0., subdivisions+1)[1:]]
            history = material_path(points)
            for stress, ep, reference, expected_ep, _, _, _ in history:
                np.testing.assert_allclose(stress, reference, rtol=2e-12, atol=2e-5)
                self.assertAlmostEqual(ep, expected_ep, delta=2e-15)
            peak, final = history[subdivisions-1], history[-1]
            self.assertAlmostEqual(peak[0][0, 0], -90e6, delta=2e-5)
            self.assertAlmostEqual(peak[1], .002714285714285714, delta=2e-15)
            self.assertAlmostEqual(final[0][0, 0], 90e6, delta=2e-5)
            self.assertAlmostEqual(final[1], .004142857142857143, delta=2e-15)
            self.assertAlmostEqual(final[5], .37285714285714283, delta=2e-12)
            self.assertAlmostEqual(final[4][1, 1], .0006428571428571428, delta=1e-17)
            self.assertGreater(final[6], 0.)

    def test_zero_stress_unload_is_not_return_to_zero_axial_displacement(self):
        h = material_path([(-.004, True), (-.002714285714285714, False)])
        np.testing.assert_allclose(h[-1][0], np.zeros((3, 3)), atol=2e-5)
        self.assertAlmostEqual(h[-1][1], .002714285714285714, delta=2e-15)
        self.assertAlmostEqual(h[-1][5], .24428571428571427, delta=2e-12)
        self.assertNotEqual(float(h[-1][4][0, 0]), 0.)

    def test_predictor_return_energy_identity_is_not_trapezoid_work(self):
        # Symmetric trace-free affine increment gives independent q_trial=210MPa.
        inc = np.diag([-.002, .001, .001])[None]
        S, ep, plastic, returned, _ = update(np.zeros((1, 3, 3)), np.zeros(1),
            inc/.05, np.array([293.15]), .05, G, PerfectPlasticBenchmark(Y, 293.15))
        dep = (210e6-90e6)/(3*35e9)
        expected_plastic = 90e6*dep
        expected_return = 1.5*35e9*dep*dep
        loss = (210e6**2-90e6**2)/(6*35e9)
        self.assertAlmostEqual(float(ep[0]), dep, delta=2e-15)
        self.assertAlmostEqual(float(plastic[0]), expected_plastic, delta=2e-6)
        self.assertAlmostEqual(float(returned[0]), expected_return, delta=2e-6)
        self.assertAlmostEqual(float(plastic[0]+returned[0]), loss, delta=2e-6)
        Wtrap = float(np.sum(.5*S[0]*inc[0]))
        Unew = 90e6**2/(6*35e9)
        self.assertGreater(abs(Wtrap-Unew-float(plastic[0])), 1e4)

    def test_fks_independent_formula_and_disclosed_rate_floor(self):
        # Direct paper-form arithmetic independent of FKS's stable exp/log implementation.
        law = FKS()
        for ep, rate, temperature in ((0., .04, 723.15), (.002, 10., 600.), (.01, 100., 293.15)):
            ts = (temperature-293.15)/(855.15-293.15)
            effective_rate = rate if rate >= 1. else 1.
            expected = (317.6e6+1e6*np.arctan(70.*ep))
            expected *= 1.+12.44*ts**5.65*np.log(effective_rate)
            expected *= 1.-(1.+np.exp(-13.*ts))**(-100.)
            self.assertAlmostEqual(float(law.yield_stress(ep, rate, temperature))/expected, 1., delta=2e-11)
        self.assertAlmostEqual(float(law.yield_stress(0., .04, 723.15)), 1517290.6902550498, delta=2e-5)
        self.assertEqual(float(law.yield_stress(0., .04, 723.15)), float(law.yield_stress(0., 1., 723.15)))

    def test_fks_invalid_temperature_and_inputs_are_rejected_not_clipped(self):
        law = FKS()
        for temperature in (293.14, 855.15, 900., np.nan, np.inf):
            with self.assertRaises(ValueError):
                law.yield_stress(0., .1, temperature)
        for ep, rate in ((-1e-6, .1), (0., -1.), (np.nan, 1.)):
            with self.assertRaises(ValueError):
                law.yield_stress(ep, rate, 723.15)


class PlasticBlockQ2(unittest.TestCase):
    def assertAcceptedEqual(self, before, b):
        after = snapshot(b)
        self.assertEqual(before[0], after[0]); self.assertEqual(before[2], after[2])
        for name in before[1]:
            np.testing.assert_array_equal(before[1][name], after[1][name])

    def affine_candidate(self, b, t=.05):
        strain, _, _ = independent_uniaxial(b.cap_displacement(t)/.01, True)
        u = b.X@strain.T
        u[b.fixed] = b.prescribed(t)[b.fixed]
        return b.trial(u, t)

    def test_three_dimensional_domain_rollers_and_free_tangential_components(self):
        b = PlasticBlock(3)
        self.assertEqual(len(b.X), 27); self.assertEqual(int(b.free.sum()), 45)
        self.assertTrue(all(len(np.unique(b.X[:, a])) == 3 for a in range(3)))
        self.assertAlmostEqual(float(b.V.sum()), 1e-6, delta=1e-20)
        self.assertAlmostEqual(float(b.pa.m.sum()), .0027, delta=1e-16)
        expected = b.X == 0.; expected[b.X[:, 0] == .01, 0] = True
        np.testing.assert_array_equal(b.fixed, expected)
        self.assertTrue(np.any(~b.fixed[b.cap, 1])); self.assertTrue(np.any(~b.fixed[b.cap, 2]))
        for t, cap in ((0., 0.), (.05, -2e-5), (.1, -4e-5), (.15, -2e-5), (.2, 0.)):
            self.assertAlmostEqual(b.cap_displacement(t), cap, delta=1e-19)
        for t in (-1., .20001, np.nan):
            with self.assertRaises(ValueError): b.prescribed(t)

    def test_elastic_tangent_has_independent_affine_energy_and_force_reference(self):
        b = PlasticBlock(3)
        u = -.0001*b.X*np.array([1., 0., 0.]); u[b.fixed] = b.prescribed(.0025)[b.fixed]
        candidate = b.trial(u, .0025); tangent = b.tangent(candidate)
        A = 1e-6*np.array([[.3, -.7, .4], [1.2, .1, -.6], [-.2, .9, .8]])
        direction = b.X@A.T
        # Existing predictor has zero accepted S, but its finite-step rotation can
        # contribute O(base strain); use symmetric perturbation to set spin exactly 0.
        A = .5*(A+A.T); direction = b.X@A.T
        independent_stress = np.broadcast_to(70e9*A, (27, 3, 3))
        force = internal_force(b.mesh.ref, independent_stress, b.mesh.i, b.mesh.j,
                               b.mesh.g, correction=b.mesh.B)
        matrix = -(tangent@direction.ravel()).reshape(27, 3)
        np.testing.assert_allclose(matrix, force, rtol=2e-8, atol=1e-7)
        expected_energy = .5*70e9*1e-6*np.sum(A*A)
        self.assertAlmostEqual(float(.5*direction.ravel()@tangent@direction.ravel()), expected_energy, delta=2e-14)

    def test_plastic_tangent_matches_actual_force_differences_at_other_scales(self):
        b = PlasticBlock(3); c = self.affine_candidate(b)
        self.assertTrue(np.all(c['ep'] > 0.))
        before = snapshot(b); tangent = b.tangent(c)
        direction = np.random.default_rng(2042).normal(size=b.X.shape)
        direction[b.fixed] = 0.; direction /= np.linalg.norm(direction)
        expected = -(tangent@direction.ravel()).reshape(b.X.shape)
        for delta_m in (2e-10, 7e-10):
            plus = b.trial(c['displacement']+delta_m*direction, .05)['internal_force']
            minus = b.trial(c['displacement']-delta_m*direction, .05)['internal_force']
            actual = (plus-minus)/(2*delta_m)
            self.assertLess(np.linalg.norm(actual-expected)/np.linalg.norm(actual), 3e-6)
        self.assertAcceptedEqual(before, b)

    def test_trial_and_tangent_are_repeatable_without_any_accepted_state_mutation(self):
        b = PlasticBlock(3); before = snapshot(b); first = self.affine_candidate(b)
        self.assertGreater(float(first['plastic_increment_J'].sum()), 0.)
        self.affine_candidate(b, .025); b.tangent(first); second = self.affine_candidate(b)
        for name in first: np.testing.assert_array_equal(first[name], second[name])
        self.assertAcceptedEqual(before, b)

    def test_physical_time_enters_rate_not_iteration_number(self):
        b = PlasticBlock(3); gradient = np.broadcast_to(np.diag([-.002, .001, .001]), (27, 3, 3)).copy()
        first = b.constitutive(gradient, .05); second = b.constitutive(gradient, .1)
        for index in (0, 1, 2, 3, 5): np.testing.assert_allclose(first[index], second[index], rtol=1e-13, atol=1e-10)
        np.testing.assert_allclose(first[4], 2.*second[4], rtol=1e-14)

    def test_unconverged_failure_and_invalid_boundary_leave_full_state_unchanged(self):
        b = PlasticBlock(3); before = snapshot(b)
        with self.assertRaises(CandidateFailure): b.solve_candidate(.05, maximum_iterations=0)
        self.assertAcceptedEqual(before, b)
        # Uniform affine n3 fields may already balance: manufacture a genuinely
        # unbalanced free interior DOF, not a presumed failure by function name.
        bad_u = self.affine_candidate(b)['displacement'].copy()
        center = np.flatnonzero(np.all(b.X == .005, axis=1))[0]
        self.assertFalse(b.fixed[center, 0]); bad_u[center, 0] += 1e-7
        unbalanced = b.trial(bad_u, .05)
        self.assertGreater(np.linalg.norm(unbalanced['residual'].ravel()[b.free])/b.force_scale, 1e-8)
        with self.assertRaisesRegex(ValueError, 'unconverged'): b.commit(unbalanced)
        self.assertAcceptedEqual(before, b)
        u = np.zeros_like(b.X)
        with self.assertRaisesRegex(ValueError, 'support/actuator'): b.trial(u, .05)
        self.assertAcceptedEqual(before, b)

    def test_commit_recomputes_candidate_and_rejects_stale_step(self):
        b = PlasticBlock(3); before = snapshot(b)
        c = b.solve_candidate(.05); self.assertAcceptedEqual(before, b)
        expected = b.trial(c['displacement'], .05)
        c['stress'][:] = 1e90; c['S'][:] = -1e90; c['ep'][:] = 1e90
        c['plastic_increment_J'][:] = 1e90; c['return_increment_J'][:] = 1e90
        c['residual'][:] = 0.
        b.commit(c)
        np.testing.assert_array_equal(b.stress, expected['stress'])
        np.testing.assert_array_equal(b.pa.S.reshape(-1, 3, 3), expected['S'])
        np.testing.assert_array_equal(b.pa.ep, expected['ep'])
        np.testing.assert_array_equal(b.pa.plastic_work, expected['plastic_increment_J'])
        np.testing.assert_array_equal(b.pa.return_dissipation, expected['return_increment_J'])
        after = snapshot(b)
        with self.assertRaisesRegex(ValueError, 'stale'): b.commit(c)
        self.assertAcceptedEqual(after, b)

    def test_single_plastic_equilibrium_has_actual_reaction_and_separate_energy_ledgers(self):
        b = PlasticBlock(3); ids = b.pa.particle_id.copy(); b.advance(.05)
        r = b.record()
        self.assertGreater(r['maximum_ep'], 0.); self.assertGreater(r['plastic_work_J'], 0.)
        self.assertGreater(r['return_dissipation_diagnostic_J'], 0.)
        self.assertLess(abs(r['work_identity_residual_J']), 1e-6)
        self.assertLess(r['free_force_relative'], 1e-8)
        self.assertNotEqual(r['work_minus_elastic_minus_plastic_J'], 0.)
        np.testing.assert_allclose(b.reaction.sum(axis=0), [0., 0., 0.], atol=1e-3)
        np.testing.assert_array_equal(b.pa.particle_id, ids)
        for name in ('u', 'v', 'w'): np.testing.assert_array_equal(getattr(b.pa, name), np.zeros(27))
        np.testing.assert_array_equal(b.pa.temperature, np.full(27, 293.15))

    def test_passive_labels_do_not_enter_plastic_solution(self):
        a, b = PlasticBlock(3), PlasticBlock(3)
        b.pa.tracer_type[:] = 53-11*b.pa.tracer_type
        a.advance(.05); b.advance(.05)
        for name in ('displacement', 'stress', 'reaction', 'rate'):
            np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
        for name in ('S', 'ep', 'plastic_work', 'return_dissipation', 'strain_rate'):
            np.testing.assert_array_equal(getattr(a.pa, name), getattr(b.pa, name))
        self.assertEqual(a.record(), b.record())

    def test_fks_invalid_temperature_trial_does_not_mutate_state(self):
        b = PlasticBlock(3, 'fks'); b.pa.temperature[:] = 855.15
        before = snapshot(b)
        with self.assertRaisesRegex(ValueError, 'temperature'): self.affine_candidate(b)
        self.assertAcceptedEqual(before, b)


if __name__ == '__main__':
    unittest.main(verbosity=2)

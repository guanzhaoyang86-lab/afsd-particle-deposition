"""Independent event-Schur algebra/linear-map baselines, not AFSD trials."""
from copy import deepcopy
from pathlib import Path
import hashlib
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.stage_event_preconditioner import StageEventPreconditioner, CLASSIFICATION
from afsd_backend.pysph.stage_continuity_preconditioner import ContinuityStagePreconditioner
from afsd_backend.pysph.stage_iteration import FIELDS
from afsd_backend.pysph.structural_mobility import FrozenElastic, ConstrainedResponse
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.updated_trial import particle_hash


def fields(n=2):
    return dict(position_m=np.zeros((n, 3)), velocity_m_s=np.zeros((n, 3)),
        density_kg_m3=np.array([5., 4.])[:n], temperature_K=np.full(n, 400.),
        mean_velocity_m_s=np.zeros((n, 3)))


def copied(y):
    return {k: v.copy() for k, v in y.items()}


def elastic_fixture(compression=False):
    # Two original DOFs share one spring, kdev=4. Compression has independent
    # literal E rows (1,-1) and (.5,-.5), V=(2,3), kappa=3.
    e = np.zeros((2, 6))
    if compression:
        e[:, [0, 3]] = [[1., -1.], [.5, -.5]]
    gradient = np.zeros((2, 3, 3, 6)); gradient[:, 0, 0] = e
    k = np.zeros((6, 6)); k[np.ix_([0, 3], [0, 3])] = [[4., -4.], [-4., 4.]]
    k += e.T@np.diag([6., 9.])@e
    return FrozenElastic(gradient.reshape(18, 6), k, [2., 3.], [1., 3.],
                         np.array([[0., 0., 0.], [1., 0., 0.]]),
                         np.array([19, 4]), 2., 3., 1)


def event_blocks():
    b = np.zeros((2, 3, 3)); c = np.zeros_like(b)
    b[:, 0, 0] = [19/54, 1/6]; c[:, 0, 0] = [7/9, 1/3]
    return b, c


def spring(compression=False, constrained=False):
    b, c = event_blocks(); free = np.ones((2, 3), bool)
    if constrained:
        free[1, 0] = False; b[1, 0] = 0.; c[1, 0] = 0.
    return StageEventPreconditioner.from_response_blocks(elastic_fixture(compression),
        b, c, 1., free, [5., 4.], {'particle_id': np.array([19, 4]),
        'scope': 'LITERAL_MANUFACTURED_LINEAR_RESPONSE_NOT_CONTACT_VALIDATION'})


class EventSchurAlgebraTests(unittest.TestCase):
    def close(self, a, b):
        np.testing.assert_allclose(a, b, rtol=2e-12, atol=2e-13)

    def test_01_literal_event_reversal_is_not_endpoint_half(self):
        # Independent analytic interval: m=1,u0=1,F=-2,muFn=1,dt=1.
        # tau=1/(1-F), v1=(F+1)(1-tau),
        # bar=tau/2+(F+1)(1-tau)^2/2. Differentiate these formulas.
        force = -2.; tau = 1/(1-force); tau_prime = 1/(1-force)**2
        endpoint_derivative = (1-tau)-(force+1)*tau_prime
        mean_derivative = tau_prime/2+(1-tau)**2/2-(force+1)*(1-tau)*tau_prime
        self.close(mean_derivative, 19/54); self.close(endpoint_derivative, 7/9)
        self.assertNotAlmostEqual(mean_derivative, endpoint_derivative/2)
        b, c = event_blocks()
        self.close(b[0, 0, 0], 19/54); self.close(c[0, 0, 0], 7/9)

    def test_02_literal_two_dof_nonsymmetric_shared_elastic_solve(self):
        candidate = fields(); mapping = copied(candidate)
        mapping['mean_velocity_m_s'][:, 0] = [1., 2.]
        pc = spring(); got = pc.apply(candidate, mapping); delta = got['correction']
        # A=[[46/27,-19/27],[-1/3,4/3]], rhs=(1,2).
        self.close(delta['mean_velocity_m_s'][:, 0], [74/55, 101/55])
        self.close(delta['position_m'][:, 0], [19/55, -9/55])
        self.close(delta['velocity_m_s'][:, 0], [42/55, -18/55])
        self.close(got['diagnostics']['delta_internal_force_N'][[0, 3]], [54/55, -54/55])
        self.assertFalse(np.allclose(pc.matrix, pc.matrix.T))
        self.assertGreater(abs(delta['velocity_m_s'][0, 0]-2*delta['position_m'][0, 0]), .01)
        self.assertEqual(got['diagnostics']['classification'], CLASSIFICATION)

    def test_03_nonzero_prescribed_axis_couples_into_other_particle(self):
        candidate = fields(); mapping = copied(candidate)
        mapping['mean_velocity_m_s'][:, 0] = [1., 2.]
        mapping['position_m'][1, 0] = .3; mapping['velocity_m_s'][1, 0] = -.7
        got = spring(constrained=True).apply(candidate, mapping); delta = got['correction']
        self.close(delta['mean_velocity_m_s'][:, 0], [65/46, 2.])
        self.close(delta['position_m'][:, 0], [19/46, .3])
        self.close(delta['velocity_m_s'][:, 0], [21/23, -.7])
        self.close(got['diagnostics']['delta_internal_force_N'][[0, 3]], [27/23, -27/23])
        self.close(got['diagnostics']['algebraic_constrained_row_defect_m_s'], [0.])
        self.close(got['diagnostics']['prescribed_delta_defect_m_s'], [0.])

    def test_04_independent_density_elastic_four_by_four_root(self):
        # Literal force path and coupled root; never evaluate with the new E/K/H.
        def physical_linear_map(y):
            bar = y['mean_velocity_m_s'][:, 0]
            r = y['density_kg_m3']-np.array([5., 4.])
            force = -2*(bar[0]-bar[1])+.6*r[0]+.5625*r[1]
            out = fields()
            out['mean_velocity_m_s'][:, 0] = [.3+(19/54)*force, -.2-force/6]
            out['velocity_m_s'][:, 0] = [.7+(7/9)*force, .4-force/3]
            out['position_m'][:, 0] = [.1, -.4]+out['mean_velocity_m_s'][:, 0]
            out['density_kg_m3'] = [5.2-5*(bar[0]-bar[1]), 3.9-2*(bar[0]-bar[1])]
            out['temperature_K'] = np.array([401., 399.])
            return out
        a = np.array([[1+19/27, -19/27, -.6*19/54, -.5625*19/54],
                      [-1/3, 1+1/3, .1, .5625/6],
                      [5., -5., 1., 0.], [2., -2., 0., 1.]])
        reference = np.linalg.solve(a, [.3, -.2, .2, -.1])
        candidate = fields(); candidate['mean_velocity_m_s'][:, 0] = [.8, -.4]
        candidate['density_kg_m3'] += [1.2, -.3]
        candidate['position_m'][:, 0] = [3., -1.]
        root = spring(compression=True).apply(candidate, physical_linear_map(candidate))['proposal']
        self.close(root['mean_velocity_m_s'][:, 0], reference[:2])
        self.close(root['density_kg_m3'], reference[2:]+[5., 4.])
        for key in FIELDS:
            self.close(root[key], physical_linear_map(root)[key])

    def test_05_free_inertial_limit_reduces_to_old_continuity_schur(self):
        elastic = elastic_fixture(compression=True); dt = .7
        rng = np.random.default_rng(718)
        candidate = fields(); mapping = {k: v+rng.normal(size=v.shape) for k, v in candidate.items()}
        for constrained in (False, True):
            free = np.ones((2, 3), bool); free[1, 0] = not constrained
            b = dt/(2*elastic.mass[:, None, None])*np.broadcast_to(np.eye(3), (2, 3, 3))
            b[~free] = 0.; c = 2*b
            old_response = ConstrainedResponse(np.repeat(elastic.mass, 3), elastic.stiffness, dt, free.ravel())
            e = np.zeros((2, 6)); e[:, [0, 3]] = [[1., -1.], [.5, -.5]]
            old = ContinuityStagePreconditioner.from_response(old_response, elastic.particle_ids,
                e, elastic.volume, [5., 4.], 3.)
            new = StageEventPreconditioner.from_response_blocks(elastic, b, c, dt, free, [5., 4.])
            left, right = old.apply(candidate, mapping), new.apply(candidate, mapping)
            for key in FIELDS:
                self.close(left['correction'][key], right['correction'][key])

    def test_06_singular_B_zero_responses_and_rigid_translation_are_valid(self):
        pc = spring(); candidate = fields(); mapping = copied(candidate)
        mapping['mean_velocity_m_s'][:, 0] = 3.
        got = pc.apply(candidate, mapping)
        self.close(got['correction']['mean_velocity_m_s'][:, 0], [3., 3.])
        self.close(got['correction']['position_m'], np.zeros((2, 3)))
        zero = np.zeros((2, 3, 3))
        pc = StageEventPreconditioner.from_response_blocks(elastic_fixture(True), zero, zero,
            1., np.ones((2, 3), bool), [5., 4.])
        mapping['mean_velocity_m_s'][:, 0] = [1., -1.]
        got = pc.apply(candidate, mapping)
        self.close(got['correction']['mean_velocity_m_s'], mapping['mean_velocity_m_s'])
        self.close(got['correction']['density_kg_m3'], [-10., -4.])

    def test_07_singular_schur_fails_without_pseudoinverse(self):
        elastic = elastic_fixture(); k = np.zeros((6, 6)); k[0, 0] = 2.
        malformed = FrozenElastic(elastic.gradient_map, k, elastic.volume, elastic.mass,
            elastic.position, elastic.particle_ids, elastic.shear, elastic.bulk, 1)
        b, c = event_blocks(); b[0, 0, 0] = -1.
        pc = StageEventPreconditioner.from_response_blocks(malformed, b, c, 1.,
            np.ones((2, 3), bool), [5., 4.])
        with self.assertRaisesRegex(ValueError, 'singular event Schur matrix'):
            pc.apply(fields(), fields())

    def test_08_full_residual_backward_error_and_matrix_hash_are_independent(self):
        pc = spring(compression=True, constrained=True); candidate = fields()
        rng = np.random.default_rng(92)
        mapping = {k: v+rng.normal(size=v.shape) for k, v in candidate.items()}
        got = pc.apply(candidate, mapping); info = got['diagnostics']
        delta = got['correction']['mean_velocity_m_s'].ravel()
        rhs = info['rhs_m_s']; residual = pc.matrix@delta-rhs
        numerator = np.linalg.norm(residual, np.inf)
        denominator = np.linalg.norm(pc.matrix, np.inf)*np.linalg.norm(delta, np.inf)+np.linalg.norm(rhs, np.inf)
        self.close(info['full_linear_residual_m_s'], residual)
        self.close(info['backward_error_inf'], numerator/denominator)
        self.assertLessEqual(info['backward_error_inf'], 1e-12)
        self.assertLessEqual(info['free_backward_error_inf'], 1e-12)
        self.close(info['mean_backsubstitution_defect_m_s'], np.zeros(6))
        self.assertEqual(info['matrix_sha256'], hashlib.sha256(pc.matrix.tobytes()).hexdigest())
        self.assertFalse(pc.matrix.flags.writeable)
        self.assertNotIn('matrix', info)

    def test_09_zero_defect_preserves_independent_event_mean(self):
        candidate = fields(); candidate['mean_velocity_m_s'][:, 0] = -1/18
        candidate['velocity_m_s'][:, 0] = -2/3
        got = spring(compression=True).apply(candidate, copied(candidate))
        for key in FIELDS:
            np.testing.assert_array_equal(got['proposal'][key], candidate[key])
            np.testing.assert_array_equal(got['correction'][key], np.zeros_like(candidate[key]))
        self.assertEqual(got['diagnostics']['backward_error_inf'], 0.)

    def test_10_input_output_provenance_purity_and_no_clipping(self):
        b, c = event_blocks(); before_b, before_c = b.copy(), c.copy()
        provenance = {'particle_id': np.array([19, 4]), 'nested': {'array': np.array([3.])}}
        original = deepcopy(provenance)
        pc = StageEventPreconditioner.from_response_blocks(elastic_fixture(), b, c, 1.,
            np.ones((2, 3), bool), [5., 4.], provenance)
        candidate = fields(); mapping = copied(candidate)
        mapping['temperature_K'][:] = -1.; mapping['density_kg_m3'][:] = -3.
        before_candidate, before_mapping = copied(candidate), copied(mapping)
        got = pc.apply(candidate, mapping)
        np.testing.assert_array_equal(got['proposal']['temperature_K'], [-1., -1.])
        np.testing.assert_array_equal(got['proposal']['density_kg_m3'], [-3., -3.])
        got['diagnostics']['event_response_provenance']['nested']['array'][:] = 99.
        got['diagnostics']['B_s_kg'][:] = 99.; got['diagnostics']['particle_id'][:] = 99
        got['diagnostics']['direct_defect']['position_m'][:] = 99.
        b[:] = 88.; c[:] = 77.; provenance['nested']['array'][:] = 66.
        np.testing.assert_array_equal(pc.B, before_b); np.testing.assert_array_equal(pc.C, before_c)
        np.testing.assert_array_equal(pc.metadata['event_response_provenance']['nested']['array'],
                                      original['nested']['array'])
        for key in FIELDS:
            np.testing.assert_array_equal(candidate[key], before_candidate[key])
            np.testing.assert_array_equal(mapping[key], before_mapping[key])
        self.assertFalse(got['diagnostics']['accepted']); self.assertFalse(got['diagnostics']['physical_map'])

    def test_11_invalid_blocks_order_constraints_and_fields_reject(self):
        elastic = elastic_fixture(); b, c = event_blocks(); free = np.ones((2, 3), bool)
        for invalid in (np.zeros((2, 3)), np.full((2, 3, 3), np.nan)):
            with self.assertRaises(ValueError):
                StageEventPreconditioner.from_response_blocks(elastic, invalid, c, 1., free, [5., 4.])
        with self.assertRaisesRegex(ValueError, 'order mismatch'):
            StageEventPreconditioner.from_response_blocks(elastic, b, c, 1., free, [5., 4.],
                                                         {'particle_id': [4, 19]})
        free[0, 0] = False
        with self.assertRaisesRegex(ValueError, 'exactly zero'):
            StageEventPreconditioner.from_response_blocks(elastic, b, c, 1., free, [5., 4.])
        for key in FIELDS:
            invalid = fields(); invalid[key].flat[0] = np.inf
            with self.assertRaises(ValueError): spring().apply(fields(), invalid)

    def test_12_offdiagonal_nonsymmetric_B_retains_force_cross_axis_response(self):
        elastic = elastic_fixture(); b, c = event_blocks()
        b[0, 1, 0] = .2; c[0, 2, 0] = -.4
        pc = StageEventPreconditioner.from_response_blocks(elastic, b, c, 1.,
            np.ones((2, 3), bool), [5., 4.])
        candidate = fields(); mapping = copied(candidate)
        mapping['mean_velocity_m_s'][:, 0] = [1., 2.]
        got = pc.apply(candidate, mapping)
        self.close(got['correction']['mean_velocity_m_s'][0, 1], .2*54/55)
        self.close(got['correction']['position_m'][0, 1], .2*54/55)
        self.close(got['correction']['velocity_m_s'][0, 2], -.4*54/55)


class AcceptedFrozenSPHTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        xyz = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                                   [-.0005, 0., .0005], indexing='ij')).reshape(3, -1).T
        cls.p = dict(dx=.0005, h_ratio=1.3, rho0=2700., width=.0015,
            marker_low=-.001, marker_high=.001, bar_length=.01, actuator_depth=.0005,
            temperature=400., traverse=.00212, feed=.00212, shear=26.3e9,
            poisson=.33, release_time=1.)
        cls.pa = make_array(xyz, np.zeros(27, int), cls.p)
        cls.dt = 2e-7
        cls.b = cls.dt/(2*cls.pa.m[:, None, None])*np.broadcast_to(np.eye(3), (27, 3, 3))
        cls.c = 2*cls.b
        cls.before = particle_hash(cls.pa)
        cls.pc = StageEventPreconditioner(cls.pa, cls.p, cls.dt, 0., cls.b, cls.c,
            {'particle_id': cls.pa.particle_id.copy(), 'scope': 'FREE_FROZEN_LINEAR_MAP_ONLY'})
        cls.candidate = dict(position_m=xyz.copy(), velocity_m_s=np.zeros_like(xyz),
            mean_velocity_m_s=np.zeros_like(xyz), density_kg_m3=cls.pa.rho.copy(),
            temperature_K=cls.pa.temperature.copy())

    def test_13_passive_labels_particle_order_and_source_history(self):
        q = self.pa.extract_particles(np.arange(27, dtype=np.uint32))
        q.tracer_type[:] = np.arange(27) % 3; q.wholebar_tag[:] = np.arange(27) % 4
        other = StageEventPreconditioner(q, self.p, self.dt, 0., self.b, self.c, {})
        np.testing.assert_array_equal(self.pc.matrix, other.matrix)
        mapping = copied(self.candidate)
        mapping['mean_velocity_m_s'] = np.sin(np.arange(81)).reshape(27, 3)
        baseline = self.pc.apply(self.candidate, mapping)
        order = np.arange(26, -1, -1, dtype=np.uint32)
        permuted = q.extract_particles(order)
        perm_pc = StageEventPreconditioner(permuted, self.p, self.dt, 0., self.b[order],
            self.c[order], {'particle_id': permuted.particle_id.copy()})
        actual = perm_pc.apply({k: v[order] for k, v in self.candidate.items()},
                               {k: v[order] for k, v in mapping.items()})
        for key in FIELDS:
            np.testing.assert_allclose(actual['correction'][key], baseline['correction'][key][order],
                                       rtol=2e-10, atol=3e-13)
        self.assertEqual(self.before, particle_hash(self.pa))
        np.testing.assert_array_equal(perm_pc.particle_ids, self.pa.particle_id[order])

    def test_14_all_six_rigid_modes_retained(self):
        xyz = self.candidate['position_m']
        for axis in np.eye(3):
            for mode in (np.broadcast_to(axis, (27, 3)), np.cross(axis, xyz)):
                mapping = copied(self.candidate); mapping['mean_velocity_m_s'] = mode.copy()
                got = self.pc.apply(self.candidate, mapping)['correction']
                self.assertLessEqual(np.linalg.norm(got['mean_velocity_m_s']-mode),
                                     2e-11*np.linalg.norm(mode))
                scale = self.dt*self.p['rho0']*np.linalg.norm(mode)/self.p['dx']
                self.assertLessEqual(np.linalg.norm(got['density_kg_m3']), 2e-11*scale)

    def test_15_original_actuator_and_release_boundary(self):
        q = self.pa.extract_particles(np.arange(27, dtype=np.uint32)); q.z0[-1] = .0098
        b, c = self.b.copy(), self.c.copy(); b[-1, 2] = 0.; c[-1, 2] = 0.
        before = particle_hash(q)
        pc = StageEventPreconditioner(q, self.p, self.dt, 0., b, c, {})
        self.assertEqual(np.count_nonzero(~pc.free), 1)
        mapping = copied(self.candidate); mapping['mean_velocity_m_s'][-1, 2] = -.004
        got = pc.apply(self.candidate, mapping)
        self.assertEqual(got['correction']['mean_velocity_m_s'][-1, 2], -.004)
        self.assertEqual(got['correction']['position_m'][-1, 2], 0.)
        self.assertEqual(got['correction']['velocity_m_s'][-1, 2], 0.)
        self.assertEqual(before, particle_hash(q))
        with self.assertRaisesRegex(ValueError, 'release must lie'):
            StageEventPreconditioner(q, self.p, .2, .9, b, c, {})


if __name__ == '__main__':
    unittest.main(verbosity=2)

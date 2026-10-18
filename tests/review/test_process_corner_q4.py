"""Q4 independent spring/continuum references; constants are not solver outputs."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import probe_process_corner
from probe_process_corner import CornerBlock, resolve_active_set, reduced_correction, check_work
from probe_process_plastic import CandidateFailure


def accepted_snapshot(block):
    arrays = {name: value.copy() for name, value in vars(block).items()
              if isinstance(value, np.ndarray)}
    arrays.update({'pa.'+name: prop.get_npy_array().copy()
                   for name, prop in block.pa.properties.items()})
    arrays.update({'mesh.'+name: value.copy() for name, value in vars(block.mesh).items()
                   if isinstance(value, np.ndarray)})
    names = ('time', 'count', 'wall_work', 'contact_work', 'support_work',
             'constraint_work', 'internal_work', 'accepted_floor_count',
             'last_wall_increment', 'last_contact_increment', 'last_residual',
             'last_force_defect', 'last_torque_defect', 'last_iterations',
             'last_event', 'events', 'last_trace')
    return arrays, {name: copy.deepcopy(getattr(block, name)) for name in names
                    if hasattr(block, name)}


def spring_branch(K, walls, active):
    """Test fixture: solve supplied linear body with the ACTUAL reduced solve."""
    dofs = np.arange(len(walls))
    unknown = np.ones(len(K), dtype=bool)
    unknown[dofs[active]] = False
    u = np.zeros(len(K)); u[dofs[active]] = walls[active]
    correction, _ = reduced_correction(K, -K@u, unknown)
    u += correction
    force = -K@u
    normal = np.zeros(len(walls)); normal[active] = force[dofs[active]]
    return dict(displacement=u, normal_force=normal, gap=walls-u[dofs], active=active.copy())


class CornerQ4Tests(unittest.TestCase):
    def assert_unchanged(self, block, before):
        after = accepted_snapshot(block)
        self.assertEqual(before[0].keys(), after[0].keys())
        for name in before[0]:
            np.testing.assert_array_equal(before[0][name], after[0][name], err_msg=name)
        self.assertEqual(before[1], after[1])

    def test_01_two_dof_active_sets_have_literal_coupled_reference(self):
        K = np.array([[2., -1.], [-1., 2.]])
        limits = dict(maximum_active_sets=8, contact_force_tolerance_N=1e-12, gap_tolerance_m=1e-12)
        cases = [((-1., -1.), (-1., -1.), (1., 1.), (True, True)),
                 ((0., -1.), (-.5, -1.), (0., 1.5), (False, True)),
                 ((-1., -.25), (-1., -.5), (1.5, 0.), (True, False)),
                 ((-1., -.75), (-1., -.75), (1.25, .5), (True, True))]
        for walls, expected_u, expected_normal, expected_active in cases:
            for initial in ((False, False), (True, True)):
                with self.subTest(walls=walls, initial=initial):
                    d = np.array(walls)
                    c = resolve_active_set(lambda active: spring_branch(K, d, active), initial, limits)
                    np.testing.assert_allclose(c['displacement'], expected_u, atol=1e-14, rtol=0.)
                    np.testing.assert_allclose(c['normal_force'], expected_normal, atol=1e-14, rtol=0.)
                    np.testing.assert_array_equal(c['active'], expected_active)
                    self.assertTrue((c['gap'] >= -1e-14).all())

    def test_02_shared_free_dof_couples_both_walls_and_release(self):
        # z=(x+y)/3; both fixed -1 gives z=-2/3 and lambda=(4/3,4/3).
        # With Y open: 2*y-z=0, 3*z-x-y=0 => y=x/5, z=2*x/5.
        K = np.array([[2., 0., -1.], [0., 2., -1.], [-1., -1., 3.]])
        limits = dict(maximum_active_sets=8, contact_force_tolerance_N=1e-12, gap_tolerance_m=1e-12)
        for walls, expected_u, expected_normal in [
                ((-1., -1.), (-1., -1., -2/3), (4/3, 4/3)),
                ((-1., 0.), (-1., -1/5, -2/5), (8/5, 0.))]:
            with self.subTest(walls=walls):
                c = resolve_active_set(lambda active: spring_branch(K, np.array(walls), active),
                                       (True, True), limits)
                np.testing.assert_allclose(c['displacement'], expected_u, atol=1e-14, rtol=0.)
                np.testing.assert_allclose(c['normal_force'], expected_normal, atol=1e-14, rtol=0.)
                np.testing.assert_allclose(-K@c['displacement']-np.r_[c['normal_force'], 0.],
                                           0., atol=1e-14, rtol=0.)

    def test_03_inventory_preserves_shared_edge_without_duplicate_dofs_or_mass(self):
        b = CornerBlock()
        self.assertEqual(len(b.X), 27); self.assertEqual(int(b.free.sum()), 54)
        self.assertEqual(len(b.nodes), 18); self.assertEqual(len(np.unique(b.nodes)), 15)
        _, counts = np.unique(b.nodes, return_counts=True)
        self.assertEqual(int(np.sum(counts == 2)), 3)
        self.assertEqual(len(np.unique(b.contact_dofs)), 18)
        both_unknown = b.free.copy(); both_unknown[b.contact_dofs] = False
        self.assertEqual(int(both_unknown.sum()), 36)
        self.assertTrue(b.free[b.contact_dofs].all())
        np.testing.assert_array_equal(b.fixed, b.X == 0.)
        self.assertAlmostEqual(float(b.pa.m.sum()), .0027, delta=1e-15)
        identities = list(zip(b.pa.particle_id[b.nodes].tolist(), b.faces.tolist()))
        self.assertEqual(len(set(identities)), 18)
        for axis in range(3): self.assertEqual(len(np.unique(b.X[:, axis])), 3)

    def test_04_trial_tangent_and_active_set_are_pure(self):
        b = CornerBlock()
        before = accepted_snapshot(b)
        c = b.solve_candidate(.0125)
        b.tangent(c)
        b._equilibrate(.0125, np.ones(18, dtype=bool))
        self.assert_unchanged(b, before)
        self.assertTrue(c['active'].all())

    def test_05_per_wall_work_cancellation_must_not_pass(self):
        zero = np.zeros(2)
        # Sum of defects is zero, but each wall violates work conjugacy.
        with self.assertRaisesRegex(ValueError, 'per-wall'):
            check_work(np.array([.01, -.01]), zero, zero, zero, 0., 0., 1e-6)
        with self.assertRaisesRegex(ValueError, 'per-wall'):
            check_work(zero, zero, np.array([.01, -.01]), zero, 0., 0., 1e-6)
        check_work(np.array([.1, -.03]), np.array([.1, -.03]), zero, zero, 0., .07, 1e-6)
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            check_work(np.array([np.nan, 0.]), zero, zero, zero, 0., 0., 1e-6)

    def test_06_stale_and_mutated_candidate_cannot_change_accepted_state(self):
        b = CornerBlock()
        c = b.solve_candidate(.0125)
        corrupt = copy.deepcopy(c)
        center = np.flatnonzero(np.all(b.X == .005, axis=1))[0]
        corrupt['displacement'][center, 0] += 1e-7
        before = accepted_snapshot(b)
        with self.assertRaisesRegex(CandidateFailure, 'equilibrium/complementarity'):
            b.commit(corrupt)
        self.assert_unchanged(b, before)
        b.commit(c)
        before = accepted_snapshot(b)
        with self.assertRaisesRegex(ValueError, 'stale'): b.commit(c)
        self.assert_unchanged(b, before)

    def test_07_candidate_ledgers_and_face_forces_recomputed(self):
        b = CornerBlock()
        c = b.solve_candidate(.0125)
        c['ep'][:] = 123.; c['plastic_increment_J'][:] = 456.
        c['normal_force'][:] = -789.; c['face_force'][:] = 987.
        b.commit(c)
        np.testing.assert_array_equal(b.pa.ep, 0.)
        np.testing.assert_array_equal(b.pa.plastic_work, 0.)
        np.testing.assert_allclose(b.face_force.sum(axis=1), [[-3500., 0., 0.], [0., -3500., 0.]],
                                   atol=1e-4, rtol=0.)
        edge = np.flatnonzero((b.X[:, 0] == .01) & (b.X[:, 1] == .01))
        self.assertTrue((b.contact_force[edge, 0] < 0).all())
        self.assertTrue((b.contact_force[edge, 1] < 0).all())
        np.testing.assert_array_equal(b.contact_force, b.face_force.sum(axis=0))

    def test_08_missing_release_metadata_rejected_before_any_commit(self):
        b = CornerBlock()
        c = b.solve_candidate(.0125); c['event_kind'] = 'RELEASE_LIMIT'
        c.pop('event_bracket_s', None)
        before = accepted_snapshot(b)
        with self.assertRaisesRegex(ValueError, 'invalid release event bracket'):
            b.commit(c)
        self.assert_unchanged(b, before)

    def test_09_unconverged_branch_rejected_without_state_changes(self):
        b = CornerBlock(); before = accepted_snapshot(b)
        with self.assertRaisesRegex(CandidateFailure, 'did not converge'):
            b.solve_candidate(.0125, maximum_iterations=0)
        self.assert_unchanged(b, before)

    def test_10_program_knots_reject_crossing_but_accept_float_noise(self):
        b = CornerBlock('elastic', 'staggered')
        b.advance(.1); before = accepted_snapshot(b)
        with self.assertRaisesRegex(ValueError, 'program corner'): b.advance(.1625)
        self.assert_unchanged(b, before)
        # linspace can represent the intended .15 exactly one ULP above its literal.
        b.advance(float(np.nextafter(.15, np.inf)))
        np.testing.assert_allclose(b.wall_offsets(b.time), [-20e-6, -40e-6], atol=1e-19, rtol=0.)
        b.advance(.1625)
        self.assertEqual(b.time, .1625)

    def test_11_passive_labels_really_changed_do_not_change_mechanics(self):
        a, b = CornerBlock(), CornerBlock()
        b.pa.tracer_type[:] = 71+np.arange(len(b.X)) % 5
        self.assertFalse(np.array_equal(a.pa.tracer_type, b.pa.tracer_type))
        a.advance(.0125); b.advance(.0125)
        for name in ('displacement', 'stress', 'normal_force', 'face_force',
                     'wall_work_by_face', 'material_work_by_face'):
            np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
        np.testing.assert_array_equal(a.pa.ep, b.pa.ep)
        np.testing.assert_array_equal(a.pa.plastic_work, b.pa.plastic_work)

    def test_12_j2_staggered_shared_state_partial_release_and_free_residual_shape(self):
        b = CornerBlock('j2', 'staggered')
        for t in np.linspace(0., .1, 9)[1:]: b.advance(float(t))
        # Independent biaxial J2: n=(-1/2,-1/2,1), ep=2*(.004-Y/E).
        np.testing.assert_allclose(b.pa.ep, .005428571428571429, atol=1e-10, rtol=0.)
        np.testing.assert_allclose(b.face_force.sum(axis=1), [[-9000., 0., 0.], [0., -9000., 0.]],
                                   atol=1e-4, rtol=0.)
        np.testing.assert_allclose(b.displacement[b.X[:, 2] == .01, 2], 54.28571428571429e-6,
                                   atol=2e-10, rtol=0.)
        self.assertAlmostEqual(float(b.pa.plastic_work.sum()), .4885714285714286, delta=1e-7)
        np.testing.assert_allclose(b.wall_work_by_face, [.3, .3], atol=1e-7, rtol=0.)
        ep, plastic = b.pa.ep.copy(), b.pa.plastic_work.copy()
        b.advance(.1125); b.advance(.125)
        before = accepted_snapshot(b)
        unsplit = b.solve_candidate(.1375)
        self.assertFalse(unsplit['active'][b.faces == 0].any())
        self.assertTrue(unsplit['active'][b.faces == 1].all())
        with self.assertRaisesRegex(CandidateFailure, 'per-wall contact work'):
            b.commit(unsplit)
        self.assert_unchanged(b, before)
        b.advance(.1375)
        self.assertFalse(b.active[b.faces == 0].any())
        self.assertTrue(b.active[b.faces == 1].all())
        np.testing.assert_allclose(b.normal_force[b.faces == 0], 0., atol=0., rtol=0.)
        self.assertAlmostEqual(float(b.normal_force[b.faces == 1].sum()), 9000., delta=1e-4)
        x_after_release = b.displacement[:, 0].copy()
        b.advance(.15)
        np.testing.assert_allclose(b.displacement[:, 0], x_after_release, atol=1e-12, rtol=0.)
        for t in (.1625, .175, .1875, .2): b.advance(t)
        self.assertEqual(len(b.events), 2)
        for event, reference, feature in zip(b.events, [.13214285714285715, .16607142857142858], [1, 2]):
            self.assertAlmostEqual(event['physical_time_s'], reference, delta=1e-8)
            self.assertEqual(event['near_zero_feature_ids'], [feature])
        self.assertFalse(b.active.any())
        np.testing.assert_array_equal(b.normal_force, 0.)
        np.testing.assert_allclose(b.pa.ep, ep, atol=1e-12, rtol=0.)
        np.testing.assert_allclose(b.pa.plastic_work, plastic, atol=1e-10, rtol=0.)
        for axis in (0, 1):
            np.testing.assert_allclose(b.displacement[b.X[:, axis] == .01, axis], -27.142857142857144e-6,
                                       atol=2e-10, rtol=0.)
        np.testing.assert_allclose(b.stress, 0., atol=1., rtol=0.)
        np.testing.assert_allclose(b.wall_work_by_face, [.24214285714285715]*2, atol=1e-7, rtol=0.)
        np.testing.assert_allclose(b.wall_work_by_face, b.material_work_by_face, atol=1e-6, rtol=0.)
        self.assertAlmostEqual(b.internal_work, .4842857142857143, delta=1e-7)
        self.assertAlmostEqual(b.wall_work-float(b.pa.plastic_work.sum()), -.0042857142857143, delta=1e-7)

    def test_13_configuration_drift_and_production_scope_rejected(self):
        original = json.loads(probe_process_corner.CONFIG.read_text(encoding='utf-8'))
        cases = [('staggered_y_return_start', .14), ('wall_minimum', -50e-6),
                 ('production_authorized', True), ('wall_axes', [0, 0])]
        for field, value in cases:
            with self.subTest(field=field):
                cfg = copy.deepcopy(original)
                if field in cfg: cfg[field] = value
                else: cfg['parameters'][field]['value'] = value
                with tempfile.TemporaryDirectory(prefix='q4_config_rejection_') as folder:
                    path = Path(folder)/'INVALID_TEST_ONLY.json'
                    path.write_text(json.dumps(cfg), encoding='utf-8')
                    with patch.object(probe_process_corner, 'CONFIG', path):
                        with self.assertRaisesRegex(ValueError, 'declared'):
                            CornerBlock()


if __name__ == '__main__': unittest.main(verbosity=2)

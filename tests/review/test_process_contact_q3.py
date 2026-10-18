"""Independent Q3 unilateral references; no tested function builds expected values."""
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
from probe_process_contact import ContactBlock, reduced_correction
import probe_process_contact
from probe_process_plastic import CandidateFailure


def accepted_snapshot(block):
    """All particle/accepted NumPy arrays, plus accepted scalar/event ledgers."""
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
    ledger = {name: copy.deepcopy(getattr(block, name)) for name in names
              if hasattr(block, name)}
    return arrays, ledger


class ContactQ3Tests(unittest.TestCase):
    def assert_unchanged(self, block, before):
        after = accepted_snapshot(block)
        self.assertEqual(before[0].keys(), after[0].keys())
        for name in before[0]:
            np.testing.assert_array_equal(before[0][name], after[0][name], err_msg=name)
        self.assertEqual(before[1], after[1])

    def test_01_coupled_two_dof_reaction_is_not_independent_diagonal(self):
        # Independent mechanics: 2*u2-u1=0 => u2=-1/2; lambda=-f1=3/2.
        K = np.array([[2., -1.], [-1., 2.]])
        u = np.array([-1., 0.])
        correction, condition = reduced_correction(K, -K@u, [False, True])
        np.testing.assert_array_equal(correction, [0., -.5])
        u += correction
        internal_force = -K@u
        self.assertEqual(float(internal_force[1]), 0.)
        self.assertEqual(float(internal_force[0]), 1.5)
        self.assertNotEqual(float(internal_force[0]), 2.)
        self.assertEqual(condition, 1.)

    def test_02_condition_applies_to_solved_subblock_without_pinv(self):
        # The complete tangent has a constrained mechanism, the solved block does not.
        K = np.diag([0., 2.])
        correction, condition = reduced_correction(K, np.array([0., 4.]), [False, True])
        np.testing.assert_array_equal(correction, [0., 2.])
        self.assertEqual(condition, 1.)
        with self.assertRaisesRegex(ValueError, 'singular/ill-conditioned'):
            reduced_correction(K, np.array([0., 4.]), [True, True])

    def test_03_three_dimensions_and_rollers_only(self):
        b = ContactBlock()
        self.assertEqual(len(b.X), 27)
        self.assertEqual(int(b.free.sum()), 54)
        self.assertEqual(len(b.nodes), 9)
        np.testing.assert_array_equal(b.fixed, b.X == 0.)
        self.assertTrue(b.free[b.contact_dofs].all())
        for axis in range(3):
            self.assertEqual(len(np.unique(b.X[:, axis])), 3)
        self.assertAlmostEqual(float(b.pa.m.sum()), .0027, delta=1e-15)
        self.assertFalse(b.active.any())

    def test_04_trial_tangent_and_active_set_do_not_mutate_accepted_state(self):
        b = ContactBlock()
        before = accepted_snapshot(b)
        c = b.solve_candidate(.0125)
        b.tangent(c)
        b._equilibrate(.0125, np.ones(len(b.nodes), dtype=bool))
        self.assert_unchanged(b, before)
        self.assertGreater(float(c['normal_force'].sum()), 0.)

    def test_05_stale_candidate_rejected(self):
        b = ContactBlock()
        c = b.solve_candidate(.0125)
        b.commit(c)
        before = accepted_snapshot(b)
        with self.assertRaisesRegex(ValueError, 'stale'):
            b.commit(c)
        self.assert_unchanged(b, before)

    def test_06_mutated_displacement_rejected_without_state_changes(self):
        b = ContactBlock()
        c = b.solve_candidate(.0125)
        center = np.flatnonzero(np.all(b.X == .005, axis=1))[0]
        c['displacement'][center, 0] += 1e-7
        before = accepted_snapshot(b)
        with self.assertRaisesRegex(CandidateFailure, 'equilibrium/complementarity'):
            b.commit(c)
        self.assert_unchanged(b, before)

    def test_07_unconverged_branch_rejected_without_state_changes(self):
        b = ContactBlock()
        before = accepted_snapshot(b)
        with self.assertRaisesRegex(CandidateFailure, 'did not converge'):
            b.solve_candidate(.0125, maximum_iterations=0)
        self.assert_unchanged(b, before)

    def test_08_candidate_material_ledger_is_recomputed_not_trusted(self):
        b = ContactBlock()
        c = b.solve_candidate(.0125)
        c['ep'][:] = 123.
        c['plastic_increment_J'][:] = 456.
        c['normal_force'][:] = -789.
        b.commit(c)
        np.testing.assert_array_equal(b.pa.ep, 0.)
        np.testing.assert_array_equal(b.pa.plastic_work, 0.)
        # Actual first increment: E*(.0005)*A = 3500 N compression.
        self.assertAlmostEqual(float(b.normal_force.sum()), 3500., delta=1e-4)

    def test_09_passive_labels_do_not_change_mechanics(self):
        a, b = ContactBlock(), ContactBlock()
        b.pa.tracer_type[:] = 71+np.arange(len(b.X)) % 5
        self.assertFalse(np.array_equal(a.pa.tracer_type, b.pa.tracer_type))
        a.advance(.0125); b.advance(.0125)
        for left, right in ((a.displacement, b.displacement), (a.stress, b.stress),
                            (a.normal_force, b.normal_force), (a.pa.ep, b.pa.ep),
                            (a.pa.plastic_work, b.pa.plastic_work)):
            np.testing.assert_array_equal(left, right)

    def test_10_elastic_load_return_has_no_plasticity_or_residual(self):
        b = ContactBlock('elastic')
        b.advance(.1)
        self.assertAlmostEqual(float(b.normal_force.sum()), 28000., delta=1e-4)
        self.assertAlmostEqual(b.wall_work, .56, delta=1e-7)
        self.assertAlmostEqual(b.contact_work, .56, delta=1e-7)
        np.testing.assert_array_equal(b.pa.ep, 0.)
        b.advance(.2)
        np.testing.assert_allclose(b.displacement, 0., atol=1e-12, rtol=0.)
        np.testing.assert_allclose(b.stress, 0., atol=.1, rtol=0.)
        np.testing.assert_array_equal(b.pa.ep, 0.)
        np.testing.assert_array_equal(b.pa.plastic_work, 0.)
        self.assertAlmostEqual(b.wall_work, 0., delta=1e-7)
        self.assertAlmostEqual(b.contact_work, 0., delta=1e-7)
        # Wall only returns to initial position: positive-gap elastic extension NOT_RUN.

    def test_11_missing_release_metadata_rejected_atomically(self):
        b = ContactBlock()
        c = b.solve_candidate(.05)
        c['event_kind'] = 'RELEASE_LIMIT'
        c.pop('event_bracket_s', None)
        before = accepted_snapshot(b)
        with self.assertRaises((KeyError, ValueError, CandidateFailure)):
            b.commit(c)
        self.assert_unchanged(b, before)

    def test_12_j2_release_event_free_shape_and_work_independent_reference(self):
        b = ContactBlock('j2')
        for t in np.linspace(0., .1, 9)[1:]:
            b.advance(float(t))
        # Literal independent continuum values from PROCESS_CONTACT_Q3_REFERENCES.
        self.assertAlmostEqual(float(b.normal_force.sum()), 9000., delta=1e-4)
        np.testing.assert_allclose(b.pa.ep, .002714285714285714, atol=1e-10, rtol=0.)
        self.assertAlmostEqual(float(b.pa.plastic_work.sum()), .2442857142857143, delta=1e-7)
        self.assertAlmostEqual(b.wall_work, .3000000000000000, delta=1e-7)
        b.advance(.1125); b.advance(.125)
        before = accepted_snapshot(b)
        unsplit = b.solve_candidate(.1375)
        self.assertFalse(unsplit['active'].any())
        # The endpoints are equilibrated, but wall travelled farther than the material.
        with self.assertRaisesRegex(CandidateFailure, 'work requires event split'):
            b.commit(unsplit)
        self.assert_unchanged(b, before)
        captured = []
        b.advance(.1375, lambda: captured.append((b.time, b.last_event)))
        self.assertEqual(len(b.events), 1)
        self.assertTrue(any(kind == 'RELEASE_LIMIT' for _, kind in captured))
        self.assertAlmostEqual(b.events[0]['physical_time_s'], .13214285714285715, delta=1e-8)
        self.assertFalse(b.active.any())
        free_displacement = b.displacement.copy()
        ep, plastic = b.pa.ep.copy(), b.pa.plastic_work.copy()
        for t in (.15, .1625, .175, .1875, .2):
            b.advance(t)
        np.testing.assert_allclose(b.displacement, free_displacement, atol=1e-12, rtol=0.)
        np.testing.assert_array_equal(b.pa.ep, ep)
        np.testing.assert_array_equal(b.pa.plastic_work, plastic)
        np.testing.assert_allclose(b.displacement[b.nodes, 0], -27.142857142857144e-6,
                                   atol=2e-10, rtol=0.)
        np.testing.assert_allclose(-b.displacement[b.nodes, 0], 27.142857142857144e-6,
                                   atol=2e-10, rtol=0.)
        np.testing.assert_allclose(b.stress, 0., atol=1., rtol=0.)
        np.testing.assert_array_equal(b.normal_force, 0.)
        self.assertAlmostEqual(b.wall_work, .2421428571428572, delta=1e-7)
        self.assertAlmostEqual(b.contact_work, .2421428571428572, delta=1e-7)
        self.assertAlmostEqual(b.internal_work, .2421428571428572, delta=1e-7)
        # Cross-yield trapezoid deficit remains visible, never added as heat.
        self.assertAlmostEqual(b.wall_work-float(b.pa.plastic_work.sum()),
                               -.0021428571428571, delta=1e-7)

    def test_13_configuration_drift_and_production_authorization_rejected(self):
        original = json.loads(probe_process_contact.CONFIG.read_text(encoding='utf-8'))
        for field, value in (('phase_duration', .2), ('temperature', 323.15),
                             ('production_authorized', True)):
            with self.subTest(field=field):
                cfg = copy.deepcopy(original)
                if field == 'production_authorized':
                    cfg[field] = value
                    message = 'production/normal scope mismatch'
                else:
                    cfg['parameters'][field]['value'] = value
                    message = 'declared/used parameter mismatch'
                with tempfile.TemporaryDirectory(prefix='q3_config_rejection_') as folder:
                    path = Path(folder)/'INVALID_TEST_ONLY.json'
                    path.write_text(json.dumps(cfg), encoding='utf-8')
                    with patch.object(probe_process_contact, 'CONFIG', path):
                        with self.assertRaisesRegex(ValueError, message):
                            ContactBlock()


if __name__ == '__main__':
    unittest.main(verbosity=2)

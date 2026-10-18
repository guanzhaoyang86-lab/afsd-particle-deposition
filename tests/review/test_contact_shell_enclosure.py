"""Analytic supplement fixtures; no AFSD equilibrium or material solve."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from batched_contact import BatchedChainRuleContact
from contact_shell_audit import CERTIFIED, FOUND, UNRESOLVED
from contact_shell_enclosure import audit_contact_shell_enclosure, CLASSIFICATION


def pose(theta=0., origin=(0., 0., 0.)):
    c, s = np.cos(theta), np.sin(theta)
    return dict(origin=np.asarray(origin, float), rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


class PlaneEnclosureTests(unittest.TestCase):
    def setUp(self):
        self.contact = BatchedChainRuleContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        self.p = pose()

    def test_parallel_near_wall_translation_certifies_without_more_geometry_queries(self):
        x0 = np.array([[.399999, .0, .5]])
        x1 = np.array([[.399999, .01, .5]])
        r = audit_contact_shell_enclosure(self.contact, x0, x1, self.p, self.p, maximum_queries=2)
        self.assertEqual(r['original_audit']['status'], UNRESOLVED)
        self.assertEqual(r['status'], CERTIFIED)
        self.assertEqual(r['geometry_queries'], 2)
        self.assertEqual(r['additional_finite_patch_queries'], 0)
        self.assertEqual(r['plane_path_evaluations'], 2)
        self.assertEqual(r['classification'], CLASSIFICATION)
        self.assertTrue(all(c['finite_patch_distance_lower_bound_m'] > .1 for c in r['planar_certificates']))

    def test_endpoint_chord_alone_would_miss_curved_wall_contact(self):
        x = np.array([[.405, 0., .5]])
        p0, p1 = pose(-np.pi/4), pose(np.pi/4)
        # Endpoint plane distances >r on the same side; midpoint is inside the
        # contact shell. Curvature subtraction must prevent endpoint-only pass.
        r = audit_contact_shell_enclosure(self.contact, x, x, p0, p1,
            maximum_queries=2, maximum_plane_evaluations=8)
        self.assertEqual(r['status'], UNRESOLVED)
        self.assertEqual(r['geometry_queries'], 2)
        self.assertTrue(r['planar_unresolved_intervals'])
        full = audit_contact_shell_enclosure(self.contact, x, x, p0, p1)
        self.assertEqual(full['status'], FOUND)
        self.assertEqual(full['plane_path_evaluations'], 0)

    def test_same_rigid_endpoints_use_world_chord_path(self):
        x0 = np.array([[.6, 0., -.04]])
        p1 = pose(np.pi/2)
        x1 = x0@p1['rotation'].T
        r = audit_contact_shell_enclosure(self.contact, x0, x1, self.p, p1)
        self.assertEqual(r['status'], FOUND)
        self.assertEqual(r['plane_path_evaluations'], 0)
        self.assertTrue(any(e['feature'] == 'BORE_X_PLUS' for e in r['hidden_contacts']))

    def test_invalid_solid_is_never_upgraded(self):
        x = np.array([[.55, 0., .5]])
        r = audit_contact_shell_enclosure(self.contact, x, x, pose(-np.pi/4), pose(np.pi/4))
        self.assertEqual(r['status'], UNRESOLVED)
        self.assertEqual(r['failure_kind'], 'INVALID_SOLID_PATH')
        self.assertEqual(r['plane_path_evaluations'], 0)

    def test_outer_cylinder_remainder_is_not_a_plane(self):
        x0 = np.array([[2.100001, 0., .5]])
        x1 = np.array([[2.100001, .01, .5]])
        r = audit_contact_shell_enclosure(self.contact, x0, x1, self.p, self.p, maximum_queries=2)
        self.assertEqual(r['status'], UNRESOLVED)
        self.assertTrue(any(c['feature'] == 'OUTER_CYLINDER' for c in r['planar_unresolved_intervals']))

    def test_base_stays_with_original_certificate(self):
        self.contact.base_z = -1.
        x0 = np.array([[0., 0., -.8]])
        x1 = np.array([[.1, 0., -.7]])
        r = audit_contact_shell_enclosure(self.contact, x0, x1, self.p, self.p)
        self.assertEqual(r['status'], CERTIFIED)
        self.assertEqual(r['original_audit']['status'], CERTIFIED)
        self.assertEqual(r['plane_path_evaluations'], 0)

    def test_roundoff_margin_prevents_ambiguous_tangency_upgrade(self):
        x0 = np.array([[.4, 0., .5]])
        x1 = np.array([[.4, .01, .5]])
        # Binary .5-.4 is slightly less than .1: use a radius equal to its
        # actually represented gap so the original endpoint test is inactive.
        self.contact.particle_radius = float(.5-.4)
        r = audit_contact_shell_enclosure(self.contact, x0, x1, self.p, self.p,
            maximum_queries=2, maximum_depth=0)
        self.assertEqual(r['status'], UNRESOLVED)
        self.assertTrue(all(c['roundoff_margin_m'] > 0. for c in r['planar_unresolved_intervals']))

    def test_accepted_pose_roundoff_defect_is_retained_and_inputs_pure(self):
        p = pose(); p['rotation'][0, 0] += 1e-13
        x0 = np.array([[.39999, 0., .5]])
        x1 = np.array([[.39999, .01, .5]])
        before = (x0.copy(), x1.copy(), p['rotation'].copy())
        r = audit_contact_shell_enclosure(self.contact, x0, x1, p, p, maximum_queries=2)
        self.assertEqual(r['status'], CERTIFIED)
        self.assertGreater(r['rotation_defect'], 0.)
        self.assertGreater(r['maximum_planar_roundoff_margin_m'], 2048*np.finfo(float).eps*2.)
        for actual, old in zip((x0, x1, p['rotation']), before):
            np.testing.assert_array_equal(actual, old)

    def test_invalid_analytic_budget_rejected(self):
        x = np.array([[0., 0., .5]])
        for budget in (1, 257, True):
            with self.assertRaises(ValueError):
                audit_contact_shell_enclosure(self.contact, x, x, self.p, self.p, maximum_plane_evaluations=budget)


if __name__ == '__main__':
    unittest.main(verbosity=2)

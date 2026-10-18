"""Generic finite-patch lower-bound fixtures, not AFSD process validation."""
from pathlib import Path
import math
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from batched_contact import BatchedChainRuleContact
from contact_shell_audit import CERTIFIED, FOUND, UNRESOLVED
from contact_shell_enclosure import audit_contact_shell_enclosure
from contact_shell_finite_enclosure import audit_contact_shell_finite_enclosure, finite_patch_box_lower_bound


def pose(theta=0., origin=(0., 0., 0.)):
    c, s = np.cos(theta), np.sin(theta)
    return dict(origin=np.asarray(origin, float), rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


class FiniteEnclosureTests(unittest.TestCase):
    def setUp(self):
        self.contact = BatchedChainRuleContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        self.p = pose()

    def bound(self, feature, lo, hi=None):
        return finite_patch_box_lower_bound(self.contact, feature, lo, lo if hi is None else hi)

    def near_lower(self, actual, expected):
        self.assertLessEqual(actual, expected)
        self.assertAlmostEqual(actual, expected, delta=2e-12)

    def test_bore_rectangle_includes_both_finite_edge_gaps(self):
        self.near_lower(self.bound(3, [.4, .6, -.1]), math.sqrt(.03))
        self.near_lower(self.bound(3, [.4, .6, -.1], [.45, .7, -.1]), .15)
        self.near_lower(self.bound(2, [-.4, .6, 1.1]), math.sqrt(.03))
        self.near_lower(self.bound(4, [.6, -.4, -.1]), math.sqrt(.03))
        self.near_lower(self.bound(5, [.6, .4, 1.1]), math.sqrt(.03))

    def test_caps_retain_square_opening_and_circle_constraint(self):
        self.near_lower(self.bound(0, [.35, .1, -.21]), math.hypot(.15, .21))
        self.near_lower(self.bound(1, [.35, .1, 1.21]), math.hypot(.15, .21))
        self.near_lower(self.bound(0, [3., 0., -.2]), math.hypot(1., .2))
        self.near_lower(self.bound(1, [3., 0., 1.2]), math.hypot(1., .2))

    def test_lip_clear_path_where_infinite_plane_certificate_cannot_pass(self):
        self.contact.particle_radius = .25
        first = np.array([[.35, .1, -.21]])
        last = np.array([[.35, .12, -.21]])
        old = audit_contact_shell_enclosure(self.contact, first, last, self.p, self.p,
            maximum_queries=2, maximum_plane_evaluations=4)
        new = audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, self.p, maximum_queries=2)
        self.assertEqual(old['status'], UNRESOLVED)
        self.assertEqual(new['original_audit']['status'], UNRESOLVED)
        self.assertEqual(new['status'], CERTIFIED)
        self.assertEqual(new['geometry_queries'], 2)
        self.assertEqual(new['additional_finite_patch_queries'], 0)
        self.assertTrue(any(c['feature'] == 'SHOULDER' for c in new['finite_patch_certificates']))

    def test_outer_radius_range_uses_entire_box_and_finite_height(self):
        self.assertEqual(self.bound(6, [-2.2, 0., .5], [2.2, 0., .5]), 0.)
        self.near_lower(self.bound(6, [2.2, 0., 1.2]), math.hypot(.2, .2))
        self.near_lower(self.bound(6, [0., 0., .5]), 2.)
        first = np.array([[2.1001, 0., .5]])
        last = np.array([[2.1001, .01, .5]])
        result = audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, self.p, maximum_queries=2)
        self.assertEqual(result['original_audit']['status'], UNRESOLVED)
        self.assertEqual(result['status'], CERTIFIED)
        self.assertTrue(any(c['feature'] == 'OUTER_CYLINDER' for c in result['finite_patch_certificates']))

    def test_hidden_outer_shell_chord_is_not_certified_from_endpoint_radii(self):
        radius, angle = 2.12, .2
        first = np.array([[radius*np.cos(angle), -radius*np.sin(angle), .5]])
        last = first.copy(); last[:, 1] *= -1
        limited = audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, self.p,
            maximum_queries=2, maximum_path_evaluations=8)
        self.assertEqual(limited['status'], UNRESOLVED)
        full = audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, self.p)
        self.assertEqual(full['status'], FOUND)
        self.assertTrue(any(e['feature'] == 'OUTER_CYLINDER' for e in full['hidden_contacts']))

    def test_rotation_curve_hidden_wall_contact_is_not_endpoint_only_certified(self):
        x = np.array([[.405, 0., .5]])
        p0, p1 = pose(-np.pi/4), pose(np.pi/4)
        limited = audit_contact_shell_finite_enclosure(self.contact, x, x, p0, p1,
            maximum_queries=2, maximum_path_evaluations=8)
        self.assertEqual(limited['status'], UNRESOLVED)
        full = audit_contact_shell_finite_enclosure(self.contact, x, x, p0, p1)
        self.assertEqual(full['status'], FOUND)
        self.assertEqual(full['analytic_path_evaluations'], 0)

    def test_rigid_motion_endpoints_still_use_world_chord(self):
        first = np.array([[.6, 0., -.04]])
        p1 = pose(np.pi/2)
        last = first@p1['rotation'].T
        result = audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, p1)
        self.assertEqual(result['status'], FOUND)
        self.assertTrue(any(e['feature'] == 'BORE_X_PLUS' for e in result['hidden_contacts']))

    def test_invalid_solid_is_never_upgraded(self):
        x = np.array([[.55, 0., .5]])
        result = audit_contact_shell_finite_enclosure(self.contact, x, x, pose(-np.pi/4), pose(np.pi/4))
        self.assertEqual(result['status'], UNRESOLVED)
        self.assertEqual(result['failure_kind'], 'INVALID_SOLID_PATH')
        self.assertEqual(result['analytic_path_evaluations'], 0)

    def test_base_stays_in_original_world_coordinate_affine_audit(self):
        self.contact.base_z = -1.
        first = np.array([[0., 0., -.8]])
        last = np.array([[.1, 0., -.7]])
        result = audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, pose(.1))
        self.assertEqual(result['status'], CERTIFIED)
        self.assertFalse(any(c['feature'] == 'BASE' for c in result['finite_patch_certificates']))
        with self.assertRaises(ValueError):
            self.bound(7, [0., 0., -.8])

    def test_scale_consistency_and_pure_inputs(self):
        for scale in (1e-3, 1., 1e3):
            contact = BatchedChainRuleContact(2*scale, scale, scale, .25*scale, 100., 200., .3, base_z=None)
            first = np.array([[.35, .1, -.21]])*scale
            last = np.array([[.35, .12, -.21]])*scale
            before = first.copy(), last.copy()
            result = audit_contact_shell_finite_enclosure(contact, first, last, self.p, self.p, maximum_queries=2)
            self.assertEqual(result['status'], CERTIFIED)
            np.testing.assert_array_equal(first, before[0]); np.testing.assert_array_equal(last, before[1])
            for entry in result['finite_patch_certificates']:
                self.assertGreater(entry['finite_patch_distance_lower_bound_m'], .25*scale)

    def test_solid_box_distance_bounds_never_exceed_sampled_patch_distances(self):
        rng = np.random.default_rng(18641)
        # Sampling is a bound sanity check only; the production certificate is
        # derived from closed-set containment, not from these finite samples.
        for _ in range(60):
            centre = rng.uniform([-.7, -.7, -.4], [.7, .7, -.15])
            width = rng.uniform(0., .03, 3)
            lo, hi = centre-width, centre+width
            samples = rng.uniform(lo, hi, (30, 3))
            distances = self.contact._patches_batch(samples, reject_active_ties=False)[0]
            for feature in range(7):
                self.assertLessEqual(self.bound(feature, lo, hi), float(distances[:, feature].min()))

    def test_tangency_roundoff_is_not_upgraded_and_invalid_bounds_reject(self):
        self.contact.particle_radius = float(.5-.4)
        first = np.array([[.4, 0., .5]])
        last = np.array([[.4, .01, .5]])
        result = audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, self.p,
            maximum_queries=2, maximum_depth=0)
        self.assertEqual(result['status'], UNRESOLVED)
        for lower, upper in (([0, 0, 0], [-1, 1, 1]), ([np.nan, 0, 0], [1, 1, 1])):
            with self.assertRaises(ValueError):
                self.bound(0, lower, upper)
        with self.assertRaises(ValueError):
            audit_contact_shell_finite_enclosure(self.contact, first, last, self.p, self.p, maximum_path_evaluations=257)


if __name__ == '__main__':
    unittest.main(verbosity=2)

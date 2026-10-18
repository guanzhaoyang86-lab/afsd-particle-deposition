"""Manufactured geometry cases for event routing, not material validation."""
import unittest
import numpy as np
from visible_support_events_20261007 import audit_visibility_events, EVENT, UNRESOLVED, NO_OBSERVED_EVENT


class VisibilityEventTests(unittest.TestCase):
    def setUp(self):
        self.p = dict(tool_radius=3., bore_width=2., layer=0., tool_top=2.)
        self.pose = dict(origin=np.zeros(3), rotation=np.eye(3))
        self.i, self.j, self.ids = np.array([0]), np.array([1]), np.array([91, 13])

    def audit(self, x0, x1=None, pose1=None, **options):
        return audit_visibility_events(np.asarray(x0, float),
            np.asarray(x0 if x1 is None else x1, float), self.i, self.j,
            self.pose, self.pose if pose1 is None else pose1, self.p, self.ids, **options)

    def test_known_release_localized_and_never_accepted(self):
        # A vertical edge at x=1.2 moves into a unit half-width bore. The
        # strict geometric release is f=.5; arithmetic policy releases just after it.
        x0 = np.array([[1.2, 0., .5], [1.2, 0., 1.5]])
        x1 = x0-np.array([.4, 0., 0.])
        r = self.audit(x0, x1)
        self.assertEqual(r['status'], EVENT)
        self.assertTrue(r['requires_resolve'])
        self.assertFalse(r['acceptance_allowed'])
        b = r['transition_brackets'][0]
        self.assertAlmostEqual(.5, .5*(b['start_fraction']+b['end_fraction']), places=8)
        self.assertLessEqual(b['bracket_width_fraction'], r['fraction_tolerance'])
        self.assertEqual(b['direction'], 'BLOCKED_TO_VISIBLE')
        self.assertEqual(r['endpoint_changed_edges'][0]['start_strict_intersection_length_m'], 1.)
        self.assertEqual(r['edge_deletions'], 0)

    def test_reverse_is_loss_not_activation(self):
        r = self.audit([[.8, 0., .5], [.8, 0., 1.5]], [[1.2, 0., .5], [1.2, 0., 1.5]])
        self.assertEqual(r['transition_brackets'][0]['direction'], 'VISIBLE_TO_BLOCKED')

    def test_visible_endpoints_can_have_two_internal_changes(self):
        x = np.array([[1.1*np.sqrt(2), 1.1*np.sqrt(2), -1.], [0., 0., 1.]])
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        r = self.audit(x, pose1=dict(origin=np.zeros(3), rotation=rotation))
        self.assertEqual(r['status'], EVENT)
        self.assertEqual(len(r['endpoint_changed_edges']), 0)
        self.assertEqual(len(r['transition_brackets']), 2)

    def test_stationary_visible_is_not_path_pass(self):
        r = self.audit([[.2, 0., .5], [.8, 0., 1.5]])
        self.assertEqual(r['status'], NO_OBSERVED_EVENT)
        self.assertFalse(r['acceptance_allowed'])
        self.assertFalse(r['requires_resolve'])

    def test_budget_exhaustion_preserves_observed_bracket_and_fails_closed(self):
        r = self.audit([[1.2, 0., .5], [1.2, 0., 1.5]],
                       [[.8, 0., .5], [.8, 0., 1.5]], maximum_evaluations=17)
        self.assertEqual(r['status'], UNRESOLVED)
        self.assertEqual(r['geometry_evaluations'], 17)
        self.assertEqual(r['transition_brackets'][0]['status'], 'EVALUATION_BUDGET_EXHAUSTED')
        self.assertTrue(r['requires_resolve'])

    def test_collapsed_interior_is_unresolved(self):
        r = self.audit([[1.5, -4., 1.], [1.5, -3., 1.]],
                       [[1.5, 4., 1.], [1.5, 3., 1.]])
        self.assertEqual(r['status'], UNRESOLVED)
        self.assertTrue(any(row['fraction'] == .5 for row in r['unresolved_samples']))

    def test_large_scale_partition_ambiguity_is_unresolved(self):
        r = self.audit([[-1e8, 2.8, 1.], [1e8, 2.8, 1.]])
        self.assertEqual(r['status'], UNRESOLVED)
        self.assertTrue(r['requires_resolve'])

    def test_input_purity_and_observational_ids(self):
        x0 = np.array([[1.2, 0., .5], [1.2, 0., 1.5]])
        x1 = x0-np.array([.4, 0., 0.])
        old = [a.copy() for a in (x0, x1, self.i, self.j, self.ids)]
        a = self.audit(x0, x1)
        for before, after in zip(old, (x0, x1, self.i, self.j, self.ids)):
            np.testing.assert_array_equal(before, after)
        self.ids = np.array([5, 8])
        b = self.audit(x0, x1)
        self.assertEqual(a['query_log'], b['query_log'])
        self.assertNotEqual(a['input_sha256']['supplied_topology'], b['input_sha256']['supplied_topology'])

    def test_invalid_options_rejected(self):
        for options in ({'sample_intervals': 0}, {'maximum_evaluations': 3},
                        {'fraction_tolerance': 0.}, {'sample_intervals': True}):
            with self.assertRaises(ValueError):
                self.audit([[.2, 0., .5], [.8, 0., 1.5]], **options)


if __name__ == '__main__':
    unittest.main()

"""Pure geometry/domain tests; no AFSD force, heat, NNPS or state advances."""
import copy
import unittest
import numpy as np

import visible_support_numerics as legacy
import visible_support_path_tight as old_path
from visible_support import filter_visible_support
from audit_solid_occluded_neighbors import _partitions
import visible_support_remote_domain_20261007 as remote


class RemoteFiniteToolDomainTests(unittest.TestCase):
    def setUp(self):
        self.p = dict(tool_radius=2.9, bore_width=2., layer=0., tool_top=2.)
        self.pose = dict(origin=np.zeros(3), rotation=np.eye(3))
        self.i, self.j, self.ids = np.array([0]), np.array([1]), np.array([31, 17])

    def check(self, start, direction, **kwargs):
        args = dict(half_bore=1., zlo=0., zhi=2.)
        args.update(kwargs)
        return remote.check_finite_partition_numerical_domain(np.array(start, float),
            np.array(direction, float), 2.9, 1., **args)

    def endpoint(self, x):
        return remote.require_visible_support_numerical_domain(np.array(x, float),
            self.i, self.j, self.pose, self.p, self.ids)

    def path(self, x0, x1=None):
        return remote.audit_visible_support_path(np.array(x0, float),
            np.array(x0 if x1 is None else x1, float), self.i, self.j,
            self.pose, self.pose, self.p, self.ids)

    def test_long_bar_above_and_far_substrate_below_tool_are_proved_remote(self):
        for start, direction in (([[.2, 0., 10.]], [[.6, 0., 1.]]),
                                  ([[.2, 0., -10.]], [[.6, 0., -1.]])):
            old = legacy.check_partition_numerical_domain(np.array(start), np.array(direction), 2.9, 1.)
            self.assertEqual(old['status'], legacy.UNRESOLVED)
            current = self.check(start, direction)
            self.assertEqual(current['status'], legacy.PASS)
            self.assertEqual(current['remote_z_edge_indices'], [0])
            self.assertEqual(current['near_edge_indices'], [])
            self.assertIsNone(current['legacy_domain'])

    def test_radial_box_lower_bound_certifies_whole_segment(self):
        r = self.check([[10., -.5, .5]], [[0., 1., 1.]])
        self.assertEqual(r['status'], legacy.PASS)
        self.assertEqual(r['remote_radial_edge_indices'], [0])
        self.assertEqual(r['remote_z_edge_indices'], [])

    def test_far_endpoints_do_not_exempt_segment_crossing_tool(self):
        r = self.check([[-1e8, 2.8, 1.]], [[2e8, 0., 0.]])
        self.assertEqual(r['remote_edge_indices'], [])
        self.assertEqual(r['status'], legacy.UNRESOLVED)
        self.assertIn('COORDINATE_EDGE_OR_RADIUS_SCALE_RATIO_EXCEEDS_EIGHT', r['reasons'])

    def test_near_tangent_discriminant_uncertainty_is_unchanged(self):
        start, direction = np.array([[-2., 2.9, 1.]]), np.array([[4., 0., 0.]])
        r = self.check(start, direction)
        self.assertEqual(r['legacy_domain'], legacy.check_partition_numerical_domain(start, direction, 2.9, 1.))
        self.assertEqual(r['status'], legacy.UNRESOLVED)
        self.assertEqual(r['uncertain_edge_indices'], [0])

    def test_mixed_rows_map_legacy_uncertainty_back_to_original_order(self):
        r = self.check([[.2, 0., 10.], [-2., 2.9, 1.]], [[.6, 0., 1.], [4., 0., 0.]])
        self.assertEqual(r['remote_edge_indices'], [0])
        self.assertEqual(r['near_edge_indices'], [1])
        self.assertEqual(r['uncertain_edge_indices'], [1])
        self.assertEqual(r['legacy_domain']['uncertain_edge_indices'], [0])

    def test_actual_inflated_query_slab_is_used(self):
        r = self.check([[.2, 0., 10.]], [[.6, 0., 1.]], zlo=-20., zhi=20., half_bore=-3.)
        self.assertEqual(r['remote_edge_indices'], [])
        self.assertEqual(r['status'], legacy.UNRESOLVED)

    def test_actual_expanded_endpoint_query_preserves_grazing_block(self):
        x = np.array([[1.5, 0., 2.+1e-15], [1.6, 0., 2.+1e-15]])
        r = self.endpoint(x)
        self.assertEqual(r['proven_visible_raw_edge_indices'], [])
        self.assertGreater(r['expanded']['queried_zhi_m'], 2.)
        g = filter_visible_support(x, self.i, self.j, np.zeros((1, 3)), self.pose, self.p, self.ids)
        self.assertEqual(g['blocked_raw_edge_indices'].tolist(), [0])
        remote.require_remote_filter_consistency(r, g)

    def test_original_filter_retains_remote_edge_and_gradient_bitwise(self):
        x = np.array([[.2, 0., 10.], [.8, 0., 11.]])
        gradient = np.array([[1.234, -4.567, 8.9]])
        r = self.endpoint(x)
        g = filter_visible_support(x, self.i, self.j, gradient, self.pose, self.p, self.ids)
        remote.require_remote_filter_consistency(r, g)
        np.testing.assert_array_equal(g['i'], self.i)
        np.testing.assert_array_equal(g['j'], self.j)
        np.testing.assert_array_equal(g['g'], gradient)
        self.assertEqual(r['proven_visible_raw_edge_indices'], [0])
        self.assertEqual(r['edge_deletions'], 0)

    def test_remote_filter_disagreement_fails_instead_of_deleting_support(self):
        x = np.array([[.2, 0., 10.], [.8, 0., 11.]])
        r = self.endpoint(x)
        bad = dict(kept_raw_edge_indices=np.array([], int), blocked_raw_edge_indices=np.array([0]))
        with self.assertRaises(ArithmeticError):
            remote.require_remote_filter_consistency(r, bad)

    def test_remote_partition_postcondition_detects_false_solid_intersection(self):
        start, direction = np.array([[.2, 0., 10.]]), np.array([[.6, 0., 1.]])
        r = self.check(start, direction)
        result = _partitions(start, direction, 2.9, 1., 0., 2.)
        remote.require_remote_partition_consistency(r, result)
        corrupt = copy.deepcopy(result)
        corrupt[-1][0] = 1e-30
        with self.assertRaises(ArithmeticError):
            remote.require_remote_partition_consistency(r, corrupt)

    def test_stationary_remote_path_is_certified_without_dropping_edges(self):
        x = [[.2, 0., 10.], [.8, 0., 11.]]
        r = self.path(x)
        self.assertEqual(r['status'], old_path.CERTIFIED)
        self.assertEqual(r['edge_count'], 1)
        self.assertEqual(r['edge_deletions'], 0)
        self.assertEqual(r['analytic_batch_queries'], 4)
        self.assertEqual(r['classification'], remote.CLASSIFICATION)

    def test_actual_enclosure_can_reenter_near_domain_and_still_fail(self):
        r = self.path([[7., 0., -1.], [7., 0., -2.]], [[-7., 0., -1.], [-7., 0., -2.]])
        self.assertEqual(r['status'], old_path.UNRESOLVED)
        self.assertEqual(r['failure_kind'], 'PARTITION_NUMERICAL_DOMAIN_UNRESOLVED')
        self.assertEqual(r['numerical_domain_failures'][0]['kind'], 'enclosure')

    def test_strict_solid_path_failure_remains_failure(self):
        r = self.path([[1.2, 0., .5], [1.2, 0., 1.5]])
        self.assertEqual(r['status'], old_path.FOUND)
        self.assertGreater(r['found_edges'][0]['total_inside_length_m'], 0.)

    def test_units_rotation_translation_and_input_purity(self):
        for factor in (1e-4, 1., 1e3):
            c, s = np.cos(.61), np.sin(.61)
            rotation = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
            origin = factor*np.array([.2, -.1, .3])
            x = factor*np.array([[.2, 0., 10.], [.8, 0., 11.]])@rotation.T+origin
            p = dict(tool_radius=2.9*factor, bore_width=2.*factor, layer=0., tool_top=2.*factor)
            before = x.copy()
            r = remote.require_visible_support_numerical_domain(x, self.i, self.j,
                dict(origin=origin, rotation=rotation), p, self.ids)
            self.assertEqual(r['status'], legacy.PASS)
            self.assertEqual(r['proven_visible_raw_edge_indices'], [0])
            np.testing.assert_array_equal(x, before)

    def test_tiny_z_clearance_is_not_mistaken_for_a_remote_proof(self):
        r = self.check([[1.2, 0., 2.+1e-15]], [[.2, 0., 0.]])
        self.assertEqual(r['remote_edge_indices'], [])
        self.assertIsNotNone(r['legacy_domain'])

    def test_extreme_remote_coordinates_do_not_license_partition_overflow(self):
        r = self.check([[0., 0., 1e100]], [[1., 0., 1.]])
        self.assertEqual(r['status'], legacy.UNRESOLVED)
        self.assertIn('UNCHANGED_PARTITION_ARITHMETIC_MAGNITUDE_UNRESOLVED', r['reasons'])

    def test_empty_rows_and_invalid_slab(self):
        r = self.check(np.empty((0, 3)), np.empty((0, 3)))
        self.assertEqual(r['status'], legacy.PASS)
        self.assertEqual(r['remote_edge_indices'], [])
        with self.assertRaises(ValueError):
            self.check([[.2, 0., 10.]], [[.6, 0., 1.]], zlo=2., zhi=2.)

    def test_path_ast_only_changes_domain_call_and_partition_postcondition(self):
        audit = remote.static_path_adaptation_audit()
        self.assertTrue(audit['passed'])
        self.assertEqual(audit['domain_calls'], 1)
        self.assertEqual(audit['consistency_assertions'], 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)

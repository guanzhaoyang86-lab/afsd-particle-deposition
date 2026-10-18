"""Independent declared-path shell checks; no material/process solve."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from batched_contact import BatchedChainRuleContact
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from contact_shell_audit import audit_contact_shell, CERTIFIED, FOUND, UNRESOLVED


def pose(theta=0., origin=(0., 0., 0.)):
    c, s = np.cos(theta), np.sin(theta)
    return dict(origin=np.array(origin), rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


class ShellAuditTests(unittest.TestCase):
    def setUp(self):
        self.contact = BatchedChainRuleContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        self.p = pose()

    def audit(self, first, last=None, **kwargs):
        first = np.asarray(first, float).reshape(-1, 3)
        last = first.copy() if last is None else np.asarray(last, float).reshape(-1, 3)
        return audit_contact_shell(self.contact, first, last, self.p, self.p, **kwargs)

    def test_stationary_and_small_separated_motion_certify(self):
        for last in (None, [[.02, .01, .51]]):
            result = self.audit([[.01, 0., .5]], last)
            self.assertEqual(result['status'], CERTIFIED)
            self.assertEqual(result['geometry_queries'], 2)
            self.assertEqual(result['endpoint_inactive_pair_count'], 7)
            self.assertIsNone(result['solid_path_valid'])

    def test_rotation_hidden_wall_shell_with_inactive_endpoints(self):
        x = np.array([[.405, 0., .5]])
        result = audit_contact_shell(self.contact, x, x, pose(-np.pi/4), pose(np.pi/4))
        self.assertEqual(result['status'], FOUND)
        event = next(e for e in result['hidden_contacts'] if e['feature'] == 'BORE_X_PLUS')
        self.assertEqual(event['fraction'], .5)
        self.assertAlmostEqual(event['gap_m'], .095)
        self.assertFalse(event['within_roundoff_margin'])
        self.assertLessEqual(result['geometry_queries'], 256)

    def test_query_budget_and_depth_exhaustion_are_unresolved(self):
        x = np.array([[.405, 0., .5]])
        for kwargs, reason in ((dict(maximum_queries=2), 'MAXIMUM_GEOMETRY_QUERIES'),
                               (dict(maximum_depth=0), 'MAXIMUM_SUBDIVISION_DEPTH')):
            result = audit_contact_shell(self.contact, x, x, pose(-np.pi/4), pose(np.pi/4), **kwargs)
            self.assertEqual(result['status'], UNRESOLVED)
            self.assertEqual(result['failure_kind'], reason)
            self.assertTrue(result['unresolved_intervals'])

    def test_narrow_event_missed_by_three_fixed_samples_is_found(self):
        x = np.array([[.4001, 0., .5]])
        p0, p1 = pose(-.125), pose(.875)
        for fraction in (0., .25, .5, .75, 1.):
            geometry = self.contact.geometry(x, pose(-.125+fraction))
            self.assertFalse(geometry['active'].any())
        result = audit_contact_shell(self.contact, x, x, p0, p1)
        self.assertEqual(result['status'], FOUND)
        event = next(e for e in result['hidden_contacts'] if e['feature'] == 'BORE_X_PLUS')
        self.assertEqual(event['fraction'], .125)
        self.assertLess(event['gap_m'], self.contact.particle_radius)

    def test_base_clear_motion_and_stationary_tangency(self):
        self.contact.base_z = -1.
        result = self.audit([[0., 0., -.8]], [[.1, 0., -.7]])
        self.assertEqual(result['status'], CERTIFIED)
        # r=.1 is exactly the base distance of z=-.9 only approximately in
        # binary arithmetic, so choose exactly representable test dimensions.
        self.contact.particle_radius = .125
        result = self.audit([[0., 0., -.875]])
        self.assertEqual(result['status'], CERTIFIED)
        self.assertIn([0, 7], result['checked_pairs'])
        bad = self.audit([[0., 0., -1.01]])
        self.assertEqual(bad['status'], UNRESOLVED)
        self.assertEqual(bad['failure_kind'], 'INVALID_SOLID_PATH')
        self.assertFalse(bad['solid_path_valid'])

    def test_common_rigid_endpoints_use_actual_chord_not_body_interpolation(self):
        first = np.array([[.6, 0., -.04]])
        last_pose = pose(np.pi/2)
        last = first @ last_pose['rotation'].T
        result = audit_contact_shell(self.contact, first, last, self.p, last_pose)
        self.assertEqual(result['status'], FOUND)
        event = next(e for e in result['hidden_contacts'] if e['feature'] == 'BORE_X_PLUS')
        self.assertEqual(event['fraction'], .5)
        expected = np.hypot(.5-.6/np.sqrt(2.), .04)
        self.assertAlmostEqual(event['gap_m'], expected)

    def test_solid_intrusion_is_explicit_failure_not_shell_pass(self):
        x = np.array([[.55, 0., .5]])
        result = audit_contact_shell(self.contact, x, x, pose(-np.pi/4), pose(np.pi/4))
        self.assertEqual(result['status'], UNRESOLVED)
        self.assertEqual(result['failure_kind'], 'INVALID_SOLID_PATH')
        self.assertFalse(result['solid_path_valid'])
        self.assertEqual(result['invalid_geometry']['fraction'], .5)
        self.assertIn('inside', result['invalid_geometry']['reason'])

    def test_scalar_batch_distance_paths_match_and_inputs_remain_pure(self):
        scalar = ChainRuleDiscreteContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        x = np.array([[.405, 0., .5]])
        before = x.copy()
        p0, p1 = pose(-np.pi/4), pose(np.pi/4)
        actual = audit_contact_shell(self.contact, x, x, p0, p1)
        expected = audit_contact_shell(scalar, x, x, p0, p1)
        self.assertEqual(actual, expected)
        np.testing.assert_array_equal(x, before)

    def test_active_endpoint_pair_is_explicitly_outside_audit_scope(self):
        result = self.audit([[.45, .1, .5]])
        self.assertEqual(result['status'], CERTIFIED)
        self.assertNotIn([0, 3], result['checked_pairs'])
        self.assertIn('active-endpoint', result['limitation'])

    def test_empty_input_and_invalid_budget(self):
        result = audit_contact_shell(self.contact, np.empty((0, 3)), np.empty((0, 3)), self.p, self.p)
        self.assertEqual(result['status'], CERTIFIED)
        self.assertEqual(result['checked_pairs'], [])
        for kwargs in (dict(maximum_queries=1), dict(maximum_depth=33), dict(maximum_depth=True)):
            with self.assertRaises(ValueError):
                self.audit([[0., 0., .5]], **kwargs)


if __name__ == '__main__':
    unittest.main(verbosity=2)

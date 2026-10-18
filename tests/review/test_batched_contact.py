"""Scalar/batch implementation equivalence only; no AFSD validation claim."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from batched_contact import BatchedChainRuleContact, IMPLEMENTATION
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from afsd_backend.pysph.quasistatic_contact import ContactGeometryError
import test_discrete_contact as original_tests

pose = original_tests.pose


class BatchedInheritedChainTests(original_tests.ChainRuleDiscreteContactTests):
    def setUp(self):
        self.contact = BatchedChainRuleContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        self.p = pose()


class BatchedEquivalenceTests(unittest.TestCase):
    def setUp(self):
        args = (2., 1., 1., .1, 100., 200., .3)
        self.scalar = ChainRuleDiscreteContact(*args, base_z=None)
        self.batch = BatchedChainRuleContact(*args, base_z=None)
        self.p = pose()

    def compare(self, actual, expected):
        self.assertEqual(actual.keys(), expected.keys())
        for key, value in expected.items():
            with self.subTest(field=key):
                if isinstance(value, np.ndarray):
                    if value.dtype == bool:
                        np.testing.assert_array_equal(actual[key], value)
                    else:
                        np.testing.assert_allclose(actual[key], value, rtol=3e-13, atol=3e-13)
                else:
                    self.assertEqual(actual[key], value)

    def matching_rejection(self, points, pattern=None):
        messages = []
        for contact in (self.scalar, self.batch):
            with self.assertRaises(ContactGeometryError) as caught:
                contact.geometry(np.asarray(points, float), self.p)
            messages.append(str(caught.exception))
        self.assertEqual(messages[0], messages[1])
        if pattern:
            self.assertIn(pattern, messages[0])

    def test_random_geometry_all_branches_and_scales(self):
        rng = np.random.default_rng(614209)
        bore = rng.uniform([-.49, -.49, -.5], [.49, .49, 1.5], (250, 3))
        caps = rng.uniform([-2.5, -2.5, -.2], [2.5, 2.5, -.005], (250, 3))
        top = caps.copy(); top[:, 2] = 1.-top[:, 2]
        theta = rng.uniform(-np.pi, np.pi, 250)
        radius = rng.uniform(2.001, 2.4, 250)
        outer = np.column_stack([radius*np.cos(theta), radius*np.sin(theta), rng.uniform(-.2, 1.2, 250)])
        points = np.vstack((bore, caps, top, outer))
        for scale in (1., .001):
            args = (2.*scale, scale, scale, .1*scale, 100., 200., .3)
            old = ChainRuleDiscreteContact(*args, base_z=None)
            new = BatchedChainRuleContact(*args, base_z=None)
            for moved in (pose(), pose(np.array([.03, -.04, .02])*scale, .3)):
                world = points*scale @ moved['rotation'].T + moved['origin']
                self.compare(new.geometry(world, moved), old.geometry(world, moved))

    def test_exact_and_near_medial_axis_tie_rejections(self):
        for point in ([.47, .47, -.04], [-.47, .47, 1.04],
                      [.47, .47+np.spacing(.47), -.04]):
            self.matching_rejection([point], 'multiple closest')
        # A tie too far from both caps is inactive; both bore walls contribute.
        points = np.array([[.45, .45, .5], [0., 0., .5], [0., 0., 3.]])
        self.compare(self.batch.geometry(points, self.p), self.scalar.geometry(points, self.p))

    def test_rejection_order_and_exact_tool_boundary(self):
        for points in (([.47, .47, -.04], [.7, .1, .5]),
                       ([.7, .1, .5], [.47, .47, -.04]),
                       ([.5, .1, .5],), ([2., 0., .5],), ([.8, .1, 0.],)):
            self.matching_rejection(points)

    def test_ulp_geometry_and_active_shell_branches_match_scalar(self):
        points = []
        def around(value):
            return (np.nextafter(value, -np.inf), value, np.nextafter(value, np.inf))
        # Non-axis-aligned circle points expose scalar-dot vs batched-reduce
        # norm roundoff; invalid scalar rows must never become batch accepted.
        for theta in np.linspace(.001, 2*np.pi-.001, 31):
            for radius in around(2.):
                points.append([radius*np.cos(theta), radius*np.sin(theta), .5])
        for side in (-1., 1.):
            for wall in around(side*.5):
                points.append([wall, .1, .5])
                points.append([.1, wall, .5])
            for coordinate in around(side*.4):
                points.append([coordinate, .1, .5])
                points.append([.1, coordinate, .5])
        for z in (*around(0.), *around(1.)):
            points.extend(([.8, .1, z], [.4, .1, z]))
        for z in (*around(-.1), *around(1.1)):
            points.append([.8, .1, z])
        # Compare one row at a time so acceptance and the exact diagnostic
        # cannot be obscured by an earlier rejected row in the batch.
        for point in points:
            with self.subTest(point=point):
                x = np.asarray([point])
                try:
                    expected = self.scalar.geometry(x, self.p)
                except ContactGeometryError as error:
                    with self.assertRaises(ContactGeometryError) as caught:
                        self.batch.geometry(x, self.p)
                    self.assertEqual(str(caught.exception), str(error))
                else:
                    actual = self.batch.geometry(x, self.p)
                    for field in expected:
                        np.testing.assert_array_equal(actual[field], expected[field])

    def test_all_work_fields_multi_contact_friction_and_lip(self):
        x0 = np.array([[.44, .44, .5], [2.04, 0., .4], [.46, .1, -.03], [.36, .1, .5]])
        x1 = np.array([[.46, .46, .51], [2.06, .003, .41], [.48, .105, -.05], [.42, .16, .5]])
        history = np.zeros((4, 8, 3)); history[0, 3, 2] = .003
        before = [v.copy() for v in (x0, x1, history)]
        self.compare(self.batch.evaluate(x0, x1, self.p, self.p, .1, history),
                     self.scalar.evaluate(x0, x1, self.p, self.p, .1, history))
        for value, original in zip((x0, x1, history), before):
            np.testing.assert_array_equal(value, original)

    def test_common_rigid_motion_full_output(self):
        x0 = np.array([[.45, .1, .5], [.1, -.44, .4]])
        last_pose = pose([.03, -.04, .02], .25)
        x1 = x0 @ last_pose['rotation'].T+last_pose['origin']
        history = np.zeros((2, 8, 3)); history[0, 3, 1] = .003
        self.compare(self.batch.evaluate(x0, x1, self.p, last_pose, .1, history),
                     self.scalar.evaluate(x0, x1, self.p, last_pose, .1, history))

    def test_inherited_sweep_rejects_tunneling_and_certifies_open_bore(self):
        self.assertIs(BatchedChainRuleContact._sweep, ChainRuleDiscreteContact._sweep)
        self.assertIs(BatchedChainRuleContact._patches, ChainRuleDiscreteContact._patches)
        self.assertIs(BatchedChainRuleContact.evaluate, ChainRuleDiscreteContact.evaluate)
        for contact in (self.scalar, self.batch):
            with self.assertRaises(ContactGeometryError):
                contact.evaluate(np.array([[0., 0., .5]]), np.array([[2.2, 0., .5]]), self.p, self.p, .1)
        first = np.array([[0., 0., 1.2]]); last = np.array([[0., 0., -.2]])
        self.compare(self.batch.evaluate(first, last, self.p, self.p, .1),
                     self.scalar.evaluate(first, last, self.p, self.p, .1))

    def test_base_is_laboratory_geometry_and_below_base_rejects(self):
        args = (2., 1., 1., .1, 100., 200., .3)
        old = ChainRuleDiscreteContact(*args, base_z=-1.)
        new = BatchedChainRuleContact(*args, base_z=-1.)
        x0 = np.array([[0., 0., -.92]]); x1 = np.array([[0., 0., -.96]])
        moved = pose([.05, -.02, 0.], .2)
        self.compare(new.evaluate(x0, x1, self.p, moved, .1), old.evaluate(x0, x1, self.p, moved, .1))
        for contact in (old, new):
            with self.assertRaisesRegex(ContactGeometryError, 'below rigid base'):
                contact.geometry(np.array([[0., 0., -1.01]]), self.p)

    def test_empty_batch_and_scalar_rejects_ties_flag(self):
        points = np.empty((0, 3))
        self.compare(self.batch.geometry(points, self.p), self.scalar.geometry(points, self.p))
        p = np.array([[.47, .47, -.04]])
        actual = self.batch._patches_batch(p, reject_active_ties=False)
        expected = self.scalar._patches(p[0], reject_active_ties=False)
        for value, oracle in zip(actual, expected):
            np.testing.assert_array_equal(value[0], oracle)
        self.assertEqual(self.batch.implementation, IMPLEMENTATION)


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Analytic geometry/spring references; component tests, never AFSD validation."""
from pathlib import Path
import inspect
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.quasistatic_contact import (
    QuasistaticContact, ContactGeometryError, FEATURES)


def pose(origin=(0., 0., 0.), theta=0.):
    c, s = np.cos(theta), np.sin(theta)
    return dict(origin=np.asarray(origin, float),
                rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


class QuasistaticContactTests(unittest.TestCase):
    def setUp(self):
        self.contact = QuasistaticContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        self.p = pose()

    def close(self, actual, expected, tol=1e-12):
        np.testing.assert_allclose(actual, expected, atol=tol, rtol=0.)

    def evaluate(self, first, last=None, history=None):
        last = first if last is None else last
        return self.contact.evaluate(np.asarray(first, float).reshape(-1, 3),
            np.asarray(last, float).reshape(-1, 3), self.p, self.p, .1, history)

    def test_open_bore_passes_both_planes_and_top_remains_open(self):
        r = self.evaluate([0., 0., 1.2], [0., 0., -.2])
        self.close(r['force_N'], [[0., 0., 0.]])
        self.assertTrue(r['sweep']['certified'])
        self.assertGreater(r['sweep']['intermediate_geometry_queries'], 0)
        self.assertFalse(np.any(r['active']))

    def test_finite_lip_is_three_four_five_and_duplicate_patches_explicit(self):
        r = self.evaluate([.47, .1, -.04])
        self.close(r['gaps_m'][0, [0, 3]], [.05, .05])
        self.close(r['normals'][0, [0, 3]], [[-.6, 0., -.8], [-.6, 0., -.8]])
        self.close(r['normal_force_by_feature_N'][0, [0, 3]], [[-3., 0., -4.], [-3., 0., -4.]])
        self.close(r['normal_potential1_J'].sum(), .25)
        self.assertIn('DUPLICATE_LIP_CONTRIBUTIONS_RETAINED', r['normal_model'])

    def test_two_bore_faces_admit_actual_multi_contact(self):
        r = self.evaluate([.45, .45, .5], [.45, .45, .52])
        self.close(r['normal_force_by_feature_N'].sum(axis=1), [[-5., -5., 0.]])
        self.close(r['tangential_force_by_feature_N'].sum(axis=1), [[0., 0., -3.]])
        self.close(r['friction_dissipation_J'].sum(), 2*.01875)
        self.assertEqual(set(np.flatnonzero(r['active'][0])), {3, 5})

    def test_top_and_outer_cylinder_are_finite(self):
        for point, active in (([.47, .1, 1.04], (1, 3)), ([2.03, 0., -.04], (0, 6)),
                              ([2.04, 0., .5], (6,))):
            r = self.evaluate(point)
            self.assertEqual(set(np.flatnonzero(r['active'][0])), set(active))
        r = self.evaluate([.8, .2, -1.])
        self.assertFalse(np.any(r['active']))

    def test_normal_force_is_negative_literal_lip_potential_gradient(self):
        point = np.array([.47, .1, -.04])
        h = 1e-7
        # Two finite patches share this lip; potential = 2*(kn/2)*(r-distance)^2.
        def reference_U(x):
            return 100*(.1-np.sqrt((x[0]-.5)**2+x[2]**2))**2
        gradient = np.array([(reference_U(point+h*e)-reference_U(point-h*e))/(2*h)
                             for e in np.eye(3)])
        r = self.evaluate(point)
        self.close(r['force_N'][0], -gradient, 3e-7)
        self.close(r['normal_potential1_J'].sum(), reference_U(point))

    def test_conservative_gradient_on_all_smooth_finite_patch_branches(self):
        h = 1e-7
        for point in ([.45, .1, .5], [.45, .45, .5], [.8, .2, -.04],
                      [.47, .1, -.04], [.47, .1, 1.04], [2.03, 0., -.04], [2.04, 0., .4]):
            point = np.asarray(point)
            gradient = []
            for e in np.eye(3):
                up = self.contact.geometry((point+h*e)[None, :], self.p)['potential_J'].sum()
                down = self.contact.geometry((point-h*e)[None, :], self.p)['potential_J'].sum()
                gradient.append((up-down)/(2*h))
            force = self.contact.geometry(point[None, :], self.p)['force_by_feature_N'].sum(axis=(0, 1))
            self.close(force, -np.asarray(gradient), 3e-7)

    def test_stick_spring_force_storage_and_BE_dissipation_are_distinct(self):
        r = self.evaluate([.45, .1, .5], [.45, .101, .5])
        self.close(r['force_N'], [[-5., -.2, 0.]])
        self.close(r['tangential_potential1_J'].sum(), .0001)
        self.close(r['algorithmic_dissipation_J'].sum(), .0001)
        self.close(r['friction_dissipation_J'].sum(), 0.)
        self.assertFalse(np.any(r['sliding']))

    def test_radial_return_matches_literal_Coulomb_slip_and_energy(self):
        r = self.evaluate([.45, .1, .5], [.45, .12, .5])
        self.close(r['force_N'], [[-5., -1.5, 0.]])
        self.close(r['history'][0, 3], [0., .0075, 0.])
        self.close(r['friction_dissipation_J'].sum(), .01875)
        self.close(r['tangential_potential1_J'].sum(), .005625)
        self.close(r['algorithmic_dissipation_J'].sum(), .005625)
        self.close(r['particle_work_J'].sum(), -.03)
        self.close(r['tangential_work_identity_error_J'], np.zeros((1, 8)))

    def test_reverse_motion_unloads_then_reverses_without_velocity_overwrite(self):
        first = self.evaluate([.45, .1, .5], [.45, .12, .5])
        unloading = self.evaluate([.45, .12, .5], [.45, .118, .5], first['history'])
        self.close(unloading['tangential_force_by_feature_N'][0, 3], [0., -1.1, 0.])
        self.close(unloading['friction_dissipation_J'].sum(), 0.)
        reverse = self.evaluate([.45, .12, .5], [.45, .08, .5], first['history'])
        self.close(reverse['tangential_force_by_feature_N'][0, 3], [0., 1.5, 0.])
        self.assertGreater(reverse['friction_dissipation_J'].sum(), 0.)

    def test_lost_contact_releases_old_spring_energy_in_separate_account(self):
        history = np.zeros((1, 8, 3)); history[0, 3, 1] = .005
        r = self.evaluate([.45, .1, .5], [.2, .1, .5], history)
        self.close(r['history'], np.zeros((1, 8, 3)))
        self.close(r['release_energy_J'].sum(), .0025)
        self.close(r['friction_dissipation_J'].sum(), 0.)
        self.close(r['algorithmic_dissipation_J'].sum(), 0.)
        self.assertEqual(r['release_energy_policy'], 'UNASSIGNED_ACCOUNT_NOT_AUTOMATICALLY_HEAT')

    def test_changing_normal_reports_projection_energy_without_clipping(self):
        history = np.zeros((1, 8, 3)); history[0, 3, 2] = .005
        r = self.evaluate([.45, .1, .01], [.47, .1, -.04], history)
        self.close(r['release_energy_J'][0, 3], .0016)
        self.close(r['tangential_work_identity_error_J'], np.zeros((1, 8)))

    def test_common_finite_tool_rotation_is_objective_and_has_zero_relative_slip(self):
        first = np.array([[.45, .1, .5]])
        final_pose = pose([.03, -.04, .02], .25)
        last = first@final_pose['rotation'].T+final_pose['origin']
        history = np.zeros((1, 8, 3)); history[0, 3, 1] = .003
        # Linear particle path and rotating tool may approach each other between
        # endpoints; this small-angle path remains outside with positive margin.
        r = self.contact.evaluate(first, last, self.p, final_pose, .1, history)
        self.close(r['history'], history)
        self.close(r['friction_dissipation_J'].sum(), 0.)
        self.close(r['normal_endpoint_work_error_J'].sum(), 0.)
        self.close(r['particle_work_J'], r['wall_work_J'])

    def test_rigid_pose_derivative_equals_normal_wall_power(self):
        point = np.array([[.47, .1, -.04]])
        translation = np.array([.3, -.2, .1]); omega = np.array([0., 0., 2.])
        h = 1e-7
        def potential(t):
            return self.contact.geometry(point, pose(t*translation, 2*t))['potential_J'].sum()
        derivative = (potential(h)-potential(-h))/(2*h)
        force = self.contact.geometry(point, self.p)['force_by_feature_N'].sum(axis=(0, 1))
        velocity = translation+np.cross(omega, point[0])
        self.close(derivative, force@velocity, 5e-7)

    def test_endpoint_normal_work_error_is_raw_and_converges_under_refinement(self):
        def accumulate(count):
            error = 0.
            for k in range(count):
                start = .44+.02*k/count; end = .44+.02*(k+1)/count
                r = self.evaluate([start, .1, .5], [end, .1, .5])
                error += r['normal_endpoint_work_error_J'].sum()
            return error
        self.close(accumulate(1), -.02)
        self.close(accumulate(2), -.01)
        self.close(accumulate(4), -.005)

    def test_lip_transition_potential_and_force_remain_continuous(self):
        eps = 1e-8
        for pair in (([.45, .1, eps], [.45, .1, -eps]),
                     ([.5-eps, .1, -.05], [.5+eps, .1, -.05])):
            first, last = [self.contact.geometry(np.array([p]), self.p) for p in pair]
            self.close(first['potential_J'], last['potential_J'], 1e-12)
            self.close(first['force_by_feature_N'], last['force_by_feature_N'], 3e-6)

    def test_inside_endpoint_and_swept_tunneling_rejected(self):
        with self.assertRaisesRegex(ContactGeometryError, 'inside'):
            self.evaluate([.7, .1, .5])
        # Both endpoints outside; the straight path crosses the annular solid.
        with self.assertRaises(ContactGeometryError):
            self.evaluate([0., 0., .5], [2.2, 0., .5])
        r = self.contact.evaluate(np.array([[0., 0., .5]]), np.array([[2.2, 0., .5]]),
                                 self.p, self.p, .1, check_path=False)
        self.assertFalse(r['sweep']['certified'])
        with self.assertRaises(ContactGeometryError):
            self.contact.evaluate(np.array([[0., 0., .5]]), np.array([[.7, 0., .5]]),
                                  self.p, self.p, .1, check_path=False)

    def test_active_cap_medial_axis_rejected_but_multi_bore_contact_allowed(self):
        with self.assertRaisesRegex(ContactGeometryError, 'multiple closest'):
            self.evaluate([.47, .47, -.04])
        self.evaluate([.45, .45, .5])

    def test_base_sphere_surface_tangency_and_compression_are_valid(self):
        c = QuasistaticContact(2., 1., 1., .1, 100., 200., .3, base_z=-1.)
        for z, force in ((-.9, 0.), (-.95, 5.), (-1., 10.)):
            x = np.array([[0., 0., z]])
            r = c.evaluate(x, x, self.p, self.p, .1)
            self.close(r['force_N'], [[0., 0., force]])
        with self.assertRaisesRegex(ContactGeometryError, 'below rigid base'):
            c.geometry(np.array([[0., 0., -1.01]]), self.p)

    def test_reactions_and_torque_balance_all_features(self):
        x0 = np.array([[.45, .45, .5], [2.04, 0., .4]])
        x1 = x0+np.array([0., 0., .01])
        r = self.contact.evaluate(x0, x1, self.p, self.p, .1)
        self.close(r['particle_reaction_on_tool_N']+r['force_N'].sum(axis=0), np.zeros(3))
        self.close(r['particle_reaction_torque_on_tool_Nm']+np.cross(x1, r['force_N']).sum(axis=0), np.zeros(3))

    def test_inputs_are_not_mutated_and_labels_absent(self):
        first = np.array([[.45, .1, .5]]); last = np.array([[.45, .12, .5]])
        history = np.zeros((1, 8, 3))
        copies = [a.copy() for a in (first, last, history)]
        self.contact.evaluate(first, last, self.p, self.p, .1, history)
        for actual, expected in zip((first, last, history), copies):
            np.testing.assert_array_equal(actual, expected)
        parameters = inspect.signature(QuasistaticContact.evaluate).parameters
        self.assertFalse(any(word in name.lower() for name in parameters
                             for word in ('label', 'tracer', 'as_rs', 'mass')))
        self.assertEqual(len(FEATURES), 8)

    def test_invalid_inputs_reject_without_hidden_repair(self):
        with self.assertRaises(ValueError):
            QuasistaticContact(2., 1., 1., .1, 100., 200., -.1)
        with self.assertRaises(ValueError):
            self.contact.evaluate(np.zeros((1, 3)), np.zeros((1, 3)), self.p, self.p, 0.)
        with self.assertRaises(ValueError):
            self.contact.geometry(np.array([[np.nan, 0., 0.]]), self.p)
        with self.assertRaises(ValueError):
            self.contact.geometry(np.zeros((1, 3)), dict(origin=np.zeros(3), rotation=-np.eye(3)))


if __name__ == '__main__':
    unittest.main(verbosity=2)

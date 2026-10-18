"""Analytic finite-contact component checks; no AFSD process certification."""
from pathlib import Path
import inspect
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.discrete_contact import (DiscreteContact, CLASSIFICATION,
    ChainRuleDiscreteContact, CHAIN_RULE_CLASSIFICATION)
from afsd_backend.pysph.quasistatic_contact import ContactGeometryError


def pose(origin=(0., 0., 0.), theta=0.):
    c, s = np.cos(theta), np.sin(theta)
    return dict(origin=np.asarray(origin, float),
                rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


class DiscreteContactTests(unittest.TestCase):
    def setUp(self):
        self.contact = DiscreteContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        self.p = pose()

    def close(self, actual, expected, tol=1e-12):
        np.testing.assert_allclose(actual, expected, atol=tol, rtol=0.)

    def evaluate(self, first, last=None, history=None):
        last = first if last is None else last
        return self.contact.evaluate(np.asarray(first, float).reshape(-1, 3),
            np.asarray(last, float).reshape(-1, 3), self.p, self.p, .1, history)

    def test_planar_compression_uses_average_force_in_actual_equilibrium(self):
        r = self.evaluate([.44, .1, .5], [.46, .1, .5])
        self.close(r['normal_force_by_feature_N'][0, 3], [-5., 0., 0.])
        self.close(r['force_N'], [[-5., 0., 0.]])
        self.close(r['normal_particle_work_J'].sum(), -.1)
        self.close((r['normal_potential1_J']-r['normal_potential0_J']).sum(), .1)
        self.close(r['dg_identity_error_J'], np.zeros((1, 8)))
        self.close(r['normal_endpoint_quadrature_error_J'].sum(), -.02)
        self.assertEqual(r['classification'], CLASSIFICATION)

    def test_activation_and_release_are_secant_forces_not_endpoint_forces(self):
        activation = self.evaluate([.38, .1, .5], [.44, .1, .5])
        # U0=0, U1=.08 J, dx=.06 m; exact secant normal is -4/3 N.
        self.close(activation['normal_force_by_feature_N'][0, 3], [-4/3, 0., 0.])
        self.close(activation['dg_identity_error_J'], np.zeros((1, 8)))
        release = self.evaluate([.44, .1, .5], [.38, .1, .5])
        self.close(release['normal_force_by_feature_N'], activation['normal_force_by_feature_N'])
        self.close(release['normal_particle_work_J'].sum(), .08)
        self.assertFalse(np.any(release['active']))

    def test_nearly_normal_activation_has_input_scaled_projection_certificate(self):
        # Midpoint is inactive; DG force is exactly parallel to displacement
        # in real arithmetic, so true tangential motion is zero.  Normalizing
        # that force and applying I-n*n^T leaves a floating residual spring.
        # Its work is O(eps*kt*|dy|^2), not O(eps) times its tiny spring energy.
        r = self.evaluate([.36, .1, .5], [.42, .11, .5])
        feature = (0, 3)
        actual = abs(r['tangential_work_identity_error_J'][feature])
        tolerance = r['tangential_work_identity_tolerance_J'][feature]
        self.assertGreater(actual, 0.)  # Raw error was not zeroed or hidden.
        output_scale = (abs(np.sum(r['tangential_force_by_feature_N'][feature]*
            r['relative_step_body_m'][feature]))+
            r['tangential_potential0_J'][feature]+r['tangential_potential1_J'][feature]+
            abs(r['friction_dissipation_J'][feature])+r['algorithmic_dissipation_J'][feature]+
            abs(r['release_energy_J'][feature]))
        self.assertGreater(actual, 256*np.finfo(float).eps*output_scale)
        self.assertLessEqual(actual, tolerance)
        # Machine certificate remains far below normal contact work.
        self.assertLess(tolerance, 1e-10*abs(r['normal_particle_work_J'][feature]))
        self.assertGreater(r['tangential_projection_roundoff_bound_J'][feature], 0.)
        self.close(r['dg_identity_error_J'], np.zeros((1, 8)))

    def test_projection_certificate_does_not_remove_actual_small_tangential_work(self):
        r = self.evaluate([.45, .1, .5], [.45001, .1000001, .5])
        # Constant flat normal direction: exact elastic tangential spring is
        # 1e-7 m and its BE storage/loss are each 1e-12 J.
        self.close(r['history'][0, 3], [0., 1e-7, 0.], 1e-16)
        self.close(r['tangential_force_by_feature_N'][0, 3], [0., -2e-5, 0.], 2e-14)
        self.close(r['tangential_potential1_J'][0, 3], 1e-12, 2e-21)
        self.close(r['algorithmic_dissipation_J'][0, 3], 1e-12, 2e-21)
        self.assertLess(r['tangential_work_identity_tolerance_J'][0, 3], 1e-8*1e-12)
        self.close(r['friction_dissipation_J'].sum(), 0.)

    def test_lip_force_has_literal_potential_work_and_endpoint_difference(self):
        first = np.array([.46, .1, -.03]); last = np.array([.48, .1, -.05])
        r = self.evaluate(first, last)
        # Bore X+ and annular shoulder share this nearest lip.
        def literal_U(x):
            return 100*(.1-np.hypot(x[0]-.5, x[2]))**2
        expected_change = literal_U(last)-literal_U(first)
        self.close((r['normal_potential1_J']-r['normal_potential0_J']).sum(), expected_change)
        force_work = np.sum(r['normal_force_by_feature_N'], axis=(0, 1))@(last-first)
        self.close(force_work, -expected_change)
        self.close(r['dg_identity_error_J'], np.zeros((1, 8)))
        self.assertGreater(abs(r['normal_endpoint_quadrature_error_J'].sum()), 1e-5)

    def test_force_symmetric_in_endpoints_for_frictionless_geometry(self):
        first, last = [.46, .1, -.03], [.48, .1, -.05]
        forward = self.evaluate(first, last)
        reverse = self.evaluate(last, first)
        self.close(forward['normal_force_by_feature_N'], reverse['normal_force_by_feature_N'])
        self.close(forward['normal_particle_work_J'], -reverse['normal_particle_work_J'])

    def test_zero_step_and_tiny_step_have_bounded_derivative_limit(self):
        first = np.array([.47, .1, -.04])
        same = self.evaluate(first)
        self.close(same['normal_force_by_feature_N'][0, [0, 3]], [[-3., 0., -4.], [-3., 0., -4.]])
        self.assertTrue(np.all(same['normal_derivative_limit_branch']))
        tiny = self.evaluate(first, first+np.array([1e-15, 0., 0.]))
        self.close(tiny['normal_force_by_feature_N'], same['normal_force_by_feature_N'], 2e-12)
        self.assertTrue(np.all(np.abs(tiny['dg_identity_error_J']) <= tiny['dg_identity_tolerance_J']))

    def test_BE_friction_uses_discrete_load_and_keeps_algorithmic_loss(self):
        # Flat force changes 4 -> 6 N, while slip is orthogonal to the normal.
        r = self.evaluate([.44, .1, .5], [.46, .12, .5])
        self.close(r['normal_load_N'][0, 3], 5.)
        self.close(r['endpoint_normal_load_N'][0, 3], 6.)
        self.close(r['history'][0, 3], [0., .0075, 0.])
        self.close(r['tangential_force_by_feature_N'][0, 3], [0., -1.5, 0.])
        self.close(r['friction_dissipation_J'].sum(), .01875)
        self.close(r['algorithmic_dissipation_J'].sum(), .005625)
        self.close(r['tangential_work_identity_error_J'], np.zeros((1, 8)))
        self.assertEqual(r['algorithmic_dissipation_policy'], 'REPORTED_NOT_AUTOMATICALLY_HEAT')

    def test_contact_loss_releases_spring_even_with_finite_interval_load(self):
        history = np.zeros((1, 8, 3)); history[0, 3, 1] = .005
        r = self.evaluate([.44, .1, .5], [.38, .1, .5], history)
        self.assertGreater(r['normal_load_N'][0, 3], 0.)
        self.assertFalse(np.any(r['friction_active']))
        self.close(r['history'], np.zeros_like(history))
        self.close(r['release_energy_J'].sum(), .0025)
        self.close(r['friction_dissipation_J'].sum(), 0.)
        self.close(r['tangential_work_identity_error_J'], np.zeros((1, 8)))

    def test_common_finite_rigid_motion_has_no_relative_work(self):
        first = np.array([[.45, .1, .5]])
        last_pose = pose([.03, -.04, .02], .25)
        last = first@last_pose['rotation'].T+last_pose['origin']
        history = np.zeros((1, 8, 3)); history[0, 3, 1] = .003
        r = self.contact.evaluate(first, last, self.p, last_pose, .1, history)
        self.close(r['history'], history)
        self.close(r['normal_force_by_feature_N'][0, 3], np.array([-5., 0., 0.])@last_pose['rotation'].T)
        self.close(r['friction_dissipation_J'].sum(), 0.)
        self.close(r['normal_particle_work_J'], r['normal_wall_work_J'])
        self.close(r['particle_work_J'], r['wall_work_J'])
        self.assertTrue(r['sweep']['certified'])

    def test_superposed_coordinate_rotation_preserves_general_finite_step(self):
        x0 = np.array([[.46, .1, -.03]]); x1 = np.array([[.48, .105, -.05]])
        plain = self.contact.evaluate(x0, x1, self.p, self.p, .1)
        other = pose([.1, -.2, .3], .6); Q = other['rotation']; c = other['origin']
        moved = self.contact.evaluate(x0@Q.T+c, x1@Q.T+c, other, other, .1)
        self.close(moved['force_N'], plain['force_N']@Q.T, 2e-12)
        self.close(moved['history'], plain['history'])
        for key in ('particle_work_J', 'normal_potential1_J', 'friction_dissipation_J',
                    'algorithmic_dissipation_J', 'release_energy_J'):
            self.close(moved[key], plain[key])

    def test_cumulative_normal_energy_has_no_per_step_defect_budget(self):
        for count in (1, 2, 4, 17):
            total_work = 0.
            for k in range(count):
                r = self.evaluate([.42+.06*k/count, .1, .5], [.42+.06*(k+1)/count, .1, .5])
                total_work += r['normal_particle_work_J'].sum()
            # Initial U=.02 J, final U=.32 J.
            self.close(total_work, -.3)

    def test_base_uses_laboratory_increment_and_does_not_follow_tool(self):
        contact = DiscreteContact(2., 1., 1., .1, 100., 200., .3, base_z=-1.)
        x0 = np.array([[0., 0., -.92]]); x1 = np.array([[0., 0., -.96]])
        r = contact.evaluate(x0, x1, self.p, pose([.05, -.02, 0.], .2), .1)
        self.close(r['normal_force_by_feature_N'][0, 7], [0., 0., 4.])
        self.close(r['wall_step_m'][0, 7], np.zeros(3))
        self.close(r['normal_particle_work_J'].sum(), -.16)
        self.close(r['dg_identity_error_J'], np.zeros((1, 8)))

    def test_swept_solid_rejection_and_open_bore_are_retained(self):
        with self.assertRaises(ContactGeometryError):
            self.evaluate([0., 0., .5], [2.2, 0., .5])
        r = self.evaluate([0., 0., 1.2], [0., 0., -.2])
        self.assertTrue(r['sweep']['certified'])
        self.close(r['force_N'], [[0., 0., 0.]])
        with self.assertRaises(ContactGeometryError):
            self.evaluate([.7, .1, .5])

    def test_swept_admissible_large_lip_step_with_attractive_DG_is_rejected(self):
        # This path stays outside the solid (left of the bore face until it
        # drops below the shoulder), but its DG force on one feature points
        # against that feature's endpoint outward normal.  Energy consistency
        # alone is insufficient to accept such a large contact increment.
        first = [0.4131859470338192, .1, 0.08424301051916916]
        last = [0.5086697285535392, .1, -0.02752729085247463]
        with self.assertRaisesRegex(ContactGeometryError, 'attracts'):
            self.evaluate(first, last)

    def test_multiple_bore_patches_and_reactions_remain_separate(self):
        first = np.array([[.44, .44, .5], [2.04, 0., .4]])
        last = np.array([[.46, .46, .51], [2.06, .003, .41]])
        r = self.contact.evaluate(first, last, self.p, self.p, .1)
        self.assertEqual(set(np.flatnonzero(r['active'][0])), {3, 5})
        self.close(r['particle_reaction_on_tool_N']+r['force_N'].sum(axis=0), np.zeros(3))
        self.close(r['particle_reaction_torque_on_tool_Nm']+
                   np.cross(last, r['force_N']).sum(axis=0), np.zeros(3))
        self.close(r['dg_identity_error_J'], np.zeros((2, 8)))
        self.close(r['tangential_work_identity_error_J'], np.zeros((2, 8)))

    def test_inputs_pure_and_no_labels_or_particle_mass_in_contact_API(self):
        x0 = np.array([[.44, .1, .5]]); x1 = np.array([[.46, .12, .5]])
        history = np.zeros((1, 8, 3))
        copies = [a.copy() for a in (x0, x1, history)]
        r = self.contact.evaluate(x0, x1, self.p, self.p, .1, history)
        for a, before in zip((x0, x1, history), copies):
            np.testing.assert_array_equal(a, before)
        self.assertFalse(np.shares_memory(r['history'], history))
        parameters = inspect.signature(DiscreteContact.evaluate).parameters
        self.assertFalse(any(word in name.lower() for name in parameters
                             for word in ('label', 'tracer', 'as_rs', 'mass')))


class ChainRuleDiscreteContactTests(unittest.TestCase):
    close = DiscreteContactTests.close
    evaluate = DiscreteContactTests.evaluate

    def setUp(self):
        self.contact = ChainRuleDiscreteContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        self.p = pose()

    def test_onset_large_tangential_motion_keeps_physical_plane_normal(self):
        old = DiscreteContact(2., 1., 1., .1, 100., 200., .3, base_z=None)
        for tangent in (.01, .06):
            first = np.array([[.36, .1, .5]]); last = np.array([[.42, .1+tangent, .5]])
            before = old.evaluate(first, last, self.p, self.p, .1)
            after = self.contact.evaluate(first, last, self.p, self.p, .1)
            # phi change=.02 J, normal displacement=.06 m.  The exact scalar
            # secant is 1/3 N, regardless of the tangential displacement.
            self.close(after['normal_force_by_feature_N'][0, 3], [-1/3, 0., 0.])
            self.assertGreater(abs(before['normal_force_by_feature_N'][0, 3, 1]), .05)
            self.close(after['normal_particle_work_J'].sum(), -.02)
            self.close(after['dg_identity_error_J'], np.zeros((1, 8)))
            self.assertGreater(after['friction_dissipation_J'].sum(), 0.)
            self.assertEqual(after['classification'], CHAIN_RULE_CLASSIFICATION)

    def test_equal_distance_uses_finite_derivative_and_tangential_spring(self):
        r = self.evaluate([.45, .1, .5], [.45, .12, .5])
        self.close(r['normal_force_by_feature_N'][0, 3], [-5., 0., 0.])
        self.close(r['tangential_force_by_feature_N'][0, 3], [0., -1.5, 0.])
        self.close(r['normal_particle_work_J'], np.zeros((1, 8)))
        self.close(r['friction_dissipation_J'].sum(), .01875)
        self.close(r['algorithmic_dissipation_J'].sum(), .005625)

    def test_activation_reversal_and_release_preserve_conservative_force(self):
        forward = self.evaluate([.36, .1, .5], [.42, .16, .5])
        reverse = self.evaluate([.42, .16, .5], [.36, .1, .5], forward['history'])
        self.close(forward['normal_force_by_feature_N'], reverse['normal_force_by_feature_N'])
        self.close(forward['normal_particle_work_J'], -reverse['normal_particle_work_J'])
        self.close(reverse['history'], np.zeros((1, 8, 3)))
        self.assertGreater(reverse['release_energy_J'].sum(), 0.)
        self.close(reverse['tangential_work_identity_error_J'], np.zeros((1, 8)))

    def test_curved_lip_chain_rule_has_literal_potential_work(self):
        first = np.array([.46, .1, -.03]); last = np.array([.48, .1, -.05])
        r = self.evaluate(first, last)
        def literal_U(x):
            return 100*(.1-np.hypot(x[0]-.5, x[2]))**2
        self.close(r['normal_particle_work_J'].sum(), -(literal_U(last)-literal_U(first)))
        self.close(r['dg_identity_error_J'], np.zeros((1, 8)))
        self.close(r['tangential_work_identity_error_J'], np.zeros((1, 8)))

    def test_frozen_v2_particle250_plane_artifact_is_removed_without_resolve(self):
        # Recorded body coordinates from v2/macro independent_initial.npz and
        # step_001/accepted.npz, row 250.  This re-evaluates a fixed contact
        # state, not a new coupled AFSD solution or a temporal convergence test.
        first = np.array([[-.0005, .0010000000000000002, .0007500000000000002]])
        last = np.array([[-.0004925028894589164, .0010040777741756582, .000749135733735804]])
        args = (.0027, .002504, .001, .00025, 657500., 657500., .3)
        old = DiscreteContact(*args, base_z=None).evaluate(first, last, self.p, self.p, .0005)
        new = ChainRuleDiscreteContact(*args, base_z=None).evaluate(first, last, self.p, self.p, .0005)
        self.close(np.linalg.norm(old['normal_force_by_feature_N'][0, 5, [0, 2]]),
                   .13487008744158116, 1e-12)
        self.close(new['normal_force_by_feature_N'][0, 5, [0, 2]], np.zeros(2), 0.)
        self.assertLess(new['normal_force_by_feature_N'][0, 5, 1], 0.)
        self.close(new['normal_particle_work_J'], old['normal_particle_work_J'], 1e-18)
        self.close(new['dg_identity_error_J'], np.zeros((1, 8)), 1e-18)

    # These general geometry/rigid-motion invariants must hold for both forms.
    test_common_finite_rigid_motion_has_no_relative_work = DiscreteContactTests.test_common_finite_rigid_motion_has_no_relative_work
    test_superposed_coordinate_rotation_preserves_general_finite_step = DiscreteContactTests.test_superposed_coordinate_rotation_preserves_general_finite_step
    test_swept_solid_rejection_and_open_bore_are_retained = DiscreteContactTests.test_swept_solid_rejection_and_open_bore_are_retained
    test_cumulative_normal_energy_has_no_per_step_defect_budget = DiscreteContactTests.test_cumulative_normal_energy_has_no_per_step_defect_budget


if __name__ == '__main__':
    unittest.main(verbosity=2)

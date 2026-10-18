"""Independent literal/Euclidean references; never an AFSD integration test.

Preregistered tolerances and limits: PROCESS_FINITE_TOOL_COMPONENT.md.
"""
from pathlib import Path
import inspect
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.finite_tool_geometry import (
    FiniteOpenSquareTool, NonsmoothContactError, normal_penalty)


class FiniteToolGeometryTests(unittest.TestCase):
    def setUp(self):
        self.tool = FiniteOpenSquareTool(2., 1., 1.)

    def close(self, actual, expected, tol=1e-12):
        np.testing.assert_allclose(actual, expected, atol=tol, rtol=0.)

    def test_01_open_bore_not_capped_top_or_bottom(self):
        for z in (-.2, 0., .5, 1., 1.2):
            q = self.tool.query([0., 0., z])
            distance = .5 if 0. <= z <= 1. else np.sqrt(.25+.04)
            self.close(q['signed_distance_m'], distance)
            self.assertFalse(q['inside_solid'])
            p = normal_penalty(q, .1, 100.)
            self.assertFalse(p['active']); self.close(p['normal_force_N'], [0., 0., 0.])

    def test_02_smooth_bore_face_literal(self):
        q = self.tool.query([.45, .1, .4])
        self.close(q['signed_distance_m'], .05)
        self.close(q['closest_points_m'], [[.5, .1, .4]])
        self.close(q['normal'], [-1., 0., 0.])
        self.assertEqual(q['boundary_features'], ('BORE_X_PLUS',))
        p = normal_penalty(q, .1, 100.)
        self.close(p['normal_force_N'], [-5., 0., 0.]); self.close(p['potential_J'], .125)

    def test_03_flat_shoulder_remains_plane_outside_square(self):
        q = self.tool.query([.8, .2, -.04])
        self.close(q['signed_distance_m'], .04)
        self.close(q['closest_points_m'], [[.8, .2, 0.]])
        self.close(q['normal'], [0., 0., -1.])

    def test_04_inner_mouth_exact_diagonal_not_boolean_sdf(self):
        q = self.tool.query([.47, .1, -.04])
        self.close(q['signed_distance_m'], .05)  # 3-4-5 triangle, not max(.03,.04)
        self.close(q['closest_points_m'], [[.5, .1, 0.]])
        self.close(q['normal'], [-.6, 0., -.8])
        self.assertEqual(q['boundary_features'], ('BORE_X_PLUS', 'SHOULDER'))
        self.assertTrue(q['boundary_is_sharp'])
        p = normal_penalty(q, .1, 100.)
        self.close(p['normal_force_N'], [-3., 0., -4.]); self.close(p['potential_J'], .125)

    def test_05_top_open_mouth_exact_distance(self):
        q = self.tool.query([.47, .1, 1.04])
        self.close(q['signed_distance_m'], .05)
        self.close(q['closest_points_m'], [[.5, .1, 1.]])
        self.close(q['normal'], [-.6, 0., .8])

    def test_06_outer_shoulder_rim_exact_diagonal(self):
        q = self.tool.query([2.03, 0., -.04])
        self.close(q['signed_distance_m'], .05)
        self.close(q['closest_points_m'], [[2., 0., 0.]])
        self.close(q['normal'], [.6, 0., -.8])
        self.assertEqual(q['boundary_features'], ('OUTER_CYLINDER', 'SHOULDER'))

    def test_07_inside_solid_signed_distance_and_outward_normal(self):
        q = self.tool.query([.55, .1, .4])
        self.assertTrue(q['inside_solid'])
        self.close(q['signed_distance_m'], -.05); self.close(q['normal'], [-1., 0., 0.])
        p = normal_penalty(q, .1, 100.)
        self.close(p['normal_force_N'], [-15., 0., 0.]); self.close(p['potential_J'], 1.125)

    def test_08_concave_bore_corner_inside_tool_has_unique_diagonal(self):
        q = self.tool.query([.53, .54, .5])
        self.close(q['signed_distance_m'], -.05)
        self.close(q['closest_points_m'], [[.5, .5, .5]])
        self.close(q['normal'], [-.6, -.8, 0.])

    def test_09_multiple_closest_points_rejected_when_active(self):
        q = self.tool.query([.45, .45, .5])
        self.assertEqual(len(q['closest_points_m']), 2)
        self.assertIsNone(q['normal'])
        with self.assertRaisesRegex(NonsmoothContactError, 'MULTIPLE_CLOSEST_POINTS'):
            normal_penalty(q, .1, 100.)
        self.assertFalse(normal_penalty(q, .01, 100.)['active'])

    def test_10_exact_sharp_edge_rejected_not_rounded(self):
        for point in ([.5, .1, 0.], [2., 0., 0.], [.5, .5, .5]):
            q = self.tool.query(point)
            self.close(q['signed_distance_m'], 0.)
            with self.assertRaisesRegex(NonsmoothContactError, 'BOUNDARY_EDGE_OR_CORNER'):
                normal_penalty(q, .1, 100.)

    def test_11_exact_smooth_boundary_normal(self):
        q = self.tool.query([.5, .1, .5])
        self.close(q['normal'], [-1., 0., 0.]); self.close(q['signed_distance_m'], 0.)
        self.assertEqual(q['normal_status'], 'UNIQUE_SMOOTH_BOUNDARY_NORMAL')

    def test_12_general_3d_rotation_covariance(self):
        # Cyclic permutation is proper rotation: body (x,y,z) -> (z,x,y).
        rotation = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
        origin = np.array([3., -2., 4.])
        q = self.tool.query(origin+rotation@[.47, .1, -.04], origin, rotation)
        self.close(q['signed_distance_m'], .05)
        self.close(q['closest_points_m'], [origin+rotation@[.5, .1, 0.]])
        self.close(q['normal'], [-.8, -.6, 0.])

    def test_13_rotating_bore_point_velocity_and_normal_power_location(self):
        q = self.tool.query([.47, .1, -.04], translation_velocity=[.3, -.2, .1],
                            angular_velocity=[0., 0., 2.])
        self.close(q['wall_velocity_at_application_m_s'], [.1, .74, .1])
        self.close(q['wall_velocity_at_closest_m_s'], [[.1, .8, .1]])
        force = np.array([-3., 0., -4.])
        self.close(force@q['wall_velocity_at_application_m_s'], -.7)
        self.close(q['wall_velocity_at_closest_m_s']@force, [-.7])

    def test_14_force_matches_independent_potential_formula_gradient(self):
        # Literal corner distance: sqrt((x-.5)^2+z^2); y does not enter.
        point = np.array([.47, .1, -.04]); h = 1e-7
        def reference_U(x):
            d = np.sqrt((x[0]-.5)**2+x[2]**2)
            return 50*(.1-d)**2
        gradient = np.array([(reference_U(point+h*e)-reference_U(point-h*e))/(2*h)
                             for e in np.eye(3)])
        p = normal_penalty(self.tool.query(point), .1, 100.)
        self.close(p['normal_force_N'], -gradient, tol=2e-7)
        self.close(p['potential_J'], reference_U(point))

    def test_15_component_potential_gradient_on_all_finite_features(self):
        h = 1e-7
        for point in ([.45, .1, .4], [.8, .2, -.04], [.47, .1, -.04],
                      [.47, .1, 1.04], [2.03, 0., -.04], [2.04, 0., .4],
                      [.55, .1, .4], [.53, .54, .5]):
            point = np.array(point)
            p = normal_penalty(self.tool.query(point), .1, 100.)
            gradient = []
            for e in np.eye(3):
                up = normal_penalty(self.tool.query(point+h*e), .1, 100.)['potential_J']
                down = normal_penalty(self.tool.query(point-h*e), .1, 100.)['potential_J']
                gradient.append((up-down)/(2*h))
            self.close(p['normal_force_N'], -np.array(gradient), tol=2e-7)

    def test_16_bore_to_mouth_to_shoulder_potential_never_disappears(self):
        # Constant-distance circular arc around real sharp lip: no rounded tool.
        for angle in np.linspace(0., np.pi/2, 33):
            point = [.5-.05*np.cos(angle), .1, -.05*np.sin(angle)]
            q = self.tool.query(point); p = normal_penalty(q, .1, 100.)
            self.close(q['signed_distance_m'], .05); self.close(p['potential_J'], .125)
            self.close(q['normal'], [-np.cos(angle), 0., -np.sin(angle)])
        # Analytic limits on both face-mask crossing coordinates.
        eps = 1e-8
        for point in ([.45, .1, eps], [.45, .1, -eps],
                      [.5-eps, .1, -.05], [.5+eps, .1, -.05]):
            q = self.tool.query(point); p = normal_penalty(q, .1, 100.)
            self.close(q['signed_distance_m'], .05, tol=1e-12)
            self.close(p['potential_J'], .125, tol=1e-12)

    def test_17_path_work_across_lip_matches_independent_end_potential(self):
        # Hand path radius r(theta)=.04+.02 theta/(pi/2); n dot x'=r'.
        # Hence exact force work = int 100*(.1-r) dr = U(.04)-U(.06)=.1 J.
        angles = np.linspace(0., np.pi/2, 1001)
        power = []
        rate = .02/(np.pi/2)
        for theta in angles:
            r = .04+rate*theta
            n = np.array([-np.cos(theta), 0., -np.sin(theta)])
            dn = np.array([np.sin(theta), 0., -np.cos(theta)])
            point = np.array([.5, .1, 0.])+r*n
            force = normal_penalty(self.tool.query(point), .1, 100.)['normal_force_N']
            power.append(force@(rate*n+r*dn))
        power = np.array(power)
        integral = .5*np.sum((power[1:]+power[:-1])*np.diff(angles))
        self.close(integral, .1, tol=2e-6)

    def test_18_rigid_time_derivative_of_potential_equals_normal_wall_power(self):
        point = np.array([.47, .1, -.04]); h = 1e-7
        def evaluated(t):
            c, s = np.cos(2*t), np.sin(2*t)
            R = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
            q = self.tool.query(point, [0.3*t, -.2*t, .1*t], R,
                                [.3, -.2, .1], [0., 0., 2.])
            return q, normal_penalty(q, .1, 100.)
        q, p = evaluated(0.)
        dUdt = (evaluated(h)[1]['potential_J']-evaluated(-h)[1]['potential_J'])/(2*h)
        self.close(dUdt, p['normal_force_N']@q['wall_velocity_at_application_m_s'], tol=2e-7)

    def test_19_millimetre_geometry_is_not_absolute_tolerance_scaled(self):
        tool = FiniteOpenSquareTool(.002, .001, .001)
        q = tool.query([.00047, .0001, -.00004])
        self.close(q['signed_distance_m'], .00005, tol=1e-15)
        self.close(q['normal'], [-.6, 0., -.8])

    def test_20_invalid_geometry_rotation_or_nonfinite_input_rejected(self):
        for args in ((0., 1., 1.), (1., 2., 1.), (2., 0., 1.),
                     (2., 1., -1.), (np.inf, 1., 1.)):
            with self.assertRaises(ValueError): FiniteOpenSquareTool(*args)
        for kwargs in ({'point': [np.nan, 0., 0.]}, {'point': [0., 0.]},
                       {'point': [0., 0., 0.], 'rotation': -np.eye(3)},
                       {'point': [0., 0., 0.], 'angular_velocity': [0., 0., np.inf]}):
            with self.assertRaises(ValueError): self.tool.query(**kwargs)
        for r, k in ((0., 100.), (.1, -1.), (.1, np.nan)):
            with self.assertRaises(ValueError): normal_penalty(self.tool.query([.4, .1, .5]), r, k)

    def test_21_labels_absent_and_queries_do_not_mutate_positions(self):
        point = np.array([.47, .1, -.04]); before = point.copy()
        self.tool.query(point)
        np.testing.assert_array_equal(point, before)
        for function in (FiniteOpenSquareTool.query, normal_penalty):
            self.assertFalse(any(word in name.lower() for name in inspect.signature(function).parameters
                                 for word in ('tracer', 'label', 'as_rs')))

    def test_22_tangent_near_lip_is_not_spurious_multiple_closest_points(self):
        # Regression for initial global-distance-tolerance false tie.  The
        # squared distance difference is only 1e-16 but the cp gap is 1e-8 m.
        for point, normal, feature in (([.45, .1, 1e-8], [-1., 0., 0.], 'BORE_X_PLUS'),
                                       ([.5+1e-8, .1, -.05], [0., 0., -1.], 'SHOULDER')):
            q = self.tool.query(point)
            self.assertEqual(q['boundary_features'], (feature,))
            self.close(q['normal'], normal)
            self.assertTrue(normal_penalty(q, .1, 100.)['active'])


if __name__ == '__main__':
    unittest.main(verbosity=2)

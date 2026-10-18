"""Independent prescribed-path contact references, never a free-body trajectory.

Preregistered before first execution: analytic impulse/work/heat errors <=1e-10
N s/J; identity and zero correction residuals <=2e-12 N s/J; transformation
comparisons <=1e-10 in those units. These are O(1) manufactured paths with
literal analytic references, not an AFSD temporal accuracy allowance. A true
circular path tests only the generic quadrature primitive; the actual contact
interface is separately tested with its exact affine/quadratic HHT paths.
"""
import copy
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import hht_curved_contact as tested
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.updated_trial import particle_hash


NORMAL_FORCE = 3.
FRICTION = .2
CAPACITY = NORMAL_FORCE*FRICTION
ABS_REFERENCE = 1e-10


def emit(name, **result):
    print('HHT_CURVED_CONTACT_EVIDENCE '+json.dumps(dict(name=name,
        classification='PRESCRIBED_PATH_COMPONENT_REFERENCE_NOT_PHYSICAL_TRAJECTORY',
        **result), sort_keys=True, allow_nan=False), flush=True)


def fixture(relative0=(0., 1., 1.), acceleration0=(0., 0., 0.),
            acceleration1=None, dt=.7, wall0=(0., 0., 0.), wall_change=(0., 0., 0.),
            rotation=None, boost=(0., 0., 0.), translated=(0., 0., 0.), bad_force=None,
            mutate_callback=False):
    """One actual ParticleArray; kinematics prescribed ONLY for component tests."""
    rotation = np.eye(3) if rotation is None else np.asarray(rotation, float)
    boost, translated = np.asarray(boost, float), np.asarray(translated, float)
    normal = rotation@np.array([1., 0., 0.])
    wall0 = rotation@np.asarray(wall0, float)+boost
    dw = rotation@np.asarray(wall_change, float)
    relative0 = rotation@np.asarray(relative0, float)
    a0 = rotation@np.asarray(acceleration0, float)
    a1 = a0.copy() if acceleration1 is None else rotation@np.asarray(acceleration1, float)
    parameters = dict(dx=.005, h_ratio=1.3, rho0=2700., width=.01,
        marker_low=-.1, marker_high=.1, bar_length=.01, actuator_depth=.001,
        temperature=400., traverse=0., feed=0., shear=35e9, bulk=70e9/3,
        cp=900., conductivity=167., room=293.15, melt=855.15,
        friction=FRICTION, plastic_heat_fraction=1., heat_partition=.5,
        contact_coulomb_relative_tolerance=.01, hht_weight=.9)
    position = rotation@np.array([0., .03, -.02])+translated
    pa = make_array(position[None, :], np.zeros(1, int), parameters)
    pa.u[:], pa.v[:], pa.w[:] = (wall0+relative0).reshape(3, 1)
    t0 = 2.3
    def boundary(private, at):
        fraction = (at-t0)/dt
        wall = wall0+fraction*dw
        velocity = np.array([private.u[0], private.v[0], private.w[0]])
        relative = velocity-wall
        tangent = relative-(relative@normal)*normal
        speed = np.linalg.norm(tangent)
        traction = np.zeros(3) if speed == 0. else -CAPACITY*tangent/speed
        if bad_force == 'assisting': traction = -traction
        if bad_force == 'undersaturated': traction *= .5
        if bad_force == 'within_old_one_percent_but_not_actual_coulomb': traction *= .995
        if bad_force == 'excess_capacity': traction *= 2.
        if bad_force == 'normal_component': traction += .1*normal
        if mutate_callback: private.x[0] += 1.
        return dict(force_definition='NONCONTACT_APPLIED_FORCE_ONLY', force_N=np.zeros((1, 3)),
            contacts=[dict(particle_index=0, contact_id='plane:1', normal=normal.copy(),
                normal_force_N=NORMAL_FORCE, tangential_force_N=traction,
                wall_velocity_m_s=wall, wall_temperature_K=410., conductance_W_K=.025,
                penetration_m=1e-6, normal_potential_J=.5*NORMAL_FORCE*1e-6,
                normal_stiffness_N_m=3e6)])
    def certificate(array, old_acc, new_acc, step, p, left, right):
        return dict(minimum_overlap_m=1e-6, maximum_overlap_m=1e-6,
                    fixed_normals=True, affine_tangential_wall_velocity=True)
    return dict(pa=pa, p=parameters, boundary=boundary, t0=t0,
                acc0=a0[None, :], acc1=a1[None, :], dt=dt,
                weight=.9, gap_certificate=certificate)


class CurvedContactReferences(unittest.TestCase):
    def run_contact(self, case):
        original = particle_hash(case['pa'])
        parameters = copy.deepcopy(case['p'])
        result = tested.integrate(**case)
        self.assertEqual(particle_hash(case['pa']), original)
        self.assertEqual(case['p'], parameters)
        for key in ('actual_impulse_Ns', 'actual_material_work_J', 'actual_wall_work_J', 'friction_raw_J'):
            self.assertTrue(np.isfinite(result[key]).all(), key)
        np.testing.assert_allclose(np.asarray(result['actual_material_work_J']),
            np.asarray(result['actual_wall_work_J'])-result['friction_raw_J'], atol=2e-12, rtol=0.)
        return result

    def assert_reference(self, result, impulse, heat, wall_work):
        np.testing.assert_allclose(result['actual_impulse_Ns'][0], impulse,
                                   atol=ABS_REFERENCE, rtol=0.)
        self.assertAlmostEqual(float(result['friction_raw_J'][0]), heat, delta=ABS_REFERENCE)
        self.assertAlmostEqual(float(result['actual_wall_work_J'][0]), wall_work, delta=ABS_REFERENCE)
        self.assertAlmostEqual(float(result['actual_material_work_J'][0]), wall_work-heat,
                               delta=ABS_REFERENCE)

    def test_01_straight_stationary_and_affine_moving_wall(self):
        h = .6; slip = np.array([0., 2., 1.])
        impulse = h*(NORMAL_FORCE*np.array([1., 0., 0.])-CAPACITY*slip/np.linalg.norm(slip))
        heat = CAPACITY*h*np.linalg.norm(slip)
        for w0, dw in ((np.zeros(3), np.zeros(3)),
                       (np.array([0., .2, -.4]), np.array([0., .3, .1]))):
            with self.subTest(wall=w0.tolist(), change=dw.tolist()):
                result = self.run_contact(fixture(relative0=slip, acceleration0=dw/h,
                                                  dt=h, wall0=w0, wall_change=dw))
                wall_work = float(impulse@(w0+.5*dw))
                self.assert_reference(result, impulse, heat, wall_work)
                np.testing.assert_allclose(result['delta_impulse_Ns'], 0., atol=2e-12, rtol=0.)
                np.testing.assert_allclose(result['delta_material_work_J'], 0., atol=2e-12, rtol=0.)
        emit('straight_and_affine_wall', expected_impulse_Ns=impulse.tolist(), expected_raw_heat_J=float(heat))

    def test_02_both_components_turn_on_exact_hht_affine_path(self):
        # u(s)=(0,1-s,1+s). Unit traction has TWO nonconstant components.
        # |u|=sqrt(2)*sqrt(1+s²); these integrals use literal log/sqrt only.
        h = .7; wall = np.array([0., .2, -.3])
        a = math.log(1.+math.sqrt(2.))
        impulse = np.array([NORMAL_FORCE*h,
            -CAPACITY*h/math.sqrt(2.)*(a-math.sqrt(2.)+1.),
            -CAPACITY*h/math.sqrt(2.)*(a+math.sqrt(2.)-1.)])
        heat = CAPACITY*h*(1.+a/math.sqrt(2.))
        result = self.run_contact(fixture(acceleration0=(0., -1./h, 1./h), dt=h, wall0=wall))
        self.assert_reference(result, impulse, heat, float(wall@impulse))
        self.assertGreater(float(np.linalg.norm(result['delta_impulse_Ns'])), 1e-3)
        emit('affine_vector_turn_analytic', expected_impulse_Ns=impulse.tolist(),
            actual_impulse_Ns=result['actual_impulse_Ns'][0].tolist(), expected_raw_heat_J=heat,
            actual_raw_heat_J=float(result['friction_raw_J'][0]))

    def test_03_true_circular_velocity_generic_quadrature_only(self):
        # This is a prescribed analytic trajectory for the generic quadrature
        # primitive, not a false claim that a circle is the polynomial HHT path.
        h, speed, first, angle = .7, 1.7, .2, 1.1
        wall = np.array([.2, -.4, .3])
        def evaluator(s):
            theta = first+angle*s
            slip = speed*np.array([0., math.cos(theta), math.sin(theta)])
            tangent_force = -CAPACITY*slip/speed
            force = np.array([NORMAL_FORCE, 0., 0.])+tangent_force
            return np.r_[force, force@(wall+slip), force@wall, -tangent_force@slip]
        result = tested.refined_gauss(evaluator, h)
        end = first+angle
        impulse = np.array([NORMAL_FORCE*h,
            -CAPACITY*h*(math.sin(end)-math.sin(first))/angle,
            CAPACITY*h*(math.cos(end)-math.cos(first))/angle])
        heat = CAPACITY*speed*h
        expected = np.r_[impulse, impulse@wall-heat, impulse@wall, heat]
        np.testing.assert_allclose(result['integral'], expected, atol=ABS_REFERENCE, rtol=0.)
        self.assertIn(result['order'], (16, 32))
        emit('true_circle_quadrature_only', expected=expected.tolist(),
             actual=result['integral'].tolist(), physical_contact_certificate=False)

    def test_04_scalar_sign_crossing_is_not_full_vector_stop(self):
        h = .8
        result = self.run_contact(fixture(relative0=(0., -.5, 1.), acceleration0=(0., 1./h, 0.), dt=h))
        asinh_half = math.log(.5+math.sqrt(1.25))
        impulse = np.array([NORMAL_FORCE*h, 0., -2*CAPACITY*h*asinh_half])
        heat = CAPACITY*h*(.5*math.sqrt(1.25)+asinh_half)
        self.assert_reference(result, impulse, heat, 0.)
        emit('sign_only_turn', minimum_analytic_vector_speed_m_s=1.,
             scalar_component_crosses_zero=True, raw_heat_J=float(result['friction_raw_J'][0]))

    def test_05_full_vector_zero_and_double_touch_reject_without_mutation(self):
        cases = [fixture(relative0=(0., -.5, 0.), acceleration0=(0., 1., 0.), dt=1.),
                 fixture(relative0=(0., .25, 0.), acceleration0=(0., 40., 0.),
                         acceleration1=(0., -80./3., 0.), dt=1.)]
        for number, case in enumerate(cases):
            with self.subTest(path=number):
                original = particle_hash(case['pa'])
                with self.assertRaises(ValueError):
                    tested.integrate(**case)
                self.assertEqual(particle_hash(case['pa']), original)

    def test_06_bad_force_negative_heat_and_callback_mutation_reject(self):
        # A 0.5% mismatch is inside the retained historical direction allowance
        # but is not the actual Coulomb law: the strengthened source-contract
        # check must reject it instead of overwriting or silently accepting it.
        for reason in ('assisting', 'undersaturated', 'excess_capacity', 'normal_component',
                       'within_old_one_percent_but_not_actual_coulomb'):
            with self.subTest(reason=reason):
                case = fixture(bad_force=reason); original = particle_hash(case['pa'])
                with self.assertRaises(ValueError):
                    tested.integrate(**case)
                self.assertEqual(particle_hash(case['pa']), original)
        case = fixture(mutate_callback=True); original = particle_hash(case['pa'])
        with self.assertRaises(ValueError):
            tested.integrate(**case)
        self.assertEqual(particle_hash(case['pa']), original)
        # The generic primitive itself must not silently clip signed inputs.
        signed = tested.refined_gauss(lambda s: np.array([-2.]), .7)
        np.testing.assert_allclose(signed['integral'], [-1.4], atol=2e-12, rtol=0.)

    def test_07_rigid_rotation_and_galilean_work_covariance(self):
        h = .7; acceleration = (0., -1./h, 1./h)
        base = self.run_contact(fixture(acceleration0=acceleration, dt=h))
        boost = np.array([.3, -.8, .4])
        changed = self.run_contact(fixture(acceleration0=acceleration, dt=h, boost=boost,
                                           translated=(.02, -.03, .07)))
        np.testing.assert_allclose(changed['actual_impulse_Ns'], base['actual_impulse_Ns'],
                                   atol=ABS_REFERENCE, rtol=0.)
        np.testing.assert_allclose(changed['friction_raw_J'], base['friction_raw_J'], atol=ABS_REFERENCE, rtol=0.)
        work_shift = float(boost@base['actual_impulse_Ns'][0])
        for key in ('actual_wall_work_J', 'actual_material_work_J'):
            np.testing.assert_allclose(changed[key], np.asarray(base[key])+work_shift,
                                       atol=ABS_REFERENCE, rtol=0.)
        rotation = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
        rotated = self.run_contact(fixture(acceleration0=acceleration, dt=h, rotation=rotation))
        np.testing.assert_allclose(rotated['actual_impulse_Ns'][0], rotation@base['actual_impulse_Ns'][0],
                                   atol=ABS_REFERENCE, rtol=0.)
        for key in ('actual_material_work_J', 'actual_wall_work_J', 'friction_raw_J'):
            np.testing.assert_allclose(rotated[key], base[key], atol=ABS_REFERENCE, rtol=0.)

    def test_08_missing_or_failed_gap_certificate_and_bad_quadrature_reject(self):
        for certificate in (None, lambda *args: dict(minimum_overlap_m=-1e-6,
                maximum_overlap_m=1e-6, fixed_normals=True, affine_tangential_wall_velocity=True)):
            case = fixture(); case['gap_certificate'] = certificate
            original = particle_hash(case['pa'])
            with self.assertRaises(ValueError):
                tested.integrate(**case)
            self.assertEqual(particle_hash(case['pa']), original)
        for h in (0., -1., float('nan')):
            with self.assertRaises(ValueError):
                tested.refined_gauss(lambda s: np.ones(3), h)
        with self.assertRaises(ValueError):
            tested.refined_gauss(lambda s: np.array([float('nan')]), .1)


if __name__ == '__main__':
    unittest.main(verbosity=2)

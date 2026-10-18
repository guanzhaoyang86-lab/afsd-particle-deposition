"""Independent references for candidate-only actual-event force derivatives."""
import copy
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import event_response as er
from afsd_backend.pysph.contact_kick import kick


def reversal(**changes):
    inputs = dict(v0=[1., 0., 0.], mass=1., other_force=[-2., 0., -1.],
        normals=[[0., 0., 1.]], normal_forces=[1.], wall_velocities=[[0., 0., 0.]],
        mu=1., dt=1., free_axes=[1., 0., 0.], prescribed_velocity=[1., 0., 0.],
        contact_ids=['face_z'])
    inputs.update(changes)
    return inputs


def smooth_two_faces():
    normals = np.array([[1., 0., 0.], [1/np.sqrt(2), 1/np.sqrt(2), 0.]])
    fn = np.array([.3, .4])
    return dict(v0=np.array([1., .4, .3]), mass=1.3,
        other_force=np.array([.2, -.1, 0.])-np.sum(fn[:, None]*normals, axis=0),
        normals=normals, normal_forces=fn, wall_velocities=np.zeros((2, 3)),
        mu=1., dt=.2, free_axes=np.array([1., 1., 0.]),
        prescribed_velocity=np.array([1., .4, .3]), contact_ids=['a', 'b'])


def independent_smooth_sensitivity(case):
    """RK4 of Coulomb ODE plus variational equation, no tested kick/response."""
    n = np.asarray(case['normals']); walls = np.asarray(case['wall_velocities'])
    projectors = np.eye(3)-n[:, :, None]*n[:, None, :]
    mobility = np.diag(case['free_axes'])/case['mass']
    total = case['other_force']+np.sum(case['normal_forces'][:, None]*n, axis=0)
    capacities = case['mu']*case['normal_forces']
    # State: velocity, derivative of velocity, integrated derivative of velocity.
    y = np.r_[case['v0'], np.zeros(18)]

    def rhs(state):
        v = state[:3]; J = state[3:12].reshape(3, 3)
        force = total.copy(); A = np.zeros((3, 3))
        for R, wall, limit in zip(projectors, walls, capacities):
            u = R@(v-wall); speed = np.linalg.norm(u)
            if speed < .05:
                raise AssertionError('independent smooth fixture approached zero slip')
            force -= limit*u/speed
            A -= limit*(R/speed-np.outer(u, u)/speed**3)
        return np.r_[mobility@force, (mobility@(A@J)+mobility).ravel(), J.ravel()]

    h = case['dt']/2048
    for _ in range(2048):
        a = rhs(y); b = rhs(y+h*a/2); c = rhs(y+h*b/2); d = rhs(y+h*c)
        y += h*(a+2*b+2*c+d)/6
    return y[12:].reshape(3, 3)/case['dt'], y[3:12].reshape(3, 3)


class ActualEventResponse(unittest.TestCase):
    def record(self, name, result):
        # External runner captures initial evidence before any repair run.
        print('EVENT_RESPONSE_EVIDENCE '+json.dumps(dict(name=name, result=er._plain(result)),
              sort_keys=True, allow_nan=False), flush=True)
        self.assertFalse(result['accepted'])
        return result

    def available(self, name, **inputs):
        result = self.record(name, er.response(**inputs))
        self.assertTrue(result['available'], result['reason'])
        for probe in result['probes']:
            self.assertTrue(probe['raw_heat_certificate_pass'])
            self.assertEqual(len(probe['result_sha256']), 64)
        return result

    def test_literal_one_dimensional_reversal(self):
        r = self.available('literal_reversal', **reversal())
        self.assertAlmostEqual(r['baseline']['mean_velocity_m_s'][0], -1/18, delta=2e-11)
        self.assertAlmostEqual(r['baseline']['velocity_m_s'][0], -2/3, delta=2e-11)
        self.assertAlmostEqual(r['baseline']['friction_raw_J'], 7/18, delta=2e-11)
        self.assertAlmostEqual(r['B'][0, 0], 19/54, delta=2e-7)
        self.assertAlmostEqual(r['C'][0, 0], 7/9, delta=2e-7)
        self.assertGreater(abs(r['B'][0, 0]-.5*r['C'][0, 0]), .03)
        self.assertEqual(r['kick_calls'], 5)
        np.testing.assert_array_equal(r['B'][1:], np.zeros((2, 3)))
        np.testing.assert_array_equal(r['C'][:, 1:], np.zeros((3, 2)))
        self.assertTrue(r['column_evidence'][0]['refinement_pass'])
        self.assertTrue(r['column_evidence'][0]['one_sided_pass'])

    def test_free_response_analytic_no_kick_with_actuator(self):
        with patch.object(er, 'kick', side_effect=AssertionError('no-contact must be analytic')):
            r = self.available('free_analytic', v0=[1., 2., -.4], mass=2.,
                other_force=[4., -8., 6.], normals=[], normal_forces=[], wall_velocities=[],
                mu=.3, dt=.25, free_axes=[1., 1., 0.], prescribed_velocity=[0., 0., -.4])
        np.testing.assert_allclose(r['B'], np.diag([.0625, .0625, 0.]), atol=1e-14, rtol=0)
        np.testing.assert_allclose(r['C'], np.diag([.125, .125, 0.]), atol=1e-14, rtol=0)
        self.assertEqual(r['kick_calls'], 0)

    def test_static_stick_zero_tangent(self):
        r = self.available('static_stick', **reversal(v0=[0., 0., 0.], mass=2.,
            other_force=[.2, .3, -1.], dt=.5, free_axes=[1., 1., 1.], prescribed_velocity=[0., 0., 0.]))
        np.testing.assert_allclose(r['B'], np.diag([0., 0., .125]), atol=1e-10, rtol=0)
        np.testing.assert_allclose(r['C'], np.diag([0., 0., .25]), atol=1e-10, rtol=0)
        self.assertEqual(r['kick_calls'], 13)

    def test_moving_wall_common_boost(self):
        reference = self.available('boost_reference', **reversal())
        boosted = self.available('boosted_reversal', **reversal(v0=[3., -1., .7],
            wall_velocities=[[2., -1., .7]], prescribed_velocity=[3., -1., .7]))
        np.testing.assert_allclose(boosted['B'], reference['B'], atol=1e-8, rtol=0)
        np.testing.assert_allclose(boosted['C'], reference['C'], atol=1e-8, rtol=0)
        self.assertAlmostEqual(boosted['baseline']['friction_raw_J'], 7/18, delta=2e-11)
        np.testing.assert_allclose(boosted['baseline']['mean_velocity_m_s'],
                                   [2-1/18, -1, .7], atol=2e-11, rtol=0)

    def test_near_zero_slip_boost_keeps_raw_heat_certificate(self):
        speed = (2.+1e-12)-2.
        r = self.available('tiny_boosted_slip', **reversal(v0=[2.+1e-12, -1., .7],
            other_force=[.2, 0., -1.], wall_velocities=[[2., -1., .7]],
            prescribed_velocity=[2.+1e-12, -1., .7]))
        self.assertAlmostEqual(r['baseline']['friction_raw_J'], speed*speed/1.6, delta=1e-38)
        for probe in r['probes']:
            self.assertGreaterEqual(probe['minimum_raw_plus_rounding_bound_J'], 0.)
            self.assertEqual(probe['reconstructed_heat_max_excess_J'], 0.)

    def test_two_face_shared_inertia_and_actuator_literal(self):
        case = reversal(v0=[0., 1., 0.], other_force=[-1., -4., -1.],
            normals=[[1., 0., 0.], [0., 0., 1.]], normal_forces=[1., 1.],
            wall_velocities=np.zeros((2, 3)), free_axes=[0., 1., 0.],
            prescribed_velocity=[0., 1., 0.], contact_ids=['x', 'z'], keep_raw_traces=True)
        r = self.available('two_face_reversal', **case)
        self.assertAlmostEqual(r['baseline']['mean_velocity_m_s'][1], -11/18, delta=2e-11)
        self.assertAlmostEqual(r['baseline']['velocity_m_s'][1], -5/3, delta=2e-11)
        self.assertAlmostEqual(r['baseline']['friction_raw_J'], 14/9, delta=2e-11)
        self.assertAlmostEqual(r['B'][1, 1], 11/27, delta=3e-7)
        self.assertAlmostEqual(r['C'][1, 1], 8/9, delta=3e-7)
        coupling = r['probes'][0]['raw_kick']['contact_mobility']
        self.assertEqual(coupling[0, 1, 1, 1], 1.)
        for field in ('normals', 'normal_forces', 'wall_velocities', 'contact_ids'):
            case[field] = case[field][::-1]
        reordered = self.available('two_face_reordered', **case)
        np.testing.assert_allclose(reordered['B'], r['B'], atol=1e-8, rtol=0)
        np.testing.assert_allclose(reordered['C'], r['C'], atol=1e-8, rtol=0)

    def test_two_free_axis_cross_coupling_independent_ode(self):
        case = smooth_two_faces()
        expected_B, expected_C = independent_smooth_sensitivity(case)
        print('INDEPENDENT_SMOOTH_REFERENCE '+json.dumps(dict(B=expected_B.tolist(), C=expected_C.tolist())))
        r = self.available('smooth_two_face', **case)
        np.testing.assert_allclose(r['B'], expected_B, atol=2e-5, rtol=0)
        np.testing.assert_allclose(r['C'], expected_C, atol=2e-5, rtol=0)
        self.assertGreater(abs(r['B'][0, 1]), 1e-5)
        self.assertGreater(abs(r['C'][0, 1]), 1e-4)
        for field in ('normals', 'normal_forces', 'wall_velocities', 'contact_ids'):
            case[field] = case[field][::-1]
        reordered = self.available('smooth_two_face_reordered', **case)
        np.testing.assert_allclose(reordered['B'], r['B'], atol=1e-8, rtol=0)
        np.testing.assert_allclose(reordered['C'], r['C'], atol=1e-8, rtol=0)

    def test_nonsmooth_cone_threshold_unavailable(self):
        r = self.record('cone_threshold', er.response(**reversal(v0=[0., 0., 0.],
            other_force=[1., 0., -1.], prescribed_velocity=[0., 0., 0.])))
        self.assertFalse(r['available'])
        self.assertEqual(r['reason'], 'ENDPOINT_ACTIVE_CONE_CHANGED')
        self.assertIsNone(r['B']); self.assertIsNone(r['C'])

    def test_journal_adapter_replay_and_detached_result(self):
        from afsd_backend.pysph.coupled_boundary import evaluate_boundary
        p = dict(dx=.2, layer=1., tool_top=2., tool_radius=2., bore_width=1., substrate_depth=1.,
            traverse=0., omega=0., base_temperature=300., tool_temperature=450., normal_stiffness=20.,
            friction=.3, contact_conductance=10., bar_length=2., actuator_depth=.25, feed=.5, release_time=10.)
        boundary = evaluate_boundary(np.array([[.45, .45, 1.5]]), np.array([[.2, .1, -.5]]),
            np.array([1.]), np.array([400.]), 0., .02, p, np.array([[1., 1., 0.]]),
            np.array([1.9]), np.array([101]))
        trace = boundary['traces'][0]; original = er._digest(trace)
        r = self.record('journal_adapter', er.response_from_trace(trace, mu=.3))
        self.assertTrue(r['available'], r['reason'])
        self.assertTrue(r['saved_replay_comparison']['passed'])
        self.assertEqual(r['saved_replay_comparison']['endpoint_max_difference_m_s'], 0.)
        r['B'][:] = 99.
        self.assertEqual(er._digest(trace), original)
        modified = copy.deepcopy(trace); modified['result']['velocity'][0] += .01
        bad = self.record('journal_mismatch', er.response_from_trace(modified, mu=.3))
        self.assertFalse(bad['available'])
        self.assertEqual(bad['reason'], 'SAVED_JOURNAL_REPLAY_MISMATCH')
        with self.assertRaisesRegex(ValueError, 'dt must equal'):
            er.response_from_trace(trace, dt=.01, mu=.3)

    def test_failure_context_is_retained_without_fallback(self):
        error = RuntimeError('manufactured old cone failure')
        error.cone_context = dict(signed_gap=-1.)
        with patch.object(er, 'kick', side_effect=error):
            r = self.record('failure_context', er.response(**reversal()))
        self.assertFalse(r['available']); self.assertEqual(r['kick_calls'], 1)
        self.assertEqual(r['probes'][0]['cone_context']['signed_gap'], -1.)
        self.assertIsNone(r['B']); self.assertIsNone(r['C'])

    def test_adaptive_interval_count_is_not_branch_signature(self):
        a = dict(events=[], endpoint_cone=[['a', 'SLIDING']],
                 endpoint_slip_vectors_m_s=np.array([[1., 0., 0.]]), interval_count=12)
        b = dict(a, interval_count=98, rejected_intervals=27)
        self.assertIsNone(er._branch_change(a, b, ['a']))


if __name__ == '__main__':
    unittest.main(verbosity=2)

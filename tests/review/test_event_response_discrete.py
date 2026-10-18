"""Actual DISCRETE map derivative checks; original continuum accuracy FAIL stays.

This suite cannot certify continuous ODE/process accuracy or accept any state.
"""
import copy
import json
import unittest

import numpy as np

import test_event_response as fixture

er = fixture.er


def independent_four_midpoint_sensitivity(case):
    """Independent implicit midpoint equations and implicit variational solve.

    Four fixed 0.05-s steps, with no tested cone, kick, adaptive trace or Jacobian
    in this reference. Discrete partition is a stated fixture assumption.
    """
    n = case['normals']; R = np.eye(3)-n[:, :, None]*n[:, None, :]
    M = np.diag(case['free_axes'])/case['mass']
    external = case['other_force']+np.sum(case['normal_forces'][:, None]*n, axis=0)
    caps = case['mu']*case['normal_forces']

    def force_jac(v):
        f = external.copy(); jac = np.zeros((3, 3))
        for projection, wall, limit in zip(R, case['wall_velocities'], caps):
            u = projection@(v-wall); speed = np.linalg.norm(u)
            if speed < .05:
                raise AssertionError('reference left smooth sliding branch')
            f -= limit*u/speed
            jac -= limit*(projection/speed-np.outer(u, u)/speed**3)
        return M@f, M@jac

    v = case['v0'].copy(); J = np.zeros((3, 3)); integrated = np.zeros((3, 3))
    h = .05
    for _ in range(4):
        w = v.copy()
        for _ in range(8):
            acceleration, A = force_jac(.5*(v+w))
            w -= np.linalg.solve(np.eye(3)-.5*h*A, w-v-h*acceleration)
        acceleration, A = force_jac(.5*(v+w))
        if np.linalg.norm(w-v-h*acceleration) > 1e-14:
            raise AssertionError('independent midpoint residual')
        following = np.linalg.solve(np.eye(3)-.5*h*A, (np.eye(3)+.5*h*A)@J+h*M)
        integrated += .5*h*(J+following)
        v = w; J = following
    return integrated/.2, J


class DiscreteEventResponse(unittest.TestCase):
    def record(self, name, result):
        print('EVENT_RESPONSE_DISCRETE_EVIDENCE '+json.dumps(dict(name=name,
            result=er._plain(result)), sort_keys=True, allow_nan=False), flush=True)
        self.assertFalse(result['accepted'])
        return result

    def test_independent_discrete_response_all_entries_and_cross_axes(self):
        case = fixture.smooth_two_faces()
        expected_B, expected_C = independent_four_midpoint_sensitivity(case)
        r = self.record('discrete_reference', er.response(**case, keep_raw_traces=True))
        self.assertTrue(r['available'], r['reason'])
        np.testing.assert_allclose(r['B'], expected_B, atol=1e-7, rtol=0)
        np.testing.assert_allclose(r['C'], expected_C, atol=1e-7, rtol=0)
        for matrix in (r['B'], r['C']):
            self.assertGreater(abs(matrix[0, 1]), 1e-4)
            self.assertGreater(abs(matrix[1, 0]), 1e-4)
            np.testing.assert_array_equal(matrix[2], np.zeros(3))
            np.testing.assert_array_equal(matrix[:, 2], np.zeros(3))
        # Explicitly check the assumed partition; never infer a general physical
        # branch from this numerical partition in the implementation.
        for probe in r['probes']:
            np.testing.assert_allclose([a['dt_s'] for a in probe['raw_kick']['intervals']],
                                       [.05]*4, atol=1e-15, rtol=0)
            self.assertTrue(probe['raw_heat_certificate_pass'])
        print('DISCRETE_REFERENCE_ERRORS '+json.dumps(dict(
            max_B_s_kg=float(np.max(np.abs(r['B']-expected_B))),
            max_C_s_kg=float(np.max(np.abs(r['C']-expected_C))),
            original_continuum_accuracy='FAIL_PRESERVED_NOT_RETESTED_HERE')))

    def test_cross_axis_response_contact_order_and_boost(self):
        case = fixture.smooth_two_faces()
        r = self.record('order_reference', er.response(**case))
        for key in ('normals', 'normal_forces', 'wall_velocities', 'contact_ids'):
            case[key] = case[key][::-1]
        boost = np.array([2., -1., .7])
        case['v0'] = case['v0']+boost
        case['prescribed_velocity'] = case['prescribed_velocity']+boost
        case['wall_velocities'] = case['wall_velocities']+boost
        changed = self.record('permuted_and_boosted', er.response(**case))
        self.assertTrue(r['available']); self.assertTrue(changed['available'], changed['reason'])
        np.testing.assert_allclose(changed['B'], r['B'], atol=1e-8, rtol=0)
        np.testing.assert_allclose(changed['C'], r['C'], atol=1e-8, rtol=0)
        self.assertAlmostEqual(changed['baseline']['friction_raw_J'],
                               r['baseline']['friction_raw_J'], delta=2e-11)

    def test_actual_trace_replay_and_force_input_purity(self):
        case = fixture.smooth_two_faces(); pristine = copy.deepcopy(case)
        physical_args = {k: v for k, v in case.items() if k != 'contact_ids'}
        saved = fixture.kick(**physical_args)
        r = self.record('discrete_replay', er.response(**case, baseline_result=saved))
        self.assertTrue(r['available'], r['reason'])
        self.assertTrue(r['saved_replay_comparison']['passed'])
        self.assertEqual(r['saved_replay_comparison']['mean_max_difference_m_s'], 0.)
        for key in case:
            np.testing.assert_array_equal(case[key], pristine[key])
        for evidence in r['column_evidence']:
            self.assertTrue(evidence['refinement_pass'])
            self.assertTrue(evidence['one_sided_pass'])

    def test_inconsistent_refinement_returns_unavailable(self):
        r = self.record('deliberately_exact_refinement_bound', er.response(**fixture.reversal(),
            derivative_rtol=0., derivative_atol_factor=0.))
        self.assertFalse(r['available'])
        self.assertEqual(r['reason'], 'CENTRAL_DIFFERENCE_REFINEMENT_INCONSISTENT')
        self.assertIsNone(r['B']); self.assertIsNone(r['C'])


if __name__ == '__main__':
    unittest.main(verbosity=2)

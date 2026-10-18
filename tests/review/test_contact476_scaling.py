"""Independent cross-scale KKT references and the complete frozen local kick."""
import itertools
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.contact_cone import solve
from afsd_backend.pysph.contact_kick import kick, ContactConvergenceError
from afsd_backend.pysph.kick_diagnostics import clean
from review_artifacts import artifact_path

CAPTURE = ROOT/'outputs/engineering_e3_20261001/dx05_cfl01/REJECTED_contact_debug.json'
REFERENCE = ROOT/'outputs/engineering_e3_20261001/contact476_reference.json'


def context():
    return json.loads(CAPTURE.read_text(encoding='utf-8'))['failed_local_context']


def rotated():
    # Fixed rational orthogonal frame, independent of any contact basis code.
    return np.array([[0., -0.8, 0.6], [1., 0., 0.], [0., 0.6, 0.8]])


class Contact476Scaling(unittest.TestCase):
    def test_frozen_decimal_reference_force_and_local_impulse(self):
        c = context()['failed_cone']
        expected = np.array(json.loads(REFERENCE.read_text())['reference_force_tangent_N'], float)
        T = np.array(c['tangent_bases']); M = np.array(c['mobility'])
        f = solve(c['base_tangent'], T, M, c['coefficient'], c['limits_N'])
        actual = np.array([T[k].T@f[k] for k in range(2)])
        error = np.linalg.norm(actual-expected)
        reference_cartesian = np.array([T[k]@expected[k] for k in range(2)])
        dv = 2*c['coefficient']*M@(f-reference_cartesian).sum(axis=0)
        report = dict(force_error_N=float(error), local_endpoint_error_m_s=float(np.linalg.norm(dv)),
                      expected_tangent_N=expected, actual_tangent_N=actual)
        artifact_path('contact476_force_reference.json').write_text(json.dumps(clean(report), indent=2)+'\n')
        self.assertLessEqual(error, 1e-10)
        self.assertLessEqual(np.linalg.norm(dv), 1e-24)

    def test_cross_scale_unique_kkt_matrix(self):
        rows = []
        T0 = np.array([[[1., 0.], [0., 0.], [0., 1.]],
                       [[0., 0.], [1., 0.], [0., 1.]]])
        for coefficient, interior, constrained, rotate, reverse in itertools.product(
                (1e-24, 1e-16, 1e-8, 1.), (False, True), (False, True),
                (False, True), (False, True)):
            Q = rotated() if rotate else np.eye(3)
            T = np.array([Q@t for t in T0])
            M = Q@np.diag([1., 0. if constrained else 1., 1.])@Q.T
            caps = np.array([.4, .7])
            expected = np.array([[.24, .32], [-.42, .56]])
            if interior:
                expected[0] *= .35
            # Construct exact KKT data BEFORE the tested solve.
            # Face 0 is uniquely determined through its positive tangent block;
            # face 1 has a strictly positive multiplier even with constrained M.
            multipliers = np.array([0. if interior else .3*coefficient, 1.])
            A = np.concatenate(list(T), axis=1)
            H = coefficient*(A.T@M@A)
            base = -(H@expected.ravel()).reshape(2, 2)-multipliers[:, None]*expected
            order = np.array([1, 0]) if reverse else np.array([0, 1])
            T, caps, base, expected = T[order], caps[order], base[order], expected[order]
            row = dict(coefficient=coefficient, interior=interior, constrained=constrained,
                       rotate=rotate, reverse=reverse)
            try:
                cartesian = solve(base, T, M, coefficient, caps)
                actual = np.array([T[k].T@cartesian[k] for k in range(2)])
                g = base+coefficient*np.array([t.T@M@cartesian.sum(axis=0) for t in T])
                fs = caps.max(); vs = max(np.linalg.norm(base), coefficient*np.linalg.norm(M, 2)*fs)
                gap = float((np.sum(g*actual)+caps@np.linalg.norm(g, axis=1))/(fs*vs))
                force_error = float(np.linalg.norm(actual-expected)/fs)
                feasible = float(np.max(np.linalg.norm(actual, axis=1)-caps)/fs)
                row.update(force_error=force_error, feasibility=feasible, scaled_gap=gap,
                           status='PASS' if force_error<=1e-10 and feasible<=2e-12 and abs(gap)<=2e-12 else 'FAIL')
            except RuntimeError as error:
                row.update(status='FAIL', error=str(error))
            rows.append(row)
        artifact_path('contact476_scale_matrix.json').write_text(json.dumps(rows, indent=2)+'\n')
        failures = [row for row in rows if row['status'] != 'PASS']
        print('SCALE_MATRIX', len(rows), 'cases', len(failures), 'failures')
        self.assertEqual(failures, [], str(failures[:6]))

    def test_frozen_full_kick_time_impulse_work_and_raw_heat(self):
        c = context(); faces = c['faces']; mu = .3
        fn = np.array(c['failed_cone']['limits_N'])/mu
        try:
            result = kick(c['velocity_before'], c['mass_kg'], c['internal_force'],
                          [a['normal'] for a in faces], fn, [a['wall_velocity'] for a in faces],
                          mu, c['dt_s'], free_axes=c['free_axes'],
                          prescribed_velocity=c['prescribed_velocity'])
        except ContactConvergenceError as error:
            artifact_path('contact476_full_kick.json').write_text(json.dumps(clean(dict(
                status='FAIL', error=str(error), cone=getattr(error, 'cone_context', None),
                interval=getattr(error, 'interval_context', None))), indent=2)+'\n')
            self.fail(str(error))
        dt = c['dt_s']; mass = c['mass_kg']
        impulse = dt*np.array(c['internal_force'])+result['normal_impulse'].sum(axis=0)+\
            result['tangential_impulse'].sum(axis=0)+result['constraint_impulse']
        np.testing.assert_allclose(mass*(result['velocity']-c['velocity_before']), impulse,
                                   atol=1e-22, rtol=1e-9)
        self.assertAlmostEqual(sum(a['dt_s'] for a in result['intervals']), dt, delta=1e-22)
        self.assertLessEqual(abs(result['kinetic_work_residual_J']), 1e-20)
        for a in result['intervals']:
            raw = -np.sum(a['tangential_impulse']*a['tangential_relative_midpoint'], axis=1)
            np.testing.assert_allclose(a['friction_heat_raw_J'], raw, atol=1e-28, rtol=1e-10)
            self.assertTrue(np.all(a['friction_heat_raw_J'] >= -a['friction_rounding_bound_J']))
        artifact_path('contact476_full_kick.json').write_text(json.dumps(clean(dict(
            status='PASS_LOCAL_IDENTITY_NOT_TRAJECTORY_REFERENCE', result=result)), indent=2)+'\n')


if __name__ == '__main__':
    unittest.main(verbosity=2)

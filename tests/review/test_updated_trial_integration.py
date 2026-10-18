"""Independent joint references for the updated trial, not AFSD trajectories.

Manufactured affine/rotation fields below are candidate probes, never prescribed
free-material motion in the real runner. References do not call the tested
kernel, material, gradient, contact, or force routine to generate answers.
"""
from pathlib import Path
import hashlib
import json
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
from run_coupled_afsd import settings, initial_particles
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.updated_trial import UpdatedTrial
from afsd_backend.pysph.coupled_boundary import evaluate_boundary
from afsd_backend.pysph.output import restore


def cloud():
    """27 genuine 3D PySPH particles, clear of walls and the axial actuator."""
    _, p, _ = settings()
    x = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                            [.0001, .0006, .0011], indexing='ij')).reshape(3, -1).T
    pa = make_array(x, np.zeros(27, dtype=int), p)
    pa.u[:] = 0.; pa.v[:] = 0.; pa.w[:] = 0.
    return pa, p, x


def independent_moment(x, h, mass, rho):
    """All ordered pairs; directly transcribed scalar radial derivative.

    This does not use NNPS, HyperbolicSpline or corrected_operators. Each
    candidate position and volume enters the independent moment matrix.
    """
    result = np.zeros((len(x), 3, 3))
    for a in range(len(x)):
        for b in range(len(x)):
            if a == b:
                continue
            displacement = x[a]-x[b]
            radius = float(np.sqrt(np.dot(displacement, displacement)))
            hab = .5*(h[a]+h[b]); q = radius/hab
            if q >= 2.:
                continue
            derivative = 3*q*q-6 if q < 1. else -3*(2-q)**2
            g = (15/(62*np.pi*hab**4))*derivative*displacement/radius
            result[a] += mass[b]/rho[b]*np.outer(x[b]-x[a], g)
    return result


def array_snapshot(pa):
    return {key: prop.get_npy_array().copy() for key, prop in pa.properties.items()}


class UpdatedTrialJointReferences(unittest.TestCase):
    def test_01_updated_affine_geometry_independent_moment_and_power(self):
        pa, p, x0 = cloud(); trial = UpdatedTrial(pa, p, 0.)
        A = np.array([[2., -3., 1.], [4., -1., .5], [.25, -.75, 1.5]])
        F = np.array([[1.012, .008, -.004], [.003, .991, .005], [0., .002, 1.006]])
        candidate = trial.candidate()
        candidate['position_m'] = x0@F.T
        xm = .5*(x0+candidate['position_m'])
        candidate['mean_velocity_m_s'] = xm@A.T
        candidate['density_kg_m3'] *= 1.+np.linspace(-.001, .001, 27)
        result = trial.evaluate(candidate, 1e-7)
        mid = result['midpoint']
        np.testing.assert_allclose(mid['velocity_gradient_per_s'],
                                   np.broadcast_to(A, (27, 3, 3)), rtol=0., atol=1e-10)
        np.testing.assert_array_equal(mid['position_m'], xm)
        rho = .5*(pa.rho+candidate['density_kg_m3'])
        moment = independent_moment(xm, pa.h, pa.m, rho)
        independent_B = np.linalg.inv(moment)
        np.testing.assert_allclose(mid['B'], independent_B, rtol=2e-13, atol=2e-13)
        original_B = np.linalg.inv(independent_moment(x0, pa.h, pa.m, pa.rho))
        self.assertGreater(float(np.max(np.abs(mid['B']-original_B))), 1e-4)
        force = result['internal_force_N']; stress = mid['stress_Pa']
        power = float(np.sum(candidate['mean_velocity_m_s']*force)+
                      np.sum((pa.m/rho)[:, None, None]*stress*A))
        scale = float(np.sum(np.linalg.norm(force, axis=1)*
                             np.linalg.norm(candidate['mean_velocity_m_s'], axis=1)))
        self.assertLessEqual(abs(power), 1e-12+1e-10*scale)
        self.assertLessEqual(np.linalg.norm(force.sum(axis=0)),
                             1e-12+1e-10*np.linalg.norm(force, axis=1).sum())
        torque = np.cross(xm, force)
        self.assertLessEqual(np.linalg.norm(torque.sum(axis=0)),
                             1e-14+1e-10*np.linalg.norm(torque, axis=1).sum())
        self.assertFalse(result['accepted']); self.assertFalse(result['restart_allowed'])

    def test_02_isotropic_continuity_independent_rational_density(self):
        pa, p, x0 = cloud(); trial = UpdatedTrial(pa, p, 0.)
        rate, dt = 10., 1e-7
        candidate = trial.candidate()
        candidate['position_m'] = x0*(1+.5*rate*dt)/(1-.5*rate*dt)
        candidate['mean_velocity_m_s'] = rate*.5*(x0+candidate['position_m'])
        expected = pa.rho*(1-1.5*rate*dt)/(1+1.5*rate*dt)
        candidate['density_kg_m3'] = expected.copy()
        result = trial.evaluate(candidate, dt)
        np.testing.assert_allclose(result['map']['density_kg_m3'], expected,
                                   rtol=2e-15, atol=0.)
        np.testing.assert_allclose(result['midpoint']['velocity_gradient_per_s'],
                                   np.broadcast_to(rate*np.eye(3), (27, 3, 3)),
                                   rtol=0., atol=1e-10)
        # Subtraction rho1-rho0 divides rounding by dt; this is an explicit
        # arithmetic budget, not a relaxed physical density tolerance.
        arithmetic = 8*np.finfo(float).eps*p['rho0']/dt
        self.assertLessEqual(np.max(np.abs(result['continuity_residual_kg_m3_s'])), arithmetic)
        np.testing.assert_array_equal(pa.rho, np.full(27, p['rho0']))

    def test_03_nonzero_stress_rotation_independent_sine_cosine_reference(self):
        pa, p, x0 = cloud()
        S0 = np.array([[900., 200., -100.], [200., -300., 150.], [-100., 150., -600.]])
        pa.S[:] = np.tile(S0.ravel(), 27)
        trial = UpdatedTrial(pa, p, 0.)
        angle, dt = .01, 1e-5
        c, s = np.cos(angle), np.sin(angle)
        geometry_Q = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
        candidate = trial.candidate(); candidate['position_m'] = x0@geometry_Q.T
        # The geometrical midpoint finite difference gives spin
        # omega=2*tan(angle/2)/dt, not exactly angle/dt.
        numerical_angle = 2*np.tan(angle/2)
        skew = np.array([[0., -numerical_angle/dt, 0.],
                         [numerical_angle/dt, 0., 0.], [0., 0., 0.]])
        candidate['mean_velocity_m_s'] = .5*(x0+candidate['position_m'])@skew.T
        # Generate the analytic skew field directly: subtracting almost equal
        # coordinates first creates a tiny symmetric strain amplified by G.
        # Keep the geometric Cayley relation independently checked in meters.
        np.testing.assert_allclose(dt*candidate['mean_velocity_m_s'],
                                   candidate['position_m']-x0, rtol=0., atol=3e-19)
        cn, sn = np.cos(numerical_angle), np.sin(numerical_angle)
        stress_Q = np.array([[cn, -sn, 0.], [sn, cn, 0.], [0., 0., 1.]])
        expected = stress_Q@S0@stress_Q.T
        result = trial.evaluate(candidate, dt)
        np.testing.assert_allclose(result['material']['deviatoric_stress_Pa'],
                                   np.broadcast_to(expected, (27, 3, 3)), rtol=0., atol=5e-7)
        np.testing.assert_array_equal(result['material']['ep'], np.zeros(27))
        np.testing.assert_allclose(result['midpoint']['velocity_gradient_per_s'],
            np.broadcast_to(skew, (27, 3, 3)), rtol=0., atol=1e-10)
        continuous = geometry_Q@S0@geometry_Q.T
        self.assertGreater(np.linalg.norm(expected-continuous), 1e-5)
        np.testing.assert_array_equal(pa.S.reshape(-1, 3, 3), np.broadcast_to(S0, (27, 3, 3)))

    def test_04_real_inventory_labels_permutation_repeat_and_rejection_are_passive(self):
        _, p, _ = settings(); a = initial_particles(p); b = initial_particles(p)
        b.tracer_type[:] = 99-b.tracer_type; b.wholebar_tag[:] = 71
        b.as_rs_tag[:] = -23; b.actuator[:] = 0
        order = np.random.default_rng(811).permutation(342).astype(np.uint32)
        b = b.extract_particles(order)
        snapshot = array_snapshot(a)
        first, second = UpdatedTrial(a, p, 0.), UpdatedTrial(b, p, 0.)
        baseline = first.accepted_hash()
        reference = first.evaluate(first.candidate(), 1e-9)
        permuted = second.evaluate(second.candidate(), 1e-9)
        sort = np.argsort(permuted['particle_id'])
        for key in reference['map']:
            np.testing.assert_allclose(reference['map'][key], permuted['map'][key][sort],
                                       rtol=1e-13, atol=1e-12)
        np.testing.assert_allclose(reference['material']['deviatoric_stress_Pa'],
                                   permuted['material']['deviatoric_stress_Pa'][sort],
                                   rtol=1e-13, atol=1e-6)
        bad = first.candidate(); bad['temperature_K'][0] = p['melt']
        with self.assertRaisesRegex(ValueError, 'thermal solid-domain'):
            first.evaluate(bad, 1e-9)
        repeated = first.evaluate(first.candidate(), 1e-9)
        for key in reference['map']:
            np.testing.assert_array_equal(reference['map'][key], repeated['map'][key])
        self.assertEqual(first.accepted_hash(), baseline)
        for key, values in snapshot.items():
            np.testing.assert_array_equal(a.properties[key].get_npy_array(), values)
        self.assertEqual(np.linalg.matrix_rank(reference['midpoint']['position_m']-
                         reference['midpoint']['position_m'].mean(axis=0)), 3)
        self.assertEqual(np.count_nonzero(a.component_id == 0), 192)
        self.assertEqual(np.count_nonzero(a.component_id == 1), 150)

    def test_05_true_reversal_boundary_average_displacement_and_heat(self):
        # Genuine old finite-base geometry: rad=.1, gap=.05, kn20 => Fn1.
        # Independent tangent reference: u=1-3t until 1/3, then -(t-1/3).
        # The added z force balances the base normal and is NOT an actuator.
        p = dict(dx=.2, layer=1., tool_top=2., tool_radius=2., bore_width=1.,
                 substrate_depth=1., traverse=0., omega=0., base_temperature=400.,
                 tool_temperature=400., normal_stiffness=20., friction=1.,
                 contact_conductance=10., bar_length=2., actuator_depth=.25,
                 feed=.5, release_time=10.)
        r = evaluate_boundary(np.array([[0., 0., -.95]]), np.array([[1., 0., 0.]]),
            np.ones(1), np.array([400.]), 0., 1., p,
            np.array([[-2., 0., -1.]]), np.zeros(1), np.array([611]))
        np.testing.assert_allclose(r['velocity'], [[-2/3, 0., 0.]], rtol=0., atol=1e-10)
        np.testing.assert_allclose(r['work_velocity'], [[-1/18, 0., 0.]], rtol=0., atol=1e-10)
        self.assertAlmostEqual(r['ledger']['friction_raw_J'], 7/18, delta=1e-10)
        self.assertAlmostEqual(r['ledger']['internal_work_J'], 1/9, delta=1e-10)
        self.assertAlmostEqual(r['ledger']['delta_kinetic_J'], -5/18, delta=1e-10)
        self.assertAlmostEqual(r['ledger']['kinetic_work_residual_J'], 0., delta=1e-10)
        self.assertAlmostEqual(r['traces'][0]['result']['events'][0]['time_s'], 1/3, delta=1e-10)
        j_t = r['traces'][0]['result']['tangential_impulse'].sum(axis=0)
        endpoint_surrogate = -float(j_t@(.5*(np.array([1., 0., 0.])+r['velocity'][0])))
        self.assertAlmostEqual(endpoint_surrogate, -1/18, delta=1e-10)
        self.assertLess(endpoint_surrogate, 0.)
        np.testing.assert_allclose(r['ledger']['actuator_impulse_Ns'], 0., rtol=0., atol=1e-15)

    def test_06_true_plastic_checkpoint_sources_identity_and_history_not_zeroed(self):
        case = ROOT/'outputs/coupled_afsd_20261002/short_07us'
        manifest = json.loads((case/'manifest.json').read_text(encoding='utf-8'))
        path = case/'plastic_accepted.npz'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         manifest['artifact_sha256']['plastic_accepted.npz'])
        for name, digest in manifest['source_sha256'].items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(), digest, name)
        pa, state = restore(path)
        self.assertEqual(state['count'], 62); self.assertTrue(state['accepted'])
        self.assertTrue(state['restart_allowed'])
        self.assertEqual(state['source_hashes'], manifest['source_sha256'])
        self.assertEqual(state['configuration_hash'], manifest['configuration_sha256'])
        _, p, _ = settings(); self.assertEqual(state['parameters'], p)
        self.assertGreater(float(pa.ep.max()), 0.); self.assertGreater(float(np.linalg.norm(pa.S)), 0.)
        self.assertGreater(float(np.linalg.norm(np.column_stack((pa.u, pa.v, pa.w)))), 0.)
        self.assertGreater(float(np.max(np.abs(pa.temperature-p['temperature']))), 0.)
        saved = array_snapshot(pa); trial = UpdatedTrial(pa, p, state['t'])
        candidate = trial.candidate()
        np.testing.assert_array_equal(candidate['velocity_m_s'], np.column_stack((pa.u, pa.v, pa.w)))
        np.testing.assert_array_equal(candidate['temperature_K'], pa.temperature)
        self.assertEqual(trial.material_scope, 'UNCHANGED_FKS_APPROXIMATE_PARAMETERS')
        for key, values in saved.items():
            np.testing.assert_array_equal(pa.properties[key].get_npy_array(), values)
        print('REAL_CHECKPOINT_REFERENCE '+json.dumps(dict(status='VERIFIED_IMPORT_ONLY_NOT_ADVANCED',
            step=state['count'], time_s=state['t'], particles=len(pa.x), ep_max=float(pa.ep.max()),
            source_count=len(manifest['source_sha256']), checkpoint_sha256=manifest['artifact_sha256'][path.name])),
            flush=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Bounded actual-PySPH transaction tests; no AFSD or coupled accuracy claim."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.mechanics import positions, velocities
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.hht_coupled import (initialize, advance, restart_at_control_event,
    restore_checkpoint, CoupledStageError, FORCE_DEFINITION, density_pressure_correction)


def fixture():
    x = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
        [-.0005, 0., .0005], indexing='ij')).reshape(3, -1).T
    p = dict(dx=.0005, h_ratio=1.3, rho0=2700., width=.0015, marker_low=-.001,
        marker_high=.001, bar_length=.01, actuator_depth=.0005, temperature=400.,
        traverse=0., feed=0., shear=35e9, bulk=70e9/3, cp=1000., conductivity=120.,
        room=293.15, melt=855.15, friction=.2, plastic_heat_fraction=.9,
        heat_partition=.5, contact_coulomb_relative_tolerance=.01)
    return make_array(x, np.zeros(27, int), p), p


def empty(pa, at):
    return dict(force_definition=FORCE_DEFINITION, force_N=np.zeros((len(pa.x), 3)), contacts=[])


def body(acceleration):
    def callback(pa, at):
        return dict(force_definition=FORCE_DEFINITION,
            force_N=pa.m[:, None]*np.asarray(acceleration), contacts=[])
    return callback


def contact(pa, at):
    rows = []
    for index in range(9):
        rows.append(dict(particle_index=index, contact_id='left:'+str(int(pa.particle_id[index])),
            normal=np.array([1., 0., 0.]), normal_force_N=.1,
            tangential_force_N=np.array([0., 0., -.02]), wall_velocity_m_s=np.zeros(3),
            wall_temperature_K=450., conductance_W_K=.01, penetration_m=1e-8,
            normal_potential_J=.5e7*1e-16, normal_stiffness_N_m=1e7))
    return dict(force_definition=FORCE_DEFINITION, force_N=np.zeros((len(pa.x), 3)), contacts=rows)


def heat_reference(pa, p, dt, qfriction, wall, bath):
    """Independent scalar all-pairs radial polynomial, no tested heat calls."""
    x = positions(pa); count = len(x); lap = np.zeros((count, count))
    for i in range(count):
        for j in range(i+1, count):
            radius = np.linalg.norm(x[i]-x[j]); smoothing = .5*(pa.h[i]+pa.h[j]); q = radius/smoothing
            if q >= 2: continue
            derivative = 3*q*q-6 if q < 1 else -3*(2-q)**2
            radial = 15/(62*np.pi*smoothing**4)*derivative
            conductance = (pa.m[i]/pa.rho[i])*(pa.m[j]/pa.rho[j])*2*p['conductivity']*(-radial/radius)
            lap[i, i] += conductance; lap[j, j] += conductance
            lap[i, j] -= conductance; lap[j, i] -= conductance
    capacity = np.diag(pa.m*p['cp']); operator = lap+np.diag(wall)
    return pa.temperature+np.linalg.solve(capacity+dt*operator,
        p['heat_partition']*qfriction+dt*(bath-operator@pa.temperature))


class CoupledTests(unittest.TestCase):
    def test_01_zero_27node_equilibrium_and_full_checkpoint(self):
        pa, p = fixture(); original = particle_hash(pa)
        state = initialize(pa, p, empty)
        after, next_state, report = advance(pa, state, p, empty, 1e-8)
        np.testing.assert_array_equal(positions(after), positions(pa))
        np.testing.assert_array_equal(velocities(after), 0.)
        np.testing.assert_array_equal(after.rho, pa.rho)
        np.testing.assert_array_equal(after.temperature, pa.temperature)
        np.testing.assert_array_equal(after.S, pa.S)
        self.assertEqual(report['evaluations'], 1)
        self.assertEqual(particle_hash(pa), original)
        restored, restored_state = restore_checkpoint(pa, next_state, p)
        self.assertEqual(particle_hash(restored), particle_hash(after))
        self.assertEqual(restored_state['state_hash'], next_state['state_hash'])
        altered = copy.deepcopy(next_state); altered['acceleration_m_s2'][0, 0] = 1.
        with self.assertRaisesRegex(ValueError, 'semantic checkpoint hash'):
            restore_checkpoint(pa, altered, p)

    def test_02_constant_body_acceleration_independent_analytic_motion(self):
        pa, p = fixture(); pa.u[:] = .03; pa.v[:] = -.02; pa.w[:] = .01
        acceleration = np.array([2., -3., 5.]); callback = body(acceleration); dt = 2e-8
        state = initialize(pa, p, callback)
        after, _, report = advance(pa, state, p, callback, dt)
        np.testing.assert_allclose(positions(after), positions(pa)+dt*velocities(pa)+.5*dt*dt*acceleration,
            atol=2e-18, rtol=0.)
        np.testing.assert_allclose(velocities(after), velocities(pa)+dt*acceleration, atol=2e-14, rtol=0.)
        np.testing.assert_allclose(after.S, 0., atol=2e-5, rtol=0.)
        self.assertLess(max(report['residuals'].values()), 1e-9)

    def test_03_actual_thermal_graph_independent_reference(self):
        pa, p = fixture(); pa.temperature[:] = np.linspace(400., 430., 27)
        dt = 1e-4; z = np.zeros(27)
        reference = heat_reference(pa, p, dt, z, z, z)
        state = initialize(pa, p, empty)
        after, _, report = advance(pa, state, p, empty, dt)
        np.testing.assert_allclose(after.temperature, reference, atol=2e-12, rtol=0.)
        self.assertLess(abs(report['thermal']['ledger']['mc_delta_T_residual_J']), 1e-12)
        np.testing.assert_array_equal(positions(after), positions(pa))

    def test_04_contact_heat_same_history_and_passive_tags(self):
        pa, p = fixture(); pa.w[:] = .02; pa.S[:] = np.tile(np.diag([1000., -500., -500.]).ravel(), 27)
        original = particle_hash(pa); state = initialize(pa, p, contact); state_hash = state['state_hash']
        after, next_state, report = advance(pa, state, p, contact, 1e-9)
        again, same_state, again_report = advance(pa, state, p, contact, 1e-9)
        self.assertEqual(particle_hash(after), particle_hash(again))
        self.assertEqual(next_state['state_hash'], same_state['state_hash'])
        self.assertEqual(particle_hash(pa), original); self.assertEqual(state['state_hash'], state_hash)
        self.assertGreater(float(report['contact']['friction_raw_J'].sum()), 0.)
        self.assertGreater(float((after.temperature-pa.temperature).max()), 0.)
        self.assertTrue(all(row['whole_path_passive'] for row in report['contact']['contacts']))
        self.assertFalse(report['algorithmic_terms_are_material_heat'])
        tagged = pa.extract_particles(np.arange(27, dtype=np.uint32))
        tagged.tracer_type[:] = np.arange(27) % 3; tagged.wholebar_tag[:] = np.arange(27) % 4
        tagged_state = initialize(tagged, p, contact)
        tagged_after, _, _ = advance(tagged, tagged_state, p, contact, 1e-9)
        for key in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'temperature', 'S', 'ep'):
            np.testing.assert_array_equal(getattr(tagged_after, key), getattr(after, key))
        self.assertLess(float(np.max(np.abs(report['kinetic']['kinetic_identity_residual_J']))), 1e-18)

    def test_05_control_reset_preserves_full_physical_state_and_restart(self):
        pa, p = fixture(); first = body([0., 0., 2.]); second = body([0., 0., -3.])
        state = initialize(pa, p, first)
        after, accepted, _ = advance(pa, state, p, first, 1e-8)
        original = particle_hash(after)
        with self.assertRaisesRegex(CoupledStageError, 'explicit event restart'):
            advance(after, accepted, p, second, 1e-8)
        restarted, event = restart_at_control_event(after, accepted, p, second)
        self.assertEqual(particle_hash(after), original)
        np.testing.assert_array_equal(restarted['acceleration_m_s2'], np.tile([0., 0., -3.], (27, 1)))
        self.assertEqual(event['physical_energy_jump_J'], 0.)
        self.assertFalse(event['full_algorithmic_energy_jump_certified'])
        self.assertFalse(event['numerical_terms_are_material_heat'])
        restored, reload_state = restore_checkpoint(pa, restarted, p)
        onward, state1, _ = advance(after, restarted, p, second, 1e-8)
        re_onward, state2, _ = advance(restored, reload_state, p, second, 1e-8)
        self.assertEqual(particle_hash(onward), particle_hash(re_onward))
        self.assertEqual(state1['state_hash'], state2['state_hash'])

    def test_06_invalid_force_definition_and_plastic_scope_reject_atomically(self):
        pa, p = fixture()
        def ambiguous(array, at): return dict(force_N=np.zeros((27, 3)), contacts=[])
        with self.assertRaisesRegex(ValueError, 'NONCONTACT_APPLIED_FORCE_ONLY'):
            initialize(pa, p, ambiguous)
        pa.S[:] = np.tile(np.diag([1e9, -5e8, -5e8]).ravel(), 27)
        original = particle_hash(pa); state = initialize(pa, p, empty)
        with self.assertRaisesRegex(CoupledStageError, 'plastic response outside'):
            advance(pa, state, p, empty, 1e-12)
        self.assertEqual(particle_hash(pa), original)

    def test_07_density_pressure_correction_matches_independent_full_block(self):
        A = np.array([[4., 1.], [1., 3.]])
        C = np.array([[2.], [-1.]])
        L = np.array([[-.1, .2]])
        RA = np.array([3., -2.]); Rrho = np.array([.5]); weight = .9
        # Independent uneliminated joint system, not the tested Schur formula.
        full = np.block([[A, -weight*C], [-L, np.eye(1)]])
        reference = np.linalg.solve(full, -np.r_[RA, Rrho])
        da, drho, certificate = density_pressure_correction(A, C, L, RA, Rrho, weight)
        np.testing.assert_allclose(np.r_[da, drho], reference, atol=2e-15, rtol=0.)
        self.assertLess(np.max(np.abs(certificate['acceleration_block_residual_N'])), 2e-15)
        # C=0: acceleration is independent, density still follows its response.
        da, drho, _ = density_pressure_correction(np.diag([2., 4.]), np.zeros((2, 1)),
            np.array([[.3, -.2]]), np.array([2., -8.]), np.array([.7]), .9)
        np.testing.assert_allclose(da, [-1., 2.], atol=0., rtol=0.)
        np.testing.assert_allclose(drho, [-1.4], atol=2e-16, rtol=0.)

    def test_08_saved_coarse_step6_converges_same_root_gates_and_budget(self):
        sys.path.insert(0, str(ROOT/'scripts'))
        from probe_hht_coupled import setup, boundary_for
        from afsd_backend.pysph.output import restore
        _, p = setup()
        source = ROOT/'outputs/process_hht_20261003/sweep_initial/grid_008/stage_005/accepted.npz'
        pa, state = restore(source)
        pa, state = restore_checkpoint(pa, state, p)
        original = particle_hash(pa); original_state = state['state_hash']
        after, next_state, report = advance(pa, state, p, boundary_for(p, True), p['end_time']/8)
        self.assertLessEqual(report['evaluations'], 24)
        self.assertLessEqual(report['residuals']['velocity_residual_m_s'], 1e-9)
        self.assertLessEqual(report['residuals']['position_residual_m'], 1e-12)
        self.assertLessEqual(report['residuals']['density_relative_residual'], 1e-10)
        self.assertLessEqual(report['residuals']['temperature_residual_K'], 5e-12)
        self.assertEqual(particle_hash(pa), original); self.assertEqual(state['state_hash'], original_state)
        self.assertAlmostEqual(next_state['time_s'], 1.5e-6, delta=1e-21)
        print('SAVED_STEP6_REPAIR1 '+repr(dict(evaluations=report['evaluations'],
            residuals=report['residuals'], elapsed_s=report['elapsed_s'])), flush=True)

    def test_09_rejected_candidates_keep_raw_residuals_and_cannot_restart(self):
        pa, p = fixture(); pa.w[:] = .02
        state = initialize(pa, p, contact); original = particle_hash(pa)
        # Failure injection only: restore the real 24-evaluation bound afterward.
        with patch('afsd_backend.pysph.hht_coupled.MAX_EVALUATIONS', 1):
            with self.assertRaises(CoupledStageError) as caught:
                advance(pa, state, p, contact, 1e-9)
        raw = caught.exception.diagnostics['last_evaluated_candidate']
        self.assertFalse(raw['accepted']); self.assertFalse(raw['restart_allowed'])
        self.assertEqual(raw['acceleration_m_s2'].shape, (27, 3))
        self.assertEqual(raw['equilibrium_residual_N'].shape, (27, 3))
        self.assertEqual(raw['density_residual_kg_m3'].shape, (27,))
        self.assertEqual(raw['temperature_residual_K'].shape, (27,))
        self.assertEqual(len(raw['physical_candidate_hash']), 64)
        self.assertEqual(particle_hash(pa), original)


if __name__ == '__main__':
    unittest.main(verbosity=2)

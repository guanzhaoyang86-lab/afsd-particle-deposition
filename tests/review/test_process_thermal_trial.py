"""Independent heat references and genuine PySPH joint-interface fixtures.

No AFSD trajectory is advanced. No tested routine generates reference values.
Manufactured velocities/material points below are tests, not free-body controls.
"""
import copy
import math
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'scripts'))
from run_coupled_afsd import settings
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.updated_trial import UpdatedTrial, particle_hash, FIELDS
from afsd_backend.pysph.process_thermal_trial import (
    ProcessThermalTrial, CLASSIFICATION, acceptance, certify_boundary_heat, certify_process_thermal)
from afsd_backend.pysph.implicit_thermal import advance_implicit_thermal
from afsd_backend.pysph.coupled_thermal import advance_thermal
from afsd_backend.pysph.coupled_boundary import evaluate_boundary


def cloud(temperature=400., shift=0.):
    _, p, _ = settings(); p = dict(p, temperature=temperature)
    x = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                            [.0001+shift, .0006+shift, .0011+shift], indexing='ij')).reshape(3, -1).T
    return make_array(x, np.zeros(27, dtype=int), p), p


def independent_heat_matrix(pa, p):
    """All pairs, scalar paper radial polynomial; no NNPS/kernel/heat calls."""
    x = np.column_stack((pa.x, pa.y, pa.z)); L = np.zeros((len(x), len(x)))
    for a in range(len(x)):
        for b in range(a+1, len(x)):
            radius = np.linalg.norm(x[a]-x[b]); h = .5*(pa.h[a]+pa.h[b]); q = radius/h
            if q >= 2.: continue
            derivative = 3*q*q-6 if q < 1 else -3*(2-q)**2
            radial = 15/(62*np.pi*h**4)*derivative
            gij = (pa.m[a]/pa.rho[a])*(pa.m[b]/pa.rho[b])*2*p['conductivity']*(-radial/radius)
            L[a, a] += gij; L[b, b] += gij; L[a, b] -= gij; L[b, a] -= gij
    return L


def reverse_boundary():
    p = dict(dx=.2, layer=1., tool_top=2., tool_radius=2., bore_width=1.,
        substrate_depth=1., traverse=0., omega=0., base_temperature=500.,
        tool_temperature=600., normal_stiffness=20., friction=1., contact_conductance=10.,
        bar_length=2., actuator_depth=.25, feed=.5, release_time=10.)
    b = evaluate_boundary(np.array([[0., 0., -.95]]), np.array([[1., 0., 0.]]),
        np.ones(1), np.array([400.]), 0., 1., p, np.array([[-2., 0., -1.]]),
        np.zeros(1), np.array([611]))
    return b, p


class ProcessThermalTrialTests(unittest.TestCase):
    def test_01_two_node_hand_reference_and_legacy_explicit_limit(self):
        # rho=m=cp=1, k=.5, separation 1, grad=[1,0,0] => gij=1 W/K.
        m = np.ones(2); rho = m.copy(); T = np.array([300., 400.]); x = np.array([[0.,0.,0.], [1.,0.,0.]])
        i = np.array([0]); j = np.array([1]); g = np.array([[1.,0.,0.]]); z = np.zeros(2)
        r = advance_implicit_thermal(m, rho, T, .5, x, i, j, g, 1., 1., z, z, z, z, z, .9, .5, 293., 855.)
        np.testing.assert_allclose(r['temperature'], [1000/3, 1100/3], rtol=0., atol=2e-12)
        self.assertEqual(r['former_explicit_dt_bound_s'], .25)
        with self.assertRaisesRegex(ValueError, 'explicit conduction time-step bound'):
            advance_thermal(m, rho, T, .5, x, i, j, g, 1., 1., z, z, z, z, .9, .5, 293., 855.)

    def test_02_single_node_wall_hand_reference(self):
        r = advance_implicit_thermal(np.ones(1), np.ones(1), np.array([300.]), 1., np.zeros((1,3)),
            np.array([], int), np.array([], int), np.zeros((0,3)), .5, 2., np.array([4.]),
            np.zeros(1), np.zeros(1), np.array([3.]), np.array([1500.]), 1., 0., 293., 855.)
        np.testing.assert_allclose(r['temperature'], [1354/3.5], rtol=0., atol=2e-12)
        self.assertAlmostEqual(r['ledger']['contact_in_J'], .5*(1500-3*1354/3.5), delta=2e-12)

    def test_03_zero_27point_equilibrium_real_full_acceptance(self):
        pa, p = cloud(); trial = ProcessThermalTrial(pa, p, 0.); c = trial.candidate()
        e = trial.evaluate(c, 1e-9)
        for key in FIELDS: np.testing.assert_array_equal(e['residual'][key], 0.)
        result = acceptance(pa, p, trial, c, e)
        self.assertEqual(result['variant'], CLASSIFICATION)
        self.assertFalse(e['accepted']); self.assertFalse(e['restart_allowed'])
        self.assertFalse(result['production_authorized'])

    def test_04_27point_nonuniform_temperature_independent_graph_solution(self):
        pa, p = cloud(); pa.temperature[:] = np.linspace(400., 450., 27)
        dt = 1e-3; L = independent_heat_matrix(pa, p); C = np.diag(pa.m*p['cp'])
        reference = pa.temperature+np.linalg.solve(C+dt*L, -dt*L@pa.temperature)
        trial = ProcessThermalTrial(pa, p, 0.); c = trial.candidate(); c['temperature_K'] = reference.copy()
        e = trial.evaluate(c, dt)
        np.testing.assert_allclose(e['map']['temperature_K'], reference, rtol=0., atol=2e-12)
        a = acceptance(pa, p, trial, c, e)
        self.assertLess(e['thermal']['former_explicit_dt_bound_s'], dt)
        self.assertLess(abs(a['thermal_residual_J']), 1e-12)
        with self.assertRaisesRegex(ValueError, 'explicit conduction time-step bound'):
            UpdatedTrial(pa, p, 0.).evaluate(c, dt)

    def test_05_same_history_repeat_and_detached_candidate(self):
        pa, p = cloud(); pa.ep[:] = .01; pa.S[:] = np.tile(np.diag([200.,-100.,-100.]).ravel(), 27)
        original = particle_hash(pa); trial = ProcessThermalTrial(pa, p, 0.); c = trial.candidate(); saved = copy.deepcopy(c)
        a = trial.evaluate(c, 1e-10); b = trial.evaluate(c, 1e-10)
        for key in FIELDS:
            np.testing.assert_array_equal(a['map'][key], b['map'][key]); np.testing.assert_array_equal(c[key], saved[key])
        self.assertEqual(trial.accepted_hash(), original); self.assertEqual(particle_hash(pa), original)
        a['material']['ep'][:] = 100.; a['map']['temperature_K'][:] = 600.
        again = trial.evaluate(c, 1e-10)
        np.testing.assert_array_equal(again['material']['ep'], .01)
        np.testing.assert_array_equal(again['map']['temperature_K'], b['map']['temperature_K'])

    def test_06_passive_tags_and_reorder_do_not_change_physics(self):
        pa, p = cloud(); pa.temperature[:] = np.linspace(400., 450., 27)
        pb = pa.extract_particles(np.arange(26,-1,-1, dtype=np.uint32)); pb.tracer_type[:] = 99
        pb.wholebar_tag[:] = -71; pb.actuator[:] = 1
        a = ProcessThermalTrial(pa, p, 0.); b = ProcessThermalTrial(pb, p, 0.)
        ea = a.evaluate(a.candidate(), 1e-6); eb = b.evaluate(b.candidate(), 1e-6); order = np.argsort(eb['particle_id'])
        for key in FIELDS: np.testing.assert_allclose(ea['map'][key], eb['map'][key][order], rtol=0., atol=2e-12)
        np.testing.assert_array_equal(ea['material']['ep'], eb['material']['ep'][order])

    def test_07_candidate_temperature_feedback_original_fks_independent_scalar_return(self):
        pa, p = cloud(723.15); S0 = np.diag([2e7,-1e7,-1e7]); pa.S[:] = np.tile(S0.ravel(), 27)
        trial = ProcessThermalTrial(pa, p, 0.); c = trial.candidate(); dt = 1e-10
        cold = trial.evaluate(c, dt); c['temperature_K'] += 10.; hot = trial.evaluate(c, dt)
        T = 728.15; ts = (T-p['room'])/(p['melt']-p['room'])
        # At zero rate, rate_floor=reference_rate=1 => lambda=1.
        theta = 1-(1+math.exp(-p['c1']*ts))**(-1/p['c2'])
        lo, hi = 0., 3e7/(3*p['shear'])
        for _ in range(90):
            dep = .5*(lo+hi); y = (p['a1']+p['a2']*math.atan(p['a3']*dep))*theta
            if 3e7-3*p['shear']*dep > y: lo = dep
            else: hi = dep
        reference_dep = .5*(lo+hi); reference_S = S0*(3e7-3*p['shear']*reference_dep)/3e7
        np.testing.assert_allclose(hot['material']['ep'], reference_dep, rtol=0., atol=1e-15)
        np.testing.assert_allclose(hot['material']['deviatoric_stress_Pa'], np.broadcast_to(reference_S,(27,3,3)), rtol=0., atol=2e-6)
        self.assertGreater(float(np.max(np.abs(hot['material']['ep']-cold['material']['ep']))), 1e-7)
        self.assertGreater(float(hot['thermal']['ledger']['plastic_heat_J']), 0.)
        np.testing.assert_array_equal(pa.ep, 0.)

    def test_08_nonthermal_map_and_event_history_equal_frozen_original(self):
        pa, p = cloud(); pa.S[:] = np.tile(np.diag([200.,-100.,-100.]).ravel(), 27)
        old = UpdatedTrial(pa, p, 0.); new = ProcessThermalTrial(pa, p, 0.)
        c = new.candidate(); a = old.evaluate(c, 1e-10); b = new.evaluate(c, 1e-10)
        for key in FIELDS: np.testing.assert_array_equal(a['map'][key], b['map'][key])
        for key in a['material']: np.testing.assert_array_equal(a['material'][key], b['material'][key])
        np.testing.assert_array_equal(a['internal_force_N'], b['internal_force_N'])
        np.testing.assert_array_equal(a['midpoint']['B'], b['midpoint']['B'])
        self.assertNotEqual(a['classification'], b['classification'])

    def test_09_real_reversal_events_and_face_reservoirs(self):
        b, p = reverse_boundary(); cert = certify_boundary_heat(b, np.array([611]), 1., p)
        self.assertAlmostEqual(cert['friction_raw_J'][0], 7/18, delta=1e-10)
        self.assertAlmostEqual(b['velocity'][0,0], -2/3, delta=1e-10)
        self.assertAlmostEqual(b['work_velocity'][0,0], -1/18, delta=1e-10)
        self.assertEqual(cert['interval_count'], 2)
        np.testing.assert_allclose(cert['wall_row_W_K'], [.4], rtol=0., atol=1e-15)
        np.testing.assert_allclose(cert['wall_weighted_temperature_W'], [200.], rtol=0., atol=1e-12)

    def test_10_inflated_roundoff_certificate_cannot_mask_negative_event(self):
        b, p = reverse_boundary(); interval = b['traces'][0]['result']['intervals'][0]
        interval['friction_rounding_bound_J'] *= 100.
        with self.assertRaisesRegex(ValueError, 'rounding certificate inflated'):
            certify_boundary_heat(b, np.array([611]), 1., p)

    def test_11_negative_event_cannot_hide_in_positive_aggregate(self):
        b, p = reverse_boundary(); b['traces'][0]['result']['intervals'][0]['friction_heat_raw_J'] *= -1.
        self.assertGreater(float(b['friction_raw_J'][0]), 0.)
        with self.assertRaisesRegex(ValueError, 'friction impulse-work identity'):
            certify_boundary_heat(b, np.array([611]), 1., p)

    def test_12_contact_face_temperature_must_match_actual_reservoir(self):
        b, p = reverse_boundary(); b['traces'][0]['contact_rows'][0]['wall_temperature_K'] += 1.
        with self.assertRaisesRegex(ValueError, 'reservoir mismatch'):
            certify_boundary_heat(b, np.array([611]), 1., p)

    def test_13_full_acceptance_does_not_bypass_old_mechanical_thresholds(self):
        pa, p = cloud(); trial = ProcessThermalTrial(pa, p, 0.); c = trial.candidate(); e = trial.evaluate(c, 1e-9)
        for key, value, message in [('power_conjugacy_residual_W', 1e-4, 'spatial power'),
                                   ('internal_net_force_N', np.ones(3), 'internal net force'),
                                   ('internal_net_torque_Nm', np.ones(3), 'internal net torque')]:
            bad = copy.deepcopy(e); bad['spatial'][key] = value
            with self.assertRaisesRegex(ValueError, message): acceptance(pa, p, trial, c, bad)
        bad = copy.deepcopy(e); bad['material_scope'] = 'MANUFACTURED_TEST_ONLY'
        with self.assertRaisesRegex(ValueError, 'manufactured material'): acceptance(pa, p, trial, c, bad)
        bad = copy.deepcopy(c); bad['mean_velocity_m_s'][:,0] = 1e-5
        with self.assertRaisesRegex(ValueError, 'global residual'): acceptance(pa, p, trial, bad, e)

    def test_14_outside_solid_temperature_rejected_not_clipped(self):
        pa, p = cloud(); trial = ProcessThermalTrial(pa, p, 0.); before = trial.accepted_hash()
        for value in [p['room']-1., p['melt']]:
            c = trial.candidate(); c['temperature_K'][0] = value
            with self.assertRaisesRegex(ValueError, 'thermal solid-domain'): trial.evaluate(c, 1e-9)
        self.assertEqual(before, trial.accepted_hash())

    def test_15_be_candidate_balance_recomputed_not_copied(self):
        pa, p = cloud(); trial = ProcessThermalTrial(pa, p, 0.); c = trial.candidate(); e = trial.evaluate(c, 1e-9)
        c['temperature_K'][0] += 1.
        with self.assertRaisesRegex(ValueError, 'candidate BE per-particle thermal balance'):
            certify_process_thermal(trial, c, e)

    def test_16_real_27point_base_heat_uses_end_temperature(self):
        # Base is z=-1.5mm, support radius=.25mm: 25um initial penalty overlap.
        # First fixture used 150um and correctly failed the unchanged 50um guard.
        pa, p = cloud(400., shift=-.001375); trial = ProcessThermalTrial(pa, p, 0.); c = trial.candidate()
        e = trial.evaluate(c, 1e-10); cert = e['thermal_boundary_certificate']; dt = e['dt_s']
        self.assertGreater(e['boundary']['ledger']['contacts'], 0)
        reference = dt*(cert['wall_weighted_temperature_W']-cert['wall_row_W_K']*e['map']['temperature_K'])
        np.testing.assert_allclose(e['thermal']['particle_energy_J']['contact_in_J'], reference, rtol=0., atol=1e-25)
        np.testing.assert_allclose(e['boundary']['contact_power_W']*dt, reference, rtol=0., atol=1e-25)
        self.assertAlmostEqual(e['boundary']['ledger']['contact_conduction_into_material_J'], float(reference.sum()), delta=1e-24)
        self.assertIn('midpoint_temperature_conduction_diagnostic_W', e['boundary'])
        self.assertFalse(e['accepted'])


if __name__ == '__main__': unittest.main(verbosity=2)

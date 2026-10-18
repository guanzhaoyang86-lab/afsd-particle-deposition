"""Pure trial regressions; manufactured paths are not accepted AFSD states."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'scripts'))
from run_coupled_afsd import settings
from afsd_backend.pysph import updated_trial as u
from afsd_backend.pysph.particles import make_array


class ConstantYield:
    def __init__(self, value): self.value = value
    def yield_stress(self, ep, rate, temperature):
        return np.full(np.broadcast(ep, rate, temperature).shape, self.value)


def cloud():
    _, p, _ = settings()
    p = dict(p, temperature=400.)
    xyz = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                               [.0001, .0006, .0011], indexing='ij')).reshape(3, -1).T
    pa = make_array(xyz, np.zeros(27, dtype=int), p)
    return pa, p


class UpdatedTrialTests(unittest.TestCase):
    def test_01_zero_equilibrium_is_not_an_accepted_commit(self):
        pa, p = cloud(); trial = u.UpdatedTrial(pa, p, 0.)
        result = trial.evaluate(trial.candidate(), 1e-9)
        for key in u.FIELDS: np.testing.assert_array_equal(result['residual'][key], 0.)
        self.assertFalse(result['accepted']); self.assertFalse(result['restart_allowed'])
        self.assertEqual(result['material_scope'], 'UNCHANGED_FKS_APPROXIMATE_PARAMETERS')
        self.assertEqual(result['boundary']['ledger']['contacts'], 0)

    def test_02_trial_and_returned_arrays_do_not_alias_or_accumulate_history(self):
        pa, p = cloud(); pa.ep[:] = .01; pa.plastic_work[:] = np.arange(27)*1e-9
        pa.S.reshape(-1, 3, 3)[:] = np.diag([2e5, -1e5, -1e5])
        before = u.particle_hash(pa); trial = u.UpdatedTrial(pa, p, 0.)
        candidate = trial.candidate(); snapshot = copy.deepcopy(candidate)
        a = trial.evaluate(candidate, 1e-10)
        a['map']['temperature_K'][:] = 500.; a['material']['ep'][:] = 100.
        b = trial.evaluate(candidate, 1e-10)
        np.testing.assert_array_equal(b['material']['ep'], .01)
        for key in u.FIELDS: np.testing.assert_array_equal(candidate[key], snapshot[key])
        self.assertEqual(before, u.particle_hash(pa)); self.assertEqual(before, trial.accepted_hash())
        pa.ep[:] = .5  # Constructor detached its accepted snapshot.
        c = trial.evaluate(candidate, 1e-10)
        np.testing.assert_array_equal(c['material']['ep'], .01)

    def test_03_candidate_shape_nonfinite_and_missing_fields_fail_without_mutation(self):
        pa, p = cloud(); trial = u.UpdatedTrial(pa, p, 0.); before = trial.accepted_hash()
        for key in u.FIELDS:
            c = trial.candidate(); c[key] = c[key][:-1]
            with self.assertRaisesRegex(ValueError, 'invalid/nonfinite'): trial.evaluate(c, 1e-9)
            c = trial.candidate(); c[key].flat[0] = np.nan
            with self.assertRaisesRegex(ValueError, 'invalid/nonfinite'): trial.evaluate(c, 1e-9)
        c = trial.candidate(); del c['mean_velocity_m_s']
        with self.assertRaisesRegex(ValueError, 'exactly five'): trial.evaluate(c, 1e-9)
        self.assertEqual(before, trial.accepted_hash())

    def test_04_domains_remain_hard_failures_not_clipping(self):
        pa, p = cloud(); trial = u.UpdatedTrial(pa, p, 0.); before = trial.accepted_hash()
        for key, value, text in [('temperature_K', p['melt'], 'thermal'),
                                 ('density_kg_m3', -1., 'nonpositive'),
                                 ('velocity_m_s', 100., 'speed')]:
            c = trial.candidate(); c[key].flat[0] = value
            with self.assertRaisesRegex(ValueError, text): trial.evaluate(c, 1e-9)
        c = trial.candidate(); c['position_m'][0, 0] += 51e-6
        with self.assertRaisesRegex(ValueError, 'displacement'): trial.evaluate(c, 1e-9)
        for dt in (0., -1., np.nan):
            with self.assertRaises(ValueError): trial.evaluate(trial.candidate(), dt)
        self.assertEqual(before, trial.accepted_hash())

    def test_05_release_discontinuity_requires_step_boundary(self):
        pa, p = cloud(); trial = u.UpdatedTrial(pa, p, p['release_time']-1e-8)
        with self.assertRaisesRegex(ValueError, 'step boundary'):
            trial.evaluate(trial.candidate(), 2e-8)

    def test_06_independent_midpoint_density_formula(self):
        pa, p = cloud(); trial = u.UpdatedTrial(pa, p, 0.); dt = 1e-9; a = 25.
        c = trial.candidate(); c['mean_velocity_m_s'] = a*c['position_m']
        c['density_kg_m3'][:] = 2700.*(1-1.5*a*dt)/(1+1.5*a*dt)
        result = trial.evaluate(c, dt)
        np.testing.assert_allclose(result['map']['density_kg_m3'], c['density_kg_m3'], rtol=0, atol=5e-12)
        # Subtracting densities in Rrho amplifies float roundoff by 1/dt.
        self.assertLess(np.max(np.abs(result['continuity_residual_kg_m3_s'])), 1e-3)
        np.testing.assert_array_equal(pa.rho, 2700.)

    def test_07_independent_j2_plastic_work_and_return_heat_separation(self):
        pa, p = cloud()
        # Explicit manufactured trial-field amplitude, never a production law.
        p = dict(p, shear=30e9, max_speed=1000., plastic_heat_fraction=.9)
        p['normal_stiffness'] = p['penalty_factor']*p['shear']*p['dx']**2/p['contact_reference_spacing']
        trial = u.UpdatedTrial(pa, p, 0., law=ConstantYield(90e6))
        c = trial.candidate(); dt = 1e-8
        gradient = np.diag([2e5, -1e5, -1e5])
        c['mean_velocity_m_s'] = c['position_m']@gradient.T
        result = trial.evaluate(c, dt)
        expected = np.broadcast_to(np.diag([60e6, -30e6, -30e6]), (27, 3, 3))
        np.testing.assert_allclose(result['material']['deviatoric_stress_Pa'], expected, rtol=2e-12, atol=1e-6)
        np.testing.assert_allclose(result['material']['ep'], .001, rtol=2e-12)
        expected_p = 90e6*.001*(.0005**3)
        expected_return = 1.5*30e9*.001**2*(.0005**3)
        np.testing.assert_allclose(result['material']['plastic_work_J'], expected_p, rtol=2e-12)
        np.testing.assert_allclose(result['material']['return_dissipation_J'], expected_return, rtol=2e-12)
        expected_t = 400.+.9*expected_p/(2700.*.0005**3*900.)
        np.testing.assert_allclose(result['map']['temperature_K'], expected_t, rtol=0, atol=1e-12)
        self.assertEqual(result['material_scope'], 'MANUFACTURED_TEST_ONLY')
        np.testing.assert_array_equal(pa.ep, 0.)

    def test_08_temperature_feedback_and_pair_conduction_layers_are_explicit(self):
        pa, p = cloud(); pa.temperature[:] = 390.+np.arange(27)
        trial = u.UpdatedTrial(pa, p, 0.); c = trial.candidate()
        c['temperature_K'] += 20.
        original_boundary = u.evaluate_boundary
        with patch.object(u, 'evaluate_boundary', wraps=original_boundary) as observer:
            result = trial.evaluate(c, 1e-9)
        np.testing.assert_array_equal(observer.call_args.args[3], pa.temperature+10.)
        np.testing.assert_array_equal(result['midpoint']['temperature_K'], pa.temperature+10.)
        # Independent scalar heat-pair accumulation from returned geometry.
        mid = result['midpoint']; power = np.zeros(27)
        for k, (i, j) in enumerate(zip(mid['i'], mid['j'])):
            r = mid['position_m'][i]-mid['position_m'][j]
            coefficient = 2*167.*(.0005**3)**2*(-float(r@mid['g'][k]))/float(r@r)
            q = coefficient*(pa.temperature[j]-pa.temperature[i])
            power[i] += q; power[j] -= q
        np.testing.assert_allclose(result['thermal']['pair_power_W'], power, rtol=2e-13, atol=1e-14)
        expected = pa.temperature+1e-9*power/(pa.m*900.)
        np.testing.assert_allclose(result['map']['temperature_K'], expected, rtol=0, atol=1e-12)

    def test_09_map_uses_event_average_not_endpoint_midpoint(self):
        pa, p = cloud(); pa.u[:] = 1.; trial = u.UpdatedTrial(pa, p, 0.)
        c = trial.candidate(); dt = 1e-9; real_boundary = u.evaluate_boundary
        def boundary_with_documented_displacement(*args, **kwargs):
            out = real_boundary(*args, **kwargs)
            # A manufactured interface probe only; physical reversal is tested
            # independently using the actual contact kick in integration suite.
            out['velocity'][:, 0] = -2./3
            out['work_velocity'][:, 0] = -1./18
            return out
        with patch.object(u, 'evaluate_boundary', side_effect=boundary_with_documented_displacement):
            result = trial.evaluate(c, dt)
        np.testing.assert_allclose(result['map']['position_m'][:, 0], pa.x-dt/18., rtol=0, atol=1e-19)
        np.testing.assert_allclose(result['map']['mean_velocity_m_s'][:, 0], -1./18, rtol=0, atol=0)
        self.assertGreater(np.max(np.abs(result['map']['position_m'][:, 0]-(pa.x+dt/6))), 1e-10)

    def test_10_thermal_failure_cannot_commit_trial_history(self):
        pa, p = cloud(); trial = u.UpdatedTrial(pa, p, 0.); before = trial.accepted_hash()
        with patch.object(u, 'advance_thermal', side_effect=ValueError('INJECTED_NO_COMMIT')):
            with self.assertRaisesRegex(ValueError, 'INJECTED_NO_COMMIT'):
                trial.evaluate(trial.candidate(), 1e-9)
        self.assertEqual(before, trial.accepted_hash()); self.assertEqual(before, u.particle_hash(pa))
        result = trial.evaluate(trial.candidate(), 1e-9)
        self.assertFalse(result['accepted'])

    def test_11_degenerate_3d_neighborhood_fails_without_pseudoinverse(self):
        pa, p = cloud()
        # Keep one distinct plane; collapsing all layers would instead invoke
        # the earlier, equally valid coincident-particle cusp rejection.
        pa = pa.extract_particles(np.flatnonzero(pa.z == .0006).astype(np.uint32))
        trial = u.UpdatedTrial(pa, p, 0.)
        with self.assertRaisesRegex(ValueError, 'degenerate/ill-conditioned'):
            trial.evaluate(trial.candidate(), 1e-9)

    def test_12_metadata_and_units_distinguish_residual_from_thermal_source(self):
        pa, p = cloud(); trial = u.UpdatedTrial(pa, p, 0.); c = trial.candidate()
        c['temperature_K'] += 5.
        result = trial.evaluate(c, 1e-9)
        np.testing.assert_array_equal(result['residual']['temperature_K'], 5.)
        self.assertEqual(result['thermal']['ledger']['total_heat_input_J'], 0.)
        self.assertEqual(result['thermal']['ledger']['mc_delta_T_residual_J'], 0.)
        self.assertFalse(result['material_energy_certified']); self.assertFalse(result['production_authorized'])


if __name__ == '__main__': unittest.main(verbosity=2)

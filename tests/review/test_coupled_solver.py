"""Short real-342-particle coupling regressions, not a deposition trajectory."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'scripts'))
from run_coupled_afsd import settings, initial_particles
from afsd_backend.pysph import coupled_solver as c
from afsd_backend.pysph.output import restore, checkpoint
from afsd_backend.pysph.coupled_boundary import BoundaryEvaluationError
from afsd_scaffold.vtk_legacy import read_point_cloud


def make():
    _, p, _ = settings()
    return c.IntegratedSolver(initial_particles(p), p,
        source_hashes={'test': 'frozen-unit-source'}, configuration_hash='frozen-unit-config')


class CoupledSolverTests(unittest.TestCase):
    def test_01_actual_deformable_components_initial_energy_momentum_and_open_tool(self):
        b = make(); p = b.pa
        self.assertEqual(len(p.x), 342)
        self.assertEqual(np.count_nonzero(p.component_id == 0), 192)
        self.assertEqual(np.count_nonzero(p.component_id == 1), 150)
        cap = p.z0 > .0025
        self.assertEqual(np.count_nonzero(cap), 25)
        np.testing.assert_array_equal(p.u, np.zeros(342))
        np.testing.assert_array_equal(p.w[cap], np.full(25, -.00212))
        np.testing.assert_array_equal(p.w[~cap], np.zeros(317))
        self.assertAlmostEqual(b.initial_kinetic_J, .5*25*3.375e-7*.00212**2, delta=1e-25)
        np.testing.assert_allclose(b.initial_linear_momentum, [0, 0, -25*3.375e-7*.00212], rtol=0, atol=1e-21)
        self.assertGreater(b.initial_normal_potential_J, 0.)
        self.assertEqual(np.linalg.matrix_rank(c.positions(p)-c.positions(p).mean(axis=0)), 3)
        self.assertFalse(b.state()['production_authorized'])

    def test_02_one_step_actual_three_dimensional_force_boundary_and_heat(self):
        b = make(); old_t = b.pa.temperature.copy(); b.advance(max_dt=1e-9)
        j = b.last_journal; traces = b.last_contact_trace
        self.assertEqual(len(traces), 342)
        self.assertEqual(j['stage'], 'ACCEPTED_COUPLED_DRIFT_KICK_DRIFT')
        self.assertGreater(np.linalg.norm(b.pa.S), 0.)
        # The first contact/material kick need not have propagated motion into
        # the initially quiescent substrate yet. Check actual free bar response
        # and cross-component neighbor connectivity, not an invented onset time.
        free_bar = (b.pa.component_id == 1)&(b.pa.z0 <= .0025)
        self.assertGreater(np.linalg.norm(c.velocities(b.pa)[free_bar]), 0.)
        self.assertTrue(np.any(b.pa.component_id[j['i']] != b.pa.component_id[j['j']]))
        self.assertGreater(np.max(np.abs(b.pa.temperature-old_t)), 0.)
        self.assertGreater(j['spatial']['pair_count'], 0)
        self.assertLess(j['spatial']['condition_max'], 1e6)
        self.assertGreater(b.ledger['friction_raw_J'], 0.)
        for trace in traces:
            np.testing.assert_array_equal(trace['free_axes'][:2], [1, 1])
            for interval in trace['result']['intervals']:
                np.testing.assert_allclose(trace['free_axes']*interval['constraint_impulse'], 0., rtol=0, atol=1e-17)

    def test_03_independent_thermal_bookkeeping_no_duplicate_friction(self):
        b = make(); t0 = b.pa.temperature.copy(); b.advance(max_dt=1e-9); j = b.last_journal
        mechanical = j['boundary']; th = j['thermal']; m = b.pa.m
        expected = (b.p['plastic_heat_fraction']*j['plastic_work_J']+
            b.p['heat_partition']*mechanical['friction_raw_J']+
            1e-9*mechanical['contact_power_W']+1e-9*th['pair_power_W'])
        measured = m*900*(b.pa.temperature-t0)
        np.testing.assert_allclose(measured, expected, rtol=0, atol=5e-17)
        self.assertAlmostEqual(float(np.sum(th['pair_power_W'])), 0., delta=1e-15)
        self.assertAlmostEqual(b.ledger['friction_heat_J'], .5*b.ledger['friction_raw_J'], delta=1e-22)
        self.assertAlmostEqual(b.ledger['wall_heat_out_J'], .5*b.ledger['friction_raw_J']-b.ledger['contact_heat_J'], delta=1e-22)
        self.assertAlmostEqual(b.ledger['thermal_balance_residual_J'], float(np.sum(measured-expected)), delta=1e-22)

    def test_04_labels_cached_actuator_and_as_rs_cannot_change_coupled_physics(self):
        a, b = make(), make()
        b.pa.tracer_type[:] = 99-b.pa.tracer_type; b.pa.wholebar_tag[:] = 17
        b.pa.as_rs_tag[:] = -9; b.pa.actuator[:] = 0
        a.advance(max_dt=1e-9); b.advance(max_dt=1e-9)
        for k in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'S', 'ep', 'temperature', 'plastic_work'):
            np.testing.assert_array_equal(getattr(a.pa, k), getattr(b.pa, k))
        self.assertEqual(a.ledger, b.ledger)

    def test_05_nonzero_stress_temperature_checkpoint_continuation(self):
        a = make(); a.advance(max_dt=1e-9)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; a.checkpoint(path)
            b = c.IntegratedSolver.from_checkpoint(path, a.p,
                source_hashes=a.source_hashes, configuration_hash=a.configuration_hash)
            a.advance(max_dt=1e-9); b.advance(max_dt=1e-9)
        for k in a.pa.properties: np.testing.assert_array_equal(getattr(a.pa, k), getattr(b.pa, k))
        self.assertEqual(a.ledger, b.ledger); self.assertEqual(a.time, b.time)

    def test_06_reorder_retains_id_material_temperature_and_ledger(self):
        a = make(); a.advance(max_dt=1e-9)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; a.checkpoint(path); pa, state = restore(path)
            order = np.random.default_rng(86).permutation(342).astype(np.int32)
            b = c.IntegratedSolver(pa.extract_particles(order), a.p, state,
                source_hashes=a.source_hashes, configuration_hash=a.configuration_hash)
            a.advance(max_dt=1e-9); b.advance(max_dt=1e-9)
        sort = np.argsort(b.pa.particle_id)
        for k in ('particle_id', 'tracer_type', 'component_id', 'm'):
            np.testing.assert_array_equal(getattr(a.pa, k), getattr(b.pa, k)[sort])
        for k in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'temperature', 'ep'):
            np.testing.assert_allclose(getattr(a.pa, k), getattr(b.pa, k)[sort], rtol=0, atol=1e-12)
        np.testing.assert_allclose(a.pa.S.reshape(-1, 9), b.pa.S.reshape(-1, 9)[sort], rtol=0, atol=1e-6)
        for k in a.ledger: self.assertAlmostEqual(a.ledger[k], b.ledger[k], delta=1e-12)

    def test_07_thermal_failure_restores_material_boundary_ledger_and_counters(self):
        b = make(); b.advance(max_dt=1e-9)
        arrays = {k: getattr(b.pa, k).copy() for k in b.pa.properties}
        old_state = copy.deepcopy(b.state())
        with patch.object(c, 'advance_thermal', side_effect=ValueError('INJECTED_THERMAL_REJECTION')):
            with self.assertRaisesRegex(ValueError, 'INJECTED_THERMAL_REJECTION'): b.advance(max_dt=1e-9)
        for k, a in arrays.items(): np.testing.assert_array_equal(getattr(b.pa, k), a)
        self.assertEqual(b.ledger, old_state['ledger']); self.assertEqual(b.time, old_state['t'])
        self.assertEqual(b.count, old_state['count']); self.assertEqual(b.floor_evaluations, old_state['floor_evaluations'])
        self.assertEqual(b.rejected_candidate['classification'], 'REJECTED')
        self.assertFalse(b.rejected_candidate['restart_allowed'])
        self.assertFalse(np.array_equal(b.rejected_candidate['properties']['S'], arrays['S']))
        self.assertEqual(b.rejected_candidate['journal']['stage'], 'CANDIDATE_MATERIAL_BOUNDARY_BEFORE_THERMAL')

    def test_08_original_boundary_error_context_preserved(self):
        b = make(); failure = BoundaryEvaluationError('INJECTED_PARTICLE_CONTACT')
        failure.contact_context = {'particle_id': 7, 'feature_ids': [2, 5], 'failed_interval': {'raw': -1.}}
        with patch.object(c, 'evaluate_boundary', side_effect=failure):
            with self.assertRaisesRegex(BoundaryEvaluationError, 'INJECTED_PARTICLE_CONTACT'): b.advance(max_dt=1e-9)
        self.assertEqual(b.rejected_candidate['boundary_context'], failure.contact_context)
        self.assertEqual(b.count, 0)

    def test_09_nonzero_feature_mask_loss_rejected_but_normal_opening_allowed(self):
        b = make(); x = c.positions(b.pa); pid = int(b.pa.particle_id[0]); cid = f'{pid}:0'
        x0 = x.copy(); x1 = x.copy()
        x0[0, 2] = -b.p['substrate_depth']+.1*b.p['dx']
        x1[0, 2] = x0[0, 2]
        before = {'potential_by_contact_id_J': {cid: 1e-6}}
        after = {'potential_by_contact_id_J': {}}
        with self.assertRaisesRegex(RuntimeError, 'UNSUPPORTED_STORED_ENERGY_FEATURE_TRANSITION'):
            b._feature_transition_guard(before, after, x0, x1, 0., 1e-9)
        x1[0, 2] = -b.p['substrate_depth']+b.p['dx']
        result = b._feature_transition_guard(before, after, x0, x1, 0., 1e-9)
        self.assertEqual(result['removed_contact_ids'], [cid])

    def test_10_checkpoint_rejects_parameter_source_and_rejected_state(self):
        a = make()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; a.checkpoint(path)
            p = dict(a.p, feed=.00106)
            with self.assertRaisesRegex(ValueError, 'parameters mismatch'):
                c.IntegratedSolver.from_checkpoint(path, p, a.source_hashes, a.configuration_hash)
            with self.assertRaisesRegex(ValueError, 'source mismatch'):
                c.IntegratedSolver.from_checkpoint(path, a.p, {'test': 'different'}, a.configuration_hash)
            pa, state = restore(path); state['classification'] = 'REJECTED'; state['restart_allowed'] = False
            checkpoint(Path(td)/'rejected.npz', pa, state)
            with self.assertRaises(ValueError): c.IntegratedSolver.from_checkpoint(Path(td)/'rejected.npz', a.p)

    def test_11_actual_vtk_has_solved_temperature_stress_density_not_fake_energy(self):
        b = make(); b.advance(max_dt=1e-9)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'state.vtk'; b.vtk(path); v = read_point_cloud(path)
        np.testing.assert_array_equal(v['scalars']['temperature_K'], b.pa.temperature)
        np.testing.assert_array_equal(v['scalars']['density_kg_per_m3'], b.pa.rho)
        self.assertIn('stress_yz_Pa', v['scalars'])
        self.assertNotIn('internal_energy_J_per_kg', v['scalars'])
        r = b.record(); np.testing.assert_array_equal(r['as_rs_tag'], b.pa.as_rs_tag)

    def test_12_stable_dt_and_initial_energy_subtraction(self):
        b = make(); p = b.p
        wave = .1*.00065/(np.sqrt((b.bulk+4*26.3e9/3)/2700)+.00212)
        normal = .1*np.sqrt(3.375e-7/b.kn)
        self.assertLessEqual(b.stable_dt(), min(wave, normal)+1e-20)
        self.assertEqual(b.diagnostics(False)['energy_residual_J'], 0.)
        self.assertGreater(b.initial_kinetic_J, 0.)
        with self.assertRaises(ValueError): b.advance(max_dt=0.)
        b.pa.temperature[0] = np.nan
        with self.assertRaisesRegex(ValueError, 'nonfinite current'): b.diagnostics(fail=True)


if __name__ == '__main__': unittest.main(verbosity=2)

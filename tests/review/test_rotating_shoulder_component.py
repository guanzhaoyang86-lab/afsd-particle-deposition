"""Short (at most two-step) component checks, independent literal references.

These tests are NOT the full 40/80-step component runs or continuum tests.
"""
from pathlib import Path
import sys
import tempfile
import unittest
import copy
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import probe_rotating_shoulder_component as component
from afsd_scaffold.vtk_legacy import read_point_cloud


class RotatingShoulderShortChecks(unittest.TestCase):
    def test_01_literal_first_step_and_normal_reaction(self):
        block = component.RotatingShoulderComponent(2e-5)
        record = block.advance()
        expected_v = np.array([[0., .004, 0.], [-.004, 0., 0.],
                               [0., -.004, 0.], [.004, 0., 0.]])
        expected_x = np.array([[.01, 4e-8, .0013], [-4e-8, .012, .0013],
                               [-.01, -4e-8, .0013], [4e-8, -.012, .0013]])
        np.testing.assert_allclose(record['velocity_m_s'], expected_v, rtol=0., atol=1e-14)
        np.testing.assert_allclose(record['position_m'], expected_x, rtol=0., atol=1e-16)
        self.assertAlmostEqual(record['kinetic_J'], 3.2e-8, delta=1e-19)
        self.assertAlmostEqual(record['cumulative']['friction_raw_J'], 3.488e-6, delta=1e-18)
        self.assertAlmostEqual(record['cumulative']['wall_work_J'], 3.52e-6, delta=1e-18)
        for trace in block.last_traces:
            result = trace['result']
            np.testing.assert_allclose(result['normal_impulse'], [[0., 0., -2e-5]], atol=1e-19)
            np.testing.assert_allclose(result['constraint_impulse'], [0., 0., 2e-5], atol=1e-19)
            self.assertAlmostEqual(result['constraint_work_J'], 0., delta=1e-18)
            self.assertEqual(sum(a['normal_work_J'].sum() for a in result['intervals']), 0.)
            self.assertEqual(trace['feature_id'], 1)
        self.assertLess(abs(record['global_energy_residual_J']), 1e-18)

    def test_02_passive_labels_do_not_change_two_step_dynamics(self):
        a, b = [component.RotatingShoulderComponent(2e-5) for _ in range(2)]
        b.pa.tracer_type[:] = [99, 0, 81, 0]; b.pa.as_rs_tag[:] = [-1, 1, 17, -12]
        for _ in range(2): a.advance(); b.advance()
        np.testing.assert_array_equal(component.positions(a.pa), component.positions(b.pa))
        np.testing.assert_array_equal(component.velocities(a.pa), component.velocities(b.pa))
        self.assertEqual(a.cumulative, b.cumulative)

    def test_03_permutation_uses_persistent_ids(self):
        a = component.RotatingShoulderComponent(2e-5)
        pa = a.pa.extract_particles(np.array([3, 1, 0, 2], dtype=np.int32))
        b = component.RotatingShoulderComponent(2e-5, pa=pa)
        for _ in range(2): a.advance(); b.advance()
        order = np.argsort(b.pa.particle_id)
        np.testing.assert_array_equal(a.pa.particle_id, b.pa.particle_id[order])
        np.testing.assert_array_equal(component.positions(a.pa), component.positions(b.pa)[order])
        np.testing.assert_array_equal(component.velocities(a.pa), component.velocities(b.pa)[order])
        np.testing.assert_array_equal(a.pa.tracer_type, b.pa.tracer_type[order])

    def test_04_real_pysph_checkpoint_restart_two_step_identity(self):
        a, b = [component.RotatingShoulderComponent(1e-5) for _ in range(2)]
        a.advance(); a.advance(); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; b.save_checkpoint(path)
            c = component.RotatingShoulderComponent.from_checkpoint(path); c.advance()
        self.assertEqual(a.step, c.step); self.assertEqual(a.time, c.time)
        self.assertEqual(a.cumulative, c.cumulative)
        for key in ('particle_id', 'tracer_type', 'as_rs_tag', 'x', 'y', 'z', 'u', 'v', 'w'):
            np.testing.assert_array_equal(getattr(a.pa, key), getattr(c.pa, key))

    def test_05_callback_failure_rolls_back_and_marks_rejected_candidate(self):
        b = component.RotatingShoulderComponent(2e-5); b.advance()
        before_x = component.positions(b.pa).copy(); before_v = component.velocities(b.pa).copy()
        before_cumulative = dict(b.cumulative)
        with patch.object(component, 'kick', side_effect=ValueError('INJECTED_CONTACT_FAILURE')):
            with self.assertRaisesRegex(ValueError, 'INJECTED_CONTACT_FAILURE'): b.advance()
        np.testing.assert_array_equal(component.positions(b.pa), before_x)
        np.testing.assert_array_equal(component.velocities(b.pa), before_v)
        self.assertEqual(b.step, 1); self.assertEqual(b.cumulative, before_cumulative)
        self.assertEqual(b.rejected_candidate['classification'], 'REJECTED')
        self.assertFalse(b.rejected_candidate['restart_allowed'])
        self.assertFalse(np.array_equal(b.rejected_candidate['position_m'], before_x))

    def test_06_vtk_omits_unsolved_material_and_thermal_fields(self):
        b = component.RotatingShoulderComponent(2e-5); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'frame.vtk'; b.vtk(path); data = read_point_cloud(path)
        self.assertEqual(set(data['scalars']), {'particle_id', 'tracer_type', 'as_rs_tag', 'mass_kg'})
        self.assertEqual(set(data['vectors']), {'velocity_m_per_s'})
        self.assertIn('STRESS_THERMAL_NOT_SOLVED', data['title'])
        np.testing.assert_array_equal(data['scalars']['particle_id'], b.pa.particle_id)
        np.testing.assert_array_equal(data['points'], component.positions(b.pa))

    def test_07_rk4_reference_independent_of_contact_code(self):
        with patch.object(component, 'kick', side_effect=AssertionError('reference called kick')):
            with patch.object(component, 'contact_rows', side_effect=AssertionError('reference called geometry')):
                r = component.reference_rk4(1e-7, final_time=2e-5)
        self.assertAlmostEqual(r['velocity_m_s'][0, 1], .004, delta=1e-11)
        self.assertAlmostEqual(r['position_m'][0, 1], 4e-8, delta=1e-16)
        self.assertLess(r['velocity_m_s'][0, 0], 0.)
        self.assertGreater(r['minimum_slip_m_s'], .19)
        np.testing.assert_allclose(r['kinetic_J']+r['friction_raw_J']-r['wall_work_J'], 0., atol=1e-18)

    def test_08_no_hidden_particle_scaling_or_extra_force(self):
        b = component.RotatingShoulderComponent(2e-5)
        np.testing.assert_array_equal(b.pa.m, [.001]*4)
        b.advance()
        self.assertTrue(all(a['result']['other_work_J'] == 0. for a in b.last_traces))
        self.assertEqual(b.config['physical_scaling'], 'NONE')
        self.assertEqual(b.config['field_status']['stress'], 'NOT_SOLVED')
        with self.assertRaises(ValueError): component.RotatingShoulderComponent(5e-6)

    def test_09_reference_rejects_configuration_drift(self):
        changed = copy.deepcopy(component.INPUT_CONFIG)
        changed['parameters']['mass']['value'] = .002
        with patch.object(component, 'INPUT_CONFIG', changed):
            with self.assertRaisesRegex(ValueError, 'reference configuration mismatch: mass'):
                component.reference_rk4(1e-7, final_time=2e-5)
        with patch.dict(component.PARAMETERS, {'omega': 40.}):
            with self.assertRaisesRegex(ValueError, 'reference configuration mismatch: omega'):
                component.reference_rk4(1e-7, final_time=2e-5)

    def test_10_reference_literal_provenance_and_hashes(self):
        r = component.reference_rk4(1e-7, final_time=2e-5)
        p = r['reference_parameters']
        self.assertEqual(p['mass_kg'], .001); self.assertEqual(p['friction_force_cap_N'], .2)
        self.assertEqual(p['omega_rad_s'], [0., 0., 20.])
        self.assertEqual(p['requested_end_s'], 2e-5)
        self.assertEqual(r['source_sha256'], component.hashes())
        self.assertIn('docs/ROTATING_CONTACT_CONTRACT.md', r['source_sha256'])
        self.assertEqual(r['configuration_sha256'], component.config_hash(r['configuration']))

    def test_11_checkpoint_restart_is_not_production_authorization(self):
        b = component.RotatingShoulderComponent(1e-5); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; b.save_checkpoint(path)
            pa, state = component.restore(path)
            self.assertFalse(state['production_authorized'])
            self.assertTrue(state['restart_allowed']); self.assertTrue(state['accepted'])
            fields = state['configuration']['field_status']
            self.assertIn('NOT_SOLVED', fields['checkpoint_rho'])
            self.assertIn('NOT_SOLVED', fields['checkpoint_e'])
            self.assertIn('NOT_SPH_RESOLUTION', fields['checkpoint_h'])
            state['restart_allowed'] = False; state['classification'] = 'REJECTED'
            component.checkpoint(Path(td)/'rejected.npz', pa, state)
            with self.assertRaisesRegex(ValueError, 'REJECTED candidate'):
                component.RotatingShoulderComponent.from_checkpoint(Path(td)/'rejected.npz')


if __name__ == '__main__': unittest.main(verbosity=2)

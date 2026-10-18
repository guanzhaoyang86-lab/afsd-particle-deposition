"""Real PySPH factory/state tests with no nonlinear solve or physical acceptance.

The nonzero prefix fixture below is manufactured bookkeeping evidence only;
it is never reported, checkpointed or reused as a solved deposition state.
"""
from pathlib import Path
import copy
import json
import math
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
sys.path.insert(0, str(ROOT/'src'))
import afsd_process_framework_20261007 as model
from afsd_backend.pysph import discrete_process
from afsd_backend.pysph.interval_acceptance import (
    initial_interval_energy_state, require_interval_energy)
from afsd_backend.pysph.output import checkpoint, restore
from afsd_backend.pysph.particles import IDENTITY
from pysph.base.particle_array import ParticleArray
from pysph.solver.utils import load


class ProcessFrameworkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original_pa, cls.original_state = model.independent_initial_state()

    def setUp(self):
        self.pa = self.original_pa.extract_particles(
            np.arange(len(self.original_pa.x), dtype=np.uint32))
        self.state = copy.deepcopy(self.original_state)
        self.p = self.state['parameters']

    def assert_state_arrays_unchanged(self, before):
        for key, value in before.items():
            with self.subTest(property=key):
                np.testing.assert_array_equal(self.pa.properties[key].get_npy_array(), value)

    def test_real_pysph_factory_has_full_1350_particle_geometry_and_mass(self):
        self.assertIsInstance(self.pa, ParticleArray)
        self.assertEqual(len(self.pa.x), 1350)
        substrate = self.pa.component_id == 0
        bar = self.pa.component_id == 1
        self.assertEqual(int(substrate.sum()), 1200)
        self.assertEqual(int(bar.sum()), 150)
        self.assertEqual(len(set(self.pa.particle_id)), 1350)
        np.testing.assert_array_equal(self.pa.particle_id, np.arange(1, 1351))
        self.assertAlmostEqual(float(self.pa.m.sum()), 2700*1350*(.0005**3), places=17)
        self.assertAlmostEqual(np.ptp(self.pa.x[substrate])+.0005, .010, places=15)
        self.assertAlmostEqual(np.ptp(self.pa.y[substrate])+.0005, .010, places=15)
        self.assertAlmostEqual(np.ptp(self.pa.z[substrate])+.0005, .0015, places=15)
        self.assertAlmostEqual(np.ptp(self.pa.z[bar])+.0005, .003, places=15)
        self.assertEqual(self.state['backend']['name'], 'PySPH')
        self.assertEqual(self.state['backend']['commit'],
                         '1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5')

    def test_free_mask_constrains_only_bottom_xyz_and_material_top_cap_z(self):
        bottom = (self.pa.component_id == 0) & (self.pa.z0 < -.001)
        cap = (self.pa.component_id == 1) & (self.pa.z0 > .0025)
        np.testing.assert_array_equal(self.pa.actuator.astype(bool), cap)
        expected = np.ones((1350, 3), dtype=bool)
        expected[bottom] = False; expected[cap, 2] = False
        np.testing.assert_array_equal(self.state['free'], expected)
        self.assertEqual(int(bottom.sum()), 400)
        self.assertEqual(int(cap.sum()), 25)
        self.assertEqual(int(self.state['free'].sum()), 2825)
        self.assertTrue(np.all(self.state['free'][cap, :2]))
        self.assertTrue(np.all(self.state['free'][(self.pa.component_id == 1) & ~cap]))

    def test_initial_material_is_at_rest_and_not_given_tool_rotation(self):
        before = {k: v.get_npy_array().copy() for k, v in self.pa.properties.items()}
        trial = model.ProcessTrial(self.pa, self.state, .001)
        np.testing.assert_array_equal(np.column_stack((self.pa.u, self.pa.v, self.pa.w)),
                                      np.zeros((1350, 3)))
        np.testing.assert_array_equal(trial.prescribed(.001), np.zeros((1350, 3)))
        self.assertGreater(trial.motion.end[0], 0.)
        self.assertFalse(self.state['process_program']['free_material_trajectory_prescribed'])
        self.assert_state_arrays_unchanged(before)

    def test_factory_preserves_material_ids_energy_and_all_contact_history_fields(self):
        expected_pa, expected_state = discrete_process.independent_initial_state(self.p)
        self.assertEqual(set(self.pa.properties), set(expected_pa.properties))
        for key in self.pa.properties:
            with self.subTest(property=key):
                np.testing.assert_array_equal(self.pa.properties[key].get_npy_array(),
                                              expected_pa.properties[key].get_npy_array())
        for key in IDENTITY:
            self.assertIn(key, self.pa.properties)
        self.assertEqual(self.state['interval_energy'], expected_state['interval_energy'])
        self.assertEqual(self.state['bore_observation'], expected_state['bore_observation'])
        np.testing.assert_array_equal(self.state['contact_history'], expected_state['contact_history'])
        self.assertEqual(self.state['contact_history'].shape, (1350, 8, 3))
        self.assertTrue(np.all(self.pa.temperature == 723.15))
        self.assertTrue(np.all(self.pa.S == 0.))
        self.assertTrue(np.all(self.pa.ep == 0.))
        self.assertEqual(self.state['crossing_ids'], [])
        self.assertEqual(self.state['crossing_mass_kg'], 0.)
        self.assertEqual(len(self.state['bore_observation']['eligible_ids']), 75)
        self.assertEqual(int(((self.pa.component_id == 1) & (self.pa.z < .0015)).sum()), 75)

    def test_factory_has_independent_unaccepted_provenance(self):
        self.assertEqual(self.state['classification'], model.CLASSIFICATION)
        self.assertEqual(self.state['time_s'], 0.)
        self.assertEqual(self.state['step'], 0)
        self.assertTrue(self.state['independent_initialization'])
        self.assertIsNone(self.state['parent_checkpoint'])
        self.assertFalse(self.state['accepted'])
        self.assertFalse(self.state['restart_allowed'])
        self.assertTrue(self.state['diagnostic_only'])
        self.assertEqual(self.state['physical_qualification'], 'NOT_RUN')
        self.assertEqual(self.state['accepted_time_advanced_s'], 0.)
        self.assertEqual(self.state['initial_particle_ids'], self.pa.particle_id.tolist())
        self.assertEqual(self.state['initial_mass_kg'], float(self.pa.m.sum()))
        self.assertEqual(self.state['machine_history'], [])

    def test_default_configuration_has_full_160ms_finite_domain_inventory_plan(self):
        config, p, program, audit = model.load_configuration()
        self.assertEqual(p['substrate_width'], .010)
        self.assertEqual(p['bore_width'], .002504)
        self.assertEqual(p['omega'], 10*math.pi)
        self.assertAlmostEqual(program.end_s, .160, places=15)
        self.assertTrue(audit['passed'])
        self.assertEqual(audit['failures'], [])
        self.assertEqual(audit['particle_count'], 1350)
        self.assertAlmostEqual(audit['programmed_feed_m'], .0002756, places=15)
        self.assertAlmostEqual(audit['programmed_traverse_m'], .0002332, places=15)
        self.assertAlmostEqual(audit['programmed_revolutions'], .75, places=14)
        self.assertGreater(audit['final_forward_edge_margin_m'], audit['required_edge_buffer_m'])
        self.assertLess(audit['programmed_feed_m'], audit['usable_axial_command_m'])
        self.assertEqual(config['bar_rotation'], 'FREE_XY_CONTACT_DRIVEN_BY_ROTATING_SQUARE_BORE')
        self.assertIn('NOT_COLD_START', config['thermal_initialization'])

    def test_frozen_physics_operators_are_inherited_without_replacement(self):
        self.assertIs(model.ProcessTrial.evaluate, model.old.RampedHotTrial.evaluate)
        self.assertIs(model.ProcessTrial.tangent, model.old.RampedHotTrial.tangent)
        self.assertIs(model.ProcessTrial.prescribed, model.old.RampedHotTrial.prescribed)
        audit = model.operator_reuse_audit()
        self.assertEqual(audit['status'], 'PASS_UNCHANGED_INHERITED_THERMOMECHANICAL_OPERATORS')
        self.assertFalse(audit['runtime_monkeypatch'])
        self.assertFalse(audit['material_trajectory_prescribed'])

    def test_old_state_or_classification_only_relabel_cannot_be_continued(self):
        _, old_state = discrete_process.independent_initial_state(self.p)
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, old_state, .001)
        old_state['classification'] = model.CLASSIFICATION
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, old_state, .001)
        for key, value in (('classification', discrete_process.CLASSIFICATION),
                           ('accepted', True), ('restart_allowed', True),
                           ('independent_initialization', False),
                           ('parent_checkpoint', 'old/accepted.npz')):
            bad = copy.deepcopy(self.state); bad[key] = value
            with self.subTest(field=key), self.assertRaises(ValueError):
                model.validate_state(bad)

    def test_changed_parameters_configuration_or_program_rejected(self):
        mutations = [lambda s: s['parameters'].update(friction=.9),
                     lambda s: s['framework_config']['phases'][0].update(duration_s=.03),
                     lambda s: s['process_program'].update(shoulder_height_m=.002),
                     lambda s: s['domain_inventory_plan'].update(programmed_feed_m=.0001)]
        for index, change in enumerate(mutations):
            bad = copy.deepcopy(self.state); change(bad)
            with self.subTest(index=index), self.assertRaises(ValueError):
                model.validate_state(bad)

    def test_added_constraints_and_wrong_interval_rejected(self):
        cap_index = np.flatnonzero(self.pa.actuator)[0]
        self.state['free'][cap_index, 0] = False
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, self.state, .001)
        self.state['free'][cap_index, 0] = True
        trial = model.ProcessTrial(self.pa, self.state, .001)
        with self.assertRaises(ValueError):
            trial.prescribed(.002)

    def test_persistent_id_loss_or_total_mass_change_rejected(self):
        self.pa.particle_id[-1] = 999999
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, self.state, .001)
        self.pa.particle_id[-1] = 1350
        self.pa.m[-1] *= 1.1
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, self.state, .001)

    def test_redistributed_particle_mass_with_unchanged_total_rejected(self):
        before = float(self.pa.m.sum())
        delta = float(self.pa.m[-1])*.1
        self.pa.m[-1] += delta; self.pa.m[-2] -= delta
        self.assertAlmostEqual(float(self.pa.m.sum()), before, places=17)
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, self.state, .001)

    def test_component_and_actuator_identity_corruption_rejected(self):
        # A bar particle away from the actuator/bottom must remain bar material.
        free_bar = np.flatnonzero((self.pa.component_id == 1) & ~self.pa.actuator.astype(bool))[0]
        self.pa.component_id[free_bar] = 0
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, self.state, .001)
        self.pa.component_id[free_bar] = 1
        # Updating the free mask along with a bad actuator flag must not turn a
        # previously free material point into an artificial kinematic driver.
        self.pa.actuator[free_bar] = 1
        self.state['free'][free_bar, 2] = False
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, self.state, .001)

    def manufactured_prefix(self):
        """Bookkeeping only: does not call or impersonate a mechanics solve."""
        state = copy.deepcopy(self.state)
        program = model.validate_state(state)
        interval = program.step(0., .001)
        state['step'] = 1; state['time_s'] = interval.t1
        state['machine_history'] = [interval.to_dict()]
        state['bore_observation']['last_time_s'] = interval.t1
        state['interval_energy'] = require_interval_energy(state['interval_energy'],
            dict(global_raw_residual_J=1e-10, absolute_work_scale_J=1e-8))
        state['contact_history'][1200, 0, 0] = 1e-8
        state['physical_qualification'] = 'MANUFACTURED_BOOKKEEPING_FIXTURE_NOT_SOLVED'
        return state

    def test_existing_energy_and_contact_prefix_is_not_reset_by_trial(self):
        state = self.manufactured_prefix()
        original_energy = copy.deepcopy(state['interval_energy'])
        original_motion = copy.deepcopy(state['machine_history'])
        original_contact = state['contact_history'].copy()
        model.ProcessTrial(self.pa, state, .001)
        self.assertEqual(state['interval_energy'], original_energy)
        self.assertEqual(state['machine_history'], original_motion)
        np.testing.assert_array_equal(state['contact_history'], original_contact)
        state['interval_energy'] = initial_interval_energy_state()
        with self.assertRaises(ValueError):
            model.validate_state(state)

    def test_broken_machine_or_observer_prefix_is_rejected(self):
        for field in ('machine', 'observer', 'time'):
            state = self.manufactured_prefix()
            if field == 'machine':
                state['machine_history'][0]['feed_displacement_m'] = -.1
            elif field == 'observer':
                state['bore_observation']['last_time_s'] = 0.
            else:
                state['time_s'] = .002
            with self.subTest(field=field), self.assertRaises(ValueError):
                model.validate_state(state)

    def test_checkpoint_keeps_fields_and_refuses_accepted_restart(self):
        with tempfile.TemporaryDirectory(prefix='framework_unit_', dir=ROOT/'tmp') as temp:
            path = Path(temp)/'independent_unaccepted.npz'
            checkpoint(path, self.pa, self.state)
            saved = load(str(path))
            other = saved['arrays']['solid']; state = saved['solver_data']
            for key in self.pa.properties:
                with self.subTest(property=key):
                    np.testing.assert_array_equal(other.properties[key].get_npy_array(),
                                                  self.pa.properties[key].get_npy_array())
            self.assertEqual(state['interval_energy'], self.state['interval_energy'])
            self.assertEqual(state['bore_observation'], self.state['bore_observation'])
            self.assertEqual(state['machine_history'], self.state['machine_history'])
            np.testing.assert_array_equal(state['contact_history'], self.state['contact_history'])
            np.testing.assert_array_equal(state['free'], self.state['free'])
            model.validate_state(state)
            with self.assertRaises(ValueError):
                restore(path)


if __name__ == '__main__':
    unittest.main()

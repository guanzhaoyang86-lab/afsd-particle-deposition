"""Independent rectangular geometry and history tests, without physics solves."""
from pathlib import Path
import copy
import math
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
sys.path.insert(0, str(ROOT/'src'))
import afsd_deposition_setup_20261007 as model
from afsd_backend.pysph.interval_acceptance import initial_interval_energy_state, require_interval_energy
from afsd_backend.pysph.particles import make_array, IDENTITY
from pysph.base.particle_array import ParticleArray


class RectangularDepositionSetupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.initial_pa, cls.initial_state = model.independent_initial_state()

    def setUp(self):
        self.pa = self.initial_pa.extract_particles(np.arange(len(self.initial_pa.x), dtype=np.uint32))
        self.state = copy.deepcopy(self.initial_state)
        self.config = self.state['framework_config']
        self.p = self.state['parameters']

    def test_true_rectangle_cell_bounds_and_nonoverlapping_world_coordinates(self):
        self.assertIsInstance(self.pa, ParticleArray)
        self.assertEqual(len(self.pa.x), 1381)
        substrate = self.pa.component_id == 0
        bar = self.pa.component_id == 1
        self.assertEqual(int(substrate.sum()), 1056)
        self.assertEqual(int(bar.sum()), 325)
        points = np.column_stack((self.pa.x, self.pa.y, self.pa.z))
        self.assertEqual(len(np.unique(points, axis=0)), 1381)
        np.testing.assert_array_equal(points, np.column_stack((self.pa.x0, self.pa.y0, self.pa.z0)))
        bounds = np.column_stack((points[substrate].min(axis=0)-.00025,
                                  points[substrate].max(axis=0)+.00025))
        np.testing.assert_allclose(bounds, [[-.004, .007], [-.004, .004], [-.0015, 0.]], atol=1e-15, rtol=0.)
        self.assertAlmostEqual(float(points[substrate, 0].mean()), .0015, places=15)
        self.assertAlmostEqual(float(points[bar, 0].mean()), 0., places=15)
        self.assertAlmostEqual(float(points[bar, 1].mean()), 0., places=15)
        self.assertAlmostEqual(float(points[bar, 2].min()), .00025, places=15)
        self.assertAlmostEqual(float(points[bar, 2].max()), .00625, places=15)
        self.assertEqual(self.p['substrate_width'], .008)
        self.assertNotEqual(self.p['substrate_width'], bounds[0, 1]-bounds[0, 0])

    def test_true_mass_volume_ids_and_upstream_particle_fields(self):
        np.testing.assert_array_equal(self.pa.particle_id, np.arange(1, 1382))
        self.assertTrue(np.all(self.pa.m == 2700*.0005**3))
        self.assertAlmostEqual(float(self.pa.m.sum()), 2700*(.011*.008*.0015+.0025**2*.0065), places=17)
        self.assertAlmostEqual(float(self.pa.m[self.pa.component_id == 1].sum()/2700), 40.625e-9, places=19)
        self.assertEqual(self.state['initial_particle_ids'], self.pa.particle_id.tolist())
        for key in IDENTITY:
            self.assertIn(key, self.pa.properties)
        expected = make_array(np.column_stack((self.pa.x, self.pa.y, self.pa.z)), self.pa.component_id, self.p)
        for key in IDENTITY:
            with self.subTest(field=key):
                np.testing.assert_array_equal(getattr(self.pa, key), getattr(expected, key))
        for key in ('plastic_work', 'return_dissipation', 'S', 'D', 'ep', 'temperature', 'rho', 'h'):
            np.testing.assert_array_equal(getattr(self.pa, key), getattr(expected, key))

    def test_only_original_cap_z_and_substrate_bottom_xyz_are_constrained(self):
        cap = (self.pa.component_id == 1) & (self.pa.z0 > .006)
        bottom = (self.pa.component_id == 0) & (self.pa.z0 < -.001)
        self.assertEqual(int(cap.sum()), 25)
        self.assertEqual(int(bottom.sum()), 352)
        np.testing.assert_array_equal(self.pa.actuator.astype(bool), cap)
        expected = np.ones((1381, 3), bool)
        expected[bottom] = False; expected[cap, 2] = False
        np.testing.assert_array_equal(self.state['free'], expected)
        self.assertEqual(int(expected.sum()), 3062)
        self.assertTrue(np.all(expected[cap, :2]))
        self.assertTrue(np.all(expected[(self.pa.component_id == 1) & ~cap]))
        above_tool = (self.pa.component_id == 1) & (self.pa.z0 > self.p['tool_top'])
        self.assertEqual(int(above_tool.sum()), 200)
        self.assertTrue(np.all(expected[above_tool, :2]))

    def test_zero_initial_material_velocity_and_simultaneous_command_entry(self):
        np.testing.assert_array_equal(np.column_stack((self.pa.u, self.pa.v, self.pa.w)), np.zeros((1381, 3)))
        trial = model.ProcessTrial(self.pa, self.state, .001)
        u = trial.prescribed(.001)
        cap = self.pa.actuator.astype(bool)
        self.assertTrue(np.all(u[cap, 2] < 0.))
        self.assertTrue(np.all(u[~cap] == 0.))
        self.assertTrue(np.all(u[:, :2] == 0.))
        self.assertGreater(trial.motion.end[0], 0.)
        self.assertEqual(trial.motion.end[2], 0.)
        self.assertFalse(self.state['process_program']['free_material_trajectory_prescribed'])

    def test_prebonded_initial_lower_material_is_excluded_from_new_supply(self):
        bar = self.pa.component_id == 1
        below = bar & (self.pa.z0 < .0015)
        self.assertEqual(int(below.sum()), 75)
        self.assertEqual(len(self.state['bore_observation']['eligible_ids']), 250)
        self.assertFalse(set(self.pa.particle_id[below]) & set(self.state['bore_observation']['eligible_ids']))
        self.assertEqual(self.state['crossing_ids'], [])
        self.assertEqual(self.state['crossing_mass_kg'], 0.)
        self.assertTrue(np.all(self.pa.temperature == 723.15))
        self.assertTrue(np.all(self.pa.ep == 0.))
        self.assertEqual(self.state['contact_history'].shape, (1381, 8, 3))
        self.assertTrue(np.all(self.state['contact_history'] == 0.))

    def test_integrated_stationary_fill_and_track_are_explicit_command_budgets(self):
        _, _, program, audit = model.load_configuration()
        fill = program.integral(.220)
        self.assertAlmostEqual(fill[1], .0004452, places=15)
        self.assertEqual(fill[2], 0.)
        for t in (.001, .020, .100, .219):
            self.assertEqual(program.velocity(t)[2], 0.)
        self.assertAlmostEqual(audit['nominal_stationary_swept_volume_m3'], 2.7825e-9, places=19)
        self.assertAlmostEqual(audit['programmed_duration_s'], 1.6450943396226414, places=14)
        self.assertAlmostEqual(audit['programmed_traverse_m'], .003, places=15)
        self.assertAlmostEqual(audit['programmed_feed_m'], .0034664, places=15)
        self.assertAlmostEqual(audit['nominal_command_swept_volume_m3'], 21.665e-9, places=19)
        self.assertAlmostEqual(audit['programmed_revolutions'], 5*(program.end_s-.010), places=13)
        self.assertEqual(audit['stationary_fill_status'], 'PLANNED_COMMAND_ONLY_NOT_PROVEN_FILLED')
        self.assertEqual(audit['deposition_status'], 'NOT_PREDICTED_BY_GEOMETRY_OR_COMMAND_INTEGRALS')

    def test_full_pass_has_kernel_buffer_and_finite_cap_reserve(self):
        audit = self.state['domain_inventory_plan']
        self.assertTrue(audit['passed'])
        self.assertEqual(audit['substrate_cell_counts'], [22, 16, 3])
        self.assertEqual(audit['bar_cell_counts'], [5, 5, 13])
        self.assertEqual(audit['actuator_particle_count'], 25)
        for margin in audit['edge_margins_m'].values():
            self.assertAlmostEqual(margin, .0013, places=15)
        self.assertAlmostEqual(audit['required_edge_buffer_m'], .0013, places=15)
        self.assertAlmostEqual(audit['usable_axial_command_m'], .00425, places=15)
        self.assertAlmostEqual(audit['remaining_axial_command_margin_m'], .0007836, places=15)
        self.assertEqual(audit['replenishment'], 'FINITE_BAR_NO_AUTOMATIC_INSERTION')
        self.assertEqual(audit['domain_convergence'], 'NOT_TESTED')

    def test_too_small_rectangle_or_bar_inventory_rejected(self):
        bad = copy.deepcopy(self.config)
        bad['geometry']['substrate_x_bounds_m'][1] = .0065
        with self.assertRaisesRegex(ValueError, 'RECTANGLE_DOES_NOT_COVER'):
            model._resolve_configuration(bad)
        bad = copy.deepcopy(self.config)
        bad['parameter_overrides']['bar_length']['value'] = .0055
        with self.assertRaisesRegex(ValueError, 'FINITE_FEEDSTOCK'):
            model._resolve_configuration(bad)

    def test_explicit_rectangle_and_intent_must_match_actual_parameters(self):
        changes = [lambda c: c['geometry'].update(substrate_width_compatibility='SQUARE'),
                   lambda c: c['geometry'].update(substrate_z_bounds_m=[-.001, 0.]),
                   lambda c: c['geometry'].update(bar_axis_xy_m=[.001, 0.]),
                   lambda c: c['geometry'].update(substrate_y_bounds_m=[-.0035, .004]),
                   lambda c: c['geometry'].update(substrate_x_bounds_m=[-.004, .0069]),
                   lambda c: c['process_intent'].update(planned_axial_feed_before_traverse_m=.001),
                   lambda c: c['process_intent'].update(traverse_requires_observed_spreading=False),
                   lambda c: c['process_intent'].update(stationary_feed_end_s=.210)]
        for index, change in enumerate(changes):
            bad = copy.deepcopy(self.config); change(bad)
            with self.subTest(index=index), self.assertRaises(ValueError):
                model._resolve_configuration(bad)

    def test_geometry_hash_or_saved_budget_tampering_rejected(self):
        changes = [lambda s: s['geometry_definition']['substrate_x_bounds_m'].__setitem__(1, .008),
                   lambda s: s.update(geometry_definition_sha256='0'*64),
                   lambda s: s['independent_geometry_audit'].update(total_mass_kg=1.),
                   lambda s: s['domain_inventory_plan'].update(particle_count=1350),
                   lambda s: s['framework_config']['geometry'].update(substrate_y_bounds_m=[-.005, .005])]
        for index, change in enumerate(changes):
            bad = copy.deepcopy(self.state); change(bad)
            with self.subTest(index=index), self.assertRaises(ValueError):
                model.validate_state(bad)

    def test_initial_reference_coordinate_mass_or_label_changes_rejected(self):
        for key in ('x0', 'm', 'component_id', 'actuator', 'wholebar_tag'):
            pa = self.initial_pa.extract_particles(np.arange(len(self.initial_pa.x), dtype=np.uint32))
            values = getattr(pa, key); values[-1] += 1 if values.dtype.kind in 'iu' else .00001
            with self.subTest(field=key), self.assertRaises(ValueError):
                model.ProcessTrial(pa, self.state, .001)

    def test_old_framework_cannot_be_relabelled_as_rectangular_origin(self):
        _, oldstate = model.parent.independent_initial_state()
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, oldstate, .001)
        oldstate['classification'] = model.CLASSIFICATION
        with self.assertRaises(ValueError):
            model.ProcessTrial(self.pa, oldstate, .001)
        self.assertIsNone(self.state['parent_checkpoint'])
        self.assertFalse(self.state['independent_geometry_audit']['old_checkpoint_read'])
        self.assertFalse(self.state['independent_geometry_audit']['material_trimmed_or_repositioned'])

    def test_frozen_physical_operators_and_interval_constraint_are_unchanged(self):
        for name in ('evaluate', 'tangent', 'prescribed', '_require_interval'):
            self.assertIs(getattr(model.ProcessTrial, name), getattr(model.parent.ProcessTrial, name))
        audit = model.operator_reuse_audit()
        self.assertEqual(audit['status'], 'PASS_UNCHANGED_INHERITED_THERMOMECHANICAL_OPERATORS')
        self.assertFalse(audit['runtime_monkeypatch'])

    def test_existing_energy_contact_and_machine_prefix_preserved(self):
        # Manufactured bookkeeping, not a physical solution or publication.
        state = copy.deepcopy(self.state)
        step = model.validate_state(state).step(0., .001)
        state.update(time_s=step.t1, step=1, machine_history=[step.to_dict()])
        state['bore_observation']['last_time_s'] = step.t1
        state['interval_energy'] = require_interval_energy(state['interval_energy'],
            dict(global_raw_residual_J=1e-10, absolute_work_scale_J=1e-8))
        state['contact_history'][1056, 0, 0] = 1e-8
        previous_energy = copy.deepcopy(state['interval_energy'])
        previous_contact = state['contact_history'].copy()
        model.ProcessTrial(self.pa, state, .001)
        self.assertEqual(state['interval_energy'], previous_energy)
        np.testing.assert_array_equal(state['contact_history'], previous_contact)
        state['interval_energy'] = initial_interval_energy_state()
        with self.assertRaises(ValueError):
            model.validate_state(state)

    def test_constructor_never_promotes_independent_state_to_accepted_restart(self):
        before = copy.deepcopy(self.state['bore_observation'])
        model.ProcessTrial(self.pa, self.state, .001)
        self.assertEqual(self.state['bore_observation'], before)
        self.assertFalse(self.state['accepted'])
        self.assertFalse(self.state['restart_allowed'])
        self.assertTrue(self.state['diagnostic_only'])
        self.assertEqual(self.state['physical_qualification'], 'NOT_RUN')
        self.assertEqual(self.state['accepted_time_advanced_s'], 0.)


if __name__ == '__main__':
    unittest.main()

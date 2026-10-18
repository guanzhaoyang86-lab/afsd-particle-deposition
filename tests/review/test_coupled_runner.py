"""Runner/configuration/output gates. No integrated physics steps are executed."""
from pathlib import Path
import contextlib
import copy
import io
import json
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_coupled_afsd as r
from afsd_backend.pysph.geometry import surfaces, as_rs
from afsd_scaffold.vtk_legacy import read_point_cloud


class CoupledRunnerTests(unittest.TestCase):
    def setUp(self):
        self.cfg, self.p, self.resolved = r.settings()

    def test_01_frozen_true_3d_particle_inventory_and_material_units(self):
        pa = r.initial_particles(self.p)
        self.assertEqual(len(pa.x), 342)
        self.assertEqual(np.count_nonzero(pa.component_id == 0), 192)
        self.assertEqual(np.count_nonzero(pa.component_id == 1), 150)
        self.assertEqual(self.p['dx'], .0005)
        self.assertEqual(self.p['layer'], .0015)
        self.assertAlmostEqual(self.p['layer']/self.p['dx'], 3., delta=1e-14)
        self.assertEqual(self.p['shear'], 26.3e9)
        self.assertEqual(self.p['rho0'], 2700.)
        self.assertEqual(self.p['normal_stiffness'], 657500.)
        np.testing.assert_allclose(pa.m, np.full(342, 3.375e-7), rtol=0, atol=1e-21)
        np.testing.assert_array_equal(pa.rho, np.full(342, 2700.))
        x = r.positions(pa)
        for component in (0, 1):
            subset = x[pa.component_id == component]
            self.assertEqual(np.linalg.matrix_rank(subset-subset.mean(axis=0)), 3)
        self.assertEqual(len(np.unique(pa.z[pa.component_id == 0])), 3)
        self.assertEqual(len(np.unique(pa.z[pa.component_id == 1])), 6)

    def test_02_only_top_25_initial_z_velocities_are_prescribed(self):
        pa = r.initial_particles(self.p)
        expected = (pa.component_id == 1) & (pa.z0 > .0025)
        self.assertEqual(np.count_nonzero(expected), 25)
        np.testing.assert_array_equal(pa.actuator, expected.astype(int))
        np.testing.assert_array_equal(pa.u, np.zeros(342))
        np.testing.assert_array_equal(pa.v, np.zeros(342))
        np.testing.assert_array_equal(pa.w, np.where(expected, -.00212, 0.))
        command = r.control(pa.z0, 0., 0., .003, .0005, .00212, 1.)
        np.testing.assert_array_equal(command['free_axes'][:, :2], np.ones((342, 2)))
        np.testing.assert_array_equal(command['free_axes'][:, 2], (~expected).astype(float))
        # No integrated solver is imported or run to manufacture a geometry test.
        self.assertTrue(np.all(pa.S == 0.)); self.assertTrue(np.all(pa.ep == 0.))

    def test_03_square_bore_center_remains_open_and_base_is_true_plane(self):
        center = np.array([[0., 0., .002]])
        self.assertEqual(sum(len(a[0]) for a in surfaces(center, 0., self.p)), 0)
        self.assertAlmostEqual(self.p['bore_width'], .002496, delta=1e-15)
        self.assertEqual(self.p['substrate_depth'], .0015)
        pa = r.initial_particles(self.p)
        bottom = pa.z.min()
        self.assertAlmostEqual(bottom+self.p['substrate_depth'], .00025, delta=1e-17)
        self.assertFalse(self.cfg['production_authorized'])

    def test_04_new_parameters_have_explicit_units_sources_and_statuses(self):
        expected = dict(r.UNITS, feed='m/s', release_time='s')
        for key, entry in self.cfg['parameters'].items():
            self.assertEqual(entry['unit'], expected[key], key)
            self.assertTrue(entry['source']); self.assertTrue(entry['status'])
            self.assertTrue(np.isfinite(entry['value']))
        self.assertIn('parent_resolved_configuration', self.resolved)
        self.assertEqual(self.resolved['derived_normal_stiffness']['unit'], 'N/m')
        self.assertEqual(self.cfg['maximum_physics_attempts'], 1)
        self.assertLessEqual(self.cfg['external_wall_seconds'], 600)

    def test_05_invalid_units_production_and_unresolved_parameter_rejected(self):
        for attack in ('unit', 'production', 'source', 'nonfinite'):
            cfg = copy.deepcopy(self.cfg)
            if attack == 'unit': cfg['parameters']['feed']['unit'] = 'mm/s'
            elif attack == 'production': cfg['production_authorized'] = True
            elif attack == 'source': cfg['parameters']['feed']['source'] = ''
            else: cfg['parameters']['feed']['value'] = float('nan')
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory)/'config.json'; path.write_text(json.dumps(cfg), encoding='utf-8')
                with patch.object(r, 'CONFIG', path):
                    with self.subTest(attack=attack):
                        with self.assertRaises(ValueError): r.settings()

    def test_06_markers_and_as_rs_are_original_position_diagnostics(self):
        pa = r.initial_particles(self.p)
        outside = (pa.z0 < .002) | (pa.z0 > .003)
        np.testing.assert_array_equal(pa.tracer_type[outside], np.zeros(np.count_nonzero(outside)))
        self.assertGreater(np.count_nonzero(pa.wholebar_tag[outside]), 0)
        self.assertTrue(np.all(pa.material_id == 1))
        expected = np.where(pa.y0 < 0, 1, np.where(pa.y0 > 0, -1, 0))
        np.testing.assert_array_equal(pa.as_rs_tag, expected)
        np.testing.assert_array_equal(pa.as_rs_tag, as_rs(pa.y0, self.p['omega'], self.p['traverse']))
        np.testing.assert_array_equal(pa.particle_id, np.arange(1, 343))

    def probe_manifest(self, slowest=.1):
        return dict(status='COMPLETE_COST_PROBE', configuration_sha256=r.fingerprint(self.resolved),
            source_sha256={name: 'HASH' for name in r.SOURCES}, artifact_sha256={'history.json': 'HASH'},
            accepted_dt_s=[1e-8]*5, step_wall_seconds=[slowest]*5,
            peak_process_working_set_bytes=123456, accepted_steps=5)

    def select(self, manifest):
        with tempfile.TemporaryDirectory(dir=ROOT/'logs/coupled_afsd_20261002') as directory:
            path = Path(directory)/'probe.json'; path.write_text(json.dumps(manifest), encoding='utf-8')
            output = Path(directory)/'budget.json'
            with patch.object(r, 'sha', return_value='HASH'), contextlib.redirect_stdout(io.StringIO()):
                status = r.select_budget(path, output)
            return status, json.loads(output.read_text(encoding='utf-8'))

    def test_07_budget_chooses_longest_permitted_horizon_not_physics_outcome(self):
        code, report = self.select(self.probe_manifest(.1))
        self.assertEqual(code, 0)
        self.assertEqual(report['choice']['end_time_s'], 4e-6)
        self.assertEqual(report['choice']['estimated_steps'], 400)
        self.assertAlmostEqual(report['choice']['forecast_wall_seconds'], 90., delta=1e-12)
        code, report = self.select(self.probe_manifest(1.))
        self.assertEqual(code, 0); self.assertEqual(report['choice']['end_time_s'], 2e-6)
        self.assertAlmostEqual(report['choice']['forecast_wall_seconds'], 330., delta=1e-12)

    def test_08_cost_over_budget_returns_not_run_without_simulator_import(self):
        code, report = self.select(self.probe_manifest(2.))
        self.assertEqual(code, 2)
        self.assertEqual(report['status'], 'NOT_RUN_COST_EXCEEDS_BUDGET')
        self.assertIsNone(report['choice'])
        self.assertTrue(all(not e['permitted'] for e in report['estimates']))

    def test_09_stale_configuration_sources_and_artifacts_rejected(self):
        for attack in ('configuration', 'source', 'artifact', 'incomplete', 'traversal'):
            manifest = self.probe_manifest()
            if attack == 'configuration': manifest['configuration_sha256'] = 'wrong'
            elif attack == 'source': manifest['source_sha256'][r.SOURCES[0]] = 'wrong'
            elif attack == 'artifact': manifest['artifact_sha256']['history.json'] = 'wrong'
            elif attack == 'incomplete': manifest['status'] = 'BLOCKED'
            else: manifest['artifact_sha256'] = {'../history.json': 'HASH'}
            with self.subTest(attack=attack):
                with self.assertRaises(ValueError): self.select(manifest)

    def test_10_budget_output_cannot_overwrite_existing_evidence(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'logs/coupled_afsd_20261002') as directory:
            path = Path(directory)/'probe.json'; path.write_text(json.dumps(self.probe_manifest()), encoding='utf-8')
            output = Path(directory)/'budget.json'; output.write_text('ORIGINAL', encoding='utf-8')
            with patch.object(r, 'sha', return_value='HASH'):
                with self.assertRaises(FileExistsError): r.select_budget(path, output)
            self.assertEqual(output.read_text(encoding='utf-8'), 'ORIGINAL')

    def test_11_vtk_real_solved_fields_no_fictitious_certified_free_energy(self):
        pa = r.initial_particles(self.p)
        bulk = 2*26.3e9*(1+.33)/(3*(1-2*.33))
        solver = types.SimpleNamespace(pa=pa, p=self.p, bulk=bulk, count=0, time=0.)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'frame.vtk'; r.write_vtk(path, solver)
            data = read_point_cloud(path)
        np.testing.assert_array_equal(data['points'], r.positions(pa))
        np.testing.assert_array_equal(data['scalars']['temperature_K'], np.full(342, 723.15))
        np.testing.assert_array_equal(data['scalars']['particle_id'], pa.particle_id)
        np.testing.assert_allclose(data['scalars']['thermal_energy_from_room_J_per_kg'], np.full(342, 387000.), rtol=0, atol=1e-8)
        self.assertIn('stress_yz_Pa', data['scalars'])
        self.assertFalse(any('free_energy' in key or 'total_energy' in key for key in data['scalars']))
        self.assertNotIn('equivalent_deviatoric_strain_rate_per_s', data['scalars'])

    def test_12_rejected_nonfinite_output_never_becomes_accepted(self):
        with self.assertRaisesRegex(ValueError, 'nonfinite accepted'): r.clean({'value': np.nan})
        rejected = r.clean({'value': np.nan}, rejected=True)
        self.assertEqual(rejected['value']['classification'], 'REJECTED_DEBUG_VALUE')
        self.assertEqual(rejected['value']['nonfinite'], 'nan')

    def test_13_budget_stale_source_rejected_before_integrated_solver_instantiation(self):
        constructor = Mock(side_effect=AssertionError('physics must not be instantiated'))
        stub = types.ModuleType('afsd_backend.pysph.coupled_solver'); stub.IntegratedSolver = constructor
        receipt = dict(status='READY_FOR_ONE_REGISTERED_TRIAL', source_sha256={'changed': 'wrong'},
            configuration_sha256=r.fingerprint(self.resolved), choice={'end_time_s': 2e-6},
            probe_manifest_path='unused_mock_probe.json')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'budget.json'; path.write_text(json.dumps(receipt), encoding='utf-8')
            with patch.dict(sys.modules, {'afsd_backend.pysph.coupled_solver': stub}), patch.object(r, 'sha', return_value='HASH'), patch.object(r, 'budget_report', return_value=receipt):
                with self.assertRaisesRegex(ValueError, 'stale'): r.run(Path(directory)/'trial', budget=path)
            constructor.assert_not_called()

    def test_14_modified_choice_cannot_bypass_recomputed_budget(self):
        constructor = Mock(side_effect=AssertionError('physics must not be instantiated'))
        stub = types.ModuleType('afsd_backend.pysph.coupled_solver'); stub.IntegratedSolver = constructor
        with tempfile.TemporaryDirectory(dir=ROOT/'logs/coupled_afsd_20261002') as directory:
            probe = Path(directory)/'probe.json'
            probe.write_text(json.dumps(self.probe_manifest(1.)), encoding='utf-8')
            with patch.object(r, 'sha', return_value='HASH'):
                receipt = r.budget_report(probe)
                self.assertEqual(receipt['choice']['end_time_s'], 2e-6)
                receipt['choice'] = dict(end_time_s=4e-6, estimated_steps=400,
                    forecast_wall_seconds=1., permitted=True)
                path = Path(directory)/'budget.json'; path.write_text(json.dumps(receipt), encoding='utf-8')
                with patch.dict(sys.modules, {'afsd_backend.pysph.coupled_solver': stub}):
                    with self.assertRaisesRegex(ValueError, 'recomputed'): r.run(Path(directory)/'trial', budget=path)
            constructor.assert_not_called()

    def test_15_run_never_overwrites_existing_directory(self):
        constructor = Mock(side_effect=AssertionError('physics must not be instantiated'))
        stub = types.ModuleType('afsd_backend.pysph.coupled_solver'); stub.IntegratedSolver = constructor
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'evidence.txt'; path.write_text('ORIGINAL', encoding='utf-8')
            with patch.dict(sys.modules, {'afsd_backend.pysph.coupled_solver': stub}), patch.object(r, 'sha', return_value='HASH'):
                with self.assertRaises(FileExistsError): r.run(directory, probe=True)
            self.assertEqual(path.read_text(encoding='utf-8'), 'ORIGINAL')
            constructor.assert_not_called()

    def test_16_external_cost_manifest_cannot_be_used_as_project_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            probe = Path(directory)/'probe.json'; probe.write_text(json.dumps(self.probe_manifest()), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'within project'): r.budget_report(probe)


if __name__ == '__main__': unittest.main(verbosity=2)

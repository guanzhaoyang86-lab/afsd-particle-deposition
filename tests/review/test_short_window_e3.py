"""Read-only acceptance of the preregistered .5us window, not paper validation."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from inspect_engineering_integration import inspect
from afsd_backend.pysph.output import restore
from afsd_scaffold.vtk_legacy import read_point_cloud

DIRECTORY = ROOT/'outputs/short_window_e3_20261001/dx05_cfl01_05us'


class ShortWindowEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = inspect(DIRECTORY)
        cls.series = json.loads((DIRECTORY/'series.json').read_text())
        cls.diagnostics = [json.loads(s) for s in
                           (DIRECTORY/'diagnostics.jsonl').read_text().splitlines() if s.strip()]

    def test_preregistered_interval_complete_not_merely_old_step_crossed(self):
        r = self.report
        self.assertEqual(r['requested_end_time_s'], 5e-7)
        self.assertEqual(r['runner_status'], 'COMPLETE')
        self.assertTrue(r['reached_requested_end'])
        self.assertEqual(r['accepted_time_s'], 5e-7)
        self.assertGreater(r['accepted_step'], 21)
        self.assertGreater(r['accepted_time_s'], 2.2029619823669545e-7)
        for name in ('REJECTED_candidate.npz', 'HARD_TIMEOUT.json'):
            self.assertFalse((DIRECTORY/name).exists())

    def test_physics_configuration_and_solver_source_unchanged(self):
        prior = ROOT/'outputs/engineering_e3_20261001/dx05_cfl01/resolved_config.json'
        self.assertEqual(json.loads((DIRECTORY/'resolved_config.json').read_text()),
                         json.loads(prior.read_text()))
        frozen = json.loads((ROOT/'docs/CONTACT476_SCALE_FREEZE.json').read_text())['files_sha256']
        sources = json.loads((DIRECTORY/'source_hashes.json').read_text())
        for path, digest in sources.items():
            self.assertEqual(digest, frozen[path.replace('\\', '/')], path)
        manifest = json.loads((DIRECTORY/'manifest.json').read_text())
        self.assertIsNone(manifest['provenance']['checkpoint_parent'])
        self.assertEqual(manifest['provenance']['contact_variant'], 'event_aware_v1')

    def test_every_accepted_frame_is_finite_and_preserves_ids_labels_mass(self):
        frames = self.series['frames']
        self.assertEqual([f['step'] for f in frames], list(range(self.report['accepted_step']+1)))
        baseline = None
        diag = {d['step']: d for d in self.diagnostics}
        for frame in frames:
            data = read_point_cloud(DIRECTORY/frame['file'])
            scalar = {k: np.asarray(v) for k, v in data['scalars'].items()}
            vector = {k: np.asarray(v) for k, v in data['vectors'].items()}
            for values in (np.asarray(data['points']), *scalar.values(), *vector.values()):
                self.assertTrue(np.isfinite(values).all(), frame['file'])
            ids = scalar['particle_id'].astype(np.int64)
            order = np.argsort(ids)
            self.assertEqual(len(np.unique(ids)), 657)
            identity = np.column_stack([scalar[k][order] for k in
                ('particle_id', 'material_id', 'component_id', 'tracer_type', 'wholebar_tag', 'actuator')])
            if baseline is None:
                baseline = identity
            np.testing.assert_array_equal(identity, baseline)
            mass = scalar['mass_kg']
            expected_mass = 657*2700*(.0005**3)
            self.assertAlmostEqual(float(mass.sum()), expected_mass, delta=expected_mass*1e-12)
            if frame['step']:
                row = diag[frame['step']]
                self.assertEqual(row['alarms'], [])
                self.assertEqual(row['time_s'], frame['time_s'])
                np.testing.assert_allclose((mass[:, None]*vector['velocity_m_per_s']).sum(axis=0),
                    row['linear_momentum_kg_m_s'], rtol=1e-12, atol=1e-20)
                self.assertAlmostEqual(float(scalar['equivalent_plastic_strain'].max()),
                                       row['ep_max'], delta=1e-14)
                self.assertLessEqual(row['density_relative_change_max'], .05)
                self.assertLessEqual(row['max_penetration_m']/.0005, .1)
                self.assertLessEqual(row['speed_max_m_s'], 20.)

    def test_tool_pose_is_from_units_not_a_fitted_trajectory(self):
        for frame in self.series['frames']:
            t = frame['time_s']
            self.assertAlmostEqual(frame['tool_rotation_rad'], 10*np.pi*t, delta=1e-18)
            np.testing.assert_allclose(frame['tool_translation_m'], [.00212*t, 0., .0015],
                                       rtol=0., atol=1e-18)

    def test_recorded_energy_ledger_and_original_limits(self):
        # Independent arithmetic on saved data, NOT a constitutive-energy oracle.
        _, state = restore(DIRECTORY/'checkpoint.npz')
        for row in self.diagnostics:
            total = sum(row[k] for k in ('kinetic_J', 'elastic_shear_J', 'elastic_volume_J',
                                        'thermal_increment_J', 'contact_potential_J'))
            expected = state['initial_energy']+row['wall_work_J']+row['actuator_work_J']-row['wall_heat_out_J']
            residual = total-expected
            self.assertAlmostEqual(total, row['total_J'], delta=1e-18)
            self.assertAlmostEqual(residual, row['energy_residual_J'], delta=1e-18)
            scale = max(abs(state['initial_energy']), abs(row['wall_work_J'])+
                        abs(row['actuator_work_J'])+abs(row['wall_heat_out_J']), 1e-15)
            self.assertLessEqual(abs(residual), 1e-7+.05*scale)
            self.assertGreaterEqual(row['temperature_min_K'], 293.15)
            self.assertLessEqual(row['temperature_max_K'], 855.15)

    def test_checkpoint_contains_actual_solved_fields_and_accepted_state(self):
        pa, state = restore(DIRECTORY/'checkpoint.npz')
        last = read_point_cloud(DIRECTORY/self.series['frames'][-1]['file'])
        self.assertEqual(state['count'], self.report['accepted_step'])
        self.assertEqual(state['t'], self.report['accepted_time_s'])
        for prop, field in (('rho', 'density_kg_per_m3'), ('temperature', 'temperature_K'),
                            ('ep', 'equivalent_plastic_strain')):
            np.testing.assert_array_equal(getattr(pa, prop), last['scalars'][field])
        np.testing.assert_array_equal(np.column_stack((pa.u, pa.v, pa.w)),
                                      last['vectors']['velocity_m_per_s'])
        sigma = pa.S.reshape(-1, 3, 3)-np.asarray(last['scalars']['pressure_Pa'])[:, None, None]*np.eye(3)
        for a, b, suffix in ((0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'),
                             (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')):
            np.testing.assert_array_equal(sigma[:, a, b], last['scalars']['stress_'+suffix+'_Pa'])

    def test_short_success_or_failure_does_not_grant_scientific_validation(self):
        for name in ('production_authorized', 'paper_reproduction', 'deposition_validated',
                     'mixing_validated'):
            self.assertIs(self.report[name], False)
        self.assertLessEqual(self.report['performance']['wall_seconds'], 300.)


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Runner safety/reference tests; no registered AFSD attempt is launched."""
import copy
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts')); sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(Path(__file__).parent))
import probe_updated_stage as runner
from test_updated_trial import cloud
from afsd_backend.pysph.updated_trial import UpdatedTrial


def save_probe(directory, elapsed=1.):
    directory = Path(directory)
    names = ('resolved_config.json', 'last_accepted_reference.json', 'uncommitted_trial.json.gz')
    for name in names:
        (directory/name).write_bytes(b'INDEPENDENT_METADATA_FIXTURE_NOT_PHYSICS')
    manifest = dict(classification=runner.CLASSIFICATION,
        status='COMPLETE_COST_PROBE_UNCOMMITTED', mode='COST_PROBE',
        accepted_progress=False, production_authorized=False, settings=copy.deepcopy(runner.SETTINGS),
        source_sha256={'fixture': 'source_hash'}, configuration_sha256='config_hash',
        complete_trial_wall_seconds=elapsed,
        artifact_sha256={name: runner.sha(directory/name) for name in names})
    path = directory/'manifest.json'
    path.write_text(json.dumps(manifest), encoding='utf-8')
    return path, manifest


def call_budget(path):
    return runner.budget(path, {'fixture': 'source_hash'}, 'config_hash')


class RunnerSafetyTests(unittest.TestCase):
    def test_01_cost_arithmetic_and_budget_boundary_are_not_relaxed(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp:
            path, m = save_probe(temp, 1.)
            self.assertEqual(call_budget(path)['forecast_wall_seconds'], 78.)
            self.assertTrue(call_budget(path)['allowed'])
            m['complete_trial_wall_seconds'] = 8.75
            path.write_text(json.dumps(m), encoding='utf-8')
            self.assertEqual(call_budget(path)['forecast_wall_seconds'], 450.)
            self.assertTrue(call_budget(path)['allowed'])
            m['complete_trial_wall_seconds'] = 8.7501
            path.write_text(json.dumps(m), encoding='utf-8')
            self.assertFalse(call_budget(path)['allowed'])

    def test_02_stale_source_configuration_status_and_elapsed_fail_closed(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp:
            path, pristine = save_probe(temp)
            mutations = [('source_sha256', {'fixture': 'changed'}),
                ('configuration_sha256', 'changed'), ('status', 'BLOCKED_REJECTED_COST_PROBE'),
                ('complete_trial_wall_seconds', 0.), ('complete_trial_wall_seconds', -1.),
                ('complete_trial_wall_seconds', float('nan'))]
            for key, value in mutations:
                with self.subTest(key=key, value=value):
                    m = copy.deepcopy(pristine); m[key] = value
                    path.write_text(json.dumps(m), encoding='utf-8')
                    with self.assertRaises(ValueError): call_budget(path)

    def test_03_artifact_tampering_omission_and_path_escape_fail(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp:
            path, pristine = save_probe(temp)
            (Path(temp)/'resolved_config.json').write_text('CHANGED', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'artifact'): call_budget(path)
            for mapping in ({}, {'../escape.json': 'hash'}):
                m = copy.deepcopy(pristine); m['artifact_sha256'] = mapping
                path.write_text(json.dumps(m), encoding='utf-8')
                with self.assertRaises(ValueError): call_budget(path)

    def test_04_forged_probe_metadata_cannot_authorize_real_attempt(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp:
            path, pristine = save_probe(temp)
            bad_settings = dict(runner.SETTINGS, dt_s=1.)
            for key, value in [('classification', 'GEOMETRY_FIXTURE'),
                               ('mode', 'SINGLE_GLOBAL_ATTEMPT'), ('accepted_progress', True),
                               ('production_authorized', True), ('settings', bad_settings)]:
                with self.subTest(key=key):
                    m = copy.deepcopy(pristine); m[key] = value
                    path.write_text(json.dumps(m), encoding='utf-8')
                    with self.assertRaises(ValueError): call_budget(path)

    def test_05_existing_output_and_nonproject_target_refused_before_loading(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp:
            sentinel = Path(temp)/'keep.txt'; sentinel.write_text('KEEP', encoding='utf-8')
            with patch.object(runner, 'inputs', side_effect=AssertionError('must not load')) as load:
                with self.assertRaisesRegex(ValueError, 'fresh project-local'): runner.run(temp, probe=True)
                with self.assertRaisesRegex(ValueError, 'fresh project-local'):
                    runner.run(ROOT.parent/'outside_registered_project', probe=True)
                load.assert_not_called()
            self.assertEqual(sentinel.read_text(encoding='utf-8'), 'KEEP')

    def test_06_actual_zero_equilibrium_is_accepted_only_as_approximate_stage(self):
        pa, p = cloud(); trial = UpdatedTrial(pa, p, 0.); c = trial.candidate()
        e = trial.evaluate(c, 1e-9)
        checks = runner.acceptance(pa, p, trial, c, e)
        self.assertEqual(checks['maximum_scaled_residual'], 0.)
        self.assertEqual(checks['thermal_residual_J'], 0.)
        self.assertEqual(checks['delta_normal_potential_J'], 0.)
        self.assertEqual(checks['delta_K_J'], 0.)
        self.assertFalse(checks['production_authorized'])
        self.assertIn('NOT_ENERGY_CERTIFIED', checks['status'])

    def test_07_acceptance_rejects_fake_metadata_identity_sources_and_heat(self):
        pa, p = cloud(); trial = UpdatedTrial(pa, p, 0.); c = trial.candidate()
        pristine = trial.evaluate(c, 1e-9)
        for key, value in [('classification', 'GEOMETRY_FIXTURE'), ('accepted', True),
                           ('restart_allowed', True), ('material_scope', 'MANUFACTURED_TEST_ONLY')]:
            with self.subTest(key=key):
                e = copy.deepcopy(pristine); e[key] = value
                with self.assertRaises(ValueError): runner.acceptance(pa, p, trial, c, e)
        e = copy.deepcopy(pristine); e['particle_id'][0] += 100
        with self.assertRaisesRegex(ValueError, 'ID'): runner.acceptance(pa, p, trial, c, e)
        e = copy.deepcopy(pristine); e['thermal']['particle_energy_J']['total_heat_input_J'][0] = 1.
        with self.assertRaisesRegex(ValueError, 'thermal balance'): runner.acceptance(pa, p, trial, c, e)
        e = copy.deepcopy(pristine); interval = e['boundary']['traces'][0]['result']['intervals'][0]
        interval['friction_heat_raw_J'] = np.array([-1.]); interval['friction_rounding_bound_J'] = np.array([0.])
        with self.assertRaisesRegex(ValueError, 'negative raw friction'): runner.acceptance(pa, p, trial, c, e)

    def test_08_provenance_failure_is_saved_rejected_without_snapshot(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp:
            output = Path(temp)/'fresh_rejected'
            with patch.object(runner, 'inputs', side_effect=ValueError('INJECTED_SOURCE_DRIFT')):
                code = runner.run(output, probe=True)
            self.assertEqual(code, 2)
            m = json.loads((output/'manifest.json').read_text(encoding='utf-8'))
            self.assertEqual(m['status'], 'BLOCKED_REJECTED_COST_PROBE')
            self.assertFalse(m['accepted_progress'])
            self.assertFalse((output/'accepted_snapshot.npz').exists())
            with gzip.open(output/'REJECTED_candidate.json.gz', 'rt', encoding='utf-8') as stream:
                rejected = json.load(stream)
            self.assertFalse(rejected['accepted']); self.assertFalse(rejected['restart_allowed'])
            self.assertEqual(rejected['classification'], 'REJECTED')
            self.assertIn('INJECTED_SOURCE_DRIFT', rejected['error'])


if __name__ == '__main__': unittest.main(verbosity=2)

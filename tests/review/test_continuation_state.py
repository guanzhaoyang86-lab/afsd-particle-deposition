"""History/storage regressions, not manufactured AFSD progress or restart dynamics.

Zero-equilibrium tests call the actual 27-point UpdatedTrial and original numerical
acceptance. The explicitly named manufactured increment test changes only stored
increments to test accounting; it is not a constitutive or energy benchmark.
"""
import copy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'scripts'))
from test_updated_trial import cloud
from probe_updated_stage import acceptance
from afsd_backend.pysph import continuation_state as storage
from afsd_backend.pysph.output import checkpoint
from afsd_backend.pysph.updated_trial import UpdatedTrial, particle_hash


SOURCES = {'MANUFACTURED_STORAGE_TEST_ONLY': 'a'*64}
CONFIG = 'b'*64


def evaluated_zero():
    pa, p = cloud()
    # Real input snapshots carry this passive property; the shared cloud does not.
    pa.add_property('as_rs_tag', type='long', data=np.arange(len(pa.x)) % 3-1)
    pa.set_output_arrays(list(pa.properties))
    trial = UpdatedTrial(pa, p, 0.)
    candidate = trial.candidate(); evaluation = trial.evaluate(candidate, 1e-9)
    checks = acceptance(pa, p, trial, candidate, evaluation)
    return pa, p, candidate, evaluation, checks


def metadata(pa, p):
    return storage.state(pa, p, 1e-9, 1, SOURCES, CONFIG,
        {'fixture': 'ZERO_EQUILIBRIUM_STORAGE_TEST_NOT_AFSD'}, {},
        {'wall_work_J': 0., 'plastic_work_J': 0.})


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ContinuationStateTests(unittest.TestCase):
    def test_real_zero_trial_acceptance_commit_does_not_mutate_old(self):
        pa, p, c, e, checks = evaluated_zero(); before = particle_hash(pa)
        committed = storage.commit(pa, p, c, e, checks)
        self.assertEqual(checks['status'], 'ACCEPTED_NUMERICAL_STAGE_APPROXIMATE_NOT_ENERGY_CERTIFIED')
        self.assertEqual(particle_hash(pa), before)
        for key in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'temperature', 'ep', 'S'):
            np.testing.assert_array_equal(getattr(committed, key), getattr(pa, key))
        committed.temperature[0] = 777.
        self.assertEqual(particle_hash(pa), before)

    def test_nonaccepted_wrong_history_and_id_correspondence_rejected(self):
        pa, p, c, e, checks = evaluated_zero(); before = particle_hash(pa)
        for bad in ({}, {'status': 'REJECTED'}, {'status': 'ROOT_CONVERGED_NOT_ACCEPTED'}):
            with self.assertRaisesRegex(ValueError, 'acceptance'): storage.commit(pa, p, c, e, bad)
        wrong = copy.deepcopy(e); wrong['accepted_state_hash'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'history'): storage.commit(pa, p, c, wrong, checks)
        wrong = copy.deepcopy(e); wrong['particle_id'] = wrong['particle_id'][::-1]
        with self.assertRaisesRegex(ValueError, 'ID'): storage.commit(pa, p, c, wrong, checks)
        self.assertEqual(particle_hash(pa), before)

    def test_manufactured_history_increments_stored_exactly_once_not_a_physics_test(self):
        pa, p, c, e, checks = evaluated_zero(); n = len(pa.x)
        pa.ep[:] = .01; pa.plastic_work[:] = np.arange(n)*1e-8
        pa.return_dissipation[:] = np.arange(n)*1e-9
        before = particle_hash(pa)
        e = copy.deepcopy(e); e['accepted_state_hash'] = before
        added = (np.arange(n)+1)*2e-9
        returned = (np.arange(n)+1)*3e-10
        e['material']['ep'] = pa.ep+.003
        e['material']['plastic_work_J'] = added
        e['material']['return_dissipation_J'] = returned
        # A storage-only injected fixture, deliberately NOT an accepted FKS run.
        committed = storage.commit(pa, p, c, e, checks)
        np.testing.assert_array_equal(committed.ep, e['material']['ep'])
        np.testing.assert_allclose(committed.ep, np.full(n, .013), rtol=1e-15)
        np.testing.assert_allclose(committed.ep_rate, np.full(n, 3e6), rtol=1e-15)
        np.testing.assert_array_equal(committed.plastic_work, pa.plastic_work+added)
        np.testing.assert_array_equal(committed.return_dissipation, pa.return_dissipation+returned)
        repeated_from_same_old = storage.commit(pa, p, c, e, checks)
        self.assertEqual(particle_hash(repeated_from_same_old), particle_hash(committed))
        with self.assertRaisesRegex(ValueError, 'history'):
            storage.commit(committed, p, c, e, checks)
        self.assertEqual(particle_hash(pa), before)

    def test_checkpoint_all_valid_fields_roundtrip_and_unsolved_energy_not_published(self):
        pa, p, c, e, checks = evaluated_zero(); committed = storage.commit(pa, p, c, e, checks)
        state = metadata(committed, p)
        with tempfile.TemporaryDirectory(prefix='afsd_storage_test_') as tmp:
            directory = Path(tmp); storage.save(directory, committed, p, state)
            restored, restored_state = storage.restore(directory/'accepted.npz', p, SOURCES, CONFIG,
                                                       file_hash(directory/'accepted.npz'))
            self.assertEqual(restored_state, state)
            for key in storage.VALID_PROPERTIES:
                np.testing.assert_array_equal(getattr(restored, key), getattr(committed, key), err_msg=key)
            self.assertEqual(particle_hash(restored), particle_hash(committed))
            vtk = (directory/'accepted.vtk').read_text(encoding='utf-8')
            self.assertNotIn('internal_energy', vtk)
            self.assertIn('NOT_DEPOSITION_VALIDATION', vtk)
            self.assertIn('NOT_YET_CERTIFIED', state['restart_scope'])

    def test_file_hash_source_config_and_parameter_mismatches_rejected(self):
        pa, p, c, e, checks = evaluated_zero(); committed = storage.commit(pa, p, c, e, checks)
        with tempfile.TemporaryDirectory(prefix='afsd_storage_test_') as tmp:
            path = Path(tmp)/'accepted.npz'; checkpoint(path, committed, metadata(committed, p)); digest = file_hash(path)
            wrong_p = dict(p, shear=p['shear']*1.001)
            cases = [(p, SOURCES, CONFIG, '0'*64), (p, {'wrong': 'a'*64}, CONFIG, digest),
                     (p, SOURCES, 'wrong-config', digest), (wrong_p, SOURCES, CONFIG, digest)]
            for parameters, sources, config, expected_hash in cases:
                with self.subTest(sources=sources, config=config, expected_hash=expected_hash):
                    with self.assertRaisesRegex(ValueError, 'mismatch'):
                        storage.restore(path, parameters, sources, config, expected_hash)

    def test_rejected_and_tampered_snapshot_metadata_cannot_restore(self):
        pa, p, c, e, checks = evaluated_zero(); committed = storage.commit(pa, p, c, e, checks)
        base = metadata(committed, p)
        changes = [('classification', 'REJECTED'), ('accepted', False), ('restart_allowed', False),
                   ('production_authorized', True), ('particle_state_sha256', '0'*64),
                   ('particle_ids', base['particle_ids'][::-1]), ('valid_properties', ['x'])]
        with tempfile.TemporaryDirectory(prefix='afsd_storage_test_') as tmp:
            for index, (key, value) in enumerate(changes):
                state = copy.deepcopy(base); state[key] = value
                path = Path(tmp)/f'invalid_{index}.npz'; checkpoint(path, committed, state)
                with self.subTest(key=key):
                    with self.assertRaises(ValueError): storage.restore(path, p, SOURCES, CONFIG, file_hash(path))

    def test_particle_permutation_preserves_id_correspondence_across_commit_snapshot(self):
        pa, p, c, e, checks = evaluated_zero()
        permutation = np.arange(len(pa.x)-1, -1, -1, dtype=np.uint32)
        permuted = pa.extract_particles(permutation)
        trial = UpdatedTrial(permuted, p, 0.); candidate = trial.candidate()
        evaluation = trial.evaluate(candidate, 1e-9)
        accepted = acceptance(permuted, p, trial, candidate, evaluation)
        committed = storage.commit(permuted, p, candidate, evaluation, accepted)
        reference = {int(pa.particle_id[i]): (pa.x[i], pa.y[i], pa.z[i], pa.tracer_type[i], pa.as_rs_tag[i])
                     for i in range(len(pa.x))}
        with tempfile.TemporaryDirectory(prefix='afsd_storage_test_') as tmp:
            path = Path(tmp)/'accepted.npz'; checkpoint(path, committed, metadata(committed, p))
            restored, state = storage.restore(path, p, SOURCES, CONFIG, file_hash(path))
            actual = {int(restored.particle_id[i]): (restored.x[i], restored.y[i], restored.z[i],
                      restored.tracer_type[i], restored.as_rs_tag[i]) for i in range(len(restored.x))}
            self.assertEqual(reference, actual)
            np.testing.assert_array_equal(restored.particle_id, pa.particle_id[::-1])
            self.assertEqual(state['particle_ids'], restored.particle_id.tolist())


if __name__ == '__main__': unittest.main(verbosity=2)

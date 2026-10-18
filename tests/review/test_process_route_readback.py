"""Literal readback references; no solver advance or reference via step_metrics."""
from pathlib import Path
import copy
import io
import json
import math
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import assess_process_route as reader


def example():
    """Two particles, two faces, two contact intervals; all values prescribed."""
    intervals = [
        dict(dt_s=.2, normal_impulse=[[.1, 0., 0.], [0., .2, 0.]],
             tangential_impulse=[[0., .03, -.01], [.02, 0., .01]],
             constraint_impulse=[0., 0., -.8]),
        dict(dt_s=.3, normal_impulse=[[.3, 0., 0.], [0., .4, 0.]],
             tangential_impulse=[[0., .07, -.03], [.08, 0., .03]],
             constraint_impulse=[0., 0., -1.2])]
    trace0 = dict(particle_index=0, particle_id=91, free_axes=[1., 1., 0.],
        prescribed_velocity=[0., 0., .3], contact_rows=[dict(feature_id=2), dict(feature_id=4)],
        result=dict(intervals=intervals, normal_impulse=[[.4, 0., 0.], [0., .6, 0.]],
                    tangential_impulse=[[0., .1, -.04], [.1, 0., .04]],
                    constraint_impulse=[0., 0., -2.]))
    trace1 = dict(particle_index=1, particle_id=17, free_axes=[1., 1., 1.],
        prescribed_velocity=[0., 0., 0.], contact_rows=[],
        result=dict(intervals=[dict(dt_s=.5, normal_impulse=[], tangential_impulse=[],
                                    constraint_impulse=[0., 0., 0.])],
                    normal_impulse=[], tangential_impulse=[], constraint_impulse=[0., 0., 0.]))
    return dict(step=1, stage='ACCEPTED_COUPLED_DRIFT_KICK_DRIFT', accepted_time_s=.5,
        dt_s=.5, mass_kg=[2., 1.], particle_id=[91, 17],
        x_before_m=[[0., 0., 0.], [1., 1., 1.]],
        x_after_m=[[.2875, -.025, .15], [.65, 1.325, .45]],
        v_before_m_s=[[.2, -.1, .3], [-.2, .4, -.1]],
        v_after_m_s=[[.95, 0., .3], [-1.2, .9, -2.1]],
        internal_force_N=[[2., -1., 4.], [-2., 1., -4.]],
        ep_before=[0., 0.], ep_after=[0., .01], traces=[trace0, trace1],
        increments=dict(wall_work_J=1.2, normal_relative_work_J=-.791,
            friction_raw_J=.14275, plastic_heat_J=.02,
            midpoint_stress_work_J=.4, material_old_work_J=.3))


class ProcessRouteReadbackTests(unittest.TestCase):
    def test_01_all_faces_intervals_and_reaction_including_contact_free_particle(self):
        step = example(); saved = copy.deepcopy(step)
        jn, jt, jr, free, prescribed, count = reader.reconstruct_impulses(step)
        np.testing.assert_allclose(jn, [[.4, .6, 0.], [0., 0., 0.]], atol=1e-15, rtol=0)
        np.testing.assert_allclose(jt, [[.1, .1, 0.], [0., 0., 0.]], atol=1e-15, rtol=0)
        np.testing.assert_allclose(jr, [[0., 0., -2.], [0., 0., 0.]], atol=1e-15, rtol=0)
        np.testing.assert_array_equal(free, [[1., 1., 0.], [1., 1., 1.]])
        np.testing.assert_array_equal(prescribed, [[0., 0., .3], [0., 0., 0.]])
        self.assertEqual(count, 3)
        self.assertEqual(step, saved)

    def test_02_inertia_balances_quasistatic_residual_not_zero_by_cancellation(self):
        r = reader.step_metrics(example())
        # Free unbalanced forces: (3,.4,0), (-2,1,-4) N. Sum is NOT the norm.
        self.assertAlmostEqual(r['quasistatic_free_residual_norm_N'], math.sqrt(30.16), delta=1e-12)
        self.assertAlmostEqual(r['inertial_free_force_norm_N'], math.sqrt(30.16), delta=1e-12)
        denominator = math.sqrt(26.)+math.sqrt(2.08)+math.sqrt(.08)
        self.assertAlmostEqual(r['force_denominator_N'], denominator, delta=1e-12)
        self.assertAlmostEqual(r['quasistatic_force_ratio'], math.sqrt(30.16)/denominator, delta=1e-12)
        self.assertLess(r['dynamic_residual_norm_N'], 1e-12)
        self.assertLess(r['inertia_static_balance_defect_N'], 1e-12)
        self.assertEqual(r['dominant_unbalanced_particle_id'], 17)
        self.assertEqual(r['new_plastic_particle_ids'], [17])
        self.assertEqual(r['free_axis_reaction_force_norm_N'], 0.)

    def test_03_event_work_is_not_whole_step_endpoint_surrogate(self):
        r = reader.step_metrics(example())
        # Contact-interval midpoint velocities (.33,-.0925,.3),(.705,-.0425,.3)
        # give .018325+.247925=.26625 J. Whole-step midpoint gives .2525 J.
        self.assertAlmostEqual(r['actual_contact_work_J'], .26625, delta=1e-12)
        self.assertAlmostEqual(r['endpoint_contact_work_surrogate_J'], .2525, delta=1e-12)
        self.assertAlmostEqual(r['endpoint_minus_actual_contact_work_J'], -.01375, delta=1e-12)
        self.assertEqual(r['surrogate_status'], 'NOT_FRICTION_HEAT_NOT_PHYSICAL_SUBINTERVAL_WORK')
        self.assertEqual(r['plastic_heat_J'], .02)

    def test_04_missing_duplicate_or_wrong_identity_trace_is_rejected(self):
        for mode in ('missing', 'duplicate', 'wrong_id', 'duplicate_ids'):
            step = example()
            if mode == 'missing': step['traces'].pop()
            elif mode == 'duplicate': step['traces'].append(copy.deepcopy(step['traces'][0]))
            elif mode == 'wrong_id': step['traces'][0]['particle_id'] = 17
            else: step['particle_id'] = [91, 91]
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError): reader.reconstruct_impulses(step)

    def test_05_bad_interval_coverage_or_face_inventory_is_rejected(self):
        for mode in ('missing_interval', 'negative_dt', 'missing_face', 'nan_impulse'):
            step = example(); trace = step['traces'][0]
            if mode == 'missing_interval': trace['result']['intervals'].pop()
            elif mode == 'negative_dt': trace['result']['intervals'][0]['dt_s'] = -.2
            elif mode == 'missing_face': trace['result']['intervals'][0]['normal_impulse'].pop()
            else: trace['result']['intervals'][0]['normal_impulse'][0][0] = float('nan')
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError): reader.reconstruct_impulses(step)

    def test_06_facewise_totals_cannot_be_swapped_even_if_resultant_is_identical(self):
        step = example(); result = step['traces'][0]['result']
        result['normal_impulse'] = list(reversed(result['normal_impulse']))
        with self.assertRaises(ValueError): reader.reconstruct_impulses(step)

    def test_07_total_impulse_and_reaction_mismatch_are_rejected(self):
        for key in ('normal_impulse', 'tangential_impulse', 'constraint_impulse'):
            step = example(); result = step['traces'][0]['result']
            if key == 'constraint_impulse': result[key][0] += .01
            else: result[key][0][0] += .01
            with self.subTest(key=key):
                with self.assertRaises(ValueError): reader.reconstruct_impulses(step)

    def test_08_bad_stage_dynamics_drift_or_constraint_cannot_pass(self):
        for mode in ('rejected', 'momentum', 'drift', 'constraint'):
            step = example()
            if mode == 'rejected': step['stage'] = 'REJECTED'
            elif mode == 'momentum': step['internal_force_N'][0][0] += .1
            elif mode == 'drift': step['x_after_m'][0][0] += .01
            else: step['traces'][0]['prescribed_velocity'][2] += .1
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError): reader.step_metrics(step)

    def test_09_zero_force_ratio_is_undefined_not_false_equilibrium_certificate(self):
        step = example()
        for key in ('internal_force_N', 'v_before_m_s', 'v_after_m_s'):
            step[key] = [[0., 0., 0.], [0., 0., 0.]]
        step['x_after_m'] = copy.deepcopy(step['x_before_m'])
        for trace in step['traces']:
            trace['prescribed_velocity'] = [0., 0., 0.]
            for value in [trace['result'], *trace['result']['intervals']]:
                value['normal_impulse'] = [[0., 0., 0.] for _ in trace['contact_rows']]
                value['tangential_impulse'] = copy.deepcopy(value['normal_impulse'])
                value['constraint_impulse'] = [0., 0., 0.]
        step['increments'] = {key: 0. for key in step['increments']}
        r = reader.step_metrics(step)
        self.assertIsNone(r['quasistatic_force_ratio'])
        self.assertEqual(r['ratio_status'], 'UNDEFINED_NO_FORCE')

    def test_10_missing_accepted_step_rejected_before_summary(self):
        first = example(); third = example(); third['step'] = 3
        manifest = dict(final_record=dict(particle_id=[91, 17]))
        stream = io.StringIO(json.dumps(first)+'\n'+json.dumps(third)+'\n')
        with patch.object(reader, 'checked_inputs', return_value=(manifest, {})), \
             patch.object(reader.gzip, 'open', return_value=stream):
            with self.assertRaisesRegex(ValueError, 'gap/duplicate'):
                reader.assess(reader.CASE)

    def test_11_changed_source_hash_rejected_before_data_analysis(self):
        case = ROOT/'outputs/manufactured_readback_test'
        freeze = ROOT/'docs/MANUFACTURED_UNUSED_FREEZE.json'
        manifest_path = case/'manifest.json'
        rel = manifest_path.relative_to(ROOT).as_posix()
        frozen = dict(files_sha256={rel: 'manifest_hash', 'src/manufactured.py': 'source_hash'})
        manifest = dict(accepted_steps=67, accepted_time_s=7e-7,
            status='COMPLETE_REGISTERED_INTERVAL_NOT_SUSTAINED_DEPOSITION', production_authorized=False,
            source_sha256={'src/manufactured.py': 'source_hash'}, artifact_sha256={})
        def text(path, **unused):
            return json.dumps(frozen if path == freeze else manifest)
        def digest(path):
            return ('manifest_hash' if path == manifest_path else
                    'freeze_hash' if path == freeze else 'CHANGED_SOURCE_HASH')
        with patch.object(reader, 'CASE', case), patch.object(reader, 'FREEZE', freeze), \
             patch.object(Path, 'read_text', text), patch.object(reader, 'sha', side_effect=digest):
            with self.assertRaisesRegex(ValueError, 'source drift'):
                reader.checked_inputs(case)

    def test_12_trace_order_does_not_change_particle_assignment(self):
        step = example(); step['traces'].reverse()
        jn, jt, reaction, _, _, count = reader.reconstruct_impulses(step)
        np.testing.assert_allclose(jn, [[.4, .6, 0.], [0., 0., 0.]], atol=1e-15, rtol=0)
        np.testing.assert_allclose(jt, [[.1, .1, 0.], [0., 0., 0.]], atol=1e-15, rtol=0)
        np.testing.assert_allclose(reaction, [[0., 0., -2.], [0., 0., 0.]], atol=1e-15, rtol=0)
        self.assertEqual(count, 3)

    def test_13_relative_cli_case_is_normalized_before_reading(self):
        relative = Path('outputs/coupled_afsd_20261002/short_07us')
        expected = relative.resolve()
        observed = []
        def checked(path):
            observed.append(path)
            self.assertTrue(path.is_absolute())
            self.assertEqual(path, expected)
            raise RuntimeError('NORMALIZED_PATH_SENTINEL_NO_DATA_READ')
        with patch.object(reader, 'checked_inputs', side_effect=checked):
            with self.assertRaisesRegex(RuntimeError, 'NORMALIZED_PATH_SENTINEL_NO_DATA_READ'):
                reader.assess(relative)
        self.assertEqual(observed, [expected])

    def _read_upstream_manifest_fixture(self, current_upstream_hash):
        case = ROOT/'outputs/manufactured_upstream_readback_test'
        freeze = ROOT/'docs/MANUFACTURED_UNUSED_UPSTREAM_FREEZE.json'
        manifest_path = case/'manifest.json'
        rel = manifest_path.relative_to(ROOT).as_posix()
        upstream = 'third_party/PySPH/pysph/sph/integrator.py'
        frozen = dict(files_sha256={rel: 'frozen_manifest_hash'})
        manifest = dict(accepted_steps=67, accepted_time_s=7e-7,
            status='COMPLETE_REGISTERED_INTERVAL_NOT_SUSTAINED_DEPOSITION', production_authorized=False,
            source_sha256={upstream: 'PINNED_UPSTREAM_HASH'}, artifact_sha256={})
        self.assertNotIn(upstream, frozen['files_sha256'])
        def text(path, **unused):
            return json.dumps(frozen if path == freeze else manifest)
        def digest(path):
            if path == manifest_path: return 'frozen_manifest_hash'
            if path == freeze: return 'upstream_fixture_freeze_hash'
            self.assertEqual(path, ROOT/upstream)
            return current_upstream_hash
        with patch.object(reader, 'CASE', case), patch.object(reader, 'FREEZE', freeze), \
             patch.object(Path, 'read_text', text), patch.object(reader, 'sha', side_effect=digest):
            return reader.checked_inputs(case)

    def test_14_upstream_hash_transitively_bound_by_frozen_manifest_is_allowed(self):
        manifest, inputs = self._read_upstream_manifest_fixture('PINNED_UPSTREAM_HASH')
        upstream = 'third_party/PySPH/pysph/sph/integrator.py'
        self.assertEqual(manifest['source_sha256'][upstream], 'PINNED_UPSTREAM_HASH')
        self.assertEqual(inputs[upstream], 'PINNED_UPSTREAM_HASH')

    def test_15_transitive_upstream_actual_hash_drift_is_still_rejected(self):
        with self.assertRaisesRegex(ValueError, 'source drift third_party/PySPH/'):
            self._read_upstream_manifest_fixture('ACTUAL_UPSTREAM_HASH_CHANGED')


if __name__ == '__main__': unittest.main(verbosity=2)

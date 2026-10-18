"""Independent administrative negative checks. No physical solves or GPU launch."""
import ast
import copy
from datetime import timedelta
import inspect
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_event_sensitivity_20261010 as subject
import check_gpu_event_sensitivity_20261010 as checker


class SensitivityTests(unittest.TestCase):
    def fixture(self):
        pa = SimpleNamespace(particle_id=np.array([11, 12, 13]), m=np.array([1., 2., 3.]),
            x=np.array([0., .01, .02]), y=np.zeros(3), z=np.ones(3), temperature=np.ones(3)*600.,
            S=np.zeros(27), ep=np.zeros(3))
        parameters = dict(dx=.0005)
        state = dict(time_s=.071, parameters=parameters,
            bore_observation=dict(crossing_ids=[11], particle_ids=[11, 12, 13]))
        common = dict(source_sha256={'frozen-physical': 'a'*64},
            screened_continuation_checkpoint='fixture/common.npz', screened_continuation_sha256='b'*64,
            screened_end_time_s=.071, maximum_candidate_wall_s=40., candidate_wall_s=[30., 40., 35.],
            event_protocol=subject.frozen.PROTOCOL, dt_s=.0004)
        parent = subject.COMMON_PARENT
        saved = {subject.ROOT/common['screened_continuation_checkpoint']: (pa, state)}
        documents = {}; stages = {}
        for scale, label in ((1., 'one'), (.5, 'half'), (.25, 'quarter')):
            folder = subject.BASE/('fixture_'+label)
            meta = dict(status=subject.frozen.CONTINUE_PASS, screened_continuation_authorized=True,
                source_sha256=dict(common['source_sha256'], **{'secant-numeric': '8'*64}),
                parent_implementation_source_sha256=copy.deepcopy(common['source_sha256']),
                parent_gpu_directory=subject.rel(parent), parent_checkpoint=common['screened_continuation_checkpoint'],
                parent_sha256=common['screened_continuation_sha256'], parent_time_s=.071, dt_s=.0004,
                screened_end_time_s=.0714, event_protocol=subject.frozen.PROTOCOL,
                event_tolerance_scale=scale, event_locator_protocol=subject.LOCATOR_PROTOCOL, implementation_transition=None,
                screened_continuation_sha256=({'one':'c','half':'d','quarter':'e'}[label])*64,
                event_records={r: [] for r in subject.ROLES}, maximum_candidate_wall_s=45.,
                candidate_wall_s=[45., 40., 42.])
            for role in subject.ROLES:
                ps = copy.deepcopy(state); ps['time_s'] = .0712 if role == 'fine_1' else .0714
                saved[folder/role/'NOT_RESTARTABLE_candidate.npz'] = (copy.deepcopy(pa), ps)
            screen = subject.compare_endpoint(saved[folder/'fine_2/NOT_RESTARTABLE_candidate.npz'],
                saved[folder/'coarse/NOT_RESTARTABLE_candidate.npz'], pa.particle_id, pa.m, parameters)['screen']
            meta['pair_screen'] = screen
            screen['fine_checkpoint_sha256'] = ['9'*64, meta['screened_continuation_sha256']]
            tol = 1e-7*scale
            lo_path = folder/'lo.npz'; hi_path = folder/'hi.npz'; location = folder/'location.json'
            low = copy.deepcopy(state); low['time_s'] = .0711-tol/2
            high = copy.deepcopy(state); high['time_s'] = .0711+tol/2
            saved[lo_path] = (copy.deepcopy(pa), low); saved[hi_path] = (copy.deepcopy(pa), high)
            documents[location] = dict(status='BRACKETED_ACTIVATION',
                lo=dict(blocked=True, checkpoint=subject.rel(lo_path)),
                hi=dict(blocked=False, checkpoint=subject.rel(hi_path)))
            meta['event_records']['coarse'].append(dict(location=subject.rel(location), location_sha256='f'*64,
                tolerance_s=tol, incoming_pid_pairs=[[11, 12]]))
            stages[scale] = (folder, meta)
        return parent, common, stages, saved, documents

    def group(self, fixture):
        parent, common, stages, saved, documents = fixture
        return subject.compare_group(common, parent, stages, loader=saved.__getitem__, reader=documents.__getitem__)

    def cost(self, value=200., branch=100., events=1):
        candidate = [branch, (value-branch)/2, (value-branch)/2]
        return dict(maximum_candidate_wall_s=branch, includes_event_probes=True,
            candidate_wall_s=candidate, candidate_total_wall_s=value, worker_wall_s=value,
            supervised_wall_s=value, measured_pair_wall_s=value, event_count=events,
            parent_sha256='b'*64, dt_s=.0004)

    def test_three_pairwise_comparisons_cover_every_macro_endpoint(self):
        result = self.group(self.fixture())
        self.assertEqual(result['status'], 'PASS_SAME_PARENT_EVENT_TOLERANCE_SENSITIVITY')
        self.assertEqual(len(result['comparisons']), 3)
        self.assertTrue(all(set(v['macro_endpoint_screens']) == set(subject.ROLES) for v in result['comparisons']))
        self.assertEqual(result['selected_scale'], .25)

    def test_failed_baseline_is_not_qualified(self):
        f = self.fixture(); f[2][1.][1]['status'] = 'FAIL_GPU_STAGE_NOT_ACCEPTED'
        with self.assertRaisesRegex(ValueError, 'PASS'):
            self.group(f)

    def test_running_baseline_is_not_qualified(self):
        f = self.fixture(); f[2][1.][1]['status'] = 'RUNNING_GPU_DIAGNOSTIC'
        with self.assertRaises(ValueError): self.group(f)

    def test_unauthed_pass_flag_is_rejected(self):
        f = self.fixture(); f[2][1.][1]['screened_continuation_authorized'] = False
        with self.assertRaises(ValueError): self.group(f)

    def test_half_cannot_start_from_tolerance_one_endpoint(self):
        f = self.fixture(); f[2][.5][1]['parent_sha256'] = f[2][1.][1]['screened_continuation_sha256']
        with self.assertRaisesRegex(ValueError, 'complete parent'): self.group(f)

    def test_quarter_cannot_start_from_half_endpoint(self):
        f = self.fixture(); f[2][.25][1]['parent_gpu_directory'] = subject.rel(f[2][.5][0])
        with self.assertRaisesRegex(ValueError, 'complete parent'): self.group(f)

    def test_relocated_checkpoint_alias_cannot_replace_common_parent(self):
        f = self.fixture(); f[2][.5][1]['parent_checkpoint'] = 'fixture/copied_parent.npz'
        with self.assertRaises(ValueError): self.group(f)

    def test_different_dt_cannot_be_compared(self):
        f = self.fixture(); f[2][.25][1]['dt_s'] /= 2
        with self.assertRaisesRegex(ValueError, 'same macro interval'): self.group(f)

    def test_wrong_scale_is_rejected(self):
        f = self.fixture(); f[2][.25][1]['event_tolerance_scale'] = 1.
        with self.assertRaises(ValueError): self.group(f)

    def test_physical_source_cannot_change(self):
        f = self.fixture(); f[2][.5][1]['source_sha256']['new-physics'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'frozen event implementation'): self.group(f)

    def test_parent_implementation_map_cannot_change(self):
        f = self.fixture(); f[2][.5][1]['parent_implementation_source_sha256'] = {}
        with self.assertRaises(ValueError): self.group(f)

    def test_missing_actual_event_is_rejected(self):
        f = self.fixture(); f[2][1.][1]['event_records']['coarse'] = []
        with self.assertRaisesRegex(ValueError, 'actual qualified activation'): self.group(f)

    def test_missing_variant_is_rejected(self):
        f = self.fixture(); del f[2][.5]
        with self.assertRaisesRegex(ValueError, 'all three'): self.group(f)

    def test_particle_inventory_cannot_change(self):
        f = self.fixture(); f[3][f[2][.25][0]/'fine_2/NOT_RESTARTABLE_candidate.npz'][0].particle_id[2] = 99
        with self.assertRaisesRegex(ValueError, 'particle IDs'): self.group(f)

    def test_particle_mass_cannot_change(self):
        f = self.fixture(); f[3][f[2][.5][0]/'coarse/NOT_RESTARTABLE_candidate.npz'][0].m[0] *= 2
        with self.assertRaisesRegex(ValueError, 'masses'): self.group(f)

    def test_crossing_ids_must_match_exactly(self):
        f = self.fixture(); f[3][f[2][.25][0]/'fine_1/NOT_RESTARTABLE_candidate.npz'][1]['bore_observation']['crossing_ids'] = [11, 12]
        with self.assertRaisesRegex(ValueError, 'PID screen'): self.group(f)

    def test_global_position_gate_remains_original(self):
        f = self.fixture(); f[3][f[2][.5][0]/'fine_2/NOT_RESTARTABLE_candidate.npz'][0].x[0] += 6e-6
        with self.assertRaisesRegex(ValueError, 'screen failed'): self.group(f)

    def test_local_stress_gate_catches_single_particle_drift(self):
        f = self.fixture(); f[3][f[2][.25][0]/'fine_2/NOT_RESTARTABLE_candidate.npz'][0].S[0] += 120000.
        with self.assertRaisesRegex(ValueError, 'screen failed'): self.group(f)

    def test_temperature_gate_remains_original(self):
        f = self.fixture(); f[3][f[2][.5][0]/'coarse/NOT_RESTARTABLE_candidate.npz'][0].temperature[1] += .51
        with self.assertRaises(ValueError): self.group(f)

    def test_endpoint_absolute_time_must_match(self):
        f = self.fixture(); f[3][f[2][.5][0]/'fine_1/NOT_RESTARTABLE_candidate.npz'][1]['time_s'] += 1e-6
        with self.assertRaisesRegex(ValueError, 'absolute time'): self.group(f)

    def test_event_incoming_pair_inventory_must_match(self):
        f = self.fixture(); f[2][.5][1]['event_records']['coarse'][0]['incoming_pid_pairs'] = [[11, 13]]
        with self.assertRaisesRegex(ValueError, 'identity'): self.group(f)

    def test_extra_event_cannot_be_hidden(self):
        f = self.fixture(); records = f[2][.25][1]['event_records']['coarse']; records.append(copy.deepcopy(records[0]))
        with self.assertRaisesRegex(ValueError, 'event count'): self.group(f)

    def test_event_migration_between_macro_roles_is_rejected(self):
        f = self.fixture(); meta = f[2][.25][1]; meta['event_records']['fine_1'] = meta['event_records']['coarse']; meta['event_records']['coarse'] = []
        with self.assertRaisesRegex(ValueError, 'event count'): self.group(f)

    def test_event_timing_shift_beyond_registered_bound_fails(self):
        f = self.fixture(); folder = f[2][.25][0]
        for name in ('lo.npz', 'hi.npz'): f[3][folder/name][1]['time_s'] += 1e-5
        with self.assertRaisesRegex(ValueError, 'event timing differs'): self.group(f)

    def test_bracket_must_have_actual_blocked_low_visible_high(self):
        f = self.fixture(); f[4][f[2][.5][0]/'location.json']['lo']['blocked'] = False
        with self.assertRaisesRegex(ValueError, 'blocked-to-visible'): self.group(f)

    def test_bracket_wider_than_verified_tolerance_fails(self):
        f = self.fixture(); f[3][f[2][.5][0]/'hi.npz'][1]['time_s'] += 1e-6
        with self.assertRaisesRegex(ValueError, 'event bracket'): self.group(f)

    def test_cost_envelope_includes_expensive_baseline_probes(self):
        rows = [self.cost(150., 50., 0), self.cost(300., 140.)]
        result = subject.cost_decision(dict(screened_end_time_s=.071), .0004, 1000.,
            [dict(start=0., end=.22, stationary=True, name='fill')], rows)
        self.assertEqual(result['action'], 'STOP_COST_REASSESSMENT')
        self.assertEqual(result['forecast_s'], 480.)

    def test_whole_pair_model_does_not_triple_count_one_long_event_macro(self):
        result = subject.measured_envelope([self.cost(200., 140.)])
        self.assertEqual(result['forecast_s'], 330.)
        self.assertGreater(30.+4.5*140., 450.)
        self.assertFalse(result['guarantees_unknown_events'])

    def test_supervision_cost_is_included_in_whole_pair(self):
        row = self.cost(); row['supervised_wall_s'] = row['measured_pair_wall_s'] = 281.
        result = subject.measured_envelope([row])
        self.assertEqual(result['forecast_s'], 451.5)

    def test_hidden_pair_cost_component_cannot_be_omitted(self):
        row = self.cost(); row['supervised_wall_s'] = 300.
        with self.assertRaisesRegex(ValueError, 'inconsistent'): subject.measured_envelope([row])

    def test_no_event_pair_cannot_qualify_event_forecast(self):
        with self.assertRaisesRegex(ValueError, 'actual secant event pair'):
            subject.measured_envelope([self.cost(events=0)])

    def test_first_qualification_explicitly_marks_event_cost_unmeasured(self):
        result = subject.measured_envelope([self.cost(120., 40., 0)], scope='FIRST_SECANT_QUALIFICATION')
        self.assertFalse(result['event_model_cost_previously_measured'])
        self.assertEqual(result['forecast_s'], 210.)

    def test_qualified_event_cost_cannot_return_to_first_trial_exemption(self):
        with self.assertRaisesRegex(ValueError, 'original no-event'):
            subject.measured_envelope([self.cost()], scope='FIRST_SECANT_QUALIFICATION')

    def test_invalid_physical_branch_cost_cannot_hide_inside_pair_total(self):
        row = self.cost(); row['candidate_wall_s'] = [100., -10., 110.]
        row['maximum_candidate_wall_s'] = 110.
        with self.assertRaisesRegex(ValueError, 'inconsistent'): subject.measured_envelope([row])

    def test_original_phase_and_wall_decisions_are_identical_for_equal_forecast(self):
        phases = [dict(start=0., end=.22, stationary=True, name='fill')]
        for t, dt, remaining in ((.071, .0004, 1000.), (.071, .0004, 350.),
                                  (.2199, .0004, 1000.), (.22, .0004, 1000.), (.21999999, .0004, 1000.)):
            with self.subTest(t=t, remaining=remaining):
                new = subject.cost_decision(dict(screened_end_time_s=t), dt, remaining, phases, [self.cost(210., 70.)])
                old = subject.admin.decision(dict(screened_end_time_s=t, maximum_candidate_wall_s=70.), dt, remaining, phases, 1e-7)
                self.assertEqual({k: v for k, v in new.items() if k != 'event_cost_envelope'}, old)

    def test_expensive_prior_pair_remains_in_cross_parent_forecast(self):
        old = self.cost(275., 100.); fresh = self.cost(150., 60.)
        fresh['parent_sha256'] = 'c'*64; fresh['dt_s'] = .0002
        result = subject.measured_envelope([old, fresh])
        self.assertEqual(result['maximum_complete_pair_wall_s'], 275.)

    def test_cost_evidence_without_probes_is_rejected(self):
        with self.assertRaises(ValueError): subject.measured_envelope([dict(maximum_candidate_wall_s=40., includes_event_probes=False)])

    def test_nonfinite_and_boolean_costs_are_rejected(self):
        for value in (True, float('nan'), float('inf'), 0., -1.):
            with self.subTest(value=value), self.assertRaises(ValueError):
                subject.measured_envelope([dict(maximum_candidate_wall_s=value, includes_event_probes=True)])

    def test_wall_reserve_remains_original(self):
        result = subject.cost_decision(dict(screened_end_time_s=.071), .0004, 330.,
            [dict(start=0., end=.22, stationary=True, name='fill')],
            [self.cost()])
        self.assertEqual(result['action'], 'STOP_WALL_RESERVE')

    def test_new_admin_protocol_preserves_physical_map_and_exact_parent(self):
        parent, common, stages, _, _ = self.fixture()
        template = dict(source_sha256=stages[1.][1]['source_sha256'], v1_source_sha256=common['source_sha256'],
                        secant_components_sha256='1'*64, event_protocol=subject.frozen.PROTOCOL)
        evidence = [self.cost(135., 45.)]
        before = copy.deepcopy(common)
        result = subject.stage_protocol(template, common, parent, .0004, {'admin':'0'*64}, {}, {}, .25, evidence)
        self.assertEqual(common, before)
        self.assertEqual(result['source_sha256'], stages[1.][1]['source_sha256'])
        self.assertTrue(all(result['source_sha256'][k] == v for k, v in before['source_sha256'].items()))
        self.assertEqual(result['parent_sha256'], before['screened_continuation_sha256'])
        self.assertEqual(result['candidate_forecast_s'], 45.)
        self.assertEqual(result['administrative_protocol'], subject.PROTOCOL)

    def test_expensive_variant_cannot_fall_back_to_cheap_parent_forecast(self):
        parent, common, _, _, _ = self.fixture()
        with self.assertRaisesRegex(ValueError, 'forecast limit'):
            subject.stage_protocol(dict(source_sha256=common['source_sha256']), common, parent,
                .0004, {}, {}, {}, .5, [self.cost(300., 100.)])

    def test_siblings_advance_policy_once(self):
        parent, common, stages, saved, _ = self.fixture()
        policy = subject.adaptive.initial_state(common['screened_continuation_sha256'], .0004)
        result, _ = subject.selected_pass(policy, stages[.25][1], saved[subject.ROOT/common['screened_continuation_checkpoint']][1])
        self.assertEqual(result['passes'], 1); self.assertEqual(result['attempts'], 1)
        self.assertEqual(result['parent_sha256'], stages[.25][1]['screened_continuation_sha256'])
        self.assertEqual(policy['passes'], 0)

    def test_absolute_stop_and_ten_hour_bound_are_inherited(self):
        now = subject.frozen.ABSOLUTE_STOP-timedelta(seconds=42)
        self.assertEqual(subject.frozen.authorized_wall_seconds(36000., now), 42.)
        with self.assertRaises(ValueError): subject.frozen.authorized_wall_seconds(36001., now)
        with self.assertRaises(ValueError): subject.frozen.authorized_wall_seconds(1., subject.frozen.ABSOLUTE_STOP)

    def test_one_initial_preflight_and_no_pool_recreation_in_continuation(self):
        tree = ast.parse(inspect.getsource(subject.run))
        calls = [ast.unparse(node.func) for node in ast.walk(tree) if isinstance(node, ast.Call)]
        self.assertEqual(calls.count('frozen._initial_preflight'), 1)
        self.assertEqual(calls.count('frozen.LockedHistory'), 1)
        code = inspect.getsource(subject.continue_selected)
        self.assertNotIn('LockedHistory(', code); self.assertNotIn('_initial_preflight(', code)
        self.assertIn('run_stage(', code); self.assertIn('frozen.validated_stage(', code)

    def test_deposition_only_observes_common_then_selected_quarter_branch(self):
        tree = ast.parse(inspect.getsource(subject.run))
        observations = [ast.unparse(node) for node in ast.walk(tree) if isinstance(node, ast.Call)
                        and ast.unparse(node.func) == 'frozen.observe']
        self.assertEqual(observations, ['frozen.observe(common, frames)', 'frozen.observe(meta, frames)'])
        code = inspect.getsource(subject.run)
        self.assertIn('parent, meta = stages[.25]', code)
        self.assertIn('if len(stages) == 3:', code)

    def test_no_runtime_global_or_provider_rebinding(self):
        tree = ast.parse(Path(subject.__file__).read_text(encoding='utf-8'))
        self.assertFalse(any(isinstance(node, (ast.Global, ast.Nonlocal)) for node in ast.walk(tree)))
        for node in ast.walk(tree):
            if isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for target in targets:
                    self.assertFalse(isinstance(target, ast.Attribute) and
                        ast.unparse(target).split('.')[0] in ('frozen', 'bridge', 'admin', 'adaptive'))

    def test_worker_and_deadlines_reuse_original_live_supervision(self):
        report = subject.static_audit()
        self.assertTrue(report['inherited_event_audit']['live_supervision_unchanged'])
        expected = inspect.getsource(subject.frozen.run_stage).replace('afsd_event_worker_20261010', 'afsd_event_secant_worker_20261010')
        self.assertEqual(ast.dump(ast.parse(expected)), ast.dump(ast.parse(inspect.getsource(subject.run_stage))))

    def test_checker_binds_explicit_new_physics_and_inherits_frozen_receipt(self):
        self.assertEqual(checker.ADMINISTRATIVE_FILES, subject.ADMINISTRATIVE_FILES)
        self.assertEqual(checker.PHYSICAL_FILES, subject.PHYSICAL_FILES)
        self.assertEqual(len(subject.ADMINISTRATIVE_FILES), 5)
        code = inspect.getsource(checker.main)
        self.assertIn('visibility_event_components_20261010_v1.json', code)
        self.assertIn('physical_source_changed=True', code)
        self.assertNotIn('run_gpu_event_campaign_20261010.py', checker.ADMINISTRATIVE_FILES)

    def test_conflicting_frozen_evidence_cannot_be_overwritten(self):
        with self.assertRaises(ValueError): checker.merge({'old':'a'*64}, {'old':'b'*64})

    def test_secant_bridge_only_accepts_registered_complete_v1_parent(self):
        parent, common, stages, _, _ = self.fixture()
        template = dict(source_sha256=stages[1.][1]['source_sha256'], v1_source_sha256=common['source_sha256'],
                        secant_components_sha256='1'*64)
        result = subject.source_transition(template, common, parent)
        self.assertEqual(result['schema'], subject.TRANSITION_SCHEMA)
        self.assertEqual(result['from_source_sha256'], common['source_sha256'])
        with self.assertRaisesRegex(ValueError, 'unregistered'):
            subject.source_transition(template, common, subject.FAILED_V1)

    def test_secant_parent_requires_locator_identity(self):
        parent, common, stages, _, _ = self.fixture(); meta = stages[.25][1]
        template = dict(source_sha256=meta['source_sha256'])
        self.assertIsNone(subject.source_transition(template, meta, stages[.25][0]))
        meta['event_locator_protocol'] = 'unregistered'
        with self.assertRaises(ValueError): subject.source_transition(template, meta, stages[.25][0])

    def test_registered_fourth_source_generation_requires_both_protocols(self):
        _, _, stages, _, _ = self.fixture(); meta = stages[1.][1]
        template = dict(source_sha256=meta['source_sha256'])
        self.assertEqual(subject.source_generation(meta, template), 3)
        meta['event_locator_protocol'] = 'unknown'
        with self.assertRaisesRegex(ValueError, 'fourth implementation'): subject.source_generation(meta, template)

    def test_all_three_variants_are_new_same_parent_solves(self):
        code = inspect.getsource(subject.run)
        self.assertIn("((1., 'one'), (.5, 'half'), (.25, 'quarter'))", code)
        self.assertIn('stages = {}', code)
        self.assertNotIn('tx.gpu_receipt(baseline', code)


if __name__ == '__main__':
    unittest.main()

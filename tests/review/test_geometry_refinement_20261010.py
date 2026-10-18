"""Read-only real failure evidence and fail-closed administrative policy tests."""
import ast
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import afsd_adaptive_dt_20261010 as adaptive
import afsd_geometry_refinement_20261010 as policy

FAILURE = ROOT/'outputs/gpu_acceleration_20261007/incremental_deposition_resume_20261010_v1/attempt_0019'
SAVED_SHA = {
    'manifest.json': '3b50e45af473788b6e547cb362f4d3f0888d41092cdfcd78b97067a5d50ebb6a',
    'supervision.json': '5977142d1bb489289f14d409e9eb1823f2a2b104b06ef0bd1fca0b8011a29f6b',
    'coarse/visible_path_audit.json': 'e4e2ede6a1816516024ef9f13266452abe75892d04ddeeda24a8598bb9e48b99',
    'coarse/continuous_shell_audit.json': '8f770c014124a22696f5e1779fbfe0ce25df40ffd70e871ff7fb97bd2184d5d1',
    'coarse/geometry_audit.json': '3336dce00f1263d7fc2d3c2ac55a0f3a8df644cdcb40cf1fe7254a5ec283310c',
    'coarse/visibility_event_routing.json': 'a9ede028a21f3cc96e92f4f8fa6182914bdd6e6bbe85f3c978dd36851225e2af',
}
FILES = dict(visible_path='visible_path_audit.json', visible_shell='continuous_shell_audit.json',
    geometry='geometry_audit.json', event_router='visibility_event_routing.json')


def load_saved():
    loaded = {}
    for relative, expected in SAVED_SHA.items():
        data = (FAILURE/relative).read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            raise AssertionError('Frozen failure evidence changed: '+relative)
        loaded[relative] = json.loads(data)
    return loaded['manifest.json'], loaded['supervision.json'], {
        key: loaded['coarse/'+name] for key, name in FILES.items()}


def passing_screen(endpoint):
    ids = [7, 11]
    return dict(global_screen=dict(position_max_m=1e-8, temperature_max_K=.001,
        stress_relative=.001, plastic_relative=.001, **{'pass': True}),
        local_screen=dict(adaptive.LOCAL_LIMITS, particle_ids=ids,
            rows=[dict(particle_id=pid, position_error_m=1e-8, temperature_error_K=.001,
                stress_relative=.001, plastic_relative=.001) for pid in ids], passed=True),
        crossing_ids_match=True, passed=True, fine_checkpoint_sha256=['f'*64, endpoint])


class GeometryRecognitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.supervision, cls.reports = load_saved()

    def recognize(self, manifest=None, supervision=None, reports=None):
        return policy.recognize(self.manifest if manifest is None else manifest,
            self.supervision if supervision is None else supervision,
            self.reports if reports is None else reports)

    def test_saved_real_failure_is_budget_only_but_remains_unaccepted(self):
        diagnosis = self.recognize()
        self.assertTrue(diagnosis['eligible'])
        self.assertEqual(diagnosis['unresolved_interval_count'], 4)
        self.assertEqual(diagnosis['original_query_limit'], 256)
        self.assertEqual(diagnosis['parent_sha256'], self.manifest['parent_sha256'])
        self.assertFalse(diagnosis['accepted'])
        self.assertFalse(diagnosis['restart_allowed'])
        self.assertFalse(diagnosis['whole_path_certified'])
        self.assertEqual(diagnosis['accepted_time_advanced_s'], 0.)

    def test_observed_intersection_degeneracy_or_numerical_domain_rejected(self):
        for field in ('found_edges', 'degenerate_edges', 'numerical_domain_failures'):
            with self.subTest(field=field):
                reports = copy.deepcopy(self.reports)
                reports['visible_path'][field] = [{'edge_particle_ids': [1347, 1371]}]
                self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_budget_flag_alone_cannot_hide_a_different_interval_failure(self):
        for reason in ('MAXIMUM_SUBDIVISION_DEPTH', 'STRICT_SOLID_INTERSECTION', None):
            with self.subTest(reason=reason):
                reports = copy.deepcopy(self.reports)
                reports['visible_path']['unresolved_intervals'][-1]['reason'] = reason
                self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_original_query_limits_and_complete_accounting_required(self):
        for field, value in [('maximum_queries', 512), ('analytic_batch_queries', 255),
                ('maximum_depth', 64), ('edge_deletions', 1), ('physical_solver_evaluations', 1),
                ('enclosure_batch_queries', 128), ('query_log', []), ('unresolved_intervals', [])]:
            with self.subTest(field=field):
                reports = copy.deepcopy(self.reports); reports['visible_path'][field] = value
                self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_malformed_unresolved_intervals_fail_closed(self):
        for field, value in [('start_fraction', float('nan')), ('end_fraction', 2.),
                ('depth', True), ('edge_indices', []), ('edge_indices', [999999]),
                ('edge_particle_ids', [[1347, 1347]])]:
            with self.subTest(field=field, value=value):
                reports = copy.deepcopy(self.reports)
                reports['visible_path']['unresolved_intervals'][0][field] = value
                self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_shell_certificate_required_even_when_visibility_only_exhausted(self):
        for field, value in [('status', 'UNRESOLVED'), ('failure_kind', 'MAXIMUM_GEOMETRY_QUERIES'),
                ('hidden_contacts', [{}]), ('unresolved_intervals', [{}]),
                ('finite_patch_unresolved_intervals', [{}])]:
            with self.subTest(field=field):
                reports = copy.deepcopy(self.reports); reports['visible_shell'][field] = value
                self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_new_top_and_shoulder_remain_original_allowed_features(self):
        for features in ([], ['TOP'], ['SHOULDER'], ['SHOULDER', 'TOP']):
            reports = copy.deepcopy(self.reports)
            reports['geometry']['new_non_bore_features'] = features
            self.assertTrue(self.recognize(reports=reports)['eligible'])
        for field, value in [('new_non_bore_features', ['BASE']),
                ('transient_contacts_missed_by_endpoints', ['TOP'])]:
            reports = copy.deepcopy(self.reports); reports['geometry'][field] = value
            self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_actual_or_unresolved_router_events_rejected(self):
        for field, value in [('requires_resolve', True), ('changed_edges', [{}]),
                ('endpoint_changed_edges', [{}]), ('transition_brackets', [{}]),
                ('unresolved_samples', [{}]), ('proposed_intervals', [{}]),
                ('start_blocked_edge_count', 1), ('edge_deletions', 1),
                ('earliest_observed_bracket', {}), ('acceptance_allowed', True)]:
            with self.subTest(field=field):
                reports = copy.deepcopy(self.reports); reports['event_router'][field] = value
                self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_router_query_and_geometry_inputs_must_match_failed_path(self):
        for key in ('world_position0', 'world_position1', 'particle_id_order',
                'supplied_topology', 'supplied_edge_pid_order', 'rotation1'):
            reports = copy.deepcopy(self.reports)
            reports['event_router']['input_sha256'][key] = '0'*64
            self.assertFalse(self.recognize(reports=reports)['eligible'])
        reports = copy.deepcopy(self.reports)
        reports['event_router']['query_log'][-1]['blocked_edge_count'] = 1
        self.assertFalse(self.recognize(reports=reports)['eligible'])

    def test_exact_geometry_exception_and_failure_flags_required(self):
        for field, value in [('status', 'PASS_GPU_STAGE'), ('error', adaptive.TEMPORAL_ERROR),
                ('accepted', True), ('restart_allowed', True), ('diagnostic_only', False),
                ('accepted_time_advanced_s', .1), ('dt_s', float('nan'))]:
            with self.subTest(field=field):
                manifest = copy.deepcopy(self.manifest); manifest[field] = value
                self.assertFalse(self.recognize(manifest=manifest)['eligible'])

    def test_interruption_source_failure_and_missing_end_grant_rejected(self):
        for field, value in [('exit_code', 0), ('source_inputs_still_match', False),
                ('source_error', 'changed'), ('interrupted', True), ('error', 'TimeoutError'),
                ('live_history_grants', self.supervision['live_history_grants'][:1])]:
            with self.subTest(field=field):
                supervision = copy.deepcopy(self.supervision); supervision[field] = value
                self.assertFalse(self.recognize(supervision=supervision)['eligible'])
        for field, value in [('phase', 'start'), ('protocol_sha256', '0'*64),
                ('controller_pid', 1), ('history_locked', False)]:
            supervision = copy.deepcopy(self.supervision)
            supervision['live_history_grants'][1][field] = value
            manifest = copy.deepcopy(self.manifest)
            manifest['live_history_grants'] = supervision['live_history_grants']
            self.assertFalse(self.recognize(manifest=manifest, supervision=supervision)['eligible'])

    def test_missing_and_malformed_reports_fail_closed(self):
        for key in policy.REPORT_NAMES:
            reports = copy.deepcopy(self.reports); reports.pop(key)
            self.assertFalse(self.recognize(reports=reports)['eligible'])
        for malformed in ({}, [], None, True):
            reports = copy.deepcopy(self.reports); reports['visible_path'] = malformed
            self.assertFalse(self.recognize(reports=reports)['eligible'])


class GeometryRefinementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.supervision, cls.reports = load_saved()
        cls.diagnosis = policy.recognize(cls.manifest, cls.supervision, cls.reports)

    def state(self, dt=None, **kwargs):
        dt = self.manifest['dt_s'] if dt is None else dt
        return adaptive.initial_state(self.manifest['parent_sha256'], dt, **kwargs)

    def refine(self, state, **changes):
        pending = adaptive.begin_attempt(state, state['nominal_dt_s'])
        diagnosis = dict(self.diagnosis, attempted_dt_s=pending['pending']['dt_s'], **changes)
        return policy.refine(pending, diagnosis)

    def test_real_rejection_halves_exact_parent_without_advancement(self):
        pending = adaptive.begin_attempt(self.state(), self.manifest['dt_s'])
        original = copy.deepcopy(pending)
        result, event = policy.refine(pending, self.diagnosis)
        self.assertEqual(pending, original)
        self.assertEqual(event['action'], 'REFINE')
        self.assertEqual(result['nominal_dt_s'], self.manifest['dt_s']/2)
        self.assertEqual(result['parent_sha256'], pending['parent_sha256'])
        self.assertEqual(result['passes'], 0)
        self.assertEqual(result['attempts'], 1)
        self.assertEqual(result['last_passed_dt_s'], None)
        self.assertEqual(result['small_error_streak'], 0)
        self.assertEqual(result['cooldown_remaining'], 8)
        self.assertEqual(result['consecutive_halvings'], 1)
        self.assertFalse(event['accepted'])
        self.assertEqual(event['accepted_time_advanced_s'], 0.)
        adaptive.validate_state(result)

    def test_diagnosis_for_another_parent_or_dt_never_authorizes_retry(self):
        pending = adaptive.begin_attempt(self.state(), self.manifest['dt_s'])
        for key, value in [('parent_sha256', '0'*64), ('attempted_dt_s', self.manifest['dt_s']/2),
                ('eligible', False), ('protocol', 'OTHER'), ('whole_path_certified', True),
                ('accepted', True), ('restart_allowed', True)]:
            diagnosis = dict(self.diagnosis); diagnosis[key] = value
            result, event = policy.refine(pending, diagnosis)
            self.assertEqual(event['action'], 'STOP')
            self.assertEqual(result['parent_sha256'], pending['parent_sha256'])
            self.assertEqual(result['passes'], 0)

    def test_original_eight_halvings_limit_is_shared_with_temporal_policy(self):
        state = self.state()
        for count in range(1, 9):
            state, event = self.refine(state)
            self.assertEqual(event['action'], 'REFINE')
            self.assertEqual(state['consecutive_halvings'], count)
            self.assertEqual(state['passes'], 0)
        stopped, event = self.refine(state)
        self.assertEqual(event['reason'], 'GEOMETRY_REFINEMENT_LIMIT')
        self.assertEqual(stopped['consecutive_halvings'], 8)
        self.assertEqual(stopped['parent_sha256'], self.manifest['parent_sha256'])
        state = self.state(max_halvings=3); state['consecutive_halvings'] = 3
        self.assertEqual(self.refine(state)[1]['action'], 'STOP')

    def test_fixed_minimum_cannot_be_bypassed_by_adaptive_configuration(self):
        state, event = self.refine(self.state(dt=2e-7))
        self.assertEqual(state['nominal_dt_s'], 1e-7)
        self.assertEqual(event['action'], 'REFINE')
        self.assertEqual(self.refine(state)[1]['reason'], 'MINIMUM_DT_REACHED')
        state = self.state(dt=1.9e-7, minimum_dt_s=1e-9)
        self.assertEqual(self.refine(state)[1]['reason'], 'MINIMUM_DT_REACHED')

    def test_phase_clipped_failure_halves_actual_interval(self):
        state = self.state()
        pending = adaptive.begin_attempt(state, state['nominal_dt_s']/2, phase_clipped=True)
        diagnosis = dict(self.diagnosis, attempted_dt_s=pending['pending']['dt_s'])
        result, event = policy.refine(pending, diagnosis)
        self.assertEqual(result['nominal_dt_s'], state['nominal_dt_s']/4)
        self.assertTrue(event['phase_clipped'])

    def test_growth_failure_resets_streak_and_never_shortens_existing_cooldown(self):
        state = self.state(); state.update(small_error_streak=1, growth_pending=True, cooldown_remaining=16)
        result, event = self.refine(state)
        self.assertEqual(result['growth_failures'], 1)
        self.assertEqual(result['small_error_streak'], 0)
        self.assertEqual(result['cooldown_remaining'], 16)
        self.assertFalse(result['growth_pending'])

    def test_subsequent_original_complete_passes_advance_and_honor_cooldown(self):
        state, event = self.refine(self.state())
        reduced = state['nominal_dt_s']
        for count in range(1, 11):
            endpoint = format(count, '064x')
            pending = adaptive.begin_attempt(state, state['nominal_dt_s'])
            state, event = adaptive.resolve_pass(pending, passing_screen(endpoint), .0005,
                [7, 11], endpoint, original_candidate_gates_passed=True)
            self.assertEqual(state['parent_sha256'], endpoint)
            self.assertEqual(state['passes'], count)
            self.assertEqual(state['consecutive_halvings'], 0)
            if count <= 8:
                self.assertEqual(state['cooldown_remaining'], 8-count)
                self.assertEqual(state['small_error_streak'], 0)
            self.assertEqual(state['nominal_dt_s'], reduced if count < 10 else reduced*2)
        self.assertEqual(event['action'], 'GROW')

    def test_no_solver_gpu_filesystem_or_new_physics_imports(self):
        tree = ast.parse((ROOT/'scripts/afsd_geometry_refinement_20261010.py').read_text())
        imports = [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names]
        imports += [node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
        self.assertEqual(set(imports), {'copy', 'afsd_adaptive_dt_20261010'})
        self.assertFalse(any(isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id in {'open', 'eval', 'exec', '__import__'} for node in ast.walk(tree)))


if __name__ == '__main__':
    unittest.main()

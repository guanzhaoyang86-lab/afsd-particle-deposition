"""Bound real multi-opening failure plus fail-closed administrative mutations."""
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
import afsd_multi_event_refinement_20261011 as policy

FAILURE = ROOT/'outputs/gpu_acceleration_20261007/event_resume_90m_20261010_v1/attempt_0014'
FILES = dict(visible_path='visible_path_audit.json', visible_shell='continuous_shell_audit.json',
    geometry='geometry_audit.json', event_router='visibility_event_routing.json',
    union_event_router='proposal/union_visibility_events.json')
SAVED_SHA = {
    'manifest.json': 'c7524386acbcd596f36ca857589fe7ce0399646267296eece153ca12ea4d57f3',
    'supervision.json': 'd299ac6f974a6418f69f22d6c2df8c3addd35ccf68604fdea7f1016594af7824',
    'coarse/visible_path_audit.json': 'b4b163d94b35255ec9867dceaa576d27a74ea051aa68c12cce1dba7250c7a860',
    'coarse/continuous_shell_audit.json': 'd3fcbf5277a13bf6f92c9df5d3092e4b9925d8ed2faf96ea42302158133bbf7a',
    'coarse/geometry_audit.json': '6115b361046687cbf96c1c5ccc11eb661ed00ae94b103800f8035df4a659fec2',
    'coarse/visibility_event_routing.json': '65da66973a7ed51b4ee587d3373e17cbd6eb831b73e285824bb2bac4fe393241',
    'coarse/proposal/union_visibility_events.json': '4469724701036f5abf159ba93972007e3da0fddc7a5744e7c52a40dda79f1756',
}


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


class RecognitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.supervision, cls.reports = load_saved()

    def recognize(self, manifest=None, supervision=None, reports=None):
        return policy.recognize(self.manifest if manifest is None else manifest,
            self.supervision if supervision is None else supervision,
            self.reports if reports is None else reports)

    def mutate(self, *path, value):
        reports = copy.deepcopy(self.reports)
        target = reports
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        result = self.recognize(reports=reports)
        self.assertFalse(result['eligible'], (path, value, result))

    def test_real_failure_matches_two_diagnostic_openings_and_remains_fail(self):
        original = copy.deepcopy((self.manifest, self.supervision, self.reports))
        result = self.recognize()
        self.assertTrue(result['eligible'], result)
        self.assertEqual(result['observed_opening_event_count'], 2)
        self.assertEqual(result['opening_particle_id_pairs'], [[1059, 1087], [1111, 1243]])
        self.assertEqual(result['parent_sha256'], self.manifest['parent_sha256'])
        self.assertEqual(result['attempted_dt_s'], .000474609375)
        self.assertEqual(result['accepted_time_advanced_s'], 0.)
        for key in ('accepted', 'restart_allowed', 'whole_path_certified',
                'event_uniqueness_certified', 'physical_substeps_solved'):
            self.assertIs(result[key], False)
        self.assertEqual((self.manifest, self.supervision, self.reports), original)

    def test_exact_failed_geometry_manifest_required(self):
        for key, value in [('status', 'PASS'), ('error', 'OTHER'), ('accepted', True),
                ('restart_allowed', True), ('diagnostic_only', False),
                ('screened_continuation_authorized', True), ('accepted_time_advanced_s', 1e-9),
                ('accepted_time_advanced_s', False), ('dt_s', float('nan')), ('parent_sha256', 'x')]:
            with self.subTest(key=key, value=value):
                manifest = dict(self.manifest, **{key: value})
                self.assertFalse(self.recognize(manifest=manifest)['eligible'])

    def test_supervision_and_both_live_grants_required(self):
        for key, value in [('exit_code', True), ('exit_code', 0), ('interrupted', True),
                ('source_inputs_still_match', False), ('source_error', 'changed'),
                ('error', 'TimeoutError'), ('manifest_sha256', 'bad'),
                ('live_history_grants', self.supervision['live_history_grants'][:1])]:
            with self.subTest(key=key):
                self.assertFalse(self.recognize(supervision=dict(self.supervision, **{key: value}))['eligible'])
        for key, value in [('phase', 'start'), ('controller_pid', 1),
                ('protocol_sha256', '0'*64), ('history_locked', False)]:
            supervision = copy.deepcopy(self.supervision)
            supervision['live_history_grants'][1][key] = value
            manifest = dict(self.manifest, live_history_grants=supervision['live_history_grants'])
            self.assertFalse(self.recognize(manifest=manifest, supervision=supervision)['eligible'])

    def test_five_original_reports_mandatory(self):
        for key in policy.REPORT_NAMES:
            reports = dict(self.reports); reports.pop(key)
            self.assertFalse(self.recognize(reports=reports)['eligible'])
            self.mutate(key, value=None)
        self.assertFalse(self.recognize(reports=dict(self.reports, extra={}))['eligible'])

    def test_actual_intersection_at_any_nonzero_fraction_rejected(self):
        for fraction in (1e-12, .5, 1., False, float('nan')):
            self.mutate('visible_path', 'found_edges', 0, 'fraction', value=fraction)
        self.mutate('visible_path', 'found_edges', 0, 'total_inside_length_m', value=0.)
        self.mutate('visible_path', 'found_edges', 0, 'strict_intervals', value=[[.5, .4]])

    def test_domain_degeneracy_unresolved_budget_or_modified_path_rejected(self):
        for key in ('degenerate_edges', 'numerical_domain_failures', 'unresolved_intervals'):
            self.mutate('visible_path', key, value=[{}])
        for key, value in [('status', 'UNRESOLVED_VISIBLE_SUPPORT_PATH'),
                ('failure_kind', 'MAXIMUM_ANALYTIC_BATCH_QUERIES'), ('maximum_queries', 512),
                ('maximum_depth', 64), ('physical_solver_evaluations', 1),
                ('edge_deletions', 1), ('prism_batch_queries', 1), ('evaluated_fractions', [0., .5, 1.]),
                ('evaluated_fractions', [False, True])]:
            self.mutate('visible_path', key, value=value)

    def test_single_or_duplicate_pair_is_not_multi_opening_scope(self):
        self.mutate('visible_path', 'found_edges', value=self.reports['visible_path']['found_edges'][:1])
        self.mutate('visible_path', 'found_edges', value=[self.reports['visible_path']['found_edges'][0]]*2)
        for key, value in [('particle_id_i', True), ('i', -1), ('j', 99999), ('edge_index', 99999)]:
            self.mutate('visible_path', 'found_edges', 0, key, value=value)

    def test_closing_unknown_and_nonlocalized_events_rejected_in_either_graph(self):
        for graph in ('event_router', 'union_event_router'):
            for key, value in [('direction', 'VISIBLE_TO_BLOCKED'), ('start_blocked', False),
                    ('end_blocked', True), ('status', 'MAXIMUM_EVALUATIONS'),
                    ('single_event_certified', True), ('physical_substeps_solved', True),
                    ('end_fraction', float('inf')), ('bracket_width_fraction', 1e-5)]:
                self.mutate(graph, 'transition_brackets', 0, key, value=value)
            self.mutate(graph, 'unresolved_samples', value=[{}])
            self.mutate(graph, 'changed_edges', 0, 'particle_id_i', value=99999)
            self.mutate(graph, 'endpoint_changed_edges', value=self.reports[graph]['endpoint_changed_edges'][:1])
            self.mutate(graph, 'transition_brackets', value=self.reports[graph]['transition_brackets']*2)

    def test_union_exposes_extra_outgoing_pair_missing_from_active_graph(self):
        reports = copy.deepcopy(self.reports)
        router = reports['union_event_router']
        extra = dict(router['endpoint_changed_edges'][0], particle_id_i=9001,
            particle_id_j=9002, edge_index=20, i=20, j=21, start_blocked=False,
            end_blocked=True, start_strict_intersection_length_m=0., end_strict_intersection_length_m=1e-6)
        router['endpoint_changed_edges'].append(extra)
        self.assertFalse(self.recognize(reports=reports)['eligible'])
        self.assertTrue(self.recognize()['eligible'])

    def test_shell_hidden_centre_and_new_geometry_failures_rejected(self):
        for key in ('hidden_contacts', 'unresolved_intervals', 'finite_patch_unresolved_intervals'):
            self.mutate('visible_shell', key, value=[{}])
        self.mutate('visible_shell', 'solid_path_valid', value=False)
        self.mutate('visible_shell', 'solid_path_valid', value=1)
        self.mutate('visible_shell', 'status', value='UNRESOLVED')
        self.mutate('geometry', 'transient_contacts_missed_by_endpoints', value=[{}])
        self.mutate('geometry', 'new_non_bore_features', value=['OUTER_CYLINDER'])

    def test_input_hashes_and_active_topology_bind_all_reports(self):
        for graph in ('event_router', 'union_event_router'):
            for key in ('world_position0', 'world_position1', 'particle_id_order',
                    'origin0', 'origin1', 'rotation0', 'rotation1'):
                self.mutate(graph, 'input_sha256', key, value='0'*64)
        for key in ('supplied_topology', 'supplied_edge_pid_order'):
            self.mutate('event_router', 'input_sha256', key, value='0'*64)
            self.mutate('union_event_router', 'input_sha256', key, value='invalid')
        # Union naturally has additional edges; requiring topology equality is incorrect.
        self.assertNotEqual(self.reports['event_router']['input_sha256']['supplied_topology'],
            self.reports['union_event_router']['input_sha256']['supplied_topology'])

    def test_router_counts_samples_effects_and_original_limits_immutable(self):
        for graph in ('event_router', 'union_event_router'):
            for key, value in [('requires_resolve', False), ('acceptance_allowed', True),
                    ('maximum_evaluations', 512), ('sample_intervals', 8),
                    ('fraction_tolerance', 1e-5), ('start_blocked_edge_count', 3),
                    ('end_blocked_edge_count', 1), ('geometry_evaluations', 70),
                    ('force_evaluations', 1), ('particle_count', 100),
                    ('earliest_observed_bracket', None), ('proposed_intervals', [])]:
                self.mutate(graph, key, value=value)
            self.mutate(graph, 'query_log', 3, 'blocked_edge_count', value=0)
            self.mutate(graph, 'query_log', 0, 'fraction', value=False)


class RefinementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.supervision, cls.reports = load_saved()
        cls.diagnosis = policy.recognize(cls.manifest, cls.supervision, cls.reports)

    def state(self, dt=None, **kwargs):
        return adaptive.initial_state(self.manifest['parent_sha256'],
            self.manifest['dt_s'] if dt is None else dt, **kwargs)

    def refine(self, state):
        pending = adaptive.begin_attempt(state, state['nominal_dt_s'])
        return policy.refine(pending, dict(self.diagnosis, attempted_dt_s=pending['pending']['dt_s']))

    def test_real_rejection_halves_same_parent_without_advancing_or_mutating(self):
        pending = adaptive.begin_attempt(self.state(), self.manifest['dt_s'])
        original = copy.deepcopy((pending, self.diagnosis))
        result, event = policy.refine(pending, self.diagnosis)
        self.assertEqual((pending, self.diagnosis), original)
        self.assertEqual(event['action'], 'REFINE')
        self.assertEqual(result['nominal_dt_s'], .0002373046875)
        self.assertEqual(result['parent_sha256'], pending['parent_sha256'])
        self.assertEqual(result['passes'], 0)
        self.assertEqual(result['attempts'], 1)
        self.assertIsNone(result['last_passed_dt_s'])
        self.assertIsNone(result['pending'])
        self.assertEqual(result['cooldown_remaining'], 8)
        self.assertFalse(event['accepted'])
        self.assertEqual(event['accepted_time_advanced_s'], 0.)
        adaptive.validate_state(result)

    def test_mismatched_parent_dt_or_invalid_diagnosis_stops(self):
        pending = adaptive.begin_attempt(self.state(), self.manifest['dt_s'])
        for key, value in [('parent_sha256', '0'*64), ('attempted_dt_s', 1e-6),
                ('eligible', False), ('protocol', 'OTHER'), ('whole_path_certified', True),
                ('accepted', True), ('restart_allowed', True), ('failed_candidate_retained', False),
                ('event_uniqueness_certified', True), ('physical_substeps_solved', True),
                ('accepted_time_advanced_s', 1e-6), ('accepted_time_advanced_s', False)]:
            state, event = policy.refine(pending, dict(self.diagnosis, **{key: value}))
            self.assertEqual(event['action'], 'STOP')
            self.assertEqual(state['parent_sha256'], pending['parent_sha256'])
            self.assertEqual(state['passes'], 0)

    def test_existing_halving_counter_and_original_eight_limit_shared(self):
        state = self.state()
        for count in range(1, 9):
            state, event = self.refine(state)
            self.assertEqual(event['action'], 'REFINE')
            self.assertEqual(state['consecutive_halvings'], count)
            self.assertEqual(state['passes'], 0)
        stopped, event = self.refine(state)
        self.assertEqual(event['reason'], 'MULTI_EVENT_REFINEMENT_LIMIT')
        self.assertEqual(stopped['consecutive_halvings'], 8)
        state = self.state(max_halvings=3); state['consecutive_halvings'] = 3
        self.assertEqual(self.refine(state)[1]['action'], 'STOP')

    def test_fixed_minimum_dt_cannot_be_relaxed_by_config(self):
        state, event = self.refine(self.state(dt=2e-7))
        self.assertEqual(event['action'], 'REFINE')
        self.assertEqual(state['nominal_dt_s'], 1e-7)
        self.assertEqual(self.refine(state)[1]['reason'], 'MINIMUM_DT_REACHED')
        self.assertEqual(self.refine(self.state(dt=1.9e-7, minimum_dt_s=1e-9))[1]['reason'],
            'MINIMUM_DT_REACHED')

    def test_actual_phase_clipped_interval_is_halved(self):
        state = self.state()
        pending = adaptive.begin_attempt(state, state['nominal_dt_s']/2, phase_clipped=True)
        result, event = policy.refine(pending, dict(self.diagnosis, attempted_dt_s=pending['pending']['dt_s']))
        self.assertEqual(result['nominal_dt_s'], state['nominal_dt_s']/4)
        self.assertTrue(event['phase_clipped'])

    def test_small_streak_reset_growth_accounted_existing_cooldown_preserved(self):
        state = self.state(); state.update(small_error_streak=1, growth_pending=True, cooldown_remaining=16)
        result, event = self.refine(state)
        self.assertEqual(result['growth_failures'], 1)
        self.assertEqual(result['small_error_streak'], 0)
        self.assertEqual(result['cooldown_remaining'], 16)
        self.assertFalse(result['growth_pending'])

    def test_frozen_adaptive_complete_pass_still_required_and_compatible(self):
        state, event = self.refine(self.state())
        smaller = state['nominal_dt_s']
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
            self.assertEqual(state['nominal_dt_s'], smaller if count < 10 else 2*smaller)
        self.assertEqual(event['action'], 'GROW')

    def test_no_physics_gpu_or_filesystem_in_policy(self):
        tree = ast.parse((ROOT/'scripts/afsd_multi_event_refinement_20261011.py').read_text())
        imports = [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names]
        imports += [node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
        self.assertEqual(set(imports), {'copy', 'afsd_adaptive_dt_20261010'})
        self.assertFalse(any(isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id in {'open', 'eval', 'exec', '__import__'} for node in ast.walk(tree)))


if __name__ == '__main__':
    unittest.main()

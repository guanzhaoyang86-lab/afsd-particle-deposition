"""Real bound failure evidence and narrowly scoped nonlinear retry tests."""
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
import afsd_nonlinear_refinement_20261011 as policy

FAILURE = ROOT/'outputs/gpu_acceleration_20261007/event_resume_8h_20261011_v1/attempt_0040'
FILES = dict(secant_failure='broyden/secant_failure.json', broyden_iterations='broyden/iterations.json',
    broyden_candidate='broyden/last_UNCOMMITTED_candidate.json',
    broyden_graph='broyden/last_UNCOMMITTED_candidate_graph.json',
    newton_iterations='newton_retry/iterations.json', newton_candidate='newton_retry/last_UNCOMMITTED_candidate.json')
SAVED_SHA = {
    'manifest.json': '1e675e7c8c114e5381d732a59d9036369085905cba1259f814183d8069ab3fe2',
    'supervision.json': '8d5095b0218e0ab016a54e4ae70f6594162dceae4082f1fcd209ce643e5e5dd6',
    'coarse/proposal/broyden/secant_failure.json': 'fda350bce9c0438a4a2885cc46dd1863a4191380e6080f7598c57e8979eceb19',
    'coarse/proposal/broyden/iterations.json': 'b7dcce7fd31075224ee717ba8dc2e2f6a2a12a1e886cea4c9073d965022a1f10',
    'coarse/proposal/broyden/last_UNCOMMITTED_candidate.json': '64dad9121e850a8eeb3764d572d97e6de36490a7f52be525d2f8192752705ad6',
    'coarse/proposal/broyden/last_UNCOMMITTED_candidate_graph.json': '1ed587b50a0ad0662f33589a21cf7122f19d1f9c06a60e9d3ff0416cedb9504c',
    'coarse/proposal/newton_retry/iterations.json': '6934f66f8ad937b57d0c650b98749ed393b5b6d4d3ec2cc3cd3582156af76cac',
    'coarse/proposal/newton_retry/last_UNCOMMITTED_candidate.json': 'a6a376f44c7f5d710f0a5d3fe350d876a0edbe9cb41b3333d662f9ed02101259',
}


def load_saved():
    loaded = {}
    for relative, expected in SAVED_SHA.items():
        data = (FAILURE/relative).read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            raise AssertionError('Frozen failure evidence changed: '+relative)
        loaded[relative] = json.loads(data)
    return loaded['manifest.json'], loaded['supervision.json'], {
        key: loaded['coarse/proposal/'+name] for key, name in FILES.items()}


def changed(tree, path, value):
    """Copy only the mutated branch; never alter the original large journals."""
    result = copy.copy(tree)
    if len(path) == 1:
        result[path[0]] = value
    else:
        result[path[0]] = changed(tree[path[0]], path[1:], value)
    return result


class RecognitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.supervision, cls.reports = load_saved()

    def recognize(self, manifest=None, supervision=None, reports=None):
        return policy.recognize(self.manifest if manifest is None else manifest,
            self.supervision if supervision is None else supervision,
            self.reports if reports is None else reports)

    def reject(self, path, value):
        result = self.recognize(reports=changed(self.reports, path, value))
        self.assertFalse(result['eligible'], (path, value, result))

    def test_real_old40_recognized_but_never_accepted_or_certified(self):
        before = hashlib.sha256(json.dumps(self.reports, sort_keys=True).encode()).hexdigest()
        result = self.recognize()
        self.assertTrue(result['eligible'], result)
        self.assertEqual(result['broyden_evaluations'], 48)
        self.assertEqual(result['newton_evaluations'], 27)
        self.assertEqual(result['scope'], 'FIRST_COARSE_PROPOSAL_ONLY')
        self.assertEqual(result['parent_sha256'], self.manifest['parent_sha256'])
        self.assertEqual(result['attempted_dt_s'], .000474609375)
        self.assertEqual(result['accepted_time_advanced_s'], 0.)
        for key in ('accepted', 'restart_allowed', 'whole_path_certified', 'nonlinear_converged'):
            self.assertIs(result[key], False)
        self.assertEqual(hashlib.sha256(json.dumps(self.reports, sort_keys=True).encode()).hexdigest(), before)

    def test_broad_value_errors_timeouts_and_physical_failures_do_not_match(self):
        for error in ["ValueError('complete residual line search/evaluation budget failed')",
                "RefinableCandidateError('registered iteration/time budget exhausted')",
                "RefinableCandidateError('registered density domain exceeded')",
                "GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')"]:
            self.assertFalse(self.recognize(manifest=dict(self.manifest, error=error))['eligible'])
        for key, value in [('accepted', True), ('restart_allowed', True), ('status', 'PASS'),
                ('diagnostic_only', False), ('accepted_time_advanced_s', 1e-9),
                ('accepted_time_advanced_s', False), ('dt_s', float('nan'))]:
            self.assertFalse(self.recognize(manifest=dict(self.manifest, **{key: value}))['eligible'])

    def test_supervision_source_and_live_grants_are_required(self):
        for key, value in [('exit_code', 0), ('interrupted', True), ('source_inputs_still_match', False),
                ('source_error', 'changed'), ('termination_reason', 'timeout'),
                ('live_history_grants', self.supervision['live_history_grants'][:1])]:
            self.assertFalse(self.recognize(supervision=dict(self.supervision, **{key: value}))['eligible'])
        supervision = copy.deepcopy(self.supervision)
        supervision['live_history_grants'][1]['controller_pid'] += 1
        manifest = dict(self.manifest, live_history_grants=supervision['live_history_grants'])
        self.assertFalse(self.recognize(manifest=manifest, supervision=supervision)['eligible'])

    def test_all_six_reports_mandatory_and_fail_closed(self):
        for name in policy.REPORT_NAMES:
            reports = dict(self.reports); reports.pop(name)
            self.assertFalse(self.recognize(reports=reports)['eligible'])
            self.reject((name,), None)
        self.assertFalse(self.recognize(reports=dict(self.reports, candidate_qualification={}))['eligible'])

    def test_broyden_exact_error_refresh_and_iteration_counts_required(self):
        for key, value in [('error', "ValueError('nonpositive BE density denominator')"),
                ('jacobian_refreshes', 1), ('jacobian_builds', 4), ('evaluations', 47),
                ('accepted', True), ('controller', 'OTHER')]:
            self.reject(('secant_failure', key), value)
        self.reject(('broyden_iterations', -1, 'terminal_error'), 'unsafe secant matrix')
        self.reject(('broyden_iterations', 0, 'jacobian_builds', 0, 'reason'), 'LINEAR_FAILURE')
        self.reject(('broyden_iterations', 2, 'secant', 'status'), 'REFRESH_REQUIRED_UNSAFE_UPDATE')

    def test_any_domain_rejection_or_nonfinite_trial_prevents_retry(self):
        for name in ('broyden_iterations', 'newton_iterations'):
            row = dict(self.reports[name][-1]['line_search'][0], rejected='contact geometry invalid')
            self.reject((name, -1, 'line_search', 0), row)
            self.reject((name, -1, 'line_search', 0, 'merit_N'), float('nan'))
            self.reject((name, 0, 'density_residual'), 1.01e-10)
            self.reject((name, 0, 'thermal_residual_K'), 1.01e-8)
            self.reject((name, 0, 'virtual_work_residual_J'), 1.)
            self.reject((name, 0, 'iteration'), False)

    def test_complete_original_backtracks_and_actual_armijo_failure_required(self):
        for name in ('broyden_iterations', 'newton_iterations'):
            self.reject((name, -1, 'line_search'), self.reports[name][-1]['line_search'][:7])
            self.reject((name, -1, 'line_search', 0, 'alpha'), .9)
            self.reject((name, -1, 'line_search', 0, 'merit_N'), 0.)
            self.reject((name, 1, 'evaluations'), 99)
        self.reject(('broyden_iterations', -1, 'accepted_alpha'), 1.)

    def test_stale_candidate_wrong_dt_force_limit_and_mixed_initial_state_rejected(self):
        for name in ('broyden_candidate', 'newton_candidate'):
            self.reject((name, 'accepted'), True)
            self.reject((name, 'restart_allowed'), True)
            self.reject((name, 'contact', 'dt_s'), self.manifest['dt_s']/2)
            self.reject((name, 'metrics', 'force_limit_N'), .001)
            self.reject((name, 'metrics', 'scaled_residual'), .5)
            self.reject((name, 'contact', 'sweep', 'certified'), False)
            self.reject((name, 'contact', 'force_N', 0, 0), float('inf'))
        self.reject(('newton_iterations', 0, 'material_endpoint_work_J'), .02)

    def test_energy_and_graph_domain_failures_not_diagnosed_as_force_only(self):
        for name in ('broyden_candidate', 'newton_candidate'):
            self.reject((name, 'energy', 'registered_energy_limit_J'), .01)
            self.reject((name, 'energy', 'global_raw_residual_J'), 1.)
        self.reject(('broyden_graph', 'graph_fingerprint'), '0'*64)
        self.reject(('broyden_graph', 'diagnostics', 'blocked_edge_count'), 1)
        self.reject(('broyden_graph', 'diagnostics', 'physical_solid_intersection_length_tolerance_m'), 1e-9)
        self.reject(('broyden_graph', 'diagnostics', 'canonical_orientation_preserved'), False)


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

    def test_real40_returns_same_complete_parent_half_step_without_pass(self):
        pending = adaptive.begin_attempt(self.state(), self.manifest['dt_s'])
        original = copy.deepcopy(pending)
        state, event = policy.refine(pending, self.diagnosis)
        self.assertEqual(pending, original)
        self.assertEqual(event['action'], 'REFINE')
        self.assertEqual(state['status'], 'READY')
        self.assertEqual(state['parent_sha256'], pending['parent_sha256'])
        self.assertEqual(state['nominal_dt_s'], .0002373046875)
        self.assertEqual(state['passes'], 0)
        self.assertEqual(state['attempts'], 1)
        self.assertEqual(state['cooldown_remaining'], 8)
        self.assertIsNone(state['pending'])
        self.assertIsNone(state['last_passed_dt_s'])
        self.assertEqual(event['accepted_time_advanced_s'], 0.)
        self.assertFalse(event['accepted'])
        adaptive.validate_state(state)

    def test_wrong_parent_dt_scope_or_diagnosis_never_authorizes_refine(self):
        pending = adaptive.begin_attempt(self.state(), self.manifest['dt_s'])
        for key, value in [('parent_sha256', '0'*64), ('attempted_dt_s', 1e-6), ('eligible', False),
                ('scope', 'FINE_1'), ('protocol', 'OTHER'), ('accepted', True), ('restart_allowed', True),
                ('whole_path_certified', True), ('nonlinear_converged', True),
                ('failed_candidate_retained', False), ('accepted_time_advanced_s', False)]:
            state, event = policy.refine(pending, dict(self.diagnosis, **{key: value}))
            self.assertEqual(event['action'], 'STOP')
            self.assertEqual(state['parent_sha256'], pending['parent_sha256'])
            self.assertEqual(state['passes'], 0)

    def test_original_eight_halvings_shared_minimum_unrelaxed(self):
        state = self.state()
        for count in range(1, 9):
            state, event = self.refine(state)
            self.assertEqual(event['action'], 'REFINE')
            self.assertEqual(state['consecutive_halvings'], count)
        self.assertEqual(self.refine(state)[1]['reason'], 'NONLINEAR_REFINEMENT_LIMIT')
        state, event = self.refine(self.state(dt=2e-7))
        self.assertEqual(state['nominal_dt_s'], 1e-7)
        self.assertEqual(self.refine(state)[1]['reason'], 'MINIMUM_DT_REACHED')
        self.assertEqual(self.refine(self.state(dt=1.9e-7, minimum_dt_s=1e-9))[1]['reason'], 'MINIMUM_DT_REACHED')

    def test_actual_clipped_interval_and_growth_cooldown(self):
        state = self.state(); state.update(small_error_streak=1, growth_pending=True, cooldown_remaining=16)
        pending = adaptive.begin_attempt(state, state['nominal_dt_s']/2, phase_clipped=True)
        result, event = policy.refine(pending, dict(self.diagnosis, attempted_dt_s=pending['pending']['dt_s']))
        self.assertEqual(result['nominal_dt_s'], state['nominal_dt_s']/4)
        self.assertEqual(result['growth_failures'], 0)  # A phase-clipped attempt is not a growth trial.
        self.assertEqual(result['small_error_streak'], 0)
        self.assertEqual(result['cooldown_remaining'], 16)
        self.assertFalse(result['growth_pending'])
        full, _ = self.refine(state)
        self.assertEqual(full['growth_failures'], 1)

    def test_frozen_adaptive_resolve_pass_remains_required_and_compatible(self):
        state, _ = self.refine(self.state())
        endpoint = 'f'*64
        row = dict(particle_id=7, position_error_m=1e-8, temperature_error_K=.001,
            stress_relative=.001, plastic_relative=.001)
        screen = dict(global_screen=dict(position_max_m=1e-8, temperature_max_K=.001,
            stress_relative=.001, plastic_relative=.001, **{'pass': True}),
            local_screen=dict(adaptive.LOCAL_LIMITS, particle_ids=[7], rows=[row], passed=True),
            crossing_ids_match=True, passed=True, fine_checkpoint_sha256=['e'*64, endpoint])
        pending = adaptive.begin_attempt(state, state['nominal_dt_s'])
        result, event = adaptive.resolve_pass(pending, screen, .0005, [7], endpoint,
            original_candidate_gates_passed=True)
        self.assertEqual(result['parent_sha256'], endpoint)
        self.assertEqual(result['passes'], 1)
        self.assertEqual(result['cooldown_remaining'], 7)
        self.assertEqual(result['nominal_dt_s'], state['nominal_dt_s'])

    def test_no_physics_gpu_filesystem_or_eval_imported(self):
        tree = ast.parse((ROOT/'scripts/afsd_nonlinear_refinement_20261011.py').read_text())
        imports = [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names]
        imports += [node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
        self.assertEqual(set(imports), {'copy', 'math', 'afsd_adaptive_dt_20261010'})
        self.assertFalse(any(isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id in {'open', 'exec', 'eval', '__import__'} for node in ast.walk(tree)))


if __name__ == '__main__':
    unittest.main()

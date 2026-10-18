"""Administrative timestep policy: fixtures and read-only saved real evidence."""
import ast
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import afsd_adaptive_dt_20261010 as policy

IDS = [7, 11]
DX = .0005
MAXIMUM = policy.DEFAULT_MAXIMUM_DT_S
INITIAL = MAXIMUM/4.


def checkpoint(number):
    return format(number, '064x')


def screen(ratio=.1, endpoint=1):
    return dict(global_screen=dict(position_max_m=ratio*.01*DX, temperature_max_K=ratio*.5,
            stress_relative=ratio*.1, plastic_relative=ratio*.2, **{'pass': ratio <= 1.}),
        local_screen=dict(policy.LOCAL_LIMITS, particle_ids=IDS[:], rows=[
            dict(particle_id=pid, position_error_m=ratio*5e-6, temperature_error_K=ratio*.5,
                stress_relative=ratio*.1, plastic_relative=ratio*.2) for pid in IDS],
            passed=ratio <= 1.), crossing_ids_match=True, passed=ratio <= 1., accepted=False,
        coarse_checkpoint_sha256=checkpoint(8000+endpoint),
        fine_checkpoint_sha256=[checkpoint(9000+endpoint), checkpoint(endpoint)])


def initial(**kwargs):
    return policy.initial_state(checkpoint(0), INITIAL, **kwargs)


def completed(state, ratio=.1, phase_dt=None, supplied=None):
    number = state['passes']+1
    value = screen(ratio, number) if supplied is None else copy.deepcopy(supplied)
    value['fine_checkpoint_sha256'][-1] = checkpoint(number)
    pending = policy.begin_attempt(state, state['nominal_dt_s'] if phase_dt is None else phase_dt,
        phase_clipped=phase_dt is not None)
    return policy.resolve_pass(pending, value, DX, IDS, checkpoint(number),
        original_candidate_gates_passed=True)


def rejected(state, *, supplied=None, manifest=None, gates=True):
    pending = policy.begin_attempt(state, state['nominal_dt_s'])
    failure = dict(status='FAIL_GPU_STAGE_NOT_ACCEPTED', error=policy.TEMPORAL_ERROR)
    return policy.resolve_failure(pending, failure if manifest is None else manifest,
        screen(1.2) if supplied is None else supplied, DX, IDS,
        original_candidate_gates_passed=gates)


class CompleteScreenTests(unittest.TestCase):
    def test_every_global_field_can_prevent_growth_even_with_small_stress(self):
        for name, threshold in dict(position_max_m=.01*DX, temperature_max_K=.5,
                stress_relative=.1, plastic_relative=.2).items():
            with self.subTest(field=name):
                value = screen(.01); value['global_screen'][name] = .8*threshold
                state, event = completed(initial(), supplied=value)
                state, event = completed(state, supplied=value)
                self.assertEqual(state['nominal_dt_s'], INITIAL)
                self.assertEqual(state['small_error_streak'], 0)
                self.assertAlmostEqual(event['normalized_screen']['ratio'], .8)

    def test_every_local_field_and_last_particle_can_prevent_growth(self):
        for name, threshold in dict(position_error_m=5e-6, temperature_error_K=.5,
                stress_relative=.1, plastic_relative=.2).items():
            with self.subTest(field=name):
                value = screen(.01); value['local_screen']['rows'][-1][name] = .8*threshold
                state, event = completed(initial(), supplied=value)
                state, event = completed(state, supplied=value)
                self.assertEqual(state['nominal_dt_s'], INITIAL)
                self.assertEqual(event['normalized_screen']['local_worst_particle_ids'][name], IDS[-1])

    def test_changed_local_limits_and_denominator_floors_rejected(self):
        for key in policy.LOCAL_LIMITS:
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'thresholds or floors'):
                value = screen(); value['local_screen'][key] *= 2.
                policy.normalized_screen(value, DX, IDS)

    def test_missing_duplicate_reordered_or_unknown_ids_rejected(self):
        for changed in ([7], [7, 7], [11, 7], [7, 12]):
            with self.subTest(ids=changed), self.assertRaisesRegex(ValueError, 'every complete-parent'):
                value = screen(); value['local_screen']['particle_ids'] = changed
                policy.normalized_screen(value, DX, IDS)
        value = screen(); value['local_screen']['rows'].reverse()
        with self.assertRaisesRegex(ValueError, 'every complete-parent'):
            policy.normalized_screen(value, DX, IDS)

    def test_invalid_values_do_not_become_growth_evidence(self):
        for invalid in (float('nan'), float('inf'), -1., True, '0', 1e308):
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                value = screen(); value['local_screen']['rows'][-1]['position_error_m'] = invalid
                policy.normalized_screen(value, DX, IDS)

    def test_declared_pass_flags_must_match_recomputed_error(self):
        for scope, key in (('global_screen', 'pass'), ('local_screen', 'passed')):
            value = screen(); value[scope][key] = False
            with self.subTest(scope=scope), self.assertRaisesRegex(ValueError, 'original limits'):
                policy.normalized_screen(value, DX, IDS)

    def test_crossing_mismatch_cannot_be_smoothed_into_a_ratio(self):
        value = screen(); value.update(crossing_ids_match=False, passed=False)
        normalized = policy.normalized_screen(value, DX, IDS)
        self.assertFalse(normalized['passed'])
        pending = policy.begin_attempt(initial(), INITIAL)
        with self.assertRaisesRegex(ValueError, 'complete original pair PASS'):
            policy.resolve_pass(pending, value, DX, IDS, checkpoint(1), original_candidate_gates_passed=True)

    def test_pair_flag_inconsistent_with_crossing_rejected(self):
        value = screen(); value['crossing_ids_match'] = False
        with self.assertRaisesRegex(ValueError, 'crossing flags'):
            policy.normalized_screen(value, DX, IDS)


class RecoveryPolicyTests(unittest.TestCase):
    def test_two_small_complete_passes_allow_one_bounded_growth(self):
        first, event = completed(initial())
        self.assertEqual(event['action'], 'HOLD'); self.assertEqual(first['nominal_dt_s'], INITIAL)
        second, event = completed(first)
        self.assertEqual(event['action'], 'GROW'); self.assertEqual(second['nominal_dt_s'], 2.*INITIAL)
        self.assertTrue(second['growth_pending'])
        self.assertEqual(second['parent_sha256'], checkpoint(2))

    def test_growth_returns_to_original_ceiling_without_exceeding_it(self):
        state = initial()
        for unused in range(10):
            state, event = completed(state)
            self.assertLessEqual(state['nominal_dt_s'], MAXIMUM)
        self.assertEqual(state['nominal_dt_s'], MAXIMUM)
        self.assertEqual(event['reason'], 'REGISTERED_MAXIMUM_DT')

    def test_marginal_pass_resets_consecutive_small_error_streak(self):
        state, unused = completed(initial())
        state, event = completed(state, .8)
        self.assertEqual(state['small_error_streak'], 0)
        state, event = completed(state)
        self.assertEqual(event['action'], 'HOLD')
        state, event = completed(state)
        self.assertEqual(event['action'], 'GROW')

    def test_phase_endpoint_does_not_permanently_shrink_nominal_dt_or_grant_growth(self):
        state, unused = completed(initial())
        state, event = completed(state, phase_dt=INITIAL/3.)
        self.assertEqual(state['nominal_dt_s'], INITIAL)
        self.assertEqual(state['last_passed_dt_s'], INITIAL/3.)
        self.assertEqual(state['small_error_streak'], 0)
        self.assertEqual(event['reason'], 'PHASE_ENDPOINT_PRESERVES_NOMINAL_DT')
        state, event = completed(state)
        self.assertEqual(state['nominal_dt_s'], INITIAL)

    def test_rounding_adjusted_phase_boundary_is_accepted(self):
        pending = policy.begin_attempt(initial(), INITIAL*(1.+2e-16), phase_clipped=True)
        state, unused = policy.resolve_pass(pending, screen(), DX, IDS, checkpoint(1),
            original_candidate_gates_passed=True)
        self.assertEqual(state['nominal_dt_s'], INITIAL)

    def test_driver_cannot_silently_shorten_or_grow_nominal_proposal(self):
        for actual, clipped in ((INITIAL/2., False), (2.*INITIAL, False), (2.*INITIAL, True)):
            with self.subTest(actual=actual, clipped=clipped), self.assertRaises(ValueError):
                policy.begin_attempt(initial(), actual, phase_clipped=clipped)

    def test_only_verified_exact_fine_two_endpoint_advances_parent(self):
        pending = policy.begin_attempt(initial(), INITIAL)
        for endpoint, gates in ((checkpoint(0), True), (checkpoint(3), True), (checkpoint(1), False),
                (checkpoint(1), 1)):
            with self.subTest(endpoint=endpoint, gates=gates), self.assertRaises(ValueError):
                policy.resolve_pass(pending, screen(), DX, IDS, endpoint,
                    original_candidate_gates_passed=gates)

    def test_state_is_json_serializable_and_calls_do_not_mutate_inputs(self):
        state = initial(); before = copy.deepcopy(state)
        pending = policy.begin_attempt(state, INITIAL)
        self.assertEqual(state, before)
        pending_before = copy.deepcopy(pending); value = screen(); value_before = copy.deepcopy(value)
        result, event = policy.resolve_pass(pending, value, DX, IDS, checkpoint(1),
            original_candidate_gates_passed=True)
        self.assertEqual(pending, pending_before); self.assertEqual(value, value_before)
        self.assertEqual(json.loads(json.dumps(result)), result)
        self.assertEqual(json.loads(json.dumps(event)), event)

    def test_resolved_or_duplicate_pending_attempt_cannot_be_reused(self):
        pending = policy.begin_attempt(initial(), INITIAL)
        with self.assertRaises(ValueError):
            policy.begin_attempt(pending, INITIAL)
        ready, unused = completed(initial())
        with self.assertRaises(ValueError):
            policy.resolve_pass(ready, screen(), DX, IDS, checkpoint(1), original_candidate_gates_passed=True)

    def test_invalid_policy_and_nonfinite_inputs_reject(self):
        for changes in (dict(maximum_dt_s=.00501), dict(maximum_dt_s=0.), dict(minimum_dt_s=1e-6),
                dict(growth_factor=2.01), dict(growth_factor=True), dict(small_error_ratio=.251),
                dict(required_small_passes=1), dict(required_small_passes=True),
                dict(cooldown_passes=0), dict(max_halvings=9)):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                initial(**changes)
        self.assertEqual(initial(maximum_dt_s=.005)['config']['maximum_dt_s'], .005)

    def test_module_imports_only_standard_library_and_no_physical_code(self):
        tree = ast.parse(Path(policy.__file__).read_text(encoding='utf-8'))
        names = {alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
            for alias in node.names}
        self.assertEqual(names, {'copy', 'math'})
        self.assertFalse(any(isinstance(node, ast.ImportFrom) for node in ast.walk(tree)))


class FailurePolicyTests(unittest.TestCase):
    def test_temporal_refinement_retries_exact_unchanged_complete_parent(self):
        state = initial(); smaller, event = rejected(state)
        self.assertEqual(event['action'], 'REFINE')
        self.assertEqual(smaller['nominal_dt_s'], INITIAL/2.)
        self.assertEqual(smaller['parent_sha256'], state['parent_sha256'])
        self.assertEqual(smaller['passes'], 0)
        self.assertGreater(smaller['cooldown_remaining'], 0)

    def test_failed_growth_cools_down_before_another_larger_trial(self):
        state, unused = completed(initial()); state, unused = completed(state)
        state, event = rejected(state)
        self.assertTrue(event['growth_trial']); self.assertEqual(state['nominal_dt_s'], INITIAL)
        self.assertEqual(state['cooldown_remaining'], 4)
        for remaining in (3, 2, 1, 0):
            state, event = completed(state)
            self.assertEqual(event['reason'], 'GROWTH_RETRY_COOLDOWN')
            self.assertEqual(state['cooldown_remaining'], remaining)
            self.assertEqual(state['nominal_dt_s'], INITIAL)
        state, event = completed(state); self.assertEqual(event['action'], 'HOLD')
        state, event = completed(state); self.assertEqual(event['action'], 'GROW')
        state, event = rejected(state)
        self.assertEqual(state['cooldown_remaining'], 8)

    def test_phase_clip_does_not_consume_cooldown(self):
        state, unused = rejected(initial())
        state, event = completed(state, phase_dt=INITIAL/4.)
        self.assertEqual(state['cooldown_remaining'], 4)
        self.assertEqual(state['nominal_dt_s'], INITIAL/2.)

    def test_physical_source_cuda_geometry_and_unrecognized_errors_stop(self):
        for error in ("ValueError('energy ledger failed')", "ValueError('geometry gate failed')",
                "RuntimeError('CUDA error')", "ValueError('source changed')", 'other temporal fail'):
            with self.subTest(error=error):
                state, event = rejected(initial(), manifest=dict(status='FAIL_GPU_STAGE_NOT_ACCEPTED', error=error))
                self.assertEqual(event['action'], 'STOP'); self.assertEqual(state['nominal_dt_s'], INITIAL)
                self.assertEqual(state['parent_sha256'], checkpoint(0))

    def test_temporal_error_without_all_candidate_gates_cannot_refine(self):
        for gates in (False, None, 1):
            with self.subTest(gates=gates):
                state, event = rejected(initial(), gates=gates)
                self.assertEqual(event['action'], 'STOP')

    def test_crossing_mismatch_stops_instead_of_halving(self):
        value = screen(1.2); value['crossing_ids_match'] = False
        state, event = rejected(initial(), supplied=value)
        self.assertEqual(event['action'], 'STOP')

    def test_missing_corrupt_or_inconsistent_temporal_screens_stop(self):
        values = [dict(passed=False, crossing_ids_match=True), screen()]
        nonfinite = screen(1.2); nonfinite['global_screen']['stress_relative'] = float('nan'); values.append(nonfinite)
        ids = screen(1.2); ids['local_screen']['particle_ids'].pop(); values.append(ids)
        for value in values:
            with self.subTest(value=list(value)):
                state, event = rejected(initial(), supplied=value)
                self.assertEqual(event['action'], 'STOP')

    def test_halving_limit_is_enforced_without_publishing_failed_endpoint(self):
        state = initial(max_halvings=1)
        state, event = rejected(state); self.assertEqual(event['action'], 'REFINE')
        state, event = rejected(state)
        self.assertEqual(event['reason'], 'TEMPORAL_REFINEMENT_LIMIT')
        self.assertEqual(state['parent_sha256'], checkpoint(0)); self.assertEqual(state['status'], 'STOPPED')

    def test_minimum_dt_is_enforced(self):
        state = policy.initial_state(checkpoint(0), 1e-7)
        state, event = rejected(state)
        self.assertEqual(event['reason'], 'MINIMUM_DT_REACHED')


class SavedEvidenceTests(unittest.TestCase):
    def read_screen(self, attempt):
        path = ROOT/f'outputs/gpu_acceleration_20261007/campaign_4h_20261009_v1/attempt_{attempt:04d}/pair_screen.json'
        return json.loads(path.read_text(encoding='utf-8'))

    def test_latest_real_1381_particle_screen_has_small_error_in_all_fields(self):
        value = self.read_screen(52)
        result = policy.normalized_screen(value, DX, list(range(1, 1382)))
        self.assertEqual(result['compared_particle_count'], 1381)
        self.assertTrue(result['passed']); self.assertLess(result['ratio'], .05)
        self.assertGreater(result['ratio'], .04)
        self.assertEqual(result['ratio'], max(row['stress_relative']/.1 for row in value['local_screen']['rows']))

    def test_old_stress_failures_remain_failures_under_new_policy(self):
        for attempt, lower, upper in ((9, 3.18, 3.20), (10, 1.68, 1.69)):
            with self.subTest(attempt=attempt):
                result = policy.normalized_screen(self.read_screen(attempt), DX, list(range(1, 1382)))
                self.assertFalse(result['passed']); self.assertGreater(result['ratio'], lower)
                self.assertLess(result['ratio'], upper)


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Scheduling regression tests; no physical runner import or subprocess solve."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import campaign as c


class SchedulerTests(unittest.TestCase):
    def meta(self):
        return dict(next_proposed_macro_dt_s=.00025, stage_target_reached=True,
            stage_end_time_s=.060, phase_end_time_s=.22,
            pair_screens=[dict(macro_index=1, dt_s=.0005, time_s=.05975, refinement_level=0,
                               maximum_original_tolerance_utilization=.4),
                          dict(macro_index=2, dt_s=.00025, time_s=.060, refinement_level=0,
                               maximum_original_tolerance_utilization=.4)])

    def attempts(self):
        return [dict(macro_index=1, level=0, passed=True), dict(macro_index=2, level=0, passed=True)]

    def test_absolute_deadline_counts_preparation_and_never_resets(self):
        self.assertEqual(c.DEADLINE_EPOCH-c.START_EPOCH, 7200.)
        mono_deadline = 50.+(c.DEADLINE_EPOCH-(c.START_EPOCH+1000.))
        self.assertEqual(c.remaining_seconds(c.START_EPOCH+2000., 1050., mono_deadline), 5200.)
        self.assertFalse(c.can_launch(c.CHILD_WRAPPER_LIMIT+c.FINAL_RESERVE+c.TERMINATION_RESERVE-1))
        self.assertTrue(c.can_launch(c.CHILD_WRAPPER_LIMIT+c.FINAL_RESERVE+c.TERMINATION_RESERVE))

    def test_wall_clock_rollback_does_not_extend(self):
        self.assertEqual(c.remaining_seconds(c.START_EPOCH-1000., 80., 100.), 20.)
        self.assertLess(c.remaining_seconds(c.DEADLINE_EPOCH+1., 1., 100.), 0.)

    def test_actual_cost_quota_and_original_forecast_limit(self):
        q = c.budget(9.68510099996638, 18.317446299999574)
        self.assertEqual(q['candidate_quota'], 15)
        self.assertEqual(q['max_macrosteps'], 5)
        self.assertLessEqual(q['forecast_s'], 450.)
        self.assertEqual(c.budget(1., 2.)['candidate_quota'], 27)
        self.assertFalse(c.budget(100., 120.)['runnable'])

    def test_sixty_ms_is_observation_not_end_of_campaign(self):
        self.assertEqual(c.endpoint(.024271484375), .060)
        self.assertEqual(c.endpoint(.060), .120439208984375)
        self.assertGreater(c.endpoint(.060), .060)

    def test_proven_endpoint_truncation_restores_passed_half_ms(self):
        result = c.resume_dt(self.meta(), self.attempts())
        self.assertEqual(result['dt_s'], .0005)
        self.assertEqual(result['reason'], 'RESTORED_ALREADY_PASSED_DT_AFTER_ENDPOINT_TRUNCATION')

    def test_temporal_halving_preserves_runner_next_dt(self):
        m = self.meta(); m['pair_screens'][-1]['refinement_level'] = 1
        self.assertEqual(c.resume_dt(m, self.attempts())['dt_s'], .00025)
        a = self.attempts(); a.append(dict(macro_index=3, level=0, passed=False))
        self.assertEqual(c.resume_dt(self.meta(), a)['dt_s'], .00025)

    def test_incomplete_endpoint_and_unpassed_growth_never_restore(self):
        m = self.meta(); m['stage_target_reached'] = False
        self.assertEqual(c.resume_dt(m, self.attempts())['dt_s'], .00025)
        m = self.meta(); m['pair_screens'][0].update(dt_s=.0003, maximum_original_tolerance_utilization=.2)
        self.assertEqual(c.resume_dt(m, self.attempts())['dt_s'], .00025)

    def test_atomic_unique_output_blocks_second_instance_and_stale_reuse(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            with patch.object(c, 'BASE', base):
                c.claim_output(base/'run')
                self.assertTrue((base/'run/instance.lock.json').is_file())
                with self.assertRaises(FileExistsError):
                    c.claim_output(base/'run')
                with self.assertRaises(ValueError):
                    c.claim_output(base.parent/'outside')

    def test_source_drift_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)/'source.py'; p.write_text('original')
            digest = c.old.sha(p); p.write_text('drift')
            with patch.object(c.old, 'inside', return_value=p):
                with self.assertRaises(ValueError):
                    c.old.verify({'source.py': digest})

    def test_full_prefix_is_valid_but_failure_no_growth_wrong_parent_are_not(self):
        parent = c.PARENT
        pm = dict(screened_continuation_sha256='parent-sha', screened_end_time_s=.024,
                  screened_end_step=10)
        good = dict(status='PASS_SCREENED_PREFIX_BUDGET_BOUND_NOT_TARGET',
            parent_stage_manifest=c.old.rel(parent/'manifest.json'), parent_stage_manifest_sha256='manifest-sha',
            parent_checkpoint_sha256='parent-sha', screened_end_time_s=.025, screened_end_step=12)
        with patch.object(c.old, 'sha', return_value='manifest-sha'), patch.object(c.old, 'snapshot', return_value={'ok': True}):
            self.assertEqual(c.validate_child(good, parent, parent, pm), {'ok': True})
            for changes in ({'status': 'FAIL'}, {'screened_end_time_s': .024},
                            {'screened_end_step': 10}, {'parent_checkpoint_sha256': 'wrong'},
                            {'screened_end_time_s': .061}):
                with self.assertRaises(ValueError):
                    c.validate_child(dict(good, **changes), parent, parent, pm)

    def test_failure_receipt_cannot_promote_even_with_pass_manifest(self):
        with patch.object(c.old, 'read', side_effect=[{'status': next(iter(c.old.PASS_STATUSES))}, {'exit_code': 2}]):
            with self.assertRaises(ValueError):
                c.old.verify_receipt(c.PARENT, c.old.PASS_STATUSES)

    def test_import_is_stdlib_and_frozen_supervisor_only(self):
        self.assertNotIn('run_afsd_deposition_stage_remote_20261007', c.sys.modules)
        self.assertFalse(any(k.startswith('afsd_backend') for k in c.sys.modules))


if __name__ == '__main__':
    unittest.main()

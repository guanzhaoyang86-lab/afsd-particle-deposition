"""Pure finite-budget/stop checks; no physical solver or runner imports."""
from pathlib import Path
import tempfile
import unittest
from supervise_afsd_deposition_20261007 import (candidate_budget, stop_after_stage,
    inside, atomic_json, read, PASS_STATUSES)


class DepositionSupervisorTests(unittest.TestCase):
    def test_quota_matches_measured_cost_and_never_exceeds_forecast(self):
        result = candidate_budget(8., 30.)
        self.assertEqual(result['candidate_quota'], 9)
        self.assertEqual(result['max_macrosteps'], 3)
        self.assertLessEqual(result['forecast_s'], 450.)

    def test_fast_candidate_quota_is_capped_and_slow_candidate_stops(self):
        self.assertEqual(candidate_budget(1., 2.)['candidate_quota'], 18)
        self.assertFalse(candidate_budget(100., 120.)['runnable'])
        self.assertTrue(candidate_budget(1., (450.-30.)/4.5)['runnable'])

    def test_failed_child_stops_even_if_it_wrote_new_time_or_shoulder_mass(self):
        self.assertEqual(stop_after_stage(.003, .004, 'FAIL', 1., .22), 'CHILD_FAILED_OR_NOT_SCREENED')

    def test_no_increase_is_not_hidden_by_a_pass_label(self):
        status = next(iter(PASS_STATUSES))
        self.assertEqual(stop_after_stage(.003, .003, status, 0., .22), 'NO_PHYSICAL_TIME_INCREASE')

    def test_first_shoulder_and_global_target_are_independent_stops(self):
        status = next(iter(PASS_STATUSES))
        self.assertEqual(stop_after_stage(.003, .004, status, 1e-9, .22), 'FIRST_SCREENED_NEW_MATERIAL_UNDER_SHOULDER')
        self.assertEqual(stop_after_stage(.219, .22, status, 0., .22), 'GLOBAL_DIAGNOSTIC_TARGET_REACHED')
        self.assertIsNone(stop_after_stage(.003, .004, status, 0., .22))

    def test_output_path_escape_is_rejected_and_atomic_progress_is_readable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(ValueError):
                inside(Path('../outside'), root)
            target = root/'progress.json'
            atomic_json(target, dict(status='RUNNING', physical_time_s=.003))
            atomic_json(target, dict(status='PASS', physical_time_s=.004))
            self.assertEqual(read(target), dict(status='PASS', physical_time_s=.004))
            self.assertEqual(list(root.glob('*.pending')), [])


if __name__ == '__main__':
    unittest.main()

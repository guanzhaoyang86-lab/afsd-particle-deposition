"""Wiring certificates; no physical AFSD advance or manufactured PASS receipt."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_process_event as base
import run_process_event_refined as refined
from afsd_backend.pysph.process_thermal_trial import ProcessThermalTrial, acceptance
from afsd_backend.pysph.stage_structural_root import solve_stage


class RefinedDriver(unittest.TestCase):
    def test_numerical_strategy_only_original_physical_interfaces(self):
        saved = base.response_from_trace, base.SOURCES, base.SETTINGS
        try:
            wired = refined.configure()
            self.assertIs(wired.ProcessThermalTrial, ProcessThermalTrial)
            self.assertIs(wired.acceptance, acceptance)
            self.assertIs(wired.solve_stage, solve_stage)
            self.assertIs(wired.response_from_trace, refined.response_from_trace)
            self.assertEqual(wired.SETTINGS['dt_s'], 5e-7)
            self.assertEqual(wired.SETTINGS['maximum_evaluations'], 32)
            self.assertEqual(wired.SETTINGS['maximum_steps'], 1)
            self.assertEqual(wired.SETTINGS['inner_limit_s'], 540.)
            self.assertEqual(wired.SETTINGS['force_difference_maximum_halvings'], 2)
            self.assertTrue(set(refined.EXTRA_SOURCES) <= set(wired.SOURCES))
            original = set(saved[1]); self.assertTrue(original <= set(wired.SOURCES))
            self.assertEqual(base.cost_budget(dict(preconditioner_build_wall_seconds=2.,
                trial_proposal_journal_wall_seconds=4.))['forecast_wall_seconds'], 225.)
            self.assertFalse(base.cost_budget(dict(preconditioner_build_wall_seconds=2.,
                trial_proposal_journal_wall_seconds=10.))['allowed'])
        finally:
            base.response_from_trace, base.SOURCES, base.SETTINGS = saved


if __name__ == '__main__':
    unittest.main(verbosity=2)

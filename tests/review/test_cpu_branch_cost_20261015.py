"""Administrative cost-model regression checks; no physical solve is executed."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import afsd_cpu_branch_cost_20261015 as costs


class BranchCostForecast(unittest.TestCase):
    def test_empty_model_cannot_authorize_a_branch(self):
        result = costs.forecast(costs.empty_model())
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "NO_COMPLETE_BRANCH_SAMPLES")

    def test_observed_slow_probe_no_longer_assumes_all_twelve_solves_are_slow(self):
        # Observed 004 coarse wall and worst probe; other values are synthetic.
        model = costs.update(costs.empty_model(), 163.56472830000075, 10,
                             [34.5509886] + [12.0] * 9)
        result = costs.forecast(model)
        self.assertLess(result["forecast_s"], 450.0)
        self.assertTrue(result["eligible"])
        self.assertAlmostEqual(result["forecast_s"], 324.4165109400014)
        self.assertEqual(model["maximum_single_solve_s"], 34.5509886)

    def test_no_event_samples_contribute_without_event_selection_bias(self):
        model = costs.update(costs.empty_model(), 18.0, 1, [12.0])
        self.assertEqual(model["maximum_projected_branch_s"], 216.0)
        self.assertEqual(costs.forecast(model)["forecast_s"], 354.0)

    def test_entire_supervised_wall_including_proof_and_output_is_charged(self):
        low = costs.update(costs.empty_model(), 12.0, 1, [12.0])
        high = costs.update(costs.empty_model(), 20.0, 1, [12.0])
        self.assertGreater(costs.forecast(high)["forecast_s"], costs.forecast(low)["forecast_s"])
        self.assertEqual(costs.forecast(high)["forecast_s"], 390.0)

    def test_actual_projection_over_limit_is_rejected(self):
        model = costs.update(costs.empty_model(), 250.0, 10, [20.0] * 10)
        result = costs.forecast(model)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["forecast_s"], 480.0)
        self.assertEqual(result["reason"], "BRANCH_COST_REASSESSMENT_REQUIRED")

    def test_forecast_limit_is_inclusive(self):
        exact = costs.update(costs.empty_model(), 280.0, 12, [20.0] * 12)
        above = costs.update(costs.empty_model(), 280.000001, 12, [20.0] * 12)
        self.assertTrue(costs.forecast(exact)["eligible"])
        self.assertFalse(costs.forecast(above)["eligible"])

    def test_single_solve_floor_is_applied_even_when_projection_is_lower(self):
        model = costs.update(costs.empty_model(), 200.0, 12, [160.0] + [1.0] * 11)
        result = costs.forecast(model)
        self.assertEqual(result["single_solve_floor_s"], 320.0)
        self.assertEqual(result["forecast_s"], 510.0)
        self.assertFalse(result["eligible"])

    def test_history_only_increases_including_single_solve_cost(self):
        first = costs.update(costs.empty_model(), 18.0, 1, [12.0])
        second = costs.update(first, 100.0, 10, [8.0] * 10)
        third = costs.update(second, 130.0, 12, [110.0] + [1.0] * 11)
        self.assertEqual(second["maximum_projected_branch_s"], 216.0)
        self.assertEqual(second["maximum_single_solve_s"], 12.0)
        self.assertEqual(third["maximum_projected_branch_s"], 216.0)
        self.assertEqual(third["maximum_single_solve_s"], 110.0)
        self.assertGreater(costs.forecast(third)["forecast_s"], costs.forecast(second)["forecast_s"])
        self.assertEqual(third["sample_count"], 3)

    def test_no_mutation_on_success_or_failure_or_forecast(self):
        original = costs.empty_model()
        before = deepcopy(original)
        sample = [10.0]
        result = costs.update(original, 11.0, 1, sample)
        self.assertIsNot(result, original)
        self.assertEqual(original, before)
        self.assertEqual(sample, [10.0])
        copied = deepcopy(result)
        costs.forecast(result)
        with self.assertRaises(ValueError):
            costs.update(result, 1.0, 1, sample)
        self.assertEqual(result, copied)
        self.assertEqual(sample, [10.0])
        original["sample_count"] = 7
        self.assertEqual(costs.empty_model()["sample_count"], 0)

    def test_invalid_wall_and_solve_costs_rejected(self):
        invalid = [0.0, -1.0, float("nan"), float("inf"), -float("inf"), True, "10", None]
        for value in invalid:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    costs.update(costs.empty_model(), value, 1, [1.0])
                with self.assertRaises(ValueError):
                    costs.update(costs.empty_model(), 10.0, 1, [value])

    def test_invalid_solve_counts_and_cost_list_lengths_rejected(self):
        for value in (True, False, 0, -1, 13, 1.0, float("nan"), "1", None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                costs.update(costs.empty_model(), 10.0, value, [1.0])
        for value in ([], [1.0, 1.0], None, "1", {0: 1.0}):
            with self.subTest(value=value), self.assertRaises(ValueError):
                costs.update(costs.empty_model(), 10.0, 1, value)

    def test_serial_solve_total_must_fit_full_wall(self):
        with self.assertRaises(ValueError):
            costs.update(costs.empty_model(), 10.0, 2, [6.0, 6.0])
        costs.update(costs.empty_model(), 12.0 - 1e-8, 2, [6.0, 6.0])

    def test_bad_or_relaxed_policy_is_rejected(self):
        changes = {"maximum_branch_solves": [0, 13, True, 12.0],
                   "forecast_limit_s": [600.0, -1, float("nan")],
                   "fixed_overhead_s": [0.0, True], "safety_factor": [1.0, float("inf")],
                   "single_solve_floor_multiplier": [1.0],
                   "internal_limit_s": [600.0], "external_limit_s": [650.0],
                   "schema": ["legacy"]}
        for key, values in changes.items():
            for value in values:
                model = costs.empty_model()
                model[key] = value
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    costs.forecast(model)

    def test_bad_history_and_schema_rejected(self):
        for field, value in (("maximum_projected_branch_s", -1.0),
                             ("maximum_single_solve_s", float("nan")),
                             ("maximum_single_solve_s", True),
                             ("sample_count", True), ("sample_count", 1.0),
                             ("sample_count", -1), ("maximum_single_solve_s", 1.0)):
            model = costs.empty_model()
            model[field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                costs.forecast(model)
        for model in (None, {}, dict(costs.empty_model(), unexpected=True)):
            with self.subTest(model=model), self.assertRaises(ValueError):
                costs.forecast(model)

    def test_forged_measured_history_cannot_qualify_with_zero_or_inconsistent_cost(self):
        for projected, single in ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (5.0, 6.0)):
            model = dict(costs.empty_model(), sample_count=1,
                         maximum_projected_branch_s=projected,
                         maximum_single_solve_s=single)
            with self.subTest(projected=projected, single=single), self.assertRaises(ValueError):
                costs.forecast(model)

    def test_numeric_overflow_rejected(self):
        with self.assertRaises(ValueError):
            costs.update(costs.empty_model(), 1e308, 1, [1.0])
        with self.assertRaises(ValueError):
            costs.update(costs.empty_model(), 1e308, 2, [1e308, 1e308])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
import unittest
from pathlib import Path

from afsd_scaffold.config import (
    ConfigurationError,
    assert_physics_ready,
    load_resolved_case,
    unresolved_parameters,
    validate_parameter_tree,
)


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "configs" / "cases"


class CaseConfigurationTests(unittest.TestCase):
    def test_three_cases_resolve_and_keep_si_values(self) -> None:
        expected_feed_speeds = {
            "starved": 0.00106,
            "copper-core": 0.00212,
            "overfed": 0.00424,
        }
        for case_id, expected_speed in expected_feed_speeds.items():
            with self.subTest(case_id=case_id):
                config = load_resolved_case(CASE_DIR / f"{case_id}.json")
                self.assertEqual(config["case_id"], case_id)
                self.assertEqual(config["calculation_units"]["system"], "SI")
                self.assertEqual(config["display_units"]["feed_speed"], "mm/s")
                speed = config["parameters"]["process"]["feedstock_axial_speed"]
                self.assertEqual(speed["unit"], "m/s")
                self.assertAlmostEqual(speed["value"], expected_speed)
                self.assertGreater(len(unresolved_parameters(config)), 0)

    def test_unknowns_are_null_not_zero(self) -> None:
        config = load_resolved_case(CASE_DIR / "starved.json")
        for path in unresolved_parameters(config):
            node = config["parameters"]
            for component in path.split(".")[1:]:
                node = node[component]
            self.assertIsNone(node["value"], path)

    def test_each_authored_file_is_valid_json(self) -> None:
        paths = [ROOT / "configs" / "base.json", *sorted(CASE_DIR.glob("*.json"))]
        for path in paths:
            with self.subTest(path=path.name):
                with path.open("r", encoding="utf-8") as stream:
                    self.assertIsInstance(json.load(stream), dict)

    def test_unresolved_parameter_rejects_fallback_zero(self) -> None:
        with self.assertRaisesRegex(ConfigurationError, "must use JSON null"):
            validate_parameter_tree(
                {
                    "bad": {
                        "value": 0,
                        "unit": "m",
                        "source": "missing paper value",
                        "status": "UNRESOLVED",
                    }
                }
            )

    def test_calculation_parameter_rejects_display_units(self) -> None:
        for forbidden_unit in ("mm", "MPa", "rpm", "degC"):
            with self.subTest(unit=forbidden_unit), self.assertRaisesRegex(
                ConfigurationError, "not an allowed SI calculation unit"
            ):
                validate_parameter_tree(
                    {
                        "bad": {
                            "value": 1.0,
                            "unit": forbidden_unit,
                            "source": "test-only invalid display unit",
                            "status": "ASSUMED",
                        }
                    }
                )

    def test_resolved_dimensional_values_are_finite_numbers(self) -> None:
        for invalid_value in (True, "1.0", float("nan"), float("inf"), -float("inf")):
            with self.subTest(value=invalid_value), self.assertRaisesRegex(
                ConfigurationError, "finite int/float"
            ):
                validate_parameter_tree(
                    {
                        "bad": {
                            "value": invalid_value,
                            "unit": "m",
                            "source": "test-only invalid dimensional value",
                            "status": "ASSUMED",
                        }
                    }
                )
        validate_parameter_tree(
            {
                "valid": {
                    "value": 1.25,
                    "unit": "m",
                    "source": "test-only finite dimensional value",
                    "status": "ASSUMED",
                }
            }
        )

    def test_dimensionless_values_allow_enums_booleans_and_finite_counts(self) -> None:
        for valid_value in ("ENUM_VALUE", False, 7, 0.25):
            with self.subTest(valid_value=valid_value):
                validate_parameter_tree(
                    {
                        "valid": {
                            "value": valid_value,
                            "unit": "1",
                            "source": "test-only dimensionless value",
                            "status": "ASSUMED",
                        }
                    }
                )
        for invalid_value in (float("nan"), float("inf"), [], {}):
            with self.subTest(invalid_value=invalid_value), self.assertRaises(
                ConfigurationError
            ):
                validate_parameter_tree(
                    {
                        "bad": {
                            "value": invalid_value,
                            "unit": "1",
                            "source": "test-only invalid dimensionless value",
                            "status": "ASSUMED",
                        }
                    }
                )

    def test_default_and_paper_comparison_gates_check_all_parameters(self) -> None:
        config = load_resolved_case(CASE_DIR / "starved.json")
        with self.assertRaisesRegex(ConfigurationError, "physics run blocked"):
            assert_physics_ready(config)
        with self.assertRaisesRegex(ConfigurationError, "full configuration"):
            assert_physics_ready(
                config,
                required_paths={"process.feedstock_axial_speed"},
                run_classification="PAPER_COMPARISON",
            )

    def test_debug_gate_can_use_explicit_nonempty_backend_path_set(self) -> None:
        config = load_resolved_case(CASE_DIR / "starved.json")
        assert_physics_ready(
            config,
            required_paths={
                "parameters.process.feedstock_axial_speed",
                "material.initial_density",
            },
            run_classification="DEBUG",
        )
        with self.assertRaisesRegex(ConfigurationError, "pressure_closure"):
            assert_physics_ready(
                config,
                required_paths={"material.pressure_closure"},
                run_classification="DEBUG",
            )
        with self.assertRaisesRegex(ConfigurationError, "may not be empty"):
            assert_physics_ready(
                config, required_paths=set(), run_classification="DEBUG"
            )

    def test_fully_resolved_paper_comparison_passes_gate(self) -> None:
        config = load_resolved_case(CASE_DIR / "starved.json")
        for path in unresolved_parameters(config):
            parameter = config["parameters"]
            for component in path.split(".")[1:]:
                parameter = parameter[component]
            parameter["value"] = "TEST_ONLY" if parameter["unit"] == "1" else 1.0
            parameter["status"] = "ASSUMED"
        assert_physics_ready(config, run_classification="PAPER_COMPARISON")


if __name__ == "__main__":
    unittest.main()


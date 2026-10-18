"""Real rejected-diagnostic classification and negative administrative checks."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import afsd_cpu_multi_event_retry_20261015 as retry

FIXTURE = ROOT / "outputs/short_bar_event_continue_20261015_v2/run/continuation_0003/coarse"


class MultiEventRetry(unittest.TestCase):
    def setUp(self):
        self.report, self.route, self.path = [json.loads((FIXTURE / name).read_text(encoding="utf-8"))
            for name in ("manifest.json", "visibility_event_routing.json", "visible_path_audit.json")]
        self.extra = {key: json.loads((FIXTURE / name).read_text(encoding="utf-8"))
            for key, name in (("union_route", "proposal/union_visibility_events.json"),
                              ("shell", "continuous_shell_audit.json"), ("geometry", "geometry_audit.json"))}

    def allowed(self):
        return retry.eligible(self.report, self.route, self.path, **self.extra)

    def test_real_separated_two_opening_failure_can_be_retried(self):
        self.assertTrue(self.allowed())
        self.assertEqual({(item["particle_id_i"], item["particle_id_j"])
                          for item in self.route["transition_brackets"]}, {(1229, 1251), (1091, 1213)})

    def test_no_mutation(self):
        before = deepcopy((self.report, self.route, self.path, self.extra))
        self.allowed()
        self.assertEqual(before, (self.report, self.route, self.path, self.extra))

    def test_single_and_excess_events_rejected(self):
        original = deepcopy(self.route["transition_brackets"])
        for changes in ([], original[:1], original * 3):
            self.route["transition_brackets"] = changes
            self.assertFalse(self.allowed())

    def test_closing_event_or_unknown_status_rejected(self):
        for field, value in (("direction", "VISIBLE_TO_BLOCKED"), ("status", "UNKNOWN"),
                             ("start_blocked", False), ("end_blocked", True),
                             ("start_blocked", 1), ("single_event_certified", True),
                             ("physical_substeps_solved", True)):
            original = self.route["transition_brackets"][0][field]
            self.route["transition_brackets"][0][field] = value
            with self.subTest(field=field, value=value):
                self.assertFalse(self.allowed())
            self.route["transition_brackets"][0][field] = original

    def test_unknown_samples_path_degeneracy_and_numerical_failures_rejected(self):
        for target, field in ((self.route, "unresolved_samples"), (self.path, "degenerate_edges"),
                              (self.path, "numerical_domain_failures"), (self.path, "unresolved_intervals")):
            target[field] = [{"unexpected": True}]
            with self.subTest(field=field):
                self.assertFalse(self.allowed())
            target[field] = []

    def test_nonfinite_diagnostics_anywhere_rejected(self):
        for target in (self.report, self.route, self.path):
            for value in (float("nan"), float("inf"), -float("inf")):
                target["additional_diagnostic"] = {"sample": [value]}
                self.assertFalse(self.allowed())
            del target["additional_diagnostic"]

    def test_invalid_fractions_and_zero_or_inconsistent_width_rejected(self):
        bracket = self.route["transition_brackets"][0]
        original = deepcopy(bracket)
        for field, value in (("start_fraction", 0), ("start_fraction", -1),
                             ("end_fraction", 1.1), ("start_fraction", True),
                             ("bracket_width_fraction", 0), ("bracket_width_fraction", 0.5)):
            bracket[field] = value
            with self.subTest(field=field, value=value):
                self.assertFalse(self.allowed())
            bracket.update(original)

    def test_overlapping_and_touching_brackets_rejected(self):
        first, second = self.route["transition_brackets"]
        for lo in (first["start_fraction"], first["end_fraction"]):
            second["start_fraction"] = lo
            second["end_fraction"] = lo + first["bracket_width_fraction"]
            self.assertFalse(self.allowed())

    def test_duplicate_edge_or_mismatched_changed_endpoint_and_found_sets_rejected(self):
        for field in ("transition_brackets", "changed_edges", "endpoint_changed_edges"):
            original = deepcopy(self.route[field])
            self.route[field][1] = deepcopy(self.route[field][0])
            with self.subTest(field=field):
                self.assertFalse(self.allowed())
            self.route[field] = original
        self.path["found_edges"] = self.path["found_edges"][:1]
        self.assertFalse(self.allowed())

    def test_repeated_same_edge_at_zero_allowed_but_later_intersection_rejected(self):
        another = deepcopy(self.path["found_edges"][0])
        self.path["found_edges"].append(another)
        self.assertTrue(self.allowed())
        another["fraction"] = 0.1
        self.assertFalse(self.allowed())

    def test_unrelated_found_edge_or_bad_intersection_rejected(self):
        item = self.path["found_edges"][0]
        original = deepcopy(item)
        for field, value in (("particle_id_i", 42), ("edge_index", 42),
                             ("fraction", -0.1), ("maximum_contiguous_inside_length_m", 0),
                             ("strict_intervals", []), ("strict_intervals", [[0.5, 0.2]])):
            item[field] = value
            self.assertFalse(self.allowed())
            item.update(original)

    def test_wrong_report_error_trace_or_acceptance_flags_rejected(self):
        original = deepcopy(self.report)
        for field, value in (("status", "PASS_CPU_EVENT_BRANCH"), ("error", "GeometryPathStop('other')"),
                             ("traceback", "only one observed activation per event solve"),
                             ("accepted", True), ("restart_allowed", True),
                             ("diagnostic_only", False), ("screened_continuation_authorized", True),
                             ("solves", 2), ("solves", True)):
            self.report[field] = value
            self.assertFalse(self.allowed())
            self.report = deepcopy(original)

    def test_wrong_route_and_path_status_rejected(self):
        for target, field, value in ((self.route, "status", "UNKNOWN"),
                                     (self.path, "status", "PASS"),
                                     (self.path, "failure_kind", "BUDGET_EXHAUSTED"),
                                     (self.route, "acceptance_allowed", True),
                                     (self.route, "force_evaluations", 1)):
            original = target[field]
            target[field] = value
            self.assertFalse(self.allowed())
            target[field] = original

    def test_missing_or_malformed_data_fail_closed(self):
        for value in (None, {}, [], "unknown", {"status": "FAIL_CPU_EVENT_BRANCH"}):
            self.assertFalse(retry.eligible(value, self.route, self.path, **self.extra))
            self.assertFalse(retry.eligible(self.report, value, self.path, **self.extra))
            self.assertFalse(retry.eligible(self.report, self.route, value, **self.extra))
            for key in self.extra:
                damaged = dict(self.extra, **{key: value})
                self.assertFalse(retry.eligible(self.report, self.route, self.path, **damaged))

    def test_union_must_match_all_opening_edges_and_brackets(self):
        union = self.extra["union_route"]
        original = deepcopy(union)
        for field, value in (("transition_brackets", union["transition_brackets"][:1]),
                             ("unresolved_samples", [{"unknown": True}]),
                             ("status", "UNKNOWN")):
            union[field] = value
            self.assertFalse(self.allowed())
            union.update(deepcopy(original))
        union["transition_brackets"][0]["direction"] = "VISIBLE_TO_BLOCKED"
        self.assertFalse(self.allowed())

    def test_shell_and_geometry_original_gates_must_pass(self):
        for key, field, value in (("shell", "status", "UNKNOWN"),
                                  ("shell", "hidden_contacts", [{"hidden": True}]),
                                  ("shell", "finite_patch_unresolved_intervals", [{"unresolved": True}]),
                                  ("geometry", "transient_contacts_missed_by_endpoints", ["hidden"]),
                                  ("geometry", "new_non_bore_features", ["OTHER"])):
            target = self.extra[key]
            original = target[field]
            target[field] = value
            self.assertFalse(self.allowed())
            target[field] = original


if __name__ == "__main__":
    unittest.main()

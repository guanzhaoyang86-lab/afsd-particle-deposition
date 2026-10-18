"""Narrow administrative permission to retry a rejected pair at half its dt.

This recognizes multiple separately observed blocked-to-visible edges in the
frozen CPU driver's declared proposal path. The diagnostic brackets are neither
physical event times nor proof of uniqueness. Returning True only permits the
caller to discard the *entire* failed pair and solve a smaller coarse/fine pair
from its unchanged accepted parent. It never accepts a failed candidate.

The caller owns source/artifact hashes, supervision, parent identity, finite
retry count, minimum dt, budgets, and all original physical acceptance gates.
No physical solver or filesystem access is used here.
"""
import math


MAXIMUM_OBSERVED_OPENINGS = 4
_ERROR = "GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')"
_CAUSE_LINE = "ValueError: only one observed activation per event solve"
_FINAL_LINE = ("run_afsd_deposition_stage_remote_20261007.GeometryPathStop: "
               "geometry/shell/visibility path failed; separate event resolution required")


def _require(condition):
    if not condition:
        raise ValueError("not an eligible separated multiple-opening diagnostic")


def _json_finite(value):
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        _require(math.isfinite(value))
        return
    if isinstance(value, list):
        for item in value:
            _json_finite(item)
        return
    if isinstance(value, dict):
        for key, item in value.items():
            _require(isinstance(key, str))
            _json_finite(item)
        return
    raise ValueError("input is not finite JSON data")


def _number(value, minimum=None, maximum=None, positive=False):
    _require(not isinstance(value, bool) and isinstance(value, (int, float)))
    _require(math.isfinite(value))
    _require(not positive or value > 0)
    _require(minimum is None or value >= minimum)
    _require(maximum is None or value <= maximum)
    return value


def _integer(value):
    _require(type(value) is int and value >= 0)
    return value


def _edge(item):
    _require(isinstance(item, dict))
    left = _integer(item["particle_id_i"])
    right = _integer(item["particle_id_j"])
    _require(left != right)
    return tuple(sorted((left, right)))


def _edge_record(item):
    key = _edge(item)
    i, j = _integer(item["i"]), _integer(item["j"])
    _require(i != j)
    index = _integer(item["edge_index"])
    # Keep the orientation in this record: all audits use one supplied topology.
    return key, (index, i, j, item["particle_id_i"], item["particle_id_j"])


def _unique_edges(items):
    _require(isinstance(items, list))
    result = {}
    for item in items:
        key, record = _edge_record(item)
        _require(key not in result)
        result[key] = record
    _require(len({record[0] for record in result.values()}) == len(result))
    return result


def _validate_route(route, *, union=False):
    _require(route["classification"] == "OBSERVED_VISIBLE_SUPPORT_REMOTE_DOMAIN_EVENT_ROUTER_20261007")
    _require(route["status"] == "VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE")
    _require(route["requires_resolve"] is True and route["acceptance_allowed"] is False)
    _require(route["unresolved_samples"] == [])
    for name in ("physical_solver_evaluations", "nnps_evaluations", "force_evaluations",
                 "material_updates", "edge_deletions"):
        _require(type(route[name]) is int and route[name] == 0)
    brackets = route["transition_brackets"]
    _require(isinstance(brackets, list) and 2 <= len(brackets) <= MAXIMUM_OBSERVED_OPENINGS)
    edges = _unique_edges(brackets)
    _require(_unique_edges(route["changed_edges"]) == edges)
    _require(_unique_edges(route["endpoint_changed_edges"]) == edges)
    start_count = _integer(route["start_blocked_edge_count"])
    end_count = _integer(route["end_blocked_edge_count"])
    _require(start_count - end_count == len(edges))
    _require(union or end_count == 0)
    tolerance = _number(route["fraction_tolerance"], positive=True, maximum=1.0)
    intervals = []
    for item in brackets:
        _require(item["status"] == "BRACKET_LOCALIZED_NOT_UNIQUE_OR_PHYSICALLY_RESOLVED")
        _require(item["direction"] == "BLOCKED_TO_VISIBLE")
        _require(item["start_blocked"] is True and item["end_blocked"] is False)
        _require(item["single_event_certified"] is False and item["physical_substeps_solved"] is False)
        lo = _number(item["start_fraction"], positive=True, maximum=1.0)
        hi = _number(item["end_fraction"], positive=True, maximum=1.0)
        width = _number(item["bracket_width_fraction"], positive=True, maximum=1.0)
        _require(lo < hi)
        roundoff = 4.0 * math.ulp(max(lo, hi))
        _require(abs(width - (hi - lo)) <= roundoff and width <= tolerance + roundoff)
        intervals.append((lo, hi))
    intervals.sort()
    _require(all(left[1] < right[0] for left, right in zip(intervals, intervals[1:])))
    for item in route["endpoint_changed_edges"]:
        _require(item["start_blocked"] is True and item["end_blocked"] is False)
        _number(item["start_strict_intersection_length_m"], positive=True)
        _require(_number(item["end_strict_intersection_length_m"], minimum=0) == 0)
    return edges


def _eligible(report, route, path, union_route, shell, geometry):
    for item in (report, route, path, union_route, shell, geometry):
        _require(isinstance(item, dict))
        _json_finite(item)
    _require(report["status"] == "FAIL_CPU_EVENT_BRANCH")
    _require(report["error"] == _ERROR)
    trace = report["traceback"]
    _require(isinstance(trace, str))
    lines = trace.splitlines()
    _require(lines.count(_CAUSE_LINE) == 1 and lines[-1] == _FINAL_LINE)
    for name in ("accepted", "restart_allowed", "screened_continuation_authorized"):
        _require(report[name] is False)
    _require(report["diagnostic_only"] is True)
    _require(type(report["solves"]) is int and report["solves"] == 1)

    edges = _validate_route(route)
    union_edges = _validate_route(union_route, union=True)
    _require(set(edges) == set(union_edges))
    # The union graph uses different edge indices, but must describe the same
    # particle pairs, point ordering, opening brackets, and endpoint evidence.
    for name in ("transition_brackets", "changed_edges", "endpoint_changed_edges"):
        retained = {_edge(item): {k: v for k, v in item.items() if k != "edge_index"}
                    for item in route[name]}
        expanded = {_edge(item): {k: v for k, v in item.items() if k != "edge_index"}
                    for item in union_route[name]}
        _require(retained == expanded)
    _require(shell["status"] == "CERTIFIED_NO_HIDDEN_SHELL_CONTACT")
    _require(shell["failure_kind"] is None)
    for name in ("hidden_contacts", "unresolved_intervals", "finite_patch_unresolved_intervals"):
        _require(shell[name] == [])
    _require(geometry["transient_contacts_missed_by_endpoints"] == [])
    features = geometry["new_non_bore_features"]
    _require(isinstance(features, list) and all(isinstance(item, str) for item in features))
    _require(set(features) <= {"SHOULDER", "TOP"})

    _require(path["status"] == "SOLID_OCCLUSION_FOUND_ON_RETAINED_EDGE_PATH")
    _require(path["failure_kind"] == "STRICT_SOLID_INTERSECTION_ON_FIXED_RETAINED_PATH")
    _require(path["classification"] == "FIXED_RETAINED_GRAPH_CURVATURE_PRISM_SAT_20261010_V1")
    for name in ("degenerate_edges", "numerical_domain_failures", "unresolved_intervals"):
        _require(path[name] == [])
    for name in ("physical_solver_evaluations", "nnps_evaluations", "edge_deletions"):
        _require(type(path[name]) is int and path[name] == 0)
    found = path["found_edges"]
    _require(isinstance(found, list) and len(found) >= len(edges))
    found_keys = set()
    for item in found:
        key, record = _edge_record(item)
        _require(key in edges and record == edges[key])
        _require(_number(item["fraction"], minimum=0.0, maximum=1.0) == 0.0)
        longest = _number(item["maximum_contiguous_inside_length_m"], positive=True)
        total = _number(item["total_inside_length_m"], positive=True)
        _require(total >= longest)
        strict = item["strict_intervals"]
        _require(isinstance(strict, list) and len(strict) > 0)
        for interval in strict:
            _require(isinstance(interval, list) and len(interval) == 2)
            lo = _number(interval[0], minimum=0.0, maximum=1.0)
            hi = _number(interval[1], minimum=0.0, maximum=1.0)
            _require(lo < hi)
        found_keys.add(key)
    _require(found_keys == set(edges))
    return True


def eligible(report, route, path_audit, *, union_route, shell, geometry):
    """Fail closed for missing, malformed, unknown, or nonfinite diagnostics."""
    try:
        return _eligible(report, route, path_audit, union_route, shell, geometry)
    except (KeyError, TypeError, ValueError, IndexError, OverflowError, RecursionError):
        return False

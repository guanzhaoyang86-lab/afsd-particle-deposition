"""Empirical scheduling costs for complete, independently verified CPU branches.

The caller must first verify each branch and its supervision receipt. This pure
module neither verifies physics nor grants continuation permission. Measurements
include the entire supervised branch wall time, including startup, event proofs,
and output. Every complete branch with 1..12 solves contributes, including
no-event branches. This is an empirical forecast, not a deterministic bound;
the independent 540/600 second child limits remain mandatory.
"""
from collections.abc import Mapping
import math
from numbers import Real


SCHEMA = "AFSD_CPU_BRANCH_COST_20261015_V1"
_POLICY = {
    "schema": SCHEMA,
    "maximum_branch_solves": 12,
    "forecast_limit_s": 450.0,
    "fixed_overhead_s": 30.0,
    "safety_factor": 1.5,
    "single_solve_floor_multiplier": 2.0,
    "internal_limit_s": 540.0,
    "external_limit_s": 600.0,
}
_HISTORY_FIELDS = {
    "maximum_projected_branch_s", "maximum_single_solve_s", "sample_count"
}


def _nonnegative(value, name):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(name + " must be a finite nonnegative real number")
    try:
        result = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(name + " must be a finite nonnegative real number") from error
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(name + " must be a finite nonnegative real number")
    return result


def _count(value, name, minimum=0, maximum=None):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(name + " must be an integer in the permitted range")
    if maximum is not None and value > maximum:
        raise ValueError(name + " must be an integer in the permitted range")
    return value


def _positive(value, name):
    result = _nonnegative(value, name)
    if result <= 0.0:
        raise ValueError(name + " must be strictly positive")
    return result


def _validate_model(model):
    if not isinstance(model, Mapping):
        raise ValueError("cost model must be a mapping")
    if set(model) != set(_POLICY) | _HISTORY_FIELDS:
        raise ValueError("cost model fields do not match its fixed schema")
    if model["schema"] != SCHEMA:
        raise ValueError("unsupported cost model schema")
    for name, expected in _POLICY.items():
        if name == "schema":
            continue
        if name == "maximum_branch_solves":
            value = _count(model[name], name, minimum=1)
        else:
            value = _nonnegative(model[name], name)
        if value != expected:
            raise ValueError("cost policy must retain fixed " + name)
    count = _count(model["sample_count"], "sample_count")
    projected = _nonnegative(model["maximum_projected_branch_s"],
                             "maximum_projected_branch_s")
    single = _nonnegative(model["maximum_single_solve_s"], "maximum_single_solve_s")
    if count == 0 and (projected != 0.0 or single != 0.0):
        raise ValueError("empty history cannot contain measured costs")
    if count > 0 and (projected <= 0.0 or single <= 0.0):
        raise ValueError("measured history must contain strictly positive costs")
    if projected < single:
        raise ValueError("projected whole-branch cost cannot be below a single solve")
    return count, projected, single


def empty_model():
    """Return a fresh model with the immutable, predeclared scheduling policy."""
    return dict(_POLICY, maximum_projected_branch_s=0.0,
                maximum_single_solve_s=0.0, sample_count=0)


def update(model, wall_seconds, solves, solve_costs):
    """Return new monotone history from one complete verified branch.

    ``solve_costs`` is the list/tuple of elapsed seconds of every physical solve
    in this branch, in any order. The number of costs must equal ``solves``.
    Whole-branch wall time must encompass those serial solve durations. The
    input model and costs are never changed, including when validation fails.
    """
    count, projected, single = _validate_model(model)
    wall = _positive(wall_seconds, "wall_seconds")
    n = _count(solves, "solves", minimum=1, maximum=_POLICY["maximum_branch_solves"])
    if not isinstance(solve_costs, (list, tuple)) or len(solve_costs) != n:
        raise ValueError("solve_costs must contain exactly one duration per solve")
    costs = [_positive(value, "solve_costs") for value in solve_costs]
    try:
        total_solve_s = math.fsum(costs)
    except OverflowError as error:
        raise ValueError("solve duration total must be finite") from error
    tolerance = max(1e-6, wall * 1e-9)
    if total_solve_s > wall + tolerance:
        raise ValueError("supervised wall time must include all physical solve durations")
    sample_projected = _nonnegative(wall * (_POLICY["maximum_branch_solves"] / n),
                                    "projected branch duration")
    result = dict(model)
    result.update(maximum_projected_branch_s=max(projected, sample_projected),
                  maximum_single_solve_s=max(single, max(costs)),
                  sample_count=count + 1)
    # Reject numeric overflow now, not only when this state is later scheduled.
    forecast(result)
    return result


def forecast(model):
    """Return the transparent empirical estimate and the unchanged 450 s gate.

    Having an eligible estimate does not supersede source/receipt validation,
    user pause, campaign deadline, physical gates, or child hard time limits.
    """
    count, projected, single = _validate_model(model)
    single_floor = _nonnegative(_POLICY["single_solve_floor_multiplier"] * single,
                                "single-solve floor duration")
    basis = max(projected, single_floor)
    predicted = _nonnegative(_POLICY["fixed_overhead_s"] +
                             _POLICY["safety_factor"] * basis, "forecast duration")
    eligible = count > 0 and predicted <= _POLICY["forecast_limit_s"]
    reason = ("NO_COMPLETE_BRANCH_SAMPLES" if count == 0 else
              "WITHIN_FORECAST_LIMIT" if eligible else "BRANCH_COST_REASSESSMENT_REQUIRED")
    return {
        "schema": SCHEMA,
        "estimate_kind": "empirical_not_deterministic_bound",
        "sample_count": count,
        "maximum_projected_branch_s": projected,
        "maximum_single_solve_s": single,
        "single_solve_floor_s": single_floor,
        "basis_s": basis,
        "fixed_overhead_s": _POLICY["fixed_overhead_s"],
        "safety_factor": _POLICY["safety_factor"],
        "maximum_branch_solves": _POLICY["maximum_branch_solves"],
        "forecast_s": predicted,
        "forecast_limit_s": _POLICY["forecast_limit_s"],
        "internal_limit_s": _POLICY["internal_limit_s"],
        "external_limit_s": _POLICY["external_limit_s"],
        "eligible": eligible,
        "reason": reason,
    }

"""Load and validate incomplete AFSD paper-case configurations.

This module deliberately does not fill unresolved values.  A validated case may
still be unsuitable for a physics run; callers must apply a backend-specific
required-parameter gate before launching a solver.
"""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping


PARAMETER_KEYS = frozenset({"value", "unit", "source", "status"})
ALLOWED_STATUSES = frozenset(
    {
        "UNRESOLVED",
        "PAPER_CONFIRMED",
        "REFERENCE_CONFIRMED",
        "DERIVED",
        "ASSUMED",
        "CALIBRATED",
    }
)
CASE_IDS = frozenset({"starved", "copper-core", "overfed"})
SI_PARAMETER_UNITS = frozenset(
    {
        "1",
        "m",
        "s",
        "kg",
        "K",
        "rad",
        "Pa",
        "m/s",
        "rad/s",
        "kg/m^3",
        "W/(m*K)",
        "J/(kg*K)",
        "W/(m^2*K)",
        "1/s",
        "J",
    }
)
EXPECTED_CALCULATION_UNITS = {
    "system": "SI",
    "length": "m",
    "time": "s",
    "mass": "kg",
    "temperature": "K",
    "angle": "rad",
    "stress": "Pa",
}
PHYSICS_RUN_CLASSIFICATIONS = frozenset({"DEBUG", "APPROXIMATE", "PAPER_COMPARISON"})


class ConfigurationError(ValueError):
    """Raised when a configuration violates the on-disk contract."""


def _is_finite_number(value: object) -> bool:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        data = json.load(stream)
    if not isinstance(data, dict):
        raise ConfigurationError(f"{path}: top-level JSON value must be an object")
    return data


def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(dict(base))
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], Mapping)
            and isinstance(value, Mapping)
        ):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _walk_parameters(
    node: Mapping[str, Any], prefix: str = "parameters"
) -> Iterator[tuple[str, Mapping[str, Any]]]:
    if set(node) == PARAMETER_KEYS:
        yield prefix, node
        return
    if not node:
        raise ConfigurationError(f"{prefix}: parameter group may not be empty")
    for key, value in node.items():
        path = f"{prefix}.{key}"
        if not isinstance(value, Mapping):
            raise ConfigurationError(
                f"{path}: expected a parameter object or nested parameter group"
            )
        yield from _walk_parameters(value, path)


def validate_parameter_tree(parameters: Mapping[str, Any]) -> None:
    """Validate parameter provenance, units, and unresolved-value semantics."""

    for path, parameter in _walk_parameters(parameters):
        keys = set(parameter)
        if keys != PARAMETER_KEYS:
            missing = sorted(PARAMETER_KEYS - keys)
            extra = sorted(keys - PARAMETER_KEYS)
            raise ConfigurationError(f"{path}: missing={missing}, extra={extra}")
        unit = parameter["unit"]
        source = parameter["source"]
        status = parameter["status"]
        if not isinstance(unit, str) or not unit.strip():
            raise ConfigurationError(f"{path}.unit: expected a non-empty string")
        if unit not in SI_PARAMETER_UNITS:
            raise ConfigurationError(
                f"{path}.unit: {unit!r} is not an allowed SI calculation unit; "
                "put presentation units under display_units"
            )
        if not isinstance(source, str) or not source.strip():
            raise ConfigurationError(f"{path}.source: expected a non-empty string")
        if status not in ALLOWED_STATUSES:
            raise ConfigurationError(f"{path}.status: unsupported value {status!r}")
        if status == "UNRESOLVED" and parameter["value"] is not None:
            raise ConfigurationError(
                f"{path}: UNRESOLVED parameters must use JSON null, never a fallback value"
            )
        if status != "UNRESOLVED" and parameter["value"] is None:
            raise ConfigurationError(
                f"{path}: null requires status UNRESOLVED and explicit provenance"
            )
        if status != "UNRESOLVED":
            value = parameter["value"]
            if unit != "1" and not _is_finite_number(value):
                raise ConfigurationError(
                    f"{path}.value: resolved dimensional parameters must be finite "
                    "int/float values (bool and strings are forbidden)"
                )
            if unit == "1":
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    if not _is_finite_number(value):
                        raise ConfigurationError(
                            f"{path}.value: resolved numeric values must be finite"
                        )
                elif isinstance(value, str):
                    if not value:
                        raise ConfigurationError(
                            f"{path}.value: dimensionless enum strings may not be empty"
                        )
                elif not isinstance(value, bool):
                    raise ConfigurationError(
                        f"{path}.value: unit '1' accepts only finite numbers, booleans, "
                        "or non-empty enum strings"
                    )


def load_resolved_case(case_path: str | Path) -> dict[str, Any]:
    """Resolve one case's base file and validate the merged configuration."""

    path = Path(case_path).resolve()
    case = _read_json(path)
    if case.get("configuration_role") != "case":
        raise ConfigurationError(f"{path}: expected configuration_role='case'")
    extends = case.get("extends")
    if not isinstance(extends, str) or not extends:
        raise ConfigurationError(f"{path}: case must name a base file in 'extends'")
    base_path = (path.parent / extends).resolve()
    base = _read_json(base_path)
    if base.get("configuration_role") != "base":
        raise ConfigurationError(f"{base_path}: expected configuration_role='base'")

    override = {key: value for key, value in case.items() if key != "extends"}
    merged = _deep_merge(base, override)
    merged["configuration_role"] = "resolved_case"
    merged["resolved_from"] = {
        "base": str(base_path),
        "case": str(path),
    }
    validate_resolved_case(merged)
    return merged


def validate_resolved_case(config: Mapping[str, Any]) -> None:
    if config.get("schema_version") != "0.1.0":
        raise ConfigurationError("unsupported configuration schema_version")
    if config.get("configuration_role") != "resolved_case":
        raise ConfigurationError("configuration is not a resolved case")
    if config.get("case_id") not in CASE_IDS:
        raise ConfigurationError(f"unsupported case_id {config.get('case_id')!r}")
    units = config.get("calculation_units")
    if not isinstance(units, Mapping):
        raise ConfigurationError("calculation_units must be an object")
    for key, expected in EXPECTED_CALCULATION_UNITS.items():
        if units.get(key) != expected:
            raise ConfigurationError(
                f"calculation_units.{key} must be canonical SI unit {expected!r}"
            )
    display_units = config.get("display_units")
    if not isinstance(display_units, Mapping):
        raise ConfigurationError("display_units must be separate from calculation_units")
    parameters = config.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ConfigurationError("parameters must be an object")
    validate_parameter_tree(parameters)


def unresolved_parameters(config: Mapping[str, Any]) -> list[str]:
    parameters = config.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ConfigurationError("parameters must be an object")
    return [
        path
        for path, parameter in _walk_parameters(parameters)
        if parameter["status"] == "UNRESOLVED"
    ]


def _parameter_at_path(
    config: Mapping[str, Any], requested_path: str
) -> tuple[str, Mapping[str, Any]]:
    if not isinstance(requested_path, str) or not requested_path.strip():
        raise ConfigurationError("required parameter paths must be non-empty strings")
    components = requested_path.split(".")
    if components[0] == "parameters":
        components = components[1:]
    if not components or any(not component for component in components):
        raise ConfigurationError(f"invalid required parameter path {requested_path!r}")
    node: Any = config.get("parameters")
    canonical = "parameters"
    for component in components:
        canonical = f"{canonical}.{component}"
        if not isinstance(node, Mapping) or component not in node:
            raise ConfigurationError(f"unknown required parameter path {requested_path!r}")
        node = node[component]
    if not isinstance(node, Mapping) or set(node) != PARAMETER_KEYS:
        raise ConfigurationError(
            f"required path {requested_path!r} names a group, not a parameter"
        )
    return canonical, node


def assert_physics_ready(
    config: Mapping[str, Any],
    required_paths: Iterable[str] | None = None,
    *,
    run_classification: str | None = None,
) -> None:
    """Fail closed when a resolved case is incomplete for a physics launch.

    With no explicit required path set, every unresolved parameter blocks every
    physics classification. DEBUG/APPROXIMATE callers may pass a non-empty,
    backend-reviewed path set to gate only capabilities used by that run.
    PAPER_COMPARISON always checks the full tree and cannot use this narrowing.
    """

    validate_resolved_case(config)
    classification = run_classification or config.get("intended_run_classification")
    if classification not in PHYSICS_RUN_CLASSIFICATIONS:
        raise ConfigurationError(
            f"unsupported physics run classification {classification!r}"
        )

    if classification == "PAPER_COMPARISON" or required_paths is None:
        selected = list(_walk_parameters(config["parameters"]))
    else:
        if isinstance(required_paths, (str, bytes)):
            raise ConfigurationError(
                "required_paths must be a non-empty collection, not one string"
            )
        requested = list(required_paths)
        if not requested:
            raise ConfigurationError(
                "explicit backend required_paths may not be empty"
            )
        selected_by_path = dict(
            _parameter_at_path(config, requested_path) for requested_path in requested
        )
        selected = sorted(selected_by_path.items())

    unresolved = [
        path for path, parameter in selected if parameter["status"] == "UNRESOLVED"
    ]
    if unresolved:
        scope = (
            "full configuration"
            if classification == "PAPER_COMPARISON" or required_paths is None
            else "backend-required paths"
        )
        raise ConfigurationError(
            f"physics run blocked: {scope} contain UNRESOLVED parameters: "
            + ", ".join(unresolved)
        )


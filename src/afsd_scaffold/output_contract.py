"""Small dependency-free checks for AFSD run manifests."""

from __future__ import annotations

from datetime import datetime
import re
from typing import Any, Mapping


RUN_CLASSIFICATIONS = frozenset(
    {"GEOMETRY_FIXTURE", "DEBUG", "APPROXIMATE", "PAPER_COMPARISON"}
)
RUN_STATUSES = frozenset({"COMPLETE", "INCOMPLETE", "FAILED", "BLOCKED"})
PHYSICS_RUN_CLASSIFICATIONS = RUN_CLASSIFICATIONS - {"GEOMETRY_FIXTURE"}
OUTPUT_FORMATS = frozenset({"VTK_LEGACY_ASCII", "VTU", "HDF5"})
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_PATTERN = re.compile(r"^[0-9a-f]{40,64}$")
REQUIRED_PARTICLE_FIELDS = frozenset(
    {
        "particle_id",
        "initial_position_m",
        "material_id",
        "tracer_type",
        "component_id",
    }
)


class ManifestError(ValueError):
    """Raised when output metadata is incomplete or self-contradictory."""


def _require_mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = parent.get(key)
    if not isinstance(value, Mapping):
        raise ManifestError(f"{key}: expected an object")
    return value


def _require_keys(mapping: Mapping[str, Any], keys: set[str], prefix: str) -> None:
    missing = keys - set(mapping)
    if missing:
        raise ManifestError(f"{prefix}: missing required keys {sorted(missing)}")


def _is_nonempty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _is_string_list(value: object) -> bool:
    return isinstance(value, list) and all(_is_nonempty_string(item) for item in value)


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    """Validate invariants needed even when a JSON Schema package is absent."""

    if manifest.get("schema_version") != "0.1.0":
        raise ManifestError("unsupported manifest schema_version")
    run_id = manifest.get("run_id")
    if not isinstance(run_id, str) or not run_id:
        raise ManifestError("run_id must be a non-empty string")
    classification = manifest.get("run_classification")
    if classification not in RUN_CLASSIFICATIONS:
        raise ManifestError(f"unsupported run_classification {classification!r}")
    if manifest.get("status") not in RUN_STATUSES:
        raise ManifestError("unsupported run status")
    if not isinstance(manifest.get("physics_claims_allowed"), bool):
        raise ManifestError("physics_claims_allowed must be boolean")
    if classification == "GEOMETRY_FIXTURE" and manifest["physics_claims_allowed"]:
        raise ManifestError("GEOMETRY_FIXTURE cannot permit physics claims")

    created = manifest.get("created_utc")
    if not isinstance(created, str):
        raise ManifestError("created_utc must be an ISO-8601 string")
    try:
        datetime.fromisoformat(created.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ManifestError("created_utc is not valid ISO-8601") from exc

    backend = _require_mapping(manifest, "backend")
    _require_keys(
        backend,
        {"name", "repository", "commit", "dirty", "compiler", "build_mode"},
        "backend",
    )
    if not _is_nonempty_string(backend.get("name")):
        raise ManifestError("backend.name must be a non-empty string")
    configuration = _require_mapping(manifest, "configuration")
    _require_keys(configuration, {"case_id", "path", "sha256"}, "configuration")
    if classification in PHYSICS_RUN_CLASSIFICATIONS:
        for key in ("repository", "compiler", "build_mode"):
            if not _is_nonempty_string(backend.get(key)):
                raise ManifestError(f"backend.{key} is required for a physics run")
        commit = backend.get("commit")
        if not isinstance(commit, str) or not GIT_COMMIT_PATTERN.fullmatch(commit):
            raise ManifestError(
                "backend.commit must be a pinned 40-64 digit lowercase hex commit for a physics run"
            )
        if not isinstance(backend.get("dirty"), bool):
            raise ManifestError("backend.dirty must be boolean for a physics run")
        for key in ("case_id", "path"):
            if not _is_nonempty_string(configuration.get(key)):
                raise ManifestError(f"configuration.{key} is required for a physics run")
        digest = configuration.get("sha256")
        if not isinstance(digest, str) or not SHA256_PATTERN.fullmatch(digest):
            raise ManifestError(
                "configuration.sha256 must be a 64 digit lowercase hex digest for a physics run"
            )
    units = _require_mapping(manifest, "units")
    expected_units = {
        "calculation_system": "SI",
        "length": "m",
        "time": "s",
        "mass": "kg",
        "temperature": "K",
        "stress": "Pa",
    }
    for key, expected in expected_units.items():
        if units.get(key) != expected:
            raise ManifestError(f"units.{key} must be {expected!r}")
    frames = _require_mapping(manifest, "coordinate_frames")
    if not frames.get("stored_frame") or not frames.get("definition"):
        raise ManifestError("coordinate frame name and definition are required")
    particles = _require_mapping(manifest, "particles")
    _require_keys(
        particles,
        {
            "file",
            "format",
            "sha256",
            "count",
            "persistent_id_field",
            "initial_position_field",
            "material_id_field",
            "tracer_type_field",
            "fields",
        },
        "particles",
    )
    if not _is_nonempty_string(particles.get("file")):
        raise ManifestError("particles.file must be a non-empty string")
    if particles.get("format") not in OUTPUT_FORMATS:
        raise ManifestError("particles.format is unsupported")
    particle_digest = particles.get("sha256")
    if not isinstance(particle_digest, str) or not SHA256_PATTERN.fullmatch(
        particle_digest
    ):
        raise ManifestError("particles.sha256 must be a 64 digit lowercase hex digest")
    if not isinstance(particles.get("count"), int) or particles["count"] < 0:
        raise ManifestError("particles.count must be a non-negative integer")
    fixed_names = {
        "persistent_id_field": "particle_id",
        "initial_position_field": "initial_position_m",
        "material_id_field": "material_id",
        "tracer_type_field": "tracer_type",
    }
    for key, expected in fixed_names.items():
        if particles.get(key) != expected:
            raise ManifestError(f"particles.{key} must be {expected!r}")
    fields = particles.get("fields")
    if not isinstance(fields, list) or not all(isinstance(item, str) for item in fields):
        raise ManifestError("particles.fields must be a list of field names")
    missing = REQUIRED_PARTICLE_FIELDS - set(fields)
    if missing:
        raise ManifestError(f"particles.fields missing {sorted(missing)}")
    provenance = _require_mapping(manifest, "provenance")
    _require_keys(
        provenance,
        {"command", "logs", "source_files", "approximations"},
        "provenance",
    )
    if not _is_nonempty_string(provenance.get("command")):
        raise ManifestError("provenance.command must be a non-empty string")
    for key in ("logs", "source_files", "approximations"):
        if not _is_string_list(provenance.get(key)):
            raise ManifestError(f"provenance.{key} must be a list of non-empty strings")
    if classification in PHYSICS_RUN_CLASSIFICATIONS:
        if not provenance["logs"]:
            raise ManifestError("provenance.logs must name at least one log for a physics run")
        if not provenance["source_files"]:
            raise ManifestError(
                "provenance.source_files must name reviewed source paths for a physics run"
            )


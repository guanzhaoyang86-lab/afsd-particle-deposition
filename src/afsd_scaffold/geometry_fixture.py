"""Synthetic 3D fixture for IDs, passive labels, and VTK I/O.

The coordinates are intentionally small and schematic.  This is not an AFSD
initial condition, a solver result, or evidence that a backend is three
dimensional.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from .output_contract import validate_manifest
from .vtk_legacy import Vector3, write_point_cloud


MATERIAL_AA6061 = 1
MATERIAL_RIGID_TOOL = 2
MATERIAL_RIGID_BASE = 3

TRACER_NONE = 0
TRACER_OXIDE_SURFACE_PROXY = 1
TRACER_COPPER_CORE_PROXY = 2

COMPONENT_SUBSTRATE = 1
COMPONENT_FEEDSTOCK = 2
COMPONENT_TOOL = 3
COMPONENT_BASE = 4


@dataclass(frozen=True)
class FixtureParticle:
    particle_id: int
    position_m: Vector3
    material_id: int
    tracer_type: int
    component_id: int


def build_geometry_fixture() -> list[FixtureParticle]:
    """Build deterministic, synthetic topology with two passive tracer types."""

    particles: list[FixtureParticle] = []

    # Two planes make the fixture genuinely 3D as data, without claiming a 3D
    # stress/neighbor solve.  Dimensions are schematic and not paper geometry.
    particle_id = 1000
    for z in (-0.001, 0.0):
        for y_index in range(-2, 3):
            for x_index in range(-3, 4):
                particles.append(
                    FixtureParticle(
                        particle_id,
                        (0.002 * x_index, 0.002 * y_index, z),
                        MATERIAL_AA6061,
                        TRACER_NONE,
                        COMPONENT_SUBSTRATE,
                    )
                )
                particle_id += 1

    particle_id = 2000
    for z_index in range(4):
        for y_index in range(-1, 2):
            for x_index in range(-1, 2):
                is_center = x_index == 0 and y_index == 0
                tracer = (
                    TRACER_COPPER_CORE_PROXY
                    if is_center
                    else TRACER_OXIDE_SURFACE_PROXY
                )
                particles.append(
                    FixtureParticle(
                        particle_id,
                        (0.0015 * x_index, 0.0015 * y_index, 0.0015 * (z_index + 1)),
                        MATERIAL_AA6061,
                        tracer,
                        COMPONENT_FEEDSTOCK,
                    )
                )
                particle_id += 1

    tool_xy = (
        (-0.005, -0.005),
        (0.0, -0.005),
        (0.005, -0.005),
        (0.005, 0.0),
        (0.005, 0.005),
        (0.0, 0.005),
        (-0.005, 0.005),
        (-0.005, 0.0),
    )
    particles.extend(
        FixtureParticle(
            3000 + index,
            (x, y, 0.0045),
            MATERIAL_RIGID_TOOL,
            TRACER_NONE,
            COMPONENT_TOOL,
        )
        for index, (x, y) in enumerate(tool_xy)
    )

    base_xy = ((-0.006, -0.004), (0.006, -0.004), (0.006, 0.004), (-0.006, 0.004))
    particles.extend(
        FixtureParticle(
            4000 + index,
            (x, y, -0.002),
            MATERIAL_RIGID_BASE,
            TRACER_NONE,
            COMPONENT_BASE,
        )
        for index, (x, y) in enumerate(base_xy)
    )
    return particles


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def generate_geometry_fixture(output_dir: str | Path) -> tuple[Path, Path]:
    """Write a VTK fixture and its self-identifying manifest."""

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    vtk_path = directory / "AFSD_GEOMETRY_FIXTURE.vtk"
    manifest_path = directory / "AFSD_GEOMETRY_FIXTURE.manifest.json"
    particles = build_geometry_fixture()
    points = [particle.position_m for particle in particles]
    scalar_fields = {
        "particle_id": [particle.particle_id for particle in particles],
        "material_id": [particle.material_id for particle in particles],
        "tracer_type": [particle.tracer_type for particle in particles],
        "component_id": [particle.component_id for particle in particles],
    }
    vector_fields = {"initial_position_m": points}
    write_point_cloud(
        vtk_path,
        points,
        scalar_fields,
        vector_fields,
        title="AFSD GEOMETRY_FIXTURE - SYNTHETIC - NO PHYSICS RESULT",
    )

    manifest = {
        "schema_version": "0.1.0",
        "run_id": "geometry-fixture-v1",
        "run_classification": "GEOMETRY_FIXTURE",
        "physics_claims_allowed": False,
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "status": "COMPLETE",
        "backend": {
            "name": "afsd_scaffold.geometry_fixture",
            "repository": None,
            "commit": None,
            "dirty": None,
            "compiler": None,
            "build_mode": None,
        },
        "configuration": {"case_id": None, "path": None, "sha256": None},
        "units": {
            "calculation_system": "SI",
            "length": "m",
            "time": "s",
            "mass": "kg",
            "temperature": "K",
            "stress": "Pa",
        },
        "coordinate_frames": {
            "stored_frame": "synthetic_fixture_frame",
            "definition": "Schematic right-handed Cartesian coordinates; no mapping to the paper figure or AS/RS is claimed.",
        },
        "particles": {
            "file": vtk_path.name,
            "format": "VTK_LEGACY_ASCII",
            "sha256": _sha256(vtk_path),
            "count": len(particles),
            "persistent_id_field": "particle_id",
            "initial_position_field": "initial_position_m",
            "material_id_field": "material_id",
            "tracer_type_field": "tracer_type",
            "fields": [
                "particle_id",
                "initial_position_m",
                "material_id",
                "tracer_type",
                "component_id",
            ],
        },
        "codebook": {
            "material_id": {
                "1": "AA6061_T6_DEFORMABLE_TAG_HOST",
                "2": "RIGID_TOOL_FIXTURE_MARKER",
                "3": "RIGID_BASE_FIXTURE_MARKER",
            },
            "tracer_type": {
                "0": "NONE",
                "1": "OXIDE_SURFACE_PROXY_PASSIVE",
                "2": "COPPER_CORE_PROXY_PASSIVE",
            },
            "component_id": {
                "1": "SUBSTRATE",
                "2": "FEEDSTOCK",
                "3": "TOOL",
                "4": "BASE",
            },
        },
        "provenance": {
            "command": "python -B scripts/generate_geometry_fixture.py",
            "logs": [],
            "source_files": [
                "src/afsd_scaffold/geometry_fixture.py",
                "src/afsd_scaffold/vtk_legacy.py",
            ],
            "approximations": [
                "Synthetic schematic coordinates; no paper dimensions are inferred"
            ],
            "geometry_source": "synthetic topology/serialization fixture",
            "warning": "GEOMETRY_FIXTURE only; not a solver run, paper geometry, or reproduced result",
            "tag_semantics": "Both oxide and copper tracer proxies retain material_id=1 and may not alter physics.",
        },
    }
    validate_manifest(manifest)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    return vtk_path, manifest_path


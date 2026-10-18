from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from afsd_scaffold.geometry_fixture import (
    MATERIAL_AA6061,
    TRACER_COPPER_CORE_PROXY,
    TRACER_OXIDE_SURFACE_PROXY,
    build_geometry_fixture,
    generate_geometry_fixture,
)
from afsd_scaffold.output_contract import validate_manifest
from afsd_scaffold.vtk_legacy import read_point_cloud


class GeometryFixtureTests(unittest.TestCase):
    def test_legacy_vtk_roundtrip_preserves_ids_and_passive_tags(self) -> None:
        expected = build_geometry_fixture()
        with tempfile.TemporaryDirectory() as temporary_directory:
            vtk_path, manifest_path = generate_geometry_fixture(temporary_directory)
            parsed = read_point_cloud(vtk_path)
            scalars = parsed["scalars"]
            vectors = parsed["vectors"]

            expected_ids = [particle.particle_id for particle in expected]
            self.assertEqual(scalars["particle_id"], expected_ids)
            self.assertEqual(len(set(scalars["particle_id"])), len(expected_ids))
            self.assertEqual(
                scalars["tracer_type"], [particle.tracer_type for particle in expected]
            )
            self.assertEqual(
                scalars["material_id"], [particle.material_id for particle in expected]
            )
            self.assertEqual(vectors["initial_position_m"], parsed["points"])
            self.assertIn("GEOMETRY_FIXTURE", parsed["title"])
            self.assertIn("NO PHYSICS RESULT", parsed["title"])

            for material_id, tracer_type in zip(
                scalars["material_id"], scalars["tracer_type"]
            ):
                if tracer_type in {
                    TRACER_OXIDE_SURFACE_PROXY,
                    TRACER_COPPER_CORE_PROXY,
                }:
                    self.assertEqual(material_id, MATERIAL_AA6061)

            with manifest_path.open("r", encoding="utf-8") as stream:
                manifest = json.load(stream)
            validate_manifest(manifest)
            self.assertEqual(manifest["run_classification"], "GEOMETRY_FIXTURE")
            self.assertFalse(manifest["physics_claims_allowed"])
            digest = hashlib.sha256(vtk_path.read_bytes()).hexdigest()
            self.assertEqual(manifest["particles"]["sha256"], digest)

    def test_fixture_contains_two_tracers_but_never_changes_host_material(self) -> None:
        particles = build_geometry_fixture()
        tagged = [p for p in particles if p.tracer_type != 0]
        self.assertGreater(len(tagged), 0)
        self.assertEqual(
            {p.tracer_type for p in tagged},
            {TRACER_OXIDE_SURFACE_PROXY, TRACER_COPPER_CORE_PROXY},
        )
        self.assertEqual({p.material_id for p in tagged}, {MATERIAL_AA6061})


if __name__ == "__main__":
    unittest.main()


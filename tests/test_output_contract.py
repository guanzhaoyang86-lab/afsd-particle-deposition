from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

from afsd_scaffold.geometry_fixture import generate_geometry_fixture
from afsd_scaffold.output_contract import ManifestError, validate_manifest


ROOT = Path(__file__).resolve().parents[1]


def fixture_manifest() -> dict[str, object]:
    with tempfile.TemporaryDirectory() as temporary_directory:
        _, manifest_path = generate_geometry_fixture(temporary_directory)
        with manifest_path.open("r", encoding="utf-8") as stream:
            return json.load(stream)


def valid_debug_manifest() -> dict[str, object]:
    manifest = fixture_manifest()
    manifest["run_id"] = "debug-contract-example"
    manifest["run_classification"] = "DEBUG"
    manifest["backend"].update(
        {
            "repository": "https://example.invalid/pinned-backend.git",
            "commit": "a" * 40,
            "dirty": False,
            "compiler": "test-compiler 1.0",
            "build_mode": "Release",
        }
    )
    manifest["configuration"].update(
        {
            "case_id": "starved",
            "path": "configs/cases/starved.resolved.json",
            "sha256": "b" * 64,
        }
    )
    manifest["provenance"].update(
        {
            "command": "backend --config configs/cases/starved.resolved.json",
            "logs": ["logs/runs/debug-contract-example.log"],
            "source_files": ["third_party/backend/actual_call_path.cpp"],
            "approximations": ["test-only contract example"],
        }
    )
    return manifest


class OutputContractTests(unittest.TestCase):
    def test_schema_documents_are_parseable_json(self) -> None:
        for name in (
            "afsd-case.schema.json",
            "run-manifest.schema.json",
            "particle-record.schema.json",
        ):
            with self.subTest(name=name):
                with (ROOT / "schemas" / name).open("r", encoding="utf-8") as stream:
                    schema = json.load(stream)
                self.assertEqual(schema["$schema"], "https://json-schema.org/draft/2020-12/schema")

    def test_machine_schemas_require_resolved_and_auditable_fields(self) -> None:
        with (ROOT / "schemas" / "afsd-case.schema.json").open(
            "r", encoding="utf-8"
        ) as stream:
            case_schema = json.load(stream)
        self.assertIn("resolved_from", case_schema["properties"])
        self.assertIn("allOf", case_schema)

        with (ROOT / "schemas" / "run-manifest.schema.json").open(
            "r", encoding="utf-8"
        ) as stream:
            run_schema = json.load(stream)
        backend_required = set(run_schema["properties"]["backend"]["required"])
        self.assertTrue(
            {"repository", "commit", "dirty", "compiler", "build_mode"}
            <= backend_required
        )
        self.assertIn("sha256", run_schema["properties"]["particles"]["required"])
        self.assertEqual(
            set(run_schema["properties"]["provenance"]["required"]),
            {"command", "logs", "source_files", "approximations"},
        )
        self.assertIn("allOf", run_schema)

    def test_geometry_fixture_cannot_enable_physics_claims(self) -> None:
        invalid = fixture_manifest()
        invalid["physics_claims_allowed"] = True
        with self.assertRaisesRegex(ManifestError, "cannot permit physics claims"):
            validate_manifest(invalid)

    def test_persistent_identity_fields_are_mandatory(self) -> None:
        manifest = fixture_manifest()
        manifest["particles"]["fields"].remove("particle_id")
        with self.assertRaisesRegex(ManifestError, "particle_id"):
            validate_manifest(manifest)

    def test_particle_artifact_hash_is_mandatory(self) -> None:
        manifest = fixture_manifest()
        manifest["particles"].pop("sha256")
        with self.assertRaisesRegex(ManifestError, "sha256"):
            validate_manifest(manifest)

    def test_valid_physics_manifest_has_pinned_build_and_log_provenance(self) -> None:
        validate_manifest(valid_debug_manifest())

    def test_physics_manifest_fails_closed_on_missing_provenance(self) -> None:
        cases = (
            (("backend", "commit"), None, "backend.commit"),
            (("backend", "compiler"), None, "backend.compiler"),
            (("configuration", "sha256"), None, "configuration.sha256"),
            (("provenance", "logs"), [], "provenance.logs"),
            (("provenance", "source_files"), [], "provenance.source_files"),
        )
        for path, invalid_value, expected_message in cases:
            with self.subTest(path=".".join(path)):
                manifest = copy.deepcopy(valid_debug_manifest())
                manifest[path[0]][path[1]] = invalid_value
                with self.assertRaisesRegex(ManifestError, expected_message):
                    validate_manifest(manifest)


if __name__ == "__main__":
    unittest.main()

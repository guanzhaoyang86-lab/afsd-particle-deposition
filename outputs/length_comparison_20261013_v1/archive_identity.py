"""Read-only checkpoint identity/reference audit; performs no physical solves.

Run with the repository's .envs/pysph/python.exe -B. The only written artifact is
archived_long_identity.json next to this independent analysis script.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from pysph.solver.utils import load


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
NEW = ROOT / "outputs/length_comparison_20261013_v1/long_6p5mm/stage_001"
OLD = ROOT / "outputs/deposition_campaign_20261007/stage_001"
CPU = ROOT / "outputs/deposition_campaign_20261007/stage_013"
GPU = ROOT / "outputs/gpu_acceleration_20261007/campaign_10h_20261008_v1/attempt_0017"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checkpoint(path, expected_sha):
    actual = sha(path)
    assert actual == expected_sha, str(path)
    saved = load(str(path))
    pa, state = saved["arrays"]["solid"], saved["solver_data"]
    identity = {
        "path": str(path.resolve()), "sha256": actual,
        "manifest_sha256_matches": True,
        "time_s": state["time_s"], "step": state["step"],
        "particle_count": len(pa.x),
        "accepted": state["accepted"],
        "restart_allowed": state["restart_allowed"],
        "diagnostic_only": state["diagnostic_only"],
    }
    return pa, state, identity


def equal(a, b):
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        x, y = np.asarray(a), np.asarray(b)
        return x.dtype == y.dtype and x.shape == y.shape and x.tobytes() == y.tobytes()
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(equal(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    return a == b


def diff_maps(a, b):
    common = sorted(a.keys() & b.keys())
    changed = {k: {"new": a[k], "old": b[k]} for k in common if a[k] != b[k]}
    return {
        "new_count": len(a), "old_count": len(b),
        "common_count": len(common), "common_identical_count": len(common) - len(changed),
        "changed_common": changed,
        "new_only": {k: a[k] for k in sorted(a.keys() - b.keys())},
        "old_only": {k: b[k] for k in sorted(b.keys() - a.keys())},
    }


def pair_identity(new_meta, old_meta, path_key, sha_key):
    a, s, new_identity = checkpoint(ROOT / new_meta[path_key], new_meta[sha_key])
    b, t, old_identity = checkpoint(ROOT / old_meta[path_key], old_meta[sha_key])
    assert a.properties.keys() == b.properties.keys()
    attrs = {
        key: equal(a.properties[key].get_npy_array(), b.properties[key].get_npy_array())
        for key in sorted(a.properties)
    }
    state_equal = {key: equal(s[key], t[key]) for key in sorted(s.keys() & t.keys())}
    assert s.keys() == t.keys()
    allowed = {
        "framework_config", "framework_config_sha256", "process_program",
        "domain_inventory_plan", "source_sha256", "implementation_source_sha256",
    }
    changed = [key for key, same in state_equal.items() if not same]
    assert set(changed) <= allowed
    assert len(attrs) == 40 and all(attrs.values())
    return {
        "new": new_identity, "old": old_identity,
        "exact_property_definition": "Same property name, dtype, shape and C-order bytes",
        "property_count": len(attrs), "exact_property_count": sum(attrs.values()),
        "property_exact": attrs, "state_recursive_exact": state_equal,
        "state_difference_keys": changed,
        "all_other_state_items_recursively_exact": True,
    }, s, t


def reference(folder, baseline_sources, kind):
    meta = read(folder / "manifest.json")
    pa, state, identity = checkpoint(
        ROOT / meta["screened_continuation_checkpoint"], meta["screened_continuation_sha256"]
    )
    screen = meta["pair_screen"] if "pair_screen" in meta else meta["pair_screens"][-1]
    screen_path = folder / "pair_screen.json" if "pair_screen" in meta else ROOT / meta["screened_history"][-1]["screen_path"]
    assert equal(read(screen_path), screen)
    gates = {
        "pair_passed": screen["passed"],
        "global_passed": screen["global_screen"]["pass"],
        "local_passed": screen["local_screen"]["passed"],
        "crossing_ids_match": screen["crossing_ids_match"],
    }
    assert all(gates.values())
    supervision = read(folder / "supervision.json")
    assert supervision["exit_code"] == 0
    assert sha(folder / "manifest.json") == supervision["manifest_sha256"]
    assert sha(folder / "run.log") == supervision["log_sha256"]
    artifacts = {
        path: digest for path, digest in meta["artifacts_sha256"].items()
        if "/fine_2/" in path or path.endswith("/pair_screen.json")
    }
    checks = {path: sha(ROOT / path) == digest for path, digest in artifacts.items()}
    assert checks and all(checks.values())
    assert all(meta["source_sha256"].get(k) == v for k, v in baseline_sources.items())
    additional = {k: v for k, v in meta["source_sha256"].items() if k not in baseline_sources}
    rows = screen["local_screen"]["rows"]
    return {
        "kind": kind, "checkpoint": identity,
        "manifest_path": str((folder / "manifest.json").resolve()),
        "manifest_sha256": sha(folder / "manifest.json"),
        "status": meta["status"], "crossed_particle_count": len(state["crossing_ids"]),
        "gate": gates, "global_screen": screen["global_screen"],
        "local_error_maxima": {
            k: max(row[k] for row in rows)
            for k in ("position_error_m", "temperature_error_K", "stress_relative", "plastic_relative")
        },
        "pair_screen_path": str(screen_path.resolve()), "pair_screen_sha256": sha(screen_path),
        "supervision_exit_code": supervision["exit_code"],
        "supervision_manifest_and_log_sha256_match": True,
        "selected_local_artifacts_checked": len(checks),
        "selected_local_artifact_matches": checks,
        "source_count": len(meta["source_sha256"]),
        "baseline_source_count": len(baseline_sources), "baseline_source_map_identical": True,
        "additional_source_count": len(additional),
        "additional_gpu_overlay_source_count": sum(k.startswith(".envs/gpu_overlay") for k in additional),
        "additional_non_overlay_sources": {k: v for k, v in additional.items() if not k.startswith(".envs/gpu_overlay")},
        "use": "EXTERNAL_READ_ONLY_REFERENCE_ONLY_NOT_NEW_CONFIG_RESTART",
        "pure_cpu_bitwise_extension_guaranteed": False,
    }


def main():
    new_meta, old_meta = read(NEW / "manifest.json"), read(OLD / "manifest.json")
    initial, _, _ = pair_identity(new_meta, old_meta, "initial_checkpoint", "initial_sha256")
    endpoint, new_state, old_state = pair_identity(
        new_meta, old_meta, "screened_continuation_checkpoint", "screened_continuation_sha256"
    )
    sources = diff_maps(new_meta["source_sha256"], old_meta["source_sha256"])
    assert sources["common_identical_count"] == 134
    assert not sources["changed_common"] and not sources["old_only"]
    assert set(sources["new_only"]) == {"outputs/flow_diagnostic_20261013_v1/length_controls/long_6p5mm.json"}
    new_program, old_program = new_state["process_program"], old_state["process_program"]
    assert new_program["phases"][:3] == old_program["phases"][:3]
    new_config, old_config = new_state["framework_config"], old_state["framework_config"]
    config_diff = {
        key: {"new": new_config.get(key), "old": old_config.get(key)}
        for key in sorted(new_config.keys() | old_config.keys())
        if not equal(new_config.get(key), old_config.get(key))
    }
    evidence = {
        "schema": "ARCHIVED_LONG_READ_ONLY_IDENTITY_20261013_V1",
        "status": "PASS_READ_ONLY_PREFIX_IDENTITY_AND_REFERENCE_CHECKS",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "readback_script": str(Path(__file__).resolve()), "readback_script_sha256": sha(Path(__file__)),
        "physical_solver_evaluations": 0,
        "existing_artifacts_modified": False,
        "initial_pair": initial, "endpoint_3ms_pair": endpoint,
        "source_map_difference": sources,
        "configuration_difference": config_diff,
        "program_comparison": {
            "first_three_phases_exact": True,
            "common_first_three_phases": new_program["phases"][:3],
            "common_phase_duration_s": sum(p["duration_s"] for p in new_program["phases"][:3]),
            "authorized_reference_scope_end_s": 0.220,
            "new_future_tail": new_program["phases"][3:],
            "old_future_tail": old_program["phases"][3:],
            "new_program_end_s": new_program["end_s"], "old_program_end_s": old_program["end_s"],
        },
        "cpu_reference": reference(CPU, old_meta["source_sha256"], "ARCHIVED_CPU_11P5625_MS"),
        "gpu_reference": reference(GPU, old_meta["source_sha256"], "ARCHIVED_GPU_BRIDGE_20_MS"),
        "limitations": [
            "The old checkpoints retain their old configuration, program and source identities; they are not restart parents for the new configuration.",
            "Bitwise equality is established only at the checked initial and 3 ms endpoints, not for uncomputed later new-config CPU states.",
            "Any comparison to newly solved short-bar states must report differing time-step histories and, for 20 ms, the explicit GPU bridge implementation.",
            "Reference checkpoints remain diagnostic_only, accepted=false and restart_allowed=false; gate PASS is not deposition validation.",
            "This audit reads selected checkpoint, source-map and local gate evidence; it does not re-audit the complete historical dependency chain.",
        ],
    }
    target = OUT / "archived_long_identity.json"
    target.write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": evidence["status"], "output": str(target),
        "initial_exact_properties": initial["exact_property_count"],
        "endpoint_exact_properties": endpoint["exact_property_count"],
        "cpu_reference_ms": evidence["cpu_reference"]["checkpoint"]["time_s"] * 1e3,
        "gpu_reference_ms": evidence["gpu_reference"]["checkpoint"]["time_s"] * 1e3,
        "physical_solver_evaluations": 0}, ensure_ascii=False))


if __name__ == "__main__":
    main()

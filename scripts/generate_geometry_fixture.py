#!/usr/bin/env python3
"""Generate the explicitly non-physical AFSD geometry/label fixture."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from afsd_scaffold.geometry_fixture import generate_geometry_fixture  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs" / "geometry_fixture",
    )
    args = parser.parse_args()
    vtk_path, manifest_path = generate_geometry_fixture(args.output_dir)
    print("classification=GEOMETRY_FIXTURE")
    print("physics_result=false")
    print(f"vtk={vtk_path}")
    print(f"manifest={manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


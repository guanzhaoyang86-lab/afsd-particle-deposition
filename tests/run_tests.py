#!/usr/bin/env python3
"""Dependency-free test entry point used by CI and recorded local runs."""

from __future__ import annotations

import platform
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))


def main() -> int:
    print(f"python={sys.version.split()[0]}")
    print(f"platform={platform.platform()}")
    print(f"repo_root={REPO_ROOT}")
    suite = unittest.defaultTestLoader.discover(str(REPO_ROOT / "tests"), "test_*.py")
    result = unittest.TextTestRunner(verbosity=2, stream=sys.stdout).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())


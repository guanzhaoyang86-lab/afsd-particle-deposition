"""Independent review tests; criteria were frozen before execution."""
import argparse
from pathlib import Path
import sys
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
parser = argparse.ArgumentParser()
parser.add_argument('--pattern', default='test_*review.py')
args = parser.parse_args()
suite = unittest.defaultTestLoader.discover(str(ROOT/'tests/review'), pattern=args.pattern)
result = unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(0 if result.wasSuccessful() and result.testsRun else 1)

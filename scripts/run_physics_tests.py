"""Run one gate or all backend tests without new test dependencies."""
import argparse
import sys
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
parser = argparse.ArgumentParser()
parser.add_argument('--gate', choices=['a', 'b', 'c', 'd', 'all'], default='all')
args = parser.parse_args()
pattern = 'test_*.py' if args.gate == 'all' else f'test_gate_{args.gate}.py'
suite = unittest.defaultTestLoader.discover(str(ROOT/'tests/backend'), pattern=pattern)
result = unittest.TextTestRunner(verbosity=2).run(suite)
sys.exit(0 if result.wasSuccessful() and result.testsRun else 1)

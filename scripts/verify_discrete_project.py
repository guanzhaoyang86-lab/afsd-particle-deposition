"""Bounded current Python regression discovery, preserving all failures."""
from pathlib import Path
import json
import subprocess
import sys
import time
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))


if __name__=='__main__':
    if '--worker' in sys.argv:
        suite=unittest.TestSuite()
        for start in ('tests','tests/review'):
            suite.addTests(unittest.defaultTestLoader.discover(str(ROOT/start),'test_*.py',top_level_dir=str(ROOT/start)))
        result=unittest.TextTestRunner(verbosity=2,stream=sys.stdout).run(suite)
        raise SystemExit(0 if result.wasSuccessful() else 1)
    directory=ROOT/'logs/discrete_contact_process_20261004/full_regression'
    directory.mkdir(parents=True,exist_ok=False)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--worker']
    start=time.perf_counter()
    with (directory/'stdout.log').open('w',encoding='utf-8') as log:
        try: code=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=540).returncode
        except subprocess.TimeoutExpired: code=124
    receipt=dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-start,
        time_limit_s=540,scope='ALL_CURRENT_PYTHON_UNITTEST_DISCOVERY_ROOT_AND_REVIEW',
        physical_validation=False,failures_retained=True)
    (directory/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt));raise SystemExit(code)

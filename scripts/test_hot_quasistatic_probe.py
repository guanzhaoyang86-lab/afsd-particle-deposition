"""Freeze sources and record the actual selected regression command/exit code."""
from pathlib import Path
import subprocess
import sys
import time
from run_hot_quasistatic_probe import ROOT, LOGS, sources, write

if __name__ == '__main__':
    LOGS.mkdir(parents=True, exist_ok=True)
    if (LOGS/'tests.json').exists() or (LOGS/'tests.log').exists():
        raise ValueError('refuse to overwrite evidence')
    frozen = sources()
    cmd = [sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'tests/review',
           '-p', 'test_quasistatic*.py', '-v']
    start = time.perf_counter()
    with (LOGS/'tests.log').open('w', encoding='utf-8') as stream:
        completed = subprocess.run(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=120.)
    if frozen != sources(): raise ValueError('source changed during tests')
    report = dict(command=cmd, exit_code=completed.returncode, source_sha256=frozen,
        wall_seconds=time.perf_counter()-start,
        scope='SELECTED_NUMERICAL_COMPONENT_TESTS_NOT_FULL_REGRESSION_OR_AFSD_VALIDATION')
    write(LOGS/'tests.json', report)
    print((LOGS/'tests.log').read_text(encoding='utf-8'))
    raise SystemExit(completed.returncode)

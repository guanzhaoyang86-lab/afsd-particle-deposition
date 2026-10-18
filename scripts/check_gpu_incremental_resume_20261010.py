"""Execute bounded-grant resume tests and freeze their administrative evidence."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/gpu_acceleration_20261007'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = BASE/'incremental_resume_components_20261010_v1.json'
    log = output.with_suffix('.log')
    if output.exists() or log.exists():
        raise ValueError('new component receipt paths required')
    original_path = BASE/'incremental_components_20261010_v1.json'
    original = json.loads(original_path.read_text(encoding='utf-8'))
    if original['status'] != 'PASS_INCREMENTAL_COMPONENTS' or original['exit_code'] != 0:
        raise ValueError('qualified original incremental components required')
    for name, digest in {**original['source_sha256'], **original['artifacts_sha256']}.items():
        if sha(ROOT/name) != digest:
            raise ValueError('original frozen evidence changed: '+name)
    sources = dict(original['source_sha256'])
    sources[original_path.relative_to(ROOT).as_posix()] = sha(original_path)
    for name in ('scripts/run_gpu_incremental_resume_20261010.py',
            'scripts/afsd_incremental_worker_resume_20261010.py',
            'scripts/check_gpu_incremental_resume_20261010.py',
            'tests/review/test_incremental_grant_resume_20261010.py',
            'tests/review/test_incremental_resume_driver_20261010.py',
            'docs/GPU_INCREMENTAL_RESUME_20261010.md'):
        sources[name] = sha(ROOT/name)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    checks = []; count = 0
    with log.open('x', encoding='utf-8') as stream:
        for script in ('tests/review/test_incremental_grant_resume_20261010.py',
                'tests/review/test_incremental_resume_driver_20261010.py'):
            command = [sys.executable, '-B', script]; began = time.perf_counter()
            result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120.)
            checks.append(dict(command=command, exit_code=result.returncode,
                wall_seconds=time.perf_counter()-began, stdout=result.stdout, stderr=result.stderr))
            stream.write(json.dumps(command)+'\n'+result.stdout+result.stderr+'\n'); stream.flush()
            match = re.search(r'Ran (\d+) tests?', result.stderr)
            count += int(match.group(1)) if match else 0
    for name, digest in sources.items():
        if sha(ROOT/name) != digest:
            raise ValueError('administrative source changed during tests: '+name)
    success = all(row['exit_code'] == 0 for row in checks) and count > 4
    artifacts = dict(original['artifacts_sha256'])
    artifacts[log.relative_to(ROOT).as_posix()] = sha(log)
    receipt = dict(status='PASS_INCREMENTAL_COMPONENTS' if success else 'FAIL_INCREMENTAL_COMPONENTS',
        exit_code=0 if success else 2, source_sha256=sources, artifacts_sha256=artifacts,
        tests=count, checks=checks, inherited_original_tests=original['tests'],
        inherited_original_passed=original['passed'], inherited_original_skipped=original['skipped'],
        physical_solver_evaluations=0, scope='BOUND_GRANT_WAIT_WITHIN_ORIGINAL_DEADLINES_ONLY',
        physical_source_changed=False, accepted=False, restart_allowed=False)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(status=receipt['status'], tests=count, exit_code=receipt['exit_code'], receipt=str(output))))
    return receipt['exit_code']


if __name__ == '__main__':
    raise SystemExit(main())

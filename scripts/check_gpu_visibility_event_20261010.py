"""Execute registered event components and bind old failure evidence once.

Run only after all new physical and administrative files are finalized. This
checker writes a new immutable receipt; it never launches a GPU campaign.
"""
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
PHYSICAL_FILES = (
    'scripts/run_gpu_visibility_event_20261010.py', 'scripts/afsd_gpu_event_trial_20261010.py',
    'scripts/afsd_event_graph_20261010.py', 'scripts/afsd_event_locator_20261010.py',
    'scripts/afsd_event_blocked_path_20261010.py',
    'scripts/check_visible_support_operators.py',
    'tests/review/test_event_graph_20261010.py', 'tests/review/test_event_locator_20261010.py',
    'tests/review/test_event_blocked_path_20261010.py',
    'tests/review/test_gpu_visibility_event_20261010.py', 'docs/GPU_VISIBILITY_EVENT_CONTRACT_20261010.md')
ADMINISTRATIVE_FILES = (
    'scripts/run_gpu_event_campaign_20261010.py', 'scripts/afsd_event_worker_20261010.py',
    'scripts/afsd_event_receipt_20261010.py', 'scripts/check_gpu_visibility_event_20261010.py',
    'tests/review/test_event_receipt_20261010.py', 'tests/review/test_gpu_event_campaign_20261010.py')
TEST_FILES = (
    'tests/review/test_event_graph_20261010.py', 'tests/review/test_event_locator_20261010.py',
    'tests/review/test_event_blocked_path_20261010.py',
    'tests/review/test_event_receipt_20261010.py', 'tests/review/test_gpu_visibility_event_20261010.py',
    'tests/review/test_gpu_event_campaign_20261010.py')


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    output = BASE/'visibility_event_components_20261010_v1.json'
    log = output.with_suffix('.log')
    if output.exists() or log.exists():
        raise ValueError('new component receipt paths required')
    original_path = BASE/'ten_hour_components_20261010_v1.json'
    original = json.loads(original_path.read_text(encoding='utf-8'))
    if original['status'] != 'PASS_INCREMENTAL_COMPONENTS' or original['exit_code'] != 0:
        raise ValueError('qualified original ten-hour components required')
    for name, digest in {**original['source_sha256'], **original['artifacts_sha256']}.items():
        if sha(ROOT/name) != digest:
            raise ValueError('original frozen evidence changed: '+name)
    sources = dict(original['source_sha256'])
    sources[original_path.relative_to(ROOT).as_posix()] = sha(original_path)
    for name in (*PHYSICAL_FILES, *ADMINISTRATIVE_FILES):
        digest = sha(ROOT/name)
        if name in sources and sources[name] != digest:
            raise ValueError('new registration would overwrite old source: '+name)
        sources[name] = digest
    physical = {name: sources[name] for name in PHYSICAL_FILES}
    artifacts = dict(original['artifacts_sha256'])
    fixture = BASE/'campaign_10h_20261010_v1/attempt_0044'
    failure = json.loads((fixture/'manifest.json').read_text(encoding='utf-8'))
    supervision = json.loads((fixture/'supervision.json').read_text(encoding='utf-8'))
    if (failure['status'] != 'FAIL_GPU_STAGE_NOT_ACCEPTED' or supervision['exit_code'] == 0
            or not supervision['source_inputs_still_match']
            or sha(fixture/'manifest.json') != supervision['manifest_sha256']
            or sha(fixture/'run.log') != supervision['log_sha256']):
        raise ValueError('exact sealed old44 failure receipt required')
    for name, digest in failure['artifacts_sha256'].items():
        if sha(ROOT/name) != digest:
            raise ValueError('old44 failure artifact changed: '+name)
        artifacts[name] = digest
    for name in ('manifest.json', 'supervision.json', 'run.log'):
        path = fixture/name
        artifacts[path.relative_to(ROOT).as_posix()] = sha(path)
    before = {**sources, **artifacts}
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    checks = []; count = 0
    with log.open('x', encoding='utf-8') as stream:
        for script in TEST_FILES:
            command = [sys.executable, '-B', script]; began = time.perf_counter()
            result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120.)
            checks.append(dict(command=command, exit_code=result.returncode,
                wall_seconds=time.perf_counter()-began, stdout=result.stdout, stderr=result.stderr))
            stream.write(json.dumps(command)+'\n'+result.stdout+result.stderr+'\n'); stream.flush()
            match = re.search(r'Ran (\d+) tests?', result.stderr)
            count += int(match.group(1)) if match else 0
    for name, digest in before.items():
        if sha(ROOT/name) != digest:
            raise ValueError('source/evidence changed during component checks: '+name)
    success = all(row['exit_code'] == 0 for row in checks) and count > 20
    artifacts[log.relative_to(ROOT).as_posix()] = sha(log)
    receipt = dict(status='PASS_INCREMENTAL_COMPONENTS' if success else 'FAIL_INCREMENTAL_COMPONENTS',
        exit_code=0 if success else 2, source_sha256=sources, physical_source_sha256=physical,
        artifacts_sha256=artifacts, tests=count, checks=checks,
        inherited_original_tests=original['tests'], inherited_original_receipt=original_path.relative_to(ROOT).as_posix(),
        physical_source_changed=True, full_physical_solver_runs=0, gpu_campaign_runs=0,
        scope='EXPLICIT_COUPLED_VISIBILITY_EVENT_MODEL_AND_ORIGINAL_LIMITS',
        accepted=False, restart_allowed=False)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(status=receipt['status'], tests=count, exit_code=receipt['exit_code'], receipt=str(output))))
    return receipt['exit_code']


if __name__ == '__main__':
    raise SystemExit(main())

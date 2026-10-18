"""Freeze the swept-prism path proof bridge and sensitivity administration; no GPU launch."""
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
ADMINISTRATIVE_FILES = ('scripts/run_gpu_event_prism_campaign_20261010.py',
    'scripts/afsd_event_prism_worker_20261010.py',
    'scripts/check_gpu_event_prism_20261010.py',
    'tests/review/test_gpu_event_prism_campaign_20261010.py', 'docs/GPU_EVENT_PRISM_CAMPAIGN_20261010.md')
PHYSICAL_FILES = ('scripts/run_gpu_visibility_event_prism_20261010.py',
    'scripts/visible_support_prism_20261010.py', 'tests/review/test_visible_support_prism_20261010.py',
    'tests/review/test_gpu_visibility_event_prism_20261010.py', 'docs/GPU_VISIBILITY_EVENT_PRISM_20261010.md')
TEST_FILES = ('tests/review/test_visible_support_prism_20261010.py',
    'tests/review/test_gpu_visibility_event_prism_20261010.py', 'tests/review/test_gpu_event_prism_campaign_20261010.py')


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def merge(*maps):
    result = {}
    for mapping in maps:
        for key, value in mapping.items():
            if key in result and result[key] != value:
                raise ValueError('conflicting frozen component: '+key)
            result[key] = value
    return result


def verify(mapping):
    for name, digest in mapping.items():
        if sha(ROOT/name) != digest:
            raise ValueError('frozen evidence changed: '+name)


def main():
    output = BASE/'event_prism_components_20261010_v1.json'; log = output.with_suffix('.log')
    if output.exists() or log.exists():
        raise ValueError('new administrative receipt and log paths required')
    inherited = BASE/'event_box_components_20261010_v1.json'
    previous = json.loads(inherited.read_text(encoding='utf-8'))
    if previous['status'] != 'PASS_EVENT_BOX_COMPONENTS' or previous['exit_code'] != 0:
        raise ValueError('executed frozen box v3 components required')
    previous_sha = sha(inherited)
    inherited_closure = merge(previous['source_sha256'], previous['artifacts_sha256'],
        {inherited.relative_to(ROOT).as_posix(): previous_sha})
    verify(inherited_closure)
    sources = merge(inherited_closure, {name: sha(ROOT/name) for name in (*PHYSICAL_FILES, *ADMINISTRATIVE_FILES)})
    failed = BASE/'event_box_campaign_20261010_v1/sensitivity_one'
    failure = json.loads((failed/'manifest.json').read_text(encoding='utf-8'))
    supervision = json.loads((failed/'supervision.json').read_text(encoding='utf-8'))
    if (failure['status'] != 'FAIL_GPU_STAGE_NOT_ACCEPTED' or supervision['exit_code'] == 0
            or not supervision['source_inputs_still_match']
            or sha(failed/'manifest.json') != supervision['manifest_sha256']
            or sha(failed/'run.log') != supervision['log_sha256']):
        raise ValueError('exact original v3 post-event active-path failure receipt required')
    artifacts = merge(failure['artifacts_sha256'], {
        (failed/name).relative_to(ROOT).as_posix(): sha(failed/name) for name in ('manifest.json', 'supervision.json', 'run.log')})
    verify(artifacts)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    checks = []; tests = 0
    with log.open('x', encoding='utf-8') as stream:
        for script in TEST_FILES:
            command = [sys.executable, '-B', script]; began = time.perf_counter()
            result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120.)
            stream.write(json.dumps(command)+'\n'+result.stdout+result.stderr+'\n'); stream.flush()
            count = re.search(r'Ran (\d+) tests?', result.stderr)
            tests += int(count.group(1)) if count else 0
            checks.append(dict(command=command, exit_code=result.returncode, wall_seconds=time.perf_counter()-began,
                               stdout=result.stdout, stderr=result.stderr))
    verify(sources); verify(artifacts)
    success = all(row['exit_code'] == 0 for row in checks) and tests >= 20
    artifacts[log.relative_to(ROOT).as_posix()] = sha(log)
    receipt = dict(status='PASS_EVENT_PRISM_COMPONENTS' if success else 'FAIL_EVENT_PRISM_COMPONENTS',
        exit_code=0 if success else 2, source_sha256=sources,
        physical_source_sha256={name: sources[name] for name in PHYSICAL_FILES},
        artifacts_sha256=artifacts, tests=tests, checks=checks,
        inherited_receipt=inherited.relative_to(ROOT).as_posix(), inherited_receipt_sha256=previous_sha,
        inherited_tests=previous['tests'], inherited_physical_source_sha256=previous['physical_source_sha256'],
        physical_source_changed=True, gpu_campaign_runs=0, full_physical_solver_runs=0,
        failed_v3_stage=failed.relative_to(ROOT).as_posix(),
        runtime_new_baseline_PASS_required=True, scope='EXPLICIT_SWEPT_PRISM_PATH_PROOF_SOURCE_BRIDGE_AND_SAME_PARENT_SENSITIVITY',
        accepted=False, restart_allowed=False)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(status=receipt['status'], tests=tests, receipt=str(output))))
    return receipt['exit_code']


if __name__ == '__main__':
    raise SystemExit(main())

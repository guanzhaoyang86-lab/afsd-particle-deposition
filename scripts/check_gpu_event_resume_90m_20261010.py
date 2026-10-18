"""Freeze a bounded administrative resume of the unchanged v5 GPU physics."""
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
OLD_SESSION = BASE/'event_delay_campaign_20261010_v1'
PARENT = OLD_SESSION/'attempt_0007'
PARENT_SHA256 = 'b1f950f4f5b947ed8d13f2b256f33b7ee7fa9f50b8dbbc5ba2f9b8541d8e9c67'
ADMINISTRATIVE_FILES = (
    'scripts/run_gpu_event_resume_90m_20261010.py',
    'scripts/afsd_event_resume_90m_worker_20261010.py',
    'scripts/check_gpu_event_resume_90m_20261010.py',
    'tests/review/test_gpu_event_resume_90m_20261010.py',
    'docs/GPU_EVENT_RESUME_90M_20261010.md')
TEST_FILES = ('tests/review/test_gpu_event_resume_90m_20261010.py',)


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


def rel(path):
    return path.relative_to(ROOT).as_posix()


def parent_evidence():
    """Component provenance only; runtime independently validates the full chain."""
    meta = json.loads((PARENT/'manifest.json').read_text(encoding='utf-8'))
    supervision = json.loads((PARENT/'supervision.json').read_text(encoding='utf-8'))
    summary = json.loads((OLD_SESSION/'final_summary.json').read_text(encoding='utf-8'))
    final_audit = json.loads((OLD_SESSION/'final_history_audit.json').read_text(encoding='utf-8'))
    if (meta['status'] != 'PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION'
            or meta['screened_continuation_authorized'] is not True
            or meta['screened_continuation_sha256'] != PARENT_SHA256
            or sha(ROOT/meta['screened_continuation_checkpoint']) != PARENT_SHA256
            or supervision['exit_code'] != 0 or not supervision['source_inputs_still_match']
            or sha(PARENT/'manifest.json') != supervision['manifest_sha256']
            or sha(PARENT/'run.log') != supervision['log_sha256']
            or summary['status'] != 'STOP_COST_REASSESSMENT' or summary['exit_code'] != 0
            or summary['finished'] is not True or summary['latest_pass_stage'] != rel(PARENT)
            or summary['checkpoint_sha256'] != PARENT_SHA256
            or final_audit['failed'] is not False or final_audit['final_audits'] != 1
            or final_audit['expected_files'] != final_audit['locked_files']):
        raise ValueError('exact complete old7 parent and normal cost-stop provenance required')
    session_files = ['protocol.json', 'final_summary.json', 'final_history_audit.json',
                     'sensitivity_comparison.json', 'deposition_progress.json']
    session_files += [f'attempt_{i:04d}_receipt.json' for i in range(1, 8)]
    artifacts = merge(meta['artifacts_sha256'],
        {rel(PARENT/name): sha(PARENT/name) for name in ('manifest.json', 'supervision.json', 'run.log')},
        {rel(OLD_SESSION/name): sha(OLD_SESSION/name) for name in session_files})
    return meta, artifacts


def main():
    output = BASE/'event_resume_90m_components_20261010_v1.json'
    log = output.with_suffix('.log')
    if output.exists() or log.exists():
        raise ValueError('new receipt and log paths required')
    inherited = BASE/'event_delay_components_20261010_v1.json'
    previous = json.loads(inherited.read_text(encoding='utf-8'))
    if previous['status'] != 'PASS_EVENT_DELAY_COMPONENTS' or previous['exit_code'] != 0:
        raise ValueError('executed frozen v5 component receipt required')
    inherited_sha = sha(inherited)
    parent, artifacts = parent_evidence()
    closure = merge(previous['source_sha256'], previous['artifacts_sha256'],
                    {rel(inherited): inherited_sha}, parent['source_sha256'])
    verify(closure)
    verify(artifacts)
    sources = merge(closure, {name: sha(ROOT/name) for name in ADMINISTRATIVE_FILES})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', PYTHONPATH='')
    checks = []
    tests = 0
    with log.open('x', encoding='utf-8') as stream:
        for script in TEST_FILES:
            command = [sys.executable, '-B', script]
            began = time.perf_counter()
            result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                                    text=True, timeout=120.)
            stream.write(json.dumps(command)+'\n'+result.stdout+result.stderr+'\n')
            stream.flush()
            count = re.search(r'Ran (\d+) tests?', result.stderr)
            tests += int(count.group(1)) if count else 0
            checks.append(dict(command=command, exit_code=result.returncode,
                wall_seconds=time.perf_counter()-began, stdout=result.stdout, stderr=result.stderr))
    verify(sources)
    verify(artifacts)
    success = bool(checks) and all(row['exit_code'] == 0 for row in checks) and tests >= 10
    artifacts[rel(log)] = sha(log)
    receipt = dict(status='PASS_EVENT_RESUME_90M_COMPONENTS' if success else 'FAIL_EVENT_RESUME_90M_COMPONENTS',
        exit_code=0 if success else 2, source_sha256=sources, artifacts_sha256=artifacts,
        physical_source_sha256=dict(parent['source_sha256']), physical_source_changed=False,
        administrative_files=list(ADMINISTRATIVE_FILES), tests=tests, checks=checks,
        inherited_receipt=rel(inherited), inherited_receipt_sha256=inherited_sha,
        inherited_tests=previous['tests'], parent_directory=rel(PARENT), parent_checkpoint_sha256=PARENT_SHA256,
        old_session_directory=rel(OLD_SESSION), absolute_stop='2026-10-11T01:14:34+08:00',
        registered_wall_limits_s=dict(forecast=1200., physical_soft=1260., worker_hard=1290., controller=1380.),
        scope='ADMINISTRATIVE_ONLY_EXACT_V5_LAST_COMPLETE_PARENT_RESUME',
        gpu_campaign_runs=0, full_physical_solver_runs=0, accepted=False, restart_allowed=False)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=receipt['status'], tests=tests, receipt=str(output))))
    return receipt['exit_code']


if __name__ == '__main__':
    raise SystemExit(main())

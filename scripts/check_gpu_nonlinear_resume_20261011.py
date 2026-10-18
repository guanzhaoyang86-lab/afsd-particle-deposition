"""Freeze administrative nonlinear recovery; no physical solver executions."""
from pathlib import Path
import json
import os
import re
import subprocess
import sys
import time
import check_gpu_event_resume_8h_20261011 as old_check

ROOT, BASE = old_check.ROOT, old_check.BASE
require, sha, read, rel = old_check.require, old_check.sha, old_check.read, old_check.rel
merge, verify, stage_evidence = old_check.merge, old_check.verify, old_check.stage_evidence
OLD_SESSION = BASE/'event_resume_8h_20261011_v1'
PARENT = OLD_SESSION/'attempt_0039'
FAILED = OLD_SESSION/'attempt_0040'
PARENT_SHA256 = 'aedf8bfed9c67f38ec12bf38e1e9ec77e37cbba39914f63064dfedbfb537c474'
ADMINISTRATIVE_FILES = (
    'scripts/run_gpu_nonlinear_resume_20261011.py',
    'scripts/afsd_nonlinear_resume_worker_20261011.py',
    'scripts/check_gpu_nonlinear_resume_20261011.py',
    'tests/review/test_gpu_nonlinear_resume_20261011.py',
    'docs/GPU_NONLINEAR_RESUME_20261011.md',
    'scripts/afsd_nonlinear_refinement_20261011.py',
    'tests/review/test_afsd_nonlinear_refinement_20261011.py')
TEST_FILES = (ADMINISTRATIVE_FILES[3], ADMINISTRATIVE_FILES[6])
FAILED_NUMBERS = frozenset((1, 2, 4, 6, 40))


def parent_evidence():
    import run_gpu_nonlinear_resume_20261011 as runner
    summary = read(OLD_SESSION/'final_summary.json')
    require(summary['finished'] is True and summary['status'] == 'STOPPED_CHILD_FAILURE'
        and summary['exit_code'] == 2 and summary['completed_stage_count'] == 35
        and summary['attempt_count'] == 40 and summary['latest_pass_stage'] == rel(PARENT)
        and summary['failed_attempt'] == rel(FAILED) and summary['checkpoint_sha256'] == PARENT_SHA256
        and summary['physical_time_s'] == .08958959960937454 and summary['screened_steps'] == 599
        and summary['adaptive_state']['status'] == 'STOPPED' and summary['adaptive_state']['pending'] is None,
        'exact stopped nonlinear session and final complete parent required')
    audit_evidence = runner.verify_old_audit()
    names = ['protocol.json', 'final_summary.json', 'final_history_audit.json', 'deposition_progress.json']
    names += [f'attempt_{i:04d}_receipt.json' for i in range(1, 41)]
    names += [f'attempt_{i:04d}_failed_cost.json' for i in sorted(FAILED_NUMBERS)]
    artifacts = {rel(OLD_SESSION/name): sha(OLD_SESSION/name) for name in names}
    initial_parent = BASE/'event_resume_90m_20261010_v1/attempt_0013'
    parent = read(initial_parent/'manifest.json'); parent_dir = initial_parent
    physical_source = dict(parent['source_sha256'])
    require(len(physical_source) == 636, 'original 636-entry physical source required')
    for number in range(1, 41):
        folder = OLD_SESSION/f'attempt_{number:04d}'; passed = number not in FAILED_NUMBERS
        meta, stage_artifacts = stage_evidence(folder, passed=passed)
        protocol = read(folder/'protocol.json'); receipt = read(OLD_SESSION/f'attempt_{number:04d}_receipt.json')
        require(meta['source_sha256'] == physical_source == protocol['source_sha256']
            and meta['parent_implementation_source_sha256'] == physical_source
            and meta.get('implementation_transition') is None and meta['event_tolerance_scale'] == .25,
            'old eight-hour physical implementation changed')
        require(meta['parent_sha256'] == parent['screened_continuation_sha256']
            and meta['parent_checkpoint'] == parent['screened_continuation_checkpoint']
            and meta['parent_gpu_directory'] == rel(parent_dir)
            and protocol['parent_sha256'] == meta['parent_sha256']
            and protocol['dt_s'] == meta['dt_s'], 'old complete-parent chain mismatch')
        require(receipt['attempt'] == number and receipt['stage'] == rel(folder)
            and receipt['parent'] == rel(parent_dir) and receipt['dt_s'] == meta['dt_s']
            and receipt['exit_code'] == (0 if passed else 2), 'old administrative receipt mismatch')
        stage_artifacts = merge(stage_artifacts, {rel(folder/'protocol.json'): sha(folder/'protocol.json')})
        if passed:
            require(receipt['verified'] is True
                and receipt['checkpoint_sha256'] == meta['screened_continuation_sha256'],
                'old verified parent receipt missing')
            parent, parent_dir = meta, folder
        else:
            cost = read(OLD_SESSION/f'attempt_{number:04d}_failed_cost.json')
            require(cost['stage'] == rel(folder) and cost['complete_pair'] is False and cost['usable_parent'] is False
                and cost['supervision_sha256'] == sha(folder/'supervision.json'), 'retained failed resource evidence differs')
            require(all(cost['retained_artifacts_sha256'].get(k) == v for k, v in stage_artifacts.items()),
                'failed cost must bind full stage evidence')
            stage_artifacts = merge(stage_artifacts, cost['retained_artifacts_sha256'])
            if number == 40:
                require(meta['error'] == "RefinableCandidateError('complete residual line search/evaluation budget failed')"
                    and receipt['adaptive_event'] == summary['last_adaptive_event']
                    and receipt['adaptive_event']['action'] == 'STOP', 'original nonlinear stop must remain a failure')
            else:
                require(meta['error'] == "GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')"
                    and receipt['adaptive_event']['action'] == 'REFINE', 'old multi-event recovery changed')
        artifacts = merge(artifacts, stage_artifacts)
    require(parent_dir == PARENT and parent['screened_continuation_sha256'] == PARENT_SHA256
        and parent['screened_end_time_s'] == .08958959960937454 and parent['screened_end_step'] == 599,
        'last complete state differs')
    return parent, artifacts, audit_evidence


def main():
    output = BASE/'nonlinear_resume_components_20261011_v1.json'; log = output.with_suffix('.log')
    require(not output.exists() and not log.exists(), 'exclusive new component receipt/log required')
    inherited = BASE/'event_resume_8h_components_20261011_v1.json'; prior = read(inherited)
    require(prior['status'] == 'PASS_EVENT_RESUME_8H_COMPONENTS' and prior['exit_code'] == 0,
        'executed frozen eight-hour component receipt required')
    inherited_sha = sha(inherited); parent, artifacts, audit_evidence = parent_evidence()
    closure = merge(prior['source_sha256'], prior['artifacts_sha256'],
        {rel(inherited): inherited_sha}, parent['source_sha256'])
    verify(closure); verify(artifacts)
    sources = merge(closure, {name: sha(ROOT/name) for name in ADMINISTRATIVE_FILES})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    checks = []; tests = 0
    with log.open('x', encoding='utf-8') as stream:
        for script in TEST_FILES:
            command = [sys.executable, '-B', script]; began = time.perf_counter()
            result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=180.)
            stream.write(json.dumps(command)+'\n'+result.stdout+result.stderr+'\n'); stream.flush()
            count = re.search(r'Ran (\d+) tests?', result.stderr)
            tests += int(count.group(1)) if count else 0
            checks.append(dict(command=command, exit_code=result.returncode,
                wall_seconds=time.perf_counter()-began, stdout=result.stdout, stderr=result.stderr))
    verify(sources); verify(artifacts)
    success = len(checks) == 2 and all(row['exit_code'] == 0 for row in checks) and tests >= 20
    artifacts[rel(log)] = sha(log)
    receipt = dict(status='PASS_NONLINEAR_RESUME_COMPONENTS' if success else 'FAIL_NONLINEAR_RESUME_COMPONENTS',
        exit_code=0 if success else 2, source_sha256=sources, artifacts_sha256=artifacts,
        physical_source_sha256=dict(parent['source_sha256']), physical_source_changed=False,
        administrative_files=list(ADMINISTRATIVE_FILES), tests=tests, checks=checks,
        inherited_receipt=rel(inherited), inherited_receipt_sha256=inherited_sha, inherited_tests=prior['tests'],
        parent_directory=rel(PARENT), parent_checkpoint_sha256=PARENT_SHA256, failed_stage=rel(FAILED),
        old_session_directory=rel(OLD_SESSION), old_final_audit_evidence=audit_evidence,
        authorized_at='2026-10-11T15:41:47+08:00', absolute_stop='2026-10-11T23:41:47+08:00',
        initial_retry_dt_s=.0002373046875, original_failed_stage_dt_s=.000474609375,
        scope='ADMINISTRATIVE_ONLY_EXACT_V5_COMPLETE_PARENT_NONLINEAR_RETRY',
        gpu_campaign_runs=0, full_physical_solver_runs=0, accepted=False, restart_allowed=False)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(status=receipt['status'], tests=tests, receipt=str(output))))
    return receipt['exit_code']


if __name__ == '__main__': raise SystemExit(main())

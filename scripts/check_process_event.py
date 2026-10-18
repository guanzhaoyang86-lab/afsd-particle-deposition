"""Bounded old/new regressions and immutable receipts for one event-root attempt.

This entry point runs no AFSD trial. Missing new sources leave every suite
NOT_RUN. A qualified numerical receipt retains one exactly identified
continuous-reference accuracy failure; it never says all tests passed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_GLOBAL_FREEZE.json'
EXPECTED_FROZEN_FILES = 3565
OLD_TESTS = ['structural_mobility', 'implicit_thermal', 'stage_anderson',
             'continuation_state', 'stage_preconditioner',
             'stage_continuity_preconditioner', 'stage_structural_root',
             'process_thermal_trial', 'process_global_runner']
NEW_TESTS = ['event_response', 'event_response_discrete',
             'stage_event_preconditioner', 'process_event_runner']
TESTS = OLD_TESTS + NEW_TESTS
RETAINED_METHOD = 'test_two_free_axis_cross_coupling_independent_ode'
RETAINED_TEST_ID = '__main__.ActualEventResponse.'+RETAINED_METHOD
INITIAL_FAILURE_LOG = 'logs/process_event_20261003/event_response_initial.log'
DOCUMENT_UPDATES = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
REQUIRED = {'scripts/run_process_event.py', 'docs/PROCESS_EVENT_ROOT_CONTRACT.md',
            'docs/PROCESS_EVENT_RESPONSE_CONTRACT.md', 'docs/PROCESS_EVENT_SCHUR_CONTRACT.md'}
REQUIRED |= {'tests/review/test_'+name+'.py' for name in TESTS}
REQUIRED.add(INITIAL_FAILURE_LOG)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    pending = path.with_suffix('.pending')
    pending.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    pending.replace(path)


def checked_source(rel):
    path = (ROOT/rel).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise ValueError('unsafe/missing source '+str(rel))
    return path


def verify_freeze():
    frozen = json.loads(PRIOR.read_text(encoding='utf-8'))
    files = frozen['files_sha256']
    if len(files) != EXPECTED_FROZEN_FILES:
        raise ValueError('unexpected PROCESS_GLOBAL_FREEZE file count')
    changed = []
    for rel, expected in files.items():
        path = checked_source(rel)
        if sha(path) != expected:
            if rel not in DOCUMENT_UPDATES:
                raise ValueError('old evidence/source changed '+rel)
            changed.append(rel)
    return len(files), sorted(changed)


def suite_result(name, output, return_code):
    """Parse unittest's result; permit only the preserved B-vs-ODE assertion.

    A failure of availability, a different assertion in the same method,
    another test, an error, skip, timeout or changed test count cannot borrow
    this narrowly scoped numerical-only exception.
    """
    counts = re.findall(r'^Ran (\d+) tests? in ', output, re.M)
    count = int(counts[0]) if len(counts) == 1 else 0
    failures = re.findall(r'^FAIL: (test_\w+) \(([\w.]+)\)\s*$', output, re.M)
    errors = re.findall(r'^ERROR: (test_\w+) \(([\w.]+)\)\s*$', output, re.M)
    failed_summaries = re.findall(r'^FAILED \(([^\r\n]*)\)\s*$', output, re.M)
    summary_counts = {}
    if len(failed_summaries) == 1:
        summary_counts = {key.strip(): int(value) for key, value in
                          re.findall(r'([a-z ]+)=(\d+)', failed_summaries[0])}
    failed = max(len(failures), summary_counts.get('failures', 0))
    errored = max(len(errors), summary_counts.get('errors', 0))
    other = sum(value for key, value in summary_counts.items()
                if key not in ('failures', 'errors'))
    passed = (return_code == 0 and len(counts) == 1 and count > 0 and
              re.search(r'^OK\s*$', output, re.M) is not None and
              not failures and not errors and not failed_summaries)
    known_failure = (
        name == 'event_response' and return_code == 1 and count == 11 and len(counts) == 1
        and failures == [(RETAINED_METHOD, RETAINED_TEST_ID)] and not errors
        and failed_summaries == ['failures=1']
        and "np.testing.assert_allclose(r['B'], expected_B, atol=2e-5, rtol=0)" in output
        and 'Not equal to tolerance rtol=0, atol=2e-05' in output
        and 'Max absolute difference among violations:' in output
        and re.search(r'^OK\s*$', output, re.M) is None)
    # This qualified gate requires reproducing the recorded accuracy failure.
    # A newly passing original suite must be reviewed, not quietly relabelled.
    if name == 'event_response':
        numerical_gate = known_failure
    else:
        numerical_gate = passed
    return dict(status=('FAIL_RETAINED_CONTINUOUS_REFERENCE' if known_failure else
                        'PASS' if passed else 'FAIL'),
        numerical_gate_pass=numerical_gate, all_tests_pass=passed,
        tests_run=count, tests_passing=(count if passed else max(0, count-failed-errored-other)),
        tests_failed=failed, tests_errors=errored, tests_other_nonpassing=other,
        failed_test_ids=[test_id for _, test_id in failures],
        error_test_ids=[test_id for _, test_id in errors],
        accuracy='FAIL_RETAINED' if known_failure else 'NOT_ASSESSED_BY_THIS_SUITE')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    out = ap.parse_args().output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('fresh local output required')
    out.mkdir(parents=True)
    started = time.perf_counter()
    code = 2
    report = dict(status='RUNNING', classification='EVENT_ROOT_PREFLIGHT_NOT_PROCESS_VALIDATION',
        tests_run=0, tests_passing=0, tests_failed=0, tests_errors=0,
        tests_other_nonpassing=0, all_tests_pass=False,
        continuous_reference_accuracy='NOT_RUN', process_accuracy='NOT_VALIDATED',
        project_complete=False, production_authorized=False, accepted_steps=0,
        planned_suites=TESTS, tests=[], suite_timeout_s=60, expected_old_tests=103,
        gate_revision='POST_FAILURE_DISCRETE_DERIVATIVE_GATE_WITH_ORIGINAL_ACCURACY_FAIL_RETAINED')
    try:
        checked, changed = verify_freeze()
        freeze_hash = sha(PRIOR)
        report.update(previous_freeze_sha256=freeze_hash, previous_files_checked=checked,
                      authorized_document_changes=changed)
        missing = sorted(rel for rel in REQUIRED if not (ROOT/rel).is_file())
        if missing:
            report.update(status='NOT_RUN_PREFLIGHT_INCOMPLETE_SOURCES', missing_sources=missing)
        else:
            # Import only after all registered test/contract files exist.  The
            # driver owns its complete input inventory; no source is inferred
            # solely from a filename here.
            from run_process_event import SOURCES
            watched = set(SOURCES) | REQUIRED | {Path(__file__).relative_to(ROOT).as_posix()}
            watched |= {p.relative_to(ROOT).as_posix()
                        for p in (ROOT/'src/afsd_backend/pysph').glob('*.py')}
            missing = sorted(rel for rel in watched if not (ROOT/rel).is_file())
            if missing:
                report.update(status='NOT_RUN_PREFLIGHT_INCOMPLETE_SOURCES', missing_sources=missing)
            else:
                sources = {rel: sha(checked_source(rel)) for rel in sorted(watched)}
                report['source_sha256'] = sources
                initial = checked_source(INITIAL_FAILURE_LOG).read_text(encoding='utf-8')
                initial_result = suite_result('event_response', initial, 1)
                if not initial_result['numerical_gate_pass']:
                    raise ValueError('initial continuous-reference failure evidence does not match registered exception')
                report['retained_failure'] = dict(test_id=RETAINED_TEST_ID,
                    original_log=INITIAL_FAILURE_LOG, original_log_sha256=sha(ROOT/INITIAL_FAILURE_LOG),
                    original_accuracy_tolerance_s_kg=2e-5,
                    original_status='FAIL_RETAINED', override_scope='NUMERICAL_DISCRETE_ROOT_ONLY')
                env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                           OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
                for name in TESTS:
                    cmd = [sys.executable, '-B', 'tests/review/test_'+name+'.py']
                    tick = time.perf_counter()
                    try:
                        completed = subprocess.run(cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, encoding='utf-8',
                            errors='replace', timeout=60)
                        output, rc = completed.stdout, completed.returncode
                    except subprocess.TimeoutExpired as error:
                        output = error.stdout or b''
                        if isinstance(output, bytes):
                            output = output.decode('utf-8', errors='replace')
                        output += '\nTIMEOUT_60S\n'
                        rc = 124
                    (out/(name+'.log')).write_text(output, encoding='utf-8')
                    result = suite_result(name, output, rc)
                    row = dict(name=name, **result, return_code=rc,
                        wall_seconds=time.perf_counter()-tick, command=cmd, log=name+'.log')
                    report['tests'].append(row)
                    for key in ('tests_run', 'tests_passing', 'tests_failed',
                                'tests_errors', 'tests_other_nonpassing'):
                        report[key] += row[key]
                    if name == 'event_response' and row['numerical_gate_pass']:
                        report['continuous_reference_accuracy'] = 'FAIL_RETAINED'
                    write(out/'report.json', report)
                    print(json.dumps(row), flush=True)
                    if not row['numerical_gate_pass']:
                        raise ValueError('preflight failed '+name)
                old_count = sum(row['tests_run'] for row in report['tests'] if row['name'] in OLD_TESTS)
                report['old_tests_run'] = old_count
                report['new_tests_run'] = report['tests_run']-old_count
                if old_count != 103:
                    raise ValueError('old regression inventory no longer totals 103')
                if (report['tests_failed'] != 1 or report['tests_errors'] != 0 or
                        report['tests_other_nonpassing'] != 0 or
                        report['tests_passing']+report['tests_failed'] != report['tests_run'] or
                        report['continuous_reference_accuracy'] != 'FAIL_RETAINED'):
                    raise ValueError('unexpected failure inventory for qualified numerical gate')
                for rel, expected in sources.items():
                    if sha(checked_source(rel)) != expected:
                        raise ValueError('source drift during tests '+rel)
                if sha(PRIOR) != freeze_hash:
                    raise ValueError('freeze manifest drift during tests')
                # Recheck artifacts as well as code; a test must not mutate
                # an earlier FAIL, rejected journal, trajectory or checkpoint.
                checked_after, changed_after = verify_freeze()
                report.update(previous_files_rechecked=checked_after,
                              authorized_document_changes=changed_after,
                              status='NUMERICAL_EVENT_PREFLIGHT_PASS_WITH_CONTINUOUS_REFERENCE_FAIL')
                code = 0
    except ModuleNotFoundError as error:
        report.update(status='NOT_RUN_PREFLIGHT_INCOMPLETE_SOURCES', error=repr(error))
    except Exception as error:
        report.update(status='FAIL_PREFLIGHT', error=repr(error))
    report.update(wall_seconds=time.perf_counter()-started,
                  artifact_sha256={p.name: sha(p) for p in out.glob('*.log')})
    write(out/'report.json', report)
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ('source_sha256', 'tests', 'artifact_sha256')}, allow_nan=False), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())

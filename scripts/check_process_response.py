"""Bounded response-component gate; no saved trial, scout or AFSD advance."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time

from check_process_event import suite_result
from check_process_stop import TESTS as OLD_TESTS, inventory as old_inventory
from probe_hht_curved import ROOT, sha, read, write, require

NEW_TESTS = ['contact_path_response', 'contact_response_stage']
TESTS = OLD_TESTS + NEW_TESTS
EXPECTED_OLD_TESTS = 233
EXPECTED_NEW_TESTS = {'contact_path_response': 7, 'contact_response_stage': 5}
EXPECTED_TOTAL_TESTS = 245
PRIOR = ROOT/'docs/PROCESS_STOP_FREEZE.json'
EXPECTED_FROZEN_FILES = 4516
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
INTERNAL_LIMIT_SECONDS = 55.


def verify_previous():
    old = read(PRIOR)
    require(len(old['files_sha256']) == EXPECTED_FROZEN_FILES, 'old inventory differs')
    changed = []
    for name, digest in old['files_sha256'].items():
        target = (ROOT/name).resolve()
        require(target.is_relative_to(ROOT) and target.is_file(), 'missing old file '+name)
        if sha(target) != digest:
            require(name in ALLOWED, 'old source/evidence changed '+name)
            changed.append(name)
    return dict(files_checked=EXPECTED_FROZEN_FILES,
                authorized_document_changes=sorted(changed), manifest_sha256=sha(PRIOR))


def inventory():
    names = set(old_inventory())
    names.update('tests/review/test_'+name+'.py' for name in NEW_TESTS)
    names.update([
        'src/afsd_backend/pysph/contact_path_response.py',
        'src/afsd_backend/pysph/contact_response_preconditioner.py',
        'src/afsd_backend/pysph/hht_response_stage.py',
        'scripts/check_process_response.py', 'scripts/run_process_response_command.py',
        'scripts/probe_hht_response.py',
        'scripts/check_process_stop.py', 'scripts/check_process_curved.py',
        'scripts/check_process_hht.py', 'scripts/check_process_event.py',
        'docs/PROCESS_RESPONSE_CONTRACT.md', 'docs/PROCESS_RESPONSE_TEST_CONTRACT.md',
        'docs/PROCESS_RESPONSE_LIMITS_CONTRACT.md', 'docs/PROCESS_STOP_FREEZE.json'])
    result = {}
    for name in sorted(names):
        path = (ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and path.is_file(), 'missing/unsafe source '+name)
        result[name] = sha(path)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    out = parser.parse_args().output.resolve()
    require(out.is_relative_to(ROOT) and not out.exists(), 'fresh project output required')
    out.mkdir(parents=True)
    started = time.perf_counter()
    deadline = started+INTERNAL_LIMIT_SECONDS
    report = dict(status='RUNNING', tests=[], planned_suites=TESTS,
        classification='FIXED_QUADRATURE_CONTACT_RESPONSE_COMPONENT_PREFLIGHT_NOT_PROCESS',
        production_authorized=False, new_afsd_steps=0, saved_trial_executed=False,
        scout_executed=False, continuous_cross_zero_jacobian_certified=False,
        old_accuracy_failures='RETAINED_NOT_RECLASSIFIED',
        expected_old_tests=EXPECTED_OLD_TESTS, expected_tests=EXPECTED_TOTAL_TESTS,
        internal_limit_seconds=INTERNAL_LIMIT_SECONDS)
    sources = None
    passed = False
    try:
        require(len(TESTS) == len(set(TESTS)), 'duplicate suite identity')
        report['prior_before'] = verify_previous()
        sources = inventory()
        report['source_sha256'] = sources
        write(out/'report.json', report)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
        for name in TESTS:
            remaining = deadline-time.perf_counter()
            require(remaining > 0, 'internal55s component budget exhausted')
            command = [sys.executable, '-B', 'tests/review/test_'+name+'.py']
            tick = time.perf_counter()
            try:
                result = subprocess.run(command, cwd=ROOT, env=env,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                    encoding='utf-8', errors='replace', timeout=remaining)
                output, code = result.stdout, result.returncode
            except subprocess.TimeoutExpired as error:
                output = error.stdout or b''
                if isinstance(output, bytes):
                    output = output.decode('utf-8', errors='replace')
                output += '\nTIMEOUT_REMAINING_INTERNAL55S_BUDGET\n'
                code = 124
            (out/(name+'.log')).write_text(output, encoding='utf-8')
            row = dict(name=name, **suite_result(name, output, code), command=command,
                       return_code=code, wall_seconds=time.perf_counter()-tick,
                       timeout_seconds=remaining, log=name+'.log')
            report['tests'].append(row)
            write(out/'report.json', report)
            print(json.dumps(row), flush=True)
            require(row['all_tests_pass'], 'component test failed '+name)
        require(sum(row['tests_run'] for row in report['tests']
                    if row['name'] in OLD_TESTS) == EXPECTED_OLD_TESTS,
                'old233 suite inventory differs')
        for name, expected in EXPECTED_NEW_TESTS.items():
            row = next(row for row in report['tests'] if row['name'] == name)
            require(row['tests_run'] == expected, 'new suite count differs '+name)
        require(sum(row['tests_run'] for row in report['tests']) == EXPECTED_TOTAL_TESTS,
                'total245 suite inventory differs')
        require(time.perf_counter() <= deadline, 'internal55s component budget exhausted')
        passed = True
    except Exception as error:
        report['error'] = repr(error)

    # Preserve post-run integrity checks even when a component fails. A failed
    # suite or old/source drift can never borrow another package's old PASS.
    integrity_errors = []
    try:
        report['prior_after'] = verify_previous()
    except Exception as error:
        integrity_errors.append('prior_after: '+repr(error))
    if sources is not None:
        try:
            require(inventory() == sources, 'source drift during preflight')
            report['source_unchanged_after'] = True
        except Exception as error:
            report['source_unchanged_after'] = False
            integrity_errors.append('source_after: '+repr(error))
    if integrity_errors:
        report['integrity_errors'] = integrity_errors
        passed = False
    report.update(
        status=('PASS_RESPONSE_COMPONENT_PREFLIGHT_NOT_PROCESS' if passed
                else 'FAIL_RESPONSE_PREFLIGHT'),
        tests_run=sum(row['tests_run'] for row in report['tests']),
        not_run_suites=[name for name in TESTS
                        if name not in [row['name'] for row in report['tests']]],
        wall_seconds=time.perf_counter()-started,
        artifact_sha256={path.name: sha(path) for path in out.glob('*.log')})
    write(out/'report.json', report)
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('source_sha256', 'artifact_sha256', 'tests')}), flush=True)
    return 0 if passed else 2


if __name__ == '__main__':
    raise SystemExit(main())

"""Actual component tests and immutable receipt for the new process runner."""
from pathlib import Path
import argparse
import json
import sys
import time
import unittest
import afsd_process_framework_20261007 as model

ROOT = model.ROOT
SUITES = (
    ('tests/review', 'test_afsd_process_motion_20261007.py'),
    ('tests/review', 'test_afsd_process_framework_20261007.py'),
    ('scripts', 'test_visible_support_events_20261007.py'),
    ('scripts', 'test_afsd_process_observer_20261007.py'),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = (ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT/'outputs') or output.exists():
        raise ValueError('new receipt inside project outputs required')
    output.parent.mkdir(parents=True, exist_ok=True)
    files = [Path(__file__).resolve(), ROOT/'docs/PROCESS_FRAMEWORK_CONTRACT_20261007.md',
        ROOT/'scripts/visible_support_events_20261007.py', ROOT/'scripts/afsd_process_observer_20261007.py']
    files += [ROOT/folder/name for folder, name in SUITES]
    sources = model.source_hashes()
    sources.update({p.relative_to(ROOT).as_posix(): model.old.expanded.sha(p) for p in files})
    checks = []; outputs = {}; start = time.perf_counter()
    for folder, name in SUITES:
        suite = unittest.TestLoader().discover(str(ROOT/folder), pattern=name, top_level_dir=str(ROOT/folder))
        log = output.with_name(output.stem+'_'+Path(name).stem+'.log')
        with log.open('x', encoding='utf-8') as stream:
            result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
        passed = result.wasSuccessful() and result.testsRun > 0 and not result.skipped
        checks.append(dict(name=name, status='PASS' if passed else 'FAIL', tests_run=result.testsRun,
            failures=len(result.failures), errors=len(result.errors), skipped=len(result.skipped),
            exit_code=0 if passed else 2, log=log.relative_to(ROOT).as_posix()))
        outputs[log.relative_to(ROOT).as_posix()] = model.old.expanded.sha(log)
    reuse = model.operator_reuse_audit()
    for name, expected in sources.items():
        if model.old.expanded.sha(ROOT/name) != expected:
            raise ValueError('source changed while checking '+name)
    passed = all(row['status'] == 'PASS' for row in checks)
    report = dict(status='PASS_AFSD_FRAMEWORK_COMPONENTS_NOT_PROCESS_VALIDATION' if passed else
        'FAIL_AFSD_FRAMEWORK_COMPONENTS', checks=checks, tests_run=sum(r['tests_run'] for r in checks),
        operator_reuse_audit=reuse, source_sha256=sources, output_sha256=outputs,
        wall_seconds=time.perf_counter()-start, physical_process_solved=False,
        accepted_steps=0, accepted_time_advanced_s=0.)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k:report[k] for k in ('status', 'tests_run', 'wall_seconds')}))
    return 0 if passed else 2


if __name__ == '__main__':
    sys.exit(main())

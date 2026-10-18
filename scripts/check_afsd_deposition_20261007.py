"""Execute new deposition components and bind unchanged upstream receipts."""
from pathlib import Path
import argparse
import json
import sys
import time
import unittest
import afsd_deposition_setup_20261007 as model

ROOT = model.ROOT
BASE_RECEIPT = ROOT/'outputs/process_framework_20261007/components_v1.json'
SUITES = (
    ('tests/review', 'test_afsd_deposition_setup_20261007.py'),
    ('scripts', 'test_afsd_deposition_stage_20261007.py'),
    ('scripts', 'test_afsd_deposition_evidence_20261007.py'),
    ('scripts', 'test_afsd_deposition_readback_20261007.py'),
)


def sha(path):
    return model.old.expanded.sha(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    output = (ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT/'outputs') or output.exists():
        raise ValueError('new component receipt inside project outputs required')
    upstream = json.loads(BASE_RECEIPT.read_text(encoding='utf-8'))
    if not upstream['status'].startswith('PASS') or not all(r['status'] == 'PASS' for r in upstream['checks']):
        raise ValueError('unchanged upstream components must have passed')
    evidence = {BASE_RECEIPT.relative_to(ROOT).as_posix(): sha(BASE_RECEIPT)}
    for key in ('source_sha256', 'output_sha256', 'evidence_sha256'):
        for name, digest in upstream.get(key, {}).items():
            if sha(ROOT/name) != digest:
                raise ValueError('upstream receipt changed: '+name)
            evidence[name] = digest
    files = [Path(__file__).resolve(), ROOT/'scripts/run_afsd_deposition_stage_20261007.py',
        ROOT/'scripts/afsd_deposition_evidence_20261007.py',
        ROOT/'scripts/readback_afsd_deposition_20261007.py',
        ROOT/'docs/DEPOSITION_STAGE_RUNNER_CONTRACT_20261007.md']
    files += [ROOT/folder/name for folder, name in SUITES]
    sources = model.source_hashes()
    sources.update({p.relative_to(ROOT).as_posix(): sha(p) for p in files})
    output.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    checks, logs = [], {}
    for folder, name in SUITES:
        suite = unittest.TestLoader().discover(str(ROOT/folder), pattern=name, top_level_dir=str(ROOT/folder))
        log = output.with_name(output.stem+'_'+Path(name).stem+'.log')
        with log.open('x', encoding='utf-8') as stream:
            result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
        passed = result.wasSuccessful() and result.testsRun > 0 and not result.skipped
        checks.append(dict(name=name, status='PASS' if passed else 'FAIL', tests_run=result.testsRun,
            failures=len(result.failures), errors=len(result.errors), skipped=len(result.skipped),
            exit_code=0 if passed else 2, log=log.relative_to(ROOT).as_posix()))
        logs[log.relative_to(ROOT).as_posix()] = sha(log)
    for name, digest in {**sources, **evidence}.items():
        if sha(ROOT/name) != digest:
            raise ValueError('component inputs changed during execution: '+name)
    passed = all(row['status'] == 'PASS' for row in checks)
    report = dict(status='PASS_DEPOSITION_COMPONENTS_NOT_PHYSICAL_VALIDATION' if passed else
        'FAIL_DEPOSITION_COMPONENTS', tests_run=sum(row['tests_run'] for row in checks),
        checks=checks, source_sha256=sources, evidence_sha256=evidence, output_sha256=logs,
        operator_reuse_audit=model.operator_reuse_audit(), upstream_tests_reexecuted=False,
        unchanged_upstream_pass_count=upstream['tests_run'], wall_seconds=time.perf_counter()-started,
        physical_solver_evaluations=0, accepted_time_advanced_s=0.)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k:report[k] for k in ('status', 'tests_run', 'wall_seconds')}))
    return 0 if passed else 2


if __name__ == '__main__':
    sys.exit(main())

"""Archive unchanged physics + measured sensitivity; integrity, not certification."""
import json
from pathlib import Path
from freeze_process_friction_q5 import sha, read, safe_file, write_new
from freeze_process_feasible_predictor import tests

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_PREDICTOR_ONSET_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_PREDICTOR_SENSITIVITY_FREEZE.json'
LOGS = 'logs/process_predictor_sensitivity_20261002'
OUTPUTS = 'outputs/process_predictor_sensitivity_20261002'
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {TARGET.relative_to(ROOT).as_posix(), LOGS+'/freeze_create.log', LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():
        raise FileExistsError('never overwrite a freeze')
    old = read(PRIOR)
    before = read(safe_file(LOGS+'/freeze_before.log'))
    if before['status'] != 'PASS' or before['files_checked'] != len(old['files_sha256']):
        raise ValueError('missing old archive verification')
    changed = [s for s, h in old['files_sha256'].items() if sha(safe_file(s)) != h]
    if set(changed)-ALLOWED:
        raise ValueError('old evidence changed: '+str(changed))
    case = ROOT/OUTPUTS/'n3_i8'
    r = read(case/'report.json')
    rb = read(ROOT/OUTPUTS/'n3_i8_readback.json')
    if rb['completion_status'] != r['status'] or rb['accepted_time_s'] != r['accepted_time_s']:
        raise ValueError('readback mismatch')
    for name, h in r['artifact_sha256'].items():
        if Path(name).name != name or sha(case/name) != h:
            raise ValueError('artifact mismatch')
    comparison_path = ROOT/OUTPUTS/'comparison.json'
    comparison = read(comparison_path) if comparison_path.exists() else None
    if r['status'] == 'J2_CONE_PATH_COMPLETE':
        if comparison is None or comparison['status'] != 'COMPARISON_COMPLETE':
            raise ValueError('complete run missing comparison')
        provenance = comparison['provenance']
        if provenance['contract_sha256'] != sha(ROOT/'docs/PROCESS_PREDICTOR_SENSITIVITY_CONTRACT.md'):
            raise ValueError('measurement contract changed')
        if provenance['comparator_sha256'] != sha(ROOT/'scripts/compare_process_predictor_sensitivity.py'):
            raise ValueError('comparator changed after reading')
    elif r['status'] != 'BLOCKED' or comparison is not None:
        raise ValueError('invalid partial comparison')
    files = {safe_file(s) for s in old['files_sha256']} | {PRIOR}
    for directory, pattern in [('docs', 'PROCESS_PREDICTOR_SENSITIVITY*'),
            ('scripts', '*process_predictor_sensitivity*.py'),
            ('tests/review', 'test_process_predictor_sensitivity*.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS, OUTPUTS):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files = {p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    result = dict(schema_version='1.0', date='2026-10-02',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False, scope='integrity only; NOT scientific acceptance',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed, old_failures_preserved=True),
        actual_tests={name: tests(LOGS+'/'+name+'.log') for name in
                     ('comparator_tests_final', 'integration_compatibility')},
        case=dict(status=r['status'], error=r['error'], accepted_time_s=r['accepted_time_s'],
                  accepted_steps=r['accepted_steps'], events=len(r['events']), readback=rb['status']),
        comparison=None if comparison is None else dict(status=comparison['status'],
            screening_status=comparison['screening_status'], warnings=comparison['warnings']),
        excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=result['preservation'],
                         case=result['case'], comparison=result['comparison']), indent=2))


if __name__ == '__main__':
    main()

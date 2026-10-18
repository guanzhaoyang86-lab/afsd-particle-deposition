"""Freeze the residual-interface/route assessment, not new AFSD acceptance."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha, read, safe_file, write_new
from freeze_process_feasible_predictor import tests

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/COUPLED_AFSD_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_ROUTE_FREEZE.json'
LOGS = 'logs/process_route_20261002'
OUT = 'outputs/process_route_20261002'
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {TARGET.relative_to(ROOT).as_posix(), LOGS+'/freeze_create.log', LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists(): raise FileExistsError('never overwrite evidence')
    old = read(PRIOR); before = read(safe_file(LOGS+'/freeze_before.log'))
    if before['status'] != 'PASS' or before['files_checked'] != len(old['files_sha256']):
        raise ValueError('missing prior integrity check')
    changed = [p for p, h in old['files_sha256'].items() if sha(safe_file(p)) != h]
    if set(changed)-ALLOWED: raise ValueError('prior evidence/source changed '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8')) > 800: raise ValueError('HANDOFF over800')
    report = read(ROOT/OUT/'assessment_complete.json')
    if (report['status'] != 'READ_ONLY_ROUTE_ASSESSMENT_COMPLETE_NOT_NEW_PHYSICAL_SOLUTION'
            or report['accepted_progress'] is not False or report['production_authorized'] is not False):
        raise ValueError('assessment scope/status mismatch')
    for name, h in dict(report['input_sha256'], **report['source_sha256']).items():
        if sha(safe_file(name)) != h: raise ValueError('assessment source/input drift '+name)
    files = {safe_file(p) for p in old['files_sha256']} | {PRIOR}
    for directory, pattern in [('docs', 'PROCESS_ROUTE*'), ('scripts', '*process_route*.py'),
          ('src/afsd_backend/pysph', 'process_residual.py'), ('tests/review', 'test_process_residual.py'),
          ('tests/review', 'test_process_route_readback.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS, OUT): files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files = {p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    logs = sorted((ROOT/LOGS).glob('*tests*.log'))
    result = dict(schema_version='1.0', date='2026-10-02',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False, accepted_physical_progress=False,
        scope='Integrity only; pure residual interface/readonly assessment, not a working global implicit AFSD solver',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed, old_failures_preserved=True),
        actual_test_logs={p.name: tests(p.relative_to(ROOT).as_posix()) for p in logs},
        assessment=dict(status=report['status'], accepted_steps_checked=report['accepted_steps_checked'],
            contact_subintervals_checked=report['contact_subintervals_checked']),
        files_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=result['preservation']), indent=2))


if __name__ == '__main__': main()

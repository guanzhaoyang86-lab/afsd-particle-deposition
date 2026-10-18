"""Archive this bounded integration stage; never scientific acceptance."""
import json
from pathlib import Path
from freeze_process_friction_q5 import sha, read, safe_file, write_new
from freeze_process_feasible_predictor import tests

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/ROTATING_CONTACT_FREEZE.json'
TARGET = ROOT/'docs/ROTATING_DEFORMABLE_FREEZE.json'
LOGS = 'logs/rotating_deformable_20261002'
OUTPUTS = 'outputs/rotating_deformable_20261002'
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {TARGET.relative_to(ROOT).as_posix(), LOGS+'/freeze_create_final.log', LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists(): raise FileExistsError('never overwrite a freeze')
    old, before = read(PRIOR), read(safe_file(LOGS+'/freeze_before.log'))
    if before['status'] != 'PASS' or before['files_checked'] != len(old['files_sha256']):
        raise ValueError('missing prior verification')
    changed = [s for s, h in old['files_sha256'].items() if sha(safe_file(s)) != h]
    if set(changed)-ALLOWED: raise ValueError('previous evidence changed: '+str(changed))
    cases = {}
    for path in sorted((ROOT/OUTPUTS).glob('*/manifest.json')):
        r = read(path)
        for src, h in r['source_sha256'].items():
            if sha(safe_file(src)) != h: raise ValueError('run source changed: '+src)
        for artifact, h in r['artifact_sha256'].items():
            if Path(artifact).name != artifact or sha(path.parent/artifact) != h:
                raise ValueError('run artifact changed: '+artifact)
        cases[path.parent.name] = {k: r.get(k) for k in
            ('status', 'accepted_steps', 'accepted_time_s', 'wall_seconds', 'production_authorized')}
    files = {safe_file(s) for s in old['files_sha256']} | {PRIOR}
    for directory, pattern in [('docs', 'ROTATING_DEFORMABLE*'),
            ('configs/benchmarks', 'rotating_deformable_block.json'),
            ('scripts', '*rotating_deformable*.py'),
            ('tests/review', 'test_rotating_deformable*.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS, OUTPUTS):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files = {p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    test_inventory = {p.name: tests(p.relative_to(ROOT).as_posix())
        for p in sorted((ROOT/LOGS).glob('*tests*.log'))}
    for name in ('old_contact_events.log', 'old_mechanics.log'):
        test_inventory[name] = tests(LOGS+'/'+name)
    result = dict(schema_version='1.0', date='2026-10-02',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False, scope='file integrity only; NOT scientific acceptance',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed, old_failures_preserved=True),
        actual_test_logs=test_inventory, cases=cases,
        excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files),
        preservation=result['preservation'], cases=cases), indent=2))


if __name__ == '__main__': main()

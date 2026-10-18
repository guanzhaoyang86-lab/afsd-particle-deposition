"""Freeze two bounded component units; archive integrity is not physics validation."""
import json
from pathlib import Path
from freeze_process_friction_q5 import sha, read, safe_file, write_new
from freeze_process_feasible_predictor import tests

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_PREDICTOR_SENSITIVITY_FREEZE.json'
TARGET = ROOT/'docs/ROTATING_CONTACT_FREEZE.json'
LOGS = 'logs/rotating_contact_20261002'
OUTPUTS = 'outputs/rotating_contact_20261002'
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {TARGET.relative_to(ROOT).as_posix(), LOGS+'/freeze_create.log', LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():
        raise FileExistsError('never overwrite a freeze')
    old, before = read(PRIOR), read(safe_file(LOGS+'/freeze_before.log'))
    if before['status'] != 'PASS' or before['files_checked'] != len(old['files_sha256']):
        raise ValueError('missing prior verification')
    changed = [s for s, h in old['files_sha256'].items() if sha(safe_file(s)) != h]
    if set(changed)-ALLOWED:
        raise ValueError('previous evidence changed: '+str(changed))
    readback = read(ROOT/OUTPUTS/'readback.json')
    cases = {}
    for name in ('dt2e5', 'dt1e5'):
        path = ROOT/OUTPUTS/name
        r = read(path/'manifest.json')
        if r['status'] != 'COMPLETE_MANUFACTURED_COMPONENT_NOT_AFSD' or r['accepted_time_s'] != .0008:
            raise ValueError('not complete at frozen endpoint')
        for src, expected in r['source_sha256'].items():
            if sha(safe_file(src)) != expected:
                raise ValueError('run source changed: '+src)
        for artifact, expected in r['artifact_sha256'].items():
            if Path(artifact).name != artifact or sha(path/artifact) != expected:
                raise ValueError('run artifact changed')
        cases[name] = dict(status=r['status'], accepted_steps=r['accepted_steps'],
                           time_s=r['accepted_time_s'], inner_wall_seconds=r['wall_seconds'])
    files = {safe_file(s) for s in old['files_sha256']} | {PRIOR}
    for directory, pattern in [('docs', 'ROTATING_CONTACT*'),
            ('configs/benchmarks', 'rotating_shoulder_component.json'),
            ('scripts', '*rotating*contact*.py'), ('scripts', '*rotating_shoulder*.py'),
            ('src/afsd_backend/pysph', 'rigid_contact_kinematics.py'),
            ('tests/review', 'test_rigid_contact_kinematics.py'),
            ('tests/review', 'test_rotating_shoulder*.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS, OUTPUTS):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files = {p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    result = dict(schema_version='1.0', date='2026-10-02', backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'], production_authorized=False,
        scope='integrity only; NOT scientific acceptance',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed, old_failures_preserved=True),
        actual_tests={n: tests(LOGS+'/'+n+'.log') for n in
            ('kinematics_final', 'component_final_gate', 'old_contact_events', 'readback_tests_final')},
        preserved_initial_test_failure=tests(LOGS+'/component_short_tests_initial.log'),
        cases=cases, readback_status=readback['status'],
        excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=result['preservation'],
                         cases=cases, readback_status=readback['status']), indent=2))


if __name__ == '__main__':
    main()

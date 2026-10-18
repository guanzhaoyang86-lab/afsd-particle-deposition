"""Freeze onset evidence. Integrity only, never grants scientific acceptance."""
import json
from pathlib import Path
from freeze_process_friction_cone_q6 import sha, read, safe_file, test_log, write_new

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_PLASTIC_FRICTION_Q7_RECOVERY_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_PLASTIC_FRICTION_ONSET_FREEZE.json'
LOGS = 'logs/process_plastic_friction_onset_20261002'
OUTPUTS = 'outputs/process_plastic_friction_onset_20261002'
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {TARGET.relative_to(ROOT).as_posix(), LOGS+'/freeze_create.log', LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():
        raise FileExistsError('never overwrite a freeze')
    old = read(PRIOR)
    before = read(safe_file(LOGS+'/freeze_before.log'))
    if before['status'] != 'PASS' or before['files_checked'] != len(old['files_sha256']):
        raise ValueError('missing initial integrity check')
    changed = [s for s,h in old['files_sha256'].items() if sha(safe_file(s)) != h]
    if set(changed)-ALLOWED:
        raise ValueError('unauthorized old evidence change: '+str(changed))
    files = {safe_file(s) for s in old['files_sha256']} | {PRIOR}
    for directory,pattern in [('docs','PROCESS_PLASTIC_FRICTION_ONSET*'),
            ('configs/benchmarks','process_plastic_friction_onset*'),
            ('scripts','*process_plastic_friction_onset*'),
            ('tests/review','test_process_plastic_friction_onset*.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS, OUTPUTS):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    cases = {}
    for case in ('n3_i4','n3_i8'):
        folder = ROOT/OUTPUTS/case
        if not (folder/'report.json').is_file():
            marker = read(folder/'NOT_RUN.json')
            if marker['status'] != 'NOT_RUN' or not marker['reason']:
                raise ValueError('missing explicit case outcome')
            cases[case] = marker
            continue
        r = read(folder/'report.json')
        rb = read(ROOT/OUTPUTS/(case+'_readback.json'))
        if r['status'] not in ('BLOCKED','J2_CONE_PATH_COMPLETE') or rb['completion_status'] != r['status']:
            raise ValueError('run/readback status mismatch')
        if rb['accepted_time_s'] != r['accepted_time_s'] or rb['onset_audit']['variant'] != 'Q7_BOUNDARY_ELIGIBILITY_V1':
            raise ValueError('run/readback variant mismatch')
        for name,digest in r['source_sha256'].items():
            if sha(safe_file(name)) != digest: raise ValueError('source changed after run')
        for name,digest in r['artifact_sha256'].items():
            if Path(name).name != name or sha(folder/name) != digest: raise ValueError('artifact mismatch')
        extra = read(folder/'onset_provenance.json')
        for name,digest in extra['source_sha256'].items():
            if sha(safe_file(name)) != digest: raise ValueError('recovery source changed after run')
        cases[case] = dict(status=r['status'], error=r['error'], accepted_steps=r['accepted_steps'],
            accepted_time_s=r['accepted_time_s'], events=len(r['events']),
            readback_status=rb['status'], report_sha256=sha(folder/'report.json'),
            readback_sha256=sha(ROOT/OUTPUTS/(case+'_readback.json')))
    files = {p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    result = dict(schema_version='1.0', date='2026-10-02',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False, scope='integrity only; NOT scientific acceptance',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed, old_failures_preserved=True),
        actual_tests={name:test_log(LOGS+'/'+name+'.log') for name in
            ('onset_tests_final','q7_compatibility','q6_compatibility')},
        cases=cases, excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)})
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),preservation=result['preservation'],
                         cases=cases,production_authorized=False),indent=2))


if __name__ == '__main__':
    main()

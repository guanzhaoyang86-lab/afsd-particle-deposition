"""Seal curved-contact evidence, not completion or physical acceptance."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha,read,safe_file,write_new
from freeze_process_event import check_report

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_SLOW_FREEZE.json'
TARGET=ROOT/'docs/PROCESS_CURVED_FREEZE.json'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={'docs/PROCESS_CURVED_FREEZE.json','logs/process_curved_20261003/freeze_create.log',
    'logs/process_curved_20261003/freeze_final.log'}


def main():
    if TARGET.exists():raise FileExistsError('never overwrite evidence')
    old=read(PRIOR)
    changed=[name for name,h in old['files_sha256'].items() if sha(safe_file(name))!=h]
    if set(changed)-ALLOWED:raise ValueError('previous evidence changed '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('HANDOFF over800')
    preflight=ROOT/'logs/process_curved_20261003/preflight_initial/report.json'
    checks=read(preflight);check_report(preflight.parent,checks)
    if checks['status']!='PASS_CURVED_COMPONENT_PREFLIGHT_NOT_PROCESS' or checks['tests_run']!=189:
        raise ValueError('wrong component receipt')
    reference=ROOT/'logs/process_slow_20261003/reference_repair1/report.json'
    check_report(reference.parent,read(reference))
    pilot=ROOT/'outputs/process_curved_20261003/pilot_initial/report.json'
    scout=ROOT/'outputs/process_curved_20261003/S0_initial/report.json'
    for path in (pilot,scout):check_report(path.parent,read(path))
    if read(scout)['run']['published_steps']!=13 or abs(read(scout)['run']['time_s']-5e-6)>1e-19:
        raise ValueError('unexpected S0 history')
    private=ROOT/'outputs/process_curved_20261003/saved_red_initial/report.json'
    check_report(private.parent,read(private))
    if read(private)['status']!='PASS_PRIVATE_SAVED_RED_REGRESSION_NOT_TRAJECTORY' or read(private)['accepted']:
        raise ValueError('private diagnostic candidate cannot be promoted')
    replay=ROOT/'logs/process_curved_20261003/replay_initial/report.json'
    replay_value=read(replay)
    if replay_value['status']!='PASS_BITWISE_PUBLISHED_INTERVAL_REPLAY_NOT_PROCESS':
        raise ValueError('one-interval replay evidence required')
    if replay_value['actual_advance_calls']!=1 or not replay_value['complete_semantic_hash_equal']:
        raise ValueError('replay did not match original complete state')
    runs=[]
    for base in ('logs','outputs'):
        for folder in (ROOT/base).glob('process_curved*'):
            for path in folder.rglob('report.json'):
                report=read(path);check_report(path.parent,report,False)
                runs.append(dict(path=path.relative_to(ROOT).as_posix(),status=report['status']))
    files={safe_file(name) for name in old['files_sha256']}|{PRIOR,Path(__file__).resolve()}
    files.update(safe_file(name) for name in checks['source_sha256'])
    files.update(safe_file(name) for name in ('logs/process_slow_20261003/freeze_final.log',
        'logs/process_slow_20261003/freeze_final_process_scope.log'))
    files.update(path for path in (ROOT/'docs').glob('PROCESS_CURVED*') if path.is_file())
    for base in ('logs','outputs'):
        for folder in (ROOT/base).glob('process_curved*'):
            files.update(path for path in folder.rglob('*') if path.is_file())
    files={path for path in files if path.relative_to(ROOT).as_posix() not in EXCLUDED}
    value=dict(schema_version='1.0',date='2026-10-03',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],project_complete=False,production_authorized=False,
        scope='FILE_INTEGRITY_ONLY_LOCAL_CURVED_CONTACT_NOT_AFSD_TRANSPORT',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        final_preflight=preflight.relative_to(ROOT).as_posix(),runs=runs,
        new_afsd_steps=0,afsd_accepted_time_s=850.393632910904e-9,
        local_scout_accepted_time_s=5e-6,independent_shared_prefix_end_s=5.9375e-6,
        files_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,value)
    print(json.dumps(dict(status='CREATED',files=len(files),project_complete=False,
        new_afsd_steps=0,preservation=value['preservation']),indent=2))


if __name__=='__main__':main()

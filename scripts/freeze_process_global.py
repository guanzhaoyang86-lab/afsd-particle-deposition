"""Chain immutable joint-root evidence; no process-completion certification."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha,read,safe_file,write_new
ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_COMPLETION_FREEZE.json'
TARGET=ROOT/'docs/PROCESS_GLOBAL_FREEZE.json'
LOGS='logs/process_global_20261003';OUT='outputs/process_global_20261003'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():raise FileExistsError('never overwrite frozen evidence')
    old=read(PRIOR);changed=[rel for rel,h in old['files_sha256'].items() if sha(safe_file(rel))!=h]
    if set(changed)-ALLOWED:raise ValueError('old evidence/source drift '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('handoff over800')
    preflight=read(ROOT/LOGS/'preflight_initial/report.json')
    attempt=read(ROOT/OUT/'root_attempt/manifest.json')
    if preflight['status']!='PASS_GLOBAL_PREFLIGHT_NOT_PROCESS_VALIDATION' or preflight['project_complete']:
        raise ValueError('preflight missing/wrong scope')
    if attempt['status']=='RUNNING_UNCOMMITTED' or attempt['production_authorized']:
        raise ValueError('unfinished or wrong-scope actual attempt')
    for report in (preflight,read(ROOT/OUT/'cost_probe/manifest.json'),attempt):
        for rel,h in dict(report.get('input_sha256',{}),**report['source_sha256']).items():
            if sha(safe_file(rel))!=h:raise ValueError('actual source/input changed '+rel)
    for directory,report in ((ROOT/LOGS/'preflight_initial',preflight),(ROOT/OUT/'root_attempt',attempt),
                            (ROOT/OUT/'cost_probe',read(ROOT/OUT/'cost_probe/manifest.json'))):
        for rel,h in report['artifact_sha256'].items():
            p=(directory/rel).resolve()
            if not p.is_relative_to(directory) or sha(p)!=h:raise ValueError('artifact drift '+rel)
    files={safe_file(rel) for rel in old['files_sha256']}|{PRIOR}
    files.update(safe_file(rel) for rel in preflight['source_sha256'])
    files.update(p for p in (ROOT/'docs').glob('PROCESS_GLOBAL*') if p.is_file())
    for name in ('freeze_process_global.py','assess_global_step_scales.py',
                 'assess_process_global_heat.py','assess_process_global_nonlinearity.py'):
        p=ROOT/'scripts'/name
        if p.exists():files.add(p)
    for directory in (LOGS,OUT,'outputs/process_global_step_audit_20261003','logs/process_global_step_audit_20261003'):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    result=dict(schema_version='1.0',date='2026-10-03',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],project_complete=False,production_authorized=False,
        scope='File integrity only; larger-step numerical roots do not certify process transport',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        actual_attempt_status=attempt['status'],new_coupled_steps=attempt['accepted_steps'],
        accepted_time_s=attempt['accepted_time_s'],preflight=LOGS+'/preflight_initial/report.json',
        actual_attempt=OUT+'/root_attempt/manifest.json',
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),project_complete=False,
        new_coupled_steps=result['new_coupled_steps'],preservation=result['preservation']),indent=2))


if __name__=='__main__':main()

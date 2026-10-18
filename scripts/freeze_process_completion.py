"""Preserve process-component evidence, never certify process completion."""
from pathlib import Path
import json
import re
from freeze_process_friction_q5 import sha,read,safe_file,write_new
ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/NONLINEAR_CONTINUATION_FREEZE.json'
TARGET=ROOT/'docs/PROCESS_COMPLETION_FREEZE.json'
LOGS='logs/process_completion_20261003';OUT='outputs/process_completion_20261003'
ALLOWED={'README.md','HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():raise FileExistsError('never overwrite frozen evidence')
    old=read(PRIOR);before=read(ROOT/LOGS/'freeze_before.log')
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):raise ValueError('initial freeze missing')
    changed=[rel for rel,h in old['files_sha256'].items() if sha(safe_file(rel))!=h]
    if set(changed)-ALLOWED:raise ValueError('old physical/evidence drift '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('handoff over800')
    preflight=read(ROOT/LOGS/'preflight/report.json')
    if (preflight['status']!='PASS_COMPONENT_PREFLIGHT_PROCESS_STILL_BLOCKED' or
            preflight['project_complete'] or preflight['new_coupled_steps']!=0 or
            preflight['total_distinct_component_tests']!=61):raise ValueError('preflight scope/count changed')
    reports=[preflight,read(ROOT/OUT/'thermal_replay.json'),read(ROOT/OUT/'response_probe/report.json')]
    for report in reports:
        if report['accepted_progress'] or report['production_authorized'] or report['new_coupled_steps']!=0:
            raise ValueError('component cannot be accepted process progress')
        for rel,h in dict(report.get('input_sha256',{}),**report['source_sha256']).items():
            if sha(safe_file(rel))!=h:raise ValueError('source/input changed '+rel)
    for rel,h in preflight['artifact_sha256'].items():
        if sha(ROOT/LOGS/'preflight'/rel)!=h:raise ValueError('preflight test log changed')
    files={safe_file(rel) for rel in old['files_sha256']}|{PRIOR,ROOT/'README.md'}
    for pattern in ('PROCESS_COMPLETION*','PROCESS_FINITE_TOOL*','PROCESS_INLET*','PROCESS_IMPLICIT_THERMAL*',
                    'PROCESS_STRUCTURAL_MOBILITY*','PROCESS_THERMAL_INDEPENDENT*','PROCESS_GEOMETRY_INDEPENDENT*',
                    'PROCESS_RESPONSE*'):
        files.update(p for p in (ROOT/'docs').glob(pattern) if p.is_file())
    for module in ('finite_tool_geometry','inlet_protocol','implicit_thermal','structural_mobility'):
        files.add(ROOT/'src/afsd_backend/pysph'/(module+'.py'))
        files.add(ROOT/'tests/review'/('test_'+module+'.py'))
    for script in ('check_process_readiness.py','replay_process_thermal.py','probe_process_response.py','freeze_process_completion.py'):
        files.add(ROOT/'scripts'/script)
    for directory in (LOGS,OUT):files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    logs={}
    for path in sorted((ROOT/LOGS).glob('*tests*.log')):
        txt=path.read_text(encoding='utf-8-sig');runs=re.findall(r'^Ran (\d+) tests? in ',txt,re.M)
        final=re.findall(r'^(OK|FAILED \([^\r\n]+\))$',txt,re.M)
        if len(runs)!=1 or len(final)!=1:raise ValueError('missing unittest terminal '+path.name)
        # A failing unittest may contain multiple failed subtests in one method.
        # Do not subtract subtest failures from the number of test methods.
        logs[path.name]=dict(tests_run=int(runs[0]),terminal=final[0],
            status='PASS' if final[0]=='OK' else 'FAIL_PRESERVED',sha256=sha(path))
    result=dict(schema_version='1.0',date='2026-10-03',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],project_complete=False,production_authorized=False,new_coupled_steps=0,
        scope='File integrity only; component PASS and offline replay do not certify macro transport',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),actual_component_logs=logs,
        preflight=LOGS+'/preflight/report.json',offline_reports=[OUT+'/thermal_replay.json',OUT+'/response_probe/report.json'],
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),project_complete=False,preservation=result['preservation']),indent=2))


if __name__=='__main__':main()

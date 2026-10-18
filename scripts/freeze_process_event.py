"""Preserve numerical event-root evidence, including every retained failure."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha, read, safe_file, write_new

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_GLOBAL_FREEZE.json'
TARGET=ROOT/'docs/PROCESS_EVENT_FREEZE.json'
LOGS='logs/process_event_20261003'; OUT='outputs/process_event_20261003'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(), LOGS+'/freeze_create.log', LOGS+'/freeze_final.log'}


def check_report(directory, report, check_sources=True):
    if report.get('status') in ('RUNNING','RUNNING_UNCOMMITTED'):
        raise ValueError('unfinished report '+str(directory))
    if report.get('project_complete') or report.get('production_authorized'):
        raise ValueError('wrong completion scope')
    if check_sources:
        for rel,h in dict(report.get('input_sha256',{}), **report['source_sha256']).items():
            if sha(safe_file(rel))!=h:raise ValueError('source/input drift '+rel)
    for rel,h in report.get('artifact_sha256',{}).items():
        p=(directory/rel).resolve()
        if not p.is_relative_to(directory) or sha(p)!=h:raise ValueError('artifact drift '+rel)


def main():
    if TARGET.exists():raise FileExistsError('never overwrite frozen evidence')
    old=read(PRIOR)
    changed=[rel for rel,h in old['files_sha256'].items() if sha(safe_file(rel))!=h]
    if set(changed)-ALLOWED:raise ValueError('old source/evidence drift '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('HANDOFF over800')
    preflight_path=ROOT/LOGS/'preflight_refined/report.json'
    final=read(preflight_path)
    if final['status']!='FAIL_PREFLIGHT' or final.get('error')!="ValueError('preflight failed event_response_refined')":
        raise ValueError('expected final failed gate required; never relabel as launch permission')
    if final['tests_failed']!=2 or final['tests_run']!=146 or final['all_tests_pass']:
        raise ValueError('original accuracy and new bounded-refinement failures not retained')
    check_report(preflight_path.parent,final)
    check_report(ROOT/LOGS/'preflight_stable',read(ROOT/LOGS/'preflight_stable/report.json'))
    # The initial checker caught its own mid-run source change. Preserve that
    # invalid receipt and logs; never certify its historical source inventory.
    initial=read(ROOT/LOGS/'preflight_initial/report.json')
    if initial['status']!='FAIL_PREFLIGHT' or 'source drift during tests' not in initial['error']:
        raise ValueError('unexpected initial orchestration evidence')
    check_report(ROOT/LOGS/'preflight_initial',initial,False)
    runs=[]
    for path in sorted((ROOT/OUT).glob('*/manifest.json')):
        report=read(path);check_report(path.parent,report)
        runs.append(dict(path=path.relative_to(ROOT).as_posix(),status=report['status'],
            mode=report['mode'],accepted_steps=report['accepted_steps'],
            accepted_time_s=report['accepted_time_s']))
    if not runs:raise ValueError('actual probe evidence required')
    files={safe_file(rel) for rel in old['files_sha256']}|{PRIOR}
    files.update(safe_file(rel) for rel in final['source_sha256'])
    files.update(p for p in (ROOT/'docs').glob('PROCESS_EVENT*') if p.is_file())
    for name in ('freeze_process_event.py','diagnose_event_response_pid1.py','assess_process_event.py'):
        path=ROOT/'scripts'/name
        if path.is_file():files.add(path)
    for directory in (LOGS,OUT):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    new_steps=sum(r['accepted_steps'] for r in runs)
    data=dict(schema_version='1.0',date='2026-10-03',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],project_complete=False,production_authorized=False,
        scope='File integrity only; event derivative/root availability is not process or transport validation',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        preflight=preflight_path.relative_to(ROOT).as_posix(),tests_passing=final['tests_passing'],
        tests_failed=final['tests_failed'],continuous_reference_accuracy='FAIL_RETAINED',
        invalid_historical_preflight=LOGS+'/preflight_initial/report.json',
        runs=runs,new_coupled_steps=new_steps,accepted_time_s=max(r['accepted_time_s'] for r in runs),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,data)
    print(json.dumps(dict(status='CREATED',files=len(files),new_coupled_steps=new_steps,
        project_complete=False,preservation=data['preservation']),indent=2))


if __name__=='__main__':main()

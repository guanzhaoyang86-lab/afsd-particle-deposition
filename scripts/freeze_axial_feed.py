"""Preserve local axial-feed evidence and all prior FAILs; not a science certificate."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha,read,safe_file,write_new
from freeze_process_feasible_predictor import tests
ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/NORMAL_CONTACT_FREEZE.json';TARGET=ROOT/'docs/AXIAL_FEED_FREEZE.json'
LOGS='logs/axial_feed_20261002';OUTPUTS='outputs/axial_feed_20261002'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():raise FileExistsError('never overwrite evidence freeze')
    old=read(PRIOR);before=read(safe_file(LOGS+'/freeze_before.log'))
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):raise ValueError('missing previous verification')
    changed=[s for s,h in old['files_sha256'].items() if sha(safe_file(s))!=h]
    if set(changed)-ALLOWED:raise ValueError('unauthorized old evidence change '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('HANDOFF over800 chars')
    cases={}
    for path in sorted((ROOT/OUTPUTS).glob('*/manifest.json')):
        m=read(path)
        for s,h in m['source_sha256'].items():
            if sha(safe_file(s))!=h:raise ValueError('run source drift '+s)
        for name,h in m['artifact_sha256'].items():
            if Path(name).name!=name or sha(path.parent/name)!=h:raise ValueError('artifact drift '+name)
        if m['production_authorized'] is not False:raise ValueError('unintended production permission')
        cases[path.parent.name]={k:m.get(k) for k in ('status','accepted_steps','accepted_time_s','wall_seconds','production_authorized')}
    reference=read(ROOT/OUTPUTS/'actuator_reference.json')
    for s,h in reference['source_sha256'].items():
        if sha(safe_file(s))!=h:raise ValueError('reference source drift '+s)
    files={safe_file(s) for s in old['files_sha256']}|{PRIOR}
    for directory,pattern in [('docs','AXIAL_FEED*'),('configs/benchmarks','axial_feed.json'),
        ('scripts','*axial*.py'),('tests/review','test_axial*.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    files.add(ROOT/'src/afsd_backend/pysph/axial_actuator.py')
    for directory in (LOGS,OUTPUTS):files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    test_paths=sorted({p for pattern in ('*unit*.log','old_*tests.log') for p in (ROOT/LOGS).glob(pattern)})
    inventory={p.name:tests(p.relative_to(ROOT).as_posix()) for p in test_paths}
    final_new=sum(inventory[name]['passing_tests'] for name in
        ('actuator_unit_initial.log','block_unit_source_final.log','readback_unit_complete.log'))
    result=dict(schema_version='1.0',date='2026-10-02',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],production_authorized=False,scope='file integrity only; NOT scientific acceptance',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        actual_test_logs=inventory,unique_final_new_short_passes=final_new,cases=cases,
        excluded_self_or_future_paths=sorted(EXCLUDED),files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)})
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),preservation=result['preservation'],
        unique_final_new_short_passes=final_new,cases=cases),indent=2))


if __name__=='__main__':main()

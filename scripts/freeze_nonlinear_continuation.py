"""Append-only evidence freeze for global continuation, not paper validation."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha,read,safe_file,write_new
from freeze_process_feasible_predictor import tests

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/UPDATED_TRIAL_FREEZE.json'
TARGET=ROOT/'docs/NONLINEAR_CONTINUATION_FREEZE.json'
LOGS='logs/nonlinear_continuation_20261003'
OUT='outputs/nonlinear_continuation_20261003'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():raise FileExistsError('never overwrite frozen evidence')
    old=read(PRIOR);before=read(ROOT/LOGS/'freeze_before.log')
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):raise ValueError('initial freeze verification missing')
    changed=[rel for rel,h in old['files_sha256'].items() if sha(safe_file(rel))!=h]
    if set(changed)-ALLOWED:raise ValueError('old physics/evidence changed '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('HANDOFF over800')
    runs=[]
    for path in sorted((ROOT/OUT).glob('*/manifest.json')):
        m=read(path)
        if m['status']=='RUNNING_UNCOMMITTED' or m['production_authorized'] is not False:raise ValueError('invalid run scope/status')
        for rel,h in dict(m['input_sha256'],**m['source_sha256']).items():
            if sha(safe_file(rel))!=h:raise ValueError('source/input drift '+rel)
        for name,h in m['artifact_sha256'].items():
            item=(path.parent/name).resolve()
            if not item.is_relative_to(path.parent) or sha(item)!=h:raise ValueError('artifact drift '+name)
        runs.append(dict(path=path.relative_to(ROOT).as_posix(),status=m['status'],accepted_steps=m['accepted_steps']))
    report=read(ROOT/OUT/'independent_readback.json')
    if report['reader_sha256']!=sha(ROOT/'scripts/assess_continuation_result.py'):
        raise ValueError('independent reader changed after readback')
    cost=read(ROOT/OUT/'process_feasibility_corrected.json')
    if cost['source_sha256']!=sha(ROOT/'scripts/assess_nonlinear_continuation_cost.py'):
        raise ValueError('corrected feasibility reader drift')
    for rel,h in cost['input_sha256'].items():
        if sha(safe_file(rel))!=h:raise ValueError('corrected feasibility input drift '+rel)
    files={safe_file(p) for p in old['files_sha256']}|{PRIOR}
    for directory,pattern in [('docs','NONLINEAR_CONTINUATION*'),('scripts','*nonlinear_continuation*.py'),
        ('scripts','assess_continuation_result.py'),('src/afsd_backend/pysph','stage_anderson.py'),
        ('src/afsd_backend/pysph','continuation_state.py'),('src/afsd_backend/pysph','transport_observer.py'),
        ('tests/review','test_stage_anderson.py'),('tests/review','test_transport_observer.py'),
        ('tests/review','test_continuation_state.py'),('tests/review','test_nonlinear_continuation_runner.py'),
        ('tests/review','test_nonlinear_feasibility.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS,OUT,'outputs/nonlinear_continuation_audit_20261003'):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    result=dict(schema_version='1.0',date='2026-10-03',backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        production_authorized=False,scope='Integrity only; short accepted continuation does not certify transport/deposition',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        actual_test_logs={p.name:tests(p.relative_to(ROOT).as_posix()) for p in sorted((ROOT/LOGS).glob('*tests*.log'))
            if not p.name.endswith(('.stdout.log','.stderr.log'))},
        runs=runs,independent_readback_status=report['status'],
        feasibility_report=OUT+'/process_feasibility_corrected.json',
        invalid_preserved_feasibility_report=OUT+'/process_feasibility.json',
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,result);print(json.dumps(dict(status='CREATED',files=len(files),runs=runs,preservation=result['preservation']),indent=2))


if __name__=='__main__':main()

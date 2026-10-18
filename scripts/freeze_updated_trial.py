"""Append-only evidence freeze, not a scientific acceptance certificate."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha, read, safe_file, write_new
from freeze_process_feasible_predictor import tests

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_ROUTE_FREEZE.json'
TARGET=ROOT/'docs/UPDATED_TRIAL_FREEZE.json'
LOGS='logs/updated_trial_20261003'
OUT='outputs/updated_trial_20261003'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists(): raise FileExistsError('never overwrite evidence')
    old=read(PRIOR); before=read(safe_file(LOGS+'/freeze_before.log'))
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):
        raise ValueError('missing initial integrity check')
    changed=[p for p,h in old['files_sha256'].items() if sha(safe_file(p))!=h]
    if set(changed)-ALLOWED: raise ValueError('old physics/evidence changed '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800: raise ValueError('HANDOFF over800')
    manifests=sorted((ROOT/OUT).glob('*/manifest.json'))
    if not manifests: raise ValueError('no real trial evidence')
    independent=read(ROOT/OUT/'independent_readback.json')
    if independent['reader_sha256']!=sha(ROOT/'scripts/assess_updated_stage.py') or independent['all_input_hashes_unchanged'] is not True:
        raise ValueError('independent reader/evidence changed')
    if independent['production_authorized'] is not False:
        raise ValueError('invalid independent scope')
    statuses=[]
    for path in manifests:
        report=read(path)
        if report['production_authorized'] is not False or report['status']=='RUNNING_UNCOMMITTED':
            raise ValueError('incomplete/unsafe trial scope')
        for rel,h in dict(report.get('input_sha256',{}),**report.get('source_sha256',{})).items():
            if sha(safe_file(rel))!=h: raise ValueError('run provenance changed '+rel)
        for name,h in report['artifact_sha256'].items():
            if Path(name).name!=name or sha(path.parent/name)!=h: raise ValueError('artifact drift '+name)
        statuses.append(dict(path=str(path.relative_to(ROOT)),status=report['status'],accepted_progress=report['accepted_progress']))
    files={safe_file(p) for p in old['files_sha256']}|{PRIOR}
    for directory,pattern in [('docs','UPDATED_TRIAL*'),('scripts','*updated*stage*.py'),
        ('scripts','freeze_updated_trial.py'),('src/afsd_backend/pysph','updated_trial.py'),
        ('src/afsd_backend/pysph','stage_iteration.py'),('tests/review','test_updated*.py'),
        ('tests/review','test_stage_iteration.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS,OUT): files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    logs=sorted((ROOT/LOGS).glob('*tests*.log'))+sorted((ROOT/LOGS).glob('trial_unit*.log'))
    result=dict(schema_version='1.0',date='2026-10-03',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],production_authorized=False,
        scope='Integrity only; updated global numerical stage is approximate, not sustained AFSD or paper validation',
        previous_freeze=str(PRIOR.relative_to(ROOT)),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        actual_test_logs={p.name:tests(str(p.relative_to(ROOT))) for p in logs},runs=statuses,
        independent_readback_status=independent['status'],
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),runs=statuses,preservation=result['preservation']),indent=2))


if __name__=='__main__':main()

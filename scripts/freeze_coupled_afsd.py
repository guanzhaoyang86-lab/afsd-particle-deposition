"""Append integrity evidence for integrated gates 1-4; not scientific acceptance."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha,read,safe_file,write_new
from freeze_process_feasible_predictor import tests

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/AXIAL_FEED_FREEZE.json'
TARGET=ROOT/'docs/COUPLED_AFSD_FREEZE.json'
LOGS='logs/coupled_afsd_20261002';OUTPUTS='outputs/coupled_afsd_20261002'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists():raise FileExistsError('never overwrite an evidence freeze')
    old=read(PRIOR);before=read(safe_file(LOGS+'/freeze_before.log'))
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):
        raise ValueError('previous integrity baseline missing')
    changed=[s for s,h in old['files_sha256'].items() if sha(safe_file(s))!=h]
    if set(changed)-ALLOWED:raise ValueError('unauthorized prior source/evidence change '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('HANDOFF over800 chars')
    cases={}
    for manifest in sorted((ROOT/OUTPUTS).glob('*/manifest.json')):
        m=read(manifest)
        for s,h in m['source_sha256'].items():
            if sha(safe_file(s))!=h:raise ValueError('runtime source drift '+s)
        for name,h in m.get('artifact_sha256',{}).items():
            if Path(name).name!=name or sha(manifest.parent/name)!=h:raise ValueError('artifact drift '+name)
        if m['production_authorized'] is not False:raise ValueError('production not permitted')
        cases[manifest.parent.name]={k:m.get(k) for k in ('status','accepted_steps','accepted_time_s','wall_seconds','error')}
    files={safe_file(s) for s in old['files_sha256']}|{PRIOR}
    for directory,pattern in [('docs','COUPLED_AFSD*'),('scripts','*coupled*.py'),
        ('tests/review','test_coupled*.py'),('src/afsd_backend/pysph','coupled_*.py'),
        ('configs/benchmarks','coupled_afsd.json')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS,OUTPUTS):files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    logs=sorted({p for pattern in ('*tests*.log','*unit*.log','old_*.log') for p in (ROOT/LOGS).glob(pattern)})
    inventory={p.name:tests(p.relative_to(ROOT).as_posix()) for p in logs}
    result=dict(schema_version='1.0',date='2026-10-02',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],production_authorized=False,
        scope='file integrity only; NOT scientific acceptance or sustained deposition',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        actual_test_logs=inventory,cases=cases,excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)})
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),preservation=result['preservation'],cases=cases),indent=2))


if __name__=='__main__':main()

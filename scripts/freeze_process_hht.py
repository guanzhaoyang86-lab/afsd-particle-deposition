"""Seal HHT local evidence without overwriting earlier failures or AFSD history."""
from pathlib import Path
import argparse
import json
from freeze_process_friction_q5 import sha,read,safe_file,write_new
from freeze_process_event import check_report

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_EVENT_FREEZE.json';TARGET=ROOT/'docs/PROCESS_HHT_FREEZE.json'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),'logs/process_hht_20261003/freeze_create.log',
          'logs/process_hht_20261003/freeze_final.log'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--preflight',type=Path,required=True)
    args=parser.parse_args();preflight=args.preflight.resolve()
    if TARGET.exists():raise FileExistsError('never overwrite frozen evidence')
    old=read(PRIOR);changed=[name for name,h in old['files_sha256'].items() if sha(safe_file(name))!=h]
    if set(changed)-ALLOWED:raise ValueError('old source/evidence changed '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('HANDOFF exceeds800')
    report=read(preflight);check_report(preflight.parent,report)
    if report['status']!='PASS_HHT_COMPONENT_PREFLIGHT_NOT_AFSD':raise ValueError('final component receipt required')
    runs=[]
    # Historical source receipts may be stale following the explicitly saved
    # repair. All corresponding artifacts and failures must still be intact.
    for base in ('logs','outputs'):
        for folder in (ROOT/base).glob('process_hht*'):
            for path in folder.rglob('report.json'):
                data=read(path);check_report(path.parent,data,False)
                runs.append(dict(path=path.relative_to(ROOT).as_posix(),status=data['status']))
    files={safe_file(name) for name in old['files_sha256']}|{PRIOR}
    files.update(safe_file(name) for name in report['source_sha256'])
    for folder,pattern in (('docs','PROCESS_HHT*'),('scripts','*hht*.py'),
            ('tests/review','test_hht*.py'),('src/afsd_backend/pysph','hht*.py')):
        files.update(path for path in (ROOT/folder).glob(pattern) if path.is_file())
    for base in ('logs','outputs'):
        for directory in (ROOT/base).glob('process_hht*'):
            files.update(path for path in directory.rglob('*') if path.is_file())
    files={path for path in files if path.relative_to(ROOT).as_posix() not in EXCLUDED}
    value=dict(schema_version='1.0',date='2026-10-03',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],project_complete=False,production_authorized=False,
        scope='FILE_INTEGRITY_ONLY_LOCAL_HHT_BENCHMARK_NOT_AFSD_TRANSPORT',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        final_preflight=preflight.relative_to(ROOT).as_posix(),runs=runs,new_afsd_steps=0,
        afsd_accepted_time_s=850.393632910904e-9,
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,value);print(json.dumps(dict(status='CREATED',files=len(files),new_afsd_steps=0,
        project_complete=False,preservation=value['preservation']),indent=2))


if __name__=='__main__':main()

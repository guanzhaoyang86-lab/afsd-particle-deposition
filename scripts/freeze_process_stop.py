"""Seal new stop evidence; integrity is NOT physical completion or authorization."""
from pathlib import Path
import argparse
import json
from check_process_stop import verify_previous,inventory
from probe_hht_curved import ROOT,sha,read,require
from freeze_process_friction_q5 import write_new

TARGET=ROOT/'docs/PROCESS_STOP_FREEZE.json'
PRIOR=ROOT/'docs/PROCESS_CURVED_FREEZE.json'
EXCLUDED={'docs/PROCESS_STOP_FREEZE.json','logs/process_stop_20261003/freeze_create.log',
    'logs/process_stop_20261003/freeze_final.log'}


def emit(value,name):
    output=json.dumps(value)+'\n';path=ROOT/'logs/process_stop_20261003'/name
    if not path.exists():path.write_text(output,encoding='utf-8')
    print(output,end='')


def verify_receipt(path):
    data=read(path)
    for name,digest in data.get('source_sha256',{}).items():require(sha(ROOT/name)==digest,'source mismatch '+name)
    for name,digest in data.get('artifact_sha256',{}).items():
        target=(path.parent/name).resolve();require(target.is_relative_to(path.parent.resolve()),'artifact escapes receipt')
        require(sha(target)==digest,'artifact mismatch '+name)
    return data


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    if args.verify:
        manifest=read(TARGET)
        for name,digest in manifest['files_sha256'].items():require(sha(ROOT/name)==digest,'sealed file changed '+name)
        result=dict(status='PASS_FILE_INTEGRITY_ONLY_NOT_PROCESS',files_checked=len(manifest['files_sha256']),
            project_complete=False,new_afsd_steps=0)
        emit(result,'freeze_final.log');return
    require(not TARGET.exists(),'never overwrite sealed evidence')
    previous=verify_previous();require(len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))<=800,'HANDOFF over800')
    checks=verify_receipt(ROOT/'logs/process_stop_20261003/preflight_initial/report.json')
    require(checks['status']=='PASS_STOP_COMPONENT_PREFLIGHT_NOT_PROCESS' and checks['tests_run']==233,'wrong component evidence')
    require(checks['source_sha256']==inventory(),'preflight stale')
    ref=verify_receipt(ROOT/'logs/process_stop_20261003/reference_endpoint2/report.json')
    require(ref['status']=='PASS_EVENT_REFERENCE_PREFIX_NEXT_BRANCH_BLOCKED' and ref['certified_common_end_time_s']==6.5e-6,'wrong reference prefix')
    run=verify_receipt(ROOT/'outputs/process_stop_20261003/scout_initial/report.json')
    require(run['status']=='STOPPED_LOCAL_EVENT_SCOUT_NOT_PROCESS' and run['time_s']==6e-6,'unexpected actual scout evidence')
    history=read(ROOT/'outputs/process_stop_20261003/scout_initial/trajectory/accepted_history.json')
    require(len(history)==2 and history[-1]['time_s']==6e-6,'wrong accepted local prefix')
    old=read(PRIOR);files={ROOT/name for name in old['files_sha256']}|{PRIOR,Path(__file__).resolve()}
    files.update(ROOT/name for name in checks['source_sha256'])
    files.update(path for path in (ROOT/'docs').glob('PROCESS_STOP*') if path.is_file())
    for base in ('logs','outputs'):
        folder=ROOT/base/'process_stop_20261003'
        files.update(path for path in folder.rglob('*') if path.is_file())
    files={path for path in files if path.relative_to(ROOT).as_posix() not in EXCLUDED}
    value=dict(schema_version='1.0',date='2026-10-03',scope='FILE_INTEGRITY_ONLY_LOCATED_STOP_DEVELOPMENT_NOT_AFSD',
        backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        project_complete=False,production_authorized=False,new_afsd_steps=0,
        afsd_accepted_time_s=850.393632910904e-9,afsd_transport_crossings=0,
        local_coupled_accepted_time_s=6e-6,independent_common_end_time_s=6.5e-6,
        coupled_stop_resolved=False,general3D_sticking='NOT_IMPLEMENTED',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=previous,old_failures_preserved=True,
        final_preflight='logs/process_stop_20261003/preflight_initial/report.json',
        files_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,value)
    emit(dict(status='CREATED_FILE_INTEGRITY_ONLY',files=len(files),project_complete=False,preservation=previous),'freeze_create.log')


if __name__=='__main__':main()

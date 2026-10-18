"""Seal response experiment evidence, never certify AFSD completion."""
from pathlib import Path
import argparse
import json
from check_process_response import verify_previous,inventory
from probe_hht_curved import ROOT,sha,read,require
from freeze_process_stop import verify_receipt
from freeze_process_friction_q5 import write_new

TARGET=ROOT/'docs/PROCESS_RESPONSE_FREEZE.json'
PRIOR=ROOT/'docs/PROCESS_STOP_FREEZE.json'
EXCLUDED={'docs/PROCESS_RESPONSE_FREEZE.json','logs/process_response_20261003/freeze_create.log',
    'logs/process_response_20261003/freeze_final.log'}


def emit(value,name):
    output=json.dumps(value)+'\n';path=ROOT/'logs/process_response_20261003'/name
    if not path.exists():path.write_text(output,encoding='utf-8')
    print(output,end='')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    if args.verify:
        manifest=read(TARGET)
        for name,digest in manifest['files_sha256'].items():
            path=(ROOT/name).resolve();require(path.is_relative_to(ROOT),'unsafe sealed path')
            require(sha(path)==digest,'sealed file changed '+name)
        emit(dict(status='PASS_FILE_INTEGRITY_ONLY_NOT_PROCESS',files_checked=len(manifest['files_sha256']),
            project_complete=False,new_afsd_steps=0),'freeze_final.log');return
    require(not TARGET.exists(),'never overwrite sealed evidence')
    previous=verify_previous();require(len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))<=800,'HANDOFF over800')
    checks=verify_receipt(ROOT/'logs/process_response_20261003/preflight_initial/report.json')
    require(checks['status']=='PASS_RESPONSE_COMPONENT_PREFLIGHT_NOT_PROCESS' and checks['tests_run']==245,'wrong component receipt')
    require(checks['source_sha256']==inventory(),'stale component sources')
    limits=verify_receipt(ROOT/'logs/process_response_20261003/limits_initial/report.json')
    require(limits['status']=='PASS_MATHEMATICAL_NONDIFFERENTIABILITY_EVIDENCE_NOT_SOLVER_VALIDATION'
        and len(limits['checks'])==130 and all(q['passed'] for q in limits['checks']),'wrong limits evidence')
    run=verify_receipt(ROOT/'outputs/process_response_20261003/saved_trial005_initial/report.json')
    require(run['status']=='PASS_PRIVATE_NUMERICAL_ROOT_NOT_EVENT_OR_PROCESS' and run['evaluations']==4,'wrong saved experiment')
    require(not run['accepted'] and not run['restart_allowed'] and run['actual_private_coupled_calls']==1,'incorrect private permissions')
    review=verify_receipt(ROOT/'logs/process_response_20261003/output_readback_initial/report.json')
    require(review['status']=='PASS_STORED_PRIVATE_READBACK_NOT_EVENT_OR_ACCURACY' and
        review['checks_run']==1218 and all(q['passed'] for q in review['checks']),'wrong independent readback')
    for name in ('REPORT','LIMITS','OUTPUT_REVIEW','END_TO_END_REVIEW','NEXT_WORK_PACKAGE'):
        require((ROOT/('docs/PROCESS_RESPONSE_'+name+'.md')).is_file(),'missing final report '+name)
    old=read(PRIOR);files={ROOT/name for name in old['files_sha256']}|{PRIOR,Path(__file__).resolve()}
    files.update(ROOT/name for name in checks['source_sha256'])
    for receipt in (limits,review):files.update(ROOT/name for name in receipt['source_sha256'])
    files.update(path for path in (ROOT/'docs').glob('PROCESS_RESPONSE*') if path.is_file())
    for name in ('logs','outputs'):
        files.update(path for path in (ROOT/name/'process_response_20261003').rglob('*') if path.is_file())
    files={path for path in files if path.relative_to(ROOT).as_posix() not in EXCLUDED}
    value=dict(schema_version='1.0',date='2026-10-03',scope='FILE_INTEGRITY_ONLY_PRIVATE_NUMERICAL_RESPONSE_NOT_AFSD',
        backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        project_complete=False,production_authorized=False,new_afsd_steps=0,new_local_published_steps=0,
        afsd_accepted_time_s=850.393632910904e-9,afsd_transport_crossings=0,
        local_coupled_accepted_time_s=6e-6,independent_common_end_time_s=6.5e-6,
        saved_private_root='PASS_NOT_EVENT_OR_PROCESS',coupled_stop_resolved=False,
        continuous_cross_zero_jacobian='NOT_DIFFERENTIABLE',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=previous,old_failures_preserved=True,
        final_preflight='logs/process_response_20261003/preflight_initial/report.json',
        files_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,value)
    emit(dict(status='CREATED_FILE_INTEGRITY_ONLY',files=len(files),project_complete=False,preservation=previous),'freeze_create.log')


if __name__=='__main__':main()

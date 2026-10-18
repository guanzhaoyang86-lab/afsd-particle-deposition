"""Seal saved real342 response evidence; never certify AFSD completion."""
from pathlib import Path
import argparse,json
from check_process_bridge import verify_previous,inventory
from probe_hht_curved import ROOT,sha,read,require
from freeze_process_stop import verify_receipt
from freeze_process_friction_q5 import write_new

TARGET=ROOT/'docs/PROCESS_BRIDGE_FREEZE.json'
PRIOR=ROOT/'docs/PROCESS_RESPONSE_FREEZE.json'
LOGS=ROOT/'logs/process_bridge_20261003'
EXCLUDED={'docs/PROCESS_BRIDGE_FREEZE.json','logs/process_bridge_20261003/freeze_create.log',
          'logs/process_bridge_20261003/freeze_final.log'}


def emit(value,name):
    path=LOGS/name
    if not path.exists():write_new(path,value)
    print(json.dumps(value),flush=True)


def supervision(path,log,source_base=None):
    r=read(path)
    require(r['return_code']==0 and r['wall_seconds']<=r['external_limit_seconds']<=600,'invalid bounded process')
    require(r['log_sha256']==sha(log),'supervised log changed')
    for name,digest in r.get('source_sha256',{}).items():
        child=(source_base/name).resolve()
        require(child.is_relative_to(source_base.resolve()) and sha(child)==digest,'supervised source changed')
    return r


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    if args.verify:
        frozen=read(TARGET)
        for name,digest in frozen['files_sha256'].items():
            path=(ROOT/name).resolve()
            require(path.is_relative_to(ROOT) and sha(path)==digest,'sealed file changed '+name)
        emit(dict(status='PASS_FILE_INTEGRITY_ONLY_NOT_PROCESS',files_checked=len(frozen['files_sha256']),
                  project_complete=False,new_afsd_steps=0),'freeze_final.log')
        return
    require(not TARGET.exists(),'never overwrite evidence')
    previous=verify_previous()
    require(len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))<=800,'HANDOFF exceeds800')
    tests=verify_receipt(LOGS/'preflight_initial/report.json')
    require(tests['status']=='PASS_BRIDGE_PREFLIGHT_NOT_PROCESS' and tests['tests_run']==273,'wrong regression receipt')
    require(tests['source_sha256']==inventory(),'test sources changed')
    supervision(LOGS/'preflight_initial_supervision.json',LOGS/'preflight_initial.log')
    runpath=ROOT/'outputs/process_bridge_20261003/real342_initial/report.json'
    run=verify_receipt(runpath)
    for name,digest in run['input_sha256'].items():
        path=(ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and sha(path)==digest,'bound input changed '+name)
    require(run['status']=='NONSMOOTH_NOT_CERTIFIED' and run['actual_evaluations']==run['actual_evaluations_attempted']==4,
            'wrong actual response experiment')
    require(run['base_reentry_exact'] and run['parent_unchanged'] and all(x['audit_passed'] for x in run['measurements']),
            'mapped accounts/history failed')
    require(all(run[k] is False for k in ('accepted','restart_allowed','can_publish','production_authorized')),
            'diagnostic gained physical authority')
    require(run['new_afsd_steps']==run['new_transport_events']==0,'unexpected physical progress')
    supervision(LOGS/'real342_initial_supervision.json',LOGS/'real342_initial.log')
    review=read(LOGS/'output_readback_repair1/report.json')
    require(review['status']=='PASS_INDEPENDENT_READBACK_NOT_PROCESS_ACCEPTANCE' and review['failure_count']==0
            and review['check_count']==471275,'independent readback incomplete')
    require(all(x['failures']==0 for x in review['checks'].values()),'readback category failed')
    require(review['checks']['EMPTY_NOT_APPLICABLE_shape_finite']['count']==3936,'empty checks misclassified')
    summary=review['summary']
    require(summary['no_new_solver_evaluations'] and summary['base_reentry_exact'] and
            summary['original_report_sha256']==sha(runpath),'readback does not bind actual outputs')
    require(summary['readback_source_sha256']==sha(LOGS/'output_readback.py') and
            summary['contract_sha256']==sha(ROOT/'docs/PROCESS_BRIDGE_OUTPUT_REVIEW_CONTRACT.md'),'readback source drift')
    supervision(LOGS/'output_readback_repair1_supervision.json',LOGS/'output_readback_repair1.log',LOGS)
    response=verify_receipt(ROOT/'logs/process_bridge_response_review_20261003/readback.json')
    require(response['status']=='COMPLETE_READBACK_NO_NEW_EVALUATIONS' and response['new_physical_evaluations']==0,
            'response readback invalid')
    require(sha(ROOT/response['full_arrays']['path'])==response['full_arrays']['sha256'],'response arrays changed')
    actiondir=ROOT/'logs/process_bridge_action_repair1_20261004'
    action=verify_receipt(actiondir/'report.json')
    require(action['status']=='COMPLETE_ALGEBRA_CHECK_NO_PHYSICAL_ACCEPTANCE' and action['new_evaluations']==0
            and action['parent_unchanged'] and action['source_unchanged'],'action incomplete')
    require([x['name'] for x in action['models']]==['structural','continuity'],'action inventory changed')
    supervision(actiondir/'supervision.json',actiondir/'run.log')
    causaldir=ROOT/'logs/process_bridge_action_review_20261004'
    causal=read(causaldir/'report.json');timing=read(causaldir/'supervision.json')
    require(causal['status']=='READ_ONLY_SIGNED_ARRAY_DECOMPOSITION_NOT_CAUSAL_PROOF' and
            causal['new_evaluations']==causal['preconditioner_matrix_builds']==causal['preconditioner_applications']==0,
            'causal scope changed')
    require(causal['input_files_checked']==172 and causal['all_watched_files_unchanged_after'],'causal incomplete')
    for name,digest in causal['watched_sha256'].items():
        path=(ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and sha(path)==digest,'causal input changed '+name)
    require(causal['reader_sha256']==sha(causaldir/'review_saved_action.py') and
            timing['report_sha256']==sha(causaldir/'report.json'),'causal output/source changed')
    require(timing['completed'] and timing['exit_code']==0 and timing['invocations']==1 and
            timing['external_seconds']<=timing['external_limit_seconds']<=600,'causal bounded execution failed')
    require(read(ROOT/'logs/process_bridge_action_20261003/report.json')['status']=='FAIL_ALGEBRA_CHECK',
            'original serialization failure lost')
    require(read(LOGS/'output_readback_initial/report.json')['status']=='FAIL_INDEPENDENT_READBACK',
            'original empty-array failure lost')
    for suffix in ('REPORT','OUTPUT_REVIEW','RESPONSE_REVIEW','ACTION_REVIEW','ACTION_CAUSAL_REVIEW'):
        require((ROOT/('docs/PROCESS_BRIDGE_'+suffix+'.md')).is_file(),'missing final review '+suffix)
    old=read(PRIOR)
    files={ROOT/name for name in old['files_sha256']}|{PRIOR,Path(__file__).resolve()}
    for receipt in (tests,run,response,action):files.update(ROOT/name for name in receipt.get('source_sha256',{}))
    files.update(p for p in (ROOT/'docs').glob('PROCESS_BRIDGE*') if p.is_file())
    for base in ('logs','outputs'):
        for folder in (ROOT/base).glob('process_bridge*'):
            if folder.is_dir():files.update(p for p in folder.rglob('*') if p.is_file())
    files={p.resolve() for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    require(all(p.is_relative_to(ROOT) and p.is_file() for p in files),'unsafe seal')
    frozen=dict(schema_version='1.0',date='2026-10-04',
        scope='FILE_INTEGRITY_ONLY_REAL342_RESPONSE_AND_OFFLINE_ACTION_NOT_AFSD',
        backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        project_complete=False,production_authorized=False,new_afsd_steps=0,new_transport_events=0,
        afsd_accepted_time_s=850.393632910904e-9,afsd_transport_crossings=0,
        diagnostic_evaluations=4,diagnostic_status=run['status'],new_accepted_state=False,
        independent_readback='PASS_ACCOUNTING_NOT_ROOT',algebra_action='COMPLETE_NOT_PHYSICAL_ACCEPTANCE',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=previous,old_failures_preserved=True,
        final_preflight='logs/process_bridge_20261003/preflight_initial/report.json',
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,frozen)
    emit(dict(status='CREATED_FILE_INTEGRITY_ONLY',files=len(files),project_complete=False,preservation=previous),
         'freeze_create.log')


if __name__=='__main__':main()

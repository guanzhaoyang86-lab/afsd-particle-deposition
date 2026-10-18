"""Seal bounded material-response evidence; NOT AFSD completion."""
from pathlib import Path
import argparse,json,re
from probe_hht_curved import ROOT,sha,read,require
from freeze_process_friction_q5 import write_new
TARGET=ROOT/'docs/MATERIAL_RESPONSE_FREEZE.json'
PRIOR=ROOT/'docs/PROCESS_BRIDGE_FREEZE.json'
LOG=ROOT/'logs/process_material_response_20261004'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={'docs/MATERIAL_RESPONSE_FREEZE.json',
          'logs/process_material_response_20261004/freeze_create.log',
          'logs/process_material_response_20261004/freeze_final.log'}

def bind(values,base):
    for name,digest in values.items():
        path=(base/name).resolve()
        require(path.is_relative_to(base.resolve()) and sha(path)==digest,'bound file changed '+name)

def tests(path,count):
    r=read(path)
    require(r['return_code']==0 and r['sources_unchanged'] and
            r['source_hashes_before']==r['source_hashes_after'],'invalid test receipt')
    require(r['wall_seconds']<=r['external_limit_seconds']<=600,'test budget')
    bind(r['source_hashes_before'],ROOT)
    logs=[p for p in path.parent.glob('*.log') if sha(p)==r['log_sha256']]
    require(len(logs)==1,'test log binding')
    text=logs[0].read_text(encoding='utf-8')
    require(re.search(r'^Ran '+str(count)+r' tests in ',text,re.M) and re.search(r'^OK\s*$',text,re.M),'test count/outcome')
    return r

def emit(value,name):
    if not (LOG/name).exists():write_new(LOG/name,value)
    print(json.dumps(value),flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    if args.verify:
        r=read(TARGET);bind(r['files_sha256'],ROOT)
        emit(dict(status='PASS_FILE_INTEGRITY_ONLY_NOT_AFSD',files_checked=len(r['files_sha256']),
                  project_complete=False,new_afsd_steps=0),'freeze_final.log');return
    require(not TARGET.exists(),'never overwrite evidence')
    old=read(PRIOR);require(len(old['files_sha256'])==4745,'wrong prior evidence inventory')
    changes=[]
    for name,digest in old['files_sha256'].items():
        path=(ROOT/name).resolve();require(path.is_relative_to(ROOT) and path.is_file(),'old file missing')
        if sha(path)!=digest:
            require(name in ALLOWED,'old evidence changed '+name);changes.append(name)
    require(len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))<=800,'HANDOFF over800')
    one=tests(LOG/'initial_supervision.json',11)
    two=tests(ROOT/'logs/material_response_20261004/stage_tests_initial_supervision.json',12)
    actiondir=ROOT/'outputs/material_response_20261004/saved_action_initial'
    action=read(actiondir/'report.json')
    require(action['status']=='COMPLETE_SAVED_ACTION_NOT_ROOT_OR_PROCESS' and action['parent_unchanged'],'action failed')
    require(action['new_evaluations']==action['new_afsd_steps']==action['new_transport_events']==0,'unexpected physical evaluation')
    require(all(action[k] is False for k in ('accepted','restart_allowed','can_publish','production_authorized')),'authority changed')
    bind(action['source_sha256'],ROOT);bind(action['artifact_sha256'],actiondir)
    s=read(LOG/'saved_action_initial_supervision.json')
    require(s['return_code']==0 and s['wall_seconds']<=s['external_limit_seconds']<=600,'action budget')
    require(s['log_sha256']==sha(LOG/'saved_action_initial.log'),'action log changed')
    require(s['supervisor_sha256']==sha(ROOT/'scripts/run_material_response_action.py'),'action supervisor changed')
    # Full review journal is required, not inferred from an existing filename.
    readback=read(LOG/'output_readback.json')
    require(readback['status']=='PASS_READBACK_ONLY' and readback['sources_unchanged'] and
            readback['source_hashes_before']==readback['source_hashes_after'],'independent action readback failed')
    bind(readback['source_hashes_before'],ROOT)
    require(readback['new_evaluations']==readback['matrix_builds']==readback['preconditioner_applications']==0,'readback ran physics/action')
    rb=read(LOG/'output_readback_supervision.json')
    require(rb['return_code']==0 and rb['wall_seconds']<=rb['external_limit_seconds']<=600,'readback budget')
    require(rb['log_sha256']==sha(LOG/'output_readback.log') and
            rb['reader_sha256']==sha(LOG/'output_readback.py') and
            rb['supervisor_sha256']==sha(LOG/'output_readback_run.py'),'readback sources/log changed')
    for name in ('REPORT','OUTPUT_REVIEW','INTEGRATION_DESIGN','CONTRACT','ACTION_CONTRACT','TEST_CONTRACT'):
        require((ROOT/('docs/MATERIAL_RESPONSE_'+name+'.md')).is_file(),'missing report '+name)
    files={ROOT/name for name in old['files_sha256']}|{PRIOR,Path(__file__).resolve()}
    for r in (one,two):files.update(ROOT/n for n in r['source_hashes_before'])
    files.update(ROOT/n for n in action['source_sha256'])
    files.update(ROOT/n for n in readback['source_hashes_before'])
    for pattern in ('MATERIAL_RESPONSE*','MATERIAL_STAGE_RESPONSE*'):
        files.update(p for p in (ROOT/'docs').glob(pattern) if p.is_file())
    for rel in ('logs/process_material_response_20261004','logs/material_response_20261004','outputs/material_response_20261004'):
        files.update(p for p in (ROOT/rel).rglob('*') if p.is_file())
    files={p.resolve() for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    value=dict(schema_version='1.0',date='2026-10-04',scope='NUMERICAL_MATERIAL_RESPONSE_FILE_INTEGRITY_NOT_AFSD',
        backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        project_complete=False,production_authorized=False,new_afsd_steps=0,new_physical_evaluations=0,
        afsd_accepted_time_s=850.393632910904e-9,afsd_transport_crossings=0,
        component_tests_passed=23,saved_action='COMPLETE_NOT_ROOT',temperature_action_improved=False,
        new_root_or_transport='NOT_RUN',old_failures_preserved=True,
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(files_checked=4745,authorized_document_changes=sorted(changes)),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)},
        excluded_self_or_future_paths=sorted(EXCLUDED))
    write_new(TARGET,value)
    emit(dict(status='CREATED_FILE_INTEGRITY_ONLY',files=len(files),project_complete=False,
              preservation=value['preservation']),'freeze_create.log')

if __name__=='__main__':main()

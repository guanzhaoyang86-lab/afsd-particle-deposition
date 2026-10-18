"""Freeze administrative recovery and eight-hour continuation; no GPU solves."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/gpu_acceleration_20261007'
OLD_SESSION=BASE/'event_resume_90m_20261010_v1'
PARENT=OLD_SESSION/'attempt_0013'
FAILED=OLD_SESSION/'attempt_0014'
PARENT_SHA256='0079328d9ea46263d62cf975299497ea6774136d22e5bdf0ffba5d6cdcb91c97'
ADMINISTRATIVE_FILES=(
    'scripts/run_gpu_event_resume_8h_20261011.py',
    'scripts/afsd_event_resume_8h_worker_20261011.py',
    'scripts/check_gpu_event_resume_8h_20261011.py',
    'tests/review/test_gpu_event_resume_8h_20261011.py',
    'docs/GPU_EVENT_RESUME_8H_20261011.md',
    'scripts/afsd_multi_event_refinement_20261011.py',
    'tests/review/test_afsd_multi_event_refinement_20261011.py')
TEST_FILES=(ADMINISTRATIVE_FILES[3],ADMINISTRATIVE_FILES[6])


def require(value,message):
    if not value: raise ValueError(message)


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def rel(path):
    return Path(path).relative_to(ROOT).as_posix()


def merge(*maps):
    result={}
    for mapping in maps:
        for key,value in mapping.items():
            require(key not in result or result[key]==value,'conflicting frozen evidence: '+key)
            result[key]=value
    return result


def verify(mapping):
    for name,digest in mapping.items():
        path=(ROOT/name).resolve()
        require(path.is_relative_to(ROOT),'evidence escaped project')
        require(sha(path)==digest,'frozen evidence changed: '+name)


def stage_evidence(folder, *, passed):
    meta=read(folder/'manifest.json'); supervision=read(folder/'supervision.json')
    require(supervision['exit_code']==(0 if passed else 2) and supervision['source_inputs_still_match'] is True
        and supervision['interrupted'] is False and supervision['manifest_sha256']==sha(folder/'manifest.json')
        and supervision['log_sha256']==sha(folder/'run.log'),'complete stage supervision/SHA required')
    require(meta['accepted'] is False and meta['restart_allowed'] is False and meta['diagnostic_only'] is True,
        'unmodified diagnostic status required')
    require(meta['status']==('PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION' if passed else 'FAIL_GPU_STAGE_NOT_ACCEPTED')
        and meta['screened_continuation_authorized'] is passed,'stage outcome must remain unchanged')
    artifacts=merge(meta['artifacts_sha256'],{rel(folder/n):sha(folder/n) for n in ('manifest.json','supervision.json','run.log')})
    if passed:
        require(sha(ROOT/meta['screened_continuation_checkpoint'])==meta['screened_continuation_sha256'],
            'complete checkpoint SHA mismatch')
        screen=read(folder/'pair_screen.json')
        require(screen==meta['pair_screen'] and screen['passed'] is True and screen['global_screen']['pass'] is True
            and screen['local_screen']['passed'] is True and screen['crossing_ids_match'] is True,
            'original complete pair screens must pass')
    return meta,artifacts


def parent_evidence():
    summary=read(OLD_SESSION/'final_summary.json'); audit=read(OLD_SESSION/'final_history_audit.json')
    require(summary['finished'] is True and summary['status']=='STOPPED_CHILD_FAILURE' and summary['exit_code']==2
        and summary['completed_stage_count']==13 and summary['attempt_count']==14
        and summary['latest_pass_stage']==rel(PARENT) and summary['checkpoint_sha256']==PARENT_SHA256
        and summary['physical_time_s']==.08134326171874973 and summary['screened_steps']==527,
        'exact old thirteen-PASS/fourteenth-FAIL session required')
    require(audit['failed'] is False and audit['final_audits']==1
        and audit['unique_files']==audit['expected_files']==audit['locked_files']>0,
        'old session final history audit must be complete')
    names=['protocol.json','final_summary.json','final_history_audit.json','deposition_progress.json']
    names += [f'attempt_{i:04d}_receipt.json' for i in range(1,15)]
    artifacts={rel(OLD_SESSION/name):sha(OLD_SESSION/name) for name in names}
    parent=None
    for i in range(1,14):
        folder=OLD_SESSION/f'attempt_{i:04d}'
        meta,stage_artifacts=stage_evidence(folder,passed=True)
        require(len(meta['source_sha256'])==636,'exact original physical source count required')
        receipt=read(OLD_SESSION/f'attempt_{i:04d}_receipt.json')
        require(receipt['verified'] is True and receipt['exit_code']==0 and receipt['stage']==rel(folder)
            and receipt['checkpoint_sha256']==meta['screened_continuation_sha256'],
            'old complete administrative receipt required')
        if parent is not None:
            require(meta['parent_sha256']==parent['screened_continuation_sha256']
                and meta['source_sha256']==parent['source_sha256'],'old selected PASS chain changed')
        artifacts=merge(artifacts,stage_artifacts); parent=meta
    require(parent['screened_continuation_sha256']==PARENT_SHA256
        and parent['screened_end_time_s']==.08134326171874973 and parent['screened_end_step']==527,
        'exact complete starting checkpoint required')
    failure,failed_artifacts=stage_evidence(FAILED,passed=False)
    require(failure['parent_sha256']==PARENT_SHA256 and failure['parent_gpu_directory']==rel(PARENT)
        and failure['source_sha256']==parent['source_sha256']
        and failure['error']=="GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')",
        'registered old failure must retain its exact full parent and original physics')
    artifacts=merge(artifacts,failed_artifacts)
    return parent,artifacts


def main():
    output=BASE/'event_resume_8h_components_20261011_v1.json'; log=output.with_suffix('.log')
    require(not output.exists() and not log.exists(),'exclusive new component receipt/log required')
    inherited=BASE/'event_resume_90m_components_20261010_v1.json'; prior=read(inherited)
    require(prior['status']=='PASS_EVENT_RESUME_90M_COMPONENTS' and prior['exit_code']==0,
        'executed inherited 90-minute component receipt required')
    inherited_sha=sha(inherited); parent,artifacts=parent_evidence()
    closure=merge(prior['source_sha256'],prior['artifacts_sha256'],{rel(inherited):inherited_sha},parent['source_sha256'])
    verify(closure); verify(artifacts)
    sources=merge(closure,{name:sha(ROOT/name) for name in ADMINISTRATIVE_FILES})
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH='')
    checks=[]; tests=0
    with log.open('x',encoding='utf-8') as stream:
        for script in TEST_FILES:
            command=[sys.executable,'-B',script]; began=time.perf_counter()
            result=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True,timeout=120.)
            stream.write(json.dumps(command)+'\n'+result.stdout+result.stderr+'\n'); stream.flush()
            count=re.search(r'Ran (\d+) tests?',result.stderr)
            tests += int(count.group(1)) if count else 0
            checks.append(dict(command=command,exit_code=result.returncode,wall_seconds=time.perf_counter()-began,
                stdout=result.stdout,stderr=result.stderr))
    verify(sources); verify(artifacts)
    success=len(checks)==2 and all(r['exit_code']==0 for r in checks) and tests>=20
    artifacts[rel(log)]=sha(log)
    receipt=dict(status='PASS_EVENT_RESUME_8H_COMPONENTS' if success else 'FAIL_EVENT_RESUME_8H_COMPONENTS',
        exit_code=0 if success else 2,source_sha256=sources,artifacts_sha256=artifacts,
        physical_source_sha256=dict(parent['source_sha256']),physical_source_changed=False,
        administrative_files=list(ADMINISTRATIVE_FILES),tests=tests,checks=checks,
        inherited_receipt=rel(inherited),inherited_receipt_sha256=inherited_sha,inherited_tests=prior['tests'],
        parent_directory=rel(PARENT),parent_checkpoint_sha256=PARENT_SHA256,failed_stage=rel(FAILED),
        old_session_directory=rel(OLD_SESSION),authorized_at='2026-10-11T15:41:47+08:00',
        absolute_stop='2026-10-11T23:41:47+08:00',registered_wall_limits_s=dict(base_physical_soft=1260.,
            maximum_physical_soft=3600.,worker_extra=30.,controller_extra=120.),
        scope='ADMINISTRATIVE_ONLY_EXACT_V5_COMPLETE_PARENT_MULTIEVENT_RETRY_8H',
        gpu_campaign_runs=0,full_physical_solver_runs=0,accepted=False,restart_allowed=False)
    with output.open('x',encoding='utf-8') as stream:
        json.dump(receipt,stream,indent=2,allow_nan=False); stream.write('\n')
    print(json.dumps(dict(status=receipt['status'],tests=tests,receipt=str(output))))
    return receipt['exit_code']


if __name__=='__main__': raise SystemExit(main())

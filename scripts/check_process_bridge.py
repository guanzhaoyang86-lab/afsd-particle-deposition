"""Old components plus real-process bridge interfaces, no AFSD candidate run."""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
from check_process_response import TESTS as OLD_TESTS, inventory as previous_inventory
from check_process_event import suite_result
from probe_hht_curved import ROOT,sha,read,write,require

TESTS=OLD_TESTS+['afsd_process_bridge','afsd_bridge_audit','afsd_bridge_driver']
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
PRIOR=ROOT/'docs/PROCESS_RESPONSE_FREEZE.json'


def verify_previous():
    old=read(PRIOR);require(len(old['files_sha256'])==4600,'wrong previous inventory');changed=[]
    for name,digest in old['files_sha256'].items():
        path=(ROOT/name).resolve();require(path.is_relative_to(ROOT) and path.is_file(),'unsafe/missing '+name)
        if sha(path)!=digest:
            require(name in ALLOWED,'old evidence changed '+name);changed.append(name)
    return dict(files_checked=4600,authorized_document_changes=sorted(changed),manifest_sha256=sha(PRIOR))


def inventory():
    names=set(previous_inventory())
    names.update(p.relative_to(ROOT).as_posix() for p in (ROOT/'src').rglob('*.py'))
    names.update(['scripts/afsd_process_bridge.py','scripts/afsd_bridge_audit.py',
        'scripts/probe_afsd_process_bridge.py','scripts/check_process_bridge.py','scripts/run_process_bridge_command.py',
        'scripts/probe_updated_stage.py','scripts/assess_process_route.py','scripts/run_coupled_afsd.py',
        'scripts/check_process_response.py','tests/review/test_afsd_process_bridge.py','tests/review/test_afsd_bridge_audit.py',
        'tests/review/test_afsd_bridge_driver.py',
        'docs/PROCESS_RESPONSE_FREEZE.json','docs/PROCESS_RESPONSE_NEXT_WORK_PACKAGE.md',
        'docs/PROCESS_BRIDGE_CONTRACT.md','docs/PROCESS_BRIDGE_AUDIT_CONTRACT.md','docs/PROCESS_BRIDGE_TEST_CONTRACT.md'])
    return {name:sha(ROOT/name) for name in sorted(names)}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    out=ap.parse_args().output.resolve();require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output')
    out.mkdir(parents=True);started=time.perf_counter();sources=None
    report=dict(status='RUNNING',tests=[],new_afsd_steps=0,actual_bridge_evaluations=0,production_authorized=False)
    code=2
    try:
        report['prior_before']=verify_previous();sources=inventory();report['source_sha256']=sources
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        for name in TESTS:
            remain=55-(time.perf_counter()-started);require(remain>0,'55s preflight budget')
            cmd=[sys.executable,'-B','tests/review/test_'+name+'.py'];tick=time.perf_counter()
            try:
                result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                    text=True,encoding='utf-8',errors='replace',timeout=remain)
                output,rc=result.stdout,result.returncode
            except subprocess.TimeoutExpired as error:
                output=error.stdout or b''
                if isinstance(output,bytes):output=output.decode('utf-8',errors='replace')
                output+='\nTIMEOUT_SHARED55S\n';rc=124
            (out/(name+'.log')).write_text(output,encoding='utf-8')
            row=dict(name=name,**suite_result(name,output,rc),command=cmd,return_code=rc,
                wall_seconds=time.perf_counter()-tick);report['tests'].append(row)
            write(out/'report.json',report);print(json.dumps(row),flush=True)
            require(row['all_tests_pass'],'component failed '+name)
        require(sum(x['tests_run'] for x in report['tests'] if x['name'] in OLD_TESTS)==245,'old245 inventory differs')
        for name,count in {'afsd_process_bridge':14,'afsd_bridge_audit':11,'afsd_bridge_driver':3}.items():
            require(next(x['tests_run'] for x in report['tests'] if x['name']==name)==count,'new suite inventory differs '+name)
        require(inventory()==sources,'source drift');report['status']='PASS_BRIDGE_PREFLIGHT_NOT_PROCESS';code=0
    except Exception as error:report.update(status='FAIL_BRIDGE_PREFLIGHT',error=repr(error))
    try:
        report['prior_after']=verify_previous()
        if sources is not None:require(inventory()==sources,'source drift after execution')
    except Exception as error:report.update(status='FAIL_BRIDGE_PREFLIGHT_INTEGRITY',integrity_error=repr(error));code=2
    report.update(wall_seconds=time.perf_counter()-started,tests_run=sum(x['tests_run'] for x in report['tests']),
        not_run_suites=[name for name in TESTS if name not in [x['name'] for x in report['tests']]],
        artifact_sha256={p.name:sha(p) for p in out.glob('*.log')})
    write(out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('tests','source_sha256','artifact_sha256')}),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

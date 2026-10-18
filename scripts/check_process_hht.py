"""Bounded new HHT component preflight, not AFSD acceptance or old FAIL waiver."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time

from probe_hht_coupled import ROOT,inventory,sha,read,write,require
from check_process_event import OLD_TESTS,suite_result

NEW_TESTS=['hht_linear','hht_path','hht_coupled','hht_reference','hht_runner','hht_publication','hht_restart']
TESTS=OLD_TESTS+NEW_TESTS
PRIOR=ROOT/'docs/PROCESS_EVENT_FREEZE.json'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}


def freeze_check():
    frozen=read(PRIOR)['files_sha256'];require(len(frozen)==3656,'unexpected frozen inventory')
    changed=[]
    for name,expected in frozen.items():
        path=(ROOT/name).resolve();require(path.is_relative_to(ROOT) and path.is_file(),'missing old evidence '+name)
        if sha(path)!=expected:
            require(name in ALLOWED,'old scientific evidence changed '+name);changed.append(name)
    return dict(files_checked=len(frozen),allowed_document_changes=changed,manifest_sha256=sha(PRIOR))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output.resolve();require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);started=time.perf_counter();code=2
    report=dict(status='RUNNING',classification='HHT_LOCAL_COMPONENT_PREFLIGHT_NOT_AFSD',tests=[],
        planned_suites=TESTS,afsd_steps=0,production_authorized=False,
        old_known_failures='PRESERVED_BY_HASH_NOT_RECLASSIFIED_NOT_RERUN_IN_THIS_PREFLIGHT',
        old_event_suite_counts=dict(total=147,passing=145,failing=2),suite_timeout_s=60)
    try:
        report['previous_freeze_before']=freeze_check();sources=inventory()
        extra=['scripts/check_process_hht.py','scripts/check_process_event.py',
            'docs/PROCESS_HHT_LINEAR_CONTRACT.md']+['tests/review/test_'+name+'.py' for name in TESTS]
        sources.update({name:sha(ROOT/name) for name in extra});report['source_sha256']=sources
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        for name in TESTS:
            command=[sys.executable,'-B','tests/review/test_'+name+'.py'];tick=time.perf_counter()
            try:
                result=subprocess.run(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                    text=True,encoding='utf-8',errors='replace',timeout=60)
                output,rc=result.stdout,result.returncode
            except subprocess.TimeoutExpired as error:
                output=error.stdout or b''
                if isinstance(output,bytes):output=output.decode('utf-8',errors='replace')
                output+='\nTIMEOUT_60S\n';rc=124
            (out/(name+'.log')).write_text(output,encoding='utf-8')
            row=dict(name=name,**suite_result(name,output,rc),command=command,return_code=rc,
                wall_seconds=time.perf_counter()-tick,log=name+'.log')
            report['tests'].append(row);write(out/'report.json',report);print(json.dumps(row),flush=True)
            require(row['all_tests_pass'],'component suite failed '+name)
        report['tests_run']=sum(row['tests_run'] for row in report['tests'])
        report['old_tests_run']=sum(row['tests_run'] for row in report['tests'] if row['name'] in OLD_TESTS)
        require(report['old_tests_run']==103,'old regression count changed')
        for name,expected in sources.items():require(sha(ROOT/name)==expected,'source drift '+name)
        report['previous_freeze_after']=freeze_check()
        report['status']='PASS_HHT_COMPONENT_PREFLIGHT_NOT_AFSD';code=0
    except Exception as exc:report.update(status='FAIL_HHT_PREFLIGHT',error=repr(exc))
    report.update(wall_seconds=time.perf_counter()-started,
        artifact_sha256={path.name:sha(path) for path in out.glob('*.log')})
    write(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','artifact_sha256','tests')}),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

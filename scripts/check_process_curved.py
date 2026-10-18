"""No-trajectory preflight for the registered slow-load diagnostic scope."""
from pathlib import Path
import argparse
import os
import subprocess
import sys
import time
import json
from check_process_event import suite_result
from check_process_hht import TESTS as HHT_TESTS
from probe_hht_coupled import inventory as original_inventory
from probe_hht_curved import ROOT,inventory,sha,read,write,require

TESTS=HHT_TESTS+['hht_slow_inputs','hht_slow_reference','hht_slow_runner',
    'hht_curved_identity','hht_curved_contact','hht_curved_coupled']
PRIOR=ROOT/'docs/PROCESS_SLOW_FREEZE.json'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}


def verify_previous():
    old=read(PRIOR);require(len(old['files_sha256'])==4221,'previous inventory differs')
    changed=[]
    for name,h in old['files_sha256'].items():
        path=(ROOT/name).resolve();require(path.is_relative_to(ROOT) and path.is_file(),'missing previous file '+name)
        if sha(path)!=h:
            require(name in ALLOWED,'previous science/code changed '+name);changed.append(name)
    return dict(files_checked=4221,authorized_document_changes=changed,manifest_sha256=sha(PRIOR))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output.resolve();require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);start=time.perf_counter();report=dict(status='RUNNING',tests=[],planned_suites=TESTS,
        classification='CURVED_STARTUP_COMPONENT_PREFLIGHT_NOT_PROCESS',production_authorized=False,
        old_accuracy_failures='RETAINED_NOT_RECLASSIFIED',new_afsd_steps=0)
    code=2
    try:
        report['prior_before']=verify_previous();sources=dict(original_inventory(),**inventory())
        names=['scripts/check_process_curved.py','scripts/check_process_hht.py','scripts/check_process_event.py',
            'scripts/probe_hht_curved_regression.py']
        names+=['tests/review/test_'+name+'.py' for name in TESTS]
        sources.update({name:sha(ROOT/name) for name in names});report['source_sha256']=sources
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
                wall_seconds=time.perf_counter()-tick,log=name+'.log');report['tests'].append(row)
            write(out/'report.json',report);print(json.dumps(row),flush=True)
            require(row['all_tests_pass'],'component test failed '+name)
        require(sum(row['tests_run'] for row in report['tests'] if row['name'] in HHT_TESTS)==145,'old HHT regression inventory differs')
        for name,h in sources.items():require(sha(ROOT/name)==h,'source drift during preflight '+name)
        report.update(status='PASS_CURVED_COMPONENT_PREFLIGHT_NOT_PROCESS',prior_after=verify_previous())
        code=0
    except Exception as error:report.update(status='FAIL_CURVED_PREFLIGHT',error=repr(error))
    report.update(tests_run=sum(row['tests_run'] for row in report['tests']),
        not_run_suites=[name for name in TESTS if name not in [row['name'] for row in report['tests']]],
        wall_seconds=time.perf_counter()-start,artifact_sha256={path.name:sha(path) for path in out.glob('*.log')})
    write(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','artifact_sha256','tests')}),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

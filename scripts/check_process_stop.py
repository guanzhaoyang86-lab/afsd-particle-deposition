"""Bounded old/new component gate. No AFSD or local trajectory is advanced."""
from pathlib import Path
import argparse
import os
import subprocess
import sys
import time
import json
from check_process_event import suite_result
from check_process_curved import TESTS as OLD_TESTS
from probe_hht_curved import ROOT,sha,read,write,require

TESTS=OLD_TESTS+['hht_zero_event','hht_event_reference','hht_stop_contact','hht_stop_stage','hht_stop_transactions']
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}


def verify_previous():
    path=ROOT/'docs/PROCESS_CURVED_FREEZE.json';old=read(path);changed=[]
    require(len(old['files_sha256'])==4365,'old inventory differs')
    for name,digest in old['files_sha256'].items():
        target=(ROOT/name).resolve();require(target.is_relative_to(ROOT) and target.is_file(),'missing old file '+name)
        if sha(target)!=digest:
            require(name in ALLOWED,'old source/evidence changed '+name);changed.append(name)
    return dict(files_checked=4365,authorized_document_changes=changed,manifest_sha256=sha(path))


def inventory():
    names={q.relative_to(ROOT).as_posix() for q in (ROOT/'src/afsd_backend/pysph').glob('*.py')}
    names.update('tests/review/test_'+test+'.py' for test in TESTS)
    names.update(['scripts/check_process_stop.py','scripts/probe_hht_stop.py','scripts/hht_stop_boundary.py',
        'scripts/run_process_stop_command.py',
        'scripts/probe_hht_curved.py','scripts/probe_hht_slow.py','scripts/probe_hht_coupled.py',
        'scripts/hht_slow_inputs.py','scripts/reference_hht_event.py','scripts/reference_hht_slow.py',
        'scripts/reference_hht_block.py','docs/PROCESS_STOP_CONTRACT.md',
        'docs/PROCESS_STOP_BENCHMARK_CONTRACT.md'])
    return {name:sha(ROOT/name) for name in sorted(names)}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    out=ap.parse_args().output.resolve();require(out.is_relative_to(ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);start=time.perf_counter();report=dict(status='RUNNING',tests=[],
        new_afsd_steps=0,production_authorized=False,planned_suites=TESTS);code=2
    try:
        report['prior_before']=verify_previous();sources=inventory();report['source_sha256']=sources
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        for name in TESTS:
            cmd=[sys.executable,'-B','tests/review/test_'+name+'.py'];tick=time.perf_counter()
            try:
                result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                    text=True,encoding='utf-8',errors='replace',timeout=60)
                output,rc=result.stdout,result.returncode
            except subprocess.TimeoutExpired as error:
                output=error.stdout or b''
                if isinstance(output,bytes):output=output.decode('utf-8',errors='replace')
                output+='\nTIMEOUT_60S\n';rc=124
            (out/(name+'.log')).write_text(output,encoding='utf-8')
            row=dict(name=name,**suite_result(name,output,rc),command=cmd,return_code=rc,
                wall_seconds=time.perf_counter()-tick,log=name+'.log');report['tests'].append(row)
            write(out/'report.json',report);print(json.dumps(row),flush=True)
            require(row['all_tests_pass'],'component test failed '+name)
        require(sum(row['tests_run'] for row in report['tests'] if row['name'] in OLD_TESTS)==189,'old189 suite inventory differs')
        require(inventory()==sources,'source drift during preflight')
        report.update(status='PASS_STOP_COMPONENT_PREFLIGHT_NOT_PROCESS',prior_after=verify_previous());code=0
    except Exception as error:report.update(status='FAIL_STOP_PREFLIGHT',error=repr(error))
    report.update(tests_run=sum(row['tests_run'] for row in report['tests']),
        not_run_suites=[name for name in TESTS if name not in [row['name'] for row in report['tests']]],
        wall_seconds=time.perf_counter()-start,artifact_sha256={q.name:sha(q) for q in out.glob('*.log')})
    write(out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','artifact_sha256','tests')}),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

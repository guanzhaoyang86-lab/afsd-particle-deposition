"""Bounded regressions and source receipts before one larger-step experiment."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_COMPLETION_FREEZE.json'
TESTS=['structural_mobility','implicit_thermal','stage_anderson','continuation_state',
       'stage_preconditioner','stage_continuity_preconditioner','stage_structural_root',
       'process_thermal_trial','process_global_runner']
DOCUMENT_UPDATES={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,data):
    pending=path.with_suffix('.pending');pending.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8');pending.replace(path)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    out=ap.parse_args().output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh local output required')
    out.mkdir(parents=True);started=time.perf_counter();code=2
    r=dict(status='RUNNING',classification='GLOBAL_ROOT_PREFLIGHT_NOT_PROCESS_VALIDATION',tests_run=0,
        project_complete=False,production_authorized=False,accepted_steps=0,tests=[])
    try:
        frozen=json.loads(PRIOR.read_text(encoding='utf-8'));changed=[]
        for rel,h in frozen['files_sha256'].items():
            path=(ROOT/rel).resolve()
            if not path.is_relative_to(ROOT) or not path.is_file():raise ValueError('unsafe/missing frozen path '+rel)
            if sha(path)!=h:
                if rel not in DOCUMENT_UPDATES:raise ValueError('old evidence/source changed '+rel)
                changed.append(rel)
        from run_process_global import SOURCES
        watched=set(SOURCES)|{Path(__file__).relative_to(ROOT).as_posix()}
        watched|={'tests/review/test_'+name+'.py' for name in TESTS}
        # Include all directly imported backend modules, including frozen helpers.
        watched|={p.relative_to(ROOT).as_posix() for p in (ROOT/'src/afsd_backend/pysph').glob('*.py')}
        sources={rel:sha(ROOT/rel) for rel in sorted(watched)}
        r.update(source_sha256=sources,previous_freeze_sha256=sha(PRIOR),
            previous_files_checked=len(frozen['files_sha256']),authorized_document_changes=changed)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        for name in TESTS:
            cmd=[sys.executable,'-B','tests/review/test_'+name+'.py'];tick=time.perf_counter()
            try:
                completed=subprocess.run(cmd,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                    text=True,encoding='utf-8',errors='replace',timeout=60)
                output=completed.stdout;rc=completed.returncode
            except subprocess.TimeoutExpired as e:
                output=e.stdout or b''
                if isinstance(output,bytes):output=output.decode('utf-8',errors='replace')
                output+='\nTIMEOUT_60S\n';rc=124
            (out/(name+'.log')).write_text(output,encoding='utf-8')
            counts=re.findall(r'^Ran (\d+) tests? in ',output,re.M)
            good=rc==0 and len(counts)==1 and re.search(r'^OK\s*$',output,re.M) is not None
            row=dict(name=name,status='PASS' if good else 'FAIL',return_code=rc,tests_run=int(counts[0]) if len(counts)==1 else 0,
                wall_seconds=time.perf_counter()-tick,command=cmd,log=name+'.log')
            r['tests'].append(row);r['tests_run']+=row['tests_run'];write(out/'report.json',r)
            print(json.dumps(row),flush=True)
            if not good:raise ValueError('preflight failed '+name)
        for rel,h in sources.items():
            if sha(ROOT/rel)!=h:raise ValueError('source drift during tests '+rel)
        r['status']='PASS_GLOBAL_PREFLIGHT_NOT_PROCESS_VALIDATION';code=0
    except Exception as error:r.update(status='FAIL_PREFLIGHT',error=repr(error))
    r.update(wall_seconds=time.perf_counter()-started,artifact_sha256={p.name:sha(p) for p in out.glob('*.log')})
    write(out/'report.json',r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('source_sha256','tests','artifact_sha256')},allow_nan=False),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

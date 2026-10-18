"""One bounded component preflight. PASS never means AFSD process completion."""
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
PRIOR=ROOT/'docs/NONLINEAR_CONTINUATION_FREEZE.json'
DOCUMENT_UPDATES={'README.md','HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
COMPONENTS=[('finite_geometry','finite_tool_geometry'),('inlet','inlet_protocol'),
            ('thermal','implicit_thermal'),('structural_response','structural_mobility')]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,data):
    pending=path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    pending.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh project-local output directory required')
    out.mkdir(parents=True);start=time.perf_counter();code=2
    result=dict(classification='PROCESS_COMPONENT_PREFLIGHT_NOT_COUPLED_ADVANCE',
        status='RUNNING',project_complete=False,production_authorized=False,
        accepted_progress=False,new_coupled_steps=0,components=[],
        process_integration='NOT_RUN',material_transport='NOT_ESTABLISHED',
        permitted_next_action='JOINT_PROCESS_STEP_DEVELOPMENT_ONLY_NOT_PRODUCTION')
    write(out/'report.json',result)
    try:
        frozen=json.loads(PRIOR.read_text(encoding='utf-8'));changed=[]
        for rel,h in frozen['files_sha256'].items():
            path=(ROOT/rel).resolve()
            if not path.is_relative_to(ROOT) or not path.is_file():raise ValueError('missing/unsafe frozen path '+rel)
            if sha(path)!=h:
                if rel not in DOCUMENT_UPDATES:raise ValueError('old physics/evidence changed '+rel)
                changed.append(rel)
        upstream=ROOT/'third_party/PySPH'
        git=['git','-c','safe.directory='+upstream.as_posix(),'-C',str(upstream)]
        pin=subprocess.check_output(git+['rev-parse','HEAD'],text=True,timeout=10).strip()
        dirty=subprocess.check_output(git+['status','--porcelain','--untracked-files=no'],text=True,timeout=10).strip()
        if pin!=frozen['backend_commit'] or dirty:raise ValueError('fixed upstream changed')
        result.update(backend_commit=pin,previous_freeze_sha256=sha(PRIOR),
            previous_files_checked=len(frozen['files_sha256']),authorized_document_changes=changed)
        watched={Path(__file__).relative_to(ROOT).as_posix():sha(__file__)}
        for name,module in COMPONENTS:
            for rel in ('src/afsd_backend/pysph/'+module+'.py','tests/review/test_'+module+'.py'):
                watched[rel]=sha(ROOT/rel)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        for name,module in COMPONENTS:
            command=[sys.executable,'-B','tests/review/test_'+module+'.py'];tick=time.perf_counter()
            try:
                completed=subprocess.run(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                    text=True,encoding='utf-8',errors='replace',timeout=60)
                output=completed.stdout;returncode=completed.returncode
            except subprocess.TimeoutExpired as error:
                output=error.stdout or b''
                if isinstance(output,bytes):output=output.decode('utf-8',errors='replace')
                output+='\nTIMEOUT_60S\n';returncode=124
            (out/(name+'.log')).write_text(output,encoding='utf-8')
            counts=re.findall(r'^Ran (\d+) tests? in ',output,re.M)
            success=returncode==0 and len(counts)==1 and re.search(r'^OK\s*$',output,re.M) is not None
            row=dict(component=name,status='PASS_COMPONENT_ONLY' if success else 'FAIL_COMPONENT',
                return_code=returncode,tests_run=int(counts[0]) if len(counts)==1 else None,
                wall_seconds=time.perf_counter()-tick,command=command,log=name+'.log')
            result['components'].append(row);write(out/'report.json',result);print(json.dumps(row),flush=True)
            if not success:raise ValueError('component failed '+name)
        for rel,h in watched.items():
            if sha(ROOT/rel)!=h:raise ValueError('component source changed during preflight '+rel)
        result.update(status='PASS_COMPONENT_PREFLIGHT_PROCESS_STILL_BLOCKED',
            source_sha256=watched,total_distinct_component_tests=sum(r['tests_run'] for r in result['components']),
            remaining_gates=['joint_inertial_event_material_thermal_step','finite_feature_and_friction_coupling',
                'accepted_history_bound_inlet_commit','true_shoulder_passage_and_sustained_flow',
                'same_end_time_sensitivity','recomputed_restart_equivalence','bounded_cost_process_demonstration'])
        code=0
    except Exception as error:result.update(status='FAIL_PREFLIGHT',error=repr(error))
    result['wall_seconds']=time.perf_counter()-start
    result['artifact_sha256']={p.name:sha(p) for p in sorted(out.glob('*.log'))}
    write(out/'report.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','components','artifact_sha256')},allow_nan=False),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

"""Re-run bounded unit/RED checks. Does not launch the 2us case or change parameters."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
CASES=[
    ('contact_class',['tests/review/test_contact_class.py'],False),
    ('contact439_open_risk',['tests/review/test_contact439_open_risk.py'],True),
    ('contact454_precision',['tests/review/test_contact454_precision.py'],False),
    ('event_contact',['tests/review/test_contact_events.py'],False),
    ('closed_elastic_cycle',['tests/review/test_coupling_benchmarks.py','--cycle'],False),
    ('spatial_probes',['tests/review/test_spatial_probes.py'],False),
    ('event_solver_smoke',['tests/review/test_event_solver_smoke.py'],False),
    ('legacy_reversal_RED',['tests/review/test_coupling_benchmarks.py','--legacy-red'],True),
    ('boundary_patch_RED',['tests/review/test_spatial_probes.py','--boundary-red'],True),
    ('new_candidate_RED',['tests/review/test_new_candidate_red.py'],True),
    ('original_scientific',['scripts/run_scientific_tests.py'],True),
    ('original_backend',['scripts/run_physics_tests.py','--gate','all'],False),
    ('original_scaffold',['tests/run_tests.py'],False),
]
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,help='new, non-existing log directory')
    args=parser.parse_args()
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    results=[]
    for name,command,known_red in CASES:
        start=time.perf_counter()
        with (out/(name+'.log')).open('w',encoding='utf-8') as log:
            try:
                p=subprocess.run([sys.executable,'-B',*command],cwd=ROOT,stdout=log,
                    stderr=subprocess.STDOUT,timeout=60,
                    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1',
                         'AFSD_REVIEW_ARTIFACT_DIR':str(out/'artifacts')})
                code=p.returncode
            except subprocess.TimeoutExpired:
                code=124
        value=(out/(name+'.log')).read_text(encoding='utf-8')
        counts=re.findall(r'Ran (\d+) tests?',value)
        failures=re.findall(r'failures=(\d+)',value)
        errors=re.findall(r'errors=(\d+)',value)
        results.append(dict(name=name,command=[sys.executable,'-B',*command],
            exit_code=code,status='PASS' if code==0 else 'FAIL',
            known_red=known_red,tests_run=int(counts[-1]) if counts else None,
            failures=int(failures[-1]) if failures else 0,
            errors=int(errors[-1]) if errors else 0,wall_seconds=time.perf_counter()-start))
        print(name,results[-1]['status'],'exit',code,flush=True)
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((ROOT/'src/afsd_backend').rglob('*.py'))}
    report=dict(overall_status='FAIL' if any(r['exit_code'] for r in results) else 'PASS',
        meaning='RED tests remain FAIL, never xfail or relabelled as green',
        results=results,tests_run=sum(r['tests_run'] or 0 for r in results),
        failures=sum(r['failures'] for r in results),errors=sum(r['errors'] for r in results),
        source_hashes=hashes,real_2us_run='NOT_RUN_BY_THIS_ENTRY',production_authorized=False)
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return 1 if report['overall_status']=='FAIL' else 0
if __name__=='__main__':raise SystemExit(main())


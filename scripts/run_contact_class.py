"""Bounded fixed contact class; never starts an AFSD simulation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    results=[]
    for name in ('test_contact_class.py','test_contact439_open_risk.py',
                 'test_contact454_precision.py','test_contact_events.py'):
        start=time.perf_counter()
        with (out/(name+'.log')).open('w',encoding='utf-8') as stream:
            try:
                p=subprocess.run([sys.executable,'-B','tests/review/'+name],cwd=ROOT,
                    stdout=stream,stderr=subprocess.STDOUT,timeout=60,
                    env={**os.environ,'AFSD_REVIEW_ARTIFACT_DIR':str(out/'artifacts')})
                code=p.returncode
            except subprocess.TimeoutExpired:code=124
        results.append(dict(name=name,exit_code=code,status='PASS' if code==0 else 'FAIL',
                            wall_seconds=time.perf_counter()-start))
        print(name,results[-1]['status'],flush=True)
    hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((ROOT/'src/afsd_backend/pysph').glob('*.py'))}
    status='FAIL' if any(x['exit_code'] for x in results) else 'PASS'
    (out/'summary.json').write_text(json.dumps(dict(status=status,results=results,
        source_hashes=hashes,real_afsd='NOT_RUN',production_authorized=False),indent=2)+'\n')
    return int(status!='PASS')


if __name__=='__main__':raise SystemExit(main())

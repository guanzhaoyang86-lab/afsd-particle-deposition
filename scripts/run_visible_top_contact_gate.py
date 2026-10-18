"""Immutable execution receipt for TOP components only; no process state loaded."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--name',required=True);args=parser.parse_args()
    if not args.name.isidentifier():raise ValueError('simple unique receipt name required')
    out=ROOT/'outputs/visible_top_components_20261004'/args.name;out.mkdir(parents=True,exist_ok=False)
    paths=[Path(__file__),ROOT/'scripts/test_visible_top_contact.py',ROOT/'docs/VISIBLE_TOP_COMPONENT_CONTRACT_20261004.md',
        ROOT/'scripts/batched_contact.py',ROOT/'scripts/contact_shell_audit.py',ROOT/'scripts/contact_shell_finite_enclosure.py',
        ROOT/'src/afsd_backend/pysph/discrete_contact.py',ROOT/'src/afsd_backend/pysph/quasistatic_contact.py']
    frozen={p.relative_to(ROOT).as_posix():sha(p) for p in paths}
    command=[sys.executable,'-B',str(ROOT/'scripts/test_visible_top_contact.py')]
    begin=time.perf_counter();timed_out=False
    with (out/'tests.log').open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=120,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124;timed_out=True
    changed=[p for p,h in frozen.items() if sha(ROOT/p)!=h]
    log=(out/'tests.log').read_text(encoding='utf-8');test_rows=[line for line in log.splitlines() if line.startswith('test_')]
    passed=code==0 and not changed and 'Ran 15 tests' in log and log.rstrip().endswith('OK') and len(test_rows)==15
    receipt=dict(status='PASS_TOP_CONTACT_COMPONENTS_NOT_PROCESS_VALIDATION' if passed else 'FAIL_TOP_CONTACT_COMPONENT_GATE',
        command=command,exit_code=code,tests_run=len(test_rows),test_results=test_rows,wall_seconds=time.perf_counter()-begin,
        timeout_seconds=120,timed_out=timed_out,source_sha256=frozen,changed_sources=changed,
        output_sha256={(out/'tests.log').relative_to(ROOT).as_posix():sha(out/'tests.log')},
        physical_solver_evaluations=0,nnps_evaluations=0,accepted_time_advanced_s=0.,
        actual_failed_candidate_shell_or_visibility_audits='NOT_RUN_BY_COMPONENT_GATE')
    (out/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'],exit_code=code,summary=str(out/'summary.json'))))
    return 0 if passed else 2


if __name__=='__main__':raise SystemExit(main())

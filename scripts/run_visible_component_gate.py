"""Bind unchanged geometry and final graph-aware search to actual test logs."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time
import run_discrete_contact_process as core


def main():
    root=core.ROOT;target=root/'outputs/visible_support_20261004/component_receipt.json'
    logs=root/'logs/visible_support_20261004/components_v1'
    if target.exists() or logs.exists():raise ValueError('refuse receipt/log overwrite')
    logs.mkdir();files=[Path(__file__).resolve(),root/'scripts/visible_support.py',
        root/'scripts/test_visible_support.py',root/'scripts/check_visibility_equilibrium.py',
        root/'scripts/visibility_equilibrium.py',root/'scripts/secant_equilibrium.py',
        root/'scripts/audit_solid_occluded_neighbors.py']
    frozen=dict(core.frozen_sources(),**{path.relative_to(root).as_posix():core.sha(path) for path in files})
    rows=[];start=time.perf_counter()
    for script,count in [('test_visible_support.py',10),('check_visibility_equilibrium.py',3)]:
        command=[sys.executable,'-B',str(root/'scripts'/script)];log=logs/(script+'.log')
        begin=time.perf_counter()
        with log.open('w',encoding='utf-8') as out:
            code=subprocess.run(command,cwd=root,stdout=out,stderr=subprocess.STDOUT,timeout=60,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        text=log.read_text(encoding='utf-8')
        passed=code==0 and f'Ran {count} tests' in text and text.rstrip().endswith('OK')
        rows.append(dict(command=command,exit_code=code,expected_test_count=count,passed=passed,
            log=log.relative_to(root).as_posix(),log_sha256=core.sha(log),wall_seconds=time.perf_counter()-begin))
    changed=[name for name,digest in frozen.items() if core.sha(root/name)!=digest]
    status='PASS_VISIBLE_GEOMETRY_AND_SEARCH_COMPONENTS_ONLY' if not changed and all(r['passed'] for r in rows) else 'FAIL_COMPONENT_GATE'
    core.write(target,dict(status=status,source_sha256=frozen,changed_sources=changed,tests=rows,
        count=sum(r['expected_test_count'] for r in rows),physical_solver_evaluations=0,
        accepted_time_advanced_s=0.,wall_seconds=time.perf_counter()-start))
    print(json.dumps(dict(status=status,tests=13,output=str(target))))
    return 0 if status.startswith('PASS') else 2


if __name__=='__main__':raise SystemExit(main())

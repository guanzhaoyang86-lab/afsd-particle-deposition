"""Record standalone endpoint/path geometry component receipts without physics."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/coupled_cost_20261004/visible_support_geometry_components.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if OUT.exists():raise RuntimeError('refuse geometry component receipt overwrite')
    rows=[]
    for filename,expected_count in [('test_visible_support.py',10),('test_visible_support_path.py',11)]:
        started=time.perf_counter()
        result=subprocess.run([sys.executable,str(ROOT/'scripts'/filename)],cwd=ROOT,text=True,capture_output=True,timeout=60)
        rows.append(dict(script='scripts/'+filename,expected_component_count=expected_count,exit_code=result.returncode,
            stdout=result.stdout,stderr=result.stderr,wall_seconds=time.perf_counter()-started,
            pass_receipt=bool(result.returncode==0 and f'Ran {expected_count} tests' in result.stderr and result.stderr.rstrip().endswith('OK'))))
    sources=[Path(__file__),ROOT/'scripts/visible_support.py',ROOT/'scripts/visible_support_path.py',
        ROOT/'scripts/test_visible_support.py',ROOT/'scripts/test_visible_support_path.py',
        ROOT/'scripts/audit_solid_occluded_neighbors.py',ROOT/'src/afsd_backend/pysph/quasistatic_contact.py']
    report=dict(status='PASS_GEOMETRY_COMPONENTS' if all(row['pass_receipt'] for row in rows) else 'FAIL_GEOMETRY_COMPONENTS',
        component_test_count=sum(row['expected_component_count'] for row in rows),rows=rows,
        physical_solver_evaluations=0,nnps_evaluations=0,
        source_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in sources},
        limitations=['These are endpoint/path geometry component tests, not process or physical validation.'])
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(output=str(OUT),status=report['status'],component_test_count=report['component_test_count'])))
    if report['status']!='PASS_GEOMETRY_COMPONENTS':raise SystemExit(1)


if __name__=='__main__':main()

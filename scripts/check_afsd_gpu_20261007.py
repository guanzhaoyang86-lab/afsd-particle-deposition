"""Execute and bind GPU component tests before any AFSD force-root stage."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output=ROOT/'outputs/gpu_acceleration_20261007/components_v1'
    output.mkdir(parents=True, exist_ok=False)
    files=list((ROOT/'scripts').glob('*afsd_gpu*20261007.py'))+[ROOT/'docs/GPU_ACCELERATION_CONTRACT_20261007.md']
    files += [p for p in (ROOT/'.envs/gpu_overlay_20261007').rglob('*') if p.is_file() and (p.suffix in ('.py','.pyd') or p.name=='METADATA')]
    mapping={p.relative_to(ROOT).as_posix():sha(p) for p in files}
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH='')
    checks=[]; began=time.perf_counter()
    commands=[('linalg',['scripts/afsd_gpu_linalg_20261007.py',str(output/'linalg.json')]),
        ('thermal',['scripts/test_afsd_gpu_thermal_20261007.py','--gpu','--output',str(output/'thermal.json')]),
        ('equilibrium',['scripts/test_afsd_gpu_equilibrium_20261007.py','-v']),
        ('runner',['scripts/test_afsd_gpu_runner_20261007.py','-v'])]
    for name, args in commands:
        cmd=[sys.executable,'-B',*args]
        with (output/(name+'.log')).open('x',encoding='utf-8') as stream:
            result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=60.)
        checks.append(dict(name=name,exit_code=result.returncode,command=cmd))
    stable=all(sha(ROOT/p)==h for p,h in mapping.items())
    passed=stable and all(row['exit_code']==0 for row in checks)
    report=dict(status='PASS_GPU_COMPONENTS' if passed else 'FAIL_GPU_COMPONENTS',checks=checks,
        source_sha256=mapping,source_inputs_still_match=stable,wall_seconds=time.perf_counter()-began,
        artifacts_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in output.iterdir() if p.is_file()},
        physical_time_advanced_s=0.,accepted=False)
    (output/'manifest.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(dict(status=report['status'],wall_s=report['wall_seconds'],checks=checks)))
    return 0 if passed else 1


if __name__=='__main__':
    raise SystemExit(main())

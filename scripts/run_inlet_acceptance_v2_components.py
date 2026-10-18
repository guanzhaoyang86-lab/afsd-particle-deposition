"""Save a bounded manufactured-component receipt; no physical run."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/inlet_acceptance_v2_20261004'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    receipt_path=OUT/'component_receipt.json'
    if OUT.exists():raise RuntimeError('refuse existing inlet component output directory')
    if (ROOT/'src/afsd_backend/pysph/inlet_acceptance_v2.py').exists():
        raise RuntimeError('new component must remain outside frozen src source scan')
    sources=[Path(__file__),ROOT/'scripts/inlet_acceptance_v2.py',ROOT/'tests/review/test_inlet_acceptance_v2.py',
        ROOT/'tests/review/test_inlet_protocol.py',ROOT/'src/afsd_backend/pysph/inlet_protocol.py',
        ROOT/'src/afsd_backend/pysph/particles.py',ROOT/'src/afsd_backend/pysph/interval_acceptance.py',
        ROOT/'src/afsd_backend/pysph/output.py',ROOT/'docs/INLET_ACCEPTANCE_V2_COMPONENT_CONTRACT_20261004.md']
    before={path.relative_to(ROOT).as_posix():sha(path) for path in sources}
    command=[sys.executable,str(ROOT/'tests/review/test_inlet_acceptance_v2.py')]
    started=time.perf_counter()
    result=subprocess.run(command,cwd=ROOT,text=True,capture_output=True,timeout=60)
    seconds=time.perf_counter()-started
    after={path.relative_to(ROOT).as_posix():sha(path) for path in sources}
    passed=result.returncode==0 and 'Ran 11 tests' in result.stderr and result.stderr.rstrip().endswith('OK') and before==after
    OUT.mkdir(parents=True)
    stdout=OUT/'component.stdout.log';stderr=OUT/'component.stderr.log'
    stdout.write_text(result.stdout,encoding='utf-8');stderr.write_text(result.stderr,encoding='utf-8')
    report=dict(status='PASS_INLET_ACCEPTANCE_V2_COMPONENT_ONLY' if passed else 'FAIL_INLET_ACCEPTANCE_V2_COMPONENT_ONLY',
        classification='MANUFACTURED_TRANSACTION_COMPONENTS_NOT_AFSD',command=command,cwd=str(ROOT),
        exit_code=result.returncode,wall_seconds=seconds,test_count=11,physical_solver_evaluations=0,
        nnps_evaluations=0,force_evaluations=0,thermal_evaluations=0,new_process_accepted_time_s=0.,
        manufactured_pa_updates=True,source_sha256=after,sources_unchanged_during_test=before==after,
        log_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in (stdout,stderr)},
        limitations=['Two insertion/restart chains use manufactured PA updates and receipt fixtures, not a real inlet benchmark.',
            'Actual process source-work and total cumulative process energy closure remain external and unimplemented here.',
            'Diagnostic scope permits future real benchmark chaining but cannot be upgraded to process scope.',
            'No existing solver, source scanner, controller, accepted checkpoint, or old failure is modified.'])
    receipt_path.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(output=str(receipt_path),status=report['status'],test_count=11,exit_code=result.returncode,wall_seconds=seconds)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()

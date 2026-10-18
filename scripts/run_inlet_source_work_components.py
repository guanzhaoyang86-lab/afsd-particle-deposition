"""Bounded receipts for manufactured source-work identities; no physics run."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/inlet_source_work_20261005'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if OUT.exists():raise RuntimeError('refuse existing source-work component output')
    sources=[Path(__file__),ROOT/'scripts/inlet_source_work.py',ROOT/'scripts/test_inlet_source_work.py',
        ROOT/'scripts/inlet_acceptance_v2.py',ROOT/'src/afsd_backend/pysph/inlet_protocol.py',
        ROOT/'src/afsd_backend/pysph/particles.py',ROOT/'src/afsd_backend/pysph/interval_acceptance.py',
        ROOT/'src/afsd_backend/pysph/mechanics.py',ROOT/'src/afsd_backend/pysph/quasistatic_process.py',
        ROOT/'src/afsd_backend/pysph/implicit_thermal.py',ROOT/'src/afsd_backend/pysph/material.py',
        ROOT/'docs/INLET_SOURCE_WORK_COMPONENT_CONTRACT_20261005.md']
    before={path.relative_to(ROOT).as_posix():sha(path) for path in sources}
    command=[sys.executable,str(ROOT/'scripts/test_inlet_source_work.py')]
    started=time.perf_counter();result=subprocess.run(command,cwd=ROOT,text=True,capture_output=True,timeout=60)
    elapsed=time.perf_counter()-started
    passed=result.returncode==0 and 'Ran 11 tests' in result.stderr and result.stderr.rstrip().endswith('OK')
    OUT.mkdir(parents=True)
    stdout=OUT/'component.stdout.log';stderr=OUT/'component.stderr.log'
    stdout.write_text(result.stdout,encoding='utf-8');stderr.write_text(result.stderr,encoding='utf-8')
    artifacts=[stdout,stderr];cases=[]
    if passed:
        from test_inlet_source_work import fixture,evaluate
        from inlet_source_work import source_work_record
        for name,options in [('heated_elastic',{}),('negative_work_unload',dict(unload=True)),('retained_raw_defect',dict(raw_defect=.012))]:
            calculated=evaluate(fixture(**options));path=OUT/(name+'.json')
            path.write_text(json.dumps(calculated,indent=2,allow_nan=False),encoding='utf-8');artifacts.append(path)
            terms=source_work_record(calculated,path)
            work_path=OUT/(name+'_v2_work.json');work_path.write_text(json.dumps(terms,indent=2),encoding='utf-8');artifacts.append(work_path)
            cases.append(dict(name=name,classification='MANUFACTURED_ENERGY_IDENTITY_NOT_COUPLED_PHYSICS',
                source_work=terms,energy=calculated['open_energy']))
    after={path.relative_to(ROOT).as_posix():sha(path) for path in sources};passed=passed and before==after
    report=dict(status='PASS_INLET_SOURCE_WORK_COMPONENT_ONLY' if passed else 'FAIL_INLET_SOURCE_WORK_COMPONENT_ONLY',
        command=command,cwd=str(ROOT),exit_code=result.returncode,wall_seconds=elapsed,test_count=11,
        physical_solver_evaluations=0,nnps_evaluations=0,force_evaluations=0,thermal_evaluations=0,
        source_sha256=after,sources_unchanged_during_test=before==after,
        artifact_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in artifacts},cases=cases,
        development_note='The initial independent fixture storage reference omitted V=m/rho; corrected from 2/3 to 1/3 J. Solver formula and acceptance thresholds were unchanged.',
        limitations=['No real insertion, equilibrium solve, thermal solve, or process acceptance was run.',
            'Stress-free reference-density inventory birth and explicitly zero exterior traction only.',
            'Per-field consistency roundoff bounds are not process energy acceptance thresholds.',
            'A cumulative open-process ledger and long-bar-reference comparison remain necessary.'])
    receipt=OUT/'component_receipt.json';receipt.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(output=str(receipt),status=report['status'],test_count=11,exit_code=result.returncode,wall_seconds=elapsed)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()

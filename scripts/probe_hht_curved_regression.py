"""One saved RED interval versus the named curved variant, never published."""
from pathlib import Path
import argparse
import gzip
import json
import os
import threading
import time
import numpy as np
import probe_hht_curved as driver
from afsd_backend.pysph import hht_coupled as old,hht_curved_coupled as new
from afsd_backend.pysph.output import restore,checkpoint
from afsd_backend.pysph.updated_trial import particle_hash


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for flag in ('output','preflight','reference'):ap.add_argument('--'+flag,type=Path,required=True)
    a=ap.parse_args();out=a.output.resolve()
    driver.require(out.is_relative_to(driver.ROOT) and not out.exists(),'fresh project output required')
    out.mkdir(parents=True);began=time.perf_counter();sources=driver.inventory()
    report=dict(status='NOT_RUN_GATES_PENDING',classification='PRIVATE_CURVED_SAVED_RED_REGRESSION',
        accepted=False,restart_allowed=False,new_afsd_steps=0,production_authorized=False,
        source_sha256=sources,script_sha256=driver.sha(Path(__file__)))
    timer=threading.Timer(55.,lambda:os._exit(124));timer.daemon=True;timer.start();code=2
    try:
        ref,receipts=driver.receipts(a.preflight,a.reference,sources)
        gate=driver.read(a.preflight)
        driver.require(gate['source_sha256'].get('scripts/probe_hht_curved_regression.py')==driver.sha(Path(__file__)),'missing regression source receipt')
        report['receipts']=receipts
        parent=driver.ROOT/'outputs/process_slow_20261003/S0_initial'
        receipt=driver.read(parent/'report.json')
        for name,digest in receipt['source_sha256'].items():driver.require(driver.sha(driver.ROOT/name)==digest,'old source drift')
        for name,digest in receipt['artifact_sha256'].items():driver.require(driver.sha(parent/name)==digest,'old evidence drift')
        history=driver.read(parent/'trajectory/accepted_history.json')
        driver.require(len(history)==8 and history[-1]['time_s']==2e-6,'original RED parent required')
        pa,state=restore(parent/'trajectory'/history[-1]['checkpoint'])
        initial,p=driver.setup();pa,state=old.restore_checkpoint(pa,state,p)
        before=particle_hash(pa);bc=driver.boundary_for(p);dt=5e-7
        report.update(parent_report_sha256=driver.sha(parent/'report.json'),parent_state_hash=state['state_hash'],dt_s=dt)
        try:old.advance(pa,state,p,bc,dt)
        except old.CoupledStageError as error:
            driver.require(str(error)=='contact turning-traction quadrature law defect','different old failure')
            report['baseline_red']=dict(error=str(error),diagnostics=error.diagnostics)
        else:raise ValueError('original RED unexpectedly passed')
        migrated,migration=new.migrate_diagnostic(pa,state,p);report['migration']=migration
        candidate,next_state,diag=new.advance(pa,migrated,p,bc,dt,driver.slow_gap_certificate)
        driver.checks(candidate,pa,initial,next_state,diag,p)
        driver.require(particle_hash(pa)==before,'accepted parent mutated')
        driver.require(driver.inventory()==sources,'source drift')
        report.update(status='PASS_PRIVATE_SAVED_RED_REGRESSION_NOT_TRAJECTORY',
            candidate_time_s=next_state['time_s'],evaluations=diag['evaluations'],residuals=diag['residuals'],
            candidate_particle_hash=particle_hash(candidate),
            maximum_coulomb_defect=diag['maximum_sampled_coulomb_relative_defect'],
            raw_friction_J=float(diag['contact']['friction_raw_J'].sum()),
            correction_impulse_Ns=diag['contact']['delta_impulse_Ns'],
            kinetic_identity_max_J=float(np.max(np.abs(diag['kinetic']['kinetic_identity_residual_J']))),
            parent_unchanged=True)
        checkpoint(out/'PRIVATE_NOT_RESTARTABLE.npz',candidate,dict(classification=report['classification'],
            accepted=False,restart_allowed=False,diagnostic_state=next_state))
        with gzip.open(out/'private_diagnostic.json.gz','wt',encoding='utf-8') as stream:
            json.dump(driver.clean(diag),stream,allow_nan=False)
        code=0
    except Exception as error:
        report.update(status='FAIL_PRIVATE_CURVED_REGRESSION',error=repr(error),
            diagnostic=getattr(error,'diagnostics',None))
    finally:timer.cancel()
    report.update(wall_seconds=time.perf_counter()-began,
        artifact_sha256={p.name:driver.sha(p) for p in out.iterdir() if p.is_file()})
    driver.write(out/'report.json',report)
    print(json.dumps(driver.clean({k:v for k,v in report.items() if k not in ('source_sha256','baseline_red','correction_impulse_Ns','diagnostic')})),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

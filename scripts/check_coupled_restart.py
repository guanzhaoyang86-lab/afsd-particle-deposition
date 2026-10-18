"""Replay two already accepted steps from real coupled checkpoints, not extend a run."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from run_coupled_afsd import ROOT,SOURCES,settings,sha,fingerprint,write_json,clean
from afsd_backend.pysph.coupled_solver import IntegratedSolver


def check(case):
    started=time.perf_counter();case=Path(case)
    manifest=json.loads((case/'manifest.json').read_text(encoding='utf-8'))
    paths=list(SOURCES);cfg,p,resolved=settings()
    if 'scripts/run_coupled_short_window.py' in manifest['source_sha256']:
        from run_coupled_short_window import short_settings,EXTRA
        cfg,p,resolved=short_settings();paths+=EXTRA
    source={s:sha(ROOT/s) for s in paths}
    if manifest['source_sha256']!=source or manifest['configuration_sha256']!=fingerprint(resolved):
        raise ValueError('source/configuration changed')
    history=json.loads((case/'history.json').read_text(encoding='utf-8'))
    if sha(case/'history.json')!=manifest['artifact_sha256']['history.json']:raise ValueError('history artifact changed')
    if len(history)!=manifest['accepted_steps']+1:raise ValueError('accepted sequence incomplete')
    cases=[]
    for name in ('nonzero_accepted.npz','plastic_accepted.npz'):
        path=case/name
        if not path.exists():
            cases.append(dict(checkpoint=name,status='NOT_RUN',reason='no such accepted checkpoint'));continue
        if sha(path)!=manifest['artifact_sha256'][name]:raise ValueError('checkpoint artifact changed')
        s=IntegratedSolver.from_checkpoint(path,p,source_hashes=source,configuration_hash=fingerprint(resolved))
        if s.count+2>manifest['accepted_steps']:
            cases.append(dict(checkpoint=name,status='NOT_RUN',reason='fewer than two stored successor steps'));continue
        maxima={};tested=[];origin=s.count
        for offset in range(3):
            ref=history[s.count];actual=clean(s.record())
            if actual['time_s']!=ref['time_s']:raise ValueError('restart time mismatch')
            for key in ('particle_id','component_id','material_id','tracer_type','wholebar_tag'):
                if actual[key]!=ref[key]:raise ValueError('restart identity '+key)
            for key,tol in [('position_m',1e-12),('velocity_m_s',1e-12),('density_kg_m3',1e-12),
                            ('temperature_K',1e-12),('deviatoric_stress_Pa',1e-6),('ep',1e-12)]:
                diff=float(np.max(np.abs(np.asarray(actual[key])-np.asarray(ref[key]))))
                if not np.isfinite(diff) or diff>tol:raise ValueError('restart state '+key)
                maxima[key]=max(maxima.get(key,0.),diff)
            for key,value in ref['cumulative'].items():
                if abs(actual['cumulative'][key]-value)>1e-15:raise ValueError('restart ledger '+key)
            tested.append(s.count)
            if offset<2:s.advance(max_dt=manifest['accepted_dt_s'][s.count])
        cases.append(dict(checkpoint=name,status='PASS_ACCEPTED_TRAJECTORY_RESTART',
            source_step=origin,tested_steps=tested,maximum_field_difference=maxima))
    return dict(status='PASS_AVAILABLE_RESTARTS_NOT_PHYSICS_REFERENCE' if any(c['status'].startswith('PASS') for c in cases) else 'NOT_RUN',
        production_authorized=False,cases=cases,wall_seconds=time.perf_counter()-started,
        manifest_sha256=sha(case/'manifest.json'),checker_sha256=sha(Path(__file__)),
        note='Uses previously accepted trajectories; no original files modified and no extension of physical end-time.')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);a=parser.parse_args()
    if a.output.exists():raise FileExistsError('fresh report path required')
    try:report=check(a.case);code=0 if report['status'].startswith('PASS') else 2
    except Exception as e:report=dict(status='FAIL_RESTART',error=repr(e),production_authorized=False);code=2
    a.output.parent.mkdir(parents=True,exist_ok=True);write_json(a.output,report)
    print(json.dumps(clean(report),allow_nan=False));return code


if __name__=='__main__':raise SystemExit(main())

"""Offline component wiring using saved real intervals; NO coupled advance."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
import time
import numpy as np
from pysph.solver.utils import load
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.implicit_thermal import advance_implicit_thermal

CASE=ROOT/'outputs/nonlinear_continuation_20261003/continuation'
SOURCES=['scripts/replay_process_thermal.py','src/afsd_backend/pysph/implicit_thermal.py',
         'src/afsd_backend/pysph/coupled_thermal.py','src/afsd_backend/pysph/thermal.py',
         'docs/PROCESS_IMPLICIT_THERMAL_COMPONENT.md']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def replay():
    start=time.perf_counter();manifest=json.loads((CASE/'manifest.json').read_text(encoding='utf-8'))
    if manifest['status']!='COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION':
        raise ValueError('completed registered source required')
    watched={CASE/'manifest.json':sha(CASE/'manifest.json')}
    for rel,h in manifest['source_sha256'].items():
        path=ROOT/rel
        if sha(path)!=h:raise ValueError('old source drift '+rel)
        watched[path]=h
    for rel,h in manifest['artifact_sha256'].items():
        path=(CASE/rel).resolve()
        if not path.is_relative_to(CASE) or sha(path)!=h:raise ValueError('old artifact drift '+rel)
        watched[path]=h
    parent=ROOT/manifest['parent_reference']['path']
    if sha(parent)!=manifest['parent_reference']['sha256']:raise ValueError('old parent drift')
    watched[parent]=sha(parent)
    p=json.loads((CASE/'resolved_config.json').read_text(encoding='utf-8'))['parameters']
    history=json.loads((CASE/'accepted_history.json').read_text(encoding='utf-8'))
    source_hashes={rel:sha(ROOT/rel) for rel in SOURCES}
    results=[];before=parent
    for entry in history:
        old=load(str(before))['arrays']['solid'];directory=(CASE/entry['snapshot']).parent
        with gzip.open(directory/'accepted_trial.json.gz','rt',encoding='utf-8') as stream:e=json.load(stream)['evaluation']
        mid=e['midpoint'];boundary=e['boundary'];n=len(old.x)
        x=np.asarray(mid['position_m']);i=np.asarray(mid['i'],int);j=np.asarray(mid['j'],int)
        g=np.asarray(mid['g']);rho=np.asarray(mid['density_kg_m3']);dt=e['dt_s']
        wall=np.zeros(n);bath=np.zeros(n);bound=np.zeros(n);qfromtraces=np.zeros(n)
        for trace in boundary['traces']:
            a=trace['particle_index']
            if int(trace['particle_id'])!=int(old.particle_id[a]):raise ValueError('trace ID correspondence')
            for face in trace['contact_rows']:
                conductance=p['contact_conductance']*p['dx']**2
                wall[a]+=conductance;bath[a]+=conductance*face['wall_temperature_K']
            for interval in trace['result']['intervals']:
                bound[a]+=sum(interval['friction_rounding_bound_J'])
                qfromtraces[a]+=sum(interval['friction_heat_raw_J'])
        raw=np.asarray(boundary['friction_raw_J'])
        if not np.allclose(qfromtraces,raw,rtol=1e-12,atol=1e-25):raise ValueError('friction source trace mismatch')
        if not np.allclose(wall,boundary['wall_row_W_K'],rtol=1e-13,atol=1e-15):raise ValueError('wall conductance mismatch')
        result=advance_implicit_thermal(old.m,rho,old.temperature,p['conductivity'],x,i,j,g,dt,p['cp'],
            e['material']['plastic_work_J'],raw,bound,wall,bath,
            p['plastic_heat_fraction'],p['heat_partition'],p['room'],p['melt'])
        # Independent pair/temperature ledger algebra, not returned ledger fields.
        rij=x[i]-x[j]
        coefficient=(old.m[i]/rho[i])*(old.m[j]/rho[j])*(2*p['conductivity'])*(-np.sum(rij*g,axis=1))/np.sum(rij**2,axis=1)
        if not np.allclose(coefficient,result['pair_conductance_W_K'],rtol=1e-13,atol=1e-20):raise ValueError('pair coefficient mismatch')
        rows=np.zeros(n);np.add.at(rows,i,coefficient);np.add.at(rows,j,coefficient)
        if not np.allclose(rows,e['thermal']['conduction_row_W_K'],rtol=1e-13,atol=1e-16):raise ValueError('old conduction row differs')
        t=result['temperature'];pair=dt*coefficient*(t[j]-t[i]);pairnode=np.zeros(n)
        np.add.at(pairnode,i,pair);np.add.at(pairnode,j,-pair)
        contact=dt*(bath-wall*t)
        sources=p['plastic_heat_fraction']*np.asarray(e['material']['plastic_work_J'])+p['heat_partition']*raw
        measured=old.m*p['cp']*(t-old.temperature);reference=sources+contact+pairnode
        limit=1e-12+1e-10*float(np.abs(reference).sum())
        if float(np.max(np.abs(measured-reference)))>limit:raise ValueError('independent per-particle energy balance')
        results.append(dict(source_stage_index=entry['stage_index'],dt_s=dt,particle_count=n,pair_count=len(i),
            former_explicit_bound_s=result['former_explicit_dt_bound_s'],
            max_difference_from_original_explicit_temperature_K=float(np.max(np.abs(t-np.asarray(e['thermal']['temperature'])))),
            independent_max_heat_residual_J=float(np.max(np.abs(measured-reference))),
            independent_total_heat_residual_J=float(np.sum(measured-sources-contact)),
            boundary_heat_into_material_J=float(contact.sum()),temperature_K=t.tolist(),
            temperature_min_K=float(t.min()),temperature_max_K=float(t.max())))
        before=CASE/entry['snapshot']
    for path,h in watched.items():
        if sha(path)!=h:raise ValueError('source input changed during replay')
    for rel,h in source_hashes.items():
        if sha(ROOT/rel)!=h:raise ValueError('replay source changed during run')
    return dict(status='PASS_OFFLINE_IMPLICIT_HEAT_REPLAY_NOT_COUPLED_ADVANCE',accepted_progress=False,
        production_authorized=False,restart_allowed=False,old_physical_end_time_s=manifest['accepted_time_s'],
        new_coupled_steps=0,stages=results,source_sha256=source_hashes,
        input_sha256={path.relative_to(ROOT).as_posix():h for path,h in watched.items()},
        wall_seconds=time.perf_counter()-start,
        interpretation='Each source interval replayed independently; new T NEVER fed to the next saved interval.')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    target=a.output.resolve()
    if not target.is_relative_to(ROOT) or target.exists():raise ValueError('fresh project-local output required')
    target.parent.mkdir(parents=True,exist_ok=True)
    try:r=replay();code=0
    except Exception as error:r=dict(status='FAIL_OFFLINE_HEAT_REPLAY',error=repr(error),accepted_progress=False);code=2
    target.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in r.items() if k not in ('input_sha256','stages')},allow_nan=False))
    if code==0:print(json.dumps([{k:v for k,v in s.items() if k!='temperature_K'} for s in r['stages']],allow_nan=False))
    return code


if __name__=='__main__':raise SystemExit(main())

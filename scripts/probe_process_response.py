"""One bounded frozen-elastic actual-cloud probe, NOT an AFSD time step."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
import threading
import time
import numpy as np
from pysph.solver.utils import load
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.structural_mobility import build_frozen_elastic
from afsd_backend.pysph.mechanics import bulk_modulus,corrected_operators,internal_force,positions
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.coupled_boundary import normal_snapshot
from afsd_backend.pysph.axial_actuator import control
from afsd_backend.pysph.application import memory_bytes

CASE=ROOT/'outputs/nonlinear_continuation_20261003/continuation'
SOURCES=['scripts/probe_process_response.py','src/afsd_backend/pysph/structural_mobility.py',
         'docs/PROCESS_RESPONSE_PROBE_CONTRACT.md']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p,data):
    pending=p.with_suffix('.pending')
    pending.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8');pending.replace(p)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    out=args.output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh project-local output directory required')
    out.mkdir(parents=True);start=time.perf_counter()
    r=dict(classification='FROZEN_ELASTIC_ACTUAL_CLOUD_COST_NOT_PROCESS_SOLVER',status='RUNNING',
        accepted_progress=False,new_coupled_steps=0,production_authorized=False,restart_allowed=False)
    write(out/'report.json',r)
    def deadline():
        try:write(out/'HARD_TIMEOUT.json',dict(status='NOT_COMPLETED_50S_NO_PHYSICAL_ADVANCE',accepted_progress=False))
        finally:os._exit(124)
    timer=threading.Timer(50.,deadline);timer.daemon=True;timer.start()
    try:
        manifest=json.loads((CASE/'manifest.json').read_text(encoding='utf-8'))
        if manifest['status']!='COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION':raise ValueError('registered source required')
        watched={CASE/'manifest.json':sha(CASE/'manifest.json')}
        for rel,h in manifest['source_sha256'].items():
            path=ROOT/rel
            if sha(path)!=h:raise ValueError('old source drift '+rel)
            watched[path]=h
        for rel,h in manifest['artifact_sha256'].items():
            path=(CASE/rel).resolve()
            if not path.is_relative_to(CASE) or sha(path)!=h:raise ValueError('old artifact drift '+rel)
            watched[path]=h
        path=CASE/manifest['last_accepted_snapshot'];data=load(str(path));pa=data['arrays']['solid'];state=data['solver_data']
        if not state['accepted'] or state['particle_state_sha256']!=particle_hash(pa):raise ValueError('accepted state mismatch')
        before=particle_hash(pa);p=state['parameters'];sources={rel:sha(ROOT/rel) for rel in SOURCES}
        bulk=bulk_modulus(p['shear'],p['poisson']);tick=time.perf_counter();op=build_frozen_elastic(pa,p['shear'],bulk)
        r.update(assembly_wall_seconds=time.perf_counter()-tick,particle_count=len(pa.x),pair_count=op.pair_count,
            gradient_shape=list(op.gradient_map.shape),gradient_bytes=op.gradient_map.nbytes,stiffness_bytes=op.stiffness.nbytes)
        command=control(pa.z0,state['t'],0.,p['bar_length'],p['actuator_depth'],p['feed'],p['release_time'])
        snapshot=normal_snapshot(positions(pa),pa.particle_id,state['t'],p);rows=snapshot['rows']
        active=np.flatnonzero(snapshot['normal_force_N']>0);gmap=np.zeros((len(active),3*len(pa.x)))
        for a,k in enumerate(active):
            at=int(rows['particle_index'][k]);gmap[a,3*at:3*at+3]=rows['normal'][k]
        r['active_normal_contact_rows']=len(active);r['responses']=[]
        for dt in (5e-8,.001):
            tick=time.perf_counter();response=op.response_operator(dt,command['free_axes'])
            setup=time.perf_counter()-tick;tick=time.perf_counter();w=response.contact_response(gmap,at='endpoint')
            ids=rows['particle_id'][active];cross=ids[:,None]!=ids[None,:]
            item=dict(dt_s=dt,matrix_setup_wall_seconds=setup,contact_response_wall_seconds=time.perf_counter()-tick,
                diagnostics=response.diagnostics,response_shape=list(w.shape),
                max_cross_particle_response_per_kg=float(np.max(np.abs(w[cross]),initial=0.)),
                response_asymmetry_per_kg=float(np.max(np.abs(w-w.T),initial=0.)),
                interpretation='UNPRESTRESSED_ELASTIC_ENDPOINT_RESPONSE_NOT_NONLINEAR_MOBILITY')
            r['responses'].append(item);write(out/'report.json',r);print(json.dumps(item),flush=True)
        # Independent existing corrected-force path with literal elastic stress.
        direction=np.sin(np.arange(3*len(pa.x))).reshape(-1,3)*1e-9
        q=pa.extract_particles(np.arange(len(pa.x),dtype=np.uint32));q.u[:],q.v[:],q.w[:]=direction.T
        i,j,g,l,_,b=corrected_operators(q);strain=.5*(l+l.swapaxes(1,2));trace=np.trace(strain,axis1=1,axis2=2)
        stress=2*p['shear']*(strain-trace[:,None,None]*np.eye(3)/3)+bulk*trace[:,None,None]*np.eye(3)
        force=internal_force(q,stress,i,j,g,correction=b);difference=float(np.linalg.norm(op.internal_force(direction)-force))
        if difference>1e-12+2e-11*np.linalg.norm(force):raise ValueError('independent frozen elastic force mismatch')
        if before!=particle_hash(pa):raise ValueError('source particle mutated')
        for source,h in watched.items():
            if sha(source)!=h:raise ValueError('source input changed during probe')
        for rel,h in sources.items():
            if sha(ROOT/rel)!=h:raise ValueError('component source changed during probe')
        r.update(status='PASS_FROZEN_RESPONSE_COMPONENT_NOT_COUPLED_ADVANCE',
            independent_force_difference_N=difference,old_physical_end_time_s=state['t'],
            source_sha256=sources,input_sha256={path.relative_to(ROOT).as_posix():h for path,h in watched.items()})
        code=0
    except Exception as error:r.update(status='FAIL_FROZEN_RESPONSE_PROBE',error=repr(error));code=2
    finally:timer.cancel()
    r.update(wall_seconds=time.perf_counter()-start,peak_working_set_bytes=memory_bytes())
    write(out/'report.json',r);print(json.dumps({k:v for k,v in r.items() if k not in ('responses','source_sha256','input_sha256')}),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

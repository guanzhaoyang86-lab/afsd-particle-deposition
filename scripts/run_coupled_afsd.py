"""One bounded integrated rotating-bore/feed/substrate/plastic/thermal adapter.

No production permission, old-source mutation, restart migration or fabricated
trajectory. Sources/geometry are fixed before cost probe and registered trial.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import os
import sys
import threading
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.config_adapter import load_settings,UNITS
from afsd_backend.pysph.particles import create_geometry
from afsd_backend.pysph.axial_actuator import control
from afsd_backend.pysph.geometry import as_rs
from afsd_backend.pysph.mechanics import positions,velocities,pressure
from afsd_backend.pysph.material import equivalent
from afsd_backend.pysph.application import memory_bytes
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG=ROOT/'configs/benchmarks/coupled_afsd.json'
CLASSIFICATION='ENGINEERING_APPROXIMATE_COUPLED_AFSD_LOCAL_DEBUG'
SOURCES=['scripts/run_coupled_afsd.py','configs/benchmarks/coupled_afsd.json',
    'configs/debug/starved.json','configs/base.json','docs/COUPLED_AFSD_CONTRACT.md',
    'src/afsd_backend/pysph/coupled_solver.py','src/afsd_backend/pysph/coupled_boundary.py',
    'src/afsd_backend/pysph/coupled_thermal.py','src/afsd_backend/pysph/axial_actuator.py',
    'src/afsd_backend/pysph/solver.py','src/afsd_backend/pysph/material.py',
    'src/afsd_backend/pysph/mechanics.py','src/afsd_backend/pysph/kernel.py',
    'src/afsd_backend/pysph/contact.py','src/afsd_backend/pysph/contact_kick.py',
    'src/afsd_backend/pysph/contact_cone.py','src/afsd_backend/pysph/geometry.py',
    'src/afsd_backend/pysph/rigid_contact_kinematics.py','src/afsd_backend/pysph/thermal.py',
    'src/afsd_backend/pysph/particles.py','src/afsd_backend/pysph/integrator.py',
    'src/afsd_backend/pysph/output.py','src/afsd_backend/pysph/config_adapter.py',
    'src/afsd_backend/pysph/application.py','src/afsd_scaffold/config.py',
    'src/afsd_scaffold/vtk_legacy.py',
    'third_party/PySPH/pysph/sph/integrator.py',
    'third_party/PySPH/pysph/sph/integrator_step.py',
    'third_party/PySPH/pysph/base/nnps.py',
    'third_party/PySPH/pysph/base/nnps_base.pyx',
    'third_party/PySPH/pysph/base/linked_list_nnps.pyx']


def clean(value,rejected=False):
    if isinstance(value,np.ndarray):return clean(value.tolist(),rejected)
    if isinstance(value,np.generic):return clean(value.item(),rejected)
    if isinstance(value,dict):return {str(k):clean(v,rejected) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [clean(v,rejected) for v in value]
    if isinstance(value,float) and not np.isfinite(value):
        if rejected:return dict(nonfinite=repr(value),classification='REJECTED_DEBUG_VALUE')
        raise ValueError('nonfinite accepted output')
    return value


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fingerprint(data):return hashlib.sha256(json.dumps(clean(data),sort_keys=True,allow_nan=False).encode()).hexdigest()


def write_json(path,value,rejected=False):
    path=Path(path);pending=path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(clean(value,rejected),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    pending.replace(path)


def settings():
    cfg=json.loads(CONFIG.read_text(encoding='utf-8'))
    if cfg['classification']!=CLASSIFICATION or cfg['production_authorized'] is not False:
        raise ValueError('named nonproduction configuration required')
    original,p,consumed=load_settings(ROOT/cfg['parent'])
    expected=dict(UNITS,feed='m/s',release_time='s')
    for k,entry in cfg['parameters'].items():
        if k not in expected or entry['unit']!=expected[k] or not entry['source'] or not entry['status']:
            raise ValueError('invalid explicit parameter '+k)
        if not np.isfinite(entry['value']):raise ValueError('unresolved parameter '+k)
        p[k]=entry['value']
    p['normal_stiffness']=p['penalty_factor']*p['shear']*p['dx']**2/p['contact_reference_spacing']
    p['maximum_displacement_m']=cfg['maximum_displacement_m']
    if p['cfl']!=cfg['registered_cfl'] or p['layer']/p['dx']<3-1e-12:
        raise ValueError('registered CFL or resolved layer changed')
    resolved=dict(classification=CLASSIFICATION,parent_resolved_configuration=original,
        explicit_configuration=cfg,parameters=p,consumed_parent_parameters=consumed,
        derived_normal_stiffness=dict(value=p['normal_stiffness'],unit='N/m',
            source='penalty_factor*G*dx^2/contact_reference_spacing',status='INHERITED_APPROXIMATE'))
    return cfg,p,resolved


def initial_particles(p):
    pa=create_geometry(p)
    # All initial free velocities are zero. Only the declared axial cap moves.
    pa.u[:]=0.;pa.v[:]=0.;pa.w[:]=0.
    command=control(pa.z0,0.,0.,p['bar_length'],p['actuator_depth'],p['feed'],p['release_time'])
    pa.w[:]=command['prescribed_velocity'][:,2]
    pa.actuator[:]=command['eligible'].astype(int)  # metadata only
    pa.add_property('as_rs_tag',type='long',data=as_rs(pa.y0,p['omega'],p['traverse']))
    pa.set_output_arrays(list(pa.properties))
    return pa


def write_vtk(path,solver):
    pa,p=solver.pa,solver.p
    S=pa.S.reshape(-1,3,3);pr=pressure(pa.rho,p['rho0'],solver.bulk);sig=S-pr[:,None,None]*np.eye(3)
    scalars={k:getattr(pa,k).tolist() for k in ('particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','as_rs_tag')}
    scalars.update(mass_kg=pa.m.tolist(),temperature_K=pa.temperature.tolist(),density_kg_per_m3=pa.rho.tolist(),
        pressure_Pa=pr.tolist(),equivalent_stress_Pa=equivalent(S).tolist(),equivalent_plastic_strain=pa.ep.tolist(),
        thermal_energy_from_room_J_per_kg=(p['cp']*(pa.temperature-p['room'])).tolist())
    for a,b,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        scalars['stress_'+name+'_Pa']=sig[:,a,b].tolist()
    if solver.count:
        scalars['equivalent_deviatoric_strain_rate_per_s']=pa.strain_rate.tolist()
        scalars['equivalent_plastic_strain_rate_per_s']=pa.ep_rate.tolist()
    write_point_cloud(path,positions(pa).tolist(),scalars,
        dict(velocity_m_per_s=velocities(pa).tolist(),initial_position_m=np.column_stack((pa.x0,pa.y0,pa.z0)).tolist()),
        title=f'{CLASSIFICATION} time_s={solver.time:.17g}; SOLVED_THERMOMECHANICAL_NOT_PAPER_VALIDATED')


def budget_report(probe_manifest):
    probe_manifest=Path(probe_manifest).resolve()
    if not probe_manifest.is_relative_to(ROOT):raise ValueError('cost evidence must be within project')
    cfg,p,resolved=settings();m=json.loads(Path(probe_manifest).read_text(encoding='utf-8'))
    if m['status']!='COMPLETE_COST_PROBE' or m['configuration_sha256']!=fingerprint(resolved):
        raise ValueError('complete same-configuration cost probe required')
    if m['source_sha256']!={s:sha(ROOT/s) for s in SOURCES}:raise ValueError('cost probe source drift')
    for name,h in m['artifact_sha256'].items():
        if Path(name).name!=name or sha(Path(probe_manifest).parent/name)!=h:raise ValueError('cost artifact drift')
    dt=min(m['accepted_dt_s']);per_step=max(m['step_wall_seconds'])
    estimates=[]
    for end in cfg['registered_end_times_s']:
        steps=int(np.ceil(end/dt))
        prediction=cfg['forecast_multiplier']*per_step*steps+cfg['forecast_overhead_seconds']
        estimates.append(dict(end_time_s=end,estimated_steps=steps,forecast_wall_seconds=prediction,
            permitted=prediction<=cfg['forecast_launch_limit_seconds']))
    allowed=[x for x in estimates if x['permitted']]
    choice=allowed[-1] if allowed else None
    report=dict(status='READY_FOR_ONE_REGISTERED_TRIAL' if choice else 'NOT_RUN_COST_EXCEEDS_BUDGET',
        production_authorized=False,probe_manifest_path=probe_manifest.relative_to(ROOT).as_posix(),
        probe_manifest_sha256=sha(probe_manifest),source_sha256=m['source_sha256'],
        configuration_sha256=m['configuration_sha256'],estimates=estimates,choice=choice,
        method='1.5 * slowest measured complete step incl evidence I/O * estimated steps +30s; no plastic outcome fitting',
        peak_process_working_set_bytes=m['peak_process_working_set_bytes'],
        internal_limit_s=cfg['internal_wall_seconds'],external_limit_s=cfg['external_wall_seconds'])
    return report


def select_budget(probe_manifest,output):
    report=budget_report(probe_manifest)
    if Path(output).exists():raise FileExistsError('fresh budget file required')
    write_json(output,report);print(json.dumps(clean(report),allow_nan=False))
    return 0 if report['choice'] else 2


def run(output,probe=False,budget=None):
    from afsd_backend.pysph.coupled_solver import IntegratedSolver
    cfg,p,resolved=settings();source={s:sha(ROOT/s) for s in SOURCES};config_hash=fingerprint(resolved)
    directory=Path(output)
    if directory.exists():raise FileExistsError('fresh run directory required')
    if probe:target_steps=cfg['cost_probe_steps'];end=None;wall_limit=160.
    else:
        if budget is None:raise ValueError('prerun budget required')
        receipt=json.loads(Path(budget).read_text(encoding='utf-8'))
        expected=budget_report(ROOT/receipt['probe_manifest_path'])
        if receipt!=expected:raise ValueError('budget receipt differs from recomputed cost evidence')
        if (receipt['status']!='READY_FOR_ONE_REGISTERED_TRIAL' or receipt['source_sha256']!=source or
            receipt['configuration_sha256']!=config_hash):raise ValueError('budget missing or stale')
        end=receipt['choice']['end_time_s']
        if end not in cfg['registered_end_times_s']:raise ValueError('unregistered horizon')
        target_steps=10000;wall_limit=cfg['internal_wall_seconds']
    directory.mkdir(parents=True);started=time.perf_counter()
    write_json(directory/'resolved_config.json',resolved)
    solver=IntegratedSolver(initial_particles(p),p,source_hashes=source,configuration_hash=config_hash)
    history=[];frames=[];times=[];dts=[];error=None;saved=dict(step=0,time_s=0.)
    manifest=dict(classification=CLASSIFICATION,production_authorized=False,status='RUNNING',
        backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5',source_sha256=source,
        configuration_sha256=config_hash,mode='COST_PROBE' if probe else 'REGISTERED_INTEGRATED_TRIAL',
        requested_end_s=end,requested_steps=target_steps,maximum_wall_seconds=wall_limit,
        particle_count=len(solver.pa.x),initial_stable_dt_s=solver.stable_dt(),
        gate_definitions='docs/COUPLED_AFSD_CONTRACT.md; integration is not sustained deposition',
        field_status=dict(temperature='SOLVED',plasticity='SOLVED_FKS_APPROXIMATE',
            substrate='DEFORMABLE',feedstock='DEFORMABLE',infinite_inflow='NOT_IMPLEMENTED',
            total_free_energy='NOT_CERTIFIED_DIAGNOSTIC_ONLY'))
    write_json(directory/'manifest.json',manifest)
    def capture(vtk=False):
        record=clean(solver.record());history.append(record)
        solver.checkpoint(directory/'last_accepted.npz')
        if solver.count==2:solver.checkpoint(directory/'nonzero_accepted.npz')
        if np.any(solver.pa.ep>0) and not (directory/'plastic_accepted.npz').exists():
            solver.checkpoint(directory/'plastic_accepted.npz')
        with (directory/'history.jsonl').open('a',encoding='utf-8') as accepted_stream:
            accepted_stream.write(json.dumps(record,allow_nan=False)+'\n')
        if vtk:
            name=f'accepted_{solver.count:05d}.vtk';write_vtk(directory/name,solver)
            frames.append(dict(step=solver.count,time_s=solver.time,file=name,sha256=sha(directory/name)))
            write_json(directory/'series.json',dict(classification=CLASSIFICATION,frames=frames))
        saved.update(step=solver.count,time_s=solver.time)
    capture(vtk=True)
    def timed_out():
        write_json(directory/'HARD_TIMEOUT.json',dict(status='INCOMPLETE_WALL_LIMIT',
            last_saved_step=saved['step'],last_saved_time_s=saved['time_s'],restart_allowed=False,
            note='in-progress candidate is not an accepted endpoint; use saved checkpoint/history only'))
        os._exit(124)
    watchdog=threading.Timer(max(.001,wall_limit-(time.perf_counter()-started)),timed_out)
    watchdog.daemon=True;watchdog.start()
    try:
        with gzip.open(directory/'steps.jsonl.gz','wt',encoding='utf-8',compresslevel=1) as stream:
            for _ in range(target_steps):
                if end is not None and solver.time>=end:break
                if time.perf_counter()-started>=wall_limit:raise TimeoutError('wall budget at accepted boundary')
                tick=time.perf_counter();dt=solver.advance(max_dt=None if end is None else end-solver.time)
                dts.append(dt)
                stream.write(json.dumps(clean(solver.last_journal),allow_nan=False)+'\n');stream.flush()
                capture(vtk=solver.count%cfg['vtk_every_steps']==0 or (end is not None and solver.time>=end))
                times.append(time.perf_counter()-tick)
                if solver.count%10==0:print(json.dumps(dict(step=solver.count,time_s=solver.time,
                    ep_max=float(solver.pa.ep.max()),temperature_min_K=float(solver.pa.temperature.min()),
                    temperature_max_K=float(solver.pa.temperature.max()),wall_s=time.perf_counter()-started)),flush=True)
        status='COMPLETE_COST_PROBE' if probe else 'COMPLETE_REGISTERED_INTERVAL_NOT_SUSTAINED_DEPOSITION'
        if end is not None and solver.time<end:status='BLOCKED_STEP_CAP'
    except Exception as exc:
        status='BLOCKED';error=repr(exc)
        rejected=getattr(solver,'rejected_candidate',dict(classification='REJECTED',error=error,restart_allowed=False))
        write_json(directory/'REJECTED_candidate.json',rejected,rejected=True)
        print(error,flush=True)
    finally:
        watchdog.cancel()
        write_json(directory/'history.json',history)
        if not frames or frames[-1]['step']!=solver.count:
            name=f'accepted_{solver.count:05d}.vtk';write_vtk(directory/name,solver)
            frames.append(dict(step=solver.count,time_s=solver.time,file=name,sha256=sha(directory/name)))
            write_json(directory/'series.json',dict(classification=CLASSIFICATION,frames=frames))
        manifest.update(status=status,error=error,accepted_steps=solver.count,accepted_time_s=solver.time,
            wall_seconds=time.perf_counter()-started,step_wall_seconds=times,accepted_dt_s=dts,
            peak_process_working_set_bytes=memory_bytes(),particle_array_bytes=sum(a.get_npy_array().nbytes for a in solver.pa.properties.values()),
            final_diagnostics=solver.diagnostics(fail=False),final_record=clean(solver.record()))
        manifest['artifact_sha256']={f.name:sha(f) for f in sorted(directory.iterdir()) if f.is_file() and f.name!='manifest.json'}
        write_json(directory/'manifest.json',manifest)
    print(json.dumps(clean({k:v for k,v in manifest.items() if k not in ('final_record','final_diagnostics','source_sha256','artifact_sha256','step_wall_seconds','accepted_dt_s')}),allow_nan=False),flush=True)
    return 0 if status.startswith('COMPLETE') else 2


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--probe',action='store_true');group.add_argument('--budget',type=Path)
    group.add_argument('--select-budget',type=Path,metavar='PROBE_MANIFEST')
    a=parser.parse_args()
    return select_budget(a.select_budget,a.output) if a.select_budget else run(a.output,a.probe,a.budget)


if __name__=='__main__':raise SystemExit(main())

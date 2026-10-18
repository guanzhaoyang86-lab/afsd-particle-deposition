"""Fixed 12/14-Newton PURE branch diagnostic; NEVER a production restart.

Budget preregistered in docs/PROCESS_PLASTIC_FRICTION_Q7_BRANCH_DIAGNOSTIC.md.
Do not widen it after failure. No commit, advance, event search, or accepted
state output is permitted here. All solved candidate files remain REJECTED.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from probe_process_plastic_friction import PlasticFrictionBlock, ROOT
from probe_process_friction_cone import feasible_cone, vector_contact_work
from probe_process_plastic import CandidateFailure
from inspect_process_plastic_friction import inspect
from pysph.solver.utils import load


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_signature(a):
    a=np.asarray(a)
    return dict(dtype=a.dtype.str,shape=list(a.shape),strides=list(a.strides),sha256=hashlib.sha256(a.tobytes()).hexdigest())


def state_signature(b):
    """Every accepted particle/reference property, all block arrays/metadata."""
    result=dict(particles={},reference={},arrays={},metadata={},law=vars(b.law).copy())
    for key,pa in (('particles',b.pa),('reference',b.mesh.ref)):
        result[key]={name:array_signature(p.get_npy_array()) for name,p in sorted(pa.properties.items())}
    for name,value in sorted(vars(b).items()):
        if isinstance(value,np.ndarray):result['arrays'][name]=array_signature(value)
        elif isinstance(value,(str,int,float,bool,list,dict)) or value is None:
            result['metadata'][name]=copy.deepcopy(value)
    return result


def restore_diagnostic(path,r,h):
    """Read-only diagnostic reconstruction, bypassing NO production API guard.

    output.restore correctly rejects this checkpoint for production use; it is
    not called or weakened. Raw checkpoint loading is explicit only here for
    pure candidate evaluation. The accepted body cannot commit or advance.
    """
    b=PlasticFrictionBlock();saved=load(str(path/'checkpoint.npz'))
    pa=saved['arrays']['solid'];meta=saved['solver_data']
    if meta['restart_allowed'] is not False or meta['step']!=6 or meta['physical_time_s']!=.05:
        raise ValueError('diagnostic is restricted to this first accepted failure checkpoint')
    if set(pa.properties)!=set(b.pa.properties):raise ValueError('particle inventory differs from constructor')
    for name in pa.properties:
        actual=pa.properties[name].get_npy_array();base=b.pa.properties[name].get_npy_array()
        if actual.shape!=base.shape or actual.dtype!=base.dtype:raise ValueError('particle schema mismatch: '+name)
    if not np.array_equal(b.X,h['reference_position_m']) or not np.array_equal(b.V,h['reference_volume_m3']):
        raise ValueError('reference geometry/quadrature reconstruction mismatch')
    if not np.array_equal(pa.particle_id,h['particle_id']):raise ValueError('checkpoint ID order mismatch')
    b.pa=pa
    mapping=dict(displacement='displacement_m',stress='stress_Pa',normal_force='normal_force_N',
        tangent_force='tangent_force_global_yz_N',modes='interval_modes',contact_force='contact_force_N',
        support_force='support_force_N',friction='friction_dissipation_J',last_friction='interval_friction_J',
        last_slip='interval_relative_slip_global_yz_m',last_impulse='tangent_impulse_global_yz_N_s',
        last_relative_rate='relative_kinematic_rate_global_yz_m_s',rate='kinematic_interval_rate_m_s',
        last_law_residual='contact_law_residual_global_yz_N',last_plastic='interval_plastic_work_J',
        last_return='interval_return_diagnostic_J')
    for attr,key in mapping.items():
        # Actual state() allocates these via zeros_like(X), whose X is Fortran
        # layout. Preserve that layout too: reduction order affects exact replay.
        order='F' if attr in ('contact_force','support_force') else 'C'
        setattr(b,attr,np.array(h[key][-1],copy=True,order=order))
    row=r['records'][-1]
    scalars=dict(time='physical_time_s',count='step',wall_work='wall_work_J',
        contact_work='contact_material_work_J',support_work='support_work_J',internal_work='stress_work_J',
        last_wall_increment='interval_wall_work_J',last_contact_increment='interval_contact_material_work_J',
        last_residual='free_force_relative',last_iterations='corrections',last_trace='iteration_trace',
        last_force_defect='internal_force_relative',last_torque_defect='internal_torque_relative',
        accepted_floor_count='accepted_rate_floor_evaluations',last_event='event_kind')
    for attr,key in scalars.items():setattr(b,attr,copy.deepcopy(row[key]))
    b.constraint_work=b.contact_work+b.support_work;b.reaction=b.contact_force+b.support_force
    b.events=copy.deepcopy(r['events'])
    reconstructed=b.record()
    if reconstructed!=row:
        mismatch={key:dict(reconstructed=reconstructed.get(key),saved=row.get(key))
            for key in sorted(set(reconstructed)|set(row)) if reconstructed.get(key)!=row.get(key)}
        raise ValueError('reconstructed accepted record differs from saved record: '+json.dumps(mismatch))
    # No accidental path progress/commit may turn this debug load into restart.
    def forbidden(*args,**kwargs):raise RuntimeError('diagnostic cannot commit/advance accepted state')
    b.commit=forbidden;b.advance=forbidden
    return b


def assess_candidate(b,c,released):
    ft=b.limits['force_tolerance_N'];st=b.limits['slip_tolerance_m']
    normal=c['normal_force'];tau=c['tangent_force'];slip=c['relative_increment']
    capacity=b.mu*normal;length=np.linalg.norm(slip,axis=1);moving=length>st
    margin=capacity-np.linalg.norm(tau,axis=1)
    unknown=b.free.copy();unknown[b.contact_dofs]=False;unknown[b.tangent_dofs[~released].ravel()]=False
    root=c['residual'].ravel().copy();root[b.tangent_dofs[released].ravel()]=c['law_residual_N'][released].ravel()
    flow=np.linalg.norm(length[:,None]*tau+capacity[:,None]*slip,axis=1)
    allowed=ft*length+capacity*st
    oldmargin=b.mu*b.normal_force-np.linalg.norm(b.tangent_force,axis=1)
    oldoff=np.abs(oldmargin)>ft;opposite=np.sum(b.tangent_force*tau,axis=1)<0
    work=vector_contact_work(b.contact_force,c['contact_force'],c['displacement']-b.displacement,
        b.wall(c['time_s'])-b.wall(b.time),b.nodes,c['dt_s'])
    forced_violation=(~released)&(margin < -ft)
    result=dict(root_relative=float(np.linalg.norm(root[unknown])/b.force_scale),
        root_threshold_relative=b.ccfg['numerical']['force_relative'],
        root_pass=bool(np.linalg.norm(root[unknown])/b.force_scale<=b.ccfg['numerical']['force_relative']),
        equilibrium_relative=c['equilibrium_relative'],maximum_contact_root_N=float(np.linalg.norm(c['law_residual_N'],axis=1).max()),
        maximum_released_root_N=float(np.linalg.norm(c['law_residual_N'][released],axis=1).max()) if released.any() else 0.,
        all_contact_root_pass=bool(np.all(np.linalg.norm(c['law_residual_N'],axis=1)<=ft)),
        forced_stick_feasible=bool(not forced_violation.any()),
        released_flow_feasible=bool(np.all(flow[released&moving]<=allowed[released&moving])),
        branch_cone_feasible=bool(feasible_cone(c,released,b.mu,b.limits)),
        old_boundary_transition_risk=bool(np.any(moving&(oldoff|opposite))),
        event_split_required=bool(forced_violation.any() or np.any(moving&(oldoff|opposite))),
        minimum_raw_candidate_friction_J=float(work['friction_increment_J'].min()),
        raw_candidate_friction_J=work['friction_increment_J'].tolist(),
        candidate_plastic_work_J=float(c['plastic_increment_J'].sum()),
        contacts=[])
    for j,node in enumerate(b.nodes):
        result['contacts'].append(dict(particle_id=int(b.pa.particle_id[node]),feature_id=1,
            released=bool(released[j]),normal_N=float(normal[j]),tau_global_yz_N=tau[j].tolist(),
            cone_margin_N=float(margin[j]),slip_global_yz_m=slip[j].tolist(),moving=bool(moving[j]),
            root_norm_N=float(np.linalg.norm(c['law_residual_N'][j])),
            flow_residual_N_m=float(flow[j]),flow_allowance_N_m=float(allowed[j]),
            forced_stick_violated=bool(forced_violation[j]),old_boundary_off=bool(oldoff[j]),
            old_new_tangent_opposite=bool(opposite[j])))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True,type=Path);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter();path=args.run
    readback=inspect(path)
    if readback['status']!='PASS_PARTIAL_READBACK':raise ValueError('requires original partial BLOCKED run')
    r=json.loads((path/'report.json').read_text(encoding='utf-8'))
    rejected=json.loads((path/'REJECTED_manifest.json').read_text(encoding='utf-8'))
    if (r['accepted_steps']!=6 or r['accepted_time_s']!=.05 or rejected['candidate_time_s']!=.0625 or
            rejected['candidate_base_step']!=6 or rejected['error']!='cone branch equilibrium/root did not converge'):
        raise ValueError('not the preregistered first failure')
    with np.load(path/'history.npz',allow_pickle=False) as a:h={k:a[k].copy() for k in a.files}
    with np.load(path/'REJECTED_candidate.npz',allow_pickle=False) as a:
        rejected_arrays={k:a[k].copy() for k in a.files};released=a['released'].copy()
    if released.shape!=(9,) or released.dtype!=np.dtype(bool):raise ValueError('invalid saved released mask')
    evidence={name:sha(path/name) for name in ('report.json','history.npz','checkpoint.npz','REJECTED_manifest.json','REJECTED_candidate.npz')}
    attempts=[]
    for budget in (12,14):
        b=restore_diagnostic(path,r,h);before=state_signature(b);tick=time.perf_counter();error=None
        try:
            c=b.branch(.0625,released,maximum_iterations=budget)
            trace=c['iteration_trace'];status='BRANCH_ROOT_RETURNED_NOT_ACCEPTED'
        except CandidateFailure as exc:
            c=exc.candidate;trace=exc.iterations;error=str(exc);status='BRANCH_ROOT_BLOCKED'
        if c is None:raise ValueError('no candidate retained for diagnostic')
        elapsed=time.perf_counter()-tick;after=state_signature(b)
        if before!=after:raise ValueError('pure branch mutated accepted state')
        replay=None
        if budget==12:
            replay={key:bool(np.array_equal(c[key],value)) for key,value in rejected_arrays.items()}
            if not all(replay.values()) or trace!=rejected['iterations']:
                raise ValueError('12-Newton branch does not reproduce original failure exactly')
        np.savez_compressed(args.output/f'REJECTED_diagnostic_budget{budget}.npz',
            classification=np.array('REJECTED_DIAGNOSTIC'),restart_allowed=np.array(False),
            accepted_time_s=np.array(b.time),accepted_step=np.array(b.count),candidate_time_s=np.array(c['time_s']),
            **{k:v for k,v in c.items() if isinstance(v,np.ndarray)})
        attempts.append(dict(maximum_iterations=budget,status=status,error=error,wall_seconds=elapsed,
            classification='REJECTED_DIAGNOSTIC',restart_allowed=False,accepted_state_unchanged=True,
            accepted_time_s=b.time,accepted_step=b.count,exact_original_candidate_replay=replay,
            iterations=trace,assessment=assess_candidate(b,c,released)))
    if any(sha(path/name)!=digest for name,digest in evidence.items()):raise ValueError('input artifacts changed during diagnostic')
    report=dict(classification='PURE_ACCEPTED_STATE_BRANCH_DIAGNOSTIC_NOT_RESTART',production_authorized=False,
        status='DIAGNOSTIC_COMPLETE_NOT_PATH_COMPLETE',accepted_time_s=.05,accepted_step=6,
        target_candidate_time_s=.0625,readback_status=readback['status'],released_mask=released.tolist(),
        input_sha256=evidence,source_sha256=r['source_sha256'],diagnostic_source_sha256=sha(Path(__file__)),
        preregistration_sha256=sha(ROOT/'docs/PROCESS_PLASTIC_FRICTION_Q7_BRANCH_DIAGNOSTIC.md'),
        attempts=attempts,wall_seconds=time.perf_counter()-started,
        warning='Root convergence alone does not imply event feasibility, commit, path completion, or production restart.',
        artifact_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output/'diagnostic.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,allow_nan=False));return 0


if __name__=='__main__':raise SystemExit(main())

"""Read-only Q3 geometry, unilateral contact, work and independent J2 checks."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud
from pysph.solver.utils import load


def inspect(directory):
    path=Path(directory)
    r=json.loads((path/'report.json').read_text())
    cfg=json.loads((path/'resolved_config.json').read_text())
    for name,digest in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('source changed: '+name)
    for name,digest in r['artifact_sha256'].items():
        if hashlib.sha256((path/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('artifact changed: '+name)
    with np.load(path/'history.npz',allow_pickle=False) as f:
        h={k:f[k].copy() for k in f.files}
    if not all(np.isfinite(a).all() for a in h.values()):
        raise ValueError('nonfinite accepted state')
    t=h['physical_time_s']; rows=r['records']; ids=h['particle_id']
    if (len(rows)!=len(t) or len(t)!=r['accepted_steps']+1 or t[0]!=0 or
            t[-1]!=r['accepted_time_s'] or np.any(np.diff(t)<=0)):
        raise ValueError('ledger/history time inventory mismatch')
    if r['production_authorized'] or cfg['production_authorized']:
        raise ValueError('Q3 cannot authorize production')
    X=h['reference_position_m']; V=h['reference_volume_m3']; u=h['displacement_m']
    nodes=np.flatnonzero(X[:,0]==.01); fixed=X==0
    if (len(np.unique(ids))!=len(ids) or not np.equal(ids,np.floor(ids)).all() or
            not np.array_equal(h['contact_particle_id'],ids[nodes]) or
            not np.array_equal(h['contact_feature_id'],np.ones(len(nodes),dtype=int))):
        raise ValueError('persistent contact identity mismatch')
    if not np.array_equal(h['support_components'],fixed) or np.any(u[:,fixed]!=0):
        raise ValueError('permanent roller boundary changed')
    if not np.array_equal(h['position_m'],X[None]+u):
        raise ValueError('position/displacement mismatch')
    if abs(V.sum()/1e-6-1)>1e-12 or not np.allclose(h['mass_kg'],2700*V,rtol=1e-12,atol=0):
        raise ValueError('physical geometry/mass changed')
    ep=h['equivalent_plastic_strain']; P=h['plastic_work_J']; R=h['return_dissipation_diagnostic_J']
    if any(np.any(a<0) or np.any(np.diff(a,axis=0)<0) for a in (ep,P,R)):
        raise ValueError('raw negative plastic/return history')
    normal=h['normal_force_N']; active=h['active_contacts']; gap=h['gap_m']
    fc=h['contact_force_N']; fs=h['support_force_N']
    expected_fc=np.zeros_like(fc); expected_fc[:,nodes,0]=-normal
    if not np.array_equal(fc,expected_fc) or np.any(fs[:,~fixed]!=0):
        raise ValueError('contact direction/support field mismatch')
    if np.any(normal[~active]!=0) or np.any(normal< -9e-5) or np.any(gap< -1e-11) or np.any(np.abs(gap[active])>1e-11):
        raise ValueError('unilateral contact/complementarity violated')
    # Independent prescribed wall program, never particle displacement input.
    wall=np.array([-4e-5*(x/.1 if x<=.1 else 2-x/.1) for x in t])
    if not np.allclose(gap,wall[:,None]-u[:,nodes,0],rtol=0,atol=1e-18):
        raise ValueError('geometric gap mismatch')
    Wwall=Wcontact=Wsupport=0.; max_p_defect=max_ret_defect=max_interval_gapwork=0.
    for i,row in enumerate(rows):
        if row['step']!=i or row['physical_time_s']!=t[i]:
            raise ValueError('ledger/history frame mismatch')
        if row['active_contacts']!=int(active[i].sum()) or abs(row['wall_force_on_material_x_N']+normal[i].sum())>1e-9:
            raise ValueError('contact ledger mismatch')
        if i:
            du=u[i]-u[i-1]
            wc=float(np.sum(.5*(fc[i]+fc[i-1])*du))
            ww=-float(.5*(normal[i]+normal[i-1]).sum())*(wall[i]-wall[i-1])
            wr=float(np.sum(.5*(fs[i]+fs[i-1])*du))
            Wcontact+=wc; Wwall+=ww; Wsupport+=wr
            max_interval_gapwork=max(max_interval_gapwork,abs(ww-wc))
            if abs(ww-wc)>1e-6 or abs(wc-row['interval_contact_material_work_J'])>1e-12 or abs(ww-row['interval_wall_work_J'])>1e-12:
                raise ValueError('unsplit event or incorrect interval contact work')
            dep=ep[i]-ep[i-1]
            expectedP=90e6*V*dep if r['selected_law']=='j2' else np.zeros_like(V)
            max_p_defect=max(max_p_defect,float(np.max(np.abs(P[i]-P[i-1]-expectedP))))
            max_ret_defect=max(max_ret_defect,float(np.max(np.abs(R[i]-R[i-1]-1.5*35e9*V*dep**2))))
        for actual,name in ((Wwall,'wall_work_J'),(Wcontact,'contact_material_work_J'),(Wsupport,'support_work_J')):
            if abs(actual-row[name])>1e-12:
                raise ValueError('actual force/displacement integral mismatch')
        if abs(Wcontact+Wsupport-row['stress_work_J'])>1e-6 or abs(Wwall-Wcontact)>1e-6:
            raise ValueError('virtual work identity mismatch')
        if abs(P[i].sum()-row['plastic_work_J'])>1e-12 or abs(R[i].sum()-row['return_dissipation_diagnostic_J'])>1e-12:
            raise ValueError('particle/global energy ledger mismatch')
        sigma=h['stress_Pa'][i]; trace=np.trace(sigma,axis1=1,axis2=2)
        dev=sigma-trace[:,None,None]*np.eye(3)/3
        Es=float(np.sum(V*(np.sum(dev*dev,axis=(1,2))/(4*35e9)+trace**2/(18*(70e9/3)))))
        if abs(Es-row['elastic_diagnostic_J'])>1e-10:
            raise ValueError('independent elastic diagnostic mismatch')
    if max_p_defect>1e-10 or max_ret_defect>1e-10:
        raise ValueError('independent material energy formula mismatch')
    vtk=read_point_cloud(path/'final.vtk'); saved=load(str(path/'checkpoint.npz'))
    pa=saved['arrays']['solid']
    if saved['solver_data']['restart_allowed'] is not False:
        raise ValueError('diagnostic snapshot must not authorize production restart')
    if (not np.array_equal(vtk['scalars']['particle_id'],ids) or not np.array_equal(pa.particle_id,ids) or
            not np.array_equal(vtk['scalars']['tracer_type'],h['tracer_type']) or
            not np.array_equal(pa.tracer_type,h['tracer_type']) or not np.array_equal(pa.ep,ep[-1]) or
            not np.array_equal(vtk['points'],h['position_m'][-1])):
        raise ValueError('ID/label/accepted output mismatch')
    if any('temperature' in k or 'velocity' in k for k in vtk['scalars']|vtk['vectors']):
        raise ValueError('unsolved field published')
    complete=r['status']=='UNILATERAL_PATH_COMPLETE' and t[-1]==.2
    result=dict(status='PASS_COMPLETE_READBACK' if complete else 'PASS_PARTIAL_READBACK',
        completion_status=r['status'],production_authorized=False,accepted_time_s=float(t[-1]),
        accepted_steps=r['accepted_steps'],law=r['selected_law'],events=r['events'],
        independent_wall_work_J=Wwall,independent_contact_work_J=Wcontact,
        maximum_interval_gapwork_J=max_interval_gapwork,maximum_particle_plastic_defect_J=max_p_defect,
        maximum_particle_return_defect_J=max_ret_defect,
        final_mean_cap_displacement_m=float(np.mean(u[-1,nodes,0])),
        final_stress_max_abs_Pa=float(np.max(np.abs(h['stress_Pa'][-1]))),
        final_volume_mean_ep=float(np.sum(V*ep[-1])/V.sum()),final=rows[-1],
        note='single frictionless platen/roller-supported cube; NOT AFSD or unconstrained free-body release')
    if r['selected_law']=='j2':
        result['independent_continuum_reference']=dict(release_time_s=.13214285714285714,
            final_axial_displacement_m=-.00002714285714285714,
            ep_final=.002714285714285714,plastic_work_J=.24428571428571427,
            reference_is_implemented_trajectory=False)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    result=inspect(a.run)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2,allow_nan=False))
    raise SystemExit(0 if result['status']=='PASS_COMPLETE_READBACK' else 2)

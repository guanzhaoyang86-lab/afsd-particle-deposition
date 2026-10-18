"""Read-only Q4 per-wall work/contact identity and independent material references."""
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
    path=Path(directory); r=json.loads((path/'report.json').read_text())
    cfg=json.loads((path/'resolved_config.json').read_text())
    for name,digest in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest: raise ValueError('source changed: '+name)
    for name,digest in r['artifact_sha256'].items():
        if hashlib.sha256((path/name).read_bytes()).hexdigest()!=digest: raise ValueError('artifact changed: '+name)
    with np.load(path/'history.npz',allow_pickle=False) as f: h={k:f[k].copy() for k in f.files}
    if not all(np.isfinite(x).all() for x in h.values()): raise ValueError('nonfinite accepted history')
    t=h['physical_time_s']; rows=r['records']; X=h['reference_position_m']; V=h['reference_volume_m3']
    ids=h['particle_id']; u=h['displacement_m']; fixed=X==0
    if (len(t)!=len(rows) or len(t)!=r['accepted_steps']+1 or t[0]!=0 or t[-1]!=r['accepted_time_s'] or
            np.any(np.diff(t)<=0)): raise ValueError('ledger/history inventory mismatch')
    if r['production_authorized'] or cfg['production_authorized']: raise ValueError('Q4 cannot authorize production')
    nodes=np.concatenate([np.flatnonzero(X[:,a]==.01) for a in (0,1)])
    face=np.repeat(np.arange(2),9); axis=face; dofs=3*nodes+axis
    pairs=np.column_stack((h['contact_particle_id'],h['contact_feature_id']))
    if (len(np.unique(ids))!=27 or not np.equal(ids,np.floor(ids)).all() or len(nodes)!=18 or
            not np.array_equal(pairs,np.column_stack((ids[nodes],face+1))) or
            len(np.unique(pairs,axis=0))!=18 or len(np.unique(h['contact_particle_id']))!=15 or
            not np.array_equal(h['contact_axis'],axis) or not np.array_equal(h['contact_normal'],-np.eye(3)[axis])):
        raise ValueError('particle/feature contact identity mismatch')
    if (not np.array_equal(h['support_components'],fixed) or np.any(u[:,fixed]!=0) or
            not np.array_equal(h['position_m'],X[None]+u)): raise ValueError('geometry/support mismatch')
    if abs(V.sum()/1e-6-1)>1e-12 or not np.allclose(h['mass_kg'],2700*V,rtol=1e-12,atol=0):
        raise ValueError('physical mass/volume mismatch')
    # Independently reconstruct rigid wall program, not free particle motion.
    wall=[]
    for time in t:
        dx=-4e-5*(time/.1 if time<=.1 else 2-time/.1)
        if r['schedule']=='symmetric' or time<=.1: dy=dx
        elif time<=.15: dy=-4e-5
        else: dy=-4e-5*(.2-time)/.05
        wall.append([dx,dy])
    wall=np.asarray(wall)
    normal=h['normal_force_N']; active=h['active_contacts']; gap=h['gap_m']
    ff=h['contact_force_by_wall_N']; support=h['support_force_N']
    expected=np.zeros_like(ff)
    for f in range(2):
        for row in np.flatnonzero(face==f): expected[:,f,nodes[row],axis[row]]-=normal[:,row]
    if not np.array_equal(ff,expected) or np.any(support[:,~fixed]!=0): raise ValueError('per-wall force assembly mismatch')
    if (np.any(normal[~active]!=0) or np.any(normal< -9e-5) or np.any(gap< -1e-11) or
            np.any(np.abs(gap[active])>1e-11)): raise ValueError('unilateral contact violated')
    if not np.allclose(gap,wall[:,face]-u.reshape(len(t),-1)[:,dofs],rtol=0,atol=1e-18):
        raise ValueError('geometric contact gap mismatch')
    ep=h['equivalent_plastic_strain']; P=h['plastic_work_J']; R=h['return_dissipation_diagnostic_J']
    if any(np.any(a<0) or np.any(np.diff(a,axis=0)<0) for a in (ep,P,R)): raise ValueError('raw negative material history')
    Ww=np.zeros(2); Wc=np.zeros(2); Wr=0.; max_gapwork=np.zeros(2); p_defect=ret_defect=0.
    for i,row in enumerate(rows):
        if row['step']!=i or row['physical_time_s']!=t[i]: raise ValueError('ledger/history frame mismatch')
        if row['active_contacts_by_wall']!=[int(active[i,face==f].sum()) for f in range(2)]:
            raise ValueError('per-wall active inventory mismatch')
        if not np.allclose(row['wall_force_on_material_N'],ff[i].sum(axis=1),rtol=0,atol=1e-9):
            raise ValueError('per-wall reaction mismatch')
        if i:
            du=u[i]-u[i-1]; mid=.5*(normal[i]+normal[i-1])
            wc=np.sum(.5*(ff[i]+ff[i-1])*du[None],axis=(1,2))
            ww=-np.array([mid[face==f].sum() for f in range(2)])*(wall[i]-wall[i-1])
            wr=float(np.sum(.5*(support[i]+support[i-1])*du))
            Ww+=ww; Wc+=wc; Wr+=wr; max_gapwork=np.maximum(max_gapwork,np.abs(ww-wc))
            if (np.any(np.abs(ww-wc)>1e-6) or
                    not np.allclose(row['interval_wall_work_J'],ww,rtol=0,atol=1e-12) or
                    not np.allclose(row['interval_contact_material_work_J'],wc,rtol=0,atol=1e-12)):
                raise ValueError('per-wall interval work mismatch or unsplit event')
            dep=ep[i]-ep[i-1]
            expected_p=90e6*V*dep if r['selected_law']=='j2' else np.zeros_like(V)
            p_defect=max(p_defect,float(np.max(np.abs(P[i]-P[i-1]-expected_p))))
            ret_defect=max(ret_defect,float(np.max(np.abs(R[i]-R[i-1]-1.5*35e9*V*dep**2))))
        if (not np.allclose(row['wall_work_J'],Ww,rtol=0,atol=1e-12) or
                not np.allclose(row['contact_material_work_J'],Wc,rtol=0,atol=1e-12) or
                abs(row['support_work_J']-Wr)>1e-12 or abs(row['stress_work_J']-Wc.sum()-Wr)>1e-6 or
                np.any(np.abs(Ww-Wc)>1e-6)): raise ValueError('per-wall cumulative/virtual work mismatch')
        if abs(P[i].sum()-row['plastic_work_J'])>1e-12 or abs(R[i].sum()-row['return_dissipation_diagnostic_J'])>1e-12:
            raise ValueError('material/global energy mismatch')
        sigma=h['stress_Pa'][i]; trace=np.trace(sigma,axis1=1,axis2=2)
        dev=sigma-trace[:,None,None]*np.eye(3)/3
        Es=float(np.sum(V*(np.sum(dev*dev,axis=(1,2))/(4*35e9)+trace**2/(18*(70e9/3)))))
        if abs(Es-row['elastic_diagnostic_J'])>1e-10: raise ValueError('elastic diagnostic mismatch')
    if p_defect>1e-10 or ret_defect>1e-10: raise ValueError('independent material energy formula mismatch')
    event_steps=[event['step'] for event in r['events']]
    ledger_event_steps=[i for i,row in enumerate(rows) if row['event_kind']=='RELEASE_LIMIT']
    if event_steps!=ledger_event_steps:
        raise ValueError('release event/ledger inventory mismatch')
    for event in r['events']:
        i=event['step']; bracket=event['bracket_s']
        if (event['kind']!='RELEASE_LIMIT' or not 0<i<len(t) or event['physical_time_s']!=t[i] or rows[i]['event_kind']!='RELEASE_LIMIT' or
                len(bracket)!=2 or bracket[0]!=t[i] or not 0<=bracket[1]-bracket[0]<=1e-12):
            raise ValueError('accepted release event mismatch')
        near=active[i] & (normal[i]<=9e-5)
        if (event['near_zero_contact_rows']!=np.flatnonzero(near).tolist() or
                event['near_zero_feature_ids']!=(np.unique(face[near])+1).tolist()):
            raise ValueError('release feature identity mismatch')
    vtk=read_point_cloud(path/'final.vtk'); saved=load(str(path/'checkpoint.npz')); pa=saved['arrays']['solid']
    if saved['solver_data']['restart_allowed'] is not False: raise ValueError('diagnostic snapshot is not restartable')
    if (not np.array_equal(vtk['scalars']['particle_id'],ids) or not np.array_equal(pa.particle_id,ids) or
            not np.array_equal(vtk['scalars']['tracer_type'],h['tracer_type']) or not np.array_equal(pa.tracer_type,h['tracer_type']) or
            not np.array_equal(pa.ep,ep[-1]) or not np.array_equal(vtk['points'],h['position_m'][-1])):
        raise ValueError('accepted ID/label/output mismatch')
    if any('temperature' in key or 'velocity' in key for key in vtk['scalars']|vtk['vectors']): raise ValueError('unsolved field published')
    complete=r['status']=='CORNER_PATH_COMPLETE' and t[-1]==.2
    if complete and any(not np.any(np.abs(t-x)<=1e-14) for x in np.linspace(0.,.2,2*r['increments_per_phase']+1)):
        raise ValueError('completed path missing requested time grid')
    result=dict(status='PASS_COMPLETE_READBACK' if complete else 'PASS_PARTIAL_READBACK',completion_status=r['status'],
        production_authorized=False,law=r['selected_law'],schedule=r['schedule'],accepted_steps=r['accepted_steps'],
        accepted_time_s=float(t[-1]),events=r['events'],independent_wall_work_J=Ww.tolist(),
        independent_material_contact_work_J=Wc.tolist(),maximum_interval_gapwork_by_wall_J=max_gapwork.tolist(),
        maximum_particle_plastic_defect_J=p_defect,maximum_particle_return_defect_J=ret_defect,
        final_surface_displacement_m=[float(u[-1,X[:,a]==.01,a].mean()) for a in range(3)],
        final_stress_max_abs_Pa=float(np.max(np.abs(h['stress_Pa'][-1]))),final=rows[-1],
        caution='Two orthogonal fixed-normal benchmark walls, NOT rotating tool, thermal AFSD, or finite-strain validation.')
    if r['selected_law']=='j2':
        result['independent_continuum_reference']=dict(ep_final=.005428571428571428,
            plastic_work_J=.48857142857142855,final_surface_displacement_m=[-27.14285714285714e-6,-27.14285714285714e-6,54.28571428571428e-6],
            release_time_by_wall_s=[.13214285714285715,.13214285714285715 if r['schedule']=='symmetric' else .16607142857142856])
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--run',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    result=inspect(a.run); a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2,allow_nan=False)); raise SystemExit(0 if result['status']=='PASS_COMPLETE_READBACK' else 2)

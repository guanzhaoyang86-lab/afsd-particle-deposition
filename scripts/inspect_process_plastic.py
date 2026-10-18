"""Independent saved-data work/ID/material checks for Q2, not AFSD validation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud
from pysph.solver.utils import load


def inspect(directory):
    path=Path(directory)
    report=json.loads((path/'report.json').read_text())
    config=json.loads((path/'resolved_config.json').read_text())
    for name,digest in report['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('source changed: '+name)
    for name,digest in report['artifact_sha256'].items():
        if hashlib.sha256((path/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('artifact changed: '+name)
    with np.load(path/'history.npz',allow_pickle=False) as f:
        h={key:f[key].copy() for key in f.files}
    if not all(np.isfinite(value).all() for value in h.values()):
        raise ValueError('nonfinite accepted history')
    t=h['physical_time_s']; V=h['reference_volume_m3']; X=h['reference_position_m']
    if t[0]!=0 or t[-1]!=report['accepted_time_s'] or len(t)!=report['accepted_steps']+1 or np.any(np.diff(t)<=0):
        raise ValueError('accepted endpoint/time mismatch')
    complete=(report['status']=='DISCRETE_PLASTIC_PATH_COMPLETE' and t[-1]==.2)
    if report['production_authorized']:
        raise ValueError('Q2 cannot authorize production')
    ids=h['particle_id']
    if len(np.unique(ids))!=len(ids) or not np.equal(ids,np.floor(ids)).all():
        raise ValueError('bad persistent ID')
    if abs(V.sum()/1e-6-1)>1e-12 or not np.allclose(h['mass_kg'],2700*V,rtol=1e-12,atol=0):
        raise ValueError('physical volume/mass changed')
    if not np.array_equal(h['position_m'],X[None]+h['displacement_m']):
        raise ValueError('position/displacement mismatch')
    fixed=X==0; cap=X[:,0]==.01; fixed[cap,0]=True
    if not np.array_equal(fixed,h['constrained_components']) or not np.array_equal(cap,h['cap_particles']):
        raise ValueError('boundary inventory changed')
    ep=h['equivalent_plastic_strain']; P=h['plastic_work_J']; R=h['return_dissipation_diagnostic_J']
    if np.any(ep<0) or np.any(np.diff(ep,axis=0)<0) or np.any(P<0) or np.any(np.diff(P,axis=0)<0):
        raise ValueError('raw negative plastic history')
    if np.any(R<0) or np.any(np.diff(R,axis=0)<0):
        raise ValueError('raw negative return history')
    work=0.; maximum_P_defect=0.; maximum_return_defect=0.
    records=report['records']
    if len(records)!=len(t):
        raise ValueError('ledger/history frame count mismatch')
    for i,row in enumerate(records):
        if row['step']!=i or row['physical_time_s']!=t[i]:
            raise ValueError('ledger/history step or time mismatch')
        prescribed=np.zeros_like(X)
        fraction=t[i]/.1 if t[i]<=.1 else 2-t[i]/.1
        prescribed[cap,0]=-.004*.01*fraction
        if not np.array_equal(h['displacement_m'][i][fixed],prescribed[fixed]):
            raise ValueError('actuator/support moved off physical program')
        if i:
            du=h['displacement_m'][i]-h['displacement_m'][i-1]
            work+=float(np.sum(.5*(h['reaction_N'][i]+h['reaction_N'][i-1])*du))
            if abs(work-row['actuator_work_J'])>1e-12:
                raise ValueError('actuator work not actual reaction-displacement integral')
            dep=ep[i]-ep[i-1]
            if report['selected_law']=='j2':
                yield_value=90e6
            else:
                f=config['fks_parameters']; temperature=config['parameters']['fks_temperature']['value']
                ts=(temperature-f['room'])/(f['melt']-f['room'])
                theta=-math.expm1(-math.log1p(math.exp(-f['c1']*ts))/f['c2'])
                rate=h['deviatoric_strain_rate_per_s'][i]
                lam=1+f['b1']*ts**f['b2']*np.log(np.maximum(rate,f['rate_floor'])/f['reference_rate'])
                yield_value=(f['a1']+f['a2']*np.arctan(f['a3']*ep[i]))*lam*theta
            expected_P=V*yield_value*dep
            expected_return=V*1.5*35e9*dep**2
            maximum_P_defect=max(maximum_P_defect,float(np.max(np.abs((P[i]-P[i-1])-expected_P))))
            maximum_return_defect=max(maximum_return_defect,float(np.max(np.abs((R[i]-R[i-1])-expected_return))))
        if abs(work-row['stress_work_J'])>1e-6:
            raise ValueError('actuator/internal virtual work mismatch')
        if abs(P[i].sum()-row['plastic_work_J'])>1e-12 or abs(R[i].sum()-row['return_dissipation_diagnostic_J'])>1e-12:
            raise ValueError('per-particle work/ledger mismatch')
    if maximum_P_defect>1e-10 or maximum_return_defect>1e-10:
        raise ValueError('plastic work/return differs from independent material reference')
    vtk=read_point_cloud(path/'final.vtk'); saved=load(str(path/'checkpoint.npz'))
    pa=saved['arrays']['solid']; state=saved['solver_data']
    if state['restart_allowed'] is not False:
        raise ValueError('diagnostic snapshot incorrectly restartable')
    for actual in (vtk['scalars']['particle_id'],pa.particle_id):
        if not np.array_equal(actual,ids): raise ValueError('persistent ID lost at output')
    if not np.array_equal(vtk['points'],h['position_m'][-1]) or not np.array_equal(pa.ep,ep[-1]):
        raise ValueError('accepted end output mismatch')
    if any('temperature' in k or 'velocity' in k for k in vtk['scalars']|vtk['vectors']):
        raise ValueError('unsolved thermal/dynamic field published')
    final=records[-1]
    result=dict(status='PASS_COMPLETE_READBACK' if complete else 'PASS_PARTIAL_READBACK',
        completion_status=report['status'],accepted_end_s=float(t[-1]),accepted_steps=report['accepted_steps'],
        selected_law=report['selected_law'],n=report['n'],production_authorized=False,
        independent_actuator_work_J=work,maximum_particle_plastic_work_defect_J=maximum_P_defect,
        maximum_particle_return_defect_J=maximum_return_defect,final=final,
        final_volume_mean_sigma_xx_Pa=float(np.sum(V*h['stress_Pa'][-1,:,0,0])/V.sum()),
        final_lateral_displacement_max_m=float(np.linalg.norm(h['displacement_m'][-1,:,1:],axis=1).max()),
        caution='Axial actuator return is NOT free unloaded residual shape; plastic work is NOT solved heat.')
    if report['selected_law']=='j2' and complete:
        result['independent_continuum_reference']=dict(ep_final=.008-3*90e6/70e9,
            sigma_xx_final_Pa=90e6,plastic_work_final_J=90e6*(.008-3*90e6/70e9)*1e-6,
            note='uniform continuum reference, not an admission threshold for coarse SPH')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    result=inspect(args.run)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2,allow_nan=False))
    raise SystemExit(0 if result['status']=='PASS_COMPLETE_READBACK' else 2)

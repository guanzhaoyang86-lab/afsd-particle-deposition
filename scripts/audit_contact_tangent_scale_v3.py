"""Saved true-parent contact-only finite-difference diagnostic; no process solve."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import threading
import time
import traceback
import numpy as np
from batched_contact import BatchedChainRuleContact

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/'outputs/expanded_hot_budget_prefix_20261005/window_001/fine/step_012/accepted.npz'
PARENT_SHA='12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b'
FOLDER=ROOT/'outputs/contact_onset_time_study_20261005/v1/eight'
EPSILONS=(1e-9,1e-10,1e-11,1e-12)
FLAGS=('active','friction_active','sliding','normal_derivative_limit_branch')


def require(value,message):
    if not bool(value):raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    return value


def write(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(clean(value),f,indent=2,allow_nan=False);f.write('\n')


def pose(p,t):
    c,s=np.cos(p['omega']*t),np.sin(p['omega']*t)
    return dict(origin=np.array([p['traverse']*t,0.,p['layer']]),rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]))


def audit(output):
    start=time.perf_counter();deadline=start+25.;evidence={};evaluations=0;code=2
    report=dict(status='FAIL_CONTACT_ONLY_TANGENT_SCALE_AUDIT',accepted=False,restart_allowed=False,
        physical_solver_evaluations=0,nnps_evaluations=0,material_evaluations=0,thermal_evaluations=0,
        process_time_advanced_s=0.,contact_path_check=False,classification='SAVED_CONTACT_COMPONENT_SCALE_DIAGNOSTIC_NOT_PROCESS')
    def bind(path,digest=None):
        path=Path(path).resolve();require(path.is_relative_to(ROOT),'input outside project')
        value=sha(path);require(digest is None or value==digest,'source/input drift '+rel(path));evidence[rel(path)]=value
        return path
    def timeout():
        try:write(output/'HARD_TIMEOUT.json',dict(status='TIMEOUT_CONTACT_ONLY_DIAGNOSTIC',limit_s=27.,accepted=False))
        finally:os._exit(124)
    timer=threading.Timer(27.,timeout);timer.daemon=True;timer.start()
    try:
        bind(Path(__file__).resolve());bind(PARENT,PARENT_SHA)
        rows=read(bind(FOLDER/'candidate_history.json'));manifest=read(bind(FOLDER/'manifest.json'))
        require(manifest['status']=='COMPLETED_CONTACT_ONSET_RESOLUTION_DIAGNOSTIC_NOT_ACCEPTED' and
            manifest['physical_trajectory_completed'] is True and len(rows)==8,'completed saved eight-step study required')
        for name,digest in manifest['source_sha256'].items():bind(ROOT/name,digest)
        row=rows[0];require(row['parent_sha256']==PARENT_SHA and row['dt_s']==9.765625e-7,'step001 is not the fixed true-parent candidate')
        for name,digest in row['artifacts'].items():bind(FOLDER/name,digest)
        with np.load(PARENT,allow_pickle=True) as z:
            state=z['solver_data'].item();pa=z['particles'].item()['solid']['arrays']
        require(state['step']==14 and state['accepted'] is True and state['restart_allowed'] is True,'parent is not actual accepted step14')
        with np.load(FOLDER/'step_001/FULL_UNCOMMITTED_candidate.npz',allow_pickle=False) as z:
            saved={k:z[k].copy() for k in z.files}
        detail=read(FOLDER/'step_001/FULL_UNCOMMITTED_candidate.json');p=state['parameters'];dt=row['dt_s'];ids=pa['particle_id']
        require(np.array_equal(ids,saved['particle_id']) and not bool(saved['accepted']) and not bool(saved['restart_allowed']),
            'saved diagnostic identity/classification differs')
        x0=np.column_stack([pa[k] for k in ('x','y','z')]);x1=saved['position_m'];history=state['contact_history']
        require(np.array_equal(x0+saved['displacement_m'],x1),'saved endpoint differs from same-parent displacement')
        before=[value.copy() for value in (x0,x1,history)]
        contact=BatchedChainRuleContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
            tool_height=p['tool_top']-p['layer'],particle_radius=p['dx']/2,normal_stiffness=p['normal_stiffness'],
            tangential_stiffness=p['normal_stiffness'],friction=p['friction'],base_z=-p['substrate_depth'])
        args=(pose(p,state['time_s']),pose(p,state['time_s']+dt),dt,history)
        def evaluate(end):
            nonlocal evaluations
            require(time.perf_counter()<deadline,'registered contact-only deadline exhausted')
            evaluations+=1
            return contact.evaluate(x0,end,*args,check_path=False)
        base=evaluate(x1);differences={}
        for key,value in base.items():
            if isinstance(value,np.ndarray):
                other=np.asarray(detail['contact'][key]);require(value.shape==other.shape,'saved contact shape differs '+key)
                require(np.allclose(value,other,rtol=1e-12,atol=1e-15,equal_nan=False),'contact-only baseline differs from actual saved contact '+key)
                finite=np.isfinite(value)&np.isfinite(other)
                differences[key]=float(np.max(np.abs(value[finite].astype(float)-other[finite].astype(float)))) if finite.any() else 0.
        require(np.array_equal(base['force_N'],saved['contact_force_N']),'baseline force not identical to actual full candidate')
        forces=np.zeros((4,3,2,len(ids),3));normal=np.zeros_like(forces);tangent=np.zeros_like(forces)
        flags=np.zeros((4,3,2,len(ids),8,4),bool);blocks=np.zeros((4,len(ids),3,3));scales=[]
        affected=set();features=list(base['features'])
        def changed_pairs(mask):
            return [dict(particle_id=int(ids[i]),feature=features[j]) for i,j in np.argwhere(mask)]
        for hindex,h in enumerate(EPSILONS):
            axisrows=[]
            for axis in range(3):
                endplus=x1.copy();endminus=x1.copy();endplus[:,axis]+=h;endminus[:,axis]-=h
                plus,minus=evaluate(endplus),evaluate(endminus)
                for sign,result in enumerate((plus,minus)):
                    forces[hindex,axis,sign]=result['force_N']
                    normal[hindex,axis,sign]=result['normal_force_by_feature_N'].sum(axis=1)
                    tangent[hindex,axis,sign]=result['tangential_force_by_feature_N'].sum(axis=1)
                    for f,key in enumerate(FLAGS):flags[hindex,axis,sign,:,:,f]=result[key]
                blocks[hindex,:,:,axis]=(plus['force_N']-minus['force_N'])/(2*h)
                oneplus=(plus['force_N']-base['force_N'])/h;oneminus=(base['force_N']-minus['force_N'])/h
                branch={}
                for key in FLAGS:
                    differing=plus[key]!=minus[key];touch=(plus[key]!=base[key])|(minus[key]!=base[key])
                    affected.update(int(v) for v in ids[np.any(touch,axis=1)])
                    branch[key]=dict(plus_count=int(plus[key].sum()),minus_count=int(minus[key].sum()),
                        base_count=int(base[key].sum()),plus_minus_changed_pairs=changed_pairs(differing),
                        either_side_changed_from_base=changed_pairs(touch))
                axisrows.append(dict(axis=axis,axis_name='xyz'[axis],branch=branch,
                    central_column_max_N_m=float(np.max(np.abs(blocks[hindex,:,:,axis]))),
                    one_sided_derivative_disagreement_max_N_m=float(np.max(np.abs(oneplus-oneminus))),
                    maximum_one_sided_disagreement_particle_id=int(ids[np.argmax(np.max(np.abs(oneplus-oneminus),axis=1))])))
            scales.append(dict(epsilon_m=h,axes=axisrows,block_frobenius_norm_N_m=float(np.linalg.norm(blocks[hindex])),
                maximum_block_frobenius_N_m=float(np.linalg.norm(blocks[hindex],axis=(1,2)).max())))
        deltas=[]
        for index in range(3):
            delta=blocks[index]-blocks[index+1];byid=np.linalg.norm(delta,axis=(1,2))
            denominator=np.maximum(np.linalg.norm(blocks[index+1],axis=(1,2)),1.)
            k=int(np.argmax(byid));relative=byid/denominator
            deltas.append(dict(coarser_epsilon_m=EPSILONS[index],finer_epsilon_m=EPSILONS[index+1],
                total_relative_block_difference=float(np.linalg.norm(delta)/max(np.linalg.norm(blocks[index+1]),1.)),
                maximum_particle_block_difference_N_m=float(byid[k]),maximum_particle_id=int(ids[k]),
                maximum_particle_relative_difference=float(relative.max()),
                maximum_relative_particle_id=int(ids[np.argmax(relative)]),relative_denominator_floor_N_m=1.))
        # Exact participant IDs, penetration/clearance and matrices allow independent reading.
        participant=[]
        close=np.argwhere(np.abs(base['gaps_m']-contact.particle_radius)<=2e-9)
        affected.update(int(ids[i]) for i,j in close)
        affected.add(int(ids[np.argmax(np.linalg.norm(blocks[0]-blocks[-1],axis=(1,2)))]))
        for pid in sorted(affected):
            i=int(np.flatnonzero(ids==pid)[0])
            participant.append(dict(particle_id=pid,component_id=int(pa['component_id'][i]),free=state['free'][i],
                saved_displacement_m=saved['displacement_m'][i],parent_history_m=history[i],
                end_clearance_to_penalty_shell_m=base['gaps_m'][i]-contact.particle_radius,
                baseline_active=base['active'][i],baseline_sliding=base['sliding'][i],
                contact_force_N=base['force_N'][i],blocks_by_epsilon_N_m=blocks[:,i]))
        np.savez_compressed(output/'contact_scale_arrays_NOT_RESTARTABLE.npz',particle_id=ids,
            x0=x0,x1=x1,history0=history,epsilon_m=np.asarray(EPSILONS),forces_N=forces,
            normal_forces_N=normal,tangential_forces_N=tangent,flags=flags,blocks_N_m=blocks,
            baseline_force_N=base['force_N'],baseline_flags=np.stack([base[k] for k in FLAGS],axis=-1),
            accepted=False,restart_allowed=False)
        require(all(np.array_equal(a,b) for a,b in zip(before,(x0,x1,history))),'contact-only calls mutated input')
        report.update(status='PASS_CONTACT_ONLY_SCALE_DIAGNOSTIC_NOT_PROCESS_QUALIFICATION',
            parent_sha256=PARENT_SHA,parent_time_s=state['time_s'],candidate_time_s=row['time_s'],dt_s=dt,
            saved_candidate_sha256=row['sha256'],particle_count=len(ids),baseline_contact_field_differences=differences,
            baseline_exact_saved_force_match=True,epsilon_m=EPSILONS,
            displacement_max_m=float(np.linalg.norm(saved['displacement_m'],axis=1).max()),
            normal_stiffness_N_m=p['normal_stiffness'],scales=scales,adjacent_epsilon_differences=deltas,
            affected_or_near_shell_particles=participant,
            limitations=['One actual endpoint and its true parent only; not all failed line-search iterates.',
                'Contact-local block derivative only; no NNPS, material/thermal tangent, full Jacobian, root or process acceptance.',
                'Branch crossing of a finite difference can alter a search approximation without proving the unique cause of nonlinear failure.',
                'Smaller h convergence or lack of convergence is observational, not a new approved epsilon or time-accuracy solution.'])
        code=0
    except Exception as error:report.update(error=repr(error),traceback=traceback.format_exc())
    finally:
        timer.cancel()
        for name,digest in evidence.items():
            if sha(ROOT/name)!=digest:report.update(status='FAIL_CONTACT_ONLY_TANGENT_SCALE_AUDIT',error='bound source/input changed');code=2
        report.update(contact_evaluations=evaluations,wall_seconds=time.perf_counter()-start,exit_code=code,input_sha256=evidence,
            output_sha256={rel(p):sha(p) for p in output.iterdir() if p.is_file()})
        write(output/'summary.json',report)
    print(json.dumps({k:report[k] for k in ('status','exit_code','wall_seconds','contact_evaluations')}));return code


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True);args=parser.parse_args()
    out=(ROOT/args.output).resolve();require(out.is_relative_to(ROOT/'outputs') and out!=ROOT/'outputs','new output directory required')
    out.mkdir(parents=True,exist_ok=False);return audit(out)


if __name__=='__main__':raise SystemExit(main())

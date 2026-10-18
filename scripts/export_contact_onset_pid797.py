"""Export actual saved PID797 states only; no physical evaluation or interpolation."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:sys.path.insert(0,str(ROOT/'src'))
from pysph.solver.utils import load

BASE=ROOT/'outputs/contact_onset_time_study_20261005/v1'
OUT=BASE/'visualization_new'

def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(value,message):
    if not bool(value):raise ValueError(message)
def rel(path):return Path(path).relative_to(ROOT).as_posix()


def main():
    require(not OUT.exists(),'visualization output already exists');OUT.mkdir()
    summary_path=BASE/'readback_eight/summary.json';summary=read(summary_path)
    require(summary['status']=='PASS_CONTACT_ONSET_TIME_STUDY_READBACK_NOT_PHYSICAL_VALIDATION' and
        summary['actual_accepted_global_step']==14 and abs(summary['actual_accepted_time_s']-62.5e-6)<1e-16,
        'registered readback/accepted parent differs')
    bindings={rel(summary_path):sha(summary_path),rel(Path(__file__)):sha(Path(__file__))}
    rows=[]
    def extract(path,digest,branch,dt_s):
        require(sha(path)==digest,'actual snapshot changed: '+rel(path));bindings[rel(path)]=digest
        saved=load(str(path));pa=saved['arrays']['solid'];state=saved['solver_data'];index=np.flatnonzero(pa.particle_id==797)
        require(len(index)==1,'persistent ID797 absent/duplicated');k=int(index[0]);S=pa.S.reshape(-1,3,3)[k]
        if branch!='shared_parent':require(state.get('accepted') is False and state.get('restart_allowed') is False,'diagnostic checkpoint restart flag differs')
        row=dict(branch=branch,time_s=float(state['time_s']),dt_s=dt_s,particle_id=797,
            q_Pa=float(np.sqrt(1.5*np.sum(S*S))),stress_frobenius_Pa=float(np.linalg.norm(S)),
            ep=float(pa.ep[k]),strain_rate_per_s=float(pa.strain_rate[k]),ep_rate_per_s=float(pa.ep_rate[k]),
            temperature_K=float(pa.temperature[k]),checkpoint=rel(path),checkpoint_sha256=digest)
        rows.append(row);return state,row
    parent=ROOT/summary['actual_accepted_parent'];state,parentrow=extract(parent,summary['actual_accepted_parent_sha256'],'shared_parent',None)
    require(state['step']==14 and state['parameters']['rate_floor']==1.,'accepted parent rate-floor/history differs')
    protocol_path=BASE/'four/protocol.json';protocol=read(protocol_path);bindings[rel(protocol_path)]=sha(protocol_path)
    directories=[('h',ROOT/protocol['reference_directory'],2,3.90625e-6),('h_over_2',BASE/'four',4,1.953125e-6),('h_over_4',BASE/'eight',8,.9765625e-6)]
    for branch,directory,count,dt in directories:
        history_path=directory/'candidate_history.json';history=read(history_path);bindings[rel(history_path)]=sha(history_path)
        require(len(history)==count,'registered resolution count differs');previous=summary['actual_accepted_parent_sha256']
        for item in history:
            require(item['parent_sha256']==previous and item['dt_s']==dt,'candidate chain or timestep differs')
            snapshot=directory/item['checkpoint'];_,row=extract(snapshot,item['sha256'],branch,dt)
            require(row['time_s']==item['time_s'],'actual checkpoint/candidate time differs');previous=item['sha256']
            if branch!='h':
                evidence=next(s for s in summary['stages'] if s['stage']==('four' if branch=='h_over_2' else 'eight'))['PID797_saved_path'][item['step']-1]
                require(row['ep']==evidence['ep'] and row['strain_rate_per_s']==evidence['strain_rate_per_s'] and
                    abs(row['stress_frobenius_Pa']-evidence['stress_norm_Pa'])<=1e-9,'PID797 disagrees with independent readback')
    four,eight=summary['stages']
    ep_floor=four['terminal_error_reduction']['denominator_floors']['plastic_strain']
    require(ep_floor==1e-5 and eight['terminal_error_reduction']['denominator_floors']['plastic_strain']==ep_floor,'ep comparison floor differs')
    errors=[four['terminal_error_reduction']['current_errors']['temperature_max_K'],eight['terminal_error_reduction']['current_errors']['temperature_max_K']]
    require(errors[1]>errors[0],'registered nonmonotone temperature evidence changed')
    metadata=dict(classification='SAVED_CONTACT_ONSET_RESOLUTION_DIAGNOSTIC_PLOT_NOT_VALIDATION',particle_id=797,
        accepted_process_time_s=summary['actual_accepted_time_s'],accepted_process_global_step=14,
        diagnostic_end_time_s=max(row['time_s'] for row in rows),rate_floor_per_s=1.,ep_comparison_denominator_floor=ep_floor,
        q_definition='sqrt(3/2 S:S) from saved deviatoric stress, not the stress Frobenius norm',
        shared_initial_point='ACTUAL_ACCEPTED_PARENT_AT_62_5_US_NOT_ORIGINAL_T0_SCRATCH_RATE',
        curve_policy='STRAIGHT_SEGMENTS_JOIN_SAVED_CANDIDATE_POINTS_ONLY_NO_SMOOTHING_OR_SUBSTEP_RECONSTRUCTION',
        terminal_global_adjacent_temperature_errors_K=errors,temperature_errors_monotone=False,
        temporal_screens_all_pass={s['stage']:s['temporal_screens_all_pass'] for s in summary['stages']},
        source_and_input_sha256=bindings,physical_evaluations=0,nnps_evaluations=0,accepted_time_advanced_s=0.,
        rows=rows,command=[sys.executable,*sys.argv])
    with (OUT/'pid797_saved_states.csv').open('x',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    metadata['csv_sha256']=sha(OUT/'pid797_saved_states.csv')
    (OUT/'data.json').write_text(json.dumps(metadata,indent=2,allow_nan=False),encoding='utf-8')
    for name,digest in bindings.items():require(sha(ROOT/name)==digest,'plot input changed during export')
    print(json.dumps(dict(status='PASS_SAVED_ONLY_PID797_EXPORT',count=len(rows),output=rel(OUT),physical_evaluations=0)))


if __name__=='__main__':main()

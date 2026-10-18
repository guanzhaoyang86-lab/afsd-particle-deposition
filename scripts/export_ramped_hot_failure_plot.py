"""Saved-only PID-specific ramp failure evidence; never evaluates physics."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:sys.path.insert(0,str(ROOT/'src'))
from pysph.solver.utils import load

BASE=ROOT/'outputs/ramped_hot_startup_20261005/v1'
OUT=BASE/'visualization'
END=.0035


def require(value,message):
    if not bool(value):raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def main():
    begin=time.perf_counter();require(not OUT.exists(),'fresh visualization output required')
    OUT.mkdir();bindings={rel(Path(__file__)):sha(Path(__file__))}
    def read(path):
        bindings[rel(path)]=sha(path)
        return json.loads(Path(path).read_text(encoding='utf-8'))
    screen=read(BASE/'chunk_004/fine/pair_001_screen.json')
    require(not screen['passed'] and abs(screen['time_s']-END)<2e-18,'registered terminal failure differs')
    worst=max(screen['local_screen']['rows'],key=lambda row:row['stress_relative'])
    pid=int(worst['particle_id']);require(pid==791,'registered largest stress discrepancy PID changed')
    initial=BASE/'initial/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
    cost=read(BASE/'cost/manifest.json');require(sha(initial)==cost['initial_sha256'],'initial snapshot drift')
    rows=[];stored={};first_loaded={};last_sha={}
    def extract(path,digest,branch,dt):
        require(sha(path)==digest,'snapshot hash differs: '+rel(path));bindings[rel(path)]=digest
        saved=load(str(path));pa=saved['arrays']['solid'];state=saved['solver_data']
        require(state['accepted'] is False and state['restart_allowed'] is False and state['diagnostic_only'] is True,
                'diagnostic snapshot was relabelled')
        require(len(pa.x)==918 and len(np.unique(pa.particle_id))==918,'inventory/IDs changed')
        found=np.flatnonzero(pa.particle_id==pid);require(len(found)==1,'selected persistent ID missing')
        k=int(found[0]);S=pa.S.reshape(-1,3,3)[k].copy()
        row=dict(branch=branch,particle_id=pid,time_s=float(state['time_s']),dt_s=dt,
            q_Pa=float(np.sqrt(1.5*np.sum(S*S))),stress_frobenius_Pa=float(np.linalg.norm(S)),
            ep=float(pa.ep[k]),strain_rate_per_s=float(pa.strain_rate[k]),ep_rate_per_s=float(pa.ep_rate[k]),
            temperature_K=float(pa.temperature[k]),rate_is_computed=branch!='independent_initial',
            checkpoint=rel(path),checkpoint_sha256=digest)
        require(row['time_s']<=END+2e-18,'unpaired coarse endpoint must not enter visualization')
        rows.append(row);stored[(branch,row['time_s'])]=(S,float(pa.ep[k]));return state,row
    state,_=extract(initial,cost['initial_sha256'],'independent_initial',None)
    require(state['time_s']==0. and state['step']==0 and state['parameters']['rate_floor']==1.,'initial scope/rate floor differs')
    parameters=state['parameters']
    for stage,dt,count in (('coarse',.0005,7),('fine',.00025,14)):
        previous=cost['initial_sha256'];selected=0
        for chunk in range(1,5):
            folder=BASE/f'chunk_{chunk:03d}'/stage
            history=read(folder/'candidate_history.json');manifest=read(folder/'manifest.json')
            require(manifest['accepted'] is False and manifest['accepted_time_advanced_s']==0.,'study published a state')
            for item in history:
                require(item['parent_sha256']==previous and item['dt_s']==dt,'branch chain or resolution differs')
                previous=item['sha256']
                if item['time_s']>END+2e-18:continue
                current,row=extract(folder/item['checkpoint'],item['sha256'],stage,dt)
                require(row['time_s']==item['time_s'] and current['step']==selected+1,'saved time/step differs')
                require(current['parameters']==parameters,'physical parameters changed')
                response=read(folder/item['checkpoint'].split('/')[0]/'response.json')
                require(response==item['response'],'response/history differs')
                if response['nonbase_loaded_ids'] and stage not in first_loaded:first_loaded[stage]=row['time_s']
                selected+=1;last_sha[stage]=item['sha256']
        require(selected==count,'unexpected saved resolution point count')
    require(last_sha['coarse']==screen['coarse_sha256'] and last_sha['fine']==screen['fine_checkpoint_sha256'][-1],
            'terminal plot endpoints do not match failed comparison')
    require(all(abs(t-.003)<2e-18 for t in first_loaded.values()) and len(first_loaded)==2,'first saved load guide differs')
    coarse=next(r for r in rows if r['branch']=='coarse' and abs(r['time_s']-END)<2e-18)
    fine=next(r for r in rows if r['branch']=='fine' and abs(r['time_s']-END)<2e-18)
    Sa,epa=stored[('fine',fine['time_s'])];Sb,epb=stored[('coarse',coarse['time_s'])]
    stress=float(np.linalg.norm(Sa-Sb)/max(np.linalg.norm(Sa),1e6))
    plastic=abs(epa-epb)/max(abs(epa),1e-5)
    require(abs(stress-worst['stress_relative'])<=2e-15 and abs(plastic-worst['plastic_relative'])<=2e-15,
            'terminal raw particle fields disagree with screened differences')
    with (OUT/'pid791_saved_states.csv').open('x',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    metadata=dict(status='PASS_SAVED_RAMP_FAILURE_PLOT_EXPORT_NOT_PHYSICAL_VALIDATION',
        classification='ENGINEERING_APPROXIMATE_INDEPENDENT_RAMP_DIAGNOSTIC_NOT_ACCEPTED',particle_id=pid,
        selection='ARGMAX_LOCAL_TENSOR_STRESS_RELATIVE_ERROR_AT_FAILED_3_5MS_SCREEN',
        terminal_local_screen=worst,terminal_global_screen=screen['global_screen'],
        terminal_time_s=END,first_saved_tool_load_time_s_by_branch=first_loaded,
        first_saved_tool_load_is_not_exact_continuous_onset=True,
        q_definition='sqrt(3/2 S:S) from saved deviatoric stress; local screen compares full tensors',
        strain_rate_floor_per_s=1.,stress_comparison_floor_Pa=1e6,plastic_comparison_floor=1e-5,
        initial_rate_policy='SCRATCH_INITIAL_RATE_IS_UNEVALUATED_AND_NOT_PLOTTED_AS_A_RATE_RESULT',
        accepted_time_advanced_s=0.,preexisting_accepted_reference_time_s=68.359375e-6,
        curve_policy='STRAIGHT_SEGMENTS_BETWEEN_SAVED_POINTS_ONLY; NO_INTERPOLATED_PHYSICAL_FIELDS',
        unpaired_coarse_4ms_excluded=True,physical_evaluations=0,nnps_evaluations=0,
        source_and_input_sha256=bindings,rows=rows,csv_sha256=sha(OUT/'pid791_saved_states.csv'),
        command=[sys.executable,*sys.argv],wall_seconds=time.perf_counter()-begin)
    for name,digest in bindings.items():require(sha(ROOT/name)==digest,'input drift during saved-only extraction')
    (OUT/'data.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=metadata['status'],particle_id=pid,points=len(rows),wall_seconds=metadata['wall_seconds'])))


if __name__=='__main__':main()

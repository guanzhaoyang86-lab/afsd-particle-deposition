"""Saved-only review of PID893 after the first increased-dt process window.

Imports only pure file/array observations from the frozen earlier review.
No root, material update, contact evaluation or NNPS is executed.
"""
import argparse
import json
import math
from pathlib import Path
import time
import numpy as np
import review_expanded_pid797_saved as pure

ROOT=pure.ROOT
WINDOW=ROOT/'outputs/expanded_hot_campaign_v2_20261005/window_002'
PREVIOUS=ROOT/'outputs/expanded_hot_campaign_v2_20261005/window_001'
PARENT=PREVIOUS/'fine/step_006/accepted.npz'
PARENT_SHA='eed12fa1410394838358297e9519d850ea84f09e0065d85dbe89584b75caf4af'
AUDIT_SHA='d67eb4f43f016d313cf3c8483385c081aec33608170be304442a31405945345b'
PID=893


def state(saved,checkpoint):
    data=pure.load(str(saved.bind(checkpoint)));pa=data['arrays']['solid'];s=data['solver_data']
    indices=np.flatnonzero(pa.particle_id==PID);assert len(indices)==1;i=int(indices[0])
    stress=pa.S.reshape(-1,3,3)[i]
    return dict(checkpoint=checkpoint.relative_to(ROOT).as_posix(),checkpoint_sha256=pure.sha(checkpoint),
        time_s=s['time_s'],global_step=s['step'],accepted=s['accepted'],restart_allowed=s['restart_allowed'],
        ep=float(pa.ep[i]),ep_rate_per_s=float(pa.ep_rate[i]),strain_rate_per_s=float(pa.strain_rate[i]),
        temperature_K=float(pa.temperature[i]),stress_Pa=stress.tolist(),stress_norm_Pa=float(np.linalg.norm(stress)),
        equivalent_stress_Pa=float(np.sqrt(1.5*np.sum(stress*stress))),
        symmetric_gradient_per_s=pa.D.reshape(-1,3,3)[i].tolist(),
        position_m=[float(pa.x[i]),float(pa.y[i]),float(pa.z[i])]),i,s


def with_journal(saved,checkpoint):
    fields,index,s=state(saved,checkpoint)
    journal=saved.read(checkpoint.parent/'FULL_UNCOMMITTED_candidate.json')
    saved.bind(checkpoint.parent/'FULL_UNCOMMITTED_candidate.npz')
    fields.update(contact=pure.contact_row(journal,index),fresh_metrics=journal['metrics'],
        raw_energy={k:v for k,v in journal['energy'].items() if not isinstance(v,(dict,list))},
        particle_material_energy={k:v[index] for k,v in journal['energy']['material_breakdown']['per_particle'].items()})
    return fields,s


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True);args=parser.parse_args()
    out=(ROOT/args.output).resolve();assert out.is_relative_to(ROOT/'outputs') and not out.exists()
    started=time.perf_counter();saved=pure.Saved();saved.bind(Path(__file__));saved.bind(Path(pure.__file__))
    audit=saved.read(WINDOW/'readback/summary.json');assert pure.sha(WINDOW/'readback/summary.json')==AUDIT_SHA
    assert audit['status']=='PASS_EXPANDED_HOT_CAMPAIGN_V2_READBACK_NOT_PHYSICAL_VALIDATION'
    assert audit['added_accepted_steps']==0 and audit['accepted_global_steps']==20 and not audit['window_completed']
    protocol=saved.read(WINDOW/'window_protocol.json');assert protocol['parent_sha256']==PARENT_SHA==pure.sha(PARENT)
    params=protocol['parameters'];material=ROOT/'src/afsd_backend/pysph/material.py';saved.bind(material)
    assert protocol['source_sha256'][material.relative_to(ROOT).as_posix()]==pure.sha(material)
    parent,parent_state=with_journal(saved,PARENT);records={}
    for name,relative in [('coarse_1','coarse/step_001'),('fine_1','fine/step_001'),('fine_2','fine/step_002')]:
        sd=WINDOW/relative;fields,s=with_journal(saved,sd/'NOT_RESTARTABLE_candidate.npz')
        assert not s['accepted'] and not s['restart_allowed']
        assert s['interval_energy']['records'][:20]==parent_state['interval_energy']['records']
        before=records['fine_1'] if name=='fine_2' else parent
        dt=fields['time_s']-before['time_s'];dep=fields['ep']-before['ep']
        assert math.isclose(dep/dt,fields['ep_rate_per_s'],rel_tol=1e-12,abs_tol=1e-14)
        fields.update(dt_s=dt,ep_increment=dep)
        fields['scalar_saved_field_FKS_interpretation']=pure.scalar_fks(fields['ep'],fields['strain_rate_per_s'],
            .5*(before['temperature_K']+fields['temperature_K']),params)
        fields['scalar_saved_field_FKS_interpretation']['limitation']='Saved-field scalar formula with mean endpoint temperature; not a constitutive return-map replay.'
        graph=saved.read(sd/'visible_graph.json');saved.bind(sd/'visible_graph.npz')
        fields['graph']=dict(fingerprint=graph['graph_fingerprint'],raw_edges=graph['diagnostics']['raw_edge_count'],
            visible_edges=graph['diagnostics']['visible_edge_count'],blocked_edges=graph['diagnostics']['blocked_edge_count'])
        fields['root_attempts']=[pure.root_log(saved,sd/'broyden')]
        if (sd/'newton_retry/iterations.json').exists():fields['root_attempts'].append(pure.root_log(saved,sd/'newton_retry'))
        records[name]=fields
    accepted_path=[]
    previous_protocol=saved.read(PREVIOUS/'window_protocol.json')
    previous_parent=ROOT/previous_protocol['parent_checkpoint'];fields,_=with_journal(saved,previous_parent);accepted_path.append(fields)
    for index in range(1,7):
        fields,_=with_journal(saved,PREVIOUS/f'fine/step_{index:03d}/accepted.npz')
        assert fields['accepted'] and fields['restart_allowed'];accepted_path.append(fields)
    current=saved.read(WINDOW/'fine/pair_001_temporal_screen.json')
    earlier=saved.read(PREVIOUS/'fine/pair_003_temporal_screen.json')
    row=next(v for v in current['local_screen']['rows'] if v['particle_id']==PID)
    coarse,fine=records['coarse_1'],records['fine_2']
    stress_error=float(np.linalg.norm(np.asarray(coarse['stress_Pa'])-np.asarray(fine['stress_Pa'])))
    ep_error=abs(coarse['ep']-fine['ep']);stress_den=max(fine['stress_norm_Pa'],1e6);ep_den=max(fine['ep'],1e-5)
    assert math.isclose(row['stress_relative'],stress_error/stress_den,rel_tol=1e-14)
    assert math.isclose(row['plastic_relative'],ep_error/ep_den,rel_tol=1e-14)
    costs=[]
    for stage in ('coarse','fine'):
        manifest=saved.read(WINDOW/stage/'manifest.json');rows=saved.read(WINDOW/stage/'candidate_history.json')
        costs.append(dict(stage=stage,status=manifest['status'],complete_candidate_count=len(rows),
            candidate_costs_s=[v['wall_seconds'] for v in rows],stage_wall_s=manifest['wall_seconds'],accepted_steps=manifest['accepted_steps']))
    report=dict(status='SAVED_ONLY_GROWTH_FAILURE_OBSERVATION_NOT_PHYSICAL_VALIDATION',particle_id=PID,
        physical_solver_evaluations=0,nnps_evaluations=0,contact_evaluations=0,constitutive_update_evaluations=0,
        accepted=False,restart_allowed=False,new_accepted_steps=0,actual_accepted_global_steps=20,
        actual_accepted_time_s=parent['time_s'],parent=parent,candidates=records,previous_accepted_path=accepted_path,
        earlier_passing_local_screen=next(v for v in earlier['local_screen']['rows'] if v['particle_id']==PID),
        current_failed_local_screen=row,current_global_screen=current['global_screen'],
        denominators=dict(stress_error_Pa=stress_error,stress_denominator_Pa=stress_den,stress_floor_active=False,
            ep_error=ep_error,ep_denominator=ep_den,ep_floor_active=fine['ep']<1e-5,
            ep_allowed_absolute_difference=.2*ep_den,ep_relative_without_floor=ep_error/fine['ep']),
        complete_accepted_energy_prefix=parent_state['interval_energy'],stage_costs=costs,
        total_two_stage_wall_s=math.fsum(v['stage_wall_s'] for v in costs),
        interpretation=[
            'The attempted increase from 1.953125/0.9765625 us to 3.90625/1.953125 us failed its first temporal comparison; this is not evidence that the previous smaller-step continuation would fail.',
            'PID893 already had sliding BORE_Y_MINUS contact in the true parent. All three stored candidate intervals retain that face and sliding status; this is not a newly created contact at the failed window start.',
            'Coarse and fine1 plastify, while fine2 has exactly zero additional ep. All saved effective rates are positive and above the explicit 1/s rate floor, with no observed effective-rate zero crossing.',
            'The plastic comparison floor is active. It enlarges the denominator relative to the actual fine ep, moderating the reported difference; removing it would increase the error percentage. The retained absolute difference still exceeds the registered limit.',
            'Unchanged saved graph and passing final force gates do not establish root uniqueness or isolate a sole causal mechanism.',
            'The previous passing PID893 screen uses a different interval and actual smaller dt; it cannot retrospectively validate the failed increase.',
            'These stages added zero accepted time. Their costs and this failed attempt must be retained in method and affordability review; no automatic smaller-step retry or next window is authorized by this report.'
        ])
    for name,value in saved.files.items():assert pure.sha(ROOT/name)==value
    report.update(input_and_source_sha256=saved.files,wall_seconds=time.perf_counter()-started)
    out.mkdir(parents=True,exist_ok=False);(out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('status','new_accepted_steps','actual_accepted_time_s','denominators','total_two_stage_wall_s','wall_seconds')},indent=2))
    return 0


if __name__=='__main__':raise SystemExit(main())

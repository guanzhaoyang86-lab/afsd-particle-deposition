"""Bounded saved-only forensic of the failed 3.5 ms smooth-start comparison."""
import argparse
import json
import math
from pathlib import Path
import time
import numpy as np
import review_expanded_pid797_saved as pure

ROOT=pure.ROOT
BASE=ROOT/'outputs/ramped_hot_startup_20261005/v1'
AUDIT_SHA='25746d994a5b2d98c68c3e6929dee005b507e3f3527f52ae23ed154b815a50be'
CLASS='INDEPENDENT_RAMPED_20MS_HOT_QS_VISIBLE_DIAGNOSTIC_V1'


def load(saved,path):
    value=pure.load(str(saved.bind(path)))
    return value['arrays']['solid'],value['solver_data']


def particle_fields(pa,state,pid):
    indices=np.flatnonzero(pa.particle_id==pid);assert len(indices)==1;i=int(indices[0])
    stress=pa.S.reshape(-1,3,3)[i];rate=float(pa.strain_rate[i])
    return dict(particle_id=pid,index=i,component_id=int(pa.component_id[i]),actuator=int(pa.actuator[i]),
        free_axes=np.asarray(state['free'])[i].tolist(),position_m=[float(pa.x[i]),float(pa.y[i]),float(pa.z[i])],
        stress_Pa=stress.tolist(),stress_norm_Pa=float(np.linalg.norm(stress)),
        equivalent_stress_Pa=float(np.sqrt(1.5*np.sum(stress**2))),ep=float(pa.ep[i]),ep_rate_per_s=float(pa.ep_rate[i]),
        strain_rate_per_s=rate,rate_floor_active=rate<=state['parameters']['rate_floor'],temperature_K=float(pa.temperature[i]),
        symmetric_gradient_per_s=pa.D.reshape(-1,3,3)[i].tolist())


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True);args=parser.parse_args()
    out=(ROOT/args.output).resolve();assert out.is_relative_to(ROOT/'outputs') and not out.exists()
    started=time.perf_counter();saved=pure.Saved();saved.bind(Path(__file__));saved.bind(Path(pure.__file__))
    audit_path=BASE/'chunk_004/readback/summary.json';audit=saved.read(audit_path);assert pure.sha(audit_path)==AUDIT_SHA
    assert audit['status']=='PASS_RAMP_CHUNK_SAVED_READBACK_NOT_ACCEPTED' and not audit['temporal_screens_all_pass']
    assert not audit['loaded_feasibility_screen_passed'] and audit['accepted_time_advanced_s']==0.
    for name,digest in audit['evidence_sha256'].items():assert pure.sha(saved.bind(ROOT/name))==digest
    screen=saved.read(BASE/'chunk_004/fine/pair_001_screen.json')
    rows=screen['local_screen']['rows'];maximum_stress=max(rows,key=lambda v:v['stress_relative'])
    maximum_ep=max(rows,key=lambda v:v['plastic_relative']);ids=sorted({maximum_stress['particle_id'],maximum_ep['particle_id']})
    initial,initial_state=load(saved,BASE/'initial/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz')
    p=initial_state['parameters'];saved.bind(ROOT/'src/afsd_backend/pysph/material.py')
    histories={};costs=[];endpoints={}
    for stage in ('coarse','fine'):
        previous,ps=initial,initial_state;history=[]
        for chunk in range(1,5):
            folder=BASE/f'chunk_{chunk:03d}'/stage
            manifest=saved.read(folder/'manifest.json');supervision=saved.read(folder/'supervision.json')
            candidates=saved.read(folder/'candidate_history.json')
            costs.append(dict(chunk=chunk,stage=stage,status=manifest['status'],candidate_count=len(candidates),
                candidate_wall_seconds=[v['wall_seconds'] for v in candidates],stage_wall_seconds=manifest['wall_seconds'],
                supervised_wall_seconds=supervision['wall_seconds'],forecast_s=manifest['forecast_s']))
            for row in candidates:
                path=folder/row['checkpoint'];pa,state=load(saved,path);assert pure.sha(path)==row['sha256']
                assert state['classification']==CLASS and not state['accepted'] and not state['restart_allowed'] and state['diagnostic_only']
                assert state['interval_energy']['records'][:-1]==ps['interval_energy']['records']
                sd=path.parent;j=saved.read(sd/'FULL_UNCOMMITTED_candidate.json');saved.bind(sd/'FULL_UNCOMMITTED_candidate.npz')
                graph=saved.read(sd/'visible_graph.json');saved.bind(sd/'visible_graph.npz')
                dt=row['dt_s'];fields=[]
                for pid in ids:
                    f=particle_fields(pa,state,pid);before=particle_fields(previous,ps,pid);i=f['index']
                    dep=f['ep']-before['ep'];assert math.isclose(dep/dt,f['ep_rate_per_s'],rel_tol=1e-12,abs_tol=1e-14)
                    scalar=pure.scalar_fks(f['ep'],f['strain_rate_per_s'],.5*(f['temperature_K']+before['temperature_K']),p)
                    scalar['limitation']='Arithmetic from saved ep/rate and mean endpoint temperature, not a constitutive update or return-map replay.'
                    f.update(ep_increment=dep,scalar_saved_field_FKS_interpretation=scalar,
                        q_minus_inferred_yield_Pa=f['equivalent_stress_Pa']-scalar['scalar_yield_Pa'],
                        rate_floor_crossed_from_previous=(before['strain_rate_per_s']<=p['rate_floor'])!=(f['strain_rate_per_s']<=p['rate_floor']),
                        contact=pure.contact_row(j,i),particle_material_energy={k:v[i] for k,v in j['energy']['material_breakdown']['per_particle'].items()})
                    fields.append(f)
                attempts=[]
                for sub in ('broyden','newton_retry'):
                    if (sd/sub/'iterations.json').exists():attempts.append(pure.root_log(saved,sd/sub))
                record=dict(time_s=state['time_s'],global_step=state['step'],dt_s=dt,checkpoint=path.relative_to(ROOT).as_posix(),
                    checkpoint_sha256=row['sha256'],particles=fields,response=row['response'],root_attempts=attempts,
                    controller=row['controller'],fallback_reason=row['fallback_reason'],candidate_wall_seconds=row['wall_seconds'],
                    fresh_metrics=j['metrics'],raw_energy={k:v for k,v in j['energy'].items() if not isinstance(v,(dict,list))},
                    complete_interval_energy=state['interval_energy'],graph_fingerprint=graph['graph_fingerprint'],
                    raw_edge_count=graph['diagnostics']['raw_edge_count'],visible_edge_count=graph['diagnostics']['visible_edge_count'],
                    blocked_edge_count=graph['diagnostics']['blocked_edge_count'],motion=saved.read(sd/'motion.json'))
                assert j['metrics']['scaled_residual']<=1 and abs(j['energy']['global_raw_residual_J'])<=j['energy']['registered_energy_limit_J']
                history.append(record)
                if abs(state['time_s']-.0035)<1e-15:endpoints[stage]=pa
                previous,ps=pa,state
        histories[stage]=history
    denominators=[]
    for pid in ids:
        c,f=(particle_fields(endpoints[s],initial_state,pid) for s in ('coarse','fine'))
        se=float(np.linalg.norm(np.asarray(c['stress_Pa'])-np.asarray(f['stress_Pa'])));ee=abs(c['ep']-f['ep'])
        sd=max(f['stress_norm_Pa'],1e6);ed=max(f['ep'],1e-5);row=next(v for v in rows if v['particle_id']==pid)
        assert math.isclose(se/sd,row['stress_relative'],rel_tol=1e-14) and math.isclose(ee/ed,row['plastic_relative'],rel_tol=1e-14)
        denominators.append(dict(particle_id=pid,stress_error_Pa=se,stress_denominator_Pa=sd,stress_floor_active=f['stress_norm_Pa']<1e6,
            ep_error=ee,ep_denominator=ed,ep_floor_active=f['ep']<1e-5,stress_relative=se/sd,plastic_relative=ee/ed,
            stress_allowed_absolute_difference_Pa=.1*sd,ep_allowed_absolute_difference=.2*ed))
    c,f=endpoints['coarse'],endpoints['fine'];norm=np.linalg.norm
    global_denominators=dict(stress_error_norm_Pa=float(norm(f.S-c.S)),fine_stress_norm_Pa=float(norm(f.S)),
        stress_denominator_Pa=float(max(norm(f.S),1e6*np.sqrt(len(f.x)))),stress_floor_active=bool(norm(f.S)<1e6*np.sqrt(len(f.x))),
        plastic_error_norm=float(norm(f.ep-c.ep)),fine_plastic_norm=float(norm(f.ep)),
        plastic_denominator=float(max(norm(f.ep),1e-5*np.sqrt(len(f.x)))),plastic_floor_active=bool(norm(f.ep)<1e-5*np.sqrt(len(f.x))))
    for name,digest in saved.files.items():assert pure.sha(ROOT/name)==digest
    report=dict(status='SAVED_ONLY_RAMP_POST_CONTACT_FAILURE_OBSERVATION_NOT_PHYSICAL_VALIDATION',
        physical_solver_evaluations=0,nnps_evaluations=0,contact_evaluations=0,constitutive_update_evaluations=0,
        accepted=False,restart_allowed=False,new_accepted_steps=0,accepted_time_advanced_s=0.,
        maximum_stress_row=maximum_stress,maximum_plastic_row=maximum_ep,failed_global_screen=screen['global_screen'],
        local_denominators=denominators,global_denominators=global_denominators,histories=histories,
        counts=dict(stress_local_failed=sum(v['stress_relative']>.1 for v in rows),plastic_local_failed=sum(v['plastic_relative']>.2 for v in rows)),
        top_stress_rows=sorted(rows,key=lambda v:v['stress_relative'],reverse=True)[:10],
        top_plastic_rows=sorted(rows,key=lambda v:v['plastic_relative'],reverse=True)[:10],stage_costs=costs,
        total_root_stage_wall_seconds=math.fsum(v['stage_wall_seconds'] for v in costs),
        total_root_stage_supervised_seconds=math.fsum(v['supervised_wall_seconds'] for v in costs),
        total_experiment_supervised_seconds=audit['total_supervised_wall_seconds'],
        limitations=['The comparison is between two continuous histories from a common independent t=0, not the same immediate 3 ms parent.',
            'Saved force convergence and raw/full-prefix energy gates are distinct from temporal accuracy.',
            'Inferred scalar yield is saved-field arithmetic; it does not replay material evolution or establish a unique failure cause.',
            'Coarse 4 ms is retained unused diagnostic evidence after the 3.5 ms fine comparison failed.',
            'No automatic refinement, controller change, acceptance, or new physical evaluation is authorized by this observation.'],
        input_and_source_sha256=saved.files,wall_seconds=time.perf_counter()-started)
    out.mkdir(parents=True,exist_ok=False);(out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('status','maximum_stress_row','maximum_plastic_row','local_denominators','counts','total_experiment_supervised_seconds','wall_seconds')},indent=2));return 0


if __name__=='__main__':raise SystemExit(main())

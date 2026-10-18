"""Read-only saved-state diagnosis of the failed seventh expanded budget pair."""
from pathlib import Path
import argparse
import json
import numpy as np
import readback_expanded_hot_budget_prefix as rb


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    output=rb.project_path(args.output)
    rb.require(not output.exists(),'refuse diagnostic evidence overwrite')
    output.mkdir(parents=True,exist_ok=False)
    reader=rb.Reader();reader.bind(Path(__file__));reader.bind(Path(rb.__file__))
    base=rb.ROOT/'outputs/expanded_hot_budget_prefix_20261005/window_001'
    summary=reader.read(base/'readback/summary.json')
    rb.require(summary['status']==rb.PASS_READBACK and summary['added_accepted_steps']==12 and
        summary['accepted_global_steps']==14 and summary['budget_prefix_completed'] is False,
        'diagnostic requires the actual audited failed-budget prefix')
    reader.sources(summary['evidence_sha256'])
    screen=reader.read(base/'fine/pair_007_temporal_screen.json')
    stress=max(screen['local_screen']['rows'],key=lambda r:r['stress_relative'])
    plastic=max(screen['local_screen']['rows'],key=lambda r:r['plastic_relative'])
    rb.require(stress['particle_id']==plastic['particle_id']==797 and
        screen['global_screen']['pass'] and not screen['local_screen']['passed'],
        'saved failure differs from the registered diagnostic target')
    cases=[('coarse_pair6',rb.COARSE_WINDOW/'coarse/step_006'),
           ('fine_true_accepted_pair6',base/'fine/step_012'),
           ('coarse_pair7_not_accepted',rb.COARSE_WINDOW/'coarse/step_007'),
           ('fine_pair7_half1_not_accepted',base/'fine/step_013'),
           ('fine_pair7_half2_not_accepted',base/'fine/step_014')]
    rows=[]
    for name,folder in cases:
        published=(folder/'accepted.npz').exists()
        checkpoint=folder/('accepted.npz' if published else 'NOT_RESTARTABLE_candidate.npz')
        pa,state=reader.checkpoint(checkpoint)
        journal=reader.npz(folder/'FULL_UNCOMMITTED_candidate.npz')
        report=reader.read(folder/'FULL_UNCOMMITTED_candidate.json')
        graph=reader.npz(folder/'visible_graph.npz')
        k=int(np.flatnonzero(pa.particle_id==797)[0]); contact=report['contact']
        incident=(graph['i']==k)|(graph['j']==k)
        neighbor=np.unique(np.r_[graph['i'][incident],graph['j'][incident]])
        neighbor=neighbor[neighbor!=k]
        contact_rows=[]
        for f,feature in enumerate(contact['features']):
            contact_rows.append(dict(feature=feature,**{key:np.asarray(contact[key])[k,f]
                for key in ('active','friction_active','sliding','gaps_m','overlap_m','penetration_m',
                    'normal_load_N','endpoint_normal_load_N','force_by_feature_N',
                    'normal_force_by_feature_N','tangential_force_by_feature_N',
                    'normal_potential0_J','normal_potential1_J','tangential_potential0_J',
                    'tangential_potential1_J','friction_dissipation_J','algorithmic_dissipation_J',
                    'release_energy_J')}))
        s=pa.S.reshape(-1,3,3)[k]
        rows.append(dict(name=name,checkpoint=checkpoint.relative_to(rb.ROOT).as_posix(),
            checkpoint_sha256=rb.digest(checkpoint),published=published,accepted=state['accepted'],
            restart_allowed=state['restart_allowed'],time_s=state['time_s'],global_step=state['step'],
            particle_id=797,actuator=pa.actuator[k],free_axes=np.asarray(state['free'])[k],
            position_m=rb.positions(pa)[k],velocity_m_s=np.array([pa.u[k],pa.v[k],pa.w[k]]),
            stress_Pa=s,stress_frobenius_Pa=np.linalg.norm(s),
            equivalent_stress_Pa=np.sqrt(1.5*np.sum(s*s)),plastic_strain=pa.ep[k],
            plastic_work_J=pa.plastic_work[k],temperature_K=pa.temperature[k],density_kg_m3=pa.rho[k],
            saved_force_N=journal['force_N'][k],saved_internal_force_N=journal['internal_force_N'][k],
            saved_contact_force_N=journal['contact_force_N'][k],complete_map_metrics=report['metrics'],
            contact_by_feature=contact_rows,retained_neighbor_ids=pa.particle_id[neighbor],
            retained_actuator_neighbor_ids=pa.particle_id[neighbor][pa.actuator[neighbor]>0],
            neighbor_plastic_strain=pa.ep[neighbor],
            global_raw_residual_J=state['energy']['global_raw_residual_J']))
    reader.finish()
    result=dict(status='SAVED_FAILED_PAIR_LOCAL_RESPONSE_DIAGNOSTIC_NOT_PHYSICAL_VALIDATION',
        physical_solver_evaluations=0,nnps_evaluations=0,accepted_time_added_s=0.,
        original_failure_preserved=True,accepted_parent=summary['accepted_history'][-1],
        global_screen=screen['global_screen'],maximum_local_stress_row=stress,
        maximum_local_plastic_row=plastic,saved_states=rows,evidence_sha256=reader.evidence,
        limitations=['Saved force components and material fields are observations, not a new force or constitutive solve.',
            'A local coarse/fine discrepancy identifies failed temporal agreement; it does not by itself establish a unique cause.',
            'Neither failed half of pair seven is accepted or restartable.'])
    (output/'summary.json').write_text(json.dumps(rb.clean(result),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],output=str(output/'summary.json'),
        sha256=rb.digest(output/'summary.json'))))


if __name__=='__main__':main()

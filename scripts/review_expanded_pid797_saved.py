"""Saved-only observations of the failed same-parent contact-onset comparison.

No force, NNPS, contact, or constitutive update is called. The scalar FKS value
below is arithmetic interpretation of saved fields, not a return-map replay.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from pysph.solver.utils import load

WINDOW = ROOT / 'outputs/expanded_hot_campaign_20261005/window_001'
PARENT = ROOT / 'outputs/expanded_hot_budget_prefix_20261005/window_001/fine/step_012/accepted.npz'
PARENT_SHA = '12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b'
AUDIT_SHA = '89f999f398a9f3de70535500121bbba3a359ab3b39cbd351a0e8b00960e88658'
PID = 797


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Saved:
    def __init__(self):
        self.files = {}

    def bind(self, path):
        path = path.resolve()
        assert path.is_relative_to(ROOT) and path.is_file()
        name = path.relative_to(ROOT).as_posix()
        digest = sha(path)
        assert self.files.get(name, digest) == digest
        self.files[name] = digest
        return path

    def read(self, path):
        return json.loads(self.bind(path).read_text(encoding='utf-8'))

    def state(self, path):
        data = load(str(self.bind(path)))
        pa, state = data['arrays']['solid'], data['solver_data']
        indices = np.flatnonzero(pa.particle_id == PID)
        assert len(indices) == 1
        i = int(indices[0])
        stress = pa.S.reshape(-1, 3, 3)[i]
        fields = dict(checkpoint=path.relative_to(ROOT).as_posix(), checkpoint_sha256=sha(path),
            particle_id=PID, particle_count=len(pa.x), particle_index=i,
            time_s=state['time_s'], global_step=state['step'],
            accepted=state['accepted'], restart_allowed=state['restart_allowed'],
            ep=float(pa.ep[i]), ep_rate_per_s=float(pa.ep_rate[i]),
            strain_rate_per_s=float(pa.strain_rate[i]), temperature_K=float(pa.temperature[i]),
            deviatoric_stress_Pa=stress.tolist(), stress_norm_Pa=float(np.linalg.norm(stress)),
            equivalent_stress_Pa=float(np.sqrt(1.5*np.sum(stress*stress))),
            symmetric_gradient_per_s=pa.D.reshape(-1, 3, 3)[i].tolist(),
            position_m=[float(pa.x[i]), float(pa.y[i]), float(pa.z[i])],
            interval_energy={k:v for k,v in state['interval_energy'].items() if k != 'records'})
        return fields, i, state


def scalar_fks(ep, rate, temperature, p):
    ts = (temperature-p['room'])/(p['melt']-p['room'])
    rate_factor = 1+p['b1']*ts**p['b2']*math.log(max(rate,p['rate_floor'])/p['reference_rate'])
    theta = -math.expm1(-math.log1p(math.exp(-p['c1']*ts))/p['c2'])
    return dict(temperature_K=temperature, normalized_temperature=ts,
        rate_floor_per_s=p['rate_floor'], rate_factor=rate_factor, thermal_factor=theta,
        scalar_yield_Pa=(p['a1']+p['a2']*math.atan(p['a3']*ep))*rate_factor*theta)


def contact_row(journal, index):
    contact = journal['contact']
    keys = ['active','sliding','friction_active','normal_load_N','endpoint_normal_load_N',
        'overlap_m','normal_potential0_J','normal_potential1_J',
        'tangential_potential0_J','tangential_potential1_J','friction_dissipation_J',
        'algorithmic_dissipation_J','release_energy_J','normal_particle_work_J',
        'normal_wall_work_J','particle_work_J','wall_work_J','relative_step_body_m',
        'force_by_feature_N','normal_force_by_feature_N','tangential_force_by_feature_N','history']
    rows = {}
    for j, feature in enumerate(contact['features']):
        rows[feature] = {key:contact[key][index][j] for key in keys}
    return rows


def root_log(saved, directory):
    rows = saved.read(directory/'iterations.json')
    fingerprints = set()
    changed = 0
    for row in rows:
        if 'graph_fingerprint' in row:
            fingerprints.add(row['graph_fingerprint'])
        for trial in row.get('line_search', []):
            if 'graph_fingerprint' in trial:
                fingerprints.add(trial['graph_fingerprint'])
            changed += int(trial.get('graph_changed', False))
    result = dict(log=(directory/'iterations.json').relative_to(ROOT).as_posix(),
        row_count=len(rows), final_iteration=rows[-1]['iteration'],
        final_evaluations=rows[-1].get('evaluations_at_failure',rows[-1]['evaluations']),
        initial_force_residual_N=rows[0]['force_residual_N'],
        final_force_residual_N=rows[-1]['force_residual_N'],
        final_scaled_residual=rows[-1]['scaled_residual'],
        final_force_gate_pass=rows[-1]['force_residual_N'] <= rows[-1]['force_limit_N'],
        recorded_graph_fingerprints=sorted(fingerprints), recorded_graph_change_count=changed,
        residual_sequence_N=[r['force_residual_N'] for r in rows])
    failure = directory/'secant_failure.json'
    if failure.exists():
        result['retained_failure'] = saved.read(failure)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    output=(ROOT/args.output).resolve()
    assert output.is_relative_to(ROOT/'outputs') and not output.exists()
    started=time.perf_counter(); saved=Saved(); saved.bind(Path(__file__))
    protocol=saved.read(WINDOW/'window_protocol.json')
    audit=saved.read(WINDOW/'readback/summary.json')
    assert sha(WINDOW/'readback/summary.json') == AUDIT_SHA
    assert audit['status']=='PASS_EXPANDED_HOT_CAMPAIGN_READBACK_NOT_PHYSICAL_VALIDATION'
    assert audit['added_accepted_steps']==0 and audit['accepted_global_steps']==14
    assert sha(PARENT)==PARENT_SHA and protocol['parent_sha256']==PARENT_SHA
    material_source=ROOT/'src/afsd_backend/pysph/material.py'
    saved.bind(material_source)
    assert protocol['source_sha256'][material_source.relative_to(ROOT).as_posix()]==sha(material_source)
    params=protocol['parameters']
    parent, parent_index, parent_state=saved.state(PARENT)
    parent_journal=saved.read(PARENT.parent/'FULL_UNCOMMITTED_candidate.json')
    parent['last_increment_contact']=contact_row(parent_journal,parent_index)
    records={}
    prior_fields={'coarse_1':parent, 'fine_1':parent}
    for name, rel in [('coarse_1','coarse/step_001'),('fine_1','fine/step_001'),('fine_2','fine/step_002')]:
        directory=WINDOW/rel
        fields,index,state=saved.state(directory/'NOT_RESTARTABLE_candidate.npz')
        assert not fields['accepted'] and not fields['restart_allowed']
        journal=saved.read(directory/'FULL_UNCOMMITTED_candidate.json')
        saved.bind(directory/'FULL_UNCOMMITTED_candidate.npz')
        graph=saved.read(directory/'visible_graph.json')
        saved.bind(directory/'visible_graph.npz')
        fields['contact']=contact_row(journal,index)
        fields['fresh_metrics']=journal['metrics']
        fields['raw_energy']={k:v for k,v in journal['energy'].items() if not isinstance(v,(dict,list))}
        fields['particle_material_energy']={k:v[index] for k,v in journal['energy']['material_breakdown']['per_particle'].items()}
        fields['graph']=dict(fingerprint=graph['graph_fingerprint'],
            raw_edges=graph['diagnostics']['raw_edge_count'],visible_edges=graph['diagnostics']['visible_edge_count'],
            blocked_edges=graph['diagnostics']['blocked_edge_count'])
        before=records['fine_1'] if name=='fine_2' else prior_fields[name]
        dt=fields['time_s']-before['time_s']; dep=fields['ep']-before['ep']
        assert math.isclose(dep/dt, fields['ep_rate_per_s'],rel_tol=1e-12,abs_tol=1e-14)
        fields['dt_s']=dt; fields['ep_increment']=dep
        fields['scalar_saved_field_FKS_interpretation']=scalar_fks(fields['ep'],fields['strain_rate_per_s'],
            .5*(before['temperature_K']+fields['temperature_K']),params)
        fields['scalar_saved_field_FKS_interpretation']['limitation']='Arithmetic from saved endpoint ep/rate and mean endpoint temperature; midpoint temperature is inferred within saved thermal fixed-point residual, not a constitutive update or exact replay.'
        fields['root_attempts']=[root_log(saved,directory/'broyden')]
        if (directory/'newton_retry/iterations.json').exists():
            fields['root_attempts'].append(root_log(saved,directory/'newton_retry'))
        records[name]=fields
    screen=saved.read(WINDOW/'fine/pair_001_temporal_screen.json')
    oldscreen=saved.read(ROOT/'outputs/expanded_hot_budget_prefix_20261005/window_001/fine/pair_007_temporal_screen.json')
    row=next(r for r in screen['local_screen']['rows'] if r['particle_id']==PID)
    coarse,fine=records['coarse_1'],records['fine_2']
    stress_error=float(np.linalg.norm(np.asarray(coarse['deviatoric_stress_Pa'])-np.asarray(fine['deviatoric_stress_Pa'])))
    ep_error=abs(coarse['ep']-fine['ep'])
    stress_den=max(fine['stress_norm_Pa'],1e6); ep_den=max(fine['ep'],1e-5)
    assert math.isclose(stress_error/stress_den,row['stress_relative'],rel_tol=1e-14)
    assert math.isclose(ep_error/ep_den,row['plastic_relative'],rel_tol=1e-14)
    report=dict(status='SAVED_ONLY_METHOD_OBSERVATION_COMPLETE_NOT_PHYSICAL_VALIDATION',
        physical_solver_evaluations=0,nnps_evaluations=0,contact_evaluations=0,constitutive_update_evaluations=0,
        accepted=False,restart_allowed=False,new_accepted_steps=0,actual_accepted_global_steps=14,
        actual_accepted_time_s=parent['time_s'],particle_id=PID,parent=parent,candidates=records,
        current_global_screen=screen['global_screen'],current_PID797_local_screen=row,
        previous_long_window_PID797_local_screen=next(r for r in oldscreen['local_screen']['rows'] if r['particle_id']==PID),
        denominator_check=dict(stress_error_Pa=stress_error,stress_denominator_Pa=stress_den,
            stress_comparison_floor_Pa=1e6,stress_floor_active=fine['stress_norm_Pa']<=1e6,
            ep_error=ep_error,ep_denominator=ep_den,ep_comparison_floor=1e-5,ep_floor_active=fine['ep']<=1e-5),
        interpretation=[
            'Same-parent first-pair failure persists after shortening the window; earlier long-window accumulated coarse/fine drift is not a sufficient explanation.',
            'All saved effective strain rates are positive. The parent is below the explicit FKS 1/s rate floor, and the trial candidates are above it; no effective-rate zero crossing is observed.',
            'Coarse PID797 plastifies; fine1 plastifies; fine2 has exactly zero additional ep. This is a different stored elastic/plastic history, not evidence that a particular mechanism is uniquely responsible.',
            'BORE_Y_PLUS is active and sliding in all three candidate journals; no reversal of the saved tangential x/z increments is observed. Saved endpoints do not prove continuous contact history.',
            'Comparison denominator floors are inactive at PID797; they do not account for this failed screen.',
            'Final roots pass the registered force gates. Fine1 exhausted Broyden and passed the registered Newton retry. These logs neither certify root uniqueness nor eliminate solver-tolerance influence.',
            'One pre-registered bounded time refinement study is justified to test resolution of the coupled contact/rate/material response. Convergence or eventual acceptance is not presumed.',
            'A 4-substep and at most one 8-substep branch from the same true parent, compared at all shared endpoints with all original gates, can distinguish improving temporal resolution from persistent/nonmonotone branch sensitivity. Stop at that bound and review the method if unresolved.'
        ],
        limitations=['No new physical calculation, publication, or accepted time.',
            'These observations do not establish contact-history, whole-history, spatial, or experimental convergence.',
            'No coefficient, floor, tolerance, history, or energy ledger has been changed.'])
    for name,digest in saved.files.items():
        assert sha(ROOT/name)==digest, 'input changed during saved-only observation'
    report['input_and_source_sha256']=saved.files
    report['wall_seconds']=time.perf_counter()-started
    output.mkdir(parents=True,exist_ok=False)
    (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['status','physical_solver_evaluations','new_accepted_steps','actual_accepted_time_s','denominator_check','wall_seconds']},indent=2))
    return 0


if __name__=='__main__':
    raise SystemExit(main())

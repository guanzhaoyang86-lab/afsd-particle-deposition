"""Read terminal saved arrays only; no full trial, material or contact calls."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
import numpy as np
from pysph.solver.utils import load

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); case = args.case.resolve(); out = args.output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('fresh project-local report path required')
    manifest = json.loads((case/'manifest.json').read_text(encoding='utf-8'))
    if 'RUNNING' in manifest['status']:
        raise ValueError('terminal manifest required; not a monitoring/restart script')
    paths = [case/'manifest.json', case/'iterations.json', case/'latest_trial.json.gz',
             case/'resolved_config.json']
    watched = {str(path): sha(path) for path in paths}
    for rel, expected in manifest['source_sha256'].items():
        if sha(ROOT/rel) != expected:
            raise ValueError('source drift '+rel)
    trace = json.loads((case/'iterations.json').read_text(encoding='utf-8'))
    saved = json.loads(gzip.decompress((case/'latest_trial.json.gz').read_bytes()))
    evaluation, candidate = saved['evaluation'], saved['candidate']
    p = json.loads((case/'resolved_config.json').read_text(encoding='utf-8'))['parameters']
    parent = ROOT/manifest['parent_reference']['path']
    if sha(parent) != manifest['parent_reference']['sha256']:
        raise ValueError('parent snapshot drift')
    watched[str(parent)] = sha(parent)
    old = load(str(parent))['arrays']['solid']
    ids = np.asarray(evaluation['particle_id'], int)
    if not np.array_equal(ids, old.particle_id):
        raise ValueError('parent ID ordering mismatch')
    residual = np.asarray(evaluation['residual']['mean_velocity_m_s'])
    rnorm = np.linalg.norm(residual, axis=1); squared = rnorm*rnorm
    records = evaluation['boundary']['traces']; dt = evaluation['dt_s']
    rows = []
    for at, record in enumerate(records):
        if record['particle_id'] != int(ids[at]):
            raise ValueError('boundary trace ordering mismatch')
        v0 = np.asarray(record['material_velocity_before_m_s'])
        v1 = np.asarray(record['result']['velocity'])
        mean = np.asarray(evaluation['map']['mean_velocity_m_s'])[at]
        events = record['result']['events']
        faces = record['contact_rows']
        rows.append(dict(particle_id=int(ids[at]),
            mean_velocity_residual_m_s=residual[at].tolist(), mean_residual_norm_m_s=float(rnorm[at]),
            contact_features=[a['feature_id'] for a in faces],
            normal_force_N=[a['normal_force_N'] for a in faces],
            normal_stiffness_dimensionless_per_face=dt*dt*p['normal_stiffness']/(4*record['mass_kg']),
            actuator_active=record['actuator_active'],
            events=events, event_count=len(events),
            latest_zero_event_fraction=(max(a['time_s'] for a in events)/dt if events else None),
            mean_minus_endpoint_average_m_s=(mean-.5*(v0+v1)).tolist(),
            strain_rate_per_s=evaluation['material']['strain_rate_per_s'][at],
            plastic_work_J=evaluation['material']['plastic_work_J'][at]))
    def share(mask):
        return float(squared[np.asarray(mask, bool)].sum()/squared.sum()) if squared.sum() else 0.
    # Independently reconstruct only the final corotated elastic trial/yield
    # margin from archived gradient and original stress. No updater invocation.
    l = np.asarray(evaluation['midpoint']['velocity_gradient_per_s'])
    d = .5*(l+l.swapaxes(1, 2)); dev = d-np.trace(d, axis1=1, axis2=2)[:, None, None]*np.eye(3)/3
    spin = .5*(l-l.swapaxes(1, 2)); w = .5*dt*spin
    angle = np.sqrt(np.sum(w*w, axis=(1, 2))/2)
    # Rodrigues coefficients use removable limits in independent scalar form.
    a = np.ones(len(ids)); b = np.full(len(ids), .5); nonzero = angle != 0
    a[nonzero] = np.sin(angle[nonzero])/angle[nonzero]
    b[nonzero] = 2*np.sin(angle[nonzero]/2)**2/angle[nonzero]**2
    half = np.eye(3)+a[:, None, None]*w+b[:, None, None]*(w@w); full = half@half
    s0 = old.S.reshape(-1, 3, 3)
    stress_trial = full@s0@full.swapaxes(1, 2)+2*p['shear']*dt*(half@dev@half.swapaxes(1, 2))
    q_trial = np.sqrt(1.5*np.sum(stress_trial*stress_trial, axis=(1, 2)))
    rate = np.sqrt((2/3)*np.sum(dev*dev, axis=(1, 2)))
    temperature = np.asarray(evaluation['midpoint']['temperature_K'])
    ts = (temperature-p['room'])/(p['melt']-p['room'])
    # rate-floor is reported unchanged; not added to or removed from the map.
    lam = 1+p['b1']*ts**p['b2']*np.log(np.maximum(rate, p['rate_floor'])/p['reference_rate'])
    theta = 1-(1+np.exp(-p['c1']*ts))**(-1/p['c2'])
    yield0 = (p['a1']+p['a2']*np.arctan(p['a3']*old.ep))*lam*theta
    margin = yield0-q_trial
    comparison = np.max(np.abs(rate-np.asarray(evaluation['material']['strain_rate_per_s'])))
    report = dict(classification='READ_ONLY_FINAL_RAW_ROOT_DIAGNOSIS_NOT_PHYSICAL_ADVANCE',
        status='EVIDENCE_RECORDED_NOT_CAUSAL_JACOBIAN_PROOF',
        source_run_status=manifest['status'], accepted_steps=manifest['accepted_steps'],
        accepted_time_s=manifest['accepted_time_s'], attempted_dt_s=dt,
        evaluations=len(trace), first=trace[0], last=trace[-1],
        best_raw_scaled_residual=min(a['maximum_scaled_residual'] for a in trace if 'maximum_scaled_residual' in a),
        last_five_raw_scaled_residuals=[a.get('maximum_scaled_residual') for a in trace[-5:]],
        residual_squared_share=dict(no_contact=share([not a['contact_features'] for a in rows]),
            multiple_contacts=share([len(a['contact_features'])>1 for a in rows]),
            zero_event_particles=share([a['event_count']>0 for a in rows])),
        largest_residual_particles=sorted(rows, key=lambda a: a['mean_residual_norm_m_s'], reverse=True)[:16],
        final_material=dict(plastic_work_positive_count=int(np.count_nonzero(np.asarray(evaluation['material']['plastic_work_J'])>0)),
            active_trial_yield_count_independent=int(np.count_nonzero(margin<0)),
            minimum_trial_yield_margin_Pa=float(margin.min()),
            minimum_margin_particle_id=int(ids[np.argmin(margin)]),
            independent_rate_difference_per_s=float(comparison),
            closest_to_rate_floor_particle_id=int(ids[np.argmin(abs(rate-p['rate_floor']))]),
            closest_to_rate_floor_distance_per_s=float(np.min(abs(rate-p['rate_floor']))),
            interpretation='FINAL_CANDIDATE_ONLY_NOT_ALL_EARLIER_BRANCHES'),
        feature_guards=evaluation['feature_guards'], spatial=evaluation['spatial'],
        boundary_ledger=evaluation['boundary']['ledger'],
        preconditioner_backward_error_inf=evaluation.get('numerical_proposal_only', {}).get('diagnostics', {}).get('backward_error_inf'),
        unresolved='No finite-difference full-trial or contact Jacobian was run. Spatial/material/contact tangent attribution is evidence-supported, not isolated proof.',
        input_sha256={str(Path(k).relative_to(ROOT)):v for k,v in watched.items()},
        analysis_source_sha256=sha(Path(__file__)))
    for path, expected in watched.items():
        if sha(Path(path)) != expected:
            raise ValueError('input changed during read-only diagnosis '+path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('status','source_run_status','evaluations','best_raw_scaled_residual',
                                          'last_five_raw_scaled_residuals','residual_squared_share','final_material')}))


if __name__ == '__main__':
    main()

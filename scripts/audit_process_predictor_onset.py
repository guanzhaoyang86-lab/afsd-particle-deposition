"""Independent saved-array audit only: standard library + NumPy, no solver."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def same_array(a, b):
    return (a.dtype == b.dtype and a.shape == b.shape
            and a.tobytes(order='C') == b.tobytes(order='C'))


def verify_artifacts(path, report):
    for name, expected in report['artifact_sha256'].items():
        require(digest(path/name) == expected, 'artifact hash mismatch: '+name)
    for name, expected in report['source_sha256'].items():
        require(digest(ROOT/name) == expected, 'baseline source changed: '+name)


def audit(run, baseline):
    report = read_json(run/'report.json')
    old_report = read_json(baseline/'report.json')
    verify_artifacts(run, report)
    verify_artifacts(baseline, old_report)
    provenance = read_json(run/'predictor_provenance.json')
    supplemental_hashes = 0
    for key, inventory in provenance.items():
        if key == 'source_sha256' or key.endswith('_source_sha256'):
            for name, expected in inventory.items():
                require(digest(ROOT/name) == expected, 'provenance source changed: '+name)
                supplemental_hashes += 1

    with np.load(run/'history.npz', allow_pickle=False) as raw:
        history = {name: raw[name] for name in raw.files}
    with np.load(baseline/'history.npz', allow_pickle=False) as raw:
        old = {name: raw[name] for name in raw.files}
    require(set(history) == set(old), 'history field inventory changed')
    n_old = len(old['physical_time_s'])
    n_new = len(history['physical_time_s'])
    static = {'particle_id', 'tracer_type', 'reference_position_m',
              'reference_volume_m3', 'mass_kg', 'support_components',
              'contact_particle_id', 'contact_feature_id'}
    comparisons = []
    for name in sorted(history):
        after = history[name] if name in static else history[name][:n_old]
        equal = same_array(old[name], after)
        comparisons.append(dict(field=name, kind='STATIC' if name in static else 'TIME_PREFIX',
                                identical_dtype_shape_bytes=equal))
        require(equal, 'physical history prefix differs: '+name)
        require(np.isfinite(history[name]).all(), 'nonfinite history: '+name)

    checkpoint_comparison = []
    record_diff_keys = set()
    for step in range(n_old):
        name = f'accepted_{step:04d}.npz'
        with np.load(run/name, allow_pickle=False) as new, np.load(baseline/name, allow_pickle=False) as previous:
            require(set(new.files) == set(previous.files), 'accepted inventory differs')
            different = [key for key in new.files if not same_array(new[key], previous[key])]
            require(set(different) <= {'record_json'}, 'accepted physical arrays differ')
            if 'record_json' in different:
                newer, older = json.loads(str(new['record_json'])), json.loads(str(previous['record_json']))
                keys = [key for key in newer.keys() | older.keys() if newer.get(key) != older.get(key)]
                require(set(keys) <= {'iteration_trace'}, 'accepted scalar physics record differs')
                record_diff_keys.update(keys)
            checkpoint_comparison.append(dict(step=step, differing_fields=different))

    particle_id = history['particle_id']
    contact_id = history['contact_particle_id']
    require(len(np.unique(particle_id)) == len(particle_id), 'duplicate particle ID')
    require(len(np.unique(contact_id)) == len(contact_id), 'duplicate contact ID')
    positions = {int(pid): row for row, pid in enumerate(particle_id)}
    indices = np.asarray([positions[int(pid)] for pid in contact_id])
    for step in range(n_new):
        with np.load(run/f'accepted_{step:04d}.npz', allow_pickle=False) as state:
            require(same_array(state['particle_id'], particle_id), 'accepted ID inventory changed')
            require(same_array(state['tracer_type'], history['tracer_type']), 'accepted tags changed')
            require(state['restart_allowed'].item() is False, 'debug accepted restart enabled')
    require(np.array_equal(history['contact_feature_id'], np.ones(len(indices), dtype=int)),
            'unexpected feature inventory')

    records = report['records']
    require(len(records) == n_new, 'record/frame mismatch')
    wall = np.asarray([r['wall_displacement_m'] for r in records])
    du = np.diff(history['displacement_m'], axis=0)
    dw = np.diff(wall, axis=0)
    tau = history['tangent_force_global_yz_N']
    midtau = .5*(tau[1:] + tau[:-1])
    slip = du[:, indices, 1:] - dw[:, None, 1:]
    q = -np.sum(midtau*slip, axis=2)
    force = history['contact_force_N']
    midforce = .5*(force[1:] + force[:-1])
    ww = np.einsum('si,si->s', midforce.sum(axis=1), dw)
    wm = np.einsum('spi,spi->s', midforce, du)
    supports = history['support_force_N']
    wr = np.einsum('spi,spi->s', .5*(supports[1:] + supports[:-1]), du)
    q_saved = history['interval_friction_J'][1:]
    q_error = float(np.max(np.abs(q-q_saved)))
    require(q_error <= 1e-14, 'signed midpoint Q reconstruction mismatch')
    contact_identity = ww-wm-q.sum(axis=1)
    cumulative_q = np.cumsum(q, axis=0)
    require(np.max(np.abs(cumulative_q-history['friction_dissipation_J'][1:])) <= 1e-14,
            'cumulative signed friction mismatch')
    require(np.max(np.abs(contact_identity)) <= 1e-12, 'contact work identity failed')
    for name, values in [('interval_wall_work_J', ww),
                         ('interval_contact_material_work_J', wm)]:
        require(np.max(np.abs(values-np.asarray([r[name] for r in records[1:]]))) <= 1e-14,
                'saved work differs from independent array reconstruction')
    raw_negative = np.argwhere(q < 0)
    material_gap = np.array([r['stress_work_J']-r['elastic_diagnostic_J']-r['plastic_work_J']
                             for r in records])
    require(np.max(np.abs(material_gap-np.asarray([r['work_minus_elastic_minus_plastic_J']
                                                   for r in records]))) <= 1e-14,
            'material warning arithmetic mismatch')
    stress_contact = np.cumsum(wm+wr)-np.asarray([r['stress_work_J'] for r in records[1:]])

    attempts = [read_json(p) for p in sorted(run.glob('attempt_*.json'))]
    roots = []
    predictor_roots = []
    accepted_trace = []
    correction_max = 0
    predictor_rows = 0
    for attempt in attempts:
        require(attempt['error'] is None, 'unexpected failed attempt')
        for row in attempt['trace']:
            if row['stage'] == 'ACCEPTED':
                accepted_trace.append(row)
            if row.get('status') != 'ROOT_RETURNED_NOT_ACCEPTED':
                continue
            iterations = row['iterations']
            used = sum('line_search_factor' in item for item in iterations)
            require(used == iterations[-1]['iteration'] <= 14, 'shared correction budget differs')
            correction_max = max(correction_max, used)
            predictors = [item for item in iterations if item['stage']=='INTERIOR_PREDICTOR']
            predictor_rows += len(predictors)
            roots.append(row)
            if predictors:
                predictor_roots.append(dict(requested_time_s=attempt['requested_time_s'],
                    base_step=attempt['base_step'], branch_stage=row['stage'],
                    branch_time_s=row['time_s'], predictor_corrections=len(predictors),
                    total_shared_corrections=used,
                    initial_minimum_normal_N=min(predictors[0]['normal_force_N']),
                    final_minimum_normal_N=min(iterations[-1]['normal_force_N'])))
    require(len(accepted_trace) == report['accepted_steps'], 'accepted trace count differs')
    require([r['step'] for r in accepted_trace] == list(range(1, n_new)), 'accepted chronology differs')
    ft = 9e-5  # Frozen actual branch force tolerance; not fitted to this output.
    new_events = []
    for event in report['events']:
        if event['step'] < n_old:
            continue
        released = np.asarray(event['released_rows'], dtype=bool)
        margin = np.asarray(event['cone_margin_N'])
        candidates = (~released) & (np.abs(margin) <= ft)
        require(candidates.any(), 'event has no newly eligible contact within frozen tolerance')
        new_events.append(dict(step=event['step'], time_s=event['physical_time_s'],
            unreleased_near_boundary_particle_id=contact_id[candidates].tolist(),
            feature_id=history['contact_feature_id'][candidates].tolist(),
            already_released_particle_id=contact_id[released].tolist(),
            bracket_width_s=event['bracket_s'][1]-event['bracket_s'][0]))
    final = records[-1]
    return dict(status='PASS_INDEPENDENT_SAVED_ARRAY_AUDIT_NOT_AFSD_VALIDATION',
        imports='standard library and NumPy only; no solver/root/material imports or executions',
        artifact_hash_counts=dict(new=len(report['artifact_sha256']), old=len(old_report['artifact_sha256']),
                                  provenance_source_entries=supplemental_hashes),
        history_prefix=dict(old_frames=n_old, new_frames=n_new, all_fields_byte_identical=True,
                            fields=comparisons),
        accepted_prefix=dict(states=n_old, physical_arrays_byte_identical=True,
            metadata_only_difference_keys=sorted(record_diff_keys), states_detail=checkpoint_comparison,
            note='record_json iteration_trace changed schema; no scalar physical record differs'),
        path=dict(accepted_steps=report['accepted_steps'], accepted_time_s=report['accepted_time_s'],
                  total_events=len(report['events']), requests=len(attempts), returned_roots=len(roots),
                  predictor_branches=len(predictor_roots), predictor_corrections=predictor_rows,
                  maximum_total_shared_corrections=correction_max, predictor_details=predictor_roots,
                  events_after_old_accepted=new_events),
        independent_work=dict(raw_Q_min_J=float(q.min()), negative_raw_components=len(raw_negative),
            negative_entries=[dict(step=int(i)+1, particle_id=int(contact_id[j]), raw_Q_J=float(q[i,j]))
                              for i,j in raw_negative],
            raw_Q_reconstruction_max_error_J=q_error,
            friction_total_J=float(q.sum()), wall_work_J=float(ww.sum()), material_contact_work_J=float(wm.sum()),
            support_work_J=float(wr.sum()),
            maximum_interval_contact_work_identity_J=float(np.max(np.abs(contact_identity))),
            maximum_stress_contact_ledger_difference_J=float(np.max(np.abs(stress_contact))),
            final_stress_work_J=final['stress_work_J'], final_elastic_diagnostic_J=final['elastic_diagnostic_J'],
            final_plastic_work_J=final['plastic_work_J'], final_material_warning_J=float(material_gap[-1]),
            maximum_abs_material_warning_J=float(np.max(np.abs(material_gap))),
            material_warning_steps=report['material_work_warning_steps'],
            note='stress_work and elastic diagnostic read from recorded scalars, NOT independently recomputed material'),
        inventory=dict(all_history_fields_finite=True, particle_count=len(particle_id),
            contact_particle_id=contact_id.tolist(), particle_ids_stable_all_accepted=True,
            tags_stable_all_accepted=True, mass_kg=float(history['mass_kg'].sum()),
            minimum_accepted_normal_excluding_initial_N=float(history['normal_force_N'][1:].min()),
            maximum_ep=float(history['equivalent_plastic_strain'].max()),
            maximum_contact_law_residual_N=float(np.linalg.norm(history['contact_law_residual_global_yz_N'],axis=2).max())),
        not_validated=['material energy closure', 'normal opening/recontact', 'thermal',
                       'FKS', 'AFSD deposition', '8-increment sensitivity', 'spatial convergence'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.run, args.baseline), indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

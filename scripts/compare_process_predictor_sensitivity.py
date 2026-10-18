"""Read-only same-problem 4/8 comparison; NumPy only, no solver imports.

Sensitivity flags are engineering screens, NOT physical/continuum certificates.
See PROCESS_PREDICTOR_SENSITIVITY_CONTRACT.md for preregistered definitions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = 'docs/PROCESS_PREDICTOR_SENSITIVITY_CONTRACT.md'
TIME_TOL = 1e-12
STATIC = ('particle_id', 'tracer_type', 'reference_position_m',
          'reference_volume_m3', 'mass_kg', 'support_components',
          'contact_particle_id', 'contact_feature_id')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def same(a, b):
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def match_times(times, requested):
    times = np.asarray(times)
    require(times.ndim == 1 and np.isfinite(times).all() and
            np.all(np.diff(times) > 0), 'invalid accepted time sequence')
    indices = []
    for t in requested:
        found = np.flatnonzero(np.abs(times-t) <= TIME_TOL)
        require(len(found) == 1, f'missing or ambiguous common time {t}')
        indices.append(int(found[0]))
    require(len(set(indices)) == len(indices), 'reused accepted frame')
    return np.asarray(indices)


def field_norm(a, rank):
    a = np.asarray(a, dtype=float)
    require(np.isfinite(a).all(), 'nonfinite comparison field')
    if rank:
        return np.sqrt(np.sum(a*a, axis=tuple(range(a.ndim-rank, a.ndim))))
    return np.abs(a)


def metric(a, b, nominal, rank=0):
    a, b = np.asarray(a), np.asarray(b)
    require(a.shape == b.shape and a.ndim >= 1, 'comparison shape mismatch')
    require(nominal > 0 and np.isfinite(nominal), 'invalid nominal scale')
    an, bn, dn = field_norm(a, rank), field_norm(b, rank), field_norm(b-a, rank)
    denominator = max(float(an.max()), float(bn.max()), .01*nominal)
    per_time = dn.reshape(len(dn), -1).max(axis=1)
    relative = float(per_time.max()/denominator)
    return dict(status='WITHIN_SCREEN' if relative <= .05 else 'SENSITIVE_WARNING',
                nominal_scale=nominal, denominator=denominator,
                maximum_absolute_difference=float(per_time.max()),
                maximum_normalized_difference=relative,
                final_absolute_difference=float(per_time[-1]),
                final_normalized_difference=float(per_time[-1]/denominator),
                per_time_absolute_difference=per_time.tolist(),
                per_time_normalized_difference=(per_time/denominator).tolist(),
                coarse_final=a[-1].tolist(), fine_final=b[-1].tolist())


def event_inventory(report, h):
    inventory, groups, sequence, counts = {}, [], [], {}
    ids, features = h['contact_particle_id'], h['contact_feature_id']
    for event in report['events']:
        t = float(event['physical_time_s'])
        require(np.isfinite(t) and 0 < t <= .2+TIME_TOL, 'invalid event time')
        phase = int(np.searchsorted(np.array([.05, .1, .15, .2]), t-TIME_TOL))+1
        margin = np.asarray(event['cone_margin_N'])
        released = np.asarray(event['released_rows'], dtype=bool)
        require(margin.shape == released.shape == ids.shape, 'event contact inventory mismatch')
        require(np.isfinite(margin).all(), 'nonfinite event margin')
        rows = np.flatnonzero((~released) & (np.abs(margin) <= 9e-5))
        require(len(rows) > 0, 'event without newly eligible contact')
        keys = []
        for row in rows:
            base = (phase, int(ids[row]), int(features[row]))
            counts[base] = counts.get(base, 0)+1
            key = ':'.join(map(str, (*base, counts[base])))
            inventory[key] = t
            keys.append(key)
        groups.append(dict(time_s=t, keys=keys))
        sequence.extend(sorted(keys))
    return dict(times_by_key=inventory, groups=groups, sequence=sequence)


def compare_events(a, b):
    ak, bk = set(a['times_by_key']), set(b['times_by_key'])
    pairs = [dict(key=k, coarse_time_s=a['times_by_key'][k],
                  fine_time_s=b['times_by_key'][k],
                  signed_difference_s=b['times_by_key'][k]-a['times_by_key'][k])
             for k in sorted(ak & bk)]
    shift = max((abs(p['signed_difference_s']) for p in pairs), default=0.)
    status = ('EVENT_INVENTORY_MISMATCH' if ak != bk else
              'WITHIN_SCREEN' if shift <= .00125 else 'SENSITIVE_WARNING')
    return dict(status=status, coarse_only=sorted(ak-bk), fine_only=sorted(bk-ak),
                maximum_paired_absolute_shift_s=shift, paired=pairs,
                sequence_equal=a['sequence'] == b['sequence'],
                coarse_groups=a['groups'], fine_groups=b['groups'],
                pairing_note='phase:particle_id:feature_id:occurrence; never raw event index',
                scope='Explicit CONE_LIMIT newly eligible boundary contacts, not all slip-mode changes')


def work_audit(h, records):
    lookup = {int(pid): i for i, pid in enumerate(h['particle_id'])}
    ci = np.asarray([lookup[int(pid)] for pid in h['contact_particle_id']])
    wall = np.asarray([r['wall_displacement_m'] for r in records])
    du, dw = np.diff(h['displacement_m'], axis=0), np.diff(wall, axis=0)
    tau = h['tangent_force_global_yz_N']
    q = -np.sum(.5*(tau[1:]+tau[:-1])*(du[:, ci, 1:]-dw[:, None, 1:]), axis=2)
    force = h['contact_force_N']
    midpoint = .5*(force[1:]+force[:-1])
    ww = np.einsum('ti,ti->t', midpoint.sum(axis=1), dw)
    wm = np.einsum('tpi,tpi->t', midpoint, du)
    residual = ww-wm-q.sum(axis=1)
    q_error = float(np.max(np.abs(q-h['interval_friction_J'][1:])))
    require(q_error <= 1e-14, 'signed Q reconstruction mismatch')
    require(np.max(np.abs(np.cumsum(q, axis=0)-h['friction_dissipation_J'][1:])) <= 1e-14,
            'cumulative signed Q mismatch')
    require(np.max(np.abs(residual)) <= 1e-12, 'contact work identity failure')
    require(float(q.min()) >= -1e-10, 'negative signed Q exceeds frozen gate')
    for key, values in [('interval_wall_work_J', ww), ('interval_contact_material_work_J', wm)]:
        require(np.max(np.abs(values-np.asarray([r[key] for r in records[1:]]))) <= 1e-14,
                'saved work reconstruction mismatch '+key)
    return dict(raw_minimum_Q_J=float(q.min()), raw_negative_components=int(np.sum(q < 0)),
                reconstructed_Q_J=float(q.sum()), reconstructed_Wwall_J=float(ww.sum()),
                reconstructed_Wmaterial_J=float(wm.sum()), maximum_Q_reconstruction_error_J=q_error,
                maximum_interval_contact_identity_J=float(np.max(np.abs(residual))))


def load_case(path, reader_path, expected_increments):
    r, cfg, rb = read(path/'report.json'), read(path/'resolved_config.json'), read(reader_path)
    require(r['status'] == 'J2_CONE_PATH_COMPLETE' and r['error'] is None,
            'BLOCKED_NOT_COMPARABLE: incomplete path')
    require(abs(r['accepted_time_s']-.2) <= TIME_TOL and r['requested_end_s'] == .2,
            'BLOCKED_NOT_COMPARABLE: unequal endpoint')
    require(r['increments_per_phase'] == cfg['increments_per_phase'] == expected_increments,
            'wrong increment grid')
    require(rb['status'] == 'PASS_COMPLETE_READBACK' and
            rb['completion_status'] == r['status'] and rb['accepted_time_s'] == r['accepted_time_s'] and
            rb['accepted_steps'] == r['accepted_steps'], 'missing matching complete readback')
    require(rb['material_energy_certified'] is False and r['production_authorized'] is False,
            'unauthorized certification')
    require(rb['final'] == r['records'][-1] and rb['events'] == r['events'] and
            rb['material_work_warning_steps'] == r['material_work_warning_steps'],
            'readback evidence belongs to a different run')
    for name, digest in r['artifact_sha256'].items():
        require(Path(name).name == name and sha(path/name) == digest, 'artifact changed: '+name)
    provenance = read(path/'predictor_provenance.json')
    audit = rb['predictor_onset_audit']
    require(audit['variant'] == provenance['variant'] and
            audit['onset_source_sha256'] == provenance['onset_source_sha256'] and
            audit['supplemental_source_sha256'] == provenance['source_sha256'],
            'readback source inventory mismatch')
    sources = dict(r['source_sha256'])
    for key, value in provenance.items():
        if key == 'source_sha256' or key.endswith('_source_sha256'):
            for name, digest in value.items():
                require(name not in sources or sources[name] == digest, 'inconsistent source hash')
                sources[name] = digest
    require(len(sources) == 23, 'incomplete source inventory')
    for name, digest in sources.items():
        require(sha(ROOT/name) == digest, 'source changed: '+name)
    with np.load(path/'history.npz', allow_pickle=False) as raw:
        h = {k: raw[k] for k in raw.files}
    for k, value in h.items():
        require(np.isfinite(value).all(), 'nonfinite field '+k)
    require(len(np.unique(h['particle_id'])) == len(h['particle_id']), 'duplicate particle IDs')
    require(len(np.unique(h['contact_particle_id'])) == len(h['contact_particle_id']), 'duplicate contact IDs')
    require(len(r['records']) == len(h['physical_time_s']) == r['accepted_steps']+1, 'frame count mismatch')
    require(np.array_equal(h['physical_time_s'], [v['physical_time_s'] for v in r['records']]),
            'record time mismatch')
    require(h['physical_time_s'][0] == 0. and abs(h['physical_time_s'][-1]-.2) <= TIME_TOL,
            'history endpoints not complete')
    return dict(path=path, r=r, cfg=cfg, rb=rb, h=h, sources=sources,
                readback_sha256=sha(reader_path), work=work_audit(h, r['records']))


def compare_loaded(coarse, fine):
    a, b = coarse, fine
    ca, cb = dict(a['cfg']), dict(b['cfg'])
    ca.pop('increments_per_phase'); cb.pop('increments_per_phase')
    require(ca == cb, 'physical/configuration mismatch')
    require(a['sources'] == b['sources'], 'source inventory mismatch')
    for name in STATIC:
        require(same(a['h'][name], b['h'][name]), 'static geometry/ID/boundary mismatch: '+name)
    requested = np.linspace(0., .2, 17)
    indices = [match_times(c['h']['physical_time_s'], requested) for c in (a, b)]
    data = []
    for c, ix in zip((a, b), indices):
        h, records = c['h'], [c['r']['records'][i] for i in ix]
        require(all(r['event_kind'] == 'REQUESTED' for r in records[1:]),
                'common time is not a requested accepted frame')
        wall = np.asarray([r['wall_displacement_m'] for r in records])
        force = h['contact_force_N'][ix].sum(axis=1)
        require(np.max(np.abs(force-np.asarray([r['wall_force_on_material_N'] for r in records]))) <= 1e-10,
                'wall resultant differs from saved contact sum')
        weights = h['reference_volume_m3']/h['reference_volume_m3'].sum()
        data.append(dict(wall=wall, displacement_m=h['displacement_m'][ix],
                         contact_resultant_N=force, contact_tangent_resultant_N=force[:, 1:],
                         stress_Pa=h['stress_Pa'][ix], ep=h['equivalent_plastic_strain'][ix],
                         volume_mean_ep=h['equivalent_plastic_strain'][ix]@weights,
                         plastic_work_J=h['plastic_work_J'][ix].sum(axis=1),
                         friction_work_J=h['friction_dissipation_J'][ix].sum(axis=1),
                         wall_work_J=np.asarray([r['wall_work_J'] for r in records]),
                         stress_work_J=np.asarray([r['stress_work_J'] for r in records])))
    require(np.max(np.abs(data[0]['wall']-data[1]['wall'])) <= 1e-14, 'wall program changed')
    definitions = dict(displacement_m=(40e-6, 1), contact_resultant_N=(9000., 1),
                       contact_tangent_resultant_N=(900., 1), stress_Pa=(90e6, 2),
                       ep=(.004, 0), volume_mean_ep=(.004, 0), plastic_work_J=(.36, 0),
                       friction_work_J=(.36, 0), wall_work_J=(.36, 0), stress_work_J=(.36, 0))
    metrics = {name: metric(data[0][name], data[1][name], *args) for name, args in definitions.items()}
    events = compare_events(event_inventory(a['r'], a['h']), event_inventory(b['r'], b['h']))
    warnings = [key for key, m in metrics.items() if m['status'] != 'WITHIN_SCREEN']
    if events['status'] != 'WITHIN_SCREEN':
        warnings.append('events')
    summaries = []
    for c in (a, b):
        r, h, last = c['r'], c['h'], c['r']['records'][-1]
        summaries.append(dict(increments_per_phase=r['increments_per_phase'],
            accepted_steps=r['accepted_steps'], events=len(r['events']), accepted_time_s=r['accepted_time_s'],
            advance_wall_s=r['advance_wall_s'], selected_array_bytes=r['selected_array_bytes'],
            mass_kg=float(h['mass_kg'].sum()), maximum_final_ep=float(h['equivalent_plastic_strain'][-1].max()),
            material_warning_steps=r['material_work_warning_steps'],
            signed_material_warning_J=last['stress_work_J']-last['elastic_diagnostic_J']-last['plastic_work_J'],
            material_energy_certified=False, independent_work=c['work']))
    return dict(status='COMPARISON_COMPLETE', screening_status='SENSITIVE_WARNING' if warnings else 'WITHIN_SCREEN',
        warnings=warnings, classification='ENGINEERING_APPROXIMATE_SENSITIVITY_NOT_AFSD_VALIDATION',
        production_authorized=False, common_time_s=requested.tolist(),
        matched_indices=[v.tolist() for v in indices], metrics=metrics, events=events, cases=summaries,
        caveats=['Two grids do not establish convergence order or continuum accuracy.',
                 'No interpolation; same accepted endpoint and unchanged physics.',
                 'Material energy warning is NOT corrected or converted to heat.',
                 'Spatial/thermal/rotating-tool/AFSD/production NOT_RUN.'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coarse', type=Path, required=True)
    p.add_argument('--fine', type=Path, required=True)
    p.add_argument('--coarse-readback', type=Path, required=True)
    p.add_argument('--fine-readback', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    require(not args.output.exists(), 'never overwrite previous comparison')
    a = load_case(args.coarse, args.coarse_readback, 4)
    b = load_case(args.fine, args.fine_readback, 8)
    result = compare_loaded(a, b)
    result['provenance'] = dict(contract_sha256=sha(ROOT/CONTRACT), comparator_sha256=sha(__file__),
        sources=a['sources'], inputs=[dict(path=str(c['path']), report_sha256=sha(c['path']/'report.json'),
                                         readback_sha256=c['readback_sha256']) for c in (a, b)])
    with args.output.open('x', encoding='utf-8') as f:
        f.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0  # warnings are reported, not an invitation to tune physics


if __name__ == '__main__':
    raise SystemExit(main())

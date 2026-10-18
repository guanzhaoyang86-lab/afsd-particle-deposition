"""Read-only common-endpoint E2 sensitivity analysis; no solver changes."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def norm(a, v):
    return float(np.sqrt(np.sum(v*np.sum(a.reshape(len(v), -1)**2, axis=1))))


def summarize(root):
    names = ['n5_cfl010', 'n5_cfl005', 'n5_cfl0025', 'n7_cfl005', 'n9_cfl005']
    data, histories = {}, {}
    for name in names:
        directory = root/name
        r = json.loads((directory/'report.json').read_text(encoding='utf-8'))
        if r['actual_steps'] != r['steps_to_common_end'] or abs(r['actual_end_s']-r['common_end_s']) > 1e-18:
            raise ValueError('incomplete common end: '+name)
        for file, expected in r['artifact_sha256'].items():
            if hashlib.sha256((directory/file).read_bytes()).hexdigest() != expected:
                raise ValueError('artifact hash mismatch: '+name+'/'+file)
        with np.load(directory/'history.npz', allow_pickle=False) as h:
            histories[name] = {key: h[key].copy() for key in h.files}
        data[name] = r
    first = data[names[0]]
    for r in data.values():
        if r['source_sha256'] != first['source_sha256'] or r['common_end_s'] != first['common_end_s']:
            raise ValueError('different problem/code or terminal time')
    differences = []
    for a, b in zip(names[:2], names[1:3]):
        ha, hb = histories[a], histories[b]
        np.testing.assert_array_equal(ha['particle_id'], hb['particle_id'])
        np.testing.assert_array_equal(ha['reference_position_m'], hb['reference_position_m'])
        scale = norm(ha['reference_final_displacement_m'], ha['volume_m3'])
        differences.append(norm(ha['position_m'][-1]-hb['position_m'][-1], ha['volume_m3'])/scale)
    ratio = differences[1]/differences[0] if differences[0] > 1e-10 else None
    temporal = differences[-1] <= .01 and (max(differences) <= 1e-10 or ratio <= .4)
    spatial = [data[k]['metrics'] for k in ('n5_cfl005', 'n7_cfl005', 'n9_cfl005')]
    space_ratios = {key: [b[key]/a[key] for a, b in zip(spatial, spatial[1:])] for key in ('displacement_relative', 'stress_relative')}
    space_trend = all(x <= .9 for values in space_ratios.values() for x in values)
    # Finest-grid physical accuracy targets; coarse-screen failures remain visible.
    fine_checks = data['n9_cfl005']['checks']
    numerical = all(all(r['checks'][k] for k in ('mass_relative', 'energy_relative', 'force', 'torque', 'constraints', 'elastic_only', 'ids')) for r in data.values())
    rows = []
    for name, r in data.items():
        h = histories[name]
        q = np.pi/(4*(r['n']-1))
        rows.append(dict(name=name, status=r['status'], steps=r['actual_steps'], dt_s=r['dt_s'],
            setup_wall_s=r['setup_wall_s'], advance_wall_s=r['advance_wall_s'],
            state_array_bytes=r['state_array_bytes'], dense_gradient_bytes=r['dense_gradient_bytes'],
            metrics=r['metrics'], final_ledger=r['records'][-1],
            actual_body_resultant_N=float(h['body_force_N'][:, 0].sum()),
            independent_trapezoid_body_resultant_N=float(70*q/np.tan(q)),
            continuum_body_resultant_N=70.))
    ok = temporal and space_trend and all(fine_checks.values()) and numerical
    return dict(status='PASS_LIMITED_SMALL_STRAIN_BENCHMARK' if ok else 'FAIL_ENGINEERING_TARGETS',
        production_authorized=False, common_end_s=first['common_end_s'],
        temporal_status='PASS' if temporal else 'FAIL', temporal_displacement_differences=differences,
        temporal_difference_ratio=ratio, spatial_trend_status='PASS' if space_trend else 'FAIL',
        spatial_error_ratios=space_ratios, all_numerical_invariants=numerical,
        finest_grid_checks=fine_checks, rows=rows,
        not_validated=['finite_strain', 'plastic_contact', 'thermal', 'AA6061', 'AFSD', 'checkpoint_restart'],
        references='docs/ENGINEERING_BLOCK_E2.md preregistered linear continuum; not experimental data')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = summarize(args.root)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(0 if result['status'].startswith('PASS') else 1)

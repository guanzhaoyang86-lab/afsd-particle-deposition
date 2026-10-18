"""Read Q1 saved states/independent balance; never advances the benchmark."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from pysph.solver.utils import load
from afsd_scaffold.vtk_legacy import read_point_cloud


def inspect(directory):
    path = Path(directory)
    r = json.loads((path/'report.json').read_text())
    if r['status'] != 'DISCRETE_EQUILIBRIUM_COMPLETE' or r['accepted_time_s'] != .1:
        raise ValueError('incomplete run cannot enter common-end comparison')
    if r['production_authorized'] or r['restart_allowed']:
        raise ValueError('unsupported production/restart claim')
    for relative, digest in r['source_sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise ValueError('changed source: '+relative)
    for name, digest in r['artifact_sha256'].items():
        if hashlib.sha256((path/name).read_bytes()).hexdigest() != digest:
            raise ValueError('changed artifact: '+name)
    with np.load(path/'history.npz', allow_pickle=False) as f:
        h = {key: f[key].copy() for key in f.files}
    if any(not np.isfinite(value).all() for value in h.values()):
        raise ValueError('nonfinite history')
    ids = h['particle_id']; x = h['reference_position_m']; volume = h['reference_volume_m3']
    if len(np.unique(ids)) != len(ids) or not np.equal(ids, np.floor(ids)).all():
        raise ValueError('invalid persistent ID')
    if not np.allclose(h['mass_kg'], 2700*volume, rtol=1e-12, atol=0):
        raise ValueError('unscaled mass/quadrature mismatch')
    if abs(volume.sum()/1e-6-1) > 1e-12:
        raise ValueError('changed physical domain')
    t = h['physical_time_s']
    if not np.array_equal(t, [row['physical_time_s'] for row in r['records']]) or np.any(np.diff(t) <= 0):
        raise ValueError('history/report time mismatch')
    if len(t) != r['increments_requested']+1 or t[0] != 0 or t[-1] != .1:
        raise ValueError('wrong common endpoint or frame inventory')
    if not np.array_equal(h['constrained_components'], x == 0):
        raise ValueError('changed roller boundary')
    if np.any(h['displacement_m'][:, h['constrained_components']] != 0):
        raise ValueError('constrained component moved')
    if not np.array_equal(h['position_m'], x[None]+h['displacement_m']):
        raise ValueError('position/displacement mismatch')
    k = np.pi/.02
    force_balance = []; external_work = 0.
    q = np.pi/(4*(r['n']-1))
    for index, row in enumerate(r['records']):
        s = t[index]/.1; factor = 3*s*s-2*s*s*s
        expected = np.zeros_like(x)
        expected[:, 0] = volume*factor*70e9*1e-5*k*np.sin(k*x[:, 0])
        if not np.allclose(h['external_force_N'][index], expected, rtol=1e-12, atol=1e-15):
            raise ValueError('force differs from independent prescribed body load')
        if abs(expected[:, 0].sum()-factor*70*q/np.tan(q)) > 1e-11:
            raise ValueError('independent force quadrature identity failed')
        defect = float(np.linalg.norm((h['external_force_N'][index]+h['support_reaction_N'][index]).sum(axis=0))/70)
        if defect > 1e-7:
            raise ValueError('external/support resultant unbalanced')
        force_balance.append(defect)
        if index:
            du = h['displacement_m'][index]-h['displacement_m'][index-1]
            external_work += float(np.sum(.5*(h['external_force_N'][index]+h['external_force_N'][index-1])*du))
            if abs(external_work-row['external_work_J']) > 1e-12*1.75e-6:
                raise ValueError('saved external work not actual force-displacement integral')
        if abs(row['work_identity_residual_J']) > 1e-7*1.75e-6:
            raise ValueError('recorded internal/external work mismatch')
    data = read_point_cloud(path/'final.vtk')
    if not np.array_equal(data['points'], h['position_m'][-1]):
        raise ValueError('VTK does not match accepted positions')
    if not np.array_equal(data['scalars']['particle_id'], ids):
        raise ValueError('VTK ID mismatch')
    saved = load(str(path/'checkpoint.npz'))
    pa = saved['arrays']['solid']
    if saved['solver_data']['restart_allowed'] is not False:
        raise ValueError('diagnostic snapshot advertised as restart')
    if not np.array_equal(pa.particle_id, ids) or not np.array_equal(np.column_stack((pa.x, pa.y, pa.z)), h['position_m'][-1]):
        raise ValueError('checkpoint primary state mismatch')
    if any('temperature' in key or 'velocity' in key for key in data['scalars'] | data['vectors']):
        raise ValueError('unsolved thermal/dynamic field published')
    return r, h, dict(status='PASS_READBACK_ONLY', n=r['n'], increments=r['accepted_steps'],
        common_end_s=float(t[-1]), maximum_resultant_defect_relative=max(force_balance),
        independent_external_work_J=external_work, continuum_accuracy=r['continuum_accuracy'],
        compute_wall_s=r['wall_seconds'], geometry='fixed cube[0,.01]^3,three normal roller planes',
        production_authorized=False, not_afsd=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    rows = []; histories = {}
    for name in ('n3_i4', 'n5_i4', 'n5_i8', 'n7_i4'):
        _, data, row = inspect(args.directory/name)
        histories[name] = data; rows.append(row)
    coarse, fine = histories['n5_i4'], histories['n5_i8']
    V = coarse['reference_volume_m3']
    if not np.array_equal(coarse['reference_position_m'], fine['reference_position_m']):
        raise ValueError('time cases changed geometry')
    def relative_difference(field):
        a, b = coarse[field][-1], fine[field][-1]
        square = lambda value: np.sum(V*np.sum(value.reshape(len(V), -1)**2, axis=1))
        return float(np.sqrt(square(a-b)/square(b)))
    result = dict(status='PASS_READBACK', scope='Q1_QUASISTATIC_ONLY_NOT_AFSD', rows=rows,
        time_increment_comparison=dict(common_end_s=.1, coarse_increments=4, fine_increments=8,
            displacement_relative=relative_difference('displacement_m'),
            stress_relative=relative_difference('stress_Pa'),
            interpretation='same elastic equilibrium, NOT dynamic/FKS/contact time convergence'),
        continuum_precision_gate_applied=False, production_authorized=False)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

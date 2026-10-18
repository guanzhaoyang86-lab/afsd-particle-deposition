"""Independent post-run assertions; reads trusted locally generated checkpoints only."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.output_contract import validate_manifest
from afsd_scaffold.vtk_legacy import read_point_cloud
from afsd_backend.pysph.output import digest, restore


def validate(folder):
    folder = Path(folder)
    manifest = json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
    validate_manifest(manifest)
    assert digest(folder/manifest['configuration']['path']) == manifest['configuration']['sha256']
    assert digest(folder/manifest['particles']['file']) == manifest['particles']['sha256']
    assert digest(folder/'checkpoint.npz') == manifest['checkpoint']['sha256']
    series = json.loads((folder/'series.json').read_text(encoding='utf-8'))['frames']
    identity = None; times = []
    for frame in series:
        assert digest(folder/frame['file']) == frame['sha256']
        state = read_point_cloud(folder/frame['file']); s, v = state['scalars'], state['vectors']
        ids = s['particle_id']; assert len(set(ids)) == len(ids) and min(ids) > 0
        now = {pid: (v['initial_position_m'][i], s['material_id'][i], s['component_id'][i], s['tracer_type'][i], s['wholebar_tag'][i]) for i, pid in enumerate(ids)}
        if identity is not None: assert now == identity
        identity = now
        for values in list(s.values())+list(v.values())+[state['points']]: assert np.isfinite(values).all()
        for field in ('density_kg_per_m3', 'temperature_K', 'mass_kg', 'smoothing_length_m'): assert np.min(s[field]) > 0
        assert np.min(s['equivalent_plastic_strain']) >= 0
        pressure = np.array(s['pressure_Pa'])
        trace = sum(np.array(s['stress_'+axis+'_Pa']) for axis in ('xx', 'yy', 'zz'))
        np.testing.assert_allclose(trace, -3*pressure, atol=1e-7)
        q2 = sum((np.array(s['stress_'+axis+'_Pa'])+pressure)**2 for axis in ('xx', 'yy', 'zz'))
        q2 += 2*sum(np.array(s['stress_'+axis+'_Pa'])**2 for axis in ('xy', 'xz', 'yz'))
        np.testing.assert_allclose(np.sqrt(1.5*q2), s['equivalent_stress_Pa'], rtol=1e-10, atol=1e-7)
        times.append(frame['time_s'])
    assert all(b > a for a, b in zip(times, times[1:]))
    pa, restart = restore(folder/'checkpoint.npz')
    assert restart['count'] == series[-1]['step']
    assert set(pa.particle_id.tolist()) == set(identity)
    final = read_point_cloud(folder/series[-1]['file'])
    np.testing.assert_array_equal(pa.temperature, final['scalars']['temperature_K'])
    np.testing.assert_array_equal(np.column_stack((pa.x, pa.y, pa.z)), final['points'])
    print(f'PASS {folder.name}: status={manifest["status"]}, frames={len(series)}, particles={len(identity)}, last_step={restart["count"]}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('folders', nargs='+')
    for folder in parser.parse_args().folders: validate(folder)

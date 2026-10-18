"""Read two complete states; measure fixed cohorts without advancing physics."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import run_gpu_nonlinear_resume_20261011 as runner

def read_state(stage):
    manifest = runner.read(ROOT / stage / 'manifest.json')
    path = ROOT / manifest['screened_continuation_checkpoint']
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    assert sha == manifest['screened_continuation_sha256']
    pa, state = runner.load(path)
    return pa, state, {'checkpoint': str(path.relative_to(ROOT)), 'sha256': sha}

a, sa, source_a = read_state('outputs/gpu_acceleration_20261007/event_resume_8h_20261011_v1/attempt_0039')
b, sb, source_b = read_state('outputs/gpu_acceleration_20261007/nonlinear_resume_20261011_v1/attempt_0086')
ids = np.asarray(a.particle_id)
assert np.array_equal(ids, b.particle_id)
assert np.array_equal(a.m, b.m)
assert sa['bore_observation']['eligible_ids'] == sb['bore_observation']['eligible_ids']
assert sorted(sa['bore_observation']['crossing_ids']) == sorted(sb['bore_observation']['crossing_ids'])
xa = np.column_stack((a.x, a.y, a.z))
xb = np.column_stack((b.x, b.y, b.z))
displacement = xb - xa
radial_change = np.linalg.norm(xb[:, :2], axis=1) - np.linalg.norm(xa[:, :2], axis=1)
elapsed = sb['time_s'] - sa['time_s']
eligible = np.isin(ids, sa['bore_observation']['eligible_ids'])
crossed = np.isin(ids, sb['bore_observation']['crossing_ids'])
masks = {
    'actuator': np.asarray(a.actuator) != 0,
    'eligible_uncrossed_at_start_z_1p5_to_2mm': eligible & ~crossed & (a.z >= .0015) & (a.z <= .002),
    'previously_crossed': crossed,
}
groups = {}
for name, mask in masks.items():
    groups[name] = {
        'count': int(mask.sum()), 'particle_ids': ids[mask].astype(int).tolist(),
        'dz_mean_um': float(displacement[mask, 2].mean() * 1e6),
        'dz_min_um': float(displacement[mask, 2].min() * 1e6),
        'dz_max_um': float(displacement[mask, 2].max() * 1e6),
        'radial_change_mean_um': float(radial_change[mask].mean() * 1e6),
        'radial_change_each_um': (radial_change[mask] * 1e6).tolist(),
    }
commanded_dz = -sa['parameters']['feed'] * elapsed
error = float(np.max(np.abs(displacement[masks['actuator'], 2] - commanded_dz)))
assert error < 1e-12
result = {
    'classification': 'READ_ONLY_PAIRED_CHECKPOINT_KINEMATIC_OBSERVATION',
    'physical_simulation_steps_executed': 0,
    'sources': [source_a, source_b],
    'start_time_s': sa['time_s'], 'end_time_s': sb['time_s'], 'interval_s': elapsed,
    'groups': groups,
    'commanded_actuator_dz_um': commanded_dz * 1e6,
    'actuator_displacement_max_error_m': error,
    'limitations': [
        'Fixed cohorts between two saved states; displacement is not mass flux.',
        'Endpoint comparison does not establish a monotonic trajectory.',
        'Radial change uses the stationary tool axis at x=y=0.',
        'No contact-force cause or spatial convergence is established.',
    ],
}
Path(__file__).with_name('readback.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
print(json.dumps({k: v for k, v in result.items() if k not in ('sources', 'limitations')}, indent=2))

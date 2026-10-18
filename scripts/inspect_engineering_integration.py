"""Read-only real-coupling evidence; never accepts a REJECTED candidate."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud
from afsd_backend.pysph.output import restore


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def accepted_endpoint(series):
    frames = series['frames']
    if not frames or any(b['step'] <= a['step'] or b['time_s'] <= a['time_s'] for a, b in zip(frames, frames[1:])):
        raise ValueError('invalid accepted frame order')
    return frames[-1]


def inspect(directory):
    directory = Path(directory)
    series = json.loads((directory/'series.json').read_text(encoding='utf-8'))
    performance = json.loads((directory/'performance.json').read_text(encoding='utf-8'))
    manifest = json.loads((directory/'manifest.json').read_text(encoding='utf-8'))
    sources = json.loads((directory/'source_hashes.json').read_text(encoding='utf-8'))
    for file, digest in sources.items():
        if sha(ROOT/file) != digest:
            raise ValueError('changed source: '+file)
    for frame in series['frames']:
        if sha(directory/frame['file']) != frame['sha256']:
            raise ValueError('changed frame: '+frame['file'])
    if sha(directory/'checkpoint.npz') != manifest['checkpoint']['sha256']:
        raise ValueError('checkpoint hash mismatch')
    if sha(directory/'resolved_config.json') != manifest['configuration']['sha256']:
        raise ValueError('resolved config hash mismatch')
    last = accepted_endpoint(series)
    first_data = read_point_cloud(directory/series['frames'][0]['file'])
    final_data = read_point_cloud(directory/last['file'])
    first_ids = np.asarray(first_data['scalars']['particle_id'], dtype=np.int64)
    final_ids = np.asarray(final_data['scalars']['particle_id'], dtype=np.int64)
    first_order, final_order = np.argsort(first_ids), np.argsort(final_ids)
    np.testing.assert_array_equal(first_ids[first_order], final_ids[final_order])
    if len(np.unique(final_ids)) != len(final_ids):
        raise ValueError('duplicate persistent ID')
    for field in ('material_id', 'component_id', 'tracer_type', 'wholebar_tag', 'actuator'):
        np.testing.assert_array_equal(np.asarray(first_data['scalars'][field])[first_order],
                                      np.asarray(final_data['scalars'][field])[final_order])
    pa, state = restore(directory/'checkpoint.npz')
    if state['count'] != last['step'] or state['t'] != last['time_s']:
        raise ValueError('checkpoint and accepted series disagree')
    np.testing.assert_array_equal(pa.particle_id, final_ids)
    np.testing.assert_array_equal(np.column_stack((pa.x, pa.y, pa.z)), final_data['points'])
    displacement = np.asarray(final_data['points'])[final_order]-np.asarray(first_data['points'])[first_order]
    scalar = {k: np.asarray(v)[final_order] for k, v in final_data['scalars'].items()}
    velocity = np.asarray(final_data['vectors']['velocity_m_per_s'])[final_order]
    free = scalar['actuator'] == 0
    diagnostics = [json.loads(row) for row in (directory/'diagnostics.jsonl').read_text(encoding='utf-8').splitlines() if row.strip()]
    accepted_steps = {f['step'] for f in series['frames']}
    accepted_diagnostics = [d for d in diagnostics if d['step'] in accepted_steps]
    if any(d['alarms'] for d in accepted_diagnostics):
        raise ValueError('alarmed state published as accepted')
    if last['step'] and (not accepted_diagnostics or accepted_diagnostics[-1]['step'] != last['step']):
        raise ValueError('missing accepted endpoint diagnostics')
    rejected_path = directory/'REJECTED_candidate.npz'
    rejection = None
    if rejected_path.exists():
        try:
            restore(rejected_path)
        except ValueError:
            pass
        else:
            raise ValueError('rejected candidate can be restored')
        debug = json.loads((directory/'REJECTED_contact_debug.json').read_text(encoding='utf-8'))
        local = debug.get('failed_local_context') or {}
        cone = local.get('failed_cone') or {}
        rejection = dict(candidate_step=last['step']+1, particle_id=local.get('particle_id'),
            features=[x['feature'] for x in local.get('faces', [])],
            time_mid_s=local.get('time_mid_s'), full_step_dt_s=local.get('dt_s'),
            coefficient=cone.get('coefficient'), scaled_residual=cone.get('scaled_residual'),
            target=cone.get('target'), signed_gap=cone.get('signed_gap'),
            restart_allowed=False, reason=debug['error'].strip().splitlines()[-1])
    requested = performance['requested_end_time_s']
    reached = requested is not None and last['time_s'] >= requested and performance['status'] == 'COMPLETE'
    config = json.loads((directory/'resolved_config.json').read_text(encoding='utf-8'))
    process = config['parameters']['process']
    feed = process['feedstock_axial_speed']['value']; traverse = process['tool_traverse_speed']['value']
    omega = process['tool_angular_speed']['value']; t = last['time_s']
    return dict(status='IMPLEMENTATION_INTERVAL_COMPLETE' if reached else 'INCOMPLETE_COUPLED_IMPLEMENTATION',
        readback_status='PASS', source_and_artifact_integrity='PASS', production_authorized=False,
        paper_reproduction=False, deposition_validated=False, mixing_validated=False,
        runner_status=performance['status'], accepted_step=last['step'], accepted_frames=len(series['frames']),
        accepted_time_s=t, requested_end_time_s=requested, reached_requested_end=reached,
        particle_count=len(final_ids), last_accepted_diagnostics=accepted_diagnostics[-1] if accepted_diagnostics else None,
        observed=dict(free_particles=int(free.sum()), free_displacement_max_m=float(np.linalg.norm(displacement[free], axis=1).max()),
            free_speed_max_m_s=float(np.linalg.norm(velocity[free], axis=1).max()),
            maximum_equivalent_plastic_strain=float(scalar['equivalent_plastic_strain'].max()),
            plastified_particle_count=int(np.count_nonzero(scalar['equivalent_plastic_strain'] > 0)),
            substrate_displacement_max_m=float(np.linalg.norm(displacement[scalar['component_id'] == 0], axis=1).max()),
            temperature_min_K=float(scalar['temperature_K'].min()), temperature_max_K=float(scalar['temperature_K'].max()),
            tool_rotation_rad=omega*t, tool_travel_m=traverse*t, prescribed_cap_feed_m=feed*t),
        rejected=rejection, performance={k:performance[k] for k in ('wall_seconds','peak_process_working_set_bytes',
            'dt_min_s','dt_max_s','physics_step_median_seconds','physical_seconds_per_wall_second')},
        evidence={name:sha(directory/name) for name in ('performance.json','series.json','manifest.json',
            'checkpoint.npz','resolved_config.json','source_hashes.json')},
        interpretation='Real force-driven hot-start response only; no macroscopic deposition or plastic flow established')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = inspect(args.run)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(report, indent=2, allow_nan=False))
    raise SystemExit(0 if report['reached_requested_end'] else 1)

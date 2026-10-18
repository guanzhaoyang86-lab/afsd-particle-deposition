"""Package 23 real checkpoint samples after the independent refined audit.

Imports only unchanged readback helpers; never edits the old reader/artifacts,
evaluates a material or contact model, interpolates, or advances physical time.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import readback_discrete_process as reader

ROOT, BASE = reader.ROOT, reader.BASE
AUDIT = ROOT/'outputs/discrete_material_options_20261004/post_refinement_audit.json'


def package(output):
    begin = time.perf_counter()
    audit = reader.read(AUDIT)
    reader.require(audit['status'] == 'PASS_READ_ONLY_POST_REFINEMENT_AUDIT_NOT_TRANSPORT',
                   'the independent refined audit has not passed')
    reader.require(audit['accepted_steps_from_independent_initial'] == 22 and
                   audit['inheritance_and_all_22_prefixes_valid'], 'not the audited 22-step chain')
    for name, digest in audit['source_input_sha256'].items():
        reader.require(reader.sha(ROOT/name) == digest, 'audited input changed: '+name)
    for name, digest in audit['old_coarse_readback_unchanged_sha256'].items():
        reader.require(reader.sha(BASE/'readback'/name) == digest, 'old readback changed')
    reader.require(not output.exists(), 'refuse to overwrite refined readback')
    initial_file = BASE/'quarters/independent_initial.npz'
    initial, initial_state = reader.restore(initial_file)
    snapshots = [(initial, initial_state)]
    field_frames = [reader.fields(initial)]
    pose_rows, sample_sources = [], [dict(time_s=0., checkpoint=str(initial_file.relative_to(ROOT)),
                                        sha256=reader.sha(initial_file))]
    for directory in (BASE/'quarters', BASE/'post_start/halves', BASE/'post_start/pilot_refinement'):
        meta = reader.read(directory/'manifest.json')
        for name, digest in meta['source_sha256'].items():
            reader.require(reader.sha(ROOT/name) == digest, 'frozen source changed: '+name)
        history = reader.read(directory/'accepted_history.json')
        for row in history:
            path = directory/row['checkpoint']
            reader.require(reader.sha(path) == row['sha256'], 'accepted snapshot hash changed')
            pa, state = reader.restore(path)
            reader.require(state['step'] == snapshots[-1][1]['step']+1, 'wrong physical chain')
            reader.require(abs(state['time_s']-snapshots[-1][1]['time_s']-row['dt_s']) < 1e-15,
                           'physical sample time gap or overlap')
            data = reader.fields(pa)
            reader.require(np.array_equal(data['particle_id'], field_frames[0]['particle_id']),
                           'particle IDs changed')
            snapshots.append((pa, state)); field_frames.append(data)
            pose_rows.append(reader.read(path.parent/'tool_pose.json'))
            sample_sources.append(dict(time_s=state['time_s'], checkpoint=str(path.relative_to(ROOT)),
                                       sha256=reader.sha(path)))
    reader.require(len(snapshots) == 23, 'expected 23 actual samples')
    final_pa, final_state = snapshots[-1]
    reader.require(final_state['interval_energy'] == audit['whole_interval_energy'],
                   'final audited cumulative interval changed')
    p = final_state['parameters']
    poses = [reader.expected_pose(p, 0.)]+pose_rows
    times = np.array([state['time_s'] for _, state in snapshots])
    arrays = {f'particle_property__{key}': np.stack([frame[key] for frame in field_frames])
              for key in field_frames[0]}
    stress = arrays['particle_property__S'].reshape(23, len(initial.x), 3, 3)
    density = arrays['particle_property__rho'][..., 0]
    bulk = 2*p['shear']*(1+p['poisson'])/(3*(1-2*p['poisson']))
    pressure = bulk*(density/p['rho0']-1.)
    arrays.update(time_s=times, particle_id=field_frames[0]['particle_id'].ravel(),
        position_m=np.stack([reader.xyz(pa) for pa, _ in snapshots]),
        velocity_m_s=np.stack([np.column_stack((pa.u, pa.v, pa.w))[np.argsort(pa.particle_id)]
                              for pa, _ in snapshots]),
        temperature_K=arrays['particle_property__temperature'][..., 0],
        density_kg_m3=density, deviatoric_stress_Pa=stress,
        cauchy_stress_Pa=stress-pressure[..., None, None]*np.eye(3),
        plastic_strain=arrays['particle_property__ep'][..., 0],
        tool_origin_m=np.array([pose['origin'] for pose in poses]),
        tool_rotation=np.array([pose['rotation'] for pose in poses]),
        contact_history_m=np.stack([state['contact_history'][np.argsort(pa.particle_id)]
                                   for pa, state in snapshots]),
        crossing_mass_kg=np.array([state.get('crossing_mass_kg', 0.) for _, state in snapshots]),
        solved_rate_sample=np.array([False]+[True]*22),
        actual_checkpoint_samples_only=np.array(True), interpolation_performed=np.array(False))
    output.mkdir(parents=True)
    trajectory = output/'accepted_trajectory.npz'
    np.savez_compressed(trajectory, **arrays)
    # Verify the exported object, not merely successful file creation.
    with np.load(trajectory, allow_pickle=False) as saved:
        reader.require(set(saved.files) == set(arrays), 'exported array names changed')
        for name, value in arrays.items():
            reader.require(np.array_equal(saved[name], value), 'exported array changed: '+name)
    summary = dict(status='PASS_REFINED_ACTUAL_CHECKPOINT_PACKAGING_NOT_NEW_PHYSICAL_VALIDATION',
        classification=final_state['classification'], physical_evaluations=0,
        samples=23, accepted_steps=22, particle_count=len(initial.x), accepted_time_s=times[-1],
        branch_chain=['independent_initial', 'quarters_4', 'post_start_halves_2', 'pilot_refinement_16'],
        alternative_coarse_pilot_excluded=True, initial_zero_rates_are_not_computed_rates=True,
        final=audit['final'], temporal_screen=audit['independently_recomputed_temporal_comparison'],
        interval_energy=final_state['interval_energy'], cost=audit['cost'],
        audited_particle_properties=sorted(field_frames[0]), sample_sources=sample_sources,
        audit_file=str(AUDIT.relative_to(ROOT)), audit_sha256=reader.sha(AUDIT),
        packaging_script_sha256=reader.sha(Path(__file__)),
        unchanged_reader_script_sha256=reader.sha(Path(reader.__file__)),
        trajectory_sha256=reader.sha(trajectory), old_readback_untouched=True,
        wall_seconds=time.perf_counter()-begin)
    (output/'summary.json').write_text(json.dumps(reader.clean(summary), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(reader.clean({key: summary[key] for key in
        ('status', 'samples', 'accepted_steps', 'particle_count', 'accepted_time_s', 'wall_seconds')})))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=BASE/'readback_refined')
    arguments = parser.parse_args()
    package(arguments.output.resolve())

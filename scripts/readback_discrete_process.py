"""Read and audit the one v3 accepted chain; never evaluates the solver.

Chain: independent initial -> quarters (4) -> post-start halves (2) -> the
published pilot prefix. Alternative macro/half branches are never concatenated.
Run after the pilot manifest has left RUNNING. No interpolation or new physics.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from pysph.solver.utils import load
from afsd_backend.pysph.output import restore
from afsd_backend.pysph.interval_acceptance import (
    validate_interval_energy_state, validate_bore_crossing_state)

BASE = ROOT/'outputs/discrete_contact_process_20261004/v3'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: clean(val) for key, val in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(val) for val in value]
    return value


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fields(pa):
    order = np.argsort(pa.particle_id)
    return {name: prop.get_npy_array().reshape(len(order), -1)[order].copy()
            for name, prop in pa.properties.items()}


def identical_fields(first, second, context):
    require(set(first) == set(second), context+': property names differ')
    for key in first:
        require(np.array_equal(first[key], second[key]), context+': '+key+' differs')


def xyz(pa):
    return np.column_stack((pa.x, pa.y, pa.z))[np.argsort(pa.particle_id)]


def expected_pose(parameters, time_s):
    angle = parameters['omega']*time_s
    c, s = np.cos(angle), np.sin(angle)
    return dict(origin=np.array([parameters['traverse']*time_s, 0., parameters['layer']]),
                rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


def audit(output):
    begin = time.perf_counter()
    stages = [('quarters', BASE/'quarters'),
              ('post_start_halves', BASE/'post_start/halves'),
              ('pilot', BASE/'post_start/pilot')]
    manifests = {name: read(path/'manifest.json') for name, path in stages}
    require(not manifests['pilot']['status'].startswith('RUNNING'),
            'pilot is still running; rerun this read-only command after completion')
    require(manifests['quarters']['status'] == 'PASS_SHORT_TIME_REFINEMENT_NOT_TRANSPORT',
            'quarter startup screen did not pass')
    require(manifests['post_start_halves']['status'] == 'PASS_POST_START_TIME_SCREEN',
            'post-start h/h2 screen did not pass')
    require(manifests['quarters']['temporal_comparison']['pass'] and
            manifests['post_start_halves']['temporal_comparison']['pass'],
            'required temporal screen absent')
    classification = manifests['quarters']['classification']
    initial_file = BASE/'quarters/independent_initial.npz'
    initial, initial_state = restore(initial_file)
    initial_fields = fields(initial)
    require(initial_state['time_s'] == 0. and initial_state['step'] == 0,
            'chain must start at the declared independent initialization')
    validate_interval_energy_state(initial_state['interval_energy'])
    validate_bore_crossing_state(initial_state['bore_observation'])
    snapshots = [(initial, initial_state)]
    field_frames = [initial_fields]
    source_digests = {str(initial_file.relative_to(ROOT)): sha(initial_file)}
    rows, raw_records, thermal_rows, pose_rows = [], [], [], []
    previous_file, previous_pa, previous_state = initial_file, initial, initial_state
    baseline = read(BASE/'cost/manifest.json')['source_sha256']
    parameters = None

    # The manifest is the authority for the outer controller; checkpoints also
    # carry the unchanged physical baseline and post-start controller metadata.
    for stage_name, directory in stages:
        meta = manifests[stage_name]
        require(meta['classification'] == classification, 'stage classification changed')
        for name, digest in meta['source_sha256'].items():
            require(sha(ROOT/name) == digest, 'current source differs from frozen stage: '+name)
        source_digests[str((directory/'manifest.json').relative_to(ROOT))] = sha(directory/'manifest.json')
        history = read(directory/'accepted_history.json')
        require(len(history) == meta['accepted_steps'], 'manifest/history accepted count differs')
        if stage_name == 'quarters':
            require(len(history) == 4, 'four accepted startup steps required')
        else:
            require(meta['parent_sha256'] == sha(previous_file), 'continuation branch parent differs')
            require(ROOT/meta['parent_checkpoint'] == previous_file, 'continuation parent path differs')
            require(meta['starting_time_s'] == previous_state['time_s'], 'continuation time differs')
            parent_pa, parent_state = restore(directory/'continuation_initial.npz')
            identical_fields(fields(parent_pa), fields(previous_pa), 'continuation initial fields')
            for key in ('time_s', 'step', 'classification', 'interval_energy',
                        'bore_observation', 'contact_history', 'free'):
                require(clean(parent_state[key]) == clean(previous_state[key]),
                        'continuation initial history differs: '+key)
            if stage_name == 'post_start_halves':
                require(len(history) == 2, 'two accepted post-start half steps required')
        for entry in history:
            path = directory/entry['checkpoint']
            require(sha(path) == entry['sha256'], 'accepted checkpoint hash differs')
            source_digests[str(path.relative_to(ROOT))] = sha(path)
            pa, state = restore(path)
            data = fields(pa)
            if parameters is None:
                parameters = state['parameters']
            require(state['parameters'] == parameters, 'physical parameters changed within chain')
            require(state['source_sha256'] == baseline, 'checkpoint physical source changed')
            require(state['parameter_provenance'] == meta['provenance'], 'parameter provenance differs')
            require(state['classification'] == classification, 'checkpoint model changed')
            require(state['step'] == previous_state['step']+1, 'global physical step is discontinuous')
            require(state['time_s'] == entry['time_s'], 'checkpoint/index time differs')
            require(abs(state['time_s']-previous_state['time_s']-entry['dt_s']) <= 1e-15,
                    'accepted chain time has a gap or overlap')
            for key in ('particle_id', 'material_id', 'component_id', 'tracer_type',
                        'wholebar_tag', 'actuator', 'm', 'x0', 'y0', 'z0'):
                require(np.array_equal(data[key], initial_fields[key]), 'persistent field changed: '+key)
            require(set(data) == set(initial_fields), 'particle property set changed')
            for name, value in data.items():
                require(np.isfinite(value).all(), 'nonfinite particle field: '+name)
            staged = load(str(path.parent/'UNPUBLISHED_candidate.npz'))
            require(staged['solver_data'].get('restart_allowed') is False,
                    'unpublished snapshot incorrectly permits restart')
            identical_fields(data, fields(staged['arrays']['solid']), 'published/staged all fields')
            for key, value in state.items():
                require(clean(staged['solver_data'][key]) == clean(value),
                        'published/staged state differs: '+key)
            validate_interval_energy_state(state['interval_energy'])
            validate_bore_crossing_state(state['bore_observation'])
            energy = read(path.parent/'energy.json')
            require(clean(state['energy']) == energy['mechanical'], 'checkpoint/mechanical ledger differs')
            require(state['interval_energy'] == energy['interval'], 'checkpoint/interval ledger differs')
            require(state['bore_observation'] == read(path.parent/'transport.json'),
                    'checkpoint/transport record differs')
            raw = state['energy']['global_raw_residual_J']
            scale = state['energy']['absolute_work_scale_J']
            raw_records.append(dict(global_raw_residual_J=raw, absolute_work_scale_J=scale))
            require(state['interval_energy']['records'] == raw_records,
                    'cumulative ledger reset or added another branch')
            require(state['interval_energy']['steps'] == state['step'], 'interval step count differs')
            require(np.all(pa.temperature >= parameters['room']) and
                    np.all(pa.temperature < parameters['melt']), 'temperature domain violation')
            require(np.all(pa.ep >= 0.) and np.all(pa.strain_rate >= 0.), 'negative solved strain/rate')
            require(np.all(data['ep'] >= field_frames[-1]['ep']), 'plastic history decreased')
            require(np.max(np.abs(pa.rho/parameters['rho0']-1)) <= parameters['density_limit'],
                    'density gate violated')
            thermal = energy['thermal']
            measured = float(np.sum(data['m'].ravel()*parameters['cp']*
                                   (data['temperature'].ravel()-field_frames[-1]['temperature'].ravel())))
            require(abs(measured-thermal['measured_thermal_increment_J']) <= 1e-12,
                    'readback mass-capacity-temperature increment differs from thermal ledger')
            require(abs(thermal['mc_delta_T_residual_J']) <= 1e-12+1e-10*abs(thermal['total_heat_input_J']),
                    'thermal sum balance exceeded reported algebraic scale')
            pose = read(path.parent/'tool_pose.json')
            expected = expected_pose(parameters, state['time_s'])
            for key in ('origin', 'rotation'):
                require(np.allclose(pose[key], expected[key], atol=1e-14, rtol=0.),
                        'recorded tool pose differs from prescribed program')
            candidate = np.load(path.parent/'last_UNCOMMITTED_candidate.npz', allow_pickle=False)
            order = np.argsort(candidate['particle_id'])
            require(np.array_equal(candidate['particle_id'][order], data['particle_id'].ravel()),
                    'candidate persistent IDs differ')
            checks = [('position_m', xyz(pa), 1e-15),
                      ('velocity_m_s', np.column_stack((pa.u, pa.v, pa.w))[np.argsort(pa.particle_id)], 1e-12),
                      ('density_kg_m3', data['rho'].ravel(), 1e-9),
                      ('stress_Pa', data['S'].reshape(-1, 3, 3), 1e-6),
                      ('temperature_K', data['temperature'].ravel(), 1e-9),
                      ('plastic_strain', data['ep'].ravel(), 1e-14)]
            for key, actual, tolerance in checks:
                require(np.allclose(candidate[key][order], actual, atol=tolerance, rtol=0.),
                        'candidate and fresh published field differ: '+key)
            journal = read(path.parent/'last_UNCOMMITTED_candidate.json')
            require(journal['metrics']['force_residual_N'] <= journal['metrics']['force_limit_N'],
                    'last full candidate did not meet force gate')
            require(journal['contact']['sweep']['certified'], 'saved swept geometry not certified')
            require(float(np.max(journal['contact']['penetration_m'])) <= .1*parameters['dx'],
                    'saved penalty penetration gate violated')
            if stage_name != 'quarters':
                geometry = read(path.parent/'geometry_audit.json')
                require(not geometry['transient_contacts_missed_by_endpoints'] and
                        not geometry['new_non_bore_features'], 'post-start contact event gate violated')
                require(state['continuation_parent_sha256'] == meta['parent_sha256'],
                        'controller parent provenance differs')
                for name, digest in state['controller_source_sha256'].items():
                    require(meta['source_sha256'][name] == digest, 'controller source provenance differs')
            require(state['crossing_ids'] == state['bore_observation']['crossing_ids'] and
                    state['crossing_mass_kg'] == state['bore_observation']['crossing_mass_kg'],
                    'crossing fields differ from verified observation')
            require(abs(float(np.sum(state['energy']['prescribed_reaction_N']*
                                     candidate['displacement_m']))-
                        state['energy']['prescribed_boundary_work_J']) <= 1e-12,
                    'prescribed reaction work readback differs')
            S = data['S'].reshape(-1, 3, 3)
            q = np.sqrt(1.5*np.sum(S*S, axis=(1, 2)))
            row = dict(stage=stage_name, step=state['step'], time_s=state['time_s'],
                dt_s=entry['dt_s'], checkpoint=str(path.relative_to(ROOT)), sha256=sha(path),
                particles=len(pa.x), force_residual_N=journal['metrics']['force_residual_N'],
                force_limit_N=journal['metrics']['force_limit_N'],
                temperature_min_K=float(pa.temperature.min()), temperature_max_K=float(pa.temperature.max()),
                q_max_Pa=float(q.max()), ep_max=float(pa.ep.max()),
                max_net_displacement_m=float(np.linalg.norm(xyz(pa)-xyz(initial), axis=1).max()),
                max_step_displacement_m=float(np.linalg.norm(xyz(pa)-xyz(previous_pa), axis=1).max()),
                global_raw_residual_J=raw, interval_limit_J=state['interval_energy']['interval_limit_J'],
                cumulative_absolute_step_residual_J=state['interval_energy']['cumulative_absolute_step_residual_J'],
                crossing_mass_kg=state['crossing_mass_kg'], published_particle_properties=len(data),
                wall_seconds=entry['wall_seconds'])
            rows.append(row); thermal_rows.append(thermal); pose_rows.append(pose)
            snapshots.append((pa, state)); field_frames.append(data)
            previous_file, previous_pa, previous_state = path, pa, state
        require(meta['accepted_time_s'] == previous_state['time_s'], 'stage final accepted time differs')
        if history:
            require(meta['interval_energy'] == previous_state['interval_energy'],
                    'stage manifest final interval differs')

    require(parameters is not None, 'no accepted chain fields')
    require(not output.exists(), 'refuse to overwrite a prior readback directory')
    output.mkdir(parents=True)
    poses = [expected_pose(parameters, 0.)]+pose_rows
    times = np.array([state['time_s'] for _, state in snapshots])
    arrays = {f'particle_property__{key}': np.stack([frame[key] for frame in field_frames])
              for key in initial_fields}
    stress = arrays['particle_property__S'].reshape(len(times), len(initial.x), 3, 3)
    density = arrays['particle_property__rho'][..., 0]
    bulk = 2*parameters['shear']*(1+parameters['poisson'])/(3*(1-2*parameters['poisson']))
    pressure = bulk*(density/parameters['rho0']-1.)
    arrays.update(time_s=times, particle_id=initial_fields['particle_id'].ravel(),
        position_m=np.stack([xyz(pa) for pa, _ in snapshots]),
        temperature_K=arrays['particle_property__temperature'][..., 0],
        density_kg_m3=density, deviatoric_stress_Pa=stress,
        cauchy_stress_Pa=stress-pressure[..., None, None]*np.eye(3),
        plastic_strain=arrays['particle_property__ep'][..., 0],
        tool_origin_m=np.array([pose['origin'] for pose in poses]),
        tool_rotation=np.array([pose['rotation'] for pose in poses]),
        contact_history_m=np.stack([state['contact_history'] for _, state in snapshots]),
        crossing_mass_kg=np.array([state.get('crossing_mass_kg', 0.) for _, state in snapshots]),
        solved_rate_sample=np.array([False]+[True]*len(rows)),
        actual_checkpoint_samples_only=np.array(True), interpolation_performed=np.array(False))
    np.savez_compressed(output/'accepted_trajectory.npz', **arrays)
    thermal_sum = {key: math.fsum(row[key] for row in thermal_rows) for key in thermal_rows[0]}
    summary = dict(status='PASS_READBACK_AUDIT_NOT_PROCESS_OR_EXPERIMENTAL_VALIDATION',
        classification=classification, physical_evaluations=0, accepted_time_s=times[-1],
        samples=len(times), accepted_steps=len(rows), particle_count=len(initial.x),
        branch_chain=['independent_initial', 'quarters', 'post_start_halves', 'pilot_published_prefix'],
        alternative_macro_branches_excluded=True, inherited_interval_never_reset=True,
        pilot_status=manifests['pilot']['status'], interval_energy=previous_state['interval_energy'],
        crossing_mass_kg=previous_state['crossing_mass_kg'], crossing_ids=previous_state['crossing_ids'],
        thermal_whole_interval=thermal_sum, source_and_input_sha256=source_digests,
        readback_script_sha256=sha(Path(__file__)), trajectory_sha256=sha(output/'accepted_trajectory.npz'),
        particle_properties_saved=sorted(initial_fields),
        initial_zero_rates_are_not_computed_rates=True,
        startup_temporal_screen=manifests['quarters']['temporal_comparison'],
        post_start_temporal_screen=manifests['post_start_halves']['temporal_comparison'],
        later_pilot_time_accuracy='NOT_REFINEMENT_VALIDATED_AT_EVERY_PILOT_ENDPOINT',
        rows=rows, wall_seconds=time.perf_counter()-begin)
    (output/'summary.json').write_text(json.dumps(clean(summary), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(clean({key: summary[key] for key in
        ('status', 'accepted_time_s', 'samples', 'accepted_steps', 'particle_count',
         'pilot_status', 'crossing_mass_kg', 'wall_seconds')})))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=BASE/'readback')
    arguments = parser.parse_args()
    audit(arguments.output.resolve())

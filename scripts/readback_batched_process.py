"""Audit/package the actual 136-step fine chain, optionally a linked event prefix.

No solver/contact/material evaluation, no interpolated trajectory frames, and
no modification of existing evidence. Pure crossing observation is replayed
on the stored segments under its existing declared observation contract.
"""
from pathlib import Path
import argparse
import json
import math
import time

import numpy as np

import readback_discrete_process as reader
from afsd_backend.pysph.interval_acceptance import (
    initial_bore_crossing_state, observe_square_bore_crossings)

ROOT, OLD = reader.ROOT, reader.BASE
BASE = ROOT/'outputs/coupled_cost_20261004/batched'
require, read, sha, clean = reader.require, reader.read, reader.sha, reader.clean


def _default_stages():
    return [
        ('quarters', OLD/'quarters', 4, 'PASS_SHORT_TIME_REFINEMENT_NOT_TRANSPORT', 'initial'),
        ('old_halves', OLD/'post_start/halves', 2, 'PASS_POST_START_TIME_SCREEN', 'old'),
        ('pilot_refinement', OLD/'post_start/pilot_refinement', 16, 'PASS_1_TO_5MS_TIME_SCREEN_NOT_TRANSPORT', 'old'),
        ('batch_halves', BASE/'halves', 2, 'PASS_5_TO_6MS_TIME_SCREEN_NOT_TRANSPORT', 'batch'),
        ('window_eighth', BASE/'window_6_20ms/refinement/eighth', 112, 'PASS_WINDOW_REFINEMENT_NOT_TRANSPORT', 'batch')]


def _matching_journal(path, pa, state):
    """Select the successful full residual journal, including fallback folders."""
    data = reader.fields(pa)
    order_pa = np.argsort(pa.particle_id)
    checks = [('position_m', reader.xyz(pa), 1e-15),
              ('velocity_m_s', np.column_stack((pa.u, pa.v, pa.w))[order_pa], 1e-12),
              ('density_kg_m3', data['rho'].ravel(), 1e-9),
              ('stress_Pa', data['S'].reshape(-1, 3, 3), 1e-6),
              ('temperature_K', data['temperature'].ravel(), 1e-9),
              ('plastic_strain', data['ep'].ravel(), 1e-14),
              ('contact_history', np.asarray(state['contact_history'])[order_pa], 1e-15)]
    matches = []
    for candidate_path in sorted(path.parent.rglob('last_UNCOMMITTED_candidate.npz')):
        journal_path = candidate_path.with_suffix('.json')
        if not journal_path.exists():
            continue
        journal = read(journal_path)
        if journal['metrics']['force_residual_N'] > journal['metrics']['force_limit_N']:
            continue
        with np.load(candidate_path, allow_pickle=False) as saved:
            order = np.argsort(saved['particle_id'])
            if not np.array_equal(saved['particle_id'][order], data['particle_id'].ravel()):
                continue
            if not all(np.allclose(saved[key][order], actual, rtol=0., atol=tolerance)
                       for key, actual, tolerance in checks):
                continue
            candidate = {key: saved[key].copy() for key in saved.files}
        matches.append((candidate_path, journal_path, candidate, journal))
    require(bool(matches), 'no converged full-residual journal matches '+str(path))
    # More than one numerically identical successful journal is harmless, but
    # every selected field still has to match the independently published data.
    return matches[-1]


def package(output, append_stages=()):
    begin = time.perf_counter()
    require(not output.exists(), 'refuse to overwrite existing readback directory')
    stages = _default_stages()
    for number, directory in enumerate(append_stages, 1):
        directory = Path(directory).resolve()
        require(directory.parent == (BASE/'first_shoulder_event').resolve() and directory.name in ('approach', 'halves'),
                'only explicit first-shoulder approach/halves accepted prefixes may be appended')
        stages.append(('event_prefix_'+str(number), directory, None, None, 'event'))
    hashes = {}

    def bind(path, expected=None):
        digest = sha(path)
        if expected is not None:
            require(digest == expected, 'file hash mismatch: '+str(path))
        hashes[path.relative_to(ROOT).as_posix()] = digest
        return digest

    cost_file = OLD/'cost/manifest.json'
    cost = read(cost_file); bind(cost_file)
    baseline, provenance = cost['source_sha256'], cost['provenance']
    bind(ROOT/provenance['parameter_source'], provenance['source_sha256'])
    initial_file = OLD/'quarters/independent_initial.npz'
    initial, initial_state = reader.restore(initial_file)
    bind(initial_file)
    initial_data = reader.fields(initial)
    require(initial_state['time_s'] == 0. and initial_state['step'] == 0,
            'chain must start at independent time zero')
    reader.validate_interval_energy_state(initial_state['interval_energy'])
    reader.validate_bore_crossing_state(initial_state['bore_observation'])
    classification = initial_state['classification']
    field_frames, snapshots = [initial_data], [(initial, initial_state)]
    previous_file, previous_pa, previous_state = initial_file, initial, initial_state
    raw_records, thermal_rows, rows, pose_rows, stage_rows = [], [], [], [], []
    contact_active_frames = [np.zeros((len(initial.x), 8), bool)]
    # No contact evaluator is called at initialization. Its zero-filled contact
    # array is explicitly masked unavailable and its count exported as -1;
    # these values must not be presented as computed contact observations.
    parameters = None

    for stage_name, directory, count, expected_status, mode in stages:
        manifest_path, history_path = directory/'manifest.json', directory/'accepted_history.json'
        meta, history = read(manifest_path), read(history_path)
        bind(manifest_path); bind(history_path)
        require('RUNNING' not in meta['status'], 'stage is still running: '+stage_name)
        if expected_status is not None:
            require(meta['status'] == expected_status, 'required successful stage missing: '+stage_name)
        elif directory.name == 'approach':
            require(meta['status'] in ('EVENT_CANDIDATE_NOT_ACCEPTED', 'BOUNDED_APPROACH_COMPLETE_NO_EVENT', 'FAILED_WITH_ACCEPTED_PREFIX'),
                    'approach has no finalized accepted prefix')
        elif directory.name == 'halves':
            require(meta['status'] in ('PASS_FIRST_SHOULDER_LOCAL_TIME_SCREEN_NOT_TRANSPORT', 'FAILED_WITH_ACCEPTED_PREFIX'),
                    'event halves have no qualified accepted prefix')
            require(meta['global_time_screen']['pass'] and meta['local_time_screen']['passed'],
                    'event halves require BOTH global and local time screens')
            screen_path = directory/'temporal_comparison.json'
            screen = read(screen_path); bind(screen_path)
            require(screen['global_screen'] == meta['global_time_screen'] and screen['local_screen'] == meta['local_time_screen'],
                    'event temporal record differs from manifest')
        if 'temporal_comparison' in meta:
            require(meta['temporal_comparison']['pass'], 'failed time screen cannot qualify appended chain')
            screen_path = directory/'temporal_comparison.json'
            if screen_path.exists():
                bind(screen_path)
        require(meta['classification'] == classification, 'stage model classification changed')
        require(meta['provenance'] == provenance, 'stage parameter provenance changed')
        require(len(history) == meta.get('accepted_steps', len(history)), 'manifest/history count differs')
        if count is not None:
            require(len(history) == count, 'unexpected stage accepted count: '+stage_name)
        require(bool(history), 'empty event stage cannot extend trajectory')
        for name, digest in meta['source_sha256'].items():
            bind(ROOT/name, digest)
        if mode != 'initial':
            require(meta['parent_sha256'] == sha(previous_file), 'stage direct branch-parent SHA differs')
            require((ROOT/meta['parent_checkpoint']).resolve() == previous_file.resolve(), 'stage parent path differs')
            require(meta['starting_time_s'] == previous_state['time_s'], 'stage starting time differs')
            parent_path = directory/'continuation_initial.npz'
            parent_pa, parent_state = reader.restore(parent_path); bind(parent_path)
            reader.identical_fields(reader.fields(parent_pa), reader.fields(previous_pa), 'continuation fields')
            for key in ('time_s', 'step', 'classification', 'interval_energy', 'bore_observation', 'contact_history', 'free'):
                require(clean(parent_state[key]) == clean(previous_state[key]), 'continuation metadata differs: '+key)
        start_step = previous_state['step']

        for entry in history:
            path = directory/entry['checkpoint']
            bind(path, entry['sha256'])
            pa, state = reader.restore(path)
            data = reader.fields(pa)
            if parameters is None:
                parameters = state['parameters']
                observed_initial = initial_bore_crossing_state(reader.xyz(initial),
                    initial_data['particle_id'].ravel(), initial_data['component_id'].ravel(),
                    initial_data['m'].ravel(), 0., reader.expected_pose(parameters, 0.), parameters['bore_width'])
                require(observed_initial == initial_state['bore_observation'], 'initial eligible labels differ')
            require(state['parameters'] == parameters, 'physical parameters changed')
            require(state['source_sha256'] == baseline, 'physical source baseline changed')
            require(state['parameter_provenance'] == provenance, 'checkpoint input provenance changed')
            require(state['classification'] == classification, 'checkpoint classification changed')
            require(state['step'] == previous_state['step']+1, 'physical global step gap')
            require(entry.get('global_step', state['step']) == state['step'], 'history global step differs')
            require(state['time_s'] == entry['time_s'], 'checkpoint/history time differs')
            require(abs(state['time_s']-previous_state['time_s']-entry['dt_s']) <= 1e-15,
                    'physical chain has a gap or overlap')
            for key in ('particle_id', 'material_id', 'component_id', 'tracer_type', 'wholebar_tag', 'actuator', 'm', 'x0', 'y0', 'z0'):
                require(np.array_equal(data[key], initial_data[key]), 'persistent field changed: '+key)
            require(set(data) == set(initial_data), 'particle properties changed')
            for name, value in data.items():
                require(np.isfinite(value).all(), 'nonfinite field: '+name)
            require(np.array_equal(state['free'], initial_state['free']), 'mechanical constraint identity changed')
            if mode in ('batch', 'event'):
                direct_hash = sha(previous_file)
                require(entry['parent_sha256'] == direct_hash, 'history direct parent SHA differs')
                require(state['continuation_parent_sha256'] == direct_hash, 'checkpoint direct parent SHA differs')
                require(state['branch_origin_sha256'] == meta.get('branch_origin_sha256', meta['parent_sha256']),
                        'checkpoint branch origin differs')
                require(state['implementation'] == 'NUMPY_BATCHED_FINITE_PATCH_EQUIVALENT', 'implementation changed')
                require(state['implementation_source_sha256'] == meta['source_sha256'], 'implementation source metadata differs')
                require(state['controller_source_sha256'] == meta['source_sha256'], 'controller source metadata differs')
            elif mode == 'old':
                require(state['continuation_parent_sha256'] == meta['parent_sha256'], 'old controller branch parent differs')
                for name, digest in state['controller_source_sha256'].items():
                    require(meta['source_sha256'][name] == digest, 'old controller source metadata differs')

            staged_path = path.parent/'UNPUBLISHED_candidate.npz'
            staged = reader.load(str(staged_path)); bind(staged_path)
            require(staged['solver_data'].get('restart_allowed') is False, 'unpublished candidate permits restart')
            reader.identical_fields(data, reader.fields(staged['arrays']['solid']), 'staged/published fields')
            for key, value in state.items():
                require(clean(staged['solver_data'][key]) == clean(value), 'staged/published metadata differs: '+key)
            reader.validate_interval_energy_state(state['interval_energy'])
            reader.validate_bore_crossing_state(state['bore_observation'])
            energy_path, transport_path, pose_path = (path.parent/name for name in ('energy.json', 'transport.json', 'tool_pose.json'))
            energy, transport, pose = read(energy_path), read(transport_path), read(pose_path)
            for file in (energy_path, transport_path, pose_path):
                bind(file)
            require(clean(state['energy']) == energy['mechanical'], 'mechanical ledger differs')
            require(state['interval_energy'] == energy['interval'], 'interval ledger differs')
            require(state['bore_observation'] == transport, 'transport ledger differs')
            raw_records.append({key: state['energy'][key] for key in ('global_raw_residual_J', 'absolute_work_scale_J')})
            require(state['interval_energy']['records'] == raw_records, 'energy history reset or branch mixed')
            require(state['interval_energy']['steps'] == state['step'], 'energy prefix/global step mismatch')
            observer = observe_square_bore_crossings(previous_state['bore_observation'], reader.xyz(previous_pa), reader.xyz(pa),
                data['particle_id'].ravel(), data['component_id'].ravel(), data['m'].ravel(),
                previous_state['time_s'], state['time_s'], lambda t: reader.expected_pose(parameters, t))
            require(observer['next_state'] == state['bore_observation'], 'passive crossing observation failed independent replay')
            require(state['crossing_ids'] == transport['crossing_ids'] and state['crossing_mass_kg'] == transport['crossing_mass_kg'],
                    'stored crossing fields differ')
            require(np.all(pa.temperature >= parameters['room']) and np.all(pa.temperature < parameters['melt']), 'temperature domain violation')
            require(np.all(pa.ep >= 0.) and np.all(pa.strain_rate >= 0.), 'negative strain/rate')
            require(np.all(data['ep'] >= field_frames[-1]['ep']), 'plastic history decreased')
            require(np.max(np.abs(pa.rho/parameters['rho0']-1.)) <= parameters['density_limit'], 'density gate violation')
            thermal = energy['thermal']
            measured = float(np.sum(data['m'].ravel()*parameters['cp']*(data['temperature']-field_frames[-1]['temperature']).ravel()))
            require(abs(measured-thermal['measured_thermal_increment_J']) <= 1e-12, 'thermal increment readback differs')
            require(abs(thermal['mc_delta_T_residual_J']) <= 1e-12+1e-10*abs(thermal['total_heat_input_J']), 'thermal sum residual too large')
            expected = reader.expected_pose(parameters, state['time_s'])
            for key in ('origin', 'rotation'):
                require(np.allclose(pose[key], expected[key], atol=1e-14, rtol=0.), 'tool pose differs')

            candidate_path, journal_path, candidate, journal = _matching_journal(path, pa, state)
            bind(candidate_path); bind(journal_path)
            contact = journal['contact']
            require(contact['sweep']['certified'], 'candidate solid sweep not certified')
            require(np.max(contact['penetration_m']) <= .1*parameters['dx'], 'penetration limit exceeded')
            require(np.max(np.linalg.norm(candidate['displacement_m'], axis=1)) <= .1*parameters['dx'], 'single-step displacement limit exceeded')
            require(abs(float(np.sum(state['energy']['prescribed_reaction_N']*candidate['displacement_m']))-
                        state['energy']['prescribed_boundary_work_J']) <= 1e-12, 'reaction work differs')
            order = np.argsort(candidate['particle_id'])
            active = np.asarray(contact['active'], bool)[order]
            geometry_path = path.parent/'geometry_audit.json'
            if geometry_path.exists():
                geometry = read(geometry_path); bind(geometry_path)
                require(not geometry['transient_contacts_missed_by_endpoints'], 'saved geometry has hidden contact')
                require(np.array_equal(geometry['endpoint_counts'], active.sum(axis=0)), 'geometry/contact counts differ')
                if mode != 'event':
                    require(not geometry['new_non_bore_features'], 'unexpected pre-event new feature')
            elif mode != 'initial':
                raise ValueError('missing contact geometry audit')
            for shell_path in path.parent.glob('*shell*.json'):
                shell = read(shell_path); bind(shell_path)
                require(shell.get('status') == 'CERTIFIED_NO_HIDDEN_SHELL_CONTACT', 'saved shell certificate did not pass')
            if mode == 'event':
                require((path.parent/'continuous_shell_audit.json').exists(), 'event prefix missing continuous shell certificate')
                events_path = path.parent/'sampled_nonbore_events.json'
                events = read(events_path); bind(events_path)
                if directory.name == 'approach':
                    require(not events and not geometry['new_non_bore_features'], 'approach committed a new contact event')
                else:
                    require(set(geometry['new_non_bore_features']) <= {'SHOULDER'}, 'event outside shoulder scope')
                    equivalent_path = path.parent/'scalar_batch_equivalence.json'
                    equivalent = read(equivalent_path); bind(equivalent_path)
                    require(equivalent['status'] == 'PASS_SAME_CANDIDATE_COMPLETE_MAP', 'event scalar/batch equivalence absent')
            ids = data['particle_id'].ravel()
            eligible = np.isin(ids, transport['eligible_ids'])
            eligible_z = reader.xyz(pa)[eligible, 2]-expected['origin'][2]
            above = eligible_z > 0.
            S = data['S'].reshape(-1, 3, 3)
            row = dict(stage=stage_name, step=state['step'], time_s=state['time_s'], dt_s=entry['dt_s'],
                checkpoint=path.relative_to(ROOT).as_posix(), sha256=sha(path),
                previous_checkpoint_sha256=sha(previous_file), force_residual_N=journal['metrics']['force_residual_N'],
                force_limit_N=journal['metrics']['force_limit_N'],
                temperature_min_K=float(pa.temperature.min()), temperature_max_K=float(pa.temperature.max()),
                equivalent_deviatoric_stress_max_Pa=float(np.sqrt(1.5*np.sum(S*S, axis=(1, 2))).max()),
                plastic_strain_max=float(pa.ep.max()),
                max_net_displacement_m=float(np.linalg.norm(reader.xyz(pa)-reader.xyz(initial), axis=1).max()),
                max_step_displacement_m=float(np.linalg.norm(reader.xyz(pa)-reader.xyz(previous_pa), axis=1).max()),
                shoulder_active_count=int(active[:, 0].sum()), shoulder_active_ids=ids[active[:, 0]].astype(int).tolist(),
                eligible_crossing_count=len(transport['crossing_ids']), crossing_mass_kg=transport['crossing_mass_kg'],
                closest_eligible_above_shoulder_m=float(eligible_z[above].min()) if above.any() else None,
                minimum_eligible_signed_shoulder_distance_m=float(eligible_z.min()) if len(eligible_z) else None,
                global_raw_residual_J=state['energy']['global_raw_residual_J'],
                absolute_work_scale_J=state['energy']['absolute_work_scale_J'],
                interval_limit_J=state['interval_energy']['interval_limit_J'],
                cumulative_absolute_step_residual_J=state['interval_energy']['cumulative_absolute_step_residual_J'],
                thermal=thermal, wall_seconds=entry['wall_seconds'])
            rows.append(row); thermal_rows.append(thermal); pose_rows.append(pose)
            field_frames.append(data); snapshots.append((pa, state)); contact_active_frames.append(active)
            previous_file, previous_pa, previous_state = path, pa, state
        require(meta['accepted_time_s'] == previous_state['time_s'], 'stage final time differs')
        require(meta['interval_energy'] == previous_state['interval_energy'], 'stage final energy prefix differs')
        stage_rows.append(dict(name=stage_name, directory=directory.relative_to(ROOT).as_posix(), status=meta['status'],
            accepted_steps=len(history), initial_global_step=start_step, final_global_step=previous_state['step'],
            temporal_comparison=meta.get('temporal_comparison'), global_event_time_screen=meta.get('global_time_screen'),
            local_event_time_screen=meta.get('local_time_screen'), extra_event_prefix=(mode == 'event')))

    require(len(rows) >= 136 and abs(snapshots[136][1]['time_s']-.020) <= 1e-15, 'complete qualified 20ms chain missing')
    require(parameters is not None, 'no accepted physical data')
    for name, digest in hashes.items():
        require(sha(ROOT/name) == digest, 'evidence changed during readback: '+name)
    times = np.array([state['time_s'] for _, state in snapshots])
    poses = [reader.expected_pose(parameters, 0.)]+pose_rows
    arrays = {f'particle_property__{key}': np.stack([frame[key] for frame in field_frames]) for key in initial_data}
    stress = arrays['particle_property__S'].reshape(len(times), len(initial.x), 3, 3)
    density = arrays['particle_property__rho'][..., 0]
    bulk = 2*parameters['shear']*(1+parameters['poisson'])/(3*(1-2*parameters['poisson']))
    pressure = bulk*(density/parameters['rho0']-1.)
    arrays.update(time_s=times, particle_id=initial_data['particle_id'].ravel(),
        position_m=np.stack([reader.xyz(pa) for pa, _ in snapshots]),
        velocity_m_s=np.stack([np.column_stack((pa.u, pa.v, pa.w))[np.argsort(pa.particle_id)] for pa, _ in snapshots]),
        temperature_K=arrays['particle_property__temperature'][..., 0], density_kg_m3=density,
        deviatoric_stress_Pa=stress, cauchy_stress_Pa=stress-pressure[..., None, None]*np.eye(3),
        plastic_strain=arrays['particle_property__ep'][..., 0],
        tool_origin_m=np.array([pose['origin'] for pose in poses]), tool_rotation=np.array([pose['rotation'] for pose in poses]),
        contact_history_m=np.stack([state['contact_history'][np.argsort(pa.particle_id)] for pa, state in snapshots]),
        contact_active=np.stack(contact_active_frames),
        crossing_mass_kg=np.array([state['bore_observation']['crossing_mass_kg'] for _, state in snapshots]),
        eligible_crossing_count=np.array([len(state['bore_observation']['crossing_ids']) for _, state in snapshots]),
        shoulder_active_count=np.array([-1]+[int(frame[:, 0].sum()) for frame in contact_active_frames[1:]]),
        contact_activity_available=np.array([False]+[True]*len(rows)),
        solved_rate_sample=np.array([False]+[True]*len(rows)), actual_checkpoint_samples_only=np.array(True),
        interpolation_performed=np.array(False))
    output.mkdir(parents=True)
    trajectory = output/'accepted_trajectory.npz'
    np.savez_compressed(trajectory, **arrays)
    with np.load(trajectory, allow_pickle=False) as saved:
        require(set(saved.files) == set(arrays), 'exported arrays differ')
        for key, value in arrays.items():
            require(np.array_equal(saved[key], value), 'exported array differs: '+key)
    summary = dict(status='PASS_ACTUAL_CHAIN_READBACK_NOT_PHYSICAL_OR_EXPERIMENTAL_VALIDATION',
        classification=classification, physical_evaluations=0, samples=len(times), accepted_steps=len(rows),
        particle_count=len(initial.x), accepted_time_s=float(times[-1]), branch_chain=stage_rows,
        alternative_branches_excluded=True, inherited_interval_never_reset=True,
        initial_zero_rates_are_not_computed_rates=True, initial_contact_activity_unavailable_and_masked=True,
        passive_crossing_observer_replayed=True, shoulder_contact_is_not_centre_crossing=True,
        final=rows[-1], crossing_ids=previous_state['bore_observation']['crossing_ids'],
        crossing_events=previous_state['bore_observation']['events'], eligible_ids=previous_state['bore_observation']['eligible_ids'],
        interval_energy=previous_state['interval_energy'],
        thermal_whole_interval={key: math.fsum(row[key] for row in thermal_rows) for key in thermal_rows[0]},
        total_mass_kg=float(np.sum(initial.m)), particle_properties_saved=sorted(initial_data),
        source_and_input_sha256=hashes, readback_script_sha256=sha(Path(__file__)),
        unchanged_reader_helper_sha256=sha(Path(reader.__file__)), trajectory_sha256=sha(trajectory),
        rows=rows, wall_seconds=time.perf_counter()-begin)
    (output/'summary.json').write_text(json.dumps(clean(summary), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(clean({key: summary[key] for key in ('status', 'samples', 'accepted_steps', 'particle_count', 'accepted_time_s', 'wall_seconds')})))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=BASE/'readback_actual_chain')
    parser.add_argument('--append-stage', type=Path, action='append', default=[])
    arguments = parser.parse_args()
    package(arguments.output.resolve(), arguments.append_stage)

"""Readback audit and trajectory packaging only; never executes a physical step."""
from pathlib import Path
import json
import sys
import numpy as np
from run_hot_quasistatic_probe import ROOT, BASE, LOGS, sha, write, sources
from afsd_backend.pysph.output import restore
from afsd_backend.pysph.quasistatic_process import pose


def main():
    directory = BASE/'readback'
    directory.mkdir(exist_ok=False)
    refinement = BASE/'refinement'
    report = json.loads((refinement/'manifest.json').read_text())
    history = json.loads((refinement/'accepted_history.json').read_text())
    attempts = json.loads((refinement/'attempts.json').read_text())
    root = json.loads((BASE/'root/manifest.json').read_text())
    cost = json.loads((BASE/'cost/manifest.json').read_text())
    tests = json.loads((LOGS/'tests.json').read_text())
    assert tests['exit_code'] == 0 and sources() == cost['source_sha256'] == root['source_sha256']
    assert report['accepted_steps'] == len(history) == sum(a['accepted'] for a in attempts)
    assert report['accepted_time_s'] == history[-1]['time_s']
    pa, state = restore(refinement/'independent_initial.npz')
    p = json.loads((BASE/'cost/resolved_config.json').read_text())['parameters']
    ids = pa.particle_id.copy()
    identity = {key: getattr(pa, key).copy() for key in
                ('material_id', 'component_id', 'tracer_type', 'wholebar_tag', 'actuator')}
    initial = np.column_stack((pa.x, pa.y, pa.z)); mass = pa.m.copy()
    times = [0.]; positions = [initial]; velocities = [np.column_stack((pa.u, pa.v, pa.w))]
    temperatures = [pa.temperature.copy()]; stresses = [pa.S.reshape(-1,3,3).copy()]
    strains = [pa.ep.copy()]; densities = [pa.rho.copy()]
    contact_history = [state['contact_history'].copy()]
    audit = []; checksums = {}
    for entry in history:
        path = refinement/entry['checkpoint']
        assert sha(path) == entry['checkpoint_sha256']
        checksums[str(path.relative_to(ROOT))] = sha(path)
        pa, state = restore(path)
        assert np.array_equal(ids, pa.particle_id) and np.array_equal(mass, pa.m)
        for key, value in identity.items(): assert np.array_equal(value, getattr(pa, key))
        assert state['time_s'] == entry['time_s'] and state['time_s'] > times[-1]
        assert state['classification'] == report['classification']
        r = json.loads((path.parent/'latest_UNCOMMITTED_candidate.json').read_text())
        energy = json.loads((path.parent/'energy.json').read_text())
        mech = energy['mechanical']; thermal = energy['thermal']
        assert r['metrics']['scaled_residual'] <= 1
        assert abs(mech['global_raw_residual_J']) <= mech['registered_energy_limit_J']
        assert r['contact']['sweep']['certified']
        assert np.max(r['contact']['penetration_m']) <= .1*p['dx']
        assert r['metrics']['density_residual'] <= 1e-10
        assert r['metrics']['thermal_residual_K'] <= 1e-8
        assert np.all(pa.ep >= strains[-1])
        times.append(state['time_s']); positions.append(np.column_stack((pa.x,pa.y,pa.z)))
        velocities.append(np.column_stack((pa.u,pa.v,pa.w)))
        temperatures.append(pa.temperature.copy()); stresses.append(pa.S.reshape(-1,3,3).copy())
        strains.append(pa.ep.copy()); densities.append(pa.rho.copy())
        contact_history.append(state['contact_history'].copy())
        audit.append(dict(time_s=state['time_s'], dt_s=entry['dt_s'],
            force_residual_N=entry['force_residual_N'], energy_raw_J=entry['energy_raw_J'],
            energy_limit_J=entry['energy_limit_J'],
            absolute_allowance_fraction_of_energy_limit=1e-7/entry['energy_limit_J'],
            maximum_overlap_m=float(np.max(r['contact']['penetration_m'])),
            active_contact_features=int(np.sum(r['contact']['active'])),
            material_discrete_identity_error_J=mech['material_breakdown']['material_discrete_identity_error_J'],
            measured_thermal_increment_J=thermal['measured_thermal_increment_J'],
            thermal_balance_residual_J=thermal['mc_delta_T_residual_J']))
    above_initial = (pa.component_id == 1) & (initial[:,2] >= p['layer'])
    center_gap = float(np.min(pa.z[above_initial]-p['layer']))
    assert center_gap > 0 and not state['crossing_ids']
    net = np.linalg.norm(positions[-1]-initial, axis=1)
    free_bar = (pa.component_id == 1) & (~pa.actuator.astype(bool))
    origins = np.array([pose(p,t)['origin'] for t in times])
    rotations = np.array([pose(p,t)['rotation'] for t in times])
    np.savez_compressed(directory/'accepted_trajectory.npz', classification=report['classification'],
        time_s=times, particle_id=ids, position_m=positions, mean_velocity_m_s=velocities,
        temperature_K=temperatures, deviatoric_stress_Pa=stresses, plastic_strain=strains,
        density_kg_m3=densities, mass_kg=mass, contact_history_m=contact_history,
        tool_origin_m=origins, tool_rotation=rotations, **identity)
    with np.load(directory/'accepted_trajectory.npz', allow_pickle=False) as traj:
        np.testing.assert_array_equal(traj['particle_id'], ids)
        np.testing.assert_array_equal(traj['position_m'][-1], positions[-1])
        assert len(traj['time_s']) == len(history)+1
    nominal = .00025/p['feed']
    summary = dict(classification=report['classification'], analysis_only_no_new_physical_evaluations=True,
        status='READBACK_PASS_STARTUP_TARGET_NOT_REACHED', selected_tests_passed=48,
        newly_added_tests=39, earlier_q1_tests_rerun=9, full_regression='NOT_RUN',
        old_dynamic_time_s=root['original_dynamic_time_s'], independent_quasistatic_time_s=times[-1],
        accepted_steps=len(history), rejected_steps=sum(not a['accepted'] for a in attempts),
        startup_target_s=report['target_time_s'], target_met=False,
        all_accepted_ids_tags_mass_and_fields_read_back=True,
        physical_duration_not_summed_with_dynamic_model=True,
        crossing_mass_kg=0., minimum_previously_above_bar_center_gap_m=center_gap,
        max_net_displacement_m=float(net.max()), max_free_bar_net_displacement_m=float(net[free_bar].max()),
        max_plastic_strain=float(pa.ep.max()), maximum_equivalent_stress_Pa=float(np.sqrt(1.5*np.sum(pa.S.reshape(-1,3,3)**2,axis=(1,2))).max()),
        temperature_range_K=[float(pa.temperature.min()), float(pa.temperature.max())],
        refinement_wall_s=report['wall_seconds'], root_wall_s=root['wall_seconds'],
        nominal_250um_feed_time_s=nominal,
        naive_startup_rate_extrapolation_to_nominal_feed_wall_s=report['wall_seconds']/times[-1]*nominal,
        extrapolation_is_not_a_production_prediction=True,
        constitutive_or_temporal_accuracy='NOT_VALIDATED',
        accepted_snapshot_sha256=checksums, accepted_step_audit=audit,
        source_sha256={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))})
    write(directory/'summary.json', summary)
    write(directory/'artifacts.json', {p.name:sha(p) for p in directory.iterdir() if p.is_file() and p.name!='artifacts.json'})
    print(json.dumps({k:v for k,v in summary.items() if k not in ('accepted_snapshot_sha256','source_sha256','accepted_step_audit')},indent=2))


if __name__ == '__main__': main()

"""Lightweight saved-field energy/deformation audit of the refined 30 ms chain."""
from pathlib import Path
import json
import math
import time
import numpy as np
import readback_discrete_process as reader
from afsd_backend.pysph.mechanics import bulk_modulus, volumetric_energy

ROOT = reader.ROOT
BASE = ROOT/'outputs/coupled_cost_20261004/batched/readback_30ms_refined'
OUTPUT = ROOT/'outputs/coupled_cost_20261004/energy30ms_audit.json'


def rigid_residual(x0, x1):
    a, b = x0-x0.mean(axis=0), x1-x1.mean(axis=0)
    left, _, right = np.linalg.svd(a.T@b)
    adjust = np.eye(3); adjust[-1, -1] = np.linalg.det(left@right)
    distance = np.linalg.norm(b-a@(left@adjust@right), axis=1)
    return dict(rms_m=float(np.sqrt(np.mean(distance**2))), maximum_m=float(distance.max()))


def main():
    begin = time.perf_counter(); reader.require(not OUTPUT.exists(), 'refuse to overwrite 30 ms audit')
    summary_path, trajectory_path = BASE/'summary.json', BASE/'accepted_trajectory.npz'
    summary = reader.read(summary_path)
    reader.require(summary['accepted_steps'] == 281 and summary['samples'] == 282 and
        abs(summary['accepted_time_s']-.030) < 1e-15 and
        reader.sha(trajectory_path) == summary['trajectory_sha256'], 'wrong refined actual chain')
    for name, digest in summary['finite_refinement_source_and_evidence_sha256'].items():
        reader.require(reader.sha(ROOT/name) == digest, 'readback source/input changed: '+name)
    hashes = {path.relative_to(ROOT).as_posix(): reader.sha(path) for path in
        (Path(__file__), summary_path, trajectory_path)}
    with np.load(trajectory_path, allow_pickle=False) as saved:
        time_s, x, rho, S, T, ep, velocity = (saved[name].copy() for name in
            ('time_s', 'position_m', 'density_kg_m3', 'deviatoric_stress_Pa', 'temperature_K', 'plastic_strain', 'velocity_m_s'))
        field = {name: saved['particle_property__'+name][..., 0].copy()
            for name in ('particle_id', 'm', 'plastic_work', 'component_id', 'actuator', 'strain_rate')}
    final_path = ROOT/summary['final']['checkpoint']
    final_pa, state = reader.restore(final_path); p = state['parameters']
    hashes[final_path.relative_to(ROOT).as_posix()] = reader.sha(final_path)
    order = np.argsort(final_pa.particle_id)
    reader.require(np.array_equal(field['particle_id'][-1], final_pa.particle_id[order]), 'trajectory final IDs differ')
    mass = field['m'][0]
    reader.require(np.array_equal(field['m'], np.broadcast_to(mass, field['m'].shape)), 'mass changed')
    K = bulk_modulus(p['shear'], p['poisson'])
    stored_dev = np.sum(mass[None, :]/rho*np.sum(S*S, axis=(2, 3))/(4*p['shear']), axis=1)
    stored_vol = np.sum(mass[None, :]*volumetric_energy(rho, p['rho0'], K), axis=1)
    deltaU = np.diff(stored_dev+stored_vol)
    increments = np.diff(field['plastic_work'], axis=0)
    reader.require(np.all(increments >= -1e-18), 'saved plastic-work accumulator decreased')
    energy_rows, raw_records = [], []
    max_energy_error = max_thermal_error = max_prefix_usage = 0.
    local_abs = local_signed = 0.
    for index, row in enumerate(summary['rows']):
        path = (ROOT/row['checkpoint']).parent/'energy.json'
        document = reader.read(path); e = document['mechanical']; thermal = document['thermal']
        hashes[path.relative_to(ROOT).as_posix()] = reader.sha(path)
        raw = e['prescribed_boundary_work_J']+e['tool_work_J']-(deltaU[index]+increments[index].sum()+
            e['normal_storage_increment_J']+e['tangential_storage_increment_J']+e['friction_sliding_J'])
        error = max(abs(deltaU[index]-e['material_storage_increment_J']),
            abs(increments[index].sum()-e['plastic_work_J']), abs(raw-e['global_raw_residual_J']))
        reader.require(error < 1e-11, 'raw energy does not reconstruct from actual fields')
        max_energy_error = max(max_energy_error, error)
        heat = float(np.sum(mass*p['cp']*(T[index+1]-T[index])))
        heat_error = abs(heat-thermal['measured_thermal_increment_J'])
        reader.require(heat_error < 1e-12, 'thermal ledger does not bind to actual fields')
        max_thermal_error = max(max_thermal_error, heat_error)
        raw_records.append({key: e[key] for key in ('global_raw_residual_J', 'absolute_work_scale_J')})
        reader.validate_interval_energy_state(document['interval'])
        reader.require(document['interval']['records'] == raw_records, 'actual accepted energy prefix reset or mixed')
        absolute = math.fsum(abs(record['global_raw_residual_J']) for record in raw_records)
        allowance = document['interval']['interval_limit_J']
        reader.require(absolute <= allowance, 'saved prefix allowance exceeded')
        max_prefix_usage = max(max_prefix_usage, absolute/allowance)
        local_abs += e['material_breakdown']['material_raw_energy_defect_absolute_sum_J']
        local_signed += e['material_breakdown']['material_raw_energy_defect_J']
        energy_rows.append(e)
    sums = {key: math.fsum(row[key] for row in energy_rows) for key, value in energy_rows[0].items()
        if isinstance(value, (int, float)) and not isinstance(value, bool) and key != 'registered_energy_limit_J'}
    Wext = sums['prescribed_boundary_work_J']+sums['tool_work_J']
    absolute = math.fsum(abs(record['global_raw_residual_J']) for record in raw_records)
    scale = math.fsum(record['absolute_work_scale_J'] for record in raw_records)
    plastic = field['plastic_work'][-1]-field['plastic_work'][0]
    floor = field['strain_rate'][1:] < p['rate_floor']
    floor_work = np.sum(np.where(floor, increments, 0.), axis=0)
    component, actuator = field['component_id'][0], field['actuator'][0].astype(bool)
    free = np.asarray(state['free'])[order]
    masks = dict(free_bar=(component == 1)&~actuator, top_bar_z_actuated=actuator,
        free_substrate=(component == 0)&free.any(axis=1), clamped_substrate=(component == 0)&~free.any(axis=1),
        initially_above_shoulder_bar=(component == 1)&(x[0, :, 2] > p['layer']),
        initially_below_shoulder_bar=(component == 1)&(x[0, :, 2] < p['layer']))
    groups = {}
    for name, mask in masks.items():
        groups[name] = dict(particles=int(mask.sum()), plastic_work_J=float(plastic[mask].sum()),
            fraction_of_total_plastic_work=float(plastic[mask].sum()/plastic.sum()),
            plastic_work_during_floor_active_steps_J=float(floor_work[mask].sum()),
            final_ep_max=float(ep[-1, mask].max()), final_ep_mean=float(ep[-1, mask].mean()),
            final_rate_mean_per_s=float(field['strain_rate'][-1, mask].mean()),
            max_net_displacement_m=float(np.linalg.norm(x[-1, mask]-x[0, mask], axis=1).max()),
            best_rigid_fit_residual=rigid_residual(x[0, mask], x[-1, mask]))
    result = dict(status='PASS_SAVED_30MS_ENERGY_AND_DEFORMATION_AUDIT_NOT_PHYSICAL_VALIDATION',
        physical_solver_evaluations=0, final_time_s=float(time_s[-1]), accepted_steps=281,
        all_energy_prefixes_reconstructed=True, old_25_25ms_branch_samples=0,
        energy=dict(cumulative_signed_raw_J=sums['global_raw_residual_J'], cumulative_absolute_step_raw_J=absolute,
            cumulative_registered_scale_J=scale, cumulative_external_work_J=Wext,
            raw_absolute_over_external_work=absolute/abs(Wext), raw_absolute_over_registered_scale=absolute/scale,
            registered_scale_over_external_work=scale/abs(Wext),
            interval_allowance_over_external_work=state['interval_energy']['interval_limit_J']/abs(Wext),
            maximum_prefix_allowance_usage=max_prefix_usage, material_local_absolute_defect_J=local_abs,
            material_local_signed_defect_J=local_signed, largest_field_energy_reconstruction_error_J=max_energy_error,
            largest_thermal_field_binding_error_J=max_thermal_error, final_deviatoric_storage_J=float(stored_dev[-1]),
            final_volumetric_storage_J=float(stored_vol[-1]), cumulative_plastic_work_J=float(plastic.sum())),
        mechanical_sums=sums, thermal_sums=summary['thermal_whole_interval'], groups=groups,
        floor=dict(rate_floor_per_s=p['rate_floor'], final_particle_count=int(floor[-1].sum()),
            cumulative_floor_active_plastic_work_J=float(floor_work.sum()),
            fraction_of_plastic_work_during_floor_active_steps=float(floor_work.sum()/plastic.sum())),
        final_state=summary['final'], source_and_input_sha256=hashes,
        interpretation=['Nonzero rigid-fit residual and free-particle plastic work show numerical nonrigid deformation, not merely prescribed boundary translation.',
            'Initially below-shoulder stock plastic work is not proof of through-aperture transport; crossing observer remains separate.',
            'The cumulative scale sums both work and storage/dissipation magnitudes, so its 5% threshold is not a 5% external-work error bound.',
            'No spatial convergence, low-rate law sensitivity, dynamic/QS comparison, or experimental validation is added.'],
        wall_seconds=time.perf_counter()-begin)
    OUTPUT.write_text(json.dumps(reader.clean(result), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(output=str(OUTPUT), energy=result['energy'], groups=groups, floor=result['floor'], wall_seconds=result['wall_seconds'])))


if __name__ == '__main__':
    main()

"""Reconstruct 431 saved accepted energy prefixes from actual particle fields."""
from pathlib import Path
import json
import math
import time
import numpy as np
import readback_discrete_process as reader
from afsd_backend.pysph.mechanics import bulk_modulus, volumetric_energy

ROOT = reader.ROOT
BASE = ROOT/'outputs/coupled_cost_20261004/batched/readback_pre_occlusion_39_375ms'
OUT = ROOT/'outputs/coupled_cost_20261004/energy_pre_occlusion_39_375ms_audit.json'


def main():
    begin = time.perf_counter(); reader.require(not OUT.exists(), 'refuse to overwrite energy audit')
    summary = reader.read(BASE/'summary.json'); trajectory = BASE/'accepted_trajectory.npz'
    reader.require(summary['accepted_steps'] == 431 and summary['samples'] == 432 and
        abs(summary['accepted_time_s']-.039375) <= 1e-15 and reader.sha(trajectory) == summary['trajectory_sha256'],
        'wrong accepted chain or changed trajectory')
    hashes = {path.relative_to(ROOT).as_posix(): reader.sha(path) for path in (Path(__file__), BASE/'summary.json', trajectory)}
    with np.load(trajectory, allow_pickle=False) as data:
        rho, S, T, x = (data[name] for name in ('density_kg_m3', 'deviatoric_stress_Pa', 'temperature_K', 'position_m'))
        m = data['particle_property__m'][0, :, 0]
        work = data['particle_property__plastic_work'][..., 0]
        component = data['particle_property__component_id'][0, :, 0]
        rate = data['particle_property__strain_rate'][..., 0]
        _, state = reader.restore(ROOT/summary['final']['checkpoint']); p = state['parameters']
        K = bulk_modulus(p['shear'], p['poisson'])
        storage = np.sum(m[None, :]/rho*np.sum(S*S, axis=(2, 3))/(4*p['shear'])+
            m[None, :]*volumetric_energy(rho, p['rho0'], K), axis=1)
        plastic = np.diff(work, axis=0)
        reader.require(np.all(plastic >= -1e-18), 'plastic work history decreased')
        heat = np.sum(m[None, :]*p['cp']*np.diff(T, axis=0), axis=1)
        groups = {}
        for label, mask in (('initially_above_shoulder_bar', (component == 1)&(x[0, :, 2] > p['layer'])),
                            ('initially_below_shoulder_bar', (component == 1)&(x[0, :, 2] < p['layer']))):
            cumulative = work[-1]-work[0]
            groups[label] = dict(particles=int(mask.sum()), plastic_work_J=float(cumulative[mask].sum()),
                fraction_of_all_plastic_work=float(cumulative[mask].sum()/cumulative.sum()))
        floor_fraction = float(np.sum(np.where(rate[1:] < p['rate_floor'], plastic, 0.))/plastic.sum())
    records, energy = [], []
    max_field_error = max_heat_error = max_prefix_usage = 0.
    for index, row in enumerate(summary['rows']):
        path = (ROOT/row['checkpoint']).parent/'energy.json'
        document = reader.read(path); e = document['mechanical']; ledger = document['interval']
        hashes[path.relative_to(ROOT).as_posix()] = reader.sha(path)
        delta_storage = storage[index+1]-storage[index]
        reconstructed_raw = e['prescribed_boundary_work_J']+e['tool_work_J']-(delta_storage+plastic[index].sum()+
            e['normal_storage_increment_J']+e['tangential_storage_increment_J']+e['friction_sliding_J'])
        discrepancy = max(abs(delta_storage-e['material_storage_increment_J']),
            abs(plastic[index].sum()-e['plastic_work_J']), abs(reconstructed_raw-e['global_raw_residual_J']))
        thermal_discrepancy = abs(heat[index]-document['thermal']['measured_thermal_increment_J'])
        reader.require(discrepancy < 1e-11 and thermal_discrepancy < 1e-12, 'saved energy does not bind to actual fields')
        max_field_error = max(max_field_error, discrepancy); max_heat_error = max(max_heat_error, thermal_discrepancy)
        records.append({name: e[name] for name in ('global_raw_residual_J', 'absolute_work_scale_J')})
        reader.validate_interval_energy_state(ledger)
        reader.require(ledger['records'] == records and ledger['steps'] == index+1, 'energy prefix reset, branch mixed, or gap')
        absolute = math.fsum(abs(record['global_raw_residual_J']) for record in records)
        reader.require(absolute <= ledger['interval_limit_J'], 'actual accepted energy prefix fails')
        max_prefix_usage = max(max_prefix_usage, absolute/ledger['interval_limit_J']); energy.append(e)
    def totals(rows):
        ext = math.fsum(row['prescribed_boundary_work_J']+row['tool_work_J'] for row in rows)
        absolute = math.fsum(abs(row['global_raw_residual_J']) for row in rows)
        scale = math.fsum(row['absolute_work_scale_J'] for row in rows)
        return dict(external_work_J=ext, signed_raw_J=math.fsum(row['global_raw_residual_J'] for row in rows),
            absolute_step_raw_J=absolute, registered_scale_J=scale, raw_absolute_over_external_work=absolute/abs(ext),
            registered_scale_over_external_work=scale/abs(ext), plastic_work_J=math.fsum(row['plastic_work_J'] for row in rows),
            material_local_absolute_defect_J=math.fsum(row['material_breakdown']['material_raw_energy_defect_absolute_sum_J'] for row in rows),
            material_local_signed_defect_J=math.fsum(row['material_breakdown']['material_raw_energy_defect_J'] for row in rows))
    result = dict(status='PASS_ALL_431_SAVED_ENERGY_PREFIXES_NOT_PHYSICAL_VALIDATION', physical_solver_evaluations=0,
        accepted_steps=431, accepted_time_s=summary['accepted_time_s'], same_parent_response_accepted_steps=32,
        whole_interval=totals(energy), accepted_two_ms_response_only=totals(energy[399:]),
        largest_field_energy_reconstruction_error_J=max_field_error, largest_thermal_field_binding_error_J=max_heat_error,
        maximum_prefix_allowance_usage=max_prefix_usage, final_interval_limit_J=state['interval_energy']['interval_limit_J'],
        final_allowance_over_external_work=state['interval_energy']['interval_limit_J']/abs(totals(energy)['external_work_J']),
        initial_material_location_groups=groups, plastic_work_fraction_during_rate_floor_active_steps=floor_fraction,
        original_four_ms_status='FAILED_NO_NEW_ACCEPTED_STATE', occlusion_model_corrected=False,
        source_and_input_sha256=hashes, wall_seconds=time.perf_counter()-begin,
        limitations=['The 5% cumulative work-plus-storage scale is not a 5% external-work error guarantee.',
            'A small raw energy defect does not validate spatial accuracy, low-rate material assumptions, or deposition.'])
    OUT.write_text(json.dumps(reader.clean(result), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

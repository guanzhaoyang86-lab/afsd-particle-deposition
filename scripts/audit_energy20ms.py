"""Independent read-only audit of the qualified 136-step, 20 ms branch.

No trial.evaluate, constitutive integration, contact solve, or time advance.
The yield-law evaluations below diagnose saved fields only.
"""
from pathlib import Path
import json
import math
import time
import numpy as np
import readback_discrete_process as reader
from afsd_backend.pysph.material import FKS, equivalent
from afsd_backend.pysph.mechanics import bulk_modulus, volumetric_energy
from afsd_backend.pysph.interval_acceptance import validate_interval_energy_state

ROOT = reader.ROOT
OLD = ROOT/'outputs/discrete_contact_process_20261004/v3'
NEW = ROOT/'outputs/coupled_cost_20261004/batched'
FINAL = NEW/'window_6_20ms/refinement/eighth'
OUTPUT = ROOT/'outputs/coupled_cost_20261004/energy20ms_audit.json'


def fit_rigid(first, last):
    a, b = first-first.mean(axis=0), last-last.mean(axis=0)
    left, _, right = np.linalg.svd(a.T@b)
    adjust = np.eye(3); adjust[-1, -1] = np.linalg.det(left@right)
    rotation = left@adjust@right
    residual = b-a@rotation
    return dict(rms_m=float(np.sqrt(np.mean(np.sum(residual**2, axis=1)))),
                maximum_m=float(np.linalg.norm(residual, axis=1).max()))


def main():
    begin = time.perf_counter()
    reader.require(not OUTPUT.exists(), 'refuse to overwrite audit')
    stages = [OLD/'quarters', OLD/'post_start/halves', OLD/'post_start/pilot_refinement',
              NEW/'halves', FINAL]
    initial, initial_state = reader.restore(OLD/'quarters/independent_initial.npz')
    x0 = reader.xyz(initial)
    initial_fields = reader.fields(initial)
    previous_pa, previous_state = initial, initial_state
    previous_path = OLD/'quarters/independent_initial.npz'
    records, rows, energy_rows, thermal_rows, sources = [], [], [], [], {}
    total_equivalent_strain = np.zeros(len(initial.x))
    plastic_work_by_floor_particle = np.zeros(len(initial.x))
    material_local_absolute = material_local_signed = 0.
    largest_field_bound_energy_error = largest_thermal_error = largest_yield_error = 0.
    maximum_prefix_usage = 0.
    for directory in stages:
        meta = reader.read(directory/'manifest.json')
        history = reader.read(directory/'accepted_history.json')
        sources[str((directory/'manifest.json').relative_to(ROOT))] = reader.sha(directory/'manifest.json')
        for path, digest in meta['source_sha256'].items():
            reader.require(reader.sha(ROOT/path) == digest, 'frozen source changed: '+path)
        if directory != stages[0]:
            reader.require(meta['parent_sha256'] == reader.sha(previous_path), 'wrong branch parent')
        for entry in history:
            path = directory/entry['checkpoint']
            reader.require(reader.sha(path) == entry['sha256'], 'checkpoint hash mismatch')
            sources[str(path.relative_to(ROOT))] = reader.sha(path)
            pa, state = reader.restore(path)
            p = state['parameters']; law = FKS(**{key: p[key] for key in FKS.__dataclass_fields__})
            dt = entry['dt_s']; data = reader.fields(pa); old = reader.fields(previous_pa)
            reader.require(np.array_equal(pa.particle_id, previous_pa.particle_id), 'raw ID order changed')
            reader.require(state['step'] == previous_state['step']+1, 'global step is discontinuous')
            reader.require(abs(state['time_s']-previous_state['time_s']-dt) < 1e-15, 'time is discontinuous')
            for key in ('particle_id', 'material_id', 'component_id', 'tracer_type', 'wholebar_tag', 'actuator', 'm'):
                reader.require(np.array_equal(data[key], initial_fields[key]), 'persistent quantity changed: '+key)
            validate_interval_energy_state(state['interval_energy'])
            document = reader.read(path.parent/'energy.json')
            e, thermal = document['mechanical'], document['thermal']
            reader.require(reader.clean(state['energy']) == e, 'saved mechanical ledger differs')
            records.append({key: e[key] for key in ('global_raw_residual_J', 'absolute_work_scale_J')})
            reader.require(state['interval_energy']['records'] == records, 'energy prefix reset or branch concatenation')
            reader.require(state['interval_energy'] == document['interval'], 'interval JSON differs')
            raw, scale = e['global_raw_residual_J'], e['absolute_work_scale_J']
            raw_abs = math.fsum(abs(row['global_raw_residual_J']) for row in records)
            registered_limit = 1e-7+.05*math.fsum(row['absolute_work_scale_J'] for row in records)
            reader.require(raw_abs <= registered_limit, 'absolute-prefix gate fails')
            maximum_prefix_usage = max(maximum_prefix_usage, raw_abs/registered_limit)
            mass, rho, rho0 = pa.m, pa.rho, previous_pa.rho
            S, S0 = pa.S.reshape(-1, 3, 3), previous_pa.S.reshape(-1, 3, 3)
            K = bulk_modulus(p['shear'], p['poisson'])
            deltaU = float(np.sum(mass/rho*np.sum(S*S, axis=(1, 2))/(4*p['shear'])-
                                 mass/rho0*np.sum(S0*S0, axis=(1, 2))/(4*p['shear']))+
                           np.sum(mass*(volumetric_energy(rho, p['rho0'], K)-volumetric_energy(rho0, p['rho0'], K))))
            wp = pa.plastic_work-previous_pa.plastic_work
            reader.require(np.all(wp >= -1e-18), 'plastic work accumulator decreased')
            raw_from_fields = e['prescribed_boundary_work_J']+e['tool_work_J']-(deltaU+float(wp.sum())+
                e['normal_storage_increment_J']+e['tangential_storage_increment_J']+e['friction_sliding_J'])
            field_error = max(abs(deltaU-e['material_storage_increment_J']),
                              abs(float(wp.sum())-e['plastic_work_J']), abs(raw_from_fields-raw))
            largest_field_bound_energy_error = max(largest_field_bound_energy_error, field_error)
            reader.require(field_error < 1e-11, 'energy is not bound to actual saved fields')
            measured_heat = float(np.sum(mass*p['cp']*(pa.temperature-previous_pa.temperature)))
            largest_thermal_error = max(largest_thermal_error, abs(measured_heat-thermal['measured_thermal_increment_J']))
            reader.require(abs(measured_heat-thermal['measured_thermal_increment_J']) < 1e-12, 'thermal field ledger differs')
            floor = pa.strain_rate < law.rate_floor
            plastic_work_by_floor_particle += np.where(floor, wp, 0.)
            total_equivalent_strain += pa.strain_rate*dt
            active = pa.ep_rate > 1e-12
            q = equivalent(S)
            yield_error = q-law.yield_stress(pa.ep, pa.strain_rate, .5*(previous_pa.temperature+pa.temperature))
            largest_yield_error = max(largest_yield_error, float(np.max(np.abs(yield_error[active]))) if active.any() else 0.)
            breakdown = e['material_breakdown']
            material_local_absolute += breakdown['material_raw_energy_defect_absolute_sum_J']
            material_local_signed += breakdown['material_raw_energy_defect_J']
            energy_rows.append(e); thermal_rows.append(thermal)
            rows.append(dict(step=state['step'], time_s=state['time_s'], dt_s=dt,
                raw_energy_J=raw, registered_step_scale_J=scale,
                raw_to_step_scale=abs(raw)/scale if scale else 0.,
                floor_particle_count=int(floor.sum()), floor_plastic_work_J=float(wp[floor].sum()),
                plastic_work_J=float(wp.sum()), rate_mean=float(pa.strain_rate.mean()),
                rate_max=float(pa.strain_rate.max()), ep_max=float(pa.ep.max()), q_max_Pa=float(q.max()),
                temperature_max_K=float(pa.temperature.max())))
            previous_pa, previous_state, previous_path = pa, state, path
    final, state = previous_pa, previous_state
    reader.require(state['step'] == len(records) == 136, 'wrong final branch length')
    reader.require(abs(state['time_s']-.020) < 1e-15, 'wrong final time')
    numeric_keys = [key for key, value in energy_rows[0].items()
                    if isinstance(value, (int, float)) and not isinstance(value, bool)]
    sums = {key: math.fsum(row[key] for row in energy_rows) for key in numeric_keys}
    # This summed per-step limit is explicitly NOT an acceptance allowance.
    sums.pop('registered_energy_limit_J', None)
    thermal_sums = {key: math.fsum(row[key] for row in thermal_rows) for key in thermal_rows[0]}
    Wext = sums['prescribed_boundary_work_J']+sums['tool_work_J']
    signed = math.fsum(row['global_raw_residual_J'] for row in records)
    absolute = math.fsum(abs(row['global_raw_residual_J']) for row in records)
    scale = math.fsum(row['absolute_work_scale_J'] for row in records)
    q = equivalent(final.S.reshape(-1, 3, 3)); x = reader.xyz(final)
    velocity = np.column_stack((final.u, final.v, final.w))
    floor = final.strain_rate < law.rate_floor
    free = state['free']; component = final.component_id
    groups = dict(free_bar=(component == 1)&~final.actuator.astype(bool),
        top_bar_z_actuated=final.actuator.astype(bool), free_substrate=(component == 0)&free.any(axis=1),
        clamped_substrate=(component == 0)&~free.any(axis=1),
        initially_above_shoulder_bar=(component == 1)&(initial.z > p['layer']),
        initially_below_shoulder_bar=(component == 1)&(initial.z < p['layer']))
    group_metrics = {}
    for name, mask in groups.items():
        group_metrics[name] = dict(particles=int(mask.sum()), floor_count=int(floor[mask].sum()),
            cumulative_plastic_work_J=float(final.plastic_work[mask].sum()),
            fraction_of_total_plastic_work=float(final.plastic_work[mask].sum()/final.plastic_work.sum()),
            cumulative_floor_active_plastic_work_J=float(plastic_work_by_floor_particle[mask].sum()),
            ep_max=float(final.ep[mask].max()), ep_mean=float(final.ep[mask].mean()),
            rate_mean_per_s=float(final.strain_rate[mask].mean()), rate_max_per_s=float(final.strain_rate[mask].max()),
            max_displacement_m=float(np.linalg.norm(x[mask]-x0[mask], axis=1).max()),
            max_downward_displacement_m=float((x0[mask, 2]-x[mask, 2]).max()),
            best_rigid_fit_residual=fit_rigid(x0[mask], x[mask]))
    # This is only a diagnostic algebraic extrapolation, never a model update.
    tstar=(final.temperature-law.room)/(law.melt-law.room)
    unregularized_lambda=1.+law.b1*tstar**law.b2*np.log(final.strain_rate/law.reference_rate)
    K=bulk_modulus(p['shear'], p['poisson'])
    stored_dev=float(np.sum(final.m/final.rho*np.sum(final.S.reshape(-1,3,3)**2,axis=(1,2))/(4*p['shear'])))
    stored_vol=float(np.sum(final.m*volumetric_energy(final.rho,p['rho0'],K)))
    kinetic_proxy=float(np.sum(.5*final.m*np.sum(velocity**2,axis=1)))
    final_dir=FINAL/'step_112';journal=reader.read(final_dir/'broyden/last_UNCOMMITTED_candidate.json')
    forces=np.asarray(journal['contact']['force_N'])
    previous_last,_=reader.restore(FINAL/'step_111/accepted.npz')
    vprev=np.column_stack((previous_last.u,previous_last.v,previous_last.w))
    inertia_proxy=final.m[:,None]*(velocity-vprev)/rows[-1]['dt_s']
    eligible=groups['initially_above_shoulder_bar']
    final_geometry=reader.read(final_dir/'geometry_audit.json')
    result=dict(status='PASS_SAVED_ENERGY_PREFIX_AUDIT_NOT_PHYSICAL_VALIDATION',
        physical_evaluations=0, classification=state['classification'], final_time_s=state['time_s'],
        accepted_steps=136, mass_kg=float(final.m.sum()), all_prefixes_pass=True,
        interval_energy=state['interval_energy'],
        energy_scale_assessment=dict(cumulative_signed_raw_J=signed, cumulative_absolute_step_raw_J=absolute,
            cumulative_registered_scale_J=scale, cumulative_external_work_J=Wext,
            raw_absolute_over_registered_scale=absolute/scale, raw_absolute_over_external_work=absolute/abs(Wext),
            registered_scale_over_external_work=scale/abs(Wext),
            interval_limit_over_external_work=state['interval_energy']['interval_limit_J']/abs(Wext),
            maximum_accepted_prefix_gate_usage=maximum_prefix_usage,
            material_local_absolute_defect_sum_J=material_local_absolute,
            material_local_signed_defect_sum_J=material_local_signed,
            largest_saved_field_energy_reconstruction_error_J=largest_field_bound_energy_error,
            largest_thermal_field_binding_error_J=largest_thermal_error),
        mechanical_sums=sums, thermal_sums=thermal_sums,
        rate_floor_assessment=dict(rate_floor_per_s=law.rate_floor, final_count=int(floor.sum()),
            final_fraction=float(floor.mean()), final_rate_min_per_s=float(final.strain_rate.min()),
            final_rate_max_per_s=float(final.strain_rate.max()),
            cumulative_plastic_work_in_floor_active_steps_J=float(plastic_work_by_floor_particle.sum()),
            cumulative_plastic_work_fraction_under_floor=float(plastic_work_by_floor_particle.sum()/final.plastic_work.sum()),
            final_plastic_work_fraction_under_floor=rows[-1]['floor_plastic_work_J']/rows[-1]['plastic_work_J'],
            final_unregularized_lambda_nonpositive_count=int(np.count_nonzero(unregularized_lambda <= 0)),
            largest_active_midtemperature_yield_consistency_error_Pa=largest_yield_error,
            actual_low_rate_sensitivity_process_run='NOT_RUN'),
        groups=group_metrics,
        final_state=dict(temperature_min_K=float(final.temperature.min()), temperature_max_K=float(final.temperature.max()),
            ep_max=float(final.ep.max()), q_max_Pa=float(q.max()), rho_reference_ratio_min=float((final.rho/p['rho0']).min()),
            rho_reference_ratio_max=float((final.rho/p['rho0']).max()), stored_deviatoric_energy_J=stored_dev,
            stored_volumetric_energy_J=stored_vol, stored_elastic_over_plastic_work=(stored_dev+stored_vol)/float(final.plastic_work.sum()),
            max_speed_m_s=float(np.linalg.norm(velocity,axis=1).max()),
            crossing_mass_kg=state['crossing_mass_kg'], crossing_ids=state['crossing_ids'],
            nearest_eligible_bar_center_above_shoulder_m=float((x[eligible,2]-p['layer']).min()),
            geometry_audit=final_geometry),
        quasistatic_scale_diagnostics=dict(mean_step_velocity_kinetic_energy_proxy_J=kinetic_proxy,
            kinetic_proxy_over_external_work=kinetic_proxy/abs(Wext),
            rho_Umax_squared_over_min_positive_yield=float(p['rho0']*np.max(np.sum(velocity**2,axis=1))/law.yield_stress(final.ep,final.strain_rate,final.temperature).min()),
            finite_difference_inertial_force_max_N=float(np.linalg.norm(inertia_proxy,axis=1).max()),
            material_contact_force_max_N=float(np.linalg.norm(forces,axis=1).max()),
            inertia_proxy_force_norm_over_contact_force_norm=float(np.linalg.norm(inertia_proxy)/np.linalg.norm(forces)),
            caveat='Stored velocities are step averages; these ratios are diagnostics, not a dynamic/quasistatic validation.'),
        tests_not_performed=['spatial convergence', 'inertial versus quasistatic comparison', 'FKS low-rate sensitivity',
            'prebond/interface and clamp sensitivity', 'experimental comparison', 'finite-strain cyclic constitutive energy validation'],
        rows=rows, source_input_sha256=sources, audit_script_sha256=reader.sha(Path(__file__)),
        wall_seconds=time.perf_counter()-begin)
    OUTPUT.write_text(json.dumps(reader.clean(result),indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(reader.clean({key: result[key] for key in ('status','final_time_s','accepted_steps',
          'energy_scale_assessment','rate_floor_assessment','final_state','quasistatic_scale_diagnostics','wall_seconds')}),indent=2))


if __name__=='__main__':
    main()

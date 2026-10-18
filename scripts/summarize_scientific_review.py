"""Read-only audit of failed review trials; never runs or modifies physics.

Writes only outputs/review/sensitivity_summary.json. No missing terminal state
is reconstructed, and different failure times are never called convergence.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud

CASES = ('coarse_dt', 'dt_half', 'dt_quarter', 'space_half')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def flatten_numbers(value, prefix=''):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from flatten_numbers(item, prefix+'.'+key if prefix else key)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from flatten_numbers(item, f'{prefix}[{index}]')
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield prefix, value


def parameter_leaves(value, prefix=''):
    if isinstance(value, dict) and 'value' in value and 'unit' in value:
        return {prefix: value}
    result = {}
    if isinstance(value, dict):
        for key, item in value.items():
            result.update(parameter_leaves(item, prefix+'.'+key if prefix else key))
    return result


def numerical_extrema(records):
    rows = {}
    for record in records:
        for key, value in flatten_numbers(record):
            rows.setdefault(key, []).append((value, record.get('step'), record.get('time_s')))
    result = {}
    for key, values in rows.items():
        finite = [item for item in values if math.isfinite(item[0])]
        low = min(finite, default=None, key=lambda item: item[0])
        high = max(finite, default=None, key=lambda item: item[0])
        result[key] = dict(min=low[0] if low else None, max=high[0] if high else None,
                           min_step=low[1] if low else None, max_step=high[1] if high else None,
                           samples=len(values), nonfinite_samples=len(values)-len(finite))
    return result


def energy_gate(records):
    required = ('total_J', 'energy_residual_J', 'energy_residual_relative',
                'wall_work_J', 'actuator_work_J', 'wall_heat_out_J')
    if not records or any(any(key not in row for key in required) for row in records):
        return dict(status='BLOCKED', maximum_allowance_ratio=None, reason='missing energy ledger fields')
    initial = [row['total_J']-row['energy_residual_J']-row['wall_work_J']-
               row['actuator_work_J']+row['wall_heat_out_J'] for row in records]
    spread = max(initial)-min(initial)
    if spread > 1e-12*max(abs(initial[0]), 1e-15)+1e-18:
        return dict(status='BLOCKED', maximum_allowance_ratio=None, reason='inconsistent recovered initial energy', initial_energy_spread_J=spread)
    details = []
    for row, energy0 in zip(records, initial):
        exchanged = abs(row['wall_work_J'])+abs(row['actuator_work_J'])+abs(row['wall_heat_out_J'])
        scale = max(abs(energy0), exchanged, 1e-15)
        relative = row['energy_residual_J']/scale
        if not math.isclose(relative, row['energy_residual_relative'], rel_tol=1e-10, abs_tol=1e-14):
            return dict(status='BLOCKED', maximum_allowance_ratio=None, reason='stored relative residual inconsistent with documented scale')
        allowance = .01*scale+1e-7
        details.append(dict(step=row['step'], time_s=row['time_s'], residual_J=row['energy_residual_J'],
                            exchanged_absolute_J=exchanged, denominator_scale_J=scale,
                            N08_allowance_J=allowance, signed_allowance_ratio=row['energy_residual_J']/allowance,
                            absolute_allowance_ratio=abs(row['energy_residual_J'])/allowance))
    worst = max(details, key=lambda row: row['absolute_allowance_ratio'])
    return dict(status='PASS_RECORDED_PREFIX_ONLY' if worst['absolute_allowance_ratio'] <= 1 else 'FAIL',
                maximum_allowance_ratio=worst['absolute_allowance_ratio'], worst=worst,
                initial_energy_recovered_J=initial[0], initial_energy_spread_J=spread,
                method='E0=total-residual-wall_work-actuator_work+wall_heat_out; checked constant and checked residual/scale against archived relative residual',
                limitation='Algebraic audit of recorded ledger, not independent verification of constitutive energy; failed partial step is absent.',
                samples=details)


def initial_geometry(vtk, p):
    scalars, points = vtk['scalars'], vtk['points']
    masses, densities = scalars['mass_kg'], scalars['density_kg_per_m3']
    volumes = [mass/rho for mass, rho in zip(masses, densities)]
    components = {}
    for component in sorted(set(scalars['component_id'])):
        ids = [i for i, value in enumerate(scalars['component_id']) if value == component]
        components[str(component)] = dict(count=len(ids), volume_m3=math.fsum(volumes[i] for i in ids),
            initial_full_cell_bounds_m=[[min(points[i][a] for i in ids)-p['dx']/2,
                                        max(points[i][a] for i in ids)+p['dx']/2] for a in range(3)])
    cap = [i for i, value in enumerate(scalars['actuator']) if value == 1]
    h = scalars['smoothing_length_m']
    return dict(particle_count=len(points), initial_mass_kg=math.fsum(masses), components=components,
                actuator_particle_count=len(cap), actuator_volume_m3=math.fsum(volumes[i] for i in cap),
                actuator_configured_depth_m=p['actuator_depth'],
                actuator_full_cell_z_bounds_m=[min(points[i][2] for i in cap)-p['dx']/2, max(points[i][2] for i in cap)+p['dx']/2] if cap else None,
                kn_N_m=p['penalty_factor']*p['shear']*p['dx']**2/p['contact_reference_spacing'],
                contact_area_per_particle_m2=p['dx']**2,
                kn_per_area_N_m3=p['penalty_factor']*p['shear']/p['contact_reference_spacing'],
                contact_reference_spacing_m=p['contact_reference_spacing'],
                h_min_m=min(h), h_max_m=max(h), expected_h_m=p['h_ratio']*p['dx'],
                h_is_normal_discretization_change=True, layer_spacings=p['layer']/p['dx'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review-root', type=Path, default=ROOT/'outputs/review')
    args = parser.parse_args()
    review_root = args.review_root.resolve()
    inputs = {}
    def remember(path):
        path = path.resolve()
        inputs[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = sha(path)
        return path
    def read_json(path):
        return json.loads(remember(path).read_text(encoding='utf-8-sig'))
    remember(Path(__file__))
    remember(ROOT/'src/afsd_scaffold/vtk_legacy.py')
    solver_path = remember(ROOT/'src/afsd_backend/pysph/solver.py')
    contract_path = remember(ROOT/'logs/review/contract_v1_pretest.md')
    configs = {}; cases = {}
    for name in CASES:
        directory = review_root/name
        performance = read_json(directory/'performance.json')
        config = read_json(directory/'resolved_config.json'); configs[name] = parameter_leaves(config['parameters'])
        manifest = read_json(directory/'manifest.json')
        sources = read_json(directory/'source_hashes.json')
        series = read_json(directory/'series.json')
        run_log = remember(directory/'run.log').read_text(encoding='utf-8-sig')
        diagnostic_path = remember(directory/'diagnostics.jsonl')
        records = [json.loads(line) for line in diagnostic_path.read_text(encoding='utf-8-sig').splitlines() if line.strip()]
        initial_frame = next(frame for frame in series['frames'] if frame['step'] == 0)
        vtk_path = remember(directory/initial_frame['file'])
        if sha(vtk_path) != initial_frame['sha256']: raise ValueError('initial VTK hash mismatch: '+name)
        vtk = read_point_cloud(vtk_path)
        p = {key: value['value'] for key, value in config['parameters']['debug'].items()}
        p['layer'] = config['parameters']['process']['layer_height']['value']
        geometry = initial_geometry(vtk, p)
        target = performance.get('requested_end_time_s')
        actual = performance['physical_seconds']
        reached = target is not None and math.isclose(actual, target, rel_tol=1e-12, abs_tol=1e-18)
        final = records[-1] if records else None
        extrema = numerical_extrema(records)
        mass_error = max((abs(row['mass_kg']-geometry['initial_mass_kg'])/geometry['initial_mass_kg'] for row in records), default=None)
        budget = energy_gate(records)
        logs_exit = re.findall(r'(?:exit_code|EXIT_CODE)\s*[=:]\s*(-?\d+)', run_log)
        failure = performance.get('error', '').strip().splitlines()[-1] if performance.get('error') else None
        cases[name] = dict(classification=performance['classification'], status=performance['status'],
            manifest_status=manifest['status'], quantitative_validation=False,
            exact_command=manifest['provenance']['command'], exit_code=int(logs_exit[-1]) if logs_exit else None,
            exit_code_note='No numerical exit code stored in run.log; do not infer one from FAILED.' if not logs_exit else None,
            target_time_s=target, target_matches_frozen_2us=target == 2e-6, reached_target=reached,
            completed_physical_time_s=actual, completed_steps=performance['completed_steps'],
            failed_step_attempt=performance['completed_steps']+1 if performance['status'] == 'FAILED' else None,
            failure=failure, negative_midpoint_friction_failure='negative midpoint friction dissipation' in (failure or ''),
            particle_count=performance['particle_count'], dt_min_s=performance.get('dt_min_s'), dt_max_s=performance.get('dt_max_s'),
            CPU_host_wall_seconds=performance['wall_seconds'], wall_cap_s=performance['wall_limit_seconds'],
            CPU_process_time_seconds=None, CPU_process_time_note='Not measured; recorded wall time is not CPU process time.',
            median_completed_step_wall_seconds=performance.get('physics_step_median_seconds'),
            completed_steps_per_wall_second=performance['end_to_end_steps_per_second'],
            peak_process_working_set_bytes=performance['peak_process_working_set_bytes'],
            diagnostics_count=len(records), diagnostic_last_matches_completed=bool(final and final['step'] == performance['completed_steps'] and final['time_s'] == actual),
            saved_checkpoint_step=performance['final_checkpoint_step'], saved_frame_steps=[frame['step'] for frame in series['frames']],
            final_completed_step_diagnostics=final, all_numeric_diagnostic_extrema=extrema,
            initial_geometry=geometry, source_solver_matches_archived=sources.get(str(solver_path.relative_to(ROOT))) == sha(solver_path),
            N08_recorded_prefix=dict(energy=budget, mass_relative_error_max=mass_error,
                mass_pass=mass_error is not None and mass_error <= 1e-12,
                max_penetration_over_dx=max((row['max_penetration_m']/p['dx'] for row in records), default=None),
                all_numeric_finite=all(value['nonfinite_samples'] == 0 for value in extrema.values()),
                alarms=sorted({alarm for row in records for alarm in row.get('alarms', [])}),
                complete_run_acceptance='BLOCKED' if not reached or performance['status'] != 'COMPLETE' else 'REQUIRES_FULL_CONTRACT_CHECK'),
            unavailable_terminal_fields=['volume_average_equivalent_plastic_strain', 'stress_extrema', 'terminal_particle_state_at_failure'],
            cost_projection=dict(status='EXTRAPOLATION_ONLY_NOT_A_RUN', projection=performance.get('paper_cost_projection'),
                warning='From a short failed prefix; optimistic host N-linear/dt-spacing scaling only. Nominal 1mm resolves layer by 1.5 cells and is not an admissible production case.'))
    base = configs['coarse_dt']; comparisons = {}
    for name in CASES[1:]:
        other = configs[name]
        value_diffs = {}; provenance_diffs = {}
        for key in sorted(set(base)|set(other)):
            a, b = base.get(key), other.get(key)
            if a is None or b is None or (a['value'], a['unit']) != (b['value'], b['unit']):
                value_diffs[key] = {'coarse': a, 'candidate': b}
            elif a != b:
                provenance_diffs[key] = {'coarse': a, 'candidate': b}
        unexpected = sorted(set(value_diffs)-{'debug.dx', 'debug.cfl'})
        g0, g1 = cases['coarse_dt']['initial_geometry'], cases[name]['initial_geometry']
        same_boundary = all(math.isclose(g0[key], g1[key], rel_tol=1e-12, abs_tol=0.) for key in ('actuator_volume_m3', 'kn_per_area_N_m3', 'actuator_configured_depth_m', 'contact_reference_spacing_m'))
        same_bounds = all(all(math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15) for aa, bb in zip(g0['components'][component]['initial_full_cell_bounds_m'], g1['components'][component]['initial_full_cell_bounds_m']) for a, b in zip(aa, bb)) for component in g0['components'])
        comparisons[name] = dict(parameter_value_unit_differences=value_diffs, provenance_only_differences=provenance_diffs,
            unexpected_physics_or_boundary_differences=unexpected, same_physical_parameter_values_except_dx_cfl=not unexpected,
            same_initial_full_cell_geometric_bounds=same_bounds, same_physical_actuator_volume_and_penalty_traction=same_boundary,
            allowable_h_change='h=h_ratio*dx, h_ratio unchanged', terminal_comparison_status='BLOCKED',
            metric_differences=None, reason='All trials stopped with negative midpoint friction before the frozen common endpoint; different failure times are not convergence data.')
    summary = dict(schema_version='1.0', created_utc=datetime.now(timezone.utc).isoformat(),
        classification='APPROXIMATE_FAILED_TRIAL_AUDIT', contract_version='v1', contract_sha256=sha(contract_path),
        run_simulation=False, target_time_s=2e-6, production_authorized=False,
        N09_time_sensitivity_status='BLOCKED', N10_space_sensitivity_status='BLOCKED',
        checks=dict(all_targets_frozen=all(case['target_matches_frozen_2us'] for case in cases.values()),
                    all_failures_negative_midpoint_friction=all(case['negative_midpoint_friction_failure'] for case in cases.values()),
                    any_target_reached=any(case['reached_target'] for case in cases.values())),
        cases=cases, comparison_preconditions=comparisons, input_sha256=inputs,
        limitations=['No final-time field differences or error percentages computed.',
                     'Partial failed kick is not a valid state and is absent from completed-step diagnostics.',
                     'Only initial VTK/checkpoint was saved; short-prefix diagnostics do not replace terminal field output.',
                     'No Fig8 area ratio, reference CSV, calibration, or qualitative validation is inferred.'])
    target_path = review_root/'sensitivity_summary.json'
    if target_path.exists(): raise FileExistsError('Do not overwrite existing review evidence: '+str(target_path))
    target_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    print('SUMMARY_ONLY_NO_PHYSICS_RUN '+str(target_path))
    for name, case in cases.items():
        print(json.dumps(dict(case=name, status=case['status'], N=case['particle_count'], steps=case['completed_steps'],
            physical_s=case['completed_physical_time_s'], target_s=case['target_time_s'], wall_s=case['CPU_host_wall_seconds'],
            peak_RAM_bytes=case['peak_process_working_set_bytes'],
            max_N08_allowance_ratio=case['N08_recorded_prefix']['energy']['maximum_allowance_ratio'],
            actuator_volume_m3=case['initial_geometry']['actuator_volume_m3'], kn_per_area_N_m3=case['initial_geometry']['kn_per_area_N_m3'],
            failure=case['failure'])))
    print('N09=BLOCKED N10=BLOCKED; no unequal-time convergence comparison.')


if __name__ == '__main__':
    main()

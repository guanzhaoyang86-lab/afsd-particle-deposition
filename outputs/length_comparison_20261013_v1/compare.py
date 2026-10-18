"""Read qualified, equal-time independent length controls; never run physics.

Cross-case cohorts are defined by each independent initial geometry. IDs are
persistent within each case, but are intentionally NOT matched across cases.
Only a new analysis JSON is written. No run, checkpoint or source is changed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import run_afsd_deposition_stage_remote_20261007 as runner
from afsd_deposition_progress_20261010 import frame_from_saved, summarize_progress


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def xyz(pa):
    return np.column_stack((pa.x, pa.y, pa.z))


def qualified_case(stage):
    stage = (ROOT / stage).resolve()
    require(stage.is_relative_to(ROOT / 'outputs'), 'stage must be a project output')
    meta = runner.checked_receipt(stage, {runner.STAGE_PASS, runner.PREFIX_PASS,
                                           runner.SHOULDER_PASS})
    cost_dir = ROOT / meta['cost_directory']
    cost, target, pa, state, references, parent = runner.continuation_parent(
        cost_dir, stage, meta['source_sha256'])
    initial_pa, initial_state = runner.load(ROOT / cost['initial_checkpoint'])
    require(initial_state['time_s'] == 0. and initial_state['step'] == 0
            and initial_state['parent_checkpoint'] is None,
            'independent zero-time initial checkpoint required')
    screen_path = ROOT / meta['screened_history'][-1]['screen_path']
    screen = runner.read(screen_path)
    coarse_path = screen_path.parent / 'coarse' / 'NOT_RESTARTABLE_candidate.npz'
    require(sha(coarse_path) == screen['coarse_checkpoint_sha256'],
            'coarse state does not match the selected pair')
    coarse_pa, coarse_state = runner.load(coarse_path)
    require(abs(state['time_s'] - coarse_state['time_s']) <= 1e-13,
            'coarse/fine endpoint times differ')
    for current in (pa, coarse_pa):
        require(np.array_equal(initial_pa.particle_id, current.particle_id)
                and np.array_equal(initial_pa.m, current.m)
                and np.array_equal(initial_pa.component_id, current.component_id),
                'within-case persistent identity/order changed')
    return dict(stage=stage, meta=meta, cost=cost, initial_pa=initial_pa,
        initial_state=initial_state, pa=pa, state=state, checkpoint=target,
        coarse_pa=coarse_pa, coarse_state=coarse_state, coarse_checkpoint=coarse_path,
        pair=screen, reference_count=len(references))


def group_metrics(initial, current, mask):
    count = int(mask.sum())
    if not count:
        return dict(particle_count=0)
    start, final = xyz(initial)[mask], xyz(current)[mask]
    weights = np.asarray(initial.m)[mask]
    weights = weights / weights.sum()
    displacement = final - start
    def mean(value):
        return float(np.sum(value * weights))
    def rms_radius(points):
        center = np.sum(points[:, :2] * weights[:, None], axis=0)
        return float(np.sqrt(mean(np.sum((points[:, :2] - center)**2, axis=1))))
    initial_r = np.linalg.norm(start[:, :2], axis=1)
    final_r = np.linalg.norm(final[:, :2], axis=1)
    return dict(particle_count=count, mass_kg=float(initial.m[mask].sum()),
        initial_z_min_m=float(start[:, 2].min()), initial_z_max_m=float(start[:, 2].max()),
        mean_dz_m=mean(displacement[:, 2]), minimum_dz_m=float(displacement[:, 2].min()),
        maximum_dz_m=float(displacement[:, 2].max()),
        mean_axis_radius_change_m=mean(final_r - initial_r),
        centroid_rms_radius_m=rms_radius(final),
        centroid_rms_radius_change_m=rms_radius(final) - rms_radius(start),
        mean_ep=mean(np.asarray(current.ep)[mask]), maximum_ep=float(current.ep[mask].max()),
        mean_temperature_K=mean(np.asarray(current.temperature)[mask]))


def endpoint_metrics(case, branch):
    initial = case['initial_pa']
    pa, state, checkpoint = ((case['pa'], case['state'], case['checkpoint'])
        if branch == 'fine' else
        (case['coarse_pa'], case['coarse_state'], case['coarse_checkpoint']))
    p = state['parameters']
    bar = np.asarray(initial.component_id) == 1
    cap = np.asarray(initial.actuator) != 0
    eligible = np.isin(initial.particle_id,
                       case['initial_state']['bore_observation']['eligible_ids'])
    groups = dict(actuator_cap=cap,
        initial_bore=bar & (initial.z >= p['layer']) & (initial.z <= p['tool_top']),
        initial_eligible_near_outlet=eligible & (initial.z >= p['layer'])
            & (initial.z <= p['layer'] + p['dx']),
        initial_lower_material=bar & (initial.z < p['layer']),
        initial_free_bar_above_tool=bar & ~cap & (initial.z > p['tool_top']))
    layers = []
    for level in np.unique(np.round(np.asarray(initial.z)[bar], 12)):
        mask = bar & np.isclose(initial.z, level, atol=1e-12, rtol=0.)
        layers.append(dict(initial_z_m=float(level),
            prescribed_z_particle_count=int(np.count_nonzero(~case['initial_state']['free'][mask, 2])),
            metrics=group_metrics(initial, pa, mask)))
    saved_observer = runner.read(checkpoint.parent / 'observer.json')
    progress = summarize_progress([frame_from_saved(pa, state, observer=saved_observer)], p)['final']
    journal = checkpoint.parent / 'FULL_UNCOMMITTED_candidate.npz'
    with np.load(journal, allow_pickle=False) as data:
        require(np.array_equal(data['particle_id'], pa.particle_id), 'journal IDs differ')
        require(np.array_equal(data['position_m'], xyz(pa)), 'journal endpoint positions differ')
        force = np.asarray(data['force_N'])
        require(force.shape == xyz(pa).shape and np.isfinite(force).all(), 'invalid journal forces')
        reaction_on_actuator_z_N = float(force[cap, 2].sum())
        actuator_on_material_z_N = -reaction_on_actuator_z_N
    program = runner.model.validate_state(state)
    commanded_feed = float(program.integral(state['time_s'])[1])
    cap_metric = group_metrics(initial, pa, cap)
    require(abs(cap_metric['mean_dz_m'] + commanded_feed) < 1e-12,
            'actual actuator displacement does not match command')
    top_free = max(row['initial_z_m'] for row in layers
                   if row['prescribed_z_particle_count'] == 0)
    bore_top_layer = min(row['initial_z_m'] for row in layers if row['initial_z_m'] > p['tool_top'])
    def layer_dz(level):
        return next(row['metrics']['mean_dz_m'] for row in layers if row['initial_z_m'] == level)
    return dict(time_s=float(state['time_s']), step=int(state['step']),
        actuator=dict(commanded_downward_feed_m=commanded_feed,
            actual_mean_downward_feed_m=-cap_metric['mean_dz_m'],
            reaction_on_actuator_z_N=reaction_on_actuator_z_N,
            actuator_on_material_z_N=actuator_on_material_z_N,
            force_convention='Journal total internal plus contact force on prescribed-z cap; external actuator load on material is its negative.'),
        groups={name: group_metrics(initial, pa, mask) for name, mask in groups.items()},
        initial_z_layers=layers,
        free_upper_segment=dict(initial_lower_z_m=bore_top_layer, initial_upper_z_m=top_free,
            initial_length_m=top_free-bore_top_layer,
            shortening_m=layer_dz(bore_top_layer)-layer_dz(top_free)),
        deposition=progress, endpoint_force_journal=runner.rel(journal),
        endpoint_force_journal_sha256=sha(journal))


FIELDS = ('mean_dz_m', 'mean_axis_radius_change_m', 'centroid_rms_radius_change_m',
          'mean_ep', 'maximum_ep', 'mean_temperature_K')


def comparison(short_fine, long_fine, short_coarse, long_coarse):
    """Signed short-minus-long effect versus same-parent local pair differences."""
    effect = float(short_fine-long_fine)
    short_error = abs(float(short_fine-short_coarse))
    long_error = abs(float(long_fine-long_coarse))
    total = short_error+long_error
    return dict(short_minus_long_fine=effect, short_abs_fine_minus_coarse=short_error,
        long_abs_fine_minus_coarse=long_error, sum_local_pair_differences=total,
        effect_over_sum_local_pair_differences=abs(effect)/total if total > 0 else None,
        effect_exceeds_sum_local_pair_differences=bool(abs(effect) > total),
        interpretation='Local final-macro indicator only; not a cumulative error bound, significance test or spatial convergence result.')


def program_and_configuration_scope(long, short, archived_long_reference):
    configurations = [runner.read(ROOT / case['meta']['configuration']) for case in (long, short)]
    for key in ('schema_version', 'classification', 'base_parameters', 'base_parameters_sha256',
                'geometry', 'thermal_initialization', 'interface_initialization', 'feed_mode', 'bar_rotation'):
        require(configurations[0].get(key) == configurations[1].get(key),
                'initial model configuration differs beyond length: ' + key)
    programs = [runner.model.validate_state(case['state']) for case in (long, short)]
    time_s = float(long['state']['time_s'])
    if not archived_long_reference:
        require(long['state']['process_program'] == short['state']['process_program'],
                'common process command required')
        return dict(mode='NEW_INDEPENDENT_PAIRED_CONTROLS', whole_program_identical=True)
    require(long['stage'].is_relative_to(ROOT / 'outputs/deposition_campaign_20261007')
            and long['meta']['configuration'] == 'configs/deposition_fill_then_track_20261007.json'
            and long['meta']['cost_directory'] == 'outputs/deposition_campaign_20261007/cost_remote_v1',
            'archived reference is restricted to the original CPU remote deposition campaign')
    require(all(case['state']['numerical_domain_protocol'] == 'CONSERVATIVE_FINITE_TOOL_REMOTE_DOMAIN_V1'
                for case in (long, short)), 'only CPU remote numerical protocol may be compared here')
    for case in (long, short):
        supervision = runner.read(case['stage'] / 'supervision.json')
        require(any(Path(arg).name == 'run_afsd_deposition_stage_remote_20261007.py'
                    for arg in supervision['command']), 'CPU remote stage supervision required')
    limits = [config['process_intent']['stationary_feed_end_s'] for config in configurations]
    declared_limit = configurations[1].get('diagnostic_control', {}).get('comparison_physical_end_limit_s')
    require(declared_limit is not None and time_s <= min(*limits, float(declared_limit)) + 1e-13,
            'archived reference restricted to the declared stationary comparison window')
    for key in ('classification', 'shoulder_height_m', 'free_material_trajectory_prescribed',
                'within_step_path', 'startup_duration_source'):
        require(programs[0].to_dict()[key] == programs[1].to_dict()[key],
                'machine program model differs: ' + key)
    used_phases = [[phase for phase in program.phases if phase['start'] < time_s]
                   for program in programs]
    require(runner.core.clean(used_phases[0]) == runner.core.clean(used_phases[1]),
            'actually used machine phase definitions differ')
    probes = {0., time_s}
    for case in (long, short):
        for interval in case['state']['machine_history']:
            probes.update((float(interval['t0_s']), float(interval['t1_s'])))
    for t in probes:
        require(np.array_equal(programs[0].integral(t), programs[1].integral(t))
                and np.array_equal(programs[0].velocity(t), programs[1].velocity(t)),
                'command prefix differs at an actual saved endpoint')
    return dict(mode='READ_ONLY_ARCHIVED_CPU_REMOTE_LONG_REFERENCE',
        whole_program_identical=False, actual_used_phase_definitions_identical=True,
        command_integral_and_velocity_endpoint_checks=len(probes),
        allowed_unused_future_difference='Original 3 mm track versus common-control 1 mm track; neither reached.',
        archived_reference_is_short_case_parent=False)


def analyze(long, short, archived_long_reference=False):
    require(abs(long['state']['time_s'] - short['state']['time_s']) <= 1e-13,
            'matched physical endpoint time required; no time interpolation permitted')
    parameter_delta = {k: [long['state']['parameters'][k], short['state']['parameters'][k]]
        for k in long['state']['parameters']
        if long['state']['parameters'][k] != short['state']['parameters'][k]}
    require(set(parameter_delta) == {'bar_length'}, 'only resolved bar_length may differ')
    comparison_scope = program_and_configuration_scope(long, short, archived_long_reference)
    source_maps = []
    allowed_configurations = {case['meta']['configuration'] for case in (long, short)}
    for key in set(long['meta']['source_sha256']) & set(short['meta']['source_sha256']):
        require(long['meta']['source_sha256'][key] == short['meta']['source_sha256'][key],
                'shared source/evidence content differs: ' + key)
    for case in (long, short):
        source_maps.append({k: v for k, v in case['meta']['source_sha256'].items()
                            if k not in allowed_configurations})
    require(source_maps[0] == source_maps[1], 'non-configuration implementation sources differ')
    positions = [{tuple(np.round(row, 12)): i for i, row in enumerate(xyz(case['initial_pa']))}
                 for case in (long, short)]
    common = sorted(set(positions[0]) & set(positions[1]))
    require(len(common) == len(short['initial_pa'].x), 'short initial geometry must be a subset of long')
    long_index = np.array([positions[0][key] for key in common])
    short_index = np.array([positions[1][key] for key in common])
    common_free = (np.asarray(common)[:, 2] < short['state']['parameters']['bar_length']
                   - short['state']['parameters']['actuator_depth'])
    require(np.array_equal(long['initial_state']['free'][long_index[common_free]],
                           short['initial_state']['free'][short_index[common_free]]),
            'common lower-region constraints differ')
    endpoint = {label: {branch: endpoint_metrics(case, branch) for branch in ('fine', 'coarse')}
                for label, case in (('long', long), ('short', short))}
    lf, lc, sf, sc = (endpoint['long']['fine'], endpoint['long']['coarse'],
                       endpoint['short']['fine'], endpoint['short']['coarse'])
    comparisons = {}
    for name in ('reaction_on_actuator_z_N', 'actual_mean_downward_feed_m'):
        comparisons['actuator.'+name] = comparison(sf['actuator'][name], lf['actuator'][name],
                                                   sc['actuator'][name], lc['actuator'][name])
    for group in ('initial_bore', 'initial_eligible_near_outlet', 'initial_lower_material'):
        for field in FIELDS:
            comparisons[group+'.'+field] = comparison(sf['groups'][group][field],
                lf['groups'][group][field], sc['groups'][group][field], lc['groups'][group][field])
    gap = 'nearest_eligible_uncrossed_gap_m'
    if all(row['deposition'][gap] is not None for row in (lf, lc, sf, sc)):
        comparisons['deposition.'+gap] = comparison(sf['deposition'][gap], lf['deposition'][gap],
            sc['deposition'][gap], lc['deposition'][gap])
    layer_maps = [{row['initial_z_m']: row for row in item['initial_z_layers']} for item in (lf, lc, sf, sc)]
    layer_comparisons = []
    for level in sorted(set(layer_maps[0]) & set(layer_maps[2])):
        rows = [mapping[level] for mapping in layer_maps]
        layer_comparisons.append(dict(initial_z_m=level,
            common_z_constraint_type=rows[0]['prescribed_z_particle_count'] == rows[2]['prescribed_z_particle_count'],
            metrics={field: comparison(rows[2]['metrics'][field], rows[0]['metrics'][field],
                rows[3]['metrics'][field], rows[1]['metrics'][field]) for field in FIELDS}))
    return dict(classification=('PASSIVE_EQUAL_TIME_LENGTH_COMPARISON_WITH_ARCHIVED_CPU_REMOTE_REFERENCE'
                                if archived_long_reference else
                                'PASSIVE_QUALIFIED_EQUAL_TIME_FREE_BAR_LENGTH_COMPARISON'),
        physical_solver_evaluations=0, time_s=long['state']['time_s'],
        comparison_scope=comparison_scope,
        parameter_delta=parameter_delta, non_configuration_source_files_matched=len(source_maps[0]),
        independent_initial_geometry_matching=dict(coordinate_decimals=12, common_particles=len(common),
            cross_case_ids_not_used=True,
            common_coordinates_with_different_ids=int(np.count_nonzero(
                long['initial_pa'].particle_id[long_index] != short['initial_pa'].particle_id[short_index]))),
        receipt_checks={label: dict(stage=runner.rel(case['stage']),
            manifest_sha256=sha(case['stage']/'manifest.json'),
            cost=case['meta']['cost_directory'], initial_checkpoint_sha256=case['cost']['initial_sha256'],
            final_checkpoint_sha256=case['meta']['screened_continuation_sha256'],
            latest_pair_dt_s=case['pair']['dt_s'], latest_pair_passed=case['pair']['passed'],
            verified_reference_files=case['reference_count'])
            for label, case in (('long', long), ('short', short))},
        endpoint_metrics=endpoint, equal_time_effect_vs_local_temporal_difference=comparisons,
        common_initial_layer_comparisons=layer_comparisons,
        limitations=['Independent hot-start length sensitivity, not continuation or replacement of the old GPU 120 ms trajectory.',
            'In archived-reference mode, the old CPU long-bar result is read-only evidence and is never imported as the short-bar continuation parent. Only the identical used command prefix is compared; unused future tracks differ.',
            'Bar length changes hot free length, finite inventory and prescribed-cap height together; this is not calibrated real machine support.',
            'Mean displacement is not mass flux. Geometric shoulder membership is not metallurgical bonding.',
            'Final-macro coarse/fine differences do not bound accumulated time error, spatial error, or boundary/model error.',
            'Different macro sizes and preceding adaptive histories must be considered when interpreting the local indicator.',
            'Initial free-upper-segment populations have different lengths and are reported separately, not treated as a matched cohort.',
            'A short early ramp interval cannot decide steady-feed deposition or diagnose the 89-120 ms continuation alone.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--long-stage', required=True, type=Path)
    parser.add_argument('--short-stage', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--archived-long-reference', action='store_true',
        help='Explicitly allow the verified original CPU remote long-bar campaign as a read-only same-prefix reference.')
    args = parser.parse_args()
    output = (ROOT / args.output).resolve()
    require(output.is_relative_to(Path(__file__).resolve().parent)
            and output.suffix == '.json' and not output.exists(),
            'new analysis JSON below this comparison output folder required')
    result = analyze(qualified_case(args.long_stage), qualified_case(args.short_stage),
                     archived_long_reference=args.archived_long_reference)
    result['analysis_script_sha256'] = sha(__file__)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=result['classification'], time_s=result['time_s'],
                         output=runner.rel(output)), indent=2))


if __name__ == '__main__':
    main()

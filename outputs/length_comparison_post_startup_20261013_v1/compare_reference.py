"""Read-only equal-time CPU length control against a qualified archived GPU reference.

This is analysis, never a restart/import of the archived long state into short.
The GPU bridge is explicit; its original CPU replay is checked, not treated as a
cumulative CPU/GPU error bound. No physical solver is called by this module.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
import run_afsd_deposition_stage_remote_20261007 as cpu
import run_afsd_gpu_20261007 as gpu

spec = importlib.util.spec_from_file_location('length_compare_original',
    ROOT / 'outputs/length_comparison_20261013_v1/compare.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
require, sha = base.require, base.sha
GPU_BASE = ROOT / 'outputs/gpu_acceleration_20261007'
REPLAY = GPU_BASE / 'replay_v1'
CPU_BASELINE = ROOT / 'outputs/deposition_campaign_20261007/stage_013'
ALLOWED_GPU = {REPLAY, GPU_BASE / 'campaign_10h_20261008_v1/attempt_0017',
               GPU_BASE / 'campaign_10h_20261008_v1/attempt_0026'}
REFERENCE_TIMES = {GPU_BASE / 'campaign_10h_20261008_v1/attempt_0017': .02,
                   GPU_BASE / 'campaign_10h_20261008_v1/attempt_0026': .024271484375}
SHORT_ROOTS = (ROOT / 'outputs/length_comparison_20261013_v1/short_4p0mm',
               OUT / 'short_4p0mm')
SHORT_COST = 'outputs/length_comparison_20261013_v1/short_4p0mm/cost'
SHORT_CONFIGURATION = 'outputs/flow_diagnostic_20261013_v1/length_controls/short_4p0mm.json'


def require_supervision(supervision, manifest_sha256, log_sha256):
    require(supervision['exit_code'] == 0 and supervision['source_inputs_still_match']
            and supervision['manifest_sha256'] == manifest_sha256
            and supervision['log_sha256'] == log_sha256,
            'GPU ancestor supervision does not bind a successful unchanged stage')


def require_short_origin(short):
    require(any(short['stage'].is_relative_to(folder) for folder in SHORT_ROOTS),
            'short result must be in an explicitly allowed independent CPU length-control directory')
    require(short['meta']['cost_directory'] == SHORT_COST
            and short['meta']['configuration'] == SHORT_CONFIGURATION
            and short['cost']['initial_checkpoint'] == SHORT_COST + '/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
            and short['meta']['initial_sha256'] == short['cost']['initial_sha256'],
            'short case must bind the original independent cost/initial state/configuration')


def qualified_gpu(stage):
    stage = (ROOT / stage).resolve()
    require(stage in ALLOWED_GPU, 'only explicit original-bridge reference endpoints are allowed')
    meta = gpu.gpu_receipt(stage, {gpu.REPLAY_PASS, gpu.CONTINUE_PASS})
    if stage in REFERENCE_TIMES:
        require(abs(meta['screened_end_time_s'] - REFERENCE_TIMES[stage]) <= 1e-13,
                'archived reference differs from its declared physical endpoint')
    replay = meta if stage == REPLAY else gpu.gpu_receipt(REPLAY, {gpu.REPLAY_PASS})
    require(meta['source_sha256'] == replay['source_sha256'],
            'GPU source extension must be exactly the frozen CPU replay bridge')
    require(meta['adaptation_audit'] == replay['adaptation_audit']
            and meta['adaptation_audit']['status'] == 'PASS_EXPLICIT_GPU_PROVIDERS_ONLY'
            and meta['adaptation_audit']['candidate_gates_ast_identical'],
            'unchanged explicit GPU adaptation and candidate gates required')
    require(meta['compute_protocol'] == 'CPU_PYSPH_WITH_EXPLICIT_CUDA_FP64_DENSE_SOLVES_V1',
            'unexpected compute provider')
    require(meta['cpu_baseline'] == cpu.rel(CPU_BASELINE), 'unexpected CPU baseline')
    require(all(replay['cpu_gpu_endpoint_comparisons'][role]['pass']
                for role in ('coarse', 'fine_1', 'fine_2')),
            'original same-parent CPU/GPU replay must pass all three branches')
    baseline = base.qualified_case(CPU_BASELINE)
    # Final gpu_receipt verifies its complete frozen reference map once. Walk the
    # bound parent receipts without repeatedly rehashing the same old files.
    references = meta['reference_sha256']
    cursor, current, chain, visited = stage, meta, [], set()
    while cursor != REPLAY:
        require(cursor not in visited, 'cyclic GPU parent chain')
        visited.add(cursor)
        require(current['status'] == gpu.CONTINUE_PASS, 'nonpassing GPU parent')
        require(current['source_sha256'] == replay['source_sha256'], 'GPU chain source drift')
        gpu.validate_gpu_screened_receipt(current, cursor)
        if current.get('parent_gpu_directory') is None:
            require(current.get('validation_directory') == cpu.rel(REPLAY),
                    'first GPU continuation must explicitly use the qualified replay')
            parent_dir = REPLAY
        else:
            parent_dir = (ROOT / current['parent_gpu_directory']).resolve()
        require(parent_dir.is_relative_to(GPU_BASE), 'GPU parent escaped archived outputs')
        require(parent_dir not in visited, 'cyclic GPU parent chain')
        parent_manifest = parent_dir / 'manifest.json'
        require(references.get(cpu.rel(parent_manifest)) == sha(parent_manifest),
                'GPU parent manifest is not bound by final reference receipt')
        parent_supervision = parent_dir / 'supervision.json'
        require(references.get(cpu.rel(parent_supervision)) == sha(parent_supervision),
                'GPU parent supervision is not bound by final reference receipt')
        parent_log_sha = references.get(cpu.rel(parent_dir / 'run.log'))
        require(parent_log_sha is not None, 'GPU ancestor log missing from verified reference closure')
        require_supervision(cpu.read(parent_supervision), sha(parent_manifest), parent_log_sha)
        parent = cpu.read(parent_manifest)
        require(current['parent_checkpoint'] == parent['screened_continuation_checkpoint']
                and current['parent_sha256'] == parent['screened_continuation_sha256']
                and abs(current['parent_time_s'] - parent['screened_end_time_s']) <= 1e-13,
                'GPU chain does not continue exact previous fine state')
        chain.append(dict(stage=cpu.rel(cursor), time_s=current['screened_end_time_s'],
            complete_fine_frames=[{key: row[key] for key in ('checkpoint', 'sha256', 'time_s')}
                                  for row in current['screened_history']]))
        cursor, current = parent_dir, parent
    require(current == replay, 'GPU chain replay differs from qualified replay')
    chronological_chain = [dict(stage=cpu.rel(REPLAY), time_s=replay['screened_end_time_s'],
        complete_fine_frames=[{key: row[key] for key in ('checkpoint', 'sha256', 'time_s')}
                              for row in replay['screened_history']])] + list(reversed(chain))
    checkpoint = ROOT / meta['screened_continuation_checkpoint']
    pa, state = cpu.load(checkpoint)
    coarse_path = stage / 'coarse/NOT_RESTARTABLE_candidate.npz'
    coarse_pa, coarse_state = cpu.load(coarse_path)
    require(abs(state['time_s'] - coarse_state['time_s']) <= 1e-13, 'GPU coarse endpoint mismatch')
    for current_pa in (pa, coarse_pa):
        for name in ('particle_id', 'm', 'component_id'):
            require(np.array_equal(getattr(current_pa, name), getattr(baseline['initial_pa'], name)),
                    'GPU persistent identity/inventory differs: ' + name)
    return dict(stage=stage, meta=meta, cost=baseline['cost'], initial_pa=baseline['initial_pa'],
        initial_state=baseline['initial_state'], pa=pa, state=state, checkpoint=checkpoint,
        coarse_pa=coarse_pa, coarse_state=coarse_state, coarse_checkpoint=coarse_path,
        pair=dict(meta['pair_screen'], dt_s=meta['dt_s']),
        configuration=baseline['meta']['configuration'], cpu_baseline=baseline,
        reference_count=len(references), gpu_chain_stages=len(chain), replay=replay,
        verified_gpu_stage_chain=chronological_chain)


def check_scope(long, short, replay_check=False):
    require(long['meta']['source_sha256'] == long['replay']['source_sha256'],
            'unknown GPU source extension outside frozen replay bridge')
    require(abs(long['state']['time_s'] - short['state']['time_s']) <= 1e-13,
            'matched physical endpoint time required; no interpolation')
    require(set(long['state']['parameters']) == set(short['state']['parameters']),
            'resolved parameter fields differ')
    delta = {key: [long['state']['parameters'][key], short['state']['parameters'][key]]
             for key in long['state']['parameters']
             if long['state']['parameters'][key] != short['state']['parameters'][key]}
    if replay_check:
        require(long['stage'] == REPLAY and short['stage'] == CPU_BASELINE and not delta,
                'adapter replay check is only original GPU replay versus its original CPU baseline')
    else:
        require(long['stage'] in REFERENCE_TIMES,
                'physical length comparisons are limited to the declared 20 ms and 24.271484375 ms references')
        require(set(delta) == {'bar_length'} and delta['bar_length'] == [0.0065, 0.004],
                'only explicit long/short bar_length may differ')
        require_short_origin(short)
    configurations = [cpu.read(ROOT / long['configuration']),
                      cpu.read(ROOT / short['meta']['configuration'])]
    for key in ('schema_version', 'classification', 'base_parameters', 'base_parameters_sha256',
                'geometry', 'thermal_initialization', 'interface_initialization', 'feed_mode', 'bar_rotation'):
        require(configurations[0].get(key) == configurations[1].get(key),
                'initial configuration differs beyond length: ' + key)
    require(all(case['state']['numerical_domain_protocol'] == 'CONSERVATIVE_FINITE_TOOL_REMOTE_DOMAIN_V1'
                for case in (long, short)), 'unexpected spatial protocol')
    short_supervision = cpu.read(short['stage'] / 'supervision.json')
    require(any(Path(arg).name == 'run_afsd_deposition_stage_remote_20261007.py'
                for arg in short_supervision['command']), 'CPU remote short provider required')
    programs = [cpu.model.validate_state(case['state']) for case in (long, short)]
    time_s = float(long['state']['time_s'])
    limits = [cfg['process_intent']['stationary_feed_end_s'] for cfg in configurations]
    if not replay_check:
        limits.append(float(configurations[1]['diagnostic_control']['comparison_physical_end_limit_s']))
    require(time_s <= min(limits) + 1e-13, 'outside stationary command comparison window')
    for key in ('classification', 'shoulder_height_m', 'free_material_trajectory_prescribed',
                'within_step_path', 'startup_duration_source'):
        require(programs[0].to_dict()[key] == programs[1].to_dict()[key], 'machine model differs: ' + key)
    used = [[phase for phase in program.phases if phase['start'] < time_s] for program in programs]
    require(cpu.core.clean(used[0]) == cpu.core.clean(used[1]), 'actually used command phases differ')
    probes = {0., time_s}
    for case in (long, short):
        for row in case['state']['machine_history']:
            probes.update((float(row['t0_s']), float(row['t1_s'])))
    for t in probes:
        require(np.array_equal(programs[0].integral(t), programs[1].integral(t))
                and np.array_equal(programs[0].velocity(t), programs[1].velocity(t)),
                'actual command prefix differs')
    baseline_sources = long['cpu_baseline']['meta']['source_sha256']
    configs = {long['configuration'], short['meta']['configuration']}
    physics_cpu = {key: value for key, value in baseline_sources.items() if key not in configs}
    physics_short = {key: value for key, value in short['meta']['source_sha256'].items() if key not in configs}
    require(physics_cpu == physics_short, 'CPU physical/evidence source closure differs')
    require(all(long['meta']['source_sha256'].get(key) == value for key, value in baseline_sources.items()),
            'GPU bridge does not preserve original CPU closure')
    # qualified_gpu has already bound all added source bytes to the original replay.
    extra = set(long['meta']['source_sha256']) - set(baseline_sources)
    return dict(parameter_delta=delta, common_non_configuration_sources=len(physics_cpu),
        gpu_bridge_additional_frozen_sources=len(extra), command_prefix_probes=len(probes),
        actual_used_phase_definitions_identical=True, stationary_window_end_s=min(limits),
        compute_providers=dict(long=long['meta']['compute_protocol'], short='ORIGINAL_CPU_REMOTE_PROVIDER'),
        archived_long_is_short_parent=False,
        unused_future_track_difference=None if replay_check else 'Original long 3 mm versus control 1 mm; neither reached.')


def check_initial_cohorts(long, short):
    positions = [{tuple(np.round(row, 12)): i for i, row in enumerate(base.xyz(case['initial_pa']))}
                 for case in (long, short)]
    require(all(len(mapping) == len(case['initial_pa'].x) for mapping, case in zip(positions, (long, short))),
            'initial coordinates are not unique')
    common = sorted(set(positions[0]) & set(positions[1]))
    require(len(common) == len(short['initial_pa'].x), 'short initial geometry not subset of long')
    indices = [np.array([mapping[key] for key in common]) for mapping in positions]
    free = np.asarray(common)[:, 2] < short['state']['parameters']['bar_length'] - short['state']['parameters']['actuator_depth']
    require(np.array_equal(long['initial_state']['free'][indices[0][free]],
                           short['initial_state']['free'][indices[1][free]]), 'lower-region constraints differ')
    for name in ('m', 'rho', 'temperature', 'component_id', 'ep'):
        require(np.array_equal(getattr(long['initial_pa'], name)[indices[0]],
                               getattr(short['initial_pa'], name)[indices[1]]),
                'matched initial cohort property differs: ' + name)
    return dict(common_particles=len(common), matching='initial xyz rounded to 12 decimals; not cross-case IDs',
        common_coordinates_with_different_ids=int(np.count_nonzero(
            long['initial_pa'].particle_id[indices[0]] != short['initial_pa'].particle_id[indices[1]])))


def analyze(long, short, replay_check=False):
    scope = check_scope(long, short, replay_check)
    cohorts = check_initial_cohorts(long, short)
    endpoints = {label: {branch: base.endpoint_metrics(case, branch) for branch in ('fine', 'coarse')}
                 for label, case in (('long', long), ('short', short))}
    lf, lc, sf, sc = (endpoints['long']['fine'], endpoints['long']['coarse'],
                      endpoints['short']['fine'], endpoints['short']['coarse'])
    effects = {}
    for key in ('reaction_on_actuator_z_N', 'actual_mean_downward_feed_m'):
        effects['actuator.' + key] = base.comparison(sf['actuator'][key], lf['actuator'][key],
                                                    sc['actuator'][key], lc['actuator'][key])
    for group in ('initial_bore', 'initial_eligible_near_outlet', 'initial_lower_material'):
        for key in base.FIELDS:
            effects[group + '.' + key] = base.comparison(sf['groups'][group][key], lf['groups'][group][key],
                                                        sc['groups'][group][key], lc['groups'][group][key])
    for key in ('nearest_eligible_uncrossed_gap_m',):
        if all(ep['deposition'][key] is not None for ep in (lf, lc, sf, sc)):
            effects['deposition.' + key] = base.comparison(sf['deposition'][key], lf['deposition'][key],
                                                          sc['deposition'][key], lc['deposition'][key])
    return dict(schema='PASSIVE_EXPLICIT_GPU_LONG_REFERENCE_20261013_V1',
        classification=('PASS_ADAPTER_SAME_PARENT_CPU_GPU_REPLAY_READBACK_ONLY' if replay_check else
                        'PASSIVE_QUALIFIED_EQUAL_TIME_LENGTH_COMPARISON_WITH_EXPLICIT_GPU_REFERENCE'),
        physical_solver_evaluations=0, time_s=lf['time_s'], length_effect_experiment=not replay_check,
        comparison_scope=scope, initial_coordinate_cohorts=cohorts, endpoint_metrics=endpoints,
        equal_time_effect_vs_local_temporal_difference=effects,
        receipt_checks={label: dict(stage=cpu.rel(case['stage']), manifest_sha256=sha(case['stage']/'manifest.json'),
            checkpoint_sha256=sha(case['checkpoint']), coarse_checkpoint_sha256=sha(case['coarse_checkpoint']),
            initial_checkpoint_sha256=case['cost']['initial_sha256'], latest_pair_dt_s=case['pair']['dt_s'],
            latest_pair_passed=case['pair']['passed'], verified_reference_files=case['reference_count'])
            for label, case in (('long', long), ('short', short))},
        gpu_bridge_evidence=dict(replay_directory=cpu.rel(REPLAY), replay_manifest_sha256=sha(REPLAY/'manifest.json'),
            replay_cpu_gpu_endpoint_comparisons=long['replay']['cpu_gpu_endpoint_comparisons'],
            qualified_parent_links=long['gpu_chain_stages'], shared_original_candidate_gates=True,
            verified_gpu_stage_chain=long['verified_gpu_stage_chain']),
        limitations=['Read-only archived long reference is never a short continuation parent.',
            'CPU versus GPU providers and independent adaptive time histories differ. Original same-parent replay agreement is not a cumulative error bound.',
            'Hot-start bar-length sensitivity does not alone diagnose the old 120 ms trajectory.',
            'Mean displacement is not mass flux; shoulder membership does not establish metallurgical bonding.',
            'Last-macro fine/coarse differences are local indicators, not accumulated or spatial error bounds.',
            'Initial lower material already outside the outlet is not newly supplied deposition.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--long-stage', required=True, type=Path)
    parser.add_argument('--short-stage', required=True, type=Path)
    parser.add_argument('--gpu-long-reference', action='store_true')
    parser.add_argument('--adapter-replay-check', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    require(args.gpu_long_reference, 'explicit --gpu-long-reference acknowledgement is required')
    output = (ROOT / args.output).resolve()
    require(output.is_relative_to(OUT) and output.suffix == '.json' and not output.exists(),
            'new JSON below post-startup comparison output required')
    result = analyze(qualified_gpu(args.long_stage), base.qualified_case(args.short_stage), args.adapter_replay_check)
    result['analysis_script_sha256'] = sha(__file__)
    result['reused_cohort_analysis_script_sha256'] = sha(base.__file__)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=result['classification'], time_s=result['time_s'], output=cpu.rel(output))))


if __name__ == '__main__':
    main()

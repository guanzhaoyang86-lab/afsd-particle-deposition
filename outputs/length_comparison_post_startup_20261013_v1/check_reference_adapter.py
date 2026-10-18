"""Actual readback plus negative scope checks; no physical evaluations."""
from pathlib import Path
import argparse
import json
import time
import compare_reference as adapter


def rejected(call):
    try:
        call()
    except ValueError as error:
        return dict(passed=True, rejection=str(error))
    raise AssertionError('invalid comparison was not rejected')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    output = (adapter.ROOT / args.output).resolve()
    adapter.require(output.is_relative_to(adapter.OUT) and not output.exists() and output.suffix == '.json',
                    'new JSON within adapter output directory required')
    began = time.perf_counter()
    baseline = adapter.base.qualified_case(adapter.CPU_BASELINE)
    replay = adapter.qualified_gpu(adapter.REPLAY)
    adapter.check_scope(replay, baseline, replay_check=True)
    bad_time = dict(baseline, state=dict(baseline['state'], time_s=.02))
    bad_parameters = dict(baseline['state']['parameters'])
    bad_parameters['bar_length'] += .0001
    bad_parameter_case = dict(baseline, state=dict(baseline['state'], parameters=bad_parameters))
    unknown_sources = dict(replay['meta']['source_sha256'], UNKNOWN_NOT_IN_FROZEN_BRIDGE='0'*64)
    unknown = dict(replay, meta=dict(replay['meta'], source_sha256=unknown_sources))
    short = adapter.base.qualified_case(adapter.ROOT / 'outputs/length_comparison_20261013_v1/short_4p0mm/stage_014')
    adapter.require_short_origin(short)
    # Directory-policy probes only; no state/receipt is presented as a new solve.
    new_root_probe = dict(short, stage=adapter.OUT / 'short_4p0mm/ORIGIN_POLICY_PROBE_NOT_A_PHYSICAL_RESULT')
    adapter.require_short_origin(new_root_probe)
    wrong_directory = dict(short, stage=adapter.ROOT / 'outputs/UNKNOWN_SHORT_ORIGIN')
    wrong_cost = dict(short, meta=dict(short['meta'], cost_directory='outputs/UNKNOWN_COST'))
    replay_supervision = adapter.cpu.read(adapter.REPLAY / 'supervision.json')
    failed_supervision = dict(replay_supervision, exit_code=2)
    changed_supervision = dict(replay_supervision, source_inputs_still_match=False)
    tests = dict(
        different_time_rejected=rejected(lambda: adapter.check_scope(replay, bad_time, True)),
        different_parameters_rejected=rejected(lambda: adapter.check_scope(replay, bad_parameter_case, True)),
        unknown_gpu_source_rejected=rejected(lambda: adapter.check_scope(unknown, baseline, True)),
        replay_cannot_be_length_effect=rejected(lambda: adapter.check_scope(replay, baseline, False)),
        unauthorized_gpu_directory_rejected=rejected(lambda: adapter.qualified_gpu(
            adapter.GPU_BASE / 'campaign_10h_20261008_fast_v1/attempt_0002')),
        unknown_short_directory_rejected=rejected(lambda: adapter.require_short_origin(wrong_directory)),
        wrong_short_initial_cost_rejected=rejected(lambda: adapter.require_short_origin(wrong_cost)),
        failed_ancestor_supervision_rejected=rejected(lambda: adapter.require_supervision(failed_supervision,
            replay_supervision['manifest_sha256'], replay_supervision['log_sha256'])),
        changed_ancestor_source_rejected=rejected(lambda: adapter.require_supervision(changed_supervision,
            replay_supervision['manifest_sha256'], replay_supervision['log_sha256'])))
    endpoints = []
    for attempt in (17, 26):
        stage = adapter.GPU_BASE / f'campaign_10h_20261008_v1/attempt_{attempt:04d}'
        case = adapter.qualified_gpu(stage)
        ep = adapter.base.endpoint_metrics(case, 'fine')
        endpoints.append(dict(stage=adapter.cpu.rel(stage), time_s=ep['time_s'], step=ep['step'],
            checkpoint=adapter.cpu.rel(case['checkpoint']), checkpoint_sha256=adapter.sha(case['checkpoint']),
            manifest_sha256=adapter.sha(stage/'manifest.json'),
            all_original_receipt_hashes_and_gates_passed=True, parent_links=case['gpu_chain_stages'],
            reference_file_count=case['reference_count'], gpu_source_count=len(case['meta']['source_sha256']),
            original_cpu_source_count=len(case['cpu_baseline']['meta']['source_sha256']),
            deposition=ep['deposition'], actuator=ep['actuator'],
            verified_gpu_stage_chain=case['verified_gpu_stage_chain']))
    result = dict(schema='POST_STARTUP_GPU_REFERENCE_READBACK_CHECK_20261013_V1',
        status='PASS_REFERENCE_READBACK_AND_ADAPTER_REJECTION_CHECKS', physical_solver_evaluations=0,
        actual_post_startup_length_comparison='NOT_RUN_SHORT_ENDPOINTS_NOT_YET_AVAILABLE',
        adapter_script_sha256=adapter.sha(adapter.__file__), test_script_sha256=adapter.sha(__file__),
        negative_checks=tests, archived_long_endpoints=endpoints, wall_seconds=time.perf_counter()-began)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'], endpoints=[row['time_s'] for row in endpoints],
                         tests=len(tests), wall_seconds=result['wall_seconds'], output=adapter.cpu.rel(output))))


if __name__ == '__main__':
    main()

"""Read actual fine states during the registered constant-command interval."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
import run_afsd_deposition_stage_remote_20261007 as runner
from afsd_deposition_progress_20261010 import frame_from_saved, summarize_progress


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(value, message):
    if not value:
        raise ValueError(message)


def checked(path, digest, evidence):
    path = path.resolve()
    require(path.is_relative_to(ROOT/'outputs'), 'Evidence path outside outputs')
    actual = sha(path)
    require(actual == digest, 'Changed evidence: ' + str(path))
    evidence[str(path.relative_to(ROOT)).replace('\\', '/')] = actual


def short_rows(start, end, evidence):
    cursor = ROOT/end['stage']
    expected = end['manifest_sha256']
    rows, visited = [], set()
    while True:
        require(cursor not in visited, 'Cyclic short chain')
        visited.add(cursor)
        checked(cursor/'manifest.json', expected, evidence)
        meta = read(cursor/'manifest.json')
        require(meta['screened_continuation_authorized'] and all(s['passed'] for s in meta['pair_screens']),
                'Incomplete short state')
        rows.extend(meta['screened_history'])
        if cursor == ROOT/start['stage']:
            require(expected == start['manifest_sha256']
                    and meta['screened_continuation_sha256'] == start['checkpoint_sha256'],
                    'The two short endpoints do not share this continuation chain')
            break
        require(meta.get('parent_stage_manifest'), 'Short start not found in parent chain')
        cursor = (ROOT/meta['parent_stage_manifest']).parent
        expected = meta['parent_stage_manifest_sha256']
    return rows


def long_rows(start, end, comparison, evidence):
    folder = ROOT/end['stage']
    checked(folder/'manifest.json', end['manifest_sha256'], evidence)
    meta = read(folder/'manifest.json')
    references = {**meta['reference_sha256'], **meta['artifacts_sha256']}
    rows = []
    for stage in comparison['gpu_bridge_evidence']['verified_gpu_stage_chain']:
        for row in stage['complete_fine_frames']:
            observer = (Path(row['checkpoint']).parent/'observer.json').as_posix()
            if row['time_s'] >= .02 - 1e-13:
                require(observer in references, 'Observer absent from qualified reference closure')
                rows.append(dict(row, observer_path=observer, observer_sha256=references[observer]))
    require(any(row['sha256'] == start['checkpoint_sha256'] and abs(row['time_s']-.02) < 1e-13
                for row in rows), 'Long start is not in the end reference chain')
    return rows


def observe(rows, t0, t1, evidence):
    rows = sorted((row for row in rows if t0-1e-13 <= row['time_s'] <= t1+1e-13),
                  key=lambda row: row['time_s'])
    require(rows and abs(rows[0]['time_s']-t0) < 1e-13
            and abs(rows[-1]['time_s']-t1) < 1e-13, 'Interval endpoints absent')
    frames, parameters = [], None
    for row in rows:
        checked(ROOT/row['checkpoint'], row['sha256'], evidence)
        checked(ROOT/row['observer_path'], row['observer_sha256'], evidence)
        pa, state = runner.load(ROOT/row['checkpoint'])
        require(abs(state['time_s']-row['time_s']) < 1e-13, 'Checkpoint time differs')
        if parameters is None:
            parameters = state['parameters']
        require(parameters == state['parameters'], 'Parameters changed within an interval')
        frames.append(frame_from_saved(pa, state, observer=read(ROOT/row['observer_path'])))
    progress = summarize_progress(frames, parameters)
    return progress, len(frames)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start-comparison', required=True, type=Path)
    parser.add_argument('--end-comparison', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    start_path, end_path, output = [(ROOT/p).resolve() for p in
                                    (args.start_comparison, args.end_comparison, args.output)]
    require(output.is_relative_to(BASE) and not output.exists(), 'New interval output required')
    start, end = read(start_path), read(end_path)
    expected = 'PASSIVE_QUALIFIED_EQUAL_TIME_LENGTH_COMPARISON_WITH_EXPLICIT_GPU_REFERENCE'
    require(start['classification'] == end['classification'] == expected
            and start['length_effect_experiment'] and end['length_effect_experiment'],
            'Actual qualified length comparisons required, not replay fixtures')
    t0, t1 = start['time_s'], end['time_s']
    require(abs(t0-.02) < 1e-13 and abs(t1-.024271484375) < 1e-13, 'Unregistered interval')
    evidence = {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in (start_path, end_path)}
    results = {}
    for case in ('long', 'short'):
        a, b = start['endpoint_metrics'][case]['fine'], end['endpoint_metrics'][case]['fine']
        ar, br = start['receipt_checks'][case], end['receipt_checks'][case]
        require(ar['initial_checkpoint_sha256'] == br['initial_checkpoint_sha256'], 'Initial state changed')
        require(a['deposition']['eligible'] == b['deposition']['eligible'], 'Supplied inventory changed')
        before = a['deposition']['supplied_crossed']
        after = b['deposition']['supplied_crossed']
        require(set(before['particle_ids']) <= set(after['particle_ids']), 'Crossing history regressed')
        rows = short_rows(ar, br, evidence) if case == 'short' else long_rows(ar, br, end, evidence)
        progress, frame_count = observe(rows, t0, t1, evidence)
        require(progress['final']['supplied_crossed'] == after
                and progress['final']['supplied_under_shoulder'] == b['deposition']['supplied_under_shoulder'],
                'Passive history and comparison endpoint disagree')
        feed = b['actuator']['actual_mean_downward_feed_m'] - a['actuator']['actual_mean_downward_feed_m']
        require(abs(feed-.00212*(t1-t0)) < 1e-12, 'Cap feed differs from constant registered command')
        down = -(b['groups']['initial_eligible_near_outlet']['mean_dz_m']
                 - a['groups']['initial_eligible_near_outlet']['mean_dz_m'])
        results[case] = dict(saved_fine_frames=frame_count, cap_feed_increment_m=feed,
            new_first_crossing_count=after['particle_count']-before['particle_count'],
            new_first_crossing_mass_kg=after['mass_kg']-before['mass_kg'],
            mean_first_crossing_mass_rate_kg_s=(after['mass_kg']-before['mass_kg'])/(t1-t0),
            shoulder_start=a['deposition']['supplied_under_shoulder'],
            shoulder_end=b['deposition']['supplied_under_shoulder'],
            fixed_outlet_cohort_downward_increment_m=down,
            fixed_outlet_cohort_mean_downward_velocity_m_s=down/(t1-t0),
            fixed_cohort_displacement_over_cap_feed=down/feed,
            observed_progress=progress)
    result = dict(classification='PASSIVE_CONSTANT_COMMAND_INTERVAL_FROM_ACTUAL_FINE_STATES',
        physical_solver_evaluations=0, start_time_s=t0, end_time_s=t1, cases=results,
        evidence_sha256=evidence, script_sha256=sha(Path(__file__)),
        limitations=['First-crossing mass rate is recorded first passage, not signed net instantaneous outlet flux.',
            'Fixed-cohort displacement ratio is not mass flow efficiency.',
            'Shoulder occupancy and residence are sampled geometry, not metallurgical bonding.',
            'Events shorter than the saved fine-state sampling interval cannot be excluded by these samples.',
            'The two cases retain different compute providers and time-step paths.',
            'Four milliseconds after startup does not establish sustained steady flow.'])
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(output), start_time_s=t0, end_time_s=t1,
        cases={case: {k: v for k, v in values.items() if k not in
               ('observed_progress', 'shoulder_start', 'shoulder_end')} for case, values in results.items()}), indent=2))


if __name__ == '__main__':
    main()

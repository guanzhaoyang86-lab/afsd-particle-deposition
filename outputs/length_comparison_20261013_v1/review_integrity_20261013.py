"""Read-only audit of saved length-control evidence; no physical evaluations."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
BASE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    report_path = BASE / 'review_integrity_20261013.json'
    require(not report_path.exists(), 'Review output must be new')
    summary = read(BASE / 'short_run_summary.json')
    checked = {}

    def check(relative, expected):
        if relative not in checked:
            checked[relative] = sha(ROOT / relative)
        require(checked[relative] == expected, 'Hash mismatch: ' + relative)

    previous = None
    pairs = set()
    states = set()
    candidate_count = 0
    for row in summary['stages']:
        check(row['manifest_path'], row['manifest_sha256'])
        check(row['supervision_path'], row['supervision_sha256'])
        manifest = read(ROOT / row['manifest_path'])
        supervision = read(ROOT / row['supervision_path'])
        require(supervision['exit_code'] == 0 and supervision['source_inputs_still_match'],
                'Supervisor failure')
        require(supervision['manifest_sha256'] == row['manifest_sha256'], 'Unbound manifest')
        check(row['directory'] + '/run.log', supervision['log_sha256'])
        check(row['parent_checkpoint'], row['parent_checkpoint_sha256'])
        check(row['screened_continuation_checkpoint'], row['screened_continuation_sha256'])
        require(manifest['candidate_attempt_count'] == row['candidate_attempt_count'], 'Candidate count differs')
        require(manifest['screened_end_time_s'] == row['screened_end_time_s'], 'Endpoint time differs')
        require(manifest['screened_end_step'] == row['screened_end_step'], 'Endpoint step differs')
        if previous:
            require(row['parent_checkpoint'] == previous['screened_continuation_checkpoint']
                    and row['parent_checkpoint_sha256'] == previous['screened_continuation_sha256']
                    and row['parent_time_s'] == previous['screened_end_time_s']
                    and row['parent_step'] == previous['screened_end_step'], 'Broken stage chain')
        stage_pairs = set()
        for state in manifest['screened_history']:
            check(state['checkpoint'], state['sha256'])
            check(state['observer_path'], state['observer_sha256'])
            check(state['screen_path'], state['screen_sha256'])
            require(read(ROOT / state['screen_path'])['passed'], 'Selected pair failed')
            stage_pairs.add(state['screen_path'])
            states.add(state['checkpoint'])
        require(len(stage_pairs) == row['complete_macro_count'], 'Macro count differs')
        pairs.update(stage_pairs)
        candidate_count += row['candidate_attempt_count']
        previous = row
    require(len(pairs) == summary['totals']['complete_macro_count'] == 52, 'Total macros differ')
    require(len(states) == summary['final_screened_step'] == 104, 'Total fine states differ')
    require(candidate_count == summary['totals']['candidate_attempt_count'] == 171, 'Total candidates differ')
    original = ('outputs/gpu_acceleration_20261007/nonlinear_resume_20261011_v1/'
                'attempt_0086/fine_2/NOT_RESTARTABLE_candidate.npz')
    check(original, '987b3919c11a02a36f0c27facb80822ab3a6a1ea1f63a46907315b29ae225df0')
    old = read(BASE / 'comparison_11p5625ms.json')
    current = read(BASE / 'review_comparison_11p5625ms_20261013.json')
    require(old == current, 'Fresh original-runner readback differs from the saved comparison')
    report = dict(status='PASS_READ_ONLY_LENGTH_REVIEW', physical_evaluations=0,
        stage_count=len(summary['stages']), complete_pairs=len(pairs), fine_states=len(states),
        candidate_count=candidate_count, unique_artifacts_rehashed=len(checked),
        contiguous_parent_chain=True, fresh_comparison_unchanged=True,
        original_120ms_checkpoint_unchanged=True, summary_sha256=sha(BASE/'short_run_summary.json'),
        readback_sha256=sha(BASE/'review_comparison_11p5625ms_20261013.json'),
        script_sha256=sha(Path(__file__)), checked_sha256=checked,
        scope='Saved evidence integrity and fresh original-runner equal-time readback; not a physical validation or new simulation.')
    with report_path.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'checked_sha256'}, indent=2))


if __name__ == '__main__':
    main()

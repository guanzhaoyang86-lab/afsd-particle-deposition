"""Freeze completed E3 continuation evidence; preserve all earlier failures.

Run only after final reports/HANDOFF are saved.  This checks byte integrity and
recorded status, not scientific validity.  It neither advances nor edits physics.
"""
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/SHORT_WINDOW_E3_FREEZE.json'
TARGET = ROOT/'docs/CONTINUATION_E3_FREEZE.json'
REPORT = ROOT/'logs/continuation_e3_20261001/preservation_report.json'
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {
    'docs/CONTINUATION_E3_FREEZE.json',
    'logs/continuation_e3_20261001/freeze_create.log',
    'logs/continuation_e3_20261001/freeze_final.log',
}
REQUIRED = (
    'docs/CONTINUATION_E3_CONTRACT.md',
    'docs/CONTINUATION_E3_REPORT.md',
    'docs/MECHANISM_OBSERVATION_CONTRACT.md',
    'docs/PROCESS_TIME_DECISION.md',
    'scripts/inspect_afsd_chain.py',
    'scripts/measure_afsd_mechanisms.py',
    'scripts/freeze_continuation_e3_evidence.py',
    'tests/review/test_afsd_chain.py',
    'tests/review/test_mechanism_observation.py',
    'tests/review/test_continuation_e3_physical_evidence.py',
    'logs/short_window_e3_20261001/freeze_final.log',
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def safe_file(relative):
    result = (ROOT/relative).resolve()
    if not result.is_relative_to(ROOT) or not result.is_file():
        raise RuntimeError('missing/unsafe evidence path: '+relative)
    return result


def successful_test_log(relative, expected_count):
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8').strip()
    matches = re.findall(r'Ran (\d+) tests? in ', text)
    if matches != [str(expected_count)] or text.splitlines()[-1].strip() != 'OK':
        raise RuntimeError('unexpected unit-test result: '+relative)
    return dict(status='PASS', tests_run=expected_count, failures=0, errors=0,
                log=relative, log_sha256=sha(path))


def check_final_evidence():
    directory = ROOT/'outputs/continuation_e3_20261001'
    chain = read(directory/'chain_report.json')
    mechanism = read(directory/'mechanism_report_v2.json')
    if (chain['readback_status'] != 'PASS' or chain['status'] != 'ACCEPTED_TARGET_REACHED'
            or chain['reached_target'] is not True or chain['target_end_time_s'] != 2e-6
            or chain['accepted_end_step'] != 191 or chain['accepted_end_time_s'] != 2e-6
            or chain['accepted_unique_frames'] != 192):
        raise RuntimeError('chain report does not certify the recorded 191-step/2us endpoint')
    for name in ('production_authorized', 'paper_reproduction', 'deposition_validated', 'mixing_validated'):
        if chain[name] is not False:
            raise RuntimeError('unsupported chain scientific claim: '+name)
    for segment in chain['segments']:
        run = Path(segment['directory']).resolve()
        if not run.is_relative_to(ROOT):
            raise RuntimeError('chain segment outside project')
        for name, digest in segment['evidence_sha256'].items():
            if sha(run/name) != digest:
                raise RuntimeError('chain evidence changed after readback: '+str(run/name))
        if sha(run/'checkpoint.npz') != segment['checkpoint_sha256']:
            raise RuntimeError('accepted checkpoint changed after readback')
    final_run = directory/'to_2us'
    performance = read(final_run/'performance.json')
    manifest = read(final_run/'manifest.json')
    if performance['status'] != 'COMPLETE' or manifest['status'] != 'COMPLETE':
        raise RuntimeError('final run not COMPLETE in both recorded statuses')
    if performance['wall_seconds'] > 540 or performance['requested_end_time_s'] != 2e-6:
        raise RuntimeError('recorded real run exceeds its frozen scope/budget')
    for path in final_run.glob('REJECTED*'):
        raise RuntimeError('unexpected rejected artifact in completed continuation: '+str(path))
    if (mechanism['accepted_step'] != 191 or mechanism['accepted_time_s'] != 2e-6
            or mechanism['chain_status'] != 'ACCEPTED_TARGET_REACHED'
            or mechanism['reached_requested_target'] is not True
            or mechanism['production_authorized'] is not False
            or mechanism['paper_comparison'] is not False):
        raise RuntimeError('mechanism report endpoint/scope mismatch')
    last = mechanism['final_observation']
    if (last['sustained_deposition_validated'] is not False
            or last['free_plastic_above_floor_count'] != 3
            or last['free_plastic_particle_ids'] != [541, 542, 543]
            or mechanism['final_cumulative_plastic_work_J'] <= 0
            or mechanism['final_ever_feed_plane_crossed_mass_kg'] != 0):
        raise RuntimeError('unexpected observed plasticity/feed/deposition evidence')
    tests = {
        'chain_evidence': successful_test_log(
            'logs/continuation_e3_20261001/chain_tests_completion_gate.log', 15),
        'mechanism_observation': successful_test_log(
            'logs/continuation_e3_20261001/mechanism_tests_final.log', 7),
        'continuation_physical_evidence': successful_test_log(
            'logs/continuation_e3_20261001/physical_evidence_tests.log', 6),
        'parent_short_window_evidence': successful_test_log(
            'logs/continuation_e3_20261001/parent_checks.log', 7),
        'full_116_suite_this_round': 'NOT_RUN; prior 111 PASS / 5 FAIL retained, not reclassified',
        'dt_or_space_convergence_this_round': 'NOT_RUN',
        'gpu_parallel_and_production': 'NOT_RUN',
    }
    return chain, mechanism, performance, tests


def write_new(path, value):
    # Exclusive mode remains fail-closed if a target appears after the initial check.
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    if TARGET.exists() or REPORT.exists():
        raise FileExistsError('never overwrite a freeze or preservation report')
    old = read(PRIOR)
    if len(old['files_sha256']) != 907:
        raise RuntimeError('unexpected parent evidence inventory; expected 907 entries')
    prior_verify = read(ROOT/'logs/short_window_e3_20261001/freeze_final.log')
    if prior_verify['status'] != 'PASS' or prior_verify['files_checked'] != 907:
        raise RuntimeError('parent final integrity verification is not PASS907')
    changes = []
    hashes = {}
    for relative, expected in old['files_sha256'].items():
        path = safe_file(relative)
        hashes[relative] = sha(path)
        if hashes[relative] != expected:
            changes.append(relative)
    unexpected = set(changes)-ALLOWED
    if unexpected:
        raise RuntimeError('out-of-scope historical changes: '+str(sorted(unexpected)))
    # Validate every required final document before creating either output file.
    files = {safe_file(relative) for relative in old['files_sha256']}
    files.add(PRIOR)
    files.update(safe_file(relative) for relative in REQUIRED)
    for name in ('logs/continuation_e3_20261001', 'outputs/continuation_e3_20261001'):
        directory = ROOT/name
        if not directory.is_dir():
            raise RuntimeError('missing continuation evidence directory: '+name)
        for path in directory.rglob('*'):
            if path.is_file() and path.relative_to(ROOT).as_posix() not in EXCLUDED:
                files.add(path)
    chain, mechanism, performance, tests = check_final_evidence()
    preservation = dict(
        status='PASS', scope='byte integrity and recorded evidence only; NOT scientific acceptance',
        previous_entries=907, previous_freeze=PRIOR.relative_to(ROOT).as_posix(),
        previous_freeze_sha256=sha(PRIOR),
        authorized_modified_previous_entries=sorted(changes),
        authorized_changes={name: {'before': old['files_sha256'][name], 'after': hashes[name]}
                            for name in sorted(changes)},
        unchanged_previous_entries=907-len(changes),
        all_old_solver_sources_configs_tests_and_failure_artifacts_preserved=True,
        old_failure_records_preserved=True, production_authorized=False,
        real_run_status=performance['status'], accepted_end_step=chain['accepted_end_step'],
        accepted_end_time_s=chain['accepted_end_time_s'], real_runner_wall_seconds=performance['wall_seconds'],
        observed_free_plastic_particles=mechanism['final_observation']['free_plastic_particle_ids'],
        sustained_deposition_validated=False)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    write_new(REPORT, preservation)
    files.add(REPORT)
    result = dict(
        schema_version='1.0', date='2026-10-01',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False,
        scientific_status='ENGINEERING_APPROXIMATE; restarted 2us implementation interval COMPLETE; local free-bar startup plasticity observed, NOT sustained feeding/deposition or quantitative validation',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=preservation, actual_tests=tests,
        current_new_test_count=28, current_parent_recheck_count=7,
        excluded_self_referential_or_future_paths=sorted(EXCLUDED),
        files_sha256={path.relative_to(ROOT).as_posix(): sha(path) for path in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=preservation,
                          production_authorized=False), indent=2))


if __name__ == '__main__':
    main()

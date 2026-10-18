"""Freeze Q1 benchmark evidence without advancing physics or changing history.

Run only after reports, tests and HANDOFF are final. Byte integrity and recorded
benchmark completion are not AFSD validation or production authorization.
"""
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/CONTINUATION_E3_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_TIME_Q1_FREEZE.json'
REPORT = ROOT/'logs/process_time_q1_20261001/preservation_report.json'
PARENT_COUNT = 1093
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {
    'docs/PROCESS_TIME_Q1_FREEZE.json',
    'logs/process_time_q1_20261001/freeze_create.log',
    'logs/process_time_q1_20261001/freeze_final.log',
}
REQUIRED = (
    'docs/PROCESS_TIME_Q1_CONTRACT.md',
    'docs/PROCESS_TIME_Q1_REPORT.md',
    'configs/benchmarks/process_time_q1.json',
    'scripts/probe_process_time.py',
    'scripts/inspect_process_time.py',
    'scripts/freeze_process_time_q1.py',
    'tests/review/test_quasistatic_block_q1.py',
    'tests/review/test_process_time_readback.py',
    'logs/continuation_e3_20261001/freeze_final.log',
    # These transitive inputs are not all listed by individual run reports.
    'configs/benchmarks/engineering_block_e2.json',
    'src/afsd_backend/pysph/kernel.py',
    'src/afsd_scaffold/vtk_legacy.py',
)
CASES = {'n3_i4': (3, 4), 'n5_i4': (5, 4), 'n5_i8': (5, 8), 'n7_i4': (7, 4)}
UNSOLVED = {
    'inertia', 'thermal', 'contact', 'FKS', 'plasticity_validation',
    'finite_strain_density', 'AFSD', 'production_restart',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def safe_file(relative):
    path = (ROOT/relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise RuntimeError('missing/unsafe evidence path: '+str(relative))
    return path


def successful_test_log(relative, count):
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8').strip()
    if (re.findall(r'Ran (\d+) tests? in ', text) != [str(count)]
            or not text or text.splitlines()[-1].strip() != 'OK'
            or re.search(r'^FAILED\b', text, re.MULTILINE)):
        raise RuntimeError('unexpected unit-test result: '+relative)
    return dict(status='PASS', tests_run=count, failures=0, errors=0,
                log=relative, log_sha256=sha(path))


def check_evidence():
    base = ROOT/'outputs/process_time_q1_20261001'
    comparison = read(base/'comparison.json')
    if (comparison['status'] != 'PASS_READBACK'
            or comparison['scope'] != 'Q1_QUASISTATIC_ONLY_NOT_AFSD'
            or comparison['production_authorized'] is not False
            or comparison['continuum_precision_gate_applied'] is not False):
        raise RuntimeError('comparison is not a completed benchmark-only readback')
    comparison_rows = comparison['rows']
    if (len(comparison_rows) != len(CASES)
            or {(r['n'], r['increments']) for r in comparison_rows} != set(CASES.values())):
        raise RuntimeError('comparison inventory differs from the four declared cases')
    compared = {(r['n'], r['increments']): r for r in comparison_rows}
    cfg = read(ROOT/'configs/benchmarks/process_time_q1.json')
    if (cfg['production_authorized'] is not False
            or set(cfg['not_solved']) != UNSOLVED
            or cfg['parameters']['ramp_duration']['value'] != .1):
        raise RuntimeError('Q1 configuration scope changed')
    runs = {}
    for name, (n, increments) in CASES.items():
        directory = base/name
        r = read(directory/'report.json')
        if (r['status'] != 'DISCRETE_EQUILIBRIUM_COMPLETE'
                or r['accepted_time_s'] != .1 or r['accepted_steps'] != increments
                or r['n'] != n or r['increments_requested'] != increments
                or r['production_authorized'] is not False or r['restart_allowed'] is not False
                or r['classification'] != cfg['classification']
                or set(r['not_solved']) != UNSOLVED or r['error'] is not None):
            raise RuntimeError('case completion/scope mismatch: '+name)
        if len(r['records']) != increments+1 or r['records'][-1]['physical_time_s'] != .1:
            raise RuntimeError('case accepted-record inventory mismatch: '+name)
        if not {'history.npz', 'checkpoint.npz', 'final.vtk', 'resolved_config.json'} <= set(r['artifact_sha256']):
            raise RuntimeError('case artifact hash inventory incomplete: '+name)
        for relative, digest in r['source_sha256'].items():
            if sha(safe_file(relative)) != digest:
                raise RuntimeError('source changed after Q1 run: '+relative)
        for relative, digest in r['artifact_sha256'].items():
            path = (directory/relative).resolve()
            if not path.is_relative_to(directory.resolve()) or not path.is_file() or sha(path) != digest:
                raise RuntimeError('missing/changed/unsafe Q1 artifact: '+name+'/'+relative)
        row = compared[(n, increments)]
        if (row['status'] != 'PASS_READBACK_ONLY' or row['common_end_s'] != .1
                or row['production_authorized'] is not False or row['not_afsd'] is not True
                or row['continuum_accuracy'] != r['continuum_accuracy']
                or row['compute_wall_s'] != r['wall_seconds']):
            raise RuntimeError('comparison does not match completed run: '+name)
        runs[name] = dict(status=r['status'], particles=r['particles'], free_dofs=r['free_dofs'],
            accepted_increments=r['accepted_steps'], benchmark_load_end_s=r['accepted_time_s'],
            compute_wall_s=r['wall_seconds'], report_sha256=sha(directory/'report.json'))
    time_comparison = comparison['time_increment_comparison']
    if (time_comparison['common_end_s'] != .1 or time_comparison['coarse_increments'] != 4
            or time_comparison['fine_increments'] != 8):
        raise RuntimeError('increment comparison does not use common 0.1s endpoint')
    tests = {
        'quasistatic_block': successful_test_log('logs/process_time_q1_20261001/q1_tests.log', 9),
        'saved_evidence_readback': successful_test_log('logs/process_time_q1_20261001/readback_tests.log', 4),
        'four_case_saved_data_comparison': 'PASS_READBACK',
        'full_prior_suite_this_round': 'NOT_RUN; old PASS/FAIL records preserved, not reclassified',
        'afsd_production_or_continuation_this_round': 'NOT_RUN',
        'finite_strain_fks_plastic_contact_thermal_quasistatic': 'NOT_RUN',
        'gpu_parallel_optimization': 'NOT_RUN',
    }
    return runs, tests


def write_new(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    if TARGET.exists() or REPORT.exists():
        raise FileExistsError('never overwrite a freeze or preservation report')
    old = read(PRIOR)
    if len(old['files_sha256']) != PARENT_COUNT:
        raise RuntimeError('unexpected parent evidence inventory; expected 1093 entries')
    parent_check = read(ROOT/'logs/continuation_e3_20261001/freeze_final.log')
    if parent_check['status'] != 'PASS' or parent_check['files_checked'] != PARENT_COUNT:
        raise RuntimeError('parent final integrity verification is not PASS1093')
    hashes = {relative: sha(safe_file(relative)) for relative in old['files_sha256']}
    changed = sorted(relative for relative, digest in old['files_sha256'].items()
                     if hashes[relative] != digest)
    if set(changed)-ALLOWED:
        raise RuntimeError('out-of-scope historical changes: '+str(sorted(set(changed)-ALLOWED)))
    files = {safe_file(relative) for relative in old['files_sha256']}
    files.add(PRIOR)
    files.update(safe_file(relative) for relative in REQUIRED)
    for name in ('logs/process_time_q1_20261001', 'outputs/process_time_q1_20261001'):
        directory = ROOT/name
        if not directory.is_dir():
            raise RuntimeError('missing Q1 evidence directory: '+name)
        for path in directory.rglob('*'):
            if path.is_file() and path.relative_to(ROOT).as_posix() not in EXCLUDED:
                files.add(safe_file(path.relative_to(ROOT).as_posix()))
    runs, tests = check_evidence()
    preservation = dict(
        status='PASS', scope='byte integrity and recorded evidence only; NOT scientific acceptance',
        previous_entries=PARENT_COUNT, previous_freeze=PRIOR.relative_to(ROOT).as_posix(),
        previous_freeze_sha256=sha(PRIOR), authorized_modified_previous_entries=changed,
        authorized_changes={name: {'before': old['files_sha256'][name], 'after': hashes[name]}
                            for name in changed},
        unchanged_previous_entries=PARENT_COUNT-len(changed),
        all_old_solver_sources_configs_tests_and_failure_artifacts_preserved=True,
        old_failure_records_preserved=True, production_authorized=False,
        benchmark_scope='elastic small-strain quasi-static equilibrium; NOT AFSD',
        completed_cases=runs, not_solved=sorted(UNSOLVED),
        timing_note='compute_wall_s excludes output serialization and interpreter/import startup')
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    write_new(REPORT, preservation)
    files.add(REPORT)
    result = dict(
        schema_version='1.0', date='2026-10-01',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False,
        scientific_status='ENGINEERING_APPROXIMATE elastic small-strain quasi-static benchmark; four discrete equilibrium cases COMPLETE at 0.1s load endpoint; NOT AFSD process time or deposition validation',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=preservation, actual_tests=tests, current_new_test_count=13,
        continuum_accuracy_is_admission_gate=False,
        excluded_self_referential_or_future_paths=sorted(EXCLUDED),
        files_sha256={path.relative_to(ROOT).as_posix(): sha(path) for path in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=preservation,
                         production_authorized=False), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

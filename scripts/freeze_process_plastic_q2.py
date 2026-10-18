"""Freeze Q2 accepted and rejected evidence; never advance physics.

Run after all reports, readbacks, tests and HANDOFF are final. A frozen BLOCKED
case stays BLOCKED; byte integrity is not AFSD or production validation.
"""
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_TIME_Q1_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_PLASTIC_Q2_FREEZE.json'
REPORT = ROOT/'logs/process_plastic_q2_20261001/preservation_report.json'
PARENT_COUNT = 1134
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {
    'docs/PROCESS_PLASTIC_Q2_FREEZE.json',
    'logs/process_plastic_q2_20261001/freeze_create.log',
    'logs/process_plastic_q2_20261001/freeze_final.log',
}
REQUIRED = (
    'docs/PROCESS_PLASTIC_Q2_CONTRACT.md',
    'docs/PROCESS_PLASTIC_Q2_REPORT.md',
    'configs/benchmarks/process_plastic_q2.json',
    'scripts/probe_process_plastic.py',
    'scripts/inspect_process_plastic.py',
    'scripts/freeze_process_plastic_q2.py',
    'tests/review/test_process_plastic_q2.py',
    'tests/review/test_process_plastic_readback.py',
    'logs/process_time_q1_20261001/freeze_final.log',
    # Transitive operator/material inputs remain frozen even when individual
    # run reports list only their direct source files.
    'configs/benchmarks/process_time_q1.json',
    'configs/benchmarks/engineering_block_e2.json',
    'src/afsd_backend/pysph/kernel.py',
    'src/afsd_scaffold/vtk_legacy.py',
)
CASES = {
    'j2_n3_i8': ('j2', 3, 8),
    'j2_n3_i16': ('j2', 3, 16),
    'j2_n5_i8': ('j2', 5, 8),
    'fks_n3_i8': ('fks', 3, 8),
}
UNSOLVED = {
    'inertia', 'thermal', 'contact', 'free_unloaded_release',
    'finite_strain_density', 'AFSD', 'production_restart',
}
COMPLETE = 'DISCRETE_PLASTIC_PATH_COMPLETE'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def safe_file(relative):
    path = (ROOT/relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise RuntimeError('missing/unsafe evidence path: '+str(relative))
    return path


def successful_test_log(relative, expected=None):
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8').strip()
    matches = re.findall(r'Ran (\d+) tests? in ', text)
    if (len(matches) != 1 or int(matches[0]) <= 0
            or (expected is not None and int(matches[0]) != expected)
            or not text or text.splitlines()[-1].strip() != 'OK'
            or re.search(r'^FAILED\b', text, re.MULTILINE)):
        raise RuntimeError('unexpected unit-test result: '+relative)
    return dict(status='PASS', tests_run=int(matches[0]), failures=0, errors=0,
                log=relative, log_sha256=sha(path))


def preserved_initial_failure():
    relative = 'logs/process_plastic_q2_20261001/q2_tests_initial.log'
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8').strip()
    if (re.findall(r'Ran (\d+) tests? in ', text) != ['15']
            or text.splitlines()[-1].strip() != 'FAILED (failures=1)'):
        raise RuntimeError('initial test failure missing or reclassified')
    return dict(status='FAIL_PRESERVED', tests_run=15, failures=1, errors=0,
                log=relative, log_sha256=sha(path),
                note='Original fixture failure retained; later fixture correction has a separate log.')


def preserved_readback_red():
    relative = 'logs/process_plastic_q2_20261001/readback_tests_red.log'
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8').strip()
    counts = re.findall(r'Ran (\d+) tests? in ', text)
    failures = re.search(r'^FAILED \(([^)]+)\)$', text, re.MULTILINE)
    if len(counts) != 1 or int(counts[0]) <= 0 or failures is None:
        raise RuntimeError('reader RED regression evidence missing or reclassified')
    return dict(status='FAIL_PRESERVED', tests_run=int(counts[0]),
                original_failure_summary=failures.group(1), log=relative,
                log_sha256=sha(path),
                note='Reader records/history inventory regression before minimal reader fix.')


def check_evidence():
    base = ROOT/'outputs/process_plastic_q2_20261001'
    cfg = read(ROOT/'configs/benchmarks/process_plastic_q2.json')
    if (cfg['production_authorized'] is not False
            or set(cfg['not_solved']) != UNSOLVED
            or cfg['parameters']['phase_duration']['value'] != .1
            or cfg['parameters']['peak_compressive_strain']['value'] != .004
            or {tuple(case) for case in cfg['cases']} != set(CASES.values())):
        raise RuntimeError('Q2 configuration scope/case inventory changed')
    runs = {}
    for name, (law, n, increments) in CASES.items():
        directory = base/name
        r = read(directory/'report.json')
        rb_path = base/(name+'_readback_v2.json')
        rb = read(rb_path)
        if (r['status'] not in (COMPLETE, 'BLOCKED')
                or r['selected_law'] != law or r['n'] != n
                or r['increments_per_phase'] != increments
                or r['production_authorized'] is not False
                or r['classification'] != cfg['classification']
                or set(r['not_solved']) != UNSOLVED
                or r['requested_end_s'] != .2):
            raise RuntimeError('case status/scope mismatch: '+name)
        records = r['records']
        if (len(records) != r['accepted_steps']+1
                or records[0]['physical_time_s'] != 0.
                or records[-1]['physical_time_s'] != r['accepted_time_s']
                or [row['step'] for row in records] != list(range(len(records)))
                or any(b['physical_time_s'] <= a['physical_time_s']
                       for a, b in zip(records, records[1:]))):
            raise RuntimeError('accepted-record inventory/time mismatch: '+name)
        if r['status'] == COMPLETE:
            if (r['accepted_steps'] != 2*increments or r['accepted_time_s'] != .2
                    or r['error'] is not None):
                raise RuntimeError('completion claimed without full accepted path: '+name)
            expected_readback = 'PASS_COMPLETE_READBACK'
        else:
            if (not r['error'] or not 0 <= r['accepted_steps'] < 2*increments
                    or not 0 <= r['accepted_time_s'] < .2):
                raise RuntimeError('BLOCKED endpoint/error mismatch: '+name)
            rejected = read(directory/'REJECTED_manifest.json')
            if (rejected['classification'] != 'REJECTED'
                    or rejected['restart_allowed'] is not False
                    or rejected['last_accepted_step'] != r['accepted_steps']
                    or rejected['error'] != r['error']):
                raise RuntimeError('rejected candidate promoted/lost: '+name)
            if 'REJECTED_manifest.json' not in r['artifact_sha256']:
                raise RuntimeError('rejected manifest lacks recorded hash: '+name)
            expected_readback = 'PASS_PARTIAL_READBACK'
        # The two required completed n3 J2 paths were observed before this
        # freeze helper was requested. Other declared cases may remain BLOCKED.
        if name in ('j2_n3_i8', 'j2_n3_i16') and r['status'] != COMPLETE:
            raise RuntimeError('previously completed n3 path changed: '+name)
        if not {'history.npz', 'checkpoint.npz', 'final.vtk', 'resolved_config.json'} <= set(r['artifact_sha256']):
            raise RuntimeError('case artifact hash inventory incomplete: '+name)
        for relative, digest in r['source_sha256'].items():
            if sha(safe_file(relative)) != digest:
                raise RuntimeError('source changed after Q2 run: '+relative)
        for relative, digest in r['artifact_sha256'].items():
            path = (directory/relative).resolve()
            if not path.is_relative_to(directory.resolve()) or not path.is_file() or sha(path) != digest:
                raise RuntimeError('missing/changed/unsafe Q2 artifact: '+name+'/'+relative)
        if (rb['status'] != expected_readback or rb['completion_status'] != r['status']
                or rb['accepted_end_s'] != r['accepted_time_s']
                or rb['accepted_steps'] != r['accepted_steps']
                or rb['selected_law'] != law or rb['n'] != n
                or rb['production_authorized'] is not False or rb['final'] != records[-1]):
            raise RuntimeError('readback does not match accepted/blocked run: '+name)
        runs[name] = dict(status=r['status'], error=r['error'], selected_law=law,
            particles=r['particles'], free_dofs=r['free_dofs'],
            accepted_increments=r['accepted_steps'], benchmark_load_return_end_s=r['accepted_time_s'],
            requested_end_s=.2, wall_before_report_s=r['wall_before_report_s'],
            readback_status=rb['status'], report_sha256=sha(directory/'report.json'),
            readback_sha256=sha(rb_path))
    tests = {
        'plastic_block_material_and_transactions': successful_test_log(
            'logs/process_plastic_q2_20261001/q2_tests_fixture_fixed.log', 15),
        'saved_evidence_readback': successful_test_log(
            'logs/process_plastic_q2_20261001/readback_tests_fixed.log'),
        'original_saved_evidence_readback': successful_test_log(
            'logs/process_plastic_q2_20261001/readback_tests.log', 6),
        'original_test_failure': preserved_initial_failure(),
        'reader_red_regression': preserved_readback_red(),
        'full_prior_suite_this_round': 'NOT_RUN; old PASS/FAIL records preserved, not reclassified',
        'afsd_production_or_continuation_this_round': 'NOT_RUN',
        'finite_strain_contact_thermal_quasistatic': 'NOT_RUN',
        'free_unloaded_release': 'NOT_RUN; axial actuator return is not free unloading',
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
        raise RuntimeError('unexpected parent evidence inventory; expected 1134 entries')
    parent_check = read(ROOT/'logs/process_time_q1_20261001/freeze_final.log')
    if parent_check['status'] != 'PASS' or parent_check['files_checked'] != PARENT_COUNT:
        raise RuntimeError('parent final integrity verification is not PASS1134')
    hashes = {relative: sha(safe_file(relative)) for relative in old['files_sha256']}
    changed = sorted(relative for relative, digest in old['files_sha256'].items()
                     if hashes[relative] != digest)
    if set(changed)-ALLOWED:
        raise RuntimeError('out-of-scope historical changes: '+str(sorted(set(changed)-ALLOWED)))
    files = {safe_file(relative) for relative in old['files_sha256']}
    files.add(PRIOR)
    files.update(safe_file(relative) for relative in REQUIRED)
    for name in ('logs/process_plastic_q2_20261001', 'outputs/process_plastic_q2_20261001'):
        directory = ROOT/name
        if not directory.is_dir():
            raise RuntimeError('missing Q2 evidence directory: '+name)
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
        benchmark_scope='ENGINEERING_APPROXIMATE small-strain displacement-controlled plastic benchmark; NOT AFSD',
        case_results=runs, not_solved=sorted(UNSOLVED),
        axial_actuator_return_is_free_unloading=False,
        plastic_work_is_solved_heat=False,
        timing_note='wall_before_report_s excludes report serialization and interpreter/import startup')
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    write_new(REPORT, preservation)
    files.add(REPORT)
    result = dict(
        schema_version='1.0', date='2026-10-01',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False,
        scientific_status='ENGINEERING_APPROXIMATE small-strain displacement-controlled plastic benchmark; see per-case COMPLETE/BLOCKED status; NOT AFSD or freely unloaded shape validation',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=preservation, actual_tests=tests,
        current_new_passing_test_count=sum(tests[key]['tests_run'] for key in
            ('plastic_block_material_and_transactions', 'saved_evidence_readback')),
        original_core_test_failures_preserved=1,
        reader_red_regression_preserved=True,
        continuum_accuracy_is_admission_gate=False,
        excluded_self_referential_or_future_paths=sorted(EXCLUDED),
        files_sha256={path.relative_to(ROOT).as_posix(): sha(path) for path in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=preservation,
                         production_authorized=False), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

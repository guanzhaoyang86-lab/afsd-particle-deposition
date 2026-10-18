"""Freeze Q3 evidence after reports/readbacks are final; never run physics.

Byte integrity is not AFSD validation. BLOCKED cases and the original failing
regression remain classified as such, and no existing freeze is overwritten.
"""
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_PLASTIC_Q2_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_CONTACT_Q3_FREEZE.json'
LOGS = 'logs/process_contact_q3_20261001'
REPORT = ROOT/LOGS/'preservation_report.json'
PARENT_COUNT = 1193
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {
    'docs/PROCESS_CONTACT_Q3_FREEZE.json',
    LOGS+'/freeze_create.log',
    LOGS+'/freeze_final.log',
}
REQUIRED = (
    'docs/PROCESS_CONTACT_Q3_CONTRACT.md',
    'docs/PROCESS_CONTACT_Q3_REFERENCES.md',
    'docs/PROCESS_CONTACT_Q3_REPORT.md',
    'configs/benchmarks/process_contact_q3.json',
    'scripts/probe_process_contact.py',
    'scripts/inspect_process_contact.py',
    'scripts/freeze_process_contact_q3.py',
    'tests/review/test_process_contact_q3.py',
    'tests/review/test_process_contact_readback.py',
    'logs/process_plastic_q2_20261001/freeze_final.log',
    # Direct reports omit some transitive inputs; preserve them as well.
    'configs/benchmarks/process_plastic_q2.json',
    'configs/benchmarks/process_time_q1.json',
    'configs/benchmarks/engineering_block_e2.json',
    'src/afsd_backend/pysph/kernel.py',
    'src/afsd_backend/pysph/output.py',
    'src/afsd_scaffold/vtk_legacy.py',
)
CASES = {'elastic_n3_i8': ('elastic', 3, 8),
         'j2_n3_i8': ('j2', 3, 8), 'j2_n3_i16': ('j2', 3, 16)}
UNSOLVED = {'friction', 'thermal', 'inertia', 'FKS_contact',
            'finite_strain_density', 'rotating_or_multiface_3D_contact',
            'AFSD', 'production_restart'}
COMPLETE = 'UNILATERAL_PATH_COMPLETE'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def safe_file(relative):
    path = (ROOT/relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise RuntimeError('missing/unsafe evidence path: '+str(relative))
    return path


def test_log(relative, expected=None, failed=False):
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8')
    counts = re.findall(r'^Ran (\d+) tests? in ', text, re.MULTILINE)
    if (len(counts) != 1 or int(counts[0]) <= 0
            or (expected is not None and int(counts[0]) != expected)):
        raise RuntimeError('unexpected test inventory: '+relative)
    # A legitimate outer timing footer may follow unittest's final status.
    successes = re.findall(r'^OK$', text, re.MULTILINE)
    failures = re.findall(r'^FAILED \(([^)]+)\)$', text, re.MULTILINE)
    if failed:
        if successes or failures != ['failures=1']:
            raise RuntimeError('initial regression failure lost/reclassified: '+relative)
    elif len(successes) != 1 or failures or re.search(r'^FAILED\b', text, re.MULTILINE):
        raise RuntimeError('test result is not PASS: '+relative)
    return dict(status='FAIL_PRESERVED' if failed else 'PASS',
                tests_run=int(counts[0]), failures=int(failed), errors=0,
                log=relative, log_sha256=sha(path))


def check_evidence():
    base = ROOT/'outputs/process_contact_q3_20261001'
    cfg = read(safe_file('configs/benchmarks/process_contact_q3.json'))
    if (cfg['production_authorized'] is not False
            or set(cfg['not_solved']) != UNSOLVED
            or cfg['parameters']['phase_duration']['value'] != .1
            or cfg['parameters']['peak_wall_displacement']['value'] != -4e-5
            or {tuple(case) for case in cfg['cases']} != set(CASES.values())):
        raise RuntimeError('Q3 scope/physical program/case inventory changed')
    runs = {}
    for name, (law, n, increments) in CASES.items():
        directory = base/name
        r = read(directory/'report.json')
        rb_path = base/(name+'_readback.json')
        rb = read(rb_path)
        if (r['status'] not in (COMPLETE, 'BLOCKED')
                or r['selected_law'] != law or r['particles'] != n**3
                or r['increments_per_phase'] != increments
                or r['production_authorized'] is not False
                or r['classification'] != cfg['classification']
                or set(r['not_solved']) != UNSOLVED
                or r['requested_end_s'] != .2):
            raise RuntimeError('case status/scope mismatch: '+name)
        rows = r['records']
        if (not rows or len(rows) != r['accepted_steps']+1
                or rows[0]['physical_time_s'] != 0.
                or rows[-1]['physical_time_s'] != r['accepted_time_s']
                or [row['step'] for row in rows] != list(range(len(rows)))
                or any(b['physical_time_s'] <= a['physical_time_s']
                       for a, b in zip(rows, rows[1:]))):
            raise RuntimeError('accepted-record inventory/time mismatch: '+name)
        if r['status'] == COMPLETE:
            # Accepted release-event substeps are additional, not replacements
            # for the requested end points; never assume exactly 2*increments.
            requested = [i*.2/(2*increments) for i in range(2*increments+1)]
            if (r['accepted_steps'] < 2*increments or r['accepted_time_s'] != .2
                    or r['error'] is not None
                    or any(not any(abs(row['physical_time_s']-t) <= 1e-14
                                   for row in rows) for t in requested)):
                raise RuntimeError('completion lacks full requested/accepted path: '+name)
            expected_readback = 'PASS_COMPLETE_READBACK'
        else:
            if not r['error'] or not 0 <= r['accepted_time_s'] <= .2:
                raise RuntimeError('BLOCKED endpoint/error mismatch: '+name)
            rejected = read(directory/'REJECTED_manifest.json')
            if (rejected['classification'] != 'REJECTED'
                    or rejected['restart_allowed'] is not False
                    or rejected['last_accepted_step'] != r['accepted_steps']
                    or rejected['error'] != r['error']
                    or 'REJECTED_manifest.json' not in r['artifact_sha256']):
                raise RuntimeError('rejected state lost/promoted: '+name)
            expected_readback = 'PASS_PARTIAL_READBACK'
        for event in r['events']:
            step = event['step']
            if (event['kind'] != 'RELEASE_LIMIT' or not 0 < step < len(rows)
                    or rows[step]['event_kind'] != 'RELEASE_LIMIT'
                    or rows[step]['physical_time_s'] != event['physical_time_s']):
                raise RuntimeError('release-event ledger mismatch: '+name)
        if not {'history.npz', 'checkpoint.npz', 'final.vtk', 'resolved_config.json'} <= set(r['artifact_sha256']):
            raise RuntimeError('artifact hash inventory incomplete: '+name)
        for relative, digest in r['source_sha256'].items():
            if sha(safe_file(relative)) != digest:
                raise RuntimeError('source changed after Q3 run: '+relative)
        for relative, digest in r['artifact_sha256'].items():
            path = (directory/relative).resolve()
            if not path.is_relative_to(directory.resolve()) or not path.is_file() or sha(path) != digest:
                raise RuntimeError('missing/changed/unsafe artifact: '+name+'/'+relative)
        if (rb['status'] != expected_readback or rb['completion_status'] != r['status']
                or rb['accepted_time_s'] != r['accepted_time_s']
                or rb['accepted_steps'] != r['accepted_steps'] or rb['law'] != law
                or rb['events'] != r['events'] or rb['production_authorized'] is not False
                or rb['final'] != rows[-1]):
            raise RuntimeError('readback does not match accepted/blocked run: '+name)
        runs[name] = dict(status=r['status'], error=r['error'], selected_law=law,
            particles=r['particles'], accepted_steps=r['accepted_steps'],
            accepted_benchmark_time_s=r['accepted_time_s'], requested_end_s=.2,
            release_events=len(r['events']), wall_before_report_s=r['wall_before_report_s'],
            readback_status=rb['status'], report_sha256=sha(directory/'report.json'),
            readback_sha256=sha(rb_path))
    tests = {
        'contact_mechanics_and_transactions': test_log(LOGS+'/q3_tests_final.log', 13),
        'saved_evidence_readback': test_log(LOGS+'/readback_tests.log'),
        'reused_q2_compatibility_not_new_tests': test_log(LOGS+'/q2_compatibility.log', 23),
        'core_before_config_guard_regression': test_log(LOGS+'/q3_tests_fixed.log', 12),
        'original_atomic_commit_failure': test_log(LOGS+'/q3_tests_initial.log', 12, failed=True),
        'full_prior_suite_this_round': 'NOT_RUN; previous PASS/FAIL evidence preserved',
        'afsd_or_three_condition_production_this_round': 'NOT_RUN',
        'contact_spatial_sensitivity_n5': 'NOT_RUN',
        'thermal_friction_FKS_or_rotating_multiface_contact': 'NOT_RUN',
        'gpu_parallel_optimization': 'NOT_RUN',
    }
    return runs, tests


def write_new(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    if TARGET.exists() or REPORT.exists():
        raise FileExistsError('never overwrite freeze/preservation evidence')
    old = read(PRIOR)
    if len(old['files_sha256']) != PARENT_COUNT:
        raise RuntimeError('unexpected parent inventory; expected 1193 entries')
    parent = read(ROOT/'logs/process_plastic_q2_20261001/freeze_final.log')
    if parent['status'] != 'PASS' or parent['files_checked'] != PARENT_COUNT:
        raise RuntimeError('parent final integrity verification is not PASS1193')
    hashes = {relative: sha(safe_file(relative)) for relative in old['files_sha256']}
    changed = sorted(relative for relative, digest in old['files_sha256'].items()
                     if hashes[relative] != digest)
    if set(changed)-ALLOWED:
        raise RuntimeError('out-of-scope historical changes: '+str(sorted(set(changed)-ALLOWED)))
    files = {safe_file(relative) for relative in old['files_sha256']}
    files.add(PRIOR)
    files.update(safe_file(relative) for relative in REQUIRED)
    for name in (LOGS, 'outputs/process_contact_q3_20261001'):
        directory = ROOT/name
        if not directory.is_dir():
            raise RuntimeError('missing Q3 evidence directory: '+name)
        for path in directory.rglob('*'):
            if path.is_file() and path.relative_to(ROOT).as_posix() not in EXCLUDED:
                files.add(safe_file(path.relative_to(ROOT).as_posix()))
    runs, tests = check_evidence()
    preservation = dict(status='PASS',
        scope='byte integrity and recorded evidence only; NOT scientific acceptance',
        previous_entries=PARENT_COUNT, previous_freeze=PRIOR.relative_to(ROOT).as_posix(),
        previous_freeze_sha256=sha(PRIOR), authorized_modified_previous_entries=changed,
        authorized_changes={name: {'before': old['files_sha256'][name], 'after': hashes[name]}
                            for name in changed},
        unchanged_previous_entries=PARENT_COUNT-len(changed),
        all_old_solver_sources_configs_tests_and_failure_artifacts_preserved=True,
        old_failure_records_preserved=True, production_authorized=False,
        benchmark_scope='ENGINEERING_APPROXIMATE single frictionless platen/roller-supported small-strain cube; NOT AFSD',
        case_results=runs, not_solved=sorted(UNSOLVED),
        separated_cap_is_not_prescribed_material_motion=True,
        full_free_body_or_rotating_multiface_contact_validated=False,
        plastic_work_is_solved_heat=False,
        timing_note='wall_before_report_s excludes report serialization and interpreter/import startup')
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    write_new(REPORT, preservation)
    files.add(REPORT)
    result = dict(schema_version='1.0', date='2026-10-01',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False,
        scientific_status='ENGINEERING_APPROXIMATE single frictionless platen benchmark; per-case COMPLETE/BLOCKED; NOT AFSD validation',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=preservation, actual_tests=tests,
        current_new_passing_test_count=sum(tests[key]['tests_run'] for key in
            ('contact_mechanics_and_transactions', 'saved_evidence_readback')),
        original_core_test_failures_preserved=1,
        continuum_accuracy_is_admission_gate=False,
        excluded_self_referential_or_future_paths=sorted(EXCLUDED),
        files_sha256={path.relative_to(ROOT).as_posix(): sha(path) for path in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=preservation,
                         production_authorized=False), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

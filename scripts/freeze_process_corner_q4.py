"""Archive Q4 evidence after reports/readbacks are final; never run physics.

This checks byte integrity and recorded status, not scientific acceptance.
BLOCKED cases and any failed/incomplete test attempts remain distinguishable.
An existing freeze or preservation report is never overwritten.
"""
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_CONTACT_Q3_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_CORNER_Q4_FREEZE.json'
LOGS = 'logs/process_corner_q4_20261001'
OUTPUTS = 'outputs/process_corner_q4_20261001'
REPORT = ROOT/LOGS/'preservation_report.json'
PARENT_COUNT = 1236
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {
    'docs/PROCESS_CORNER_Q4_FREEZE.json',
    LOGS+'/freeze_create.log',
    LOGS+'/freeze_final.log',
}
REQUIRED = (
    'docs/PROCESS_CORNER_Q4_CONTRACT.md',
    'docs/PROCESS_CORNER_Q4_REFERENCES.md',
    'docs/PROCESS_CORNER_Q4_REPORT.md',
    'configs/benchmarks/process_corner_q4.json',
    'scripts/probe_process_corner.py',
    'scripts/inspect_process_corner.py',
    'scripts/freeze_process_corner_q4.py',
    'tests/review/test_process_corner_q4.py',
    'tests/review/test_process_corner_readback.py',
    LOGS+'/COMMANDS.md',
    LOGS+'/q4_tests_initial.log',
    'logs/process_contact_q3_20261001/freeze_final.log',
    # Direct run reports omit some transitive inputs; preserve them as well.
    'configs/benchmarks/process_contact_q3.json',
    'configs/benchmarks/process_plastic_q2.json',
    'configs/benchmarks/process_time_q1.json',
    'configs/benchmarks/engineering_block_e2.json',
    'src/afsd_backend/pysph/kernel.py',
    'src/afsd_backend/pysph/output.py',
    'src/afsd_scaffold/vtk_legacy.py',
)
CASES = {
    'elastic_staggered_i8': ('elastic', 'staggered', 8),
    'j2_symmetric_i8': ('j2', 'symmetric', 8),
    'j2_staggered_i8': ('j2', 'staggered', 8),
    'j2_staggered_i16': ('j2', 'staggered', 16),
}
UNSOLVED = {'friction', 'thermal', 'inertia', 'FKS_contact',
            'finite_strain_density', 'arbitrary_rotating_tool_contact',
            'AFSD', 'production_restart'}
COMPLETE = 'CORNER_PATH_COMPLETE'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def safe_file(relative):
    path = (ROOT/relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise RuntimeError('missing/unsafe evidence path: '+str(relative))
    return path


def test_log(relative, expected=None, require_pass=True):
    """Count actual tests; do not invent or require an initial failure."""
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8-sig')
    counts = re.findall(r'^Ran (\d+) tests? in ', text, re.MULTILINE)
    successes = re.findall(r'^OK$', text, re.MULTILINE)
    failures = re.findall(r'^FAILED \(([^)]+)\)$', text, re.MULTILINE)
    count = int(counts[0]) if len(counts) == 1 else None
    if require_pass:
        if (count is None or count <= 0 or len(successes) != 1 or failures
                or re.search(r'^FAILED\b', text, re.MULTILINE)
                or (expected is not None and count != expected)):
            raise RuntimeError('test inventory/result is not required PASS: '+relative)
        status = 'PASS'
    elif count and len(successes) == 1 and not failures:
        status = 'PASS'
    elif count and not successes and len(failures) == 1:
        status = 'FAIL_PRESERVED'
    else:
        status = 'INCOMPLETE_OR_UNPARSED_PRESERVED'
    return dict(status=status, tests_run=count,
                passing_tests=count if status == 'PASS' else None,
                unittest_failure_detail=failures,
                log=relative, log_sha256=sha(path))


def core_test_log():
    """Preserve the actual known strict reference failure, not an all-green gate."""
    relative = LOGS+'/q4_tests_final.log'
    result = test_log(relative, require_pass=False)
    text = safe_file(relative).read_text(encoding='utf-8-sig')
    name = 'test_12_j2_staggered_shared_state_partial_release_and_free_residual_shape'
    failed_names = re.findall(r'^FAIL: (\S+) ', text, re.MULTILINE)
    successes = re.findall(r'^test_\S+ .* \.\.\. ok$', text, re.MULTILINE)
    if (result['status'] != 'FAIL_PRESERVED' or result['tests_run'] != 13
            or result['unittest_failure_detail'] != ['failures=1']
            or failed_names != [name] or len(successes) != 12
            or re.search(r'^ERROR:', text, re.MULTILINE)
            or 'np.testing.assert_allclose(b.pa.ep, ep, atol=1e-12, rtol=0.)' not in text):
        raise RuntimeError('core outcome differs from explicit known strict failure; inspect without reclassifying')
    result.update(passing_tests=12, failing_tests=1, errors=0,
                  retained_strict_failure=name,
                  failure_scope='post-unloading plastic-strain constancy at absolute 1e-12; not waived or loosened',
                  later_assertions_in_failed_test='NOT_RUN after first failing assertion',
                  status='FAIL_PRESERVED_NOT_SCIENTIFIC_ACCEPTANCE')
    return result


def check_evidence():
    base = ROOT/OUTPUTS
    cfg = read(safe_file('configs/benchmarks/process_corner_q4.json'))
    expected_parameters = dict(length=.01, wall_minimum=-4e-5,
                               peak_time=.1, end_time=.2,
                               staggered_y_return_start=.15)
    if (cfg['production_authorized'] is not False
            or set(cfg['not_solved']) != UNSOLVED
            or cfg['wall_axes'] != [0, 1]
            or cfg['wall_normals'] != [[-1, 0, 0], [0, -1, 0]]
            or {tuple(case) for case in cfg['cases']} != set(CASES.values())
            or any(cfg['parameters'][key]['value'] != value
                   for key, value in expected_parameters.items())):
        raise RuntimeError('Q4 scope/physical program/case inventory changed')
    runs = {}
    for name, (law, schedule, increments) in CASES.items():
        directory = base/name
        r = read(safe_file(OUTPUTS+'/'+name+'/report.json'))
        rb_path = safe_file(OUTPUTS+'/'+name+'_readback.json')
        rb = read(rb_path)
        if (r['status'] not in (COMPLETE, 'BLOCKED')
                or r['selected_law'] != law or r['schedule'] != schedule
                or r['particles'] != 27 or r['contact_rows'] != 18
                or r['permanent_free_dofs'] != 54
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
            # Accepted release substeps are additional, not replacements for
            # requested end points. Do not demand exactly 2*increments steps.
            requested = [i*.2/(2*increments) for i in range(2*increments+1)]
            if (r['accepted_steps'] < 2*increments or r['accepted_time_s'] != .2
                    or r['error'] is not None
                    or any(not any(abs(row['physical_time_s']-t) <= 1e-14
                                   for row in rows) for t in requested)):
                raise RuntimeError('completion lacks requested/accepted path: '+name)
            expected_readback = 'PASS_COMPLETE_READBACK'
        else:
            if not r['error'] or not 0 <= r['accepted_time_s'] <= .2:
                raise RuntimeError('BLOCKED endpoint/error mismatch: '+name)
            rejected = read(safe_file(OUTPUTS+'/'+name+'/REJECTED_manifest.json'))
            if (rejected['classification'] != 'REJECTED'
                    or rejected['restart_allowed'] is not False
                    or rejected['last_accepted_step'] != r['accepted_steps']
                    or rejected['error'] != r['error']
                    or 'REJECTED_manifest.json' not in r['artifact_sha256']):
                raise RuntimeError('rejected state lost/promoted: '+name)
            expected_readback = 'PASS_PARTIAL_READBACK'
        for event in r['events']:
            step = event['step']
            bracket = event['bracket_s']
            contact_rows = event['near_zero_contact_rows']
            features = event['near_zero_feature_ids']
            if (event['kind'] != 'RELEASE_LIMIT' or not 0 < step < len(rows)
                    or rows[step]['event_kind'] != 'RELEASE_LIMIT'
                    or rows[step]['physical_time_s'] != event['physical_time_s']
                    or len(bracket) != 2 or bracket[0] != event['physical_time_s']
                    or not 0 <= bracket[1]-bracket[0] <= 1e-12
                    or any(not 0 <= contact_row < 18 for contact_row in contact_rows)
                    or features != sorted({contact_row//9+1 for contact_row in contact_rows})):
                raise RuntimeError('release-event identity/ledger mismatch: '+name)
        required_artifacts = {'history.npz', 'checkpoint.npz', 'final.vtk', 'resolved_config.json'}
        if not required_artifacts <= set(r['artifact_sha256']):
            raise RuntimeError('artifact hash inventory incomplete: '+name)
        for relative, digest in r['source_sha256'].items():
            if sha(safe_file(relative)) != digest:
                raise RuntimeError('source changed after Q4 run: '+relative)
        for relative, digest in r['artifact_sha256'].items():
            path = (directory/relative).resolve()
            if not path.is_relative_to(directory.resolve()) or not path.is_file() or sha(path) != digest:
                raise RuntimeError('missing/changed/unsafe artifact: '+name+'/'+relative)
        if (rb['status'] != expected_readback or rb['completion_status'] != r['status']
                or rb['accepted_time_s'] != r['accepted_time_s']
                or rb['accepted_steps'] != r['accepted_steps'] or rb['law'] != law
                or rb['schedule'] != schedule or rb['events'] != r['events']
                or rb['production_authorized'] is not False or rb['final'] != rows[-1]):
            raise RuntimeError('readback does not match accepted/blocked run: '+name)
        runs[name] = dict(status=r['status'], error=r['error'], selected_law=law,
            schedule=schedule, particles=r['particles'], contact_rows=r['contact_rows'],
            accepted_steps=r['accepted_steps'], accepted_benchmark_time_s=r['accepted_time_s'],
            requested_end_s=.2, release_events=len(r['events']),
            wall_before_report_s=r['wall_before_report_s'], readback_status=rb['status'],
            report_sha256=sha(directory/'report.json'), readback_sha256=sha(rb_path))
    prior_attempts = {
        path.name: test_log(path.relative_to(ROOT).as_posix(), require_pass=False)
        for path in sorted((ROOT/LOGS).glob('q4_tests*.log'))
        if path.name != 'q4_tests_final.log'
    }
    tests = {
        'corner_mechanics_and_transactions': core_test_log(),
        'saved_evidence_readback': test_log(LOGS+'/readback_tests.log'),
        'initial_readback_failure_preserved': test_log(LOGS+'/readback_tests_initial.log', require_pass=False),
        'reused_q3_compatibility_not_new_tests': test_log(LOGS+'/q3_compatibility.log', 19),
        'prior_q4_test_attempts_preserved': prior_attempts,
        'full_prior_suite_this_round': 'NOT_RUN; previous PASS/FAIL evidence preserved',
        'afsd_or_three_condition_production_this_round': 'NOT_RUN',
        'contact_spatial_sensitivity_n5': 'NOT_RUN',
        'thermal_friction_FKS_or_rotating_tool_contact': 'NOT_RUN',
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
        raise RuntimeError('unexpected parent inventory; expected 1236 entries')
    parent = read(safe_file('logs/process_contact_q3_20261001/freeze_final.log'))
    if parent['status'] != 'PASS' or parent['files_checked'] != PARENT_COUNT:
        raise RuntimeError('parent final integrity verification is not PASS1236')
    hashes = {relative: sha(safe_file(relative)) for relative in old['files_sha256']}
    changed = sorted(relative for relative, digest in old['files_sha256'].items()
                     if hashes[relative] != digest)
    if set(changed)-ALLOWED:
        raise RuntimeError('out-of-scope historical changes: '+str(sorted(set(changed)-ALLOWED)))
    files = {safe_file(relative) for relative in old['files_sha256']}
    files.add(PRIOR)
    files.update(safe_file(relative) for relative in REQUIRED)
    for name in (LOGS, OUTPUTS):
        directory = ROOT/name
        if not directory.is_dir():
            raise RuntimeError('missing Q4 evidence directory: '+name)
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
        benchmark_scope='ENGINEERING_APPROXIMATE two orthogonal frictionless walls/roller-supported small-strain cube; NOT AFSD',
        case_results=runs, not_solved=sorted(UNSOLVED),
        shared_material_commit_and_per_wall_ledgers=True,
        separated_material_motion_not_prescribed=True,
        arbitrary_rotating_tool_contact_validated=False,
        plastic_work_is_solved_heat=False,
        timing_note='wall_before_report_s excludes report serialization and interpreter/import startup')
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    write_new(REPORT, preservation)
    files.add(REPORT)
    result = dict(schema_version='1.0', date='2026-10-01',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        production_authorized=False,
        scientific_status='ENGINEERING_APPROXIMATE orthogonal frictionless contact benchmark; per-case COMPLETE/BLOCKED; NOT AFSD validation',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=preservation, actual_tests=tests,
        current_new_passing_test_count=sum(tests[key]['passing_tests'] for key in
            ('corner_mechanics_and_transactions', 'saved_evidence_readback')),
        current_new_strict_failing_test_count=tests['corner_mechanics_and_transactions']['failing_tests'],
        failed_or_incomplete_attempts_not_reclassified=True,
        continuum_accuracy_is_admission_gate=False,
        excluded_self_referential_or_future_paths=sorted(EXCLUDED),
        files_sha256={path.relative_to(ROOT).as_posix(): sha(path) for path in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=preservation,
                         production_authorized=False), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

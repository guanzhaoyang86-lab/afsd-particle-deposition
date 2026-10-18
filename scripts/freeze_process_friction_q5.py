"""Archive Q5 evidence without running physics or granting scientific acceptance.

Run only after reports are final. COMPLETE, BLOCKED, explicit NOT_RUN, and
failed/incomplete test attempts remain distinct. Neither output is overwritten.
Extra companion evidence outside the Q5 directories can be named with --include.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_CORNER_Q4_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_FRICTION_Q5_FREEZE.json'
LOGS = 'logs/process_friction_q5_20261001'
OUTPUTS = 'outputs/process_friction_q5_20261001'
REPORT = ROOT/LOGS/'preservation_report.json'
PARENT_COUNT = 1289
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {'docs/PROCESS_FRICTION_Q5_FREEZE.json',
            LOGS+'/freeze_create.log', LOGS+'/freeze_final.log'}
REQUIRED = (
    'docs/PROCESS_FRICTION_Q5_CONTRACT.md',
    'docs/PROCESS_FRICTION_Q5_REFERENCES.md',
    'docs/PROCESS_FRICTION_Q5_REPORT.md',
    'configs/benchmarks/process_friction_q5.json',
    'scripts/probe_process_friction.py',
    'scripts/inspect_process_friction.py',
    'scripts/freeze_process_friction_q5.py',
    'tests/review/test_process_friction_q5.py',
    'tests/review/test_process_friction_readback.py',
    'tests/review/test_process_friction_guards.py',
    LOGS+'/COMMANDS.md',
    'logs/process_corner_q4_20261001/freeze_final.log',
    # Transitive inputs not all present in the direct run's source manifest.
    'configs/benchmarks/process_contact_q3.json',
    'configs/benchmarks/process_plastic_q2.json',
    'configs/benchmarks/process_time_q1.json',
    'configs/benchmarks/engineering_block_e2.json',
    'src/afsd_backend/pysph/kernel.py',
    'src/afsd_backend/pysph/output.py',
    'src/afsd_scaffold/vtk_legacy.py',
)
CASES = {'n3_i4': 4, 'n3_i8': 8}
UNSOLVED = {'isotropic_two_tangent_friction', 'normal_release_recontact',
            'plasticity_with_friction', 'thermal', 'inertia', 'rotating_tool',
            'AFSD', 'production_restart', 'spatial_convergence'}
COMPLETE = 'DIRECTIONAL_FRICTION_PATH_COMPLETE'
KNOWN_Q4_FAILURE = 'test_12_j2_staggered_shared_state_partial_release_and_free_residual_shape'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def safe_file(relative):
    path = (ROOT/relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise RuntimeError('missing/unsafe evidence path: '+str(relative))
    return path


def not_run(relative):
    """Absence alone is not evidence that a test or case was NOT_RUN."""
    path = safe_file(relative)
    record = read(path)
    if (record.get('status') != 'NOT_RUN' or not record.get('reason')
            or record.get('production_authorized') is not False):
        raise RuntimeError('invalid explicit NOT_RUN evidence: '+relative)
    return dict(status='NOT_RUN', reason=record['reason'],
                evidence=relative, evidence_sha256=sha(path),
                production_authorized=False)


def test_log(relative):
    """Parse actual unittest outcomes; never convert a failed run to PASS."""
    path = ROOT/relative
    if not path.is_file():
        result = not_run(str(Path(relative).with_suffix('.NOT_RUN.json')).replace('\\', '/'))
        result.update(tests_run=0, passing_tests=0, failing_tests=0, errors=0)
        return result
    path = safe_file(relative)
    text = path.read_text(encoding='utf-8-sig')
    runs = re.findall(r'^Ran (\d+) tests? in ', text, re.MULTILINE)
    terminals = re.findall(r'^(OK(?: \([^\r\n]+\))?|FAILED \([^\r\n]+\))$', text, re.MULTILINE)
    verbose = re.findall(r'^test_\S+ .* \.\.\. (ok|FAIL|ERROR|skipped[^\r\n]*|expected failure|unexpected success)$',
                         text, re.MULTILINE)
    count = int(runs[0]) if len(runs) == 1 else None
    complete = count is not None and len(terminals) == 1
    passed = verbose.count('ok')
    failed = verbose.count('FAIL')
    errors = verbose.count('ERROR')
    skipped = sum(line.startswith('skipped') for line in verbose)
    expected = verbose.count('expected failure')
    unexpected = verbose.count('unexpected success')
    detailed = count is not None and len(verbose) == count
    if complete and terminals[0].startswith('OK') and not failed and not errors and not unexpected:
        status = 'PASS' if not skipped and not expected else 'COMPLETED_WITH_SKIP_OR_EXPECTED_FAILURE'
    elif complete and terminals[0].startswith('FAILED'):
        status = 'FAIL_PRESERVED'
    else:
        status = 'INCOMPLETE_OR_UNPARSED_PRESERVED'
    return dict(status=status, tests_run=count,
                passing_tests=passed if detailed else None,
                failing_tests=failed if detailed else None, errors=errors if detailed else None,
                skipped_tests=skipped if detailed else None,
                expected_failures=expected if detailed else None,
                unexpected_successes=unexpected if detailed else None,
                verbose_inventory_complete=detailed,
                terminal_lines=terminals,
                failed_test_names=re.findall(r'^FAIL: (\S+) ', text, re.MULTILINE),
                log=relative, log_sha256=sha(path))


def compatibility_log():
    relative = LOGS+'/q4_compatibility.log'
    result = test_log(relative)
    if result['status'] == 'NOT_RUN':
        return result
    text = safe_file(relative).read_text(encoding='utf-8-sig')
    expected = (result['status'] == 'FAIL_PRESERVED' and result['tests_run'] == 19
                and result['passing_tests'] == 18 and result['failing_tests'] == 1
                and result['errors'] == 0 and result['failed_test_names'] == [KNOWN_Q4_FAILURE]
                and 'np.testing.assert_allclose(b.pa.ep, ep, atol=1e-12, rtol=0.)' in text)
    result.update(matches_previous_18_pass_1_strict_failure=expected,
                  is_new_test_inventory=False,
                  later_assertions_in_failed_test='NOT_RUN after the failing assertion')
    # A changed result is evidence too, but cannot quietly be labelled unchanged.
    if not expected:
        result['compatibility_caution'] = 'DIFFERS_FROM_PREVIOUS_RESULT; inspect before scientific conclusions'
    return result


def check_case(name, increments, cfg):
    relative = OUTPUTS+'/'+name
    directory = ROOT/relative
    if not (directory/'report.json').is_file():
        if any((directory/name).exists() for name in
               ('history.npz', 'checkpoint.npz', 'final.vtk', 'REJECTED_manifest.json')):
            raise RuntimeError('partial run artifacts cannot be classified NOT_RUN: '+name)
        return not_run(relative+'/NOT_RUN.json')
    if (directory/'NOT_RUN.json').exists():
        raise RuntimeError('run report conflicts with NOT_RUN marker: '+name)
    r = read(safe_file(relative+'/report.json'))
    if (r['status'] not in (COMPLETE, 'BLOCKED') or r['particles'] != 27
            or r['contact_rows'] != 9 or r['permanent_free_dofs'] != 54
            or r['increments_per_phase'] != increments or r['requested_end_s'] != .2
            or r['production_authorized'] is not False
            or r['classification'] != cfg['classification'] or set(r['not_solved']) != UNSOLVED):
        raise RuntimeError('Q5 run scope/inventory mismatch: '+name)
    rows = r['records']
    if (not rows or len(rows) != r['accepted_steps']+1
            or [row['step'] for row in rows] != list(range(len(rows)))
            or rows[0]['physical_time_s'] != 0.
            or rows[-1]['physical_time_s'] != r['accepted_time_s']
            or not 0 <= r['accepted_time_s'] <= .2
            or any(b['physical_time_s'] <= a['physical_time_s'] for a, b in zip(rows, rows[1:]))):
        raise RuntimeError('accepted state ledger mismatch: '+name)
    if r['status'] == COMPLETE:
        if (r['error'] is not None or r['accepted_time_s'] != .2
                or r['accepted_steps'] < 4*increments
                or any(not any(abs(row['physical_time_s']-i*.2/(4*increments)) <= 1e-14
                               for row in rows) for i in range(4*increments+1))):
            raise RuntimeError('completion lacks full requested path: '+name)
        expected_readback = 'PASS_COMPLETE_READBACK'
    else:
        rejected = read(safe_file(relative+'/REJECTED_manifest.json'))
        if (not r['error'] or rejected['classification'] != 'REJECTED'
                or rejected['restart_allowed'] is not False or rejected['error'] != r['error']
                or rejected['last_accepted_step'] != r['accepted_steps']
                or 'REJECTED_manifest.json' not in r['artifact_sha256']):
            raise RuntimeError('rejected state lost/promoted: '+name)
        expected_readback = 'PASS_PARTIAL_READBACK'
    events = r['events']
    if [event['step'] for event in events] != [i for i,row in enumerate(rows) if row['event_kind']=='FRICTION_LIMIT']:
        raise RuntimeError('event/accepted ledger inventory mismatch: '+name)
    for event in events:
        step, bracket = event['step'], event['bracket_s']
        if (event['kind'] != 'FRICTION_LIMIT' or not 0 < step < len(rows)
                or event['physical_time_s'] != rows[step]['physical_time_s']
                or len(bracket) != 2 or bracket[0] != event['physical_time_s']
                or not 0 <= bracket[1]-bracket[0] <= 1e-12
                or len(event['ending_interval_modes']) != 9
                or len(event['cone_margin_N']) != 9):
            raise RuntimeError('friction event identity mismatch: '+name)
    if not {'history.npz', 'checkpoint.npz', 'final.vtk', 'resolved_config.json'} <= set(r['artifact_sha256']):
        raise RuntimeError('missing accepted artifact inventory: '+name)
    for source, digest in r['source_sha256'].items():
        if sha(safe_file(source)) != digest:
            raise RuntimeError('source changed after Q5 run: '+source)
    for artifact, digest in r['artifact_sha256'].items():
        path = (directory/artifact).resolve()
        if not path.is_relative_to(directory.resolve()) or not path.is_file() or sha(path) != digest:
            raise RuntimeError('missing/changed/unsafe accepted artifact: '+name+'/'+artifact)
    rb_path = ROOT/OUTPUTS/(name+'_readback.json')
    if rb_path.is_file():
        rb = read(rb_path)
        if (rb['status'] != expected_readback or rb['completion_status'] != r['status']
                or rb['accepted_steps'] != r['accepted_steps'] or rb['accepted_time_s'] != r['accepted_time_s']
                or rb['events'] != events or rb['final'] != rows[-1] or rb['production_authorized'] is not False):
            raise RuntimeError('readback mismatch: '+name)
        readback = dict(status=rb['status'], evidence=rb_path.relative_to(ROOT).as_posix(), sha256=sha(rb_path))
    else:
        # Do not call absent readback PASS or silently infer NOT_RUN. Explicit
        # failed/blocked/not-run readback evidence may instead be archived.
        marker = safe_file(OUTPUTS+'/'+name+'_readback_status.json')
        rb = read(marker)
        if rb.get('status') not in ('FAIL', 'BLOCKED', 'NOT_RUN') or not rb.get('reason'):
            raise RuntimeError('readback needs an explicit outcome/reason: '+name)
        evidence = safe_file(rb['evidence']) if rb.get('evidence') else None
        if rb['status'] != 'NOT_RUN' and evidence is None:
            raise RuntimeError('failed readback needs its actual evidence path: '+name)
        readback = dict(status=rb['status'], reason=rb['reason'], evidence=marker.relative_to(ROOT).as_posix(),
                        sha256=sha(marker), failure_evidence=None if evidence is None else evidence.relative_to(ROOT).as_posix())
    return dict(status=r['status'], error=r['error'], particles=r['particles'], contact_rows=r['contact_rows'],
        accepted_steps=r['accepted_steps'], accepted_benchmark_time_s=r['accepted_time_s'], requested_end_s=.2,
        increments_per_phase=increments, friction_events=len(events), wall_before_report_s=r['wall_before_report_s'],
        selected_array_bytes=r['selected_array_bytes'], memory_note=r['memory_note'],
        report_sha256=sha(directory/'report.json'), readback=readback)


def check_evidence():
    cfg = read(safe_file('configs/benchmarks/process_friction_q5.json'))
    values = dict(mu=.1, normal_displacement=-4e-5, tangent_amplitude=2e-5,
                  phase_duration=.05, temperature=293.15)
    if (cfg['production_authorized'] is not False or set(cfg['not_solved']) != UNSOLVED
            or cfg['cases'] != [4, 8] or cfg['run_timeout_s'] != 60
            or any(cfg['parameters'][key]['value'] != value for key,value in values.items())):
        raise RuntimeError('Q5 declared scope/parameter inventory changed')
    runs = {name: check_case(name, increments, cfg) for name,increments in CASES.items()}
    final_logs = {'q5_tests_final.log', 'readback_tests.log', 'guard_tests.log', 'q4_compatibility.log'}
    attempts = {path.name: test_log(path.relative_to(ROOT).as_posix())
                for path in sorted((ROOT/LOGS).glob('*.log'))
                if path.name not in final_logs and ('test' in path.stem or 'independent' in path.stem)}
    tests = dict(friction_mechanics_and_transactions=test_log(LOGS+'/q5_tests_final.log'),
                 saved_evidence_readback=test_log(LOGS+'/readback_tests.log'),
                 actual_commit_guards=test_log(LOGS+'/guard_tests.log'),
                 reused_q4_compatibility_not_new_tests=compatibility_log(),
                 prior_q5_attempts_preserved=attempts,
                 full_prior_suite_this_round='NOT_RUN; all previous PASS/FAIL evidence preserved',
                 afsd_or_three_condition_production_this_round='NOT_RUN',
                 spatial_sensitivity_n5='NOT_RUN',
                 isotropic_friction_thermal_plastic_coupling_rotating_tool='NOT_RUN',
                 gpu_parallel_optimization='NOT_RUN')
    return runs, tests


def write_new(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include', action='append', default=[], help='additional repository-relative evidence file')
    args = parser.parse_args()
    if TARGET.exists() or REPORT.exists():
        raise FileExistsError('never overwrite freeze/preservation evidence')
    old = read(PRIOR)
    if len(old['files_sha256']) != PARENT_COUNT:
        raise RuntimeError('unexpected parent inventory; expected 1289 entries')
    final = read(safe_file('logs/process_corner_q4_20261001/freeze_final.log'))
    if final['status'] != 'PASS' or final['files_checked'] != PARENT_COUNT:
        raise RuntimeError('parent final integrity verification is not PASS1289')
    hashes = {name:sha(safe_file(name)) for name in old['files_sha256']}
    changed = sorted(name for name,digest in old['files_sha256'].items() if hashes[name] != digest)
    if set(changed)-ALLOWED:
        raise RuntimeError('out-of-scope historical changes: '+str(sorted(set(changed)-ALLOWED)))
    files = {safe_file(name) for name in old['files_sha256']}
    files.add(PRIOR)
    files.update(safe_file(name) for name in REQUIRED)
    for name in args.include:
        if name in EXCLUDED:
            raise RuntimeError('cannot include self/future output: '+name)
        files.add(safe_file(name))
    for name in (LOGS, OUTPUTS):
        directory = ROOT/name
        if not directory.is_dir():
            raise RuntimeError('missing Q5 evidence directory: '+name)
        for path in directory.rglob('*'):
            if path.is_file() and path.relative_to(ROOT).as_posix() not in EXCLUDED:
                files.add(safe_file(path.relative_to(ROOT).as_posix()))
    runs, tests = check_evidence()
    preservation = dict(status='PASS', scope='byte integrity and recorded evidence only; NOT scientific acceptance',
        previous_entries=PARENT_COUNT, previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        authorized_modified_previous_entries=changed,
        authorized_changes={name:dict(before=old['files_sha256'][name], after=hashes[name]) for name in changed},
        unchanged_previous_entries=PARENT_COUNT-len(changed),
        all_old_solver_sources_configs_tests_and_failure_artifacts_preserved=True,
        old_failure_records_preserved=True, production_authorized=False,
        benchmark_scope='ENGINEERING_APPROXIMATE single fixed-normal directional-Y Coulomb wall and clamped small-strain elastic cube; NOT AFSD',
        case_results=runs, not_solved=sorted(UNSOLVED), shared_material_commit=True,
        dissipation_is_solved_temperature=False, arbitrary_rotating_tool_contact_validated=False,
        timing_note='wall_before_report_s excludes report serialization and interpreter/import startup')
    counters = {}
    for field in ('passing_tests', 'failing_tests', 'errors'):
        counts = [tests[key].get(field) for key in
                  ('friction_mechanics_and_transactions', 'saved_evidence_readback', 'actual_commit_guards')]
        counters['current_new_'+field] = sum(counts) if all(isinstance(count, int) for count in counts) else None
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    write_new(REPORT, preservation)
    files.add(REPORT)
    result = dict(schema_version='1.0', date='2026-10-01',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'], production_authorized=False,
        scientific_status='ENGINEERING_APPROXIMATE directional friction benchmark; per-case COMPLETE/BLOCKED/NOT_RUN; NOT AFSD validation',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        preservation=preservation, actual_tests=tests, **counters,
        failed_or_incomplete_attempts_not_reclassified=True, continuum_accuracy_is_admission_gate=False,
        excluded_self_referential_or_future_paths=sorted(EXCLUDED),
        files_sha256={path.relative_to(ROOT).as_posix():sha(path) for path in sorted(files)})
    write_new(TARGET, result)
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=preservation,
                         actual_tests=tests, production_authorized=False), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

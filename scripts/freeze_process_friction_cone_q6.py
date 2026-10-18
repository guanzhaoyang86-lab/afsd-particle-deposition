"""Archive Q6 evidence, not scientific acceptance, without running physics.

The Q5 archive and every historical failure remain immutable except the three
handoff/index documents. COMPLETE, BLOCKED, explicit NOT_RUN and failed test
attempts are separate outcomes. Run only after final reports/readbacks exist.
Use --include for additional repository-relative evidence outside Q6 folders.
"""
import argparse
import json
from pathlib import Path
import re

from freeze_process_friction_q5 import sha, read, safe_file, not_run, test_log, write_new


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'docs/PROCESS_FRICTION_Q5_FREEZE.json'
TARGET = ROOT/'docs/PROCESS_FRICTION_CONE_Q6_FREEZE.json'
LOGS = 'logs/process_friction_cone_q6_20261001'
OUTPUTS = 'outputs/process_friction_cone_q6_20261001'
REPORT = ROOT/LOGS/'preservation_report.json'
PARENT_COUNT = 1327
ALLOWED = {'HANDOFF.md', 'docs/IMPLEMENTATION_PLAN.md', 'docs/OUTPUT_SCHEMA.md'}
EXCLUDED = {'docs/PROCESS_FRICTION_CONE_Q6_FREEZE.json',
            LOGS+'/freeze_create.log', LOGS+'/freeze_final.log'}
REQUIRED = (
    'docs/PROCESS_FRICTION_CONE_Q6_CONTRACT.md',
    'docs/PROCESS_FRICTION_CONE_Q6_REFERENCES.md',
    'docs/PROCESS_FRICTION_CONE_Q6_REPORT.md',
    'configs/benchmarks/process_friction_cone_q6.json',
    'scripts/probe_process_friction_cone.py',
    'scripts/inspect_process_friction_cone.py',
    'scripts/freeze_process_friction_cone_q6.py',
    'tests/review/test_process_friction_cone_q6.py',
    'tests/review/test_process_friction_cone_q6_guards.py',
    'tests/review/test_process_friction_cone_q6_readback.py',
    LOGS+'/COMMANDS.md',
    'logs/process_friction_q5_20261001/freeze_final.log',
)
CASES = {'oblique_i4': (4, 0.), 'oblique_i8': (8, 0.), 'basis37_i4': (4, 37.)}
UNSOLVED = {'normal_release_recontact', 'plasticity_with_friction', 'thermal',
            'inertia', 'rotating_normal_or_tool', 'arbitrary_curved_slip_path',
            'AFSD', 'production_restart', 'spatial_convergence'}
COMPLETE = 'TWO_TANGENT_CONE_PATH_COMPLETE'


def readback_result(name, report, expected_status):
    path = ROOT/OUTPUTS/(name+'_readback.json')
    if path.is_file():
        rb = read(path)
        if (rb['status'] != expected_status or rb['completion_status'] != report['status']
                or rb['accepted_steps'] != report['accepted_steps']
                or rb['accepted_time_s'] != report['accepted_time_s']
                or rb['events'] != report['events'] or rb['final'] != report['records'][-1]
                or rb['production_authorized'] is not False):
            raise RuntimeError('Q6 readback mismatch: '+name)
        return dict(status=rb['status'], evidence=path.relative_to(ROOT).as_posix(), sha256=sha(path))
    # Missing success evidence is not automatically NOT_RUN. A failed attempt
    # needs its actual log; an explicit NOT_RUN marker needs a reason.
    marker = safe_file(OUTPUTS+'/'+name+'_readback_status.json')
    rb = read(marker)
    if rb.get('status') not in ('FAIL', 'BLOCKED', 'NOT_RUN') or not rb.get('reason'):
        raise RuntimeError('readback needs explicit outcome/reason: '+name)
    evidence = safe_file(rb['evidence']) if rb.get('evidence') else None
    if rb['status'] != 'NOT_RUN' and evidence is None:
        raise RuntimeError('failed readback lacks actual evidence: '+name)
    return dict(status=rb['status'], reason=rb['reason'],
                evidence=marker.relative_to(ROOT).as_posix(), sha256=sha(marker),
                failure_evidence=None if evidence is None else evidence.relative_to(ROOT).as_posix())


def check_case(name, increments, angle, cfg):
    relative = OUTPUTS+'/'+name
    directory = ROOT/relative
    if not (directory/'report.json').is_file():
        if directory.exists() and any(p.name != 'NOT_RUN.json' for p in directory.iterdir()):
            raise RuntimeError('partial artifacts cannot be classified NOT_RUN: '+name)
        return not_run(relative+'/NOT_RUN.json')
    if (directory/'NOT_RUN.json').exists():
        raise RuntimeError('report conflicts with NOT_RUN marker: '+name)
    r = read(safe_file(relative+'/report.json'))
    if (r['status'] not in (COMPLETE, 'BLOCKED') or r['particles'] != 27
            or r['contact_rows'] != 9 or r['permanent_free_dofs'] != 54
            or r['increments_per_phase'] != increments or r['basis_angle_deg'] != angle
            or r['requested_end_s'] != .2 or r['production_authorized'] is not False
            or r['classification'] != cfg['classification'] or set(r['not_solved']) != UNSOLVED):
        raise RuntimeError('Q6 run scope/inventory mismatch: '+name)
    rows = r['records']
    if (not rows or len(rows) != r['accepted_steps']+1
            or [row['step'] for row in rows] != list(range(len(rows)))
            or rows[0]['physical_time_s'] != 0.
            or rows[-1]['physical_time_s'] != r['accepted_time_s']
            or not 0 <= r['accepted_time_s'] <= .2
            or any(b['physical_time_s'] <= a['physical_time_s'] for a,b in zip(rows, rows[1:]))):
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
    if [e['step'] for e in events] != [i for i,row in enumerate(rows) if row['event_kind']=='CONE_LIMIT']:
        raise RuntimeError('event/accepted ledger inventory mismatch: '+name)
    for event in events:
        step, bracket = event['step'], event['bracket_s']
        if (event['kind'] != 'CONE_LIMIT' or not 0 < step < len(rows)
                or event['physical_time_s'] != rows[step]['physical_time_s']
                or len(bracket) != 2 or bracket[0] != event['physical_time_s']
                or not 0 <= bracket[1]-bracket[0] <= 1e-12
                or len(event['ending_interval_modes']) != 9
                or any(mode not in (0, 1) for mode in event['ending_interval_modes'])
                or len(event['released_rows']) != 9
                or any(type(flag) is not bool for flag in event['released_rows'])
                or len(event['cone_margin_N']) != 9):
            raise RuntimeError('cone event identity mismatch: '+name)
    if not {'history.npz', 'checkpoint.npz', 'final.vtk', 'resolved_config.json'} <= set(r['artifact_sha256']):
        raise RuntimeError('missing accepted artifact inventory: '+name)
    for source, digest in r['source_sha256'].items():
        if sha(safe_file(source)) != digest:
            raise RuntimeError('source changed after Q6 run: '+source)
    for artifact, digest in r['artifact_sha256'].items():
        path = (directory/artifact).resolve()
        if not path.is_relative_to(directory.resolve()) or not path.is_file() or sha(path) != digest:
            raise RuntimeError('missing/changed/unsafe accepted artifact: '+name+'/'+artifact)
    return dict(status=r['status'], error=r['error'], particles=r['particles'], contact_rows=r['contact_rows'],
                accepted_steps=r['accepted_steps'], accepted_benchmark_time_s=r['accepted_time_s'],
                requested_end_s=.2, increments_per_phase=increments, basis_angle_deg=angle,
                cone_events=len(events), wall_before_report_s=r['wall_before_report_s'],
                selected_array_bytes=r['selected_array_bytes'], memory_note=r['memory_note'],
                report_sha256=sha(directory/'report.json'),
                readback=readback_result(name, r, expected_readback))


def check_evidence(old):
    cfg = read(safe_file('configs/benchmarks/process_friction_cone_q6.json'))
    values = dict(mu=.1, normal_displacement=-4e-5, tangent_amplitude=2e-5,
                  tangent_direction_global_yz=[.6, .8], phase_duration=.05, kappa=7e8)
    if (cfg['production_authorized'] is not False or set(cfg['not_solved']) != UNSOLVED
            or cfg['cases'] != [[4, 0], [8, 0], [4, 37]] or cfg['run_timeout_s'] != 60
            or any(cfg['parameters'][key]['value'] != value for key,value in values.items())):
        raise RuntimeError('Q6 declared scope/parameter inventory changed')
    runs = {name:check_case(name, increments, angle, cfg) for name,(increments,angle) in CASES.items()}
    # The final aggregate contains core, guards and reader tests. Do not add
    # earlier/focused attempts a second time to the new-test inventory.
    final = test_log(LOGS+'/q6_tests_all_final.log')
    compatibility = test_log(LOGS+'/q5_compatibility.log')
    attempts = {}
    for path in sorted((ROOT/LOGS).glob('*.log')):
        if path.name in ('q6_tests_all_final.log', 'q5_compatibility.log'):
            continue
        text = path.read_text(encoding='utf-8-sig')
        if (re.search(r'^Ran \d+ tests? in ', text, re.MULTILINE)
                or any(word in path.stem for word in ('test', 'core', 'guard', 'independent'))):
            attempts[path.name] = test_log(path.relative_to(ROOT).as_posix())
    prior_q4 = old['actual_tests']['reused_q4_compatibility_not_new_tests']
    evidence = safe_file(prior_q4['log'])
    if sha(evidence) != prior_q4['log_sha256']:
        raise RuntimeError('historical Q4 failure evidence changed')
    tests = dict(final_new_q6_inventory=final,
                 reused_q5_compatibility_not_new_tests=compatibility,
                 focused_or_initial_attempts_preserved_not_added_again=attempts,
                 historical_q4_result_not_rerun=dict(**prior_q4, rerun_this_round='NOT_RUN'),
                 full_prior_suite_this_round='NOT_RUN; previous PASS/FAIL evidence preserved',
                 afsd_or_three_condition_production_this_round='NOT_RUN',
                 spatial_sensitivity_n5='NOT_RUN',
                 friction_thermal_plastic_coupling_rotating_tool='NOT_RUN',
                 gpu_parallel_optimization='NOT_RUN')
    return runs, tests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include', action='append', default=[], help='additional repository-relative evidence file')
    args = parser.parse_args()
    if TARGET.exists() or REPORT.exists():
        raise FileExistsError('never overwrite freeze/preservation evidence')
    old = read(PRIOR)
    if len(old['files_sha256']) != PARENT_COUNT:
        raise RuntimeError('unexpected parent inventory; expected 1327 entries')
    final = read(safe_file('logs/process_friction_q5_20261001/freeze_final.log'))
    if final['status'] != 'PASS' or final['files_checked'] != PARENT_COUNT:
        raise RuntimeError('parent final integrity verification is not PASS1327')
    hashes = {name:sha(safe_file(name)) for name in old['files_sha256']}
    changed = sorted(name for name,digest in old['files_sha256'].items() if hashes[name] != digest)
    if set(changed)-ALLOWED:
        raise RuntimeError('out-of-scope historical changes: '+str(sorted(set(changed)-ALLOWED)))
    files = {safe_file(name) for name in old['files_sha256']}
    files.add(PRIOR)
    files.update(safe_file(name) for name in REQUIRED)
    for name in args.include:
        path = safe_file(name)
        if path.relative_to(ROOT).as_posix() in EXCLUDED:
            raise RuntimeError('cannot include self/future output: '+name)
        files.add(path)
    for name in (LOGS, OUTPUTS):
        directory = ROOT/name
        if not directory.is_dir():
            raise RuntimeError('missing Q6 evidence directory: '+name)
        for path in directory.rglob('*'):
            if path.is_file() and path.relative_to(ROOT).as_posix() not in EXCLUDED:
                files.add(safe_file(path.relative_to(ROOT).as_posix()))
    runs, tests = check_evidence(old)
    preservation = dict(status='PASS', scope='byte integrity and recorded evidence only; NOT scientific acceptance',
        previous_entries=PARENT_COUNT, previous_freeze=PRIOR.relative_to(ROOT).as_posix(), previous_freeze_sha256=sha(PRIOR),
        authorized_modified_previous_entries=changed,
        authorized_changes={name:dict(before=old['files_sha256'][name], after=hashes[name]) for name in changed},
        unchanged_previous_entries=PARENT_COUNT-len(changed),
        all_old_solver_sources_configs_tests_and_failure_artifacts_preserved=True,
        old_failure_records_preserved=True, production_authorized=False,
        benchmark_scope='ENGINEERING_APPROXIMATE single fixed-normal circular two-tangent Coulomb wall and clamped small-strain elastic cube; NOT AFSD',
        case_results=runs, not_solved=sorted(UNSOLVED), shared_material_commit=True,
        dissipation_is_solved_temperature=False, arbitrary_rotating_tool_contact_validated=False,
        timing_note='wall_before_report_s excludes report serialization and interpreter/import startup')
    counters = {'current_new_'+field:tests['final_new_q6_inventory'].get(field)
                for field in ('passing_tests', 'failing_tests', 'errors')}
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    write_new(REPORT, preservation)
    files.add(REPORT)
    result = dict(schema_version='1.0', date='2026-10-01',
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'], production_authorized=False,
        scientific_status='ENGINEERING_APPROXIMATE two-tangent Coulomb cone benchmark; per-case COMPLETE/BLOCKED/NOT_RUN; NOT AFSD validation',
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

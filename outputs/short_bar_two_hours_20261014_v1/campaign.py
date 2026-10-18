"""Finite CPU diagnostic continuation using only the frozen runner's gates.

This administrative module imports no physical solver. Authorization includes
preparation: UTC 2026-10-14 13:06:13 through 15:06:13, never a fresh two hours.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import json
import math
import os
import subprocess
import sys
import time
import traceback

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import supervise_afsd_deposition_20261007 as old

START_UTC = '2026-10-14T13:06:13Z'
DEADLINE_UTC = '2026-10-14T15:06:13Z'
START_EPOCH = datetime.fromisoformat(START_UTC.replace('Z', '+00:00')).timestamp()
DEADLINE_EPOCH = datetime.fromisoformat(DEADLINE_UTC.replace('Z', '+00:00')).timestamp()
PARENT = ROOT/'outputs/length_comparison_post_startup_20261013_v1/short_4p0mm/stage_007'
COST = ROOT/'outputs/length_comparison_20261013_v1/short_4p0mm/cost'
CONFIG = 'outputs/flow_diagnostic_20261013_v1/length_controls/short_4p0mm.json'
INITIAL_PARENT_SHA = 'a1843b8115e1d5ee182085b98c008b6d3a925ba462c4d760c035892b13d6ba45'
PYTHON = ROOT/'.envs/pysph/python.exe'
OUTPUT = BASE/'run'
OBSERVATION_TIME = .060
TARGET_TIME = .120439208984375
MAX_DT = .0005
MIN_DT = 1e-7
CHILD_WRAPPER_LIMIT = 650.
FINAL_RESERVE = 90.
TERMINATION_RESERVE = 35.
MAX_STAGES = 128
SCHEMA = 'FINITE_CPU_SHORT_BAR_TWO_HOURS_20261014_V1'


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def budget(initial_forecast, parent_actual):
    measured = max(float(initial_forecast), float(parent_actual))
    old.require(math.isfinite(measured) and measured > 0, 'finite positive measured candidate cost required')
    quota = min(27, math.floor(420./(1.5*measured)))
    while quota and 30.+1.5*quota*measured > 450.:
        quota -= 1
    return dict(measured_candidate_s=measured, candidate_quota=quota,
                max_macrosteps=quota//3, forecast_s=30.+1.5*quota*measured, runnable=quota >= 3)


def remaining_seconds(epoch_now, monotonic_now, monotonic_deadline):
    """Clock rollback cannot extend the authorization, and forward jumps stop it."""
    return min(DEADLINE_EPOCH-epoch_now, monotonic_deadline-monotonic_now)


def can_launch(remaining):
    return remaining >= CHILD_WRAPPER_LIMIT+FINAL_RESERVE+TERMINATION_RESERVE


def endpoint(current_time):
    return OBSERVATION_TIME if current_time < OBSERVATION_TIME-1e-13 else TARGET_TIME


def resume_dt(meta, attempts):
    """Undo only an evidenced endpoint truncation; preserve temporal refinement.

    The old runner stores a truncated last interval as its next proposal. Restore
    the already passed incoming proposal only when the whole stage has no failed
    attempts/refinement, the last pair is shorter solely due to its endpoint,
    and a preceding pair actually passed that incoming interval. Otherwise use
    exactly the runner's next proposal (including a failed final budget trial).
    """
    proposed = float(meta['next_proposed_macro_dt_s'])
    old.require(math.isfinite(proposed) and MIN_DT <= proposed <= MAX_DT, 'invalid runner next dt')
    screens = meta['pair_screens']
    if (len(screens) < 2 or not attempts or any(not a['passed'] or a['level'] != 0 for a in attempts)
            or any(s['refinement_level'] != 0 for s in screens)
            or not meta.get('stage_target_reached', False)):
        return dict(dt_s=proposed, reason='RUNNER_NEXT_DT_PRESERVED')
    last, prior = screens[-1], screens[-2]
    incoming = min(MAX_DT, float(prior['dt_s'])*(1.5 if prior['maximum_original_tolerance_utilization'] <= .25 else 1.))
    available = min(meta['stage_end_time_s'], meta['phase_end_time_s'])-prior['time_s']
    already_passed = any(math.isclose(s['dt_s'], incoming, rel_tol=0., abs_tol=1e-13) for s in screens[:-1])
    if (already_passed and incoming > last['dt_s']+1e-13
            and math.isclose(last['dt_s'], available, rel_tol=0., abs_tol=1e-13)
            and math.isclose(last['time_s'], meta['stage_end_time_s'], rel_tol=0., abs_tol=1e-13)
            and attempts[-1]['macro_index'] == last['macro_index']):
        return dict(dt_s=incoming, reason='RESTORED_ALREADY_PASSED_DT_AFTER_ENDPOINT_TRUNCATION',
                    runner_next_dt_s=proposed, truncated_last_dt_s=last['dt_s'])
    return dict(dt_s=proposed, reason='RUNNER_NEXT_DT_PRESERVED')


def verify_selected(folder, cost):
    meta = old.verify_parent(folder, cost)
    selected = old.inside(Path(meta['screened_continuation_checkpoint']))
    qualification = selected.parent/'candidate_qualification.json'
    old.require(old.rel(qualification) in meta['artifacts_sha256'], 'selected physical gate receipt must be hash bound')
    q = old.read(qualification)
    old.require(q['status'] == 'PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN'
                and q['accepted'] is False and q['restart_allowed'] is False, 'selected candidate physical gates failed')
    old.require(old.sha(selected) == meta['screened_continuation_sha256'], 'selected checkpoint changed')
    return meta


def verify_chain(parent, cost):
    """Audit every complete ancestor and independent origin, without a solve."""
    initial = old.inside(Path(cost['initial_checkpoint']))
    old.require(old.sha(initial) == cost['initial_sha256'], 'independent initial checkpoint changed')
    seen, lineage, first = set(), [], None
    while parent is not None:
        parent = parent.resolve()
        old.require(parent not in seen, 'cycle in parent chain')
        seen.add(parent)
        meta = verify_selected(parent, cost)
        if first is None:
            first = meta
        lineage.append(dict(stage=old.rel(parent), manifest_sha256=old.sha(parent/'manifest.json'),
                            supervision_sha256=old.sha(parent/'supervision.json'),
                            checkpoint_sha256=meta['screened_continuation_sha256']))
        prior = meta.get('parent_stage_manifest')
        if prior is None:
            old.require(meta['parent_checkpoint_sha256'] == cost['initial_sha256'], 'origin mismatch')
            parent = None
        else:
            path = old.inside(Path(prior))
            old.require(old.sha(path) == meta['parent_stage_manifest_sha256'], 'ancestor manifest changed')
            prior_meta = old.read(path)
            old.require(meta['parent_checkpoint_sha256'] == prior_meta['screened_continuation_sha256'], 'skipped/mixed ancestor')
            parent = path.parent
    return first, lineage


def validate_child(meta, child, parent, parent_meta):
    old.require(meta['status'] in old.PASS_STATUSES, 'failed/incomplete child cannot continue')
    old.require(meta['parent_stage_manifest'] == old.rel(parent/'manifest.json')
                and meta['parent_stage_manifest_sha256'] == old.sha(parent/'manifest.json')
                and meta['parent_checkpoint_sha256'] == parent_meta['screened_continuation_sha256'], 'child parent mismatch')
    old.require(meta['screened_end_time_s'] > parent_meta['screened_end_time_s']
                and meta['screened_end_step'] > parent_meta['screened_end_step'], 'no screened advancement')
    old.require(meta['screened_end_time_s'] <= endpoint(parent_meta['screened_end_time_s'])+1e-13, 'observation endpoint skipped')
    return old.snapshot(meta, child)


def claim_output(output):
    old.require(output.resolve().is_relative_to(BASE.resolve()) and output != BASE, 'new campaign output required')
    output.mkdir(parents=False, exist_ok=False)  # Atomic single-instance ownership, retained after exit.
    with (output/'instance.lock.json').open('x', encoding='utf-8') as stream:
        json.dump(dict(pid=os.getpid(), acquired_utc=utc_now(), deadline_utc=DEADLINE_UTC,
                       policy='NO_REUSE_OR_AUTOMATIC_RESUME'), stream, indent=2)


def plan():
    old.require(not OUTPUT.exists(), 'run directory already exists; duplicate or stale resume refused')
    old.require(PYTHON.is_file(), 'existing PySPH Python interpreter required')
    cost = old.verify_receipt(COST, {old.COST_PASS})
    old.require(cost['configuration'] == CONFIG, 'short configuration mismatch')
    meta, lineage = verify_chain(PARENT, cost)
    old.require(meta['screened_continuation_sha256'] == INITIAL_PARENT_SHA
                and meta['screened_end_time_s'] == .024271484375 and meta['screened_end_step'] == 158,
                'registered full initial parent mismatch')
    sources = dict(cost['source_sha256'])
    for path in (Path(__file__), BASE/'test_campaign.py', ROOT/'scripts/supervise_afsd_deposition_20261007.py', old.RUNNER):
        sources[old.rel(path)] = old.sha(path)
    inputs = {old.rel(folder/name): old.sha(folder/name) for folder in (PARENT, COST)
              for name in ('manifest.json', 'supervision.json')}
    inputs[cost['initial_checkpoint']] = cost['initial_sha256']
    inputs[meta['screened_continuation_checkpoint']] = INITIAL_PARENT_SHA
    protocol = dict(schema=SCHEMA, authorization_start_utc=START_UTC, deadline_utc=DEADLINE_UTC,
        authorization_seconds=7200., preparation_counts=True, parent=old.rel(PARENT), cost=old.rel(COST),
        configuration=CONFIG, component_receipt=cost['component_receipt'], python=str(PYTHON),
        source_sha256=sources, input_sha256=inputs, initial_lineage=lineage,
        initial_snapshot=old.snapshot(meta, PARENT), initial_dt=resume_dt(meta, old.read(PARENT/'attempt_history.json')),
        initial_budget=budget(cost['candidate_forecast_s'], meta['maximum_candidate_wall_s']),
        observation_time_s=OBSERVATION_TIME, diagnostic_upper_time_s=TARGET_TIME,
        max_stages=MAX_STAGES, max_dt_s=MAX_DT, min_dt_s=MIN_DT, growth_factor=1.5,
        forecast_limit_s=450., runner_internal_s=540., runner_external_s=600.,
        administrative_wrapper_limit_s=CHILD_WRAPPER_LIMIT, final_audit_reserve_s=FINAL_RESERVE,
        termination_reserve_s=TERMINATION_RESERVE, stop_on_first_shoulder=False,
        physical_acceptance=False, accepted_time_advanced_s=0., no_physical_module_imported=True,
        retry_policy='ONLY_ORIGINAL_RUNNER_TEMPORAL_REFINEMENT_NO_SCHEDULER_RETRY')
    return protocol, cost, meta


def command_for(protocol, child, parent, meta, cost):
    quota = budget(cost['candidate_forecast_s'], meta['maximum_candidate_wall_s'])
    dt = resume_dt(meta, old.read(parent/'attempt_history.json'))
    command = [protocol['python'], '-B', str(old.RUNNER), '--stage', 'stage', '--output', old.rel(child),
        '--config', CONFIG, '--checks', protocol['component_receipt'], '--cost', old.rel(COST), '--parent', old.rel(parent),
        '--end-time', str(endpoint(meta['screened_end_time_s'])), '--dt', str(dt['dt_s']), '--max-dt', str(MAX_DT),
        '--min-dt', str(MIN_DT), '--growth-factor', '1.5', '--max-halvings', '8',
        '--max-candidate-attempts', str(quota['candidate_quota']), '--max-macrosteps', str(quota['max_macrosteps'])]
    return command, quota, dt


def supervise(protocol, cost, meta):
    claim_output(OUTPUT)
    old.atomic_json(OUTPUT/'protocol.json', protocol)
    epoch = time.time(); mono = time.monotonic()
    monotonic_deadline = mono+max(0., DEADLINE_EPOCH-epoch)
    stages_directory = OUTPUT/'short_4p0mm'; stages_directory.mkdir()
    parent = PARENT; stages = []; process = None
    progress = dict(schema=SCHEMA, status='STARTING', supervisor_pid=os.getpid(), finished=False,
        authorization_start_utc=START_UTC, deadline_utc=DEADLINE_UTC, completed_stage_count=0,
        stages=stages, physical_acceptance=False, accepted_time_advanced_s=0., observation_60ms=None,
        **protocol['initial_snapshot'])
    def update(**values):
        progress.update(values)
        progress.update(updated_utc=utc_now(), elapsed_from_authorization_s=time.time()-START_EPOCH,
            remaining_wall_seconds=remaining_seconds(time.time(), time.monotonic(), monotonic_deadline))
        old.atomic_json(OUTPUT/'progress.json', progress)
    update()
    code = 0
    try:
        for index in range(8, 8+MAX_STAGES):
            old.verify(protocol['source_sha256']); old.verify(protocol['input_sha256'])
            verify_selected(parent, cost)
            if meta['screened_end_time_s'] >= TARGET_TIME-1e-13:
                update(status='STOPPED_GLOBAL_DIAGNOSTIC_UPPER_BOUND'); break
            remaining = remaining_seconds(time.time(), time.monotonic(), monotonic_deadline)
            if not can_launch(remaining):
                update(status='STOPPED_ABSOLUTE_DEADLINE_RESERVE'); break
            child = stages_directory/f'stage_{index:03d}'
            old.require(not child.exists(), 'unique child output occupied; no overwrite or retry')
            command, quota, dt = command_for(protocol, child, parent, meta, cost)
            if not quota['runnable']:
                update(status='STOPPED_COST_REASSESSMENT_REQUIRED', next_budget=quota); break
            update(status='RUNNING_CHILD', current_stage=old.rel(child), current_command=command,
                   next_budget=quota, dt_policy=dt, child_pid=None)
            launch_log = OUTPUT/f'stage_{index:03d}_launcher.log'
            launched = time.monotonic(); timeout = False; termination = None
            env = dict(os.environ, **{k: '1' for k in old.THREADS}, PYTHONPATH='')
            # Additional BLAS/runtime thread switches retain the one-thread intent.
            env.update(NUMEXPR_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', BLIS_NUM_THREADS='1')
            with launch_log.open('x', encoding='utf-8') as stream:
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream,
                    stderr=subprocess.STDOUT, **old.launch_options())
                update(child_pid=process.pid)
                while process.poll() is None:
                    remaining = remaining_seconds(time.time(), time.monotonic(), monotonic_deadline)
                    if time.monotonic()-launched >= CHILD_WRAPPER_LIMIT or remaining <= FINAL_RESERVE+TERMINATION_RESERVE:
                        timeout = True; termination = old.terminate_owned_tree(process, stream)
                        process.wait(timeout=20.); break
                    try:
                        process.wait(timeout=min(5., max(.05, remaining-FINAL_RESERVE-TERMINATION_RESERVE)))
                    except subprocess.TimeoutExpired:
                        live = {}
                        if (child/'manifest.json').is_file():
                            try:
                                m = old.read(child/'manifest.json')
                                live = dict(child_manifest_status=m.get('status'),
                                            child_unverified_screened_time_s=m.get('screened_end_time_s'))
                            except (OSError, ValueError):
                                pass
                        update(**live)
                exit_code = process.returncode
            row = dict(stage=old.rel(child), parent=old.rel(parent), command=command, child_pid=process.pid,
                exit_code=exit_code, timeout=timeout, termination=termination, wall_seconds=time.monotonic()-launched,
                launcher_log=old.rel(launch_log), launcher_log_sha256=old.sha(launch_log), verified=False)
            process = None
            stages.append(row); receipt = OUTPUT/f'stage_{index:03d}_launch_receipt.json'
            old.atomic_json(receipt, row)
            if timeout or exit_code != 0:
                update(status='STOPPED_CHILD_TIMEOUT' if timeout else 'STOPPED_CHILD_FAILURE', child_exit_code=exit_code)
                code = 2; break
            old.verify(protocol['source_sha256']); old.verify(protocol['input_sha256'])
            next_meta = verify_selected(child, cost)
            selected = validate_child(next_meta, child, parent, meta)
            row.update(verified=True, status=next_meta['status'], manifest_sha256=old.sha(child/'manifest.json'),
                       supervision_sha256=old.sha(child/'supervision.json'), **selected)
            old.atomic_json(receipt, row)
            parent, meta = child, next_meta
            update(status='VERIFIED_PASS_READY_FOR_NEXT_STAGE', completed_stage_count=progress['completed_stage_count']+1,
                   child_pid=None, **selected)
            if math.isclose(meta['screened_end_time_s'], OBSERVATION_TIME, rel_tol=0., abs_tol=1e-13):
                observation = dict(recorded_utc=utc_now(), **selected)
                old.atomic_json(OUTPUT/'observation_60ms.json', observation)
                update(observation_60ms=observation)
        else:
            update(status='STOPPED_STAGE_COUNT_BOUND')
    except Exception as error:
        code = 2
        if process is not None and process.poll() is None:
            with (OUTPUT/'emergency_termination.log').open('a', encoding='utf-8') as stream:
                old.terminate_owned_tree(process, stream)
                process.wait(timeout=20.)
        update(status='STOPPED_SUPERVISION_VERIFICATION_FAILURE', error=repr(error), traceback=traceback.format_exc())
    try:
        old.verify(protocol['source_sha256']); old.verify(protocol['input_sha256'])
        final_meta, lineage = verify_chain(parent, cost)
        old.require(final_meta['screened_continuation_sha256'] == progress['last_screened_checkpoint_sha256'], 'final parent mismatch')
        old.atomic_json(OUTPUT/'final_audit.json', dict(status='PASS_READ_ONLY_FULL_PARENT_CHAIN_AUDIT',
            finished_utc=utc_now(), lineage=lineage, source_sha256=protocol['source_sha256']))
        update(final_audit='PASS_READ_ONLY_FULL_PARENT_CHAIN_AUDIT')
    except Exception as error:
        code = 2
        old.atomic_json(OUTPUT/'final_audit.json', dict(status='FAIL_READ_ONLY_FULL_PARENT_CHAIN_AUDIT',
            finished_utc=utc_now(), error=repr(error), traceback=traceback.format_exc()))
        update(status='STOPPED_FINAL_AUDIT_FAILURE', final_audit='FAIL', final_audit_error=repr(error))
    update(finished=True, supervisor_exit_code=code, child_pid=None)
    print(json.dumps({k: progress[k] for k in ('status', 'completed_stage_count', 'physical_time_s', 'deadline_utc', 'supervisor_exit_code')}), flush=True)
    return code


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--preflight-only', action='store_true')
    choice.add_argument('--run', action='store_true')
    args = parser.parse_args(argv)
    protocol, cost, meta = plan()
    if not args.run:
        print(json.dumps(protocol, indent=2)); return 0
    return supervise(protocol, cost, meta)


if __name__ == '__main__':
    raise SystemExit(main())

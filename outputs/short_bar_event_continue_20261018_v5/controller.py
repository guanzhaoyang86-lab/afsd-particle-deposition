"""Resume a verified local interval with an explicit event-path proof extension."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import importlib.util
import json
import math
import os
import subprocess
import sys
import threading
import time
import traceback

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import numpy as np
import run_cpu_visibility_event_20261018 as physical
import supervise_afsd_deposition_20261007 as old
import verify_cpu_event_20261018 as verification
import verify_cpu_event_20261015 as legacy_verification
import afsd_cpu_branch_cost_20261015 as costs
import afsd_cpu_multi_event_retry_20261015 as retry_policy
import afsd_cpu_local_interval_scheduler_20261015 as intervals

LATEST_RUN = ROOT/'outputs/short_bar_event_continue_20261015_v4/run'
PARENT = LATEST_RUN/'piece_000001_attempt_000001'
PARENT_SHA = 'e93edb7ca310b551f7a76c2f3636ce5ef8b3172bc64cf7506622001bbf094793'
PARENT_MANIFEST_SHA = '340e3303bac2eeb4ad6358794c8cd9ba98a1bad467d5b2abf42e451d455829a2'
ANCHOR_SHA = '03fefe0f56537c3793be81a1d11b500375b16eec78d5482175295f7b248c2c60'
RESUME_PINS = {
    'source_bridge.json': 'c75ab8d26312fd697150b1da67284241799a6ba99627089407ebec8a4d534772',
    'recent_parent_readback.json': 'dd43dedcd48464df9d12ca3fd1ceab2886b08c4b9ab9cfa16a16ef2472f78dc6',
    'initial_local_intervals.json': '87d81252afbe7effbc51b84964a0092100d45bb045cc97fd4c5049f6b9120536',
    'scheduler.json': '6aaa326cb6e34de6edc9baf252995e31492d32f89bed98493085720d68960194',
}
CONFIGURATION = None
CONFIGURATION_SHA = None
DEADLINE = None
MONOTONIC_DEADLINE = None
ROLES = ('coarse', 'fine_1', 'fine_2')
RUN = BASE/'run'
PYTHON = ROOT/'.envs/pysph/python.exe'
require, read, sha, rel = old.require, old.read, old.sha, old.rel


def now():
    return datetime.now(timezone.utc).isoformat()


def configure(path, expected_sha):
    """Read one immutable launch deadline; workers never receive a fresh budget."""
    global CONFIGURATION, CONFIGURATION_SHA, DEADLINE, MONOTONIC_DEADLINE
    target = Path(path).resolve()
    require(target == BASE/'run.configuration.json', 'registered launch configuration required')
    require(isinstance(expected_sha, str) and sha(target) == expected_sha, 'launch configuration changed')
    config = read(target)
    require(config['schema'] == 'CPU_LOCAL_INTERVAL_RUN_CONFIGURATION_V1', 'unknown launch configuration')
    duration = config['wall_seconds']
    require(type(duration) in (int, float) and math.isfinite(duration) and duration >= 650.,
            'finite launch duration of at least 650 seconds required')
    start = datetime.fromisoformat(config['started_utc'])
    end = datetime.fromisoformat(config['deadline_utc'])
    require(start.utcoffset() is not None and end.utcoffset() is not None, 'timezone-aware deadline required')
    deadline = end.timestamp()
    require(math.isfinite(deadline) and abs((end-start).total_seconds()-duration) <= 1e-6,
            'launch start/duration/deadline disagree')
    CONFIGURATION, CONFIGURATION_SHA = target, expected_sha
    DEADLINE = deadline
    MONOTONIC_DEADLINE = time.monotonic()+max(0., DEADLINE-time.time())


def remaining():
    require(DEADLINE is not None and MONOTONIC_DEADLINE is not None, 'launch configuration required')
    return min(DEADLINE-time.time(), MONOTONIC_DEADLINE-time.monotonic())


def write(path, obj):
    old.atomic_json(path, obj)


def artifacts(folder):
    return {rel(p): sha(p) for p in sorted(folder.rglob('*')) if p.is_file()
            and p != folder/'manifest.json' and p != folder/'supervision.json'}


def source_closure(legacy):
    result = dict(legacy)
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if not filename:
            continue
        path = Path(filename).resolve()
        if path.is_relative_to(ROOT) and path.suffix == '.py' and '.envs' not in path.parts:
            key = rel(path)
            require(key not in result or result[key] == sha(path), 'changed legacy source '+key)
            result[key] = sha(path)
    for name in ('scripts/afsd_gpu_event_trial_20261010.py',
                 'scripts/run_gpu_visibility_event_delay_20261010.py',
                 'tests/review/test_cpu_event_adapter_20261015.py',
                 'docs/CPU_SHORT_BAR_EVENT_REPAIR_20261015.md',
                 'tests/review/test_cpu_branch_cost_20261015.py',
                 'tests/review/test_cpu_event_continue_admin_20261015.py',
                 'docs/CPU_SHORT_BAR_BRANCH_COST_20261015.md',
                 'tests/review/test_cpu_multi_event_retry_20261015.py',
                 'tests/review/test_cpu_multi_event_admin_20261015.py',
                 'docs/CPU_SHORT_BAR_MULTI_EVENT_SPLIT_20261015.md',
                 'tests/review/test_cpu_local_interval_scheduler_20261015.py',
                 'tests/review/test_cpu_local_interval_admin_20261015.py',
                 'tests/review/test_cpu_local_interval_launch_20261015.py',
                 'docs/CPU_SHORT_BAR_LOCAL_INTERVAL_20261015.md',
                 'outputs/short_bar_event_continue_20261015_v4/launch.py',
                 'scripts/run_cpu_visibility_event_20261018.py',
                 'scripts/verify_cpu_event_20261018.py',
                 'docs/CPU_SHORT_BAR_PATH_RECOVERY_20261018.md',
                 'docs/VISIBLE_SUPPORT_ANCHORED_CONTRACT_20261018.md',
                 'tests/review/test_cpu_path_resume_admin_20261018.py',
                 'tests/review/test_cpu_path_provider_20261018.py',
                 'tests/review/test_visible_support_anchored_20261018.py',
                 'tests/review/test_cpu_path_resume_launch_20261018.py',
                 'outputs/short_bar_event_continue_20261018_v5/launch.py',
                 'outputs/short_bar_event_continue_20261018_v5/replay_geometry.py'):
        path = ROOT/name
        result[name] = sha(path)
    result[rel(Path(__file__))] = sha(Path(__file__))
    require(CONFIGURATION is not None and sha(CONFIGURATION) == CONFIGURATION_SHA,
            'immutable launch configuration required')
    result[rel(CONFIGURATION)] = CONFIGURATION_SHA
    return result


def worker(request_path):
    p = read(request_path)
    folder = ROOT/p['output']
    began = time.perf_counter()
    deadline = began+min(525., remaining()-15.)
    require(deadline > began+10., 'absolute deadline exhausted')
    report = dict(p, status='FAIL_CPU_EVENT_BRANCH', accepted=False, restart_allowed=False,
                  diagnostic_only=True, screened_continuation_authorized=False)
    budget = dict(solves=0, maximum_solves=p['maximum_solves'], maximum_events_per_branch=1)
    timer = threading.Timer(min(540., remaining()), lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    code = 2
    try:
        require(all(os.environ.get(k) == '1' for k in old.THREADS), 'single thread required')
        old.verify(p['source_sha256'])
        require(sha(ROOT/p['role_parent_checkpoint']) == p['role_parent_sha256'], 'role parent changed')
        pa, state = physical.load(ROOT/p['role_parent_checkpoint'])
        program = physical.model.validate_state(state)
        phase = program.phase(state['time_s'])
        require(p['event_tolerance_scale'] in (1., .5, .25), 'registered tolerance scale required')
        require(state['time_s']+p['role_dt_s'] <= phase['end']+1e-15
                and phase['target'][2] <= 0 and phase['initial'][2] <= 0, 'stationary phase only')
        item = physical.complete_branch(pa, state, p['role_dt_s'], folder, deadline,
            p['candidate_forecast_s'], p['source_sha256'], budget,
            ROOT/p['role_parent_checkpoint'], p['role_parent_sha256'], p['event_tolerance_scale'])
        require(sha(ROOT/p['role_parent_checkpoint']) == p['role_parent_sha256'], 'parent modified')
        report.update(status='PASS_CPU_EVENT_BRANCH_DIAGNOSTIC_ONLY', checkpoint=rel(item['checkpoint']),
            sha256=sha(item['checkpoint']), state_time_s=item['state']['time_s'], state_step=item['state']['step'],
            leaves=item['leaves'], event_records=item['event_records'], observer_path=rel(item['observer_path']),
            physical_branch_wall_s=item['wall_seconds'])
        code = 0
    except BaseException as error:
        report.update(error=repr(error), traceback=traceback.format_exc())
    finally:
        folder.mkdir(parents=True, exist_ok=True)
        try:
            old.verify(p['source_sha256'])
        except Exception as error:
            code = 2
            report.update(status='FAIL_CPU_EVENT_BRANCH', source_error=repr(error))
        timer.cancel()
        write(folder/'protocol.json', p)
        report.update(solves=budget['solves'], wall_seconds=time.perf_counter()-began,
                      artifacts_sha256=artifacts(folder))
        write(folder/'manifest.json', report)
    return code


def progress(**values):
    path = RUN/'progress.json'
    data = read(path) if path.exists() else {}
    data.update(values, updated_utc=now(), supervisor_pid=os.getpid(),
                deadline_utc=datetime.fromtimestamp(DEADLINE, timezone.utc).isoformat())
    write(path, data)


def run_branch(pair, protocol, role, parent_checkpoint, parent_sha, role_dt, remaining_solves, cost):
    require(not (RUN/'PAUSE_REQUEST.json').exists(), 'USER_PAUSE_REQUESTED')
    require(remaining() >= 650., 'ABSOLUTE_DEADLINE_RESERVE')
    maximum = min(12, remaining_solves)
    estimate = costs.forecast(cost)
    forecast = estimate['forecast_s']
    require(maximum > 0 and estimate['eligible'], 'BRANCH_COST_REASSESSMENT_REQUIRED')
    folder = pair/role
    request = dict(protocol, role=role, role_parent_checkpoint=parent_checkpoint,
        role_parent_sha256=parent_sha, role_dt_s=role_dt, output=rel(folder),
        candidate_forecast_s=cost['maximum_single_solve_s'], maximum_solves=maximum, maximum_events_per_branch=1,
        cost_model=cost, branch_cost_estimate=estimate,
        forecast_s=forecast, internal_limit_s=540., external_limit_s=600.)
    request_path = pair/(role+'.request.json')
    write(request_path, request)
    old.verify(protocol['source_sha256'])
    command = [str(PYTHON), '-B', str(Path(__file__)), '--worker', str(request_path),
        '--configuration', str(CONFIGURATION), '--expected-configuration-sha', CONFIGURATION_SHA]
    began = time.perf_counter()
    timeout = False
    env = dict(os.environ, **{k:'1' for k in old.THREADS}, PYTHONPATH='')
    with (pair/(role+'.worker.log')).open('x', encoding='utf-8') as log:
        proc = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                                **old.launch_options())
        progress(status='RUNNING_BRANCH', current_pair=rel(pair), current_role=role, child_pid=proc.pid)
        try:
            code = proc.wait(timeout=min(600., remaining()-35.))
        except subprocess.TimeoutExpired:
            timeout = True
            old.terminate_owned_tree(proc, log)
            proc.wait(timeout=20.)
            code = proc.returncode
    folder.mkdir(exist_ok=True)
    source_ok = True
    try:
        old.verify(protocol['source_sha256'])
    except Exception:
        source_ok = False
    receipt = dict(exit_code=code, timed_out=timeout, source_inputs_still_match=source_ok,
        manifest_sha256=sha(folder/'manifest.json') if (folder/'manifest.json').exists() else None,
        protocol_sha256=sha(folder/'protocol.json') if (folder/'protocol.json').exists() else sha(request_path), wall_seconds=time.perf_counter()-began,
        command=command, request_path=rel(request_path), internal_limit_s=540., external_limit_s=600.)
    write(folder/'supervision.json', receipt)
    require(code == 0 and not timeout and source_ok, 'BRANCH_FAILED: '+rel(folder))
    result = read(folder/'manifest.json')
    require(result['status'] == 'PASS_CPU_EVENT_BRANCH_DIAGNOSTIC_ONLY', 'unqualified branch')
    old.verify(result['artifacts_sha256'])
    require(result['sha256'] == sha(ROOT/result['checkpoint']), 'branch checkpoint changed')
    return result


def run_pair(name, parent_checkpoint, parent_sha, dt, scale, sources, cost, purpose):
    pair = RUN/name
    pair.mkdir()
    pa, state = physical.load(ROOT/parent_checkpoint)
    protocol = dict(schema='BOUNDED_CPU_SHORT_BAR_EVENT_PAIR_20261015_V1',
        parent_checkpoint=parent_checkpoint, parent_sha256=parent_sha,
        parent_implementation_source_sha256=state['implementation_source_sha256'],
        parent_time_s=state['time_s'], parent_fingerprint=physical.cpu.parent_fingerprint(pa, state),
        source_sha256=sources, dt_s=dt, event_tolerance_scale=scale, purpose=purpose,
        event_protocol=physical.PROTOCOL, event_locator_protocol=physical.LOCATOR_PROTOCOL,
        blocked_prefix_protocol=physical.BLOCKED_PREFIX_PROTOCOL, visible_path_protocol=physical.VISIBLE_PATH_PROTOCOL,
        accepted=False, restart_allowed=False, diagnostic_only=True)
    write(pair/'protocol.json', protocol)
    rows = []
    began = time.perf_counter()
    for role in ROLES:
        parent = rows[1] if role == 'fine_2' else None
        row = run_branch(pair, protocol, role,
            parent['checkpoint'] if parent else parent_checkpoint,
            parent['sha256'] if parent else parent_sha,
            dt if role == 'coarse' else dt/2, 25-sum(r['solves'] for r in rows), cost)
        rows.append(row)
        supervised = read(pair/role/'supervision.json')
        observed = [read(p)['wall_seconds'] for p in sorted((pair/role).rglob('physical_solve_cost.json'))]
        cost = costs.update(cost, supervised['wall_seconds'], row['solves'], observed)
        write(RUN/'cost_model.json', cost)
    states = [physical.load(ROOT/r['checkpoint']) for r in rows]
    global_screen = physical.post.compare(states[2][0], states[0][0], state['parameters'])
    local_screen = physical.event.local_compare(states[2][0], states[0][0], sorted(int(v) for v in pa.particle_id))
    matching = states[2][1]['bore_observation']['crossing_ids'] == states[0][1]['bore_observation']['crossing_ids']
    screen = dict(global_screen=global_screen, local_screen=local_screen, crossing_ids_match=matching,
        coarse_checkpoint_sha256=rows[0]['sha256'], fine_checkpoint_sha256=[r['sha256'] for r in rows[1:]],
        passed=bool(global_screen['pass'] and local_screen['passed'] and matching), accepted=False)
    write(pair/'pair_screen.json', screen)
    meta = dict(protocol, status=physical.CONTINUE_PASS if screen['passed'] else 'FAIL_CPU_PAIR_SCREEN',
        pair_screen=screen, screened_continuation_authorized=screen['passed'],
        screened_continuation_checkpoint=rows[2]['checkpoint'], screened_continuation_sha256=rows[2]['sha256'],
        screened_end_time_s=states[2][1]['time_s'], screened_end_step=states[2][1]['step'],
        screened_history=[dict(checkpoint=r['checkpoint'], sha256=r['sha256'], time_s=r['state_time_s'],
            observer_path=r['observer_path']) for r in rows[1:]],
        event_substeps=dict(schema='EXPLICIT_COUPLED_EVENT_SUBSTEPS_V1', **{k:r['leaves'] for k,r in zip(ROLES,rows)}),
        event_records={k:r['event_records'] for k,r in zip(ROLES,rows)},
        physical_solve_count=sum(r['solves'] for r in rows), parent_prefix_unchanged=True,
        maximum_candidate_wall_s=max(r['physical_branch_wall_s'] for r in rows),
        candidate_wall_s=[r['physical_branch_wall_s'] for r in rows],
        wall_seconds=time.perf_counter()-began, observer=read(ROOT/rows[2]['observer_path']),
        full_program_completed=False)
    meta['artifacts_sha256'] = artifacts(pair)
    old.verify(sources)
    write(pair/'manifest.json', meta)
    write(pair/'supervision.json', dict(exit_code=0 if screen['passed'] else 2, timed_out=False,
        source_inputs_still_match=True, manifest_sha256=sha(pair/'manifest.json'),
        protocol_sha256=sha(pair/'protocol.json'), wall_seconds=meta['wall_seconds'],
        scope='AGGREGATE_OF_THREE_INDEPENDENTLY_BOUNDED_CHILDREN', children=[rel(pair/r/'supervision.json') for r in ROLES]))
    require(screen['passed'], 'PAIR_SCREEN_FAILED: '+name)
    verification.verify_pair(pair, expected_source=sources)
    return meta, cost


def pinned_resume(name):
    target = LATEST_RUN/name
    require(name in RESUME_PINS and sha(target) == RESUME_PINS[name],
            'registered resume evidence changed: '+name)
    return read(target)


def restore_parent():
    """Read the complete immediate parent; reuse immutable historical cost evidence."""
    require(sha(PARENT/'manifest.json') == PARENT_MANIFEST_SHA, 'registered parent manifest changed')
    result = legacy_verification.verify_pair(PARENT)
    meta = result['meta']
    require(meta['screened_continuation_sha256'] == PARENT_SHA and meta['parent_sha256'] == ANCHOR_SHA,
            'registered complete continuation parent mismatch')
    bridge = pinned_resume('source_bridge.json')
    require(bridge['to_source_sha256'] == meta['source_sha256']
        and bridge['parent_sha256'] == meta['parent_sha256']
        and bridge['parent_checkpoint'] == meta['parent_checkpoint'], 'saved source bridge mismatch')
    saved = pinned_resume('recent_parent_readback.json')
    require(bridge['recent_parent_readback_sha256'] == RESUME_PINS['recent_parent_readback.json']
        and saved['status'] == 'PASS_COMPLETE_RECENT_PARENT'
        and saved['parent_sha256'] == ANCHOR_SHA, 'saved cost evidence anchor mismatch')
    model = saved['model']
    costs.forecast(model)
    records = []
    for role in ROLES:
        branch = PARENT/role
        record, supervised = read(branch/'manifest.json'), read(branch/'supervision.json')
        observed = [read(p)['wall_seconds'] for p in sorted(branch.rglob('physical_solve_cost.json'))]
        model = costs.update(model, supervised['wall_seconds'], record['solves'], observed)
        records.append(dict(branch=rel(branch), manifest_sha256=sha(branch/'manifest.json'),
            supervision_sha256=sha(branch/'supervision.json'), solves=record['solves'],
            supervised_wall_seconds=supervised['wall_seconds']))
    write(RUN/'parent_readback.json', dict(status='PASS_COMPLETE_IMMEDIATE_PARENT',
        parent_manifest_sha256=PARENT_MANIFEST_SHA, parent_sha256=PARENT_SHA,
        parent_receipt_status=result['status'], preserved_source_count=len(meta['source_sha256']),
        source_bridge_sha256=RESUME_PINS['source_bridge.json'], physical_solves_started=0))
    write(RUN/'cost_calibration.json', dict(status='PASS_PINNED_HISTORY_AND_COMPLETE_PARENT_COST',
        cached_history_sha256=RESUME_PINS['recent_parent_readback.json'],
        cached_history_samples=saved['model']['sample_count'], records=records, model=model,
        estimate=costs.forecast(model), incomplete_pair_samples_excluded=True,
        empirical_estimate_not_deterministic_bound=True, physical_solves_started=0))
    write(RUN/'cost_model.json', model)
    return meta, model


def retry_evidence(pair, expected_parent, expected_sha):
    """Allow only the completely evidenced root-proposal multiple-opening case."""
    protocol = read(pair/'protocol.json')
    require(protocol['parent_checkpoint'] == expected_parent and protocol['parent_sha256'] == expected_sha,
            'failed pair does not bind last complete parent')
    require(sha(ROOT/expected_parent) == expected_sha, 'last complete parent changed')
    old.verify(protocol['source_sha256'])
    failed = []
    for role in ROLES:
        branch = pair/role
        if not (branch/'manifest.json').exists():
            continue
        record = read(branch/'manifest.json')
        if record['status'] != 'FAIL_CPU_EVENT_BRANCH':
            continue
        supervised = read(branch/'supervision.json')
        require(supervised['exit_code'] == 2 and supervised['timed_out'] is False
            and supervised['source_inputs_still_match'] is True
            and supervised['manifest_sha256'] == sha(branch/'manifest.json')
            and supervised['protocol_sha256'] == sha(branch/'protocol.json'),
            'failed branch must be supervised, complete and unchanged')
        require(record['solves'] == 1, 'only a root proposal may request multi-event splitting')
        verification.verify_map(record['artifacts_sha256'], branch)
        own = read(branch/'protocol.json')
        require(own['parent_checkpoint'] == expected_parent and own['parent_sha256'] == expected_sha
            and own['source_sha256'] == protocol['source_sha256'], 'failed branch identity mismatch')
        role_parent, role_sha = expected_parent, expected_sha
        if role == 'fine_2':
            first = pair/'fine_1'
            first_meta, first_supervised = read(first/'manifest.json'), read(first/'supervision.json')
            require(first_meta['status'] == 'PASS_CPU_EVENT_BRANCH_DIAGNOSTIC_ONLY'
                and first_supervised['exit_code'] == 0 and first_supervised['timed_out'] is False
                and first_supervised['source_inputs_still_match'] is True
                and first_supervised['manifest_sha256'] == sha(first/'manifest.json')
                and first_supervised['protocol_sha256'] == sha(first/'protocol.json')
                and first_meta['source_sha256'] == protocol['source_sha256'], 'fine2 parent branch incomplete')
            verification.verify_map(first_meta['artifacts_sha256'], first)
            role_parent, role_sha = rel(first/'NOT_RESTARTABLE_candidate.npz'), first_meta['sha256']
            require(first_meta['checkpoint'] == role_parent, 'fine1 alias mismatch')
        identity = dict(role=role, role_dt_s=protocol['dt_s'] if role == 'coarse' else protocol['dt_s']/2.,
            role_parent_checkpoint=role_parent, role_parent_sha256=role_sha,
            parent_checkpoint=expected_parent, parent_sha256=expected_sha,
            dt_s=protocol['dt_s'], event_tolerance_scale=protocol['event_tolerance_scale'],
            source_sha256=protocol['source_sha256'])
        require(all(own.get(k) == value and record.get(k) == value for k,value in identity.items()),
                'failed role/parent/dt/source differs from its supervised request')
        require(sha(verification.path(role_parent)) == role_sha, 'actual failed role parent changed')
        paths = {name:branch/'proposal'/name for name in ('visibility_event_routing.json',
            'union_visibility_events.json', 'visible_path_audit.json', 'continuous_shell_audit.json',
            'geometry_audit.json', 'physical_solve_cost.json')}
        require(all(record['artifacts_sha256'].get(rel(p)) == sha(p) for p in paths.values()),
                'failed root proposal evidence not bound')
        journal = read(paths['physical_solve_cost.json'])
        target = branch/'proposal/NOT_RESTARTABLE_candidate.npz'
        require(journal['status'] == 'MEASURED_COMPLETED_PHYSICAL_SOLVE'
            and journal['accepted'] is False and journal['incoming_pid_pairs'] == []
            and journal['dt_s'] == identity['role_dt_s'] and journal['checkpoint'] == rel(target)
            and journal['checkpoint_sha256'] == sha(target)
            and record['artifacts_sha256'].get(rel(target)) == sha(target)
            and verification.positive(journal['wall_seconds'])
            and journal['scope'] == 'NONLINEAR_SOLVE_FRESH_ADVANCE_AND_SAVED_STATE_EXCLUDES_LATER_GEOMETRY_QUALIFICATION',
            'unbound or invalid failed-proposal cost journal')
        parent_pa, parent_state = physical.load(verification.path(role_parent))
        require(journal['parent_fingerprint'] == physical.cpu.parent_fingerprint(parent_pa, parent_state),
                'failed proposal cost uses another parent')
        ok = retry_policy.eligible(record, read(paths['visibility_event_routing.json']),
            read(paths['visible_path_audit.json']), union_route=read(paths['union_visibility_events.json']),
            shell=read(paths['continuous_shell_audit.json']), geometry=read(paths['geometry_audit.json']))
        require(ok, 'failure outside registered multiple-opening split scope')
        failed.append(dict(role=role, branch=rel(branch), manifest_sha256=sha(branch/'manifest.json'),
            supervision_sha256=sha(branch/'supervision.json'), evidence_sha256={rel(p):sha(p) for p in paths.values()},
            measured_single_solve_s=journal['wall_seconds']))
    require(len(failed) == 1, 'exactly one eligible failed branch required')
    return dict(status='ELIGIBLE_FOR_SAME_PARENT_REAL_RECOMPUTATION_ONLY', failed_pair=rel(pair),
        parent_checkpoint=expected_parent, parent_sha256=expected_sha, branches=failed,
        accepted=False, failed_candidate_used_as_parent=False)


def seed_pending(meta, model):
    """Restore the exact remaining interval; never inherit the last tiny dt as nominal."""
    initial = pinned_resume('initial_local_intervals.json')
    saved = pinned_resume('scheduler.json')
    bridge = pinned_resume('source_bridge.json')
    require(bridge['initial_local_intervals_sha256'] == RESUME_PINS['initial_local_intervals.json']
        and initial['parent_sha256'] == meta['parent_sha256'] == ANCHOR_SHA,
        'saved pending state anchor mismatch')
    require(initial['scheduler']['accepted_time_s'] == meta['parent_time_s']
        and initial['nominal_dt_s'] == initial['scheduler']['nominal_dt'],
        'initial pending interval identity mismatch')
    require(math.isclose(intervals.current_dt(initial['scheduler'], meta['parent_time_s']),
        meta['dt_s'], rel_tol=0., abs_tol=8.*math.ulp(meta['parent_time_s'])),
        'complete parent did not cover its first pending interval')
    pending = intervals.accept(initial['scheduler'], meta['screened_end_time_s'])
    require(dict(pending, pending_targets=list(pending['pending_targets'])) == saved,
            'saved pending remainder differs from independently reconstructed state')
    require(pending['target_time_s'] == .22, 'registered stationary phase target required')
    intervals.current_dt(pending, meta['screened_end_time_s'])
    write(RUN/'initial_local_intervals.json', dict(status='RESTORED_FROM_COMPLETE_PARENT_ONLY',
        parent_sha256=PARENT_SHA, nominal_dt_s=pending['nominal_dt'], scheduler=pending,
        previous_scheduler_sha256=RESUME_PINS['scheduler.json'],
        previous_initial_scheduler_sha256=RESUME_PINS['initial_local_intervals.json'],
        physical_solves_started=0, failed_candidate_used_as_parent=False))
    return pending, model


def failed_solve_floor(model, evidence):
    # A completed solve within a failed branch can raise the single-solve floor;
    # its incomplete whole-branch wall is never a projected calibration sample.
    result = dict(model)
    require(bool(evidence['branches']) and all(type(row['measured_single_solve_s']) in (int, float)
        and math.isfinite(row['measured_single_solve_s']) and row['measured_single_solve_s'] > 0
        for row in evidence['branches']), 'failed solve costs must be finite and positive')
    result['maximum_single_solve_s'] = max(result['maximum_single_solve_s'],
        *(row['measured_single_solve_s'] for row in evidence['branches']))
    result['maximum_projected_branch_s'] = max(result['maximum_projected_branch_s'],
                                             result['maximum_single_solve_s'])
    costs.forecast(result)
    return result


def stopped_status(error):
    message = str(error)
    if 'USER_PAUSE_REQUESTED' in message:
        return 'STOPPED_USER_PAUSE'
    if 'ABSOLUTE_DEADLINE_RESERVE' in message:
        return 'STOPPED_WALL_RESERVE'
    if 'BRANCH_COST_REASSESSMENT_REQUIRED' in message:
        return 'STOPPED_COST_REASSESSMENT'
    if 'PAIR_SCREEN_FAILED' in message:
        return 'STOPPED_PAIR_SCREEN_FAILURE'
    if ('NUMERICAL_TIME_RESOLUTION_EXHAUSTED' in message
            or 'FLOATING_POINT_FINE_INTERVAL_NOT_REPRESENTABLE' in message):
        return 'STOPPED_NUMERICAL_TIME_RESOLUTION'
    return 'STOPPED_VALIDATION_OR_CHILD_FAILURE'


def run():
    require(CONFIGURATION is not None and remaining() >= 650., 'ABSOLUTE_DEADLINE_RESERVE')
    require(not RUN.exists(), 'fresh run directory required')
    RUN.mkdir()
    write(RUN/'instance.lock.json', dict(pid=os.getpid(), started_utc=now()))
    progress(status='INITIAL_VERIFICATION', finished=False, last_complete_checkpoint=None)
    sources, meta = None, None
    code = 2
    try:
        checks = read(BASE/'component_checks.json')
        require(checks['exit_code'] == 0, 'component checks required')
        old.verify(checks['source_sha256'])
        meta, cost = restore_parent()
        pending, cost = seed_pending(meta, cost)
        sources = source_closure(meta['source_sha256'])
        old.verify(sources)
        write(RUN/'source_bridge.json', dict(status='EXPLICIT_AFTER_EVENT_PATH_PROOF_EXTENSION',
            parent_manifest=rel(PARENT/'manifest.json'), parent_manifest_sha256=sha(PARENT/'manifest.json'),
            parent_checkpoint=meta['screened_continuation_checkpoint'], parent_sha256=PARENT_SHA,
            from_source_sha256=meta['source_sha256'], to_source_sha256=sources,
            component_checks=rel(BASE/'component_checks.json'), component_sha256=sha(BASE/'component_checks.json'),
            cost_calibration_sha256=sha(RUN/'cost_calibration.json'),
            parent_readback_sha256=sha(RUN/'parent_readback.json'),
            initial_local_intervals_sha256=sha(RUN/'initial_local_intervals.json'),
            configuration_sha256=CONFIGURATION_SHA,
            adaptation_audit=physical.adaptation_audit()))
        nominal_dt = pending['nominal_dt']
        write(RUN/'scheduler.json', pending)
        progress(status='READY_TO_CONTINUE', last_complete_checkpoint=meta['screened_continuation_checkpoint'],
            physical_time_s=meta['screened_end_time_s'], diagnostic_steps=meta['screened_end_step'],
            validated_new_pairs=0, observer=meta['observer']['supplied_material'],
            nominal_dt_s=nominal_dt, proposed_piece_dt_s=intervals.current_dt(pending, meta['screened_end_time_s']),
            pending_boundaries=list(pending['pending_targets']), total_refinements=3)
        index = 0
        macro_index = 1
        total_refinements = 3  # Diagnostic counter only; never a stopping quota.
        while True:
            if (RUN/'PAUSE_REQUEST.json').exists():
                status = 'STOPPED_USER_PAUSE'; break
            if remaining() < 650.:
                status = 'STOPPED_WALL_RESERVE'; break
            if meta['screened_end_time_s'] >= .22-1e-13:
                status = 'STOPPED_STATIONARY_TARGET_REACHED'; break
            if not pending['pending_targets']:
                pending = intervals.start(meta['screened_end_time_s'], nominal_dt, target=.22)
                macro_index += 1
                write(RUN/'scheduler.json', pending)
            index += 1
            attempt = 0
            while True:
                if (RUN/'PAUSE_REQUEST.json').exists():
                    raise ValueError('USER_PAUSE_REQUESTED')
                require(remaining() >= 650., 'ABSOLUTE_DEADLINE_RESERVE')
                dt = intervals.current_dt(pending, meta['screened_end_time_s'])
                attempt += 1
                name = 'piece_%06d_attempt_%06d'%(index, attempt)
                try:
                    next_meta, cost = run_pair(name, meta['screened_continuation_checkpoint'],
                        meta['screened_continuation_sha256'], dt,
                        .25, sources, cost, 'SAME_PHYSICS_LOCALLY_SUBDIVIDED_INTERVAL')
                    next_pending = intervals.accept(pending, next_meta['screened_end_time_s'])
                    meta = next_meta
                    pending = next_pending
                    write(RUN/'scheduler.json', pending)
                    break
                except Exception as error:
                    if 'BRANCH_FAILED:' not in str(error):
                        raise
                    evidence = retry_evidence(RUN/name, meta['screened_continuation_checkpoint'],
                                              meta['screened_continuation_sha256'])
                    updated = intervals.split(pending, meta['screened_end_time_s'])
                    total_refinements += 1
                    evidence.update(old_dt_s=dt,
                        new_dt_s=intervals.current_dt(updated, meta['screened_end_time_s']),
                        total_refinements=total_refinements, pending_targets=list(updated['pending_targets']))
                    pending = updated
                    cost = failed_solve_floor(read(RUN/'cost_model.json'), evidence)
                    write(RUN/(name+'.retry.json'), evidence)
                    write(RUN/'cost_model.json', cost)
                    write(RUN/'scheduler.json', pending)
                    progress(status='RETRYING_SAME_COMPLETE_PARENT_LOCAL_INTERVAL',
                        proposed_piece_dt_s=evidence['new_dt_s'], total_refinements=total_refinements,
                        pending_boundaries=list(pending['pending_targets']), rejected_pair=rel(RUN/name), child_pid=None)
            progress(status='CONTINUATION_VERIFIED', last_complete_checkpoint=meta['screened_continuation_checkpoint'],
                physical_time_s=meta['screened_end_time_s'], diagnostic_steps=meta['screened_end_step'],
                validated_new_pairs=index, observer=meta['observer']['supplied_material'],
                last_complete_pair=rel(RUN/name), last_completed_piece_dt_s=dt,
                nominal_dt_s=nominal_dt, macro_interval=macro_index,
                pending_boundaries=list(pending['pending_targets']), total_refinements=total_refinements)
        progress(status=status, finished=True, child_pid=None)
        code = 0
    except BaseException as error:
        status = stopped_status(error)
        code = 0 if status in ('STOPPED_USER_PAUSE', 'STOPPED_WALL_RESERVE') else 2
        progress(status=status, error=repr(error), traceback=traceback.format_exc(), finished=True, child_pid=None)
    finally:
        try:
            if sources:
                old.verify(sources)
            require(sha(ROOT/read(PARENT/'manifest.json')['screened_continuation_checkpoint']) == PARENT_SHA,
                    'original parent checkpoint changed')
        except BaseException as error:
            code = 2
            progress(status='STOPPED_SOURCE_OR_PARENT_INTEGRITY', error=repr(error), finished=True, child_pid=None)
        write(RUN/'result.json', dict(exit_code=code, finished_utc=now(), pid=os.getpid()))
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker', type=Path)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--configuration', type=Path)
    parser.add_argument('--expected-configuration-sha')
    args = parser.parse_args()
    if args.worker or args.run:
        require(args.configuration is not None and args.expected_configuration_sha,
                'explicit immutable launch configuration required')
        configure(args.configuration, args.expected_configuration_sha)
    raise SystemExit(worker(args.worker) if args.worker else run() if args.run else 0)

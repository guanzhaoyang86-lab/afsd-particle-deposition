"""One locked-history session: same-parent event tolerance study, then .25 continuation."""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
import ast
import copy
import inspect
import itertools
import json
import math
import multiprocessing
import os
import sys
import time
import traceback
import numpy as np
import run_gpu_event_campaign_20261010 as frozen
import run_gpu_event_sensitivity_20261010 as previous_admin

ROOT, BASE, CPU_BASE = frozen.ROOT, frozen.BASE, frozen.CPU_BASE
bridge, admin, adaptive = frozen.bridge, frozen.admin, frozen.adaptive
require, read, load, rel = frozen.require, frozen.read, frozen.load, frozen.rel
PROTOCOL = 'SAME_PARENT_BOX_EVENT_SENSITIVITY_20261010_V1'
SESSION_PROTOCOL = PROTOCOL
LOCATOR_PROTOCOL = 'REAL_BLOCKED_LENGTH_SECANT_20261010_V1'
BLOCKED_PREFIX_PROTOCOL = 'FIXED_WITNESS_BOX_CONTINUOUS_INCOMING_BLOCKED_PREFIX_20261010_V1'
FORECAST_LIMIT_S = 900.
PHYSICAL_SOFT_LIMIT_S = 960.
WORKER_HARD_LIMIT_S = 990.
CONTROLLER_LIMIT_S = 1080.
QUALIFICATION_WALL_RESERVE_S = 1100.
TRANSITION_SCHEMA = 'EXPLICIT_BOUNDED_EVENT_V1_TO_FIXED_BOX_V3_20261010'
authorize_grant = frozen.authorize_grant
COMPONENTS = BASE/'event_box_components_20261010_v1.json'
CONTRACT = ROOT/'docs/GPU_EVENT_BOX_CAMPAIGN_20261010.md'
COMMON_PARENT = BASE/'event_validation_20261010_v1/attempt_0001'
FAILED_V1 = BASE/'event_validation_20261010_v1/attempt_0002'
FAILED_V2 = BASE/'event_secant_sensitivity_20261010_v1/sensitivity_one'
PHYSICAL_FILES = ('scripts/run_gpu_visibility_event_box_20261010.py',
    'scripts/afsd_event_blocked_box_20261010.py', 'tests/review/test_event_blocked_box_20261010.py',
    'tests/review/test_gpu_visibility_event_box_20261010.py', 'docs/GPU_VISIBILITY_EVENT_BOX_20261010.md')
ADMINISTRATIVE_FILES = ('scripts/run_gpu_event_box_campaign_20261010.py',
    'scripts/afsd_event_box_worker_20261010.py',
    'scripts/check_gpu_event_box_20261010.py',
    'tests/review/test_gpu_event_box_campaign_20261010.py', 'docs/GPU_EVENT_BOX_CAMPAIGN_20261010.md')
SCALES = (1., .5, .25)
ROLES = ('coarse', 'fine_1', 'fine_2')


def checks(pool):
    previous = previous_admin.checks(pool)
    receipt = read(COMPONENTS)
    require(receipt['status'] == 'PASS_EVENT_BOX_COMPONENTS' and receipt['exit_code'] == 0,
            'executed sensitivity administrative checks required')
    require(receipt['physical_source_changed'] is True and receipt['gpu_campaign_runs'] == 0
        and set(receipt['physical_source_sha256']) == set(PHYSICAL_FILES), 'explicit registered box source bridge required')
    require(receipt['inherited_receipt'] == rel(previous_admin.COMPONENTS)
        and receipt['inherited_receipt_sha256'] == pool.sha(previous_admin.COMPONENTS),
            'frozen secant v2 receipt required')
    require(receipt.get('failed_v2_stage') == rel(FAILED_V2)
        and all(rel(FAILED_V2/name) in receipt['artifacts_sha256'] for name in ('manifest.json', 'supervision.json', 'run.log')),
            'complete sealed v2 failure evidence registration required')
    require(set(ADMINISTRATIVE_FILES).issubset(receipt['source_sha256'])
        and all(receipt['source_sha256'].get(k) == v for k, v in previous.items())
        and all(receipt['source_sha256'].get(k) == v for k, v in receipt['physical_source_sha256'].items()),
            'new administrative closure must preserve every frozen event component')
    pool.seal(receipt['source_sha256']); pool.seal(receipt['artifacts_sha256'])
    return frozen.merge(previous, receipt['source_sha256'], receipt['artifacts_sha256'],
                        {rel(COMPONENTS): pool.sha(COMPONENTS)})


def box_sources(pool, v2):
    receipt = read(COMPONENTS)
    require(set(receipt['physical_source_sha256']) == set(PHYSICAL_FILES), 'exact registered box physical files required')
    pool.seal(receipt['source_sha256']); pool.seal(receipt['artifacts_sha256'])
    return frozen.merge(v2, receipt['physical_source_sha256'], {rel(COMPONENTS): pool.sha(COMPONENTS)})


def source_transition(template, meta, parent):
    if meta['source_sha256'] == template['source_sha256']:
        require(meta.get('event_locator_protocol') == LOCATOR_PROTOCOL
            and meta.get('blocked_prefix_protocol') == BLOCKED_PREFIX_PROTOCOL,
                'box parent lacks explicit locator and prefix protocol')
        return None
    require(Path(parent).resolve() == COMMON_PARENT.resolve()
        and meta['source_sha256'] == template['v1_source_sha256']
        and meta.get('event_protocol') == frozen.PROTOCOL, 'unregistered v1-to-box source boundary')
    return dict(schema=TRANSITION_SCHEMA, parent_directory=rel(parent),
        parent_checkpoint=meta['screened_continuation_checkpoint'], parent_sha256=meta['screened_continuation_sha256'],
        from_source_sha256=meta['source_sha256'], to_source_sha256=template['source_sha256'],
        components=rel(COMPONENTS), components_sha256=template['box_components_sha256'])


def source_generation(meta, template):
    if meta['source_sha256'] == template['source_sha256']:
        require(meta.get('event_locator_protocol') == LOCATOR_PROTOCOL and meta.get('event_protocol') == frozen.PROTOCOL
            and meta.get('blocked_prefix_protocol') == BLOCKED_PREFIX_PROTOCOL,
                'fifth implementation generation lacks explicit box event identity')
        return 4
    if meta['source_sha256'] == template['v2_source_sha256']:
        require(meta.get('event_locator_protocol') == LOCATOR_PROTOCOL and meta.get('event_protocol') == frozen.PROTOCOL,
                'fourth implementation generation lacks explicit secant identity')
        return 3
    return frozen.source_generation(meta, template['legacy_source_sha256'],
        template['search_source_sha256'], template['v1_source_sha256'])


def static_audit():
    import run_gpu_visibility_event_box_20261010 as box
    expected = inspect.getsource(frozen.run_stage).replace('afsd_event_worker_20261010', 'afsd_event_box_worker_20261010').replace('started+600.', 'started+1080.')
    require(ast.dump(ast.parse(expected)) == ast.dump(ast.parse(inspect.getsource(run_stage))),
            'live supervisor changed beyond explicit box worker name and 1080-second bound')
    old_worker = (ROOT/'scripts/afsd_event_worker_20261010.py').read_text(encoding='utf-8')
    new_worker = (ROOT/'scripts/afsd_event_box_worker_20261010.py').read_text(encoding='utf-8')
    require(ast.dump(ast.parse(old_worker.replace('import run_gpu_visibility_event_20261010 as physical',
        'import run_gpu_visibility_event_box_20261010 as physical').replace('threading.Timer(540.,timeout)',
        'threading.Timer(990.,timeout)').replace('began+525.', 'began+960.'))) == ast.dump(ast.parse(new_worker)),
            'box worker changed beyond explicit physical provider and registered 960/990-second bounds')
    require(box.BLOCKED_PREFIX_PROTOCOL == BLOCKED_PREFIX_PROTOCOL, 'physical prefix protocol differs')
    return dict(status='PASS_EXPLICIT_BOX_SOURCE_BRIDGE_AND_UNCHANGED_LIVE_SUPERVISION',
        runtime_global_rebinding=False, inherited_event_audit=frozen.static_audit(),
        box_physical_audit=box.adaptation_audit(), registered_wall_limits_s=dict(forecast=FORECAST_LIMIT_S,
            physical_soft=PHYSICAL_SOFT_LIMIT_S, worker_hard=WORKER_HARD_LIMIT_S, controller=CONTROLLER_LIMIT_S),
        absolute_stop=frozen.ABSOLUTE_STOP.isoformat(), grant_timeout_unchanged_s=180.)


def sibling_identity(meta, parent_meta, parent, scale, dt, source=None):
    """A comparison stage is a sibling, never a continuation of another variant."""
    require(meta.get('status') == frozen.CONTINUE_PASS and meta.get('screened_continuation_authorized') is True,
            'complete qualified PASS stage required')
    require(meta.get('event_protocol') == frozen.PROTOCOL and meta.get('event_tolerance_scale') == scale
        and meta.get('event_locator_protocol') == LOCATOR_PROTOCOL
        and meta.get('blocked_prefix_protocol') == BLOCKED_PREFIX_PROTOCOL,
            'explicit registered event tolerance scale required')
    target = source if source is not None else meta['source_sha256']
    require(meta['source_sha256'] == target and all(meta['source_sha256'].get(k) == v for k, v in parent_meta['source_sha256'].items())
        and meta['parent_implementation_source_sha256'] == parent_meta['source_sha256'],
            'sensitivity must retain the exact frozen event implementation')
    require(meta['parent_gpu_directory'] == rel(parent)
        and meta['parent_checkpoint'] == parent_meta['screened_continuation_checkpoint']
        and meta['parent_sha256'] == parent_meta['screened_continuation_sha256'],
            'all tolerance variants must share the exact complete parent')
    require(type(dt) in (float, int) and math.isfinite(dt) and 1e-7 <= dt <= .005
        and meta['dt_s'] == dt and meta['parent_time_s'] == parent_meta['screened_end_time_s']
        and abs(meta['screened_end_time_s']-parent_meta['screened_end_time_s']-dt) <= 1e-14,
            'all tolerance variants must cover the same macro interval')
    require(meta['pair_screen']['passed'] is True, 'original coarse/fine PASS required')
    require(set(meta['event_records']) == set(ROLES)
        and any(meta['event_records'][role] for role in ROLES), 'actual qualified activation event evidence required')
    return True


def cost_evidence(pool, folder, meta):
    measured = meta['maximum_candidate_wall_s']
    require(type(measured) in (int, float) and math.isfinite(measured) and measured > 0,
            'positive measured whole-branch event cost required')
    require(len(meta['candidate_wall_s']) == 3 and all(type(v) in (int, float) and math.isfinite(v)
        and v > 0 for v in meta['candidate_wall_s']) and max(meta['candidate_wall_s']) == measured,
            'cost must include all three whole macro branches and their probes')
    supervision = read(folder/'supervision.json')
    require(supervision['exit_code'] == 0 and supervision['manifest_sha256'] == pool.sha(folder/'manifest.json'),
            'whole-pair cost requires its complete successful supervision receipt')
    worker_wall = meta['wall_seconds']; supervised_wall = supervision['wall_seconds']
    require(all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in (worker_wall, supervised_wall)),
            'finite measured worker and supervision wall times required')
    single = single_solve_evidence(pool, folder, meta)
    return dict(stage=rel(folder), manifest_sha256=pool.sha(folder/'manifest.json'),
        supervision_sha256=pool.sha(folder/'supervision.json'), maximum_candidate_wall_s=measured,
        candidate_wall_s=list(meta['candidate_wall_s']), candidate_total_wall_s=sum(meta['candidate_wall_s']),
        worker_wall_s=worker_wall, supervised_wall_s=supervised_wall,
        measured_pair_wall_s=max(sum(meta['candidate_wall_s']), worker_wall, supervised_wall),
        parent_sha256=meta['parent_sha256'], dt_s=meta['dt_s'], event_tolerance_scale=meta['event_tolerance_scale'],
        event_count=sum(len(meta['event_records'][role]) for role in ROLES), includes_event_probes=True, **single)


def single_solve_evidence(pool, folder, meta):
    names = sorted(name for name in meta['artifacts_sha256'] if name.endswith('/physical_solve_cost.json'))
    if Path(folder).resolve() == COMMON_PARENT.resolve():
        require(not names and meta['physical_solve_count'] == 3
            and not any(meta['event_records'][role] for role in ROLES), 'only the exact no-event common parent may use a startup proxy')
        return dict(maximum_single_solve_wall_s=meta['maximum_candidate_wall_s'], single_solve_cost_is_proxy=True,
                    single_solve_cost_basis='QUALIFIED_OLD_NO_EVENT_MACRO_INCLUDES_GEOMETRY_NOT_SINGLE_SOLVE_MEASUREMENT',
                    single_solve_journals=[])
    require(meta.get('blocked_prefix_protocol') == BLOCKED_PREFIX_PROTOCOL
        and type(meta['physical_solve_count']) is int and 0 < meta['physical_solve_count'] <= 25
        and len(names) == meta['physical_solve_count'], 'every completed box physical solve requires a bound cost journal')
    rows = []; seen = set()
    for name in names:
        path = (ROOT/name).resolve(); digest = pool.sha(path); row = read(path)
        require(path.is_relative_to(Path(folder).resolve()) and digest == meta['artifacts_sha256'][name],
                'physical solve cost journal path or SHA differs')
        target = (ROOT/row['checkpoint']).resolve()
        require(target == path.parent/'NOT_RESTARTABLE_candidate.npz' and target not in seen
            and meta['artifacts_sha256'].get(row['checkpoint']) == row['checkpoint_sha256']
            and pool.sha(target) == row['checkpoint_sha256'], 'physical solve cost checkpoint path, uniqueness or SHA differs')
        seen.add(target)
        require(row['status'] == 'MEASURED_COMPLETED_PHYSICAL_SOLVE' and row.get('accepted') is False
            and row['scope'] == 'NONLINEAR_SOLVE_FRESH_ADVANCE_AND_SAVED_STATE_EXCLUDES_LATER_GEOMETRY_QUALIFICATION',
                'exact completed single-solve measurement scope required')
        value, dt = row['wall_seconds'], row['dt_s']
        require(all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in (value, dt))
            and value <= meta['wall_seconds'] and dt <= meta['dt_s']*(1.+1e-12), 'invalid single-solve cost or physical interval')
        adaptive.digest(row['parent_fingerprint'])
        rows.append(dict(journal=name, sha256=digest, **row))
    return dict(maximum_single_solve_wall_s=max(row['wall_seconds'] for row in rows),
        single_solve_cost_is_proxy=False, single_solve_cost_basis='ALL_SHA_BOUND_COMPLETED_PHYSICAL_SOLVE_JOURNALS',
        single_solve_journals=rows)


def failed_qualification_basis(pool):
    meta = read(FAILED_V2/'manifest.json'); supervision = read(FAILED_V2/'supervision.json')
    require(meta['status'] == 'FAIL_GPU_STAGE_NOT_ACCEPTED' and supervision['exit_code'] != 0
        and pool.sha(FAILED_V2/'manifest.json') == supervision['manifest_sha256']
        and pool.sha(FAILED_V2/'run.log') == supervision['log_sha256'], 'sealed failed-v2 qualification cost basis required')
    return dict(stage=rel(FAILED_V2), manifest_sha256=pool.sha(FAILED_V2/'manifest.json'),
        supervision_sha256=pool.sha(FAILED_V2/'supervision.json'), partial_worker_wall_s=meta['wall_seconds'],
        partial_supervised_wall_s=supervision['wall_seconds'], complete_pair_measured=False,
        usable_parent=False, first_box_qualification_budget_s=FORECAST_LIMIT_S,
        reason='REAL_V2_EVENT_LOCALIZED_BUT_BLOCKED_PREFIX_CERTIFICATE_EXHAUSTED_QUERY_BUDGET')


def measured_envelope(rows, *, scope='EMPIRICAL_COMPLETE_EVENT_PAIR', parent_sha256=None, dt_s=None):
    require(bool(rows), 'measured event cost evidence required')
    for row in rows:
        value = row['maximum_candidate_wall_s']
        require(type(value) in (int, float) and math.isfinite(value) and value > 0
            and row.get('includes_event_probes') is True, 'invalid measured whole-branch cost evidence')
        single = row['maximum_single_solve_wall_s']
        require(type(single) in (int, float) and math.isfinite(single) and single > 0
            and type(row['single_solve_cost_is_proxy']) is bool, 'positive explicitly scoped single-solve cost required')
    require(scope in ('FIRST_BOX_QUALIFICATION', 'EMPIRICAL_COMPLETE_EVENT_PAIR'), 'registered cost scope required')
    used = rows
    if scope == 'EMPIRICAL_COMPLETE_EVENT_PAIR':
        require(any(type(r.get('event_count')) is int and r['event_count'] > 0 for r in used),
                'qualified actual secant event pair required before empirical continuation forecast')
        for row in used:
            values = (row['candidate_total_wall_s'], row['worker_wall_s'], row['supervised_wall_s'])
            require(all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in values)
                and len(row['candidate_wall_s']) == 3
                and all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in row['candidate_wall_s'])
                and sum(row['candidate_wall_s']) == row['candidate_total_wall_s']
                and max(row['candidate_wall_s']) == row['maximum_candidate_wall_s']
                and row['measured_pair_wall_s'] == max(values), 'whole event pair cost evidence is inconsistent')
        pair = max(r['measured_pair_wall_s'] for r in used)
        formula = '30 + 1.5 * maximum_measured_complete_event_pair_seconds'
        forecast = 30.+1.5*pair
    else:
        require(len(rows) == 1 and rows[0].get('event_count') == 0,
                'first unmeasured qualification may use only its original no-event common parent')
        pair = None; formula = 'FIRST_QUALIFICATION_ONLY: REGISTERED_900_SECOND_UNMEASURED_COMPLETE_EVENT_PAIR_BUDGET'
        forecast = FORECAST_LIMIT_S
    measured_single = [r['maximum_single_solve_wall_s'] for r in used if not r['single_solve_cost_is_proxy']]
    require(measured_single or scope == 'FIRST_BOX_QUALIFICATION', 'completed box event pair needs actual single-solve measurements')
    single = max(measured_single) if measured_single else max(r['maximum_single_solve_wall_s'] for r in used)
    return dict(scope=scope, maximum_candidate_wall_s=max(row['maximum_candidate_wall_s'] for row in used),
        maximum_single_solve_wall_s=single, single_solve_cost_is_proxy=not bool(measured_single),
        maximum_complete_pair_wall_s=pair, forecast_s=forecast, formula=formula, evidence=copy.deepcopy(used),
        empirical_scheduling_estimate_only=True, guarantees_unknown_events=False,
        event_model_cost_previously_measured=scope != 'FIRST_BOX_QUALIFICATION')


def cost_decision(meta, dt, remaining, phases, rows, *, scope='EMPIRICAL_COMPLETE_EVENT_PAIR', qualification=False):
    envelope = measured_envelope(rows, scope=scope, parent_sha256=meta.get('screened_continuation_sha256'), dt_s=dt)
    t = float(meta['screened_end_time_s']); minimum_dt = 1e-7; forecast = envelope['forecast_s']
    require(all(math.isfinite(v) for v in (t, dt, remaining)) and t >= 0 and minimum_dt <= dt <= .005,
            'finite bounded cost/time inputs required')
    result = dict(forecast_s=forecast, child_timeout_s=max(0., min(CONTROLLER_LIMIT_S, remaining)), event_cost_envelope=envelope)
    if (qualification or scope == 'FIRST_BOX_QUALIFICATION') and remaining < QUALIFICATION_WALL_RESERVE_S:
        return dict(result, action='STOP_QUALIFICATION_WALL_RESERVE')
    if remaining < forecast+20.: return dict(result, action='STOP_WALL_RESERVE')
    if forecast > FORECAST_LIMIT_S: return dict(result, action='STOP_COST_REASSESSMENT')
    phase = next((p for p in phases if p['start'] <= t < p['end']), None)
    if phase is None or not phase['stationary']: return dict(result, action='STOP_REGISTERED_STATIONARY_FILL_END')
    try:
        bounded = admin.bounded_dt(t, dt, phase['end'])
    except ValueError as error:
        return dict(result, action='STOP_UNREPRESENTABLE_PHASE_BOUNDARY', reason=str(error))
    if bounded < minimum_dt: return dict(result, action='STOP_MINIMUM_DT_OR_BOUNDARY')
    return dict(result, action='RUN', dt_s=bounded, phase=phase['name'], phase_end_s=phase['end'])


def stage_protocol(template, meta, parent, dt, administrative, session_record, refs, scale, costs,
                   *, scope='EMPIRICAL_COMPLETE_EVENT_PAIR'):
    require(scale in SCALES, 'registered tolerance scale required')
    if scope == 'FIRST_BOX_QUALIFICATION':
        require(scale == 1. and Path(parent).resolve() == COMMON_PARENT.resolve()
            and meta['source_sha256'] == template['v1_source_sha256'],
                'unmeasured qualification is restricted to the first exact v1-parent tolerance-one solve')
    envelope = measured_envelope(costs, scope=scope, parent_sha256=meta['screened_continuation_sha256'], dt_s=dt)
    require(envelope['forecast_s'] <= FORECAST_LIMIT_S, 'registered box event cost forecast limit exceeded')
    transition = source_transition(template, meta, parent)
    frozen.merge(template['source_sha256'], refs, administrative, session_record)
    result = copy.deepcopy(template)
    for name in ('search_tensor_reference_stage', 'administrative_preflight_statistics'):
        result.pop(name, None)
    result.update(mode='continue', parent_checkpoint=meta['screened_continuation_checkpoint'],
        parent_sha256=meta['screened_continuation_sha256'], parent_gpu_directory=rel(parent),
        parent_implementation_source_sha256=meta['source_sha256'], dt_s=dt,
        reference_sha256=frozen.merge(refs, administrative, session_record), administrative_source_sha256=administrative,
        implementation_transition=transition, event_protocol=frozen.PROTOCOL, event_locator_protocol=LOCATOR_PROTOCOL,
        blocked_prefix_protocol=BLOCKED_PREFIX_PROTOCOL,
        administrative_protocol=PROTOCOL, event_tolerance_scale=scale,
        candidate_forecast_s=envelope['maximum_single_solve_wall_s'], forecast_s=envelope['forecast_s'],
        event_cost_envelope=envelope, sensitivity_contract=rel(CONTRACT))
    result.update(forecast_limit_s=FORECAST_LIMIT_S, internal_limit_s=WORKER_HARD_LIMIT_S,
        external_limit_s=CONTROLLER_LIMIT_S, physical_soft_deadline_s=PHYSICAL_SOFT_LIMIT_S, worker_deadline_s=PHYSICAL_SOFT_LIMIT_S,
        qualification_wall_reserve_s=QUALIFICATION_WALL_RESERVE_S)
    return result


def event_timing(meta, reader=read, loader=load):
    """Read already fully receipt-validated, SHA-bound physical brackets."""
    result = {}
    for role in ROLES:
        result[role] = []
        for event in meta['event_records'][role]:
            location = reader(ROOT/event['location'])
            require(location['status'] == 'BRACKETED_ACTIVATION'
                and location['lo']['blocked'] is True and location['hi']['blocked'] is False,
                    'real blocked-to-visible bracket required')
            _, low = loader(ROOT/location['lo']['checkpoint'])
            _, high = loader(ROOT/location['hi']['checkpoint'])
            lo, hi, tolerance = low['time_s'], high['time_s'], event['tolerance_s']
            require(all(type(v) in (float, int) and math.isfinite(v) for v in (lo, hi, tolerance))
                and 0 < hi-lo <= tolerance+1e-14 and tolerance > 0., 'invalid physical event bracket')
            pairs = event['incoming_pid_pairs']
            require(pairs and all(len(p) == 2 and type(p[0]) is int and type(p[1]) is int
                and p[0] < p[1] for p in pairs) and pairs == sorted(pairs), 'canonical event PID pairs required')
            result[role].append(dict(incoming_pid_pairs=pairs, direction='BLOCKED_TO_VISIBLE',
                low_time_s=lo, high_time_s=hi, tolerance_s=tolerance,
                low_checkpoint=location['lo']['checkpoint'], high_checkpoint=location['hi']['checkpoint'],
                location=event['location'], location_sha256=event['location_sha256']))
    return result


def compare_event_timing(left, right):
    rows = []
    for role in ROLES:
        require(len(left[role]) == len(right[role]), 'event count differs across tolerance variants')
        for number, (a, b) in enumerate(zip(left[role], right[role])):
            require(a['incoming_pid_pairs'] == b['incoming_pid_pairs']
                and a['direction'] == b['direction'] == 'BLOCKED_TO_VISIBLE',
                    'event PID identity, direction or order changed across tolerance variants')
            gap = max(0., max(a['low_time_s'], b['low_time_s'])-min(a['high_time_s'], b['high_time_s']))
            bound = a['tolerance_s']+b['tolerance_s']
            require(math.isfinite(gap) and math.isfinite(bound) and bound > 0 and gap <= bound+1e-14,
                    'event timing differs beyond summed verified activation tolerances')
            rows.append(dict(role=role, event=number, bracket_gap_s=gap, allowed_gap_s=bound,
                             incoming_pid_pairs=a['incoming_pid_pairs'], passed=True))
    return rows


def compare_endpoint(left, right, expected_ids, expected_mass, parameters):
    a, sa = left; b, sb = right
    require(np.array_equal(a.particle_id, expected_ids) and np.array_equal(b.particle_id, expected_ids),
            'persistent particle IDs or order differ from the complete common parent')
    require(np.array_equal(a.m, expected_mass) and np.array_equal(b.m, expected_mass),
            'particle masses differ from the complete common parent')
    require(sa['parameters'] == parameters and sb['parameters'] == parameters
        and abs(sa['time_s']-sb['time_s']) <= 1e-14, 'endpoint physics or absolute time differ')
    ids = sorted(int(v) for v in expected_ids)
    glob = bridge.post.compare(a, b, parameters)
    local = bridge.event.local_compare(a, b, ids)
    matching = sa['bore_observation']['crossing_ids'] == sb['bore_observation']['crossing_ids']
    screen = dict(global_screen=glob, local_screen=local, crossing_ids_match=matching,
                  passed=bool(glob['pass'] and local['passed'] and matching), accepted=False)
    normalized = adaptive.normalized_screen(screen, parameters['dx'], ids)
    require(normalized['passed'], 'cross-tolerance original global/local/PID screen failed')
    return dict(screen=screen, normalized=normalized, exact_parent_masses=True)


def compare_group(parent_meta, parent, stages, loader=load, reader=read):
    require(set(stages) == set(SCALES), 'all three registered tolerance stages required')
    baseline_dt = stages[1.][1]['dt_s']
    parent_pa, parent_state = loader(ROOT/parent_meta['screened_continuation_checkpoint'])
    endpoints = {}; timing = {}; records = []
    for scale in SCALES:
        folder, meta = stages[scale]
        sibling_identity(meta, parent_meta, parent, scale, baseline_dt, stages[1.][1]['source_sha256'])
        endpoints[scale] = {role: loader(folder/role/'NOT_RESTARTABLE_candidate.npz') for role in ROLES}
        timing[scale] = event_timing(meta, reader=reader, loader=loader)
    for a, b in itertools.combinations(SCALES, 2):
        screens = {role: compare_endpoint(endpoints[b][role], endpoints[a][role], parent_pa.particle_id,
                   parent_pa.m, parent_state['parameters']) for role in ROLES}
        events = compare_event_timing(timing[a], timing[b])
        records.append(dict(scales=[a, b], macro_endpoint_screens=screens, event_timing=events))
    return dict(status='PASS_SAME_PARENT_EVENT_TOLERANCE_SENSITIVITY', schema=PROTOCOL,
        common_parent=rel(parent), common_parent_checkpoint=parent_meta['screened_continuation_checkpoint'],
        common_parent_sha256=parent_meta['screened_continuation_sha256'], dt_s=baseline_dt,
        stages=[dict(scale=scale, directory=rel(stages[scale][0]),
            checkpoint_sha256=stages[scale][1]['screened_continuation_sha256']) for scale in SCALES],
        selected_scale=.25, selected_stage=rel(stages[.25][0]), comparisons=records,
        event_timing=timing, accepted=False, physical_model_unchanged=True)


def selected_pass(policy, selected, state):
    """Three sibling experiments advance the adaptive policy exactly once."""
    policy = adaptive.begin_attempt(policy, selected['dt_s'])
    return adaptive.resolve_pass(policy, selected['pair_screen'], state['parameters']['dx'],
        state['bore_observation']['particle_ids'], selected['screened_continuation_sha256'],
        original_candidate_gates_passed=True)


def keep_failure(pool, child, source, parent_meta, refs):
    """Seal a failed experiment too; it can never become a continuation parent."""
    if not (child/'manifest.json').exists():
        return refs, {}
    failure = read(child/'manifest.json'); supervision = read(child/'supervision.json')
    require(supervision['exit_code'] != 0 and pool.sha(child/'manifest.json') == supervision['manifest_sha256']
        and pool.sha(child/'run.log') == supervision['log_sha256'], 'failed stage supervision mismatch')
    require(failure['source_sha256'] == source and failure['parent_sha256'] == parent_meta['screened_continuation_sha256']
        and failure['parent_implementation_source_sha256'] == parent_meta['source_sha256'],
            'failed stage source or parent differs')
    pool.seal(failure['source_sha256']); pool.seal(failure['reference_sha256'])
    return frozen.bind_stage(pool, child, failure, refs), failure


def run(args):
    output = (ROOT/args.output).resolve(); parent = COMMON_PARENT.resolve()
    require(output.is_relative_to(BASE) and not output.exists(), 'new project-local output required')
    require(0 < args.wall_seconds <= 36000 and 0 <= args.max_stages <= 2000,
            'bounded original wall authorization and continuation stage count required')
    require(0 < args.maximum_dt <= .005 and 0 < args.target_time <= .22
        and 0 < args.shoulder_observe_ms <= 20., 'original continuation bounds required')
    require(read(parent/'manifest.json')['status'] == frozen.CONTINUE_PASS,
            'complete verified v1 no-event replay parent must actually PASS first')
    started = datetime.now(timezone.utc); wall = frozen.authorized_wall_seconds(args.wall_seconds, started)
    output.mkdir(parents=True); began = time.perf_counter(); deadline = began+wall
    pool = frozen.LockedHistory(ROOT); kernel = mutex = None; code = 0; frames = []
    progress = dict(status='INITIAL_FULL_HISTORY_VERIFICATION', finished=False, supervisor_pid=os.getpid(),
        child_pid=None, completed_stage_count=0, validation_variants_completed=0, attempt_count=0,
        started_at=started.astimezone(admin.LOCAL).isoformat(),
        latest_stop_at=(started+timedelta(seconds=wall)).astimezone(admin.LOCAL).isoformat(),
        accepted=False, formal_accepted_time_advanced_s=0.)
    def update(**fields):
        progress.update(fields, updated_at=admin.wall_now(), wall_seconds=time.perf_counter()-began,
            remaining_wall_seconds=max(0., deadline-time.perf_counter()))
        admin.atomic_json(output/'progress.json', progress)
    update()
    try:
        kernel, mutex = admin.singleton_and_awake()
        administrative = checks(pool); audit = static_audit()
        startup = argparse.Namespace(mode='continue', parent=rel(parent), dt=None,
            baseline=rel(CPU_BASE/'stage_013'), cost=rel(BASE/'benchmark_v1'),
            validation=rel(BASE/'replay_v1'), reference_stage=None)
        template = frozen._initial_preflight(startup, pool)
        v1 = dict(template['source_sha256'])
        v2 = previous_admin.secant_sources(pool, v1)
        template.update(v1_source_sha256=v1, v2_source_sha256=v2, source_sha256=box_sources(pool, v2),
            box_components_sha256=pool.sha(COMPONENTS), event_locator_protocol=LOCATOR_PROTOCOL,
            blocked_prefix_protocol=BLOCKED_PREFIX_PROTOCOL)
        template['unmeasured_box_qualification_cost_basis'] = failed_qualification_basis(pool)
        tx = frozen.LockedPreflight(pool)
        common = tx.gpu_receipt(parent, {frozen.CONTINUE_PASS})
        tx.finish()
        require(source_generation(common, template) == 2 and common['source_sha256'] == v1 and common['incremental_replay_particle_arrays_bitwise_equal'] ==
            dict(coarse=True, fine_1=True, fine_2=True) and not any(common['event_records'][r] for r in ROLES),
                'qualified exact v1 no-event replay common parent required')
        dt = common['dt_s']
        require(dt <= args.maximum_dt and common['screened_end_time_s']+dt <= args.target_time,
                'sensitivity interval exceeds registered recovery/target bounds')
        refs = frozen.bind_stage(pool, parent, common, dict(template['reference_sha256']))
        costs = [cost_evidence(pool, parent, common)]
        phases = admin.phase_table(read(ROOT/'configs/deposition_fill_then_track_20261007.json'))
        session = dict(schema=PROTOCOL, common_parent=rel(parent), failed_v1_evidence=rel(FAILED_V1),
            failed_v2_evidence=rel(FAILED_V2), scales=list(SCALES),
            selected_scale=.25, source_sha256=template['source_sha256'], administrative_source_sha256=administrative,
            initial_parent_source_sha256=v1, inherited_secant_source_sha256=v2, event_locator_protocol=LOCATOR_PROTOCOL,
            blocked_prefix_protocol=BLOCKED_PREFIX_PROTOCOL,
            numeric_change='SAME_SECANT_PHYSICAL_SOLVES_WITH_EXACT_OUTER_BOX_PREFIX_AND_NO_INEFFECTIVE_CONTRACTION_FALLBACK',
            initial_preflight_statistics=template['administrative_preflight_statistics'], adapter_audit=audit,
            wall_limit_s=wall, requested_wall_limit_s=args.wall_seconds, maximum_dt_s=args.maximum_dt,
            target_time_s=args.target_time, shoulder_observe_ms=args.shoulder_observe_ms,
            maximum_continuation_stages=args.max_stages, started_at=progress['started_at'],
            latest_stop_at=progress['latest_stop_at'], inherited_absolute_stop=frozen.ABSOLUTE_STOP.isoformat(),
            initial_cost_evidence=costs, accepted=False)
        session.update(registered_wall_limits_s=audit['registered_wall_limits_s'],
            qualification_wall_reserve_s=QUALIFICATION_WALL_RESERVE_S,
            unmeasured_box_qualification_cost_basis=template['unmeasured_box_qualification_cost_basis'])
        bridge.shared.write_new(output/'protocol.json', session)
        session_record = {rel(output/'protocol.json'): pool.sha(output/'protocol.json')}
        stages = {}
        update(status='INITIAL_VERIFICATION_COMPLETE', initial_preflight_statistics=template['administrative_preflight_statistics'])
        for scale, label in ((1., 'one'), (.5, 'half'), (.25, 'quarter')):
            if (output/'STOP_REQUESTED').exists() or (output/'STOP_AFTER_STAGE').exists():
                update(status='STOPPED_USER_REQUEST'); break
            if deadline-time.perf_counter() < QUALIFICATION_WALL_RESERVE_S:
                update(status='STOP_QUALIFICATION_WALL_RESERVE', required_wall_reserve_s=QUALIFICATION_WALL_RESERVE_S); break
            scope = 'FIRST_BOX_QUALIFICATION' if scale == 1. else 'EMPIRICAL_COMPLETE_EVENT_PAIR'
            action = cost_decision(common, dt, max(0., deadline-time.perf_counter()), phases, costs, scope=scope, qualification=True)
            if action['action'] != 'RUN':
                update(status=action['action'], next_decision=action); break
            require(action['dt_s'] == dt, 'sensitivity variants cannot change common macro interval')
            child = output/('sensitivity_'+label)
            protocol = stage_protocol(template, common, parent, dt, administrative, session_record, refs, scale, costs,
                                      scope=scope)
            update(status='RUNNING_SAME_PARENT_SENSITIVITY', current_attempt=rel(child), current_parent=rel(parent),
                   current_dt_s=dt, event_tolerance_scale=scale)
            exit_code, elapsed = run_stage(pool, child, protocol, deadline, update, output/'STOP_REQUESTED')
            if exit_code:
                refs, failure = keep_failure(pool, child, template['source_sha256'], common, refs)
                update(status='STOPPED_SENSITIVITY_CHILD_FAILURE', failed_attempt=rel(child), error=failure.get('error'))
                code = 2; break
            meta = frozen.validated_stage(pool, child, template['source_sha256'], common)
            require(source_generation(meta, template) == 4
                and meta.get('implementation_transition') == source_transition(template, common, parent),
                    'new secant result source bridge differs from the exact complete parent')
            sibling_identity(meta, common, parent, scale, dt, template['source_sha256'])
            refs = frozen.bind_stage(pool, child, meta, refs)
            costs.append(cost_evidence(pool, child, meta)); stages[scale] = (child, meta)
            bridge.shared.write_new(output/(label+'_receipt.json'), dict(stage=rel(child), scale=scale,
                common_parent_sha256=common['screened_continuation_sha256'], wall_seconds=elapsed,
                verified=True, selected_for_continuation=False))
            update(status='SENSITIVITY_VARIANT_VERIFIED', validation_variants_completed=len(stages), child_pid=None)
        if len(stages) == 3:
            try:
                report = compare_group(common, parent, stages)
            except Exception as error:
                bridge.shared.write_new(output/'sensitivity_comparison.json', dict(
                    status='FAIL_SAME_PARENT_EVENT_TOLERANCE_SENSITIVITY', schema=PROTOCOL, error=repr(error),
                    stages=[dict(scale=scale, directory=rel(folder),
                        manifest_sha256=pool.sha(folder/'manifest.json'), event_records=meta['event_records'])
                        for scale, (folder, meta) in stages.items()], selected_stage=None, accepted=False))
                pool.sha(output/'sensitivity_comparison.json')
                raise
            bridge.shared.write_new(output/'sensitivity_comparison.json', report)
            refs = frozen.merge(refs, {rel(output/'sensitivity_comparison.json'): pool.sha(output/'sensitivity_comparison.json')})
            observation, state = frozen.observe(common, frames)
            policy = adaptive.initial_state(common['screened_continuation_sha256'], dt, maximum_dt_s=args.maximum_dt)
            parent, meta = stages[.25]; policy, event = selected_pass(policy, meta, state)
            observation, state = frozen.observe(meta, frames)
            admin.atomic_json(output/'deposition_progress.json', observation)
            update(status='SENSITIVITY_PASS_SELECTED_QUARTER', selected_parent=rel(parent), adaptive_state=policy,
                   selected_physical_interval_count=1, last_adaptive_event=event, **admin.snapshot(meta, parent))
            code = continue_selected(args, output, pool, deadline, update, template, administrative,
                session_record, refs, parent, meta, costs, policy, observation, state, frames, phases)
    except Exception as error:
        code = 2; update(status='STOPPED_SESSION_ERROR', error=repr(error), traceback=traceback.format_exc())
    finally:
        try:
            bridge.shared.write_new(output/'final_history_audit.json', pool.final_audit())
        except Exception as error:
            code = 2; update(status='STOPPED_FINAL_HISTORY_VERIFICATION_FAILURE', final_audit_error=repr(error))
        pool.close()
        if kernel is not None:
            kernel.SetThreadExecutionState(0x80000000); kernel.CloseHandle(mutex)
        update(finished=True, exit_code=code, child_pid=None)
        bridge.shared.write_new(output/'final_summary.json', progress)
    print(json.dumps(progress), flush=True)
    return code


def continue_selected(args, output, pool, deadline, update, template, administrative, session_record,
                      refs, parent, meta, costs, policy, observation, state, frames, phases):
    """Continue the selected physical history with the frozen original retry rules."""
    completed = attempt = 0; shoulder_started = None; shoulder_start_completed = None
    while completed < args.max_stages:
        if (output/'STOP_REQUESTED').exists() or (output/'STOP_AFTER_STAGE').exists():
            update(status='STOPPED_USER_REQUEST'); return 0
        if meta['screened_end_time_s'] >= args.target_time-1e-14:
            update(status='STOPPED_STATIONARY_TARGET_REACHED'); return 0
        first_shoulder = observation['first_observed_supplied_under_shoulder_time_s']
        if first_shoulder is not None and shoulder_started is None:
            shoulder_started = first_shoulder; shoulder_start_completed = completed
        if shoulder_started is not None and meta['screened_end_time_s']-shoulder_started >= args.shoulder_observe_ms/1000 and completed-shoulder_start_completed >= 5:
            update(status='STOPPED_SHOULDER_OBSERVATION_WINDOW_COMPLETE'); return 0
        desired = policy['nominal_dt_s']
        action = cost_decision(meta, desired, max(0., deadline-time.perf_counter()), phases, costs)
        if action['action'] != 'RUN':
            update(status=action['action'], next_decision=action); return 0
        dt = min(action['dt_s'], args.target_time-meta['screened_end_time_s'])
        policy = adaptive.begin_attempt(policy, dt, phase_clipped=dt < desired*(1.-1e-12))
        attempt += 1; child = output/f'attempt_{attempt:04d}'
        protocol = stage_protocol(template, meta, parent, dt, administrative, session_record, refs, .25, costs)
        update(status='RUNNING_GPU_CHILD', attempt_count=attempt, current_attempt=rel(child), current_dt_s=dt,
               current_parent=rel(parent), adaptive_state=policy, event_tolerance_scale=.25)
        exit_code, elapsed = run_stage(pool, child, protocol, deadline, update, output/'STOP_REQUESTED')
        row = dict(attempt=attempt, stage=rel(child), parent=rel(parent), dt_s=dt,
                   exit_code=exit_code, wall_seconds=elapsed)
        if exit_code:
            refs, failure = keep_failure(pool, child, template['source_sha256'], meta, refs)
            screen = read(child/'pair_screen.json') if (child/'pair_screen.json').exists() else None
            gates = False
            if admin.classify_failure(failure, screen) == 'REFINE':
                gates = all(read(child/role/'candidate_qualification.json')['status'] ==
                    'PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN' for role in ROLES)
            recovery = frozen.geometry_budget_recovery(pool, child, failure, meta, policy,
                template['source_sha256'], protocol.get('implementation_transition'))
            if recovery is None:
                policy, event = adaptive.resolve_failure(policy, failure, screen, state['parameters']['dx'],
                    state['bore_observation']['particle_ids'], original_candidate_gates_passed=gates)
            else:
                policy, event = recovery
            row['adaptive_event'] = event
            bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json', row)
            if event['action'] == 'REFINE':
                update(status='RETRY_SMALLER_DT_FROM_SAME_COMPLETE_PARENT', adaptive_state=policy, last_adaptive_event=event)
                continue
            update(status='STOPPED_CHILD_FAILURE', failed_attempt=rel(child), error=failure.get('error', 'worker failed'))
            return 2
        next_meta = frozen.validated_stage(pool, child, template['source_sha256'], meta)
        require(source_generation(next_meta, template) == 4 and next_meta['event_tolerance_scale'] == .25
            and next_meta.get('event_locator_protocol') == LOCATOR_PROTOCOL
            and next_meta.get('implementation_transition') == source_transition(template, meta, parent),
                'selected tolerance and registered secant implementation must persist during continuation')
        refs = frozen.bind_stage(pool, child, next_meta, refs); costs.append(cost_evidence(pool, child, next_meta))
        policy, event = adaptive.resolve_pass(policy, next_meta['pair_screen'], state['parameters']['dx'],
            state['bore_observation']['particle_ids'], next_meta['screened_continuation_sha256'], original_candidate_gates_passed=True)
        parent, meta = child, next_meta; completed += 1
        observation, state = frozen.observe(meta, frames)
        admin.atomic_json(output/'deposition_progress.json', observation)
        row.update(verified=True, adaptive_event=event, locked_history_statistics=pool.stats(), **admin.snapshot(meta, parent))
        bridge.shared.write_new(output/f'attempt_{attempt:04d}_receipt.json', row)
        update(status='VERIFIED_PASS_READY_FOR_NEXT_STAGE', completed_stage_count=completed, child_pid=None,
               adaptive_state=policy, last_adaptive_event=event, **admin.snapshot(meta, parent))
    update(status='STOPPED_STAGE_LIMIT')
    return 0


def run_stage(pool, output, protocol, controller_deadline, update, stop_file):
    from afsd_event_box_worker_20261010 import stage_worker
    output.mkdir(); bridge.shared.write_new(output/'protocol.json', protocol)
    protocol_sha = pool.sha(output/'protocol.json')
    pool.seal(protocol['source_sha256']); pool.seal(protocol['reference_sha256'])
    context = multiprocessing.get_context('spawn'); server, child_channel = context.Pipe(duplex=True)
    child = context.Process(target=stage_worker, args=(str(output), child_channel, protocol_sha), daemon=False)
    started = time.perf_counter(); history_before = pool.stats()
    child.start(); child_channel.close(); update(child_pid=child.pid)
    grants = []; interrupted = False; error = None; last_update = started
    try:
        while child.is_alive():
            if time.perf_counter() >= min(started+1080., controller_deadline) or stop_file.exists():
                interrupted = True
                with (output/'termination.log').open('a', encoding='utf-8') as stream:
                    import subprocess
                    subprocess.run(['taskkill','/PID',str(child.pid),'/T','/F'], stdout=stream, stderr=subprocess.STDOUT, timeout=20.)
                child.join(timeout=20.); break
            if len(grants) < 2 and server.poll(.25):
                request = server.recv()
                phase = ('start','end')[len(grants)] if len(grants) < 2 else None
                response = authorize_grant(request, protocol_sha, phase, pool, protocol)
                grants.append(response); server.send(response)
            elif len(grants) == 2:
                child.join(timeout=.25)
            if time.perf_counter()-last_update >= 5.:
                update(); last_update = time.perf_counter()
        child.join(timeout=5.)
        code = child.exitcode if child.exitcode is not None else 124
        require(not child.is_alive(), 'worker could not be terminated')
        pool.verify(protocol['source_sha256']); pool.verify(protocol['reference_sha256'])
        require(code != 0 or len(grants) == 2, 'successful worker lacks both live grants')
    except Exception as exc:
        error = repr(exc); code = 2
        if child.is_alive():
            child.terminate(); child.join(timeout=15.)
    finally:
        server.close()
    manifest = output/'manifest.json'; log = output/'run.log'
    bridge.shared.write_new(output/'supervision.json', dict(exit_code=code, wall_seconds=time.perf_counter()-started,
        source_inputs_still_match=error is None and len(grants)==2, source_error=error, live_history_grants=grants,
        manifest_sha256=bridge.core.sha(manifest) if manifest.exists() else None,
        log_sha256=bridge.core.sha(log) if log.exists() else None, interrupted=interrupted,
        locked_history_before=history_before, locked_history_after=pool.stats(),
        accepted=False, restart_allowed=False, administrative_protocol=SESSION_PROTOCOL))
    return code, time.perf_counter()-started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--wall-seconds', type=float, default=36000.)
    parser.add_argument('--max-stages', type=int, default=1000, help='continuation stages after the sensitivity group; 0 validates only')
    parser.add_argument('--maximum-dt', type=float, default=.000474609375)
    parser.add_argument('--target-time', type=float, default=.22)
    parser.add_argument('--shoulder-observe-ms', type=float, default=5.)
    parser.add_argument('--detach', action='store_true'); parser.add_argument('--detached-child', action='store_true')
    args = parser.parse_args()
    if args.detach or args.detached_child:
        import subprocess
        from start_gpu_resume_detached_20261008 import job_policy, write_new
        output = (ROOT/args.output).resolve()
        require(output.is_relative_to(BASE) and not output.exists(), 'new detached output required')
        if args.detached_child:
            policy = job_policy()
            write_new(output.with_name(output.name+'.bootstrap.json'), dict(**policy, supervisor=rel(Path(__file__)),
                status='PASS_DETACHED_JOB_POLICY' if not policy['kill_on_job_close'] else 'FAIL_JOB_POLICY'))
            require(not policy['kill_on_job_close'], 'detached session retains kill-on-close')
        else:
            require(os.name == 'nt', 'Windows detached launch required')
            forwarded = [v for v in sys.argv[1:] if v != '--detach']
            command = [sys.executable, '-B', str(Path(__file__).resolve()), *forwarded, '--detached-child']
            startup = subprocess.STARTUPINFO(); startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW; startup.wShowWindow = 0
            flags = subprocess.CREATE_BREAKAWAY_FROM_JOB|subprocess.DETACHED_PROCESS|subprocess.CREATE_NEW_PROCESS_GROUP
            env = dict(os.environ, **{k:'1' for k in bridge.shared.THREADS}); env['PYTHONPATH'] = ''
            with output.with_name(output.name+'.stdout.log').open('xb', buffering=0) as stream:
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL, stdout=stream,
                    stderr=subprocess.STDOUT, startupinfo=startup, creationflags=flags, close_fds=True)
            receipt = dict(status='STARTED_SENSITIVITY_SESSION_PENDING_BOOTSTRAP', pid=process.pid, command=command,
                creationflags=flags, source_sha256={rel(Path(__file__)): bridge.core.sha(Path(__file__))})
            write_new(output.with_name(output.name+'.launch.json'), receipt); print(json.dumps(receipt)); return 0
    return run(args)


if __name__ == '__main__':
    raise SystemExit(main())

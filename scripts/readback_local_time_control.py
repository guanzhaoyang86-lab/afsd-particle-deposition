"""Saved-only audit of the bounded local timestep diagnostic and its replay.

No process driver, force/contact/material map or NNPS is imported or called.
The frozen ramp reader audits complete candidates; temporal comparisons are
recomputed independently from saved arrays. A readback PASS can faithfully
record a failed study and never publishes or makes a checkpoint restartable.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np

import bounded_time_control as policy
import readback_ramped_hot_startup as ramp

ROOT = ramp.ROOT
Reader, require, clean, digest = ramp.Reader, ramp.require, ramp.clean, ramp.digest
old = ramp.old
PASS = 'PASS_SAVED_LOCAL_TIME_READBACK_NOT_ACCEPTED'
FAIL = 'FAIL_SAVED_LOCAL_TIME_READBACK_NOT_ACCEPTED'
PRIMARY_PASS = 'PASS_LOCAL_PAIR_DIAGNOSTIC_NOT_ACCEPTED'
PRIMARY_REJECT = 'TEMPORAL_REJECT_DIAGNOSTIC_NOT_ACCEPTED'
REPLAY_PASS = 'PASS_REPLAY_STAGE_DIAGNOSTIC_NOT_ACCEPTED'
REPLAY_REJECT = 'TEMPORAL_REJECT_REPLAY_NOT_QUALIFIED'
PHYSICAL_FAIL = 'FAIL_LOCAL_TIME_PHYSICAL_NOT_ACCEPTED'
REPLAY_SCOPE = 'COMMON_3MS_ORIGIN_INDEPENDENT_HALF_PRIMARY_FINE_GRID'
ORIGIN = ROOT/'outputs/ramped_hot_startup_20261005/v1/chunk_003/fine/step_004/NOT_RESTARTABLE_candidate.npz'
ORIGIN_SHA = '43b2957db88cb3db23f41dadd4a0b6a9c17f3b73aa92ebb4741f51b4ff03775f'
ORIGIN_REVIEW = ROOT/'outputs/ramped_hot_startup_20261005/v1/chunk_003/readback/summary.json'
ORIGIN_REVIEW_SHA = 'ef2beb0c45691ae755c342191a7077c7f3c4c9328eaa5d173fe8694979dad859'
DIAGNOSTIC = dict(accepted=False, restart_allowed=False, accepted_time_advanced_s=0.)


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def folder_for(base, phase, index):
    return base/f"{'attempt' if phase == 'primary' else 'replay'}_{index:03d}"


def project_path(value):
    return ramp.project_path(value)


def diagnostic(record, name):
    require(record.get('accepted') is False and record.get('restart_allowed') is False,
            name+' became accepted or restartable')
    if 'accepted_steps' in record:
        require(record['accepted_steps'] == 0, name+' reports accepted steps')
    if 'accepted_time_advanced_s' in record:
        require(record['accepted_time_advanced_s'] == 0., name+' advances accepted time')


def origin_audit(reader):
    reader.bind(ORIGIN_REVIEW, ORIGIN_REVIEW_SHA)
    review = reader.read(ORIGIN_REVIEW)
    require(review['status'] == ramp.PASS_READBACK and review['temporal_screens_all_pass'],
            'original 3 ms diagnostic prefix is not qualified')
    reader.sources(review['evidence_sha256'])
    require(review['evidence_sha256'].get(rel(ORIGIN)) == ORIGIN_SHA,
            'original readback does not bind the registered 3 ms origin')
    pa, state = reader.checkpoint(reader.bind(ORIGIN, ORIGIN_SHA))
    diagnostic(state, 'origin')
    require(len(pa.x) == 918 and abs(state['time_s']-.003) <= 1e-15 and state['step'] == 12 and
            state['classification'] == ramp.CLASS and state['model_scope'] == ramp.SCOPE,
            'registered independent ramp origin differs')
    old.validate_interval_energy_state(state['interval_energy'])
    old.validate_bore_crossing_state(state['bore_observation'])
    return pa, state


def previous_review(reader, folder):
    path = folder/'readback/summary.json'
    summary = reader.read(path)
    require(summary['status'] == PASS and summary['exit_code'] == 0,
            'previous stage lacks a passing saved-only readback')
    diagnostic(summary, 'previous readback')
    reader.sources(summary['evidence_sha256'])
    return summary


def primary_chain(reader, base, count):
    """Trust only hashed prior readback arithmetic, not unreviewed candidates."""
    snapshots = []
    current = policy.LocalTimeController(ORIGIN_SHA, reader.checkpoint(ORIGIN)[1]['time_s'])
    checkpoint_map = {ORIGIN_SHA: rel(ORIGIN)}
    for index in range(1, count+1):
        folder = folder_for(base, 'primary', index)
        summary = previous_review(reader, folder)
        require(summary['phase'] == 'primary' and summary['index'] == index and
                summary['controller_before'] == current.snapshot(),
                'prior primary controller chain is broken')
        current = policy.LocalTimeController.from_snapshot(summary['controller_after'])
        require(summary['physical_stage_completed'], 'continued after a physical failure')
        for row in reader.read(folder/'candidate_history.json'):
            path = folder/row['role']/row['checkpoint']
            reader.bind(path, row['sha256'])
            checkpoint_map[row['sha256']] = rel(path)
        require(summary['next_parent_sha256'] == current.parent_sha256 and
                summary['next_parent_checkpoint'] == checkpoint_map[current.parent_sha256],
                'prior selected parent path/hash differs')
        require(index == count or current.status == 'READY_DIAGNOSTIC_ONLY',
                'primary stage followed a terminal controller')
        snapshots.append(summary)
    return current, checkpoint_map, snapshots


def prior_context(reader, base, phase, index):
    if phase == 'primary':
        controller, paths, reviews = primary_chain(reader, base, index-1)
        return controller, paths, reviews, 0
    primary_folders = sorted(base.glob('attempt_[0-9][0-9][0-9]'))
    require(1 <= len(primary_folders) <= policy.MAX_PAIR_ATTEMPTS and
            [p.name for p in primary_folders] == [f'attempt_{k:03d}' for k in range(1, len(primary_folders)+1)],
            'replay requires a contiguous bounded primary record')
    controller, paths, reviews = primary_chain(reader, base, len(primary_folders))
    require(controller.status == 'PRIMARY_COMPLETE_INDEPENDENT_REPLAY_REQUIRED',
            'replay started before a complete primary diagnostic')
    offset = 0
    for number in range(1, index):
        summary = previous_review(reader, folder_for(base, 'replay', number))
        require(summary['phase'] == 'replay' and summary['index'] == number and
                summary['controller_before'] == controller.snapshot() == summary['controller_after'] and
                summary['physical_stage_completed'] and summary['temporal_screens_all_pass'] and
                summary['replay_offset'] == offset and not summary['study_qualified'],
                'replay continued beyond failure/completion or lost ancestry')
        offset += summary['candidate_count']
        reader.bind(project_path(summary['last_checkpoint']), summary['last_sha256'])
        paths[summary['last_sha256']] = summary['last_checkpoint']
        reviews.append(summary)
    return controller, paths, reviews, offset


def stage_envelope(reader, base, phase, index, reviews):
    folder = folder_for(base, phase, index)
    protocol = reader.read(folder/'protocol.json')
    manifest = reader.read(folder/'manifest.json')
    supervision = reader.read(folder/'supervision.json')
    require(protocol['phase'] == phase and protocol['index'] == index,
            'stage protocol phase/index differs')
    for name in ('source_sha256', 'reference_sha256'):
        reader.sources(protocol[name])
    require(all(manifest.get(key) == value for key, value in protocol.items()),
            'manifest changed an immutable protocol field')
    for path in (Path(__file__), Path(policy.__file__), Path(ramp.__file__), Path(ramp.prior.__file__), Path(old.__file__)):
        reader.bind(path)
        require(protocol['source_sha256'].get(rel(path)) == reader.evidence[rel(path)],
                'protocol does not bind the actual saved-only reader/policy source: '+rel(path))
    for path in (ORIGIN, ORIGIN_REVIEW):
        require(protocol['reference_sha256'].get(rel(path)) == reader.evidence[rel(path)],
                'protocol lacks registered origin binding')
    if reviews:
        previous = reviews[-1]
        previous_path = folder_for(base, previous['phase'], previous['index'])/'readback/summary.json'
        require(protocol['reference_sha256'].get(rel(previous_path)) == digest(previous_path),
                'protocol does not bind the immediately preceding independent readback')
    reader.bind(folder/'manifest.json', supervision['manifest_sha256'])
    reader.bind(folder/'run.log', supervision['log_sha256'])
    require(supervision['source_inputs_still_match'] is True and supervision['source_error'] is None,
            'source or evidence drift during supervised stage')
    for item, name in ((protocol, 'protocol'), (manifest, 'manifest'), (supervision, 'supervision')):
        diagnostic(item, name)
    count = protocol['candidate_count']
    require(count == (3 if phase == 'primary' else protocol['step_count']),
            'stage candidate reservation differs')
    elapsed = sum(v['supervised_wall_seconds'] for v in reviews)
    previous_count = sum(v['candidate_count'] for v in reviews)
    measured = max([policy.COST_FLOOR_S]+[v['measured_candidate_cost_s'] for v in reviews])
    expected = policy.stage_budget(count, measured, elapsed, previous_count)
    require(all(protocol.get(key) == value for key, value in expected.items()),
            'stage cost forecast/shared budget does not reconstruct')
    require(supervision['internal_limit_s'] == expected['internal_limit_s'] and
            supervision['external_limit_s'] == expected['external_limit_s'],
            'supervision deadline differs from registered budget')
    require(supervision['exit_code'] in (0, 2), 'unfinished or timed-out stage is not qualified by this reader')
    require(math.isfinite(manifest['wall_seconds']) and manifest['wall_seconds'] >= 0. and
            math.isfinite(supervision['wall_seconds']) and supervision['wall_seconds'] >= 0.,
            'nonfinite or negative measured cost')
    require(manifest['wall_seconds'] < expected['internal_limit_s'] and
            supervision['wall_seconds'] < expected['external_limit_s'] and
            elapsed+supervision['wall_seconds'] <= policy.TOTAL_WALL_LIMIT_S,
            'measured stage/shared wall budget exceeded')
    actual = {rel(path): digest(path) for path in folder.rglob('*') if path.is_file() and
              not path.is_relative_to(folder/'readback') and path.name not in ('manifest.json', 'run.log', 'supervision.json')}
    require(manifest['artifacts_sha256'] == actual, 'saved artifact inventory is incomplete or changed')
    reader.sources(actual)
    require(not list(folder.rglob('accepted.npz')) and not list(folder.rglob('UNPUBLISHED_candidate.npz')),
            'diagnostic stage contains publication files')
    return folder, protocol, manifest, supervision, expected


def audit_one(reader, folder, row, previous, previous_state, parent_sha, initial, frozen, dt):
    require(row['step'] == 1 and row['parent_sha256'] == parent_sha and row['dt_s'] == dt,
            'candidate role does not bind its actual immediate parent/duration')
    role = row['role']
    require(isinstance(role, str) and '/' not in role and '\\' not in role,
            'candidate role escaped its stage')
    context = dict(initial=initial, parent=previous_state, p=previous_state['parameters'],
                   prefix=previous_state['interval_energy']['records'], dt=dt)
    case = ramp.candidate_audit(reader, role, folder/role, row,
                                previous, previous_state, context, frozen)
    # The inherited reader checks SHOULDER/TOP endpoint records. Include all
    # bore/base features here so a newly active bore pair cannot be omitted.
    journal = reader.read(folder/role/'step_001/FULL_UNCOMMITTED_candidate.json')
    contact = journal['contact']; features = contact['features']
    endpoint = np.asarray(contact['endpoint_normal_load_N'])
    recorded = {(int(e['particle_id']), e['feature']) for e in row['events'] if e['fraction'] == 1.}
    actual = {(int(case['pa'].particle_id[i]), features[j]) for i, j in np.argwhere(endpoint > 0.)}
    require(recorded == actual, 'all-feature endpoint contact observation differs')
    case['checkpoint'] = rel(folder/role/row['checkpoint'])
    return case


def screen_audit(reader, saved, fine, coarse, parameters, scope):
    diagnostic(saved, 'screen')
    require(saved['comparison_scope'] == scope and saved['time_s'] == fine['state']['time_s'] and
            abs(saved['time_s']-coarse['state']['time_s']) <= 1e-15,
            'comparison scope/shared time differs')
    require(saved['coarse_checkpoint'] == coarse['checkpoint'] and saved['coarse_sha256'] == coarse['row']['sha256'] and
            saved['fine_checkpoint'] == fine['checkpoint'] and saved['fine_sha256'] == fine['row']['sha256'],
            'screen checkpoint/hash binding differs')
    glob, local = old.compare_arrays(fine['pa'], coarse['pa'], sorted(int(v) for v in fine['pa'].particle_id), parameters)
    same = fine['state']['bore_observation']['crossing_ids'] == coarse['state']['bore_observation']['crossing_ids']
    passed = bool(glob['pass'] and local['passed'] and same)
    normalized = policy.normalized_screen_error(glob, local, parameters['dx'])
    require(clean(glob) == saved['global_screen'] and clean(local) == saved['local_screen'] and
            same == saved['crossing_ids_match'] and passed == saved['passed'] and
            normalized == saved['normalized_error'], 'saved temporal comparison does not independently replay')
    return dict(time_s=saved['time_s'], comparison_scope=scope, passed=passed,
                normalized_error=normalized, global_screen=glob, local_particle_count=len(local['rows']),
                local_max_stress_relative=max(r['stress_relative'] for r in local['rows']),
                local_max_plastic_relative=max(r['plastic_relative'] for r in local['rows']),
                crossing_ids_match=same, coarse_checkpoint=coarse['checkpoint'], coarse_sha256=coarse['row']['sha256'],
                fine_checkpoint=fine['checkpoint'], fine_sha256=fine['row']['sha256'])


def primary_audit(reader, folder, protocol, manifest, supervision, controller, paths, initial):
    require(protocol['controller_before'] == controller.snapshot(), 'controller parent snapshot differs')
    require(protocol['parent_sha256'] == controller.parent_sha256 and
            protocol['parent_checkpoint'] == paths[controller.parent_sha256], 'selected primary parent differs')
    proposal = controller.propose()
    require(proposal == protocol['proposal'], 'pure timestep proposal does not replay')
    parent_path = reader.bind(project_path(protocol['parent_checkpoint']), protocol['parent_sha256'])
    parent, parent_state = reader.checkpoint(parent_path)
    rows = reader.read(folder/'candidate_history.json')
    roles = ('coarse', 'fine_001', 'fine_002')
    require([r['role'] for r in rows] == list(roles[:len(rows)]) and len(rows) <= 3,
            'primary role order/count differs')
    cases = {}
    for row in rows:
        role = row['role']
        previous = cases['fine_001'] if role == 'fine_002' else dict(pa=parent, state=parent_state, row=dict(sha256=protocol['parent_sha256']))
        dt = proposal['coarse_dt_s'] if role == 'coarse' else proposal['fine_dt_s']
        cases[role] = audit_one(reader, folder, row, previous['pa'], previous['state'],
            previous['row']['sha256'], initial, protocol['source_sha256'], dt)
        paths[row['sha256']] = cases[role]['checkpoint']
    complete = supervision['exit_code'] == 0
    screens = []
    if complete:
        require(len(rows) == 3, 'successful primary attempt lacks all three full candidates')
        screen = manifest['screen']
        require(manifest['screens'] == [screen] and reader.read(folder/'screen.json') == screen,
                'primary saved screen ledgers differ')
        screens.append(screen_audit(reader, screen, cases['fine_002'], cases['coarse'], parent_state['parameters'], policy.COMPARISON_SCOPE))
        require(manifest['status'] == (PRIMARY_PASS if screen['passed'] else PRIMARY_REJECT),
                'primary outcome relabelled a temporal result')
        controller.resolve(proposal, physical_passed=True, temporal_passed=screen['passed'],
            normalized_error=screen['normalized_error'], fine2_sha256=rows[-1]['sha256'] if screen['passed'] else None)
    else:
        require(manifest['status'] == PHYSICAL_FAIL and manifest.get('error') and not manifest.get('screen'),
                'physical failure evidence/status differs')
        controller.resolve(proposal, physical_passed=False, temporal_passed=False)
    require(manifest['controller_after'] == controller.snapshot(), 'controller resolution does not replay')
    return cases, screens, complete, controller, dict(next_parent_checkpoint=paths[controller.parent_sha256],
        next_parent_sha256=controller.parent_sha256, replay_offset=None, replay_completed_steps=0, study_qualified=False)


def replay_plan_audit(reader, base, controller, plan):
    grid = controller.replay_grid()
    require(len(plan) == len(grid), 'replay plan length differs from qualified primary fine grid')
    primary_rows = []
    for interval in controller.qualified_intervals:
        folder = folder_for(base, 'primary', interval['attempt'])
        rows = reader.read(folder/'candidate_history.json')
        require(len(rows) == 3 and rows[-1]['sha256'] == interval['fine2_sha256'],
                'replay reference is not a qualified primary pair')
        primary_rows.append({r['role']: (folder/r['role']/r['checkpoint'], r) for r in rows})
    for expected, actual in zip(grid, plan):
        require(all(actual.get(k) == v for k, v in expected.items()), 'replay grid alters time/history or diagnostic flags')
        extras = set(actual)-set(expected)
        path, row = primary_rows[expected['primary_interval']-1][f"fine_{expected['primary_fine_step']:03d}"]
        require(actual['reference_checkpoint'] == rel(path) and actual['reference_sha256'] == row['sha256'],
                'replay comparison reference is not the actual qualified primary fine endpoint')
        reader.bind(path, row['sha256'])
        require(extras == {'reference_checkpoint', 'reference_sha256'}, 'unregistered replay plan fields')
    return grid


def replay_audit(reader, base, folder, protocol, manifest, supervision, controller, reviews, offset, initial):
    require(protocol['controller_before'] == controller.snapshot() == manifest['controller_after'],
            'replay changed the primary controller')
    plan = protocol['replay_plan']; replay_plan_audit(reader, base, controller, plan)
    require(protocol['replay_offset'] == offset and protocol['step_count'] == min(4, len(plan)-offset),
            'replay reservation skips/repeats the registered grid')
    previous = reviews[-1] if protocol['index'] > 1 else None
    parent_name = previous['last_checkpoint'] if previous else rel(ORIGIN)
    parent_hash = previous['last_sha256'] if previous else ORIGIN_SHA
    require(protocol['parent_checkpoint'] == parent_name and protocol['parent_sha256'] == parent_hash,
            'independent replay does not start from its own direct parent')
    pa, state = reader.checkpoint(reader.bind(project_path(parent_name), parent_hash))
    rows = reader.read(folder/'candidate_history.json')
    require(len(rows) <= protocol['step_count'] and [r['role'] for r in rows] ==
            [f'replay_{k:03d}' for k in range(1, len(rows)+1)], 'replay roles/count differ')
    cases = {}; screens = []; saved_screens = manifest.get('screens', [])
    for ordinal, row in enumerate(rows):
        planned = plan[offset+ordinal]
        require(abs(state['time_s']-planned['t0_s']) <= 1e-15, 'replay actual parent time is off grid')
        case = audit_one(reader, folder, row, pa, state, parent_hash, initial,
                         protocol['source_sha256'], planned['dt_s'])
        cases[row['role']] = case; pa, state, parent_hash = case['pa'], case['state'], row['sha256']
        require(abs(state['time_s']-planned['t1_s']) <= 1e-15, 'replay endpoint is off grid')
        if planned['shared_primary_fine_endpoint']:
            ref_path = reader.bind(project_path(planned['reference_checkpoint']), planned['reference_sha256'])
            ref_pa, ref_state = reader.checkpoint(ref_path)
            reference = dict(pa=ref_pa, state=ref_state, checkpoint=rel(ref_path), row=dict(sha256=planned['reference_sha256']))
            require(len(saved_screens) > len(screens), 'matching replay endpoint lacks temporal screen')
            require(reader.read(folder/f'screen_{ordinal+1:03d}.json') == saved_screens[len(screens)],
                    'replay saved screen ledgers differ')
            screen = screen_audit(reader, saved_screens[len(screens)], case, reference, state['parameters'], REPLAY_SCOPE)
            screens.append(screen)
            require(screen['passed'] or ordinal == len(rows)-1, 'replay advanced beyond a rejected comparison')
    require(len(saved_screens) == len(screens), 'replay screen ledger contains unbound entries')
    complete = supervision['exit_code'] == 0
    passed = all(s['passed'] for s in screens)
    if complete:
        if passed:
            require(manifest['status'] == REPLAY_PASS and len(rows) == protocol['step_count'],
                    'successful replay stage lacks its complete registered trajectory')
        else:
            require(manifest['status'] == REPLAY_REJECT and screens and not screens[-1]['passed'],
                    'replay temporal rejection was relabelled')
    else:
        require(manifest['status'] == PHYSICAL_FAIL and manifest.get('error') and passed,
                'replay physical failure lacks consistent evidence')
    finished = offset+len(rows)
    return cases, screens, complete, controller, dict(next_parent_checkpoint=None, next_parent_sha256=None,
        replay_offset=offset, replay_completed_steps=finished,
        study_qualified=bool(complete and passed and finished == len(plan)))


def audit(reader, base, phase, index):
    initial, origin_state = origin_audit(reader)
    controller, paths, reviews, offset = prior_context(reader, base, phase, index)
    folder, protocol, manifest, supervision, budget = stage_envelope(reader, base, phase, index, reviews)
    before = controller.snapshot()
    if phase == 'primary':
        cases, screens, complete, controller, extra = primary_audit(reader, folder, protocol, manifest,
            supervision, controller, paths, initial)
    else:
        cases, screens, complete, controller, extra = replay_audit(reader, base, folder, protocol, manifest,
            supervision, controller, reviews, offset, initial)
    rows = reader.read(folder/'candidate_history.json')
    require(manifest['candidate_steps'] == len(rows), 'manifest completed candidate count differs')
    candidate_times = [float(r['wall_seconds']) for r in rows]
    require(all(math.isfinite(v) and v >= 0. for v in candidate_times), 'invalid completed candidate cost')
    measured = max([budget['measured_candidate_cost_s']]+candidate_times+
                   ([supervision['wall_seconds']/len(rows)] if rows else []))
    terminal = list(cases.values())[-1] if cases else None
    reader.finish()
    return dict(status=PASS, phase=phase, index=index, stage_status=manifest['status'],
        physical_stage_completed=complete, physical_stage_exit_code=supervision['exit_code'],
        temporal_screens_all_pass=bool(complete and screens and all(s['passed'] for s in screens)),
        controller_before=before, controller_after=controller.snapshot(), screens=screens,
        candidate_count=len(rows), total_candidate_count=budget['previous_candidate_count']+len(rows),
        candidate_summaries=[dict(role=role, **case['summary']) for role, case in cases.items()],
        last_checkpoint=terminal['checkpoint'] if terminal else protocol['parent_checkpoint'],
        last_sha256=terminal['row']['sha256'] if terminal else protocol['parent_sha256'],
        supervised_wall_seconds=supervision['wall_seconds'], physical_stage_wall_seconds=manifest['wall_seconds'],
        total_supervised_wall_seconds=budget['previous_supervised_wall_s']+supervision['wall_seconds'],
        measured_candidate_cost_s=measured, protocol_sha256=digest(folder/'protocol.json'),
        manifest_sha256=digest(folder/'manifest.json'), source_entry_checks_including_repeats=reader.source_entries,
        **extra)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', required=True)
    parser.add_argument('--phase', choices=('primary', 'replay'), required=True)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--output')
    args = parser.parse_args()
    base = project_path(args.base)
    require(base.is_relative_to(ROOT/'outputs') and base != ROOT/'outputs' and base.is_dir(),
            'existing project output study required')
    require(1 <= args.index <= (policy.MAX_PAIR_ATTEMPTS if args.phase == 'primary' else 8),
            'stage index exceeds registered bound')
    folder = folder_for(base, args.phase, args.index)
    output = project_path(args.output) if args.output else folder/'readback'
    require(output == folder/'readback' and not output.exists(), 'fresh registered saved-only readback folder required')
    output.mkdir(parents=True, exist_ok=False)
    reader = Reader(); started = time.perf_counter()
    report = dict(status='READBACK_RUNNING', base=rel(base), phase=args.phase, index=args.index,
        physical_solver_evaluations=0, nnps_evaluations=0, material_evaluations=0, contact_evaluations=0,
        accepted_steps=0, production_authorized=False, command=[sys.executable, *sys.argv], **DIAGNOSTIC)
    try:
        report.update(audit(reader, base, args.phase, args.index))
        code = 0
    except Exception:
        report.update(status=FAIL, error=traceback.format_exc(), temporal_screens_all_pass=False, study_qualified=False)
        code = 2
    report.update(exit_code=code, wall_seconds=time.perf_counter()-started, evidence_sha256=reader.evidence,
        limitations=['Readback PASS certifies saved evidence consistency and can retain a failed study.',
            'Every checkpoint remains diagnostic, unaccepted and nonrestartable; accepted time advance is zero.',
            'Primary pairs share one immediate parent. The independent replay shares only the 3 ms origin.',
            'Local temporal screening is not whole-history, spatial, full-speed or experimental qualification.'])
    (output/'summary.json').write_text(json.dumps(clean(report), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: report.get(k) for k in ('status', 'phase', 'index', 'stage_status',
        'temporal_screens_all_pass', 'study_qualified', 'wall_seconds', 'error')}, indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())

"""One explicitly licensed local diagnostic window from the 3 ms v1 prefix.

The parent remains NOT_RESTARTABLE. Both local branches preserve its complete
history; none of these diagnostic states become formally accepted or join the
old 342-particle process. The failed v2 branch supplies wall cost only.
"""
from pathlib import Path
import argparse
import json
import math
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import run_afsd_process_framework_20261007 as inherited

ROOT, model, core = inherited.ROOT, inherited.model, inherited.core
require, read, load, rel, verify = inherited.require, inherited.read, inherited.load, inherited.rel, inherited.verify
CONTRACT = ROOT/'docs/PROCESS_CONTACT_WINDOW_CONTRACT_20261007.md'
PARENT_PROBE = ROOT/'outputs/process_framework_20261007/run_v1/probe'
COST_PROBE = ROOT/'outputs/process_framework_20261007/run_v2/probe'
PARENT_SHA = '911f472f82bcf9b6809228fbd46aefbc28b31f6ef4d3988fcefb3041b825af26'
PASS = 'PASS_LOCAL_CONTACT_WINDOW_TIME_SCREEN_DIAGNOSTIC_ONLY'


def checked_stage(folder, expected_exit):
    meta = read(folder/'manifest.json'); supervision = read(folder/'supervision.json')
    require(supervision['exit_code'] == expected_exit and supervision['source_inputs_still_match']
        and supervision['manifest_sha256'] == core.sha(folder/'manifest.json')
        and supervision['log_sha256'] == core.sha(folder/'run.log'), 'stage supervision changed: '+rel(folder))
    verify(meta['source_sha256']); verify(meta['artifacts_sha256'])
    return meta


def source_and_parent():
    parent = checked_stage(PARENT_PROBE, 0)
    require(parent['status'] == inherited.PROBE_PASS and parent['temporal_screens_all_pass']
        and len(parent['pair_screens']) == 6 and all(r['passed'] for r in parent['pair_screens'])
        and parent['accepted'] is False and parent['restart_allowed'] is False
        and parent['accepted_time_advanced_s'] == 0., 'six successful diagnostic parent comparisons required')
    checkpoint = ROOT/parent['candidate_end_checkpoint']
    require(checkpoint == PARENT_PROBE/'fine/step_012/NOT_RESTARTABLE_candidate.npz'
        and core.sha(checkpoint) == parent['candidate_end_sha256'] == PARENT_SHA,
        'only the registered v1 3 ms diagnostic parent is licensed')
    pa, state = load(checkpoint)
    program = model.validate_state(state)
    require(state['step'] == 12 and abs(state['time_s']-.003) <= 1e-14 and len(pa.x) == 1350
        and state['diagnostic_only'] and state['accepted'] is False and state['restart_allowed'] is False,
        'parent diagnostic identity, population or time differs')
    require(state['implementation_source_sha256'] == parent['source_sha256'],
            'parent checkpoint implementation binding differs')
    for k, v in parent['source_sha256'].items():
        require(core.sha(ROOT/k) == v, 'parent frozen source changed: '+k)
    fresh = inherited.source_inputs(ROOT/parent['configuration'], ROOT/parent['component_receipt'])
    require(all(fresh.get(k) == v for k, v in parent['source_sha256'].items()),
            'current physical implementation differs from diagnostic parent')
    fresh.update({rel(Path(__file__)): core.sha(Path(__file__)), rel(CONTRACT): core.sha(CONTRACT)})
    references = dict(parent['artifacts_sha256'])
    for name in ('manifest.json', 'supervision.json', 'run.log'):
        path = PARENT_PROBE/name; references[rel(path)] = core.sha(path)
    cost = checked_stage(COST_PROBE, 2)
    require(cost['status'] == 'FAIL_PROCESS_FRAMEWORK_DIAGNOSTIC_NOT_ACCEPTED'
        and any(not screen['passed'] for screen in cost['pair_screens'])
        and cost['accepted'] is False and cost['restart_allowed'] is False,
        'registered v2 failure must remain failed and unpublished')
    require(cost['source_sha256'] == parent['source_sha256'], 'cost evidence used a different physical implementation')
    rows = read(COST_PROBE/'candidate_history.json')
    measured = max(float(row['wall_seconds']) for row in rows)
    require(math.isfinite(measured) and measured > 0., 'actual finite candidate wall times required')
    # Bind failure evidence too; none of these checkpoints can become a parent.
    references.update(cost['artifacts_sha256'])
    for name in ('manifest.json', 'supervision.json', 'run.log'):
        path = COST_PROBE/name; references[rel(path)] = core.sha(path)
    verify(fresh); verify(references)
    return checkpoint, pa, state, program, fresh, references, measured


def preflight(output, pairs, coarse_dt):
    require(not output.exists(), 'new local-window output required; old evidence cannot be overwritten')
    require(isinstance(pairs, int) and 1 <= pairs <= 21 and 3*pairs <= 64,
            '1..21 pairs and at most64 new candidates required')
    require(math.isfinite(coarse_dt) and 1e-9 <= coarse_dt <= .0005,
            'positive bounded coarse increment in [1 ns, 500 us] required')
    checkpoint, pa, state, program, frozen, references, measured = source_and_parent()
    start = state['time_s']; end = start+pairs*coarse_dt
    require(end <= program.end_s and end > start, 'finite local endpoint within machine program required')
    dts = dict(coarse=coarse_dt, fine=coarse_dt/2)
    counts = dict(coarse=pairs, fine=2*pairs)
    for role in dts:
        t = start
        for _ in range(counts[role]):
            interval = program.step(t, dts[role]); t = interval.t1
    forecast = 30.+1.5*(3*pairs)*measured
    require(forecast <= 450., 'local stage predicted to exceed450s; no automatic reduced alternative')
    return dict(schema='LOCAL_CONTACT_WINDOW_DIAGNOSTIC_LICENSE_20261007',
        source_sha256=frozen, reference_sha256=references,
        diagnostic_parent_checkpoint=rel(checkpoint), diagnostic_parent_sha256=PARENT_SHA,
        diagnostic_parent_time_s=start, diagnostic_parent_step=state['step'], end_time_s=end,
        parent_formal_restart_allowed=False,
        license_scope='ONLY_LOCAL_SAME_PARENT_DIAGNOSTIC_WITH_UNCHANGED_CLASSIFICATION_AND_FULL_PREFIX',
        counts=counts, dts_s=dts, pairs=pairs, particle_count=len(pa.x),
        measured_candidate_wall_s=measured, cost_source=rel(COST_PROBE/'candidate_history.json'),
        failed_v2_state_used_as_parent=False, forecast_s=forecast,
        forecast_formula='30 + 1.5 * 3*pairs * maximum_v2_measured_candidate_wall_s',
        internal_limit_s=540., external_limit_s=600., forecast_limit_s=450.,
        accepted=False, restart_allowed=False, accepted_steps=0, accepted_time_advanced_s=0.,
        classification=model.CLASSIFICATION, old_formal_process_connected=False,
        parent_machine_prefix=state['machine_history'],
        parent_energy_prefix=state['interval_energy']['records'],
        parent_crossing_ids=state['bore_observation']['crossing_ids'])


def verify_prefix(pa, state, protocol):
    model.validate_state(state)
    n = protocol['diagnostic_parent_step']
    require(state['machine_history'][:n] == protocol['parent_machine_prefix']
        and state['interval_energy']['records'][:n] == protocol['parent_energy_prefix'],
        'diagnostic parent machine or energy prefix changed')
    require(set(protocol['parent_crossing_ids']) <= set(state['bore_observation']['crossing_ids']),
            'parent passage history was lost')
    require(len(pa.x) == protocol['particle_count'] and state['classification'] == protocol['classification']
        and state['accepted'] is False and state['restart_allowed'] is False and state['diagnostic_only'],
        'local diagnostic attempted to change parent status or population')


def worker(output):
    protocol = read(output/'protocol.json')
    started = time.perf_counter(); rows, screens, coarse_rows = [], [], []; code = 2
    report = dict(protocol, status='RUNNING_LOCAL_CONTACT_WINDOW_DIAGNOSTIC', candidate_steps=0)
    def timeout():
        try:
            core.write(output/'HARD_TIMEOUT.json', dict(status='TIMEOUT_NOT_ACCEPTED', completed_candidates=len(rows)))
        finally:
            os._exit(124)
    timer = threading.Timer(540., timeout); timer.daemon = True; timer.start()
    deadline = started+525.
    try:
        require(all(os.environ.get(k) == '1' for k in inherited.THREADS), 'single-thread CPU required')
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
        parent = ROOT/protocol['diagnostic_parent_checkpoint']
        require(core.sha(parent) == PARENT_SHA == protocol['diagnostic_parent_sha256'], 'diagnostic parent changed')
        core.write(output/'manifest.json', report)
        measured = protocol['measured_candidate_wall_s']
        for role in ('coarse', 'fine'):
            folder = output/role; folder.mkdir()
            pa, state = load(parent); verify_prefix(pa, state, protocol)
            require(state['time_s'] == protocol['diagnostic_parent_time_s']
                and state['step'] == protocol['diagnostic_parent_step'], 'both paths must load the exact same parent')
            parent_sha = PARENT_SHA
            for number in range(1, protocol['counts'][role]+1):
                require(time.perf_counter()+1.5*measured+20 < deadline, 'remaining measured stage budget exhausted')
                step = folder/f'step_{number:03d}'; step.mkdir()
                item = inherited.complete_candidate(pa, state, protocol['dts_s'][role], step,
                    deadline, measured, protocol['source_sha256'])
                pa, state = item['pa'], item['state']; verify_prefix(pa, state, protocol)
                row = dict(role=role, step=number, diagnostic_global_step=state['step'], time_s=state['time_s'],
                    dt_s=protocol['dts_s'][role], checkpoint=rel(item['checkpoint']),
                    sha256=core.sha(item['checkpoint']), parent_sha256=parent_sha,
                    wall_seconds=item['wall_seconds'], metrics=item['metrics'], accepted=False, restart_allowed=False)
                rows.append(row); parent_sha = row['sha256']
                core.write(output/'candidate_history.json', rows)
                if role == 'coarse':
                    coarse_rows.append(row)
                elif number % 2 == 0:
                    coarse_row = coarse_rows[number//2-1]
                    cp, cs = load(ROOT/coarse_row['checkpoint']); verify_prefix(cp, cs, protocol)
                    require(abs(cs['time_s']-state['time_s']) <= 1e-14, 'local coarse/fine time mismatch')
                    global_screen = inherited.post.compare(pa, cp, state['parameters'])
                    local_screen = inherited.event.local_compare(pa, cp, sorted(int(v) for v in pa.particle_id))
                    matching_ids = state['bore_observation']['crossing_ids'] == cs['bore_observation']['crossing_ids']
                    screen = dict(pair=number//2, time_s=state['time_s'],
                        scope='IDENTICAL_3MS_DIAGNOSTIC_PARENT_TWO_LOCAL_RESOLUTIONS_ALL_PARTICLE_IDS',
                        diagnostic_parent_sha256=PARENT_SHA, coarse_checkpoint_sha256=coarse_row['sha256'],
                        fine_checkpoint_sha256=parent_sha, global_screen=global_screen, local_screen=local_screen,
                        crossing_ids_match=matching_ids,
                        passed=bool(global_screen['pass'] and local_screen['passed'] and matching_ids), accepted=False)
                    screens.append(screen); core.write(output/f'pair_{number//2:03d}_screen.json', screen)
                    require(screen['passed'], 'local contact global/local/passage time screen failed; thresholds unchanged')
                measured = max(measured, item['wall_seconds'])
                report.update(candidate_steps=len(rows), current_role=role, candidate_time_s=state['time_s'], pair_screens=screens)
                core.write(output/'manifest.json', report)
                verify(protocol['source_sha256'])
                print(json.dumps(dict(role=role, step=number, time_s=state['time_s'],
                                      wall_seconds=item['wall_seconds'], accepted=False)), flush=True)
            require(abs(state['time_s']-protocol['end_time_s']) <= 1e-14, 'local diagnostic endpoint not reached')
        report.update(status=PASS, candidate_end_checkpoint=row['checkpoint'], candidate_end_sha256=row['sha256'],
            temporal_screens_all_pass=len(screens) == protocol['pairs'] and all(s['passed'] for s in screens),
            diagnostic_endpoint_s=state['time_s'], interval_energy=state['interval_energy'],
            crossing_ids=state['bore_observation']['crossing_ids'], parent_prefix_unchanged=True,
            physical_accuracy_scope='LOCAL_WINDOW_ONLY_PARENT_PREFIX_AND_COMPLETE_PROCESS_NOT_CONVERGENCE_CERTIFIED')
        verify(protocol['source_sha256']); verify(protocol['reference_sha256']); code = 0
    except Exception as error:
        report.update(status='FAIL_LOCAL_CONTACT_WINDOW_DIAGNOSTIC_NOT_ACCEPTED', error=repr(error),
            traceback=traceback.format_exc(), unresolved='CONTACT_START_TIME_ACCURACY_OR_PHYSICAL_GATE_REQUIRES_FURTHER_WORK')
    finally:
        timer.cancel(); core.write(output/'candidate_history.json', rows)
        report.update(wall_seconds=time.perf_counter()-started, candidate_steps=len(rows), pair_screens=screens,
            peak_working_set_bytes=core.memory_bytes(), accepted=False, restart_allowed=False,
            accepted_steps=0, accepted_time_advanced_s=0., old_formal_process_connected=False)
        report['artifacts_sha256'] = {rel(p): core.sha(p) for p in output.rglob('*') if p.is_file()
            and p.name not in ('manifest.json', 'run.log', 'supervision.json')}
        core.write(output/'manifest.json', report)
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--pairs', type=int, default=4)
    parser.add_argument('--coarse-dt', type=float, default=.000125)
    parser.add_argument('--preflight-only', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args(); output = (ROOT/args.output).resolve()
    require(output.is_relative_to(ROOT/'outputs') and output != ROOT/'outputs', 'new project output directory required')
    if args.worker:
        return worker(output)
    protocol = preflight(output, args.pairs, args.coarse_dt)
    if args.preflight_only:
        print(json.dumps({k: protocol[k] for k in ('diagnostic_parent_checkpoint', 'diagnostic_parent_time_s',
            'end_time_s', 'counts', 'dts_s', 'forecast_s', 'measured_candidate_wall_s')}, indent=2))
        return 0
    output.mkdir(parents=True, exist_ok=False); inherited.write_new(output/'protocol.json', protocol)
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--output', rel(output), '--worker']
    environment = dict(os.environ, **{k: '1' for k in inherited.THREADS}); environment['PYTHONPATH'] = ''
    began = time.perf_counter(); code = 124
    with (output/'run.log').open('x', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=ROOT, env=environment, stdout=stream,
                                  stderr=subprocess.STDOUT, timeout=600.).returncode
        except subprocess.TimeoutExpired:
            stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    source_error = None
    try:
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
    except Exception as error:
        source_error = repr(error); code = 2
    manifest = output/'manifest.json'
    inherited.write_new(output/'supervision.json', dict(command=command, exit_code=code,
        wall_seconds=time.perf_counter()-began, source_inputs_still_match=source_error is None,
        source_error=source_error, manifest_sha256=core.sha(manifest) if manifest.exists() else None,
        log_sha256=core.sha(output/'run.log'), internal_limit_s=540., external_limit_s=600.,
        accepted=False, restart_allowed=False, accepted_time_advanced_s=0.))
    print(json.dumps(dict(exit_code=code, status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST',
                         output=rel(output))))
    return code


if __name__ == '__main__':
    raise SystemExit(main())

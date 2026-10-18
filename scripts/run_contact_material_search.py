"""Supervised bounded contact/material search experiment; no production."""
from pathlib import Path
import argparse
import gzip
import json
import os
import re
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from probe_updated_stage import sha, require, compressed
from run_coupled_afsd import write_json, clean, fingerprint
from run_process_event import setup as parent_setup
from run_nonlinear_continuation import accumulate
from check_process_event import suite_result
from afsd_backend.pysph.contact_force_secant import response_from_trace
from afsd_backend.pysph.contact_material_response import BoundContactMaterialSearch, CLASSIFICATION
from afsd_backend.pysph.process_thermal_trial import ProcessThermalTrial, acceptance
from afsd_backend.pysph.stage_structural_root import solve_stage, METHOD
from afsd_backend.pysph.stage_iteration import initial_candidate, residual_norms, StageFailure
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.process_global_state import commit, state, save, restore
from afsd_backend.pysph.transport_observer import observe
from afsd_backend.pysph.application import memory_bytes

BASE = ROOT/'outputs/contact_material_search_20261004'
LOGS = ROOT/'logs/contact_material_search_20261004'
TEST_GATE = 'tests_review1'
SETTINGS = dict(dt_s=5e-7, maximum_evaluations=32, inner_limit_s=540.,
    contact_force_step_relative=1e-3, contact_response='FROZEN_REJECTED_ANCHOR_SECANT',
    material_response='REFRESHED_CANDIDATE_JOURNAL', beta=.75, history_depth=6)
TESTS = ['test_contact_material_search', 'test_material_response', 'test_material_stage_response',
         'test_process_thermal_trial', 'test_stage_structural_root', 'test_event_response']
NEW_SOURCES = ['scripts/run_contact_material_search.py',
    'scripts/check_process_event.py',
    'src/afsd_backend/pysph/contact_force_secant.py',
    'src/afsd_backend/pysph/contact_material_response.py',
    'docs/CONTACT_MATERIAL_SEARCH_CONTRACT.md']


def sources_now():
    paths = [*NEW_SOURCES]+['tests/review/'+name+'.py' for name in TESTS]
    paths += [p.relative_to(ROOT).as_posix() for p in (ROOT/'src').rglob('*.py')]
    return {s: sha(ROOT/s) for s in sorted(set(paths))}


def artifacts(directory):
    return {f.relative_to(directory).as_posix(): sha(f) for f in directory.rglob('*')
            if f.is_file() and f.name != 'manifest.json'}


def gate(name, sources, status):
    directory = BASE/name
    report = json.loads((directory/'manifest.json').read_text(encoding='utf-8'))
    require(report['status'] == status, 'missing successful '+name+' gate')
    require(report['source_sha256'] == sources, 'stale '+name+' source gate')
    for rel, h in report['artifact_sha256'].items():
        path = (directory/rel).resolve()
        require(path.is_relative_to(directory) and sha(path) == h, 'changed '+name+' artifact '+rel)
    return report


def worker(mode, label):
    directory = BASE/label
    require(not directory.exists(), 'refuse output overwrite')
    directory.mkdir(parents=True)
    start = time.perf_counter(); sources = sources_now(); published = False
    report = dict(classification=CLASSIFICATION, status='RUNNING', mode=mode,
        source_sha256=sources, settings=SETTINGS, accepted_steps=0,
        accepted_progress=False, production_authorized=False,
        backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5')
    timer = None; candidate = evaluation = None; trace = []
    try:
        if mode == 'tests':
            counts = []
            for name in TESTS:
                tick = time.perf_counter()
                command = [sys.executable, '-B', 'tests/review/'+name+'.py']
                result = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=60)
                content = result.stdout+result.stderr
                (directory/(name+'.log')).write_bytes(content)
                count = re.search(rb'Ran (\d+) tests?', content)
                row = dict(test=name, command=command, exit_code=result.returncode,
                    tests=int(count.group(1)) if count else None, wall_seconds=time.perf_counter()-tick)
                parsed = suite_result(name.removeprefix('test_'), content.decode('utf-8', errors='replace'), result.returncode)
                row.update(parsed)
                counts.append(row); print(json.dumps(row), flush=True)
                require(parsed['numerical_gate_pass'], 'component gate failed '+name)
            report.update(status='PASS_NUMERICAL_PREFLIGHT_WITH_RETAINED_ACCURACY_FAIL', tests=counts,
                tests_run=sum(row['tests'] for row in counts),
                tests_passing=sum(row['tests_passing'] for row in counts),
                tests_failed=sum(row['tests_failed'] for row in counts),
                continuous_contact_accuracy='FAIL_RETAINED', all_tests_pass=False)
        else:
            gate(TEST_GATE, sources, 'PASS_NUMERICAL_PREFLIGHT_WITH_RETAINED_ACCURACY_FAIL')
            pa, old, p, watched, inherited, _, parent, anchor, anchor_ref = parent_setup()
            before = particle_hash(pa)
            report.update(parent_reference=parent, numerical_anchor=anchor_ref,
                inherited_source_sha256=inherited, input_sha256=watched,
                accepted_time_s=old['t'], source_time_s=old['t'], particle_count=len(pa.x),
                requested_end_s=old['t']+SETTINGS['dt_s'])
            write_json(directory/'manifest.json', report)
            write_json(directory/'resolved_config.json', dict(parameters=p, settings=SETTINGS))
            write_json(directory/'accepted_history.json', [])
            if mode == 'response':
                B = []; C = []; rows = []
                for tr in anchor['evaluation']['boundary']['traces']:
                    r = response_from_trace(tr, p['friction'])
                    compressed(directory/f"contact_{tr['particle_id']:06d}.json.gz", r)
                    B.append(r['B']); C.append(r['C'])
                    row = dict(particle_id=tr['particle_id'], contact_count=len(tr['contact_rows']),
                        kick_calls=r['kick_calls'], branch_changes=r['branch_changes'], wall_seconds=r['wall_seconds'])
                    rows.append(row)
                    write_json(directory/'response_progress.json', rows)
                    if tr['contact_rows']: print(json.dumps(row), flush=True)
                np.savez_compressed(directory/'response.npz', particle_id=pa.particle_id,
                                    B=np.asarray(B), C=np.asarray(C))
                report.update(status='COMPLETE_NUMERICAL_CONTACT_SECANTS',
                    completed_particles=len(rows), kick_calls=sum(r['kick_calls'] for r in rows),
                    branch_change_particles=sum(bool(r['branch_changes']) for r in rows))
            else:
                response_report = gate('response', sources, 'COMPLETE_NUMERICAL_CONTACT_SECANTS')
                require(response_report['parent_reference'] == parent and response_report['numerical_anchor'] == anchor_ref,
                        'response parent or anchor changed')
                with np.load(BASE/'response/response.npz', allow_pickle=False) as data:
                    require(np.array_equal(data['particle_id'], pa.particle_id), 'response ID ordering mismatch')
                    B = data['B'].copy(); C = data['C'].copy()
                trial = ProcessThermalTrial(pa, p, old['t'])
                parent_arrays = {k: prop.get_npy_array().copy() for k, prop in pa.properties.items()}
                bound = BoundContactMaterialSearch(trial, parent_arrays, old, p, B, C)
                if mode == 'cost':
                    tick = time.perf_counter()
                    candidate = initial_candidate(bound, SETTINGS['dt_s'])
                    evaluation = bound.evaluate(candidate, SETTINGS['dt_s'])
                    proposal = bound.apply(candidate, evaluation['map'])
                    norms, scaled = residual_norms(candidate, evaluation['map'], pa.rho)
                    compressed(directory/'cost_trial.json.gz', dict(candidate=candidate, evaluation=evaluation,
                        numerical_proposal_only=proposal))
                    write_json(directory/'cost_residual.json', dict(norms=norms, scaled=scaled))
                    cost = time.perf_counter()-tick
                    forecast = 1.5*32*cost+30
                    report.update(status='COMPLETE_UNCOMMITTED_COST_PROBE',
                        trial_proposal_journal_wall_seconds=cost, forecast_root_wall_seconds=forecast,
                        root_allowed=forecast <= 450., maximum_scaled_residual=scaled,
                        build_and_apply_wall_seconds=proposal['diagnostics']['build_and_apply_wall_seconds'])
                else:
                    cost = gate('cost', sources, 'COMPLETE_UNCOMMITTED_COST_PROBE')
                    require(cost['parent_reference'] == parent and cost['root_allowed'], 'NOT_RUN_COST_GATE')
                    report['cost_forecast_wall_seconds'] = cost['forecast_root_wall_seconds']
                    def hard_stop():
                        write_json(directory/'HARD_TIMEOUT.json', dict(accepted=False, restart_allowed=False,
                            status='WALL_BUDGET_STOP'))
                        os._exit(124)
                    timer = threading.Timer(540.-(time.perf_counter()-start), hard_stop)
                    timer.daemon = True; timer.start()
                    def journal(row, c, e):
                        nonlocal candidate, evaluation
                        candidate, evaluation = c, e
                        trace.append(row)
                        compressed(directory/'latest_trial.json.gz', dict(candidate=c, evaluation=e))
                        compressed(directory/f'candidate_{len(trace):03d}.json.gz', dict(candidate=c, row=row))
                        write_json(directory/'iterations.json', trace)
                        print(json.dumps(clean(row), allow_nan=False), flush=True)
                    solved = solve_stage(bound, bound, SETTINGS['dt_s'], 32,
                        max(1., 535.-(time.perf_counter()-start)), journal)
                    candidate, evaluation = solved['candidate'], solved['evaluation']
                    checks = acceptance(pa, p, trial, candidate, evaluation)
                    after = commit(pa, p, candidate, evaluation, checks)
                    transport = observe(pa, after, old['t'], evaluation['time_after_s'], p, old['transport_cumulative'])
                    ledger = dict(boundary={}, thermal={}, energy_diagnostic_residual_J=0.)
                    accumulate(ledger, evaluation, checks)
                    all_sources = dict(inherited, **sources)
                    config = fingerprint(dict(parameters=p, settings=SETTINGS, classification=CLASSIFICATION))
                    s = state(after, p, evaluation['time_after_s'], 1, all_sources, config,
                        parent, transport['cumulative'], ledger)
                    save(directory, after, p, s)
                    loaded, _ = restore(directory/'accepted.npz', p, all_sources, config, sha(directory/'accepted.npz'))
                    require(particle_hash(loaded) == particle_hash(after), 'checkpoint roundtrip mismatch')
                    write_json(directory/'acceptance.json', checks)
                    write_json(directory/'transport.json', transport)
                    write_json(directory/'tool_pose.json', s['tool_pose'])
                    compressed(directory/'accepted_trial.json.gz', dict(candidate=candidate, evaluation=evaluation))
                    require(particle_hash(pa) == before, 'parent changed')
                    for rel, h in dict(watched, **all_sources).items():
                        require(sha(ROOT/rel) == h, 'input/source drift before acceptance '+rel)
                    entry = dict(time_s=s['t'], dt_s=SETTINGS['dt_s'], evaluations=len(trace),
                        snapshot='accepted.npz', snapshot_sha256=sha(directory/'accepted.npz'),
                        maximum_scaled_residual=checks['maximum_scaled_residual'])
                    write_json(directory/'accepted_history.json', [entry]); published = True
                    report.update(status='COMPLETE_ONE_ROOT_NOT_PROCESS_VALIDATION', accepted_steps=1,
                        accepted_progress=True, accepted_time_s=s['t'], root_controller=METHOD,
                        last_accepted_snapshot='accepted.npz', total_evaluations=len(trace),
                        temporal_accuracy='NOT_VALIDATED', transport=transport['crossing'])
            require(particle_hash(pa) == before, 'accepted parent mutation')
            for rel, h in dict(watched, **inherited).items():
                require(sha(ROOT/rel) == h, 'old evidence drift '+rel)
        require(sources_now() == sources, 'new source drift')
        code = 0
    except Exception as error:
        traceback.print_exc()
        if isinstance(error, StageFailure):
            candidate, evaluation, trace = error.candidate, error.evaluation, error.trace
        compressed(directory/'REJECTED_candidate.json.gz', dict(accepted=False, restart_allowed=False,
            error=repr(error), candidate=candidate, last_evaluation=evaluation, iterations=trace))
        report.update(status='FAILED_AFTER_PUBLICATION' if published else 'BLOCKED_NO_NEW_ACCEPTED_STATE',
                      error=repr(error), total_evaluations=len(trace)); code = 2
    finally:
        if timer: timer.cancel()
    report.update(wall_seconds=time.perf_counter()-start, peak_process_working_set_bytes=memory_bytes(),
                  artifact_sha256=artifacts(directory))
    write_json(directory/'manifest.json', report, rejected=code != 0)
    print(json.dumps(clean({k: v for k, v in report.items() if 'sha256' not in k}, rejected=code != 0)), flush=True)
    return code


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode', choices=['tests', 'response', 'cost', 'solve'])
    ap.add_argument('--label'); ap.add_argument('--worker', action='store_true')
    a = ap.parse_args(); label = a.label or a.mode
    require(re.fullmatch(r'[a-z0-9_]+', label) is not None, 'simple unique output label required')
    if a.worker: return worker(a.mode, label)
    LOGS.mkdir(parents=True, exist_ok=True)
    log = LOGS/(label+'.log'); receipt = LOGS/(label+'_supervision.json')
    require(not log.exists() and not receipt.exists(), 'refuse log overwrite')
    command = [sys.executable, '-B', str(Path(__file__).resolve()), a.mode, '--label', label, '--worker']
    timeout = dict(tests=240, response=240, cost=90, solve=600)[a.mode]
    tick = time.perf_counter()
    with log.open('wb') as stream:
        try:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = 124
    row = dict(command=command, exit_code=code, wall_seconds=time.perf_counter()-tick,
               timeout_s=timeout, log=log.relative_to(ROOT).as_posix(), log_sha256=sha(log))
    write_json(receipt, row); print(json.dumps(row), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())

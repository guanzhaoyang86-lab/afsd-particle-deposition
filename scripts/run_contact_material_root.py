"""One cost-limited root consuming frozen contact/material preflight evidence.

The original 32-evaluation cost gate failed (465 s > 450 s). This separate
driver reduces the number of evaluations to fit that same 450 s forecast.
No existing source/evidence is revised, no physical tolerance is relaxed.
"""
from pathlib import Path
import json
import math
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import run_contact_material_search as package
from run_contact_material_search import (
    ROOT, BASE, LOGS, SETTINGS, TEST_GATE, CLASSIFICATION, sources_now, gate,
    parent_setup, sha, require, compressed, write_json, clean, fingerprint,
    ProcessThermalTrial, BoundContactMaterialSearch, solve_stage, StageFailure,
    acceptance, commit, state, save, restore, observe, accumulate,
    particle_hash, memory_bytes, artifacts, METHOD)


def evaluation_budget(cost):
    if not math.isfinite(cost) or cost <= 0:
        raise ValueError('positive finite measured cost required')
    cap = min(32, math.floor((450.-30.)/(1.5*cost)))
    if cap < 2:
        raise ValueError('NOT_RUN_COST_NO_USEFUL_BOUNDED_ROOT')
    forecast = 1.5*cap*cost+30.
    return dict(maximum_evaluations=cap, forecast_wall_seconds=forecast,
                launch_limit_s=450., inner_limit_s=540., external_limit_s=600.)


def worker():
    directory = BASE/'root'
    require(not directory.exists(), 'fresh root output required')
    directory.mkdir(parents=True)
    started = time.perf_counter(); timer = None; trace = []
    candidate = evaluation = None; published = False
    report = dict(classification=CLASSIFICATION, accepted_steps=0, accepted_progress=False,
        production_authorized=False, status='RUNNING_UNCOMMITTED')
    try:
        original_sources = sources_now()
        gate(TEST_GATE, original_sources, 'PASS_NUMERICAL_PREFLIGHT_WITH_RETAINED_ACCURACY_FAIL')
        response = gate('response', original_sources, 'COMPLETE_NUMERICAL_CONTACT_SECANTS')
        cost = gate('cost', original_sources, 'COMPLETE_UNCOMMITTED_COST_PROBE')
        budget = evaluation_budget(cost['trial_proposal_journal_wall_seconds'])
        pa, old, p, watched, inherited, _, parent, anchor, anchor_ref = parent_setup()
        require(response['parent_reference'] == cost['parent_reference'] == parent and
                response['numerical_anchor'] == cost['numerical_anchor'] == anchor_ref, 'gate parent/anchor mismatch')
        sources = dict(inherited, **original_sources)
        extra = ['scripts/run_contact_material_root.py', 'tests/review/test_contact_material_budget.py',
                 'docs/CONTACT_MATERIAL_BUDGET_ADDENDUM.md']
        sources.update({rel: sha(ROOT/rel) for rel in extra})
        budget_test = json.loads((LOGS/'budget_test.json').read_text(encoding='utf-8'))
        require(budget_test['exit_code'] == 0 and budget_test['source_sha256'] ==
                {rel: sources[rel] for rel in extra}, 'unverified budget driver')
        settings = dict(SETTINGS, maximum_evaluations=budget['maximum_evaluations'])
        before = particle_hash(pa)
        report.update(settings=settings, budget=budget, source_sha256=sources, input_sha256=watched,
            parent_reference=parent, numerical_anchor=anchor_ref, particle_count=len(pa.x),
            source_time_s=old['t'], accepted_time_s=old['t'], root_controller=METHOD,
            process_accuracy='NOT_VALIDATED; OLD_CONTINUOUS_CONTACT_FAIL_RETAINED',
            backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5')
        write_json(directory/'manifest.json', report)
        write_json(directory/'accepted_history.json', [])
        write_json(directory/'resolved_config.json', dict(parameters=p, settings=settings))
        def timeout():
            write_json(directory/'HARD_TIMEOUT.json', dict(accepted=False, restart_allowed=False))
            os._exit(124)
        timer = threading.Timer(max(1., 540.-(time.perf_counter()-started)), timeout)
        timer.daemon = True; timer.start()
        with np.load(BASE/'response/response.npz', allow_pickle=False) as data:
            require(np.array_equal(data['particle_id'], pa.particle_id), 'contact response ID order')
            B, C = data['B'].copy(), data['C'].copy()
        trial = ProcessThermalTrial(pa, p, old['t'])
        arrays = {k: v.get_npy_array().copy() for k, v in pa.properties.items()}
        bound = BoundContactMaterialSearch(trial, arrays, old, p, B, C)
        def journal(row, c, e):
            nonlocal candidate, evaluation
            candidate, evaluation = c, e
            trace.append(row)
            compressed(directory/'latest_trial.json.gz', dict(candidate=c, evaluation=e))
            compressed(directory/f'candidate_{len(trace):03d}.json.gz', dict(candidate=c, row=row))
            write_json(directory/'iterations.json', trace)
            print(json.dumps(clean(row), allow_nan=False), flush=True)
        solved = solve_stage(bound, bound, settings['dt_s'], budget['maximum_evaluations'],
            max(1., 535.-(time.perf_counter()-started)), journal)
        candidate, evaluation = solved['candidate'], solved['evaluation']
        checks = acceptance(pa, p, trial, candidate, evaluation)
        after = commit(pa, p, candidate, evaluation, checks)
        transport = observe(pa, after, old['t'], evaluation['time_after_s'], p, old['transport_cumulative'])
        ledger = dict(boundary={}, thermal={}, energy_diagnostic_residual_J=0.)
        accumulate(ledger, evaluation, checks)
        config = fingerprint(dict(parameters=p, settings=settings, classification=CLASSIFICATION))
        s = state(after, p, evaluation['time_after_s'], 1, sources, config, parent, transport['cumulative'], ledger)
        save(directory, after, p, s)
        loaded, _ = restore(directory/'accepted.npz', p, sources, config, sha(directory/'accepted.npz'))
        require(particle_hash(loaded) == particle_hash(after), 'checkpoint roundtrip mismatch')
        write_json(directory/'acceptance.json', checks)
        write_json(directory/'transport.json', transport)
        write_json(directory/'tool_pose.json', s['tool_pose'])
        compressed(directory/'accepted_trial.json.gz', dict(candidate=candidate, evaluation=evaluation))
        require(particle_hash(pa) == before, 'parent mutation')
        for rel, h in dict(watched, **sources).items():
            require(sha(ROOT/rel) == h, 'source/input drift before publication '+rel)
        entry = dict(time_s=s['t'], dt_s=settings['dt_s'], evaluations=len(trace), snapshot='accepted.npz',
            snapshot_sha256=sha(directory/'accepted.npz'), maximum_scaled_residual=checks['maximum_scaled_residual'])
        write_json(directory/'accepted_history.json', [entry]); published = True
        report.update(status='COMPLETE_ONE_ROOT_NOT_PROCESS_VALIDATION', accepted_steps=1,
            accepted_progress=True, accepted_time_s=s['t'], last_accepted_snapshot='accepted.npz',
            transport=transport['crossing'], temporal_accuracy='NOT_VALIDATED')
        code = 0
    except Exception as error:
        traceback.print_exc()
        if isinstance(error, StageFailure):
            candidate, evaluation, trace = error.candidate, error.evaluation, error.trace
        compressed(directory/'REJECTED_candidate.json.gz', dict(accepted=False, restart_allowed=False,
            error=repr(error), candidate=candidate, evaluation=evaluation, iterations=trace))
        report.update(status='FAILED_AFTER_PUBLICATION' if published else 'BLOCKED_NO_NEW_ACCEPTED_STATE', error=repr(error))
        code = 2
    finally:
        if timer: timer.cancel()
    report.update(total_evaluations=len(trace), wall_seconds=time.perf_counter()-started,
                  peak_process_working_set_bytes=memory_bytes(), artifact_sha256=artifacts(directory))
    write_json(directory/'manifest.json', report, rejected=code != 0)
    print(json.dumps(clean({k: v for k, v in report.items() if 'sha256' not in k}, rejected=code != 0)), flush=True)
    return code


def main():
    if '--worker' in sys.argv:
        return worker()
    log = LOGS/'root.log'; receipt = LOGS/'root_supervision.json'
    require(not log.exists() and not receipt.exists(), 'refuse root evidence overwrite')
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker']
    start = time.perf_counter()
    with log.open('wb') as stream:
        try:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=600)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = 124
    result = dict(command=command, exit_code=code, wall_seconds=time.perf_counter()-start,
                  timeout_s=600, log_sha256=sha(log))
    write_json(receipt, result); print(json.dumps(result), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())

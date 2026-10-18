"""Single bounded 250+250 ns branch; same parent and end as the 500 ns root."""
from pathlib import Path
import copy
import gzip
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
from run_contact_material_search import (
    ROOT, BASE, LOGS, TEST_GATE, CLASSIFICATION, sources_now, gate, parent_setup,
    sha, require, compressed, write_json, clean, fingerprint, ProcessThermalTrial,
    BoundContactMaterialSearch, solve_stage, StageFailure, acceptance, commit,
    state, save, restore, observe, accumulate, particle_hash, artifacts, memory_bytes)
from afsd_backend.pysph.contact_force_secant import response
from afsd_backend.pysph.axial_actuator import control


def worker():
    out = BASE/'time_check'; require(not out.exists(), 'fresh branch output required')
    out.mkdir(parents=True); start = time.perf_counter(); timer = None
    report = dict(status='RUNNING', accepted_steps=0, production_authorized=False,
        classification='SAME_ENDPOINT_250NS_BRANCH_APPROXIMATE_NOT_PROCESS_VALIDATION')
    candidate = evaluation = None; trace = []; history = []; stage = None
    try:
        original = sources_now()
        gate(TEST_GATE, original, 'PASS_NUMERICAL_PREFLIGHT_WITH_RETAINED_ACCURACY_FAIL')
        cost = gate('cost', original, 'COMPLETE_UNCOMMITTED_COST_PROBE')
        response_cost = gate('response', original, 'COMPLETE_NUMERICAL_CONTACT_SECANTS')['wall_seconds']
        forecast = 1.5*(28*cost['trial_proposal_journal_wall_seconds']+2*response_cost)+30
        require(forecast <= 540, 'NOT_RUN_INITIAL_TIME_CHECK_COST')
        root = json.loads((BASE/'root/manifest.json').read_text(encoding='utf-8'))
        require(root['status'] == 'COMPLETE_ONE_ROOT_NOT_PROCESS_VALIDATION', 'published full step required')
        for rel, h in root['artifact_sha256'].items(): require(sha(BASE/'root'/rel) == h, 'full step drift')
        for rel, h in root['source_sha256'].items(): require(sha(ROOT/rel) == h, 'full step source drift')
        with gzip.open(BASE/'root/accepted_trial.json.gz', 'rt', encoding='utf-8') as stream:
            anchor = json.load(stream)['evaluation']
        pa, old, p, watched, inherited, _, parent, _, _ = parent_setup()
        initial_pa_hash = particle_hash(pa); beginning = old['t']; dt = 2.5e-7
        sources = dict(inherited, **original)
        for rel in ['scripts/run_contact_material_time_check.py', 'docs/CONTACT_MATERIAL_TIME_CHECK_CONTRACT.md']:
            sources[rel] = sha(ROOT/rel)
        config = fingerprint(dict(parameters=p, dt_s=dt, stages=2, maximum_evaluations_per_stage=14,
                                  classification=report['classification']))
        report.update(source_sha256=sources, input_sha256=watched, parent_reference=parent,
            source_time_s=beginning, accepted_time_s=beginning, requested_end_s=beginning+2*dt,
            particle_count=len(pa.x), launch_forecast_seconds=forecast,
            maximum_evaluations_per_stage=14, maximum_total_evaluations=28,
            numeric_anchor_path='outputs/contact_material_search_20261004/root/accepted_trial.json.gz',
            numeric_anchor_sha256=sha(BASE/'root/accepted_trial.json.gz'))
        write_json(out/'manifest.json', report)
        write_json(out/'accepted_history.json', [])
        write_json(out/'resolved_config.json', dict(parameters=p, dt_s=dt, stages=2))
        def hard_stop():
            write_json(out/'HARD_TIMEOUT.json', dict(accepted=False, restart_allowed=False))
            os._exit(124)
        timer = threading.Timer(max(1., 570.-(time.perf_counter()-start)), hard_stop)
        timer.daemon = True; timer.start()
        ledger = dict(boundary={}, thermal={}, energy_diagnostic_residual_J=0.)
        transport_cumulative = old['transport_cumulative']; now = beginning
        for index in (1, 2):
            stage = out/f'stage_{index:03d}'; stage.mkdir()
            tick = time.perf_counter(); B = []; C = []; response_rows = []
            commands = control(pa.z0, now+.5*dt, 0., p['bar_length'], p['actuator_depth'], p['feed'], p['release_time'])
            for tr in anchor['boundary']['traces']:
                i = tr['particle_index']; faces = tr['contact_rows']
                require(pa.particle_id[i] == tr['particle_id'], 'anchor ID mismatch')
                r = response([pa.u[i], pa.v[i], pa.w[i]], pa.m[i], tr['internal_force_N'],
                    [f['normal'] for f in faces], [f['normal_force_N'] for f in faces],
                    [f['wall_velocity'] for f in faces], p['friction'], dt,
                    commands['free_axes'][i], commands['prescribed_velocity'][i],
                    contact_ids=[f['contact_id'] for f in faces])
                B.append(r['B']); C.append(r['C'])
                response_rows.append(dict(particle_id=int(pa.particle_id[i]),
                    branch_changes=r['branch_changes'], kick_calls=r['kick_calls']))
            build = time.perf_counter()-tick
            np.savez_compressed(stage/'response.npz', particle_id=pa.particle_id, B=np.asarray(B), C=np.asarray(C))
            write_json(stage/'response_summary.json', dict(build_seconds=build, rows=response_rows))
            remaining = 2-index
            needed = 1.5*((1+remaining)*14*cost['trial_proposal_journal_wall_seconds']+
                          remaining*max(build, response_cost))+15
            require(time.perf_counter()-start+needed < 570, 'NOT_RUN_REMAINING_TIME_CHECK_COST')
            before = particle_hash(pa); trial = ProcessThermalTrial(pa, p, now)
            arrays = {k: v.get_npy_array().copy() for k, v in pa.properties.items()}
            bound = BoundContactMaterialSearch(trial, arrays, dict(t=now, parameters=p), p, np.asarray(B), np.asarray(C))
            trace = []
            def journal(row, c, e):
                nonlocal candidate, evaluation
                candidate, evaluation = c, e; trace.append(row)
                compressed(stage/'latest_trial.json.gz', dict(candidate=c, evaluation=e))
                compressed(stage/f'candidate_{len(trace):03d}.json.gz', dict(candidate=c, row=row))
                write_json(stage/'iterations.json', trace)
                print(json.dumps(dict(stage=index, **clean(row)), allow_nan=False), flush=True)
            solved = solve_stage(bound, bound, dt, 14, max(1., 565.-(time.perf_counter()-start)), journal)
            candidate, evaluation = solved['candidate'], solved['evaluation']
            checks = acceptance(pa, p, trial, candidate, evaluation)
            after = commit(pa, p, candidate, evaluation, checks)
            transport = observe(pa, after, now, evaluation['time_after_s'], p, transport_cumulative)
            transport_cumulative = transport['cumulative']; accumulate(ledger, evaluation, checks)
            s = state(after, p, evaluation['time_after_s'], index, sources, config, parent, transport_cumulative, ledger)
            save(stage, after, p, s)
            loaded, _ = restore(stage/'accepted.npz', p, sources, config, sha(stage/'accepted.npz'))
            require(particle_hash(loaded) == particle_hash(after), 'time branch checkpoint roundtrip')
            write_json(stage/'acceptance.json', checks); write_json(stage/'transport.json', transport)
            write_json(stage/'tool_pose.json', s['tool_pose'])
            compressed(stage/'accepted_trial.json.gz', dict(candidate=candidate, evaluation=evaluation))
            require(particle_hash(pa) == before, 'time branch parent changed')
            for rel, h in dict(watched, **sources).items(): require(sha(ROOT/rel) == h, 'time branch source/input drift')
            entry = dict(stage=index, time_s=s['t'], dt_s=dt, evaluations=len(trace),
                snapshot=f'stage_{index:03d}/accepted.npz', snapshot_sha256=sha(stage/'accepted.npz'),
                maximum_scaled_residual=checks['maximum_scaled_residual'])
            history.append(entry); write_json(out/'accepted_history.json', history)
            pa = after; now = s['t']
            report.update(accepted_steps=index, accepted_time_s=now, last_accepted_snapshot=entry['snapshot'])
            write_json(out/'manifest.json', report)
        check_end = abs(now-root['accepted_time_s'])
        require(check_end <= 4*np.finfo(float).eps*now, 'different branch endpoint')
        report.update(status='COMPLETE_SAME_ENDPOINT_TWO_STAGES', initial_particle_hash=initial_pa_hash,
            total_root_evaluations=sum(h['evaluations'] for h in history), cumulative_ledger=ledger,
            comparison_status='PENDING_INDEPENDENT_READBACK')
        code = 0
    except Exception as error:
        traceback.print_exc()
        if isinstance(error, StageFailure): candidate, evaluation, trace = error.candidate, error.evaluation, error.trace
        compressed(out/'REJECTED_candidate.json.gz', dict(accepted=False, restart_allowed=False,
            error=repr(error), candidate=candidate, evaluation=evaluation, trace=trace))
        report.update(status='INCOMPLETE_SAME_ENDPOINT_BRANCH', error=repr(error)); code = 2
    finally:
        if timer: timer.cancel()
    report.update(wall_seconds=time.perf_counter()-start, peak_process_working_set_bytes=memory_bytes(),
                  artifact_sha256=artifacts(out))
    write_json(out/'manifest.json', report, rejected=code != 0)
    print(json.dumps(clean({k: v for k, v in report.items() if 'sha256' not in k and k != 'cumulative_ledger'}, rejected=code != 0)), flush=True)
    return code


def main():
    if '--worker' in sys.argv: return worker()
    log = LOGS/'time_check.log'; receipt = LOGS/'time_check_supervision.json'
    require(not log.exists() and not receipt.exists(), 'fresh time branch evidence required')
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker']; start = time.perf_counter()
    with log.open('wb') as stream:
        try:
            r = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=600); code = r.returncode
        except subprocess.TimeoutExpired: code = 124
    row = dict(command=command, exit_code=code, wall_seconds=time.perf_counter()-start, timeout_s=600, log_sha256=sha(log))
    write_json(receipt, row); print(json.dumps(row), flush=True); return code


if __name__ == '__main__': raise SystemExit(main())

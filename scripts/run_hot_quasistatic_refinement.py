"""One registered adaptive startup experiment, preserving every acceptance."""
from pathlib import Path
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
from run_hot_quasistatic_probe import (
    ROOT, LOGS, BASE, CLASSIFICATION, parameters, sources, test_gate, write,
    sha, clean, summary, journal, verify_snapshot, memory_bytes, pose,
    independent_initial_state, QuasistaticTrial, newton_direction,
    energy_diagnostics, accepted_copy, checkpoint, vtk)

TARGET = .001
INITIAL_DT = .0005/32
MAX_ATTEMPTS = 9
DIRECTORY = BASE/'refinement'


def frozen_sources():
    paths = [Path(__file__).resolve(), ROOT/'docs/HOT_QUASISTATIC_REFINEMENT_ADDENDUM_20261004.md']
    return dict(sources(), **{str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in paths})


def solve(trial, dt, initial, directory, deadline):
    r = trial.evaluate(initial, dt); evaluations = 1; trace = []
    for iteration in range(25):
        row = dict(iteration=iteration, evaluations=evaluations, **summary(r))
        trace.append(row); write(directory/'iterations.json', trace)
        journal(directory, trial, r, 'latest_UNCOMMITTED_candidate')
        if r['scaled_residual'] <= 1:
            return r, dict(iterations=iteration, evaluations=evaluations)
        if iteration == 24: raise ValueError('Newton iteration budget exhausted')
        if time.perf_counter() > deadline: raise ValueError('startup series time budget exhausted')
        delta, linear = newton_direction(trial, r)
        merit = np.linalg.norm(r['force'][trial.state['free']]); found = False
        for backtrack in range(8):
            alpha = .5**backtrack
            try:
                candidate = trial.evaluate(r['displacement']+alpha*delta, dt)
                evaluations += 1
                candidate_merit = np.linalg.norm(candidate['force'][trial.state['free']])
                row.setdefault('line_search', []).append(dict(alpha=alpha, merit_N=float(candidate_merit)))
                if candidate_merit < (1-1e-4*alpha)*merit:
                    r = candidate; found = True; break
            except (ValueError, ArithmeticError) as error:
                evaluations += 1
                row.setdefault('line_search', []).append(dict(alpha=alpha, rejected=str(error)))
            if evaluations >= 100 or time.perf_counter() > deadline: break
        row['linear'] = linear; write(directory/'iterations.json', trace)
        if not found: raise ValueError('complete residual line search failed')
        if evaluations >= 100: raise ValueError('residual budget exhausted')


def worker():
    DIRECTORY.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter(); accepted_history = []; attempts = []
    def hard_stop():
        write(DIRECTORY/'HARD_TIMEOUT.json', dict(accepted_history=accepted_history, restart_only_published=True))
        os._exit(124)
    timer = threading.Timer(540, hard_stop); timer.daemon = True; timer.start()
    report = dict(classification=CLASSIFICATION, target_time_s=TARGET, accepted_time_s=0.,
        original_dynamic_time_s=1.3503936329109039e-6, accepted_steps=0, production_authorized=False,
        temporal_accuracy='NOT_VALIDATED_BY_ACCEPTANCE_ALONE')
    try:
        original = sources(); test_gate(original)
        cost = json.loads((BASE/'cost/manifest.json').read_text())
        root = json.loads((BASE/'root/manifest.json').read_text())
        if original != cost['source_sha256'] or original != root['source_sha256']:
            raise ValueError('frozen base version differs from measured evidence')
        root_cost = root['wall_seconds']
        forecast = 1.5*MAX_ATTEMPTS*root_cost+30
        if forecast > 450.: raise ValueError('startup sequence cost forecast exceeds registration')
        p, provenance = parameters()
        if provenance != cost['provenance']: raise ValueError('input drift')
        frozen = frozen_sources()
        report.update(source_sha256=frozen, provenance=provenance, forecast_wall_s=forecast,
                      status='RUNNING_INDEPENDENT_STARTUP', maximum_attempts=MAX_ATTEMPTS)
        write(DIRECTORY/'manifest.json', report); write(DIRECTORY/'accepted_history.json', [])
        pa, state = independent_initial_state(p)
        checkpoint(DIRECTORY/'independent_initial.npz', pa, state)
        verify_snapshot(DIRECTORY/'independent_initial.npz', pa, state)
        dt = INITIAL_DT
        for index in range(MAX_ATTEMPTS):
            if time.perf_counter()-begin+1.5*root_cost+15 > 525.:
                report['stop_reason'] = 'NO_REMAINING_BOUNDED_ROOT_BUDGET'; break
            dt = min(dt, TARGET-state['time_s'])
            if dt < 1e-6: report['stop_reason'] = 'REGISTERED_MINIMUM_DT'; break
            directory = DIRECTORY/f'attempt_{index+1:02d}'
            directory.mkdir()
            trial = QuasistaticTrial(pa, state, p)
            guess = np.column_stack((pa.u, pa.v, pa.w))*dt
            fixed = ~state['free']; guess[fixed] = trial.prescribed(dt)[fixed]
            started = time.perf_counter()
            entry = dict(attempt=index+1, start_time_s=state['time_s'], dt_s=dt, accepted=False)
            try:
                r, metrics = solve(trial, dt, guess, directory, begin+525)
                entry.update(metrics, force_residual_N=r['force_residual_N'])
                energy = energy_diagnostics(trial, r)
                entry.update(energy_raw_J=energy['global_raw_residual_J'],
                    energy_limit_J=energy['registered_energy_limit_J'])
                newpa, newstate = accepted_copy(trial, r)
                newstate.update(source_sha256=frozen, parameter_provenance=provenance, parameters=p)
                staged = dict(newstate, restart_allowed=False)
                checkpoint(directory/'UNPUBLISHED_candidate.npz', newpa, staged)
                verify_snapshot(directory/'UNPUBLISHED_candidate.npz', newpa, staged)
                vtk(directory/'accepted.vtk', newpa, p, trial.bulk, newstate['time_s'], CLASSIFICATION)
                write(directory/'energy.json', dict(mechanical=energy, thermal=r['thermal']['ledger']))
                write(directory/'tool_pose.json', pose(p, newstate['time_s']))
                if frozen_sources() != frozen or parameters() != (p, provenance):
                    raise RuntimeError('source/input drift before publication')
                checkpoint(directory/'accepted.npz', newpa, newstate)
                # Immediately preserve the accepted fact, even if later output fails.
                pa, state = newpa, newstate
                entry.update(accepted=True, time_s=state['time_s'],
                    checkpoint=str((directory/'accepted.npz').relative_to(DIRECTORY)),
                    checkpoint_sha256=sha(directory/'accepted.npz'))
                accepted_history.append(entry.copy())
                write(DIRECTORY/'accepted_history.json', accepted_history)
                fraction = abs(entry['energy_raw_J'])/entry['energy_limit_J']
                dt = min(.0005, dt*(2 if fraction < .25 else 1.25))
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                entry.update(error=repr(error)); dt *= .5
            entry['wall_seconds'] = time.perf_counter()-started
            attempts.append(entry); write(DIRECTORY/'attempts.json', attempts)
            print(json.dumps(clean(entry)), flush=True)
            report.update(accepted_steps=len(accepted_history), accepted_time_s=state['time_s'])
            write(DIRECTORY/'manifest.json', report)
            if state['time_s'] >= TARGET-1e-15:
                report['stop_reason'] = 'TARGET_REACHED'; break
        report.update(status=('COMPLETE_REGISTERED_STARTUP_NOT_TRANSPORT_VALIDATION'
                               if state['time_s'] >= TARGET-1e-15 else 'BOUNDED_STARTUP_TARGET_NOT_REACHED'),
            accepted_steps=len(accepted_history), accepted_time_s=state['time_s'],
            crossing_ids=state['crossing_ids'],
            crossing_mass_kg=sum(float(m) for ident,m in zip(pa.particle_id,pa.m) if int(ident) in state['crossing_ids']),
            max_plastic_strain=float(pa.ep.max()), temperature_min_K=float(pa.temperature.min()),
            temperature_max_K=float(pa.temperature.max()),
            max_net_displacement_m=float(np.max(np.linalg.norm(np.column_stack((pa.x-pa.x0,pa.y-pa.y0,pa.z-pa.z0)),axis=1))),
            estimated_wall_s_per_physical_s=((time.perf_counter()-begin)/state['time_s'] if state['time_s'] else None))
        code = 0 if state['time_s'] >= TARGET-1e-15 else 3
    except Exception as error:
        traceback.print_exc(); report.update(status='FAILED_WITH_PUBLISHED_ACCEPTANCE_RETAINED', error=repr(error))
        code = 2
    finally:
        timer.cancel()
    report.update(accepted_steps=len(accepted_history),
        accepted_time_s=(accepted_history[-1]['time_s'] if accepted_history else 0.),
        wall_seconds=time.perf_counter()-begin, peak_process_working_set_bytes=memory_bytes())
    write(DIRECTORY/'manifest.json', report)
    return code


if __name__ == '__main__':
    if '--worker' in sys.argv: raise SystemExit(worker())
    log = LOGS/'refinement.log'
    if log.exists() or DIRECTORY.exists(): raise ValueError('refuse evidence overwrite')
    cmd = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker']
    start = time.perf_counter()
    with log.open('w', encoding='utf-8') as output:
        try:
            code = subprocess.run(cmd, cwd=ROOT, stdout=output, stderr=subprocess.STDOUT, timeout=600.).returncode
        except subprocess.TimeoutExpired: code = 124
    write(LOGS/'refinement_supervision.json', dict(command=cmd, exit_code=code,
        external_wall_s=time.perf_counter()-start, external_limit_s=600.))
    print(json.dumps(dict(exit_code=code, log=str(log))))
    raise SystemExit(code)

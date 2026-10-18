"""Registered 16 x 0.25 ms replay from the same 1 ms pilot parent.

No physical equations change. This is one bounded temporal comparison, not
permission to continue a long deposition process. --validate-only is read-only.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import run_discrete_contact_process as core
import run_discrete_post_start as post

BASE = core.ROOT/'outputs/discrete_contact_process_20261004/v3'
DIRECTORY = BASE/'post_start/pilot_refinement'
CONTRACT = core.ROOT/'docs/DISCRETE_PILOT_REFINEMENT_CONTRACT_20261004.md'
PARENT_HASH = 'c4b07b8211b09e3e539000b4c7a75a8e471813723d1f1c7abb67da0be9327f8d'
COUNT, DT = 16, .00025


def require(condition, message):
    if not condition:
        raise ValueError(message)


def preflight():
    """Validate actual accepted metadata and cost; no material/contact evaluation."""
    baseline = core.frozen_sources()
    parameters, provenance = core.parameters()
    cost = post.read(BASE/'cost/manifest.json')
    tests = post.read(BASE/'tests/manifest.json')
    require(cost['status'] == 'PASS_MEASURED_COST_NOT_PHYSICAL_VALIDATION', 'cost gate absent')
    require(tests['status'] == 'PASS_SELECTED_COMPONENT_TESTS_NOT_FULL_REGRESSION', 'tests absent')
    for meta in (cost, tests):
        require(meta['source_sha256'] == baseline and meta['provenance'] == provenance,
                'component/cost source or inputs differ')
    parent_dir, coarse_dir = BASE/'post_start/halves', BASE/'post_start/pilot'
    parent_meta, coarse_meta = post.read(parent_dir/'manifest.json'), post.read(coarse_dir/'manifest.json')
    require(parent_meta['status'] == 'PASS_POST_START_TIME_SCREEN' and
            parent_meta['temporal_comparison']['pass'], 'parent time screen absent')
    require(coarse_meta['status'] == 'PASS_BOUNDED_5MS_PILOT_NOT_TRANSPORT_VALIDATION',
            'coarse pilot has not successfully finished')
    for meta in (parent_meta, coarse_meta):
        require(meta['baseline_source_sha256'] == baseline and meta['provenance'] == provenance,
                'parent/coarse source or inputs differ')
        post.check_sources(meta['source_sha256'])
    parent_history, coarse_history = post.read(parent_dir/'accepted_history.json'), post.read(coarse_dir/'accepted_history.json')
    require(len(parent_history) == parent_meta['accepted_steps'] == 2, 'parent local steps differ')
    require(len(coarse_history) == coarse_meta['accepted_steps'] == 8, 'coarse pilot steps differ')
    parent_file = parent_dir/parent_history[-1]['checkpoint']
    coarse_file = coarse_dir/coarse_history[-1]['checkpoint']
    parent_sha, coarse_sha = core.sha(parent_file), core.sha(coarse_file)
    require(parent_sha == parent_history[-1]['sha256'] == coarse_meta['parent_sha256'] == PARENT_HASH,
            'pilot branches do not share the registered parent')
    require(core.ROOT/coarse_meta['parent_checkpoint'] == parent_file, 'coarse parent path differs')
    require(coarse_sha == coarse_history[-1]['sha256'], 'coarse endpoint hash differs')
    pa, state = core.restore(parent_file)
    coarse, coarse_state = core.restore(coarse_file)
    for array, saved in ((pa, state), (coarse, coarse_state)):
        require(saved['classification'] == core.CLASSIFICATION, 'model classification differs')
        require(saved['source_sha256'] == baseline and saved['parameters'] == parameters and
                saved['parameter_provenance'] == provenance, 'checkpoint source or input mismatch')
        core.DiscreteContactTrial(array, saved, parameters)  # validates serialized histories only
    require(state['step'] == 6 and abs(state['time_s']-.001) <= 1e-15 and
            abs(parent_meta['accepted_time_s']-.001) <= 1e-15 and
            parent_history[-1]['time_s'] == state['time_s'], 'parent time/global step mismatch')
    require(coarse_state['step'] == 14 and abs(coarse_state['time_s']-.005) <= 1e-15 and
            abs(coarse_meta['accepted_time_s']-.005) <= 1e-15 and
            abs(coarse_meta['starting_time_s']-.001) <= 1e-15, 'coarse time/global step mismatch')
    require(coarse_state['interval_energy']['records'][:6] == state['interval_energy']['records'],
            'coarse pilot did not inherit the same cumulative energy prefix')
    require(np.array_equal(pa.particle_id, coarse.particle_id) and np.array_equal(pa.m, coarse.m),
            'parent/coarse particle identity or mass differs')
    require(state['interval_energy'] == parent_meta['interval_energy'] and
            coarse_state['interval_energy'] == coarse_meta['interval_energy'],
            'manifest/checkpoint interval differs')
    per_root = max(float(row['wall_seconds']) for row in parent_history)
    require(all(abs(row['dt_s']-DT) <= 1e-15 for row in parent_history), 'cost came from another dt')
    forecast = 30.+1.5*COUNT*per_root
    require(forecast <= 450., 'measured refinement forecast exceeds 450 seconds')
    extra_paths = (Path(__file__).resolve(), CONTRACT, Path(post.__file__).resolve(), post.CONTRACT)
    extra = {str(path.relative_to(core.ROOT)).replace('\\', '/'): core.sha(path) for path in extra_paths}
    frozen = dict(baseline, **extra)
    post.check_sources(frozen)
    return dict(baseline=baseline, parameters=parameters, provenance=provenance, cost=cost,
                parent_file=parent_file, parent_sha=parent_sha, pa=pa, state=state,
                coarse_file=coarse_file, coarse_sha=coarse_sha, coarse=coarse,
                coarse_state=coarse_state, per_root=per_root, forecast=forecast,
                extra=extra, frozen=frozen)


def worker():
    DIRECTORY.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter()
    history = []
    report = dict(classification=core.CLASSIFICATION, stage='pilot_refinement',
                  status='RUNNING_PILOT_TIME_REFINEMENT_NOT_TRANSPORT',
                  accepted_steps=0, accepted_time_s=None, production_authorized=False)
    def stop():
        core.write(DIRECTORY/'HARD_TIMEOUT.json', dict(accepted_history=history,
            published_checkpoint_paths=[str(p.relative_to(DIRECTORY)) for p in DIRECTORY.glob('step_*/accepted.npz')]))
        os._exit(124)
    timer = threading.Timer(540., stop)
    timer.daemon = True
    timer.start()
    try:
        context = preflight()
        p, provenance = context['parameters'], context['provenance']
        baseline, frozen, extra = context['baseline'], context['frozen'], context['extra']
        pa, state = context['pa'], context['state']
        per_root = context['per_root']
        report.update(source_sha256=frozen, baseline_source_sha256=baseline,
            provenance=provenance, parent_checkpoint=str(context['parent_file'].relative_to(core.ROOT)),
            parent_sha256=context['parent_sha'], coarse_endpoint=str(context['coarse_file'].relative_to(core.ROOT)),
            coarse_sha256=context['coarse_sha'], starting_time_s=state['time_s'],
            accepted_time_s=state['time_s'], starting_global_step=state['step'],
            requested_time_s=.005, requested_dts_s=[DT]*COUNT, particle_count=len(pa.x),
            forecast_wall_s=context['forecast'], forecast_root_cost_s=per_root,
            numerical_thread_environment={key: os.environ.get(key) for key in
                ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS')})
        core.write(DIRECTORY/'manifest.json', report)
        core.write(DIRECTORY/'accepted_history.json', history)
        core.checkpoint(DIRECTORY/'continuation_initial.npz', pa, state)
        initial, initial_state = core.restore(BASE/'quarters/independent_initial.npz')
        core.checkpoint(DIRECTORY/'independent_initial.npz', initial, initial_state)
        for number in range(1, COUNT+1):
            require(time.perf_counter()-begin+1.5*per_root+15 <= 525., 'remaining cost budget exhausted')
            stepdir = DIRECTORY/f'step_{number:03d}'
            stepdir.mkdir()
            trial = core.DiscreteContactTrial(pa, state, p)
            guess = np.column_stack((pa.u, pa.v, pa.w))*DT
            fixed = ~state['free']
            guess[fixed] = trial.prescribed(DT)[fixed]
            started = time.perf_counter()
            result, metrics = core.solve(trial, DT, guess, stepdir, begin+525., context['cost']['maximum_iterations'])
            geometry = post.geometry_audit(trial, result)
            core.write(stepdir/'geometry_audit.json', geometry)
            require(not geometry['transient_contacts_missed_by_endpoints'], 'unresolved interior penalty-shell contact')
            require(not geometry['new_non_bore_features'], 'new finite feature requires separate local time screen')
            next_pa, next_state, fresh = core.accepted_copy(trial, result)
            next_state.update(controller_source_sha256=extra,
                              continuation_parent_sha256=context['parent_sha'])
            post.check_sources(frozen)
            require(core.sha(context['parent_file']) == context['parent_sha'], 'parent checkpoint drift')
            require(core.sha(context['coarse_file']) == context['coarse_sha'], 'coarse endpoint drift')
            core.publish(stepdir, next_pa, next_state, trial, fresh, baseline, provenance)
            # Record the already-published fact before post-publication audits.
            # A later source/readback error must not erase this accepted prefix.
            pa, state = next_pa, next_state
            row = dict(step=number, global_step=state['step'], time_s=state['time_s'], dt_s=DT,
                checkpoint=f'step_{number:03d}/accepted.npz', **metrics,
                sha256=core.sha(stepdir/'accepted.npz'), wall_seconds=time.perf_counter()-started,
                force_residual_N=fresh['force_residual_N'], energy_raw_J=state['energy']['global_raw_residual_J'])
            history.append(row)
            core.write(DIRECTORY/'accepted_history.json', history)
            report.update(accepted_steps=len(history), accepted_time_s=state['time_s'],
                interval_energy=state['interval_energy'], crossing_ids=state['crossing_ids'],
                crossing_mass_kg=state['crossing_mass_kg'], max_plastic_strain=float(pa.ep.max()))
            core.write(DIRECTORY/'manifest.json', report)
            post.check_sources(frozen)
            _, saved = core.restore(stepdir/'accepted.npz')
            require(saved['controller_source_sha256'] == extra and
                    saved['continuation_parent_sha256'] == context['parent_sha'], 'controller provenance lost')
            require(saved['interval_energy'] == state['interval_energy'], 'published cumulative ledger differs')
            print(json.dumps(core.clean(row)), flush=True)
            per_root = max(per_root, row['wall_seconds'])
        require(state['step'] == 22 and abs(state['time_s']-.005) <= 1e-15, 'final time/global step mismatch')
        require(abs(context['coarse_state']['time_s']-state['time_s']) <= 1e-15, 'comparison time differs')
        comparison = post.compare(pa, context['coarse'], p)
        core.write(DIRECTORY/'temporal_comparison.json', comparison)
        report.update(temporal_comparison=comparison,
            status='PASS_1_TO_5MS_TIME_SCREEN_NOT_TRANSPORT' if comparison['pass'] else 'FAIL_1_TO_5MS_TIME_SCREEN',
            physical_accuracy='ENGINEERING_TWO_STEP_SIZE_SCREEN_NOT_EXPERIMENTAL_VALIDATION')
        post.check_sources(frozen)
        require(core.sha(context['coarse_file']) == context['coarse_sha'], 'coarse endpoint changed during comparison')
        code = 0 if comparison['pass'] else 3
    except Exception as error:
        traceback.print_exc()
        report.update(status='FAILED_WITH_ACCEPTED_FACTS_RETAINED' if history else 'REJECTED_NO_NEW_ACCEPTED_STATE',
                      error=repr(error))
        code = 2
    finally:
        timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin, peak_process_working_set_bytes=core.memory_bytes(),
        published_checkpoint_paths=[str(path.relative_to(DIRECTORY)) for path in DIRECTORY.glob('step_*/accepted.npz')])
    core.write(DIRECTORY/'manifest.json', report)
    print(json.dumps(core.clean({key: report[key] for key in
        ('status', 'accepted_steps', 'accepted_time_s', 'wall_seconds')})), flush=True)
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--validate-only', action='store_true')
    arguments = parser.parse_args()
    if arguments.validate_only:
        result = preflight()
        print(json.dumps(dict(status='PASS_READ_ONLY_PREFLIGHT_NO_PHYSICAL_EVALUATION',
            parent_sha256=result['parent_sha'], starting_time_s=result['state']['time_s'],
            inherited_energy_steps=result['state']['interval_energy']['steps'],
            coarse_time_s=result['coarse_state']['time_s'], forecast_wall_s=result['forecast'])))
        raise SystemExit(0)
    if arguments.worker:
        raise SystemExit(worker())
    logs = core.ROOT/'logs/discrete_contact_process_20261004/v3/post_start'
    logs.mkdir(parents=True, exist_ok=True)
    log = logs/'pilot_refinement.log'
    require(not log.exists() and not DIRECTORY.exists(), 'refuse evidence overwrite')
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker']
    child_env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
    begin = time.perf_counter()
    with log.open('w', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=core.ROOT, stdout=stream,
                stderr=subprocess.STDOUT, timeout=600, env=child_env).returncode
        except subprocess.TimeoutExpired:
            code = 124
    core.write(logs/'pilot_refinement_supervision.json', dict(command=command, exit_code=code,
        external_wall_s=time.perf_counter()-begin, limit_s=600))
    print(json.dumps(dict(stage='pilot_refinement', exit_code=code, log=str(log))))
    raise SystemExit(code)

"""Strictly gated recomputation of the local HHT fine grid's published second half."""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

import numpy as np

import probe_hht_coupled as driver

ROOT = driver.ROOT
EXTRA_SOURCES = ('scripts/restart_hht_probe.py', 'tests/review/test_hht_restart.py',
                 'docs/PROCESS_HHT_RESTART_CONTRACT.md')
PASS_SWEEP = 'PASS_LOCAL_HHT_TEMPORAL_BENCHMARK_NOT_AFSD'
PASS_PREFLIGHT = 'PASS_HHT_COMPONENT_PREFLIGHT_NOT_AFSD'
INTEGRAL_SHAPES = dict(friction_raw=(2,), reservoir_heat=(2,), wall_work=(2,),
                       impulse=(2, 3), algorithmic_filter_work_J=(), root_residual_work_J=())


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--sweep', type=Path, required=True)
    result.add_argument('--preflight', type=Path, required=True)
    result.add_argument('--output', type=Path, required=True)
    result.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    return result


def checked_path(base, relative):
    """Manifest paths must be relative and resolve inside their stated root."""
    driver.require(isinstance(relative, str) and bool(relative), 'invalid manifest path')
    part = Path(relative)
    driver.require(not part.is_absolute() and '..' not in part.parts, 'escaping manifest path')
    target = (Path(base)/part).resolve()
    driver.require(target.is_relative_to(Path(base).resolve()), 'escaping manifest path')
    return target


def verify_manifest(base, manifest, description):
    driver.require(isinstance(manifest, dict) and bool(manifest), 'missing '+description)
    for relative, expected in manifest.items():
        target = checked_path(base, relative)
        driver.require(target.is_file() and driver.sha(target) == expected,
                       description+' drift: '+relative)


def own_inventory():
    return dict(driver.inventory(), **{name: driver.sha(ROOT/name) for name in EXTRA_SOURCES})


def validate_receipt_payloads(sweep, preflight, current_sources, preflight_hash):
    """Pure checks: no filesystem or solver calls, also used by rejection tests."""
    driver.require(sweep.get('classification') == driver.CLASS and sweep.get('status') == PASS_SWEEP
                   and sweep.get('mode') == 'sweep', 'successful local three-grid sweep required')
    driver.require(sweep.get('production_authorized') is False and sweep.get('afsd_progress') is False
                   and sweep.get('new_afsd_steps') == 0, 'non-AFSD receipt required')
    driver.require(preflight.get('status') == PASS_PREFLIGHT, 'successful component preflight required')
    driver.require(sweep.get('receipts', {}).get('preflight_sha256') == preflight_hash,
                   'sweep does not bind supplied preflight')
    tested = preflight.get('source_sha256', {})
    for name, expected in current_sources.items():
        driver.require(tested.get(name) == expected, 'stale/missing preflight source '+name)
    driver.require(sweep.get('source_sha256') == current_sources, 'sweep source inventory differs')
    driver.require(sweep.get('comparison', {}).get('passed') is True, 'temporal accuracy sweep failed')
    grids = sweep.get('completed_grids', [])
    driver.require(len(grids) == 3 and [row.get('steps') for row in grids] == [8, 16, 32],
                   'complete 8/16/32 grids required')
    for row in grids:
        duration = row.get('time_s', float('nan'))
        cost = row.get('wall_seconds', float('nan'))
        driver.require(np.isfinite(duration) and abs(duration-2e-6) < 1e-19
                       and np.isfinite(cost) and cost > 0, 'invalid grid time/cost')
    fine_whole_step = grids[-1]['wall_seconds']/32.
    forecast = 2.*16.*fine_whole_step+5.
    driver.require(forecast <= 50., 'NOT_RUN_RESTART_COST_EXCEEDS_REGISTERED_LIMIT')
    return fine_whole_step, forecast


def validate_history(history, artifacts, dt):
    driver.require(isinstance(history, list) and len(history) == 32,
                   'exactly 32 published fine-grid states required')
    for number, entry in enumerate(history, 1):
        expected = f'stage_{number:03d}/accepted.npz'
        driver.require(entry.get('checkpoint') == expected, 'unpublished or unexpected checkpoint path')
        driver.require(abs(entry.get('time_s', float('nan'))-number*dt) < 1e-19,
                       'published step time mismatch')
        expected_hash = artifacts.get('grid_032/'+expected)
        driver.require(isinstance(expected_hash, str) and entry.get('checkpoint_sha256') == expected_hash,
                       'checkpoint not bound into sweep artifacts')
        driver.require(isinstance(entry.get('state_hash'), str) and len(entry['state_hash']) == 64,
                       'published semantic state hash missing')
        driver.require('grid_032/'+f'stage_{number:03d}/integrals.json' in artifacts,
                       'published integral artifact missing')


def parse_integrals(raw):
    driver.require(isinstance(raw, dict) and set(raw) == set(INTEGRAL_SHAPES),
                   'exact cumulative integral inventory required')
    result = {}
    for name, shape in INTEGRAL_SHAPES.items():
        value = np.asarray(raw[name])
        driver.require(value.dtype.kind in 'iuf' and value.shape == shape and np.isfinite(value).all(),
                       'invalid cumulative integral '+name)
        result[name] = float(value) if shape == () else np.array(value, dtype=float, copy=True)
    return result


def integral_comparison(actual, expected):
    first, second = parse_integrals(actual), parse_integrals(expected)
    differences = {name: float(np.max(np.abs(np.asarray(first[name])-second[name])))
                   for name in INTEGRAL_SHAPES}
    return dict(bitwise_equal=all(np.array_equal(first[name], second[name])
                                for name in INTEGRAL_SHAPES), maximum_absolute_differences=differences)


def load_published(grid, entry, parameters):
    path = checked_path(grid, entry['checkpoint'])
    driver.require(driver.sha(path) == entry['checkpoint_sha256'], 'published checkpoint changed')
    pa, state = driver.restore(path)
    pa, state = driver.restore_checkpoint(pa, state, parameters)
    driver.require(state['state_hash'] == entry['state_hash']
                   and state['time_s'] == entry['time_s'], 'checkpoint differs from publication')
    driver.require(len(pa.x) == 27 and np.array_equal(pa.particle_id, np.arange(1, 28)),
                   'local 27-particle fixture only')
    return pa, state


def compare_candidate(pa, state, integrals, expected_pa, expected_state, expected_integrals):
    driver.require(set(pa.properties) == set(expected_pa.properties), 'particle property inventory mismatch')
    differences = {}
    for name in sorted(pa.properties):
        actual = pa.properties[name].get_npy_array()
        expected = expected_pa.properties[name].get_npy_array()
        driver.require(actual.shape == expected.shape and actual.dtype == expected.dtype,
                       'particle property layout mismatch: '+name)
        differences[name] = float(np.max(np.abs(actual-expected))) if actual.size else 0.
    ledger = integral_comparison(integrals, expected_integrals)
    result = dict(particle_hash_equal=driver.particle_hash(pa) == driver.particle_hash(expected_pa),
        semantic_state_hash_equal=state['state_hash'] == expected_state['state_hash'],
        maximum_property_differences=differences,
        maximum_acceleration_difference_m_s2=float(np.max(np.abs(
            state['acceleration_m_s2']-expected_state['acceleration_m_s2']))),
        maximum_force_difference_N=float(np.max(np.abs(state['force_N']-expected_state['force_N']))),
        integral_comparison=ledger)
    result['passed'] = result['particle_hash_equal'] and result['semantic_state_hash_equal'] and ledger['bitwise_equal']
    return result


def run_worker(arguments):
    out = arguments.output.resolve()
    driver.require(out.is_relative_to(ROOT) and not out.exists(), 'fresh project output required')
    out.mkdir(parents=True)
    started = time.perf_counter(); history = []; state = None; comparisons = []
    recomputation_started = False
    report = dict(status='NOT_RUN_GATES_PENDING', classification='LOCAL_HHT_RECOMPUTED_RESTART_NOT_AFSD',
        production_authorized=False, afsd_progress=False, new_afsd_steps=0)
    def timeout():
        driver.write(out/'HARD_TIMEOUT.json', dict(status='WALL_BUDGET_STOP', accepted=False,
                     restart_allowed=False, publication='accepted_history.json is sole authority'))
        os._exit(124)
    timer = threading.Timer(55., timeout); timer.daemon = True; timer.start()
    try:
        sweep_path, preflight_path = arguments.sweep.resolve(), arguments.preflight.resolve()
        driver.require(sweep_path.is_relative_to(ROOT) and preflight_path.is_relative_to(ROOT),
                       'local project receipts required')
        sweep, preflight = driver.read(sweep_path), driver.read(preflight_path)
        sources = own_inventory(); preflight_hash = driver.sha(preflight_path)
        unit_cost, forecast = validate_receipt_payloads(sweep, preflight, sources, preflight_hash)
        verify_manifest(ROOT, preflight['source_sha256'], 'preflight source')
        verify_manifest(preflight_path.parent, preflight['artifact_sha256'], 'preflight artifact')
        artifacts = sweep['artifact_sha256']; sweep_dir = sweep_path.parent
        verify_manifest(sweep_dir, artifacts, 'sweep artifact')
        driver.require(not any(sweep_dir.rglob('REJECTED.json'))
                       and not any(sweep_dir.rglob('HARD_TIMEOUT.json')), 'failed sweep artifact present')
        config_path = sweep_dir/'resolved_config.json'
        driver.require(driver.sha(config_path) == sweep['config_sha256']
                       == artifacts.get('resolved_config.json'), 'configuration hash mismatch')
        initial, parameters = driver.setup()
        driver.require(driver.clean(parameters) == driver.read(config_path), 'registered configuration differs')
        dt = parameters['end_time']/32.
        grid = sweep_dir/'grid_032'; publication = grid/'accepted_history.json'
        driver.require(artifacts.get('grid_032/accepted_history.json') == driver.sha(publication),
                       'fine-grid publication not bound into sweep')
        original = driver.read(publication); validate_history(original, artifacts, dt)
        pa, state = load_published(grid, original[15], parameters)
        driver.require(state['last_event']['kind'] == 'EXPLICIT_RIGHT_LIMIT_CONTROL_FORCE_ACCELERATION_INITIALIZATION',
                       'halfway right-limit control history missing')
        driver.require(state['last_dt_s'] == dt and all(row['wall_velocity_m_s'][2] == -.2
            for row in state['boundary_snapshot']['contacts']), 'incorrect halfway time/control layer')
        integrals = parse_integrals(driver.read(grid/'stage_016/integrals.json'))
        report.update(status='RUNNING_LOCAL_RECOMPUTATION', source_sha256=sources,
            sweep_sha256=driver.sha(sweep_path), preflight_sha256=preflight_hash,
            configuration_sha256=driver.sha(config_path), fine_whole_step_seconds=unit_cost,
            forecast_seconds=forecast, resumed_from=original[15], start_time_s=state['time_s'])
        driver.write(out/'report.json', report)
        recomputation_started = True
        for number in range(17, 33):
            expected_pa, expected_state = load_published(grid, original[number-1], parameters)
            expected_integrals = driver.read(grid/f'stage_{number:03d}/integrals.json')
            driver.require(expected_state['parent_state_hash'] == state['state_hash'],
                           'original published semantic parent chain mismatch')
            candidate, candidate_state, diagnostic = driver.advance(pa, state, parameters,
                driver.boundary_for(parameters, True), dt)
            driver.physical_checks(candidate, pa, initial, candidate_state, diagnostic, parameters)
            candidate_integrals = copy.deepcopy(integrals)
            driver.add_integrals(candidate_integrals, state, candidate_state, diagnostic, parameters)
            comparison = compare_candidate(candidate, candidate_state, candidate_integrals,
                                           expected_pa, expected_state, expected_integrals)
            comparison['original_stage'] = number
            comparisons.append(comparison); driver.write(out/'comparisons.json', comparisons)
            driver.require(comparison['passed'], 'recomputed complete state/integrals differ')
            driver.require(own_inventory() == sources, 'source drift before restart candidate save')
            entry = driver.save_stage(out/f'stage_{number:03d}', candidate, candidate_state,
                                      parameters, diagnostic, candidate_integrals)
            driver.require(own_inventory() == sources, 'source drift before restart publication')
            next_history = history+[entry]
            driver.write(out/'accepted_history.json', next_history)
            pa, state, integrals, history = candidate, candidate_state, candidate_integrals, next_history
        verify_manifest(sweep_dir, artifacts, 'sweep artifact after recomputation')
        driver.require(driver.sha(sweep_path) == report['sweep_sha256']
                       and driver.sha(preflight_path) == preflight_hash
                       and own_inventory() == sources, 'receipt/source drift after recomputation')
        driver.require(len(history) == 16 and abs(state['time_s']-2e-6) < 1e-19,
                       'incomplete common-end-time restart')
        report.update(status='PASS_LOCAL_RECOMPUTED_RESTART_NOT_AFSD', recomputed_steps=16,
                      complete_state_and_integrals_bitwise_equal=True)
        code = 0
    except Exception as error:
        report.update(status='FAIL_LOCAL_RECOMPUTED_RESTART' if recomputation_started else 'NOT_RUN_RESTART_GATE_REJECTED',
                      error=repr(error), published_steps=len(history))
        driver.write(out/'REJECTED.json', dict(accepted=False, restart_allowed=False, error=repr(error),
            published_steps=len(history), last_published_time_s=None if state is None else state['time_s'],
            last_published_state_hash=None if state is None else state['state_hash']))
        code = 2
    finally:
        timer.cancel()
    report.update(wall_seconds=time.perf_counter()-started, published_steps=len(history),
        artifact_sha256={path.relative_to(out).as_posix(): driver.sha(path)
                         for path in out.rglob('*') if path.is_file() and path.name != 'report.json'})
    driver.write(out/'report.json', report)
    print(json.dumps(driver.clean(report), allow_nan=False), flush=True)
    return code


def main(argv=None):
    arguments = parser().parse_args(argv)
    if arguments.worker:
        return run_worker(arguments)
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker',
               '--sweep', str(arguments.sweep), '--preflight', str(arguments.preflight),
               '--output', str(arguments.output)]
    try:
        result = subprocess.run(command, cwd=ROOT, timeout=60,
            creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        return result.returncode
    except subprocess.TimeoutExpired:
        print('RESTART_EXTERNAL_TIMEOUT_60S', flush=True)
        return 124


if __name__ == '__main__':
    raise SystemExit(main())

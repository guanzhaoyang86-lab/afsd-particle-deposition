"""Explicit auxiliary-tangent domain bridge; physical maps and gates unchanged."""
from pathlib import Path
import argparse
import ast
import copy
import inspect
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import run_afsd_gpu_20261007 as original
import run_gpu_cached_preflight_20261008 as cache
import afsd_gpu_search_trial_20261009 as model

cpu, core, post, event, batch, old, shared = (original.cpu, original.core, original.post,
    original.event, original.batch, original.old, original.shared)
gpu, solver = original.gpu, original.solver
ROOT, BASE, CPU_BASE = original.ROOT, original.BASE, original.CPU_BASE
require, read, load, rel, verify = original.require, original.read, original.load, original.rel, original.verify
GeometryPathStop, RefinableCandidateError = original.GeometryPathStop, original.RefinableCandidateError
BENCH_PASS, REPLAY_PASS, CONTINUE_PASS = original.BENCH_PASS, original.REPLAY_PASS, original.CONTINUE_PASS
compare_cpu_endpoint, benchmark = original.compare_cpu_endpoint, original.benchmark
remote_audit_visible_support_path = original.remote_audit_visible_support_path
remote_audit_visibility_events = original.remote_audit_visibility_events
CONTRACT = ROOT/'docs/GPU_SEARCH_TENSOR_CONTRACT_20261009.md'
CHECKS = BASE/'search_tensor_components_20261009_v1.json'
TEST = ROOT/'tests/review/test_gpu_search_tensor_20261009.py'
PROTOCOL = 'EXPLICIT_RAW_SYMMETRIC_STRESS_AUXILIARY_DERIVATIVE_DOMAIN_V1'


def adaptation_audit():
    for name in ('scalar_equivalence', 'complete_candidate'):
        left = ast.parse(inspect.getsource(getattr(original, name)))
        right = ast.parse(inspect.getsource(globals()[name]))
        require(ast.dump(left) == ast.dump(right), 'physical candidate gate changed: '+name)
    left = ast.parse(inspect.getsource(original.run_pair)); left.body[0].name = '_run_pair_original'
    right = ast.parse(inspect.getsource(_run_pair_original))
    require(ast.dump(left) == ast.dump(right), 'original complete pair gates changed')
    require(model.ProcessTrial.evaluate is original.model.ProcessTrial.evaluate
        and model.advance_candidate is original.model.advance_candidate, 'physical map or advance changed')
    return dict(status='PASS_SEARCH_TANGENT_DOMAIN_ONLY', inherited_gpu=original.adaptation_audit(),
        search=model.operator_reuse_audit(), complete_candidate_ast_identical=True,
        pair_gates_ast_identical=True, physical_evaluate_and_advance_identical=True)


def frozen_sources(transaction, legacy):
    checks = transaction.read(CHECKS)
    require(checks['status'] == 'PASS_SEARCH_TENSOR_COMPONENTS' and checks['exit_code'] == 0,
        'executed search tensor components required')
    transaction.verify(checks['source_sha256']); transaction.verify(checks['artifacts_sha256'])
    required = [Path(__file__).resolve(), CONTRACT, TEST,
        ROOT/'scripts/afsd_search_tensor_20261009.py', ROOT/'scripts/afsd_gpu_search_trial_20261009.py']
    for path in required:
        require(checks['source_sha256'].get(rel(path)) == transaction.sha(path),
            'component evidence does not bind '+rel(path))
    result = dict(legacy)
    for mapping in (checks['source_sha256'], checks['artifacts_sha256'], {rel(CHECKS):transaction.sha(CHECKS)}):
        for name, digest in mapping.items():
            require(name not in result or result[name] == digest, 'frozen implementation source conflict')
            result[name] = digest
    transaction.verify(result)
    return result


def preflight(args):
    require(args.mode == 'continue' and args.parent, 'explicit complete parent required')
    transaction = cache.CachedPreflight()
    # Reuse every original CPU, GPU benchmark and GPU replay verification.
    base_args = copy.copy(args); base_args.parent = None; base_args.dt = None
    protocol = transaction.preflight(base_args)
    legacy = dict(protocol['source_sha256']); current = frozen_sources(transaction, legacy)
    valid_dir = (ROOT/args.validation).resolve()
    valid = transaction.gpu_receipt(valid_dir, {REPLAY_PASS})
    refs = protocol['reference_sha256']; folder = (ROOT/args.parent).resolve()
    chosen = None; visited = set()
    while True:
        require(folder not in visited, 'cycle in explicit search-domain continuation')
        visited.add(folder)
        meta = transaction.gpu_receipt(folder, {CONTINUE_PASS})
        require(meta['source_sha256'] in (legacy, current)
            and meta['validation_directory'] == rel(valid_dir), 'unregistered implementation or replay in parent chain')
        if meta['source_sha256'] == current:
            require(meta.get('search_tensor_protocol') == PROTOCOL, 'new implementation lacks explicit bridge identity')
        if chosen is None:
            chosen = meta
        transaction.references(folder, meta, refs)
        previous = meta.get('parent_gpu_directory')
        earlier = transaction.gpu_receipt(ROOT/previous, {CONTINUE_PASS}) if previous else valid
        require(meta['parent_sha256'] == earlier['screened_continuation_sha256'], 'explicit bridge parent chain mismatch')
        if not previous:
            break
        folder = (ROOT/previous).resolve()
    dt = args.dt if args.dt is not None else chosen['dt_s']
    require(np.isfinite(dt) and 0 < dt <= .005, 'bounded finite positive macro interval required')
    protocol.update(source_sha256=current, legacy_source_sha256=legacy,
        parent_checkpoint=chosen['screened_continuation_checkpoint'], parent_sha256=chosen['screened_continuation_sha256'],
        parent_gpu_directory=rel((ROOT/args.parent).resolve()), parent_implementation_source_sha256=chosen['source_sha256'],
        dt_s=dt, candidate_forecast_s=chosen['maximum_candidate_wall_s'],
        forecast_s=30.+4.5*chosen['maximum_candidate_wall_s'], search_tensor_protocol=PROTOCOL,
        search_tensor_components=rel(CHECKS), search_tensor_adaptation_audit=adaptation_audit())
    require(protocol['forecast_s'] <= 450., 'original measured cost forecast exceeds450s')
    if args.reference_stage:
        reference_dir = (ROOT/args.reference_stage).resolve()
        reference = transaction.gpu_receipt(reference_dir, {CONTINUE_PASS})
        require(reference['source_sha256'] == legacy and reference['parent_sha256'] == protocol['parent_sha256']
            and reference['dt_s'] == dt, 'replay must match the exact original complete parent and interval')
        transaction.references(reference_dir, reference, refs)
        protocol['search_tensor_reference_stage'] = rel(reference_dir)
    cache.bind_administration(transaction, protocol, cache.DEFAULT_CHECKS)
    transaction.verify(refs)
    protocol['administrative_preflight_statistics'] = transaction.finish()
    return protocol


def run_pair(output, protocol, report, deadline):
    _run_pair_original(output, protocol, report, deadline)
    if protocol.get('search_tensor_reference_stage'):
        reference = ROOT/protocol['search_tensor_reference_stage']; comparisons = {}
        exact = {}
        for role in ('coarse', 'fine_1', 'fine_2'):
            pa, state = load(output/role/'NOT_RESTARTABLE_candidate.npz')
            other, unused = load(reference/role/'NOT_RESTARTABLE_candidate.npz')
            comparisons[role] = compare_cpu_endpoint(dict(pa=pa, state=state), reference/role/'NOT_RESTARTABLE_candidate.npz')
            exact[role] = all(np.array_equal(prop.get_npy_array(), other.properties[name].get_npy_array())
                for name, prop in pa.properties.items())
        report.update(search_tensor_endpoint_comparisons=comparisons,
            search_tensor_particle_arrays_bitwise_equal=exact)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', default='continue', choices=('continue',))
    parser.add_argument('--output', required=True)
    parser.add_argument('--baseline', default=rel(CPU_BASE/'stage_013'))
    parser.add_argument('--cost', default=rel(BASE/'benchmark_v1'))
    parser.add_argument('--validation', default=rel(BASE/'replay_v1'))
    parser.add_argument('--parent'); parser.add_argument('--dt', type=float)
    parser.add_argument('--reference-stage'); parser.add_argument('--worker', action='store_true')
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args(); output = (ROOT/args.output).resolve()
    require(output.is_relative_to(BASE) and output != BASE, 'separate project GPU output required')
    if args.worker:
        require(read(output/'protocol.json')['mode'] == 'continue', 'only registered continuation worker allowed')
        return worker(output)
    require(not output.exists(), 'new output version required')
    protocol = preflight(args)
    if args.preflight_only:
        print(json.dumps(dict(status='PASS_SEARCH_TENSOR_PREFLIGHT', parent=protocol['parent_checkpoint'],
            forecast_s=protocol['forecast_s'], statistics=protocol['administrative_preflight_statistics']))); return 0
    output.mkdir(parents=True); shared.write_new(output/'protocol.json', protocol)
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--mode', 'continue', '--output', rel(output), '--worker']
    env = dict(os.environ, **{k:'1' for k in shared.THREADS}); env['PYTHONPATH'] = ''
    began = time.perf_counter(); code = 124
    with (output/'run.log').open('x', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=600.).returncode
        except subprocess.TimeoutExpired:
            stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    error = None
    try:
        after = cache.CachedPreflight(); after.verify(protocol['source_sha256']); after.verify(protocol['reference_sha256']); after.finish()
    except Exception as exc:
        code = 2; error = repr(exc)
    manifest = output/'manifest.json'
    shared.write_new(output/'supervision.json', dict(exit_code=code, command=command, wall_seconds=time.perf_counter()-began,
        source_inputs_still_match=error is None, source_error=error,
        manifest_sha256=core.sha(manifest) if manifest.exists() else None, log_sha256=core.sha(output/'run.log'),
        accepted=False, restart_allowed=False, search_tensor_protocol=PROTOCOL))
    print(json.dumps(dict(exit_code=code, output=rel(output), status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    return code


# STATIC_ORIGINAL_CANDIDATE_PAIR_AND_WORKER_FUNCTIONS_FOLLOW


def scalar_equivalence(trial, result, directory):
    scalar = model.ProcessTrial(trial.pa, trial.state, result['dt_s']); p = trial.p
    scalar.contact = old.ChainRuleDiscreteContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2, normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'], friction=p['friction'], base_z=-p['substrate_depth'])
    sr = scalar.evaluate(result['displacement'], result['dt_s'])
    differences = batch.candidate_equivalence(dict(sr, energy_comparison=core.energy_diagnostics(scalar, sr)),
        dict(result, energy_comparison=core.energy_diagnostics(trial, result)))
    require(sr['graph_fingerprint'] == result['graph_fingerprint'] and np.array_equal(sr['g'], result['g'])
        and np.array_equal(sr['thermal']['pair_conductance_W_K'], result['thermal']['pair_conductance_W_K']),
        'scalar/batch exact common graph and heat differ')
    core.write(directory/'scalar_batch_equivalence.json', dict(status='PASS_COMPLETE_MAP_EQUIVALENCE', differences=differences))
    return differences


def complete_candidate(pa, state, dt, directory, deadline, per_candidate, frozen):
    """Explicit new-model candidate; no rebinding of any inherited global."""
    began = time.perf_counter(); trial = model.ProcessTrial(pa, state, dt)
    attempt = directory/'broyden'; attempt.mkdir(); guess = batch.guess_for(trial, dt)
    try:
        try:
            result, metrics = solver.solve(trial, dt, guess, attempt, deadline, 24)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            require(time.perf_counter()+1.5*per_candidate+20 < deadline, 'no budget for one unchanged Newton retry')
            attempt = directory/'newton_retry'; attempt.mkdir()
            result, metrics = solver.newton_solve(trial, dt, guess, attempt, deadline, 24)
            metrics['fallback_reason'] = repr(error)
        core.journal(directory, trial, result, 'SOLVED_UNCOMMITTED_candidate'); old.save_graph(directory, trial, result)
        nextpa, nextstate, fresh = model.advance_candidate(trial, result)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        raise RefinableCandidateError(str(error)) from error
    require(np.array_equal(pa.particle_id, nextpa.particle_id) and np.array_equal(pa.m, nextpa.m)
        and model.identity_digest(pa) == model.identity_digest(nextpa), 'immutable material inventory changed')
    target, _ = event.save_candidate(directory, trial, nextpa, nextstate, fresh, frozen)
    old.save_graph(directory, trial, fresh)
    observer = shared.observe_process(pa, nextstate, fresh, trial.motion, trial.p)
    core.write(directory/'observer.json', observer); core.write(directory/'motion.json', trial.motion.to_dict())
    core.write(directory/'transport.json', nextstate['bore_observation'])
    scalar_equivalence(trial, fresh, directory)
    geometry, events = shared.inherited.geometry_audit(trial, fresh)
    core.write(directory/'geometry_audit.json', geometry); core.write(directory/'sampled_events.json', events)
    pose0, pose1 = trial.motion.pose(trial.motion.t0), trial.motion.pose(trial.motion.t1)
    shell = old.audit_contact_shell_finite_enclosure(trial.contact, trial.x0, fresh['position'], pose0, pose1)
    core.write(directory/'continuous_shell_audit.json', shell)
    path = remote_audit_visible_support_path(trial.x0, fresh['position'], fresh['i'], fresh['j'],
                                             pose0, pose1, trial.p, pa.particle_id)
    core.write(directory/'visible_path_audit.json', path)
    if (geometry['transient_contacts_missed_by_endpoints'] or set(geometry['new_non_bore_features'])-{'SHOULDER', 'TOP'}
            or shell['status'] != event.CERTIFIED or path['status'] != old.CERTIFIED):
        try:
            routing = remote_audit_visibility_events(trial.x0, fresh['position'], fresh['i'], fresh['j'],
                                                      pose0, pose1, trial.p, pa.particle_id)
        except Exception as error:
            routing = dict(status='UNRESOLVED_EVENT_DIAGNOSTIC', error=repr(error), requires_resolve=True, acceptance_allowed=False)
        core.write(directory/'visibility_event_routing.json', routing)
        raise GeometryPathStop('geometry/shell/visibility path failed; separate event resolution required')
    core.write(directory/'candidate_qualification.json', dict(status='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN',
        accepted=False, restart_allowed=False, metrics=metrics))
    core.vtk(directory/'fields.vtk', nextpa, trial.p, trial.bulk, nextstate['time_s'], model.CLASSIFICATION)
    return dict(pa=nextpa, state=nextstate, checkpoint=target, observer_path=directory/'observer.json', observer=observer,
                metrics=metrics, wall_seconds=time.perf_counter()-began)


def _run_pair_original(output, protocol, report, deadline):
    pa, state = load(ROOT/protocol['parent_checkpoint']); dt = protocol['dt_s']
    fingerprint = cpu.parent_fingerprint(pa, state)
    program = model.validate_state(state); phase = program.phase(state['time_s'])
    require(state['time_s']+dt <= phase['end']+1e-15, 'GPU interval crosses machine phase')
    if phase['target'][2] > 0 or phase['initial'][2] > 0:
        raise GeometryPathStop('GPU initial integration is bounded to stationary fill; traverse requires separate spreading review')
    items = []
    for role in ('coarse', 'fine_1', 'fine_2'):
        directory = output/role; directory.mkdir()
        pp, ss = (items[1]['pa'], items[1]['state']) if role == 'fine_2' else (pa, state)
        item = complete_candidate(pp, ss, dt if role == 'coarse' else dt/2, directory,
            deadline, protocol['candidate_forecast_s'], protocol['source_sha256'])
        items.append(item)
        require(cpu.parent_fingerprint(pa, state) == fingerprint, 'original parent mutated')
        print(json.dumps(dict(role=role, wall_s=item['wall_seconds'], time_s=item['state']['time_s'])), flush=True)
    coarse, first, last = items
    global_screen = post.compare(last['pa'], coarse['pa'], state['parameters'])
    local_screen = event.local_compare(last['pa'], coarse['pa'], sorted(int(v) for v in pa.particle_id))
    matching = last['state']['bore_observation']['crossing_ids'] == coarse['state']['bore_observation']['crossing_ids']
    screen = dict(global_screen=global_screen, local_screen=local_screen, crossing_ids_match=matching,
        coarse_checkpoint_sha256=core.sha(coarse['checkpoint']),
        fine_checkpoint_sha256=[core.sha(first['checkpoint']), core.sha(last['checkpoint'])],
        passed=bool(global_screen['pass'] and local_screen['passed'] and matching), accepted=False)
    core.write(output/'pair_screen.json', screen)
    require(screen['passed'], 'original coarse/fine screen failed; no GPU continuation')
    comparisons = {}
    if protocol['mode'] == 'replay':
        baseline = ROOT/protocol['cpu_baseline']
        for role, item in zip(('coarse', 'fine_1', 'fine_2'), items):
            comparisons[role] = compare_cpu_endpoint(item, baseline/'macro_0001/attempt_00'/role/'NOT_RESTARTABLE_candidate.npz')
    prefix_n = len(state['machine_history'])
    require(last['state']['machine_history'][:prefix_n] == state['machine_history']
        and last['state']['interval_energy']['records'][:prefix_n] == state['interval_energy']['records'], 'history prefix changed')
    report.update(status=REPLAY_PASS if protocol['mode']=='replay' else CONTINUE_PASS,
        pair_screen=screen, cpu_gpu_endpoint_comparisons=comparisons,
        screened_continuation_checkpoint=rel(last['checkpoint']), screened_continuation_sha256=core.sha(last['checkpoint']),
        screened_end_time_s=last['state']['time_s'], screened_end_step=last['state']['step'],
        parent_time_s=state['time_s'], parent_fingerprint=fingerprint,
        screened_continuation_authorized=True, parent_prefix_unchanged=True,
        maximum_candidate_wall_s=max(v['wall_seconds'] for v in items),
        candidate_wall_s=[v['wall_seconds'] for v in items],
        screened_history=[dict(checkpoint=rel(v['checkpoint']), sha256=core.sha(v['checkpoint']),
            time_s=v['state']['time_s'], observer_path=rel(v['observer_path'])) for v in items[1:]],
        observer=last['observer'], full_program_completed=False)


def worker(output):
    protocol = read(output/'protocol.json'); report = dict(protocol, status='RUNNING_GPU_DIAGNOSTIC', screened_continuation_authorized=False)
    began = time.perf_counter(); code = 2
    def timeout():
        core.write(output/'HARD_TIMEOUT.json', dict(status='TIMEOUT_NOT_ACCEPTED'))
        os._exit(124)
    timer = threading.Timer(540., timeout); timer.daemon=True; timer.start()
    try:
        require(all(os.environ.get(k)=='1' for k in shared.THREADS), 'single CPU thread baseline required')
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
        require(core.sha(ROOT/protocol['parent_checkpoint']) == protocol['parent_sha256'], 'parent changed')
        report['adaptation_audit'] = adaptation_audit(); report['gpu_environment'] = gpu.environment()
        if protocol['mode'] == 'benchmark':
            benchmark(output, protocol, report)
        else:
            run_pair(output, protocol, report, began+525.)
        verify(protocol['source_sha256']); verify(protocol['reference_sha256']); code=0
    except Exception as error:
        report.update(status='FAIL_GPU_STAGE_NOT_ACCEPTED', error=repr(error), traceback=traceback.format_exc(), screened_continuation_authorized=False)
    finally:
        timer.cancel()
        report.update(wall_seconds=time.perf_counter()-began, peak_working_set_bytes=core.memory_bytes(), gpu_statistics=gpu.statistics(),
            accepted=False, restart_allowed=False, accepted_time_advanced_s=0.)
        if code == 0 and protocol['mode'] == 'replay':
            report['cpu_gpu_stage_speedup'] = protocol['cpu_stage_wall_s']/report['wall_seconds']
        report['artifacts_sha256'] = {rel(p):core.sha(p) for p in output.rglob('*') if p.is_file() and p.name not in ('manifest.json','run.log','supervision.json')}
        core.write(output/'manifest.json', report)
    return code


if __name__ == "__main__":
    raise SystemExit(main())

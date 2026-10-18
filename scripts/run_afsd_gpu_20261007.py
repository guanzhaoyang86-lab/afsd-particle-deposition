"""Supervised FP64 GPU comparison and screened diagnostic continuation.

The original CPU history is verified unchanged. Switching compute provider is
an explicit receipt-bound boundary, never an edited historical source hash.
Every continuation still requires a coarse/two-half physical comparison.
"""
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
import run_afsd_deposition_stage_remote_20261007 as cpu
import afsd_gpu_trial_20261007 as model
import afsd_gpu_linalg_20261007 as gpu
import afsd_gpu_equilibrium_20261007 as solver
from test_afsd_gpu_equilibrium_20261007 import audit_gpu_solver_adaptation
from test_afsd_gpu_runner_20261007 import audit_gpu_candidate_adaptation
from test_afsd_gpu_thermal_20261007 import arguments_from_trial_result, compare_fixture
from visible_support_remote_domain_20261007 import audit_visible_support_path as remote_audit_visible_support_path
from visible_support_remote_events_20261007 import audit_visibility_events as remote_audit_visibility_events

core, post, event, batch, old, shared = cpu.core, cpu.post, cpu.event, cpu.batch, cpu.old, cpu.shared
ROOT, require, read, load, rel, verify = cpu.ROOT, cpu.require, cpu.read, cpu.load, cpu.rel, cpu.verify
GeometryPathStop, RefinableCandidateError = cpu.GeometryPathStop, cpu.RefinableCandidateError
BASE = ROOT/'outputs/gpu_acceleration_20261007'
CPU_BASE = ROOT/'outputs/deposition_campaign_20261007'
CONTRACT = ROOT/'docs/GPU_ACCELERATION_CONTRACT_20261007.md'
BENCH_PASS = 'PASS_GPU_MAP_AND_LINEAR_EQUIVALENCE_MEASURED_COST'
REPLAY_PASS = 'PASS_GPU_IDENTICAL_PARENT_CPU_REPLAY_AND_TIME_SCREEN'
CONTINUE_PASS = 'PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION'

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


def adaptation_audit():
    for name in ('scalar_equivalence', 'complete_candidate'):
        reference = ast.parse(inspect.getsource(getattr(cpu, name))).body[0]
        actual = ast.parse(inspect.getsource(globals()[name])).body[0]
        for node in ast.walk(actual):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
                if isinstance(node.func.value, ast.Name) and node.func.value.id == 'solver' and node.func.attr == 'newton_solve':
                    node.func.value.id, node.func.attr = 'core', 'solve'
        require(ast.dump(reference, include_attributes=False) == ast.dump(actual, include_attributes=False),
                'candidate gates or scalar comparison changed')
    return dict(status='PASS_EXPLICIT_GPU_PROVIDERS_ONLY', search=audit_gpu_solver_adaptation(),
                model=model.operator_reuse_audit(), candidate=audit_gpu_candidate_adaptation(), candidate_gates_ast_identical=True)


def gpu_receipt(folder, statuses):
    meta, supervision = read(folder/'manifest.json'), read(folder/'supervision.json')
    require(meta['status'] in statuses and supervision['exit_code'] == 0 and supervision['source_inputs_still_match']
            and supervision['manifest_sha256'] == core.sha(folder/'manifest.json')
            and supervision['log_sha256'] == core.sha(folder/'run.log'), 'invalid GPU receipt')
    verify(meta['source_sha256']); verify(meta['reference_sha256']); verify(meta['artifacts_sha256'])
    require(not meta['accepted'] and not meta['restart_allowed'] and meta['diagnostic_only']
        and meta['classification'] == model.CLASSIFICATION and meta['accepted_time_advanced_s'] == 0., 'diagnostic GPU states only')
    if meta['status'] in {REPLAY_PASS, CONTINUE_PASS}:
        validate_gpu_screened_receipt(meta, folder)
    return meta


def validate_gpu_screened_receipt(meta, folder):
    history, screen = meta['screened_history'], meta['pair_screen']
    require(meta['screened_continuation_authorized'] and len(history) == 2 and screen['passed']
        and screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match'],
        'complete passing GPU fine pair required')
    require(read(folder/'pair_screen.json') == screen, 'GPU manifest and pair screen differ')
    require(screen['fine_checkpoint_sha256'] == [r['sha256'] for r in history]
        and meta['screened_continuation_checkpoint'] == history[-1]['checkpoint']
        and meta['screened_continuation_sha256'] == history[-1]['sha256'], 'GPU continuation must bind exact last fine state')
    parent_pa, parent_state = load(ROOT/meta['parent_checkpoint'])
    require(core.sha(ROOT/meta['parent_checkpoint']) == meta['parent_sha256']
        and cpu.parent_fingerprint(parent_pa, parent_state) == meta['parent_fingerprint'], 'GPU full parent identity changed')
    prefix_n=len(parent_state['machine_history'])
    for half, (role, row) in enumerate(zip(('fine_1', 'fine_2'), history), 1):
        require((ROOT/row['checkpoint']).resolve() == (folder/role/'NOT_RESTARTABLE_candidate.npz').resolve()
            and core.sha(ROOT/row['checkpoint']) == row['sha256'], 'GPU fine checkpoint path/hash mismatch')
        _, state = load(ROOT/row['checkpoint']); model.validate_state(state)
        require(state['diagnostic_only'] and not state['accepted'] and not state['restart_allowed']
            and state['time_s'] == row['time_s'], 'GPU fine metadata mismatch')
        require(abs(state['time_s']-meta['parent_time_s']-half*meta['dt_s']/2) <= 1e-13
            and state['step'] == parent_state['step']+half, 'GPU half-step time/history discontinuity')
        require(state['machine_history'][:prefix_n] == parent_state['machine_history']
            and state['interval_energy']['records'][:prefix_n] == parent_state['interval_energy']['records'], 'GPU parent history prefix differs')
    require(state['time_s'] == meta['screened_end_time_s'] and state['step'] == meta['screened_end_step']
        and abs(state['time_s']-meta['parent_time_s']-meta['dt_s']) <= 1e-13, 'GPU interval endpoint mismatch')
    require(core.sha(folder/'coarse/NOT_RESTARTABLE_candidate.npz') == screen['coarse_checkpoint_sha256'], 'GPU coarse screen binding mismatch')


def freeze(cpu_meta):
    mapping = dict(cpu_meta['source_sha256'])
    components_path = BASE/'components_v1/manifest.json'
    components = read(components_path)
    require(components['status']=='PASS_GPU_COMPONENTS' and all(r['exit_code']==0 for r in components['checks']), 'executed GPU components required')
    verify(components['source_sha256']); verify(components['artifacts_sha256'])
    mapping.update(components['source_sha256']); mapping.update(components['artifacts_sha256'])
    files = [CONTRACT, Path(__file__).resolve(), components_path]
    files += list((ROOT/'scripts').glob('*afsd_gpu*20261007.py'))
    files += list((gpu.OVERLAY).rglob('*.pyd')) + list(gpu.OVERLAY.rglob('METADATA'))
    files += list(gpu.OVERLAY.rglob('*.py'))
    for path in files:
        mapping[rel(path)] = core.sha(path)
    verify(mapping)
    require(read(BASE/'components_v1/linalg.json')['status'] == 'PASS_CUDA_FP64_COMPONENTS', 'GPU linalg checks required')
    require(read(BASE/'components_v1/thermal.json')['status'] == 'PASS'
            and read(BASE/'components_v1/thermal.json')['gpu_requested'], 'executed GPU thermal checks required')
    adaptation_audit()
    return mapping


def references(folder, meta, mapping):
    mapping.update(meta['artifacts_sha256'])
    for name in ('manifest.json', 'supervision.json', 'run.log'):
        mapping[rel(folder/name)] = core.sha(folder/name)


def preflight(args):
    baseline = (ROOT/args.baseline).resolve()
    require(baseline == (CPU_BASE/'stage_013').resolve(), 'this registered replay is specifically CPU stage013')
    original = cpu.checked_receipt(baseline, {cpu.STAGE_PASS, cpu.PREFIX_PASS, cpu.SHOULDER_PASS})
    require(len(original['pair_screens']) == 1 and original['pair_screens'][0]['refinement_level'] == 0
        and original['screened_continuation_sha256'] == 'c3a9e4e5b194de04a5528b06d7fce6815c97a4f8b9d823aa53206dc773ff73e2',
        'registered exact CPU one-pair reference changed')
    cost_folder = ROOT/original['cost_directory']
    # Full original history, pair bindings, persistent identity and source hashes.
    _, cpu_target, latest_pa, latest_state, refs, _ = cpu.continuation_parent(cost_folder, baseline, original['source_sha256'])
    references(baseline, original, refs)
    frozen = freeze(original)
    dt = original['pair_screens'][0]['dt_s']
    replay_parent = ROOT/original['parent_checkpoint']
    protocol = dict(schema='EXPLICIT_GPU_FP64_BRIDGE_20261007', mode=args.mode,
        source_sha256=frozen, reference_sha256=refs, classification=model.CLASSIFICATION,
        cpu_baseline=rel(baseline), cpu_stage_wall_s=original['wall_seconds'],
        cpu_candidate_wall_s=original['maximum_candidate_wall_s'], dt_s=dt,
        parent_checkpoint=rel(replay_parent), parent_sha256=core.sha(replay_parent),
        cpu_latest_checkpoint=rel(cpu_target), cpu_latest_sha256=core.sha(cpu_target),
        particle_count=len(latest_pa.x), free_dof=int(latest_state['free'].sum()),
        accepted=False, restart_allowed=False, diagnostic_only=True, accepted_time_advanced_s=0.,
        internal_limit_s=540., external_limit_s=600., forecast_limit_s=450.,
        compute_protocol='CPU_PYSPH_WITH_EXPLICIT_CUDA_FP64_DENSE_SOLVES_V1')
    if args.mode == 'benchmark':
        protocol['forecast_s'] = 90.
    else:
        require(args.cost, '--cost required')
        cost_dir = (ROOT/args.cost).resolve(); cost = gpu_receipt(cost_dir, {BENCH_PASS})
        require(cost['source_sha256'] == frozen and cost['cpu_baseline'] == rel(baseline), 'GPU benchmark binding changed')
        references(cost_dir, cost, refs)
        protocol['cost_directory'] = rel(cost_dir)
        if args.mode == 'replay':
            # Predeclared conservative CPU measured stage, unchanged search budgets.
            protocol['forecast_s'] = 30.+1.5*original['wall_seconds']
            protocol['candidate_forecast_s'] = original['maximum_candidate_wall_s']
        else:
            require(args.validation, '--validation required')
            valid_dir = (ROOT/args.validation).resolve(); valid = gpu_receipt(valid_dir, {REPLAY_PASS})
            require(valid['source_sha256'] == frozen and valid['cpu_baseline'] == rel(baseline), 'GPU replay binding changed')
            references(valid_dir, valid, refs)
            chosen = valid
            if args.parent:
                parent_dir = (ROOT/args.parent).resolve(); visited = set()
                while True:
                    require(parent_dir not in visited, 'GPU continuation cycle'); visited.add(parent_dir)
                    prior = gpu_receipt(parent_dir, {CONTINUE_PASS})
                    require(prior['source_sha256'] == frozen and prior['validation_directory'] == rel(valid_dir), 'mixed GPU continuation')
                    references(parent_dir, prior, refs)
                    if chosen is valid:
                        chosen = prior
                    previous = prior.get('parent_gpu_directory')
                    earlier = gpu_receipt(ROOT/previous, {CONTINUE_PASS}) if previous else valid
                    require(prior['parent_sha256'] == earlier['screened_continuation_sha256'], 'GPU parent chain mismatch')
                    if not previous:
                        break
                    parent_dir = ROOT/previous
            protocol.update(parent_checkpoint=chosen['screened_continuation_checkpoint'], parent_sha256=chosen['screened_continuation_sha256'],
                validation_directory=rel(valid_dir), parent_gpu_directory=rel(ROOT/args.parent) if args.parent else None,
                dt_s=args.dt if args.dt is not None else chosen['dt_s'],
                candidate_forecast_s=chosen['maximum_candidate_wall_s'])
            require(0 < protocol['dt_s'] <= .005, 'bounded positive GPU macro interval required')
            protocol['forecast_s'] = 30.+4.5*protocol['candidate_forecast_s']
    require(protocol['forecast_s'] <= 450., 'GPU stage measured forecast exceeds450s')
    verify(refs)
    return protocol


def timed(function, *args):
    began = time.perf_counter(); value = function(*args)
    return value, time.perf_counter()-began


def benchmark(output, protocol, report):
    pa, state = load(ROOT/protocol['parent_checkpoint']); dt = protocol['dt_s']
    cpu_trial, gpu_trial = cpu.model.ProcessTrial(pa, state, dt), model.ProcessTrial(pa, state, dt)
    guess = batch.guess_for(cpu_trial, dt)
    gpu.solve(np.eye(8), np.ones(8)); gpu.statistics(reset=True)
    cpu_rows, gpu_rows, diffs = [], [], []
    for unused in range(3):
        cr, ct = timed(cpu_trial.evaluate, guess, dt); gr, gt = timed(gpu_trial.evaluate, guess, dt)
        diffs.append(batch.candidate_equivalence(dict(cr, energy_comparison=core.energy_diagnostics(cpu_trial, cr)),
            dict(gr, energy_comparison=core.energy_diagnostics(gpu_trial, gr))))
        for key in ('i', 'j', 'g'):
            require(np.array_equal(cr[key], gr[key]), 'CPU/GPU support changed')
        require(cr['graph_fingerprint'] == gr['graph_fingerprint'], 'CPU/GPU graph identity changed')
        cpu_rows.append(ct); gpu_rows.append(gt)
    thermal = compare_fixture(arguments_from_trial_result(cpu_trial, cr))
    tangent, tangent_wall = timed(cpu_trial.tangent, cr)
    free = state['free'].ravel(); matrix = tangent[np.ix_(free, free)]; rhs = -cr['force'][state['free']]
    rows = []
    for kind, cpu_fn, gpu_fn, inputs in (('solve', np.linalg.solve, gpu.solve, (matrix, rhs)), ('inv', np.linalg.inv, gpu.inv, (matrix,))):
        for repeat in range(3):
            ca, ct = timed(cpu_fn, *inputs); ga, gt = timed(gpu_fn, *inputs)
            require(np.isfinite(ga).all(), 'nonfinite GPU linear result')
            error = float(np.linalg.norm(ga-ca, np.inf)/max(np.linalg.norm(ca, np.inf), 1e-30))
            require(error <= 1e-8, 'CPU/GPU actual mechanical linear mismatch')
            if kind == 'solve':
                residual = float(np.linalg.norm(matrix@ga-rhs, np.inf)/max(1., np.linalg.norm(rhs, np.inf)))
                require(residual <= 1e-9, 'original mechanical linear residual failed')
            else:
                residual = float(np.linalg.norm(matrix@ga-np.eye(len(matrix)), np.inf))
                require(residual <= 1e-8, 'GPU inverse posterior residual failed')
                require(float(np.linalg.norm(matrix, np.inf)*np.linalg.norm(ga, np.inf)) <= 1e12, 'original condition gate failed')
            rows.append(dict(operation=kind, repeat=repeat, cpu_wall_s=ct, gpu_wall_s=gt, cpu_gpu_relative_error=error, residual=residual))
    report.update(status=BENCH_PASS, cpu_map_wall_s=cpu_rows, gpu_map_wall_s=gpu_rows,
        map_speedup=float(np.median(cpu_rows)/np.median(gpu_rows)), full_map_differences=diffs,
        thermal_component=thermal, tangent_cpu_wall_s=tangent_wall, linear_benchmarks=rows,
        linear_speedups={kind:float(np.median([r['cpu_wall_s'] for r in rows if r['operation']==kind])/
            np.median([r['gpu_wall_s'] for r in rows if r['operation']==kind])) for kind in ('solve', 'inv')},
        force_root_cost_measured=False, stable_timestep_s='NOT_ESTABLISHED_BY_UNSOLVED_MAP')


def compare_cpu_endpoint(item, path):
    other, state = load(path)
    require(np.array_equal(item['pa'].particle_id, other.particle_id) and np.array_equal(item['pa'].m, other.m), 'CPU/GPU inventory mismatch')
    require(item['state']['machine_history'] == state['machine_history'], 'CPU/GPU prescribed history mismatch')
    global_diff = post.compare(item['pa'], other, item['state']['parameters'])
    require(global_diff['position_max_m'] <= .001*item['state']['parameters']['dx']
        and global_diff['temperature_max_K'] <= 1e-4 and global_diff['stress_relative'] <= 1e-3
        and global_diff['plastic_relative'] <= 1e-3, 'predeclared stricter CPU/GPU endpoint comparison failed')
    require(item['state']['bore_observation']['crossing_ids'] == state['bore_observation']['crossing_ids'], 'CPU/GPU crossing identity mismatch')
    return global_diff


def run_pair(output, protocol, report, deadline):
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


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', required=True, choices=('benchmark','replay','continue')); p.add_argument('--output',required=True)
    p.add_argument('--baseline', default=rel(CPU_BASE/'stage_013')); p.add_argument('--cost'); p.add_argument('--validation')
    p.add_argument('--parent'); p.add_argument('--dt',type=float); p.add_argument('--worker',action='store_true')
    args=p.parse_args(); output=(ROOT/args.output).resolve()
    require(output.is_relative_to(BASE) and output != BASE, 'separate GPU output directory required')
    if args.worker:
        return worker(output)
    require(not output.exists(), 'new output version required')
    protocol=preflight(args); output.mkdir(parents=True); shared.write_new(output/'protocol.json',protocol)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--mode',args.mode,'--output',rel(output),'--worker']
    env=dict(os.environ, **{k:'1' for k in shared.THREADS}); env['PYTHONPATH']=''
    began=time.perf_counter(); code=124
    with (output/'run.log').open('x',encoding='utf-8') as stream:
        try:
            code=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600.).returncode
        except subprocess.TimeoutExpired:
            stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    error=None
    try:
        verify(protocol['source_sha256']); verify(protocol['reference_sha256'])
    except Exception as exc:
        code=2; error=repr(exc)
    manifest=output/'manifest.json'
    shared.write_new(output/'supervision.json',dict(exit_code=code,command=command,wall_seconds=time.perf_counter()-began,
        source_inputs_still_match=error is None, source_error=error, manifest_sha256=core.sha(manifest) if manifest.exists() else None,
        log_sha256=core.sha(output/'run.log'),accepted=False,restart_allowed=False))
    print(json.dumps(dict(exit_code=code,output=rel(output),status=read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    return code


if __name__=='__main__':
    raise SystemExit(main())

"""Supervised independent QS cost gate and one bounded, auditable root.

Use --stage cost first, then --stage root. Never overwrites existing evidence.
Physical time is independent of the old dynamic history, even if accepted.
"""
from pathlib import Path
import argparse
import ctypes
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.quasistatic_process import (
    CLASSIFICATION, QuasistaticTrial, independent_initial_state, pose,
    newton_direction, energy_diagnostics, accepted_copy)
from afsd_backend.pysph.mechanics import positions, operators
from afsd_backend.pysph.quasistatic_tangent import (
    gradient_map, material_tangent, assemble_frozen_force_tangent)
from afsd_backend.pysph.output import checkpoint, restore, vtk
from pysph.solver.utils import load

BASE = ROOT/'outputs/hot_quasistatic_probe_20261004'
LOGS = ROOT/'logs/hot_quasistatic_probe_20261004'
DT = .0005


def clean(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    return value


def write(path, value):
    pending = path.with_suffix('.pending.json')
    pending.write_text(json.dumps(clean(value), indent=2, allow_nan=False), encoding='utf-8')
    pending.replace(path)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def memory_bytes():
    class Counters(ctypes.Structure):
        _fields_ = [('cb', ctypes.c_ulong), ('PageFaultCount', ctypes.c_ulong)] + [
            (name, ctypes.c_size_t) for name in ('PeakWorkingSetSize', 'WorkingSetSize',
                'QuotaPeakPagedPoolUsage', 'QuotaPagedPoolUsage', 'QuotaPeakNonPagedPoolUsage',
                'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage')]
    counter = Counters(); counter.cb = ctypes.sizeof(counter)
    kernel = ctypes.windll.kernel32
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    fn = ctypes.windll.psapi.GetProcessMemoryInfo
    fn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_ulong]
    if not fn(kernel.GetCurrentProcess(), ctypes.byref(counter), counter.cb):
        raise OSError('working set measurement failed')
    return int(counter.PeakWorkingSetSize)


def sources():
    files = list((ROOT/'src/afsd_backend/pysph').glob('*.py'))
    files += [Path(__file__).resolve(), ROOT/'scripts/test_hot_quasistatic_probe.py',
              ROOT/'docs/HOT_QUASISTATIC_PROCESS_CONTRACT_20261004.md']
    files += list((ROOT/'tests/review').glob('test_quasistatic*.py'))
    return {str(path.relative_to(ROOT)).replace('\\', '/'): sha(path) for path in sorted(files)}


def parameters():
    path = ROOT/'outputs/contact_material_search_20261004/root/resolved_config.json'
    p = json.loads(path.read_text(encoding='utf-8'))['parameters']
    p['bore_width'] = .002504
    return p, dict(parameter_source=str(path.relative_to(ROOT)), source_sha256=sha(path),
        override=dict(bore_width_m=.002504), old_state_read=False,
        stress_free_initial_gap_per_side_m=.000002)


def summary(r):
    return {k: r[k] for k in ('force_residual_N', 'force_limit_N', 'scaled_residual',
        'density_residual', 'thermal_residual_K', 'virtual_work_residual_J',
        'material_endpoint_work_J', 'density_iterations', 'thermal_iterations', 'wall_seconds')}


def verify_snapshot(path, pa, state):
    data = load(str(path)); other = data['arrays']['solid']; saved = data['solver_data']
    for key, prop in pa.properties.items():
        if key not in other.properties or not np.array_equal(prop.get_npy_array(), other.properties[key].get_npy_array()):
            raise ValueError('checkpoint particle property mismatch: '+key)
    for key in ('time_s', 'step', 'classification', 'crossing_ids'):
        if saved[key] != state[key]: raise ValueError('checkpoint state mismatch: '+key)
    for key in ('contact_history', 'free'):
        if not np.array_equal(saved[key], state[key]): raise ValueError('checkpoint history mismatch: '+key)


def journal(directory, trial, r, name):
    np.savez_compressed(directory/(name+'.npz'), accepted=False, restart_allowed=False,
        particle_id=trial.pa.particle_id, position_m=r['position'],
        displacement_m=r['displacement'], velocity_m_s=r['displacement']/r['dt_s'],
        density_kg_m3=r['density'], stress_Pa=r['stress'], temperature_K=r['temperature'],
        plastic_strain=r['ep'], force_N=r['force'], internal_force_N=r['internal_force'],
        contact_force_N=r['contact']['force_N'], contact_history=r['contact']['history'])
    write(directory/(name+'.json'), dict(accepted=False, restart_allowed=False,
        metrics=summary(r), contact=r['contact'], energy=energy_diagnostics(trial, r),
        thermal=r['thermal']['ledger']))


def initial_audit(trial):
    pa, p = trial.pa, trial.p
    x = positions(pa)
    contact = trial.contact.evaluate(x, x, pose(p, 0.), pose(p, 0.), 1.,
                                     trial.state['contact_history'])
    force = float(np.max(np.abs(contact['force_N'])))
    potential = float(np.sum(contact['normal_potential1_J']))
    if force > 1e-10 or potential > 1e-20 or np.any(pa.S != 0) or np.any(pa.ep != 0):
        raise ValueError('independent unloaded initial equilibrium not verified')
    i, j, g, _, _ = operators(pa)
    from afsd_backend.pysph.quasistatic_process import _geometry
    H, B = _geometry(x, np.zeros_like(x), pa.m, pa.rho, i, j, g)
    V = pa.m/pa.rho
    D = gradient_map(V, i, j, g, B)
    C = material_tangent(trial.s0, pa.ep, H, pa.temperature, DT, p['shear'], trial.law)
    J = assemble_frozen_force_tangent(D, V, C['tangent_Pa_s'], pa.rho,
                                      p['rho0'], trial.bulk, DT, np.zeros_like(x))
    free = trial.state['free'].ravel()
    K = -J[np.ix_(free, free)]
    asym = float(np.max(np.abs(K-K.T)))
    ev = np.linalg.eigvalsh(K)
    rank_ratio = float(ev[0]/ev[-1])
    if ev[0] <= 1e-12*ev[-1]:
        raise ValueError('initial free elastic stiffness has unresolved null/negative modes')
    return dict(status='PASS_INDEPENDENT_MECHANICAL_EQUILIBRIUM_NOT_THERMAL_EQUILIBRIUM',
        maximum_force_N=force, normal_contact_energy_J=potential,
        elastic_stiffness_min_eigenvalue_N_m=float(ev[0]),
        elastic_stiffness_max_eigenvalue_N_m=float(ev[-1]), relative_min_eigenvalue=rank_ratio,
        stiffness_max_asymmetry_N_m=asym, free_dofs=int(free.sum()),
        fixed_dofs=int((~free).sum()), particle_count=len(x), mass_kg=float(pa.m.sum()),
        initially_below_shoulder_bar_particles=int(np.count_nonzero((pa.component_id == 1)&(pa.z < p['layer']))),
        tool_pose=pose(p, 0.), temperature_K=float(pa.temperature[0]),
        elapsed_model_time_s=0., reused_dynamic_checkpoint=False)


def test_gate(current):
    data = json.loads((LOGS/'tests.json').read_text(encoding='utf-8'))
    if data['exit_code'] != 0 or data['source_sha256'] != current:
        raise ValueError('missing matching passed numerical/component tests')


def worker(stage):
    directory = BASE/stage
    directory.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    def timeout():
        write(directory/'HARD_TIMEOUT.json', dict(accepted=False, stage=stage))
        os._exit(124)
    timer = threading.Timer(540., timeout); timer.daemon = True; timer.start()
    published = False
    report = dict(classification=CLASSIFICATION, stage=stage, accepted_steps=0,
        accepted_time_s=0., production_authorized=False, temporal_accuracy='NOT_VALIDATED',
        original_dynamic_time_s=1.3503936329109039e-6,
        backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5',
        dt_s=DT, mass_time_modulus_scaling=False, omitted_inertia=True,
        stable_timestep_s='NOT_ESTABLISHED', requested_physical_duration_s=DT,
        particle_count=342)
    try:
        frozen = sources(); test_gate(frozen)
        p, provenance = parameters()
        report.update(source_sha256=frozen, provenance=provenance)
        write(directory/'manifest.json', dict(report, status='RUNNING_UNCOMMITTED'))
        write(directory/'resolved_config.json', dict(parameters=p, classification=CLASSIFICATION,
            provenance=provenance, dt_s=DT, boundary='BOTTOM_LAYER_CLAMP_TOP_BAR_Z_ONLY'))
        write(directory/'accepted_history.json', [])
        pa, state = independent_initial_state(p)
        trial = QuasistaticTrial(pa, state, p)
        if stage == 'cost':
            begin = time.perf_counter()
            audit = initial_audit(trial)
            audit['wall_seconds'] = time.perf_counter()-begin
            write(directory/'initial_audit.json', audit)
            checkpoint(directory/'independent_initial.npz', pa, state)
            vtk(directory/'independent_initial.vtk', pa, p, trial.bulk, 0., CLASSIFICATION)
            initial_hash = sha(directory/'independent_initial.npz')
            restored, restored_state = restore(directory/'independent_initial.npz')
            if not np.array_equal(restored.particle_id, pa.particle_id) or restored_state['time_s'] != 0.:
                raise ValueError('initial checkpoint readback failed')
            begin = time.perf_counter()
            r = trial.evaluate(trial.prescribed(DT), DT)
            eval_time = time.perf_counter()-begin
            journal(directory, trial, r, 'UNCOMMITTED_cost_candidate')
            begin = time.perf_counter()
            delta, diagnostics = newton_direction(trial, r)
            tangent_time = time.perf_counter()-begin
            np.savez_compressed(directory/'search_direction.npz', displacement_correction_m=delta,
                                accepted=False, restart_allowed=False)
            # All eight permitted line-search residuals per Newton iteration.
            iteration_cost = tangent_time+8*eval_time
            iterations = min(24, int((450.-30.)/(1.5*iteration_cost)))
            if iterations < 2:
                status = 'NOT_RUN_COST_NO_USEFUL_ROOT_BUDGET'
            else:
                status = 'PASS_COST_FOR_BOUNDED_ROOT_NOT_PROCESS_VALIDATION'
            report.update(status=status, initial_audit=audit,
                initial_snapshot_sha256=initial_hash, candidate=summary(r),
                complete_trial_wall_s=eval_time, tangent_and_linear_solve_wall_s=tangent_time,
                iteration_forecast_s=iteration_cost, maximum_iterations=iterations,
                root_forecast_wall_s=1.5*iterations*iteration_cost+30., linear=diagnostics,
                predicted_0_118s_wall_s=(.118/DT)*(1.5*iteration_cost+0.),
                transport_prediction='COST_ONLY_NOT_EVIDENCE_OF_CROSSING')
        elif stage == 'root':
            cost = json.loads((BASE/'cost/manifest.json').read_text(encoding='utf-8'))
            if cost['status'] != 'PASS_COST_FOR_BOUNDED_ROOT_NOT_PROCESS_VALIDATION' or cost['source_sha256'] != frozen:
                raise ValueError('no matching feasible full candidate cost gate')
            cost_parameters = json.loads((BASE/'cost/resolved_config.json').read_text(encoding='utf-8'))['parameters']
            if cost['provenance'] != provenance or cost_parameters != p:
                raise ValueError('cost parameter provenance or resolved values changed')
            maximum = cost['maximum_iterations']
            r = trial.evaluate(trial.prescribed(DT), DT)
            trace = []; evaluations = 1
            for iteration in range(maximum+1):
                entry = dict(iteration=iteration, evaluations=evaluations, **summary(r))
                trace.append(entry)
                write(directory/'iterations.json', trace)
                journal(directory, trial, r, 'latest_UNCOMMITTED_candidate')
                print(json.dumps(clean(entry)), flush=True)
                if r['scaled_residual'] <= 1.:
                    accepted, next_state = accepted_copy(trial, r)
                    if sources() != frozen: raise ValueError('source drift before commit')
                    next_state.update(source_sha256=frozen, parameter_provenance=provenance, parameters=p)
                    staged_state = dict(next_state, restart_allowed=False)
                    checkpoint(directory/'UNPUBLISHED_candidate.npz', accepted, staged_state)
                    verify_snapshot(directory/'UNPUBLISHED_candidate.npz', accepted, staged_state)
                    vtk(directory/'accepted.vtk', accepted, p, trial.bulk, DT, CLASSIFICATION)
                    write(directory/'energy.json', dict(mechanical=next_state['energy'], thermal=r['thermal']['ledger']))
                    write(directory/'tool_pose.json', pose(p, DT))
                    write(directory/'transport.json', dict(crossing_ids=next_state['crossing_ids'],
                        crossing_mass_kg=sum(float(m) for ident, m in zip(accepted.particle_id, accepted.m)
                            if int(ident) in next_state['crossing_ids']),
                        max_particle_displacement_m=float(np.linalg.norm(r['displacement'], axis=1).max()),
                        max_new_plastic_strain=float(np.max(accepted.ep-pa.ep))))
                    if sources() != frozen or parameters() != (p, provenance):
                        raise ValueError('source/input drift before final publication')
                    checkpoint(directory/'accepted.npz', accepted, next_state)
                    published = True
                    report.update(accepted_time_s=DT, accepted_steps=1)
                    write(directory/'accepted_history.json', [dict(time_s=DT, dt_s=DT,
                        snapshot='accepted.npz', sha256=sha(directory/'accepted.npz'))])
                    report.update(status='ACCEPTED_ONE_QS_STEP_NOT_PROCESS_VALIDATION',
                        accepted_time_s=DT, accepted_steps=1, evaluations=evaluations)
                    break
                if iteration == maximum: raise ValueError('bounded Newton iteration budget exhausted')
                if time.perf_counter()-start+1.5*cost['iteration_forecast_s'] > 525.:
                    raise ValueError('remaining cost budget cannot support next iteration')
                delta, linear = newton_direction(trial, r)
                merit = float(np.linalg.norm(r['force'][state['free']]))
                accepted_search = False
                for backtrack in range(8):
                    alpha = .5**backtrack
                    try:
                        candidate = trial.evaluate(r['displacement']+alpha*delta, DT)
                        evaluations += 1
                        new_merit = float(np.linalg.norm(candidate['force'][state['free']]))
                        entry.setdefault('line_search', []).append(dict(alpha=alpha, merit_N=new_merit))
                        if new_merit < (1-1e-4*alpha)*merit:
                            r = candidate; accepted_search = True; break
                    except (ValueError, ArithmeticError) as error:
                        evaluations += 1
                        entry.setdefault('line_search', []).append(dict(alpha=alpha, rejected=str(error)))
                    if evaluations >= 100: break
                entry['linear'] = linear
                write(directory/'iterations.json', trace)
                if not accepted_search: raise ValueError('complete-residual line search failed')
                if evaluations >= 100: raise ValueError('complete residual evaluation budget exhausted')
            report['evaluations'] = evaluations
        if sources() != frozen: raise ValueError('source changed during bounded probe')
        code = 0
    except Exception as error:
        traceback.print_exc()
        report.update(status=('FAILED_AFTER_ACCEPTED_CHECKPOINT_PUBLICATION' if published
                              else 'REJECTED_NO_NEW_ACCEPTED_STATE'), error=repr(error))
        code = 2
    finally:
        timer.cancel()
    report['wall_seconds'] = time.perf_counter()-start
    report['peak_process_working_set_bytes'] = memory_bytes()
    write(directory/'manifest.json', report)
    print(json.dumps(clean(report)), flush=True)
    return code


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['cost', 'root'], required=True)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.worker: return worker(args.stage)
    LOGS.mkdir(parents=True, exist_ok=True)
    log = LOGS/(args.stage+'.log')
    if log.exists() or (BASE/args.stage).exists():
        raise ValueError('refuse to overwrite evidence')
    cmd = [sys.executable, '-B', str(Path(__file__).resolve()), '--stage', args.stage, '--worker']
    start = time.perf_counter()
    with log.open('w', encoding='utf-8') as stream:
        try:
            result = subprocess.run(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=600.)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = 124
    write(LOGS/(args.stage+'_supervision.json'), dict(command=cmd, exit_code=code,
        external_wall_s=time.perf_counter()-start, external_limit_s=600., log=str(log)))
    print(json.dumps(dict(stage=args.stage, exit_code=code, log=str(log))))
    return code


if __name__ == '__main__':
    raise SystemExit(main())

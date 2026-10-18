"""Read-only inertia/macro-step assessment; NO solver advance or checkpoint write."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import math
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.process_residual import midpoint_residual

CASE = ROOT/'outputs/coupled_afsd_20261002/short_07us'
FREEZE = ROOT/'docs/COUPLED_AFSD_FREEZE.json'
SOURCES = ['scripts/assess_process_route.py', 'src/afsd_backend/pysph/process_residual.py',
           'docs/PROCESS_ROUTE_CONTRACT.md']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def require(condition, message):
    if not condition: raise ValueError(message)


def vec(value, n, label):
    value = np.asarray(value, float)
    require(value.shape == (n, 3) and np.isfinite(value).all(), 'invalid '+label)
    return value


def checked_inputs(case):
    case = Path(case).resolve()
    require(case.resolve() == CASE.resolve(), 'only frozen actual short-window case is registered')
    frozen = json.loads(FREEZE.read_text(encoding='utf-8'))['files_sha256']
    manifest_path = case/'manifest.json'
    require(sha(manifest_path) == frozen[manifest_path.relative_to(ROOT).as_posix()], 'manifest freeze drift')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    require(manifest['accepted_steps'] == 67 and manifest['accepted_time_s'] == 7e-7,
            'not the registered complete trajectory')
    require(manifest['status'] == 'COMPLETE_REGISTERED_INTERVAL_NOT_SUSTAINED_DEPOSITION'
            and manifest['production_authorized'] is False, 'wrong source status')
    inputs = {FREEZE.relative_to(ROOT).as_posix(): sha(FREEZE),
              manifest_path.relative_to(ROOT).as_posix(): sha(manifest_path)}
    for rel, digest in manifest['source_sha256'].items():
        # Upstream files are pinned by git and by the already frozen run manifest,
        # not duplicated in the project-file index. Never skip their actual hash.
        directly_indexed = rel in frozen and frozen[rel] == digest
        upstream_transitive = rel not in frozen and rel.startswith('third_party/PySPH/')
        require((ROOT/rel).resolve().is_relative_to(ROOT) and
                (directly_indexed or upstream_transitive) and sha(ROOT/rel) == digest, 'source drift '+rel)
        inputs[rel] = digest
    for name, digest in manifest['artifact_sha256'].items():
        require(Path(name).name == name, 'unsafe artifact')
        rel = (case/name).relative_to(ROOT).as_posix()
        require(frozen.get(rel) == digest and sha(case/name) == digest, 'artifact drift '+name)
        inputs[rel] = digest
    return manifest, inputs


def reconstruct_impulses(step):
    """Every trace and every subinterval, including zero-contact particles."""
    ids = np.asarray(step['particle_id']); n = len(ids)
    require(len(np.unique(ids)) == n, 'duplicate IDs')
    normal = np.zeros((n, 3)); tangent = normal.copy(); reaction = normal.copy()
    free = normal.copy(); prescribed = normal.copy(); seen = set(); intervals = 0
    eps = np.finfo(float).eps
    for trace in step['traces']:
        i = trace['particle_index']
        require(0 <= i < n and i not in seen and ids[i] == trace['particle_id'], 'trace ID/index mismatch')
        seen.add(i); free[i] = trace['free_axes']; prescribed[i] = trace['prescribed_velocity']
        result = trace['result']; count_faces = len(trace['contact_rows']); dt_sum = 0.
        face_totals = {key: np.zeros((count_faces, 3)) for key in ('normal_impulse', 'tangential_impulse')}
        for interval in result['intervals']:
            dt = float(interval['dt_s']); require(dt > 0 and math.isfinite(dt), 'invalid subinterval dt')
            dt_sum += dt; intervals += 1
            for key, destination in [('normal_impulse', normal), ('tangential_impulse', tangent)]:
                values = np.asarray(interval[key], float).reshape(-1, 3)
                require(values.shape == (count_faces, 3) and np.isfinite(values).all(), 'face impulse shape')
                face_totals[key] += values
                destination[i] += values.sum(axis=0)
            value = np.asarray(interval['constraint_impulse'], float)
            require(value.shape == (3,) and np.isfinite(value).all(), 'constraint impulse')
            reaction[i] += value
        require(abs(dt_sum-step['dt_s']) <= 128*eps*step['dt_s'], 'subinterval coverage')
        for key, value in [('normal_impulse', normal[i]), ('tangential_impulse', tangent[i])]:
            stored_faces = np.asarray(result[key], float).reshape(-1, 3)
            require(stored_faces.shape == (count_faces, 3) and np.isfinite(stored_faces).all(),
                    'stored face impulse shape')
            face_bound = 1e-17+128*eps*np.maximum(np.linalg.norm(stored_faces, axis=1),
                                               np.linalg.norm(face_totals[key], axis=1))
            require(np.all(np.linalg.norm(stored_faces-face_totals[key], axis=1) <= face_bound),
                    'stored per-face contact impulse mismatch')
            total = stored_faces.sum(axis=0)
            bound = 1e-17+128*eps*max(np.linalg.norm(total), np.linalg.norm(value))
            require(np.linalg.norm(total-value) <= bound, 'stored contact impulse mismatch')
        total = np.asarray(result['constraint_impulse'], float)
        require(np.linalg.norm(total-reaction[i]) <= 1e-17+128*eps*np.linalg.norm(total), 'reaction mismatch')
    require(len(seen) == n, 'missing particle trace')
    return normal, tangent, reaction, free, prescribed, intervals


def step_metrics(step):
    require(step['stage'] == 'ACCEPTED_COUPLED_DRIFT_KICK_DRIFT', 'rejected step in accepted stream')
    mass = np.asarray(step['mass_kg'], float); n = len(mass); dt = float(step['dt_s'])
    x0, x1 = (vec(step[k], n, k) for k in ('x_before_m', 'x_after_m'))
    v0, v1 = (vec(step[k], n, k) for k in ('v_before_m_s', 'v_after_m_s'))
    fint = vec(step['internal_force_N'], n, 'internal force')
    jn, jt, jr, free, prescribed, intervals = reconstruct_impulses(step)
    fn, ft, reaction = jn/dt, jt/dt, jr/dt
    residual = midpoint_residual(mass, x0, x1, v0, v1, dt, fint, fn+ft,
                                np.zeros_like(fint), reaction, free, prescribed)
    components = [np.linalg.norm(free*f) for f in (fint, fn, ft)]
    denominator = float(sum(components))
    static = float(np.linalg.norm(residual['quasistatic_free_residual_N']))
    dynamic = float(np.linalg.norm(residual['dynamic_residual_N']))
    force_scale = float(sum(np.linalg.norm(f) for f in (fint, fn, ft, reaction)))
    kinetic = float(.5*np.sum(mass[:, None]*v1*v1))
    kinematic = float(np.max(np.abs(residual['kinematic_residual_m'])))
    kinematic_bound = 1e-15+128*np.finfo(float).eps*max(float(np.max(np.abs(x0))), float(np.max(np.abs(x1))))
    require(dynamic <= 1e-12+1e-10*force_scale, 'dynamic momentum identity failed')
    require(kinematic <= kinematic_bound, 'saved drift identity failed')
    require(float(np.max(np.abs(residual['constraint_velocity_residual_m_s']))) <= 1e-12,
            'axial command not respected')
    inertia_norm = float(np.linalg.norm(free*residual['inertial_force_N']))
    inertial_matching = float(np.linalg.norm(free*residual['inertial_force_N']+residual['quasistatic_free_residual_N']))
    actual_work = (step['increments']['wall_work_J']+step['increments']['normal_relative_work_J']
                   -step['increments']['friction_raw_J'])
    endpoint_surrogate = float(np.sum((jn+jt)*(.5*(v0+v1))))
    ep0, ep1 = (np.asarray(step[k]) for k in ('ep_before', 'ep_after'))
    particle_norm = np.linalg.norm(residual['quasistatic_free_residual_N'], axis=1)
    dominant = int(np.argmax(particle_norm))
    return dict(step=step['step'], time_s=step['accepted_time_s'], dt_s=dt,
        contact_subintervals=intervals, kinetic_J=kinetic,
        new_plastic_particle_ids=np.asarray(step['particle_id'])[ep1 > ep0].tolist(),
        internal_free_force_norm_N=float(components[0]), normal_free_force_norm_N=float(components[1]),
        tangent_free_force_norm_N=float(components[2]), force_denominator_N=denominator,
        quasistatic_free_residual_norm_N=static,
        quasistatic_force_ratio=None if denominator == 0 else static/denominator,
        ratio_status='UNDEFINED_NO_FORCE' if denominator == 0 else 'DIAGNOSTIC_NOT_ACCEPTANCE_THRESHOLD',
        dynamic_residual_norm_N=dynamic, inertial_free_force_norm_N=inertia_norm,
        inertia_static_balance_defect_N=inertial_matching,
        maximum_kinematic_residual_m=kinematic,
        free_axis_reaction_force_norm_N=float(np.linalg.norm(residual['free_reaction_force_N'])),
        dominant_unbalanced_particle_id=int(step['particle_id'][dominant]),
        dominant_free_force_norm_N=float(particle_norm[dominant]),
        actual_contact_work_J=float(actual_work), endpoint_contact_work_surrogate_J=endpoint_surrogate,
        endpoint_minus_actual_contact_work_J=endpoint_surrogate-actual_work,
        surrogate_status='NOT_FRICTION_HEAT_NOT_PHYSICAL_SUBINTERVAL_WORK',
        plastic_heat_J=float(step['increments']['plastic_heat_J']),
        recorded_midpoint_minus_material_work_J=float(step['increments']['midpoint_stress_work_J']-
                                                     step['increments']['material_old_work_J']))


def assess(case):
    started = time.perf_counter(); case = Path(case).resolve()
    manifest, inputs = checked_inputs(case)
    rows = []
    with gzip.open(case/'steps.jsonl.gz', 'rt', encoding='utf-8') as stream:
        for line in stream:
            step = json.loads(line)
            require(step['step'] == len(rows)+1, 'gap/duplicate in accepted steps')
            require(step['particle_id'] == manifest['final_record']['particle_id'], 'changed ID order in this trajectory')
            rows.append(step_metrics(step))
    require(len(rows) == 67 and rows[-1]['time_s'] == 7e-7, 'incomplete accepted coverage')
    require(sum(r['contact_subintervals'] for r in rows) == 189382, 'contact coverage differs from complete readback')
    for rel, digest in inputs.items(): require(sha(ROOT/rel) == digest, 'input changed during read-only assessment '+rel)
    first_plastic = next((r for r in rows if r['new_plastic_particle_ids']), None)
    defined = [r['quasistatic_force_ratio'] for r in rows if r['quasistatic_force_ratio'] is not None]
    nonequilibrium_steps = [r['step'] for r in rows if r['quasistatic_free_residual_norm_N'] >
                            1e-12+1e-10*r['force_denominator_N']]
    final = manifest['final_diagnostics']
    p = json.loads((case/'resolved_config.json').read_text(encoding='utf-8'))['parameters']
    bulk = 2*p['shear']*(1+p['poisson'])/(3*(1-2*p['poisson']))
    wave = math.sqrt((bulk+4*p['shear']/3)/p['rho0'])
    return dict(status='READ_ONLY_ROUTE_ASSESSMENT_COMPLETE_NOT_NEW_PHYSICAL_SOLUTION',
        production_authorized=False, accepted_progress=False, source_case=str(case.relative_to(ROOT)),
        input_sha256=inputs, source_sha256={s:sha(ROOT/s) for s in SOURCES},
        numpy_version=np.__version__, classification='UNSCALED_INERTIAL_RESIDUAL_INTERFACE_ONLY',
        accepted_steps_checked=len(rows), contact_subintervals_checked=sum(r['contact_subintervals'] for r in rows),
        first=rows[0], first_plastic=first_plastic, last=rows[-1], all_steps=rows,
        summary=dict(quasistatic_ratio_min=min(defined), quasistatic_ratio_max=max(defined),
            steps_not_in_force_equilibrium=nonequilibrium_steps,
            maximum_dynamic_residual_norm_N=max(r['dynamic_residual_norm_N'] for r in rows),
            maximum_kinematic_residual_m=max(r['maximum_kinematic_residual_m'] for r in rows),
            cumulative_endpoint_minus_actual_contact_work_J=math.fsum(r['endpoint_minus_actual_contact_work_J'] for r in rows),
            initial_contact_potential_J=final['initial_normal_potential_J'], final_kinetic_J=final['kinetic_J']),
        nominal_scales=dict(sound_speed_m_s=wave, tool_rim_speed_m_s=p['omega']*p['tool_radius'],
            nominal_rim_mach=p['omega']*p['tool_radius']/wave, axial_feed_m_s=p['feed'],
            warning='Small nominal Mach number does not certify that the preloaded transient is equilibrated.'),
        decision=dict(direct_checkpoint_to_quasistatic=(
            'NOT_JUSTIFIED_TRANSIENT_INERTIA_MUST_NOT_BE_DISCARDED' if nonequilibrium_steps
            else 'NO_IMBALANCE_OBSERVED_BUT_WINDOW_INSUFFICIENT_FOR_PROCESS_CERTIFICATION'),
            next='Pure updated-configuration trial + coupled force/constraint residual, then bounded global implicit solve; no production switch.',
            friction_event_warning='Whole-step midpoint momentum alone does not integrate event-resolved motion/work/heat.',
            solved_global_nonlinear_step='NOT_RUN', process_time_afsd='NOT_RUN',
            q7_reuse='Bounded nonlinear trial/solve/commit pattern only, not fixed reference geometry or manufactured material.'),
        wall_seconds=time.perf_counter()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); require(not args.output.exists(), 'fresh output only')
    try: result = assess(args.case); code = 0
    except Exception as error:
        result = dict(status='FAIL_ROUTE_ASSESSMENT', error=repr(error), accepted_progress=False,
                      production_authorized=False, source_sha256={s:sha(ROOT/s) for s in SOURCES}); code = 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream: json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps({k:v for k,v in result.items() if k not in ('all_steps', 'input_sha256')}, allow_nan=False))
    return code


if __name__ == '__main__': raise SystemExit(main())

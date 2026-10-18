"""Independent readback of saved arrays only; never calls a physical solver."""
from pathlib import Path
import gzip
import hashlib
import json
import math
import sys
import time
import numpy as np
from pysph.solver.utils import load
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/contact_material_search_20261004'


def read(path):
    with (gzip.open(path, 'rt', encoding='utf-8') if path.suffix == '.gz'
          else path.open(encoding='utf-8')) as f:
        return json.load(f)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    started = time.perf_counter(); watched = {}
    for name in ('tests_review1', 'response', 'cost', 'root'):
        directory = BASE/name; manifest = read(directory/'manifest.json')
        for rel, h in manifest.get('source_sha256', {}).items():
            check(digest(ROOT/rel) == h, 'source drift '+rel)
            watched[ROOT/rel] = h
        for rel, h in manifest.get('input_sha256', {}).items():
            check(digest(ROOT/rel) == h, 'input drift '+rel)
            watched[ROOT/rel] = h
        for rel, h in manifest.get('artifact_sha256', {}).items():
            path = (directory/rel).resolve()
            check(path.is_relative_to(directory), 'unsafe relative artifact')
            check(digest(path) == h, 'artifact drift '+rel)
            watched[path] = h
        watched[directory/'manifest.json'] = digest(directory/'manifest.json')
    manifest = read(BASE/'root/manifest.json')
    parent = load(str(ROOT/manifest['parent_reference']['path']))['arrays']['solid']
    response = read(BASE/'response/manifest.json')
    probes = 0; maximum_block_difference = 0.; branches = 0
    with np.load(BASE/'response/response.npz', allow_pickle=False) as arrays:
        for i, pid in enumerate(arrays['particle_id']):
            r = read(BASE/f'response/contact_{pid:06d}.json.gz')
            probes += r['kick_calls']; branches += bool(r['branch_changes'])
            B = np.zeros((3, 3)); C = B.copy()
            if not r['probes']:
                check(np.array_equal(np.asarray(r['C']), 2*np.asarray(r['B'])), 'analytic free relation')
                B = np.asarray(r['B']); C = np.asarray(r['C'])
            else:
                for axis in range(3):
                    q = [p for p in r['probes'] if p['axis'] == axis]
                    if not q: continue
                    check(len(q) == 2 and q[0]['step_N'] > 0 and q[1]['step_N'] < 0, 'signed probe pair')
                    step = q[0]['step_N']-q[1]['step_N']
                    B[:, axis] = (np.array(q[0]['mean_velocity_m_s'])-q[1]['mean_velocity_m_s'])/step
                    C[:, axis] = (np.array(q[0]['velocity_m_s'])-q[1]['velocity_m_s'])/step
                    check(all(p['raw_heat_certificate_pass'] for p in q), 'probe signed heat')
            delta = max(np.max(abs(B-arrays['B'][i])), np.max(abs(C-arrays['C'][i])))
            maximum_block_difference = max(maximum_block_difference, float(delta))
            check(delta == 0., 'saved secant differs from raw signed differences')
    check(probes == response['kick_calls'] and branches == response['branch_change_particles'], 'response counts')
    saved = read(BASE/'root/latest_trial.json.gz')
    c, e = saved['candidate'], saved['evaluation']
    check(e is not None, 'no completed full trial to inspect')
    tolerances = dict(position_m=1e-12, velocity_m_s=1e-9, density_kg_m3=1e-10,
                      temperature_K=5e-12, mean_velocity_m_s=1e-9)
    norms = {}
    for k in tolerances:
        residual = np.asarray(c[k])-e['map'][k]
        check(np.array_equal(residual, e['residual'][k]), 'raw residual differs '+k)
        if k == 'density_kg_m3': residual = residual/parent.rho
        norms[k] = float(np.max(abs(residual)))
    scaled = max(norms[k]/tolerances[k] for k in norms)
    iterations = read(BASE/'root/iterations.json')
    row = iterations[-1]
    if row['evaluation_status'] == 'EVALUATED':
        check(norms == row['residual_norms'] and scaled == row['maximum_scaled_residual'], 'root residual report')
    mass = parent.m; before_velocity = np.column_stack((parent.u, parent.v, parent.w))
    momentum = []; mean_error = []; heat = 0.; negative_raw = 0; interval_count = 0
    for tr in e['boundary']['traces']:
        i = tr['particle_index']; r = tr['result']; h = e['dt_s']
        check(tr['particle_id'] == parent.particle_id[i], 'trace ID')
        impulse = (np.asarray(r['normal_impulse']).reshape(-1, 3).sum(axis=0)+
                   np.asarray(r['tangential_impulse']).reshape(-1, 3).sum(axis=0)+r['constraint_impulse'])
        momentum.append(mass[i]*(np.array(c['velocity_m_s'][i])-before_velocity[i])-
                        h*np.asarray(e['internal_force_N'][i])-impulse)
        mean = np.zeros(3); coverage = 0.
        for interval in r['intervals']:
            interval_count += 1; dt = interval['dt_s']; coverage += dt
            mean += dt*.5*(np.asarray(interval['v_before'])+interval['v_after'])/h
            q = np.asarray(interval['friction_heat_raw_J']); bound = np.asarray(interval['friction_rounding_bound_J'])
            reconstructed = -np.sum(np.asarray(interval['tangential_impulse']).reshape(-1, 3)*
                np.asarray(interval['tangential_relative_midpoint']).reshape(-1, 3), axis=1)
            check(np.allclose(q, reconstructed, atol=1e-25, rtol=1e-12), 'raw friction work identity')
            check(np.all(q >= -bound), 'negative heat outside original roundoff')
            heat += q.sum(); negative_raw += int(np.count_nonzero(q < 0))
        check(abs(coverage-h) <= 128*np.finfo(float).eps*h, 'contact coverage')
        mean_error.append(mean-np.asarray(c['mean_velocity_m_s'][i]))
    dep = np.asarray(e['material']['ep'])-parent.ep
    check(np.all(dep >= 0), 'decreasing plastic history')
    report = dict(status='PASS_SAVED_ARRAY_READBACK_NOT_PROCESS_VALIDATION',
        root_status=manifest['status'], accepted_steps=manifest['accepted_steps'],
        accepted_time_s=manifest['accepted_time_s'], full_root_evaluations=len(iterations),
        contact_kick_probes=probes, branch_change_particles=branches,
        secant_reconstruction_max_difference_s_kg=maximum_block_difference,
        final_raw_residual_norms=norms, final_maximum_scaled_residual=scaled,
        final_candidate_scope='ACCEPTED_NUMERICAL_ROOT' if manifest['accepted_steps'] else 'REJECTED_NOT_PHYSICAL_PROGRESS',
        candidate_momentum_residual_max_Ns=float(np.max(abs(np.asarray(momentum)))),
        candidate_mean_velocity_defect_max_m_s=float(np.max(abs(np.asarray(mean_error)))),
        raw_friction_J=float(heat), negative_roundoff_interval_entries=negative_raw,
        contact_intervals=interval_count, new_plastic_particle_count=int(np.count_nonzero(dep > 0)),
        new_plastic_strain_max=float(np.max(dep)), root_wall_seconds=manifest['wall_seconds'],
        memory_peak_bytes=manifest['peak_process_working_set_bytes'])
    if manifest['accepted_steps']:
        check(scaled <= 1., 'accepted unconverged root')
        accepted = load(str(BASE/'root/accepted.npz'))['arrays']['solid']
        check(np.array_equal(accepted.particle_id, parent.particle_id) and np.array_equal(accepted.m, mass), 'identity/mass changed')
        p = read(BASE/'root/resolved_config.json')['parameters']
        bar = parent.component_id == 1
        down = bar & (parent.z >= p['layer']) & (accepted.z < p['layer'])
        report.update(downward_shoulder_crossing_mass_kg=float(mass[down].sum()),
            free_bar_step_displacement_max_m=float(np.max(np.linalg.norm(
                np.column_stack((accepted.x-parent.x, accepted.y-parent.y, accepted.z-parent.z))
                [bar & (parent.actuator == 0)], axis=1), initial=0.)))
        nominal_duration = .00025/p['feed']
        report['nominal_feed_length_s'] = nominal_duration
        report['same_step_cost_illustration_not_prediction_days'] = (
            math.ceil(nominal_duration/e['dt_s'])*manifest['wall_seconds']/86400)
    old = read(ROOT/'outputs/process_global_20261003/root_attempt/iterations.json')
    report['old_500ns_final_scaled_residual'] = old[-1]['maximum_scaled_residual']
    report['old_500ns_evaluations'] = len(old)
    if len(iterations) <= len(old):
        report['old_500ns_same_evaluation_scaled_residual'] = old[len(iterations)-1]['maximum_scaled_residual']
    for path, h in watched.items(): check(digest(path) == h, 'reader changed evidence '+str(path))
    report.update(watched_files_unchanged=len(watched), wall_seconds=time.perf_counter()-started)
    output = BASE/'readback.json'
    check(not output.exists(), 'refuse readback overwrite')
    output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(report, allow_nan=False))


if __name__ == '__main__':
    main()

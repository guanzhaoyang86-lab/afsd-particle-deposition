"""Independent read-only event-average stage audit; never advances physics.

References are explicit algebra on saved particles and contact intervals. This
reader does not import UpdatedTrial, the acceptance routine, contact, material,
mechanical force, or thermal update implementations.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
import time
import numpy as np
from pysph.solver.utils import load

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CASE = ROOT/'outputs/updated_trial_20261003/global_attempt'
CLASSIFICATION = 'UPDATED_CONFIGURATION_EVENT_AVERAGE_STAGE_IMPLICIT_APPROXIMATE'
FIELDS = ('position_m', 'velocity_m_s', 'density_kg_m3', 'temperature_K', 'mean_velocity_m_s')
TOL = dict(position_m=1e-12, velocity_m_s=1e-9, density_kg_m3=1e-10,
           temperature_K=5e-12, mean_velocity_m_s=1e-9)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def require(value, message):
    if not value:
        raise ValueError(message)


def arr(value, shape=None):
    out = np.asarray(value, dtype=float)
    require(np.isfinite(out).all(), 'nonfinite readback array')
    if shape is not None:
        require(out.shape == shape, 'readback shape '+str((out.shape, shape)))
    return out


def close(a, b, label, atol=1e-18, rtol=1e-12):
    a, b = arr(a), arr(b)
    require(a.shape == b.shape and np.all(np.abs(a-b) <= atol+rtol*np.abs(b)), label)


def clean(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [clean(v) for v in value]
    return value


def verify_files(case, manifest):
    watched = {str((case/'manifest.json').relative_to(ROOT)): sha(case/'manifest.json')}
    for section in ('source_sha256', 'input_sha256'):
        require(bool(manifest.get(section)), 'missing source/input provenance')
        for name, digest in manifest[section].items():
            path = (ROOT/name).resolve()
            require(path.is_relative_to(ROOT) and sha(path) == digest, section+' drift '+name)
            watched[name] = digest
    for name, digest in manifest['artifact_sha256'].items():
        require(Path(name).name == name, 'unsafe artifact path')
        path = case/name
        require(sha(path) == digest, 'artifact drift '+name)
        watched[str(path.relative_to(ROOT))] = digest
    return watched


def audit(case):
    case = Path(case).resolve()
    require(case.is_relative_to(ROOT), 'project-local case required')
    manifest = json.loads((case/'manifest.json').read_text(encoding='utf-8'))
    require(manifest['classification'] == CLASSIFICATION, 'wrong stage classification')
    require(manifest['production_authorized'] is False, 'production not permitted')
    status = manifest['status']
    accepted = status == 'ACCEPTED_ONE_APPROXIMATE_GLOBAL_STAGE_NOT_SUSTAINED_DEPOSITION'
    require(accepted or status == 'BLOCKED_REJECTED_GLOBAL_STAGE', 'case is not a completed global attempt')
    watched = verify_files(case, manifest)
    ref = json.loads((case/'last_accepted_reference.json').read_text(encoding='utf-8'))
    source = (ROOT/ref['path']).resolve()
    require(source.is_relative_to(ROOT) and sha(source) == ref['sha256'], 'initial accepted source drift')
    original = load(str(source)); pa = original['arrays']['solid']; state = original['solver_data']
    resolved = json.loads((case/'resolved_config.json').read_text(encoding='utf-8'))
    p = resolved['parameters']
    require(p == state['parameters'] and state['count'] == 62 and state['accepted'] is True,
            'not the unchanged accepted plastic source')
    require(np.any(pa.ep > 0) and np.any(pa.S != 0), 'nonzero accepted plastic/stress history missing')
    config_payload = dict(parameters=p, stage_settings=resolved['settings'], classification=CLASSIFICATION)
    config_sha = hashlib.sha256(json.dumps(config_payload, sort_keys=True, allow_nan=False).encode()).hexdigest()
    require(config_sha == manifest['configuration_sha256'], 'configuration digest mismatch')
    report = dict(status='PASS_ACCEPTED_SINGLE_STAGE_READBACK' if accepted else 'PASS_BLOCKED_EVIDENCE_ONLY_NOT_ACCEPTED',
                  source_status=status, accepted_physical_stage=accepted, production_authorized=False,
                  source_step=state['count'], source_time_s=state['t'], source_checkpoint_sha256=ref['sha256'],
                  particle_count=len(pa.x), source_files_checked=len(manifest['source_sha256']),
                  no_solver_or_contact_or_material_calls=True)
    if not accepted:
        with gzip.open(case/'REJECTED_candidate.json.gz', 'rt', encoding='utf-8') as stream:
            rejected = json.load(stream)
        require(rejected['classification'] == 'REJECTED' and rejected['accepted'] is False and
                rejected['restart_allowed'] is False, 'rejected candidate falsely accepted')
        report.update(rejected_error=rejected['error'], rejected_candidate_scope='DEBUG_ONLY_NEVER_PHYSICAL_END_STATE')
    journal = case/'latest_uncommitted_trial.json.gz'
    if not journal.exists():
        require(not accepted, 'accepted case lacks final evaluated candidate')
        report['completed_candidate_checks'] = 'NOT_RUN_NO_COMPLETED_EVALUATION'
        return report, watched
    with gzip.open(journal, 'rt', encoding='utf-8') as stream:
        payload = json.load(stream)
    c, e = payload['candidate'], payload['evaluation']
    require(e['classification'] == CLASSIFICATION and e['accepted'] is False and e['restart_allowed'] is False,
            'trial archive is not uncommitted')
    require(e['material_scope'] == 'UNCHANGED_FKS_APPROXIMATE_PARAMETERS', 'manufactured material in real candidate')
    n = len(pa.x); m = pa.m.copy(); x0 = np.column_stack((pa.x, pa.y, pa.z)); v0 = np.column_stack((pa.u, pa.v, pa.w))
    ids = np.asarray(e['particle_id'], dtype=np.int64)
    require(np.array_equal(ids, pa.particle_id) and len(np.unique(ids)) == n, 'particle ID correspondence')
    candidate = {k: arr(c[k], (n,) if k in ('density_kg_m3', 'temperature_K') else (n, 3)) for k in FIELDS}
    mapping = {k: arr(e['map'][k], candidate[k].shape) for k in FIELDS}
    norms = {k: float(np.max(np.abs(candidate[k]-mapping[k])/(pa.rho if k == 'density_kg_m3' else 1.))) for k in FIELDS}
    scaled = max(norms[k]/TOL[k] for k in FIELDS)
    if accepted: require(scaled <= 1., 'accepted root exceeds frozen residual tolerance')
    dt = float(e['dt_s']); require(dt == 5e-8 and e['time_before_s'] == state['t'], 'wrong true dt/start state')
    close(e['time_after_s'], state['t']+dt, 'wrong end time', atol=0., rtol=0.)
    mid = e['midpoint']; xm = arr(mid['position_m'], (n, 3)); rho = arr(mid['density_kg_m3'], (n,))
    mean = candidate['mean_velocity_m_s']; force = arr(e['internal_force_N'], (n, 3))
    close(xm, .5*(x0+candidate['position_m']), 'midpoint positions')
    close(rho, .5*(pa.rho+candidate['density_kg_m3']), 'midpoint densities', atol=1e-12)
    close(mid['temperature_K'], .5*(pa.temperature+candidate['temperature_K']), 'midpoint temperature', atol=1e-12)
    close(mid['mass_kg'], m, 'unscaled masses', atol=0., rtol=0.)
    i, j = np.asarray(mid['i'], int), np.asarray(mid['j'], int)
    require(i.shape == j.shape and np.all(i < j) and np.all(i >= 0) and np.all(j < n), 'invalid unique pair indices')
    require(len(np.unique(np.column_stack((i, j)), axis=0)) == len(i), 'duplicate pair')
    g = arr(mid['g'], (len(i), 3)); B = arr(mid['B'], (n, 3, 3)); L = arr(mid['velocity_gradient_per_s'], (n, 3, 3))
    rij = xm[i]-xm[j]; radius = np.linalg.norm(rij, axis=1); h = .5*(pa.h[i]+pa.h[j]); q = radius/h
    df = np.where(q < 1., 3*q*q-6., np.where(q < 2., -3*(2-q)**2, 0.))
    independent_g = (15/(62*np.pi*h**4)*df/radius)[:, None]*rij
    close(g, independent_g, 'independent raw 3D radial gradient', atol=1e-3, rtol=2e-13)
    volume = m/rho; raw = np.zeros((n, 3, 3)); moment = raw.copy()
    pair = (mean[j]-mean[i])[:, :, None]*g[:, None, :]
    moment_pair = (-rij)[:, :, None]*g[:, None, :]
    for a, b in ((i, j), (j, i)):
        np.add.at(raw, a, volume[b, None, None]*pair)
        np.add.at(moment, a, volume[b, None, None]*moment_pair)
    close(B, np.linalg.inv(moment), 'independent moment inverse', atol=1e-12, rtol=2e-12)
    close(L, raw@B, 'independent current velocity gradient', atol=1e-10, rtol=1e-12)
    sigma = arr(mid['stress_Pa'], (n, 3, 3))
    spatial_power = float(np.sum(mean*force)+np.sum(volume[:, None, None]*sigma*L))
    power_scale = float(np.sum(np.linalg.norm(mean, axis=1)*np.linalg.norm(force, axis=1)))
    require(abs(spatial_power) <= 1e-12+1e-10*power_scale, 'spatial power conjugacy')
    net_force = force.sum(axis=0); net_torque = np.cross(xm, force).sum(axis=0)
    require(np.linalg.norm(net_force) <= 1e-12+1e-10*np.linalg.norm(force, axis=1).sum(), 'net internal force')
    require(np.linalg.norm(net_torque) <= 1e-14+1e-10*np.linalg.norm(np.cross(xm, force), axis=1).sum(), 'net internal torque')
    divergence = np.trace(L, axis1=1, axis2=2)
    close(mapping['density_kg_m3'], pa.rho*(1-.5*dt*divergence)/(1+.5*dt*divergence), 'continuity map', atol=1e-12)
    ep = arr(e['material']['ep'], (n,)); S = arr(e['material']['deviatoric_stress_Pa'], (n, 3, 3))
    dep = ep-pa.ep; require(np.all(dep >= 0.), 'plastic history decreased')
    plastic = np.sqrt(1.5*np.sum(S*S, axis=(1, 2)))*dep*volume
    returned = 1.5*p['shear']*dep*dep*volume
    close(e['material']['plastic_work_J'], plastic, 'same-history physical plastic work', atol=1e-18)
    close(e['material']['return_dissipation_J'], returned, 'same-history numerical return work', atol=1e-18)
    b = e['boundary']; qfr = np.zeros(n); contact_flux = np.zeros(n); event_mean = np.zeros((n, 3))
    endpoint = event_mean.copy(); external_impulse = event_mean.copy(); external_angular = np.zeros(3)
    wall_work = normal_relative_work = constraint_work = internal_work = 0.; intervals = 0; seen = set()
    max_interval_momentum = 0.; min_raw_q = None
    for trace in b['traces']:
        k = int(trace['particle_index']); require(k not in seen and ids[k] == trace['particle_id'], 'trace identity')
        seen.add(k); rows = trace['contact_rows']; nc = len(rows)
        normals = np.array([r['normal'] for r in rows], dtype=float).reshape(nc, 3)
        walls = np.array([r['wall_velocity'] for r in rows], dtype=float).reshape(nc, 3)
        fn = np.array([r['normal_force_N'] for r in rows], dtype=float)
        free = arr(trace['free_axes'], (3,)); prescribed = arr(trace['prescribed_velocity'], (3,))
        expected_cap = bool(pa.z0[k] > p['bar_length']-p['actuator_depth'] and pa.z0[k] <= p['bar_length'] and e['time_mid_s'] < p['release_time'])
        close(free, np.array([1., 1., 0. if expected_cap else 1.]), 'original-material-only actuator', atol=0., rtol=0.)
        for row in rows:
            require(row['particle_id'] == ids[k] and row['contact_id'] == str(ids[k])+':'+str(row['feature_id']), 'contact persistent identity')
            contact_flux[k] += p['contact_conductance']*p['dx']**2*(row['wall_temperature_K']-mid['temperature_K'][k])
        clock = 0.; previous = v0[k].copy(); total_jn = np.zeros((nc, 3)); total_jt = total_jn.copy(); reaction = np.zeros(3)
        frame = walls[0] if nc else v0[k]
        speed_scale = max(np.linalg.norm(v0[k]-frame), np.linalg.norm(walls-frame) if nc else 0.,
                          dt*np.linalg.norm(free/m[k]*(force[k]+(fn[:, None]*normals).sum(axis=0))), np.finfo(float).tiny)
        for interval in trace['result']['intervals']:
            intervals += 1; ih = float(interval['dt_s']); require(ih > 0, 'negative interval dt')
            close(interval['start_s'], clock, 'interval coverage start', atol=128*np.finfo(float).eps*dt)
            vb, va = arr(interval['v_before'], (3,)), arr(interval['v_after'], (3,)); vm = .5*(vb+va)
            close(vb, previous, 'interval velocity continuity', atol=1e-12)
            jn = np.asarray(interval['normal_impulse'], float).reshape(nc, 3)
            jt = np.asarray(interval['tangential_impulse'], float).reshape(nc, 3)
            jr = arr(interval['constraint_impulse'], (3,))
            close(jn, ih*fn[:, None]*normals, 'normal penalty impulse', atol=1e-18)
            relative = vm-walls
            tangential = relative-np.sum(relative*normals, axis=1)[:, None]*normals
            raw_q = -np.sum(jt*tangential, axis=1)
            rounding = 256*np.finfo(float).eps*np.linalg.norm(jt, axis=1)*(speed_scale+np.linalg.norm(vb-frame)+np.linalg.norm(va-frame)+2*np.linalg.norm(walls-frame, axis=1))
            close(interval['friction_rounding_bound_J'], rounding, 'independent roundoff certificate', atol=1e-29, rtol=1e-10)
            stored_q = np.asarray(interval['friction_heat_raw_J'], float)
            close(stored_q, raw_q, 'independent signed friction heat', atol=1e-24, rtol=1e-10)
            require(np.all(stored_q >= -rounding-1e-29), 'negative friction beyond raw rounding certificate')
            if nc: min_raw_q = float(stored_q.min()) if min_raw_q is None else min(min_raw_q, float(stored_q.min()))
            defect = m[k]*(va-vb)-ih*force[k]-jn.sum(axis=0)-jt.sum(axis=0)-jr
            max_interval_momentum = max(max_interval_momentum, float(np.max(np.abs(defect))))
            require(np.max(np.abs(defect)) <= 1e-17, 'subinterval momentum')
            require(np.max(np.abs(free*jr)) <= 1e-17, 'reaction on free axis')
            event_mean[k] += ih*vm/dt; qfr[k] += stored_q.sum()
            wall_work += float(np.sum((jn+jt)*walls)); constraint_work += float(jr@vm)
            normal_relative_work += float(np.sum(jn*(vm-walls))); internal_work += float(ih*force[k]@vm)
            total_jn += jn; total_jt += jt; reaction += jr
            previous = va; clock += ih
        close(clock, dt, 'contact interval duration', atol=128*np.finfo(float).eps*dt)
        close(np.asarray(trace['result']['normal_impulse']).reshape(nc, 3), total_jn, 'sum face normal impulses', atol=1e-18)
        close(np.asarray(trace['result']['tangential_impulse']).reshape(nc, 3), total_jt, 'sum face tangent impulses', atol=1e-18)
        close(trace['result']['constraint_impulse'], reaction, 'sum reaction impulses', atol=1e-18)
        endpoint[k] = previous; external_impulse[k] = total_jn.sum(axis=0)+total_jt.sum(axis=0)+reaction
        external_angular += np.cross(xm[k], external_impulse[k])
        if accepted:
            require(np.max(np.abs((1-free)*(candidate['velocity_m_s'][k]-prescribed))) <= 1e-12, 'saved axial constraint')
    require(len(seen) == n, 'missing particle traces')
    close(b['friction_raw_J'], qfr, 'raw friction per-particle totals', atol=1e-18)
    close(b['work_velocity'], event_mean, 'event-averaged actual velocities', atol=1e-12)
    close(mapping['velocity_m_s'], endpoint, 'actual endpoint map', atol=1e-12)
    close(mapping['mean_velocity_m_s'], event_mean, 'actual mean map', atol=1e-12)
    close(mapping['position_m'], x0+dt*event_mean, 'actual event-integrated position map', atol=1e-18)
    for key, value in [('friction_raw_J', qfr.sum()), ('wall_work_J', wall_work),
                       ('actuator_work_J', constraint_work), ('internal_work_J', internal_work),
                       ('normal_relative_work_J', normal_relative_work)]:
        close(b['ledger'][key], value, 'independent work '+key, atol=1e-18, rtol=1e-10)
    close(b['contact_power_W'], contact_flux, 'contact-only signed thermal flux', atol=1e-14)
    coefficient = volume[i]*volume[j]*(2*p['conductivity'])*(-np.sum(rij*g, axis=1))/np.sum(rij*rij, axis=1)
    require(np.all(coefficient >= 0), 'anti-diffusive conduction')
    pair = coefficient*(pa.temperature[j]-pa.temperature[i]); heat_power = np.zeros(n)
    np.add.at(heat_power, i, pair); np.add.at(heat_power, j, -pair)
    source_heat = p['plastic_heat_fraction']*plastic+p['heat_partition']*qfr+dt*(contact_flux+heat_power)
    close(e['thermal']['pair_power_W'], heat_power, 'independent signed pair conduction', atol=1e-14)
    close(e['thermal']['particle_energy_J']['total_heat_input_J'], source_heat, 'independent heat without residual/return', atol=1e-18)
    close(mapping['temperature_K'], pa.temperature+source_heat/(m*p['cp']), 'temperature source map', atol=1e-12)
    measured_heat = m*p['cp']*(candidate['temperature_K']-pa.temperature)
    thermal_residual = float(np.sum(measured_heat-source_heat))
    if accepted: require(abs(thermal_residual) <= 1e-12+1e-9*np.abs(source_heat).sum(), 'saved candidate thermal balance')
    x1, v1 = candidate['position_m'], candidate['velocity_m_s']
    angular_delta = np.sum(m[:, None]*(np.cross(x1, v1)-np.cross(x0, v0)), axis=0)
    extra = np.sum(m[:, None]*np.cross(x1-x0, .5*(v0+v1)), axis=0)
    angular_identity_defect = angular_delta-external_angular-dt*net_torque-extra
    dK = float(.5*np.sum(m[:, None]*(v1*v1-v0*v0)))
    map_dK = float(.5*np.sum(m[:, None]*(endpoint*endpoint-v0*v0)))
    nonlinear_work = dt*float(np.sum(force*(mean-event_mean)))
    features = {int(trace['particle_id']): [int(row['feature_id']) for row in trace['contact_rows']]
                for trace in b['traces']}
    dominant = {}
    for key in FIELDS:
        defect = candidate[key]-mapping[key]
        ranked = np.argsort(np.abs(defect).ravel())[-5:][::-1]
        entries = []
        for flat in ranked:
            index = np.unravel_index(int(flat), defect.shape)
            k = index[0]
            entries.append(dict(particle_id=int(ids[k]), particle_index=int(k),
                axis='xyz'[index[1]] if len(index) == 2 else None,
                component_id=int(pa.component_id[k]), contact_feature_ids=features[int(ids[k])],
                candidate_value=float(candidate[key][index]), mapped_value=float(mapping[key][index]),
                signed_residual=float(defect[index]),
                scaled_residual=float(abs(defect[index])/(pa.rho[k] if key == 'density_kg_m3' else 1.)/TOL[key])))
        dominant[key] = entries
    report.update(evaluated_candidate_scope='FINAL_ACCEPTED_CANDIDATE' if accepted else 'LAST_COMPLETED_UNCOMMITTED_EVALUATION_ONLY',
        dt_s=dt, evaluated_end_s=e['time_after_s'], residual_norms=norms, maximum_scaled_residual=scaled,
        dominant_residuals=dominant,
        neighbor_pairs=len(i), condition_max=float(np.linalg.cond(moment).max()), spatial_power_residual_W=spatial_power,
        net_internal_force_N=net_force, net_internal_torque_Nm=net_torque, mass_kg=float(m.sum()),
        contact_subintervals=intervals, minimum_raw_friction_J=min_raw_q, raw_friction_J=float(qfr.sum()),
        interval_momentum_max_Ns=max_interval_momentum, thermal_residual_J=thermal_residual,
        physical_plastic_J=float(plastic.sum()), numerical_return_J=float(returned.sum()),
        plastic_increment_particle_ids=ids[dep > 0], delta_K_J=dK, candidate_minus_contact_map_delta_K_J=dK-map_dK,
        nonlinear_internal_work_defect_J=nonlinear_work, angular_momentum_change_kg_m2_s=angular_delta,
        external_midpoint_angular_impulse_Nms=external_angular, angular_residual_Nms=angular_delta-external_angular,
        event_average_kinematic_angular_term_Nms=extra, angular_identity_defect_Nms=angular_identity_defect,
        angular_status='APPROXIMATE_FROZEN_POSITION_MOMENT_NOT_CORRECTED', total_free_energy_certification='NOT_RUN',
        sustained_deposition='NOT_RUN', same_end_time_dt_sensitivity='NOT_RUN')
    if accepted:
        saved = load(str(case/'accepted_snapshot.npz')); after = saved['arrays']['solid']; newstate = saved['solver_data']
        require(newstate['accepted'] is True and newstate['restart_allowed'] is False, 'snapshot scope mislabelled')
        require(newstate['source_sha256'] == manifest['source_sha256'] and newstate['configuration_sha256'] == manifest['configuration_sha256'], 'snapshot source mismatch')
        for key in ('particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','as_rs_tag','m','h','x0','y0','z0'):
            close(getattr(after, key), getattr(pa, key), 'immutable identity/mass '+key, atol=0., rtol=0.)
        for name, expected in [('x',x1[:,0]),('y',x1[:,1]),('z',x1[:,2]),('u',v1[:,0]),('v',v1[:,1]),('w',v1[:,2]),
                               ('rho',candidate['density_kg_m3']),('temperature',candidate['temperature_K']),('ep',ep),('S',S.ravel())]:
            close(getattr(after, name), expected, 'saved candidate '+name, atol=0., rtol=0.)
        checks = json.loads((case/'acceptance.json').read_text(encoding='utf-8'))
        for key in ('delta_K_J','candidate_minus_contact_map_delta_K_J','nonlinear_internal_work_defect_J',
                    'angular_momentum_change_kg_m2_s','external_midpoint_angular_impulse_Nms','angular_residual_Nms',
                    'event_average_kinematic_angular_term_Nms','angular_identity_defect_Nms'):
            close(checks[key], report[key], 'independent saved diagnostic '+key, atol=1e-18, rtol=1e-9)
        report['snapshot'] = 'PASS_REAL_CANDIDATE_IDENTITY_FIELDS_NOT_RESTART_VALIDATED'
    return report, watched


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=DEFAULT_CASE)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); output = args.output.resolve()
    require(output.is_relative_to(ROOT) and not output.exists(), 'fresh project-local audit output required')
    started = time.perf_counter(); code = 0
    try:
        report, watched = audit(args.case)
        for name, digest in watched.items(): require(sha(ROOT/name) == digest, 'read-only input drift '+name)
        report['all_input_hashes_unchanged'] = True
    except Exception as error:
        report = dict(status='FAIL_INDEPENDENT_READBACK', error=repr(error), production_authorized=False)
        code = 2
    report.update(reader_sha256=sha(__file__), wall_seconds=time.perf_counter()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(clean(report), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(clean(report), allow_nan=False), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())

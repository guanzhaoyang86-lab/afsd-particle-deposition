"""Independent saved-array audit: planar point masses, NOT a 3D SPH solid.

Does not import the run, contact, geometry, material, integrator, or RK4 code.
Literal SI inputs and tolerances were frozen before the physical runs.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud

CLASSIFICATION = 'ENGINEERING_APPROXIMATE_CONTACT_COMPONENT_NOT_AFSD'
IDS = [101, 102, 103, 104]
X0 = [[.01, 0., .0013], [0., .012, .0013], [-.01, 0., .0013], [0., -.012, .0013]]
PARAMETERS = dict(initial_position=X0, initial_velocity=[0., 0., 0.], mass=.001,
    normal_force=1., mu=.2, omega=20., traverse=0., layer=.0015, tool_top=.02,
    tool_radius=.02, bore_width=.008, dx=.0005, substrate_depth=.01,
    tool_temperature=293.15, base_temperature=293.15, end_time=.0008)
SOURCE_PATHS = {'scripts/probe_rotating_shoulder_component.py',
    'configs/benchmarks/rotating_shoulder_component.json', 'docs/ROTATING_CONTACT_CONTRACT.md',
    'src/afsd_backend/pysph/rigid_contact_kinematics.py', 'src/afsd_backend/pysph/geometry.py',
    'src/afsd_backend/pysph/contact_kick.py', 'src/afsd_backend/pysph/contact_cone.py',
    'src/afsd_backend/pysph/integrator.py', 'src/afsd_backend/pysph/output.py',
    'src/afsd_scaffold/vtk_legacy.py', 'third_party/PySPH/pysph/sph/integrator.py',
    'third_party/PySPH/pysph/sph/integrator_step.py'}
ENERGY_TOL = 5e-12


def require(condition, message):
    if not condition: raise ValueError(message)


def arr(value, shape, name):
    result = np.asarray(value, dtype=float)
    require(result.shape == shape and np.isfinite(result).all(), 'shape/nonfinite '+name)
    return result


def close(actual, expected, tolerance, name):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), 'shape/nonfinite '+name)
    error = float(np.max(np.abs(a-b))) if a.size else 0.
    require(error <= tolerance, name+': absolute defect '+str(error)+' > '+str(tolerance))
    return error


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def file_hash(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    def reject(value): raise ValueError('nonfinite JSON constant '+value)
    return json.loads(Path(path).read_text(encoding='utf-8'), parse_constant=reject)


def audit_sources(sources):
    require(set(sources) == SOURCE_PATHS, 'source inventory changed')
    for name, digest in sources.items(): require(file_hash(ROOT/name) == digest, 'source hash '+name)


def audit_configuration(config, dt):
    require(config['classification'] == CLASSIFICATION, 'component classification')
    require(config['parameters'] == PARAMETERS, 'frozen physical parameters changed')
    require(config['dt_s'] == dt and dt in (2e-5, 1e-5), 'unregistered dt')
    require(config['free_axes'] == [1, 1, 0] and config['ids'] == IDS, 'free axes/IDs changed')
    require(config['initial_positions_m'] == X0 and config['initial_velocity_m_s'] == [0., 0., 0.], 'initial state changed')
    require(config['physical_scaling'] == 'NONE' and config['tags_are_observation_only'] is True, 'scaling/passive tags')
    require(config['normal_law'] == 'PRESCRIBED 1 N; NOT PENALTY OR DEFORMABLE SOLID', 'normal law')
    for name in ('temperature', 'density', 'stress', 'strain'):
        require(config['field_status'][name] == 'NOT_SOLVED', 'unsolved field misclassified '+name)
    source = config['source_configuration']
    require(source['classification'] == CLASSIFICATION and source['production_authorized'] is False and source['units'] == 'SI', 'source config policy')
    require({k: a['value'] for k, a in source['parameters'].items()} == PARAMETERS, 'source physical values')
    require(all(a.get('unit') and a.get('status') for a in source['parameters'].values()), 'parameter unit/status missing')
    return source


def audit_interval(interval, before, wall, dt, start):
    """Standalone arithmetic audit; synthetic tests use independently written values."""
    h = float(interval['dt_s'])
    require(np.isfinite(h) and h > 0 and h <= dt, 'invalid subinterval duration')
    close(interval['start_s'], start, 1e-15, 'subinterval chronology')
    vb = arr(interval['v_before'], (3,), 'interval before')
    va = arr(interval['v_after'], (3,), 'interval after')
    close(vb, before, 1e-12, 'velocity interval continuity')
    jn = arr(interval['normal_impulse'], (1, 3), 'normal impulse')[0]
    jt = arr(interval['tangential_impulse'], (1, 3), 'tangent impulse')[0]
    jc = arr(interval['constraint_impulse'], (3,), 'constraint impulse')
    close(jn, [0., 0., -h], 1e-14, 'prescribed normal load impulse')
    close(jc, [0., 0., h], 1e-14, 'normal actuator reaction')
    close(jt[2], 0., 1e-14, 'tangent impulse has normal component')
    close(.001*(va-vb), jn+jt+jc, 1e-14, 'impulse momentum')
    midpoint = .5*(vb+va)
    relative = midpoint-wall; relative[2] = 0.
    close(interval['tangential_relative_midpoint'], [relative], 1e-12, 'relative midpoint')
    q = -float(jt@relative)
    raw_q = float(arr(interval['friction_heat_raw_J'], (1,), 'raw friction')[0])
    close(raw_q, q, 1e-14, 'signed friction impulse work')
    # Sign must not be lost even when its magnitude falls below an energy tolerance.
    if q < 0: require(raw_q < 0, 'negative raw friction was clipped or made absolute')
    rounding = float(arr(interval['friction_rounding_bound_J'], (1,), 'rounding bound')[0])
    require(rounding >= 0 and raw_q >= -rounding, 'negative friction exceeds original rejection bound')
    ww, wc, wn = float((jn+jt)@wall), float(jc@midpoint), float(jn@midpoint)
    dk = .0005*float(va@va-vb@vb)
    for name, expected in [('wall_work_J', ww), ('constraint_work_J', wc), ('other_work_J', 0.),
                           ('delta_kinetic_J', dk), ('kinetic_work_residual_J', 0.)]:
        close(interval[name], expected, 1e-14, name)
    close(interval['normal_work_J'], [wn], 1e-14, 'normal work')
    close(dk-ww-wc+q, 0., ENERGY_TOL, 'interval work identity')
    return dict(v_after=va, dt=h, q=q, wall=ww, constraint=wc, normal=wn, dk=dk,
                jn=jn, jt=jt, reaction=jc, negative_raw=int(raw_q < 0), raw=raw_q)


def audit_sliding_friction(interval, kick_initial, wall):
    """Additional algebraic law audit, not a changed physical acceptance gate.

    Single shoulder, xy-free, no other projected force: normal motion is constrained.
    Recompute the existing roundoff expression instead of trusting a saved bound.
    """
    vb, va = np.asarray(interval['v_before']), np.asarray(interval['v_after'])
    jt = np.asarray(interval['tangential_impulse'])[0]
    u = .5*(vb+va)-wall; u[2] = 0.
    speed = float(np.linalg.norm(u))
    require(speed >= .03, 'component left independently verified sliding-only domain')
    expected = -.2*interval['dt_s']*u/speed
    close(jt, expected, 1e-14, 'Coulomb sliding impulse direction/cap')
    scale = float(np.linalg.norm(kick_initial-wall))
    close(interval['relative_speed_scale_m_s'], scale, 1e-13, 'kick relative speed scale')
    close(interval['reference_frame_velocity'], wall, 1e-13, 'kick Galilean frame')
    rounding = 256*np.finfo(float).eps*np.linalg.norm(jt)*(scale+np.linalg.norm(vb-wall)+np.linalg.norm(va-wall))
    saved = float(interval['friction_rounding_bound_J'][0])
    close(saved, float(rounding), max(1e-35, 1e-12*float(rounding)), 'original roundoff bound')


def audit_trace(trace, xbefore, vbefore, xafter, vafter, dt, step, pid):
    require(trace['particle_id'] == pid and trace['feature_id'] == 1 and trace['step'] == step, 'trace identity')
    require(trace['stage'] == 'MIDPOINT_FROZEN_GEOMETRY_KICK', 'trace stage')
    close(trace['time_mid_s'], (step-.5)*dt, 1e-15, 'midpoint time')
    xmid = arr(trace['position_m'], (3,), 'application point')
    close(xmid, xbefore+.5*dt*vbefore, 1e-15, 'first drift')
    close(xafter, xmid+.5*dt*vafter, 1e-15, 'second drift')
    row = trace['contact_row']
    require(row['particle_id'] == pid and row['feature_id'] == 1 and row['contact_id'] == str(pid)+':1', 'contact row identity')
    close(row['application_point'], xmid, 1e-15, 'force application')
    close(row['contact_point'], [xmid[0], xmid[1], .0015], 1e-15, 'shoulder projection')
    close(row['normal'], [0., 0., -1.], 1e-15, 'shoulder normal')
    close(row['body_origin'], [0., 0., .0015], 1e-15, 'rigid origin')
    close(row['translation_velocity'], [0., 0., 0.], 1e-15, 'rigid translation')
    close(row['omega'], [0., 0., 20.], 1e-15, 'rigid omega')
    close(row['gap'], .0015-xmid[2], 1e-15, 'shoulder gap')
    close(row['penetration'], .00025-(.0015-xmid[2]), 1e-15, 'shoulder overlap')
    angle = (step-.5)*dt*20.
    c, s = np.cos(angle), np.sin(angle)
    local = np.array([c*xmid[0]+s*xmid[1], -s*xmid[0]+c*xmid[1]])
    require(np.max(np.abs(local)) >= .004 and np.linalg.norm(local) <= .02, 'shoulder material outside finite annulus')
    wall = np.array([-20*xmid[1], 20*xmid[0], 0.])
    close(row['wall_velocity'], wall, 1e-13, 'independent rotating wall velocity')
    close(row['wall_velocity_at_face'], wall, 1e-13, 'shoulder face velocity')
    close(trace['material_velocity_before_m_s'], vbefore, 1e-12, 'trace starting velocity')
    result = trace['result']; close(result['velocity'], vafter, 1e-12, 'final kick velocity')
    close(result['constraint_projector'], np.diag([1., 1., 0.]), 1e-15, 'constraint mobility projection')
    sums = dict(q=0., wall=0., constraint=0., normal=0., dk=0., negative_raw=0)
    jn, jt, reaction = np.zeros(3), np.zeros(3), np.zeros(3)
    elapsed = 0.; previous = vbefore
    raw_min = float('inf')
    require(bool(result['intervals']), 'empty accepted kick')
    for interval in result['intervals']:
        a = audit_interval(interval, previous, wall, dt, elapsed)
        audit_sliding_friction(interval, vbefore, wall)
        previous = a['v_after']; elapsed += a['dt']; raw_min = min(raw_min, a['raw'])
        for name in sums: sums[name] += a[name]
        jn += a['jn']; jt += a['jt']; reaction += a['reaction']
    close(elapsed, dt, 1e-15, 'kick duration')
    close(previous, vafter, 1e-12, 'kick end velocity')
    for key, expected in [('friction_heat_J', [sums['q']]), ('normal_impulse', [jn]),
                          ('tangential_impulse', [jt]), ('constraint_impulse', reaction),
                          ('wall_work_J', sums['wall']), ('constraint_work_J', sums['constraint']),
                          ('other_work_J', 0.), ('delta_kinetic_J', sums['dk']),
                          ('kinetic_work_residual_J', 0.)]: close(result[key], expected, 1e-14, 'kick sum '+key)
    sums.update(torque_impulse_z=float(np.cross(xmid, jn+jt+reaction)[2]), min_raw=raw_min)
    return sums


def audit_run(path, dt):
    path = Path(path); manifest = load(path/'manifest.json'); config = manifest['configuration']
    require(manifest['status'] == 'COMPLETE_MANUFACTURED_COMPONENT_NOT_AFSD', 'unequal/failed endpoint cannot be compared')
    require(manifest['classification'] == CLASSIFICATION and manifest['production_authorized'] is False, 'manifest policy')
    audit_configuration(config, dt); audit_sources(manifest['source_sha256'])
    require(manifest['configuration_sha256'] == canonical_hash(config), 'configuration hash')
    files = {str(p.relative_to(path)).replace('\\', '/') for p in path.rglob('*') if p.is_file() and p.name != 'manifest.json'}
    require(files == set(manifest['artifact_sha256']), 'artifact inventory')
    for name, digest in manifest['artifact_sha256'].items(): require(file_hash(path/name) == digest, 'artifact hash '+name)
    steps = round(.0008/dt)
    require(manifest['accepted_steps'] == steps, 'accepted step count')
    close(manifest['accepted_time_s'], .0008, 1e-15, 'common physical endpoint')
    require(0 < manifest['wall_seconds'] <= 60 and manifest['maximum_wall_seconds'] <= 60, 'wall budget')
    history = load(path/'history.json')
    traces = [json.loads(line) for line in (path/'contact_intervals.jsonl').read_text(encoding='utf-8').splitlines()]
    require(len(history) == steps+1 and len(traces) == steps, 'history/journal length')
    require('last_accepted.npz' in files and not any('reject' in name for name in files), 'accepted checkpoint inventory')
    totals = dict(q=0., wall=0., constraint=0., normal=0., dk=0., negative_raw=0)
    torque = 0.; minimum_q = float('inf'); residual_max = 0.; previous = None
    for step, record in enumerate(history):
        require(record['step'] == step and record['particle_id'] == IDS, 'persistent history ID/order')
        require(record['tracer_type'] == [0, 1, 2, 1] and record['as_rs_tag'] == [0]*4, 'passive label history')
        close(record['time_s'], step*dt, 1e-15, 'history time')
        x = arr(record['position_m'], (4, 3), 'history position')
        v = arr(record['velocity_m_s'], (4, 3), 'history velocity')
        close(x[:, 2], [.0013]*4, 1e-14, 'fixed normal geometry')
        close(v[:, 2], [0.]*4, 1e-14, 'fixed normal velocity')
        close(record['tool_origin_m'], [0., 0., .0015], 1e-15, 'tool origin')
        close(record['tool_angle_rad'], 20*step*dt, 1e-15, 'tool angle')
        kinetic = .0005*float(np.sum(v*v)); close(record['kinetic_J'], kinetic, 1e-14, 'actual kinetic energy')
        if step == 0:
            close(x, X0, 0., 'literal initial geometry'); close(v, np.zeros((4, 3)), 0., 'initial rest')
        else:
            batch = traces[step-1]
            require(len(batch) == 4 and {r['particle_id'] for r in batch} == set(IDS), 'journal ID inventory')
            for trace in batch:
                i = IDS.index(trace['particle_id'])
                a = audit_trace(trace, previous[0][i], previous[1][i], x[i], v[i], dt, step, IDS[i])
                for key in totals: totals[key] += a[key]
                torque += a['torque_impulse_z']; minimum_q = min(minimum_q, a['min_raw'])
        cumulative = record['cumulative']
        for key, expected in [('friction_raw_J', totals['q']), ('wall_work_J', totals['wall']),
            ('constraint_work_J', totals['constraint']), ('normal_work_J', totals['normal']),
            ('other_work_J', 0.), ('delta_kinetic_J', totals['dk']), ('contact_kinetic_residual_J', 0.)]:
            close(cumulative[key], expected, 1e-13, 'cumulative '+key)
        residual = kinetic+totals['q']-totals['wall']-totals['constraint']
        close(record['global_energy_residual_J'], residual, 1e-13, 'saved energy residual')
        residual_max = max(residual_max, close(residual, 0., ENERGY_TOL, 'global energy residual'))
        lz = .001*float(np.sum(np.cross(x, v)[:, 2]))
        close(lz, torque, 2.5e-13, 'angular impulse conservation')
        close(totals['wall'], 20*lz, ENERGY_TOL, 'rotating work equals omega delta Lz')
        vtk = read_point_cloud(path/f'frame_{step:04d}.vtk')
        require(set(vtk['scalars']) == {'particle_id', 'tracer_type', 'as_rs_tag', 'mass_kg'} and set(vtk['vectors']) == {'velocity_m_per_s'}, 'VTK fabricated/omitted field inventory')
        require(CLASSIFICATION in vtk['title'] and 'STRESS_THERMAL_NOT_SOLVED' in vtk['title'], 'VTK classification')
        close(vtk['points'], x, 0., 'VTK positions'); close(vtk['vectors']['velocity_m_per_s'], v, 0., 'VTK velocity')
        for key in ('particle_id', 'tracer_type', 'as_rs_tag'): require(vtk['scalars'][key] == record[key], 'VTK '+key)
        close(vtk['scalars']['mass_kg'], [.001]*4, 0., 'unscaled prescribed particle masses')
        previous = (x, v)
    require(manifest['final_record'] == history[-1], 'final manifest record differs from history')
    return dict(path=str(path), dt_s=dt, accepted_steps=steps, endpoint_s=.0008,
        wall_seconds=manifest['wall_seconds'], configuration=config, source_sha256=manifest['source_sha256'],
        final_position_m=previous[0].tolist(), final_velocity_m_s=previous[1].tolist(),
        friction_raw_J=totals['q'], wall_work_J=totals['wall'], kinetic_J=kinetic,
        total_mass_kg=.004, maximum_energy_residual_J=residual_max,
        angular_momentum_z_kg_m2_s=lz, torque_impulse_z_N_m_s=torque,
        negative_raw_intervals=totals['negative_raw'], minimum_raw_interval_friction_J=minimum_q,
        checkpoint='HASH_VERIFIED; LOADER/ID_RESTART COVERED BY SEPARATE SHORT TEST, NOT RELOADED HERE')


def audit_reference(path, dt, sources, source_config):
    result = load(path)
    require(result['method'] == 'INDEPENDENT_RK4_CONTINUOUS_SLIP_ODE_NOT_CONTACT_KICK', 'reference method')
    require(result['dt_s'] == dt and result['steps'] == round(.0008/dt), 'reference step count')
    close(result['time_s'], .0008, 1e-15, 'reference endpoint')
    require(result['source_sha256'] == sources, 'reference source mismatch')
    require(result['configuration_sha256'] == canonical_hash(result['configuration']), 'reference config hash')
    require(result['configuration']['source_configuration'] == source_config, 'reference physical configuration')
    require(result['configuration']['reference_dt_s'] == dt and result['configuration']['reference_parameters'] == result['reference_parameters'], 'reference metadata disagreement')
    expected = dict(mass_kg=.001, omega_rad_s=[0., 0., 20.], friction_force_cap_N=.2,
        normal_force_N=1., mu=.2, traverse_m_s=0., initial_position_m=X0,
        initial_velocity_m_s=[0., 0., 0.], full_interval_end_s=.0008,
        requested_end_s=.0008, source='INDEPENDENT_FROZEN_LITERALS')
    require(result['reference_parameters'] == expected, 'independent reference literal parameters')
    require(np.isfinite(result['minimum_slip_m_s']) and result['minimum_slip_m_s'] >= .03, 'reference left continuous-slip domain')
    require(0 < result['wall_seconds'] <= 60, 'reference wall budget')
    x = arr(result['position_m'], (4, 3), 'reference position')
    v = arr(result['velocity_m_s'], (4, 3), 'reference velocity')
    q = arr(result['friction_raw_J'], (4,), 'reference friction')
    w = arr(result['wall_work_J'], (4,), 'reference wall work')
    close(x[:, 2], [.0013]*4, 1e-14, 'reference z'); close(v[:, 2], [0.]*4, 1e-14, 'reference vz')
    close(result['kinetic_J'], .0005*np.sum(v*v, axis=1), 1e-14, 'reference K')
    close(.0005*np.sum(v*v, axis=1)+q-w, np.zeros(4), ENERGY_TOL, 'reference energy')
    return result


def reference_refinement(a, b):
    differences = {}
    for key, bound in [('position_m', 1e-11), ('velocity_m_s', 1e-8),
                       ('friction_raw_J', 1e-11), ('wall_work_J', 1e-11)]:
        delta = np.asarray(a[key])-np.asarray(b[key])
        error = float(np.linalg.norm(delta, axis=1).max()) if delta.ndim == 2 else float(np.abs(delta).max())
        require(np.isfinite(error) and error <= bound, 'reference refinement '+key)
        differences[key] = dict(maximum_difference=error, tolerance=bound)
    return differences


def compare_to_reference(run, reference):
    result = {}
    for key, field, bound in [('position_m', 'final_position_m', 1e-8), ('velocity_m_s', 'final_velocity_m_s', 1e-5)]:
        error = float(np.linalg.norm(np.asarray(run[field])-np.asarray(reference[key]), axis=1).max())
        require(np.isfinite(error) and error <= bound, 'component/reference '+key)
        result[key] = dict(maximum_vector_error=error, tolerance=bound)
    for key in ('friction_raw_J', 'wall_work_J'):
        expected = float(np.sum(reference[key])); error = abs(run[key]-expected)
        bound = 1e-10+2e-4*abs(expected)
        require(np.isfinite(error) and error <= bound, 'component/reference '+key)
        result[key] = dict(component=run[key], reference=expected, absolute_error=error, tolerance=bound)
    return result


def inspect(coarse, fine, reference_coarse, reference_fine):
    a, b = audit_run(coarse, 2e-5), audit_run(fine, 1e-5)
    ca, cb = dict(a['configuration']), dict(b['configuration']); ca.pop('dt_s'); cb.pop('dt_s')
    require(ca == cb and a['source_sha256'] == b['source_sha256'], 'physical/source differences between dt cases')
    r1 = audit_reference(reference_coarse, 1e-7, a['source_sha256'], ca['source_configuration'])
    r2 = audit_reference(reference_fine, 5e-8, a['source_sha256'], ca['source_configuration'])
    refinement = reference_refinement(r1, r2)
    comparison = {'coarse': compare_to_reference(a, r2), 'fine': compare_to_reference(b, r2)}
    for run in (a, b): run.pop('configuration'); run.pop('source_sha256')
    return dict(status='PASS_COMPONENT_READBACK_NOT_AFSD',
        classification='PLANAR_POINT_MASS_DYNAMICS_WITH_3D_CONTACT_INTERFACE_NOT_SPH_SOLID',
        production_authorized=False, coarse=a, fine=b, reference_refinement=refinement,
        component_reference=comparison, reference_minimum_slip_m_s=r2['minimum_slip_m_s'],
        claims_excluded=['3D_SPH_SOLID', 'NORMAL_OPENING', 'ROTATING_BORE_FORCE_TRANSFER',
            'THERMAL', 'FKS', 'FEEDING', 'AFSD_DEPOSITION', 'CONVERGENCE_ORDER'],
        input_file_sha256={str(p): file_hash(p) for p in [Path(coarse)/'manifest.json', Path(fine)/'manifest.json', Path(reference_coarse), Path(reference_fine)]})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('coarse', 'fine', 'reference-coarse', 'reference-fine', 'output'): p.add_argument('--'+name, type=Path, required=True)
    args = p.parse_args()
    require(not args.output.exists(), 'fresh output JSON required')
    try:
        result = inspect(args.coarse, args.fine, args.reference_coarse, args.reference_fine); code = 0
    except (ValueError, KeyError, OSError, TypeError, IndexError) as error:
        result = dict(status='FAIL_COMPONENT_READBACK_NOT_AFSD', error=str(error), production_authorized=False); code = 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, allow_nan=False), flush=True)
    return code


if __name__ == '__main__': raise SystemExit(main())

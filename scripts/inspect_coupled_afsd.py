"""Independent saved-array coupled AFSD audit; never runs a physical solver.

References are literal geometry, radial calculus, balance identities and array
sums. No solver/material/contact/thermal implementation generates an answer.
An accepted short interval is NOT evidence of sustained deposition.
"""
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
from afsd_scaffold.vtk_legacy import read_point_cloud

CLASSIFICATION = 'ENGINEERING_APPROXIMATE_COUPLED_AFSD_LOCAL_DEBUG'
MASS = 2700*.0005**3
G = 26.3e9
BULK = 2*G*1.33/(3*.34)
OMEGA = 10*np.pi
SOURCE_PATHS = [
    'scripts/run_coupled_afsd.py', 'configs/benchmarks/coupled_afsd.json',
    'configs/debug/starved.json', 'configs/base.json', 'docs/COUPLED_AFSD_CONTRACT.md',
    *['src/afsd_backend/pysph/'+name+'.py' for name in (
        'coupled_solver', 'coupled_boundary', 'coupled_thermal', 'axial_actuator',
        'solver', 'material', 'mechanics', 'kernel', 'contact', 'contact_kick',
        'contact_cone', 'geometry', 'rigid_contact_kinematics', 'thermal',
        'particles', 'integrator', 'output', 'config_adapter', 'application')],
    'src/afsd_scaffold/config.py', 'src/afsd_scaffold/vtk_legacy.py',
    'third_party/PySPH/pysph/sph/integrator.py',
    'third_party/PySPH/pysph/sph/integrator_step.py',
    'third_party/PySPH/pysph/base/nnps.py',
    'third_party/PySPH/pysph/base/nnps_base.pyx',
    'third_party/PySPH/pysph/base/linked_list_nnps.pyx']
SHORT_WINDOW_SOURCES = ['scripts/run_coupled_short_window.py', 'docs/COUPLED_AFSD_SHORT_WINDOW.md']
ORIGINAL_PROBE = ROOT/'outputs/coupled_afsd_20261002/cost_probe'


def require(condition, message):
    if not condition: raise ValueError(message)


def array(value, shape, name):
    a = np.asarray(value, float)
    require(a.shape == shape and np.isfinite(a).all(), 'shape/nonfinite '+name)
    return a


def close(actual, expected, atol, name):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), 'shape/nonfinite '+name)
    defect = float(np.max(np.abs(a-b))) if a.size else 0.
    require(defect <= atol, name+': defect '+str(defect)+' > '+str(atol))
    return defect


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'),
        parse_constant=lambda v: (_ for _ in ()).throw(ValueError('nonfinite JSON '+v)))


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def initial_geometry():
    substrate = np.array([[x, y, z] for x in (np.arange(8)+.5)*.0005-.002
        for y in (np.arange(8)+.5)*.0005-.002 for z in (np.arange(3)+.5)*.0005-.0015])
    bar = np.array([[x, y, z] for x in (np.arange(5)+.5)*.0005-.00125
        for y in (np.arange(5)+.5)*.0005-.00125 for z in (np.arange(6)+.5)*.0005])
    return np.vstack((substrate, bar))


def independent_short_choice(probe):
    """Enumerate a preregistered 100 ns grid using literal unchanged bounds."""
    dt = np.asarray(probe['accepted_dt_s'], float)
    wall = np.asarray(probe['step_wall_seconds'], float)
    require(dt.shape == wall.shape == (5,) and np.isfinite(dt).all() and np.isfinite(wall).all()
        and np.all(dt > 0) and np.all(wall > 0), 'five actual positive finite cost steps')
    minimum_dt, slowest = float(dt.min()), float(wall.max())
    allowed = None
    # Evaluate the cost directly at successive quantized horizons; do not use
    # the wrapper's capacity inversion or import its choice implementation.
    for index in range(1, 100001):
        end = index*1e-7
        steps = math.ceil(end/minimum_dt)
        prediction = 1.5*slowest*steps+30.
        if prediction > 450.: break
        allowed = dict(end_time_s=end, estimated_steps=steps, forecast_wall_seconds=prediction, permitted=True)
    else: raise ValueError('short horizon unexpectedly outside bounded cost audit domain')
    require(allowed is not None, 'no affordable 100 ns window')
    return allowed


def audit_short_configuration(resolved, original, probe, probe_sha):
    require(resolved['parameters'] == original['parameters'], 'short window must not change any physical parameter')
    cfg = original['explicit_configuration']
    require(resolved['original_explicit_configuration'] == cfg, 'retained original configuration')
    override = resolved['short_window_budget_override']
    require(set(override) == {'classification','contract','choice','origin_probe_manifest_sha256','physical_parameters_unchanged'}, 'exact short override schema')
    require(override['classification'] == 'AFFORDABLE_SHORT_WINDOW_NOT_ORIGINAL_2US_COMPLETION'
        and override['contract'] == 'docs/COUPLED_AFSD_SHORT_WINDOW.md'
        and override['physical_parameters_unchanged'] is True, 'declared nonproduction short horizon')
    require(override['origin_probe_manifest_sha256'] == probe_sha, 'actual original cost receipt hash')
    require(probe['status'] == 'COMPLETE_COST_PROBE' and probe['accepted_steps'] == 5, 'complete original five-step cost evidence')
    require(probe['configuration_sha256'] == fingerprint(original), 'original cost configuration provenance')
    choice = independent_short_choice(probe)
    require(override['choice'] == choice, 'longest affordable quantized horizon, no outcome-based shortening')
    require(resolved['explicit_configuration'] == dict(cfg, registered_end_times_s=[choice['end_time_s']]), 'only scheduled endtime may differ')
    expected = dict(original, explicit_configuration=resolved['explicit_configuration'],
        original_explicit_configuration=cfg, short_window_budget_override=override)
    require(resolved == expected, 'no undeclared short-window resolved changes')
    for name, value in [('forecast_launch_limit_seconds',450), ('forecast_multiplier',1.5),
        ('forecast_overhead_seconds',30), ('internal_wall_seconds',540), ('external_wall_seconds',570)]:
        close(cfg[name], value, 0., 'unchanged cost safety '+name)
    require(cfg['registered_end_times_s'] == [2e-6,4e-6], 'original 2/4us schedule retained')
    return choice


def audit_configuration(resolved):
    require(resolved['classification'] == CLASSIFICATION, 'configuration classification')
    p = resolved['parameters']
    values = dict(dx=.0005, h_ratio=1.3, substrate_width=.004, substrate_depth=.0015,
        width=.0025, bar_length=.003, layer=.0015, tool_top=.0025, tool_radius=.0027,
        bore_width=.002496, actuator_depth=.0005, feed=.00212, traverse=.00212,
        omega=OMEGA, release_time=1., temperature=723.15, tool_temperature=733.15,
        base_temperature=703.15, rho0=2700., shear=G, poisson=.33, conductivity=167.,
        cp=900., friction=.3, contact_conductance=1000., heat_partition=.5,
        plastic_heat_fraction=1., normal_stiffness=657500., cfl=.1,
        maximum_displacement_m=.00005, room=293.15, melt=855.15,
        marker_low=.002, marker_high=.003, lift_speed=0., dwell_duration=0.,
        a1=317.6e6, a2=1e6, a3=70., b1=12.44, b2=5.65, c1=13., c2=.01,
        reference_rate=1., rate_floor=1., penalty_factor=.05, contact_reference_spacing=.0005,
        density_limit=.05, penetration_limit=.1, max_speed=20., energy_absolute_limit=1e-7,
        energy_relative_limit=.05)
    for key, expected in values.items(): close(p[key], expected, 2e-14*max(abs(expected), 1e-10), 'frozen '+key)
    require(p['spatial_operator'] == 'corrected_work_conjugate', 'unpaired spatial variant')
    if 'short_window_budget_override' in resolved:
        original = load(ORIGINAL_PROBE/'resolved_config.json'); probe = load(ORIGINAL_PROBE/'manifest.json')
        require(original['explicit_configuration'] == load(ROOT/'configs/benchmarks/coupled_afsd.json'), 'original explicit configuration differs')
        audit_short_configuration(resolved, original, probe, sha(ORIGINAL_PROBE/'manifest.json'))
        require(set(probe['source_sha256']) == set(SOURCE_PATHS), 'original exact runtime source inventory')
        for name, expected in probe['source_sha256'].items(): require(sha(ROOT/name) == expected, 'original source drift '+name)
        original_files = {f.name for f in ORIGINAL_PROBE.iterdir() if f.is_file() and f.name != 'manifest.json'}
        require(set(probe['artifact_sha256']) == original_files, 'original cost artifact inventory')
        for name, expected in probe['artifact_sha256'].items():
            require(Path(name).name == name and sha(ORIGINAL_PROBE/name) == expected, 'original cost artifact drift '+name)
    else:
        require(resolved['explicit_configuration'] == load(ROOT/'configs/benchmarks/coupled_afsd.json'), 'explicit configuration differs')
    return p


def gradient(v, i, j, g, volume, correction):
    outer = (v[j]-v[i])[:, :, None]*g[:, None, :]
    raw = np.zeros((len(v), 3, 3))
    np.add.at(raw, i, volume[j, None, None]*outer)
    np.add.at(raw, j, volume[i, None, None]*outer)
    return raw@correction


def audit_spatial(journal):
    x = np.asarray(journal['x_mid_m'], float); n = len(x)
    array(x, (n, 3), 'spatial positions')
    i = np.asarray(journal['i'], int); j = np.asarray(journal['j'], int)
    require(i.shape == j.shape and i.ndim == 1, 'pair shape')
    require(np.all((i >= 0)&(j < n)&(i < j)), 'canonical pair indices')
    keys = i*n+j
    require(len(np.unique(keys)) == len(keys), 'duplicate NNPS pair')
    distance2 = np.sum((x[:, None]-x[None, :])**2, axis=2)
    ai, aj = np.where(np.triu((distance2 < (.0013)**2)&(distance2 > 0), 1))
    require(np.array_equal(np.sort(keys), ai*n+aj), 'incomplete or non-3D NNPS pair inventory')
    delta = x[i]-x[j]; radius = np.linalg.norm(delta, axis=1); q = radius/.00065
    derivative = np.where(q < 1, 3*q*q-6, np.where(q < 2, -3*(2-q)**2, 0.))
    g = (15/(62*np.pi*.00065**4)*derivative/radius)[:, None]*delta
    close(journal['g'], g, 1e-11*max(float(np.max(np.abs(g), initial=0)), 1.), 'independent hyperbolic gradient')
    m = array(journal['mass_kg'], (n,), 'mass'); rho = array(journal['rho_before_kg_m3'], (n,), 'rho')
    close(m, np.full(n, MASS), 1e-20, 'unscaled particle mass')
    close(journal['force_density_kg_m3'], rho, 0., 'force density time layer')
    volume = m/rho; pair = (-delta)[:, :, None]*g[:, None, :]
    moment = np.zeros((n, 3, 3))
    np.add.at(moment, i, volume[j, None, None]*pair)
    np.add.at(moment, j, volume[i, None, None]*pair)
    eigen = np.linalg.eigvalsh(moment)
    require(np.all(eigen[:, 0] > 0) and np.all(eigen[:, -1] <= 1e6*eigen[:, 0]), 'non-3D/ill-conditioned neighborhood')
    B = np.linalg.inv(moment)
    close(journal['B'], B, 1e-10*max(float(np.max(np.abs(B))), 1.), 'moment inverse')
    before = array(journal['v_before_m_s'], (n, 3), 'before velocity')
    actual_work_v = array(journal['work_velocity_m_s'], (n, 3), 'actual impulse-interval work velocity')
    L = gradient(before, i, j, g, volume, B)
    Lw = gradient(actual_work_v, i, j, g, volume, B)
    close(journal['L_before_per_s'], L, 1e-9*max(float(np.max(np.abs(L))), 1.), 'actual pre-kick gradient')
    close(journal['L_midpoint_per_s'], Lw, 1e-9*max(float(np.max(np.abs(Lw))), 1.), 'actual work gradient')
    drho = -rho*np.trace(L, axis1=1, axis2=2)
    close(journal['drho_kg_m3_s'], drho, 1e-9*max(float(np.max(np.abs(drho))), 1.), 'continuity')
    sigma = array(journal['sigma_mid_Pa'], (n, 3, 3), 'midpoint stress')
    gi = np.einsum('nji,nj->ni', B[i], g); gj = np.einsum('nji,nj->ni', B[j], g)
    pairforce = (volume[i]*volume[j])[:, None]*(np.einsum('nij,nj->ni', sigma[i], gi)+np.einsum('nij,nj->ni', sigma[j], gj))
    force = np.zeros((n, 3)); np.add.at(force, i, pairforce); np.add.at(force, j, -pairforce)
    fscale = np.linalg.norm(force, axis=1).sum()
    close(journal['internal_force_N'], force, 1e-12+1e-10*fscale, 'adjoint force')
    torque = np.cross(x, force); tscale = np.linalg.norm(torque, axis=1).sum()
    require(np.linalg.norm(force.sum(axis=0)) <= 1e-12+1e-10*fscale, 'net internal force')
    require(np.linalg.norm(torque.sum(axis=0)) <= 1e-14+1e-10*tscale, 'net internal torque')
    power = float(np.sum(force*actual_work_v)); stress_power = float(np.sum(volume[:, None, None]*sigma*Lw))
    pscale = float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(actual_work_v, axis=1)))
    require(abs(power+stress_power) <= 1e-12+1e-10*pscale, 'discrete power conjugacy')
    dt = journal['dt_s']
    close(journal['increments']['material_old_work_J'], dt*np.sum(volume[:, None, None]*sigma*L), 1e-15+dt*1e-10*pscale, 'old-L work separated')
    close(journal['increments']['midpoint_stress_work_J'], dt*stress_power, 1e-15+dt*1e-10*pscale, 'true work-adjoint stress work')
    return dict(condition_max=float(np.max(eigen[:, -1]/eigen[:, 0])),
        minimum_neighborhood_rank=3, power_residual_W=power+stress_power, pair_count=len(i))


def audit_plastic(journal):
    n = len(journal['mass_kg']); dt = journal['dt_s']
    ep0 = array(journal['ep_before'], (n,), 'ep before'); ep1 = array(journal['ep_after'], (n,), 'ep after')
    dep = ep1-ep0
    require(np.all(dep >= 0), 'negative plastic increment, never clipped')
    S = array(journal['S_after_Pa'], (n, 3, 3), 'new deviatoric stress')
    close(S, S.swapaxes(1, 2), 1e-6, 'symmetric stress')
    close(np.trace(S, axis1=1, axis2=2), np.zeros(n), 1e-5, 'deviatoric trace')
    qnew = np.sqrt(1.5*np.sum(S*S, axis=(1, 2)))
    volume = np.asarray(journal['mass_kg'])/np.asarray(journal['rho_before_kg_m3'])
    plastic = qnew*dep*volume; returned = 1.5*G*dep*dep*volume
    close(journal['plastic_work_J'], plastic, 1e-15+1e-9*float(np.sum(plastic)), 'qnew dep V plastic work')
    close(journal['return_dissipation_J'], returned, 1e-15+1e-9*float(np.sum(returned)), 'separate return dissipation')
    rho_mid = np.asarray(journal['rho_before_kg_m3'])+.5*dt*np.asarray(journal['drho_kg_m3_s'])
    sigma = .5*(np.asarray(journal['S_before_Pa'])+S)-BULK*(rho_mid/2700-1)[:, None, None]*np.eye(3)
    close(journal['sigma_mid_Pa'], sigma, 1e-5, 'pressure closure and midpoint stress')
    return dep, plastic, returned


def audit_thermal(journal, plastic, friction, conduction, wall_row):
    n = len(plastic); dt = journal['dt_s']; m = np.asarray(journal['mass_kg'])
    rho = np.asarray(journal['rho_before_kg_m3']); x = np.asarray(journal['x_mid_m'])
    old = array(journal['temperature_before_K'], (n,), 'temperature before')
    new = array(journal['temperature_after_K'], (n,), 'temperature after')
    require(np.all((old >= 293.15)&(old < 855.15)&(new >= 293.15)&(new < 855.15)), 'temperature supported range')
    i, j = np.asarray(journal['i'], int), np.asarray(journal['j'], int); g = np.asarray(journal['g']).reshape(-1, 3)
    delta = x[i]-x[j]; vi = m/rho
    coefficient = vi[i]*vi[j]*334*(-np.sum(delta*g, axis=1))/np.sum(delta*delta, axis=1)
    require(np.all(coefficient >= 0), 'heat conduction must be down temperature gradient')
    pair = coefficient*(old[j]-old[i]); power = np.zeros(n); row = np.zeros(n)
    np.add.at(power, i, pair); np.add.at(power, j, -pair)
    np.add.at(row, i, coefficient); np.add.at(row, j, coefficient)
    require(abs(float(power.sum())) <= 1e-12+1e-10*float(np.sum(np.abs(pair))), 'internal heat must conserve')
    bound = .25*np.min(m*900/(row+wall_row))
    require(dt <= bound*(1+1e-12), 'thermal stability bound')
    sources = dict(plastic_work_J=plastic, plastic_heat_J=plastic,
        plastic_unthermalized_J=np.zeros(n), friction_raw_J=friction,
        friction_material_heat_J=.5*friction, friction_wall_heat_J=.5*friction,
        contact_in_J=dt*conduction, internal_conduction_J=dt*power)
    expected = plastic+.5*friction+dt*(conduction+power)
    measured = m*900*(new-old); scale = float(np.sum(np.abs(plastic)+.5*np.abs(friction)+dt*(np.abs(conduction)+np.abs(power))))
    tolerance = 1e-12+1e-9*scale
    close(measured, expected, tolerance, 'thermal mcDeltaT independent source closure')
    sources.update(total_heat_input_J=expected, measured_thermal_increment_J=measured, mc_delta_T_residual_J=measured-expected)
    saved = journal['thermal']
    for name, value in sources.items():
        close(saved['particle_energy_J'][name], value, tolerance, 'thermal particle '+name)
        close(saved['ledger'][name], value.sum(), tolerance, 'thermal ledger '+name)
    close(saved['pair_power_W'], power, 1e-12+1e-9*float(np.sum(np.abs(power))), 'independent pair heat')
    close(saved['conduction_row_W_K'], row, 1e-12+1e-9*float(row.sum()), 'conductance row')
    close(saved['ledger']['wall_heat_out_J'], np.sum(.5*friction-dt*conduction), tolerance, 'wall heat exported')
    return dict(max_mc_delta_T_residual_J=float(np.max(np.abs(measured-expected), initial=0)),
        plastic_heat_J=float(plastic.sum()), friction_material_heat_J=float(.5*friction.sum()),
        contact_in_J=float(dt*conduction.sum()), internal_conduction_net_J=float(dt*power.sum()))


def reference_faces(x, ids, t):
    """Literal finite planes; face-selection limitations remain explicit."""
    angle = OMEGA*t; c, s = np.cos(angle), np.sin(angle)
    local = (x[:, :2]-[.00212*t, 0.])@np.array([[c, -s], [s, c]])
    result = {}
    def add(index, feature, gap, normal):
        result[f'{int(ids[index])}:{feature}'] = dict(index=int(index), feature=feature,
            gap=float(gap), penetration=float(.00025-gap), normal=np.asarray(normal, float))
    for k, point in enumerate(x):
        gap = point[2]+.0015
        if gap < .00025: add(k, 0, gap, [0, 0, 1])
        q = local[k]; gap = .0015-point[2]
        if np.max(np.abs(q)) >= .002496/2 and np.linalg.norm(q) <= .0027 and -.00025 < gap < .00025:
            add(k, 1, gap, [0, 0, -1])
        if .0015 <= point[2] <= .0025:
            for axis in range(2):
                for side, sign in enumerate((-1, 1)):
                    gap = .002496/2-sign*q[axis]
                    if abs(q[1-axis]) < .002496/2+.00025 and -.0005 < gap < .00025:
                        local_n = np.zeros(2); local_n[axis] = -sign
                        normal = np.r_[np.array([[c, -s], [s, c]])@local_n, 0.]
                        add(k, 2+2*axis+side, gap, normal)
    return result


def audit_snapshot(snapshot, x, ids, t):
    faces = reference_faces(x, ids, t); rows = snapshot['rows']
    require(set(rows['contact_id']) == set(faces), 'finite contact feature inventory')
    require(len(rows['contact_id']) == len(faces), 'duplicate contact feature')
    for a, cid in enumerate(rows['contact_id']):
        f = faces[cid]; k = f['index']; feature = f['feature']
        require(rows['particle_index'][a] == k and rows['particle_id'][a] == ids[k] and rows['feature_id'][a] == feature, 'contact identity')
        close(rows['gap'][a], f['gap'], 1e-16, 'literal finite-plane gap')
        close(rows['penetration'][a], f['penetration'], 1e-16, 'literal overlap')
        close(rows['normal'][a], f['normal'], 1e-14, 'rotating normal')
        wall = np.zeros(3) if feature == 0 else np.array([.00212-OMEGA*x[k, 1], OMEGA*(x[k, 0]-.00212*t), 0.])
        close(rows['application_point'][a], x[k], 0., 'force application point')
        close(rows['wall_velocity'][a], wall, 1e-14, 'centre power-conjugate wall velocity')
    U = sum(.5*657500*f['penetration']**2 for f in faces.values())
    close(snapshot['normal_potential_J'], U, 1e-15+1e-10*U, 'sum current feature potential')
    penetration = max((f['penetration'] for f in faces.values()), default=0.)
    require(penetration <= .00005, 'contact penetration limit')
    return faces, U


def audit_interval(interval, before, mass, force, normals, fn, walls, free):
    h = interval['dt_s']; require(h > 0, 'positive contact subinterval')
    va = array(interval['v_before'], (3,), 'kick before'); vb = array(interval['v_after'], (3,), 'kick after')
    close(va, before, 1e-12, 'subinterval velocity continuity'); mid = .5*(va+vb)
    count = len(fn)
    jn = np.asarray(interval['normal_impulse'], float).reshape(count, 3)
    jt = np.asarray(interval['tangential_impulse'], float).reshape(count, 3)
    reaction = array(interval['constraint_impulse'], (3,), 'reaction')
    close(jn, h*fn[:, None]*normals, 1e-17+1e-10*h*float(np.sum(fn)), 'normal impulse')
    close(np.sum(jt*normals, axis=1), np.zeros(count), 1e-17, 'tangential impulse orthogonal')
    require(np.all(np.linalg.norm(jt, axis=1) <= .3*fn*h+1e-17), 'Coulomb impulse bound')
    momentum_roundoff = 64*np.finfo(float).eps*mass*(np.linalg.norm(va)+np.linalg.norm(vb))
    close(mass*(vb-va), h*force+jn.sum(axis=0)+jt.sum(axis=0)+reaction, 1e-17+momentum_roundoff, 'kick momentum')
    require(np.max(np.abs(reaction*free), initial=0) <= 1e-17, 'reaction on free axis')
    u = mid-walls; ut = u-np.sum(u*normals, axis=1)[:, None]*normals
    raw = -np.sum(jt*ut, axis=1)
    close(interval['friction_heat_raw_J'], raw, 1e-18+1e-10*float(np.sum(np.abs(raw))), 'signed actual interval friction')
    # The existing algorithm certifies roundoff in a translating frame. Rebuild
    # its scale from saved physical fields, never accept an arbitrary saved bound.
    frame = array(interval['reference_frame_velocity'], (3,), 'roundoff frame')
    scale = float(interval['relative_speed_scale_m_s'])
    rounding = 256*np.finfo(float).eps*np.linalg.norm(jt, axis=1)*(scale+np.linalg.norm(va-frame)+np.linalg.norm(vb-frame)+2*np.linalg.norm(walls-frame, axis=1))
    close(interval['friction_rounding_bound_J'], rounding, 1e-28+1e-9*float(rounding.sum()), 'friction rounding bound')
    require(np.all(raw >= -rounding-1e-28), 'negative raw friction beyond roundoff')
    work = dict(internal_work_J=float(h*force@mid), normal_material_work_J=float(np.sum(jn*mid)),
        normal_wall_work_J=float(np.sum(jn*walls)), wall_work_J=float(np.sum((jn+jt)*walls)),
        actuator_work_J=float(reaction@mid), delta_kinetic_J=float(.5*mass*(vb@vb-va@va)), friction_raw_J=float(raw.sum()))
    residual = work['delta_kinetic_J']-work['internal_work_J']-work['actuator_work_J']-float(np.sum((jn+jt)*mid))
    require(abs(residual) <= 1e-15+1e-9*sum(abs(v) for v in work.values()), 'contact kinetic work identity')
    for saved, key in [('other_work_J','internal_work_J'), ('constraint_work_J','actuator_work_J'), ('wall_work_J','wall_work_J'), ('delta_kinetic_J','delta_kinetic_J')]:
        close(interval[saved], work[key], 1e-15+1e-9*abs(work[key]), 'interval '+saved)
    close(interval['normal_work_J'], np.sum(jn*mid, axis=1), 1e-15, 'normal material work')
    return vb, h*mid, jn, jt, reaction, raw, work


def audit_boundary(journal, original):
    x = np.asarray(journal['x_mid_m']); ids = np.asarray(journal['particle_id'], int)
    n = len(x); dt = journal['dt_s']; t = journal['time_mid_s']
    faces, _ = audit_snapshot(journal['normal_midpoint'], x, ids, t)
    by_particle = [[] for _ in range(n)]
    for cid, value in faces.items():
        if value['penetration'] > 0: by_particle[value['index']].append(cid)
    b = journal['boundary']; traces = journal['traces']
    require(len(traces) == n and b['traces'] == traces, 'complete all-particle boundary trace')
    controlled = original[:, 2] > .0025
    require(np.count_nonzero(controlled) == 25, '25 original-material actuator particles')
    free = np.ones((n, 3)); free[controlled, 2] = 0
    commanded = np.zeros((n, 3)); commanded[controlled, 2] = -.00212
    close(b['command']['free_axes'], free, 0., 'only z constrained')
    close(b['command']['prescribed_velocity'], commanded, 0., 'actuator independent initial-position command')
    require(np.array_equal(b['command']['active'], controlled), 'actuator selection')
    close(b['command']['commanded_displacement_m'], -.00212*t, 1e-20, 'actuator command, not free trajectory')
    require(t < 1 and b['command']['phase'] == 'DRIVEN', 'registered actuator phase')
    impulses = {body: np.zeros((n, 3)) for body in ('base', 'tool', 'actuator')}
    totals = dict.fromkeys(('friction_raw_J', 'normal_material_work_J', 'normal_wall_work_J',
        'wall_work_J', 'actuator_work_J', 'internal_work_J', 'delta_kinetic_J'), 0.)
    q = np.zeros(n); flux = np.zeros(n); wallrow = np.zeros(n); used = set(); minimum = None; intervals = 0
    for trace in traces:
        k = trace['particle_index']; require(k not in used and 0 <= k < n, 'unique particle trace'); used.add(k)
        require(trace['particle_id'] == ids[k], 'trace particle ID')
        close(trace['position_m'], x[k], 0., 'trace application point')
        close(trace['free_axes'], free[k], 0., 'trace only z constrained')
        close(trace['prescribed_velocity'], commanded[k], 0., 'trace axial command')
        require(bool(trace['actuator_active']) == bool(controlled[k]), 'trace active role')
        frows = trace['contact_rows']; require(set(f['contact_id'] for f in frows) == set(by_particle[k]), 'all active faces passed to joint kick')
        normals = np.asarray([f['normal'] for f in frows], float).reshape(-1, 3)
        walls = np.asarray([f['wall_velocity'] for f in frows], float).reshape(-1, 3)
        fn = np.array([657500*faces[f['contact_id']]['penetration'] for f in frows])
        for f, force in zip(frows, fn):
            close(f['normal_force_N'], force, 1e-11, 'penalty force literal')
        temperature = journal['temperature_before_K'][k]
        reservoirs = np.array([703.15 if f['feature_id'] == 0 else 733.15 for f in frows])
        flux[k] = .00025*np.sum(reservoirs-temperature); wallrow[k] = len(frows)*.00025
        result = trace['result']; current = np.asarray(journal['v_before_m_s'][k]); work_v = np.zeros(3); elapsed = 0.
        accum_n = np.zeros((len(fn), 3)); accum_t = accum_n.copy(); accum_c = np.zeros(3)
        for a in result['intervals']:
            close(a['start_s'], elapsed, 1e-20, 'contact interval time continuity')
            frame = walls[0] if len(walls) else np.asarray(journal['v_before_m_s'][k])
            vinit = np.asarray(journal['v_before_m_s'][k])
            total_other = np.asarray(journal['internal_force_N'][k])+np.sum(fn[:, None]*normals, axis=0)
            independent_scale = max(np.linalg.norm(vinit-frame), np.linalg.norm(walls-frame) if len(walls) else 0.,
                dt*np.linalg.norm(free[k]*total_other/journal['mass_kg'][k]), np.finfo(float).tiny)
            close(a['reference_frame_velocity'], frame, 0., 'physical translating roundoff frame')
            close(a['relative_speed_scale_m_s'], independent_scale, 1e-12*max(independent_scale, 1e-12), 'independent friction roundoff scale')
            current, weighted, jn, jt, reaction, raw, works = audit_interval(a, current,
                journal['mass_kg'][k], np.asarray(journal['internal_force_N'][k]), normals, fn, walls, free[k])
            elapsed += a['dt_s']; work_v += weighted; accum_n += jn; accum_t += jt; accum_c += reaction
            q[k] += raw.sum(); intervals += 1
            if len(raw): minimum = float(raw.min()) if minimum is None else min(minimum, float(raw.min()))
            for name, value in works.items(): totals[name] += value
        close(elapsed, dt, 1e-20, 'full candidate kick covered')
        close(current, journal['v_after_m_s'][k], 1e-12, 'actual final kick velocity')
        close(work_v/dt, journal['work_velocity_m_s'][k], 1e-12, 'actual integrated work velocity')
        if controlled[k]: close(current[2], -.00212, 1e-12, 'z command preserved')
        close(np.asarray(result['normal_impulse']).reshape(-1, 3), accum_n, 1e-17, 'aggregate normal impulse')
        close(np.asarray(result['tangential_impulse']).reshape(-1, 3), accum_t, 1e-17, 'aggregate tangent impulse')
        close(result['constraint_impulse'], accum_c, 1e-17, 'aggregate actuator impulse')
        close(trace['actuator_force_on_material_N'], accum_c/dt, 1e-10, 'actuator reaction force')
        for a, f in enumerate(frows): impulses['base' if f['feature_id'] == 0 else 'tool'][k] += accum_n[a]+accum_t[a]
        impulses['actuator'][k] = accum_c
    for key in ('friction_raw_J', 'contact_power_W', 'wall_row_W_K'):
        reference = {'friction_raw_J': q, 'contact_power_W': flux, 'wall_row_W_K': wallrow}[key]
        close(b[key], reference, 1e-14, 'boundary '+key)
    totals['normal_relative_work_J'] = totals['normal_material_work_J']-totals['normal_wall_work_J']
    total_impulse = sum(impulses.values(), np.zeros((n, 3)))
    close(np.asarray(journal['mass_kg'])[:, None]*(np.asarray(journal['v_after_m_s'])-np.asarray(journal['v_before_m_s'])),
        dt*np.asarray(journal['internal_force_N'])+total_impulse, 1e-16, 'all internal/contact/actuator impulse')
    for body, values in impulses.items():
        vector = values.sum(axis=0); angular = np.cross(x, values).sum(axis=0)
        close(b['ledger'][body+'_impulse_Ns'], vector, 1e-16, 'body impulse '+body)
        close(b['ledger'][body+'_angular_impulse_Nms'], angular, 1e-18, 'body angular impulse '+body)
        for axis, a in enumerate('xyz'):
            totals[body+'_impulse_'+a+'_Ns'] = float(vector[axis])
            totals[body+'_angular_impulse_'+a+'_Nms'] = float(angular[axis])
    tool = impulses['tool']; tr = tool.sum(axis=0); torque = np.cross(x-np.array([.00212*t, 0, .0015]), tool).sum(axis=0)
    close(totals['wall_work_J'], .00212*tr[0]+OMEGA*torque[2], 1e-15+1e-9*abs(totals['wall_work_J']), 'translated-rotating tool work identity')
    for key, value in totals.items(): close(journal['increments'][key], value, 1e-14+1e-9*abs(value), 'integrated increment '+key)
    return q, flux, wallrow, totals, dict(minimum_raw_friction_J=minimum, intervals=intervals,
        features=sorted(set(f['feature'] for f in faces.values())),
        external_impulse_Ns=total_impulse.sum(axis=0), external_angular_impulse_Nms=np.cross(x, total_impulse).sum(axis=0))


def audit_record(record, original, identity):
    n = len(original)
    ids = np.asarray(record['particle_id'], int)
    require(np.array_equal(ids, identity['particle_id']), 'persistent ordered ID inventory')
    for key in ('component_id', 'material_id', 'tracer_type', 'wholebar_tag'):
        require(np.array_equal(record[key], identity[key]), 'passive identity '+key)
    close(record['actuator'], (original[:, 2] > .0025).astype(int), 0., 'original actuator metadata')
    close(record['as_rs_tag'], np.sign(-original[:, 1]), 0., 'original AS/RS metadata')
    x = array(record['position_m'], (n, 3), 'record position'); v = array(record['velocity_m_s'], (n, 3), 'record velocity')
    rho = array(record['density_kg_m3'], (n,), 'record density'); T = array(record['temperature_K'], (n,), 'record temperature')
    mass = array(record['particle_mass_kg'], (n,), 'record mass')
    close(mass, np.full(n, MASS), 1e-20, 'constant mass per particle')
    require(np.max(np.abs(rho/2700-1)) <= .05, 'accepted density limit')
    require(np.max(np.linalg.norm(v, axis=1)) <= 20, 'accepted speed limit')
    require(np.max(np.linalg.norm(x-original, axis=1)) <= .00005, 'finite-feature displacement limit')
    require(np.all((T >= 293.15)&(T < 855.15)), 'accepted temperature range')
    S = array(record['deviatoric_stress_Pa'], (n, 3, 3), 'record S'); ep = array(record['ep'], (n,), 'record ep')
    require(np.all(ep >= 0), 'nonnegative accumulated ep')
    pr = BULK*(rho/2700-1); sigma = S-pr[:, None, None]*np.eye(3)
    close(record['pressure_Pa'], pr, 1e-5, 'record pressure')
    close(record['stress_Pa'], sigma, 1e-5, 'record Cauchy stress')
    close(record['equivalent_stress_Pa'], np.sqrt(1.5*np.sum(S*S, axis=(1, 2))), 1e-5, 'record Mises')
    K = float(.5*np.sum(mass[:, None]*v*v)); P = np.sum(mass[:, None]*v, axis=0); L = np.sum(mass[:, None]*np.cross(x, v), axis=0)
    close(record['kinetic_J'], K, 1e-15, 'record kinetic')
    close(record['linear_momentum_kg_m_s'], P, 1e-16, 'record linear momentum')
    close(record['angular_momentum_kg_m2_s'], L, 1e-18, 'record angular momentum')
    close(record['mass_kg'], n*MASS, n*MASS*1e-12, 'total mass')
    close(record['thermal_increment_J'], np.sum(mass*900*(T-723.15)), 1e-12, 'thermal state increment')
    close(record['tool_origin_m'], [.00212*record['time_s'], 0., .0015], 1e-18, 'tool translation')
    close(record['tool_angle_rad'], OMEGA*record['time_s'], 1e-18, 'tool rotation')
    controlled = original[:, 2] > .0025
    close(v[controlled, 2], np.full(25, -.00212), 1e-12, 'top z control')
    close(x[controlled, 2], original[controlled, 2]-.00212*record['time_s'], 1e-15, 'actuator command stroke')
    return dict(x=x, v=v, rho=rho, T=T, S=S, ep=ep, K=K, P=P, L=L)


def audit_vtk(path, record, original, previous=None):
    data = read_point_cloud(path); scalar = data['scalars']; vector = data['vectors']
    require('SOLVED_THERMOMECHANICAL_NOT_PAPER_VALIDATED' in data['title'], 'VTK classification')
    require('internal_energy_J_per_kg' not in scalar, 'uncertified free-energy field must not be faked')
    close(data['points'], record['position_m'], 0., 'VTK position')
    close(vector['velocity_m_per_s'], record['velocity_m_s'], 0., 'VTK velocity')
    close(vector['initial_position_m'], original, 1e-18, 'VTK persistent original position')
    for name in ('particle_id', 'component_id', 'material_id', 'tracer_type', 'wholebar_tag'):
        require(scalar[name] == record[name], 'VTK identity '+name)
    pairs = dict(mass_kg='particle_mass_kg', temperature_K='temperature_K', density_kg_per_m3='density_kg_m3',
        pressure_Pa='pressure_Pa', equivalent_stress_Pa='equivalent_stress_Pa', equivalent_plastic_strain='ep')
    for name, key in pairs.items(): close(scalar[name], record[key], 1e-5 if 'stress' in name or 'pressure' in name else 1e-12, 'VTK '+name)
    close(scalar['actuator'], (original[:, 2] > .0025).astype(int), 0., 'VTK actuator metadata')
    close(scalar['as_rs_tag'], np.sign(-original[:, 1]), 0., 'AS y<0, never swapped for fit')
    close(scalar['thermal_energy_from_room_J_per_kg'], 900*(np.asarray(record['temperature_K'])-293.15), 1e-9, 'VTK solved thermal field')
    for a, b, name in ((0,0,'xx'), (1,1,'yy'), (2,2,'zz'), (0,1,'xy'), (0,2,'xz'), (1,2,'yz')):
        close(scalar['stress_'+name+'_Pa'], np.asarray(record['stress_Pa'])[:, a, b], 1e-5, 'VTK 3D stress '+name)
    if record['step']:
        close(scalar['equivalent_deviatoric_strain_rate_per_s'], record['strain_rate_per_s'], 1e-9, 'VTK solved rate')
        require(previous is not None, 'previous accepted state for actual plastic rate')
        close(scalar['equivalent_plastic_strain_rate_per_s'],
            (np.asarray(record['ep'])-np.asarray(previous['ep']))/(record['time_s']-previous['time_s']),
            1e-9, 'VTK solved plastic increment rate')
    else: require(not any('rate' in k for k in scalar), 'initial unevaluated rates not fabricated')


def provenance(directory):
    directory = Path(directory).resolve(); require(directory.is_relative_to(ROOT), 'audit only local project artifacts')
    manifest = load(directory/'manifest.json'); resolved = load(directory/'resolved_config.json')
    require(manifest['classification'] == CLASSIFICATION and manifest['production_authorized'] is False, 'nonproduction manifest')
    require(manifest['backend_commit'] == '1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5', 'upstream pin')
    short = 'short_window_budget_override' in resolved
    require(set(manifest['source_sha256']) == set(SOURCE_PATHS+(SHORT_WINDOW_SOURCES if short else [])), 'exact consumed source inventory')
    if short:
        require(manifest['mode'] == 'REGISTERED_INTEGRATED_TRIAL', 'short wrapper is not another cost probe')
        require(manifest['requested_end_s'] == resolved['short_window_budget_override']['choice']['end_time_s'],
            'short manifest request agrees with cost-selected horizon even for a blocked prefix')
    for name, expected in manifest['source_sha256'].items(): require(sha(ROOT/name) == expected, 'source hash '+name)
    require(fingerprint(resolved) == manifest['configuration_sha256'], 'resolved configuration hash')
    audit_configuration(resolved)
    actual_files = {p.name for p in directory.iterdir() if p.is_file() and p.name != 'manifest.json'}
    require(set(manifest['artifact_sha256']) == actual_files, 'exact artifact inventory')
    for name, expected in manifest['artifact_sha256'].items():
        require(Path(name).name == name and sha(directory/name) == expected, 'artifact hash '+name)
    return manifest, resolved


def audit_checkpoint(path, manifest, history, original):
    # PySPH's own NPZ schema is pickled dictionaries. Only a hash-verified file
    # under this project's run directory is read; never use this on uploads.
    with np.load(path, allow_pickle=True) as saved:
        require(int(saved['version']) == 2, 'supported upstream checkpoint schema')
        state = saved['solver_data'].item(); arrays = saved['particles'].item()['solid']['arrays']
    require(state['classification'] == CLASSIFICATION and state['accepted'] is True and state['restart_allowed'] is True, 'checkpoint accepted, never REJECTED')
    require(state['source_hashes'] == manifest['source_sha256'] and state['configuration_hash'] == manifest['configuration_sha256'], 'checkpoint provenance')
    record = history[int(state['count'])]
    close(state['t'], record['time_s'], 0., 'checkpoint time')
    require(np.array_equal(arrays['particle_id'], record['particle_id']), 'checkpoint IDs')
    for key, field in [('rho','density_kg_m3'), ('temperature','temperature_K'), ('ep','ep'), ('m','particle_mass_kg'),
                       ('plastic_work','plastic_work_J'), ('return_dissipation','return_dissipation_J')]:
        close(arrays[key], record[field], 0., 'checkpoint '+key)
    close(np.column_stack([arrays[k] for k in 'xyz']), record['position_m'], 0., 'checkpoint position')
    close(np.column_stack([arrays[k] for k in ('u','v','w')]), record['velocity_m_s'], 0., 'checkpoint velocity')
    close(np.asarray(arrays['S']).reshape(-1, 3, 3), record['deviatoric_stress_Pa'], 0., 'checkpoint stress')
    close(np.column_stack([arrays[k] for k in ('x0','y0','z0')]), original, 1e-18, 'checkpoint original material positions')
    for key in ('component_id', 'material_id', 'tracer_type', 'wholebar_tag'):
        require(np.array_equal(arrays[key], record[key]), 'checkpoint passive '+key)
    require(state['ledger'] == record['cumulative'], 'checkpoint complete ledger')
    return dict(file=Path(path).name, step=state['count'], time_s=state['t'], status='SAVED_ARRAYS_CHECKED_NOT_CONTINUATION_TEST')


def audit_case(directory):
    started = time.perf_counter(); directory = Path(directory).resolve()
    manifest, resolved = provenance(directory)
    history = load(directory/'history.json'); n = 342
    require(manifest['particle_count'] == n and len(history) == manifest['accepted_steps']+1, 'all accepted states retained')
    with (directory/'history.jsonl').open(encoding='utf-8') as stream:
        count = 0
        for count, line in enumerate(stream, 1): require(json.loads(line) == history[count-1], 'incremental/final history mismatch')
    require(count == len(history), 'incremental accepted history count')
    original = initial_geometry(); initial = history[0]
    expected_component = np.r_[np.zeros(192, int), np.ones(150, int)]
    require(np.array_equal(initial['particle_id'], np.arange(1, n+1)), '342 unique original IDs')
    require(np.array_equal(initial['component_id'], expected_component), '192 substrate plus 150 bar')
    close(initial['material_id'], np.ones(n), 0., 'all AA6061 material, tracers passive')
    outer = (expected_component == 1)&(np.max(np.abs(original[:, :2]), axis=1) >= .00125-.0005*.51)
    core = (expected_component == 1)&(np.linalg.norm(original[:, :2], axis=1) <= .0005*.51)
    whole = np.where(outer, 1, np.where(core, 2, 0))
    close(initial['wholebar_tag'], whole, 0., 'wholebar initial tags')
    close(initial['tracer_type'], np.where((original[:, 2] >= .002)&(original[:, 2] <= .003), whole, 0), 0., 'local 2-3mm tags distinct from wholebar')
    state = audit_record(initial, original, initial)
    close(state['x'], original, 1e-18, 'literal 3D initial geometry')
    v0 = np.zeros((n, 3)); v0[original[:, 2] > .0025, 2] = -.00212
    close(state['v'], v0, 0., 'initial free material at rest')
    close(state['K'], .5*25*MASS*.00212**2, 1e-24, 'initial actuator K0')
    close(state['P'], [0., 0., -25*MASS*.00212], 1e-22, 'initial actuator P0')
    close(state['S'], np.zeros((n, 3, 3)), 0., 'no manufactured initial stress')
    close(state['ep'], np.zeros(n), 0., 'no initial manufactured plasticity')
    close(state['T'], np.full(n, 723.15), 0., 'declared hot start, not solved preheat')
    sums = {key: 0. for key in initial['cumulative']}; P0 = state['P']; L0 = state['L']; K0 = state['K']
    stats = dict(accepted_steps=0, active_plastic_steps=0, longest_consecutive_plastic_steps=0,
        maximum_condition_number=0., maximum_spatial_power_residual_W=0., maximum_thermal_particle_residual_J=0.,
        minimum_neighbor_rank=None, contact_subintervals=0, minimum_raw_friction_J=None, energy_warning_frames=0)
    streak = 0; features = set(); free_bar_plastic = set(); substrate_plastic = set(); substep_thermal = []
    external = np.zeros(3); angular = np.zeros(3); first_plastic = None
    first_features = {key: [0., 0.] for key in reference_faces(original, initial['particle_id'], 0.)}; last_features = set()
    thermal_alias = dict(plastic_heat_J='plastic_heat_J', nonthermal_plastic_J='plastic_unthermalized_J',
        friction_heat_J='friction_material_heat_J', friction_heat_to_wall_J='friction_wall_heat_J',
        contact_heat_J='contact_in_J', wall_heat_out_J='wall_heat_out_J',
        internal_conduction_net_J='internal_conduction_J', thermal_balance_residual_J='mc_delta_T_residual_J')
    with gzip.open(directory/'steps.jsonl.gz', 'rt', encoding='utf-8') as stream:
        for step, line in enumerate(stream, 1):
            journal = json.loads(line); require(step < len(history), 'journal cannot exceed accepted history')
            record = history[step]; previous = history[step-1]
            require(journal['stage'] == 'ACCEPTED_COUPLED_DRIFT_KICK_DRIFT' and journal['step'] == step == record['step'], 'journal accepted status and step')
            require(np.array_equal(journal['particle_id'], initial['particle_id']), 'journal persistent IDs')
            dt = journal['dt_s']; require(dt > 0, 'positive physical dt')
            close(record['time_s'], previous['time_s']+dt, 1e-20, 'physical time progression')
            close(journal['time_mid_s'], previous['time_s']+.5*dt, 1e-20, 'force evaluation time')
            close(journal['accepted_time_s'], record['time_s'], 0., 'accepted endpoint time')
            close(journal['v_before_m_s'], state['v'], 0., 'previous accepted velocity')
            close(journal['x_before_m'], state['x'], 0., 'previous accepted position')
            close(journal['x_mid_m'], state['x']+.5*dt*state['v'], 1e-18, 'true first half drift')
            close(journal['rho_before_kg_m3'], state['rho'], 0., 'previous accepted density')
            close(journal['temperature_before_K'], state['T'], 0., 'previous accepted temperature')
            close(journal['S_before_Pa'], state['S'], 0., 'previous accepted stress')
            close(journal['ep_before'], state['ep'], 0., 'previous accepted plastic strain')
            sound = np.sqrt((BULK+4*G/3)/2700)
            acoustic = .1*.00065/(sound+np.linalg.norm(state['v'], axis=1).max())
            require(dt <= min(acoustic, .1*np.sqrt(MASS/657500))*(1+1e-12), 'actual acoustic/normal stable timestep')
            spatial = audit_spatial(journal); dep, plastic, returned = audit_plastic(journal)
            q, contact, row, mechanical, contact_stats = audit_boundary(journal, original)
            thermal = audit_thermal(journal, plastic, q, contact, row); substep_thermal.append(thermal)
            next_state = audit_record(record, original, initial)
            raw_L = np.asarray(journal['L_before_per_s']); D = .5*(raw_L+raw_L.swapaxes(1,2))
            dev = D-np.trace(D, axis1=1, axis2=2)[:, None, None]*np.eye(3)/3
            close(record['strain_rate_per_s'], np.sqrt(2/3*np.sum(dev*dev, axis=(1,2))), 1e-9, 'actual 3D strain-rate invariant')
            close(next_state['v'], journal['v_after_m_s'], 1e-12, 'accepted endpoint equals actual kick')
            close(next_state['x'], np.asarray(journal['x_mid_m'])+.5*dt*next_state['v'], 1e-18, 'true second half drift')
            close(next_state['x'], journal['x_after_m'], 0., 'saved endpoint geometry')
            close(next_state['rho'], state['rho']+dt*np.asarray(journal['drho_kg_m3_s']), 1e-10, 'actual continuity endpoint')
            close(next_state['rho'], journal['rho_after_kg_m3'], 0., 'saved endpoint density')
            for field, key in [('T','temperature_after_K'), ('S','S_after_Pa'), ('ep','ep_after')]:
                close(next_state[field], journal[key], 1e-6 if field == 'S' else 1e-12, 'accepted '+field+' equals candidate')
            close(record['plastic_work_J'], np.asarray(previous['plastic_work_J'])+plastic, 1e-15, 'stored accumulated plastic work')
            close(record['return_dissipation_J'], np.asarray(previous['return_dissipation_J'])+returned, 1e-15, 'stored separate return work')
            _, ub = audit_snapshot(journal['normal_before'], state['x'], np.asarray(initial['particle_id']), previous['time_s'])
            fa, ua = audit_snapshot(journal['normal_after'], next_state['x'], np.asarray(initial['particle_id']), record['time_s'])
            mechanical['normal_temporal_residual_J'] = ua-ub+mechanical['normal_relative_work_J']
            close(record['contact_potential_J'], ua, 1e-15, 'accepted normal potential')
            for cid, face in fa.items():
                if cid not in first_features: first_features[cid] = [previous['time_s'], record['time_s']]
            last_features = set(fa)
            mechanical['plastic_work_J'] = float(plastic.sum()); mechanical['numerical_return_dissipation_J'] = float(returned.sum())
            for key, savedkey in thermal_alias.items(): mechanical[key] = journal['thermal']['ledger'][savedkey]
            mechanical['material_old_work_J'] = journal['increments']['material_old_work_J']
            mechanical['midpoint_stress_work_J'] = journal['increments']['midpoint_stress_work_J']
            mechanical['external_mass_kg'] = 0.; mechanical['external_feed_enthalpy_J'] = 0.
            require(set(mechanical) == set(sums) == set(journal['increments']), 'complete independent ledger inventory')
            for key, value in mechanical.items():
                close(journal['increments'][key], value, 1e-12+1e-9*abs(value), 'independent increment '+key)
                sums[key] += value
                close(record['cumulative'][key], sums[key], 1e-12+1e-9*abs(sums[key]), 'accumulated raw ledger '+key)
            external += contact_stats['external_impulse_Ns']; angular += contact_stats['external_angular_impulse_Nms']
            close(next_state['P']-P0, external, 1e-15, 'global external linear impulse')
            close(next_state['L']-L0, angular, 1e-17, 'global external angular impulse')
            kinetic_identity = next_state['K']-K0+sums['friction_raw_J']-sums['wall_work_J']-sums['actuator_work_J']-sums['internal_work_J']-sums['normal_relative_work_J']
            require(abs(kinetic_identity) <= 1e-15+1e-9*sum(abs(sums[k]) for k in ('friction_raw_J','wall_work_J','actuator_work_J','internal_work_J','normal_relative_work_J')), 'global K/contact/constraint identity')
            rho = next_state['rho']; S = next_state['S']
            shear_energy = np.sum(MASS/rho*np.sum(S*S, axis=(1, 2))/(4*G))
            ratio = rho/2700; volume_energy = np.sum(MASS*BULK/2700*(np.log(ratio)+1/ratio-1))
            energy = next_state['K']+shear_energy+volume_energy+ua+np.sum(MASS*900*(next_state['T']-723.15))
            residual = energy-initial['total_J']-sums['wall_work_J']-sums['actuator_work_J']+sums['wall_heat_out_J']
            close(record['energy_residual_J'], residual, 1e-12, 'diagnostic residual, not physical heat')
            scale = max(abs(initial['total_J']), abs(sums['wall_work_J'])+abs(sums['actuator_work_J'])+abs(sums['wall_heat_out_J']), 1e-15)
            warning = abs(residual) > 1e-7+.05*scale
            require((record['energy_status'] == 'ENERGY_WARNING') == warning, 'honest diagnostic warning')
            stats['energy_warning_frames'] += int(warning)
            active = bool(np.any(dep > 0)); streak = streak+1 if active else 0
            stats['active_plastic_steps'] += int(active); stats['longest_consecutive_plastic_steps'] = max(stats['longest_consecutive_plastic_steps'], streak)
            if active and first_plastic is None: first_plastic = dict(step=step, time_s=record['time_s'])
            ids = np.asarray(initial['particle_id']); free_bar = (expected_component == 1)&(original[:, 2] <= .0025)
            free_bar_plastic.update(int(a) for a in ids[(dep > 0)&free_bar]); substrate_plastic.update(int(a) for a in ids[(dep > 0)&(expected_component == 0)])
            features.update(contact_stats['features']); stats['contact_subintervals'] += contact_stats['intervals']
            stats['maximum_condition_number'] = max(stats['maximum_condition_number'], spatial['condition_max'])
            stats['maximum_spatial_power_residual_W'] = max(stats['maximum_spatial_power_residual_W'], abs(spatial['power_residual_W']))
            stats['maximum_thermal_particle_residual_J'] = max(stats['maximum_thermal_particle_residual_J'], thermal['max_mc_delta_T_residual_J'])
            stats['minimum_neighbor_rank'] = 3
            minimum = contact_stats['minimum_raw_friction_J']
            if minimum is not None: stats['minimum_raw_friction_J'] = minimum if stats['minimum_raw_friction_J'] is None else min(minimum, stats['minimum_raw_friction_J'])
            stats['accepted_steps'] = step; state = next_state
    require(stats['accepted_steps'] == manifest['accepted_steps'], 'all accepted steps independently audited')
    require(manifest['final_record'] == history[-1], 'manifest final accepted record')
    close(manifest['accepted_time_s'], history[-1]['time_s'], 0., 'manifest accepted time')
    complete = manifest['status'].startswith('COMPLETE')
    short = 'short_window_budget_override' in resolved
    if complete and manifest['mode'] == 'COST_PROBE':
        require(manifest['status'] == 'COMPLETE_COST_PROBE' and manifest['accepted_steps'] == 5 and manifest['requested_end_s'] is None,
            'cost probe cannot masquerade as registered full interval')
    if complete and manifest['mode'] == 'REGISTERED_INTEGRATED_TRIAL':
        allowed = [resolved['short_window_budget_override']['choice']['end_time_s']] if short else (2e-6, 4e-6)
        require(manifest['requested_end_s'] in allowed, 'explicitly registered endtime')
        close(manifest['accepted_time_s'], manifest['requested_end_s'], 1e-20, 'common registered endpoint reached')
    if not complete:
        require(manifest['status'] in ('BLOCKED', 'BLOCKED_STEP_CAP'), 'known incomplete status')
        if (directory/'REJECTED_candidate.json').exists():
            rejected = load(directory/'REJECTED_candidate.json')
            require(rejected['classification'] == 'REJECTED' and rejected['restart_allowed'] is False and rejected['accepted'] is False, 'rejected state must never be endpoint')
            require(rejected['step'] == manifest['accepted_steps']+1, 'first rejected candidate follows accepted prefix')
            require(rejected['source_hashes'] == manifest['source_sha256'] and rejected['configuration_hash'] == manifest['configuration_sha256'], 'rejected provenance')
    frames = load(directory/'series.json')['frames']; expected_frames = sorted(set([0, manifest['accepted_steps'], *range(10, manifest['accepted_steps']+1, 10)]))
    require([f['step'] for f in frames] == expected_frames, 'all promised VTK cadence and final accepted frame')
    for frame in frames:
        require(sha(directory/frame['file']) == frame['sha256'], 'VTK hash')
        audit_vtk(directory/frame['file'], history[frame['step']], original, history[frame['step']-1] if frame['step'] else None)
    checkpoints = [audit_checkpoint(path, manifest, history, original) for path in sorted(directory.glob('*accepted.npz'))]
    require(any(q['file'] == 'last_accepted.npz' and q['step'] == manifest['accepted_steps'] for q in checkpoints), 'last accepted checkpoint')
    free_bar = (expected_component == 1)&(original[:, 2] <= .0025); substrate = expected_component == 0
    d = np.linalg.norm(state['x']-original, axis=1)
    final_ep = state['ep']; epsummary = {name: dict(nonzero_count=int(np.count_nonzero(final_ep[mask])), maximum=float(np.max(final_ep[mask], initial=0)), ids=[int(i) for i in np.asarray(initial['particle_id'])[mask&(final_ep > 0)]]) for name, mask in [('free_bar',free_bar), ('substrate',substrate), ('actuator',~(free_bar|substrate))]}
    thermal_change = float(np.max(np.abs(state['T']-723.15)))
    mechanisms = dict(rotating_open_bore_and_axial_feed='ACTIVE' if features-{0} and sums['wall_work_J'] != 0 else 'CODE_CONNECTED_NOT_YET_OBSERVED',
        deformable_substrate_and_feedstock='ACTIVE' if np.max(d[substrate]) > 0 and np.max(d[free_bar]) > 0 else 'CODE_CONNECTED_NOT_YET_OBSERVED',
        plasticity='ACTIVE_FINITE_WINDOW' if stats['active_plastic_steps'] else 'IMPLEMENTED_NO_PLASTIC_INCREMENT_OBSERVED',
        thermal='ACTIVE_SOLVED' if thermal_change > 0 else 'IMPLEMENTED_NO_THERMAL_CHANGE_OBSERVED',
        sustained_deposition='NOT_DEMONSTRATED', three_case_production='NOT_RUN', spatial_time_convergence='NOT_RUN')
    readback_status = ('PASS_COST_PROBE_READBACK_NOT_REGISTERED_TRIAL' if manifest['mode'] == 'COST_PROBE'
        else 'PASS_SHORT_WINDOW_NOT_2US_NOT_DEPOSITION' if short
        else 'PASS_REGISTERED_INTERVAL_READBACK_NOT_SUSTAINED_DEPOSITION') if complete else 'PASS_ACCEPTED_PREFIX_ONLY_PHYSICAL_RUN_BLOCKED'
    return dict(status=readback_status,
        mode=manifest['mode'], physical_run_status=manifest['status'], classification=CLASSIFICATION, production_authorized=False,
        short_window_override=resolved.get('short_window_budget_override'),
        original_2us_4us_status='NOT_RUN_COST_EXCEEDS_BUDGET' if short else 'REFER_TO_PREREGISTERED_COST_EVIDENCE',
        case=str(directory), source_sha256=manifest['source_sha256'], manifest_sha256=sha(directory/'manifest.json'),
        reader_sha256=sha(Path(__file__)), accepted_end_s=manifest['accepted_time_s'], requested_end_s=manifest['requested_end_s'],
        particle_count=n, substrate_particles=192, bar_particles=150, axial_controlled_particles=25,
        independent_accepted_states=len(history), vtk_frames_checked=len(frames), checkpoint_arrays=checkpoints,
        statistics=stats, first_plastic=first_plastic, free_bar_plastic_particle_ids=sorted(free_bar_plastic),
        free_substrate_plastic_particle_ids=sorted(substrate_plastic), ep_distribution=epsummary,
        feature_ids_observed=sorted(features), first_sampled_contact_presence_brackets_s=first_features,
        final_contact_ids=sorted(last_features), cumulative=sums, mechanisms=mechanisms,
        substrate_displacement_max_m=float(np.max(d[substrate])), free_bar_displacement_max_m=float(np.max(d[free_bar])),
        maximum_temperature_change_K=thermal_change, final_temperature_range_K=[float(state['T'].min()),float(state['T'].max())],
        final_stress_xz_max_Pa=float(np.max(np.abs(state['S'][:,0,2]))), final_stress_yz_max_Pa=float(np.max(np.abs(state['S'][:,1,2]))),
        total_energy_certified=False, total_energy_note='Es is a hypoelastic diagnostic; warnings retained, return/residual never added as heat',
        limitations=['No exact full-trajectory solution','No continuous inlet, bonding evolution or sustained deposition certification',
            'No new spatial/time convergence test','Contact finite edges/exterior remain incomplete',
            'Checkpoint saved arrays verified; restart execution is a separate test'],
        reader_wall_seconds=time.perf_counter()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); require(not args.output.exists(), 'fresh output required')
    try: result = audit_case(args.case); code = 0
    except Exception as error:
        result = dict(status='FAIL_COUPLED_READBACK', error=str(error), reader_sha256=sha(Path(__file__)), production_authorized=False); code = 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, allow_nan=False)); return code


if __name__ == '__main__': raise SystemExit(main())

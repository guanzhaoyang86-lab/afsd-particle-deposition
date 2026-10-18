"""Bounded, detached HHT / current-geometry elastic FKS / BE benchmark stage.

Physical masses, paired SPH operators and the existing material/heat routines
are reused. This is NOT an AFSD continuation schema or a contact-event driver.
Callbacks supply noncontact applied force and separately identified contacts.
Interior events reject; the caller must locate/split and explicitly restart.
"""
import copy
import hashlib
import time
import numpy as np

from . import hht_path
from .material import FKS, update
from .mechanics import positions, velocities, corrected_operators, internal_force, pressure
from .implicit_thermal import advance_implicit_thermal
from .structural_mobility import build_frozen_elastic
from .updated_trial import particle_hash

CLASSIFICATION = 'PROCESS_HHT_EVENT_THERMAL_APPROXIMATE_COUPLED_ELASTIC_BENCHMARK'
SCHEMA = 'HHT_CURRENT_GEOMETRY_ELASTIC_BE_STATE_V1'
FORCE_DEFINITION = 'NONCONTACT_APPLIED_FORCE_ONLY'
MAX_EVALUATIONS = 24


class CoupledStageError(RuntimeError):
    """Rejected private candidate; the accepted state remains unchanged."""
    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = {} if diagnostics is None else diagnostics


def _hash(value):
    digest = hashlib.sha256()
    def visit(v):
        if isinstance(v, dict):
            digest.update(b'dict')
            for key in sorted(v):
                visit(str(key)); visit(v[key])
        elif isinstance(v, np.ndarray):
            digest.update(b'array'+v.dtype.str.encode()+str(v.shape).encode()+v.tobytes())
        elif isinstance(v, (list, tuple)):
            digest.update(b'list'+str(len(v)).encode())
            for item in v: visit(item)
        elif isinstance(v, np.generic): visit(v.item())
        elif isinstance(v, (str, int, float, bool, type(None))):
            digest.update((type(v).__name__+':'+repr(v)+';').encode())
        else:
            raise ValueError('unsupported checkpoint/parameter value '+type(v).__name__)
    visit(value)
    return digest.hexdigest()


def _detach(pa):
    return pa.extract_particles(np.arange(len(pa.x), dtype=np.uint32))


def _array(value, shape, name):
    result = np.array(value, dtype=float, copy=True)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError('finite matching '+name+' required')
    return result


def _parameters(p):
    for key in ('rho0', 'shear', 'bulk', 'cp', 'conductivity', 'room', 'melt'):
        if not np.isfinite(p[key]) or p[key] <= 0:
            raise ValueError('positive finite '+key+' required')
    for key in ('plastic_heat_fraction', 'heat_partition'):
        if not np.isfinite(p[key]) or not 0 <= p[key] <= 1:
            raise ValueError('heat fraction outside [0,1]')
    if p['melt'] <= p['room']:
        raise ValueError('ordered solid temperature limits required')
    if not np.isfinite(p['friction']) or p['friction'] < 0:
        raise ValueError('nonnegative friction required')
    weight = float(p.get('hht_weight', .9))
    gamma, beta = hht_path.parameters(weight)
    direction_limit = float(p.get('contact_coulomb_relative_tolerance', 1e-3))
    if not np.isfinite(direction_limit) or not 0 <= direction_limit <= .01:
        raise ValueError('registered Coulomb quadrature defect bound in [0,.01] required')
    return weight, gamma, beta, direction_limit


def _valid_particles(pa, p):
    if len(pa.x) == 0 or len(np.unique(pa.particle_id)) != len(pa.x):
        raise ValueError('nonempty unique particle inventory required')
    for key, prop in pa.properties.items():
        if not np.isfinite(prop.get_npy_array()).all():
            raise ValueError('nonfinite particle property '+key)
    if np.any(pa.m <= 0) or np.any(pa.rho <= 0) or np.any(pa.h <= 0) or np.any(pa.ep < 0):
        raise ValueError('invalid physical particle state')
    if np.any(pa.temperature < p['room']) or np.any(pa.temperature >= p['melt']):
        raise ValueError('particle outside supported solid temperature domain')


def _boundary(pa, callback, at, p):
    private = _detach(pa)
    before = particle_hash(private)
    raw = callback(private, float(at))
    if particle_hash(private) != before:
        raise ValueError('boundary callback mutated its particle argument')
    if raw.get('force_definition') != FORCE_DEFINITION:
        raise ValueError('explicit NONCONTACT_APPLIED_FORCE_ONLY boundary definition required')
    count = len(pa.x)
    applied = _array(raw['force_N'], (count, 3), 'noncontact applied force')
    total = applied.copy(); contacts = []
    seen = set()
    for source in raw.get('contacts', []):
        row = copy.deepcopy(source)
        index = int(row['particle_index']); cid = str(row['contact_id'])
        if not 0 <= index < count or cid in seen:
            raise ValueError('invalid/duplicate contact identity')
        seen.add(cid)
        row['particle_index'] = index; row['contact_id'] = cid
        row['particle_id'] = int(pa.particle_id[index])
        for key in ('normal', 'tangential_force_N', 'wall_velocity_m_s'):
            row[key] = _array(row[key], (3,), key)
        if abs(np.linalg.norm(row['normal'])-1) > 1e-12:
            raise ValueError('unit contact normal required')
        for key in ('normal_force_N', 'conductance_W_K', 'penetration_m', 'normal_potential_J',
                    'normal_stiffness_N_m', 'wall_temperature_K'):
            row[key] = float(row[key])
            if not np.isfinite(row[key]) or row[key] < 0:
                raise ValueError('finite nonnegative contact '+key+' required')
        if row['normal_stiffness_N_m'] <= 0 or row['wall_temperature_K'] <= 0:
            raise ValueError('positive contact stiffness/reservoir temperature required')
        kn = row['normal_stiffness_N_m']; pen = row['penetration_m']
        if abs(row['normal_force_N']-kn*pen) > 1e-14+1e-12*abs(kn*pen):
            raise ValueError('contact normal penalty force inconsistent with penetration')
        if abs(row['normal_potential_J']-.5*kn*pen*pen) > 1e-20+1e-12*.5*kn*pen*pen:
            raise ValueError('contact normal potential inconsistent with penetration')
        cap = p['friction']*row['normal_force_N']; f = row['tangential_force_N']
        scale = max(cap, np.linalg.norm(f), np.finfo(float).tiny)
        if abs(f@row['normal']) > 256*np.finfo(float).eps*scale or np.linalg.norm(f) > cap+256*np.finfo(float).eps*scale:
            raise ValueError('contact endpoint force outside tangent friction cone')
        row['capacity_N'] = cap
        total[index] += row['normal_force_N']*row['normal']+f
        contacts.append(row)
    contacts.sort(key=lambda row: row['contact_id'])
    return dict(time_s=float(at), force_definition=FORCE_DEFINITION,
        noncontact_force_N=applied, total_boundary_force_N=total, contacts=contacts)


def _force(pa, p):
    i, j, g, gradient, _, correction = corrected_operators(pa)
    sigma = pa.S.reshape(-1, 3, 3)-pressure(pa.rho, p['rho0'], p['bulk'])[:, None, None]*np.eye(3)
    force = internal_force(pa, sigma, i, j, g, correction=correction)
    volume = pa.m/pa.rho
    material_power = float(np.sum(volume[:, None, None]*sigma*gradient))
    kinetic_power = float(np.sum(force*velocities(pa)))
    scale = float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(velocities(pa), axis=1)))
    defect = kinetic_power+material_power
    if abs(defect) > 1e-12+1e-10*scale:
        raise ValueError('current-geometry paired spatial power failed')
    return force, dict(power_conjugacy_residual_W=defect, power_absolute_scale_W=scale,
        internal_power_W=kinetic_power, stress_power_W=material_power,
        internal_net_force_N=force.sum(axis=0),
        internal_net_torque_Nm=np.cross(positions(pa), force).sum(axis=0), pair_count=len(i))


def _snapshot(pa):
    return {key: prop.get_npy_array().copy() for key, prop in pa.properties.items()}


def _state(pa, p, at, acc, force, boundary, parent=None, dt=None, event=None):
    value = dict(schema=SCHEMA, classification=CLASSIFICATION, time_s=float(at),
        acceleration_m_s2=np.array(acc, copy=True), force_N=np.array(force, copy=True),
        boundary_snapshot=copy.deepcopy(boundary), particle_id=pa.particle_id.copy(),
        particle_hash=particle_hash(pa), parameter_hash=_hash(p),
        particle_snapshot=_snapshot(pa), parent_state_hash=parent, last_dt_s=dt,
        last_event=copy.deepcopy(event), production_authorized=False,
        material_energy_certified=False, restart_allowed=True)
    value['state_hash'] = _hash(value)
    return value


def _validate_state(pa, state, p):
    _parameters(p); _valid_particles(pa, p)
    if state.get('schema') != SCHEMA or state.get('classification') != CLASSIFICATION:
        raise ValueError('registered HHT coupled state schema required')
    payload = {key: value for key, value in state.items() if key != 'state_hash'}
    if state.get('state_hash') != _hash(payload):
        raise ValueError('HHT accepted semantic checkpoint hash mismatch')
    if state['particle_hash'] != particle_hash(pa) or state['parameter_hash'] != _hash(p):
        raise ValueError('accepted particle/parameter history changed')
    if not np.array_equal(state['particle_id'], pa.particle_id):
        raise ValueError('accepted particle ID correspondence lost')
    for key, values in state['particle_snapshot'].items():
        if key not in pa.properties or not np.array_equal(values, pa.properties[key].get_npy_array()):
            raise ValueError('full physical checkpoint differs from supplied particles')


def initialize(pa, p, boundary, time=0.):
    """Initialize actual force/physical acceleration without changing the state."""
    _parameters(p); _valid_particles(pa, p)
    if not np.isfinite(time) or time < 0:
        raise ValueError('finite nonnegative initial time required')
    before = particle_hash(pa); private = _detach(pa)
    internal, _ = _force(private, p)
    contact = _boundary(private, boundary, time, p)
    force = internal+contact['total_boundary_force_N']
    if particle_hash(pa) != before:
        raise RuntimeError('initialization mutated input')
    return _state(pa, p, time, force/pa.m[:, None], force, contact)


def restore_checkpoint(template_pa, state, p):
    """Restore every property using a matching ParticleArray layout, then verify."""
    restored = _detach(template_pa)
    if set(restored.properties) != set(state['particle_snapshot']):
        raise ValueError('checkpoint/template property inventory mismatch')
    for key, values in state['particle_snapshot'].items():
        target = restored.properties[key].get_npy_array()
        if target.shape != values.shape or target.dtype != values.dtype:
            raise ValueError('checkpoint/template property layout mismatch')
        target[:] = values
    copied = copy.deepcopy(state)
    _validate_state(restored, copied, p)
    return restored, copied


def restart_at_control_event(pa, state, p, new_boundary):
    """Explicit right-limit force/acceleration reset; no physical field reset."""
    _validate_state(pa, state, p)
    before = particle_hash(pa)
    internal, spatial = _force(_detach(pa), p)
    boundary = _boundary(pa, new_boundary, state['time_s'], p)
    force = internal+boundary['total_boundary_force_N']; acc = force/pa.m[:, None]
    weight, gamma, beta, _ = _parameters(p)
    delta = beta-gamma/2; last_dt = state['last_dt_s']
    storage = (None if last_dt is None else .5*last_dt**2*delta*float(np.sum(
        pa.m[:, None]*(acc*acc-state['acceleration_m_s2']**2))))
    event = dict(kind='EXPLICIT_RIGHT_LIMIT_CONTROL_FORCE_ACCELERATION_INITIALIZATION',
        time_s=state['time_s'], old_acceleration_m_s2=state['acceleration_m_s2'].copy(),
        new_acceleration_m_s2=acc.copy(), old_force_N=state['force_N'].copy(), new_force_N=force.copy(),
        old_boundary_snapshot=copy.deepcopy(state['boundary_snapshot']),
        new_boundary_snapshot=copy.deepcopy(boundary), physical_state_hash_before=before,
        physical_state_hash_after=particle_hash(pa), physical_energy_jump_J=0.,
        last_dt_acceleration_storage_jump_J=storage,
        acceleration_storage_definition='0.5*h_last^2*(beta-gamma/2)*sum(m*A^2)',
        full_algorithmic_energy_jump_certified=False, numerical_terms_are_material_heat=False,
        spatial=spatial)
    result = _state(pa, p, state['time_s'], acc, force, boundary,
        parent=state['state_hash'], dt=last_dt, event=event)
    return result, copy.deepcopy(event)


def _contact_work(old, new, pa, v1, acc0, acc1, dt, weight, direction_limit):
    left = old['contacts']; right = new['contacts']
    if [row['contact_id'] for row in left] != [row['contact_id'] for row in right]:
        raise CoupledStageError('contact identity switch requires a located event')
    count = len(pa.x); q = np.zeros(count); bounds = np.zeros(count)
    wall = np.zeros(count); bath = np.zeros(count); audits = []
    physical_wall_work = 0.; normal_relative_work = 0.; normal_material_work = 0.
    v0 = velocities(pa)
    for row0, row1 in zip(left, right):
        index = row0['particle_index']
        if (row1['particle_index'] != index or row1['particle_id'] != row0['particle_id'] or
                not np.array_equal(row1['normal'], row0['normal']) or
                row1['normal_stiffness_N_m'] != row0['normal_stiffness_N_m']):
            raise CoupledStageError('fixed-normal/stable-contact scope violated')
        audit = hht_path.contact_audit(row0['normal'], row0['tangential_force_N'],
            row1['tangential_force_N'], row0['wall_velocity_m_s'], row1['wall_velocity_m_s'],
            v0[index], v1[index], acc0[index], acc1[index], dt,
            row0['capacity_N'], row1['capacity_N'], weight)
        audit.update(contact_id=row0['contact_id'], particle_id=row0['particle_id'], particle_index=index)
        if not audit['whole_path_passive'] or audit['requires_event_split']:
            raise CoupledStageError('contact path needs passivity/event split', dict(contact=audit))
        # Registered quadrature-law error, separate from the exact power guard.
        fractions = sorted(set([0., .25, .5, .75, 1.]+
            [row['fraction'] for row in audit['quadrature']]+audit['power_stationary_fractions']))
        defects = []
        for fraction in fractions:
            u = audit['slip_coefficients_m_s']
            slip = u[0]+fraction*u[1]+fraction*fraction*u[2]
            speed = np.linalg.norm(slip)
            force = row0['tangential_force_N']+fraction*(row1['tangential_force_N']-row0['tangential_force_N'])
            cap = row0['capacity_N']+fraction*(row1['capacity_N']-row0['capacity_N'])
            defect = (0. if speed <= audit['zero_speed_bound_m_s'] else
                float(np.linalg.norm(force+cap*slip/speed)/max(cap, np.finfo(float).tiny)))
            defects.append(defect)
        audit['sampled_coulomb_defect_fractions'] = fractions
        audit['sampled_coulomb_relative_defects'] = defects
        audit['maximum_sampled_coulomb_relative_defect'] = max(defects)
        audit['coulomb_relative_error_bound'] = direction_limit
        audit['interior_coulomb_direction'] = 'SAMPLED_QUADRATURE_DEFECT_APPROXIMATE_NOT_CONTINUOUS_LAW_PROOF'
        if max(defects) > direction_limit+512*np.finfo(float).eps:
            raise CoupledStageError('contact turning-traction quadrature law defect', dict(contact=audit))
        if any(row['friction_heat_raw_J'] < -.5*audit['friction_rounding_bound_J']
               for row in audit['quadrature']):
            raise CoupledStageError('negative individual Gauss friction work', dict(contact=audit))
        q[index] += audit['friction_heat_raw_J']; bounds[index] += audit['friction_rounding_bound_J']
        wall[index] += row1['conductance_W_K']
        bath[index] += row1['conductance_W_K']*row1['wall_temperature_K']
        n0 = row0['normal_force_N']*row0['normal']; n1 = row1['normal_force_N']*row1['normal']
        normal = hht_path.linear_work(n0, n1, v0[index], v1[index], acc0[index], acc1[index], dt, weight)
        normal_material_work += float(normal['work_J'])
        normal_wall = 0.; contact_wall = 0.
        for quadrature in normal['quadrature']:
            fraction = quadrature['fraction']
            wallv = row0['wall_velocity_m_s']+fraction*(row1['wall_velocity_m_s']-row0['wall_velocity_m_s'])
            ft = row0['tangential_force_N']+fraction*(row1['tangential_force_N']-row0['tangential_force_N'])
            normal_wall += float(quadrature['impulse_Ns']@wallv)
            contact_wall += float((quadrature['impulse_Ns']+.5*dt*ft)@wallv)
        physical_wall_work += contact_wall
        normal_relative_work += float(normal['work_J'])-normal_wall
        audits.append(audit)
    return dict(friction_raw_J=q, friction_rounding_bound_J=bounds,
        wall_row_W_K=wall, wall_weighted_temperature_W=bath, contacts=audits,
        wall_work_J=physical_wall_work, normal_relative_work_J=normal_relative_work,
        normal_material_work_J=normal_material_work,
        normal_potential_change_J=sum(row['normal_potential_J'] for row in right)-sum(row['normal_potential_J'] for row in left))


def density_pressure_correction(acceleration_matrix, pressure_map, density_response,
                                force_residual, density_residual, weight):
    """Exact solve of a FROZEN numerical block, never a new physical map.

    A da - weight*C drho = -RA; -L da + drho = -Rrho.
    Keeping the pressure residual RHS and L*da density response is necessary:
    bulk stiffness in A's Schur complement alone is not this block solve.
    """
    A = np.asarray(acceleration_matrix, float); C = np.asarray(pressure_map, float)
    L = np.asarray(density_response, float); RA = np.asarray(force_residual, float)
    Rrho = np.asarray(density_residual, float)
    na, nr = RA.size, Rrho.size
    if (RA.shape != (na,) or Rrho.shape != (nr,) or A.shape != (na, na) or
            C.shape != (na, nr) or L.shape != (nr, na) or
            not all(np.isfinite(v).all() for v in (A, C, L, RA, Rrho))):
        raise ValueError('finite compatible density-pressure block required')
    matrix = A-weight*C@L
    da = np.linalg.solve(matrix, -RA-weight*C@Rrho)
    drho = L@da-Rrho
    return da, drho, dict(
        acceleration_block_residual_N=A@da-weight*C@drho+RA,
        density_block_residual_kg_m3=-L@da+drho+Rrho,
        omitted_pressure_residual_force_norm_N=float(np.linalg.norm(weight*C@Rrho)),
        density_response_correction_norm_kg_m3=float(np.linalg.norm(L@da)))


def _rejected_candidate(evaluated):
    """Inspectable failed candidate, never an accepted/restartable state."""
    acc, density, temperature = evaluated['candidate']
    return dict(classification='REJECTED_HHT_NUMERICAL_CANDIDATE', accepted=False,
        restart_allowed=False, physical_candidate_hash=particle_hash(evaluated['pa']),
        acceleration_m_s2=acc.copy(), density_kg_m3=density.copy(), temperature_K=temperature.copy(),
        position_m=positions(evaluated['pa']), velocity_m_s=velocities(evaluated['pa']),
        endpoint_force_N=evaluated['force1'].copy(),
        acceleration_residual_m_s2=evaluated['residual'].copy(),
        equilibrium_residual_N=evaluated['pa'].m[:, None]*evaluated['residual'],
        density_residual_kg_m3=density-evaluated['mapped_density'],
        temperature_residual_K=temperature-evaluated['thermal']['temperature'],
        scaled_merit=evaluated['merit'], residual_metrics=copy.deepcopy(evaluated['metrics']))


def advance(pa, state, p, boundary, dt):
    """Solve one event-free joint candidate and return a detached accepted stage.

All physical DOFs are free. Contact/control discontinuities need caller splits.
At most 24 full residual evaluations, including line-search candidates, occur.
"""
    started = time.perf_counter(); _validate_state(pa, state, p)
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('positive finite HHT stage dt required')
    weight, gamma, beta, direction_limit = _parameters(p)
    accepted_hash = particle_hash(pa); before_state_hash = state['state_hash']
    old = _detach(pa); x0 = positions(old); v0 = velocities(old)
    acc0 = state['acceleration_m_s2']; force0 = state['force_N']
    start = state['time_s']; end = start+dt
    left_boundary = _boundary(old, boundary, start, p)
    if _hash(left_boundary) != _hash(state['boundary_snapshot']):
        raise CoupledStageError('left boundary force/control history changed; explicit event restart required')
    frozen = build_frozen_elastic(old, p['shear'], p['bulk'])
    gradient_map = frozen.gradient_map.reshape(len(old.x), 3, 3, 3*len(old.x))
    divergence_map = np.einsum('niik->nk', gradient_map)
    pressure_map = divergence_map.T*(frozen.volume*p['bulk']/p['rho0'])[None, :]
    bulk_stiffness = divergence_map.T@((frozen.volume*p['bulk'])[:, None]*divergence_map)
    deviatoric_stiffness = frozen.stiffness-bulk_stiffness
    midpoint_acceleration_response = gamma/2+1.5*(beta-gamma/2)
    acceleration_matrix = (np.diag(np.repeat(old.m, 3))+
        weight*midpoint_acceleration_response*dt*dt*deviatoric_stiffness)
    # This is only a numerical preconditioner; actual endpoint forces below
    # are always recomputed from current geometry, density and material state.
    law = FKS(**{key: p[key] for key in FKS.__dataclass_fields__ if key in p})
    evaluations = 0; attempts = []; current = None; last_evaluated = None

    def evaluate(candidate):
        nonlocal evaluations, last_evaluated
        if evaluations >= MAX_EVALUATIONS:
            raise CoupledStageError('HHT joint evaluation budget exceeded', dict(evaluations=evaluations, attempts=attempts))
        evaluations += 1
        acc1, density1, temperature1 = candidate
        acc1 = _array(acc1, v0.shape, 'candidate acceleration')
        density1 = _array(density1, old.rho.shape, 'candidate density')
        temperature1 = _array(temperature1, old.temperature.shape, 'candidate temperature')
        endpoint = _detach(old)
        x1, v1 = hht_path.endpoint(x0, v0, acc0, acc1, dt, weight)
        endpoint.x[:], endpoint.y[:], endpoint.z[:] = x1.T
        endpoint.u[:], endpoint.v[:], endpoint.w[:] = v1.T
        endpoint.rho[:] = density1; endpoint.temperature[:] = temperature1
        _valid_particles(endpoint, p)
        middle = _detach(old)
        xm, vm = hht_path.sample(x0, v0, acc0, acc1, dt, .5, weight)
        middle.x[:], middle.y[:], middle.z[:] = xm.T
        middle.u[:], middle.v[:], middle.w[:] = vm.T
        middle.rho[:] = .5*(old.rho+density1)
        middle.temperature[:] = .5*(old.temperature+temperature1)
        i, j, raw_gradient, gradient, _, correction = corrected_operators(middle)
        stress, ep, plastic_density, return_density, rate = update(old.S.reshape(-1, 3, 3).copy(),
            old.ep.copy(), gradient, middle.temperature, dt, p['shear'], law)
        if np.any(ep != old.ep) or np.any(plastic_density != 0) or np.any(return_density != 0):
            raise CoupledStageError('plastic response outside registered elastic coupled benchmark scope')
        endpoint.S[:] = stress.ravel(); endpoint.ep[:] = ep
        endpoint.strain_rate[:] = rate
        endpoint.D[:] = (.5*(gradient+gradient.swapaxes(1, 2))).ravel()
        divergence = np.trace(gradient, axis1=1, axis2=2)
        numerator = 1-.5*dt*divergence; denominator = 1+.5*dt*divergence
        if np.any(numerator <= 0) or np.any(denominator <= 0):
            raise CoupledStageError('nonpositive implicit continuity density factor')
        mapped_density = old.rho*numerator/denominator
        internal, spatial = _force(endpoint, p)
        right_boundary = _boundary(endpoint, boundary, end, p)
        force1 = internal+right_boundary['total_boundary_force_N']
        contact = _contact_work(left_boundary, right_boundary, old, v1, acc0, acc1, dt, weight, direction_limit)
        thermal = advance_implicit_thermal(old.m, middle.rho, old.temperature, p['conductivity'],
            xm, i, j, raw_gradient, dt, p['cp'], np.zeros(len(old.x)), contact['friction_raw_J'],
            contact['friction_rounding_bound_J'], contact['wall_row_W_K'],
            contact['wall_weighted_temperature_W'], p['plastic_heat_fraction'], p['heat_partition'],
            p['room'], p['melt'])
        mapped_acc = (weight*force1+(1-weight)*force0)/old.m[:, None]
        residual = acc1-mapped_acc
        metrics = dict(velocity_residual_m_s=float(dt*np.max(np.abs(residual))),
            position_residual_m=float(beta*dt*dt*np.max(np.abs(residual))),
            density_relative_residual=float(np.max(np.abs(density1-mapped_density))/p['rho0']),
            temperature_residual_K=float(np.max(np.abs(temperature1-thermal['temperature']))))
        merit = max(metrics['velocity_residual_m_s']/1e-9,
            metrics['position_residual_m']/1e-12, metrics['density_relative_residual']/1e-10,
            metrics['temperature_residual_K']/5e-12)
        last_evaluated = dict(pa=endpoint, candidate=(acc1.copy(), density1.copy(), temperature1.copy()),
            residual=residual, mapped_density=mapped_density, thermal=thermal, force1=force1,
            boundary=right_boundary, contact=contact, spatial=spatial, metrics=metrics, merit=merit,
            middle_position_m=xm, middle_velocity_m_s=vm, gradient_per_s=gradient)
        return last_evaluated

    candidate = (acc0.copy(), old.rho.copy(), old.temperature.copy())
    try:
        current = evaluate(candidate)
        while current['merit'] > 1:
            acc1, rho1, temp1 = current['candidate']
            divergence = np.trace(current['gradient_per_s'], axis1=1, axis2=2)
            density_response = (-(old.rho*dt*dt*midpoint_acceleration_response/
                (1+.5*dt*divergence)**2)[:, None]*divergence_map)
            da, density_correction, block_diagnostic = density_pressure_correction(
                acceleration_matrix, pressure_map, density_response,
                (old.m[:, None]*current['residual']).ravel(), rho1-current['mapped_density'], weight)
            correction = da.reshape(v0.shape)
            improved = None; last_error = None
            for backtrack in range(7):
                fraction = 2.**(-backtrack)
                trial = (acc1+fraction*correction,
                    rho1+fraction*density_correction,
                    temp1+fraction*(current['thermal']['temperature']-temp1))
                try:
                    result = evaluate(trial)
                except (ValueError, CoupledStageError) as error:
                    last_error = error
                    attempts.append(dict(evaluation=evaluations, step_fraction=fraction, rejected=str(error)))
                    if evaluations >= MAX_EVALUATIONS: raise
                    continue
                attempts.append(dict(evaluation=evaluations, step_fraction=fraction,
                    merit=result['merit'], density_pressure_block=copy.deepcopy(block_diagnostic)))
                if result['merit'] <= 1 or result['merit'] < current['merit']*(1-1e-4*fraction):
                    improved = result; break
            if improved is None:
                raise CoupledStageError('bounded HHT line search failed', dict(evaluations=evaluations,
                    last_metrics=current['metrics'], attempts=attempts, last_error=str(last_error)))
            current = improved
        acc1, _, _ = current['candidate']; result_pa = current['pa']
        # Persist the exact candidate that satisfied the coupled residuals.
        result_pa.ax[:], result_pa.ay[:], result_pa.az[:] = acc1.T
        result_pa.au[:], result_pa.av[:], result_pa.aw[:] = acc1.T
        result_pa.arho[:] = (result_pa.rho-old.rho)/dt
        result_pa.ep_rate[:] = (result_pa.ep-old.ep)/dt
        result_pa.plastic_work[:] = old.plastic_work
        result_pa.return_dissipation[:] = old.return_dissipation
        kinetic = hht_path.kinetic_identity(old.m, v0, velocities(result_pa), acc0, acc1,
            force0, current['force1'], dt, weight)
        applied_work = hht_path.linear_work(left_boundary['noncontact_force_N'],
            current['boundary']['noncontact_force_N'], v0, velocities(result_pa), acc0, acc1, dt, weight)
        if particle_hash(pa) != accepted_hash or state['state_hash'] != before_state_hash:
            raise RuntimeError('joint trial mutated accepted state')
        new_state = _state(result_pa, p, end, acc1, current['force1'], current['boundary'],
            parent=before_state_hash, dt=float(dt))
        diagnostic = dict(classification=CLASSIFICATION, status='ACCEPTED_ELASTIC_COUPLED_BENCHMARK_STAGE',
            accepted=True, production_authorized=False, material_energy_certified=False,
            time_before_s=start, time_after_s=end, dt_s=float(dt), evaluations=evaluations,
            elapsed_s=time.perf_counter()-started, residuals=current['metrics'], attempts=attempts,
            contact=current['contact'], thermal=current['thermal'], kinetic=kinetic, spatial=current['spatial'],
            applied_work_J=float(np.sum(applied_work['work_J'])),
            plastic_work_J=0., return_dissipation_J=0.,
            maximum_sampled_coulomb_relative_defect=max((row['maximum_sampled_coulomb_relative_defect']
                for row in current['contact']['contacts']), default=0.),
            normal_temporal_residual_J=current['contact']['normal_potential_change_J']+current['contact']['normal_relative_work_J'],
            numerical_preconditioner='FROZEN_DENSITY_PRESSURE_SCHUR_BLOCK_WITH_DENSE_MIDPOINT_RESPONSE; ACTUAL_FORCE_RESIDUAL_ALWAYS_REEVALUATED',
            geometry_time_layer='TRUE_DENSE_PATH_MIDPOINT_AND_TRUE_NEWMARK_ENDPOINT',
            material_time_layer='UNCHANGED_FKS_ELASTIC_BRANCH_WITH_MIDPOINT_L_T',
            thermal_time_layer='BE_T1_WITH_DENSE_MIDPOINT_RADIAL_PAIRS_AND_RAW_GAUSS_CONTACT_Q',
            algorithmic_terms_are_material_heat=False, hypoelastic_energy_status='DIAGNOSTIC_NOT_CERTIFIED_FREE_ENERGY',
            accepted_parent_hash=before_state_hash, accepted_state_hash=new_state['state_hash'])
        return result_pa, new_state, diagnostic
    except Exception as error:
        if particle_hash(pa) != accepted_hash or state['state_hash'] != before_state_hash:
            raise RuntimeError('failed joint trial mutated accepted state') from error
        if isinstance(error, CoupledStageError):
            error.diagnostics.update(accepted=False, restart_allowed=False, evaluations=evaluations,
                elapsed_s=time.perf_counter()-started, accepted_parent_hash=before_state_hash)
            if current is not None:
                error.diagnostics['last_line_search_base_candidate'] = _rejected_candidate(current)
            if last_evaluated is not None:
                error.diagnostics['last_evaluated_candidate'] = _rejected_candidate(last_evaluated)
            raise
        raise CoupledStageError(str(error), dict(accepted=False, restart_allowed=False,
            evaluations=evaluations, elapsed_s=time.perf_counter()-started,
            accepted_parent_hash=before_state_hash)) from error

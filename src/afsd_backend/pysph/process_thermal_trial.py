"""Named same-history joint trial with conservative BE thermal feedback.

Only the thermal time discretization differs from frozen UpdatedTrial.
Mechanical statements in evaluate retain their original order and arguments;
see PROCESS_GLOBAL_THERMAL_VARIANT.md for the complete difference inventory.
No accepted state, geometry, contact model or material model is changed here.
"""
import numpy as np

from .updated_trial import UpdatedTrial, FIELDS, particle_hash
from .mechanics import positions, velocities, corrected_operators, internal_force, pressure
from .material import update
from .axial_actuator import validate_interval
from .coupled_boundary import evaluate_boundary, normal_snapshot
from .implicit_thermal import advance_implicit_thermal, CLASSIFICATION as COMPONENT_CLASSIFICATION
from .stage_iteration import residual_norms
from .mechanics import volumetric_energy

CLASSIFICATION = 'PROCESS_GLOBAL_EVENT_AVERAGE_BE_THERMAL_APPROXIMATE'
THERMAL_CLASSIFICATION = 'PROCESS_GLOBAL_BE_THERMAL_UNCOMMITTED_MAP'
TEMPERATURE_LAYERS = 'FKS_TMID; PAIR_AND_WALL_CONDUCTION_T1_BE; INCREMENT_FROM_T0'


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _close(a, b, atol, rtol, message):
    a, b = np.asarray(a, float), np.asarray(b, float)
    _require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
             and np.all(np.abs(a-b) <= atol+rtol*np.abs(b)), message)


def certify_boundary_heat(boundary, particle_ids, dt, p):
    """Verify EACH actual accepted contact interval before summing heat.

Independent lab-vector work and frame-relative roundoff reconstruction.
This checks the existing kick's certificate, never enlarges it or clips heat.
It does not replace global momentum/constraint/material acceptance.
"""
    ids = np.asarray(particle_ids); n = len(ids)
    wall = np.zeros(n); bath = np.zeros(n); bound = np.zeros(n); heat = np.zeros(n)
    seen = set(); interval_count = 0; minimum = None
    eps = np.finfo(float).eps
    _require(np.isfinite(dt) and dt > 0, 'positive certificate interval required')
    hface = p['contact_conductance']*p['dx']**2
    _require(np.isfinite(hface) and hface >= 0, 'invalid wall conductance')
    for trace in boundary['traces']:
        idx = trace['particle_index']
        _require(0 <= idx < n and idx not in seen and ids[idx] == trace['particle_id'], 'trace identity mismatch')
        seen.add(idx); faces = trace['contact_rows']; result = trace['result']
        normals = np.asarray([f['normal'] for f in faces], float).reshape(-1, 3)
        walls = np.asarray([f['wall_velocity'] for f in faces], float).reshape(-1, 3)
        _require(len({f['contact_id'] for f in faces}) == len(faces), 'duplicate contact face')
        for face in faces:
            expected = p['base_temperature'] if int(face['feature_id']) == 0 else p['tool_temperature']
            _require(face['normal_force_N'] > 0 and face['wall_temperature_K'] == expected,
                     'active face/reservoir mismatch')
            wall[idx] += hface; bath[idx] += hface*expected
        nc = len(faces); covered = 0.; previous = np.asarray(trace['material_velocity_before_m_s'])
        expected_frame = walls[0] if nc else previous
        normal_force = np.asarray([f['normal_force_N'] for f in faces])[:, None]*normals
        free = np.asarray(trace['free_axes'])
        _require(np.all((free == 0) | (free == 1)) and trace['mass_kg'] > 0, 'invalid contact mobility')
        other = np.asarray(trace['internal_force_N'])+normal_force.sum(axis=0)
        expected_speed = max(np.linalg.norm(previous-expected_frame),
            float(np.linalg.norm(walls-expected_frame)) if nc else 0.,
            dt*np.linalg.norm(free*other/trace['mass_kg']), np.finfo(float).tiny)
        for interval in result['intervals']:
            interval_count += 1
            h = float(interval['dt_s']); start = float(interval['start_s'])
            _require(np.isfinite(h) and h > 0 and abs(start-covered) <= 128*eps*dt, 'contact interval coverage')
            before = np.asarray(interval['v_before']); after = np.asarray(interval['v_after'])
            _close(before, previous, 1e-15, 1e-13, 'contact velocity history mismatch')
            jt = np.asarray(interval['tangential_impulse']).reshape(nc, 3)
            raw = np.asarray(interval['friction_heat_raw_J']); qround = np.asarray(interval['friction_rounding_bound_J'])
            _require(raw.shape == (nc,) and qround.shape == (nc,) and np.isfinite(raw).all()
                     and np.isfinite(qround).all() and np.all(qround >= 0), 'invalid raw friction certificate')
            frame = np.asarray(interval['reference_frame_velocity'])
            _close(frame, expected_frame, 0., 0., 'contact numerical frame mismatch')
            rel0 = before-frame; rel1 = after-frame; wall_relative = walls-frame
            speed = float(interval['relative_speed_scale_m_s'])
            _require(np.isfinite(speed) and speed > 0, 'invalid slip speed scale')
            _close(np.asarray(speed), np.asarray(expected_speed), 0., 1e-13, 'contact slip scale inflated or changed')
            # Relative frame before averaging avoids cancellation under a lab boost.
            relative = .5*(rel0+rel1)-wall_relative
            tangent = relative-normals*np.sum(normals*relative, axis=1)[:, None]
            supplied_tangent = np.asarray(interval['tangential_relative_midpoint']).reshape(nc, 3)
            velocity_budget = 128*eps*(speed+np.linalg.norm(rel0)+np.linalg.norm(rel1))
            _close(supplied_tangent, tangent, velocity_budget, 0., 'relative midpoint slip mismatch')
            reconstructed = -np.sum(jt*supplied_tangent, axis=1)
            _close(raw, reconstructed, 1e-25, 1e-12, 'friction impulse-work identity')
            expected_bound = 256*eps*np.linalg.norm(jt, axis=1)*(speed+np.linalg.norm(rel0)+
                np.linalg.norm(rel1)+2*np.linalg.norm(wall_relative, axis=1))
            _close(qround, expected_bound, 1e-30, 1e-12, 'friction rounding certificate inflated or changed')
            _require(np.all(raw >= -qround), 'negative individual friction event outside roundoff')
            heat[idx] += float(raw.sum()); bound[idx] += float(qround.sum())
            if nc:
                local = float(raw.min()); minimum = local if minimum is None else min(minimum, local)
            covered += h; previous = after
        _require(abs(covered-dt) <= 128*eps*dt, 'incomplete contact interval')
        _close(previous, result['velocity'], 1e-15, 1e-13, 'contact endpoint history mismatch')
        _close(np.asarray([heat[idx]]), np.asarray([np.sum(result['friction_heat_J'])]), 1e-25, 1e-12,
               'contact accumulated friction mismatch')
    _require(len(seen) == n, 'missing particle contact trace')
    _close(heat, boundary['friction_raw_J'], 1e-25, 1e-12, 'friction source trace mismatch')
    _close(wall, boundary['wall_row_W_K'], 1e-15, 1e-13, 'wall row/face mismatch')
    return dict(status='EACH_CONTACT_SUBINTERVAL_CERTIFIED_BEFORE_SUM',
        wall_row_W_K=wall, wall_weighted_temperature_W=bath,
        friction_rounding_bound_J=bound, friction_raw_J=heat,
        interval_count=interval_count, minimum_raw_friction_J=minimum)


def _end_temperature_boundary_ledger(boundary, certificate, thermal, dt):
    """Thermal journal only; no contact force, impulse or motion is changed."""
    flux = certificate['wall_weighted_temperature_W']-certificate['wall_row_W_K']*thermal['temperature']
    boundary['midpoint_temperature_conduction_diagnostic_W'] = boundary['contact_power_W'].copy()
    boundary['contact_power_W'] = flux.copy(); boundary['contact_conduction_W'] = flux.copy()
    boundary['ledger']['contact_conduction_into_material_J'] = float(dt*flux.sum())
    boundary['ledger']['wall_conduction_out_J'] = -float(dt*flux.sum())
    boundary['contact_power_definition'] = 'BE_T1_SIGNED_CONTACT_CONDUCTION_ONLY; RAW_FRICTION_SEPARATE'
    for trace in boundary['traces']:
        trace['midpoint_temperature_conduction_diagnostic_W'] = trace['contact_conduction_W']
        trace['contact_conduction_W'] = float(flux[trace['particle_index']])
        trace['thermal_flux_time_layer'] = 'BE_T1'


class ProcessThermalTrial(UpdatedTrial):
    """Same accepted history / five-field mechanical map, named BE heat variant."""

    def evaluate(self, candidate, dt):
        if self.accepted_hash() != self._accepted_hash:
            raise RuntimeError('accepted private snapshot changed')
        dt = float(dt)
        validate_interval(self.time, dt, self.p['release_time'])
        c = self._validate_arrays(candidate)
        p, old = self.p, self._accepted
        end_time, mid_time = self.time+dt, self.time+.5*dt
        normal_after = self._check_state(c, end_time)
        normal_before = normal_snapshot(positions(old), old.particle_id, self.time, p)
        pa = old.extract_particles(self._order)
        x0, v0 = positions(old), velocities(old)
        xmid = .5*(x0+c['position_m'])
        rhomid = .5*(old.rho+c['density_kg_m3'])
        tmid = .5*(old.temperature+c['temperature_K'])
        pa.x[:], pa.y[:], pa.z[:] = xmid.T
        pa.u[:], pa.v[:], pa.w[:] = c['mean_velocity_m_s'].T
        pa.rho[:] = rhomid
        pa.temperature[:] = tmid
        i, j, g, gradient, drho, correction = corrected_operators(pa)
        stress0 = old.S.reshape(-1, 3, 3).copy()
        stress1, ep1, plastic_density, return_density, rate = update(
            stress0, old.ep.copy(), gradient, tmid, dt, p['shear'], self.law)
        sigma = .5*(stress0+stress1)-pressure(rhomid, p['rho0'], self.bulk)[:, None, None]*np.eye(3)
        force = internal_force(pa, sigma, i, j, g, correction=correction)
        boundary = evaluate_boundary(xmid, v0, old.m, tmid, mid_time, dt,
                                     p, force, old.z0, old.particle_id)
        normal_mid = boundary['normal_snapshot']
        guard = self._feature._feature_transition_guard
        guards = dict(before_to_mid=guard(normal_before, normal_mid, x0, xmid, self.time, mid_time),
                      mid_to_after=guard(normal_mid, normal_after, xmid, c['position_m'], mid_time, end_time))
        volume = old.m/rhomid
        plastic_j, return_j = volume*plastic_density, volume*return_density
        certificate = certify_boundary_heat(boundary, old.particle_id, dt, p)
        thermal = advance_implicit_thermal(old.m, rhomid, old.temperature, p['conductivity'],
            xmid, i, j, g, dt, p['cp'], plastic_j, boundary['friction_raw_J'],
            certificate['friction_rounding_bound_J'], certificate['wall_row_W_K'],
            certificate['wall_weighted_temperature_W'], p['plastic_heat_fraction'], p['heat_partition'],
            p['room'], p['melt'])
        thermal.update(component_classification=COMPONENT_CLASSIFICATION,
            classification=THERMAL_CLASSIFICATION, status='UNCOMMITTED_BE_COUPLED_MAP')
        _end_temperature_boundary_ledger(boundary, certificate, thermal, dt)
        divergence = np.trace(gradient, axis1=1, axis2=2)
        denominator = 1+.5*dt*divergence
        numerator = 1-.5*dt*divergence
        if np.any(denominator <= 0) or np.any(numerator <= 0):
            raise ValueError('nonpositive midpoint continuity density factor')
        mapped = dict(position_m=x0+dt*boundary['work_velocity'],
            velocity_m_s=boundary['velocity'].copy(), mean_velocity_m_s=boundary['work_velocity'].copy(),
            density_kg_m3=old.rho*numerator/denominator, temperature_K=thermal['temperature'].copy())
        self._check_state(mapped, end_time)
        internal_power = float(np.sum(force*c['mean_velocity_m_s']))
        stress_power = float(np.sum(volume[:, None, None]*sigma*gradient))
        condition = np.linalg.cond(np.linalg.inv(correction))
        spatial = dict(pair_count=len(i), condition_max=float(condition.max()),
            internal_power_W=internal_power, midpoint_stress_power_W=stress_power,
            power_conjugacy_residual_W=internal_power+stress_power,
            internal_net_force_N=force.sum(axis=0), internal_net_torque_Nm=np.cross(xmid, force).sum(axis=0),
            force_absolute_scale_N=float(np.linalg.norm(force, axis=1).sum()),
            torque_absolute_scale_Nm=float(np.linalg.norm(np.cross(xmid, force), axis=1).sum()),
            power_absolute_scale_W=float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(c['mean_velocity_m_s'], axis=1))))
        if self.accepted_hash() != self._accepted_hash:
            raise RuntimeError('trial mutated accepted history')
        return dict(classification=CLASSIFICATION, accepted=False, restart_allowed=False,
            material_scope=self.material_scope, accepted_state_hash=self._accepted_hash,
            time_before_s=self.time, time_after_s=end_time, time_mid_s=mid_time, dt_s=dt,
            particle_id=old.particle_id.copy(), map=mapped,
            residual={key: c[key]-mapped[key] for key in FIELDS},
            continuity_residual_kg_m3_s=(c['density_kg_m3']-old.rho)/dt+rhomid*divergence,
            midpoint=dict(position_m=xmid, density_kg_m3=rhomid, temperature_K=tmid,
                mean_velocity_m_s=c['mean_velocity_m_s'], velocity_gradient_per_s=gradient,
                stress_Pa=sigma, i=i, j=j, g=g, B=correction, mass_kg=old.m.copy()),
            material=dict(deviatoric_stress_Pa=stress1, ep=ep1, plastic_work_J=plastic_j,
                return_dissipation_J=return_j, strain_rate_per_s=rate),
            internal_force_N=force, boundary=boundary, thermal=thermal, spatial=spatial,
            thermal_boundary_certificate=certificate,
            normal_before=normal_before, normal_after=normal_after, feature_guards=guards,
            temperature_layers=TEMPERATURE_LAYERS,
            material_energy_certified=False, production_authorized=False)


def certify_process_thermal(trial, candidate, evaluation):
    """Candidate T1 BE balances, NOT a replacement for full-stage acceptance."""
    e = evaluation; p = trial.p; old = trial._accepted; c = trial._validate_arrays(candidate)
    _require(e['classification'] == CLASSIFICATION and e['accepted'] is False and
             e['restart_allowed'] is False and e['temperature_layers'] == TEMPERATURE_LAYERS,
             'named uncommitted BE trial metadata required')
    _require(e['thermal']['classification'] == THERMAL_CLASSIFICATION, 'BE thermal metadata required')
    _require(trial.accepted_hash() == e['accepted_state_hash'] == particle_hash(old), 'accepted history changed')
    _require(np.array_equal(e['particle_id'], old.particle_id), 'BE particle identity mismatch')
    trial._check_state(c, e['time_after_s'])
    trial._check_state(e['map'], e['time_after_s'])
    cert = certify_boundary_heat(e['boundary'], old.particle_id, e['dt_s'], p)
    th = e['thermal']; mid = e['midpoint']; first = mid['i']; second = mid['j']; dt = e['dt_s']
    t = c['temperature_K']; qpair = dt*th['pair_conductance_W_K']*(t[second]-t[first])
    node = np.zeros(len(t)); np.add.at(node, first, qpair); np.add.at(node, second, -qpair)
    qwall = dt*(cert['wall_weighted_temperature_W']-cert['wall_row_W_K']*t)
    total = p['plastic_heat_fraction']*e['material']['plastic_work_J']+p['heat_partition']*cert['friction_raw_J']+node+qwall
    measured = old.m*p['cp']*(t-old.temperature); residual = measured-total
    scale = float(np.abs(total).sum())
    _require(float(np.max(np.abs(residual))) <= 1e-12+1e-10*scale, 'candidate BE per-particle thermal balance')
    _require(abs(float(residual.sum())) <= 1e-12+1e-9*scale, 'candidate BE total thermal balance')
    _require(abs(float(node.sum())) <= 1e-15+1e-10*float(np.abs(qpair).sum()), 'candidate BE pair balance')
    # Algebraic map journals must independently agree with the mapped end T.
    mapped = e['map']['temperature_K']; pair_map = dt*th['pair_conductance_W_K']*(mapped[second]-mapped[first])
    wall_map = dt*(cert['wall_weighted_temperature_W']-cert['wall_row_W_K']*mapped)
    _close(th['temperature'], mapped, 0., 0., 'map/thermal temperature mismatch')
    _close(th['pair_heat_into_i_J'], pair_map, 1e-25, 1e-12, 'BE mapped pair heat mismatch')
    _close(th['particle_energy_J']['contact_in_J'], wall_map, 1e-25, 1e-12, 'BE mapped wall heat mismatch')
    _close(e['boundary']['contact_power_W']*dt, wall_map, 1e-25, 1e-12, 'BE boundary wall ledger time-layer mismatch')
    return dict(status='CANDIDATE_BE_HEAT_CHECK_ONLY_NOT_GLOBAL_ACCEPTANCE',
        thermal_residual_J=float(residual.sum()), maximum_particle_thermal_residual_J=float(np.max(np.abs(residual))),
        pair_net_J=float(node.sum()), candidate_wall_heat_into_material_J=float(qwall.sum()),
        contact_subintervals=cert['interval_count'], former_explicit_dt_bound_s=th['former_explicit_dt_bound_s'],
        former_explicit_bound_is_not_a_MECHANICAL_macrostep_certificate=True,
        accepted=False, production_authorized=False)


def acceptance(pa, p, trial, c, e):
    """Complete old mechanical acceptance with explicitly named BE heat checks.

The old acceptance routine is not called with forged classification/fields.
Mechanical thresholds and energy *diagnostic* definitions are retained.
Only explicit thermal-bound/pair-power checks are replaced by actual BE checks.
"""
    require = _require
    require(e['classification'] == CLASSIFICATION and e['accepted'] is False and e['restart_allowed'] is False,
            'uncommitted registered BE trial metadata required')
    norms, scaled = residual_norms(c, e['map'], pa.rho)
    require(scaled <= 1, 'global residual not converged')
    require(trial.accepted_hash() == particle_hash(pa) == e['accepted_state_hash'], 'accepted history changed')
    require(np.array_equal(e['particle_id'], pa.particle_id), 'ID correspondence lost')
    require(e['material_scope'] == 'UNCHANGED_FKS_APPROXIMATE_PARAMETERS', 'manufactured material cannot be accepted')
    s = e['spatial']; b = e['boundary']; ledger = b['ledger']; th = e['thermal']
    require(abs(s['power_conjugacy_residual_W']) <= 1e-12+1e-10*s['power_absolute_scale_W'], 'spatial power')
    require(np.linalg.norm(s['internal_net_force_N']) <= 1e-12+1e-10*s['force_absolute_scale_N'], 'internal net force')
    require(np.linalg.norm(s['internal_net_torque_Nm']) <= 1e-14+1e-10*s['torque_absolute_scale_Nm'], 'internal net torque')
    heat_checks = certify_process_thermal(trial, c, e)
    thermal_residual = heat_checks['thermal_residual_J']
    work_scale = sum(abs(ledger[k]) for k in ('delta_kinetic_J','wall_work_J','actuator_work_J','internal_work_J','normal_relative_work_J','friction_raw_J'))
    for key in ('kinetic_work_residual_J', 'contact_kinetic_residual_J'):
        require(abs(ledger[key]) <= 1e-15+1e-9*work_scale, key)
    seen = set(); max_momentum = 0.; intervals = 0
    for trace in b['traces']:
        idx = trace['particle_index']; require(idx not in seen and trace['particle_id'] == pa.particle_id[idx], 'trace IDs')
        seen.add(idx); result = trace['result']; free = np.asarray(trace['free_axes'])
        reaction = np.asarray(result['constraint_impulse'])
        require(np.max(np.abs(free*reaction)) <= 1e-17, 'reaction on free axis')
        require(np.max(np.abs((1-free)*(c['velocity_m_s'][idx]-trace['prescribed_velocity']))) <= 1e-12, 'axial endpoint constraint')
        impulse = np.asarray(result['normal_impulse']).reshape(-1,3).sum(axis=0)+np.asarray(result['tangential_impulse']).reshape(-1,3).sum(axis=0)+reaction
        defect = pa.m[idx]*(c['velocity_m_s'][idx]-velocities(pa)[idx])-e['dt_s']*e['internal_force_N'][idx]-impulse
        max_momentum = max(max_momentum, float(np.max(np.abs(defect))))
        require(np.max(np.abs(defect)) <= pa.m[idx]*1e-9+1e-17, 'candidate momentum identity')
        mean = np.zeros(3); covered = 0.
        for interval in result['intervals']:
            intervals += 1; h = interval['dt_s']; covered += h
            require(h > 0, 'nonpositive contact subinterval')
            mean += h*.5*(np.asarray(interval['v_before'])+interval['v_after'])/e['dt_s']
            require(np.all(np.asarray(interval['friction_heat_raw_J']) >= -np.asarray(interval['friction_rounding_bound_J'])), 'negative raw friction beyond roundoff')
        require(abs(covered-e['dt_s']) <= 128*np.finfo(float).eps*e['dt_s'], 'incomplete contact interval')
        require(np.max(np.abs(mean-c['mean_velocity_m_s'][idx])) <= 1e-9, 'event-average velocity mismatch')
    require(len(seen) == len(pa.x), 'missing particle traces')
    require(np.all(e['material']['ep'] >= pa.ep), 'plastic history decreased')
    v0 = velocities(pa); v1 = c['velocity_m_s']; rho1 = c['density_kg_m3']; S0 = pa.S.reshape(-1,3,3); S1 = e['material']['deviatoric_stress_Pa']
    dK = float(.5*np.sum(pa.m[:,None]*(v1*v1-v0*v0)))
    dEs = float(np.sum(pa.m/rho1*np.sum(S1*S1,axis=(1,2))/(4*p['shear'])-pa.m/pa.rho*np.sum(S0*S0,axis=(1,2))/(4*p['shear'])))
    dEv = float(np.sum(pa.m*(volumetric_energy(rho1,p['rho0'],trial.bulk)-volumetric_energy(pa.rho,p['rho0'],trial.bulk))))
    dU = e['normal_after']['normal_potential_J']-e['normal_before']['normal_potential_J']
    dEth = float(np.sum(pa.m*p['cp']*(c['temperature_K']-pa.temperature)))
    diagnostic = dK+dEs+dEv+dU+dEth+th['ledger']['plastic_unthermalized_J']-ledger['wall_work_J']-ledger['actuator_work_J']+th['ledger']['wall_heat_out_J']
    x0=positions(pa); x1=c['position_m']; xmid=e['midpoint']['position_m']
    angular_delta=np.sum(pa.m[:,None]*(np.cross(x1,v1)-np.cross(x0,v0)),axis=0)
    external_angular=sum((np.asarray(ledger[k+'_angular_impulse_Nms']) for k in ('base','tool','actuator')),np.zeros(3))
    extra_kinematic=np.sum(pa.m[:,None]*np.cross(x1-x0,.5*(v0+v1)),axis=0)
    midpoint_internal_torque=e['dt_s']*np.cross(xmid,e['internal_force_N']).sum(axis=0)
    nonlinear_work=e['dt_s']*float(np.sum(e['internal_force_N']*(c['mean_velocity_m_s']-b['work_velocity'])))
    return dict(status='ACCEPTED_NUMERICAL_STAGE_APPROXIMATE_NOT_ENERGY_CERTIFIED', variant=CLASSIFICATION,
        residual_norms=norms, maximum_scaled_residual=scaled, momentum_residual_max_Ns=max_momentum,
        thermal_residual_J=thermal_residual, process_thermal_checks=heat_checks,
        contact_subintervals=intervals, spatial=s, contact_ledger=ledger, thermal_ledger=th['ledger'],
        delta_K_J=dK, delta_Es_diagnostic_J=dEs, delta_volumetric_diagnostic_J=dEv,
        candidate_minus_contact_map_delta_K_J=dK-ledger['delta_kinetic_J'],
        delta_normal_potential_J=dU, delta_thermal_J=dEth, total_energy_diagnostic_residual_J=diagnostic,
        normal_temporal_residual_J=dU+ledger['normal_relative_work_J'],
        nonlinear_internal_work_defect_J=nonlinear_work,
        angular_momentum_change_kg_m2_s=angular_delta,external_midpoint_angular_impulse_Nms=external_angular,
        angular_residual_Nms=angular_delta-external_angular,
        event_average_kinematic_angular_term_Nms=extra_kinematic,
        angular_identity_defect_Nms=angular_delta-external_angular-midpoint_internal_torque-extra_kinematic,
        angular_status='FROZEN_MIDPOINT_MOMENT_APPROXIMATE_EXTRA_KINEMATIC_TERM_REPORTED_NOT_CORRECTED',
        total_energy_status='DIAGNOSTIC_HYPOELASTIC_NOT_FREE_ENERGY_CERTIFIED_NO_RESIDUAL_HEAT',
        plastic_particle_ids=pa.particle_id[e['material']['ep']>pa.ep], production_authorized=False)

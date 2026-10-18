"""Private fixed-quadrature response experiment; old physics/root kept verbatim.

Only the numerical search matrix differs from hht_stop_stage.advance_private.
This function is NEVER a publication or physical event implementation.
"""
import copy
import time
import numpy as np
from . import hht_path
from .hht_stop_stage import (
    CoupledStageError,_validate_state,_parameters,_detach,_hash,_array,
    _valid_particles,_boundary,_force,_rejected_candidate,_state,
    density_pressure_correction,MAX_EVALUATIONS,FKS,update,positions,velocities,
    corrected_operators,advance_implicit_thermal,build_frozen_elastic,
    particle_hash,integrate_contact,curved_kinetic_identity,CLASSIFICATION,INTEGRATION)
from .contact_response_preconditioner import assemble,STRATEGY

def advance_private_response(pa, state, p, boundary, dt, gap_certificate, event_context=None):
    """Solve a joint candidate; NEVER accepted or restartable before certification.

All physical DOFs are free. Contact/control discontinuities need caller splits.
At most 24 full residual evaluations, including line-search candidates, occur.
"""
    started = time.perf_counter(); _validate_state(pa, state, p)
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('positive finite HHT stage dt required')
    if np.any(pa.actuator != 0):
        raise ValueError('actual actuator constraints outside registered numerical experiment')
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
    evaluations = 0; attempts = []; current = None; last_evaluated = None; pending_candidate = None
    response_history = []

    def evaluate(candidate):
        nonlocal evaluations, last_evaluated, pending_candidate
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
        pending_candidate = dict(classification='REJECTED_UNEQUILIBRATED_CURVED_CANDIDATE',
            accepted=False,restart_allowed=False,time_s=end,evaluation=evaluations,
            position_m=x1.copy(),velocity_m_s=v1.copy(),acceleration_m_s2=acc1.copy(),
            density_kg_m3=density1.copy(),temperature_K=temperature1.copy(),
            deviatoric_stress_Pa=stress.copy(),particle_id=endpoint.particle_id.copy(),
            endpoint_physical_force_N=force1.copy(),particle_hash=particle_hash(endpoint))
        contact = integrate_contact(old, p, boundary, start, acc0, acc1, dt, weight,
            old_boundary=left_boundary, new_boundary=right_boundary, gap_certificate=gap_certificate,
            event_context=event_context)
        thermal = advance_implicit_thermal(old.m, middle.rho, old.temperature, p['conductivity'],
            xm, i, j, raw_gradient, dt, p['cp'], np.zeros(len(old.x)), contact['friction_raw_J'],
            contact['friction_rounding_bound_J'], contact['wall_row_W_K'],
            contact['wall_weighted_temperature_W'], p['plastic_heat_fraction'], p['heat_partition'],
            p['room'], p['melt'])
        mapped_acc = (weight*force1+(1-weight)*force0+
            contact['delta_impulse_Ns']/(gamma*dt))/old.m[:, None]
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
            contact_matrix, response_record = assemble(current['contact'],current['boundary'],
                velocities(current['pa']),dt,p)
            response_matrix=acceleration_matrix+contact_matrix
            response_record.update(at_evaluation=evaluations,raw_merit=current['merit'],
                numerical_matrix_condition=float(np.linalg.cond(response_matrix)))
            response_history.append(response_record)
            da, density_correction, block_diagnostic = density_pressure_correction(
                response_matrix, pressure_map, density_response,
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
                    attempts.append(dict(evaluation=evaluations, step_fraction=fraction,rejected=str(error),
                        rejection_diagnostic=copy.deepcopy(getattr(error,'diagnostics',{}))))
                    if evaluations >= MAX_EVALUATIONS: raise
                    continue
                attempts.append(dict(evaluation=evaluations, step_fraction=fraction,
                    merit=result['merit'], density_pressure_block=copy.deepcopy(block_diagnostic)))
                if result['merit'] <= 1 or result['merit'] < current['merit']*(1-1e-4*fraction):
                    improved = result; break
            if improved is None:
                raise CoupledStageError('bounded HHT line search failed', dict(evaluations=evaluations,
                    last_metrics=current['metrics'], attempts=attempts, last_error=str(last_error),
                    last_error_diagnostic=copy.deepcopy(getattr(last_error,'diagnostics',{}))))
            current = improved
        acc1, _, _ = current['candidate']; result_pa = current['pa']
        # Persist the exact candidate that satisfied the coupled residuals.
        result_pa.ax[:], result_pa.ay[:], result_pa.az[:] = acc1.T
        result_pa.au[:], result_pa.av[:], result_pa.aw[:] = acc1.T
        result_pa.arho[:] = (result_pa.rho-old.rho)/dt
        result_pa.ep_rate[:] = (result_pa.ep-old.ep)/dt
        result_pa.plastic_work[:] = old.plastic_work
        result_pa.return_dissipation[:] = old.return_dissipation
        kinetic = curved_kinetic_identity(old.m, v0, velocities(result_pa), acc0, acc1,
            force0, current['force1'], dt, current['contact']['delta_impulse_Ns'],
            current['contact']['delta_material_work_J'], weight)
        applied_work = hht_path.linear_work(left_boundary['noncontact_force_N'],
            current['boundary']['noncontact_force_N'], v0, velocities(result_pa), acc0, acc1, dt, weight)
        if particle_hash(pa) != accepted_hash or state['state_hash'] != before_state_hash:
            raise RuntimeError('joint trial mutated accepted state')
        new_state = _state(result_pa, p, end, acc1, current['force1'], current['boundary'],
            parent=before_state_hash, dt=float(dt),interval=dict(
                integration_definition=INTEGRATION,
                delta_impulse_Ns=current['contact']['delta_impulse_Ns'],
                delta_material_work_J=current['contact']['delta_material_work_J'],
                actual_contact_impulse_Ns=current['contact']['actual_impulse_Ns'],
                actual_contact_material_work_J=current['contact']['actual_material_work_J'],
                raw_friction_J=current['contact']['friction_raw_J'],
                quadrature=current['contact']['quadrature'],contacts=current['contact']['contacts'],
                kinetic=kinetic,algorithmic_work_is_not_heat=True))
        diagnostic = dict(classification=CLASSIFICATION, status='PRIVATE_EQUILIBRATED_CANDIDATE_REQUIRES_PATH_CERTIFICATE',
            accepted=False,restart_allowed=False, production_authorized=False, material_energy_certified=False,
            time_before_s=start, time_after_s=end, dt_s=float(dt), evaluations=evaluations,
            elapsed_s=time.perf_counter()-started, residuals=current['metrics'], attempts=attempts,
            contact=current['contact'], thermal=current['thermal'], kinetic=kinetic, spatial=current['spatial'],
            applied_work_J=float(np.sum(applied_work['work_J'])),
            plastic_work_J=0., return_dissipation_J=0.,
            maximum_sampled_coulomb_relative_defect=max((row['maximum_sampled_coulomb_relative_defect']
                for row in current['contact']['contacts']), default=0.),
            normal_temporal_residual_J=current['contact']['normal_potential_change_J']+current['contact']['normal_relative_work_J'],
            numerical_preconditioner=STRATEGY,contact_response_history=response_history,
            continuous_cross_zero_jacobian_certified=False,can_publish=False,
            geometry_time_layer='TRUE_DENSE_PATH_MIDPOINT_AND_TRUE_NEWMARK_ENDPOINT',
            material_time_layer='UNCHANGED_FKS_ELASTIC_BRANCH_WITH_MIDPOINT_L_T',
            thermal_time_layer='BE_T1_WITH_DENSE_MIDPOINT_RADIAL_PAIRS_AND_RAW_CURVED_CONTACT_Q',
            integration_definition=INTEGRATION,
            classical_hht_stability_certified=False,
            algorithmic_terms_are_material_heat=False, hypoelastic_energy_status='DIAGNOSTIC_NOT_CERTIFIED_FREE_ENERGY',
            accepted_parent_hash=before_state_hash, accepted_state_hash=new_state['state_hash'])
        new_state.update(accepted=False,restart_allowed=False)
        new_state.pop('state_hash');new_state['state_hash']=_hash(new_state)
        diagnostic['candidate_state_hash']=new_state['state_hash']
        diagnostic.pop('accepted_state_hash')
        return result_pa, new_state, diagnostic
    except Exception as error:
        if particle_hash(pa) != accepted_hash or state['state_hash'] != before_state_hash:
            raise RuntimeError('failed joint trial mutated accepted state') from error
        if isinstance(error, CoupledStageError):
            error.diagnostics.update(accepted=False, restart_allowed=False, evaluations=evaluations,
                elapsed_s=time.perf_counter()-started, accepted_parent_hash=before_state_hash,
                numerical_preconditioner=STRATEGY,contact_response_history=response_history)
            error.diagnostics['pending_unequilibrated_candidate'] = pending_candidate
            if current is not None:
                error.diagnostics['last_line_search_base_candidate'] = _rejected_candidate(current)
            if last_evaluated is not None:
                error.diagnostics['last_evaluated_candidate'] = _rejected_candidate(last_evaluated)
            raise
        detail=copy.deepcopy(getattr(error,'diagnostics',{}))
        detail.update(accepted=False,restart_allowed=False,evaluations=evaluations,
            elapsed_s=time.perf_counter()-started,accepted_parent_hash=before_state_hash,
            pending_unequilibrated_candidate=pending_candidate,
            numerical_preconditioner=STRATEGY,contact_response_history=response_history)
        raise CoupledStageError(str(error),detail) from error

"""Located-event transactions; no velocity edits and no hidden heat correction."""
import copy
import numpy as np

from . import hht_path, hht_stop_stage as stage
from .hht_stop_contact import integrate
from .hht_zero_event import zero_event_limits, acceleration_storage_jump
from .hht_coupled import _hash, _force, _boundary, _detach
from .mechanics import velocities
from .updated_trial import particle_hash


def certify_interval(old, parent, candidate, trial_state, diagnostic, p,
                     boundary, gap_certificate, event_context=None):
    """Re-audit a converged private path under strict physical endpoint rules.

    A PRIVATE continuation context cannot certify anything. Exact integral
    equality prevents promotion after changing force/heat rules. Publication
    still additionally requires driver geometry/reference/receipt checks.
    """
    stage._validate_state(old, parent, p)
    if trial_state.get('accepted') is not False or trial_state.get('restart_allowed') is not False:
        raise ValueError('only a private equilibrated candidate can be certified')
    payload={k:v for k,v in trial_state.items() if k!='state_hash'}
    if trial_state.get('state_hash')!=_hash(payload) or trial_state['particle_hash']!=particle_hash(candidate):
        raise ValueError('private candidate semantic/particle hash changed')
    if trial_state['parent_state_hash']!=parent['state_hash']:
        raise ValueError('candidate belongs to another accepted parent')
    if (diagnostic.get('candidate_state_hash')!=trial_state['state_hash'] or
            diagnostic.get('accepted_parent_hash')!=parent['state_hash'] or
            diagnostic.get('time_before_s')!=parent['time_s'] or
            diagnostic.get('time_after_s')!=trial_state['time_s'] or
            diagnostic.get('dt_s')!=trial_state['last_dt_s'] or
            trial_state['time_s']!=parent['time_s']+trial_state['last_dt_s']):
        raise ValueError('coupled residual receipt does not bind this candidate')
    context=copy.deepcopy(event_context or {})
    if any(row['endpoint']=='private_continuation' for row in context.values()):
        raise ValueError('private branch continuation is never a physical certificate')
    for key,limit in [('velocity_residual_m_s',1e-9),('position_residual_m',1e-12),
                      ('density_relative_residual',1e-10),('temperature_residual_K',5e-12)]:
        value=diagnostic['residuals'][key]
        if not np.isfinite(value) or value>limit:raise ValueError('raw coupled residual failed '+key)
    contact=integrate(old,p,boundary,parent['time_s'],parent['acceleration_m_s2'],
        trial_state['acceleration_m_s2'],trial_state['last_dt_s'],p.get('hht_weight',.9),
        old_boundary=parent['boundary_snapshot'],new_boundary=trial_state['boundary_snapshot'],
        gap_certificate=gap_certificate,event_context=context)
    if not contact['physical_contact_path_certified']:
        raise ValueError('whole physical interval not certified')
    for key in ('actual_impulse_Ns','actual_material_work_J','friction_raw_J',
                'delta_impulse_Ns','delta_material_work_J'):
        if not np.array_equal(np.asarray(contact[key]),np.asarray(diagnostic['contact'][key])):
            raise ValueError('strict contact audit changed solved physical quantity '+key)
    for saved,key in [('actual_contact_impulse_Ns','actual_impulse_Ns'),
                      ('actual_contact_material_work_J','actual_material_work_J'),
                      ('raw_friction_J','friction_raw_J'),('delta_impulse_Ns','delta_impulse_Ns'),
                      ('delta_material_work_J','delta_material_work_J')]:
        if not np.array_equal(np.asarray(trial_state['last_interval'][saved]),np.asarray(contact[key])):
            raise ValueError('saved candidate interval differs from strict physical integral '+saved)
    result=copy.deepcopy(trial_state);result.update(accepted=True,restart_allowed=True)
    result['last_interval']['contacts']=copy.deepcopy(contact['contacts'])
    result['last_interval']['quadrature']=copy.deepcopy(contact['quadrature'])
    result['last_interval']['event_context']=context
    result.pop('state_hash');result['state_hash']=_hash(result)
    stage._validate_state(candidate,result,p)
    detail=copy.deepcopy(diagnostic);detail.update(accepted=True,restart_allowed=True,
        status='CERTIFIED_LOCAL_EVENT_CAPABLE_INTERVAL_NOT_PROCESS',contact=contact,
        accepted_state_hash=result['state_hash'])
    return result,detail


def discover_events(old, state, p, boundary, dt, acceleration1):
    """Discover candidate full-vector roots from an actual failed dense path.

    This is a predictor ONLY: every trial duration later re-solves the complete
    coupled root. No tracer/AS/RS field or expected reference time is consulted.
    """
    rows0=state['boundary_snapshot']['contacts'];temporary=_detach(old)
    x1,v1=hht_path.endpoint(np.column_stack((old.x,old.y,old.z)),velocities(old),
        state['acceleration_m_s2'],acceleration1,dt,p['hht_weight'])
    temporary.x[:],temporary.y[:],temporary.z[:]=x1.T
    temporary.u[:],temporary.v[:],temporary.w[:]=v1.T
    rows1={row['contact_id']:row for row in _boundary(temporary,boundary,state['time_s']+dt,p)['contacts']}
    found=[]
    for row in rows0:
        idx=row['particle_index'];other=rows1[row['contact_id']]
        audit=hht_path.contact_audit(row['normal'],np.zeros(3),np.zeros(3),
            row['wall_velocity_m_s'],other['wall_velocity_m_s'],velocities(old)[idx],v1[idx],
            state['acceleration_m_s2'][idx],acceleration1[idx],dt,0.,0.,p['hht_weight'])
        if audit['interior_zero_fractions']:
            incoming=audit['slip_coefficients_m_s'][0]
            if np.linalg.norm(incoming)==0:raise ValueError('already-zero parent needs explicit right restart')
            found.append((min(audit['interior_zero_fractions']),row['contact_id'],incoming/np.linalg.norm(incoming)))
    if not found:raise ValueError('failed path has no resolvable full-vector root predictor')
    earliest=min(item[0] for item in found)
    # Candidate grouping is only a search convenience. ALL full-vector endpoint
    # residuals must independently meet roundoff; this never imposes simultaneity.
    return {cid:dict(endpoint='private_continuation',incoming_direction=direction,
        predictor_dt_s=fraction*dt)
        for fraction,cid,direction in found if abs(fraction-earliest)*dt<=1e-12}


def locate_duration(evaluate, old, parent, context, upper_dt, maximum_trials=12):
    """Bounded safeguarded duration root, each call a complete coupled solve."""
    if not context or maximum_trials<2 or maximum_trials>12:raise ValueError('bounded nonempty event root required')
    rows={row['contact_id']:row for row in parent['boundary_snapshot']['contacts']}
    if set(context)-set(rows):raise ValueError('event contact identity missing')
    def values(pa,state):
        current={r['contact_id']:r for r in state['boundary_snapshot']['contacts']};result=[]
        for cid,entry in context.items():
            row=current[cid];n=row['normal'];u=velocities(pa)[row['particle_index']]-row['wall_velocity_m_s']
            u=u-(u@n)*n;result.append(u)
        signed=[u@entry['incoming_direction'] for u,entry in zip(result,context.values())]
        return np.array(result),float(np.mean(signed))
    initial,g0=values(old,parent)
    if g0<=0:raise ValueError('strict incoming sliding parent required')
    low=0.;high=float(upper_dt);glo=g0;ghi=None;duration=high;attempts=[]
    for iteration in range(maximum_trials):
        pa,state,diag=evaluate(duration)
        raw,g=values(pa,state)
        scale=max(1.,float(np.max(np.linalg.norm(velocities(pa),axis=1))),
            float(np.max(np.linalg.norm(velocities(old),axis=1))),
            float(duration*np.max(np.linalg.norm(state['acceleration_m_s2'],axis=1))),
            float(duration*np.max(np.linalg.norm(parent['acceleration_m_s2'],axis=1))),
            max(float(np.linalg.norm(rows[cid]['wall_velocity_m_s'])) for cid in context))
        bound=512*np.finfo(float).eps*scale
        attempts.append(dict(iteration=iteration+1,dt_s=duration,signed_root_m_s=g,
            full_tangent_m_s=raw.copy(),zero_bound_m_s=bound,
            coupled_residuals=copy.deepcopy(diag['residuals']),inner_evaluations=diag['evaluations'],
            candidate_state_hash=state['state_hash'],accepted=False,restart_allowed=False))
        if np.max(np.linalg.norm(raw,axis=1))<=bound:
            return pa,state,diag,dict(status='LOCATED_OWN_TRAJECTORY_FULL_VECTOR_ROOT',
                attempts=attempts,time_s=state['time_s'],dt_s=duration,
                bracket_s=[low,high],zero_bound_m_s=bound,raw_full_tangent_m_s=raw,
                physical_state_was_projected=False,reference_time_was_used=False)
        if g>0:low,glo=duration,g
        else:high,ghi=duration,g
        if ghi is None:
            error=ValueError('upper coupled trial did not bracket the event')
            error.diagnostics=dict(attempts=attempts);raise error
        trial=(low*ghi-high*glo)/(ghi-glo)
        if not low<trial<high or trial==duration:trial=.5*(low+high)
        duration=trial
    error=ValueError('bounded duration root did not certify every full tangent')
    error.diagnostics=dict(attempts=attempts,bracket_s=[low,high]);raise error


def restart_sliding(pa,left,p,incoming_context,right_boundary,free_axes,wall_acceleration):
    """Explicit zero-duration force/history change; all particle arrays preserved."""
    stage._validate_state(pa,left,p);before=particle_hash(pa)
    if not incoming_context:raise ValueError('no located event: acceleration reset prohibited')
    internal,spatial=_force(_detach(pa),p);rows=left['boundary_snapshot']['contacts']
    limits={};outgoing={};vv=velocities(pa)
    axes=np.asarray(free_axes)
    if axes.shape!=(len(pa.x),3):raise ValueError('actual complete free-axis mask required')
    for cid,entry in incoming_context.items():
        row=next(r for r in rows if r['contact_id']==cid);idx=row['particle_index'];n=row['normal']
        relative=vv[idx]-row['wall_velocity_m_s'];slip=relative-(relative@n)*n
        other=internal[idx]+left['boundary_snapshot']['noncontact_force_N'][idx]
        value=zero_event_limits(pa.m[idx],n,row['normal_force_N'],p['friction'],other,
            wall_acceleration[cid],slip,max(1.,np.linalg.norm(vv[idx]),
            np.linalg.norm(row['wall_velocity_m_s']),left['last_dt_s']*np.linalg.norm(left['acceleration_m_s2'][idx])),
            entry['incoming_direction'],other_contact_count=sum(r['particle_index']==idx for r in rows)-1,
            free_axes=axes[idx])
        if value['stick_feasible'] or value['feasibility_margin_is_roundoff_marginal']:
            raise ValueError('general3D sticking/marginal feasibility not implemented')
        limits[cid]=value;outgoing[cid]=dict(endpoint='right_zero',
            incoming_direction=entry['incoming_direction'],outgoing_direction=value['outgoing_direction'])
    callback=right_boundary(outgoing)
    right=_boundary(pa,callback,left['time_s'],p)
    if [r['contact_id'] for r in rows]!=[r['contact_id'] for r in right['contacts']]:
        raise ValueError('unregistered contact topology change at zero event')
    if not np.array_equal(left['boundary_snapshot']['noncontact_force_N'],right['noncontact_force_N']):
        raise ValueError('noncontact force change at slip event')
    for before_row,after_row in zip(rows,right['contacts']):
        for key in before_row:
            if key=='tangential_force_N' and before_row['contact_id'] in limits:continue
            if not np.array_equal(before_row[key],after_row[key]):
                raise ValueError('unregistered boundary change '+key)
    for row in right['contacts']:
        cid=row['contact_id']
        if cid in limits and not np.array_equal(row['tangential_force_N'],limits[cid]['right_tangential_force_N']):
            raise ValueError('right boundary differs from static-feasibility force limit')
    force=internal+right['total_boundary_force_N'];acc=force/pa.m[:,None]
    storage=acceleration_storage_jump(left['acceleration_m_s2'],acc,pa.m,left['last_dt_s'],p.get('hht_weight',.9))
    event=dict(kind='LOCATED_ZERO_SLIP_RIGHT_FORCE_LIMIT_RESTART',time_s=left['time_s'],
        limits=limits,old_canonical_acceleration_m_s2=left['acceleration_m_s2'].copy(),
        new_canonical_acceleration_m_s2=acc.copy(),old_force_N=left['force_N'].copy(),new_force_N=force.copy(),
        left_boundary_snapshot=copy.deepcopy(left['boundary_snapshot']),right_boundary_snapshot=right,
        physical_state_hash_before=before,physical_state_hash_after=particle_hash(pa),
        physical_impulse_Ns=np.zeros_like(force),physical_work_J=0.,physical_heat_J=0.,
        signed_acceleration_storage=storage,algorithmic_terms_are_material_heat=False,spatial=spatial)
    state=stage._state(pa,p,left['time_s'],acc,force,right,parent=left['state_hash'],
        dt=left['last_dt_s'],event=event,interval=left['last_interval'])
    if before!=particle_hash(pa):raise RuntimeError('event restart mutated physical fields')
    return state,outgoing,event

"""Detached curved contact with explicitly located one-sided zero endpoints.

This new module leaves the frozen event-free helper untouched. It reuses its
quadrature rule, penalty/boundary validation and physical HHT dense path.
A PRIVATE incoming-branch extension is mathematical root-search data only:
negative local raw work is retained and no such result can be published.
No particle velocity, coordinate, material field or heat is corrected here.
"""
import copy
import numpy as np

from . import hht_path
from .hht_coupled import _boundary, _detach, _hash
from .mechanics import positions, velocities
from .updated_trial import particle_hash
from .hht_curved_contact import CurvedContactError, refined_gauss
from .hht_zero_event import one_sided_slip_direction

MODES = ('private_continuation','left_zero','right_zero')


def _event_direction(value, normal, name):
    direction = np.array(value,dtype=float,copy=True)
    eps = np.finfo(float).eps
    if (direction.shape!=(3,) or not np.isfinite(direction).all() or
            abs(np.linalg.norm(direction)-1)>512*eps or abs(direction@normal)>512*eps):
        raise ValueError('unit tangential '+name+' required')
    return direction


def _event_path_certificate(audit,context,normal,velocity_scale):
    """Full polynomial check; no endpoint or coefficient is projected.

left_zero means the interval ENDS at an incoming zero (s=1).
right_zero means it STARTS at an outgoing zero (s=0) and requires the separate
outgoing_direction. The incoming_direction field always retains its meaning.
Private extension is allowed only for naturally collinear target coefficients;
this is not a frozen-direction model for a face or ordinary turning contact.
"""
    if not isinstance(context,dict) or context.get('endpoint') not in MODES:
        raise ValueError('explicit registered event endpoint mode required')
    if 'incoming_direction' not in context:
        raise ValueError('incoming one-sided direction required')
    if not np.isfinite(velocity_scale):
        raise ValueError('finite kinematic event scale required')
    mode = context['endpoint']
    incoming = _event_direction(context['incoming_direction'],normal,'incoming_direction')
    coefficients = np.array(audit['slip_coefficients_m_s'],copy=True)
    coefficient_scale = float(np.linalg.norm(coefficients,axis=1).sum())
    if not np.isfinite(coefficient_scale):
        raise ValueError('finite polynomial event scale required')
    eps = np.finfo(float).eps
    bound = 512*eps*max(1.,velocity_scale,coefficient_scale)
    result = dict(mode=mode,incoming_direction=incoming,zero_slip_bound_m_s=bound,
        raw_slip_coefficients_m_s=coefficients,velocity_was_modified=False,
        whole_interval_no_other_zero_certified=False,accepted=False,restart_allowed=False)
    if mode=='private_continuation':
        transverse=coefficients-np.outer(coefficients@incoming,incoming)
        if np.max(np.linalg.norm(transverse,axis=1))>bound:
            raise ValueError('private event extension requires naturally collinear full-vector slip')
        if coefficients[0]@incoming<=bound:
            raise ValueError('private incoming continuation must start on resolved incoming branch')
        result.update(private_continuation=True,physical_contact_path_certified=False,
            transverse_coefficients_m_s=transverse,
            maximum_transverse_coefficient_m_s=float(np.max(np.linalg.norm(transverse,axis=1))))
        return result
    endpoint=1 if mode=='left_zero' else 0
    if mode=='right_zero' and 'outgoing_direction' not in context:
        raise ValueError('right_zero requires separate force-selected outgoing_direction')
    selected=incoming if endpoint==1 else _event_direction(context['outgoing_direction'],normal,'outgoing_direction')
    # The standalone algebra guard may be stricter than the velocity-scaled
    # caller bound; it is never enlarged to obtain a desired endpoint.
    direction=one_sided_slip_direction(coefficients,endpoint,
        tolerance=min(bound,512*eps*max(1.,coefficient_scale)))
    if np.linalg.norm(direction['direction']-selected)>512*eps:
        raise ValueError('actual dense one-sided direction differs from selected physical branch')
    if audit['all_interval_stick']:
        raise ValueError('whole-interval stick needs a separate constrained solve')
    other=coefficients[0] if endpoint==1 else coefficients.sum(axis=0)
    if np.linalg.norm(other)<=bound:
        raise ValueError('both endpoints at zero require an additional located event')
    leading_norm=float(np.linalg.norm(direction['leading_coefficient']))
    # Root uncertainty from the raw endpoint roundoff. This is only an alias
    # enclosure on THIS polynomial, not an accuracy claim about the ODE event.
    alias=(4*bound/leading_norm)**(1./direction['leading_order'])
    if not np.isfinite(alias) or alias>=1.:
        raise ValueError('event endpoint and entire interval are numerically unresolved')
    u=coefficients
    squared=np.array([u[0]@u[0],2*u[0]@u[1],u[1]@u[1]+2*u[0]@u[2],
                      2*u[1]@u[2],u[2]@u[2]])
    candidates=list(audit['interior_zero_fractions'])
    candidates+=hht_path._roots([squared[1],2*squared[2],3*squared[3],4*squared[4]])
    for axis in range(3):
        candidates+=hht_path._roots(u[:,axis])
        if u[2,axis]!=0:
            vertex=-u[1,axis]/(2*u[2,axis])
            if 0<vertex<1:candidates.append(float(vertex))
    aliases=[]
    for fraction in sorted(set(candidates)):
        raw=u[0]+fraction*u[1]+fraction*fraction*u[2]
        if np.linalg.norm(raw)<=bound:
            if abs(fraction-endpoint)>alias:
                raise ValueError('resolvable earlier/interior full-vector zero requires another event')
            aliases.append(dict(fraction=float(fraction),raw_slip_m_s=raw,
                classification='ROUNDING_ENDPOINT_ALIAS_NOT_A_RESOLVED_SECOND_EVENT'))
    result.update(private_continuation=False,physical_contact_path_certified=True,
        one_sided_direction=direction['direction'],one_sided_limit=direction,
        event_endpoint_fraction=endpoint,endpoint_alias_bound_fraction=float(alias),
        endpoint_roundoff_aliases=aliases,whole_interval_no_other_zero_certified=True,
        outgoing_direction=None if endpoint==1 else selected)
    return result


def integrate(pa, p, boundary, t0, acc0, acc1, dt, weight=.9,
              old_boundary=None, new_boundary=None, gap_certificate=None, event_context=None):
    """Detached interval impulse/work/Q, with explicit event_context exceptions.

Context: contact_id -> {incoming_direction, endpoint}. Allowed endpoint values
are private_continuation, left_zero (finish at zero), right_zero (start at zero).
right_zero additionally requires outgoing_direction selected by actual force
feasibility. The first mode is REJECTED root-search data, never a checkpoint.
Negative local private raw Q is preserved; negative integrated Q still rejects.
Strict endpoint modes use one-sided force only at the registered zero endpoint;
ordinary quadrature nodes retain full-vector pointwise Coulomb unchanged.

gap_certificate(pa,a0,a1,dt,p,left,right) must return global positive
minimum_overlap_m, finite maximum_overlap_m, fixed_normals=True and
affine_tangential_wall_velocity=True. It certifies the whole interval using
the explicitly registered wall programme. The calling driver additionally
enforces its physical maximum-penetration domain; no universal5um limit is
imposed on independently manufactured SI unit benchmarks here.

All actual/delta impulse and work arrays are per particle. Deltas compare
the SAME curved contact quantities with total endpoint-linear contact, for
the coupled caller's impulse-consistent residual and energy identity.
"""
    if not np.isfinite(t0) or not np.isfinite(dt) or t0 < 0 or dt <= 0:
        raise CurvedContactError('finite nonnegative start and positive interval required')
    hht_path.parameters(weight)
    if event_context is not None and not isinstance(event_context, dict):
        raise CurvedContactError('event_context must be a contact-ID mapping')
    context_before = _hash(event_context); contexts = copy.deepcopy(event_context or {})
    event_records = {}
    x0 = positions(pa); v0 = velocities(pa)
    a0 = np.array(acc0, float, copy=True); a1 = np.array(acc1, float, copy=True)
    if (a0.shape != v0.shape or a1.shape != v0.shape or
            not np.isfinite(a0).all() or not np.isfinite(a1).all()):
        raise CurvedContactError('matching finite candidate accelerations required')
    x1, v1 = hht_path.endpoint(x0, v0, a0, a1, dt, weight)
    count = len(pa.x); before = particle_hash(pa); parameter_before = _hash(p)
    cache = {}; identities = None; left = None; right = None; certificate = None
    direction_limit = float(p.get('contact_coulomb_relative_tolerance', .01))
    if not np.isfinite(direction_limit) or not 0 <= direction_limit <= .01:
        raise CurvedContactError('unchanged Coulomb defect acceptance bound at most0.01 required')

    def sample_boundary(s):
        if s in cache: return cache[s]
        sample = _detach(pa)
        if s == 0.: x, v = x0.copy(), v0.copy()
        elif s == 1.: x, v = x1.copy(), v1.copy()
        else: x, v = hht_path.sample(x0, v0, a0, a1, dt, s, weight)
        sample.x[:], sample.y[:], sample.z[:] = x.T
        sample.u[:], sample.v[:], sample.w[:] = v.T
        measured = _boundary(sample, boundary, t0+s*dt, p)
        result = dict(boundary=measured, position_m=x, velocity_m_s=v)
        cache[s] = result
        return result

    def validate_rows(result, s):
        rows = result['boundary']['contacts']
        if [row['contact_id'] for row in rows] != identities:
            raise CurvedContactError('contact identity/feature switch requires located event')
        defects = []; computed_tractions = []
        for index, row in enumerate(rows):
            ref = left['contacts'][index]; end = right['contacts'][index]
            j = row['particle_index']; n = row['normal']
            if (j != ref['particle_index'] or row['particle_id'] != ref['particle_id'] or
                    not np.array_equal(n, ref['normal']) or
                    row['normal_stiffness_N_m'] != ref['normal_stiffness_N_m']):
                raise CurvedContactError('fixed-normal stable-identity penalty scope violated')
            pen = row['penetration_m']
            geometry_scale = max(float(np.max(np.abs(result['position_m']))), 1e-12)
            geometric_roundoff = 512*np.finfo(float).eps*geometry_scale
            if pen <= 0 or pen < certificate['minimum_overlap_m']-geometric_roundoff or pen > certificate['maximum_overlap_m']+geometric_roundoff:
                raise CurvedContactError('actual overlap outside positive whole-path certificate')
            project = lambda vector: vector-np.dot(vector, n)*n
            w = row['wall_velocity_m_s']; expected = (1-s)*project(ref['wall_velocity_m_s'])+s*project(end['wall_velocity_m_s'])
            wallscale = max(np.linalg.norm(w), np.linalg.norm(expected), 1e-12)
            if np.linalg.norm(project(w)-expected) > 512*np.finfo(float).eps*wallscale:
                raise CurvedContactError('nonaffine tangential wall programme outside zero-slip certificate')
            u = project(result['velocity_m_s'][j]-w); speed = float(np.linalg.norm(u))
            event = event_records.get(row['contact_id'])
            mode = None if event is None else event['mode']
            at_limit = ((mode == 'left_zero' and s == 1.) or
                        (mode == 'right_zero' and s == 0.))
            cap = p['friction']*row['normal_force_N']; ft = row['tangential_force_N']
            if mode == 'private_continuation':
                if speed <= event['zero_slip_bound_m_s']:
                    expected_force = -cap*event['incoming_direction']
                elif u@event['incoming_direction'] >= 0.:
                    expected_force = -cap*u/speed
                else:
                    expected_force = cap*u/speed
            elif at_limit:
                if speed > event['zero_slip_bound_m_s']:
                    raise CurvedContactError('one-sided force requires actual full-vector endpoint zero',
                        dict(contact_id=row['contact_id'],raw_slip_m_s=u.copy(),accepted=False))
                expected_force = -cap*event['one_sided_direction']
            elif speed == 0:
                raise CurvedContactError('unregistered full-vector zero slip requires located event')
            else:
                expected_force = -cap*u/speed
            defect = float(np.linalg.norm(ft-expected_force)/max(cap, np.finfo(float).tiny))
            if defect > direction_limit+512*np.finfo(float).eps:
                raise CurvedContactError('actual curved force violates original Coulomb law bound',
                    dict(contact_id=row['contact_id'], fraction=s, defect=defect, bound=direction_limit))
            force_roundoff = 512*np.finfo(float).eps*max(cap,np.linalg.norm(ft),np.finfo(float).tiny)
            if np.linalg.norm(ft-expected_force) > force_roundoff:
                raise CurvedContactError('callback not roundoff-equal to declared pointwise Coulomb law',
                    dict(contact_id=row['contact_id'], fraction=s, callback_force_N=ft.copy(),
                         pointwise_force_N=expected_force.copy(),force_roundoff_bound_N=force_roundoff))
            # The actual interval force is the declared pointwise law, NOT a
            # sampled opaque callback interpolation. Inconsistent callbacks
            # reject above; they are never silently overwritten to make a pass.
            raw_power = float(expected_force@u)
            power_bound = 256*np.finfo(float).eps*max(cap, np.linalg.norm(expected_force))*speed
            if at_limit:
                # A measure-zero endpoint force is a one-sided limit. Preserve
                # its raw residual power; do not turn that rounding into heat.
                power_bound += cap*event['zero_slip_bound_m_s']
            if raw_power > power_bound and mode != 'private_continuation':
                raise CurvedContactError('negative raw curved friction dissipation',
                    dict(contact_id=row['contact_id'], fraction=s, friction_power_W=raw_power,
                         power_rounding_bound_W=power_bound))
            defects.append(defect)
            computed_tractions.append(expected_force)
        result['coulomb_defects'] = defects
        result['pointwise_coulomb_forces_N'] = computed_tractions
        return rows

    try:
        left = sample_boundary(0.)['boundary']; right = sample_boundary(1.)['boundary']
        if old_boundary is not None and _hash(old_boundary) != _hash(left):
            raise CurvedContactError('curved interval changed accepted left boundary history')
        if new_boundary is not None and _hash(new_boundary) != _hash(right):
            raise CurvedContactError('curved interval endpoint differs from candidate boundary')
        identities = [row['contact_id'] for row in left['contacts']]
        if [row['contact_id'] for row in right['contacts']] != identities:
            raise CurvedContactError('contact endpoint topology requires located event')
        if identities:
            if gap_certificate is None:
                raise CurvedContactError('source-bound whole-path gap/kinematic certificate required')
            private = _detach(pa); private_hash = particle_hash(private)
            certificate = gap_certificate(private, a0.copy(), a1.copy(), dt,
                copy.deepcopy(p), copy.deepcopy(left), copy.deepcopy(right))
            if particle_hash(private) != private_hash:
                raise CurvedContactError('gap certificate mutated particle state')
            required = ('minimum_overlap_m', 'maximum_overlap_m', 'fixed_normals', 'affine_tangential_wall_velocity')
            if not isinstance(certificate, dict) or any(key not in certificate for key in required):
                raise CurvedContactError('explicit whole-path certificate fields required')
            low, high = certificate['minimum_overlap_m'], certificate['maximum_overlap_m']
            if (certificate['fixed_normals'] is not True or certificate['affine_tangential_wall_velocity'] is not True or
                    not np.isfinite(low) or not np.isfinite(high) or not 0 < low <= high):
                raise CurvedContactError('unsupported or nonpositive whole-path contact certificate')
        else:
            certificate = dict(minimum_overlap_m=None, maximum_overlap_m=None,
                fixed_normals=True, affine_tangential_wall_velocity=True, empty_contact_set=True)
        if any(not isinstance(cid, str) or cid not in identities for cid in contexts):
            raise CurvedContactError('event context names an absent or invalid contact identity')
        path_audits = []
        for row0, row1 in zip(left['contacts'], right['contacts']):
            j = row0['particle_index']
            # Only the polynomial full-vector zero-slip certificate is reused.
            # Dummy zero traction prevents the obsolete endpoint-linear force
            # path from vetoing valid curved sliding. No zero-force physics is
            # returned or used in any integral.
            audit = hht_path.contact_audit(row0['normal'], np.zeros(3), np.zeros(3),
                row0['wall_velocity_m_s'], row1['wall_velocity_m_s'],v0[j],v1[j],a0[j],a1[j],dt,0.,0.,weight)
            context = contexts.get(row0['contact_id'])
            if context is None:
                if (audit['requires_event_split'] or audit['all_interval_stick'] or
                        audit['minimum_checked_slip_m_s'] <= audit['zero_speed_bound_m_s']):
                    raise CurvedContactError('dense full-vector zero slip requires located event',
                        dict(contact_id=row0['contact_id'],zero_slip_audit=audit,accepted=False))
            else:
                scale = max(np.linalg.norm(v0[j]),np.linalg.norm(v1[j]),
                    np.linalg.norm(row0['wall_velocity_m_s']),np.linalg.norm(row1['wall_velocity_m_s']),
                    dt*np.linalg.norm(a0[j]),dt*np.linalg.norm(a1[j]))
                try:
                    event_records[row0['contact_id']] = _event_path_certificate(
                        audit,context,row0['normal'],float(scale))
                except ValueError as error:
                    raise CurvedContactError(str(error),dict(contact_id=row0['contact_id'],
                        zero_slip_audit=audit,event_context=copy.deepcopy(context),accepted=False)) from error
            path_audits.append({key:copy.deepcopy(audit[key]) for key in (
                'requires_event_split','interior_zero_fractions','all_interval_stick',
                'minimum_checked_slip_m_s','zero_speed_bound_m_s','slip_coefficients_m_s')})

        validate_rows(sample_boundary(0.), 0.); validate_rows(sample_boundary(1.), 1.)
        is_private = any(row['mode']=='private_continuation' for row in event_records.values())

        def integrand(s):
            result = sample_boundary(s); rows = validate_rows(result, s); values = []
            for ci,row in enumerate(rows):
                j = row['particle_index']; material = result['velocity_m_s'][j]
                wall = row['wall_velocity_m_s']; fn = row['normal_force_N']*row['normal']
                ft = result['pointwise_coulomb_forces_N'][ci]
                relative = material-wall
                tangent = relative-np.dot(relative,row['normal'])*row['normal']
                values.append(np.r_[fn, ft, fn@material, ft@material, fn@wall, ft@wall,
                                    -ft@tangent])
            return np.asarray(values, float).reshape(len(identities), 11)

        quadrature = refined_gauss(integrand, dt)
        integral = quadrature['integral']; contacts = []
        Jn = np.zeros((count,3)); Jt = np.zeros_like(Jn)
        Wn = np.zeros(count); Wt = np.zeros(count); Wnw = np.zeros(count); Wtw = np.zeros(count)
        Q = np.zeros(count); roundbound = np.zeros(count); wallrow = np.zeros(count); bath = np.zeros(count)
        for ci, (row0, row1) in enumerate(zip(left['contacts'], right['contacts'])):
            j = row0['particle_index']; value = integral[ci]
            Jn[j] += value[:3]; Jt[j] += value[3:6]
            Wn[j] += value[6]; Wt[j] += value[7]; Wnw[j] += value[8]; Wtw[j] += value[9]; Q[j] += value[10]
            wallrow[j] += row1['conductance_W_K']; bath[j] += row1['conductance_W_K']*row1['wall_temperature_K']
            nodes = []; local_rounding = 0.
            for s, w, values in zip(quadrature['nodes'],quadrature['weights'],quadrature['values']):
                sample = cache[float(s)]; row = sample['boundary']['contacts'][ci]
                material = sample['velocity_m_s'][j]; wall = row['wall_velocity_m_s']; n = row['normal']
                u = material-wall; tangent = u-np.dot(u,n)*n; cap=p['friction']*row['normal_force_N']
                actual_ft = sample['pointwise_coulomb_forces_N'][ci]
                rounding = dt*w*256*np.finfo(float).eps*max(cap,np.linalg.norm(actual_ft))*np.linalg.norm(tangent)
                local_rounding += rounding
                nv = dt*w*values[ci]
                nodes.append(dict(fraction=float(s),weight=float(w),normal_impulse_Ns=nv[:3],
                    tangential_impulse_Ns=nv[3:6],impulse_Ns=nv[3:6],actual_total_impulse_Ns=nv[:3]+nv[3:6],
                    velocity_m_s=material.copy(),wall_velocity_m_s=wall.copy(),relative_velocity_m_s=u,
                    tangential_relative_velocity_m_s=tangent,normal_material_work_J=float(nv[6]),
                    tangent_material_work_J=float(nv[7]),normal_wall_work_J=float(nv[8]),
                    tangent_wall_work_J=float(nv[9]),friction_heat_raw_J=float(nv[10]),
                    friction_rounding_bound_J=float(rounding),penetration_m=row['penetration_m'],
                    normal_force_N=row['normal_force_N'],tangential_force_N=actual_ft.copy(),
                    callback_tangential_force_N=row['tangential_force_N'].copy(),
                    callback_force_roundoff_defect_N=float(np.linalg.norm(actual_ft-row['tangential_force_N']))))
            roundbound[j] += local_rounding
            if value[10] < -local_rounding:
                raise CurvedContactError('negative integrated raw friction heat',
                    dict(contact_id=row0['contact_id'],friction_heat_raw_J=float(value[10]),
                         friction_rounding_bound_J=float(local_rounding),quadrature=nodes,
                         private_continuation=is_private,accepted=False,restart_allowed=False))
            maximum_defect = max(sample['coulomb_defects'][ci] for sample in cache.values())
            contacts.append(dict(path_audits[ci],contact_id=row0['contact_id'],particle_id=row0['particle_id'],
                particle_index=j,normal=row0['normal'].copy(),actual_normal_impulse_Ns=value[:3].copy(),
                actual_tangential_impulse_Ns=value[3:6].copy(),total_actual_impulse_Ns=(value[:3]+value[3:6]).copy(),
                actual_material_work_J=float(value[6]+value[7]),actual_wall_work_J=float(value[8]+value[9]),
                normal_material_work_J=float(value[6]),normal_relative_work_J=float(value[6]-value[8]),
                friction_heat_raw_J=float(value[10]),friction_rounding_bound_J=float(local_rounding),
                maximum_sampled_coulomb_relative_defect=maximum_defect,coulomb_relative_error_bound=direction_limit,
                whole_path_passive=not (row0['contact_id'] in event_records and
                    event_records[row0['contact_id']]['mode']=='private_continuation'),
                passivity_scope=('PRIVATE_BRANCH_EXTENSION_NOT_PHYSICAL_OR_PUBLISHABLE' if
                    row0['contact_id'] in event_records and event_records[row0['contact_id']]['mode']=='private_continuation'
                    else 'POINTWISE_COULOMB_OPEN_INTERVAL_WITH_CERTIFIED_ONE_SIDED_ZERO_ENDPOINT_IF_PRESENT'),
                event_path_certificate=copy.deepcopy(event_records.get(row0['contact_id'])),
                event_mode=contexts.get(row0['contact_id'],{}).get('endpoint'),
                restart_allowed=False,
                interior_coulomb_direction=('PRIVATE_INCOMING_BRANCH_NOT_PHYSICAL' if
                    row0['contact_id'] in event_records and event_records[row0['contact_id']]['mode']=='private_continuation'
                    else 'EXACT_POINTWISE_COULOMB_ON_OPEN_DENSE_PATH_CALLBACK_ROUNDOFF_CONSISTENT'),
                quadrature=nodes,quadrature_order=quadrature['order'],accepted=False))
        impulse = Jn+Jt; material_work = Wn+Wt; wall_work = Wnw+Wtw
        endpoint_forces = []
        for at in (left,right):
            force=np.zeros((count,3))
            for row in at['contacts']:
                force[row['particle_index']] += row['normal_force_N']*row['normal']+row['tangential_force_N']
            endpoint_forces.append(force)
        linear=hht_path.linear_work(endpoint_forces[0],endpoint_forces[1],v0,v1,a0,a1,dt,weight)
        return dict(classification=('REJECTED_PRIVATE_INCOMING_BRANCH_CONTACT_QUADRATURE' if is_private
                else 'CURVED_CONTACT_QUADRATURE_WITH_CERTIFIED_ZERO_ENDPOINTS'),
            accepted=False,restart_allowed=False,can_publish=False,private_continuation=is_private,
            physical_contact_path_certified=not is_private,whole_path_passive=not is_private,
            event_context=copy.deepcopy(contexts),event_path_certificates=copy.deepcopy(event_records),
            actual_impulse_Ns=impulse,actual_normal_impulse_Ns=Jn,actual_tangential_impulse_Ns=Jt,
            actual_material_work_J=material_work,actual_wall_work_J=wall_work,
            actual_normal_material_work_J=Wn,actual_tangent_material_work_J=Wt,
            actual_normal_wall_work_J=Wnw,actual_tangent_wall_work_J=Wtw,
            delta_impulse_Ns=impulse-linear['impulse_Ns'],delta_material_work_J=material_work-linear['work_J'],
            endpoint_linear_contact_impulse_Ns=linear['impulse_Ns'],endpoint_linear_contact_material_work_J=linear['work_J'],
            friction_raw_J=Q,friction_rounding_bound_J=roundbound,
            wall_row_W_K=wallrow,wall_weighted_temperature_W=bath,contacts=contacts,
            wall_work_J=float(wall_work.sum()),normal_material_work_J=float(Wn.sum()),
            normal_relative_work_J=float((Wn-Wnw).sum()),
            normal_potential_change_J=sum(r['normal_potential_J'] for r in right['contacts'])-sum(r['normal_potential_J'] for r in left['contacts']),
            contact_work_identity_residual_J=material_work-wall_work+Q-(Wn-Wnw),
            maximum_sampled_coulomb_relative_defect=max((c['maximum_sampled_coulomb_relative_defect'] for c in contacts),default=0.),
            quadrature=dict(order=quadrature['order'],attempted_orders=quadrature['attempted_orders'],
                difference=quadrature['difference'],tolerance=quadrature['tolerance'],error_status=quadrature['error_status']),
            boundary_evaluations=len(cache),gap_certificate=copy.deepcopy(certificate),
            wall_thermal_policy='BE_ENDPOINT_RESERVOIR_FIELDS',algorithmic_terms_are_material_heat=False)
    finally:
        if particle_hash(pa)!=before or _hash(p)!=parameter_before or _hash(event_context)!=context_before:
            raise CurvedContactError('stop contact callback/certificate changed accepted particles, parameters or event context')



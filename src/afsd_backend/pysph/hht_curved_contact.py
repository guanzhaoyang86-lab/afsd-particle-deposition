"""Event-free curved Coulomb contact quadrature on the existing HHT path.

Only fixed-normal, stable-identity contacts with certified affine tangential
wall velocity are supported. The caller supplies a source-bound whole-path
gap certificate; quadrature samples alone do not certify contact topology.
No state is committed, velocity projected, negative heat clipped, or old
endpoint-linear friction work used in place of the actual curved work.
"""
import copy
import numpy as np

from . import hht_path
from .hht_coupled import _boundary, _detach, _hash
from .mechanics import positions, velocities
from .updated_trial import particle_hash

ABSOLUTE_QUADRATURE_BOUND = 1e-14
RELATIVE_QUADRATURE_BOUND = 1e-8
QUADRATURE_ORDERS = (8, 16, 32)


class CurvedContactError(ValueError):
    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = {} if diagnostics is None else diagnostics


def refined_gauss(evaluator, dt):
    """Integrate a finite vector integrand in physical time; no physics implied.

Independent 8/16 and, if needed, 16/32 Gauss-Legendre differences must meet
1e-14 absolute plus1e-8 relative PER COMPONENT. These differences are an
estimated quadrature error, not a universal rigorous bound for arbitrary
integrands. This generic primitive deliberately does not claim an arbitrary
prescribed circular path is an HHT path or a free-material trajectory.
"""
    if not np.isfinite(dt) or dt <= 0:
        raise CurvedContactError('positive finite quadrature interval required')
    previous = None; expected_shape = None; attempts = []
    for order in QUADRATURE_ORDERS:
        abscissae, weights = np.polynomial.legendre.leggauss(order)
        nodes = .5*(abscissae+1.); weights = .5*weights
        values = []
        for fraction in nodes:
            value = np.asarray(evaluator(float(fraction)), float)
            if expected_shape is None: expected_shape = value.shape
            if value.shape != expected_shape or not np.isfinite(value).all():
                raise CurvedContactError('finite fixed-shape quadrature integrand required')
            values.append(value.copy())
        values = np.asarray(values)
        integral = dt*np.tensordot(weights, values, axes=(0, 0))
        attempts.append(order)
        if previous is not None:
            difference = integral-previous
            tolerance = ABSOLUTE_QUADRATURE_BOUND+RELATIVE_QUADRATURE_BOUND*np.maximum(
                np.abs(integral), np.abs(previous))
            if np.all(np.abs(difference) <= tolerance):
                return dict(integral=integral, order=order, nodes=nodes, weights=weights,
                    values=values, difference=difference, tolerance=tolerance,
                    attempted_orders=attempts, error_status='REFINEMENT_ESTIMATE_NOT_RIGOROUS_ARBITRARY_FUNCTION_BOUND')
        previous = integral
    raise CurvedContactError('curved contact quadrature did not converge by32nodes',
        dict(attempted_orders=attempts, difference=difference, tolerance=tolerance,
             last_integral=integral, accepted=False))


def integrate(pa, p, boundary, t0, acc0, acc1, dt, weight=.9,
              old_boundary=None, new_boundary=None, gap_certificate=None):
    """Actual contact impulse/work/Q on a detached HHT candidate dense path.

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
            if speed == 0:
                raise CurvedContactError('full-vector zero slip requires located event')
            cap = p['friction']*row['normal_force_N']; ft = row['tangential_force_N']
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
            if raw_power > power_bound:
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
        validate_rows(sample_boundary(0.), 0.); validate_rows(sample_boundary(1.), 1.)
        path_audits = []
        for row0, row1 in zip(left['contacts'], right['contacts']):
            j = row0['particle_index']
            # Only the polynomial full-vector zero-slip certificate is reused.
            # Dummy zero traction prevents the obsolete endpoint-linear force
            # path from vetoing valid curved sliding. No zero-force physics is
            # returned or used in any integral.
            audit = hht_path.contact_audit(row0['normal'], np.zeros(3), np.zeros(3),
                row0['wall_velocity_m_s'], row1['wall_velocity_m_s'],v0[j],v1[j],a0[j],a1[j],dt,0.,0.,weight)
            if (audit['requires_event_split'] or audit['all_interval_stick'] or
                    audit['minimum_checked_slip_m_s'] <= audit['zero_speed_bound_m_s']):
                raise CurvedContactError('dense full-vector zero slip requires located event',
                    dict(contact_id=row0['contact_id'], zero_slip_audit=audit, accepted=False))
            path_audits.append({key:copy.deepcopy(audit[key]) for key in (
                'requires_event_split','interior_zero_fractions','all_interval_stick',
                'minimum_checked_slip_m_s','zero_speed_bound_m_s','slip_coefficients_m_s')})

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
                raise CurvedContactError('negative integrated raw friction heat')
            maximum_defect = max(sample['coulomb_defects'][ci] for sample in cache.values())
            contacts.append(dict(path_audits[ci],contact_id=row0['contact_id'],particle_id=row0['particle_id'],
                particle_index=j,normal=row0['normal'].copy(),actual_normal_impulse_Ns=value[:3].copy(),
                actual_tangential_impulse_Ns=value[3:6].copy(),total_actual_impulse_Ns=(value[:3]+value[3:6]).copy(),
                actual_material_work_J=float(value[6]+value[7]),actual_wall_work_J=float(value[8]+value[9]),
                normal_material_work_J=float(value[6]),normal_relative_work_J=float(value[6]-value[8]),
                friction_heat_raw_J=float(value[10]),friction_rounding_bound_J=float(local_rounding),
                maximum_sampled_coulomb_relative_defect=maximum_defect,coulomb_relative_error_bound=direction_limit,
                whole_path_passive=True,passivity_scope='ANALYTIC_POINTWISE_COULOMB_WITH_SOURCE_BOUND_POSITIVE_GAP_AND_POLYNOMIAL_NONZERO_SLIP',
                interior_coulomb_direction='EXACT_POINTWISE_COULOMB_ON_DENSE_PATH_CALLBACK_ROUNDOFF_CONSISTENT',
                quadrature=nodes,quadrature_order=quadrature['order'],accepted=False))
        impulse = Jn+Jt; material_work = Wn+Wt; wall_work = Wnw+Wtw
        endpoint_forces = []
        for at in (left,right):
            force=np.zeros((count,3))
            for row in at['contacts']:
                force[row['particle_index']] += row['normal_force_N']*row['normal']+row['tangential_force_N']
            endpoint_forces.append(force)
        linear=hht_path.linear_work(endpoint_forces[0],endpoint_forces[1],v0,v1,a0,a1,dt,weight)
        return dict(classification='CURVED_EVENT_FREE_FIXED_NORMAL_CONTACT_QUADRATURE',accepted=False,
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
        if particle_hash(pa)!=before or _hash(p)!=parameter_before:
            raise CurvedContactError('curved contact callback/certificate changed accepted particles or parameters')

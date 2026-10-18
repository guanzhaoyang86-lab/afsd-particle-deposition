"""Event-response elastic/continuity correction of unaccepted stage defects.

B=d(vbar_actual)/dF_internal and C=d(v1_actual)/dF_internal are supplied
particle-local force responses.  No contact kick, material map, thermal map,
accepted history or raw acceptance threshold is replaced by this component.
"""
from copy import deepcopy
import hashlib
import numpy as np

from .stage_iteration import FIELDS
from .stage_preconditioner import _fields, _norm_inf, MAXIMUM_BACKWARD_ERROR
from .structural_mobility import FrozenElastic, build_frozen_elastic, _real, _readonly

CLASSIFICATION = 'EVENT_RESPONSE_SCHUR_CANDIDATE_PROPOSAL_NOT_PHYSICAL_MAP'


class StageEventPreconditioner:
    """Accepted-geometry frozen elastic block with separate B and C responses.

The original actuator mask determines prescribed rows.  Contact sticking is
encoded in the supplied responses, never by deleting a material degree of
freedom.  B may be singular, nonsymmetric, or zero: it is never inverted.
"""
    def __init__(self, accepted_pa, p, dt, time, B, C, provenance):
        from .axial_actuator import control, validate_interval
        from .mechanics import bulk_modulus
        from .updated_trial import particle_hash
        validate_interval(time, dt, p['release_time'])
        before = particle_hash(accepted_pa)
        elastic = build_frozen_elastic(accepted_pa, p['shear'],
                                       bulk_modulus(p['shear'], p['poisson']))
        program = control(accepted_pa.z0, time+.5*dt, 0., p['bar_length'],
                          p['actuator_depth'], p['feed'], p['release_time'])
        self._install(elastic, B, C, dt, program['free_axes'],
                      np.full(len(elastic.mass), p['rho0']), provenance)
        self.metadata.update(frozen_geometry='ACCEPTED_STATE',
            accepted_state_sha256=before, accepted_time_s=float(time),
            constraint_time_s=float(time+.5*dt), constraint_phase=program['phase'],
            pair_count=elastic.pair_count)
        if before != particle_hash(accepted_pa):
            raise RuntimeError('event preconditioner changed accepted particles')

    @classmethod
    def from_response_blocks(cls, elastic, B, C, dt, free_axes, rho_reference,
                             provenance=None):
        """Detached frozen block hook for independently manufactured baselines.

The caller can reuse an accepted-geometry FrozenElastic with new event B/C;
it must supply responses in exactly the FrozenElastic particle-ID order.
"""
        obj = cls.__new__(cls)
        obj._install(elastic, B, C, dt, free_axes, rho_reference, provenance)
        return obj

    def _install(self, elastic, B, C, dt, free_axes, rho_reference, provenance):
        if not isinstance(elastic, FrozenElastic):
            raise ValueError('explicit FrozenElastic required')
        ids = np.asarray(elastic.particle_ids)
        if (ids.ndim != 1 or not ids.size or ids.dtype.kind not in 'iu' or
                len(np.unique(ids)) != len(ids)):
            raise ValueError('unique integer particle IDs required')
        n, ndof = len(ids), 3*len(ids)
        step = _real(dt, 'dt', ())
        if step <= 0:
            raise ValueError('positive dt required')
        mask = np.asarray(free_axes)
        if (mask.shape != (n, 3) or mask.dtype.kind not in 'biuf' or
                not np.isfinite(mask).all() or not np.all((mask == 0) | (mask == 1))):
            raise ValueError('free_axes must contain exact boolean/zero-one values')
        b = _real(B, 'B', (n, 3, 3)); c = _real(C, 'C', (n, 3, 3))
        free = mask.astype(bool)
        if np.any(b[~free] != 0.) or np.any(c[~free] != 0.):
            raise ValueError('prescribed B/C output rows must be exactly zero')
        volume = _real(elastic.volume, 'volume', (n,))
        rho = _real(rho_reference, 'rho_reference', (n,))
        mass = _real(elastic.mass, 'mass', (n,))
        bulk = _real(elastic.bulk, 'bulk', ())
        if np.any(volume <= 0) or np.any(rho <= 0) or np.any(mass <= 0) or bulk <= 0:
            raise ValueError('positive volume, reference density, mass and bulk required')
        gradient = _real(elastic.gradient_map, 'gradient_map', (9*n, ndof))
        stiffness = _real(elastic.stiffness, 'stiffness', (ndof, ndof))
        e = np.trace(gradient.reshape(n, 3, 3, ndof), axis1=1, axis2=2)
        h = .5*e.T*(volume*float(bulk)/rho)[None, :]
        a = np.eye(ndof)+.5*float(step)*np.einsum(
            'nij,njk->nik', b, stiffness.reshape(n, 3, ndof)).reshape(ndof, ndof)
        if not np.isfinite(h).all() or not np.isfinite(a).all():
            raise ValueError('nonfinite event Schur operator')
        if provenance is not None and not isinstance(provenance, dict):
            raise ValueError('event response provenance must be a dictionary')
        origin = deepcopy(provenance or {})
        if 'particle_id' in origin and not np.array_equal(origin['particle_id'], ids):
            raise ValueError('event response particle-ID order mismatch')
        self.count, self.dt = n, float(step)
        self.particle_ids = _readonly(ids)
        self.B, self.C = _readonly(b), _readonly(c)
        self.free = _readonly(free.ravel())
        self.mass_per_dof = _readonly(np.repeat(mass, 3))
        self.divergence, self.pressure_force_derivative = _readonly(e), _readonly(h)
        self.stiffness, self.matrix = _readonly(stiffness), _readonly(a)
        self.rho_reference, self.volume = _readonly(rho), _readonly(volume)
        self.bulk = float(bulk)
        self.metadata = dict(classification=CLASSIFICATION, accepted=False,
            physical_map=False, restart_allowed=False, full_tangent=False,
            frozen_geometry='MANUFACTURED_OR_EXPLICIT_FROZEN_BLOCK_INPUT', dt_s=self.dt,
            acceptance='UNCHANGED_RAW_FIVE_FIELD_RESIDUAL_REQUIRED',
            mean_velocity_convention='INDEPENDENT_EVENT_AVERAGE_NOT_ENDPOINT_AVERAGE',
            density_linearization='RHO_APPROX_REFERENCE_DENSITY_FOR_PRECONDITIONER_ONLY',
            event_response_provenance=origin, event_B_equals_half_C_assumed=False,
            matrix_sha256=hashlib.sha256(a.tobytes(order='C')).hexdigest(),
            matrix_shape=list(a.shape),
            matrix_inf_norm=float(np.max(np.sum(np.abs(a), axis=1))),
            matrix_artifact_scope='READONLY_MATRIX_PROPERTY_SAVE_ONCE_PER_CONSTRUCTION',
            omitted_tangents=['geometry', 'prestress', 'objective_rotation', 'plasticity',
                             'thermal', 'contact_geometry', 'event_response_variation'])

    def apply(self, candidate, mapping):
        """Return candidate+correction, retaining the five true-map defects."""
        candidate, mapped = _fields(candidate, self.count), _fields(mapping, self.count)
        defect = {key: mapped[key]-candidate[key] for key in FIELDS}
        dbar = defect['mean_velocity_m_s'].ravel()
        pressure_force = self.pressure_force_derivative@defect['density_kg_m3']
        pressure_rhs = np.einsum('nij,nj->ni', self.B,
                                pressure_force.reshape(self.count, 3)).ravel()
        rhs = dbar+pressure_rhs
        free, constrained = self.free, ~self.free
        delta_bar = dbar.copy()
        aff = self.matrix[np.ix_(free, free)]
        shifted = rhs[free]-self.matrix[np.ix_(free, constrained)]@delta_bar[constrained]
        if np.any(free):
            try:
                delta_bar[free] = np.linalg.solve(aff, shifted)
            except np.linalg.LinAlgError as error:
                raise ValueError('singular event Schur matrix; no inverse, pinv or regularization') from error
        residual = self.matrix@delta_bar-rhs
        numerator = _norm_inf(residual)
        denominator = self.metadata['matrix_inf_norm']*_norm_inf(delta_bar)+_norm_inf(rhs)
        backward_error = numerator/denominator if denominator else 0.
        # Also audit the solved free block so large prescribed RHS values cannot
        # conceal an inaccurate free solve in a full-system normalized residual.
        free_residual = aff@delta_bar[free]-shifted
        free_matrix_norm = float(np.max(np.sum(np.abs(aff), axis=1), initial=0.))
        free_denominator = free_matrix_norm*_norm_inf(delta_bar[free])+_norm_inf(shifted)
        free_error = _norm_inf(free_residual)/free_denominator if free_denominator else 0.
        if (not np.isfinite(delta_bar).all() or not np.isfinite([backward_error, free_error]).all()
                or max(backward_error, free_error) > MAXIMUM_BACKWARD_ERROR):
            raise ValueError('event Schur linear backward error exceeds 1e-12')
        delta_force = pressure_force-.5*self.dt*(self.stiffness@delta_bar)
        force = delta_force.reshape(self.count, 3)
        mean_feedback = np.einsum('nij,nj->ni', self.B, force)
        endpoint_feedback = np.einsum('nij,nj->ni', self.C, force)
        density_feedback = -self.dt*self.rho_reference*(self.divergence@delta_bar)
        correction = {key: value.copy() for key, value in defect.items()}
        correction['mean_velocity_m_s'] = delta_bar.reshape(self.count, 3).copy()
        correction['position_m'] += self.dt*mean_feedback
        correction['velocity_m_s'] += endpoint_feedback
        correction['density_kg_m3'] += density_feedback
        proposal = {key: candidate[key]+correction[key] for key in FIELDS}
        if not all(np.isfinite(value).all() for value in proposal.values()):
            raise ValueError('nonfinite event Schur proposal; no clipping')
        info = deepcopy(self.metadata)
        info.update(particle_id=self.particle_ids.copy(),
            free_axes=free.reshape(self.count, 3).copy(),
            mass_per_dof_kg=self.mass_per_dof.copy(),
            density_reference_kg_m3=self.rho_reference.copy(), direct_defect=defect,
            B_s_kg=self.B.copy(), C_s_kg=self.C.copy(),
            pressure_density_force_rhs_N=pressure_force.copy(),
            pressure_density_rhs_m_s=pressure_rhs.copy(), rhs_m_s=rhs.copy(),
            prescribed_delta_mean_velocity_m_s=np.where(constrained, dbar, 0.),
            shifted_free_rhs_m_s=shifted.copy(), delta_mean_velocity_m_s=delta_bar.copy(),
            delta_internal_force_N=delta_force.copy(),
            density_backsubstitution_kg_m3=density_feedback.copy(),
            event_mean_feedback_m_s=mean_feedback.copy(),
            event_endpoint_feedback_m_s=endpoint_feedback.copy(),
            full_linear_residual_m_s=residual.copy(), free_linear_residual_m_s=free_residual.copy(),
            algebraic_constrained_row_defect_m_s=residual[constrained].copy(),
            prescribed_delta_defect_m_s=(delta_bar-dbar)[constrained].copy(),
            mean_backsubstitution_defect_m_s=(delta_bar-dbar-mean_feedback.ravel()),
            constrained_defect_scope='PRECONDITIONER_ALGEBRA_NOT_PHYSICAL_REACTION',
            backward_error_numerator_m_s=numerator, backward_error_denominator_m_s=denominator,
            backward_error_inf=backward_error, free_backward_error_inf=free_error,
            maximum_backward_error=MAXIMUM_BACKWARD_ERROR)
        return dict(proposal=proposal, correction=correction, diagnostics=info)

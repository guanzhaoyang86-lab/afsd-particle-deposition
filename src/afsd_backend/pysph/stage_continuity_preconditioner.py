"""Frozen elastic/continuity Schur correction of UNACCEPTED five-field defects.

The reference-density linearization is only a preconditioning model.  Actual
midpoint continuity, material, event contact and thermal maps stay unchanged.
"""
import numpy as np

from .stage_iteration import FIELDS
from .stage_preconditioner import (StagePreconditioner, _fields, _norm_inf,
                                   MAXIMUM_BACKWARD_ERROR)
from .structural_mobility import build_frozen_elastic

CLASSIFICATION = 'CONTINUITY_SCHUR_PRECONDITIONED_CANDIDATE_PROPOSAL_NOT_PHYSICAL_MAP'


class ContinuityStagePreconditioner(StagePreconditioner):
    """Freeze accepted geometry; eliminate the approximate density increment."""
    def __init__(self, accepted_pa, p, dt, time):
        from .axial_actuator import control, validate_interval
        from .mechanics import bulk_modulus
        from .updated_trial import particle_hash
        validate_interval(time, dt, p['release_time'])
        before = particle_hash(accepted_pa)
        bulk = bulk_modulus(p['shear'], p['poisson'])
        elastic = build_frozen_elastic(accepted_pa, p['shear'], bulk)
        program = control(accepted_pa.z0, time+.5*dt, 0., p['bar_length'],
                          p['actuator_depth'], p['feed'], p['release_time'])
        response = elastic.response_operator(dt, program['free_axes'])
        self._install(response, elastic.particle_ids)
        divergence = np.trace(elastic.gradient_map.reshape(self.count, 3, 3, 3*self.count),
                              axis1=1, axis2=2)
        self._continuity(divergence, elastic.volume, np.full(self.count, p['rho0']), bulk)
        self.metadata.update(frozen_geometry='ACCEPTED_STATE', accepted_state_sha256=before,
            accepted_time_s=float(time), constraint_time_s=float(time+.5*dt),
            constraint_phase=program['phase'], pair_count=elastic.pair_count)
        if before != particle_hash(accepted_pa):
            raise RuntimeError('continuity preconditioner changed accepted particles')

    @classmethod
    def from_response(cls, response, particle_ids, divergence, volume, rho_reference, bulk):
        """Explicit independent manufactured block, not a real SPH setup."""
        obj = cls.__new__(cls)
        obj._install(response, particle_ids)
        obj._continuity(divergence, volume, rho_reference, bulk)
        return obj

    def _continuity(self, divergence, volume, rho_reference, bulk):
        required = ((divergence, (self.count, 3*self.count), 'divergence'),
                    (volume, (self.count,), 'volume'),
                    (rho_reference, (self.count,), 'rho_reference'),
                    (bulk, (), 'bulk'))
        values = {}
        for value, shape, name in required:
            a = np.asarray(value)
            if (a.shape != shape or a.dtype.kind not in 'iuf' or not np.isfinite(a).all() or
                    (name != 'divergence' and np.any(a <= 0))):
                raise ValueError('invalid continuity '+name)
            values[name] = np.array(a, dtype=float, copy=True)
            values[name].flags.writeable = False
        self.divergence = values['divergence']
        self.volume = values['volume']
        self.rho_reference = values['rho_reference']
        self.bulk = float(values['bulk'])
        self.metadata.update(classification=CLASSIFICATION,
            density_linearization='RHO_APPROX_REFERENCE_DENSITY_FOR_PRECONDITIONER_ONLY',
            density_reference_kg_m3=self.rho_reference.copy())

    def apply(self, candidate, mapping):
        c, mapped = _fields(candidate, self.count), _fields(mapping, self.count)
        defect = {key: mapped[key]-c[key] for key in FIELDS}
        op, e = self.response, self.divergence
        dbar = defect['mean_velocity_m_s'].ravel()
        pressure_rhs = .25*op.dt*(e.T @ (self.volume*self.bulk/self.rho_reference *
                                       defect['density_kg_m3']))
        inertial_rhs = op.mass_per_dof*dbar
        rhs = inertial_rhs+pressure_rhs
        prescribed = np.zeros_like(dbar)
        prescribed[~op.free] = dbar[~op.free]
        solved = op.response(rhs, prescribed)
        delta_bar = solved['delta_velocity_m_s']
        free, constrained = op.free, ~op.free
        aff = op.matrix[np.ix_(free, free)]
        shifted = rhs[free]-op.matrix[np.ix_(free, constrained)]@prescribed[constrained]
        residual = aff@delta_bar[free]-shifted
        matrix_norm = float(np.max(np.sum(np.abs(aff), axis=1), initial=0.))
        numerator = _norm_inf(residual)
        denominator = matrix_norm*_norm_inf(delta_bar[free])+_norm_inf(shifted)
        error = numerator/denominator if denominator else 0.
        if not np.isfinite(error) or error > MAXIMUM_BACKWARD_ERROR:
            raise ValueError('continuity preconditioner linear backward error exceeds 1e-12')
        change = (delta_bar-dbar).reshape(self.count, 3)
        density_feedback = -op.dt*self.rho_reference*(e@delta_bar)
        correction = {key: value.copy() for key, value in defect.items()}
        correction['mean_velocity_m_s'] = delta_bar.reshape(self.count, 3).copy()
        correction['position_m'] += op.dt*change
        correction['velocity_m_s'] += 2.*change
        correction['density_kg_m3'] += density_feedback
        proposal = {key: c[key]+correction[key] for key in FIELDS}
        if not all(np.isfinite(value).all() for value in proposal.values()):
            raise ValueError('nonfinite continuity proposal; no clipping')
        info = dict(self.metadata,
            particle_id=self.particle_ids.copy(), free_axes=free.reshape(self.count, 3).copy(),
            mass_per_dof_kg=op.mass_per_dof.copy(), direct_defect=defect,
            inertial_rhs_Ns=inertial_rhs, pressure_density_rhs_Ns=pressure_rhs,
            rhs_Ns=rhs, prescribed_delta_mean_velocity_m_s=prescribed,
            shifted_free_rhs_Ns=shifted, delta_mean_velocity_m_s=delta_bar.copy(),
            density_backsubstitution_kg_m3=density_feedback,
            algebraic_reaction_impulse_Ns=solved['reaction_impulse_Ns'].copy(),
            reaction_scope='PRECONDITIONER_ALGEBRA_NOT_PHYSICAL_EXTERNAL_IMPULSE',
            free_linear_residual_Ns=residual,
            free_matrix_inf_norm_kg=matrix_norm, backward_error_numerator_Ns=numerator,
            backward_error_denominator_Ns=denominator, backward_error_inf=error,
            maximum_backward_error=MAXIMUM_BACKWARD_ERROR,
            mean_velocity_convention='INDEPENDENT_EVENT_AVERAGE_NOT_ENDPOINT_AVERAGE')
        # Detach array-valued metadata too; output journals cannot mutate P.
        info['density_reference_kg_m3'] = self.rho_reference.copy()
        return dict(proposal=proposal, correction=correction, diagnostics=info)

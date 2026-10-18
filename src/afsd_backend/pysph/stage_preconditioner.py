"""Physical-unit five-field proposal preconditioner; never a physical map.

The original event-average velocity remains an independent unknown.  This
module only changes unaccepted candidate corrections; no thermal/contact/
material residual or acceptance tolerance is replaced by the elastic model.
"""
import numpy as np

from .stage_iteration import FIELDS
from .structural_mobility import ConstrainedResponse, build_frozen_elastic

CLASSIFICATION = 'STRUCTURALLY_PRECONDITIONED_CANDIDATE_PROPOSAL_NOT_PHYSICAL_MAP'
MAXIMUM_BACKWARD_ERROR = 1e-12


def _norm_inf(value):
    return float(np.max(np.abs(value), initial=0.))


def _fields(value, count):
    if not isinstance(value, dict) or set(value) != set(FIELDS):
        raise ValueError('exact five-field candidate/map required')
    result = {}
    for key in FIELDS:
        raw = np.asarray(value[key])
        shape = (count,) if key in ('density_kg_m3', 'temperature_K') else (count, 3)
        if raw.shape != shape or raw.dtype.kind not in 'iuf' or not np.isfinite(raw).all():
            raise ValueError('invalid shape/type/nonfinite field '+key)
        result[key] = np.array(raw, dtype=float, copy=True)
    return result


class StagePreconditioner:
    """Accepted-geometry frozen elastic correction of a real five-field defect.

`time` is the accepted state's time; constraints match the actual kick midpoint.
Full process parameters are read, never changed.  No accepted ParticleArray is
stored or modified.  The outer root solver must re-evaluate the physical map
and judge its ORIGINAL residual after using this proposal.
"""
    def __init__(self, accepted_pa, p, dt, time):
        from .axial_actuator import control, validate_interval
        from .mechanics import bulk_modulus
        from .updated_trial import particle_hash
        validate_interval(time, dt, p['release_time'])
        before = particle_hash(accepted_pa)
        elastic = build_frozen_elastic(accepted_pa, p['shear'],
                                        bulk_modulus(p['shear'], p['poisson']))
        program = control(accepted_pa.z0, time+.5*dt, 0., p['bar_length'],
                          p['actuator_depth'], p['feed'], p['release_time'])
        response = elastic.response_operator(dt, program['free_axes'])
        self._install(response, elastic.particle_ids)
        self.metadata = dict(self.metadata, frozen_geometry='ACCEPTED_STATE',
            accepted_state_sha256=before, accepted_time_s=float(time),
            constraint_time_s=float(time+.5*dt), constraint_phase=program['phase'],
            pair_count=elastic.pair_count)
        if before != particle_hash(accepted_pa):
            raise RuntimeError('preconditioner construction changed accepted particles')

    @classmethod
    def from_response(cls, response, particle_ids):
        """Explicit manufactured/independent component hook, not AFSD setup."""
        obj = cls.__new__(cls)
        obj._install(response, particle_ids)
        return obj

    def _install(self, response, particle_ids):
        if not isinstance(response, ConstrainedResponse):
            raise ValueError('explicit ConstrainedResponse required')
        ids = np.asarray(particle_ids)
        if (ids.ndim != 1 or not ids.size or ids.dtype.kind not in 'iu' or
                len(np.unique(ids)) != len(ids) or len(response.mass_per_dof) != 3*len(ids)):
            raise ValueError('unique integer particle IDs and three axes required')
        # The structural component itself has detached read-only arrays.
        self.response = response
        self.particle_ids = np.array(ids, copy=True)
        self.particle_ids.flags.writeable = False
        self.count = len(ids)
        self.metadata = dict(classification=CLASSIFICATION, accepted=False,
            physical_map=False, restart_allowed=False, full_tangent=False,
            frozen_geometry='MANUFACTURED_COMPONENT_INPUT', dt_s=response.dt,
            acceptance='UNCHANGED_RAW_FIVE_FIELD_RESIDUAL_REQUIRED',
            structural_diagnostics=dict(response.diagnostics))

    def apply(self, candidate, mapping):
        """Return y+P(F(y)-y), P defect, and complete algebra diagnostics.

`mapping` must be the real five-field F(y), never an evaluation residual or
another preconditioned proposal. rho/T retain their direct defects. Values
outside the physical domain are NOT clipped: the outer evaluator rejects.
"""
        c, mapped = _fields(candidate, self.count), _fields(mapping, self.count)
        defect = {key: mapped[key]-c[key] for key in FIELDS}
        op = self.response
        dbar = defect['mean_velocity_m_s'].ravel()
        rhs = op.mass_per_dof*dbar
        prescribed = np.zeros_like(dbar)
        prescribed[~op.free] = dbar[~op.free]
        solved = op.response(rhs, prescribed)
        delta_bar = solved['delta_velocity_m_s']
        free, constrained = op.free, ~op.free
        aff = op.matrix[np.ix_(free, free)]
        shifted = rhs[free]-op.matrix[np.ix_(free, constrained)]@prescribed[constrained]
        residual = aff@delta_bar[free]-shifted
        matrix_norm = float(np.max(np.sum(np.abs(aff), axis=1), initial=0.))
        denominator = matrix_norm*_norm_inf(delta_bar[free])+_norm_inf(shifted)
        numerator = _norm_inf(residual)
        error = numerator/denominator if denominator else 0.
        if not np.isfinite(error) or error > MAXIMUM_BACKWARD_ERROR:
            raise ValueError('preconditioner free linear backward error exceeds 1e-12')
        change = (delta_bar-dbar).reshape(self.count, 3)
        correction = {key: value.copy() for key, value in defect.items()}
        correction['mean_velocity_m_s'] = delta_bar.reshape(self.count, 3).copy()
        correction['position_m'] += op.dt*change
        correction['velocity_m_s'] += 2.*change
        proposal = {key: c[key]+correction[key] for key in FIELDS}
        if not all(np.isfinite(value).all() for value in proposal.values()):
            raise ValueError('nonfinite preconditioner proposal; no clipping')
        diagnostics = dict(self.metadata,
            particle_id=self.particle_ids.copy(), free_axes=free.reshape(self.count, 3).copy(),
            mass_per_dof_kg=op.mass_per_dof.copy(), direct_defect=defect,
            rhs_Ns=rhs.copy(), prescribed_delta_mean_velocity_m_s=prescribed.copy(),
            shifted_free_rhs_Ns=shifted, delta_mean_velocity_m_s=delta_bar.copy(),
            algebraic_reaction_impulse_Ns=solved['reaction_impulse_Ns'].copy(),
            reaction_scope='PRECONDITIONER_ALGEBRA_NOT_PHYSICAL_EXTERNAL_IMPULSE',
            free_linear_residual_Ns=residual,
            free_matrix_inf_norm_kg=matrix_norm, backward_error_numerator_Ns=numerator,
            backward_error_denominator_Ns=denominator, backward_error_inf=error,
            maximum_backward_error=MAXIMUM_BACKWARD_ERROR,
            density_temperature_correction='DIRECT_TRUE_MAP_DEFECT',
            mean_velocity_convention='INDEPENDENT_EVENT_AVERAGE_NOT_ENDPOINT_AVERAGE')
        return dict(proposal=proposal, correction=correction, diagnostics=diagnostics)

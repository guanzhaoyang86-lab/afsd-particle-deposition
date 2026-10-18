"""Candidate-only material/continuity search with distinct contact B and C."""
import copy
import time
import numpy as np
from .material_stage_response import MaterialStageResponse, _array, _fingerprint
from .stage_preconditioner import _fields, MAXIMUM_BACKWARD_ERROR
from .stage_iteration import FIELDS

CLASSIFICATION = 'REFRESHED_MATERIAL_FROZEN_CONTACT_SECANT_SEARCH_APPROXIMATE'


class ContactMaterialResponse:
    def __init__(self, material, B, C):
        if not isinstance(material, MaterialStageResponse):
            raise ValueError('explicit material response required')
        self.material = material
        n = material.count
        self.B = _array(B, (n, 3, 3), 'contact mean secant')
        self.C = _array(C, (n, 3, 3), 'contact endpoint secant')
        free = material.free.reshape(n, 3)
        for a in (self.B, self.C):
            if np.any(a[~free] != 0.) or np.any(a.transpose(0, 2, 1)[~free] != 0.):
                raise ValueError('prescribed contact response rows/columns must be zero')
            a.flags.writeable = False
        H = material.Jv+(material.Jr*material.kappa[None, :])@material.E
        self.matrix = np.eye(3*n)-self._left(self.B, H)
        if not np.isfinite(self.matrix).all():
            raise ValueError('nonfinite contact material matrix')
        self.matrix.flags.writeable = False

    @staticmethod
    def _left(blocks, matrix):
        n = len(blocks)
        return np.einsum('nij,njk->nik', blocks, matrix.reshape(n, 3, -1)).reshape(3*n, -1)

    def apply(self, candidate, mapping):
        p = self.material; n = p.count
        c = _fields(candidate, n); mapped = _fields(mapping, n)
        d = {k: mapped[k]-c[k] for k in FIELDS}
        fixed_force = p.Jr@d['density_kg_m3']+p.JT@d['temperature_K']
        rhs = d['mean_velocity_m_s'].ravel()+self._left(self.B, fixed_force[:, None]).ravel()
        delta = d['mean_velocity_m_s'].ravel().copy(); f = p.free
        A = self.matrix; aff = A[np.ix_(f, f)]
        shifted = rhs[f]-A[np.ix_(f, ~f)]@delta[~f]
        if np.any(f):
            delta[f] = np.linalg.solve(aff, shifted)
        residual = aff@delta[f]-shifted
        norm = lambda x: float(np.max(abs(x), initial=0.))
        denominator = float(np.max(np.sum(abs(aff), axis=1), initial=0.))*norm(delta[f])+norm(shifted)
        backward = norm(residual)/denominator if denominator else 0.
        if not np.isfinite(delta).all() or not np.isfinite(backward) or backward > MAXIMUM_BACKWARD_ERROR:
            raise ValueError('contact material solve failed original backward error')
        correction = {k: v.copy() for k, v in d.items()}
        correction['mean_velocity_m_s'] = delta.reshape(n, 3).copy()
        correction['density_kg_m3'] += p.kappa*(p.E@delta)
        force = p.Jv@delta+p.Jr@correction['density_kg_m3']+p.JT@d['temperature_K']
        mean = self._left(self.B, force[:, None]).reshape(n, 3)
        endpoint = self._left(self.C, force[:, None]).reshape(n, 3)
        correction['position_m'] += p.dt*mean
        correction['velocity_m_s'] += endpoint
        proposal = {k: c[k]+correction[k] for k in FIELDS}
        if not all(np.isfinite(v).all() for v in proposal.values()):
            raise ValueError('nonfinite contact material proposal')
        return dict(proposal=proposal, correction=correction, diagnostics=dict(
            classification=CLASSIFICATION, accepted=False, restart_allowed=False,
            physical_map=False, full_tangent=False, backward_error_inf=backward,
            material_metadata=copy.deepcopy(p.metadata),
            mean_backsubstitution_defect_m_s=delta.reshape(n, 3)-d['mean_velocity_m_s']-mean,
            algebraic_row_defect_m_s=A@delta-rhs,
            delta_internal_force_N=force, contact_mean_feedback_m_s=mean,
            contact_endpoint_feedback_m_s=endpoint,
            reaction_scope='SEARCH_ALGEBRA_NOT_PHYSICAL_REACTION',
            omissions=['geometry', 'normal_force_and_wall_motion', 'BE_heat_feedback',
                       'contact_secant_refresh_and_classical_derivative_certificate']))


class BoundContactMaterialSearch:
    """Delegated one-call trial plus matching lazy proposal; no accepted writes."""
    def __init__(self, trial, parent, state, parameters, B, C):
        self.trial = trial
        self.parent = copy.deepcopy(parent); self.state = copy.deepcopy(state)
        self.parameters = copy.deepcopy(parameters)
        self.B = np.array(B, float, copy=True); self.C = np.array(C, float, copy=True)
        self.parent_hash = trial.accepted_hash()
        self.context = None
        self.calls = 0

    def candidate(self):
        return self.trial.candidate()

    def evaluate(self, candidate, dt):
        self.context = None
        if self.trial.accepted_hash() != self.parent_hash:
            raise ValueError('search parent changed')
        evaluation = self.trial.evaluate(candidate, dt)
        self.calls += 1
        if (evaluation['accepted_state_hash'] != self.parent_hash or
                evaluation['dt_s'] != dt or evaluation['time_before_s'] != self.state['t']):
            raise ValueError('search evaluation parent/time/dt mismatch')
        self.context = (_fingerprint(candidate), _fingerprint(evaluation['map']),
                        _fingerprint(evaluation), evaluation)
        return evaluation

    def apply(self, candidate, mapping):
        if self.context is None:
            raise ValueError('search has no completed matching trial')
        c_hash, m_hash, e_hash, evaluation = self.context
        # Consume once, including on an assembly error; no stale reuse.
        self.context = None
        if (self.trial.accepted_hash() != self.parent_hash or
                _fingerprint(candidate) != c_hash or _fingerprint(mapping) != m_hash or
                _fingerprint(evaluation) != e_hash):
            raise ValueError('stale candidate/map/evaluation/parent search context')
        tick = time.perf_counter()
        material = MaterialStageResponse.from_journal(self.parent, self.state,
            self.parameters, candidate, evaluation)
        block = ContactMaterialResponse(material, self.B, self.C)
        answer = block.apply(candidate, mapping)
        answer['diagnostics']['build_and_apply_wall_seconds'] = time.perf_counter()-tick
        return answer

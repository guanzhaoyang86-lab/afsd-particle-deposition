"""Bounded interior-equilibrium initial guess; diagnostic companion, not AFSD."""
import numpy as np

from probe_process_contact import reduced_correction
from probe_process_friction_cone import cone_rows, release_inventory
from probe_process_plastic import CandidateFailure


def solve_with_feasible_predictor(owner, t, released, maximum_iterations=14, observe=None):
    """Pure branch: interior predictor + original cone Newton share one budget.

    Interior correction is candidate initialization, not a physical substep.
    Hold every contact component only during initialization, never project N.
    No commit/advance, clipping, pseudoinverse or mutation of accepted history.
    """
    if type(maximum_iterations) is not int or not 0 <= maximum_iterations <= 14:
        raise ValueError('integer shared correction budget in [0,14] required')
    if not np.isfinite(t) or t <= owner.time:
        raise ValueError('target time must advance')
    released = release_inventory(released, len(owner.nodes)).copy()
    stick = ~released
    unknown = owner.free.copy()
    unknown[owner.contact_dofs] = False
    unknown[owner.tangent_dofs[stick].ravel()] = False
    interior = owner.free.copy()
    interior[owner.contact_dofs] = False
    interior[owner.tangent_dofs.ravel()] = False
    u = owner.displacement.copy()
    u[owner.nodes, 0] = owner.wall(t)[0]
    u[owner.nodes[stick], 1:] = (owner.displacement[owner.nodes[stick], 1:]
                                + (owner.wall(t)-owner.wall(owner.time))[1:])
    lim = owner.ccfg['numerical']
    halves = owner.qcfg['numerical']['line_search_halvings']
    trace = []
    used = 0
    predictor_corrections = 0

    def raw(displacement, stage):
        c = owner.trial(displacement, t)
        if not np.isfinite(c['internal_force']).all():
            raise CandidateFailure('nonfinite raw predictor force', c, trace)
        if observe is not None:
            observe(stage, c)
        return c

    c = raw(u, 'INITIAL_RAW')
    while np.any(c['internal_force'][owner.nodes, 0] < 0):
        force = c['internal_force'].ravel()
        norm = float(np.linalg.norm(force[interior])/owner.force_scale)
        row = dict(stage='INTERIOR_PREDICTOR', corrections_used=used,
                   interior_force_relative=norm,
                   normal_force_N=c['internal_force'][owner.nodes, 0].tolist())
        trace.append(row)
        if norm <= lim['force_relative']:
            raise CandidateFailure('negative normal at interior equilibrium; predictor blocked, not opening proof', c, trace)
        if used == maximum_iterations:
            raise CandidateFailure('shared correction budget exhausted in predictor', c, trace)
        try:
            step, cond = reduced_correction(owner.tangent(c), force, interior,
                                           lim['tangent_condition_limit'])
        except (ValueError, np.linalg.LinAlgError) as exc:
            raise CandidateFailure(str(exc), c, trace) from exc
        row['tangent_condition'] = cond
        for half in range(halves+1):
            trial_u = u + 2.**(-half)*step.reshape(u.shape)
            probe = raw(trial_u, 'PREDICTOR_LINE_SEARCH')
            after = float(np.linalg.norm(probe['internal_force'].ravel()[interior])/owner.force_scale)
            if after < norm:
                u, c = trial_u, probe
                row['line_search_factor'] = 2.**(-half)
                used += 1
                predictor_corrections += 1
                break
        else:
            raise CandidateFailure('bounded interior predictor line search failed', c, trace)

    # All raw normals are nonnegative BEFORE calling the unchanged contact law.
    while True:
        c = owner.state(c, released)
        residual = c['residual'].ravel().copy()
        residual[owner.tangent_dofs[released].ravel()] = c['law_residual_N'][released].ravel()
        norm = float(np.linalg.norm(residual[unknown])/owner.force_scale)
        row = dict(stage='CONE_ROOT', iteration=used, root_relative=norm,
                   normal_force_N=c['normal_force'].tolist())
        trace.append(row)
        if norm <= lim['force_relative']:
            c.update(iterations=used, predictor_corrections=predictor_corrections,
                     cone_corrections=used-predictor_corrections, iteration_trace=trace)
            return c
        if used == maximum_iterations:
            raise CandidateFailure('shared correction budget exhausted in cone root', c, trace)
        K, r = cone_rows(owner.tangent(c), c['residual'].ravel(),
                         owner.contact_dofs, owner.tangent_dofs,
                         c['response'], owner.mu, owner.kappa, released)
        try:
            step, cond = reduced_correction(K, r, unknown, lim['tangent_condition_limit'])
        except (ValueError, np.linalg.LinAlgError) as exc:
            raise CandidateFailure(str(exc), c, trace) from exc
        row['tangent_condition'] = cond
        for half in range(halves+1):
            trial_u = u + 2.**(-half)*step.reshape(u.shape)
            raw_probe = raw(trial_u, 'CONE_LINE_SEARCH')
            if np.any(raw_probe['internal_force'][owner.nodes, 0] < 0):
                continue  # Same original domain guard; retain raw trace, never clip.
            probe = owner.state(raw_probe, released)
            rr = probe['residual'].ravel().copy()
            rr[owner.tangent_dofs[released].ravel()] = probe['law_residual_N'][released].ravel()
            if np.linalg.norm(rr[unknown])/owner.force_scale < norm:
                u, c = trial_u, probe
                row['line_search_factor'] = 2.**(-half)
                used += 1
                break
        else:
            raise CandidateFailure('cone branch line search failed', c, trace)

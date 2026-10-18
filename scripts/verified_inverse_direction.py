"""Checked use of an existing Broyden inverse; no physics or matrix mutation.

The explicit search matrix remains authoritative. An inverse application is
used only when its actual A @ direction - rhs satisfies the original absolute
backward-error gate. At most one improving correction is attempted before one
ordinary dense solve. The caller retains graph refresh and nonlinear budgets.
"""
import time
import numpy as np

CLASSIFICATION = 'VERIFIED_BROYDEN_INVERSE_ONE_CORRECTION_DENSE_FALLBACK_V1'
CONDITION_LIMIT = 1e12
BACKWARD_FACTOR = 1e-9


class LinearDirectionBudgetError(RuntimeError):
    """Terminal elapsed-time failure; do not spend a tangent refresh on it."""


def verified_inverse_direction(matrix, inverse, rhs, *, deadline=None):
    """Return (direction, diagnostics) with the frozen solver's linear gate.

    Inputs are never written. matrix/inverse must represent the same current
    graph, free-axis ordering and Broyden version; the graph-aware caller owns
    that binding. This function neither repairs a condition failure nor changes
    the physical residual tolerance. The condition estimate is the existing
    ||A||inf ||H||inf policy, not a claim that H is an exact inverse.
    """
    begin = time.perf_counter()
    if deadline is not None and not np.isfinite(deadline):
        raise ValueError('finite absolute deadline required')

    def budget():
        if deadline is not None and time.perf_counter() > deadline:
            raise LinearDirectionBudgetError('verified linear-direction deadline exhausted')

    budget()
    matrix, inverse, rhs = np.asarray(matrix, float), np.asarray(inverse, float), np.asarray(rhs, float)
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not len(matrix) or
            inverse.shape != matrix.shape or rhs.shape != (len(matrix),) or
            not all(np.isfinite(v).all() for v in (matrix, inverse, rhs))):
        raise ValueError('finite matching nonempty square A/H and vector rhs required')
    condition = float(np.linalg.norm(matrix, np.inf)*np.linalg.norm(inverse, np.inf))
    if not np.isfinite(condition) or condition > CONDITION_LIMIT:
        raise ValueError('near-singular secant search matrix; no regularization')
    limit = BACKWARD_FACTOR*max(1., float(np.max(np.abs(rhs))))
    diagnostics = dict(classification=CLASSIFICATION, matrix_size=len(rhs),
        condition_infinity_estimate=condition, linear_backward_limit_N=limit,
        correction_attempted=False, correction_beneficial=False, dense_fallback_used=False,
        inverse_initial_backward_error_N=None, correction_backward_error_N=None,
        inverse_application_wall_s=0., correction_wall_s=0., dense_fallback_wall_s=0.)

    def residual(direction):
        if not np.isfinite(direction).all(): return None, np.inf
        with np.errstate(over='ignore', invalid='ignore'):
            value = matrix@direction-rhs
        if not np.isfinite(value).all(): return None, np.inf
        return value, float(np.max(np.abs(value)))

    start = time.perf_counter()
    with np.errstate(over='ignore', invalid='ignore'):
        direction = inverse@rhs
    defect, backward = residual(direction)
    diagnostics['inverse_application_wall_s'] = time.perf_counter()-start
    diagnostics['inverse_initial_backward_error_N'] = backward if np.isfinite(backward) else None
    budget()
    method = 'VERIFIED_INVERSE_APPLICATION'
    if backward > limit and defect is not None:
        diagnostics['correction_attempted'] = True
        start = time.perf_counter()
        with np.errstate(over='ignore', invalid='ignore'):
            corrected = direction-inverse@defect
        _, corrected_backward = residual(corrected)
        diagnostics['correction_wall_s'] = time.perf_counter()-start
        diagnostics['correction_backward_error_N'] = corrected_backward if np.isfinite(corrected_backward) else None
        budget()
        if corrected_backward < backward:
            direction, backward = corrected, corrected_backward
            diagnostics['correction_beneficial'] = True
            method = 'VERIFIED_INVERSE_ONE_BENEFICIAL_CORRECTION'
    if backward > limit:
        budget()
        diagnostics['dense_fallback_used'] = True
        start = time.perf_counter()
        direction = np.linalg.solve(matrix, rhs)
        _, backward = residual(direction)
        diagnostics['dense_fallback_wall_s'] = time.perf_counter()-start
        budget()
        method = 'VERIFIED_ORIGINAL_DENSE_FALLBACK'
    if not np.isfinite(direction).all() or not np.isfinite(backward) or backward > limit:
        raise ValueError('secant linear solve backward error')
    diagnostics.update(method=method, linear_backward_error_N=backward,
                       direction_wall_seconds=time.perf_counter()-begin)
    return direction, diagnostics

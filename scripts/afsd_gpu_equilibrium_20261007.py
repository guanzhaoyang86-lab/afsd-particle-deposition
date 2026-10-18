"""Static GPU dense-linear-algebra adaptation of the frozen CPU searches.

Only inverse/solve dispatch changes. Full physical maps, tangent assembly,
Broyden rank-one updates, condition estimates, backward checks, Armijo searches,
graph-change handling and all original budgets stay on their CPU code paths.
The backend returns NumPy float64 and does not silently fall back to the CPU.
No module globals are rebound. This module alone is not physical validation.
"""
from pathlib import Path
import time
import numpy as np
import afsd_gpu_linalg_20261007 as gpu_linalg
from secant_equilibrium import (MAX_REFRESHES, MAX_EVALUATIONS, SearchBudgetError,
    _condition, _broyden, summary, write)
from visibility_equilibrium import CONTROLLER, journal
from run_discrete_contact_process import journal as newton_journal


def _full_matrix(trial, result, free):
    tangent = np.asarray(trial.tangent(result), float)
    if tangent.shape != (len(free), len(free)) or not np.isfinite(tangent).all():
        raise ValueError('finite full tangent of physical displacement required')
    matrix = tangent[np.ix_(free, free)].copy()
    # The inverse is used ONLY for condition/update diagnostics; actual
    # directions are independently solved with np.linalg.solve each time.
    inverse = gpu_linalg.inv(matrix)
    condition = _condition(matrix, inverse)
    return matrix, inverse, condition


def solve(trial, dt, guess, directory, deadline, maximum=24):
    """Return full trial, metrics using <=24 accepted updates/100 evaluations.

    At most two full Jacobian refreshes follow the initial build.
    Any accepted actual support graph change mandates a rebuild before another
    search direction; no secant difference is taken across that change.  Refresh is
    requested after two slow accepted reductions (>0.8), a heavily damped
    step (alpha<0.25), a line-search failure, or an unsafe secant update.
    Neither residual thresholds nor physical acceptance rules are changed.
    """
    if not isinstance(maximum, (int, np.integer)) or not 0 <= maximum <= 24:
        raise ValueError('maximum must be an integer in [0,24]')
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(deadline):
        raise ValueError('positive dt and finite absolute deadline required')
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    free_shape = np.asarray(trial.state['free'], bool)
    free = free_shape.ravel()
    if free_shape.ndim != 2 or free_shape.shape[1] != 3 or not np.any(free):
        raise ValueError('nonempty physical free displacement mask required')
    initial = np.array(guess, float, copy=True)
    if initial.shape != free_shape.shape or not np.isfinite(initial).all():
        raise ValueError('finite guess must match physical displacement shape')
    fixed_values = initial[~free_shape].copy()
    rows = []; result = None; count = 0; builds = 0; refreshes = 0
    rank_updates = 0; slow_count = 0; matrix = inverse = None
    refresh_reason = None; unsafe_reason = None

    def budget():
        if count >= MAX_EVALUATIONS or time.perf_counter() > deadline:
            raise SearchBudgetError('registered secant evaluation/time budget exhausted')

    def evaluate(displacement):
        nonlocal count
        budget(); count += 1  # failed domain candidates also consume budget
        value = trial.evaluate(displacement, dt)
        if not isinstance(value.get('graph_fingerprint'), str):
            raise ValueError('every complete trial must expose its actual graph fingerprint')
        force = np.asarray(value['force'], float)
        candidate_u = np.asarray(value['displacement'], float)
        if (force.shape != free_shape.shape or candidate_u.shape != free_shape.shape or
                not np.isfinite(force).all() or not np.isfinite(candidate_u).all()):
            raise ValueError('complete candidate force/displacement must be finite')
        if not np.array_equal(candidate_u[~free_shape], fixed_values):
            raise ValueError('candidate changed prescribed physical axes')
        if not np.isfinite(value['force_limit_N']) or value['force_limit_N'] <= 0:
            raise ValueError('invalid physical force limit')
        return value

    try:
        result = evaluate(initial)
        for iteration in range(maximum+1):
            actual = float(np.max(np.abs(result['force'][free_shape])))
            row = dict(iteration=iteration, evaluations=count, controller=CONTROLLER,
                actual_free_force_residual_N=actual, graph_fingerprint=result['graph_fingerprint'], **summary(result))
            rows.append(row); write(directory/'iterations.json', rows)
            # Never return merely because an approximate/search residual or
            # stale scaled_residual claims convergence.
            if actual <= result['force_limit_N'] and result['scaled_residual'] <= 1:
                journal(directory, trial, result, 'last_UNCOMMITTED_candidate')
                return result, dict(iterations=iteration, evaluations=count,
                    jacobian_builds=builds, jacobian_refreshes=refreshes,
                    broyden_updates=rank_updates, controller=CONTROLLER)
            if iteration == maximum:
                raise SearchBudgetError('registered secant iteration budget exhausted')
            budget()
            if unsafe_reason is not None and refreshes >= MAX_REFRESHES:
                raise ValueError('unsafe secant matrix and refresh budget exhausted: '+unsafe_reason)
            if matrix is None or (refresh_reason is not None and refreshes < MAX_REFRESHES):
                reason = 'INITIAL_FULL_TANGENT' if matrix is None else refresh_reason
                if matrix is not None: refreshes += 1
                builds += 1
                matrix, inverse, condition = _full_matrix(trial, result, free)
                row.setdefault('jacobian_builds', []).append(dict(reason=reason,
                    condition_infinity=condition, build_number=builds))
                refresh_reason = unsafe_reason = None; slow_count = 0
            found = False
            while not found:
                budget()
                rhs = -np.asarray(result['force'])[free_shape]
                try:
                    direction = gpu_linalg.solve(matrix, rhs)
                    backward = float(np.max(np.abs(matrix@direction-rhs)))
                    if not np.isfinite(direction).all() or backward > 1e-9*max(1., float(np.max(np.abs(rhs)))):
                        raise ValueError('secant linear solve backward error')
                except (ValueError, np.linalg.LinAlgError):
                    if refreshes >= MAX_REFRESHES: raise
                    refreshes += 1
                    builds += 1
                    matrix, inverse, condition = _full_matrix(trial, result, free)
                    row.setdefault('jacobian_builds', []).append(dict(reason='LINEAR_FAILURE',
                        condition_infinity=condition, build_number=builds))
                    continue
                delta = np.zeros(free_shape.shape); delta[free_shape] = direction
                merit = float(np.linalg.norm(result['force'][free_shape]))
                row.setdefault('linear', []).append(dict(linear_backward_error_N=backward,
                    nnz=int(np.count_nonzero(matrix)), fixed_axis_increment_zero=True,
                    jacobian_build_number=builds,
                    condition_policy='INITIAL_INVERSE_NORM_THEN_SHERMAN_MORRISON_GUARDED_ESTIMATE'))
                for backtrack in range(8):
                    budget(); alpha = .5**backtrack
                    try:
                        candidate = evaluate(result['displacement']+alpha*delta)
                        next_merit = float(np.linalg.norm(candidate['force'][free_shape]))
                        row.setdefault('line_search', []).append(dict(alpha=alpha, merit_N=next_merit,
                            jacobian_build_number=builds, graph_fingerprint=candidate['graph_fingerprint'],
                            graph_changed=candidate['graph_fingerprint']!=result['graph_fingerprint']))
                        if next_merit < (1.-1e-4*alpha)*merit:
                            found = True; break
                    except SearchBudgetError:
                        raise
                    except (ValueError, ArithmeticError) as error:
                        row.setdefault('line_search', []).append(dict(alpha=alpha, rejected=str(error),
                            jacobian_build_number=builds))
                    write(directory/'iterations.json', rows)
                if found: break
                if refreshes >= MAX_REFRESHES:
                    raise ValueError('complete residual Armijo search failed; refresh budget exhausted')
                budget(); refreshes += 1
                builds += 1
                matrix, inverse, condition = _full_matrix(trial, result, free)
                row.setdefault('jacobian_builds', []).append(dict(reason='ARMIJO_FAILURE',
                    condition_infinity=condition, build_number=builds))
            previous = result; result = candidate
            step = (result['displacement']-previous['displacement'])[free_shape]
            change = (result['force']-previous['force'])[free_shape]
            scale = max(float(np.linalg.norm(result['displacement'][free_shape])),
                        float(np.linalg.norm(previous['displacement'][free_shape])))
            if result['graph_fingerprint'] != previous['graph_fingerprint']:
                unsafe_reason = refresh_reason = 'ACTUAL_SUPPORT_GRAPH_CHANGED'
                row['secant'] = dict(status='DISCARDED_ACROSS_GRAPH_CHANGE',
                    previous_graph=previous['graph_fingerprint'], current_graph=result['graph_fingerprint'])
            else:
                try:
                    matrix, inverse, secant = _broyden(matrix, inverse, step, change, scale)
                    row['secant'] = secant
                    if secant['status'] == 'GOOD_BROYDEN_RANK_ONE': rank_updates += 1
                except (ValueError, ArithmeticError) as error:
                    unsafe_reason = refresh_reason = str(error)
                    row['secant'] = dict(status='REFRESH_REQUIRED_UNSAFE_UPDATE', reason=str(error))
            ratio = next_merit/merit
            slow_count = slow_count+1 if ratio > .8 else 0
            if refresh_reason is None and (slow_count >= 2 or alpha < .25):
                refresh_reason = 'POOR_PROGRESS' if slow_count >= 2 else 'HEAVY_DAMPING'
            row.update(accepted_alpha=alpha, merit_ratio=ratio, evaluations_after_search=count)
            write(directory/'iterations.json', rows)
    except Exception as error:
        if rows:
            rows[-1]['terminal_error'] = repr(error)
            rows[-1]['evaluations_at_failure'] = count
        write(directory/'iterations.json', rows)
        write(directory/'secant_failure.json', dict(controller=CONTROLLER, error=repr(error),
            evaluations=count, jacobian_builds=builds, jacobian_refreshes=refreshes,
            accepted=False, restart_allowed=False))
        if result is not None:
            journal(directory, trial, result, 'last_UNCOMMITTED_candidate')
        raise


def newton_direction(trial, result):
    """No mass shift, pseudoinverse, or artificial constraints on free modes."""
    free = trial.state['free'].ravel()
    J = trial.tangent(result)
    A = J[np.ix_(free, free)]
    rhs = -result['force'].ravel()[free]
    delta = np.zeros(free.shape)
    delta[free] = gpu_linalg.solve(A, rhs)
    backward_error = float(np.max(np.abs(A@delta[free]-rhs)))
    if backward_error > 1e-9*max(1., float(np.max(np.abs(rhs)))):
        raise ValueError('quasi-static linear solve backward error')
    return delta.reshape(trial.x0.shape), dict(
        linear_backward_error_N=backward_error, nnz=int(np.count_nonzero(A)),
        coercivity='NOT_CERTIFIED_BY_LU; INITIAL_ELASTIC_RANK_AUDIT_REQUIRED', inertia=False)


def newton_solve(trial, dt, guess, directory, deadline, maximum=24):
    r = trial.evaluate(guess, dt); count = 1; rows = []
    for iteration in range(maximum+1):
        row = dict(iteration=iteration, evaluations=count, **summary(r)); rows.append(row)
        write(directory/'iterations.json', rows)
        if r['scaled_residual'] <= 1:
            newton_journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            return r, dict(iterations=iteration, evaluations=count)
        if iteration == maximum or time.perf_counter() > deadline:
            newton_journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            raise ValueError('registered iteration/time budget exhausted')
        try:
            delta, linear = newton_direction(trial, r)
        except Exception:
            newton_journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            raise
        merit = np.linalg.norm(r['force'][trial.state['free']]); found = False
        for backtrack in range(8):
            alpha = .5**backtrack
            try:
                next_r = trial.evaluate(r['displacement']+alpha*delta, dt); count += 1
                next_merit = float(np.linalg.norm(next_r['force'][trial.state['free']]))
                row.setdefault('line_search', []).append(dict(alpha=alpha, merit_N=next_merit))
                if next_merit < (1-1e-4*alpha)*merit:
                    r = next_r; found = True; break
            except (ValueError, ArithmeticError) as error:
                count += 1
                row.setdefault('line_search', []).append(dict(alpha=alpha, rejected=str(error)))
            if count >= 100 or time.perf_counter() > deadline: break
        row['linear'] = linear; write(directory/'iterations.json', rows)
        if not found or count >= 100:
            newton_journal(directory, trial, r, 'last_UNCOMMITTED_candidate')
            raise ValueError('complete residual line search/evaluation budget failed')

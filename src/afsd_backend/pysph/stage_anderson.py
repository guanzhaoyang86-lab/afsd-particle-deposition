"""Bounded type-II Anderson iteration of the unchanged five-field stage map.

Only unaccepted candidates are mixed.  This module commits no particle state,
changes no physical time step, and adds no damping force or thermal source.
"""
import time
import numpy as np

from .stage_iteration import (FIELDS, TOLERANCES, StageFailure,
                              initial_candidate, residual_norms)

BETA = .75
DEPTH = 6
MAXIMUM_QR_CONDITION = 1e10
METHOD = 'TYPE_II_ANDERSON_CANDIDATES_ONLY_NOT_PHYSICAL_DAMPING'


class _Coordinates:
    """Dimensionless increments relative to a fixed, copied reference state."""
    def __init__(self, reference):
        self.reference = {k: np.asarray(reference[k], float).copy() for k in FIELDS}
        density = self.reference['density_kg_m3']
        if np.any(density <= 0) or not np.isfinite(density).all():
            raise ValueError('nonpositive/nonfinite initial density')
        self.scales = {}
        for key in FIELDS:
            value = self.reference[key]
            if not value.size or not np.isfinite(value).all():
                raise ValueError('invalid initial candidate '+key)
            self.scales[key] = (TOLERANCES[key]*density if key == 'density_kg_m3'
                                else np.full(value.shape, TOLERANCES[key]))

    def pack(self, candidate):
        return np.concatenate([((candidate[k]-self.reference[k])/self.scales[k]).ravel()
                               for k in FIELDS])

    def difference(self, mapping, candidate):
        return np.concatenate([((np.asarray(mapping[k])-candidate[k])/self.scales[k]).ravel()
                               for k in FIELDS])

    def unpack(self, coordinates):
        result = {}; offset = 0
        for key in FIELDS:
            old = self.reference[key]; count = old.size
            result[key] = old + self.scales[key]*coordinates[offset:offset+count].reshape(old.shape)
            offset += count
        if offset != coordinates.size:
            raise ValueError('packed candidate dimension mismatch')
        return result


def _accelerate(x, f, dx, df):
    """Thin QR least squares; explicitly discard ill-conditioned oldest secants.

    No pseudoinverse, singular-value truncation, or regularized physical map is
    used.  Returned metadata describes every discarded history column.
    """
    dx = np.asarray(dx, float); df = np.asarray(df, float)
    if dx.shape != df.shape or dx.ndim != 2 or dx.shape[0] != x.size:
        raise ValueError('invalid Anderson secant shape')
    dropped = []
    while df.shape[1]:
        count = df.shape[1]
        q, r = np.linalg.qr(df, mode='reduced')
        condition = (float(np.linalg.cond(r)) if r.shape == (count, count)
                     else float('inf'))
        if np.isfinite(condition) and condition <= MAXIMUM_QR_CONDITION:
            gamma = np.linalg.solve(r, q.T @ f)
            next_x = x + BETA*f - (dx+BETA*df) @ gamma
            if not np.isfinite(next_x).all():
                raise ValueError('nonfinite Anderson candidate')
            return next_x, dict(method='ANDERSON', secant_depth=count,
                                conditioning=condition, discarded_history=dropped,
                                fallback=('DROPPED_OLDEST_ILL_CONDITIONED_SECANTS' if dropped else None))
        # A singular QR diagnostic is not a physical field.  Preserve its
        # meaning explicitly while keeping journals strict JSON (no Infinity).
        dropped.append(dict(secant_depth=count,
                            conditioning=condition if np.isfinite(condition) else None,
                            conditioning_status=('FINITE_ILL_CONDITIONED' if np.isfinite(condition)
                                                 else 'NONFINITE_SINGULAR')))
        dx, df = dx[:, 1:], df[:, 1:]
    return x+BETA*f, dict(method='PICARD', secant_depth=0, conditioning=None,
                         discarded_history=dropped,
                         fallback=('NO_WELL_CONDITIONED_SECANT' if dropped else 'INSUFFICIENT_HISTORY'))


def solve_stage(trial, dt, maximum_evaluations=32, wall_seconds=540., on_evaluation=None):
    """Return a converged *unaccepted* root, or a traceable StageFailure.

    Callback receives each evaluated candidate, including an
    accelerated candidate rejected by the residual-growth safeguard.  An
    evaluation/domain exception is retained as a final trace row and in
    StageFailure, with callback evaluation=None, never replaced with a prior
    candidate's evaluation.
    """
    if (type(maximum_evaluations) is not int or maximum_evaluations < 1 or
            not np.isfinite(wall_seconds) or wall_seconds <= 0):
        raise ValueError('invalid bounded Anderson settings')
    c = initial_candidate(trial, dt)
    coordinates = _Coordinates(c)
    density = c['density_kg_m3'].copy()
    start = time.perf_counter(); trace = []; history = []; best = None; last = None
    origin = dict(method='INITIAL', secant_depth=0, conditioning=None,
                  discarded_history=[], fallback=None)
    for index in range(maximum_evaluations):
        if time.perf_counter()-start >= wall_seconds:
            raise StageFailure('global Anderson wall budget', c, trace, None)
        tick = time.perf_counter(); last = None
        try:
            last = trial.evaluate(c, dt)
            norms, scaled = residual_norms(c, last['map'], density)
            x = coordinates.pack(c); f = coordinates.difference(last['map'], c)
        except Exception as exc:
            row = dict(evaluation=index+1, **origin, evaluation_status='FAILED_TRIAL',
                       candidate_disposition='REJECTED_TRIAL_ERROR', error=str(exc),
                       evaluation_wall_seconds=time.perf_counter()-tick)
            trace.append(row)
            if on_evaluation is not None:
                on_evaluation(row, c, None)
            raise StageFailure('global Anderson trial rejected: '+str(exc),
                               c, trace, last, exc) from exc
        reject_growth = (origin['method'] == 'ANDERSON' and best is not None and
                         scaled > 2.*best['scaled'])
        row = dict(evaluation=index+1, **origin, evaluation_status='EVALUATED', residual_norms=norms,
                   maximum_scaled_residual=scaled,
                   candidate_disposition=('REJECTED_ACCELERATED_CANDIDATE' if reject_growth
                                          else 'ROOT_CANDIDATE_NOT_PHYSICAL_STATE'),
                   evaluation_wall_seconds=time.perf_counter()-tick)
        trace.append(row)
        if on_evaluation is not None:
            on_evaluation(row, c, last)
        if time.perf_counter()-start >= wall_seconds:
            raise StageFailure('global Anderson wall budget after evaluation', c, trace, last)
        if not reject_growth and scaled <= 1.:
            return dict(status='ROOT_CONVERGED_NOT_ACCEPTED', candidate=c, evaluation=last,
                        trace=trace, wall_seconds=time.perf_counter()-start, method=METHOD)
        if index+1 == maximum_evaluations:
            break
        if reject_growth:
            # Discard all old secants, but anchor the newly restarted history
            # at the already evaluated best point.  The next fresh secant must
            # span best -> fallback; dropping this anchor can erase a crossed
            # constitutive kink and repeat the identical failed extrapolation.
            history = [(best['x'].copy(), best['f'].copy())]
            next_x = best['x']+BETA*best['f']
            origin = dict(method='PICARD', secant_depth=0, conditioning=None,
                          discarded_history=[], fallback='GROWTH_GT_2_BEST_REJECTED_RESET_HISTORY')
        else:
            if best is None or scaled < best['scaled']:
                best = dict(x=x.copy(), f=f.copy(), scaled=scaled)
            history.append((x.copy(), f.copy()))
            history = history[-(DEPTH+1):]
            if len(history) > 1:
                dx = np.column_stack([b[0]-a[0] for a, b in zip(history[:-1], history[1:])])
                df = np.column_stack([b[1]-a[1] for a, b in zip(history[:-1], history[1:])])
            else:
                dx = df = np.empty((x.size, 0))
            try:
                next_x, origin = _accelerate(x, f, dx, df)
            except Exception as exc:
                raise StageFailure('global Anderson algebra rejected: '+str(exc),
                                   c, trace, last, exc) from exc
            if origin['discarded_history']:
                history = history[-(origin['secant_depth']+1):]
        c = coordinates.unpack(next_x)
    raise StageFailure('global Anderson evaluation limit', c, trace, last)

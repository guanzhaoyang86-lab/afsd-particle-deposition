"""Raw-residual-controlled Anderson with an algebraic structural proposal.

This module never commits state, changes dt, or replaces physical impulses.
The preconditioner is not a convergence certificate.
"""
import time
import numpy as np
from .stage_iteration import StageFailure, initial_candidate, residual_norms
from .stage_anderson import _Coordinates, _accelerate, BETA, DEPTH

METHOD = 'FROZEN_ELASTIC_PRECONDITIONED_ANDERSON_RAW_FIVE_FIELD_ACCEPTANCE'


def solve_stage(trial, preconditioner, dt, maximum_evaluations=32,
                wall_seconds=540., on_evaluation=None):
    if (type(maximum_evaluations) is not int or maximum_evaluations < 1 or
            not np.isfinite(wall_seconds) or wall_seconds <= 0):
        raise ValueError('invalid bounded structural root settings')
    c = initial_candidate(trial, dt)
    coordinates = _Coordinates(c); density = c['density_kg_m3'].copy()
    start = time.perf_counter(); trace = []; history = []; best = None; last = None
    origin = dict(method='INITIAL', secant_depth=0, conditioning=None,
                  discarded_history=[], fallback=None)
    for index in range(maximum_evaluations):
        if time.perf_counter()-start >= wall_seconds:
            raise StageFailure('structural root wall budget', c, trace)
        tick = time.perf_counter(); last = None
        try:
            last = trial.evaluate(c, dt)
            norms, scaled = residual_norms(c, last['map'], density)
            x = coordinates.pack(c)
            reject_growth = (origin['method'] == 'ANDERSON' and best is not None
                             and scaled > 2.*best['scaled'])
            # No algebraic proposal can override or relax raw convergence.
            if scaled > 1. and not reject_growth:
                proposal = preconditioner.apply(c, last['map'])
                f = coordinates.difference(proposal['proposal'], c)
                if not np.isfinite(f).all():
                    raise ValueError('nonfinite preconditioned proposal')
                last['numerical_proposal_only'] = proposal
            else:
                f = None
        except Exception as exc:
            row = dict(evaluation=index+1, **origin, evaluation_status='FAILED_TRIAL_OR_PROPOSAL',
                       candidate_disposition='REJECTED_TRIAL_ERROR', error=str(exc),
                       evaluation_wall_seconds=time.perf_counter()-tick)
            trace.append(row)
            if on_evaluation is not None: on_evaluation(row, c, last)
            raise StageFailure('structural root rejected: '+str(exc), c, trace, last, exc) from exc
        row = dict(evaluation=index+1, **origin, evaluation_status='EVALUATED',
                   residual_norms=norms, maximum_scaled_residual=scaled,
                   maximum_scaled_proposal=(float(np.max(np.abs(f))) if f is not None else None),
                   candidate_disposition=('REJECTED_ACCELERATED_CANDIDATE' if reject_growth
                                          else 'ROOT_CANDIDATE_NOT_PHYSICAL_STATE'),
                   evaluation_wall_seconds=time.perf_counter()-tick)
        trace.append(row)
        if on_evaluation is not None: on_evaluation(row, c, last)
        if time.perf_counter()-start >= wall_seconds:
            raise StageFailure('structural root wall budget after evaluation', c, trace, last)
        if not reject_growth and scaled <= 1.:
            return dict(status='ROOT_CONVERGED_NOT_ACCEPTED', candidate=c, evaluation=last,
                        trace=trace, wall_seconds=time.perf_counter()-start, method=METHOD)
        if index+1 == maximum_evaluations: break
        if reject_growth:
            history = [(best['x'].copy(), best['f'].copy())]
            next_x = best['x']+BETA*best['f']
            origin = dict(method='PRECONDITIONED_PICARD', secant_depth=0, conditioning=None,
                          discarded_history=[], fallback='RAW_GROWTH_GT_2_BEST_REJECTED_RESET_HISTORY')
        else:
            if best is None or scaled < best['scaled']:
                best = dict(x=x.copy(), f=f.copy(), scaled=scaled)
            history.append((x.copy(), f.copy())); history = history[-(DEPTH+1):]
            if len(history)>1:
                dx = np.column_stack([b[0]-a[0] for a,b in zip(history[:-1],history[1:])])
                df = np.column_stack([b[1]-a[1] for a,b in zip(history[:-1],history[1:])])
            else:
                dx = df = np.empty((x.size, 0))
            try:
                next_x, origin = _accelerate(x, f, dx, df)
            except Exception as exc:
                raise StageFailure('structural root algebra rejected: '+str(exc), c, trace, last, exc) from exc
            if origin['discarded_history']:
                history = history[-(origin['secant_depth']+1):]
        c = coordinates.unpack(next_x)
    raise StageFailure('structural root evaluation limit', c, trace, last)

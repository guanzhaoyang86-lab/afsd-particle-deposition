"""Bounded global candidate fixed-point solve; no physical commit or clipping."""
import time
import numpy as np

FIELDS = ('position_m', 'velocity_m_s', 'density_kg_m3', 'temperature_K', 'mean_velocity_m_s')
TOLERANCES = dict(position_m=1e-12, velocity_m_s=1e-9, density_kg_m3=1e-10,
                  temperature_K=5e-12, mean_velocity_m_s=1e-9)


class StageFailure(RuntimeError):
    def __init__(self, message, candidate, trace, evaluation=None, original=None):
        super().__init__(message)
        self.candidate = {k: np.asarray(v).copy() for k, v in candidate.items()}
        self.trace = list(trace); self.evaluation = evaluation
        self.original_context = getattr(original, 'contact_context', getattr(original, 'feature_context', None))


def initial_candidate(trial, dt):
    if not np.isfinite(dt) or dt <= 0: raise ValueError('positive physical dt required')
    c = {k: np.asarray(v, float).copy() for k, v in trial.candidate().items()}
    if set(c) != set(FIELDS): raise ValueError('candidate interface mismatch')
    c['position_m'] += dt*c['mean_velocity_m_s']
    return c


def residual_norms(candidate, mapping, reference_density):
    if set(candidate) != set(FIELDS) or set(mapping) != set(FIELDS):
        raise ValueError('candidate/map field mismatch')
    norms = {}
    for key in FIELDS:
        a, b = np.asarray(candidate[key]), np.asarray(mapping[key])
        if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError('invalid candidate/map '+key)
        residual = a-b
        if key == 'density_kg_m3': residual = residual/reference_density
        norms[key] = float(np.max(np.abs(residual)))
    return norms, max(norms[k]/TOLERANCES[k] for k in FIELDS)


def solve_stage(trial, dt, maximum_evaluations=32, relaxation=.75, wall_seconds=540., on_evaluation=None):
    if (type(maximum_evaluations) is not int or maximum_evaluations < 1 or
            not np.isfinite(relaxation) or not 0 < relaxation <= 1 or
            not np.isfinite(wall_seconds) or wall_seconds <= 0):
        raise ValueError('invalid bounded iteration settings')
    c = initial_candidate(trial, dt); density = c['density_kg_m3'].copy()
    if np.any(density <= 0): raise ValueError('nonpositive initial density')
    start = time.perf_counter(); trace = []; last = None
    for index in range(maximum_evaluations):
        if time.perf_counter()-start >= wall_seconds:
            raise StageFailure('global candidate wall budget', c, trace, last)
        tick = time.perf_counter()
        try:
            last = trial.evaluate(c, dt)
            norms, scaled = residual_norms(c, last['map'], density)
        except Exception as exc:
            raise StageFailure('global trial rejected: '+str(exc), c, trace, last, exc) from exc
        row = dict(evaluation=index+1, residual_norms=norms, maximum_scaled_residual=scaled,
                   evaluation_wall_seconds=time.perf_counter()-tick)
        trace.append(row)
        if on_evaluation is not None: on_evaluation(row, c, last)
        if time.perf_counter()-start >= wall_seconds:
            raise StageFailure('global candidate wall budget after evaluation', c, trace, last)
        if scaled <= 1.:
            return dict(status='ROOT_CONVERGED_NOT_ACCEPTED', candidate=c, evaluation=last,
                        trace=trace, wall_seconds=time.perf_counter()-start,
                        method='DAMPED_PICARD_CANDIDATES_ONLY_NOT_PHYSICAL_DAMPING')
        if index+1 < maximum_evaluations:
            c = {k: c[k]+relaxation*(np.asarray(last['map'][k])-c[k]) for k in FIELDS}
    raise StageFailure('global candidate evaluation limit', c, trace, last)

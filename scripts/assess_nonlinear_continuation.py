"""Read-only saved-array audit; never evaluate a new physical candidate.

The directional rate diagnostic keeps the saved midpoint geometry/density
fixed. It is NOT a global Jacobian, accepted advance, or rerun of material.py.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assess():
    start = time.perf_counter()
    run = ROOT/'outputs/updated_trial_20261003/global_attempt'
    paths = [run/'latest_uncommitted_trial.json.gz', run/'iterations.json',
             run/'resolved_config.json', Path(__file__)]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    with gzip.open(paths[0], 'rt', encoding='utf-8') as handle:
        saved = json.load(handle)
    e, c = saved['evaluation'], saved['candidate']
    p = json.loads(paths[2].read_text(encoding='utf-8'))['parameters']
    trace = json.loads(paths[1].read_text(encoding='utf-8'))
    mid = e['midpoint']
    ids = np.asarray(e['particle_id'])
    i, j = np.asarray(mid['i']), np.asarray(mid['j'])
    B, g = np.asarray(mid['B']), np.asarray(mid['g'])
    volume = np.asarray(mid['mass_kg'])/np.asarray(mid['density_kg_m3'])
    mean = np.asarray(c['mean_velocity_m_s'])
    delta = np.asarray(e['map']['mean_velocity_m_s'])-mean

    def gradient(value):
        raw = np.zeros_like(B)
        pair = (value[j]-value[i])[:, :, None]*g[:, None, :]
        np.add.at(raw, i, volume[j, None, None]*pair)
        np.add.at(raw, j, volume[i, None, None]*pair)
        return raw@B

    actual_L = gradient(mean)
    reference_L = np.asarray(mid['velocity_gradient_per_s'])
    reconstruction_error = float(np.max(np.abs(actual_L-reference_L)))
    if reconstruction_error > 1e-9:
        raise ValueError('saved-array gradient reconstruction mismatch')
    directional_L = gradient(delta)
    points = []
    for pid in (266, 268, 269, 270):
        k = int(np.flatnonzero(ids == pid)[0])
        temperature = mid['temperature_K'][k]
        ts = (temperature-p['room'])/(p['melt']-p['room'])
        logarithmic_coefficient = p['b1']*ts**p['b2']
        # Independent direct expression, not the production expm1/logaddexp.
        theta = 1-(1+np.exp(-p['c1']*ts))**(-1/p['c2'])
        ep = e['material']['ep'][k]
        base = (p['a1']+p['a2']*np.arctan(p['a3']*ep))*theta
        rows = []
        for alpha in (0., .1, .25, .5, .75, 1.):
            L = actual_L[k]+alpha*directional_L[k]
            D = .5*(L+L.T)
            dev = D-np.trace(D)*np.eye(3)/3
            rate = float(np.sqrt(2*np.sum(dev**2)/3))
            # This max IS the documented original low-rate regularization;
            # no heat, stress, temperature or root defect is clipped here.
            law_rate = max(rate, p['rate_floor'])
            Y = base*(1+logarithmic_coefficient*np.log(law_rate/p['reference_rate']))
            rows.append(dict(candidate_direction_fraction=alpha,
                             frozen_geometry_rate_per_s=rate,
                             fixed_ep_yield_Pa=float(Y)))
        contact = next(t for t in e['boundary']['traces'] if t['particle_id'] == pid)
        points.append(dict(particle_id=pid, ep=ep,
            contact_features=[r['feature_id'] for r in contact['contact_rows']],
            logarithmic_rate_coefficient=float(logarithmic_coefficient),
            yield_at_floor_Pa=float(base),
            right_derivative_at_floor_Pa_s=float(base*logarithmic_coefficient/p['rate_floor']),
            hypothetical_candidate_directions=rows))
    residual = np.asarray(c['velocity_m_s'])-np.asarray(e['map']['velocity_m_s'])
    k, axis = np.unravel_index(np.argmax(np.abs(residual)), residual.shape)
    result = dict(status='PASS_READ_ONLY_DIAGNOSTIC_NOT_PHYSICAL_ADVANCE',
        source_hashes=hashes, evaluated_new_physical_candidates=0,
        accepted=False, restart_allowed=False,
        approximation='FROZEN_SAVED_GEOMETRY_DENSITY_AND_FINAL_EP_DIRECTION_ONLY',
        saved_gradient_reconstruction_max_per_s=reconstruction_error,
        largest_velocity_defect=dict(particle_id=int(ids[k]), axis=int(axis),
                                    value_m_s=float(residual[k, axis])),
        direct_neighbor_268_269=bool(np.any(((ids[i] == 268) & (ids[j] == 269)) |
                                          ((ids[i] == 269) & (ids[j] == 268)))),
        prior_iteration_count=len(trace),
        prior_minimum_velocity_defect_m_s=min(t['residual_norms']['velocity_m_s'] for t in trace),
        particles=points, wall_seconds=time.perf_counter()-start)
    if any(sha(ROOT/key) != value for key, value in hashes.items()):
        raise ValueError('read-only evidence changed during audit')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = assess()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))

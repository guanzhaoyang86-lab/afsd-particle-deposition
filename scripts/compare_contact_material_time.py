"""Read-only, same-parent/same-endpoint comparison; no solver calls or pass threshold."""
from pathlib import Path
import gzip
import hashlib
import json
import time
import numpy as np
from pysph.solver.utils import load
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/contact_material_search_20261004'


def read(p):
    with (gzip.open(p, 'rt', encoding='utf-8') if p.suffix == '.gz' else p.open(encoding='utf-8')) as f:
        return json.load(f)


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    start = time.perf_counter(); watched = {}
    for folder in ('root', 'time_check'):
        path = BASE/folder; m = read(path/'manifest.json')
        for rel, h in m['artifact_sha256'].items():
            p = (path/rel).resolve()
            if not p.is_relative_to(path) or sha(p) != h: raise ValueError('artifact drift '+rel)
            watched[p] = h
        for rel, h in m['source_sha256'].items():
            if sha(ROOT/rel) != h: raise ValueError('source drift '+rel)
            watched[ROOT/rel] = h
        watched[path/'manifest.json'] = sha(path/'manifest.json')
    coarse = read(BASE/'root/manifest.json'); fine = read(BASE/'time_check/manifest.json')
    if fine['status'] != 'COMPLETE_SAME_ENDPOINT_TWO_STAGES':
        raise ValueError('time branch incomplete: do not compare different endpoints')
    if coarse['parent_reference'] != fine['parent_reference']:
        raise ValueError('different physical parent')
    if abs(coarse['accepted_time_s']-fine['accepted_time_s']) > 4*np.finfo(float).eps*coarse['accepted_time_s']:
        raise ValueError('different endpoint')
    a = load(str(BASE/'root/accepted.npz'))['arrays']['solid']
    b = load(str(BASE/'time_check/stage_002/accepted.npz'))['arrays']['solid']
    parent = load(str(ROOT/coarse['parent_reference']['path']))['arrays']['solid']
    order_a, order_b, order_parent = [np.argsort(q.particle_id) for q in (a, b, parent)]
    for k in ('particle_id', 'm', 'x0', 'y0', 'z0', 'tracer_type', 'wholebar_tag', 'as_rs_tag'):
        if not (np.array_equal(getattr(a, k)[order_a], getattr(b, k)[order_b]) and
                np.array_equal(getattr(a, k)[order_a], getattr(parent, k)[order_parent])):
            raise ValueError('identity/mass/label mismatch '+k)
    p = read(BASE/'root/resolved_config.json')['parameters']
    def field(pa, key, order):
        if key == 'position_m': return np.column_stack((pa.x, pa.y, pa.z))[order]
        if key == 'velocity_m_s': return np.column_stack((pa.u, pa.v, pa.w))[order]
        if key == 'deviatoric_stress_Pa': return pa.S.reshape(-1, 9)[order]
        return getattr(pa, key)[order]
    measures = {}
    for key, characteristic in [('position_m', p['dx']), ('velocity_m_s', p['feed']),
        ('rho', p['rho0']), ('temperature', 1.), ('deviatoric_stress_Pa', p['shear']),
        ('ep', 1.), ('plastic_work', None), ('return_dissipation', None)]:
        x = field(a, key, order_a); y = field(b, key, order_b)
        diff = y-x; old = field(parent, key, order_parent)
        norm = lambda q: np.linalg.norm(q, axis=1) if q.ndim == 2 else np.abs(q)
        maximum = float(np.max(norm(diff)))
        measures[key] = dict(maximum_absolute_difference=maximum,
            rms_difference=float(np.sqrt(np.mean(norm(diff)**2))),
            coarse_maximum_change_from_common_parent=float(np.max(norm(x-old))),
            fine_maximum_change_from_common_parent=float(np.max(norm(y-old))),
            characteristic_scale=characteristic,
            difference_over_characteristic=(maximum/characteristic if characteristic else None))
    full = read(BASE/'root/acceptance.json')
    parts = [read(BASE/f'time_check/stage_{i:03d}/acceptance.json') for i in (1, 2)]
    ledgers = {}
    for area, keys in [('contact_ledger', ['friction_raw_J', 'wall_work_J', 'actuator_work_J', 'internal_work_J']),
                       ('thermal_ledger', ['wall_heat_out_J'])]:
        for key in keys:
            x = full[area][key]; y = sum(v[area][key] for v in parts)
            ledgers[key] = dict(coarse_J=x, fine_J=y, fine_minus_coarse_J=y-x)
    for key in ('total_energy_diagnostic_residual_J', 'normal_temporal_residual_J'):
        x = full[key]; y = sum(v[key] for v in parts)
        ledgers[key] = dict(coarse_J=x, fine_J=y, fine_minus_coarse_J=y-x)
    ea = read(BASE/'root/accepted_trial.json.gz')['evaluation']
    eb = read(BASE/'time_check/stage_002/accepted_trial.json.gz')['evaluation']
    ledgers['final_normal_potential_J'] = dict(coarse_J=ea['normal_after']['normal_potential_J'],
        fine_J=eb['normal_after']['normal_potential_J'],
        fine_minus_coarse_J=eb['normal_after']['normal_potential_J']-ea['normal_after']['normal_potential_J'])
    history = read(BASE/'time_check/accepted_history.json')
    transport = read(BASE/'time_check/stage_002/transport.json')
    result = dict(status='COMPLETE_SAME_ENDPOINT_SENSITIVITY_NOT_ACCURACY_VALIDATION',
        parent_time_s=coarse['source_time_s'], endpoint_s=coarse['accepted_time_s'],
        intervals_s=dict(coarse=[5e-7], fine=[2.5e-7, 2.5e-7]),
        fine_stage_evaluations=[x['evaluations'] for x in history],
        fine_stage_scaled_residuals=[x['maximum_scaled_residual'] for x in history],
        endpoint_fields=measures, integrated_ledgers=ledgers,
        fine_shoulder_crossing_mass_kg=transport['cumulative']['downward_event_mass_kg'],
        fine_branch_wall_seconds=fine['wall_seconds'],
        interpretation='Two time grids, no fitted tolerance, no order or continuum/experiment validation')
    for path, h in watched.items():
        if sha(path) != h: raise ValueError('comparison changed source/evidence')
    result.update(verified_unchanged_files=len(watched), wall_seconds=time.perf_counter()-start)
    output = BASE/'time_comparison.json'
    if output.exists(): raise ValueError('refuse overwrite')
    output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__': main()

"""E1 manufactured material-point history, NOT an AFSD/AA6061 simulation."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from afsd_backend.pysph.material import update

CONFIG = ROOT / 'configs/benchmarks/engineering_j2_cycle.json'


class PerfectPlasticBenchmark:
    """Test-only isothermal law. Never selected by a production configuration."""
    def __init__(self, yield_stress, temperature):
        if not np.isfinite([yield_stress, temperature]).all() or yield_stress <= 0:
            raise ValueError('invalid benchmark constants')
        self.value, self.temperature = yield_stress, temperature

    def yield_stress(self, ep, rate, temperature):
        ep, rate, temperature = np.broadcast_arrays(ep, rate, temperature)
        if not all(np.isfinite(x).all() for x in (ep, rate, temperature)):
            raise ValueError('non-finite benchmark input')
        if np.any(ep < 0) or np.any(rate < 0):
            raise ValueError('negative strain/rate')
        if np.any(temperature != self.temperature):
            raise ValueError('benchmark temperature must remain prescribed')
        return np.full(ep.shape, self.value)


def load_config(path=CONFIG):
    cfg = json.loads(Path(path).read_text(encoding='utf-8'))
    validate_config(cfg)
    return cfg


def validate_config(cfg):
    """This oracle is only for the preregistered manufactured problem."""
    if cfg['production_authorized'] is not False:
        raise ValueError('material-point benchmark cannot authorize production')
    expected = {'shear': (30e9, 'Pa'), 'yield_stress': (90e6, 'Pa'),
                'reference_volume': (1e-6, 'm^3'), 'temperature': (293.15, 'K'),
                'peak_log_strain': (.004, '1'), 'segment_duration': (1., 's'),
                'directions': ([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1/3, 2/3, 2/3]], '1'),
                'aligned_steps': ([16, 32, 64], '1'), 'unaligned_steps': ([17, 33, 65], '1')}
    for key, (value, unit) in expected.items():
        entry = cfg['parameters'][key]
        if entry['value'] != value or entry['unit'] != unit:
            raise ValueError('changed E1 parameter requires new oracle/preregistration: '+key)


def analytic_history(strain, loading, shear, yield_stress, peak):
    """Independent scalar piecewise solution, never calls material.update."""
    ey = yield_stress / (3 * shear)
    r, ep = np.empty(len(strain)), np.empty(len(strain))
    for k, (e, first) in enumerate(zip(strain, loading)):
        if first:
            if e <= ey:
                r[k], ep[k] = 3 * shear * e, 0.
            else:
                r[k], ep[k] = yield_stress, e - ey
        elif e >= peak - 2 * ey:
            r[k], ep[k] = yield_stress + 3 * shear * (e - peak), peak - ey
        else:
            r[k], ep[k] = -yield_stress, 2 * peak - 3 * ey - e
    return r, ep


def run_history(cfg, steps):
    validate_config(cfg)
    if steps not in (16, 32, 64, 17, 33, 65):
        raise ValueError('unregistered E1 time grid')
    p = {k: v['value'] for k, v in cfg['parameters'].items()}
    shear, ys, volume = p['shear'], p['yield_stress'], p['reference_volume']
    peak, duration = p['peak_log_strain'], p['segment_duration']
    n = np.asarray(p['directions'], dtype=float)
    if steps <= 0 or not np.allclose(np.sum(n*n, axis=1), 1., rtol=0., atol=1e-14):
        raise ValueError('invalid steps or orientation')
    nn = n[:, :, None] * n[:, None, :]
    a = 1.5 * nn - .5 * np.eye(3)
    dt, count = duration / steps, len(n)
    strain = np.concatenate((np.linspace(0, peak, steps+1), np.linspace(peak, 0, steps+1)[1:]))
    loading = np.arange(len(strain)) <= steps
    stress, ep = np.zeros((count, 3, 3)), np.zeros(count)
    temperature = np.full(count, p['temperature'])
    law = PerfectPlasticBenchmark(ys, p['temperature'])
    arrays = {'time_s': np.arange(2*steps+1)*dt, 'prescribed_log_strain': strain,
              'stress_Pa': [stress.copy()], 'plastic_strain': [ep.copy()],
              'cumulative_stress_work_J': [np.zeros(count)],
              'cumulative_plastic_work_J': [np.zeros(count)],
              'cumulative_algorithmic_return_J': [np.zeros(count)],
              'plastic_work_increment_J': [], 'algorithmic_return_increment_J': [],
              'equivalent_strain_rate_per_s': []}
    for k in range(1, len(strain)):
        gradient = ((strain[k]-strain[k-1])/dt)*a
        old = stress.copy()
        stress, ep, plastic, algorithmic, rate = update(stress, ep, gradient, temperature, dt, shear, law)
        work = volume * np.sum(.5*(old+stress)*gradient, axis=(1, 2))*dt
        arrays['stress_Pa'].append(stress.copy())
        arrays['plastic_strain'].append(ep.copy())
        for name, increment in [('stress_work', work), ('plastic_work', plastic*volume),
                                ('algorithmic_return', algorithmic*volume)]:
            key = 'cumulative_'+name+'_J'
            arrays[key].append(arrays[key][-1]+increment)
        arrays['plastic_work_increment_J'].append(plastic*volume)
        arrays['algorithmic_return_increment_J'].append(algorithmic*volume)
        arrays['equivalent_strain_rate_per_s'].append(rate)
    arrays = {k: np.asarray(v) for k, v in arrays.items()}
    e = strain[:, None, None, None]
    arrays['prescribed_F'] = np.exp(-e/2)*np.eye(3)+(np.exp(e)-np.exp(-e/2))*nn[None]
    arrays['elastic_shear_diagnostic_J'] = volume*np.sum(arrays['stress_Pa']**2, axis=(2, 3))/(4*shear)
    r_ref, ep_ref = analytic_history(strain, loading, shear, ys, peak)
    s_ref = (2/3)*r_ref[:, None, None, None]*a[None]
    arrays['reference_stress_Pa'] = s_ref
    arrays['reference_plastic_strain'] = np.broadcast_to(ep_ref[:, None], arrays['plastic_strain'].shape).copy()
    return arrays


def assess(cfg, steps, a):
    """The fixed E1 oracle is NOT computed from the tested updater."""
    validate_config(cfg)
    if steps not in (16, 32, 64, 17, 33, 65):
        raise ValueError('unregistered E1 time grid')
    t = cfg['tolerances']
    p = {k: v['value'] for k, v in cfg['parameters'].items()}
    expected_work = .495 if steps % 4 == 0 else .495-.315/steps**2
    dep = np.diff(a['plastic_strain'], axis=0)
    errors = {
        'stress_relative': float(np.max(np.abs(a['stress_Pa']-a['reference_stress_Pa']))/p['yield_stress']),
        'ep_absolute': float(np.max(np.abs(a['plastic_strain']-a['reference_plastic_strain']))),
        'work_reference_J': float(np.max(np.abs(a['cumulative_stress_work_J'][-1]-expected_work))),
        'plastic_work_reference_J': float(np.max(np.abs(a['cumulative_plastic_work_J'][-1]-.45))),
        'elastic_reference_J': float(np.max(np.abs(a['elastic_shear_diagnostic_J'][-1]-.045))),
        'plastic_increment_identity_J': float(np.max(np.abs(a['plastic_work_increment_J']-p['yield_stress']*dep*p['reference_volume']))),
        'return_increment_identity_J': float(np.max(np.abs(a['algorithmic_return_increment_J']-1.5*p['shear']*dep**2*p['reference_volume']))),
        'trace_relative': float(np.max(np.abs(np.trace(a['stress_Pa'], axis1=2, axis2=3)))/p['yield_stress']),
    }
    checks = {
        'finite': all(np.isfinite(x).all() for x in a.values()),
        'stress': errors['stress_relative'] <= t['stress_relative'],
        'ep': errors['ep_absolute'] <= t['ep_absolute'],
        'trace': errors['trace_relative'] <= t['stress_relative'],
        'nonnegative_raw_plastic_work': bool(np.all(a['plastic_work_increment_J'] >= 0)),
        'nonnegative_raw_algorithmic_return': bool(np.all(a['algorithmic_return_increment_J'] >= 0)),
        'common_end_time': bool(a['time_s'][-1] == 2.),
    }
    checks.update({k: v <= t['energy_absolute_J'] for k, v in errors.items() if k.endswith('_J')})
    return {'steps_per_segment': steps, 'status': 'PASS' if all(checks.values()) else 'FAIL',
            'checks': checks, 'errors': errors, 'expected_discrete_work_J': expected_work,
            'final_stress_work_J': a['cumulative_stress_work_J'][-1].tolist(),
            'final_plastic_work_J': a['cumulative_plastic_work_J'][-1].tolist(),
            'final_elastic_diagnostic_J': a['elastic_shear_diagnostic_J'][-1].tolist(),
            'algorithmic_return_J_not_heat': a['cumulative_algorithmic_return_J'][-1].tolist(),
            'continuous_work_error_J': float(np.max(np.abs(a['cumulative_stress_work_J'][-1]-.495)))}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    cfg = load_config()
    runs, hashes = [], {}
    for family in ('aligned_steps', 'unaligned_steps'):
        for steps in cfg['parameters'][family]['value']:
            arrays = run_history(cfg, steps)
            result = assess(cfg, steps, arrays)
            result['family'] = family
            artifact = args.output / ('history_N%d.npz' % steps)
            np.savez_compressed(artifact, **arrays)
            hashes[artifact.name] = sha(artifact)
            runs.append(result)
    unaligned = [r for r in runs if r['family'] == 'unaligned_steps']
    ratios = [b['continuous_work_error_J']/a['continuous_work_error_J'] for a, b in zip(unaligned, unaligned[1:])]
    refinement = all(r < cfg['tolerances']['continuous_work_error_ratio'] for r in ratios)
    refinement = refinement and unaligned[-1]['continuous_work_error_J']/.495 < cfg['tolerances']['unaligned_finest_relative']
    status = 'PASS' if refinement and all(r['status'] == 'PASS' for r in runs) else 'FAIL'
    sources = ['src/afsd_backend/pysph/material.py', 'scripts/probe_engineering_material.py',
               'configs/benchmarks/engineering_j2_cycle.json']
    report = {'classification': cfg['classification'], 'material_variant': cfg['material_variant'],
              'production_authorized': False, 'status': status, 'backend_commit': cfg['backend_commit'],
              'not_solved': cfg['not_solved'], 'wall_time_s': time.perf_counter()-started,
              'source_sha256': {p: sha(ROOT/p) for p in sources}, 'artifact_sha256': hashes,
              'unaligned_continuous_error_ratios': ratios, 'refinement_status': 'PASS' if refinement else 'FAIL',
              'runs': runs,
              'field_notes': {'prescribed_F': 'input kinematics, not solved particle deformation',
                              'elastic_shear_diagnostic_J': 'limited isochoric diagnostic, not certified general free energy',
                              'algorithmic_return': 'not added to physical heat; thermal NOT_RUN',
                              'reference_stress_Pa': 'independent piecewise scalar solution'}}
    (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if status == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())

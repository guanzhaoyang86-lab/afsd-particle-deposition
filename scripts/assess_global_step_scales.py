"""Read archived arrays and derive scales; NEVER evaluate a physical trial."""
from pathlib import Path
import gzip
import hashlib
import json
import math
import sys
import time
import numpy as np
from pysph.solver.utils import load

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'outputs/nonlinear_continuation_20261003/continuation'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.perf_counter()
    output = Path(sys.argv[1]).resolve()
    if not output.is_relative_to(ROOT) or output.exists():
        raise ValueError('fresh project-local output required')
    manifest = json.loads((CASE/'manifest.json').read_text(encoding='utf-8'))
    files = ['resolved_config.json', 'stage_004/accepted.npz',
             'stage_004/accepted_trial.json.gz', 'stage_004/acceptance.json']
    watched = {CASE/name: manifest['artifact_sha256'][name] for name in files}
    response_path = ROOT/'outputs/process_completion_20261003/response_probe/report.json'
    watched[response_path] = sha(response_path)
    for path, digest in watched.items():
        if sha(path) != digest:
            raise ValueError('archived input drift '+str(path))
    data = load(str(CASE/'stage_004/accepted.npz'))
    pa, state = data['arrays']['solid'], data['solver_data']
    if state['accepted'] is not True or state['restart_allowed'] is not True:
        raise ValueError('accepted checkpoint required')
    p = state['parameters']
    if p != json.loads((CASE/'resolved_config.json').read_text(encoding='utf-8'))['parameters']:
        raise ValueError('configuration mismatch')
    with gzip.open(CASE/'stage_004/accepted_trial.json.gz', 'rt', encoding='utf-8') as stream:
        archive = json.load(stream)
    e = archive['evaluation']
    x = np.column_stack((pa.x, pa.y, pa.z))
    x0 = np.column_stack((pa.x0, pa.y0, pa.z0))
    v = np.column_stack((pa.u, pa.v, pa.w))
    mass = pa.m.copy()
    force = np.asarray(e['internal_force_N'])
    bulk = 2*p['shear']*(1+p['poisson'])/(3*(1-2*p['poisson']))
    sound = math.sqrt((bulk+4*p['shear']/3)/p['rho0'])
    contact_omega = np.sqrt(p['normal_stiffness']/mass)
    normal = e['normal_after']
    traces = e['boundary']['traces']
    maximum_intervals = max(len(trace['result']['intervals']) for trace in traces)
    force_predictor = force/mass[:, None]
    response = json.loads(response_path.read_text(encoding='utf-8'))['responses'][0]
    # Uniform mass and an unconstrained rigid translation make lambda_min(A)=m.
    # The measured condition thus yields the max constrained elastic omega.
    if not np.all(mass == mass[0]):
        raise ValueError('condition-to-frequency shortcut requires uniform physical mass')
    elastic_omega = math.sqrt(4*(response['diagnostics']['free_matrix_condition']-1))/response['dt_s']
    result = dict(status='READ_ONLY_SCALES_NOT_A_TIME_STEP', new_coupled_steps=0,
        accepted_progress=False, production_authorized=False,
        accepted_time_s=state['t'], particle_count=len(x),
        current=dict(max_speed_m_s=float(np.linalg.norm(v, axis=1).max()),
            max_original_displacement_m=float(np.linalg.norm(x-x0, axis=1).max()),
            maximum_penetration_m=normal['max_penetration_m'],
            maximum_density_relative_change=float(np.abs(pa.rho/p['rho0']-1).max()),
            temperature_range_K=[float(pa.temperature.min()), float(pa.temperature.max())],
            maximum_internal_acceleration_m_s2=float(np.linalg.norm(force_predictor, axis=1).max()),
            thermal_bound_s=e['thermal']['thermal_dt_bound_s'],
            contact_intervals_total=e['boundary']['ledger']['contact_subintervals'],
            contact_intervals_max_per_particle=maximum_intervals,
            contact_intervals_max_particle_ids=[int(tr['particle_id']) for tr in traces
                if len(tr['result']['intervals']) == maximum_intervals]),
        fixed_scales=dict(bulk_Pa=bulk, longitudinal_speed_m_s=sound,
            acoustic_dx_time_s=p['dx']/sound, acoustic_h_time_s=p['dx']*p['h_ratio']/sound,
            old_explicit_cfl_h_time_s=p['cfl']*p['dx']*p['h_ratio']/sound,
            contact_omega_rad_s=float(contact_omega.max()),
            contact_period_s=float(2*np.pi/contact_omega.max()),
            frozen_elastic_max_omega_from_condition_rad_s=elastic_omega,
            frozen_elastic_shortest_period_s=2*np.pi/elastic_omega,
            frequency_inference='UNPRESTRESSED_FROZEN_ELASTIC; NOT_FULL_NONLINEAR_JACOBIAN'),
        dt_comparisons=[])
    for dt in (5e-7, 5e-6):
        x_predictor = x+dt*v
        endpoint_internal_only = v+dt*force_predictor
        mean_internal_only = v+.5*dt*force_predictor
        theta = elastic_omega*dt
        result['dt_comparisons'].append(dict(dt_s=dt,
            tool_rotation_rad=p['omega']*dt, tool_rotation_deg=p['omega']*dt*180/np.pi,
            axial_command_increment_m=p['feed']*dt,
            ballistic_max_original_displacement_m=float(np.linalg.norm(x_predictor-x0, axis=1).max()),
            old_force_only_endpoint_speed_max_m_s=float(np.linalg.norm(endpoint_internal_only, axis=1).max()),
            old_force_only_mean_speed_max_m_s=float(np.linalg.norm(mean_internal_only, axis=1).max()),
            old_force_only_prediction_status='NOT_ACTUAL_MAP; CONTACT_NORMAL_FRICTION_CONSTRAINT_AND_MATERIAL_FEEDBACK_OMITTED',
            contact_omega_dt=float(contact_omega.max()*dt),
            inferred_elastic_omega_dt=theta,
            midpoint_high_mode_phase_ratio=2*math.atan(theta/2)/theta,
            dt_over_archived_thermal_bound=dt/e['thermal']['thermal_dt_bound_s']))
    result['input_sha256'] = {str(path.relative_to(ROOT)).replace('\\','/'): digest for path,digest in watched.items()}
    result['script_sha256'] = sha(Path(__file__).resolve())
    result['wall_seconds'] = time.perf_counter()-started
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

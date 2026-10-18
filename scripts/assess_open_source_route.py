"""Read-only source freeze and dimensional assessment; never advances particles."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from afsd_backend.pysph.material import FKS

PIN = '940885734daae8e85c1fcd00c56e95a56cc2aa58'
UPSTREAM = ROOT / 'third_party/route_audit_20261004/FE-SPH-GPU'
CASE = ROOT / 'outputs/contact_material_search_20261004'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(
        ['git', '-c', 'safe.directory=' + UPSTREAM.as_posix(),
         '-C', str(UPSTREAM), *args], timeout=20)


def freeze():
    commit = git('rev-parse', 'HEAD').decode().strip()
    if commit != PIN or git('status', '--porcelain').strip():
        raise ValueError('upstream pin or clean checkout mismatch')
    names = git('ls-files', '-z').decode().split('\0')
    entries = {}
    for name in filter(None, names):
        path = UPSTREAM / name
        if not path.resolve().is_relative_to(UPSTREAM.resolve()) or path.is_symlink():
            raise ValueError('unexpected tracked path: ' + name)
        entries[name] = dict(bytes=path.stat().st_size, sha256=sha(path))
    return dict(repository='https://github.com/SPH-SSMP/FE-SPH-GPU',
                commit=commit, files=entries,
                license_file='LICENSE', compiled=False, executed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to((ROOT / 'outputs').resolve()) or output.exists():
        raise ValueError('new output directory inside project outputs required')
    start = time.perf_counter()
    upstream = freeze()
    protected = sorted((ROOT / 'src').rglob('*.py'))
    protected += [CASE / 'root' / name for name in (
        'resolved_config.json', 'accepted.npz', 'acceptance.json',
        'accepted_history.json', 'manifest.json')]
    protected += [CASE / 'time_comparison.json', Path(__file__).resolve()]
    before = {p.relative_to(ROOT).as_posix(): sha(p) for p in protected}
    config = json.loads((CASE / 'root/resolved_config.json').read_text())
    comparison = json.loads((CASE / 'time_comparison.json').read_text())
    manifest = json.loads((CASE / 'root/manifest.json').read_text())
    p = config['parameters']
    law = FKS(**{name: p[name] for name in FKS.__dataclass_fields__})
    length, speed = p['layer'], p['omega'] * p['tool_radius']
    characteristic_rate = speed / length / math.sqrt(3.)
    alpha = p['conductivity'] / (p['rho0'] * p['cp'])
    bulk = 2 * p['shear'] * (1 + p['poisson']) / (3 * (1 - 2*p['poisson']))
    acoustic = math.sqrt((bulk + 4*p['shear']/3) / p['rho0'])
    rate_rows = []
    for rate in (0., .01, .1, 1., 10., characteristic_rate, 100.):
        stress = float(law.yield_stress(0., rate, p['temperature']))
        rate_rows.append(dict(equivalent_rate_s_inv=rate, yield_stress_Pa=stress,
                             secant_viscosity_Pa_s=stress/(3*rate) if rate > 0 else None))
    nominal_gap = .00025  # Documented original nearest centre-to-shoulder gap.
    nominal_time = nominal_gap / p['feed']
    yield_characteristic = float(law.yield_stress(0., characteristic_rate, p['temperature']))
    target_rows = []
    for dt in (.0001, .001, .01):
        steps = math.ceil(nominal_time / dt)
        target_rows.append(dict(dt_s=dt, nominal_steps=steps,
                               rotation_deg_per_step=math.degrees(p['omega']*dt),
                               nominal_feed_m_per_step=p['feed']*dt,
                               wall_budget_s_per_step_excluding_overhead=600/steps))
    readback = {key: sha(ROOT/key) for key in before}
    if before != readback or upstream != freeze():
        raise RuntimeError('audit changed protected source or evidence')
    result = dict(
        classification='READ_ONLY_SOURCE_AND_DIMENSIONAL_ASSESSMENT_NOT_SIMULATION',
        created_utc=datetime.now(timezone.utc).isoformat(),
        upstream=upstream, input_source_sha256=before,
        protected_files_unchanged=True,
        FS_source_audit='NOT_RUN_DOWNLOAD_AUTO_REVIEW_COMPATIBILITY_FAILURE',
        initial_temperature_K=p['temperature'],
        initial_temperature_C=p['temperature']-273.15,
        initial_state_note='Already preheated; hot alternative requires mechanical equilibration, not merely warming.',
        FKS_note='Current DEBUG dimensionless b1 and rate floor, not calibrated AA6061 validation.',
        FKS_ep_zero_rate_samples=rate_rows,
        scales=dict(shoulder_peripheral_speed_m_s=speed,
                    assumed_simple_shear_equivalent_rate_s_inv=characteristic_rate,
                    rho_U2_over_yield=p['rho0']*speed**2/yield_characteristic,
                    nominal_Mach=speed/acoustic,
                    longitudinal_wave_speed_m_s=acoustic,
                    dx_wave_crossing_time_s=p['dx']/acoustic,
                    thermal_diffusivity_m2_s=alpha,
                    layer_diffusion_time_s=length**2/alpha,
                    nominal_rotational_Peclet=speed*length/alpha,
                    nominal_feed_Peclet=p['feed']*length/alpha),
        scales_note='Dimensional estimates using commanded speeds and layer length; not actual gradients, local force balance, or permission to remove inertia.',
        nominal_gap_m=nominal_gap,
        nominal_gap_provenance='docs/PROCESS_COMPLETION_TIME_ASSESSMENT.md (existing geometry assessment, not recomputed crossing)',
        old_displacement_guard_m=p['maximum_displacement_m'],
        nominal_transport_time_s=nominal_time,
        budget_requirements_not_performance_predictions=target_rows,
        naive_repeated_current_step_wall_days=(
            math.ceil(nominal_time/config['settings']['dt_s'])*manifest['wall_seconds']/86400),
        baseline=dict(accepted_time_s=comparison['endpoint_s'],
                      coarse_dt_s=config['settings']['dt_s'],
                      coarse_internal_wall_s=manifest['wall_seconds'],
                      velocity_endpoint_difference_m_s=comparison['endpoint_fields']['velocity_m_s']['maximum_absolute_difference'],
                      new_plastic_strain_change=comparison['endpoint_fields']['ep']['coarse_maximum_change_from_common_parent']),
        checks=dict(no_new_physical_time=True, upstream_not_executed=True,
                    new_backend_not_activated=True, quasi_static_validated=False,
                    process_transport_validated=False),
        wall_seconds=time.perf_counter()-start)
    output.mkdir(parents=True)
    (output / 'assessment.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(output), upstream_files=len(upstream['files']),
                         protected_files=len(before), status=result['classification'],
                         wall_seconds=result['wall_seconds']), ensure_ascii=False))


if __name__ == '__main__':
    main()

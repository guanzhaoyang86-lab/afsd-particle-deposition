"""Registered manufactured operator checks; never calls an AFSD trial.

Use the actual PySPH endpoint NNPS and one visible i/j/g subset for existing
correction, force and heat functions. This is not loaded-process validation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.mechanics import operators, internal_force
from afsd_backend.pysph.quasistatic_process import _geometry
from afsd_backend.pysph.implicit_thermal import pair_conductance, advance_implicit_thermal
from audit_solid_occluded_neighbors import audit_solid_edge_visibility
from visible_support import filter_visible_support


CONTRACT = 'docs/VISIBLE_GRAPH_BENCHMARK_CONTRACT_20261004.md'
CONTRACT_SHA256 = 'e7d594b54e678ad8f8675241fd79f574d259dd035bd24a7a58a7093658218d8c'
SOURCE_PATHS = [CONTRACT, 'scripts/check_visible_support_operators.py',
    'scripts/visible_support.py', 'scripts/audit_solid_occluded_neighbors.py',
    'src/afsd_backend/pysph/mechanics.py', 'src/afsd_backend/pysph/kernel.py',
    'src/afsd_backend/pysph/quasistatic_process.py',
    'src/afsd_backend/pysph/implicit_thermal.py',
    'src/afsd_backend/pysph/coupled_thermal.py']
PARAMETERS = dict(tool_radius=.00145, bore_width=.002, layer=0., tool_top=.0024)
POSE = dict(origin=np.zeros(3), rotation=np.eye(3))
AFFINE = np.array([[.7, -.4, .3], [.2, -.6, .8], [-.9, .5, .1]])
RHO = 2700.
CP = 900.


def source_hashes():
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in SOURCE_PATHS}


def native(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=native,
                               allow_nan=False) + '\n', encoding='utf-8')


def check(name, passed, **metrics):
    return dict(name=name, status='PASS' if bool(passed) else 'FAIL', **metrics)


def maxnorm(a):
    a = np.asarray(a)
    return float(np.linalg.norm(a.reshape(len(a), -1), axis=1).max(initial=0.))


def particle_array(x, volume, h, velocity=None):
    velocity = np.zeros_like(x) if velocity is None else velocity
    return get_particle_array(name='manufactured_visible_support',
        x=x[:, 0], y=x[:, 1], z=x[:, 2], h=np.full(len(x), h),
        m=RHO*volume, rho=np.full(len(x), RHO),
        u=velocity[:, 0], v=velocity[:, 1], w=velocity[:, 2])


def fixture(n):
    length = .0008
    spacing = length/(n-1)
    axes = [np.linspace(.0001, .0009, n), np.linspace(-.0004, .0004, n),
            np.linspace(.0008, .0016, n)]
    mesh = np.meshgrid(*axes, indexing='ij')
    first = np.stack([a.ravel() for a in mesh], axis=1)
    second = first.copy(); second[:, 0] += .0014
    x = np.concatenate((first, second))
    q = np.ones(n); q[[0, -1]] = .5
    wmesh = np.meshgrid(q, q, q, indexing='ij')
    weights = np.prod(np.stack(wmesh), axis=0).ravel()*spacing**3
    volume = np.tile(weights, 2)
    block = np.repeat(np.arange(2), n**3)
    integer_mesh = np.stack(np.meshgrid(*[np.arange(n)]*3, indexing='ij'), axis=-1).reshape(-1, 3)
    boundary = np.tile(np.any((integer_mesh == 0) | (integer_mesh == n-1), axis=1), 2)
    centers = np.array([[.0005, 0., .0012], [.0019, 0., .0012]])
    return x, volume, 1.7*spacing, block, boundary, x-centers[block]


def graph(pa, pose=POSE):
    i, j, g, _, _ = operators(pa)
    x = np.column_stack((pa.x, pa.y, pa.z))
    visible = filter_visible_support(x, i, j, g, pose, PARAMETERS,
                                     np.arange(len(x), dtype=np.int64))
    return (i, j, g), visible


def geometry(x, v, pa, graph_data):
    if isinstance(graph_data, dict):
        i, j, g = (graph_data[k] for k in ('i', 'j', 'g'))
    else:
        i, j, g = graph_data
    return _geometry(x, v, pa.m, pa.rho, i, j, g)


def force(pa, stress, visible, B):
    return internal_force(pa, stress, visible['i'], visible['j'], visible['g'], correction=B)


def heat(pa, x, graph_data, temperature, plastic=None):
    i, j, g = (graph_data[k] for k in ('i', 'j', 'g')) if isinstance(graph_data, dict) else graph_data
    n = len(x)
    return advance_implicit_thermal(pa.m, pa.rho, temperature, 150., x, i, j, g,
        .0002, CP, np.zeros(n) if plastic is None else plastic,
        np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), .9, .5, 273., 900.)


def stress_field(x):
    base = np.array([[2., .4, -.2], [.4, -1., .3], [-.2, .3, .7]])
    factor = 1.+.3*np.sin(x[:, 0]/.0008)+.2*np.cos(x[:, 2]/.0008)
    stress = 1e6*factor[:, None, None]*base
    stress[:, 0, 0] += 1e5*np.sin(x[:, 1]/.0004)
    stress[:, 2, 2] += 2e5*np.cos(x[:, 0]/.0007)
    return stress


def pairs_sorted(i, j):
    a = np.column_stack((i, j))
    return a[np.lexsort((a[:, 1], a[:, 0]))]


def resolution_case(n):
    start = time.perf_counter()
    x, volume, h, block, boundary, xi = fixture(n)
    pa = particle_array(x, volume, h)
    raw, visible = graph(pa)
    i, j, g = raw
    keep = visible['kept_raw_edge_indices']
    expected = block[i] == block[j]
    subset_ok = (np.array_equal(keep, np.flatnonzero(expected)) and
        all(np.array_equal(visible[k], a[keep]) for k, a in zip(('i', 'j', 'g'), raw)))
    geometric = audit_solid_edge_visibility(x, visible['i'], visible['j'], POSE, PARAMETERS)
    checks = [check('one_common_visible_subset', subset_ok and np.any(~expected) and
                    geometric['occluded_edge_count'] == 0,
                    raw_edge_count=len(i), raw_cross_block_count=int((~expected).sum()),
                    visible_edge_count=len(keep), blocked_edge_count=len(i)-len(keep),
                    roundoff_only_blocked_edge_count=visible['diagnostics']['roundoff_only_blocked_edge_count'],
                    expected_mask_exact=subset_ok,
                    retained_old_guard_maximum_length_m=geometric['max_occluded_interval_length_m'])]
    v = x@AFFINE.T+np.array([.03, -.02, .01])
    L, B = geometry(x, v, pa, visible)
    affine_error = maxnorm(L-AFFINE)
    translation, _ = geometry(x, np.tile([.03, -.02, .01], (len(x), 1)), pa, visible)
    translation_error = maxnorm(translation)
    omega = np.array([.3, -.5, .7])
    spin = np.array([[0., -omega[2], omega[1]], [omega[2], 0., -omega[0]], [-omega[1], omega[0], 0.]])
    Lspin, _ = geometry(x, np.cross(omega, x), pa, visible)
    spin_error = maxnorm(Lspin-spin)
    strain_error = maxnorm(.5*(Lspin+np.swapaxes(Lspin, 1, 2)))
    checks.append(check('affine_translation_rigid_rotation', affine_error <= 1e-10 and
        translation_error <= 1e-12 and spin_error <= 1e-10 and strain_error <= 1e-10,
        affine_max_frobenius_s_inv=affine_error, translation_max_frobenius_s_inv=translation_error,
        spin_max_frobenius_s_inv=spin_error, rigid_strain_max_frobenius_s_inv=strain_error))

    stress = stress_field(x)
    velocity = v+.001*np.column_stack((np.sin(x[:, 1]/.0004),
        np.cos(x[:, 0]/.0008), np.sin(x[:, 2]/.0005)))
    Lv, B = geometry(x, velocity, pa, visible)
    F = force(pa, stress, visible, B)
    sum_force = np.sum(F, axis=0)
    torque_terms = np.cross(x, F)
    sum_torque = np.sum(torque_terms, axis=0)
    particle_power = np.einsum('ni,ni->n', velocity, F)
    stress_power = volume*np.einsum('nij,nij->n', stress, Lv)
    power_residual = math.fsum(particle_power)+math.fsum(stress_power)
    force_limit = 1e-12+1e-10*float(np.linalg.norm(F, axis=1).sum())
    torque_limit = 1e-14+1e-10*float(np.linalg.norm(torque_terms, axis=1).sum())
    power_limit = 1e-12+1e-10*float(np.abs(particle_power).sum()+np.abs(stress_power).sum())
    checks.append(check('momentum_torque_and_virtual_power', np.linalg.norm(sum_force) <= force_limit and
        np.linalg.norm(sum_torque) <= torque_limit and abs(power_residual) <= power_limit,
        total_force_N=sum_force, force_residual_norm_N=float(np.linalg.norm(sum_force)), force_limit_N=force_limit,
        total_torque_N_m=sum_torque, torque_residual_norm_N_m=float(np.linalg.norm(sum_torque)), torque_limit_N_m=torque_limit,
        virtual_power_residual_W=power_residual, virtual_power_limit_W=power_limit,
        particle_virtual_power_W=math.fsum(particle_power), stress_virtual_power_W=math.fsum(stress_power)))

    imposed_velocity = velocity.copy(); imposed_velocity[block == 1] = 0.
    isolated_L, isolated_B = geometry(x, imposed_velocity, pa, visible)
    raw_L, raw_B = geometry(x, imposed_velocity, pa, raw)
    imposed_stress = stress.copy(); imposed_stress[block == 1] = 0.
    isolated_F = force(pa, imposed_stress, visible, isolated_B)
    raw_F = internal_force(pa, imposed_stress, *raw, correction=raw_B)
    isolated_gradient = maxnorm(isolated_L[block == 1])
    isolated_force = maxnorm(isolated_F[block == 1])
    raw_gradient = maxnorm(raw_L[block == 1]); raw_force = maxnorm(raw_F[block == 1])
    coefficient = pair_conductance(pa.m, pa.rho, 150., x, visible['i'], visible['j'], visible['g'])
    cross_heat_count = int(np.count_nonzero((block[visible['i']] != block[visible['j']]) & (coefficient != 0)))
    block_temperature = np.where(block == 0, 400., 500.)
    constant_heat = heat(pa, x, visible, block_temperature)
    raw_heat = heat(pa, x, raw, block_temperature)
    constant_error = float(np.max(np.abs(constant_heat['temperature']-block_temperature)))
    raw_temperature_change = float(np.max(np.abs(raw_heat['temperature']-block_temperature)))
    checks.append(check('no_direct_cross_wall_mechanical_or_heat_coupling',
        isolated_gradient == 0. and isolated_force == 0. and raw_gradient > 0. and raw_force > 0. and
        cross_heat_count == 0 and constant_error <= 1e-10 and raw_temperature_change > 1e-5,
        visible_other_block_gradient_max_s_inv=isolated_gradient,
        visible_other_block_force_max_N=isolated_force, raw_other_block_gradient_max_s_inv=raw_gradient,
        raw_other_block_force_max_N=raw_force, cross_heat_nonzero_coefficient_count=cross_heat_count,
        visible_block_uniform_temperature_max_error_K=constant_error,
        raw_temperature_max_change_K=raw_temperature_change))

    T0 = 450.+25.*np.sin(xi[:, 0]/.0004)+15.*np.cos(xi[:, 1]/.0004)+7.*np.sin(xi[:, 2]/.0004)
    plastic = 2000.*volume
    heated = heat(pa, x, visible, T0, plastic)
    signed = pa.m*CP*(heated['temperature']-T0)-.9*plastic
    rows = []
    for label, mask in [('inside_block', block == 0), ('outside_block', block == 1), ('all', np.ones(len(x), bool))]:
        residual = math.fsum(signed[mask]); supplied = math.fsum(.9*plastic[mask])
        limit = 1e-12+1e-10*abs(supplied)
        rows.append(dict(scope=label, residual_J=residual, source_heat_J=supplied,
                         limit_J=limit, passed=abs(residual) <= limit))
    checks.append(check('implicit_heat_capacity_balance_by_disconnected_block', all(a['passed'] for a in rows),
        scopes=rows, temperature_min_K=float(heated['temperature'].min()),
        temperature_max_K=float(heated['temperature'].max()), ledger=heated['ledger'],
        algebraic_residual_max_J=float(np.max(np.abs(heated['linear_system_residual_J'])))))

    axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
    skew = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
    Q = np.eye(3)+np.sin(.37)*skew+(1-np.cos(.37))*(skew@skew)
    shift = np.array([.008, -.004, .003])
    rotated_x = x@Q.T+shift
    rotated_velocity = velocity@Q.T
    rotated_stress = Q@stress@Q.T
    rotated_pa = particle_array(rotated_x, volume, h, rotated_velocity)
    rotated_raw, rotated_visible = graph(rotated_pa, dict(origin=shift, rotation=Q))
    rotated_L, rotated_B = geometry(rotated_x, rotated_velocity, rotated_pa, rotated_visible)
    rotated_F = force(rotated_pa, rotated_stress, rotated_visible, rotated_B)
    same_raw_pairs = np.array_equal(pairs_sorted(*raw[:2]), pairs_sorted(*rotated_raw[:2]))
    same_visible_pairs = np.array_equal(pairs_sorted(visible['i'], visible['j']),
                                        pairs_sorted(rotated_visible['i'], rotated_visible['j']))
    covariant_L_error = maxnorm(rotated_L-Q@Lv@Q.T)
    covariant_F_error = maxnorm(rotated_F-F@Q.T)
    covariant_L_limit = 1e-12+1e-10*maxnorm(Lv)
    covariant_F_limit = 1e-12+1e-10*maxnorm(F)
    checks.append(check('rigid_coordinate_transform_covariance', same_raw_pairs and same_visible_pairs and
        covariant_L_error <= covariant_L_limit and covariant_F_error <= covariant_F_limit,
        raw_pair_set_identical=same_raw_pairs, visible_pair_set_identical=same_visible_pairs,
        gradient_max_frobenius_error_s_inv=covariant_L_error, gradient_limit_s_inv=covariant_L_limit,
        force_max_error_N=covariant_F_error, force_limit_N=covariant_F_limit,
        rotation=Q, translation_m=shift,
        geometric_roundoff_expansion_before_m=visible['diagnostics']['roundoff_solid_expansion_m'],
        geometric_roundoff_expansion_after_m=rotated_visible['diagnostics']['roundoff_solid_expansion_m']))

    quadratic = 1000.*xi**2
    quadratic_L, _ = geometry(x, quadratic, pa, visible)
    exact = np.zeros_like(quadratic_L)
    exact[:, np.arange(3), np.arange(3)] = 2000.*xi
    square_error = np.sum((quadratic_L-exact)**2, axis=(1, 2))
    errors = {}
    for label, mask in [('all', np.ones(len(x), bool)), ('boundary', boundary), ('interior', ~boundary)]:
        errors[label+'_volume_weighted_RMS_s_inv'] = float(np.sqrt(np.sum(volume[mask]*square_error[mask])/np.sum(volume[mask])))
        errors[label+'_max_frobenius_s_inv'] = float(np.sqrt(square_error[mask].max(initial=0.)))
        errors[label+'_particle_count'] = int(mask.sum())
    return dict(n_per_axis=n, particle_count=len(x), physical_cube_volume_m3=.0008**3,
        total_volume_m3=float(volume.sum()), total_mass_kg=float(pa.m.sum()), h_m=h,
        wall_seconds=time.perf_counter()-start, checks=checks, quadratic_gradient_errors=errors,
        graph_fingerprints=visible['fingerprints'], production_authorized=False)


def degeneracy_case():
    yz = np.stack(np.meshgrid(np.linspace(-.0002, .0002, 3),
                             np.linspace(.0009, .0013, 3), indexing='ij'), axis=-1).reshape(-1, 2)
    x = np.concatenate([np.column_stack((np.full(len(yz), a), yz)) for a in (.0009, .0015)])
    pa = particle_array(x, np.full(len(x), .0002**3), .0006)
    raw, visible = graph(pa)
    raw_error = None; filtered_error = None
    raw_affine_error = None
    try:
        L, _ = geometry(x, x@AFFINE.T, pa, raw)
        raw_affine_error = maxnorm(L-AFFINE)
    except Exception as exc:
        raw_error = f'{type(exc).__name__}: {exc}'
    try:
        geometry(x, x@AFFINE.T, pa, visible)
    except ValueError as exc:
        filtered_error = str(exc)
    return check('filtered_support_degeneracy_rejected_explicitly', raw_error is None and
        raw_affine_error <= 1e-10 and filtered_error is not None and
        'singular/ill-conditioned' in filtered_error,
        raw_edge_count=len(raw[0]), visible_edge_count=len(visible['i']),
        raw_geometry_error=raw_error, raw_affine_max_error_s_inv=raw_affine_error,
        filtered_exception=filtered_error)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='outputs/visible_graph_20261004/operators_v1')
    args = parser.parse_args()
    output = (ROOT / args.output).resolve()
    if not output.is_relative_to(ROOT):
        raise ValueError('output must remain in project')
    output.mkdir(parents=True, exist_ok=False)
    before = source_hashes()
    if before[CONTRACT] != CONTRACT_SHA256:
        raise ValueError('registered benchmark contract changed')
    plan = dict(classification='MANUFACTURED_VISIBLE_GRAPH_OPERATOR_BENCHMARK_NOT_AFSD_VALIDATION',
        command=[sys.executable, *sys.argv], source_sha256=before,
        registered_resolutions=[3, 5], parameters=PARAMETERS, physical_solver_evaluations=0,
        tests_frozen_before_run=True, production_authorized=False)
    write_json(output / 'execution_plan.json', plan)
    start = time.perf_counter()
    cases = []; exceptions = []
    for n in (3, 5):
        try:
            case = resolution_case(n); cases.append(case)
            write_json(output / f'resolution_{n}.json', case)
        except Exception:
            exceptions.append(dict(case=f'resolution_{n}', exception=traceback.format_exc()))
    top_checks = []
    if len(cases) == 2:
        coarse = cases[0]['quadratic_gradient_errors']['all_volume_weighted_RMS_s_inv']
        fine = cases[1]['quadratic_gradient_errors']['all_volume_weighted_RMS_s_inv']
        ratio = fine/coarse
        top_checks.append(check('same_domain_quadratic_gradient_error_reduction', fine < .8*coarse,
            coarse_RMS_s_inv=coarse, fine_RMS_s_inv=fine, fine_to_coarse_ratio=ratio,
            required_ratio_strictly_below=.8,
            total_mass_difference_kg=cases[1]['total_mass_kg']-cases[0]['total_mass_kg'],
            total_volume_difference_m3=cases[1]['total_volume_m3']-cases[0]['total_volume_m3']))
    try:
        top_checks.append(degeneracy_case())
    except Exception:
        exceptions.append(dict(case='degeneracy', exception=traceback.format_exc()))
    after = source_hashes()
    top_checks.append(check('all_registered_sources_unchanged', before == after,
                             changed=[p for p in before if before[p] != after[p]]))
    all_checks = top_checks + [item for case in cases for item in case['checks']]
    passed = not exceptions and len(cases) == 2 and all(c['status'] == 'PASS' for c in all_checks)
    exit_code = 0 if passed else 2
    result = dict(**plan, status='PASS_REGISTERED_COMPONENT_BENCHMARKS_ONLY' if passed else 'FAIL_COMPONENT_BENCHMARKS',
        exit_code=exit_code, wall_seconds=time.perf_counter()-start, resolutions=cases,
        checks=top_checks, exceptions=exceptions, check_count=len(all_checks),
        failed_checks=[c['name'] for c in all_checks if c['status'] != 'PASS'], source_sha256_after=after,
        limitations=['Manufactured fields and fixed graphs, not any loaded AFSD simulation.',
            'Derivative-error refinement does not verify material flow, finite-step energy or deposition.',
            'Moving graph discontinuity, free-surface traction consistency and constitutive assumptions remain unverified.'])
    write_json(output / 'summary.json', result)
    print(json.dumps(dict(status=result['status'], exit_code=exit_code, output=str(output),
        wall_seconds=result['wall_seconds'], check_count=len(all_checks), failed_checks=result['failed_checks'],
        exceptions=exceptions), indent=2))
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())

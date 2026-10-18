"""One preregistered interior-point perturbation; no AFSD or production edit."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

from probe_boundary_sbp import fixture, extract_original, project, adjoint_force, STRESS

ROOT = Path(__file__).resolve().parents[1]
ELL_M = 1.
MOMENT_TOL_M2 = 1e-10
SHIFT_M = np.array([.025, -.015, .010])
CONFIG = dict(classification='STATIC_SAMPLING_FIXTURE_NOT_PHYSICAL_STATE',
              bounds_m=[-1., 1.], nodes_per_axis=5, dx_m=.5, h_m=.65,
              rho_kg_m3=1., moved_original_position_m=[0., 0., 0.],
              displacement_m=SHIFT_M.tolist(), weights_frozen=True,
              ell_m=ELL_M, moment_tolerance_m2=MOMENT_TOL_M2,
              production_authorized=False)


class IncompatibleMoments(ValueError):
    """No operator can meet the declared frozen H/E and affine conditions."""


def moment_defect(x, volume, boundary):
    # Dimensionless P; every entry of C has units m^2.
    P = np.column_stack((np.ones(len(x)), x/ELL_M))
    C = np.zeros((3, 4, 4))
    for k in range(3):
        dP = np.zeros_like(P)
        dP[:, k+1] = 1./ELL_M
        C[k] = (P.T@(volume[:, None]*dP) + dP.T@(volume[:, None]*P)
                - P.T@(boundary[k, :, None]*P))
    return C


def checked_project(x, volume, boundary, D0, i, j):
    C = moment_defect(x, volume, boundary)
    if not np.isfinite(C).all():
        raise IncompatibleMoments('nonfinite volume-boundary moment defect')
    error = float(np.max(np.abs(C)))
    if error > MOMENT_TOL_M2:
        raise IncompatibleMoments(
            f'frozen volume-boundary moments incompatible: {error:.17g} m^2 > '
            f'{MOMENT_TOL_M2:.17g} m^2; projection NOT_RUN')
    return project(x, volume, boundary, D0, i, j)


def run_case(perturb):
    pa, x, volume, boundary = fixture()
    original = x.copy()
    center = np.flatnonzero(np.all(x == 0., axis=1))
    if len(center) != 1:
        raise ValueError('expected exactly one original center')
    if perturb:
        x[center[0]] += SHIFT_M
        pa.x[:], pa.y[:], pa.z[:] = x.T
    # Necessary moments first; extraction is read-only condition/affine diagnosis.
    C = moment_defect(x, volume, boundary)
    D0, (i, j, _, _), condition = extract_original(pa)
    target = np.zeros((3, len(x), 4))
    for k in range(3):
        target[k, :, k+1] = 1./ELL_M
    P = np.column_stack((np.ones(len(x)), x/ELL_M))
    max_defect = float(np.max(np.abs(C)))
    r = dict(perturbed=perturb, particle_count=len(x), shifted_index=int(center[0]),
             volume_m3=float(volume.sum()),
             first_volume_moment_m4=(volume@x).tolist(),
             quadrature_centroid_m=(volume@x/volume.sum()).tolist(),
             moment_defect_m2=C.tolist(), max_moment_defect_m2=max_defect,
             moment_tolerance_m2=MOMENT_TOL_M2,
             original_max_moment_condition=condition, nnps_edges=len(i),
             original_affine_max=float(np.max(np.abs(D0@P-target))),
             zero_mode_test='NOT_RUN', production_authorized=False)
    states = dict(original_position_m=original, position_m=x, volume_m3=volume,
                  mass_kg=pa.m.copy(), rho_kg_m3=pa.rho.copy(), h_m=pa.h.copy(),
                  signed_face_area_m2=boundary, moment_defect_m2=C,
                  original_D_per_m=D0, neighbor_i=i, neighbor_j=j)
    try:
        D, stats, singular = checked_project(x, volume, boundary, D0, i, j)
    except IncompatibleMoments as exc:
        r.update(method_status='BLOCKED_BY_MOMENT_COMPATIBILITY',
                 projection_status='NOT_RUN', rejection_reason=str(exc),
                 accepted_operator=False)
        states['operator_status'] = np.array('REJECTED_FOR_OPERATOR_CONSTRUCTION')
    else:
        Q = volume[None, :, None]*D
        E = np.array([np.diag(e) for e in boundary])
        ext = np.einsum('ki,ak->ia', boundary, STRESS)
        fi = adjoint_force(D, volume, np.broadcast_to(STRESS, (len(x), 3, 3)))
        checks = dict(sbp_max=float(np.max(np.abs(Q+Q.transpose(0, 2, 1)-E))),
                      affine_max=float(np.max(np.abs(D@P-target))),
                      constant_patch_relative=float(np.linalg.norm(fi+ext)/np.linalg.norm(ext)),
                      constraint_compatibility_max=stats['compatibility_max'],
                      constraint_residual_max=stats['linear_constraint_max'])
        accepted = all(value <= 1e-10 for value in checks.values())
        r.update(method_status='FIXTURE_PASS' if accepted else 'FIXTURE_FAIL',
                 projection_status='RUN', checks=checks, projection=stats,
                 accepted_operator=accepted)
        states['candidate_D_per_m'] = D
        states['constraint_singular_values'] = singular
        states['operator_status'] = np.array('ACCEPTED_STATIC_FIXTURE' if accepted else 'REJECTED')
    return r, states


def evaluate():
    start = time.perf_counter()
    control, control_state = run_case(False)
    perturbed, perturbed_state = run_case(True)
    report = dict(status=perturbed['method_status'],
                  classification=CONFIG['classification'], diagnostic_completed=True,
                  control=control, perturbed=perturbed,
                  wall_seconds=time.perf_counter()-start,
                  production_authorized=False, afsd='NOT_RUN',
                  physically_accepted_terminal_state='NONE',
                  actual_call_order=['fixture', 'moment_defect', 'extract_original',
                                     'checked_project: reject incompatible or call step2 project'])
    return report, control_state, perturbed_state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    out = Path(parser.parse_args().output)
    out.mkdir(parents=True, exist_ok=False)
    report, control, perturbed = evaluate()
    configuration = json.dumps(CONFIG, sort_keys=True, indent=2)+'\n'
    report['configuration_sha256'] = hashlib.sha256(configuration.encode()).hexdigest()
    files = [Path(__file__), ROOT/'scripts/probe_boundary_sbp.py',
             ROOT/'scripts/probe_boundary_closure.py',
             ROOT/'tests/review/test_boundary_closure_step3.py',
             ROOT/'logs/boundary_step3_20261001/contract_pretest.md']
    files.extend(sorted((ROOT/'src/afsd_backend/pysph').glob('*.py')))
    report['source_and_contract_sha256'] = {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    np.savez_compressed(out/'control_sampling_fixture.npz', **control)
    np.savez_compressed(out/'REJECTED_sampling_fixture.npz', **perturbed)
    # Match the hashed UTF-8 bytes exactly, including on Windows.
    (out/'fixture_config.json').write_text(configuration, encoding='utf-8', newline='\n')
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return 0 if report['control']['accepted_operator'] and report['perturbed']['accepted_operator'] else 1


if __name__ == '__main__':
    raise SystemExit(main())


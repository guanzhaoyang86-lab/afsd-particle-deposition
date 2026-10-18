"""Fixed 5^3 SBP feasibility prototype; NEVER imported by the AFSD solver."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.mechanics import corrected_operators, internal_force
from probe_boundary_closure import bubble, STRESS

ALGEBRA_TOL = 1e-10
RANK_RTOL = 1e-12
NULL_RTOL = 1e-10


def fixture():
    """Independent volume/face quadrature, not inferred from internal forces."""
    index = np.array(list(itertools.product(range(5), repeat=3)))
    x = np.linspace(-1., 1., 5)[index]
    weights = .5*np.where((index == 0) | (index == 4), .5, 1.)
    volume = np.prod(weights, axis=1)
    boundary = np.zeros((3, len(x)))
    for axis in range(3):
        area = np.prod(weights[:, [k for k in range(3) if k != axis]], axis=1)
        boundary[axis, index[:, axis] == 0] = -area[index[:, axis] == 0]
        boundary[axis, index[:, axis] == 4] = area[index[:, axis] == 4]
    pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], m=volume,
                            rho=np.ones(len(x)), h=np.full(len(x), .65))
    return pa, x, volume, boundary


def extract_original(pa):
    """D[k] v[a] is L[a,k]; preserve actual NNPS edges and moment guard."""
    i, j, g, _, _, B = corrected_operators(pa)
    order = np.lexsort((j, i))
    i, j, g = i[order], j[order], g[order]
    volume = pa.m/pa.rho
    gi = np.einsum('nji,nj->ni', B[i], g)
    gj = np.einsum('nji,nj->ni', B[j], g)
    D = np.zeros((3, len(volume), len(volume)))
    for k in range(3):
        np.add.at(D[k], (i, j), volume[j]*gi[:, k])
        np.add.at(D[k], (i, i), -volume[j]*gi[:, k])
        np.add.at(D[k], (j, i), -volume[i]*gj[:, k])
        np.add.at(D[k], (j, j), volume[i]*gj[:, k])
    eigen = np.linalg.eigvalsh(np.linalg.inv(B))
    return D, (i, j, g, B), float(np.max(eigen[:, -1]/eigen[:, 0]))


def gradient(D, velocity):
    return np.einsum('kij,ja->iak', D, velocity)


def adjoint_force(D, volume, stress):
    return -np.einsum('kij,i,iak->ja', D, volume, stress)


def project(x, volume, boundary, D0, i, j):
    """One explicit rank-revealed projection on the ORIGINAL neighbor graph.

    SVD resolves known global affine/skew constraint dependencies. It is not a
    pseudoinverse of particle moment matrices, and compatibility is tested.
    """
    count = len(x)
    P = np.column_stack((np.ones(count), x))
    A = np.zeros((4*count, len(i)))
    edge = np.arange(len(i))
    for a in range(4):
        A[4*i+a, edge] = P[j, a]
        A[4*j+a, edge] = -P[i, a]
    Q0 = volume[None, :, None]*D0
    S0 = .5*(Q0 - Q0.transpose(0, 2, 1))
    E = np.array([np.diag(row) for row in boundary])
    base = S0 + .5*E
    target = np.zeros((3, count, 4))
    for k in range(3):
        target[k, :, k+1] = volume
    rhs = (target - base@P).reshape(3, -1).T
    U, singular, VT = np.linalg.svd(A, full_matrices=False)
    rank = int(np.sum(singular > RANK_RTOL*singular[0]))
    incompatible = rhs - U[:, :rank]@(U[:, :rank].T@rhs)
    delta = VT[:rank].T@((U[:, :rank].T@rhs)/singular[:rank, None])
    Q = base.copy()
    for k in range(3):
        Q[k, i, j] += delta[:, k]
        Q[k, j, i] -= delta[:, k]
    D = Q/volume[None, :, None]
    stats = dict(constraint_rows=len(A), edge_unknowns=len(i), rank=rank,
                 discarded_left_dimensions=len(A)-rank,
                 singular_threshold=float(RANK_RTOL*singular[0]),
                 resolved_condition=float(singular[0]/singular[rank-1]),
                 compatibility_max=float(np.max(np.abs(incompatible))),
                 linear_constraint_max=float(np.max(np.abs(A@delta-rhs))),
                 relative_Q_change=float(np.linalg.norm(Q-Q0)/np.linalg.norm(Q0)),
                 relative_D_change=float(np.linalg.norm(D-D0)/np.linalg.norm(D0)))
    return D, stats, singular


def mode_spectra(D):
    count = D.shape[1]
    scalar = np.linalg.svd(D.reshape(3*count, count), compute_uv=False)
    # Frobenius norm of symmetric gradient: off-diagonal entries count twice.
    strain = np.zeros((6*count, 3*count))
    for a in range(3):
        strain[a*count:(a+1)*count, a*count:(a+1)*count] = D[a]
    for row, (a, b) in enumerate(((0, 1), (0, 2), (1, 2)), start=3):
        strain[row*count:(row+1)*count, a*count:(a+1)*count] = D[b]/np.sqrt(2.)
        strain[row*count:(row+1)*count, b*count:(b+1)*count] = D[a]/np.sqrt(2.)
    vector = np.linalg.svd(strain, compute_uv=False)
    return scalar, vector


def evaluate():
    start = time.perf_counter()
    pa, x, volume, boundary = fixture()
    D0, (i, j, g, B), moment_condition = extract_original(pa)
    D, projection, singular = project(x, volume, boundary, D0, i, j)
    E = np.array([np.diag(row) for row in boundary])
    P = np.column_stack((np.ones(len(x)), x))
    target = np.zeros((3, len(x), 4))
    for k in range(3):
        target[k, :, k+1] = 1.
    rng = np.random.default_rng(20260930)
    v = rng.normal(size=x.shape)
    stress = rng.normal(size=(len(x), 3, 3))
    stress = .5*(stress+stress.transpose(0, 2, 1))
    pa.u[:], pa.v[:], pa.w[:] = v.T
    _, _, _, actual_L, _, _ = corrected_operators(pa)
    actual_force = internal_force(pa, stress, i, j, g, correction=B)
    constant = np.broadcast_to(STRESS, stress.shape)
    ext = np.einsum('ki,ak->ia', boundary, STRESS)
    omega = np.array([.3, -.2, .5])
    rotation = np.cross(omega, x)
    f, q = rng.normal(size=(2, len(x)))
    support = np.eye(len(x), dtype=bool)
    support[i, j] = support[j, i] = True
    records = {}
    states = {}
    for name, operator in (('original', D0), ('prototype', D)):
        Q = volume[None, :, None]*operator
        L = gradient(operator, v)
        force = adjoint_force(operator, volume, stress)
        patch_force = adjoint_force(operator, volume, constant)
        Lrotation = gradient(operator, rotation)
        scalar, vector = mode_spectra(operator)
        powers = {}
        for field, velocity, exact in (('affine', x, 48.), ('bubble', bubble(x)[0], 0.)):
            power = float(np.einsum('i,iak,ak', volume, gradient(operator, velocity), STRESS))
            powers[field] = dict(volume_W=power, surface_W=float(np.sum(ext*velocity)),
                                 independent_exact_W=exact, error_W=power-exact)
        records[name] = dict(
            sbp_max=float(np.max(np.abs(Q+Q.transpose(0, 2, 1)-E))),
            affine_max=float(np.max(np.abs(operator@P-target))),
            constant_patch_relative=float(np.linalg.norm(patch_force+ext)/np.linalg.norm(ext)),
            power_conjugacy_W=float(np.einsum('i,iak,iak', volume, stress, L)+np.sum(v*force)),
            arbitrary_ibp_max=float(max(abs(f@(volume*(d@q))+(d@f)@(volume*q)-f@(e*q))
                                        for d, e in zip(operator, boundary))),
            total_internal_force_N=float(np.max(np.abs(force.sum(axis=0)))),
            total_internal_torque_Nm=float(np.max(np.abs(np.cross(x, force).sum(axis=0)))),
            rotation_symmetric_gradient_max=float(np.max(np.abs(.5*(Lrotation+Lrotation.transpose(0, 2, 1))))),
            rotation_nonzero_stress_power_W=float(np.einsum('i,iak,iak', volume, stress, Lrotation)),
            outside_graph_max=float(np.max(np.abs(operator[:, ~support]))),
            scalar_nullity=int(np.sum(scalar <= NULL_RTOL*scalar[0])),
            strain_nullity=int(np.sum(vector <= NULL_RTOL*vector[0])),
            powers=powers)
        states[name+'_D_per_m'] = operator
        states[name+'_constant_stress_internal_force_N'] = patch_force
        states[name+'_scalar_gradient_singular_values'] = scalar
        states[name+'_symmetric_gradient_singular_values'] = vector
    bridge = dict(gradient_max=float(np.max(np.abs(gradient(D0, v)-actual_L))),
                  force_max_N=float(np.max(np.abs(adjoint_force(D0, volume, stress)-actual_force))))
    r = records['prototype']
    checks = {
        'original_interface': max(bridge.values()) <= ALGEBRA_TOL,
        'moment_guard': moment_condition <= 1e6,
        'constraint_compatibility': projection['compatibility_max'] <= ALGEBRA_TOL,
        'linear_constraints': projection['linear_constraint_max'] <= ALGEBRA_TOL,
        'sbp': r['sbp_max'] <= ALGEBRA_TOL,
        'affine': r['affine_max'] <= ALGEBRA_TOL,
        'constant_patch': r['constant_patch_relative'] <= ALGEBRA_TOL,
        'power': abs(r['power_conjugacy_W']) <= ALGEBRA_TOL,
        'arbitrary_ibp': r['arbitrary_ibp_max'] <= ALGEBRA_TOL,
        'total_force': r['total_internal_force_N'] <= ALGEBRA_TOL,
        'total_torque': r['total_internal_torque_Nm'] <= ALGEBRA_TOL,
        'rotation': max(r['rotation_symmetric_gradient_max'], abs(r['rotation_nonzero_stress_power_W'])) <= ALGEBRA_TOL,
        'support': r['outside_graph_max'] == 0.,
        'scalar_nullity': r['scalar_nullity'] == 1,
        'strain_nullity': r['strain_nullity'] == 6,
        'virtual_powers': all(abs(row['error_W']) <= ALGEBRA_TOL for row in r['powers'].values()),
    }
    states.update(position_m=x, volume_m3=volume, signed_face_area_m2=boundary,
                  external_constant_stress_traction_N=ext, neighbor_i=i, neighbor_j=j,
                  constraint_singular_values=singular)
    report = dict(status='FIXED_CUBE_PROTOTYPE_PASS' if all(checks.values()) else 'FIXED_CUBE_PROTOTYPE_FAIL',
                  classification='VIRTUAL_WORK_FIXTURE / APPROXIMATE_INDEPENDENT_PROTOTYPE',
                  particles=len(x), physical_bounds_m=[-1., 1.], dx_m=.5, h_m=.65,
                  volume_m3=float(volume.sum()), original_moment_max_condition=moment_condition,
                  bridge=bridge, projection=projection, operators=records, checks=checks,
                  wall_seconds=time.perf_counter()-start, production_authorized=False,
                  dynamic_afsd='NOT_RUN', distorted_or_moving_geometry='NOT_RUN',
                  material_change='NOT_IMPLEMENTED', production_spatial_change='NOT_IMPLEMENTED')
    return report, states


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    out = Path(parser.parse_args().output)
    out.mkdir(parents=True, exist_ok=False)
    report, state = evaluate()
    files = [Path(__file__), ROOT/'scripts/probe_boundary_closure.py',
             ROOT/'tests/review/test_boundary_closure_step2.py',
             ROOT/'logs/boundary_step2_20260930/contract_pretest.md']
    files.extend(sorted((ROOT/'src/afsd_backend/pysph').glob('*.py')))
    report['source_and_contract_sha256'] = {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    np.savez_compressed(out/'fixed_cube_operators.npz', **state)
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return 0 if all(report['checks'].values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())


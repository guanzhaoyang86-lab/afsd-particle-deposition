"""P0 Table 1 pair operators, using real PySPH 3D linked-list neighbors."""
import numpy as np
from cyarray.carray import UIntArray
from pysph.base.nnps import LinkedListNNPS
from .kernel import HyperbolicSpline


def positions(pa):
    return np.column_stack((pa.x, pa.y, pa.z))


def velocities(pa):
    return np.column_stack((pa.u, pa.v, pa.w))


def pairs(pa):
    nnps = LinkedListNNPS(dim=3, particles=[pa], radius_scale=2.0)
    neighbors = UIntArray()
    ii, jj = [], []
    for i in range(pa.get_number_of_particles()):
        nnps.get_nearest_particles(0, 0, i, neighbors)
        js = neighbors.get_npy_array()
        js = js[js > i]
        ii.extend([i]*len(js))
        jj.extend(js.tolist())
    return np.asarray(ii, dtype=int), np.asarray(jj, dtype=int)


def operators(pa):
    i, j = pairs(pa)
    x = positions(pa)
    g = HyperbolicSpline().gradient(x[i]-x[j], (pa.h[i]+pa.h[j])/2)
    dv = velocities(pa)[j]-velocities(pa)[i]
    outer = dv[:, :, None]*g[:, None, :]
    grad = np.zeros((len(x), 3, 3))
    np.add.at(grad, i, pa.m[j, None, None]/pa.rho[j, None, None]*outer)
    np.add.at(grad, j, pa.m[i, None, None]/pa.rho[i, None, None]*outer)
    # Exactly the Table 1 rho_i volume_j continuity form.
    drho = -pa.rho*np.trace(grad, axis1=1, axis2=2)
    return i, j, g, grad, drho


def corrected_operators(pa, condition_limit=1e6):
    """Explicit APPROXIMATE first-order-complete, work-conjugate variant.

    The original Table 1 operators above remain available for audit. For
    r_ij=x_j-x_i and radial grad_i W_ij, M_i=sum(V_j r_ij tensor grad_i W_ij)
    gives L_raw=A M_i in an affine velocity field. L=L_raw inv(M_i) therefore
    reproduces all three-dimensional affine fields, including rigid rotation.
    The returned B must also be passed to internal_force; correcting only L
    would destroy work conjugacy. The returned g stays radial for conduction.
    Singular/free isolated neighborhoods are unsupported, and fail closed.
    """
    i, j, g, raw_gradient, _ = operators(pa)
    x = positions(pa)
    moment_pair = (x[j]-x[i])[:, :, None]*g[:, None, :]
    moment = np.zeros_like(raw_gradient)
    np.add.at(moment, i, pa.m[j, None, None]/pa.rho[j, None, None]*moment_pair)
    np.add.at(moment, j, pa.m[i, None, None]/pa.rho[i, None, None]*moment_pair)
    if not np.isfinite(moment).all():
        raise ValueError('non-finite correction moment matrix')
    eigenvalues = np.linalg.eigvalsh(moment)
    smallest, largest = eigenvalues[:, 0], eigenvalues[:, -1]
    bad = (smallest <= 0) | (largest > condition_limit*smallest)
    if np.any(bad):
        indices = np.flatnonzero(bad)[:8].tolist()
        raise ValueError(f'degenerate/ill-conditioned 3D SPH neighborhoods at indices {indices}; correction condition limit={condition_limit:g}')
    B = np.linalg.inv(moment)
    gradient = raw_gradient@B
    drho = -pa.rho*np.trace(gradient, axis1=1, axis2=2)
    return i, j, g, gradient, drho, B


def internal_force(pa, stress, i, j, g, correction=None):
    if correction is None:
        # Original paper Table 1 discretization, retained without alteration.
        pair = pa.m[i, None]*pa.m[j, None]*np.einsum(
            'nij,nj->ni', stress[i]/pa.rho[i, None, None]**2+stress[j]/pa.rho[j, None, None]**2, g)
    else:
        # F_ij=V_i V_j (sigma_i B_i^T + sigma_j B_j^T) grad_i W_ij.
        # Sum_i v_i.F_i + sum_i V_i sigma_i:L_i = 0 to roundoff, even
        # at unequal density. Symmetric stresses + complete M also conserve
        # total angular momentum. This is not the raw Table 1 force formula.
        gi = np.einsum('nji,nj->ni', correction[i], g)
        gj = np.einsum('nji,nj->ni', correction[j], g)
        pair = (pa.m[i]/pa.rho[i]*pa.m[j]/pa.rho[j])[:, None]*(
            np.einsum('nij,nj->ni', stress[i], gi)+np.einsum('nij,nj->ni', stress[j], gj))
    force = np.zeros((pa.get_number_of_particles(), 3))
    np.add.at(force, i, pair)
    np.add.at(force, j, -pair)
    return force


def bulk_modulus(shear, poisson):
    return 2*shear*(1+poisson)/(3*(1-2*poisson))


def pressure(rho, rho0, bulk):
    if np.any(rho <= 0):
        raise ValueError('nonpositive density')
    return bulk*(rho/rho0-1)


def volumetric_energy(rho, rho0, bulk):
    ratio = rho/rho0
    return bulk/rho0*(np.log(ratio)+1/ratio-1)

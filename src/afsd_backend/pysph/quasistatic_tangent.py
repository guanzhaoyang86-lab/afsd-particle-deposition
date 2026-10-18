"""Frozen-geometry force derivative for an independently initialized QS model.

This is a modified-Newton search matrix, not a physical update or full tangent.
There is no mass matrix, inertial regularizer, pseudoinverse or symmetrization.
The physical FKS/J2 update and its disclosed low-rate floor are unchanged.
"""
import numpy as np

from .material import rotation, update
from .material_response import (NonDifferentiableMaterial, fks_partials,
                                radial_direction, trial_direction)

CLASSIFICATION = 'FROZEN_BE_QUASISTATIC_MODIFIED_NEWTON_NOT_FULL_TANGENT'


def _real(value, shape, name):
    a = np.asarray(value)
    if a.dtype.kind not in 'iuf' or a.shape != shape or not np.isfinite(a).all():
        raise ValueError('invalid finite real '+name)
    return np.array(a, dtype=float, copy=True)


def gradient_map(volume, i, j, g, correction):
    """Dense-storage D [1/m], with the actual sparse neighborhood entries.

    Accept the current evaluation's actual NNPS pairs and correction rather
    than rebuilding its geometry. The same D is conjugate to -D.T @ (V*sigma).
    """
    raw = np.asarray(volume)
    if raw.ndim != 1 or not len(raw):
        raise ValueError('nonempty volume vector required')
    n = len(raw)
    V = _real(volume, (n,), 'volume')
    if np.any(V <= 0):
        raise ValueError('positive volumes required')
    i, j = np.asarray(i), np.asarray(j)
    if (i.ndim != 1 or j.shape != i.shape or i.dtype.kind not in 'iu'
            or j.dtype.kind not in 'iu' or np.any(i < 0) or np.any(j < 0)
            or np.any(i >= n) or np.any(j >= n) or np.any(i == j)):
        raise ValueError('valid pair indices required')
    g = _real(g, (len(i), 3), 'pair gradients')
    B = _real(correction, (n, 3, 3), 'correction')
    gi = np.einsum('nji,nj->ni', B[i], g)
    gj = np.einsum('nji,nj->ni', B[j], g)
    D = np.zeros((9*n, 3*n))
    for component in range(3):
        for axis in range(3):
            a, b = V[j]*gi[:, axis], V[i]*gj[:, axis]
            np.add.at(D, (9*i+3*component+axis, 3*j+component), a)
            np.add.at(D, (9*i+3*component+axis, 3*i+component), -a)
            np.add.at(D, (9*j+3*component+axis, 3*j+component), b)
            np.add.at(D, (9*j+3*component+axis, 3*i+component), -b)
    return D


def material_tangent(stress0, ep0, gradient, temperature, dt, shear, law,
                     material=None):
    """Return dS1/dL [Pa s] and dS1/dT [Pa/K] for this same-history trial.

    ``material`` optionally supplies the exact five-tuple from material.update.
    Strictly elastic zero stress and zero strain rate are supported: the
    composite stress map is smooth even though an intermediate norm is not.
    At zero rate the actual FKS floor makes yield locally rate-independent;
    no artificial nonzero rate is substituted. Plastic yield/floor switching
    points remain unsupported classical derivatives and fail explicitly.
    """
    raw = np.asarray(stress0)
    if raw.ndim != 3 or raw.shape[1:] != (3, 3) or len(raw) == 0:
        raise ValueError('nonempty particle stress tensors required')
    n = len(raw)
    S0 = _real(stress0, (n, 3, 3), 'parent stress')
    L = _real(gradient, (n, 3, 3), 'gradient')
    ep = _real(ep0, (n,), 'parent plastic history')
    T = _real(np.broadcast_to(temperature, (n,)), (n,), 'temperature')
    if not np.isfinite([dt, shear]).all() or dt <= 0 or shear <= 0:
        raise ValueError('positive finite dt and shear required')
    base = trial_direction(S0, L, np.zeros_like(L), dt, shear)
    trial, rate = base['trial_stress_Pa'], base['rate_per_s']
    qtr = np.sqrt(1.5*np.sum(trial*trial, axis=(1, 2)))
    initial_yield = law.yield_stress(ep, rate, T)
    if np.any(qtr == initial_yield):
        raise NonDifferentiableMaterial('trial at yield switching surface')
    active = qtr > initial_yield
    result = update(S0, ep, L, T, dt, shear, law) if material is None else material
    if not isinstance(result, (tuple, list)) or len(result) != 5:
        raise ValueError('exact material.update five-tuple required')
    S1 = _real(result[0], (n, 3, 3), 'returned stress')
    ep1 = _real(result[1], (n,), 'returned plastic history')
    saved_rate = _real(result[4], (n,), 'returned rate')
    dep = ep1-ep
    if np.any(dep < 0) or not np.array_equal(active, dep > 0):
        raise ValueError('inconsistent returned plastic branch')
    if not np.allclose(saved_rate, rate, rtol=1e-12, atol=1e-12):
        raise ValueError('inconsistent returned material rate')
    factor = np.divide(qtr-3*shear*dep, qtr, out=np.ones(n), where=qtr > 0)
    if not np.allclose(S1, trial*factor[:, None, None], rtol=1e-12, atol=1e-6):
        raise ValueError('inconsistent returned material stress')
    partials = fks_partials(law, ep1[active], rate[active], T[active]) if active.any() else None
    return_residual = np.zeros(n)
    if active.any():
        return_residual[active] = (qtr[active]-3*shear*dep[active]
                                   -partials['yield_stress_Pa'])
        scale = np.maximum.reduce((qtr[active], 3*shear*np.abs(ep1[active]),
                                   partials['yield_stress_Pa']))
        # Covers the unchanged 48-bisection error and ep1-ep0 cancellation.
        if np.any(np.abs(return_residual[active]) > 256*np.finfo(float).eps*scale):
            raise ValueError('returned history does not close same-temperature FKS law')
    zero = rate == 0
    Rzero = rotation(.25*dt*(L[zero]-L[zero].swapaxes(1, 2)))
    C = np.zeros((n, 9, 9))
    for a in range(3):
        for b in range(3):
            basis = np.zeros_like(L)
            basis[:, a, b] = 1.
            dst, drate = np.zeros_like(L), np.zeros(n)
            if (~zero).any():
                tr = trial_direction(S0[~zero], L[~zero], basis[~zero], dt, shear)
                dst[~zero], drate[~zero] = tr['delta_trial_stress_Pa'], tr['delta_rate_per_s']
            if zero.any():
                direction = basis[zero]
                skew = .5*(direction-direction.swapaxes(1, 2))
                dev = (.5*(direction+direction.swapaxes(1, 2))
                       -np.trace(direction, axis1=1, axis2=2)[:, None, None]*np.eye(3)/3)
                spin_part = trial_direction(S0[zero], L[zero], skew, dt, shear)
                dst[zero] = (spin_part['delta_trial_stress_Pa']
                             +2*shear*dt*(Rzero@dev@Rzero.swapaxes(1, 2)))
            derivative = dst.copy()
            if active.any():
                response = radial_direction(trial[active], dep[active], dst[active],
                    drate[active], np.zeros(active.sum()), shear, partials,
                    np.ones(active.sum(), dtype=bool))
                derivative[active] = response['delta_stress_Pa']
            C[:, :, 3*a+b] = derivative.reshape(n, 9)
    ST = np.zeros((n, 9))
    if active.any():
        response = radial_direction(trial[active], dep[active], np.zeros_like(trial[active]),
            np.zeros(active.sum()), np.ones(active.sum()), shear, partials,
            np.ones(active.sum(), dtype=bool))
        ST[active] = response['delta_stress_Pa'].reshape(-1, 9)
    return dict(tangent_Pa_s=C, temperature_tangent_Pa_K=ST, stress_Pa=S1,
                ep=ep1, rate_per_s=rate, plastic_active=active,
                base_return_residual_Pa=return_residual,
                classification=CLASSIFICATION, physical_update=False,
                full_tangent=False)


def assemble_frozen_force_tangent(D, volume, C, rho0, rho_reference, bulk,
                                  dt, displacement):
    """Frozen d f_int / d u [N/m] for endpoint BE stress and continuity.

    L=D u/dt; rho1=rho0/(1+tr(Du)); sigma=S1-p(rho1)I.
    Excludes derivatives of D,V, thermal feedback and contact. The full force
    residual must be recomputed for every line-search candidate. No inertia.
    """
    raw = np.asarray(volume)
    if raw.ndim != 1 or not len(raw):
        raise ValueError('nonempty volume vector required')
    n = len(raw)
    V = _real(volume, (n,), 'volume')
    C = _real(C, (n, 9, 9), 'material tangent')
    oldrho = _real(rho0, (n,), 'parent density')
    ref = _real(np.broadcast_to(rho_reference, (n,)), (n,), 'reference density')
    u = _real(displacement, (n, 3), 'displacement')
    if (np.any(V <= 0) or np.any(oldrho <= 0) or np.any(ref <= 0)
            or not np.isfinite([bulk, dt]).all() or bulk <= 0 or dt <= 0):
        raise ValueError('positive finite volumes, densities, bulk and dt required')
    D = _real(D, (9*n, 3*n), 'gradient')
    E = D[9*np.arange(n)]+D[9*np.arange(n)+4]+D[9*np.arange(n)+8]
    denominator = 1+E@u.ravel()
    if np.any(denominator <= 0):
        raise ValueError('BE continuity outside positive density domain')
    pressure = -V*bulk/ref*oldrho/denominator**2
    J = np.zeros((3*n, 3*n))
    # SciPy is not an installed pinned dependency. Work only on each node's
    # actual nonzero columns, avoiding the full dense 3N-by9N-by3N product.
    for k in range(n):
        local = D[9*k:9*k+9]
        columns = np.flatnonzero(np.any(local != 0, axis=0))
        small = local[:, columns]
        divergence = E[k, columns]
        block = (-V[k]/dt*(small.T@C[k]@small)
                 +pressure[k]*np.outer(divergence, divergence))
        J[np.ix_(columns, columns)] += block
    if not np.isfinite(J).all():
        raise ValueError('nonfinite frozen force tangent')
    return J

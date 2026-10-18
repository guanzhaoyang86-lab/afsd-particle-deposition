"""Detached constant-K HHT component; no contact/material/physical acceptance.

The independently derived energy identity and preregistration are in
docs/PROCESS_HHT_LINEAR_CONTRACT.md. Numerical energy loss is never heat.
"""
from dataclasses import dataclass

import numpy as np


CLASSIFICATION = 'PROCESS_HHT_EVENT_THERMAL_APPROXIMATE_SMOOTH_LINEAR_COMPONENT'


def _real(value, name, shape=None):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or (shape is not None and raw.shape != shape):
        raise ValueError('invalid real shape/type: '+name)
    out = np.array(raw, dtype=float, copy=True)
    if not np.isfinite(out).all():
        raise ValueError('nonfinite '+name)
    return out


def _parameters(dt, weight):
    h = float(_real(dt, 'dt', ()))
    w = float(_real(weight, 'weight', ()))
    if h <= 0 or not 2./3. <= w <= 1.:
        raise ValueError('positive dt and future-force weight in [2/3,1] required')
    d = 1.-w
    return h, w, .5+d, .25*(1.+d)**2


def _system(mass, stiffness):
    m = _real(mass, 'mass')
    if m.ndim != 1 or not len(m) or np.any(m <= 0):
        raise ValueError('positive nonempty diagonal mass vector required')
    k = _real(stiffness, 'stiffness', (len(m), len(m)))
    scale = max(float(np.linalg.norm(k, 2)), np.finfo(float).tiny)
    if np.linalg.norm(k-k.T, 2) > 1e-12*scale:
        raise ValueError('nonsymmetric stiffness')
    minimum = float(np.linalg.eigvalsh(k)[0])
    if minimum < -1e-12*scale:
        raise ValueError('negative elastic stiffness; no clipping')
    return m, k, minimum


def kinematic_prediction(x0, v0, acc0, dt, weight=.9):
    """Return copies of Newmark predictors; acceleration remains explicit."""
    h, _, gamma, beta = _parameters(dt, weight)
    x = _real(x0, 'x0')
    v = _real(v0, 'v0', x.shape)
    z = _real(acc0, 'acc0', x.shape)
    return dict(x_predictor=x+h*v+h*h*(.5-beta)*z,
                v_predictor=v+h*(1.-gamma)*z, beta=beta, gamma=gamma)


def initialize_acceleration(mass, stiffness, x, force):
    """Smooth, unconstrained start only; invalid as a contact initializer."""
    m, k, _ = _system(mass, stiffness)
    u = _real(x, 'x', m.shape)
    f = _real(force, 'force', m.shape)
    return (f-k@u)/m


@dataclass(frozen=True)
class LinearStep:
    x: np.ndarray
    v: np.ndarray
    acc: np.ndarray
    displacement_increment: np.ndarray
    load_increment: np.ndarray
    ledger: dict


def step(mass, stiffness, x0, v0, acc0, f0, f1, dt, weight=.9,
         previous_displacement_increment=None, previous_load_increment=None):
    """One smooth HHT step, with saved history and raw signed energy terms.

    None histories mean zero virtual increments for initialization, never an
    inferred restart. A nonzero force-history defect is exposed in the ledger.
    The caller owns acceptance, persistent IDs, contact, material and thermal
    transactions. This function does not perform those operations.
    """
    h, w, gamma, beta = _parameters(dt, weight)
    m, k, minimum = _system(mass, stiffness)
    x = _real(x0, 'x0', m.shape)
    v = _real(v0, 'v0', m.shape)
    z = _real(acc0, 'acc0', m.shape)
    fa = _real(f0, 'f0', m.shape)
    fb = _real(f1, 'f1', m.shape)
    b = (np.zeros_like(x) if previous_displacement_increment is None else
         _real(previous_displacement_increment, 'previous displacement increment', m.shape))
    g = (np.zeros_like(x) if previous_load_increment is None else
         _real(previous_load_increment, 'previous load increment', m.shape))
    prediction = kinematic_prediction(x, v, z, h, w)
    matrix = np.diag(m)+w*beta*h*h*k
    rhs = w*(fb-k@prediction['x_predictor'])+(1.-w)*(fa-k@x)
    z1 = np.linalg.solve(matrix, rhs)
    x1 = prediction['x_predictor']+beta*h*h*z1
    v1 = prediction['v_predictor']+gamma*h*z1
    if not all(np.isfinite(a).all() for a in (matrix, rhs, z1, x1, v1)):
        raise ValueError('nonfinite HHT solution; no clipping')
    dx, dz, df = x1-x, z1-z, fb-fa
    d = 1.-w
    c = d*d/4.
    q = d*(.5-d)
    history_defect = m*z-(fa-k@x)+d*(g-k@b)
    physical0 = .5*float(v@(m*v)+x@k@x)
    physical1 = .5*float(v1@(m*v1)+x1@k@x1)
    acc_energy0 = .5*h*h*c*float(z@(m*z))
    acc_energy1 = .5*h*h*c*float(z1@(m*z1))
    increment_energy0 = .5*q*float(b@k@b)
    increment_energy1 = .5*q*float(dx@k@dx)
    algorithmic0 = physical0+acc_energy0+increment_energy0
    algorithmic1 = physical1+acc_energy1+increment_energy1
    loss_acc = h*h*d*c*float(dz@(m*dz))
    loss_increment = .5*q*float((dx-b)@k@(dx-b))
    work = .5*float(dx@(fa+fb))
    quadrature_difference = h*float(df@(v1-v))/12.
    load_work = q*float(dx@(df-g))
    history_work = (.5-d)*float(dx@history_defect)
    residual = algorithmic1-algorithmic0+loss_acc+loss_increment-work-load_work-history_work
    equilibrium = m*z1-w*(fb-k@x1)-d*(fa-k@x)
    ledger = dict(classification=CLASSIFICATION, accepted=False,
        physical_energy0_J=physical0, physical_energy1_J=physical1,
        acceleration_energy0_J=acc_energy0, acceleration_energy1_J=acc_energy1,
        increment_energy0_J=increment_energy0, increment_energy1_J=increment_energy1,
        algorithmic_energy0_J=algorithmic0, algorithmic_energy1_J=algorithmic1,
        algorithmic_acceleration_loss_J=loss_acc,
        algorithmic_increment_loss_J=loss_increment,
        external_trapezoidal_work_J=work,
        exact_linear_load_work_J=work+quadrature_difference,
        load_quadrature_difference_J=quadrature_difference,
        work_path='CUBIC_HERMITE_DISPLACEMENT_LINEAR_FORCE',
        algorithmic_load_work_J=load_work,
        acceleration_history_work_J=history_work,
        acceleration_history_defect_N=history_defect.copy(),
        acceleration_history_defect_max_N=float(np.max(np.abs(history_defect))),
        identity_residual_J=residual,
        equilibrium_residual_max_N=float(np.max(np.abs(equilibrium))),
        stiffness_min_eigenvalue_N_m=minimum,
        stiffness_psd_nonnegative_to_computed_precision=bool(minimum >= 0),
        dt_s=h, weight=w, gamma=gamma, beta=beta,
        algorithmic_terms_are_material_heat=False,
        physical_energy_monotonicity_claimed=False)
    return LinearStep(x1.copy(), v1.copy(), z1.copy(), dx.copy(), df.copy(), ledger)

"""Explicit auxiliary derivative domain for the unchanged raw-stress update.

Q S Q.T and its product derivative are defined for any finite symmetric S;
zero trace is not needed. The original physical material.update already uses
this formula without projecting S. This search-only provider retains the
original symmetry tolerance, rejects nonfinite tensors, and copies the input.
It removes only the initial trace precondition of trial_direction. The active
radial-return derivative retains every original deviatoric/branch check.

Observed trigger: saved particle 1216 unloaded from 64700.9 to 2954.2 Pa while
its inherited trace remained about -3.6e-10 Pa, exceeding the old local trace
limit of 2.52e-10 Pa. No stress, residual, acceptance gate or material update is
changed here. A full physical comparison remains required by the new bridge.
"""
import ast
import copy
import inspect
import textwrap
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:
    sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import material_response as original_response
from afsd_backend.pysph import quasistatic_tangent as original_tangent
from afsd_backend.pysph.material import rotation, update
from afsd_backend.pysph.material_response import (
    NonDifferentiableMaterial, fks_partials, radial_direction, _tensor)

# The two source modules use different _real validators; preserve each binding.
_real = original_response._real
_tangent_real = original_tangent._real
CLASSIFICATION = 'EXPLICIT_FINITE_SYMMETRIC_STRESS_SEARCH_DERIVATIVE_20261009'
SEARCH_PROTOCOL = CLASSIFICATION


def _symmetric_tensor(value, name):
    a = _tensor(value, name)
    scale = np.maximum(1., np.max(np.abs(a), axis=(-2, -1)))
    tolerance = 128*np.finfo(float).eps*scale
    if np.any(np.max(np.abs(a-a.swapaxes(-2, -1)), axis=(-2, -1)) > tolerance):
        raise ValueError(name+' must be symmetric; not projected')
    return a


def trial_direction(stress, gradient, delta_gradient, dt, shear):
    """Analytic Frechet derivative of the unchanged corotational trial stress.

    Rotations use the same Rodrigues baseline as material.update. Series only
    avoid cancellation in derivatives near zero spin; no physical field is
    projected or clipped. The returned rate derivative is undefined at zero
    deviatoric rate for nonzero deviatoric directions and explicitly rejects.
    """
    s = _symmetric_tensor(stress,'stress')
    L = _tensor(gradient,'gradient'); dL = _tensor(delta_gradient,'delta_gradient')
    h = _real(dt,'dt'); G = _real(shear,'shear')
    if h.shape or G.shape or h <= 0 or G <= 0:
        raise ValueError('positive scalar dt/shear required')
    shape = np.broadcast_shapes(s.shape[:-2],L.shape[:-2],dL.shape[:-2])
    s,L,dL = [np.broadcast_to(a,shape+(3,3)) for a in (s,L,dL)]
    eye = np.eye(3)
    def dev(a):
        return .5*(a+a.swapaxes(-2,-1))-np.trace(a,axis1=-2,axis2=-1)[...,None,None]*eye/3
    D,dD = dev(L),dev(dL)
    A = .25*float(h)*(L-L.swapaxes(-2,-1))
    dA = .25*float(h)*(dL-dL.swapaxes(-2,-1))
    squared = .5*np.sum(A*A,axis=(-2,-1)); angle = np.sqrt(squared)
    d_squared = np.sum(A*dA,axis=(-2,-1))
    small = squared < 1e-6
    # Evaluate the non-small denominator at 1 for series nodes, avoiding 0/0.
    angle_safe = np.where(small,1.,angle)
    da_large = (angle_safe*np.cos(angle_safe)-np.sin(angle_safe))/(2*angle_safe**3)
    db_large = (angle_safe*np.sin(angle_safe)-4*np.sin(angle_safe/2)**2)/(2*angle_safe**4)
    da = np.where(small,-1/6+squared/60-squared**2/1680+squared**3/90720,da_large)
    db = np.where(small,-1/24+squared/360-squared**2/13440+squared**3/907200,db_large)
    a = np.sinc(angle/np.pi); b = .5*np.sinc(angle/(2*np.pi))**2
    R = rotation(A)
    dR = ((da*d_squared)[...,None,None]*A+a[...,None,None]*dA
          +(db*d_squared)[...,None,None]*(A@A)+b[...,None,None]*(dA@A+A@dA))
    Q = R@R; dQ = dR@R+R@dR
    Rt=R.swapaxes(-2,-1); dRt=dR.swapaxes(-2,-1)
    Qt=Q.swapaxes(-2,-1); dQt=dQ.swapaxes(-2,-1)
    trial = Q@s@Qt+2*float(G*h)*(R@D@Rt)
    delta = dQ@s@Qt+Q@s@dQt+2*float(G*h)*(dR@D@Rt+R@dD@Rt+R@D@dRt)
    rate = np.sqrt(2/3*np.sum(D*D,axis=(-2,-1)))
    zero = rate == 0
    if np.any(zero & np.any(dD != 0,axis=(-2,-1))):
        raise NonDifferentiableMaterial('rate norm at zero with nonzero deviatoric direction')
    dr = np.divide((2/3)*np.sum(D*dD,axis=(-2,-1)),rate,
                   out=np.zeros(shape),where=~zero)
    if not all(np.isfinite(a).all() for a in (trial,delta,rate,dr)):
        raise ValueError('nonfinite corotational response')
    return dict(trial_stress_Pa=trial,delta_trial_stress_Pa=delta,
                rate_per_s=rate,delta_rate_per_s=dr,classification=CLASSIFICATION)


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
    S0 = _tangent_real(stress0, (n, 3, 3), 'parent stress')
    L = _tangent_real(gradient, (n, 3, 3), 'gradient')
    ep = _tangent_real(ep0, (n,), 'parent plastic history')
    T = _tangent_real(np.broadcast_to(temperature, (n,)), (n,), 'temperature')
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
    S1 = _tangent_real(result[0], (n, 3, 3), 'returned stress')
    ep1 = _tangent_real(result[1], (n,), 'returned plastic history')
    saved_rate = _tangent_real(result[4], (n,), 'returned rate')
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

def static_adaptation_audit(actual_tree=None):
    """Validate exact static copies; optional module AST supports negative tests."""
    actual = (ast.parse(Path(__file__).read_text(encoding='utf-8'))
              if actual_tree is None else copy.deepcopy(actual_tree))
    def function(tree, name):
        return next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    symmetric = ast.parse(textwrap.dedent(inspect.getsource(original_response._deviatoric))).body[0]
    symmetric.name = '_symmetric_tensor'
    guard = next(n for n in symmetric.body if isinstance(n, ast.If))
    if not isinstance(guard.test, ast.BoolOp) or not isinstance(guard.test.op, ast.Or):
        raise ValueError('unexpected original tensor domain check')
    guard.test = guard.test.values[0]
    for node in ast.walk(symmetric):
        if isinstance(node, ast.Constant) and node.value == ' must be symmetric deviatoric; not projected':
            node.value = ' must be symmetric; not projected'
    if ast.dump(symmetric, include_attributes=False) != ast.dump(function(actual, '_symmetric_tensor'), include_attributes=False):
        raise ValueError('finite or symmetry check changed')
    for original, name in ((original_response.trial_direction, 'trial_direction'),
                           (original_tangent.material_tangent, 'material_tangent')):
        expected = ast.parse(textwrap.dedent(inspect.getsource(original)))
        class Provider(ast.NodeTransformer):
            def visit_Name(self, node):
                if name == 'trial_direction' and node.id == '_deviatoric':
                    node.id = '_symmetric_tensor'
                if name == 'material_tangent' and node.id == '_real':
                    node.id = '_tangent_real'
                return node
        wanted = Provider().visit(expected).body[0]
        if ast.dump(wanted, include_attributes=False) != ast.dump(function(actual, name), include_attributes=False):
            raise ValueError('search derivative body outside explicit provider adaptation: '+name)
    if (_real is not original_response._real or _tensor is not original_response._tensor
            or _tangent_real is not original_tangent._real
            or rotation is not original_response.rotation or update is not original_tangent.update
            or radial_direction is not original_response.radial_direction
            or fks_partials is not original_response.fks_partials):
        raise ValueError('original numerical provider changed')
    return dict(status='PASS_EXPLICIT_SEARCH_DERIVATIVE_DOMAIN_EXTENSION',
                protocol=SEARCH_PROTOCOL, finite_and_symmetry_checks_preserved=True,
                parent_trace_precondition_removed_from_search_only=True,
                active_radial_deviatoric_guard_unchanged=True,
                physical_update_unchanged=True, stress_projected=False,
                physical_gates_changed=False, original_globals_modified=False)

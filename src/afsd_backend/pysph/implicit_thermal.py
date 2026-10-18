"""Conservative BE heat graph; named component, NOT integrated AFSD physics.

Same radial SPH pair coefficients as thermal.py. Temperatures at the end of
the interval set conductive fluxes. Geometry/coefficients stay frozen inside
this component; mechanical source energies must come from actual intervals.
"""
import math
import numpy as np
from .coupled_thermal import _vector, _scalar

CLASSIFICATION = 'ENGINEERING_APPROXIMATE_IMPLICIT_THERMAL_COMPONENT_NOT_INTEGRATED'


def pair_conductance(mass, rho, conductivity, x, i, j, gradient):
    """Undirected conductance [W/K], never a corrected mechanical gradient."""
    m = np.asarray(mass, float)
    if m.ndim != 1 or not len(m) or not np.isfinite(m).all() or np.any(m <= 0):
        raise ValueError('positive finite mass vector required')
    n = len(m)
    density = _vector(rho, n, 'rho')
    k = _vector(conductivity, n, 'conductivity', broadcast=True)
    if np.any(density <= 0) or np.any(k <= 0):
        raise ValueError('positive density and conductivity required')
    points = np.asarray(x, float); first = np.asarray(i); second = np.asarray(j)
    g = np.asarray(gradient, float)
    if points.shape != (n, 3) or not np.isfinite(points).all():
        raise ValueError('finite Nx3 positions required')
    if (first.ndim != 1 or second.shape != first.shape or
            first.dtype.kind not in 'iu' or second.dtype.kind not in 'iu'):
        raise ValueError('integer matching edge indices required')
    if g.shape != (len(first), 3) or not np.isfinite(g).all():
        raise ValueError('finite edge gradient required')
    if np.any(first < 0) or np.any(second < 0) or np.any(first >= n) or np.any(second >= n):
        raise ValueError('edge index outside particle inventory')
    if np.any(first == second):
        raise ValueError('self conduction edge')
    edges = np.sort(np.column_stack((first, second)), axis=1)
    if len(np.unique(edges, axis=0)) != len(first):
        raise ValueError('duplicate undirected conduction edge')
    rij = points[first]-points[second]
    radius2 = np.sum(rij*rij, axis=1)
    if np.any(radius2 == 0):
        raise ValueError('coincident pair')
    if np.any(np.linalg.norm(np.cross(rij, g), axis=1) >
              128*np.finfo(float).eps*np.linalg.norm(rij, axis=1)*np.linalg.norm(g, axis=1)):
        raise ValueError('raw radial gradient required')
    coefficient = (m[first]*m[second]/(density[first]*density[second]) *
                   (4*k[first]*k[second]/(k[first]+k[second])) *
                   (-np.sum(rij*g, axis=1))/radius2)
    if not np.isfinite(coefficient).all() or np.any(coefficient < 0):
        raise ValueError('negative/nonfinite pair conductance; not clipped')
    return coefficient


def advance_implicit_thermal(mass, rho, T, conductivity, x, i, j, gradient, dt, cp,
                             plastic_work_J, friction_raw_J, friction_rounding_bound_J,
                             wall_row_W_K, wall_weighted_temperature_W,
                             chi, partition, room, melt):
    """Pure frozen-coefficient heat solve and unmodified signed source ledger.

wall_weighted_temperature_W = sum(h_wall * T_wall), per particle. Supplying
an old contact-power value here is WRONG; caller must retain reservoir data.
Friction certificate belongs to the mechanical caller, not to this solver.
No clipping, artificial residual heat, mass scaling or mechanical return heat.
"""
    t0 = np.asarray(T, float)
    if t0.ndim != 1 or not len(t0) or not np.isfinite(t0).all():
        raise ValueError('nonempty finite temperature vector required')
    n = len(t0); m = _vector(mass, n, 'mass'); cpv = _vector(cp, n, 'cp', broadcast=True)
    qplastic = _vector(plastic_work_J, n, 'plastic_work_J')
    qfriction = _vector(friction_raw_J, n, 'friction_raw_J')
    qround = _vector(friction_rounding_bound_J, n, 'friction_rounding_bound_J')
    wall = _vector(wall_row_W_K, n, 'wall_row_W_K')
    bath = _vector(wall_weighted_temperature_W, n, 'wall_weighted_temperature_W')
    h = _scalar(dt, 'dt'); fraction = _scalar(chi, 'chi'); split = _scalar(partition, 'partition')
    lower = _scalar(room, 'room'); upper = _scalar(melt, 'melt')
    if h <= 0 or np.any(m <= 0) or np.any(cpv <= 0) or lower <= 0 or upper <= lower:
        raise ValueError('positive dt/capacity and ordered positive thermal domain required')
    if not 0 <= fraction <= 1 or not 0 <= split <= 1:
        raise ValueError('heat fractions outside [0,1]')
    if np.any(t0 < lower) or np.any(t0 >= upper):
        raise ValueError('initial temperature outside solid domain')
    if np.any(qplastic < 0) or np.any(qround < 0) or np.any(qfriction < -qround):
        raise ValueError('negative plastic/friction outside mechanical roundoff certificate; not clipped')
    if np.any(wall < 0) or np.any(bath < 0) or np.any((wall == 0) & (bath != 0)):
        raise ValueError('invalid nonnegative reservoir conductance/weighted temperature')
    coefficient = pair_conductance(m, rho, conductivity, x, i, j, gradient)
    first = np.asarray(i); second = np.asarray(j)
    laplacian = np.zeros((n, n))
    np.add.at(laplacian, (first, first), coefficient)
    np.add.at(laplacian, (second, second), coefficient)
    np.add.at(laplacian, (first, second), -coefficient)
    np.add.at(laplacian, (second, first), -coefficient)
    capacity = m*cpv
    operator = laplacian+np.diag(wall)
    matrix = np.diag(capacity)+h*operator
    rhs = fraction*qplastic+split*qfriction+h*(bath-operator@t0)
    if not np.isfinite(matrix).all() or not np.isfinite(rhs).all():
        raise ValueError('nonfinite heat system')
    increment = np.linalg.solve(matrix, rhs)
    t1 = t0+increment
    if not np.isfinite(t1).all() or np.any(t1 < lower) or np.any(t1 >= upper):
        raise ValueError('thermal solid-domain violation; temperature was NOT clipped')
    pair_heat = h*coefficient*(t1[second]-t1[first])
    internal = np.zeros(n)
    np.add.at(internal, first, pair_heat); np.add.at(internal, second, -pair_heat)
    particle = dict(plastic_work_J=qplastic.copy(), plastic_heat_J=fraction*qplastic,
        plastic_unthermalized_J=(1-fraction)*qplastic, friction_raw_J=qfriction.copy(),
        friction_material_heat_J=split*qfriction, friction_wall_heat_J=(1-split)*qfriction,
        contact_in_J=h*(bath-wall*t1), internal_conduction_J=internal)
    total = particle['plastic_heat_J']+particle['friction_material_heat_J']+particle['contact_in_J']+internal
    measured = capacity*(t1-t0)
    particle.update(total_heat_input_J=total, measured_thermal_increment_J=measured,
                    mc_delta_T_residual_J=measured-total)
    ledger = {key: math.fsum(float(v) for v in value) for key, value in particle.items()}
    ledger['wall_heat_out_J'] = ledger['friction_wall_heat_J']-ledger['contact_in_J']
    balance_limit = 1e-12+1e-10*float(np.abs(total).sum())
    if float(np.max(np.abs(measured-total))) > balance_limit:
        raise ValueError('implicit heat algebraic balance failed; no residual heat added')
    row = laplacian.diagonal()+wall
    old_bound = .25*float(np.min(np.divide(capacity, row, out=np.full(n, np.inf), where=row > 0)))
    return dict(classification=CLASSIFICATION, status='COMPONENT_SOLVED_NOT_INTEGRATED',
        temperature=t1, pair_conductance_W_K=coefficient, pair_heat_into_i_J=pair_heat,
        particle_energy_J=particle, ledger=ledger, linear_system_residual_J=matrix@increment-rhs,
        former_explicit_dt_bound_s=old_bound,
        temporal_accuracy='FIRST_ORDER_BE_FROZEN_COEFFICIENTS_NOT_MACROSTEP_CERTIFICATION',
        production_authorized=False)

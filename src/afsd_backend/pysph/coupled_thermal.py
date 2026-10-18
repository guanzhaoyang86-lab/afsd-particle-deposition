"""Metered explicit thermal companion, reusing the unchanged thermal operators.

Mechanical corrected gradients MUST NOT be substituted for the radial pair
gradient supplied here. Friction work is already contact-certified and remains
signed; neither return dissipation nor a mechanical residual is a heat input.
This pure function does not mutate a ParticleArray or any supplied array.
"""
import math
import numpy as np

from .thermal import conduction, temperature_step


def _vector(value, count, name, broadcast=False):
    result = np.asarray(value, dtype=float)
    if broadcast:
        try:
            result = np.broadcast_to(result, (count,))
        except ValueError as error:
            raise ValueError(name+' must be scalar or an N vector') from error
    if result.shape != (count,) or not np.isfinite(result).all():
        raise ValueError(name+' must be a finite N vector')
    return result


def _scalar(value, name):
    result = np.asarray(value, dtype=float)
    if result.shape != () or not np.isfinite(result):
        raise ValueError(name+' must be a finite scalar')
    return float(result)


def advance_thermal(mass, rho, T, conductivity, x, i, j, g, dt, cp,
                    plastic_work_J, friction_raw_J, contact_power_W,
                    wall_row_W_K, chi, partition, room, melt):
    """Return temperature, pair-conduction power and separately metered sources.

All units are SI. N-vector sources are energy during the actual mechanical
interval (plastic/friction [J]) and signed wall-to-material power [W]. The
normal penalty spring and numerical return work are NOT thermal sources.

``pair_power_W`` is the N-vector of net internal pair-conduction power, not
an edge-indexed array. ``thermal_dt_bound_s`` is the enforced quarter-capacity
bound, or infinity if all conductance rows vanish. cp and conductivity accept
a scalar or N vector; all other particle-valued inputs require N vectors.
The caller owns the contact algorithm's raw-friction acceptance certificate.
"""
    temperature = np.asarray(T, dtype=float)
    if temperature.ndim != 1 or len(temperature) == 0 or not np.isfinite(temperature).all():
        raise ValueError('T must be a nonempty finite N vector')
    count = len(temperature)
    m = _vector(mass, count, 'mass')
    density = _vector(rho, count, 'rho')
    k = _vector(conductivity, count, 'conductivity', broadcast=True)
    capacity_specific = _vector(cp, count, 'cp', broadcast=True)
    plastic = _vector(plastic_work_J, count, 'plastic_work_J')
    friction = _vector(friction_raw_J, count, 'friction_raw_J')
    contact = _vector(contact_power_W, count, 'contact_power_W')
    wall_row = _vector(wall_row_W_K, count, 'wall_row_W_K')
    if np.any(m <= 0) or np.any(density <= 0) or np.any(k <= 0) or np.any(capacity_specific <= 0):
        raise ValueError('mass/rho/conductivity/cp must be positive')
    if np.any(plastic < 0):
        raise ValueError('negative plastic work was NOT clipped')
    if np.any(wall_row < 0):
        raise ValueError('wall conductance row must be nonnegative')
    h = _scalar(dt, 'dt'); fraction = _scalar(chi, 'chi')
    split = _scalar(partition, 'partition')
    lower = _scalar(room, 'room'); upper = _scalar(melt, 'melt')
    if h <= 0 or lower <= 0 or upper <= lower:
        raise ValueError('positive dt and ordered positive temperature bounds required')
    if not 0 <= fraction <= 1 or not 0 <= split <= 1:
        raise ValueError('chi/partition must lie in [0,1]')
    if np.any(temperature < lower) or np.any(temperature >= upper):
        raise ValueError('initial temperature outside supported solid domain')
    points = np.asarray(x, dtype=float)
    if points.shape != (count, 3) or not np.isfinite(points).all():
        raise ValueError('finite N by 3 positions required')
    first = np.asarray(i); second = np.asarray(j)
    if (first.ndim != 1 or second.shape != first.shape or
            first.dtype.kind not in 'iu' or second.dtype.kind not in 'iu'):
        raise ValueError('matching integer pair-index vectors required')
    gradient = np.asarray(g, dtype=float)
    if gradient.shape != (len(first), 3) or not np.isfinite(gradient).all():
        raise ValueError('finite pair by 3 raw radial gradients required')
    if len(first):
        if (np.any(first >= count) or np.any(second >= count) or
                np.any(first < 0) or np.any(second < 0)):
            raise ValueError('pair index outside particle inventory')
        if np.any(first == second):
            raise ValueError('self conduction pair is invalid')
        canonical = np.sort(np.column_stack((first, second)), axis=1)
        if len(np.unique(canonical, axis=0)) != len(first):
            raise ValueError('duplicate conduction pair would double count heat')
        separation = points[first]-points[second]
        if np.any(np.sum(separation*separation, axis=1) == 0):
            raise ValueError('coincident conduction pair is invalid')
        # The older thermal operator assumes a radial gradient. Fail explicitly
        # if a mechanically corrected, nonradial operator is passed by mistake.
        cross = np.linalg.norm(np.cross(separation, gradient), axis=1)
        scale = np.linalg.norm(separation, axis=1)*np.linalg.norm(gradient, axis=1)
        if np.any(cross > 128*np.finfo(float).eps*scale):
            raise ValueError('thermal conduction requires raw radial pair gradient')
    power, row = conduction(m, density, temperature, k, points, first, second, gradient)
    heat_capacity = m*capacity_specific
    total_row = row+wall_row
    if (not np.isfinite(power).all() or not np.isfinite(total_row).all() or
            not np.isfinite(heat_capacity).all() or np.any(heat_capacity <= 0)):
        raise ValueError('nonfinite/invalid thermal coefficient or capacity')
    bound = .25*float(np.min(np.divide(heat_capacity, total_row,
        out=np.full(count, np.inf), where=total_row > 0)))
    if h > bound:
        raise ValueError('explicit conduction time-step bound exceeded')
    particle = dict(plastic_work_J=plastic.copy(), plastic_heat_J=fraction*plastic,
        plastic_unthermalized_J=(1-fraction)*plastic, friction_raw_J=friction.copy(),
        friction_material_heat_J=split*friction, friction_wall_heat_J=(1-split)*friction,
        contact_in_J=h*contact, internal_conduction_J=h*power)
    increment = (particle['plastic_heat_J']+particle['friction_material_heat_J']+
                 particle['contact_in_J']+particle['internal_conduction_J'])
    if not np.isfinite(increment).all():
        raise ValueError('nonfinite thermal source increment')
    updated = temperature_step(temperature, increment, m, capacity_specific, lower, upper)
    measured = heat_capacity*(updated-temperature)
    particle.update(total_heat_input_J=increment,
        measured_thermal_increment_J=measured, mc_delta_T_residual_J=measured-increment)
    ledger = {key: math.fsum(float(v) for v in value) for key, value in particle.items()}
    ledger['wall_heat_out_J'] = ledger['friction_wall_heat_J']-ledger['contact_in_J']
    return dict(temperature=updated, pair_power_W=power, conduction_row_W_K=row,
        total_row_W_K=total_row, thermal_dt_bound_s=bound,
        particle_energy_J=particle, ledger=ledger,
        status='RAW_THERMAL_SOURCE_ACCOUNTING_NOT_MATERIAL_ENERGY_CERTIFICATION')

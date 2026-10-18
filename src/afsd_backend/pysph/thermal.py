"""Conservative P0 Table 1 conduction and explicit, separately metered sources."""
import numpy as np


def conduction(mass, rho, temperature, conductivity, x, i, j, gradient):
    k = np.broadcast_to(conductivity, temperature.shape)
    if np.any(k <= 0):
        raise ValueError('positive conductivity required')
    rij = x[i]-x[j]
    coefficient = mass[i]*mass[j]/(rho[i]*rho[j]) * (4*k[i]*k[j]/(k[i]+k[j])) * (-np.sum(rij*gradient, axis=1))/np.sum(rij*rij, axis=1)
    if np.any(coefficient < 0):
        raise ValueError('heat operator has wrong sign')
    pair_power = coefficient*(temperature[j]-temperature[i])
    power = np.zeros_like(temperature)
    np.add.at(power, i, pair_power); np.add.at(power, j, -pair_power)
    row = np.zeros_like(temperature)
    np.add.at(row, i, coefficient); np.add.at(row, j, coefficient)
    return power, row


def temperature_step(temperature, energy_increment, mass, cp, room, melt):
    result = temperature+energy_increment/(mass*cp)
    if not np.isfinite(result).all() or np.any(result < room) or np.any(result >= melt):
        raise ValueError('thermal domain violation; temperature was NOT clipped')
    return result


def effusivity_partition(k1, cp1, rho1, k2, cp2, rho2):
    """P44 Eq3-110, semi-infinite perfect-contact heat split."""
    a = np.sqrt(k1*cp1*rho1); b = np.sqrt(k2*cp2*rho2)
    return a/(a+b)

"""Raw endpoint-work diagnostics for the independent quasistatic model.

These identities do not certify a finite-strain free energy. Algorithmic terms
are reported, never made heat or subtracted from the process acceptance screen.
The radial-return diagnostic and BE elastic dissipation are different accounts;
adding both to endpoint stress work would double count the return operation.
"""
import numpy as np

from .material import rotation
from .mechanics import volumetric_energy


def material_energy_breakdown(S0, S1, rho0, rho1, mass, L, dt, G, K,
                              rho_reference, plastic_work_density,
                              return_dissipation_density):
    """Return sums and local signed diagnostics; inputs and thresholds unchanged.

    Stresses/G/K are Pa, L is 1/s, densities kg/m^3, mass kg, dt s, and
    plastic/return work densities Pa. Every ``*_J`` field is in joules.
    ``rho0`` is the accepted step-start density, not the reference density.
    The ``per_particle`` dictionary preserves local cancellation information.
    """
    S0 = np.asarray(S0, float); S1 = np.asarray(S1, float)
    L = np.asarray(L, float)
    if (S0.ndim != 3 or S0.shape[1:] != (3, 3) or len(S0) == 0 or
            S1.shape != S0.shape or L.shape != S0.shape):
        raise ValueError('equal nonempty (N,3,3) stresses and gradient required')
    count = len(S0)
    vectors = []
    for name, value in (('rho0', rho0), ('rho1', rho1), ('mass', mass),
                        ('rho_reference', rho_reference),
                        ('plastic_work_density', plastic_work_density),
                        ('return_dissipation_density', return_dissipation_density)):
        value = np.asarray(value, float)
        if value.ndim == 0:
            value = np.full(count, float(value))
        if value.shape != (count,) or not np.isfinite(value).all():
            raise ValueError(f'finite scalar or (N,) {name} required')
        vectors.append(value)
    rho0, rho1, mass, rho_reference, plastic_density, return_density = vectors
    if (not np.isfinite([dt, G, K]).all() or min(dt, G, K) <= 0. or
            not all(np.isfinite(a).all() for a in (S0, S1, L)) or
            any(np.any(a <= 0.) for a in (rho0, rho1, mass, rho_reference))):
        raise ValueError('finite fields and positive dt, moduli, densities, mass required')
    # Negative input work is preserved in diagnostics, never silently clipped.
    H = dt*L
    dilation = np.trace(H, axis1=1, axis2=2)
    dev_increment = .5*(H+H.swapaxes(1, 2))-dilation[:, None, None]*np.eye(3)/3.
    half = rotation(.25*(H-H.swapaxes(1, 2)))
    full = half@half
    rotated_old = full@S0@full.swapaxes(1, 2)
    A = half@dev_increment@half.swapaxes(1, 2)
    trial = rotated_old+2.*G*A
    V0 = mass/rho0; V1 = mass/rho1
    e0 = np.sum(S0*S0, axis=(1, 2))/(4.*G)
    e1 = np.sum(S1*S1, axis=(1, 2))/(4.*G)
    etrial = np.sum(trial*trial, axis=(1, 2))/(4.*G)
    p1 = K*(rho1/rho_reference-1.)

    dUdev = V1*e1-V0*e0
    Wdev = V1*np.sum(S1*dev_increment, axis=(1, 2))
    plastic = V1*plastic_density
    return_diagnostic = V1*return_density
    BEdev = V1*np.sum((S1-rotated_old)**2, axis=(1, 2))/(4.*G)
    Espin = V1*np.sum(S1*(dev_increment-A), axis=(1, 2))
    Emeasure = -(V1-V0)*e0
    dUvol = mass*(volumetric_energy(rho1, rho_reference, K)-
                  volumetric_energy(rho0, rho_reference, K))
    Wvol = -V1*p1*dilation
    BEvol = -p1*(V1-V0)-dUvol
    # The actual signed defect is meaningful even on a rejected density trial.
    Erho = Wvol+p1*(V1-V0)
    Erho_expected = -p1*(V1-V0)*dilation
    raw = Wdev+Wvol-dUdev-dUvol-plastic
    explained = BEdev+BEvol+Espin+Emeasure+Erho
    local = dict(
        deviatoric_endpoint_work_J=Wdev,
        volumetric_endpoint_work_J=Wvol,
        deviatoric_storage_change_J=dUdev,
        volumetric_storage_change_J=dUvol,
        plastic_work_J=plastic,
        return_map_diagnostic_J=return_diagnostic,
        be_elastic_deviatoric_dissipation_J=BEdev,
        be_volumetric_dissipation_J=BEvol,
        deviatoric_spin_defect_J=Espin,
        deviatoric_volume_measure_defect_J=Emeasure,
        volumetric_density_be_defect_J=Erho,
        volumetric_density_be_expected_defect_J=Erho_expected,
        density_be_work_identity_error_J=Erho-Erho_expected,
        material_raw_energy_defect_J=raw,
        material_discrete_identity_error_J=raw-explained,
        deviatoric_discrete_identity_error_J=Wdev-dUdev-plastic-BEdev-Espin-Emeasure,
        return_map_identity_error_J=V1*(etrial-e1)-plastic-return_diagnostic,
    )
    result = {key: float(np.sum(value)) for key, value in local.items()}
    result.update(
        per_particle={key: value.copy() for key, value in local.items()},
        material_raw_energy_defect_absolute_sum_J=float(np.sum(np.abs(raw))),
        material_discrete_identity_error_absolute_sum_J=float(np.sum(np.abs(raw-explained))),
        density_closure_relative_max=float(np.max(np.abs(rho1*(1.+dilation)-rho0)/rho0)),
        minimum_be_volumetric_dissipation_J=float(np.min(BEvol)),
        minimum_input_plastic_work_density_Pa=float(np.min(plastic_density)),
        minimum_input_return_dissipation_density_Pa=float(np.min(return_density)),
        material_energy_certified=False,
        algorithmic_terms_are_heat=False,
        raw_acceptance_threshold_modified=False,
        return_diagnostic_is_additional_endpoint_dissipation=False,
        classification='QUASISTATIC_ENDPOINT_ENERGY_ALGEBRA_NOT_FREE_ENERGY_CERTIFICATION')
    return result

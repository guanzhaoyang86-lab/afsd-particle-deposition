"""Exact bookkeeping for the named impulse-corrected HHT contact variant.

This module neither integrates contact nor commits a state. Physical endpoint
forces remain physical forces; numerical acceleration/history is not relabeled.
"""
import numpy as np
from . import hht_path


def kinetic_identity(mass,v0,v1,acc0,acc1,force0,force1,dt,
                     delta_impulse_Ns,delta_material_work_J,weight=.9):
    """Return per-particle balances, including signed curved-force filter work.

Correction work is the actual MATERIAL-velocity work minus the endpoint chord
work, not heat or wall-relative friction work. Never broadcast an aggregate
energy over particles. Inconsistent returned velocities remain visible in the
kinematic, momentum and energy residuals rather than silently being replaced.
"""
    old=hht_path.kinetic_identity(mass,v0,v1,acc0,acc1,force0,force1,dt,weight)
    m=np.asarray(mass,float);v0=np.asarray(v0,float);v1=np.asarray(v1,float)
    a0=np.asarray(acc0,float);a1=np.asarray(acc1,float)
    dj=np.asarray(delta_impulse_Ns,float);dw=np.asarray(delta_material_work_J,float)
    if dj.shape!=v0.shape or not np.isfinite(dj).all():
        raise ValueError('finite per-particle correction impulse matching velocity required')
    if dw.shape!=m.shape or not np.isfinite(dw).all():
        raise ValueError('finite per-particle material-work correction matching mass required')
    mean=.5*(v0+v1)
    baseline_impulse=old['physical_impulse_Ns']
    physical_impulse=baseline_impulse+dj
    physical_work=old['physical_linear_force_work_J']+dw
    correction_filter_work=np.sum(dj*mean,axis=-1)-dw
    filter_work=old['algorithmic_filter_work_J']+correction_filter_work
    residual_impulse=old['root_residual_impulse_Ns']-dj
    root_work=np.sum(residual_impulse*mean,axis=-1)
    gamma,_=hht_path.parameters(weight)
    velocity_defect=v1-v0-dt*((1-gamma)*a0+gamma*a1)
    momentum_defect=(m[...,None]*(v1-v0)-physical_impulse-
        old['algorithmic_filter_impulse_Ns']-residual_impulse)
    energy_defect=old['delta_kinetic_J']-physical_work-filter_work-root_work
    return dict(integration_definition='HHT_CONTACT_IMPULSE_CORRECTION_V1',
        physical_impulse_Ns=physical_impulse,
        baseline_linear_force_impulse_Ns=baseline_impulse.copy(),
        contact_impulse_correction_Ns=dj.copy(),
        physical_force_work_J=physical_work,
        baseline_linear_force_work_J=old['physical_linear_force_work_J'].copy(),
        contact_material_work_correction_J=dw.copy(),
        algorithmic_filter_impulse_Ns=old['algorithmic_filter_impulse_Ns'].copy(),
        baseline_linear_algorithmic_filter_work_J=old['algorithmic_filter_work_J'].copy(),
        curved_algorithmic_filter_work_J=correction_filter_work,
        algorithmic_filter_work_J=filter_work,
        root_residual_impulse_Ns=residual_impulse,
        corrected_root_force_residual_N=residual_impulse/(gamma*dt),
        root_residual_work_J=root_work,
        kinematic_velocity_residual_m_s=velocity_defect,
        momentum_identity_residual_Ns=momentum_defect,
        delta_kinetic_J=old['delta_kinetic_J'].copy(),
        kinetic_identity_residual_J=energy_defect,
        algorithmic_work_is_not_heat=True,algorithmic_work_sign_not_assumed=True,
        physical_force_history_definition='TRUE_ENDPOINT_FORCE_NOT_CORRECTION_EFFECTIVE_FORCE',
        acceleration_history_definition='NUMERICAL_HHT_NOT_ASSUMED_FORCE_DIVIDED_BY_MASS')

"""Uncommitted updated-configuration, event-average coupled trial.

This named stage-implicit APPROXIMATE map reuses the unchanged spatial,
Jaumann/J2/FKS, event contact and explicit thermal implementations. It is not
a full implicit constitutive free-energy model. No accepted history advances
here, and the actual event-average velocity is a separate unknown from v1.
"""
import copy
import hashlib
import numpy as np

from .mechanics import (positions, velocities, corrected_operators,
                        internal_force, pressure, bulk_modulus)
from .material import FKS, update
from .axial_actuator import validate_interval
from .coupled_boundary import evaluate_boundary, normal_snapshot
from .coupled_thermal import advance_thermal
from .coupled_solver import IntegratedSolver

CLASSIFICATION = 'UPDATED_CONFIGURATION_EVENT_AVERAGE_STAGE_IMPLICIT_APPROXIMATE'
FIELDS = ('position_m', 'velocity_m_s', 'density_kg_m3', 'temperature_K',
          'mean_velocity_m_s')


def particle_hash(pa):
    """Exact inventory/property hash, including history and passive labels."""
    result = hashlib.sha256()
    for key in sorted(pa.properties):
        array = pa.properties[key].get_npy_array()
        result.update(key.encode()); result.update(array.dtype.str.encode())
        result.update(str(array.shape).encode()); result.update(array.tobytes())
    return result.hexdigest()


class _FeatureGuard:
    # Reuse the existing finite-feature loss test verbatim, without constructing
    # an Integrator or computing another material/contact kick.
    _continued_gap = IntegratedSolver._continued_gap
    _feature_transition_guard = IntegratedSolver._feature_transition_guard

    def __init__(self, pa, p):
        self.pa, self.p, self.kn = pa, p, p['normal_stiffness']


class UpdatedTrial:
    """Read-only accepted snapshot plus pure same-history candidate evaluation.

``law`` is ONLY an explicit manufactured-material test hook. Omitting it uses
the original FKS with the unchanged p values. Returned outputs never alias the
accepted ParticleArray, its property history, or the caller's candidate.
"""
    def __init__(self, accepted_pa, p, time, law=None):
        self.p = copy.deepcopy(p)
        if not np.isfinite(time) or time < 0:
            raise ValueError('finite nonnegative accepted time required')
        self.time = float(time)
        if self.p['spatial_operator'] != 'corrected_work_conjugate':
            raise ValueError('paired corrected spatial operators required')
        if self.p['lift_speed'] != 0 or self.p['dwell_duration'] != 0:
            raise ValueError('registered fixed-height hot-start programme required')
        expected = (p['penalty_factor']*p['shear']*p['dx']**2/
                    p['contact_reference_spacing'])
        if p['normal_stiffness'] != expected:
            raise ValueError('undeclared normal stiffness change')
        count = accepted_pa.get_number_of_particles()
        if not count:
            raise ValueError('nonempty accepted particle inventory required')
        self._accepted = accepted_pa.extract_particles(np.arange(count, dtype=np.uint32))
        self._order = np.arange(count, dtype=np.uint32)
        self.count = count
        self.bulk = bulk_modulus(p['shear'], p['poisson'])
        self.law = law if law is not None else FKS(**{k: p[k] for k in FKS.__dataclass_fields__})
        self.material_scope = ('MANUFACTURED_TEST_ONLY' if law is not None else
                               'UNCHANGED_FKS_APPROXIMATE_PARAMETERS')
        self._feature = _FeatureGuard(self._accepted, self.p)
        for key, prop in self._accepted.properties.items():
            if not np.isfinite(prop.get_npy_array()).all():
                raise ValueError('nonfinite accepted property '+key)
        pa = self._accepted
        if (np.any(pa.m <= 0) or np.any(pa.h <= 0) or np.any(pa.ep < 0) or
                len(np.unique(pa.particle_id)) != count):
            raise ValueError('invalid mass/smoothing length/history/identity')
        self._check_state(self.candidate(), self.time)
        self._accepted_hash = particle_hash(pa)

    def accepted_hash(self):
        return particle_hash(self._accepted)

    def candidate(self):
        pa = self._accepted
        return dict(position_m=positions(pa), velocity_m_s=velocities(pa),
                    density_kg_m3=pa.rho.copy(), temperature_K=pa.temperature.copy(),
                    mean_velocity_m_s=velocities(pa))

    def _validate_arrays(self, candidate):
        if not isinstance(candidate, dict) or set(candidate) != set(FIELDS):
            raise ValueError('candidate requires exactly five registered fields')
        result = {}
        for key in FIELDS:
            shape = (self.count,) if key in ('density_kg_m3', 'temperature_K') else (self.count, 3)
            value = np.asarray(candidate[key], dtype=float)
            if value.shape != shape or not np.isfinite(value).all():
                raise ValueError('invalid/nonfinite candidate '+key)
            result[key] = value.copy()
        return result

    def _check_state(self, state, time):
        p, pa = self.p, self._accepted
        rho, temperature = state['density_kg_m3'], state['temperature_K']
        if np.any(rho <= 0):
            raise ValueError('nonpositive candidate density')
        if np.max(np.abs(rho/p['rho0']-1)) > p['density_limit']:
            raise ValueError('density change limit')
        if np.any(temperature < p['room']) or np.any(temperature >= p['melt']):
            raise ValueError('thermal solid-domain limit')
        for key in ('velocity_m_s', 'mean_velocity_m_s'):
            if np.linalg.norm(state[key], axis=1).max() > p['max_speed']:
                raise ValueError('speed limit')
        original = np.column_stack((pa.x0, pa.y0, pa.z0))
        if np.linalg.norm(state['position_m']-original, axis=1).max() > p['maximum_displacement_m']:
            raise ValueError('maximum displacement limit')
        snapshot = normal_snapshot(state['position_m'], pa.particle_id, time, p)
        if snapshot['max_penetration_m'] > p['penetration_limit']*p['dx']:
            raise ValueError('contact penetration limit')
        return snapshot

    def evaluate(self, candidate, dt):
        """Return same-unit map defects and all uncommitted physics journals.

Pair conduction is explicitly T0, FKS/contact flux are candidate Tmid. Both
mechanical volume and thermal conductance use candidate rhomid. The resulting
T map always starts at T0; no iteration can accumulate a heat/plastic source.
"""
        if self.accepted_hash() != self._accepted_hash:
            raise RuntimeError('accepted private snapshot changed')
        dt = float(dt)
        validate_interval(self.time, dt, self.p['release_time'])
        c = self._validate_arrays(candidate)
        p, old = self.p, self._accepted
        end_time, mid_time = self.time+dt, self.time+.5*dt
        normal_after = self._check_state(c, end_time)
        normal_before = normal_snapshot(positions(old), old.particle_id, self.time, p)
        pa = old.extract_particles(self._order)
        x0, v0 = positions(old), velocities(old)
        xmid = .5*(x0+c['position_m'])
        rhomid = .5*(old.rho+c['density_kg_m3'])
        tmid = .5*(old.temperature+c['temperature_K'])
        pa.x[:], pa.y[:], pa.z[:] = xmid.T
        pa.u[:], pa.v[:], pa.w[:] = c['mean_velocity_m_s'].T
        pa.rho[:] = rhomid
        pa.temperature[:] = tmid
        i, j, g, gradient, drho, correction = corrected_operators(pa)
        stress0 = old.S.reshape(-1, 3, 3).copy()
        stress1, ep1, plastic_density, return_density, rate = update(
            stress0, old.ep.copy(), gradient, tmid, dt, p['shear'], self.law)
        sigma = .5*(stress0+stress1)-pressure(rhomid, p['rho0'], self.bulk)[:, None, None]*np.eye(3)
        force = internal_force(pa, sigma, i, j, g, correction=correction)
        boundary = evaluate_boundary(xmid, v0, old.m, tmid, mid_time, dt,
                                     p, force, old.z0, old.particle_id)
        normal_mid = boundary['normal_snapshot']
        guard = self._feature._feature_transition_guard
        guards = dict(before_to_mid=guard(normal_before, normal_mid, x0, xmid, self.time, mid_time),
                      mid_to_after=guard(normal_mid, normal_after, xmid, c['position_m'], mid_time, end_time))
        volume = old.m/rhomid
        plastic_j, return_j = volume*plastic_density, volume*return_density
        thermal = advance_thermal(old.m, rhomid, old.temperature, p['conductivity'],
            xmid, i, j, g, dt, p['cp'], plastic_j, boundary['friction_raw_J'],
            boundary['contact_power_W'], boundary['wall_row_W_K'],
            p['plastic_heat_fraction'], p['heat_partition'], p['room'], p['melt'])
        divergence = np.trace(gradient, axis1=1, axis2=2)
        denominator = 1+.5*dt*divergence
        numerator = 1-.5*dt*divergence
        if np.any(denominator <= 0) or np.any(numerator <= 0):
            raise ValueError('nonpositive midpoint continuity density factor')
        mapped = dict(position_m=x0+dt*boundary['work_velocity'],
            velocity_m_s=boundary['velocity'].copy(), mean_velocity_m_s=boundary['work_velocity'].copy(),
            density_kg_m3=old.rho*numerator/denominator, temperature_K=thermal['temperature'].copy())
        # A fixed-point map outside the registered domain is not clipped into a
        # plausible candidate. Root callers must report a rejected attempt.
        self._check_state(mapped, end_time)
        internal_power = float(np.sum(force*c['mean_velocity_m_s']))
        stress_power = float(np.sum(volume[:, None, None]*sigma*gradient))
        condition = np.linalg.cond(np.linalg.inv(correction))
        spatial = dict(pair_count=len(i), condition_max=float(condition.max()),
            internal_power_W=internal_power, midpoint_stress_power_W=stress_power,
            power_conjugacy_residual_W=internal_power+stress_power,
            internal_net_force_N=force.sum(axis=0), internal_net_torque_Nm=np.cross(xmid, force).sum(axis=0),
            force_absolute_scale_N=float(np.linalg.norm(force, axis=1).sum()),
            torque_absolute_scale_Nm=float(np.linalg.norm(np.cross(xmid, force), axis=1).sum()),
            power_absolute_scale_W=float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(c['mean_velocity_m_s'], axis=1))))
        if self.accepted_hash() != self._accepted_hash:
            raise RuntimeError('trial mutated accepted history')
        return dict(classification=CLASSIFICATION, accepted=False, restart_allowed=False,
            material_scope=self.material_scope, accepted_state_hash=self._accepted_hash,
            time_before_s=self.time, time_after_s=end_time, time_mid_s=mid_time, dt_s=dt,
            particle_id=old.particle_id.copy(), map=mapped,
            residual={key: c[key]-mapped[key] for key in FIELDS},
            continuity_residual_kg_m3_s=(c['density_kg_m3']-old.rho)/dt+rhomid*divergence,
            midpoint=dict(position_m=xmid, density_kg_m3=rhomid, temperature_K=tmid,
                mean_velocity_m_s=c['mean_velocity_m_s'], velocity_gradient_per_s=gradient,
                stress_Pa=sigma, i=i, j=j, g=g, B=correction, mass_kg=old.m.copy()),
            material=dict(deviatoric_stress_Pa=stress1, ep=ep1, plastic_work_J=plastic_j,
                return_dissipation_J=return_j, strain_rate_per_s=rate),
            internal_force_N=force, boundary=boundary, thermal=thermal, spatial=spatial,
            normal_before=normal_before, normal_after=normal_after, feature_guards=guards,
            temperature_layers='FKS_AND_CONTACT_TMID; PAIR_CONDUCTION_T0; INCREMENT_FROM_T0',
            material_energy_certified=False, production_authorized=False)

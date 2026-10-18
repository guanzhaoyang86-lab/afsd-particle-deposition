"""Bounded 3D thermo-elastoplastic DEBUG physics, with no label-dependent forces."""
import numpy as np
from .integrator import HostLeapfrog
from .material import FKS, update
from .mechanics import positions, velocities, operators, corrected_operators, internal_force, pressure, bulk_modulus, volumetric_energy
from .thermal import conduction, temperature_step
from .contact import WallContact


class SolidSolver:
    def __init__(self, pa, p, state=None, diagnostic_recorder=None, contact_variant='legacy'):
        self.pa, self.p = pa, p
        self.diagnostic_recorder = diagnostic_recorder
        if contact_variant not in ('legacy', 'event_aware_v1'):
            raise ValueError('unknown contact variant')
        if diagnostic_recorder is not None and contact_variant != 'legacy':
            raise ValueError('legacy KickRecorder cannot observe another contact law')
        if state and state.get('contact_variant', 'legacy') != contact_variant:
            raise ValueError('checkpoint contact variant mismatch')
        self.contact_variant = contact_variant
        self.last_contact_trace = []
        self.bulk = bulk_modulus(p['shear'], p['poisson'])
        self.sound = np.sqrt((self.bulk+4*p['shear']/3)/p['rho0'])
        # Hold traction stiffness fixed when changing numerical dx.
        self.kn = p['penalty_factor']*p['shear']*p['dx']**2/p['contact_reference_spacing']
        self.contact = WallContact(self.kn, p['friction'])
        self.law = FKS(**{k: p[k] for k in FKS.__dataclass_fields__})
        self.integrator = HostLeapfrog(pa, self.evaluate)
        self.time = 0.; self.count = 0
        self.ledger = dict(wall_work_J=0., actuator_work_J=0., wall_heat_out_J=0., plastic_heat_J=0., friction_heat_J=0., contact_heat_J=0., numerical_return_dissipation_J=0., external_mass_kg=0., external_feed_enthalpy_J=0.)
        self.last = {}; self.floor_evaluations = 0
        self.next_id = int(pa.particle_id.max())+1
        if state:
            self.time, self.count = state['t'], state['count']
            self.ledger = state['ledger']; self.initial_energy = state['initial_energy']
            self.initial_mass = state['initial_mass']; self.floor_evaluations = state['floor_evaluations']
            self.next_id = state['next_id']
            self.last = state['last']
        else:
            self.initial_mass = float(pa.m.sum())
            self.initial_energy = self.energy()['total_J']

    def stable_dt(self):
        p = self.p
        wave = p['cfl']*p['dx']*p['h_ratio']/(self.sound+np.linalg.norm(velocities(self.pa), axis=1).max())
        contact = .1*np.sqrt(self.pa.m.min()/self.kn)
        thermal = .1*p['rho0']*p['cp']*p['dx']**2/p['conductivity']
        wall_row = self.contact.conductance_row(positions(self.pa), self.time, p)
        wall_thermal = np.min(np.divide(.25*self.pa.m*p['cp'], wall_row, out=np.full_like(wall_row, np.inf), where=wall_row > 0))
        return float(min(wave, contact, thermal, wall_thermal))

    def evaluate(self, time, dt):
        pa, p = self.pa, self.p
        correction = None
        if p['spatial_operator'] == 'corrected_work_conjugate':
            i, j, g, L, drho, correction = corrected_operators(pa)
        elif p['spatial_operator'] == 'paper_table1_raw_audit':
            i, j, g, L, drho = operators(pa)
        else:
            raise ValueError('spatial operator not explicitly selected')
        old_s = pa.S.reshape(-1, 3, 3).copy()
        new_s, ep, plastic_density, return_density, rate = update(old_s, pa.ep, L, pa.temperature, dt, p['shear'], self.law)
        rho_mid = pa.rho+dt/2*drho
        sigma = .5*(old_s+new_s)-pressure(rho_mid, p['rho0'], self.bulk)[:, None, None]*np.eye(3)
        force = internal_force(pa, sigma, i, j, g, correction=correction)
        diagnostic_internal = force.copy() if self.diagnostic_recorder is not None else None
        if self.contact_variant == 'event_aware_v1':
            from .contact_kick import domain_kick
            acceleration, contact_heat, ledger, actuator_work, trace = domain_kick(
                self.contact, positions(pa), velocities(pa), pa.m, pa.temperature,
                time, dt, p, force, pa.actuator, pa.particle_id)
            self.ledger['actuator_work_J'] += actuator_work
            self.last_contact_trace = trace
        else:
            contact_force, contact_heat, ledger = self.contact.evaluate(positions(pa), velocities(pa), pa.m, pa.temperature, time, dt, p)
            force += contact_force
            acceleration = force/pa.m[:, None]
            # Only top-cap actuator x,z DOFs are imposed. All y and all free-body DOFs respond to force.
            cap = pa.actuator == 1
            old_v = velocities(pa)
            unconstrained = old_v+dt*acceleration
            for axis, prescribed in ((0, p['traverse']), (2, -p['feed'])):
                change = prescribed-unconstrained[cap, axis]
                midpoint = .5*(old_v[cap, axis]+prescribed)
                self.ledger['actuator_work_J'] += float(np.sum(pa.m[cap]*change*midpoint))
                acceleration[cap, axis] = (prescribed-old_v[cap, axis])/dt
            # Heat uses work done during this actual kick, including all forces and
            # actuator constraints. Do not turn old-velocity power into excess heat.
            if self.diagnostic_recorder is not None:
                self.diagnostic_recorder.record(self, time, dt, diagnostic_internal,
                    contact_force, old_v, acceleration, L, new_s, ep, drho)
            _, contact_heat, ledger = self.contact.evaluate(
                positions(pa), old_v, pa.m, pa.temperature, time, dt, p,
                heat_velocity=old_v+.5*dt*acceleration)
        pa.au[:], pa.av[:], pa.aw[:] = acceleration.T
        pa.arho[:] = drho
        pa.ep_rate[:] = (ep-pa.ep)/dt
        pa.S[:] = new_s.ravel(); pa.ep[:] = ep; pa.strain_rate[:] = rate
        pa.D[:] = (.5*(L+L.swapaxes(1, 2))).ravel()
        plastic_j = pa.m/pa.rho*plastic_density
        return_j = pa.m/pa.rho*return_density
        conductive_power, row = conduction(pa.m, pa.rho, pa.temperature, p['conductivity'], positions(pa), i, j, g)
        row += self.contact.conductance_row(positions(pa), time, p)
        thermal_dt = np.min(np.divide(pa.m*p['cp'], row, out=np.full_like(row, np.inf), where=row > 0))
        if dt > .25*thermal_dt: raise RuntimeError('explicit conduction time-step bound exceeded')
        pa.temperature[:] = temperature_step(pa.temperature, dt*(conductive_power+contact_heat)+p['plastic_heat_fraction']*plastic_j, pa.m, p['cp'], p['room'], p['melt'])
        pa.plastic_work[:] += plastic_j; pa.return_dissipation[:] += return_j
        self.ledger['wall_work_J'] += dt*ledger['wall_work_rate']
        self.ledger['wall_heat_out_J'] += dt*ledger['wall_heat_rate']
        self.ledger['plastic_heat_J'] += float(plastic_j.sum())
        self.ledger['numerical_return_dissipation_J'] += float(return_j.sum())
        self.ledger['friction_heat_J'] += dt*ledger['friction_power']*p['heat_partition']
        self.ledger['contact_heat_J'] += dt*ledger['contact_heat_rate']
        self.floor_evaluations += int(np.count_nonzero(rate < p['rate_floor']))
        self.last = {**ledger, 'pair_count': len(i), 'conduction_net_W': float(conductive_power.sum()), 'dt_s': dt}

    def advance(self, max_dt=None):
        dt = self.stable_dt()
        if max_dt is not None:
            if not np.isfinite(max_dt) or max_dt <= 0: raise ValueError('positive finite max_dt required')
            dt = min(dt, max_dt)
        self.integrator.advance(self.time, dt)
        self.time += dt; self.count += 1
        crossed = (self.pa.component_id == 1) & (self.pa.z0 >= self.p['layer']) & (self.pa.z < self.p['layer'])
        self.pa.feed_crossed[crossed] = 1
        return dt

    def energy(self):
        pa, p = self.pa, self.p
        s = pa.S.reshape(-1, 3, 3)
        _, _, contact = self.contact.evaluate(positions(pa), velocities(pa), pa.m, pa.temperature, self.time, self.stable_dt(), p)
        result = dict(kinetic_J=float(np.sum(.5*pa.m*np.sum(velocities(pa)**2, axis=1))), elastic_shear_J=float(np.sum(pa.m/pa.rho*np.sum(s*s, axis=(1, 2))/(4*p['shear']))), elastic_volume_J=float(np.sum(pa.m*volumetric_energy(pa.rho, p['rho0'], self.bulk))), thermal_increment_J=float(np.sum(pa.m*p['cp']*(pa.temperature-p['temperature']))), contact_potential_J=contact['contact_potential'])
        result['total_J'] = sum(result.values())
        result['max_penetration_m'] = contact['max_penetration']
        return result

    def diagnostics(self, fail=True):
        pa, p = self.pa, self.p
        energy = self.energy()
        expected = self.initial_energy+self.ledger['wall_work_J']+self.ledger['actuator_work_J']-self.ledger['wall_heat_out_J']
        residual = energy['total_J']-expected
        scale = max(abs(self.initial_energy), abs(self.ledger['wall_work_J'])+abs(self.ledger['actuator_work_J'])+abs(self.ledger['wall_heat_out_J']), 1e-15)
        density_change = float(np.max(np.abs(pa.rho/p['rho0']-1)))
        errors = []
        for name, array in pa.properties.items():
            if not np.isfinite(array.get_npy_array()).all(): errors.append('NaN/Inf '+name)
        if abs(pa.m.sum()-self.initial_mass) > self.initial_mass*1e-12: errors.append('mass conservation')
        if density_change > p['density_limit']: errors.append('density limit')
        if energy['max_penetration_m']/p['dx'] > p['penetration_limit']: errors.append('contact penetration limit')
        if np.max(np.linalg.norm(velocities(pa), axis=1)) > p['max_speed']: errors.append('speed limit')
        if abs(residual) > p['energy_absolute_limit']+p['energy_relative_limit']*scale: errors.append('energy residual limit')
        result = {**energy, **self.ledger, **self.last, 'time_s': self.time, 'step': self.count, 'mass_kg': float(pa.m.sum()), 'energy_residual_J': residual, 'energy_residual_relative': residual/scale, 'density_relative_change_max': density_change, 'temperature_min_K': float(pa.temperature.min()), 'temperature_max_K': float(pa.temperature.max()), 'ep_max': float(pa.ep.max()), 'speed_max_m_s': float(np.linalg.norm(velocities(pa), axis=1).max()), 'rate_floor_evaluations': self.floor_evaluations, 'feed_plane_crossed_mass_kg': float(pa.m[pa.feed_crossed == 1].sum()), 'linear_momentum_kg_m_s': np.sum(pa.m[:, None]*velocities(pa), axis=0).tolist(), 'angular_momentum_kg_m2_s': np.sum(pa.m[:, None]*np.cross(positions(pa), velocities(pa)), axis=0).tolist(), 'alarms': errors}
        if fail and errors: raise RuntimeError('; '.join(errors)+f'; diagnostics={result}')
        return result

    def state(self):
        return dict(t=self.time, count=self.count, ledger=self.ledger, initial_energy=self.initial_energy, initial_mass=self.initial_mass, floor_evaluations=self.floor_evaluations, next_id=self.next_id, last=self.last, contact_variant=self.contact_variant)

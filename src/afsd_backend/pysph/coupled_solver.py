"""One bounded coupled AFSD adapter; old numerical/physical sources unchanged.

Updated 3D neighborhoods, original Jaumann/J2/FKS, centre-collocated multi-face
contact, original-material z actuator, and conservative thermal split. This is
an ENGINEERING_APPROXIMATE hot-start model, not validated sustained deposition.
"""
import copy
import hashlib
import json
import numpy as np

from .solver import SolidSolver
from .integrator import HostLeapfrog
from .mechanics import (positions, velocities, corrected_operators, internal_force,
                        pressure, volumetric_energy)
from .material import update, equivalent
from .geometry import local_xy, tool_pose
from .axial_actuator import validate_interval, control
from .coupled_boundary import evaluate_boundary, normal_snapshot
from .coupled_thermal import advance_thermal
from .output import checkpoint as save_checkpoint, restore
from afsd_scaffold.vtk_legacy import write_point_cloud

CLASSIFICATION = 'ENGINEERING_APPROXIMATE_COUPLED_AFSD_LOCAL_DEBUG'


def _clean(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: _clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [_clean(v) for v in value]
    return value


def parameter_hash(p):
    return hashlib.sha256(json.dumps(_clean(p), sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def frozen_gradient(v, mass, rho, i, j, g, correction):
    raw = np.zeros((len(v), 3, 3))
    pair = (v[j]-v[i])[:, :, None]*g[:, None, :]
    np.add.at(raw, i, (mass[j]/rho[j])[:, None, None]*pair)
    np.add.at(raw, j, (mass[i]/rho[i])[:, None, None]*pair)
    return raw@correction


class IntegratedSolver(SolidSolver):
    def __init__(self, pa, p, state=None, source_hashes=None, configuration_hash=None):
        p = copy.deepcopy(p)
        if p['spatial_operator'] != 'corrected_work_conjugate':
            raise ValueError('coupled adapter requires paired corrected spatial operators')
        if p['lift_speed'] != 0 or p['dwell_duration'] != 0:
            raise ValueError('this registered coupled adapter is fixed-height hot start')
        expected_kn = p['penalty_factor']*p['shear']*p['dx']**2/p['contact_reference_spacing']
        if p['normal_stiffness'] != expected_kn:
            raise ValueError('undeclared normal stiffness change')
        if p['release_time'] <= 0 or p['maximum_displacement_m'] <= 0:
            raise ValueError('explicit positive actuator programme and displacement bound required')
        self.source_hashes = dict(source_hashes or {})
        self.parameter_fingerprint = parameter_hash(p)
        self.configuration_hash = configuration_hash or self.parameter_fingerprint
        # Constructs the old physical constants, contact and pinned HostLeapfrog.
        # energy() below contains no legacy friction evaluation or extra kick.
        super().__init__(pa, p, state=None, contact_variant='event_aware_v1')
        self.integrator = HostLeapfrog(pa, self.evaluate)
        self.ledger = dict(wall_work_J=0., actuator_work_J=0., internal_work_J=0.,
            delta_kinetic_J=0., friction_raw_J=0., normal_material_work_J=0.,
            normal_wall_work_J=0., normal_relative_work_J=0., normal_temporal_residual_J=0.,
            plastic_work_J=0., plastic_heat_J=0., nonthermal_plastic_J=0.,
            friction_heat_J=0., friction_heat_to_wall_J=0., contact_heat_J=0.,
            wall_heat_out_J=0., numerical_return_dissipation_J=0.,
            internal_conduction_net_J=0., thermal_balance_residual_J=0.,
            material_old_work_J=0., midpoint_stress_work_J=0.,
            external_mass_kg=0., external_feed_enthalpy_J=0.)
        for body in ('base', 'tool', 'actuator'):
            for axis in 'xyz':
                self.ledger[f'{body}_impulse_{axis}_Ns'] = 0.
                self.ledger[f'{body}_angular_impulse_{axis}_Nms'] = 0.
        self.last_journal = None; self.last_contact_trace = []; self.last_thermal = None
        self.last = {}; self.floor_evaluations = 0
        order = np.argsort(pa.particle_id)
        initial_v = velocities(pa); initial_x = positions(pa)
        self.initial_kinetic_J = float(.5*np.sum((pa.m[:, None]*initial_v**2)[order]))
        self.initial_linear_momentum = np.sum((pa.m[:, None]*initial_v)[order], axis=0)
        self.initial_angular_momentum = np.sum((pa.m[:, None]*np.cross(initial_x, initial_v))[order], axis=0)
        self.initial_normal_potential_J = self.energy()['contact_potential_J']
        self.initial_particle_ids = sorted(int(q) for q in pa.particle_id)
        self.initial_energy = self.energy()['total_J']
        if state is not None:
            if (state.get('classification') != CLASSIFICATION or not state.get('accepted') or
                    not state.get('restart_allowed')):
                raise ValueError('only accepted coupled-adapter checkpoint may restart')
            if state['parameter_fingerprint'] != self.parameter_fingerprint or state['parameters'] != p:
                raise ValueError('checkpoint physical parameters mismatch')
            if state['configuration_hash'] != self.configuration_hash or state['source_hashes'] != self.source_hashes:
                raise ValueError('checkpoint configuration/source mismatch')
            if state['initial_particle_ids'] != self.initial_particle_ids:
                raise ValueError('checkpoint persistent particle inventory mismatch')
            self.time = float(state['t']); self.count = int(state['count'])
            self.ledger = copy.deepcopy(state['ledger']); self.last = copy.deepcopy(state['last'])
            self.initial_energy = state['initial_energy']; self.initial_mass = state['initial_mass']
            self.initial_kinetic_J = state['initial_kinetic_J']
            self.initial_linear_momentum = np.asarray(state['initial_linear_momentum_kg_m_s'], float)
            self.initial_angular_momentum = np.asarray(state['initial_angular_momentum_kg_m2_s'], float)
            self.initial_normal_potential_J = state['initial_normal_potential_J']
            self.floor_evaluations = state['floor_evaluations']; self.next_id = state['next_id']

    def _normal(self, time=None):
        return normal_snapshot(positions(self.pa), self.pa.particle_id,
                               self.time if time is None else time, self.p)

    def _continued_gap(self, feature, xyz, time):
        p = self.p
        if feature == 0: return float(xyz[2]+p['substrate_depth'])
        if feature == 1: return float(p['layer']-xyz[2])
        if feature not in (2, 3, 4, 5): raise ValueError('unknown finite contact feature')
        local, _ = local_xy(np.asarray(xyz).reshape(1, 3), time, p['traverse'], p['omega'])
        axis, sign = (feature-2)//2, (-1 if feature % 2 == 0 else 1)
        return float(p['bore_width']/2-sign*local[0, axis])

    def _feature_transition_guard(self, before, after, x_before, x_after, t_before, t_after):
        old = before['potential_by_contact_id_J']; new = after['potential_by_contact_id_J']
        lookup = {int(pid): k for k, pid in enumerate(self.pa.particle_id)}
        removed = sorted(set(old)-set(new)); added = sorted(set(new)-set(old))
        for ids, state, opposite_x, opposite_t, direction in (
                (removed, old, x_after, t_after, 'DISAPPEARED'),
                (added, new, x_before, t_before, 'APPEARED')):
            for cid in ids:
                if state[cid] <= 0.: continue
                pid, feature = (int(s) for s in str(cid).split(':'))
                delta = self.p['dx']/2-self._continued_gap(feature, opposite_x[lookup[pid]], opposite_t)
                if delta > 0.:
                    error = RuntimeError('UNSUPPORTED_STORED_ENERGY_FEATURE_TRANSITION '+direction+' '+str(cid))
                    error.feature_context = dict(contact_id=cid, transition=direction,
                        represented_potential_J=state[cid], opposite_extended_penetration_m=delta,
                        opposite_extended_potential_J=.5*self.kn*delta*delta)
                    raise error
        return dict(added_contact_ids=added, removed_contact_ids=removed,
                    status='FINITE_FEATURE_MASK_CHECKED; NATURAL_ZERO_POTENTIAL_OPENING_ALLOWED')

    def evaluate(self, time, dt):
        pa, p = self.pa, self.p
        x = positions(pa); v = velocities(pa); old_rho = pa.rho.copy(); old_t = pa.temperature.copy()
        i, j, g, L, drho, B = corrected_operators(pa)
        old_s = pa.S.reshape(-1, 3, 3).copy(); old_ep = pa.ep.copy()
        new_s, ep, plastic_density, returned_density, rate = update(
            old_s, old_ep, L, old_t, dt, p['shear'], self.law)
        rho_mid = old_rho+.5*dt*drho
        sigma = .5*(old_s+new_s)-pressure(rho_mid, p['rho0'], self.bulk)[:, None, None]*np.eye(3)
        force = internal_force(pa, sigma, i, j, g, correction=B)
        self.last_journal = dict(stage='CANDIDATE_MATERIAL_BEFORE_BOUNDARY',
            step=self.count+1, time_mid_s=time, dt_s=dt, particle_id=pa.particle_id.copy(),
            x_mid_m=x.copy(), v_before_m_s=v.copy(), S_before_Pa=old_s, S_after_Pa=new_s.copy(),
            ep_before=old_ep, ep_after=ep.copy(), internal_force_N=force,
            rho_before_kg_m3=old_rho, temperature_before_K=old_t, i=i, j=j, g=g, B=B,
            L_before_per_s=L, sigma_mid_Pa=sigma, normal_before=self._candidate_normal_before)
        boundary = evaluate_boundary(x, v, pa.m, old_t, time, dt, p, force, pa.z0, pa.particle_id)
        self.last_contact_trace = boundary['traces']
        volume = pa.m/old_rho; plastic_j = volume*plastic_density; return_j = volume*returned_density
        Lmid = frozen_gradient(boundary['work_velocity'], pa.m, old_rho, i, j, g, B)
        force_power = float(np.sum(force*boundary['work_velocity']))
        stress_power = float(np.sum(volume[:, None, None]*sigma*Lmid))
        condition = np.linalg.cond(np.linalg.inv(B))
        spatial = dict(pair_count=len(i), condition_max=float(condition.max()),
            internal_power_W=force_power, midpoint_stress_power_W=stress_power,
            power_conjugacy_residual_W=force_power+stress_power,
            internal_net_force_N=force.sum(axis=0), internal_net_torque_Nm=np.cross(x, force).sum(axis=0),
            force_absolute_scale_N=float(np.linalg.norm(force, axis=1).sum()),
            torque_absolute_scale_Nm=float(np.linalg.norm(np.cross(x, force), axis=1).sum()),
            power_absolute_scale_W=float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(boundary['work_velocity'], axis=1))))
        # Commit into a candidate ParticleArray; advance() rolls back everything
        # if thermal or endpoint guards reject. Temperature is never a force fit.
        pa.au[:], pa.av[:], pa.aw[:] = ((boundary['velocity']-v)/dt).T
        pa.arho[:] = drho; pa.ep_rate[:] = (ep-old_ep)/dt
        pa.S[:] = new_s.ravel(); pa.ep[:] = ep; pa.strain_rate[:] = rate
        pa.D[:] = (.5*(L+L.swapaxes(1, 2))).ravel()
        pa.plastic_work[:] += plastic_j; pa.return_dissipation[:] += return_j
        self.last_journal.update(stage='CANDIDATE_MATERIAL_BOUNDARY_BEFORE_THERMAL',
            boundary=boundary, traces=boundary['traces'], v_after_m_s=boundary['velocity'],
            work_velocity_m_s=boundary['work_velocity'], L_midpoint_per_s=Lmid,
            force_density_kg_m3=old_rho.copy(), mass_kg=pa.m.copy(), smoothing_length_m=pa.h.copy(),
            drho_kg_m3_s=drho, rho_mid_kg_m3=rho_mid, moment_condition_number=condition,
            plastic_work_J=plastic_j, return_dissipation_J=return_j,
            normal_midpoint=normal_snapshot(x, pa.particle_id, time, p), spatial=spatial)
        thermal = advance_thermal(pa.m, old_rho, old_t, p['conductivity'], x, i, j, g,
            dt, p['cp'], plastic_j, boundary['friction_raw_J'], boundary['contact_power_W'],
            boundary['wall_row_W_K'], p['plastic_heat_fraction'], p['heat_partition'], p['room'], p['melt'])
        pa.temperature[:] = thermal['temperature']
        self.last_thermal = thermal
        self.last_journal.update(stage='CANDIDATE_COUPLED_CALLBACK_COMPLETE', thermal=thermal,
                                 temperature_after_K=pa.temperature.copy())
        self._candidate_increment = self._increments(boundary['ledger'], thermal, plastic_j, return_j,
            float(dt*np.sum(volume[:, None, None]*sigma*L)), dt*stress_power)
        self.floor_evaluations += int(np.count_nonzero(rate < p['rate_floor']))

    def _increments(self, boundary, thermal, plastic, returned, old_work, midpoint_work):
        increment = dict.fromkeys(self.ledger, 0.)
        for key in ('wall_work_J', 'actuator_work_J', 'internal_work_J', 'delta_kinetic_J',
                    'friction_raw_J', 'normal_material_work_J', 'normal_wall_work_J', 'normal_relative_work_J'):
            increment[key] = float(boundary[key])
        for body in ('base', 'tool', 'actuator'):
            for k, axis in enumerate('xyz'):
                increment[f'{body}_impulse_{axis}_Ns'] = float(boundary[f'{body}_impulse_Ns'][k])
                increment[f'{body}_angular_impulse_{axis}_Nms'] = float(boundary[f'{body}_angular_impulse_Nms'][k])
        increment.update(plastic_work_J=float(plastic.sum()), numerical_return_dissipation_J=float(returned.sum()),
            material_old_work_J=old_work, midpoint_stress_work_J=midpoint_work)
        thermal_ledger = thermal['ledger']
        mapping = dict(plastic_heat_J='plastic_heat_J', nonthermal_plastic_J='plastic_unthermalized_J',
            friction_heat_J='friction_material_heat_J', friction_heat_to_wall_J='friction_wall_heat_J',
            contact_heat_J='contact_in_J', wall_heat_out_J='wall_heat_out_J',
            internal_conduction_net_J='internal_conduction_J', thermal_balance_residual_J='mc_delta_T_residual_J')
        for target, source in mapping.items(): increment[target] = float(thermal_ledger[source])
        return increment

    def _validate_candidate(self, normal_after, increment):
        pa, p = self.pa, self.p
        for key, array in pa.properties.items():
            if not np.isfinite(array.get_npy_array()).all(): raise ValueError('nonfinite coupled candidate '+key)
        if len(np.unique(pa.particle_id)) != len(pa.particle_id): raise ValueError('duplicate particle ID')
        if abs(float(pa.m.sum())-self.initial_mass) > 1e-12*self.initial_mass: raise ValueError('mass conservation limit')
        if np.max(np.abs(pa.rho/p['rho0']-1)) > p['density_limit']: raise ValueError('density change limit')
        if np.linalg.norm(velocities(pa), axis=1).max() > p['max_speed']: raise ValueError('speed limit')
        if np.any(pa.temperature < p['room']) or np.any(pa.temperature >= p['melt']): raise ValueError('thermal solid-domain limit')
        original = np.column_stack((pa.x0, pa.y0, pa.z0))
        if np.linalg.norm(positions(pa)-original, axis=1).max() > p['maximum_displacement_m']:
            raise ValueError('finite-feature displacement limit')
        if max(normal_after['max_penetration_m'], self.last_journal['normal_midpoint']['max_penetration_m']) > p['penetration_limit']*p['dx']:
            raise ValueError('contact penetration limit')
        d = self.last_journal['spatial']
        if d['condition_max'] > 1e6: raise ValueError('spatial condition limit')
        if abs(d['power_conjugacy_residual_W']) > 1e-12+1e-10*d['power_absolute_scale_W']:
            raise ValueError('spatial power conjugacy limit')
        if np.linalg.norm(d['internal_net_force_N']) > 1e-12+1e-10*d['force_absolute_scale_N']:
            raise ValueError('net internal force limit')
        if np.linalg.norm(d['internal_net_torque_Nm']) > 1e-14+1e-10*d['torque_absolute_scale_Nm']:
            raise ValueError('net internal torque limit')
        command = self.last_journal['boundary']['command']
        for trace in self.last_contact_trace:
            for interval in trace['result']['intervals']:
                if np.max(np.abs(trace['free_axes']*interval['constraint_impulse'])) > 1e-17:
                    raise ValueError('reaction impulse on a free degree of freedom')
        if np.any(command['active']) and np.max(np.abs(pa.w[command['active']]+p['feed'])) > 1e-12:
            raise ValueError('controlled axial velocity limit')
        residual = (increment['delta_kinetic_J']+increment['friction_raw_J']-increment['wall_work_J']-
            increment['actuator_work_J']-increment['internal_work_J']-increment['normal_relative_work_J'])
        scale = sum(abs(increment[k]) for k in ('delta_kinetic_J', 'friction_raw_J', 'wall_work_J',
                                               'actuator_work_J', 'internal_work_J', 'normal_relative_work_J'))
        if abs(residual) > 1e-15+1e-9*scale: raise ValueError('coupled kinetic/contact work identity limit')
        thermal = self.last_thermal['particle_energy_J']
        thermal_scale = sum(float(np.sum(np.abs(thermal[k]))) for k in
            ('plastic_heat_J', 'friction_material_heat_J', 'contact_in_J', 'internal_conduction_J'))
        if abs(increment['thermal_balance_residual_J']) > 1e-12+1e-9*thermal_scale:
            raise ValueError('coupled thermal source/temperature balance limit')

    def advance(self, max_dt=None):
        dt = self.stable_dt()
        if max_dt is not None:
            if not np.isfinite(max_dt) or max_dt <= 0: raise ValueError('positive finite max_dt required')
            dt = min(dt, max_dt)
        if not np.isfinite(dt) or dt <= 0: raise ValueError('invalid stable timestep')
        validate_interval(self.time, dt, self.p['release_time'])
        arrays = {k: a.get_npy_array().copy() for k, a in self.pa.properties.items()}
        keys = ('time', 'count', 'ledger', 'last', 'floor_evaluations', 'next_id',
                'last_contact_trace', 'last_journal', 'last_thermal')
        accepted = {key: copy.deepcopy(getattr(self, key)) for key in keys}
        before_x = positions(self.pa); before_normal = self._normal()
        self._candidate_normal_before = before_normal
        self.last_contact_trace = []; self.last_journal = None; self.last_thermal = None
        try:
            self.integrator.advance(self.time, dt)
            end = self.time+dt; after_normal = self._normal(end)
            middle = self.last_journal['normal_midpoint']; middle_x = self.last_journal['x_mid_m']
            transitions = dict(before_to_midpoint=self._feature_transition_guard(before_normal, middle,
                before_x, middle_x, self.time, self.time+dt/2),
                midpoint_to_after=self._feature_transition_guard(middle, after_normal, middle_x,
                    positions(self.pa), self.time+dt/2, end))
            increment = self._candidate_increment
            increment['normal_temporal_residual_J'] = (after_normal['normal_potential_J']-
                before_normal['normal_potential_J']+increment['normal_relative_work_J'])
            self.last_journal.update(normal_after=after_normal, feature_transitions=transitions,
                normal_temporal_residual_J=increment['normal_temporal_residual_J'],
                x_before_m=before_x, x_after_m=positions(self.pa), rho_after_kg_m3=self.pa.rho.copy())
            self._validate_candidate(after_normal, increment)
            for key in self.ledger: self.ledger[key] += increment[key]
            self.time = end; self.count += 1
            crossed = ((self.pa.component_id == 1) & (self.pa.z0 >= self.p['layer']) & (self.pa.z < self.p['layer']))
            self.pa.feed_crossed[crossed] = 1.
            self.last = dict(dt_s=dt, **self.last_journal['spatial'],
                contacts=len(after_normal['active_contact_ids']),
                plastic_increment_J=increment['plastic_work_J'],
                normal_temporal_residual_J=increment['normal_temporal_residual_J'])
            self.last_journal.update(stage='ACCEPTED_COUPLED_DRIFT_KICK_DRIFT', increments=increment,
                                     accepted_time_s=self.time)
        except Exception as error:
            self.rejected_candidate = dict(classification='REJECTED', accepted=False, restart_allowed=False,
                error=str(error), time_s=accepted['time']+dt, step=accepted['count']+1,
                source_hashes=self.source_hashes, configuration_hash=self.configuration_hash,
                properties={k: a.get_npy_array().copy() for k, a in self.pa.properties.items()},
                journal=copy.deepcopy(self.last_journal), traces=copy.deepcopy(self.last_contact_trace),
                boundary_context=copy.deepcopy(getattr(error, 'contact_context', None)),
                feature_context=copy.deepcopy(getattr(error, 'feature_context', None)),
                last_accepted_properties=arrays, last_accepted_state=accepted)
            for key, value in arrays.items(): getattr(self.pa, key)[:] = value
            for key, value in accepted.items(): setattr(self, key, value)
            raise
        return dt

    def energy(self):
        pa, p = self.pa, self.p; s = pa.S.reshape(-1, 3, 3)
        normal = normal_snapshot(positions(pa), pa.particle_id, self.time, p)
        result = dict(kinetic_J=float(.5*np.sum(pa.m[:, None]*velocities(pa)**2)),
            elastic_shear_J=float(np.sum(pa.m/pa.rho*np.sum(s*s, axis=(1, 2))/(4*p['shear']))),
            elastic_volume_J=float(np.sum(pa.m*volumetric_energy(pa.rho, p['rho0'], self.bulk))),
            thermal_increment_J=float(np.sum(pa.m*p['cp']*(pa.temperature-p['temperature']))),
            contact_potential_J=normal['normal_potential_J'],
            nonthermal_plastic_J=float((1-p['plastic_heat_fraction'])*pa.plastic_work.sum()))
        result['total_J'] = sum(result.values())
        return result

    def diagnostics(self, fail=False):
        """Read-only diagnostics; optional state safety check, energy only warns.

        The full coupled candidate checks always run inside atomic advance().
        fail=True additionally checks the current state without another kick,
        material update or NNPS evaluation.
        """
        pa, p = self.pa, self.p
        if fail:
            if any(not np.isfinite(a.get_npy_array()).all() for a in pa.properties.values()):
                raise ValueError('nonfinite current coupled state')
            if abs(float(pa.m.sum())-self.initial_mass) > 1e-12*self.initial_mass: raise ValueError('mass conservation limit')
            if np.max(np.abs(pa.rho/p['rho0']-1)) > p['density_limit']: raise ValueError('density change limit')
            if np.linalg.norm(velocities(pa), axis=1).max() > p['max_speed']: raise ValueError('speed limit')
            if np.any(pa.temperature < p['room']) or np.any(pa.temperature >= p['melt']): raise ValueError('thermal solid-domain limit')
            if self._normal()['max_penetration_m'] > p['penetration_limit']*p['dx']: raise ValueError('contact penetration limit')
        energy = self.energy(); ledger = self.ledger
        residual = energy['total_J']-self.initial_energy-ledger['wall_work_J']-ledger['actuator_work_J']+ledger['wall_heat_out_J']
        scale = max(abs(self.initial_energy), abs(ledger['wall_work_J'])+abs(ledger['actuator_work_J'])+abs(ledger['wall_heat_out_J']), 1e-15)
        warning = abs(residual) > p['energy_absolute_limit']+p['energy_relative_limit']*scale
        normal = self._normal(); command = control(pa.z0, self.time, 0., p['bar_length'], p['actuator_depth'], p['feed'], p['release_time'])
        actual_p = np.sum(pa.m[:, None]*velocities(pa), axis=0)
        actual_l = np.sum(pa.m[:, None]*np.cross(positions(pa), velocities(pa)), axis=0)
        impulse = sum((np.array([ledger[f'{body}_impulse_{a}_Ns'] for a in 'xyz']) for body in ('base', 'tool', 'actuator')), np.zeros(3))
        angular = sum((np.array([ledger[f'{body}_angular_impulse_{a}_Nms'] for a in 'xyz']) for body in ('base', 'tool', 'actuator')), np.zeros(3))
        original = np.column_stack((pa.x0, pa.y0, pa.z0)); displacement = np.linalg.norm(positions(pa)-original, axis=1)
        bar = pa.component_id == 1; substrate = pa.component_id == 0
        return dict(**energy, cumulative=copy.deepcopy(ledger), time_s=self.time, step=self.count,
            initial_kinetic_J=self.initial_kinetic_J, initial_normal_potential_J=self.initial_normal_potential_J,
            initial_energy_J=self.initial_energy, initial_linear_momentum_kg_m_s=self.initial_linear_momentum.copy(),
            initial_angular_momentum_kg_m2_s=self.initial_angular_momentum.copy(),
            linear_momentum_kg_m_s=actual_p, angular_momentum_kg_m2_s=actual_l,
            linear_momentum_residual_Ns=actual_p-self.initial_linear_momentum-impulse,
            angular_momentum_residual_Nms=actual_l-self.initial_angular_momentum-angular,
            mass_kg=float(pa.m.sum()), density_relative_change_max=float(np.max(np.abs(pa.rho/p['rho0']-1))),
            speed_max_m_s=float(np.linalg.norm(velocities(pa), axis=1).max()),
            displacement_max_m=float(displacement.max()), temperature_min_K=float(pa.temperature.min()),
            temperature_max_K=float(pa.temperature.max()), ep_max=float(pa.ep.max()),
            plastic_particles=int(np.count_nonzero(pa.ep > 0)),
            substrate_plastic_particles=int(np.count_nonzero((pa.ep > 0)&substrate)),
            bar_plastic_particles=int(np.count_nonzero((pa.ep > 0)&bar)),
            substrate_displacement_max_m=float(np.max(displacement[substrate], initial=0.)),
            feed_plane_crossed_mass_kg=float(pa.m[pa.feed_crossed == 1].sum()),
            energy_residual_J=residual, energy_residual_relative=residual/scale,
            energy_status='ENERGY_WARNING' if warning else 'DIAGNOSTIC_WITHIN_SCREEN_NOT_FREE_ENERGY_CERTIFICATION',
            max_penetration_m=normal['max_penetration_m'], contacts=len(normal['active_contact_ids']),
            contact_feature_ids=sorted(set(int(a) for a in normal['rows']['feature_id'])),
            actuator_active_count=int(np.count_nonzero(command['active'])),
            rate_floor_evaluations=self.floor_evaluations, classification=CLASSIFICATION,
            production_authorized=False, material_energy_certified=False,
            mechanism_scope='FINITE_PREBONDED_HOT_START; NO_CONTINUOUS_INFLOW_OR_DEPOSITION_CERTIFICATION')

    def record(self):
        pa = self.pa; s = pa.S.reshape(-1, 3, 3)
        pr = pressure(pa.rho, self.p['rho0'], self.bulk)
        origin, angle = tool_pose(self.time, self.p['traverse'], self.p['omega'], self.p['layer'])
        return dict(**self.diagnostics(False), particle_id=pa.particle_id.copy(),
            component_id=pa.component_id.copy(), material_id=pa.material_id.copy(),
            tracer_type=pa.tracer_type.copy(), wholebar_tag=pa.wholebar_tag.copy(),
            actuator=pa.actuator.copy(), as_rs_tag=pa.as_rs_tag.copy() if 'as_rs_tag' in pa.properties else None,
            position_m=positions(pa), velocity_m_s=velocities(pa), temperature_K=pa.temperature.copy(),
            density_kg_m3=pa.rho.copy(), pressure_Pa=pr, deviatoric_stress_Pa=s.copy(),
            stress_Pa=s-pr[:, None, None]*np.eye(3), equivalent_stress_Pa=equivalent(s),
            ep=pa.ep.copy(), strain_rate_per_s=pa.strain_rate.copy() if self.count else None,
            plastic_work_J=pa.plastic_work.copy(), return_dissipation_J=pa.return_dissipation.copy(),
            particle_mass_kg=pa.m.copy(), tool_origin_m=origin, tool_angle_rad=angle,
            spatial=self.last_journal['spatial'] if self.last_journal is not None else None)

    def state(self):
        return dict(classification=CLASSIFICATION, accepted=True, restart_allowed=True,
            production_authorized=False, t=self.time, count=self.count,
            ledger=copy.deepcopy(self.ledger), last=copy.deepcopy(self.last),
            initial_energy=self.initial_energy, initial_mass=self.initial_mass,
            initial_kinetic_J=self.initial_kinetic_J, initial_normal_potential_J=self.initial_normal_potential_J,
            initial_linear_momentum_kg_m_s=self.initial_linear_momentum.copy(),
            initial_angular_momentum_kg_m2_s=self.initial_angular_momentum.copy(),
            initial_particle_ids=self.initial_particle_ids,
            floor_evaluations=self.floor_evaluations, next_id=self.next_id,
            parameters=copy.deepcopy(self.p), parameter_fingerprint=self.parameter_fingerprint,
            configuration_hash=self.configuration_hash, source_hashes=self.source_hashes)

    def checkpoint(self, path):
        save_checkpoint(path, self.pa, self.state())

    @classmethod
    def from_checkpoint(cls, path, p, source_hashes=None, configuration_hash=None):
        pa, state = restore(path)
        return cls(pa, p, state=state, source_hashes=source_hashes, configuration_hash=configuration_hash)

    def vtk(self, path):
        p = self.pa; r = self.record()
        scalar = {k: getattr(p, k).tolist() for k in ('particle_id', 'material_id', 'component_id', 'tracer_type', 'wholebar_tag')}
        scalar.update(temperature_K=p.temperature.tolist(), density_kg_per_m3=p.rho.tolist(),
            pressure_Pa=r['pressure_Pa'].tolist(), equivalent_stress_Pa=r['equivalent_stress_Pa'].tolist(),
            equivalent_plastic_strain=p.ep.tolist(), mass_kg=p.m.tolist(),
            plastic_work_J=p.plastic_work.tolist(), numerical_return_dissipation_J=p.return_dissipation.tolist())
        for a, b, name in ((0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'), (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')):
            scalar['stress_'+name+'_Pa'] = r['stress_Pa'][:, a, b].tolist()
        if self.count:
            scalar['equivalent_deviatoric_strain_rate_per_s'] = p.strain_rate.tolist()
            scalar['equivalent_plastic_strain_rate_per_s'] = p.ep_rate.tolist()
        write_point_cloud(path, positions(p).tolist(), scalar,
            dict(velocity_m_per_s=velocities(p).tolist(), initial_position_m=np.column_stack((p.x0, p.y0, p.z0)).tolist()),
            title=f'{CLASSIFICATION} time_s={self.time:.17g}; COUPLED_FIELDS_NOT_DEPOSITION_VALIDATION')

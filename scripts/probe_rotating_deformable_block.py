"""Bounded isothermal 3D SPH rotating-shoulder block, NOT AFSD deposition.

Companion only: all existing material, corrected spatial, contact and upstream
integrator implementations remain unchanged. Prescribed normal-load/clamp
fixture; no normal-opening model and no solved temperature field.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.integrator import HostLeapfrog
from afsd_backend.pysph.mechanics import (positions, velocities, corrected_operators,
    internal_force, pressure, bulk_modulus, volumetric_energy)
from afsd_backend.pysph.material import FKS, update, equivalent
from afsd_backend.pysph.contact_kick import kick
from afsd_backend.pysph.rigid_contact_kinematics import contact_rows
from afsd_backend.pysph.output import checkpoint as save_checkpoint, restore
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG_PATH = ROOT/'configs/benchmarks/rotating_deformable_block.json'
INPUT_CONFIG = json.loads(CONFIG_PATH.read_text(encoding='utf-8'))
PARAMETERS = {k: v['value'] for k, v in INPUT_CONFIG['parameters'].items()}
CLASSIFICATION = INPUT_CONFIG['classification']
SOURCES = ['scripts/probe_rotating_deformable_block.py',
    'configs/benchmarks/rotating_deformable_block.json',
    'docs/ROTATING_DEFORMABLE_CONTRACT.md',
    'src/afsd_backend/pysph/particles.py', 'src/afsd_backend/pysph/kernel.py',
    'src/afsd_backend/pysph/mechanics.py', 'src/afsd_backend/pysph/material.py',
    'src/afsd_backend/pysph/geometry.py', 'src/afsd_backend/pysph/rigid_contact_kinematics.py',
    'src/afsd_backend/pysph/contact_kick.py', 'src/afsd_backend/pysph/contact_cone.py',
    'src/afsd_backend/pysph/integrator.py', 'src/afsd_backend/pysph/output.py',
    'src/afsd_scaffold/vtk_legacy.py',
    'third_party/PySPH/pysph/sph/integrator.py',
    'third_party/PySPH/pysph/sph/integrator_step.py',
    'third_party/PySPH/pysph/base/nnps.py',
    'third_party/PySPH/pysph/base/nnps_base.pyx',
    'third_party/PySPH/pysph/base/linked_list_nnps.pyx']


def jsonable(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [jsonable(v) for v in value]
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(jsonable(value), indent=2, allow_nan=False)+'\n', encoding='utf-8')


def rejected_jsonable(value):
    """Preserve nonfinite rejection evidence with explicit tokens, never clip it."""
    if isinstance(value, np.ndarray): return rejected_jsonable(value.tolist())
    if isinstance(value, np.generic): return rejected_jsonable(value.item())
    if isinstance(value, dict): return {k: rejected_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [rejected_jsonable(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return {'nonfinite_float': 'NaN' if np.isnan(value) else 'Infinity' if value > 0 else '-Infinity'}
    return value


def hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def config_hash(config):
    return hashlib.sha256(json.dumps(config, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def initial_positions():
    p = PARAMETERS
    return np.array(np.meshgrid(.01+(np.arange(4)-1.5)*p['dx'],
        (np.arange(4)-1.5)*p['dx'], (np.arange(4)+.5)*p['dx'], indexing='ij')).reshape(3, -1).T


def configuration(dt):
    return dict(classification=CLASSIFICATION, source_configuration=INPUT_CONFIG,
        parameters=PARAMETERS, dt_s=dt, tags_are_observation_only=True,
        physical_scaling='NONE', spatial_operator='corrected_work_conjugate',
        normal_law='PRESCRIBED 1 N PER TOP PARTICLE; NORMAL DOF CLAMPED; NOT PENALTY',
        field_status=dict(position='SOLVED', velocity='SOLVED', density='SOLVED_CONTINUITY',
            stress='SOLVED_EXISTING_JAUMANN_J2_FKS_APPROXIMATE', ep='SOLVED',
            strain_rate='SOLVED_PRECONTACT_TIME_LAYER', temperature='PRESCRIBED_INPUT_NOT_SOLVED',
            friction_dissipation='SOLVED_RAW_NOT_DEPOSITED_AS_HEAT',
            Es='DIAGNOSTIC_ONLY_NOT_CERTIFIED_FREE_ENERGY',
            checkpoint_e='UPSTREAM_SCRATCH_NOT_SOLVED_NOT_PHYSICAL_ENERGY'),
        production_authorized=False)


def frozen_gradient(v, mass, rho, i, j, g, B):
    """Diagnostic only; same frozen spatial operator, arbitrary supplied velocity."""
    raw = np.zeros((len(v), 3, 3))
    outer = (v[j]-v[i])[:, :, None]*g[:, None, :]
    np.add.at(raw, i, (mass[j]/rho[j])[:, None, None]*outer)
    np.add.at(raw, j, (mass[i]/rho[i])[:, None, None]*outer)
    return raw@B


class Block:
    def __init__(self, dt, pa=None, state=None):
        if dt not in INPUT_CONFIG['dt_cases_s']:
            raise ValueError('only preregistered dt=5e-9 or 2.5e-9 seconds')
        self.dt = float(dt); self.config = configuration(self.dt); self.source_hashes = hashes()
        self.p = dict(PARAMETERS)
        self.bulk = bulk_modulus(self.p['shear'], self.p['poisson'])
        self.sound = float(np.sqrt((self.bulk+4*self.p['shear']/3)/self.p['rho0']))
        self.law = FKS(**{k: self.p[k] for k in FKS.__dataclass_fields__})
        if pa is None:
            x = initial_positions()
            pa = make_array(x, np.zeros(64, dtype=int), self.p)
            role = np.where(x[:, 2] == x[:, 2].min(), 1,
                            np.where(x[:, 2] == x[:, 2].max(), 2, 0))
            pa.add_property('boundary_role', type='long', data=role)
            pa.add_property('as_rs_tag', type='long', data=np.zeros(64, dtype=int))
            pa.tracer_type[:] = pa.particle_id % 3
            pa.set_output_arrays(list(pa.properties))
        self.pa = pa; self.time = 0.; self.step = 0
        self.initial_mass = float(64*self.p['rho0']*self.p['dx']**3)
        self.last_traces = []; self.last_journal = None
        self.cumulative = dict(friction_raw_J=0., wall_work_J=0., constraint_work_J=0.,
            internal_work_J=0., normal_work_J=0., delta_kinetic_J=0.,
            contact_kinetic_residual_J=0., material_old_work_J=0.,
            midpoint_stress_work_J=0., plastic_work_J=0., return_dissipation_J=0.)
        if state is not None:
            if not state.get('accepted') or not state.get('restart_allowed'):
                raise ValueError('only accepted component checkpoint can be restored')
            if state['configuration_sha256'] != config_hash(self.config):
                raise ValueError('checkpoint configuration mismatch')
            if state['source_sha256'] != self.source_hashes:
                raise ValueError('checkpoint source mismatch')
            if state['classification'] != CLASSIFICATION:
                raise ValueError('checkpoint classification mismatch')
            self.time = float(state['time_s']); self.step = int(state['step'])
            self.cumulative = dict(state['cumulative'])
        self.integrator = HostLeapfrog(self.pa, self._callback)

    @classmethod
    def from_checkpoint(cls, path):
        pa, state = restore(path)
        return cls(state['dt_s'], pa, state)

    def stable_dt(self):
        return float(self.p['cfl']*self.p['h_ratio']*self.p['dx']/
                     (self.sound+np.linalg.norm(velocities(self.pa), axis=1).max()))

    def _callback(self, midpoint_time, dt):
        p = self.pa; x = positions(p); v = velocities(p); rho = p.rho.copy()
        i, j, g, L, drho, B = corrected_operators(p)
        old_s = p.S.reshape(-1, 3, 3).copy()
        new_s, ep, plastic_density, return_density, rate = update(
            old_s, p.ep.copy(), L, p.temperature, dt, self.p['shear'], self.law)
        rho_mid = rho+.5*dt*drho
        sigma = .5*(old_s+new_s)-pressure(rho_mid, self.p['rho0'], self.bulk)[:, None, None]*np.eye(3)
        force = internal_force(p, sigma, i, j, g, correction=B)
        rows = contact_rows(x, p.particle_id, midpoint_time, self.p)
        top = p.boundary_role == 2
        if (len(rows['particle_id']) != 16 or np.any(rows['feature_id'] != 1) or
                set(rows['particle_id']) != set(p.particle_id[top])):
            raise ValueError('fixture requires exactly 16 persistent top shoulder contacts')
        after = np.empty_like(v); work_v = np.empty_like(v); traces = []
        self.last_traces = traces
        self.last_journal = dict(stage='CANDIDATE_MATERIAL_BEFORE_CONTACT',
            step=self.step+1, time_mid_s=midpoint_time, dt_s=dt,
            particle_id=p.particle_id.copy(), x_mid_m=x.copy(), v_before_m_s=v.copy(),
            rho_before_kg_m3=rho, rho_mid_kg_m3=rho_mid, S_before_Pa=old_s,
            S_after_Pa=new_s.copy(), ep_before=p.ep.copy(), ep_after=ep.copy(),
            L_before_per_s=L, sigma_mid_Pa=sigma, internal_force_N=force,
            i=i, j=j, g=g, B=B, traces=traces)
        for idx, pid in enumerate(p.particle_id):
            selected = np.flatnonzero(rows['particle_index'] == idx)
            faces = [{key: value if isinstance(value, str) else np.asarray(value)[r].copy()
                      for key, value in rows.items()} for r in selected]
            free = ([0., 0., 0.] if p.boundary_role[idx] == 1 else
                    [1., 1., 0.] if p.boundary_role[idx] == 2 else [1., 1., 1.])
            result = kick(v[idx], p.m[idx], force[idx], rows['normal'][selected],
                np.full(len(selected), self.p['normal_force']), rows['wall_velocity'][selected],
                self.p['friction'], dt, free_axes=free, prescribed_velocity=np.zeros(3))
            after[idx] = result['velocity']
            # This is the actual impulse-work average, not an assumed endpoint midpoint.
            work_v[idx] = sum(a['dt_s']*.5*(a['v_before']+a['v_after'])
                             for a in result['intervals'])/dt
            traces.append(dict(particle_id=int(pid), particle_index=idx,
                stage='MIDPOINT_FROZEN_GEOMETRY_KICK', step=self.step+1,
                time_mid_s=midpoint_time, mass_kg=float(p.m[idx]), free_axes=free,
                position_m=x[idx].copy(), material_velocity_before_m_s=v[idx].copy(),
                contact_rows=faces, internal_force_N=force[idx].copy(), result=result))
        Lmid = frozen_gradient(work_v, p.m, rho, i, j, g, B)
        volume = p.m/rho
        power_force = float(np.sum(force*work_v))
        power_stress = float(np.sum(volume[:, None, None]*sigma*Lmid))
        stress_old = float(dt*np.sum(volume[:, None, None]*sigma*L))
        net_force = force.sum(axis=0); net_torque = np.cross(x, force).sum(axis=0)
        condition = np.linalg.cond(np.linalg.inv(B))
        diagnostics = dict(pair_count=len(i), condition_max=float(condition.max()),
            internal_power_W=power_force, midpoint_stress_power_W=power_stress,
            power_conjugacy_residual_W=power_force+power_stress,
            internal_net_force_N=net_force, internal_net_torque_Nm=net_torque,
            force_absolute_scale_N=float(np.linalg.norm(force, axis=1).sum()),
            torque_absolute_scale_Nm=float(np.linalg.norm(np.cross(x, force), axis=1).sum()),
            power_absolute_scale_W=float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(work_v, axis=1))),
            material_old_work_J=stress_old, midpoint_stress_work_J=dt*power_stress,
            temporal_stress_work_split_J=dt*power_stress-stress_old)
        self.last_journal = dict(step=self.step+1, time_mid_s=midpoint_time, dt_s=dt,
            stage='CANDIDATE_MATERIAL_AND_CONTACT_COMPLETE',
            particle_id=p.particle_id.copy(), boundary_role=p.boundary_role.copy(),
            x_mid_m=x.copy(), v_before_m_s=v.copy(), v_after_m_s=after.copy(),
            work_velocity_m_s=work_v, endpoint_midpoint_velocity_m_s=.5*(v+after),
            mass_kg=p.m.copy(), rho_before_kg_m3=rho, force_density_kg_m3=rho.copy(),
            smoothing_length_m=p.h.copy(), moment_condition_number=condition,
            i=i, j=j, g=g, B=B,
            L_before_per_s=L, L_midpoint_per_s=Lmid, sigma_mid_Pa=sigma,
            S_before_Pa=old_s, S_after_Pa=new_s.copy(), ep_before=p.ep.copy(), ep_after=ep.copy(),
            internal_force_N=force, drho_kg_m3_s=drho, rho_mid_kg_m3=rho_mid,
            plastic_work_J=volume*plastic_density, return_dissipation_J=volume*return_density,
            traces=traces, diagnostics=diagnostics)
        p.au[:], p.av[:], p.aw[:] = ((after-v)/dt).T
        p.arho[:] = drho
        p.ep_rate[:] = (ep-p.ep)/dt; p.S[:] = new_s.ravel(); p.ep[:] = ep
        p.strain_rate[:] = rate; p.D[:] = (.5*(L+L.swapaxes(1, 2))).ravel()
        p.plastic_work[:] += volume*plastic_density
        p.return_dissipation[:] += volume*return_density

    def _validate_candidate(self, increment):
        p = self.pa
        for key, array in p.properties.items():
            if not np.isfinite(array.get_npy_array()).all():
                raise ValueError('nonfinite candidate property '+key)
        if abs(float(p.m.sum())-self.initial_mass) > 1e-12*self.initial_mass:
            raise ValueError('mass conservation limit')
        if np.max(np.abs(p.rho/self.p['rho0']-1)) > .01: raise ValueError('density change limit')
        if np.linalg.norm(velocities(p), axis=1).max() > 1.: raise ValueError('speed limit')
        original = np.column_stack((p.x0, p.y0, p.z0))
        if np.linalg.norm(positions(p)-original, axis=1).max() > 1e-5:
            raise ValueError('displacement limit')
        fixed = np.zeros((len(p.x), 3), dtype=bool)
        fixed[p.boundary_role == 1] = True; fixed[p.boundary_role == 2, 2] = True
        if (np.max(np.abs(velocities(p)[fixed])) > 1e-12 or
                np.max(np.abs((positions(p)-original)[fixed])) > 2e-17):
            raise ValueError('fixed-DOF motion limit')
        diag = self.last_journal['diagnostics']
        if diag['condition_max'] > 1e6: raise ValueError('condition limit')
        if abs(diag['power_conjugacy_residual_W']) > 1e-12+1e-10*diag['power_absolute_scale_W']:
            raise ValueError('spatial power conjugacy limit')
        if np.linalg.norm(diag['internal_net_force_N']) > 1e-12+1e-10*diag['force_absolute_scale_N']:
            raise ValueError('net internal force limit')
        if np.linalg.norm(diag['internal_net_torque_Nm']) > 1e-14+1e-10*diag['torque_absolute_scale_Nm']:
            raise ValueError('net internal torque limit')
        residual = (increment['delta_kinetic_J']+increment['friction_raw_J']-
            increment['wall_work_J']-increment['constraint_work_J']-increment['internal_work_J'])
        scale = sum(abs(increment[k]) for k in ('delta_kinetic_J', 'friction_raw_J',
                    'wall_work_J', 'constraint_work_J', 'internal_work_J'))
        if abs(residual) > 1e-15+1e-9*scale:
            raise ValueError('kinetic/contact work identity limit')

    def advance(self):
        if self.time+self.dt > self.p['end_time']+1e-20: raise ValueError('beyond frozen end time')
        if self.dt > self.stable_dt(): raise ValueError('preregistered acoustic timestep bound exceeded')
        before = {k: a.get_npy_array().copy() for k, a in self.pa.properties.items()}
        old_traces, old_journal = self.last_traces, self.last_journal
        old_cumulative = dict(self.cumulative)
        self.last_traces = []; self.last_journal = None
        try:
            self.integrator.advance(self.time, self.dt)
            increment = dict.fromkeys(self.cumulative, 0.)
            for trace in self.last_traces:
                result = trace['result']
                increment['friction_raw_J'] += float(result['friction_heat_J'].sum())
                for key in ('wall_work_J', 'constraint_work_J', 'delta_kinetic_J'):
                    increment[key] += result[key]
                increment['internal_work_J'] += result['other_work_J']
                increment['normal_work_J'] += sum(float(a['normal_work_J'].sum()) for a in result['intervals'])
                increment['contact_kinetic_residual_J'] += result['kinetic_work_residual_J']
            for key in ('material_old_work_J', 'midpoint_stress_work_J'):
                increment[key] = self.last_journal['diagnostics'][key]
            increment['plastic_work_J'] = float(self.last_journal['plastic_work_J'].sum())
            increment['return_dissipation_J'] = float(self.last_journal['return_dissipation_J'].sum())
            self._validate_candidate(increment)
        except Exception as error:
            self.rejected_candidate = dict(classification='REJECTED', accepted=False,
                restart_allowed=False, production_authorized=False, error=str(error),
                source_sha256=self.source_hashes, configuration_sha256=config_hash(self.config),
                time_s=self.time+self.dt, step=self.step+1,
                properties={k: a.get_npy_array().copy() for k, a in self.pa.properties.items()},
                journal=copy.deepcopy(self.last_journal), traces=copy.deepcopy(self.last_traces),
                last_accepted_time_s=self.time, last_accepted_step=self.step,
                last_accepted_properties=before, last_accepted_cumulative=old_cumulative)
            for key, values in before.items(): getattr(self.pa, key)[:] = values
            self.cumulative = old_cumulative; self.last_traces = old_traces; self.last_journal = old_journal
            raise
        for key in self.cumulative: self.cumulative[key] += increment[key]
        self.last_journal['increments'] = increment
        self.last_journal['stage'] = 'ACCEPTED_DRIFT_KICK_DRIFT_COMPLETE'
        self.last_journal['x_before_m'] = np.column_stack((before['x'], before['y'], before['z']))
        self.last_journal['x_after_m'] = positions(self.pa)
        self.last_journal['rho_after_kg_m3'] = self.pa.rho.copy()
        self.step += 1; self.time = self.step*self.dt
        return self.record()

    def record(self):
        p = self.pa; s = p.S.reshape(-1, 3, 3)
        pr = pressure(p.rho, self.p['rho0'], self.bulk)
        kinetic = float(.5*np.sum(p.m[:, None]*velocities(p)**2))
        shear = float(np.sum(p.m/p.rho*np.sum(s*s, axis=(1, 2))/(4*self.p['shear'])))
        volume = float(np.sum(p.m*volumetric_energy(p.rho, self.p['rho0'], self.bulk)))
        c = self.cumulative
        residual = (kinetic+shear+volume+c['plastic_work_J']+c['friction_raw_J']-
                    c['wall_work_J']-c['constraint_work_J'])
        warning = abs(residual) > 1e-12+.01*(abs(c['wall_work_J'])+abs(c['constraint_work_J'])+c['friction_raw_J'])
        return dict(step=self.step, time_s=self.time, particle_id=p.particle_id.copy(),
            boundary_role=p.boundary_role.copy(), tracer_type=p.tracer_type.copy(), as_rs_tag=p.as_rs_tag.copy(),
            position_m=positions(p), velocity_m_s=velocities(p), mass_kg=p.m.copy(),
            density_kg_m3=p.rho.copy(), pressure_Pa=pr, deviatoric_stress_Pa=s.copy(),
            stress_Pa=s-pr[:, None, None]*np.eye(3), ep=p.ep.copy(),
            strain_rate_per_s=p.strain_rate.copy() if self.step else None,
            ep_rate_per_s=p.ep_rate.copy() if self.step else None,
            rate_field_status='PRECONTACT_TIME_LAYER' if self.step else 'NOT_EVALUATED',
            equivalent_stress_Pa=equivalent(s), kinetic_J=kinetic,
            Es_diagnostic_J=shear+volume, shear_energy_diagnostic_J=shear,
            volume_energy_diagnostic_J=volume, total_mass_kg=float(p.m.sum()),
            density_relative_change_max=float(np.max(np.abs(p.rho/self.p['rho0']-1))),
            speed_max_m_s=float(np.linalg.norm(velocities(p), axis=1).max()),
            tool_origin_m=[self.p['traverse']*self.time, 0., self.p['layer']],
            tool_angle_rad=self.p['omega']*self.time, cumulative=dict(c),
            kinetic_contact_residual_J=kinetic+c['friction_raw_J']-c['wall_work_J']-
                c['constraint_work_J']-c['internal_work_J'],
            energy_diagnostic_residual_J=residual,
            energy_status='ENERGY_WARNING' if warning else 'DIAGNOSTIC_WITHIN_SCREEN_NOT_CERTIFIED_FREE_ENERGY',
            diagnostic_status='OLD_JAUMANN_CLOSED_CYCLE_FAIL_REMAINS; ISOTHERMAL_NO_HEAT_SOLUTION',
            spatial=self.last_journal['diagnostics'] if self.last_journal is not None else None)

    def checkpoint(self, path):
        save_checkpoint(path, self.pa, dict(classification=CLASSIFICATION, accepted=True,
            restart_allowed=True, production_authorized=False, dt_s=self.dt, step=self.step,
            time_s=self.time, cumulative=self.cumulative, configuration=self.config,
            configuration_sha256=config_hash(self.config), source_sha256=self.source_hashes))

    def vtk(self, path):
        p = self.pa; r = self.record()
        scalar = {k: getattr(p, k).tolist() for k in ('particle_id', 'material_id', 'component_id',
            'tracer_type', 'wholebar_tag', 'actuator', 'boundary_role', 'as_rs_tag')}
        scalar.update(mass_kg=p.m.tolist(), density_kg_per_m3=p.rho.tolist(),
            pressure_Pa=r['pressure_Pa'].tolist(), equivalent_stress_Pa=r['equivalent_stress_Pa'].tolist(),
            equivalent_plastic_strain=p.ep.tolist(), smoothing_length_m=p.h.tolist())
        for a, b, name in ((0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'),
                           (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')):
            scalar['stress_'+name+'_Pa'] = r['stress_Pa'][:, a, b].tolist()
        if self.step:
            scalar['equivalent_deviatoric_strain_rate_per_s'] = p.strain_rate.tolist()
            scalar['equivalent_plastic_strain_rate_per_s'] = p.ep_rate.tolist()
        write_point_cloud(path, positions(p).tolist(), scalar,
            dict(velocity_m_per_s=velocities(p).tolist(),
                 initial_position_m=np.column_stack((p.x0, p.y0, p.z0)).tolist()),
            title=f'{CLASSIFICATION} time_s={self.time:.17g}; ISOTHERMAL_INPUT_ONLY_NO_SOLVED_TEMPERATURE')


def run(dt, output, probe_steps=None, timeout_s=None):
    output = Path(output)
    if output.exists(): raise ValueError('fresh output directory required; no historical overwrite')
    if probe_steps is not None and probe_steps != 2: raise ValueError('only two-step cost probe is registered')
    timeout_s = INPUT_CONFIG['run_timeout_s'] if timeout_s is None else timeout_s
    if not 0 < timeout_s <= INPUT_CONFIG['run_timeout_s']: raise ValueError('bounded positive wall budget required')
    output.mkdir(parents=True)
    start = time.perf_counter(); block = Block(dt)
    target_steps = 2 if probe_steps is not None else int(round(PARAMETERS['end_time']/dt))
    manifest = dict(classification=CLASSIFICATION, status='RUNNING', configuration=block.config,
        configuration_sha256=config_hash(block.config), source_sha256=block.source_hashes,
        maximum_wall_seconds=timeout_s, production_authorized=False, target_steps=target_steps,
        requested_end_time_s=target_steps*dt, mode='TWO_STEP_COST_PROBE' if probe_steps else 'FULL_REGISTERED_COMPONENT',
        initial_stable_dt_s=block.stable_dt(), sound_speed_m_s=block.sound,
        energy_note='Es is a diagnostic, not thermodynamic energy certification; raw Q not thermal deposition.')
    write_json(output/'manifest.json', manifest)
    records = [block.record()]; block.checkpoint(output/'last_accepted.npz'); block.vtk(output/'frame_0000.vtk')
    step_times = []
    try:
        with (output/'steps.jsonl').open('w', encoding='utf-8') as journal:
            while block.step < target_steps:
                if time.perf_counter()-start > timeout_s: raise TimeoutError('component wall budget exceeded')
                step_start = time.perf_counter(); records.append(block.advance())
                step_times.append(time.perf_counter()-step_start)
                journal.write(json.dumps(jsonable(block.last_journal), allow_nan=False)+'\n'); journal.flush()
                block.checkpoint(output/'last_accepted.npz'); block.vtk(output/f'frame_{block.step:04d}.vtk')
        manifest['status'] = 'COMPLETE_COST_PROBE_NOT_AFSD' if probe_steps else 'COMPLETE_DEFORMABLE_COMPONENT_NOT_AFSD'
    except Exception as error:
        manifest.update(status='BLOCKED', error=str(error))
        if hasattr(block, 'rejected_candidate'):
            write_json(output/'rejected_candidate.json', rejected_jsonable(block.rejected_candidate))
    finally:
        write_json(output/'history.json', records)
        manifest.update(wall_seconds=time.perf_counter()-start, accepted_steps=block.step,
            accepted_time_s=block.time, final_record=block.record(), step_compute_seconds=step_times,
            particle_array_bytes=sum(a.get_npy_array().nbytes for a in block.pa.properties.values()),
            memory_note='ParticleArray property bytes only; NOT peak process RSS.')
        manifest['artifact_sha256'] = {str(path.relative_to(output)).replace('\\', '/'):
            hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(output.rglob('*'))
            if path.is_file() and path != output/'manifest.json'}
        write_json(output/'manifest.json', manifest)
    print(json.dumps(jsonable({k: v for k, v in manifest.items() if k not in
        ('final_record', 'configuration', 'artifact_sha256')}), allow_nan=False), flush=True)
    return 0 if manifest['status'].startswith('COMPLETE') else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dt', type=float, choices=INPUT_CONFIG['dt_cases_s'], default=5e-9)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--probe-steps', type=int, choices=[2])
    args = parser.parse_args()
    return run(args.dt, args.output, args.probe_steps)


if __name__ == '__main__': raise SystemExit(main())

"""Free-normal 3D SPH block with unilateral shoulder penalty and tool lift.

Named companion, NOT AFSD deposition. Reuses prior atomic drift/kick/drift and
unchanged material/spatial/contact implementations; no global monkey-patching.
Normal potential and its temporal residual are separate from raw friction heat.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import probe_rotating_deformable_block as base
from probe_rotating_deformable_block import (positions, velocities, corrected_operators,
    internal_force, pressure, bulk_modulus, FKS, update, HostLeapfrog, kick,
    save_checkpoint, restore, frozen_gradient, config_hash, jsonable, write_json,
    rejected_jsonable, write_point_cloud)
from afsd_backend.pysph.normal_penalty import normal_state, moving_tool_rows, tool_state

CONFIG_PATH = ROOT/'configs/benchmarks/normal_deformable_block.json'
INPUT_CONFIG = json.loads(CONFIG_PATH.read_text(encoding='utf-8'))
PARAMETERS = {**base.PARAMETERS,
    **{name: item['value'] for name, item in INPUT_CONFIG['parameters'].items()}}
CLASSIFICATION = INPUT_CONFIG['classification']
SOURCES = list(dict.fromkeys(base.SOURCES+[
    'scripts/probe_normal_deformable_block.py',
    'src/afsd_backend/pysph/normal_penalty.py',
    'configs/benchmarks/normal_deformable_block.json',
    'docs/NORMAL_CONTACT_CONTRACT.md',
    'third_party/PySPH/pysph/sph/rigid_body.py']))


def hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def configuration(dt):
    result = base.configuration(dt)
    result.update(classification=CLASSIFICATION, source_configuration=INPUT_CONFIG,
        inherited_source_configuration=base.INPUT_CONFIG, parameters=PARAMETERS,
        normal_law='UNILATERAL LINEAR PENALTY; TOP/INTERIOR XYZ FREE; BOTTOM XYZ FIXED',
        fixture='PRESCRIBED RIGID TOOL HOLD THEN LIFT; NO FREE-MATERIAL TRAJECTORY')
    result['field_status'] = dict(result['field_status'],
        normal_force='SOLVED_FROM_CURRENT_OVERLAP',
        normal_potential='SOLVED_PENALTY_POTENTIAL_NOT_HEAT',
        normal_temporal_residual='RAW_DIAGNOSTIC_NOT_HEAT',
        top_normal_displacement='SOLVED_FREE_DOF')
    return result


class Block(base.Block):
    def __init__(self, dt, pa=None, state=None):
        if dt not in INPUT_CONFIG['dt_cases_s']:
            raise ValueError('only preregistered normal-component timesteps permitted')
        # Base constructor creates the same unscaled 64-particle block. State is
        # validated against this companion below, not against the old variant.
        super().__init__(dt, pa=pa)
        self.p = dict(PARAMETERS); self.config = configuration(float(dt)); self.source_hashes = hashes()
        self.bulk = bulk_modulus(self.p['shear'], self.p['poisson'])
        self.sound = float(np.sqrt((self.bulk+4*self.p['shear']/3)/self.p['rho0']))
        self.law = FKS(**{k: self.p[k] for k in FKS.__dataclass_fields__})
        self.kn = float(self.p['normal_stiffness'])
        initial = np.column_stack((self.pa.x0, self.pa.y0, self.pa.z0))
        self.initial_normal_potential_J = float(np.sum(self.normal_snapshot(initial, 0.)['potential_J']))
        self.cumulative.update(normal_wall_work_J=0., normal_relative_work_J=0.,
            normal_temporal_residual_J=0.)
        if state is not None:
            if not state.get('accepted') or not state.get('restart_allowed'):
                raise ValueError('only accepted normal-component checkpoint may restart')
            if state['classification'] != CLASSIFICATION:
                raise ValueError('checkpoint classification mismatch')
            if state['configuration_sha256'] != config_hash(self.config):
                raise ValueError('checkpoint configuration mismatch')
            if state['source_sha256'] != self.source_hashes:
                raise ValueError('checkpoint source mismatch')
            if state['initial_normal_potential_J'] != self.initial_normal_potential_J:
                raise ValueError('checkpoint initial normal energy mismatch')
            self.time = float(state['time_s']); self.step = int(state['step'])
            self.cumulative = dict(state['cumulative'])
        self.integrator = HostLeapfrog(self.pa, self._callback)

    def stable_dt(self):
        return float(min(super().stable_dt(), .1*np.sqrt(self.pa.m.min()/self.kn)))

    def normal_snapshot(self, x, time_s):
        # Initial surface geometry selects the diagnostic cohort. Role 2 is
        # observational and cannot constrain motion or set a contact force.
        top = np.flatnonzero(self.pa.z0 == self.pa.z0.max())
        tool = tool_state(time_s, self.p, self.p['lift_start'], self.p['lift_speed'])
        gap = tool['layer']-np.asarray(x)[top, 2]
        values = normal_state(gap, self.p['dx']/2, self.kn)
        return dict(time_s=time_s, particle_id=self.pa.particle_id[top].copy(),
            particle_index=top, gap_m=gap, radius_m=self.p['dx']/2, stiffness_N_m=self.kn,
            signed_penetration_m=values['signed_penetration'], penetration_m=values['penetration'],
            active=values['active'], normal_force_N=values['normal_force_N'],
            potential_J=values['potential_J'], tool=tool)

    def advance(self):
        # This input snapshot is diagnostic only. Base advance atomically rolls
        # back ParticleArray, cumulative ledger and last accepted journal/traces.
        self._normal_before = self.normal_snapshot(positions(self.pa), self.time)
        return super().advance()

    def _callback(self, midpoint_time, dt):
        p = self.pa; x = positions(p); v = velocities(p); rho = p.rho.copy()
        i, j, g, L, drho, B = corrected_operators(p)
        old_s = p.S.reshape(-1, 3, 3).copy()
        new_s, ep, plastic_density, return_density, rate = update(
            old_s, p.ep.copy(), L, p.temperature, dt, self.p['shear'], self.law)
        rho_mid = rho+.5*dt*drho
        sigma = .5*(old_s+new_s)-pressure(rho_mid, self.p['rho0'], self.bulk)[:, None, None]*np.eye(3)
        force = internal_force(p, sigma, i, j, g, correction=B)
        rows = moving_tool_rows(x, p.particle_id, midpoint_time, self.p,
                               self.p['lift_start'], self.p['lift_speed'])
        top = p.z0 == p.z0.max()
        if (np.any(rows['feature_id'] != 1) or
                not set(rows['particle_id']).issubset(set(p.particle_id[top]))):
            raise ValueError('normal component supports top shoulder contacts only')
        normal_mid = self.normal_snapshot(x, midpoint_time)
        normal_forces = normal_state(rows['gap'], self.p['dx']/2, self.kn)['normal_force_N']
        after = np.empty_like(v); work_v = np.empty_like(v); traces = []
        self.last_traces = traces
        self.last_journal = dict(stage='CANDIDATE_MATERIAL_BEFORE_CONTACT',
            step=self.step+1, time_mid_s=midpoint_time, dt_s=dt,
            particle_id=p.particle_id.copy(), x_mid_m=x.copy(), v_before_m_s=v.copy(),
            rho_before_kg_m3=rho, rho_mid_kg_m3=rho_mid, S_before_Pa=old_s,
            S_after_Pa=new_s.copy(), ep_before=p.ep.copy(), ep_after=ep.copy(),
            L_before_per_s=L, sigma_mid_Pa=sigma, internal_force_N=force,
            i=i, j=j, g=g, B=B, traces=traces,
            normal_before=self._normal_before, normal_midpoint=normal_mid)
        for idx, pid in enumerate(p.particle_id):
            selected = np.flatnonzero(rows['particle_index'] == idx)
            faces = [{key: value if isinstance(value, str) else np.asarray(value)[r].copy()
                      for key, value in rows.items()} for r in selected]
            for face, row in zip(faces, selected):
                face['normal_force_N'] = float(normal_forces[row])
            free = [0., 0., 0.] if p.boundary_role[idx] == 1 else [1., 1., 1.]
            result = kick(v[idx], p.m[idx], force[idx], rows['normal'][selected],
                normal_forces[selected], rows['wall_velocity'][selected], self.p['friction'], dt,
                free_axes=free, prescribed_velocity=np.zeros(3))
            after[idx] = result['velocity']
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
        old_work = float(dt*np.sum(volume[:, None, None]*sigma*L))
        condition = np.linalg.cond(np.linalg.inv(B))
        diagnostics = dict(pair_count=len(i), condition_max=float(condition.max()),
            internal_power_W=power_force, midpoint_stress_power_W=power_stress,
            power_conjugacy_residual_W=power_force+power_stress,
            internal_net_force_N=force.sum(axis=0), internal_net_torque_Nm=np.cross(x, force).sum(axis=0),
            force_absolute_scale_N=float(np.linalg.norm(force, axis=1).sum()),
            torque_absolute_scale_Nm=float(np.linalg.norm(np.cross(x, force), axis=1).sum()),
            power_absolute_scale_W=float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(work_v, axis=1))),
            material_old_work_J=old_work, midpoint_stress_work_J=dt*power_stress,
            temporal_stress_work_split_J=dt*power_stress-old_work)
        self.last_journal = dict(step=self.step+1, time_mid_s=midpoint_time, dt_s=dt,
            stage='CANDIDATE_MATERIAL_AND_CONTACT_COMPLETE',
            particle_id=p.particle_id.copy(), boundary_role=p.boundary_role.copy(),
            x_mid_m=x.copy(), v_before_m_s=v.copy(), v_after_m_s=after.copy(),
            work_velocity_m_s=work_v, endpoint_midpoint_velocity_m_s=.5*(v+after),
            mass_kg=p.m.copy(), rho_before_kg_m3=rho, force_density_kg_m3=rho.copy(),
            smoothing_length_m=p.h.copy(), moment_condition_number=condition,
            i=i, j=j, g=g, B=B, L_before_per_s=L, L_midpoint_per_s=Lmid,
            sigma_mid_Pa=sigma, S_before_Pa=old_s, S_after_Pa=new_s.copy(),
            ep_before=p.ep.copy(), ep_after=ep.copy(), internal_force_N=force,
            drho_kg_m3_s=drho, rho_mid_kg_m3=rho_mid,
            plastic_work_J=volume*plastic_density, return_dissipation_J=volume*return_density,
            traces=traces, diagnostics=diagnostics,
            normal_before=self._normal_before, normal_midpoint=normal_mid)
        p.au[:], p.av[:], p.aw[:] = ((after-v)/dt).T
        p.arho[:] = drho
        p.ep_rate[:] = (ep-p.ep)/dt; p.S[:] = new_s.ravel(); p.ep[:] = ep
        p.strain_rate[:] = rate; p.D[:] = (.5*(L+L.swapaxes(1, 2))).ravel()
        p.plastic_work[:] += volume*plastic_density; p.return_dissipation[:] += volume*return_density

    def _validate_candidate(self, increment):
        p = self.pa; journal = self.last_journal
        normal_after = self.normal_snapshot(positions(p), self.time+self.dt)
        journal['normal_after'] = normal_after
        normal_wall = 0.
        for trace in self.last_traces:
            walls = np.asarray([row['wall_velocity'] for row in trace['contact_rows']]).reshape(-1, 3)
            for interval in trace['result']['intervals']:
                free_reaction = np.asarray(trace['free_axes'])*interval['constraint_impulse']
                if np.max(np.abs(free_reaction)) > 1e-17:
                    raise ValueError('reaction impulse on a free degree of freedom')
                normal_wall += float(np.sum(interval['normal_impulse']*walls))
        delta_u = float(np.sum(normal_after['potential_J'])-np.sum(journal['normal_before']['potential_J']))
        increment['normal_wall_work_J'] = normal_wall
        increment['normal_relative_work_J'] = increment['normal_work_J']-normal_wall
        increment['normal_temporal_residual_J'] = delta_u+increment['normal_relative_work_J']
        journal['normal_energy'] = dict(potential_change_J=delta_u,
            material_work_J=increment['normal_work_J'], wall_work_J=normal_wall,
            temporal_residual_J=increment['normal_temporal_residual_J'],
            status='RAW_TEMPORAL_DISCRETIZATION_DIAGNOSTIC_NOT_HEAT')
        for key, array in p.properties.items():
            if not np.isfinite(array.get_npy_array()).all(): raise ValueError('nonfinite candidate '+key)
        if abs(float(p.m.sum())-self.initial_mass) > 1e-12*self.initial_mass:
            raise ValueError('mass conservation limit')
        if np.max(np.abs(p.rho/self.p['rho0']-1)) > .01: raise ValueError('density change limit')
        if np.linalg.norm(velocities(p), axis=1).max() > 1.: raise ValueError('speed limit')
        original = np.column_stack((p.x0, p.y0, p.z0))
        if np.linalg.norm(positions(p)-original, axis=1).max() > 1e-5:
            raise ValueError('displacement limit')
        fixed = p.boundary_role == 1
        if (np.max(np.abs(velocities(p)[fixed])) > 1e-12 or
                np.max(np.abs((positions(p)-original)[fixed])) > 2e-17):
            raise ValueError('bottom fixed-DOF motion limit')
        for state in (journal['normal_before'], journal['normal_midpoint'], normal_after):
            if np.max(state['penetration_m']) > .1*self.p['dx']:
                raise ValueError('normal penetration limit')
        d = journal['diagnostics']
        if d['condition_max'] > 1e6: raise ValueError('condition limit')
        if abs(d['power_conjugacy_residual_W']) > 1e-12+1e-10*d['power_absolute_scale_W']:
            raise ValueError('spatial power conjugacy limit')
        if np.linalg.norm(d['internal_net_force_N']) > 1e-12+1e-10*d['force_absolute_scale_N']:
            raise ValueError('net internal force limit')
        if np.linalg.norm(d['internal_net_torque_Nm']) > 1e-14+1e-10*d['torque_absolute_scale_Nm']:
            raise ValueError('net internal torque limit')
        residual = (increment['delta_kinetic_J']+increment['friction_raw_J']-
            increment['wall_work_J']-increment['constraint_work_J']-increment['internal_work_J']-
            increment['normal_relative_work_J'])
        scale = sum(abs(increment[k]) for k in ('delta_kinetic_J', 'friction_raw_J',
                    'wall_work_J', 'constraint_work_J', 'internal_work_J', 'normal_relative_work_J'))
        if abs(residual) > 1e-15+1e-9*scale:
            raise ValueError('kinetic/contact/normal work identity limit')

    def record(self):
        record = super().record()
        normal = self.normal_snapshot(positions(self.pa), self.time)
        potential = float(np.sum(normal['potential_J']))
        c = self.cumulative
        full_residual = record['energy_diagnostic_residual_J']+potential-self.initial_normal_potential_J
        scale = (self.initial_normal_potential_J+abs(c['wall_work_J'])+
                 abs(c['constraint_work_J'])+abs(c['friction_raw_J']))
        record.update(normal_contact=normal, normal_potential_J=potential,
            initial_normal_potential_J=self.initial_normal_potential_J,
            contact_count=int(np.count_nonzero(normal['active'])),
            max_penetration_m=float(np.max(normal['penetration_m'])),
            tool_origin_m=normal['tool']['origin'],
            tool_translation_velocity_m_s=normal['tool']['translation_velocity'],
            kinetic_contact_residual_J=record['kinetic_contact_residual_J']-c['normal_relative_work_J'],
            energy_diagnostic_residual_J=full_residual,
            energy_status='ENERGY_WARNING' if abs(full_residual) > 1e-12+.01*scale else
                'DIAGNOSTIC_WITHIN_SCREEN_NOT_CERTIFIED_FREE_ENERGY',
            diagnostic_status='NORMAL_TEMPORAL_RESIDUAL_REPORTED_NOT_HEAT; OLD_JAUMANN_FAIL_REMAINS')
        return record

    def checkpoint(self, path):
        save_checkpoint(path, self.pa, dict(classification=CLASSIFICATION, accepted=True,
            restart_allowed=True, production_authorized=False, dt_s=self.dt, step=self.step,
            time_s=self.time, cumulative=self.cumulative, configuration=self.config,
            configuration_sha256=config_hash(self.config), source_sha256=self.source_hashes,
            initial_normal_potential_J=self.initial_normal_potential_J))

    def vtk(self, path):
        p = self.pa; r = self.record()
        scalar = {k: getattr(p, k).tolist() for k in ('particle_id', 'material_id', 'component_id',
            'tracer_type', 'wholebar_tag', 'actuator', 'boundary_role', 'as_rs_tag')}
        scalar.update(mass_kg=p.m.tolist(), density_kg_per_m3=p.rho.tolist(),
            pressure_Pa=r['pressure_Pa'].tolist(), equivalent_stress_Pa=r['equivalent_stress_Pa'].tolist(),
            equivalent_plastic_strain=p.ep.tolist(), smoothing_length_m=p.h.tolist())
        gap = r['normal_contact']['tool']['layer']-p.z
        normal = normal_state(gap, self.p['dx']/2, self.kn)
        scalar.update(signed_gap_m=gap.tolist(),
            normal_signed_penetration_m=normal['signed_penetration'].tolist(),
            normal_penetration_m=normal['penetration'].tolist(),
            normal_force_N=normal['normal_force_N'].tolist(),
            normal_potential_J=normal['potential_J'].tolist(),
            contact_active=normal['active'].astype(int).tolist())
        for a, b, name in ((0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'),
                           (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')):
            scalar['stress_'+name+'_Pa'] = r['stress_Pa'][:, a, b].tolist()
        if self.step:
            scalar['equivalent_deviatoric_strain_rate_per_s'] = p.strain_rate.tolist()
            scalar['equivalent_plastic_strain_rate_per_s'] = p.ep_rate.tolist()
        write_point_cloud(path, positions(p).tolist(), scalar,
            dict(velocity_m_per_s=velocities(p).tolist(),
                 initial_position_m=np.column_stack((p.x0, p.y0, p.z0)).tolist()),
            title=f'{CLASSIFICATION} time_s={self.time:.17g}; FREE_NORMAL_CONTACT_ISOTHERMAL_NO_THERMAL_SOLUTION')


def run(dt, output, probe_steps=None, timeout_s=None):
    output = Path(output)
    if output.exists(): raise ValueError('fresh output directory required')
    if probe_steps is not None and probe_steps != 2: raise ValueError('only two-step cost probe allowed')
    timeout_s = INPUT_CONFIG['run_timeout_s'] if timeout_s is None else timeout_s
    if not 0 < timeout_s <= INPUT_CONFIG['run_timeout_s']: raise ValueError('bounded wall budget required')
    output.mkdir(parents=True)
    start = time.perf_counter(); block = Block(dt)
    target_steps = 2 if probe_steps else int(round(PARAMETERS['end_time']/dt))
    manifest = dict(classification=CLASSIFICATION, status='RUNNING', configuration=block.config,
        configuration_sha256=config_hash(block.config), source_sha256=block.source_hashes,
        maximum_wall_seconds=timeout_s, production_authorized=False, target_steps=target_steps,
        requested_end_time_s=target_steps*dt,
        mode='TWO_STEP_COST_PROBE' if probe_steps else 'FULL_REGISTERED_NORMAL_COMPONENT',
        initial_stable_dt_s=block.stable_dt(), sound_speed_m_s=block.sound,
        contact_dt_bound_s=float(.1*np.sqrt(block.pa.m.min()/block.kn)),
        initial_normal_potential_J=block.initial_normal_potential_J,
        energy_note='Es uncertified diagnostic; normal potential separate; no residual heat or thermal solve.')
    write_json(output/'manifest.json', manifest)
    records = [block.record()]; block.checkpoint(output/'last_accepted.npz'); block.vtk(output/'frame_0000.vtk')
    step_times = []
    try:
        with (output/'steps.jsonl').open('w', encoding='utf-8') as journal:
            while block.step < target_steps:
                if time.perf_counter()-start > timeout_s: raise TimeoutError('normal component wall budget exceeded')
                step_start = time.perf_counter(); records.append(block.advance())
                step_times.append(time.perf_counter()-step_start)
                journal.write(json.dumps(jsonable(block.last_journal), allow_nan=False)+'\n'); journal.flush()
                block.checkpoint(output/'last_accepted.npz')
                if block.step == int(round(PARAMETERS['lift_start']/dt)):
                    block.checkpoint(output/'prelift_accepted.npz')
                block.vtk(output/f'frame_{block.step:04d}.vtk')
        manifest['status'] = 'COMPLETE_COST_PROBE_NOT_AFSD' if probe_steps else 'COMPLETE_NORMAL_DEFORMABLE_COMPONENT_NOT_AFSD'
    except Exception as error:
        manifest.update(status='BLOCKED', error=str(error))
        if hasattr(block, 'rejected_candidate'):
            write_json(output/'rejected_candidate.json', rejected_jsonable(block.rejected_candidate))
    finally:
        write_json(output/'history.json', records)
        manifest.update(wall_seconds=time.perf_counter()-start, accepted_steps=block.step,
            accepted_time_s=block.time, final_record=block.record(), step_compute_seconds=step_times,
            particle_array_bytes=sum(a.get_npy_array().nbytes for a in block.pa.properties.values()),
            memory_note='ParticleArray bytes only, NOT peak process RSS.')
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

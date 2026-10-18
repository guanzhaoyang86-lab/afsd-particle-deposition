"""Manufactured free-mass rotating-shoulder component, NOT an AFSD solid.

Uses the existing PySPH drift/kick/drift and event-aware contact kick unchanged.
The normal load and normal constraint are prescribed benchmark inputs; x/y
motion is solved. No material, NNPS, thermal or finite-strain field is solved.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.integrator import HostLeapfrog
from afsd_backend.pysph.contact_kick import kick
from afsd_backend.pysph.output import checkpoint, restore
from afsd_scaffold.vtk_legacy import write_point_cloud
from afsd_backend.pysph.rigid_contact_kinematics import contact_rows

CONFIG_PATH = ROOT/'configs/benchmarks/rotating_shoulder_component.json'
INPUT_CONFIG = json.loads(CONFIG_PATH.read_text(encoding='utf-8'))
CLASSIFICATION = INPUT_CONFIG['classification']
PARAMETERS = {name: entry['value'] for name, entry in INPUT_CONFIG['parameters'].items()}
INITIAL_POSITIONS = np.array(PARAMETERS['initial_position'], dtype=float)
SOURCES = ['scripts/probe_rotating_shoulder_component.py',
           'configs/benchmarks/rotating_shoulder_component.json',
           'docs/ROTATING_CONTACT_CONTRACT.md',
           'src/afsd_backend/pysph/rigid_contact_kinematics.py',
           'src/afsd_backend/pysph/geometry.py',
           'src/afsd_backend/pysph/contact_kick.py',
           'src/afsd_backend/pysph/contact_cone.py',
           'src/afsd_backend/pysph/integrator.py',
           'src/afsd_backend/pysph/output.py',
           'src/afsd_scaffold/vtk_legacy.py',
           'third_party/PySPH/pysph/sph/integrator.py',
           'third_party/PySPH/pysph/sph/integrator_step.py']


def jsonable(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [jsonable(v) for v in value]
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(jsonable(value), indent=2, allow_nan=False)+'\n', encoding='utf-8')


def hashes():
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}


def configuration(dt):
    return dict(classification=CLASSIFICATION, source_configuration=INPUT_CONFIG, parameters=PARAMETERS,
                initial_positions_m=INITIAL_POSITIONS.tolist(), initial_velocity_m_s=[0., 0., 0.],
                free_axes=[1, 1, 0], dt_s=float(dt), ids=[101, 102, 103, 104],
                tags_are_observation_only=True, physical_scaling='NONE',
                normal_law='PRESCRIBED 1 N; NOT PENALTY OR DEFORMABLE SOLID',
                field_status={'position': 'SOLVED', 'velocity': 'SOLVED',
                    'mass': 'PRESCRIBED_MANUFACTURED_INPUT',
                    'friction_dissipation': 'SOLVED_NOT_DEPOSITED_AS_TEMPERATURE',
                    'stress': 'NOT_SOLVED', 'strain': 'NOT_SOLVED',
                    'temperature': 'NOT_SOLVED', 'density': 'NOT_SOLVED',
                    'normal_constraint_reaction': 'SOLVED',
                    'checkpoint_rho': 'UPSTREAM_SCRATCH_NOT_SOLVED_NOT_PHYSICAL_DENSITY',
                    'checkpoint_e': 'UPSTREAM_SCRATCH_NOT_SOLVED_NOT_PHYSICAL_ENERGY',
                    'checkpoint_h': 'CONTACT_DETECTION_METADATA_NOT_SPH_RESOLUTION'})


def config_hash(config):
    return hashlib.sha256(json.dumps(config, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def positions(pa): return np.column_stack((pa.x, pa.y, pa.z))
def velocities(pa): return np.column_stack((pa.u, pa.v, pa.w))


def reference_rk4(dt, final_time=PARAMETERS['end_time']):
    """Independent continuous ODE; never calls any production contact function.

    State columns are r(3), v(3), Q, W_wall. Applicable only before zero slip.
    This is not a reference for stick/reversal, solved normal contact, or solids.
    """
    literal_parameters = dict(mass=.001, omega=20., mu=.2, normal_force=1.,
        traverse=0., initial_position=[[.01, 0., .0013], [0., .012, .0013],
                                      [-.01, 0., .0013], [0., -.012, .0013]],
        initial_velocity=[0., 0., 0.], end_time=.0008)
    for key, expected in literal_parameters.items():
        if (INPUT_CONFIG['parameters'][key]['value'] != expected or PARAMETERS[key] != expected):
            raise ValueError('independent literal reference configuration mismatch: '+key)
    if not 0. < final_time <= literal_parameters['end_time']:
        raise ValueError('reference time outside frozen interval')
    reference_parameters = dict(mass_kg=.001, omega_rad_s=[0., 0., 20.],
        friction_force_cap_N=.2, normal_force_N=1., mu=.2, traverse_m_s=0.,
        initial_position_m=literal_parameters['initial_position'], initial_velocity_m_s=[0., 0., 0.],
        full_interval_end_s=.0008, requested_end_s=final_time, source='INDEPENDENT_FROZEN_LITERALS')
    reference_config = dict(source_configuration=INPUT_CONFIG, reference_parameters=reference_parameters,
                            reference_dt_s=dt)
    count = int(round(final_time/dt))
    if count < 1 or abs(count*dt-final_time) > 1e-14:
        raise ValueError('reference step must divide final time')
    y = np.zeros((4, 8)); y[:, :3] = literal_parameters['initial_position']
    minimum_slip = float('inf')
    def rhs(state):
        nonlocal minimum_slip
        r, v = state[:, :3], state[:, 3:6]
        wall = np.column_stack((-20.*r[:, 1], 20.*r[:, 0], np.zeros(4)))
        slip = wall-v
        speed = np.linalg.norm(slip, axis=1)
        minimum_slip = min(minimum_slip, float(speed.min()))
        if speed.min() < .03:
            raise ValueError('independent ODE left preregistered continuous-slip domain')
        force = .2*slip/speed[:, None]
        out = np.zeros_like(state)
        out[:, :3] = v; out[:, 3:6] = force/.001
        out[:, 6] = np.sum(force*slip, axis=1)
        out[:, 7] = np.sum(force*wall, axis=1)
        return out
    for _ in range(count):
        a = rhs(y); b = rhs(y+.5*dt*a); c = rhs(y+.5*dt*b); d = rhs(y+dt*c)
        y += dt*(a+2*b+2*c+d)/6
    return dict(method='INDEPENDENT_RK4_CONTINUOUS_SLIP_ODE_NOT_CONTACT_KICK', dt_s=dt,
                reference_parameters=reference_parameters, source_sha256=hashes(),
                configuration=reference_config, configuration_sha256=config_hash(reference_config),
                steps=count, time_s=final_time, position_m=y[:, :3], velocity_m_s=y[:, 3:6],
                friction_raw_J=y[:, 6], wall_work_J=y[:, 7],
                kinetic_J=.0005*np.sum(y[:, 3:6]**2, axis=1), minimum_slip_m_s=minimum_slip)


class RotatingShoulderComponent:
    def __init__(self, dt, pa=None, state=None):
        if dt not in INPUT_CONFIG['dt_cases_s']: raise ValueError('only preregistered dt = 2e-5 or 1e-5 s')
        self.dt = dt; self.config = configuration(dt); self.source_hashes = hashes()
        if pa is None:
            x = INITIAL_POSITIONS
            pa = get_particle_array(name='solid', x=x[:, 0], y=x[:, 1], z=x[:, 2],
                                    m=np.full(4, PARAMETERS['mass']), h=np.full(4, PARAMETERS['dx']))
            for name in ('ax', 'ay', 'az', 'au', 'av', 'aw', 'arho', 'ae', 'e'):
                if name not in pa.properties: pa.add_property(name)
            for name, values in [('particle_id', [101, 102, 103, 104]),
                                  ('tracer_type', [0, 1, 2, 1]), ('as_rs_tag', [0, 0, 0, 0])]:
                pa.add_property(name, type='long', data=np.array(values, dtype=np.int64))
            for k, name in enumerate(('x0', 'y0', 'z0')): pa.add_property(name, data=x[:, k])
            pa.set_output_arrays(list(pa.properties))
        self.pa = pa; self.time = 0.; self.step = 0; self.last_traces = []
        self.cumulative = dict(friction_raw_J=0., wall_work_J=0., constraint_work_J=0.,
                              other_work_J=0., normal_work_J=0., delta_kinetic_J=0.,
                              contact_kinetic_residual_J=0.)
        if state is not None:
            if state['configuration_sha256'] != config_hash(self.config): raise ValueError('checkpoint config mismatch')
            if state['source_sha256'] != self.source_hashes: raise ValueError('checkpoint source mismatch')
            if state['classification'] != CLASSIFICATION: raise ValueError('checkpoint classification mismatch')
            self.time = state['time_s']; self.step = state['step']; self.cumulative = dict(state['cumulative'])
        self.integrator = HostLeapfrog(self.pa, self._callback)

    @classmethod
    def from_checkpoint(cls, path):
        pa, state = restore(path)
        return cls(state['dt_s'], pa=pa, state=state)

    def _callback(self, midpoint_time, dt):
        p = self.pa; x = positions(p); v = velocities(p)
        rows = contact_rows(x, p.particle_id, midpoint_time, PARAMETERS)
        traces = []; acceleration = np.zeros((4, 3))
        for i, pid in enumerate(p.particle_id):
            selected = np.flatnonzero(np.asarray(rows['particle_id']) == pid)
            if len(selected) != 1 or int(rows['feature_id'][selected[0]]) != 1:
                raise ValueError('manufactured benchmark requires exactly one shoulder row per particle')
            row = int(selected[0])
            result = kick(v[i], p.m[i], np.zeros(3), [rows['normal'][row]], [PARAMETERS['normal_force']],
                          [rows['wall_velocity'][row]], PARAMETERS['mu'], dt,
                          free_axes=[1., 1., 0.], prescribed_velocity=[v[i, 0], v[i, 1], 0.])
            acceleration[i] = (result['velocity']-v[i])/dt
            traces.append(dict(particle_id=int(pid), stage='MIDPOINT_FROZEN_GEOMETRY_KICK',
                time_mid_s=midpoint_time, step=self.step+1, position_m=x[i].copy(),
                feature_id=1, contact_row={k: a if isinstance(a, str) else np.asarray(a)[row].copy()
                                          for k, a in rows.items()},
                material_velocity_before_m_s=v[i].copy(), result=result))
        p.au[:] = acceleration[:, 0]; p.av[:] = acceleration[:, 1]; p.aw[:] = acceleration[:, 2]
        self.last_traces = traces

    def advance(self):
        """Atomically accept one unmodified upstream drift/kick/drift step."""
        if self.time+self.dt > PARAMETERS['end_time']+1e-14: raise ValueError('beyond frozen end time')
        before = {k: a.get_npy_array().copy() for k, a in self.pa.properties.items()}
        try:
            self.integrator.advance(self.time, self.dt)
            if not np.all(np.isfinite(positions(self.pa))) or not np.all(np.isfinite(velocities(self.pa))):
                raise ValueError('nonfinite candidate')
            if not np.array_equal(self.pa.z, INITIAL_POSITIONS[:, 2]) or np.any(self.pa.w != 0.):
                raise ValueError('normal constraint changed z position/velocity')
            increment = dict.fromkeys(self.cumulative, 0.)
            for trace in self.last_traces:
                result = trace['result']
                increment['friction_raw_J'] += float(result['friction_heat_J'].sum())
                for key in ('wall_work_J', 'constraint_work_J', 'other_work_J', 'delta_kinetic_J'):
                    increment[key] += result[key]
                increment['normal_work_J'] += sum(float(a['normal_work_J'].sum()) for a in result['intervals'])
                increment['contact_kinetic_residual_J'] += result['kinetic_work_residual_J']
            if abs(increment['delta_kinetic_J']+increment['friction_raw_J']-increment['wall_work_J']-
                   increment['constraint_work_J']-increment['other_work_J']) > 5e-12:
                raise ValueError('accepted-candidate work identity failed')
        except Exception:
            self.rejected_candidate = dict(classification='REJECTED', restart_allowed=False,
                time_s=self.time+self.dt, position_m=positions(self.pa).copy(), velocity_m_s=velocities(self.pa).copy())
            for key, a in before.items(): getattr(self.pa, key)[:] = a
            raise
        for key in self.cumulative: self.cumulative[key] += increment[key]
        self.step += 1; self.time = self.step*self.dt
        return self.record()

    def record(self):
        p = self.pa
        kinetic = float(.5*np.sum(p.m[:, None]*velocities(p)**2))
        return dict(step=self.step, time_s=self.time, particle_id=p.particle_id.copy(),
                    tool_origin_m=[PARAMETERS['traverse']*self.time, 0., PARAMETERS['layer']],
                    tool_angle_rad=PARAMETERS['omega']*self.time,
                    tracer_type=p.tracer_type.copy(), as_rs_tag=p.as_rs_tag.copy(),
                    position_m=positions(p), velocity_m_s=velocities(p), kinetic_J=kinetic,
                    cumulative=dict(self.cumulative), global_energy_residual_J=kinetic+
                    self.cumulative['friction_raw_J']-self.cumulative['wall_work_J']-
                    self.cumulative['constraint_work_J']-self.cumulative['other_work_J'])

    def save_checkpoint(self, path):
        checkpoint(path, self.pa, dict(classification=CLASSIFICATION, accepted=True, restart_allowed=True,
            production_authorized=False, dt_s=self.dt, step=self.step, time_s=self.time,
            cumulative=self.cumulative, configuration=self.config, configuration_sha256=config_hash(self.config),
            source_sha256=self.source_hashes))

    def vtk(self, path):
        p = self.pa
        write_point_cloud(path, positions(p).tolist(),
            dict(particle_id=p.particle_id.tolist(), tracer_type=p.tracer_type.tolist(),
                 as_rs_tag=p.as_rs_tag.tolist(), mass_kg=p.m.tolist()),
            dict(velocity_m_per_s=velocities(p).tolist()),
            title=f'{CLASSIFICATION} time_s={self.time:.17g}; STRESS_THERMAL_NOT_SOLVED')


def run(dt, output, timeout_s=55.):
    output = Path(output)
    if output.exists(): raise ValueError('fresh output directory required; no historical overwrite')
    output.mkdir(parents=True)
    start = time.perf_counter(); block = RotatingShoulderComponent(dt)
    manifest = dict(classification=CLASSIFICATION, status='RUNNING', configuration=block.config,
                    configuration_sha256=config_hash(block.config), source_sha256=block.source_hashes,
                    maximum_wall_seconds=timeout_s, production_authorized=False,
                    references='Independent RK4 is a separate reference mode, not source of solver trajectories.')
    write_json(output/'manifest.json', manifest)
    records = [block.record()]; block.save_checkpoint(output/'last_accepted.npz'); block.vtk(output/'frame_0000.vtk')
    try:
        with (output/'contact_intervals.jsonl').open('w', encoding='utf-8') as journal:
            while block.step < round(PARAMETERS['end_time']/dt):
                if time.perf_counter()-start > timeout_s: raise TimeoutError('component wall budget exceeded')
                records.append(block.advance())
                journal.write(json.dumps(jsonable(block.last_traces), allow_nan=False)+'\n'); journal.flush()
                block.save_checkpoint(output/'last_accepted.npz')
                block.vtk(output/f'frame_{block.step:04d}.vtk')
        manifest['status'] = 'COMPLETE_MANUFACTURED_COMPONENT_NOT_AFSD'
    except Exception as error:
        manifest.update(status='BLOCKED', error=str(error))
        if hasattr(block, 'rejected_candidate'): write_json(output/'rejected_candidate.json', block.rejected_candidate)
    finally:
        write_json(output/'history.json', records)
        manifest.update(wall_seconds=time.perf_counter()-start, accepted_steps=block.step,
                        accepted_time_s=block.time, final_record=block.record(),
                        particle_array_bytes=sum(a.get_npy_array().nbytes for a in block.pa.properties.values()),
                        memory_note='ParticleArray property bytes only, NOT peak process RSS.')
        manifest['artifact_sha256'] = {str(path.relative_to(output)).replace('\\', '/'):
            hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(output.rglob('*'))
            if path.is_file() and path != output/'manifest.json'}
        write_json(output/'manifest.json', manifest)
    print(json.dumps(jsonable(manifest), allow_nan=False), flush=True)
    return 0 if manifest['status'].startswith('COMPLETE') else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dt', type=float, choices=[2e-5, 1e-5], default=2e-5)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-dt', type=float, choices=[1e-7, 5e-8])
    args = parser.parse_args()
    if args.reference_dt is not None:
        if args.output.exists(): raise ValueError('fresh output file required')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        start = time.perf_counter(); result = reference_rk4(args.reference_dt)
        result['wall_seconds'] = time.perf_counter()-start
        write_json(args.output, result); print(json.dumps(jsonable(result)), flush=True); return 0
    return run(args.dt, args.output)


if __name__ == '__main__': raise SystemExit(main())

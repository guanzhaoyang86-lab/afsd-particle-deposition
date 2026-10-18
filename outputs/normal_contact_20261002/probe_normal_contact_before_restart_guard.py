"""Independent-reference normal impact benchmark; manufactured mass, NOT AFSD.

Unchanged PySPH HostLeapfrog drift/kick/drift and event-aware contact kick.
Normal force is frozen at drifted midpoint geometry. Opening is unilateral;
normal spring work defects are reported, never deposited as friction heat.
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
from afsd_backend.pysph.normal_penalty import normal_state
from afsd_backend.pysph.output import checkpoint, restore
from afsd_scaffold.vtk_legacy import write_point_cloud

CLASSIFICATION = 'MANUFACTURED_NORMAL_IMPACT_NOT_AFSD'
CONFIG_PATH = ROOT/'configs/benchmarks/normal_contact.json'
SOURCES = ['scripts/probe_normal_contact.py', 'src/afsd_backend/pysph/normal_penalty.py',
           'configs/benchmarks/normal_contact.json', 'docs/NORMAL_CONTACT_CONTRACT.md',
           'src/afsd_backend/pysph/rigid_contact_kinematics.py', 'src/afsd_backend/pysph/geometry.py',
           'src/afsd_backend/pysph/contact.py', 'src/afsd_backend/pysph/contact_kick.py',
           'src/afsd_backend/pysph/contact_cone.py', 'src/afsd_backend/pysph/integrator.py',
           'src/afsd_backend/pysph/output.py', 'src/afsd_scaffold/vtk_legacy.py',
           'third_party/PySPH/pysph/sph/rigid_body.py',
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


def config_hash(config):
    return hashlib.sha256(json.dumps(config, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def positions(pa): return np.column_stack((pa.x, pa.y, pa.z))
def velocities(pa): return np.column_stack((pa.u, pa.v, pa.w))


def configuration(dt, wall_speed):
    inputs = json.loads(CONFIG_PATH.read_text(encoding='utf-8'))
    literal = dict(mass=1., stiffness=4., radius=1., initial_gap=1.,
                   initial_relative_normal_velocity=-1., friction=0., end_time=2.)
    for name, value in literal.items():
        if inputs['parameters'][name]['value'] != value:
            raise ValueError('frozen literal normal benchmark configuration mismatch: '+name)
    if inputs['dt_cases_s'] != [.005, .0025] or inputs['wall_speed_cases_m_s'] != [0., .2]:
        raise ValueError('frozen normal benchmark case inventory mismatch')
    return dict(classification=CLASSIFICATION, source_configuration=inputs,
        mass_kg=1., radius_m=1., stiffness_N_m=4., mu=0., normal=[0., 0., 1.],
        initial_position_m=[0., 0., 1.], initial_velocity_m_s=[0., 0., wall_speed-1.],
        wall_translation_velocity_m_s=[0., 0., wall_speed], dt_s=dt, end_time_s=2.,
        free_axes=[1., 1., 1.], physical_scaling='NONE_MANUFACTURED_PARAMETERS_NOT_ALUMINUM',
        normal_law='Fn=k*positive_part(radius-gap); U=k*positive_part(radius-gap)^2/2; en=1',
        field_status=dict(position='SOLVED', velocity='SOLVED', normal_potential='SOLVED',
            friction_heat='SOLVED_ZERO_MU', temperature='NOT_SOLVED', stress='NOT_SOLVED',
            density='NOT_SOLVED', material='NOT_SOLVED',
            checkpoint_rho='UPSTREAM_SCRATCH_NOT_SOLVED', checkpoint_e='UPSTREAM_SCRATCH_NOT_SOLVED'))


class NormalImpact:
    def __init__(self, dt, wall_speed=0., pa=None, state=None):
        if dt not in (.005, .0025) or wall_speed not in (0., .2):
            raise ValueError('only preregistered dt .005/.0025 s and wall speed 0/.2 m/s')
        self.dt, self.wall_speed = float(dt), float(wall_speed)
        self.config = configuration(self.dt, self.wall_speed); self.source_hashes = hashes()
        if pa is None:
            pa = get_particle_array(name='solid', x=[0.], y=[0.], z=[1.],
                                    w=[wall_speed-1.], m=[1.], h=[1.])
            for name in ('ax', 'ay', 'az', 'au', 'av', 'aw', 'arho', 'ae', 'e'):
                if name not in pa.properties: pa.add_property(name)
            for name, value in (('particle_id', 101), ('tracer_type', 1), ('as_rs_tag', -1)):
                pa.add_property(name, type='long', data=np.array([value], dtype=np.int64))
            pa.set_output_arrays(list(pa.properties))
        self.pa = pa; self.time = 0.; self.step = 0; self.last_trace = None
        self.cumulative = dict(normal_work_J=0., wall_work_J=0., friction_raw_J=0.,
            constraint_work_J=0., other_work_J=0., normal_impulse_N_s=0.,
            normal_potential_temporal_defect_J=0., contact_kinetic_residual_J=0.)
        self.first_sampled_open_time_s = None; self.ever_compressed = False
        self.initial_kinetic_J = .5*(wall_speed-1.)**2
        if state is not None:
            if state['configuration_sha256'] != config_hash(self.config): raise ValueError('checkpoint config mismatch')
            if state['source_sha256'] != self.source_hashes: raise ValueError('checkpoint source mismatch')
            if state['classification'] != CLASSIFICATION: raise ValueError('checkpoint classification mismatch')
            self.time = state['time_s']; self.step = state['step']
            self.cumulative = dict(state['cumulative'])
            self.ever_compressed = state['ever_compressed']
            self.first_sampled_open_time_s = state['first_sampled_open_time_s']
        self.integrator = HostLeapfrog(self.pa, self._callback)

    @classmethod
    def from_checkpoint(cls, path):
        pa, state = restore(path)
        return cls(state['dt_s'], state['wall_speed_m_s'], pa, state)

    def normal(self, t=None):
        if t is None: t = self.time
        return normal_state(self.pa.z[0]-self.wall_speed*t, 1., 4.)

    def _callback(self, midpoint_time, dt):
        before = velocities(self.pa)[0].copy()
        state = self.normal(midpoint_time)
        wall = np.array([0., 0., self.wall_speed])
        result = kick(before, 1., np.zeros(3), [[0., 0., 1.]],
                      [float(state['normal_force_N'])], [wall], 0., dt)
        acceleration = (result['velocity']-before)/dt
        self.pa.au[:] = acceleration[0]; self.pa.av[:] = acceleration[1]; self.pa.aw[:] = acceleration[2]
        self.last_trace = dict(step=self.step+1, stage='MIDPOINT_FROZEN_NORMAL_KICK',
            time_mid_s=midpoint_time, particle_id=int(self.pa.particle_id[0]),
            contact_id=f'{int(self.pa.particle_id[0])}:NORMAL_PLANE',
            feature_id='NORMAL_PLANE', position_mid_m=positions(self.pa)[0].copy(),
            wall_position_mid_m=[0., 0., self.wall_speed*midpoint_time],
            normal=[0., 0., 1.], wall_velocity_m_s=wall, normal_state=state,
            material_velocity_before_m_s=before, result=result)

    def advance(self):
        if self.time+self.dt > 2.+1e-14: raise ValueError('beyond frozen end time')
        snapshot = {k: a.get_npy_array().copy() for k, a in self.pa.properties.items()}
        prior = self.record(); prior_trace = self.last_trace
        try:
            self.integrator.advance(self.time, self.dt)
            if not np.isfinite(positions(self.pa)).all() or not np.isfinite(velocities(self.pa)).all():
                raise ValueError('nonfinite candidate')
            r = self.last_trace['result']
            next_time = (self.step+1)*self.dt
            after_normal = self.normal(next_time)
            normal_work = sum(float(a['normal_work_J'].sum()) for a in r['intervals'])
            delta_U = float(after_normal['potential_J'])-prior['normal_potential_J']
            defect = normal_work+delta_U-r['wall_work_J']
            delta_K = .5*float(np.sum(velocities(self.pa)**2))-prior['kinetic_J']
            if abs(delta_K-normal_work) > 1e-12:
                raise ValueError('normal impulse kinetic work identity failed')
            if not np.array_equal(positions(self.pa)[0, :2], [0., 0.]):
                raise ValueError('manufactured normal impact moved tangentially')
            self.last_trace.update(position_before_m=prior['position_m'],
                position_after_m=positions(self.pa)[0].copy(),
                normal_state_before=dict(gap_m=prior['gap_m'], active=prior['active'],
                    potential_J=prior['normal_potential_J']), normal_state_after=after_normal,
                delta_potential_J=delta_U, delta_kinetic_J=delta_K,
                normal_potential_temporal_defect_J=defect)
        except Exception:
            self.rejected_candidate = dict(classification='REJECTED', restart_allowed=False,
                time_s=self.time+self.dt, position_m=positions(self.pa), velocity_m_s=velocities(self.pa),
                trace=self.last_trace)
            for k, value in snapshot.items(): getattr(self.pa, k)[:] = value
            self.last_trace = prior_trace
            raise
        for key in ('wall_work_J', 'constraint_work_J', 'other_work_J'):
            self.cumulative[key] += r[key]
        self.cumulative['normal_work_J'] += normal_work
        self.cumulative['normal_potential_temporal_defect_J'] += defect
        self.cumulative['contact_kinetic_residual_J'] += r['kinetic_work_residual_J']
        self.cumulative['friction_raw_J'] += float(np.sum(r['friction_heat_J']))
        self.cumulative['normal_impulse_N_s'] += float(r['normal_impulse'][:, 2].sum())
        self.time = next_time; self.step += 1
        if bool(after_normal['active']): self.ever_compressed = True
        elif self.ever_compressed and self.first_sampled_open_time_s is None:
            self.first_sampled_open_time_s = self.time
        return self.record()

    def record(self):
        state = self.normal(); kinetic = .5*float(np.sum(velocities(self.pa)**2))
        return dict(step=self.step, time_s=self.time, particle_id=int(self.pa.particle_id[0]),
            tracer_type=int(self.pa.tracer_type[0]), as_rs_tag=int(self.pa.as_rs_tag[0]),
            position_m=positions(self.pa)[0].copy(), velocity_m_s=velocities(self.pa)[0].copy(),
            wall_position_m=[0., 0., self.wall_speed*self.time],
            wall_velocity_m_s=[0., 0., self.wall_speed], gap_m=float(self.pa.z[0]-self.wall_speed*self.time),
            active=bool(state['active']), penetration_m=float(state['penetration']),
            normal_force_N=float(state['normal_force_N']), normal_potential_J=float(state['potential_J']),
            kinetic_J=kinetic, cumulative=dict(self.cumulative),
            global_energy_residual_J=kinetic+float(state['potential_J'])-self.initial_kinetic_J-
                self.cumulative['wall_work_J'], first_sampled_open_time_s=self.first_sampled_open_time_s)

    def save_checkpoint(self, path):
        checkpoint(path, self.pa, dict(classification=CLASSIFICATION, accepted=True, restart_allowed=True,
            production_authorized=False, dt_s=self.dt, wall_speed_m_s=self.wall_speed,
            step=self.step, time_s=self.time, cumulative=self.cumulative,
            first_sampled_open_time_s=self.first_sampled_open_time_s, ever_compressed=self.ever_compressed,
            configuration=self.config, configuration_sha256=config_hash(self.config), source_sha256=self.source_hashes))

    def vtk(self, path):
        write_point_cloud(path, positions(self.pa).tolist(),
            dict(particle_id=self.pa.particle_id.tolist(), tracer_type=self.pa.tracer_type.tolist(),
                 as_rs_tag=self.pa.as_rs_tag.tolist(), mass_kg=self.pa.m.tolist()),
            dict(velocity_m_per_s=velocities(self.pa).tolist()),
            title=f'{CLASSIFICATION} time_s={self.time:.17g}; STRESS_THERMAL_NOT_SOLVED')


def run(dt, wall_speed, output, probe_steps=None, timeout_s=50.):
    output = Path(output)
    if output.exists(): raise ValueError('fresh output directory required')
    if probe_steps not in (None, 2): raise ValueError('cost probe must be exactly two steps')
    output.mkdir(parents=True)
    start = time.perf_counter(); block = NormalImpact(dt, wall_speed)
    manifest = dict(classification=CLASSIFICATION, status='RUNNING', source_sha256=block.source_hashes,
        configuration=block.config, configuration_sha256=config_hash(block.config),
        production_authorized=False, maximum_wall_seconds=timeout_s,
        reference_status='Independent harmonic impact values are checked by separate tests/readback, not solver driving')
    records = [block.record()]; traces = []
    block.save_checkpoint(output/'last_accepted.npz'); block.vtk(output/'initial.vtk')
    write_json(output/'manifest.json', manifest)
    count = probe_steps if probe_steps is not None else int(round(2./dt))
    try:
        while block.step < count:
            if time.perf_counter()-start > timeout_s: raise TimeoutError('normal impact wall budget exceeded')
            records.append(block.advance()); traces.append(block.last_trace)
            block.save_checkpoint(output/'last_accepted.npz')
        manifest['status'] = 'PROBE_COMPLETE_NOT_FULL_BENCHMARK' if probe_steps else 'COMPLETE_MANUFACTURED_NORMAL_IMPACT_NOT_AFSD'
    except Exception as error:
        manifest.update(status='BLOCKED', error=str(error))
        if hasattr(block, 'rejected_candidate'): write_json(output/'rejected_candidate.json', block.rejected_candidate)
    finally:
        block.vtk(output/'final.vtk')
        write_json(output/'history.json', records); write_json(output/'normal_intervals.json', traces)
        manifest.update(wall_seconds=time.perf_counter()-start, accepted_steps=block.step,
            accepted_time_s=block.time, final_record=block.record(),
            particle_array_bytes=sum(a.get_npy_array().nbytes for a in block.pa.properties.values()),
            memory_note='ParticleArray property bytes only, NOT peak process RSS')
        manifest['artifact_sha256'] = {str(p.relative_to(output)).replace('\\', '/'):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.rglob('*'))
            if p.is_file() and p != output/'manifest.json'}
        write_json(output/'manifest.json', manifest)
    print(json.dumps(jsonable(manifest), allow_nan=False), flush=True)
    return 0 if manifest['status'] != 'BLOCKED' else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dt', type=float, choices=[.005, .0025], default=.005)
    parser.add_argument('--wall-speed', type=float, choices=[0., .2], default=0.)
    parser.add_argument('--probe-steps', type=int, choices=[2])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    return run(args.dt, args.wall_speed, args.output, args.probe_steps)


if __name__ == '__main__': raise SystemExit(main())

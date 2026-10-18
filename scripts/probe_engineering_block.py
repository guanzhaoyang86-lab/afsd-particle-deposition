"""Small-strain reference-configuration 3D load test, NOT the AFSD solver."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.integrator import HostLeapfrog
from afsd_backend.pysph.material import update
from afsd_backend.pysph.mechanics import positions, velocities, internal_force
from afsd_backend.pysph.output import checkpoint, restore
from afsd_scaffold.vtk_legacy import write_point_cloud
from probe_boundary_sbp import extract_original
from probe_engineering_material import PerfectPlasticBenchmark

CONFIG = ROOT/'configs/benchmarks/engineering_block_e2.json'


def load_config():
    cfg = json.loads(CONFIG.read_text(encoding='utf-8'))
    if cfg['production_authorized'] is not False or cfg['parameters']['poisson']['value'] != 0:
        raise ValueError('only the declared nu=0 benchmark is supported')
    expected = dict(length=.01, rho0=2700, young=70e9, poisson=0, alpha=1e-5,
                    temperature=293.15, yield_stress=1e12, h_ratio=1.3)
    if parameters(cfg) != expected:
        raise ValueError('changed manufactured problem requires new preregistration/reference')
    return cfg


def parameters(cfg):
    return {k: v['value'] for k, v in cfg['parameters'].items()}


def reference(cfg, x, t):
    p = parameters(cfg)
    ell, rho, young, alpha = (p[k] for k in ('length', 'rho0', 'young', 'alpha'))
    c = np.sqrt(young/rho); k = np.pi/(2*ell); omega = c*k
    f = 1-np.cos(omega*t)
    u, v, s = np.zeros_like(x), np.zeros_like(x), np.zeros((len(x), 3, 3))
    u[:, 0] = alpha/k*f*np.sin(k*x[:, 0])
    v[:, 0] = alpha*c*np.sin(omega*t)*np.sin(k*x[:, 0])
    s[:, 0, 0] = young*alpha*f*np.cos(k*x[:, 0])
    energy = young*ell**3*alpha**2/4
    return dict(displacement_m=u, velocity_m_s=v, stress_Pa=s,
                reaction_N=-young*alpha*ell**2*f, work_J=2*energy*f,
                kinetic_J=energy*np.sin(omega*t)**2, elastic_J=energy*f*f)


def make_problem(cfg, n, cfl):
    return Block(cfg, n, cfl)


class Block:
    def __init__(self, cfg, n, cfl):
        if n < 3 or cfl <= 0:
            raise ValueError('n>=3 and positive CFL required')
        self.cfg, self.p, self.n, self.cfl = cfg, parameters(cfg), n, cfl
        p = self.p
        self.dx = p['length']/(n-1)
        index = np.indices((n, n, n)).reshape(3, -1).T
        self.X = self.dx*index.astype(float)
        self.V = self.dx**3*np.prod(np.where((index == 0) | (index == n-1), .5, 1.), axis=1)
        base = dict(x=self.X[:, 0], y=self.X[:, 1], z=self.X[:, 2],
                    m=p['rho0']*self.V, rho=np.full(n**3, p['rho0']), h=np.full(n**3, p['h_ratio']*self.dx))
        self.pa = get_particle_array(name='solid', **base)
        self.ref = get_particle_array(name='reference', **base)
        for key in ('ax', 'ay', 'az', 'au', 'av', 'aw', 'arho', 'ae', 'e', 'ep'):
            if key not in self.pa.properties:
                self.pa.add_property(key)
        self.pa.add_property('particle_id', type='long', data=np.arange(1, n**3+1))
        self.pa.add_property('tracer_type', type='long', data=index[:, 2] % 3)
        self.pa.add_property('S', stride=9)
        self.pa.add_property('temperature', data=np.full(n**3, p['temperature']))
        for axis, key in enumerate(('x0', 'y0', 'z0')):
            self.pa.add_property(key, data=self.X[:, axis])
        self.pa.set_output_arrays(list(self.pa.properties))
        self.D, (self.i, self.j, self.g, self.B), self.condition = extract_original(self.ref)
        # ||sym grad|| <= ||grad||; nu=0 gives elastic modulus E in this bound.
        weighted = self.D*np.sqrt(self.V)[None, :, None]/np.sqrt(self.V)[None, None, :]
        matrix = weighted.reshape(3*n**3, n**3)
        norm_bound = np.sqrt(np.max(np.sum(np.abs(matrix), axis=0))*np.max(np.sum(np.abs(matrix), axis=1)))
        self.c = np.sqrt(p['young']/p['rho0'])
        self.omega_bound = self.c*norm_bound
        self.tf = p['length']/self.c
        requested = min(cfl*self.dx/self.c, 1/self.omega_bound)
        self.steps = int(np.ceil(self.tf/requested))
        self.dt = self.tf/self.steps
        self.fixed = index[:, 0] == 0
        self.body = np.zeros_like(self.X)
        k = np.pi/(2*p['length'])
        self.body[:, 0] = self.pa.m*p['alpha']*self.c**2*k*np.sin(k*self.X[:, 0])
        self.law = PerfectPlasticBenchmark(p['yield_stress'], p['temperature'])
        self.sigma = np.zeros((n**3, 3, 3))
        self.last_displacement = np.zeros_like(self.X)
        self.stress_time = self.time = 0.
        self.count = 0
        self.constraint_impulse = np.zeros(3)
        self.last_midpoint_reaction = np.zeros(3)
        self.constraint_work = self.internal_work = 0.
        self.max_force_defect = self.max_torque_defect = 0.
        self.integrator = HostLeapfrog(self.pa, self.evaluate)

    def gradient(self, field):
        return np.einsum('kij,ja->iak', self.D, field)

    def sync(self, t):
        displacement = positions(self.pa)-self.X
        interval = t-self.stress_time
        if interval <= 0:
            if interval < 0:
                raise ValueError('material time reversal')
            return
        grad = self.gradient(displacement-self.last_displacement)/interval
        stress, ep, plastic, algorithmic, rate = update(
            self.pa.S.reshape(-1, 3, 3), self.pa.ep, grad, self.pa.temperature,
            interval, self.p['young']/2, self.law)
        if np.any(ep != 0) or np.any(plastic != 0) or np.any(algorithmic != 0):
            raise RuntimeError('declared elastic benchmark yielded')
        self.pa.S[:] = stress.ravel(); self.pa.ep[:] = ep
        total = self.gradient(displacement)
        trace = np.trace(total, axis1=1, axis2=2)
        self.sigma = stress+self.p['young']/3*trace[:, None, None]*np.eye(3)
        self.last_displacement = displacement.copy()
        self.stress_time = t

    def force(self):
        return internal_force(self.ref, self.sigma, self.i, self.j, self.g, correction=self.B)

    def evaluate(self, t, dt):
        self.sync(t)
        force = self.force()
        old_v = velocities(self.pa)
        total = force+self.body
        reaction = np.zeros_like(total)
        reaction[self.fixed, 0] = -total[self.fixed, 0]-self.pa.m[self.fixed]*old_v[self.fixed, 0]/dt
        acceleration = (total+reaction)/self.pa.m[:, None]
        mid_v = old_v+.5*dt*acceleration
        self.constraint_work += float(dt*np.sum(reaction*mid_v))
        self.constraint_impulse += dt*reaction.sum(axis=0)
        self.last_midpoint_reaction = reaction.sum(axis=0)
        self.internal_work += float(dt*np.sum(force*mid_v))
        scale = self.p['young']*self.p['alpha']*self.p['length']**2
        self.max_force_defect = max(self.max_force_defect, float(np.linalg.norm(force.sum(axis=0))/scale))
        self.max_torque_defect = max(self.max_torque_defect, float(np.linalg.norm(np.cross(self.X, force).sum(axis=0))/(scale*self.p['length'])))
        self.pa.au[:], self.pa.av[:], self.pa.aw[:] = acceleration.T
        # rho is REFERENCE density in this explicitly linearized test variant.
        self.pa.arho[:] = 0.; self.pa.ae[:] = 0.

    def advance(self):
        self.integrator.advance(self.time, self.dt)
        self.count += 1
        self.time = self.count*self.dt
        self.sync(self.time)
        if not all(np.isfinite(x).all() for x in (positions(self.pa), velocities(self.pa), self.sigma)):
            raise RuntimeError('non-finite candidate')
        if np.max(np.linalg.norm(positions(self.pa)-self.X, axis=1))/self.p['length'] >= 1e-4:
            raise RuntimeError('left preregistered small displacement domain')

    def record(self):
        displacement = positions(self.pa)-self.X
        v = velocities(self.pa)
        strain = self.gradient(displacement)
        tr = np.trace(strain, axis1=1, axis2=2)
        kinetic = float(np.sum(.5*self.pa.m*np.sum(v*v, axis=1)))
        shear = self.pa.S.reshape(-1, 3, 3)
        elastic = float(np.sum(self.V*(np.sum(shear*shear, axis=(1, 2))/(2*self.p['young'])+self.p['young']/6*tr*tr)))
        work = float(np.sum(self.body*displacement))
        reaction = float(-np.sum((self.force()+self.body)[self.fixed, 0]))
        return dict(time_s=self.time, step=self.count, kinetic_J=kinetic, elastic_J=elastic,
                    body_work_J=work, constraint_work_J=self.constraint_work,
                    internal_force_work_J=self.internal_work,
                    energy_residual_J=kinetic+elastic-work-self.constraint_work,
                    end_state_reaction_x_N=reaction, constraint_impulse_N_s=self.constraint_impulse.tolist(),
                    last_kick_midpoint_reaction_N=self.last_midpoint_reaction.tolist(),
                    linear_momentum_kg_m_s=np.sum(self.pa.m[:, None]*v, axis=0).tolist(),
                    body_impulse_N_s=(self.time*self.body.sum(axis=0)).tolist(),
                    mass_kg=float(self.pa.m.sum()), max_displacement_m=float(np.max(np.linalg.norm(displacement, axis=1))))


def weighted_norm(field, volume):
    return float(np.sqrt(np.sum(volume*np.sum(field.reshape(len(volume), -1)**2, axis=1))))


def metrics(block):
    b = block; ref = reference(b.cfg, b.X, b.time)
    u = positions(b.pa)-b.X
    r = b.record()
    u_scale = weighted_norm(ref['displacement_m'], b.V)
    stress_scale = weighted_norm(ref['stress_Pa'], b.V)
    return dict(displacement_relative=weighted_norm(u-ref['displacement_m'], b.V)/u_scale,
                stress_relative=weighted_norm(b.sigma-ref['stress_Pa'], b.V)/stress_scale,
                reaction_relative=abs(r['end_state_reaction_x_N']/ref['reaction_N']-1),
                transverse_relative=weighted_norm(u[:, 1:], b.V)/u_scale,
                energy_relative=abs(r['energy_residual_J'])/ref['work_J'],
                continuum_work_relative=abs(r['body_work_J']/ref['work_J']-1),
                force_relative=b.max_force_defect, torque_relative=b.max_torque_defect,
                mass_relative=abs(b.pa.m.sum()/(b.p['rho0']*b.p['length']**3)-1),
                fixed_displacement_max_m=float(np.max(np.abs(u[b.fixed, 0]))),
                equivalent_plastic_strain_max=float(b.pa.ep.max()))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--cfl', type=float, required=True)
    parser.add_argument('--pilot-steps', type=int)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    cfg = load_config(); start = time.perf_counter()
    b = make_problem(cfg, args.n, args.cfl)
    setup_s = time.perf_counter()-start
    limit = b.steps if args.pilot_steps is None else min(b.steps, args.pilot_steps)
    records = [b.record()]; frames = [dict(time_s=0., position=positions(b.pa).copy(), velocity=velocities(b.pa), stress=b.sigma.copy())]
    run_start = time.perf_counter()
    for _ in range(limit):
        b.advance(); records.append(b.record())
        if b.count % max(1, b.steps//4) == 0 or b.count == limit:
            frames.append(dict(time_s=b.time, position=positions(b.pa).copy(), velocity=velocities(b.pa), stress=b.sigma.copy()))
    advance_s = time.perf_counter()-run_start
    error = metrics(b); tol = cfg['tolerances']
    checks = {k: error[k] <= tol[k] for k in ('mass_relative', 'energy_relative', 'displacement_relative', 'stress_relative', 'reaction_relative', 'transverse_relative')}
    checks.update(force=error['force_relative'] <= tol['force_torque_relative'],
                  torque=error['torque_relative'] <= tol['force_torque_relative'],
                  constraints=error['fixed_displacement_max_m'] == 0 and b.constraint_work == 0,
                  elastic_only=error['equivalent_plastic_strain_max'] == 0,
                  ids=bool(np.array_equal(b.pa.particle_id, np.arange(1, args.n**3+1))))
    checks = {k: bool(v) for k, v in checks.items()}
    status = 'PILOT_NOT_ACCEPTANCE' if args.pilot_steps is not None else ('PASS' if all(checks.values()) else 'FAIL')
    ref = reference(cfg, b.X, b.time)
    np.savez_compressed(args.output/'history.npz', particle_id=b.pa.particle_id.copy(), reference_position_m=b.X,
        volume_m3=b.V, mass_kg=b.pa.m.copy(), time_s=np.array([f['time_s'] for f in frames]),
        position_m=np.array([f['position'] for f in frames]), velocity_m_s=np.array([f['velocity'] for f in frames]),
        stress_Pa=np.array([f['stress'] for f in frames]), equivalent_plastic_strain=b.pa.ep.copy(),
        body_force_N=b.body, reference_final_displacement_m=ref['displacement_m'], reference_final_stress_Pa=ref['stress_Pa'])
    checkpoint(args.output/'checkpoint.npz', b.pa, dict(classification=cfg['classification'], restart_allowed=False,
        reason='diagnostic checkpoint; reference operator/material staging reconstruction not implemented', time_s=b.time,
        production_authorized=False, material_variant='E2_SYNTHETIC_ELASTIC', not_solved=cfg['not_solved']))
    scalar = {'particle_id': b.pa.particle_id.tolist(), 'stress_xx_Pa': b.sigma[:, 0, 0].tolist(),
              'equivalent_plastic_strain': b.pa.ep.tolist(), 'reference_density_kg_m3': b.pa.rho.tolist()}
    write_point_cloud(args.output/'final.vtk', positions(b.pa).tolist(), scalar,
        {'velocity_m_s': velocities(b.pa).tolist(), 'reference_position_m': b.X.tolist()},
        title='E2_SMALL_STRAIN_BENCHMARK; thermal_NOT_SOLVED; NOT_AFSD')
    source_paths = ['scripts/probe_engineering_block.py', 'scripts/probe_boundary_sbp.py',
        'scripts/probe_engineering_material.py', 'src/afsd_backend/pysph/mechanics.py',
        'src/afsd_backend/pysph/material.py', 'src/afsd_backend/pysph/integrator.py',
        'configs/benchmarks/engineering_block_e2.json']
    report = dict(status=status, classification=cfg['classification'], production_authorized=False,
        n=args.n, particles=args.n**3, cfl=args.cfl, dt_s=b.dt, steps_to_common_end=b.steps, actual_steps=b.count,
        common_end_s=b.tf, actual_end_s=b.time, setup_wall_s=setup_s, advance_wall_s=advance_s,
        seconds_per_step=advance_s/limit, predicted_full_wall_s=setup_s+b.steps*advance_s/limit,
        state_array_bytes=sum(p.get_npy_array().nbytes for p in b.pa.properties.values()),
        dense_gradient_bytes=b.D.nbytes, memory_note='array bytes, NOT process peak RSS',
        condition_max=b.condition, stability_omega_upper_per_s=b.omega_bound,
        stability_dt_omega=b.dt*b.omega_bound, metrics=error, checks=checks, records=records,
        not_solved=cfg['not_solved'], checkpoint_restart='NOT_IMPLEMENTED_DIAGNOSTIC_ONLY',
        source_sha256={p: sha(ROOT/p) for p in source_paths},
        artifact_sha256={p.name: sha(p) for p in args.output.iterdir()})
    (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'records'}, indent=2))
    return 0 if status in ('PASS', 'PILOT_NOT_ACCEPTANCE') else 1


if __name__ == '__main__':
    sys.exit(main())

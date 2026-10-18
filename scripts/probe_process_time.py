"""Q1 isolated quasi-static 3D benchmark, NOT a production AFSD integrator."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from probe_engineering_block import Block, load_config as load_e2
from afsd_backend.pysph.material import update
from afsd_backend.pysph.mechanics import internal_force
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG = ROOT/'configs/benchmarks/process_time_q1.json'


def load_config():
    cfg = json.loads(CONFIG.read_text(encoding='utf-8'))
    parent = load_e2()
    if cfg['production_authorized'] is not False:
        raise ValueError('benchmark cannot authorize production')
    for key, entry in parent['parameters'].items():
        if (cfg['parameters'][key]['value'], cfg['parameters'][key]['unit']) != (entry['value'], entry['unit']):
            raise ValueError('Q1 geometry/material changed from preregistered E2: '+key)
    if cfg['parameters']['ramp_duration']['value'] != .1:
        raise ValueError('changed physical load schedule needs explicit new contract')
    return cfg


def continuum_reference(x, factor):
    """Independent PDE reference; never used in the residual or imposed motion."""
    length, young, alpha = .01, 70e9, 1e-5
    k = np.pi/(2*length)
    u = np.zeros_like(x); stress = np.zeros((len(x), 3, 3))
    u[:, 0] = factor*alpha/k*np.sin(k*x[:, 0])
    stress[:, 0, 0] = factor*young*alpha*np.cos(k*x[:, 0])
    return dict(displacement_m=u, stress_Pa=stress, reaction_N=[-70.*factor, 0., 0.],
                elastic_J=1.75e-6*factor*factor)


class QuasiStaticBlock:
    def __init__(self, n=3, cfg=None):
        self.cfg = cfg or load_config()
        if n not in (3, 5, 7):
            raise ValueError('bounded Q1 only supports preregistered n=3/5/7')
        # Reuse existing particles, 3D NNPS, corrected pair and volume quadrature.
        # The E2 HostLeapfrog object is constructed but NEVER advanced here.
        self.mesh = Block(load_e2(), n, .1)
        self.pa, self.X, self.V = self.mesh.pa, self.mesh.X, self.mesh.V
        self.n = n
        self.fixed = self.X == 0.  # corresponding component only, not a rigid block
        self.free = ~self.fixed.reshape(-1)
        self.young = self.mesh.p['young']
        self.force_scale = self.young*self.mesh.p['alpha']*.01**2
        self.K = self.elastic_matrix()
        self.Kff = self.K[np.ix_(self.free, self.free)]
        self.symmetry_relative = float(np.linalg.norm(self.K-self.K.T)/np.linalg.norm(self.K))
        if self.symmetry_relative > 1e-11:
            raise RuntimeError('elastic tangent is not symmetric')
        eig = np.linalg.eigvalsh(self.Kff)
        self.eigen_min, self.eigen_max = float(eig[0]), float(eig[-1])
        if eig[0] <= self.cfg['numerical_solve']['stiffness_positive_ratio']*eig[-1]:
            raise RuntimeError('constrained stiffness singular/indefinite; no pinv or added constraints')
        self.displacement = np.zeros_like(self.X)
        self.stress = np.zeros((len(self.X), 3, 3))
        self.external = np.zeros_like(self.X)
        self.reaction = np.zeros_like(self.X)
        self.rate = np.zeros_like(self.X)
        self.time = 0.; self.count = 0
        self.external_work = self.internal_work = self.constraint_work = 0.
        self.last_iterations = 0; self.last_residual = 0.
        self.last_force_defect = self.last_torque_defect = 0.
        self.last_trace = []

    def elastic_matrix(self):
        D, V = self.mesh.D, self.V
        weighted = [V[:, None]*D[k] for k in range(3)]
        H = sum(D[k].T@weighted[k] for k in range(3))
        count = len(V)
        K = np.empty((count, 3, count, 3))
        for a in range(3):
            for b in range(3):
                K[:, a, :, b] = .5*self.young*((H if a == b else 0.) + D[b].T@weighted[a])
        return K.reshape(3*count, 3*count)

    def load_factor(self, t):
        duration = self.cfg['parameters']['ramp_duration']['value']
        if not np.isfinite(t) or t < 0 or t > duration:
            raise ValueError('time outside declared physical ramp')
        s = t/duration
        return 3*s*s-2*s*s*s

    def trial(self, displacement, t):
        """Pure candidate from accepted material; no mutation/iteration accumulation."""
        u = np.asarray(displacement, dtype=float)
        if u.shape != self.X.shape or not np.isfinite(u).all():
            raise ValueError('invalid candidate displacement')
        if np.any(u[self.fixed] != 0.):
            raise ValueError('support component moved')
        if t <= self.time:
            raise ValueError('candidate must advance physical time')
        factor = self.load_factor(t)
        dt = t-self.time
        gradient = self.mesh.gradient(u-self.displacement)/dt
        S, ep, plastic, returned, strain_rate = update(
            self.pa.S.reshape(-1, 3, 3), self.pa.ep, gradient,
            self.pa.temperature, dt, self.young/2, self.mesh.law)
        if np.any(ep != 0) or np.any(plastic != 0) or np.any(returned != 0):
            raise RuntimeError('Q1 is declared elastic; yielded candidate rejected')
        trace = np.trace(self.mesh.gradient(u), axis1=1, axis2=2)
        stress = S+self.young/3*trace[:, None, None]*np.eye(3)
        force = internal_force(self.mesh.ref, stress, self.mesh.i, self.mesh.j,
                               self.mesh.g, correction=self.mesh.B)
        external = factor*self.mesh.body
        residual = force+external
        return dict(base_step=self.count, time_s=float(t), dt_s=dt, factor=factor,
                    displacement=u.copy(), S=S, ep=ep, stress=stress, internal_force=force,
                    external=external, residual=residual, strain_rate=strain_rate)

    def solve_candidate(self, t, maximum_iterations=None):
        limit = (self.cfg['numerical_solve']['maximum_iterations']
                 if maximum_iterations is None else maximum_iterations)
        if limit < 0:
            raise ValueError('negative iteration limit')
        u = self.displacement.copy(); trace = []
        for iteration in range(limit+1):
            candidate = self.trial(u, t)
            relative = float(np.linalg.norm(candidate['residual'].ravel()[self.free])/self.force_scale)
            if not np.isfinite(relative):
                raise RuntimeError('nonfinite equilibrium residual')
            trace.append(dict(iteration=iteration, free_force_relative=relative))
            if relative <= self.cfg['numerical_solve']['force_relative']:
                candidate.update(iterations=iteration, iteration_trace=trace,
                                 free_force_relative=relative)
                return candidate
            if iteration != limit:
                step = np.linalg.solve(self.Kff, candidate['residual'].ravel()[self.free])
                u.ravel()[self.free] += step
        raise RuntimeError('bounded equilibrium did not converge; accepted state unchanged; '+str(trace))

    def commit(self, candidate):
        if candidate['base_step'] != self.count:
            raise ValueError('stale candidate from different accepted state')
        # Re-evaluate at commit: mutable candidate fields cannot bypass physics.
        trial = self.trial(candidate['displacement'], candidate['time_s'])
        residual = float(np.linalg.norm(trial['residual'].ravel()[self.free])/self.force_scale)
        if residual > self.cfg['numerical_solve']['force_relative'] or not np.isfinite(residual):
            raise ValueError('unconverged candidate cannot be accepted')
        if np.max(np.linalg.norm(trial['displacement'], axis=1))/.01 >= self.cfg['numerical_solve']['small_displacement_ratio']:
            raise ValueError('left small-strain benchmark domain')
        du = trial['displacement']-self.displacement
        reaction = np.zeros_like(self.X)
        reaction[self.fixed] = -trial['residual'][self.fixed]
        external_work = self.external_work+float(np.sum(.5*(self.external+trial['external'])*du))
        internal_work = self.internal_work+float(np.einsum('i,iab,iab->', self.V,
                                           .5*(self.stress+trial['stress']), self.mesh.gradient(du)))
        constraint_work = self.constraint_work+float(np.sum(.5*(self.reaction+reaction)*du))
        if abs(external_work+constraint_work-internal_work) > 1e-7*1.75e-6:
            raise RuntimeError('discrete work identity failed; no residual heat')
        # All checks completed. Commit one time; iterations did not touch these.
        self.displacement = trial['displacement'].copy()
        self.stress, self.external, self.reaction = trial['stress'], trial['external'], reaction
        self.pa.x[:], self.pa.y[:], self.pa.z[:] = (self.X+self.displacement).T
        self.pa.S[:] = trial['S'].ravel(); self.pa.ep[:] = trial['ep']
        self.rate = du/trial['dt_s']  # kinematic interval average; not dynamic velocity
        self.time = trial['time_s']; self.count += 1
        self.external_work, self.internal_work, self.constraint_work = external_work, internal_work, constraint_work
        self.last_iterations = candidate.get('iterations', None)
        self.last_trace = candidate.get('iteration_trace', [])
        self.last_residual = residual
        force = trial['internal_force']
        self.last_force_defect = float(np.linalg.norm(force.sum(axis=0))/self.force_scale)
        self.last_torque_defect = float(np.linalg.norm(np.cross(self.X, force).sum(axis=0))/(self.force_scale*.01))

    def advance(self, t):
        self.commit(self.solve_candidate(t))

    def record(self):
        grad = self.mesh.gradient(self.displacement)
        strain = .5*(grad+grad.swapaxes(1, 2))
        elastic = float(.5*self.young*np.einsum('i,iab,iab->', self.V, strain, strain))
        return dict(step=self.count, physical_time_s=self.time, load_factor=self.load_factor(self.time),
            corrections=self.last_iterations, iteration_trace=self.last_trace,
            free_force_relative=self.last_residual, external_force_N=self.external.sum(axis=0).tolist(),
            reaction_N=self.reaction.sum(axis=0).tolist(), external_work_J=self.external_work,
            stress_work_J=self.internal_work, constraint_work_J=self.constraint_work,
            elastic_small_strain_diagnostic_J=elastic,
            work_identity_residual_J=self.external_work+self.constraint_work-self.internal_work,
            elastic_diagnostic_residual_J=self.external_work+self.constraint_work-elastic,
            internal_force_relative=self.last_force_defect, internal_torque_relative=self.last_torque_defect,
            mass_kg=float(self.pa.m.sum()), maximum_ep=float(self.pa.ep.max()))


def accuracy(block):
    b = block; ref = continuum_reference(b.X, b.load_factor(b.time))
    def norm(x):
        return np.sqrt(np.sum(b.V*np.sum(x.reshape(len(b.V), -1)**2, axis=1)))
    return dict(displacement_relative=float(norm(b.displacement-ref['displacement_m'])/norm(ref['displacement_m'])),
                stress_relative=float(norm(b.stress-ref['stress_Pa'])/norm(ref['stress_Pa'])),
                reaction_relative=float(abs(b.reaction[:, 0].sum()/ref['reaction_N'][0]-1)),
                elastic_energy_relative=float(abs(b.record()['elastic_small_strain_diagnostic_J']/ref['elastic_J']-1)),
                policy='REPORTED_NOT_A_PRECISION_GATE')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--increments', type=int, choices=(4, 8), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    cfg = load_config()
    config_path = args.output/'resolved_config.json'
    config_path.write_text(json.dumps(cfg, indent=2)+'\n', encoding='utf-8')
    b = None; error = None; records = []; frames = []
    def capture():
        records.append(b.record())
        frames.append(dict(position=b.X+b.displacement, displacement=b.displacement.copy(),
                           stress=b.stress.copy(), reaction=b.reaction.copy(), external=b.external.copy(),
                           kinematic_rate=b.rate.copy()))
    try:
        b = QuasiStaticBlock(args.n, cfg)
        setup_s = time.perf_counter()-started
        capture(); tick = time.perf_counter()
        for t in np.linspace(0., .1, args.increments+1)[1:]:
            b.advance(float(t)); capture()
        advance_s = time.perf_counter()-tick
        status = 'DISCRETE_EQUILIBRIUM_COMPLETE'
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        error = str(exc); status = 'BLOCKED'
        setup_s = time.perf_counter()-started; advance_s = None
    report = dict(status=status, classification=cfg['classification'], production_authorized=False,
                  error=error, n=args.n, increments_requested=args.increments,
                  setup_wall_s=setup_s, advance_wall_s=advance_s,
                  wall_seconds=time.perf_counter()-started, not_solved=cfg['not_solved'],
                  accepted_steps=0 if b is None else b.count, accepted_time_s=0 if b is None else b.time,
                  records=records, restart_allowed=False)
    if b is not None:
        report.update(particles=len(b.X), free_dofs=int(b.free.sum()), condition_M_max=b.mesh.condition,
            tangent_symmetry_relative=b.symmetry_relative, constrained_eigenvalue_min_N_per_m=b.eigen_min,
            constrained_eigenvalue_max_N_per_m=b.eigen_max,
            array_bytes=sum(prop.get_npy_array().nbytes for prop in b.pa.properties.values())+b.K.nbytes+b.Kff.nbytes+b.mesh.D.nbytes,
            memory_note='selected arrays only, NOT process peak RSS',
            continuum_accuracy=accuracy(b) if b.time else None)
        if frames:
            np.savez_compressed(args.output/'history.npz', particle_id=b.pa.particle_id.copy(),
                tracer_type=b.pa.tracer_type.copy(), reference_position_m=b.X, mass_kg=b.pa.m.copy(),
                reference_volume_m3=b.V, constrained_components=b.fixed,
                physical_time_s=np.array([r['physical_time_s'] for r in records]),
                position_m=np.array([f['position'] for f in frames]),
                displacement_m=np.array([f['displacement'] for f in frames]),
                stress_Pa=np.array([f['stress'] for f in frames]),
                external_force_N=np.array([f['external'] for f in frames]),
                support_reaction_N=np.array([f['reaction'] for f in frames]),
                kinematic_interval_rate_m_s=np.array([f['kinematic_rate'] for f in frames]),
                equivalent_plastic_strain=b.pa.ep.copy())
            checkpoint(args.output/'checkpoint.npz', b.pa, dict(classification=cfg['classification'],
                restart_allowed=False, reason='diagnostic only; quasi-static accepted reference state not reconstructed',
                physical_time_s=b.time, step=b.count, production_authorized=False, not_solved=cfg['not_solved']))
            scalars = dict(particle_id=b.pa.particle_id.tolist(), tracer_type=b.pa.tracer_type.tolist(),
                           equivalent_plastic_strain=b.pa.ep.tolist(), reference_density_kg_m3=b.pa.rho.tolist())
            for a, c, name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
                scalars['stress_'+name+'_Pa'] = b.stress[:, a, c].tolist()
            write_point_cloud(args.output/'final.vtk', (b.X+b.displacement).tolist(), scalars,
                dict(reference_position_m=b.X.tolist(), displacement_m=b.displacement.tolist()),
                title='Q1_QUASISTATIC_BENCHMARK; NOT_AFSD; inertia_thermal_contact_NOT_SOLVED')
    files = ['scripts/probe_process_time.py', 'configs/benchmarks/process_time_q1.json',
             'scripts/probe_engineering_block.py', 'scripts/probe_boundary_sbp.py',
             'scripts/probe_engineering_material.py', 'src/afsd_backend/pysph/mechanics.py',
             'src/afsd_backend/pysph/material.py', 'src/afsd_backend/pysph/output.py']
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    report['source_sha256'] = {f: sha(ROOT/f) for f in files}
    report['artifact_sha256'] = {f.name: sha(f) for f in args.output.iterdir() if f.is_file()}
    (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if status == 'DISCRETE_EQUILIBRIUM_COMPLETE' else 2


if __name__ == '__main__':
    raise SystemExit(main())

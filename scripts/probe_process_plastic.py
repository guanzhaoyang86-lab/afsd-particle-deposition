"""Q2 displacement-controlled plastic compression/return; NOT production AFSD."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from probe_process_time import QuasiStaticBlock
from probe_engineering_material import PerfectPlasticBenchmark
from afsd_backend.pysph.material import update, FKS
from afsd_backend.pysph.mechanics import internal_force
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG = ROOT/'configs/benchmarks/process_plastic_q2.json'


class CandidateFailure(RuntimeError):
    def __init__(self, message, candidate=None, iterations=None):
        super().__init__(message)
        self.candidate = candidate
        self.iterations = iterations or []


class PlasticBlock(QuasiStaticBlock):
    def __init__(self, n=3, law='j2'):
        if n not in (3, 5) or law not in ('j2', 'fks'):
            raise ValueError('bounded Q2 supports n3/5 and explicit j2/fks only')
        super().__init__(n)
        self.qcfg = json.loads(CONFIG.read_text(encoding='utf-8'))
        if self.qcfg['production_authorized']:
            raise ValueError('Q2 never authorizes production')
        self.law_name = law
        p = self.qcfg['parameters']
        self.pa.temperature[:] = p['j2_temperature' if law == 'j2' else 'fks_temperature']['value']
        self.law = (PerfectPlasticBenchmark(p['j2_yield']['value'], p['j2_temperature']['value'])
                    if law == 'j2' else FKS())
        self.cap = self.X[:, 0] == .01
        self.fixed[self.cap, 0] = True
        self.free = ~self.fixed.ravel()
        self.force_scale = 90e6*.01**2
        self.shear, self.bulk = self.young/2, self.young/3
        for name in ('plastic_work', 'return_dissipation', 'strain_rate'):
            self.pa.add_property(name)
        self.pa.set_output_arrays(list(self.pa.properties))
        # Full tensor gradient map; no engineering-shear factor ambiguity.
        count = len(self.X)
        self.Lmap = np.zeros((count, 3, 3, count, 3))
        for a in range(3):
            for k in range(3):
                self.Lmap[:, a, k, :, a] = self.mesh.D[k]
        self.Lmap = self.Lmap.reshape(count, 9, 3*count)
        self.Lflat = self.Lmap.reshape(9*count, 3*count)
        self.accepted_floor_count = 0

    def cap_displacement(self, t):
        duration = self.qcfg['parameters']['phase_duration']['value']
        if not np.isfinite(t) or t < 0 or t > 2*duration:
            raise ValueError('outside declared physical load-return schedule')
        fraction = t/duration if t <= duration else 2-t/duration
        return -.01*self.qcfg['parameters']['peak_compressive_strain']['value']*fraction

    def prescribed(self, t):
        boundary = np.zeros_like(self.X)
        boundary[self.cap, 0] = self.cap_displacement(t)
        return boundary

    def constitutive(self, increment_gradient, dt):
        S, ep, plastic, returned, rate = update(
            self.pa.S.reshape(-1, 3, 3), self.pa.ep,
            increment_gradient/dt, self.pa.temperature, dt, self.shear, self.law)
        total_trace = np.trace(self.mesh.gradient(self.displacement)+increment_gradient, axis1=1, axis2=2)
        stress = S+self.bulk*total_trace[:, None, None]*np.eye(3)
        arrays = (S, ep, plastic, returned, rate, stress)
        if not all(np.isfinite(a).all() for a in arrays):
            raise ValueError('nonfinite material candidate')
        if np.any(ep < self.pa.ep) or np.any(plastic < 0) or np.any(returned < 0):
            raise ValueError('raw negative plastic increment/work; no clipping')
        return S, ep, plastic, returned, rate, stress

    def trial(self, displacement, t):
        u = np.asarray(displacement, dtype=float)
        prescribed = self.prescribed(t)
        if u.shape != self.X.shape or not np.isfinite(u).all():
            raise ValueError('invalid candidate displacement')
        if np.any(u[self.fixed] != prescribed[self.fixed]):
            raise ValueError('candidate violates support/actuator displacement')
        if not t > self.time:
            raise ValueError('candidate physical time must advance')
        gradient = self.mesh.gradient(u-self.displacement)
        S, ep, plastic, returned, rate, stress = self.constitutive(gradient, t-self.time)
        force = internal_force(self.mesh.ref, stress, self.mesh.i, self.mesh.j, self.mesh.g, correction=self.mesh.B)
        if not np.isfinite(force).all():
            raise ValueError('nonfinite internal force')
        return dict(base_step=self.count, time_s=float(t), dt_s=t-self.time,
            displacement=u.copy(), increment_gradient=gradient, S=S, ep=ep,
            plastic_increment_J=self.V*plastic, return_increment_J=self.V*returned,
            strain_rate=rate, stress=stress, internal_force=force, residual=force)

    def tangent(self, candidate):
        """Local material finite differences; actual paired force remains residual."""
        grad, dt = candidate['increment_gradient'], candidate['dt_s']
        delta = self.qcfg['numerical']['local_gradient_difference']
        C = np.empty((len(self.X), 9, 9))
        for column in range(9):
            a, k = divmod(column, 3)
            plus, minus = grad.copy(), grad.copy()
            plus[:, a, k] += delta; minus[:, a, k] -= delta
            stress_plus = self.constitutive(plus, dt)[-1]
            stress_minus = self.constitutive(minus, dt)[-1]
            C[:, :, column] = ((stress_plus-stress_minus)/(2*delta)).reshape(len(self.X), 9)
        tangent = self.Lflat.T@(self.V[:, None, None]*(C@self.Lmap)).reshape(9*len(self.X), -1)
        if not np.isfinite(tangent).all():
            raise ValueError('nonfinite tangent')
        return tangent

    def solve_candidate(self, t, maximum_iterations=None):
        limits = self.qcfg['numerical']
        maximum = limits['maximum_iterations'] if maximum_iterations is None else maximum_iterations
        if maximum < 0:
            raise ValueError('negative iteration count')
        u = self.displacement.copy()
        u[self.fixed] = self.prescribed(t)[self.fixed]
        trace = []
        for iteration in range(maximum+1):
            candidate = self.trial(u, t)
            residual = candidate['residual'].ravel()[self.free]
            norm = float(np.linalg.norm(residual)/self.force_scale)
            trace.append(dict(iteration=iteration, free_force_relative=norm))
            if norm <= limits['force_relative']:
                candidate.update(iterations=iteration, iteration_trace=trace, free_force_relative=norm)
                return candidate
            if iteration == maximum:
                break
            Kff = self.tangent(candidate)[np.ix_(self.free, self.free)]
            condition = float(np.linalg.cond(Kff))
            trace[-1]['tangent_condition'] = condition
            if not np.isfinite(condition) or condition > limits['tangent_condition_limit']:
                raise CandidateFailure('singular/ill-conditioned plastic tangent; no pinv', candidate, trace)
            try:
                step = np.linalg.solve(Kff, residual)
            except np.linalg.LinAlgError as exc:
                raise CandidateFailure('plastic tangent solve failed', candidate, trace) from exc
            accepted_search = False
            for halving in range(limits['line_search_halvings']+1):
                test_u = u.copy(); test_u.ravel()[self.free] += (2.**(-halving))*step
                probe = self.trial(test_u, t)
                after = float(np.linalg.norm(probe['residual'].ravel()[self.free])/self.force_scale)
                if after < norm:
                    u = test_u; trace[-1]['line_search_factor'] = 2.**(-halving)
                    accepted_search = True; break
            if not accepted_search:
                raise CandidateFailure('bounded line search failed', candidate, trace)
        raise CandidateFailure('bounded plastic equilibrium did not converge', candidate, trace)

    def commit(self, candidate):
        if candidate['base_step'] != self.count:
            raise ValueError('stale candidate')
        c = self.trial(candidate['displacement'], candidate['time_s'])
        residual = float(np.linalg.norm(c['residual'].ravel()[self.free])/self.force_scale)
        if residual > self.qcfg['numerical']['force_relative']:
            raise ValueError('unconverged candidate rejected')
        if np.max(np.linalg.norm(c['displacement'], axis=1))/.01 >= self.qcfg['numerical']['displacement_over_length_limit']:
            raise ValueError('left declared small-strain benchmark domain')
        reaction = np.zeros_like(self.X); reaction[self.fixed] = -c['internal_force'][self.fixed]
        du = c['displacement']-self.displacement
        W = self.constraint_work+float(np.sum(.5*(self.reaction+reaction)*du))
        Ws = self.internal_work+float(np.einsum('i,iab,iab->', self.V,
                          .5*(self.stress+c['stress']), self.mesh.gradient(du)))
        if not np.isfinite([W, Ws]).all() or abs(W-Ws) > 1e-6:
            raise CandidateFailure('actuator/internal work identity failed; no residual heat', c)
        # One commit after checks: rejected iterations cannot accumulate plasticity.
        self.displacement = c['displacement'].copy(); self.stress = c['stress']; self.reaction = reaction
        self.pa.x[:], self.pa.y[:], self.pa.z[:] = (self.X+self.displacement).T
        self.pa.S[:] = c['S'].ravel(); self.pa.ep[:] = c['ep']
        self.pa.plastic_work[:] += c['plastic_increment_J']
        self.pa.return_dissipation[:] += c['return_increment_J']
        self.pa.strain_rate[:] = c['strain_rate']
        self.rate = du/c['dt_s']
        self.time = c['time_s']; self.count += 1
        self.constraint_work, self.internal_work = W, Ws
        self.last_iterations = candidate.get('iterations', None)
        self.last_trace = candidate.get('iteration_trace', [])
        self.last_residual = residual
        force = c['internal_force']
        self.last_force_defect = float(np.linalg.norm(force.sum(axis=0))/self.force_scale)
        self.last_torque_defect = float(np.linalg.norm(np.cross(self.X, force).sum(axis=0))/(self.force_scale*.01))
        if self.law_name == 'fks':
            self.accepted_floor_count += int(np.count_nonzero(c['strain_rate'] < self.law.rate_floor))

    def record(self):
        S = self.pa.S.reshape(-1, 3, 3)
        tr = np.trace(self.mesh.gradient(self.displacement), axis1=1, axis2=2)
        elastic = float(np.sum(self.V*(np.sum(S*S, axis=(1,2))/(4*self.shear)+.5*self.bulk*tr*tr)))
        P = float(self.pa.plastic_work.sum())
        return dict(step=self.count, physical_time_s=self.time, cap_displacement_m=self.cap_displacement(self.time),
            cap_reaction_x_N=float(self.reaction[self.cap, 0].sum()), reaction_N=self.reaction.sum(axis=0).tolist(),
            actuator_work_J=self.constraint_work, stress_work_J=self.internal_work,
            plastic_work_J=P, return_dissipation_diagnostic_J=float(self.pa.return_dissipation.sum()),
            elastic_diagnostic_J=elastic, work_identity_residual_J=self.constraint_work-self.internal_work,
            work_minus_elastic_minus_plastic_J=self.internal_work-elastic-P,
            maximum_ep=float(self.pa.ep.max()), volume_mean_ep=float(np.sum(self.V*self.pa.ep)/self.V.sum()),
            maximum_strain_rate_per_s=float(self.pa.strain_rate.max()),
            free_force_relative=self.last_residual, internal_force_relative=self.last_force_defect,
            internal_torque_relative=self.last_torque_defect, corrections=self.last_iterations,
            iteration_trace=self.last_trace, mass_kg=float(self.pa.m.sum()),
            accepted_rate_floor_evaluations=self.accepted_floor_count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--law', choices=('j2','fks'), required=True)
    parser.add_argument('--n', type=int, choices=(3,5), required=True)
    parser.add_argument('--increments-per-phase', type=int, choices=(8,16), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter(); b = PlasticBlock(args.n, args.law)
    setup = time.perf_counter()-start
    records, frames = [], []
    def capture():
        records.append(b.record())
        frames.append(dict(position=b.X+b.displacement, u=b.displacement.copy(), stress=b.stress.copy(),
            ep=b.pa.ep.copy(), plastic=b.pa.plastic_work.copy(), returned=b.pa.return_dissipation.copy(),
            reaction=b.reaction.copy(), rate=b.rate.copy(), strain_rate=b.pa.strain_rate.copy()))
    capture(); status='DISCRETE_PLASTIC_PATH_COMPLETE'; error=None; rejected=None
    tick=time.perf_counter()
    try:
        duration=b.qcfg['parameters']['phase_duration']['value']
        for t in np.linspace(0.,2*duration,2*args.increments_per_phase+1)[1:]:
            b.advance(float(t)); capture()
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        status='BLOCKED'; error=str(exc)
        rejected=getattr(exc,'candidate',None)
        detail=dict(classification='REJECTED',restart_allowed=False,error=error,
                    attempted_time_s=float(t),last_accepted_step=b.count,
                    iterations=getattr(exc,'iterations',[]))
        if rejected is not None:
            np.savez_compressed(args.output/'REJECTED_candidate.npz',
                **{k:v for k,v in rejected.items() if isinstance(v,np.ndarray)})
        (args.output/'REJECTED_manifest.json').write_text(json.dumps(detail,indent=2)+'\n',encoding='utf-8')
    advance=time.perf_counter()-tick
    cfg={**b.qcfg,'selected_law':args.law,'n':args.n,'increments_per_phase':args.increments_per_phase,
         'elastic_geometry_inputs':b.mesh.cfg,'fks_parameters':vars(b.law) if args.law=='fks' else None}
    (args.output/'resolved_config.json').write_text(json.dumps(cfg,indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(args.output/'history.npz', particle_id=b.pa.particle_id.copy(),tracer_type=b.pa.tracer_type.copy(),
        reference_position_m=b.X,reference_volume_m3=b.V,mass_kg=b.pa.m.copy(),constrained_components=b.fixed,
        cap_particles=b.cap,physical_time_s=np.array([r['physical_time_s'] for r in records]),
        position_m=np.array([f['position'] for f in frames]), displacement_m=np.array([f['u'] for f in frames]),
        stress_Pa=np.array([f['stress'] for f in frames]),equivalent_plastic_strain=np.array([f['ep'] for f in frames]),
        plastic_work_J=np.array([f['plastic'] for f in frames]),return_dissipation_diagnostic_J=np.array([f['returned'] for f in frames]),
        reaction_N=np.array([f['reaction'] for f in frames]),kinematic_interval_rate_m_s=np.array([f['rate'] for f in frames]),
        deviatoric_strain_rate_per_s=np.array([f['strain_rate'] for f in frames]))
    checkpoint(args.output/'checkpoint.npz',b.pa,dict(classification=b.qcfg['classification'],restart_allowed=False,
        physical_time_s=b.time,step=b.count,selected_law=args.law,production_authorized=False,not_solved=b.qcfg['not_solved']))
    scalar=dict(particle_id=b.pa.particle_id.tolist(),tracer_type=b.pa.tracer_type.tolist(),
                equivalent_plastic_strain=b.pa.ep.tolist(),plastic_work_J=b.pa.plastic_work.tolist(),
                return_dissipation_diagnostic_J=b.pa.return_dissipation.tolist(),reference_density_kg_m3=b.pa.rho.tolist())
    for a,c,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        scalar['stress_'+name+'_Pa']=b.stress[:,a,c].tolist()
    write_point_cloud(args.output/'final.vtk',(b.X+b.displacement).tolist(),scalar,
        dict(reference_position_m=b.X.tolist(),displacement_m=b.displacement.tolist()),
        title='Q2_ACCEPTED_PLASTIC_BENCHMARK; NOT_AFSD; thermal_inertia_contact_NOT_SOLVED')
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    source=['scripts/probe_process_plastic.py','configs/benchmarks/process_plastic_q2.json',
            'scripts/probe_process_time.py','scripts/probe_engineering_material.py',
            'scripts/probe_engineering_block.py','scripts/probe_boundary_sbp.py',
            'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']
    report=dict(status=status,error=error,classification=b.qcfg['classification'],production_authorized=False,
        selected_law=args.law,n=args.n,particles=len(b.X),free_dofs=int(b.free.sum()),
        increments_per_phase=args.increments_per_phase,accepted_steps=b.count,accepted_time_s=b.time,
        requested_end_s=.2,setup_wall_s=setup,advance_wall_s=advance,wall_before_report_s=time.perf_counter()-start,
        selected_array_bytes=sum(prop.get_npy_array().nbytes for prop in b.pa.properties.values())+b.mesh.D.nbytes+b.Lmap.nbytes+b.K.nbytes,
        memory_note='selected accepted/operator arrays only; excludes transient tangent/solver buffers and process RSS',
        condition_M_max=b.mesh.condition,records=records,not_solved=b.qcfg['not_solved'],
        source_sha256={name:sha(ROOT/name) for name in source},
        artifact_sha256={f.name:sha(f) for f in args.output.iterdir() if f.is_file()})
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2,allow_nan=False))
    print(json.dumps(records[-1],indent=2,allow_nan=False))
    return 0 if status=='DISCRETE_PLASTIC_PATH_COMPLETE' else 2


if __name__=='__main__':raise SystemExit(main())

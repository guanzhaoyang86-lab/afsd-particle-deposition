"""Q3 single frictionless platen; isolated small-strain benchmark, NOT AFSD."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from probe_process_plastic import PlasticBlock, CandidateFailure, ROOT
from probe_engineering_material import PerfectPlasticBenchmark
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG = ROOT/'configs/benchmarks/process_contact_q3.json'


def reduced_correction(tangent, force, unknown, condition_limit=1e12):
    """Whole constrained stiffness response, never a scalar particle mobility."""
    unknown = np.asarray(unknown, dtype=bool)
    step = np.zeros_like(force, dtype=float)
    if not unknown.any():
        return step, 1.
    sub = tangent[np.ix_(unknown, unknown)]
    condition = float(np.linalg.cond(sub))
    if not np.isfinite(condition) or condition > condition_limit:
        raise ValueError('singular/ill-conditioned constrained tangent; no pinv')
    step[unknown] = np.linalg.solve(sub, force[unknown])
    return step, condition


class ContactBlock(PlasticBlock):
    def __init__(self, law='j2', n=3):
        if law not in ('elastic', 'j2') or n != 3:
            raise ValueError('Q3 bounded single-face n3 elastic/j2 only')
        super().__init__(n, 'j2')
        self.ccfg = json.loads(CONFIG.read_text(encoding='utf-8'))
        used = dict(length=self.mesh.p['length'], peak_wall_displacement=self.cap_displacement(.1),
            phase_duration=self.qcfg['parameters']['phase_duration']['value'],
            j2_yield=self.qcfg['parameters']['j2_yield']['value'], temperature=float(self.pa.temperature[0]))
        if self.ccfg['production_authorized'] is not False or self.ccfg['normal'] != [-1,0,0]:
            raise ValueError('Q3 production/normal scope mismatch')
        for name, value in used.items():
            if self.ccfg['parameters'][name]['value'] != value:
                raise ValueError('Q3 declared/used parameter mismatch: '+name)
        self.law_name = law
        if law == 'elastic':
            self.law = PerfectPlasticBenchmark(self.ccfg['parameters']['elastic_yield']['value'], 293.15)
        # Only permanent rollers; contact is NOT a permanent displacement BC.
        self.fixed = self.X == 0.
        self.free = ~self.fixed.ravel()
        self.nodes = np.flatnonzero(self.cap)
        self.contact_dofs = 3*self.nodes
        self.active = np.zeros(len(self.nodes), dtype=bool)
        self.normal_force = np.zeros(len(self.nodes))
        self.contact_force = np.zeros_like(self.X)
        self.support_force = np.zeros_like(self.X)
        self.wall_work = self.contact_work = self.support_work = 0.
        self.last_wall_increment = self.last_contact_increment = 0.
        self.events = []
        self.last_event = 'INITIAL'

    def prescribed(self, t):
        self.cap_displacement(t)  # validate the same unchanged physical program
        return np.zeros_like(self.X)

    def _equilibrate(self, t, active, maximum_iterations=None):
        """Pure branch solve; tensile active reactions are rejected by outer KKT."""
        active = np.asarray(active, dtype=bool)
        if active.shape != self.active.shape:
            raise ValueError('bad contact inventory')
        unknown = self.free.copy(); unknown[self.contact_dofs[active]] = False
        u = self.displacement.copy()
        u[self.nodes[active], 0] = self.cap_displacement(t)
        maximum = self.qcfg['numerical']['maximum_iterations'] if maximum_iterations is None else maximum_iterations
        if maximum < 0:
            raise ValueError('negative iteration limit')
        trace = []
        for iteration in range(maximum+1):
            c = self.trial(u, t)
            force = c['internal_force'].ravel()
            norm = float(np.linalg.norm(force[unknown])/self.force_scale)
            trace.append(dict(iteration=iteration, free_force_relative=norm))
            if norm <= self.ccfg['numerical']['force_relative']:
                c.update(active=active.copy(), iterations=iteration, iteration_trace=trace)
                return self.contact_state(c)
            if iteration == maximum:
                break
            try:
                step, condition = reduced_correction(self.tangent(c), force, unknown,
                    self.ccfg['numerical']['tangent_condition_limit'])
            except (ValueError, np.linalg.LinAlgError) as exc:
                raise CandidateFailure(str(exc), c, trace) from exc
            trace[-1]['tangent_condition'] = condition
            for halving in range(self.qcfg['numerical']['line_search_halvings']+1):
                test_u = u+(2.**(-halving))*step.reshape(u.shape)
                probe = self.trial(test_u, t)
                after = np.linalg.norm(probe['internal_force'].ravel()[unknown])/self.force_scale
                if after < norm:
                    u = test_u; trace[-1]['line_search_factor'] = 2.**(-halving)
                    break
            else:
                raise CandidateFailure('contact branch line search failed', c, trace)
        raise CandidateFailure('contact branch equilibrium did not converge', c, trace)

    def contact_state(self, c):
        active = np.asarray(c['active'], dtype=bool)
        if active.shape != self.active.shape:
            raise ValueError('bad contact inventory')
        normal = np.zeros(len(self.nodes))
        normal[active] = c['internal_force'][self.nodes[active], 0]
        contact = np.zeros_like(self.X); contact[self.nodes, 0] = -normal
        support = np.zeros_like(self.X)
        support[self.fixed] = -(c['internal_force']+contact)[self.fixed]
        gap = self.cap_displacement(c['time_s'])-c['displacement'][self.nodes, 0]
        residual = c['internal_force']+contact+support
        c.update(normal_force=normal, contact_force=contact, support_force=support, gap=gap,
                 equilibrium_relative=float(np.linalg.norm(residual)/self.force_scale))
        return c

    def solve_candidate(self, t, maximum_iterations=None):
        lim = self.ccfg['numerical']
        gap_guess = self.cap_displacement(t)-self.displacement[self.nodes, 0]
        active = ((self.active & (self.normal_force > lim['contact_force_tolerance_N'])) |
                  (gap_guess < -lim['gap_tolerance_m']))
        inventories = []
        for outer in range(lim['maximum_active_sets']):
            c = self._equilibrate(t, active, maximum_iterations)
            inventories.append(dict(active=active.tolist(), iterations=c['iteration_trace']))
            tensile = active & (c['normal_force'] < -lim['contact_force_tolerance_N'])
            penetrating = (~active) & (c['gap'] < -lim['gap_tolerance_m'])
            if not tensile.any() and not penetrating.any():
                c.update(active_set_trace=inventories, event_kind='REQUESTED')
                return c
            active = (active & ~tensile) | penetrating
        raise CandidateFailure('bounded contact active sets exhausted', c, inventories)

    def _release_event(self, t, loaded):
        """Locate first zero reaction on the still-contacting branch, without commit."""
        lo, hi = self.time, t
        active = self.active.copy()
        endpoint = self._equilibrate(hi, active)
        if np.min(endpoint['normal_force'][loaded]) >= 0:
            raise CandidateFailure('unbracketed release event', endpoint)
        lower = None
        for _ in range(self.ccfg['numerical']['event_bisections']):
            mid = .5*(lo+hi)
            branch = self._equilibrate(mid, active)
            if np.min(branch['normal_force'][loaded]) >= 0:
                lo, lower = mid, branch
            else:
                hi = mid
            if hi-lo <= self.ccfg['numerical']['event_time_tolerance_s']:
                break
        if lower is None or hi-lo > self.ccfg['numerical']['event_time_tolerance_s']:
            raise CandidateFailure('release event location failed', endpoint)
        lower.update(event_kind='RELEASE_LIMIT', event_bracket_s=[lo, hi])
        return lower

    def advance(self, t, on_accept=None):
        if self.time < .1 < t:
            raise ValueError('requested interval crosses prescribed program corner at .1s')
        for _ in range(self.ccfg['numerical']['maximum_events_per_requested_interval']):
            c = self.solve_candidate(t)
            loaded = self.active & (self.normal_force > self.ccfg['numerical']['contact_force_tolerance_N'])
            if np.any(loaded & ~c['active']):
                c = self._release_event(t, loaded)
            self.commit(c)
            if on_accept is not None:
                on_accept()
            if self.time == t:
                return
        raise CandidateFailure('bounded release events exhausted')

    def commit(self, candidate):
        if candidate['base_step'] != self.count:
            raise ValueError('stale candidate')
        # Reconstruct material/forces/KKT, rather than trust candidate ledger fields.
        c = self.trial(candidate['displacement'], candidate['time_s'])
        c['active'] = np.asarray(candidate['active'], dtype=bool).copy()
        c = self.contact_state(c); lim = self.ccfg['numerical']
        if (c['equilibrium_relative'] > lim['force_relative'] or
                np.any(c['normal_force'] < -lim['contact_force_tolerance_N']) or
                np.any(c['gap'] < -lim['gap_tolerance_m']) or
                np.any(np.abs(c['gap'][c['active']]) > lim['gap_tolerance_m'])):
            raise CandidateFailure('contact equilibrium/complementarity rejected', c)
        if np.max(np.linalg.norm(c['displacement'], axis=1))/.01 >= self.qcfg['numerical']['displacement_over_length_limit']:
            raise CandidateFailure('left declared small-strain domain', c)
        du = c['displacement']-self.displacement
        dd = self.cap_displacement(c['time_s'])-self.cap_displacement(self.time)
        wc = float(np.sum(.5*(self.contact_force+c['contact_force'])*du))
        ww = -float(.5*(self.normal_force+c['normal_force']).sum())*dd
        wr = float(np.sum(.5*(self.support_force+c['support_force'])*du))
        ws = float(np.einsum('i,iab,iab->', self.V, .5*(self.stress+c['stress']), self.mesh.gradient(du)))
        if not np.isfinite([wc, ww, wr, ws]).all():
            raise CandidateFailure('nonfinite contact work', c)
        if (abs(wc+wr-ws) > lim['work_tolerance_J'] or abs(ww-wc) > lim['work_tolerance_J'] or
                abs(self.wall_work+ww-self.contact_work-wc) > lim['work_tolerance_J']):
            raise CandidateFailure('contact work requires event split; no residual heat', c)
        event_kind = candidate.get('event_kind', 'DIRECT_COMMIT')
        if event_kind not in ('REQUESTED', 'RELEASE_LIMIT', 'DIRECT_COMMIT'):
            raise ValueError('unknown event metadata')
        event = None
        if event_kind == 'RELEASE_LIMIT':
            bracket = np.asarray(candidate.get('event_bracket_s', []), dtype=float)
            if (bracket.shape != (2,) or not np.isfinite(bracket).all() or
                    bracket[0] != c['time_s'] or bracket[1] < bracket[0] or bracket[1] > .2 or
                    bracket[1]-bracket[0] > lim['event_time_tolerance_s']):
                raise ValueError('invalid release event bracket metadata')
            event = dict(step=self.count+1, physical_time_s=c['time_s'],
                         bracket_s=bracket.tolist(), kind=event_kind)
        # No accepted-state mutation before the preceding checks.
        self.displacement = c['displacement'].copy(); self.stress = c['stress']
        self.pa.x[:], self.pa.y[:], self.pa.z[:] = (self.X+self.displacement).T
        self.pa.S[:] = c['S'].ravel(); self.pa.ep[:] = c['ep']
        self.pa.plastic_work[:] += c['plastic_increment_J']
        self.pa.return_dissipation[:] += c['return_increment_J']; self.pa.strain_rate[:] = c['strain_rate']
        self.rate = du/c['dt_s']; self.time = c['time_s']; self.count += 1
        self.active = c['active']; self.normal_force = c['normal_force']
        self.contact_force = c['contact_force']; self.support_force = c['support_force']
        self.reaction = self.support_force+self.contact_force
        self.wall_work += ww; self.contact_work += wc; self.support_work += wr
        self.internal_work += ws; self.constraint_work = self.contact_work+self.support_work
        self.last_wall_increment = ww; self.last_contact_increment = wc
        self.last_residual = c['equilibrium_relative']
        self.last_iterations = candidate.get('iterations', 0)
        self.last_trace = candidate.get('active_set_trace', candidate.get('iteration_trace', []))
        force = c['internal_force']
        self.last_force_defect = float(np.linalg.norm(force.sum(axis=0))/self.force_scale)
        self.last_torque_defect = float(np.linalg.norm(np.cross(self.X, force).sum(axis=0))/(self.force_scale*.01))
        self.last_event = event_kind
        if event is not None:
            self.events.append(event)

    def record(self):
        row = super().record()
        # Parent's actuator_work means sum of forces at MATERIAL coordinates.
        row.pop('actuator_work_J'); row.pop('cap_reaction_x_N')
        gap = self.cap_displacement(self.time)-self.displacement[self.nodes, 0]
        row.update(wall_displacement_m=self.cap_displacement(self.time),
            wall_force_on_material_x_N=-float(self.normal_force.sum()),
            wall_work_J=self.wall_work, contact_material_work_J=self.contact_work,
            support_work_J=self.support_work, wall_minus_material_work_J=self.wall_work-self.contact_work,
            interval_wall_work_J=self.last_wall_increment, interval_contact_material_work_J=self.last_contact_increment,
            active_contacts=int(self.active.sum()), minimum_gap_m=float(gap.min()),
            maximum_gap_m=float(gap.max()), minimum_normal_force_N=float(self.normal_force.min()),
            maximum_abs_complementarity_J=float(np.max(np.abs(self.normal_force*gap))),
            event_kind=self.last_event)
        return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--law', choices=('elastic','j2'), required=True)
    parser.add_argument('--increments-per-phase', type=int, choices=(8,16), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter(); b = ContactBlock(args.law); setup = time.perf_counter()-start
    records, frames = [], []
    def capture():
        records.append(b.record())
        frames.append(dict(position=b.X+b.displacement, displacement=b.displacement.copy(),
            stress=b.stress.copy(), ep=b.pa.ep.copy(), plastic=b.pa.plastic_work.copy(),
            returned=b.pa.return_dissipation.copy(), contact=b.contact_force.copy(), support=b.support_force.copy(),
            normal=b.normal_force.copy(), active=b.active.copy(),
            gap=b.cap_displacement(b.time)-b.displacement[b.nodes,0], rate=b.rate.copy()))
    capture(); status='UNILATERAL_PATH_COMPLETE'; error=None; tick=time.perf_counter()
    try:
        for t in np.linspace(0., .2, 2*args.increments_per_phase+1)[1:]:
            b.advance(float(t), capture)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        status='BLOCKED'; error=str(exc); c=getattr(exc,'candidate',None)
        if c is not None:
            np.savez_compressed(args.output/'REJECTED_candidate.npz',
                **{k:v for k,v in c.items() if isinstance(v,np.ndarray)})
        detail=dict(classification='REJECTED',restart_allowed=False,error=error,
            attempted_time_s=float(t),last_accepted_step=b.count,candidate_available=c is not None,
            candidate_time_s=None if c is None else c.get('time_s'),
            candidate_dt_s=None if c is None else c.get('dt_s'),
            candidate_base_step=None if c is None else c.get('base_step'),
            iterations=getattr(exc,'iterations',[]))
        (args.output/'REJECTED_manifest.json').write_text(json.dumps(detail,indent=2)+'\n',encoding='utf-8')
    advance=time.perf_counter()-tick
    config=dict(**b.ccfg,selected_law=args.law,increments_per_phase=args.increments_per_phase,
                geometry_elastic_inputs=b.mesh.cfg,plastic_parent=b.qcfg)
    (args.output/'resolved_config.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(args.output/'history.npz',particle_id=b.pa.particle_id.copy(),tracer_type=b.pa.tracer_type.copy(),
        reference_position_m=b.X,reference_volume_m3=b.V,mass_kg=b.pa.m.copy(),support_components=b.fixed,
        contact_particle_id=b.pa.particle_id[b.nodes],contact_feature_id=np.ones(len(b.nodes),dtype=int),
        physical_time_s=np.array([r['physical_time_s'] for r in records]),
        position_m=np.array([f['position'] for f in frames]),displacement_m=np.array([f['displacement'] for f in frames]),
        stress_Pa=np.array([f['stress'] for f in frames]),equivalent_plastic_strain=np.array([f['ep'] for f in frames]),
        plastic_work_J=np.array([f['plastic'] for f in frames]),return_dissipation_diagnostic_J=np.array([f['returned'] for f in frames]),
        contact_force_N=np.array([f['contact'] for f in frames]),support_force_N=np.array([f['support'] for f in frames]),
        normal_force_N=np.array([f['normal'] for f in frames]),active_contacts=np.array([f['active'] for f in frames]),
        gap_m=np.array([f['gap'] for f in frames]),kinematic_interval_rate_m_s=np.array([f['rate'] for f in frames]))
    checkpoint(args.output/'checkpoint.npz',b.pa,dict(classification=b.ccfg['classification'],restart_allowed=False,
        physical_time_s=b.time,step=b.count,production_authorized=False,not_solved=b.ccfg['not_solved']))
    scalars=dict(particle_id=b.pa.particle_id.tolist(),tracer_type=b.pa.tracer_type.tolist(),
        equivalent_plastic_strain=b.pa.ep.tolist(),plastic_work_J=b.pa.plastic_work.tolist())
    for a,c,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        scalars['stress_'+name+'_Pa']=b.stress[:,a,c].tolist()
    write_point_cloud(args.output/'final.vtk',(b.X+b.displacement).tolist(),scalars,
        dict(reference_position_m=b.X.tolist(),displacement_m=b.displacement.tolist()),
        title='Q3_ACCEPTED_UNILATERAL_BENCHMARK_NOT_AFSD_THERMAL_NOT_SOLVED')
    def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
    source=['scripts/probe_process_contact.py','configs/benchmarks/process_contact_q3.json',
        'scripts/probe_process_plastic.py','scripts/probe_process_time.py','scripts/probe_engineering_material.py',
        'scripts/probe_engineering_block.py','scripts/probe_boundary_sbp.py',
        'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']
    report=dict(status=status,error=error,classification=b.ccfg['classification'],production_authorized=False,
        selected_law=args.law,particles=len(b.X),permanent_free_dofs=int(b.free.sum()),
        increments_per_phase=args.increments_per_phase,accepted_steps=b.count,accepted_time_s=b.time,
        requested_end_s=.2,setup_wall_s=setup,advance_wall_s=advance,wall_before_report_s=time.perf_counter()-start,
        selected_array_bytes=sum(p.get_npy_array().nbytes for p in b.pa.properties.values())+b.mesh.D.nbytes+b.Lmap.nbytes+b.K.nbytes,
        memory_note='selected retained arrays only; transient tangent buffers and process RSS excluded',
        condition_M_max=b.mesh.condition,events=b.events,records=records,not_solved=b.ccfg['not_solved'],
        source_sha256={name:sha(ROOT/name) for name in source},
        artifact_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','source_sha256','artifact_sha256')},indent=2))
    print(json.dumps(records[-1],indent=2,allow_nan=False))
    return 0 if status=='UNILATERAL_PATH_COMPLETE' else 2


if __name__=='__main__': raise SystemExit(main())

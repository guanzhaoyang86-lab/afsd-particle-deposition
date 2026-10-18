"""Q5 directional-Y friction companion; NOT isotropic contact or production AFSD."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from probe_process_contact import ContactBlock, reduced_correction
from probe_process_plastic import CandidateFailure, PlasticBlock, ROOT
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG=ROOT/'configs/benchmarks/process_friction_q5.json'


def feasible(c, modes, mu, limits, strict=False):
    if not np.isfinite(mu) or mu<=0: raise ValueError('positive finite friction coefficient required')
    for key in ('normal_force','tangent_force','relative_increment'):
        if np.shape(c[key])!=np.shape(modes) or not np.isfinite(c[key]).all():
            raise CandidateFailure('invalid/nonfinite contact inventory',c)
    if np.any(c['normal_force'] < -limits['force_tolerance_N']):
        raise CandidateFailure('normal detachment outside Q5 scope', c)
    ft=0. if strict else limits['force_tolerance_N']
    st=0. if strict else limits['slip_tolerance_m']
    stick=modes==0
    return (not np.any(np.abs(c['tangent_force'][stick])-mu*c['normal_force'][stick]>ft)
            and not np.any(modes[~stick]*c['relative_increment'][~stick] < -st))


def solve_friction_modes(branch, count, mu, limits):
    """Incremental stick feasibility, followed by a WHOLE coupled branch solve."""
    modes=np.zeros(count,dtype=int); trace=[]
    for _ in range(limits['maximum_active_sets']):
        c=branch(modes); trace.append(modes.tolist())
        if feasible(c,modes,mu,limits):
            c.update(modes=modes.copy(),mode_trace=trace)
            return c
        stick=modes==0
        yields=stick & (np.abs(c['tangent_force'])-mu*c['normal_force']>limits['force_tolerance_N'])
        reverses=(~stick) & (modes*c['relative_increment'] < -limits['slip_tolerance_m'])
        modes[yields]=-np.sign(c['tangent_force'][yields]).astype(int)
        modes[reverses]=0
    raise CandidateFailure('bounded friction modes exhausted',c,trace)


def advance_contact(owner, t, on_accept=None):
    """Bounded first transition on a declared short monotone loading branch.

    owner supplies time,count,nodes,mu,limits,normal_force,tangent_force,
    branch(t,modes),commit(candidate).
    The tiny candidate is a mode probe only, never a committed physical state.
    Program corners must be split by caller. No arbitrary event-detection claim.
    """
    if not np.isfinite(t) or t<=owner.time: raise ValueError('time must advance')
    lim=owner.limits
    for _ in range(lim['maximum_events_per_interval']):
        probe=owner.time+(t-owner.time)*lim['mode_probe_fraction']
        if not owner.time<probe<=t: raise CandidateFailure('unresolved mode probe time')
        start=solve_friction_modes(lambda m:owner.branch(probe,m),len(owner.nodes),owner.mu,lim)
        modes=start['modes'].copy()
        # A finite mode probe may cross a nearby yield point. Never start a
        # sliding branch unless the ACCEPTED force is on that very boundary.
        off_boundary=np.abs(owner.tangent_force+owner.mu*modes*owner.normal_force)>lim['force_tolerance_N']
        modes[(modes!=0) & off_boundary]=0
        c=owner.branch(t,modes)
        if feasible(c,modes,owner.mu,lim):
            c.update(event_kind='REQUESTED',modes=modes,mode_trace=start['mode_trace'])
        else:
            lo,hi=owner.time,t; lower=None
            for _ in range(lim['event_bisections']):
                mid=.5*(lo+hi); b=owner.branch(mid,modes)
                if feasible(b,modes,owner.mu,lim,strict=True): lo,lower=mid,b
                else: hi=mid
                if hi-lo<=lim['event_time_tolerance_s']: break
            if lower is None or hi-lo>lim['event_time_tolerance_s']:
                raise CandidateFailure('stick/slip event not bracketed',c)
            c=lower
            c.update(event_kind='FRICTION_LIMIT',event_bracket_s=[lo,hi],
                     modes=modes,mode_trace=start['mode_trace'])
        owner.commit(c)
        if on_accept is not None: on_accept()
        if owner.time==t: return
    raise CandidateFailure('bounded friction events exhausted')


def contact_work(old_force, force, du, dw, nodes, dt):
    mid=.5*(old_force+force)
    slip=du[nodes,1]-dw[1]
    impulse=mid[nodes,1]*dt
    rate=slip/dt
    raw=-impulse*rate
    return dict(friction_increment_J=raw,relative_increment=slip,
                tangent_impulse_N_s=impulse,relative_kinematic_rate_m_s=rate,
                material_work_J=float(np.sum(mid*du)),wall_work_J=float(np.sum(mid,axis=0)@dw))


class FrictionBlock(ContactBlock):
    def __init__(self):
        super().__init__('elastic',3)
        self.fcfg=json.loads(CONFIG.read_text(encoding='utf-8'))
        expected=dict(mu=.1,normal_displacement=-4e-5,tangent_amplitude=2e-5,phase_duration=.05,temperature=293.15)
        if self.fcfg['production_authorized'] is not False or any(
                self.fcfg['parameters'][k]['value']!=v for k,v in expected.items()):
            raise ValueError('Q5 declared physics/scope mismatch')
        self.mu=self.fcfg['parameters']['mu']['value']; self.limits=self.fcfg['numerical']
        self.fixed=np.broadcast_to((self.X[:,0]==0)[:,None],self.X.shape).copy()
        self.free=~self.fixed.ravel(); self.tangent_dofs=3*self.nodes+1
        self.modes=np.zeros(len(self.nodes),dtype=int)
        self.tangent_force=np.zeros(len(self.nodes)); self.friction=np.zeros(len(self.nodes))
        self.last_friction=np.zeros(len(self.nodes)); self.last_slip=np.zeros(len(self.nodes))
        self.last_impulse=np.zeros(len(self.nodes)); self.last_relative_rate=np.zeros(len(self.nodes))

    def wall(self,t):
        if not np.isfinite(t) or t<0 or t>.2: raise ValueError('outside Q5 schedule')
        dx=-4e-5*t/.05 if t<=.05 else -4e-5
        if t<=.05: y=0.
        elif t<=.1: y=2e-5*(t-.05)/.05
        elif t<=.15: y=2e-5-4e-5*(t-.1)/.05
        else: y=-2e-5+2e-5*(t-.15)/.05
        return np.array([dx,y,0.])

    def cap_displacement(self,t): return self.wall(t)[0]

    def prescribed(self,t):
        self.wall(t)
        return np.zeros_like(self.X)

    def state(self,c,modes):
        modes=np.asarray(modes,dtype=int)
        if modes.shape!=(len(self.nodes),) or not np.isin(modes,[-1,0,1]).all():
            raise ValueError('invalid friction mode inventory')
        force=c['internal_force']; normal=force[self.nodes,0]
        tau=-force[self.nodes,1].copy(); sliding=modes!=0
        tau[sliding]=-self.mu*modes[sliding]*normal[sliding]
        contact=np.zeros_like(self.X); contact[self.nodes,0]=-normal; contact[self.nodes,1]=tau
        support=np.zeros_like(self.X); support[self.fixed]=-(force+contact)[self.fixed]
        gap=self.wall(c['time_s'])[0]-c['displacement'][self.nodes,0]
        slip=c['displacement'][self.nodes,1]-self.displacement[self.nodes,1]-(self.wall(c['time_s'])-self.wall(self.time))[1]
        residual=force+contact+support
        c.update(normal_force=normal.copy(),tangent_force=tau,contact_force=contact,support_force=support,
            relative_increment=slip,gap=gap,modes=modes.copy(),residual=residual,
            equilibrium_relative=float(np.linalg.norm(residual)/self.force_scale))
        return c

    def branch(self,t,modes,maximum_iterations=None):
        modes=np.asarray(modes,dtype=int); stick=modes==0
        unknown=self.free.copy(); unknown[self.contact_dofs]=False; unknown[self.tangent_dofs[stick]]=False
        u=self.displacement.copy(); u[self.nodes,0]=self.wall(t)[0]
        u[self.nodes[stick],1]=self.displacement[self.nodes[stick],1]+(self.wall(t)-self.wall(self.time))[1]
        maximum=self.qcfg['numerical']['maximum_iterations'] if maximum_iterations is None else maximum_iterations
        trace=[]
        for iteration in range(maximum+1):
            c=self.state(self.trial(u,t),modes)
            residual=c['residual'].ravel(); norm=float(np.linalg.norm(residual[unknown])/self.force_scale)
            trace.append(dict(iteration=iteration,free_force_relative=norm))
            if norm<=self.ccfg['numerical']['force_relative']:
                c.update(iteration_trace=trace,iterations=iteration)
                return c
            if iteration==maximum: break
            K=self.tangent(c)
            for row in np.flatnonzero(~stick):
                K[self.tangent_dofs[row]]-=self.mu*modes[row]*K[self.contact_dofs[row]]
            try: step,cond=reduced_correction(K,residual,unknown,self.ccfg['numerical']['tangent_condition_limit'])
            except (ValueError,np.linalg.LinAlgError) as exc: raise CandidateFailure(str(exc),c,trace) from exc
            trace[-1]['tangent_condition']=cond
            for half in range(self.qcfg['numerical']['line_search_halvings']+1):
                test=u+2.**(-half)*step.reshape(u.shape)
                probe=self.state(self.trial(test,t),modes)
                if np.linalg.norm(probe['residual'].ravel()[unknown])/self.force_scale < norm:
                    u=test; break
            else: raise CandidateFailure('friction branch line search failed',c,trace)
        raise CandidateFailure('friction branch equilibrium did not converge',c,trace)

    def advance(self,t,on_accept=None):
        if any(self.time<k-1e-14 and t>k+1e-14 for k in (.05,.1,.15)):
            raise ValueError('requested interval crosses Q5 program corner')
        return advance_contact(self,t,on_accept)

    def commit(self,candidate):
        if candidate['base_step']!=self.count: raise ValueError('stale candidate')
        c=self.state(self.trial(candidate['displacement'],candidate['time_s']),candidate['modes'])
        if not feasible(c,c['modes'],self.mu,self.limits): raise CandidateFailure('friction mode infeasible',c)
        if (c['equilibrium_relative']>self.ccfg['numerical']['force_relative'] or
                np.any(np.abs(c['gap'])>self.limits['slip_tolerance_m']) or
                np.any(np.abs(c['relative_increment'][c['modes']==0])>self.limits['slip_tolerance_m'])):
            raise CandidateFailure('friction equilibrium/constraint rejected',c)
        if np.any(c['ep']!=0) or np.any(c['plastic_increment_J']!=0): raise CandidateFailure('Q5 elastic-only yielded',c)
        sliding=c['modes']!=0
        old_limit_error=np.abs(self.tangent_force+self.mu*c['modes']*self.normal_force)
        if np.any(sliding & (np.abs(c['relative_increment'])>self.limits['slip_tolerance_m']) &
                  (old_limit_error>self.limits['force_tolerance_N'])):
            raise CandidateFailure('unsplit stick/slip transition',c)
        if np.max(np.linalg.norm(c['displacement'],axis=1))/.01>=.02: raise CandidateFailure('left small strain domain',c)
        du=c['displacement']-self.displacement; dw=self.wall(c['time_s'])-self.wall(self.time)
        work=contact_work(self.contact_force,c['contact_force'],du,dw,self.nodes,c['dt_s'])
        q=work['friction_increment_J']; wc=work['material_work_J']; ww=work['wall_work_J']
        wr=float(np.sum(.5*(self.support_force+c['support_force'])*du))
        ws=float(np.einsum('i,iab,iab->',self.V,.5*(self.stress+c['stress']),self.mesh.gradient(du)))
        if not np.isfinite(np.r_[q,wc,ww,wr,ws]).all(): raise CandidateFailure('nonfinite friction work',c)
        if np.any(q < -self.limits['negative_dissipation_tolerance_J']):
            raise CandidateFailure('negative raw friction dissipation',c)
        if max(abs(ww-wc-q.sum()),abs(wc+wr-ws),
               abs(self.wall_work+ww-self.contact_work-wc-(self.friction.sum()+q.sum())),
               abs(self.contact_work+wc+self.support_work+wr-self.internal_work-ws))>self.limits['work_tolerance_J']:
            raise CandidateFailure('friction work identity failed',c)
        kind=candidate.get('event_kind','DIRECT_COMMIT'); event=None
        if kind not in ('REQUESTED','FRICTION_LIMIT','DIRECT_COMMIT'): raise ValueError('unknown event kind')
        if kind=='FRICTION_LIMIT':
            bracket=np.asarray(candidate.get('event_bracket_s',[]))
            if (bracket.shape!=(2,) or not np.isfinite(bracket).all() or bracket[0]!=c['time_s'] or
                    not 0<=bracket[1]-bracket[0]<=self.limits['event_time_tolerance_s']):
                raise ValueError('invalid friction event bracket')
            event=dict(step=self.count+1,physical_time_s=c['time_s'],kind=kind,bracket_s=bracket.tolist(),
                       ending_interval_modes=c['modes'].tolist(),cone_margin_N=(self.mu*c['normal_force']-np.abs(c['tangent_force'])).tolist())
        # All checks precede mutation; material commits ONCE for all contacts.
        self.displacement=c['displacement'].copy(); self.stress=c['stress']; self.pa.S[:]=c['S'].ravel()
        self.pa.ep[:]=c['ep']; self.pa.strain_rate[:]=c['strain_rate']
        self.pa.x[:],self.pa.y[:],self.pa.z[:]=(self.X+self.displacement).T
        self.rate=du/c['dt_s']; self.time=c['time_s']; self.count+=1
        self.modes=c['modes']; self.normal_force=c['normal_force']; self.tangent_force=c['tangent_force']
        self.contact_force=c['contact_force']; self.support_force=c['support_force']; self.reaction=self.contact_force+self.support_force
        self.friction+=q; self.last_friction=q; self.last_slip=work['relative_increment']
        self.last_impulse=work['tangent_impulse_N_s']; self.last_relative_rate=work['relative_kinematic_rate_m_s']
        self.wall_work+=ww; self.contact_work+=wc; self.support_work+=wr; self.internal_work+=ws
        self.constraint_work=self.contact_work+self.support_work
        self.last_wall_increment=ww; self.last_contact_increment=wc; self.last_event=kind
        self.last_residual=c['equilibrium_relative']; self.last_iterations=candidate.get('iterations',0)
        self.last_trace=candidate.get('iteration_trace',[])
        force=c['internal_force']
        self.last_force_defect=float(np.linalg.norm(force.sum(axis=0))/self.force_scale)
        self.last_torque_defect=float(np.linalg.norm(np.cross(self.X,force).sum(axis=0))/(self.force_scale*.01))
        if event is not None: self.events.append(event)

    def record(self):
        row=PlasticBlock.record(self)
        for name in ('actuator_work_J','cap_reaction_x_N','cap_displacement_m'): row.pop(name)
        row.update(wall_displacement_m=self.wall(self.time).tolist(),wall_force_on_material_N=self.contact_force.sum(axis=0).tolist(),
            wall_work_J=self.wall_work,contact_material_work_J=self.contact_work,support_work_J=self.support_work,
            friction_dissipation_J=float(self.friction.sum()),interval_friction_J=float(self.last_friction.sum()),
            minimum_raw_contact_dissipation_J=float(self.last_friction.min()),
            interval_wall_work_J=self.last_wall_increment,interval_contact_material_work_J=self.last_contact_increment,
            friction_work_identity_J=self.wall_work-self.contact_work-float(self.friction.sum()),
            mode_counts={str(m):int((self.modes==m).sum()) for m in (-1,0,1)},
            minimum_normal_force_N=float(self.normal_force.min()),event_kind=self.last_event)
        return row


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--increments-per-phase',type=int,choices=(4,8),required=True)
    p.add_argument('--output',type=Path,required=True); a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter(); b=FrictionBlock(); setup=time.perf_counter()-start; frames=[]; records=[]
    def capture():
        records.append(b.record()); frames.append(dict(position=b.X+b.displacement,u=b.displacement.copy(),stress=b.stress.copy(),
            normal=b.normal_force.copy(),tau=b.tangent_force.copy(),modes=b.modes.copy(),force=b.contact_force.copy(),
            support=b.support_force.copy(),q=b.friction.copy(),dq=b.last_friction.copy(),slip=b.last_slip.copy(),
            impulse=b.last_impulse.copy(),relrate=b.last_relative_rate.copy(),rate=b.rate.copy()))
    capture(); tick=time.perf_counter(); status='DIRECTIONAL_FRICTION_PATH_COMPLETE'; error=None
    try:
        for t in np.linspace(0.,.2,4*a.increments_per_phase+1)[1:]: b.advance(float(t),capture)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        status='BLOCKED'; error=str(exc); c=getattr(exc,'candidate',None)
        if c is not None: np.savez_compressed(a.output/'REJECTED_candidate.npz',**{k:v for k,v in c.items() if isinstance(v,np.ndarray)})
        detail=dict(classification='REJECTED',restart_allowed=False,error=error,attempted_time_s=float(t),last_accepted_step=b.count,
            candidate_available=c is not None,candidate_time_s=None if c is None else c.get('time_s'),candidate_base_step=None if c is None else c.get('base_step'))
        (a.output/'REJECTED_manifest.json').write_text(json.dumps(detail,indent=2)+'\n',encoding='utf-8')
    advance=time.perf_counter()-tick
    cfg=dict(**b.fcfg,increments_per_phase=a.increments_per_phase,geometry_inputs=b.mesh.cfg,plastic_parent=b.qcfg,contact_parent=b.ccfg)
    (a.output/'resolved_config.json').write_text(json.dumps(cfg,indent=2)+'\n',encoding='utf-8')
    data=dict(particle_id=b.pa.particle_id.copy(),tracer_type=b.pa.tracer_type.copy(),reference_position_m=b.X,reference_volume_m3=b.V,
        mass_kg=b.pa.m.copy(),support_components=b.fixed,contact_particle_id=b.pa.particle_id[b.nodes],contact_feature_id=np.ones(len(b.nodes),dtype=int),
        physical_time_s=np.array([r['physical_time_s'] for r in records]),equivalent_plastic_strain=b.pa.ep.copy())
    mapping=dict(position_m='position',displacement_m='u',stress_Pa='stress',normal_force_N='normal',tangent_force_N='tau',
        interval_modes='modes',contact_force_N='force',support_force_N='support',friction_dissipation_J='q',interval_friction_J='dq',
        interval_relative_slip_m='slip',tangent_impulse_N_s='impulse',relative_kinematic_rate_m_s='relrate',kinematic_interval_rate_m_s='rate')
    data.update({key:np.array([f[field] for f in frames]) for key,field in mapping.items()})
    np.savez_compressed(a.output/'history.npz',**data)
    checkpoint(a.output/'checkpoint.npz',b.pa,dict(restart_allowed=False,classification=b.fcfg['classification'],physical_time_s=b.time,step=b.count,production_authorized=False))
    scalars=dict(particle_id=b.pa.particle_id.tolist(),tracer_type=b.pa.tracer_type.tolist(),equivalent_plastic_strain=b.pa.ep.tolist())
    for i,j,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')): scalars['stress_'+name+'_Pa']=b.stress[:,i,j].tolist()
    write_point_cloud(a.output/'final.vtk',(b.X+b.displacement).tolist(),scalars,dict(reference_position_m=b.X.tolist(),displacement_m=b.displacement.tolist()),title='Q5_DIRECTIONAL_ACCEPTED_NOT_AFSD_THERMAL_NOT_SOLVED')
    source=['scripts/probe_process_friction.py','configs/benchmarks/process_friction_q5.json','scripts/probe_process_contact.py','scripts/probe_process_plastic.py',
        'scripts/probe_process_time.py','scripts/probe_engineering_material.py','scripts/probe_engineering_block.py','scripts/probe_boundary_sbp.py',
        'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    report=dict(status=status,error=error,production_authorized=False,classification=b.fcfg['classification'],particles=len(b.X),contact_rows=len(b.nodes),
        permanent_free_dofs=int(b.free.sum()),increments_per_phase=a.increments_per_phase,accepted_steps=b.count,accepted_time_s=b.time,requested_end_s=.2,
        setup_wall_s=setup,advance_wall_s=advance,wall_before_report_s=time.perf_counter()-start,events=b.events,records=records,not_solved=b.fcfg['not_solved'],
        selected_array_bytes=sum(prop.get_npy_array().nbytes for prop in b.pa.properties.values())+b.K.nbytes+b.Lmap.nbytes+b.mesh.D.nbytes,
        memory_note='selected retained arrays, NOT peak RSS',condition_M_max=b.mesh.condition,
        source_sha256={s:sha(ROOT/s) for s in source},artifact_sha256={f.name:sha(f) for f in a.output.iterdir() if f.is_file()})
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','source_sha256','artifact_sha256')},indent=2)); print(json.dumps(records[-1],indent=2))
    return 0 if status=='DIRECTIONAL_FRICTION_PATH_COMPLETE' else 2


if __name__=='__main__':raise SystemExit(main())

"""Q6 full tangent disk, fixed normal, bounded quasi-static benchmark ONLY."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from probe_process_friction import FrictionBlock
from probe_process_contact import reduced_correction
from probe_process_plastic import CandidateFailure, ROOT
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG=ROOT/'configs/benchmarks/process_friction_cone_q6.json'


def release_inventory(released, count):
    mask=np.asarray(released)
    if mask.shape!=(count,) or mask.dtype!=np.dtype(bool):
        raise ValueError('invalid released cone inventory: boolean mask required')
    return mask


def cone_response(tau, normal, slip, mu, kappa, basis=None):
    """Projection defines a ROOT EQUATION; actual forces are never clipped.

    tau/slip use global YZ coordinates. Returned r,Pq,Pc also use global YZ.
    basis columns are local tangent axes in global YZ (representation only).
    """
    tau,normal,slip=map(lambda a:np.asarray(a,dtype=float),(tau,normal,slip))
    if (tau.ndim!=2 or tau.shape[1]!=2 or slip.shape!=tau.shape or normal.shape!=(len(tau),)
            or not all(np.isfinite(a).all() for a in (tau,normal,slip))):
        raise ValueError('invalid/nonfinite cone inventory')
    if not np.isfinite([mu,kappa]).all() or mu<=0 or kappa<=0: raise ValueError('positive cone numerical inputs required')
    if np.any(normal<0): raise ValueError('negative normal outside fixed-contact cone scope')
    B=np.eye(2) if basis is None else np.asarray(basis,dtype=float)
    if B.shape!=(2,2) or not np.isfinite(B).all() or not np.allclose(B.T@B,np.eye(2),atol=1e-13,rtol=0):
        raise ValueError('orthonormal tangent basis required')
    q=(tau-kappa*slip)@B; norm=np.linalg.norm(q,axis=1); C=mu*normal
    outside=norm>C; projected=q.copy(); D=np.broadcast_to(np.eye(2),(len(tau),2,2)).copy(); pc=np.zeros_like(tau)
    # At q=C=0 the inside derivative is finite; no division by zero/pinv.
    for i in np.flatnonzero(outside):
        n=q[i]/norm[i]; projected[i]=C[i]*n
        D[i]=(C[i]/norm[i])*(np.eye(2)-np.outer(n,n)); pc[i]=n
    return dict(residual=projected@B.T-tau,Pq=B@D@B.T,Pc=pc@B.T,outside=outside)


def cone_rows(K, residual, normal_dofs, tangent_dofs, response, mu, kappa, released):
    """Replace contact rows by -Jacobian of the coupled natural-map residual."""
    matrix=K.copy(); force=residual.copy()
    for i in np.flatnonzero(released):
        rows=tangent_dofs[i]; D=response['Pq'][i]; pc=response['Pc'][i]
        block=(np.eye(2)-D)@K[rows]+mu*np.outer(pc,K[normal_dofs[i]])
        block[:,rows]+=kappa*D
        matrix[rows]=block; force[rows]=response['residual'][i]
    return matrix,force


def feasible_cone(c, released, mu, limits, strict=False):
    normal=np.asarray(c['normal_force']); tau=np.asarray(c['tangent_force']); slip=np.asarray(c['relative_increment'])
    n=len(released); released=release_inventory(released,n)
    if (normal.shape!=(n,) or tau.shape!=(n,2) or slip.shape!=(n,2) or
            not all(np.isfinite(a).all() for a in (normal,tau,slip))):
        raise CandidateFailure('invalid/nonfinite cone inventory',c)
    if np.any(normal<0): raise CandidateFailure('normal detachment outside Q6 scope',c)
    ft=limits['force_tolerance_N']; threshold=0. if strict else ft
    margin=mu*normal-np.linalg.norm(tau,axis=1)
    if np.any(margin[~released]<-threshold): return False
    if np.any(margin[released]<-ft): return False
    length=np.linalg.norm(slip,axis=1); moving=released & (length>limits['slip_tolerance_m'])
    if moving.any():
        # Do not divide displacement subtraction noise by a nearly zero slip.
        # Same force/slip tolerances, propagated to this N*m vector relation.
        # The independent cone bound and natural-map root are still required.
        capacity=mu*normal[moving]; distance=length[moving]
        defect=distance[:,None]*tau[moving]+capacity[:,None]*slip[moving]
        bound=ft*distance+capacity*limits['slip_tolerance_m']
        if np.any(np.linalg.norm(defect,axis=1)>bound): return False
    return True


def solve_cone_modes(branch,count,mu,limits):
    released=np.zeros(count,dtype=bool); trace=[]
    for _ in range(limits['maximum_active_sets']):
        c=branch(released); trace.append(released.tolist())
        if feasible_cone(c,released,mu,limits):
            c.update(released=released.copy(),mode_trace=trace)
            return c
        over=np.linalg.norm(c['tangent_force'],axis=1)-mu*c['normal_force']>limits['force_tolerance_N']
        added=(~released)&over
        if not added.any(): raise CandidateFailure('released cone root/flow infeasible',c,trace)
        released|=added
    raise CandidateFailure('bounded cone modes exhausted',c,trace)


def advance_cone(owner,t,on_accept=None):
    if not np.isfinite(t) or t<=owner.time: raise ValueError('time must advance')
    lim=owner.limits
    for _ in range(lim['maximum_events_per_interval']):
        probe=owner.time+(t-owner.time)*lim['mode_probe_fraction']
        if not owner.time<probe<=t: raise CandidateFailure('unresolved cone probe time')
        start=solve_cone_modes(lambda m:owner.branch(probe,m),len(owner.nodes),owner.mu,lim)
        released=start['released'].copy()
        oldnorm=np.linalg.norm(owner.tangent_force,axis=1)
        off=np.abs(oldnorm-owner.mu*owner.normal_force)>lim['force_tolerance_N']
        opposite=np.sum(owner.tangent_force*start['tangent_force'],axis=1)<=0.
        released[off|opposite]=False
        c=owner.branch(t,released)
        if feasible_cone(c,released,owner.mu,lim): c.update(event_kind='REQUESTED')
        else:
            lo,hi=owner.time,t; lower=None
            for _ in range(lim['event_bisections']):
                mid=.5*(lo+hi); b=owner.branch(mid,released)
                if feasible_cone(b,released,owner.mu,lim,strict=True): lo,lower=mid,b
                else: hi=mid
                if hi-lo<=lim['event_time_tolerance_s']: break
            if lower is None or hi-lo>lim['event_time_tolerance_s']:
                raise CandidateFailure('cone event not bracketed',c)
            c=lower; c.update(event_kind='CONE_LIMIT',event_bracket_s=[lo,hi])
        c.update(released=released,mode_trace=start['mode_trace'])
        owner.commit(c)
        if on_accept is not None: on_accept()
        if owner.time==t:return
    raise CandidateFailure('bounded cone events exhausted')


def vector_contact_work(old_force,force,du,dw,nodes,dt):
    if not np.isfinite(dt) or dt<=0: raise ValueError('positive interval required')
    mid=.5*(old_force+force); slip=du[nodes,1:]-dw[1:]
    impulse=mid[nodes,1:]*dt; rate=slip/dt; raw=-np.sum(impulse*rate,axis=1)
    return dict(friction_increment_J=raw,relative_increment=slip,tangent_impulse_N_s=impulse,
        relative_kinematic_rate_m_s=rate,material_work_J=float(np.sum(mid*du)),wall_work_J=float(mid.sum(axis=0)@dw))


class ConeBlock(FrictionBlock):
    def __init__(self,basis_angle_deg=0.):
        super().__init__()
        self.conecfg=json.loads(CONFIG.read_text(encoding='utf-8'))
        p=self.conecfg['parameters']
        expected=dict(mu=.1,normal_displacement=-4e-5,tangent_amplitude=2e-5,phase_duration=.05,kappa=7e8)
        if (self.conecfg['production_authorized'] is not False or
                any(p[k]['value']!=v for k,v in expected.items()) or p['tangent_direction_global_yz']['value']!=[.6,.8]):
            raise ValueError('Q6 declared parameter/scope mismatch')
        if not np.isfinite(basis_angle_deg): raise ValueError('invalid basis angle')
        self.basis_angle_deg=float(basis_angle_deg); a=np.deg2rad(basis_angle_deg)
        self.basis=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        self.kappa=p['kappa']['value']; self.tangent_dofs=3*self.nodes[:,None]+np.array([1,2])
        self.tangent_force=np.zeros((len(self.nodes),2)); self.last_slip=np.zeros_like(self.tangent_force)
        self.last_impulse=np.zeros_like(self.tangent_force); self.last_relative_rate=np.zeros_like(self.tangent_force)
        self.last_law_residual=np.zeros_like(self.tangent_force)

    def wall(self,t):
        parent=FrictionBlock.wall(self,t)
        return np.r_[parent[0],parent[1]*np.array([.6,.8])]

    def state(self,c,released):
        released=release_inventory(released,len(self.nodes))
        F=c['internal_force']; normal=F[self.nodes,0].copy(); tau=-F[self.nodes,1:].copy()
        slip=c['displacement'][self.nodes,1:]-self.displacement[self.nodes,1:]-(self.wall(c['time_s'])-self.wall(self.time))[1:]
        response=cone_response(tau,normal,slip,self.mu,self.kappa,self.basis)
        contact=np.zeros_like(self.X); contact[self.nodes,0]=-normal; contact[self.nodes,1:]=tau
        support=np.zeros_like(self.X); support[self.fixed]=-(F+contact)[self.fixed]
        residual=F+contact+support
        c.update(normal_force=normal,tangent_force=tau,relative_increment=slip,
            contact_force=contact,support_force=support,residual=residual,response=response,
            gap=self.wall(c['time_s'])[0]-c['displacement'][self.nodes,0],released=released.copy(),
            modes=(np.linalg.norm(slip,axis=1)>self.limits['slip_tolerance_m']).astype(int),
            equilibrium_relative=float(np.linalg.norm(residual)/self.force_scale),
            law_residual_N=response['residual'].copy())
        return c

    def branch(self,t,released,maximum_iterations=None):
        released=release_inventory(released,len(self.nodes)); stick=~released
        unknown=self.free.copy(); unknown[self.contact_dofs]=False; unknown[self.tangent_dofs[stick].ravel()]=False
        u=self.displacement.copy(); u[self.nodes,0]=self.wall(t)[0]
        u[self.nodes[stick],1:]=self.displacement[self.nodes[stick],1:]+(self.wall(t)-self.wall(self.time))[1:]
        maximum=self.qcfg['numerical']['maximum_iterations'] if maximum_iterations is None else maximum_iterations
        if maximum<0:raise ValueError('negative iteration budget')
        trace=[]
        for iteration in range(maximum+1):
            c=self.state(self.trial(u,t),released)
            residual=c['residual'].ravel().copy()
            residual[self.tangent_dofs[released].ravel()]=c['law_residual_N'][released].ravel()
            norm=float(np.linalg.norm(residual[unknown])/self.force_scale)
            trace.append(dict(iteration=iteration,root_relative=norm))
            if norm<=self.ccfg['numerical']['force_relative']:
                c.update(iterations=iteration,iteration_trace=trace); return c
            if iteration==maximum: break
            K,r=cone_rows(self.tangent(c),c['residual'].ravel(),self.contact_dofs,self.tangent_dofs,c['response'],self.mu,self.kappa,released)
            try: step,cond=reduced_correction(K,r,unknown,self.ccfg['numerical']['tangent_condition_limit'])
            except (ValueError,np.linalg.LinAlgError) as exc: raise CandidateFailure(str(exc),c,trace) from exc
            trace[-1]['tangent_condition']=cond
            for half in range(self.qcfg['numerical']['line_search_halvings']+1):
                trial_u=u+2.**(-half)*step.reshape(u.shape)
                try: probe=self.state(self.trial(trial_u,t),released)
                except ValueError as exc:
                    if 'negative normal' in str(exc): continue
                    raise
                rr=probe['residual'].ravel().copy(); rr[self.tangent_dofs[released].ravel()]=probe['law_residual_N'][released].ravel()
                if np.linalg.norm(rr[unknown])/self.force_scale<norm:
                    u=trial_u; trace[-1]['line_search_factor']=2.**(-half); break
            else: raise CandidateFailure('cone branch line search failed',c,trace)
        raise CandidateFailure('cone branch equilibrium/root did not converge',c,trace)

    def advance(self,t,on_accept=None):
        if any(self.time<k-1e-14 and t>k+1e-14 for k in (.05,.1,.15)):
            raise ValueError('requested interval crosses Q6 program corner')
        return advance_cone(self,t,on_accept)

    def commit(self,candidate):
        if candidate['base_step']!=self.count:raise ValueError('stale candidate')
        c=self.state(self.trial(candidate['displacement'],candidate['time_s']),candidate['released'])
        if (not feasible_cone(c,np.ones(len(self.nodes),dtype=bool),self.mu,self.limits) or
                np.max(np.linalg.norm(c['law_residual_N'],axis=1))>self.limits['force_tolerance_N']):
            raise CandidateFailure('cone contact law/root rejected',c)
        if (c['equilibrium_relative']>self.ccfg['numerical']['force_relative'] or
                np.any(np.abs(c['gap'])>self.limits['slip_tolerance_m'])):
            raise CandidateFailure('cone equilibrium/normal constraint rejected',c)
        if np.any(c['ep']!=0) or np.any(c['plastic_increment_J']!=0):raise CandidateFailure('Q6 elastic-only yielded',c)
        moving=c['modes']!=0
        off=np.abs(np.linalg.norm(self.tangent_force,axis=1)-self.mu*self.normal_force)>self.limits['force_tolerance_N']
        opposite=np.sum(self.tangent_force*c['tangent_force'],axis=1)<0
        if np.any(moving & (off|opposite)):raise CandidateFailure('unsplit cone stick/slip transition',c)
        if np.max(np.linalg.norm(c['displacement'],axis=1))/.01>=.02:raise CandidateFailure('left small strain domain',c)
        du=c['displacement']-self.displacement; dw=self.wall(c['time_s'])-self.wall(self.time)
        work=vector_contact_work(self.contact_force,c['contact_force'],du,dw,self.nodes,c['dt_s'])
        q=work['friction_increment_J']; wc=work['material_work_J']; ww=work['wall_work_J']
        wr=float(np.sum(.5*(self.support_force+c['support_force'])*du))
        ws=float(np.einsum('i,iab,iab->',self.V,.5*(self.stress+c['stress']),self.mesh.gradient(du)))
        if not np.isfinite(np.r_[q,wc,ww,wr,ws]).all():raise CandidateFailure('nonfinite cone work',c)
        if np.any(q < -self.limits['negative_dissipation_tolerance_J']):raise CandidateFailure('negative raw cone friction dissipation',c)
        if max(abs(ww-wc-q.sum()),abs(wc+wr-ws),
               abs(self.wall_work+ww-self.contact_work-wc-(self.friction.sum()+q.sum())),
               abs(self.contact_work+wc+self.support_work+wr-self.internal_work-ws))>self.limits['work_tolerance_J']:
            raise CandidateFailure('cone contact work identity failed',c)
        kind=candidate.get('event_kind','DIRECT_COMMIT'); event=None
        if kind not in ('REQUESTED','CONE_LIMIT','DIRECT_COMMIT'):raise ValueError('invalid cone event kind')
        if kind=='CONE_LIMIT':
            bracket=np.asarray(candidate.get('event_bracket_s',[]),dtype=float)
            if (bracket.shape!=(2,) or not np.isfinite(bracket).all() or bracket[0]!=c['time_s'] or
                    not 0<=bracket[1]-bracket[0]<=self.limits['event_time_tolerance_s']):
                raise ValueError('invalid cone event bracket')
            event=dict(step=self.count+1,physical_time_s=c['time_s'],kind=kind,bracket_s=bracket.tolist(),
                ending_interval_modes=c['modes'].tolist(),released_rows=c['released'].tolist(),
                cone_margin_N=(self.mu*c['normal_force']-np.linalg.norm(c['tangent_force'],axis=1)).tolist())
        # No state changes above. Contact and every material point commit once.
        self.displacement=c['displacement'].copy(); self.stress=c['stress']; self.pa.S[:]=c['S'].ravel(); self.pa.ep[:]=c['ep']
        self.pa.strain_rate[:]=c['strain_rate']; self.pa.x[:],self.pa.y[:],self.pa.z[:]=(self.X+self.displacement).T
        self.rate=du/c['dt_s']; self.time=c['time_s']; self.count+=1; self.modes=c['modes']
        self.normal_force=c['normal_force']; self.tangent_force=c['tangent_force']; self.last_law_residual=c['law_residual_N']
        self.contact_force=c['contact_force']; self.support_force=c['support_force']; self.reaction=self.contact_force+self.support_force
        self.friction+=q; self.last_friction=q; self.last_slip=work['relative_increment']
        self.last_impulse=work['tangent_impulse_N_s']; self.last_relative_rate=work['relative_kinematic_rate_m_s']
        self.wall_work+=ww; self.contact_work+=wc; self.support_work+=wr; self.internal_work+=ws
        self.constraint_work=self.contact_work+self.support_work; self.last_wall_increment=ww; self.last_contact_increment=wc
        self.last_event=kind; self.last_residual=c['equilibrium_relative']; self.last_iterations=candidate.get('iterations',0)
        self.last_trace=candidate.get('iteration_trace',[]); F=c['internal_force']
        self.last_force_defect=float(np.linalg.norm(F.sum(axis=0))/self.force_scale)
        self.last_torque_defect=float(np.linalg.norm(np.cross(self.X,F).sum(axis=0))/(self.force_scale*.01))
        if event is not None:self.events.append(event)

    def record(self):
        row=super().record()
        row.update(mode_counts={str(m):int((self.modes==m).sum()) for m in (0,1)},
            maximum_contact_law_residual_N=float(np.linalg.norm(self.last_law_residual,axis=1).max()),
            minimum_cone_margin_N=float((self.mu*self.normal_force-np.linalg.norm(self.tangent_force,axis=1)).min()))
        return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--increments-per-phase',type=int,choices=(4,8),required=True)
    p.add_argument('--basis-angle-deg',type=float,choices=(0.,37.),default=0.)
    p.add_argument('--output',type=Path,required=True); a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter(); b=ConeBlock(a.basis_angle_deg); setup=time.perf_counter()-start
    frames=[]; records=[]
    def capture():
        records.append(b.record()); frames.append(dict(position=b.X+b.displacement,u=b.displacement.copy(),stress=b.stress.copy(),
            normal=b.normal_force.copy(),tau=b.tangent_force.copy(),modes=b.modes.copy(),force=b.contact_force.copy(),support=b.support_force.copy(),
            q=b.friction.copy(),dq=b.last_friction.copy(),slip=b.last_slip.copy(),impulse=b.last_impulse.copy(),relrate=b.last_relative_rate.copy(),
            rate=b.rate.copy(),lawres=b.last_law_residual.copy()))
    capture(); tick=time.perf_counter(); status='TWO_TANGENT_CONE_PATH_COMPLETE'; error=None
    try:
        for t in np.linspace(0.,.2,4*a.increments_per_phase+1)[1:]:b.advance(float(t),capture)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        status='BLOCKED'; error=str(exc); c=getattr(exc,'candidate',None)
        if c is not None:np.savez_compressed(a.output/'REJECTED_candidate.npz',**{k:v for k,v in c.items() if isinstance(v,np.ndarray)})
        detail=dict(classification='REJECTED',restart_allowed=False,error=error,attempted_time_s=float(t),last_accepted_step=b.count,
            candidate_available=c is not None,candidate_time_s=None if c is None else c.get('time_s'),candidate_base_step=None if c is None else c.get('base_step'))
        (a.output/'REJECTED_manifest.json').write_text(json.dumps(detail,indent=2)+'\n',encoding='utf-8')
    advance=time.perf_counter()-tick
    config=dict(**b.conecfg,increments_per_phase=a.increments_per_phase,basis_angle_deg=a.basis_angle_deg,tangent_basis_global_yz=b.basis.tolist(),
        friction_parent=b.fcfg,geometry_inputs=b.mesh.cfg,plastic_parent=b.qcfg,contact_parent=b.ccfg)
    (a.output/'resolved_config.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    data=dict(particle_id=b.pa.particle_id.copy(),tracer_type=b.pa.tracer_type.copy(),reference_position_m=b.X,reference_volume_m3=b.V,
        mass_kg=b.pa.m.copy(),support_components=b.fixed,contact_particle_id=b.pa.particle_id[b.nodes],contact_feature_id=np.ones(len(b.nodes),dtype=int),
        physical_time_s=np.array([r['physical_time_s'] for r in records]),equivalent_plastic_strain=b.pa.ep.copy())
    mapping=dict(position_m='position',displacement_m='u',stress_Pa='stress',normal_force_N='normal',tangent_force_global_yz_N='tau',
        interval_modes='modes',contact_force_N='force',support_force_N='support',friction_dissipation_J='q',interval_friction_J='dq',
        interval_relative_slip_global_yz_m='slip',tangent_impulse_global_yz_N_s='impulse',relative_kinematic_rate_global_yz_m_s='relrate',
        kinematic_interval_rate_m_s='rate',contact_law_residual_global_yz_N='lawres')
    data.update({key:np.array([f[field] for f in frames]) for key,field in mapping.items()}); np.savez_compressed(a.output/'history.npz',**data)
    checkpoint(a.output/'checkpoint.npz',b.pa,dict(restart_allowed=False,classification=b.conecfg['classification'],physical_time_s=b.time,step=b.count,production_authorized=False))
    scalar=dict(particle_id=b.pa.particle_id.tolist(),tracer_type=b.pa.tracer_type.tolist(),equivalent_plastic_strain=b.pa.ep.tolist())
    for i,j,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):scalar['stress_'+name+'_Pa']=b.stress[:,i,j].tolist()
    write_point_cloud(a.output/'final.vtk',(b.X+b.displacement).tolist(),scalar,dict(reference_position_m=b.X.tolist(),displacement_m=b.displacement.tolist()),title='Q6_ACCEPTED_VECTOR_CONE_NOT_AFSD_THERMAL_NOT_SOLVED')
    source=['scripts/probe_process_friction_cone.py','configs/benchmarks/process_friction_cone_q6.json',
        'scripts/probe_process_friction.py','configs/benchmarks/process_friction_q5.json','scripts/probe_process_contact.py','scripts/probe_process_plastic.py',
        'scripts/probe_process_time.py','scripts/probe_engineering_material.py','scripts/probe_engineering_block.py','scripts/probe_boundary_sbp.py',
        'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    report=dict(status=status,error=error,classification=b.conecfg['classification'],production_authorized=False,particles=len(b.X),contact_rows=len(b.nodes),
        permanent_free_dofs=int(b.free.sum()),increments_per_phase=a.increments_per_phase,basis_angle_deg=a.basis_angle_deg,
        accepted_steps=b.count,accepted_time_s=b.time,requested_end_s=.2,setup_wall_s=setup,advance_wall_s=advance,
        wall_before_report_s=time.perf_counter()-start,events=b.events,records=records,not_solved=b.conecfg['not_solved'],condition_M_max=b.mesh.condition,
        selected_array_bytes=sum(prop.get_npy_array().nbytes for prop in b.pa.properties.values())+b.K.nbytes+b.Lmap.nbytes+b.mesh.D.nbytes,
        memory_note='selected arrays NOT process peak RSS',source_sha256={s:sha(ROOT/s) for s in source},
        artifact_sha256={f.name:sha(f) for f in a.output.iterdir() if f.is_file()})
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','source_sha256','artifact_sha256')},indent=2)); print(json.dumps(records[-1],indent=2))
    return 0 if status=='TWO_TANGENT_CONE_PATH_COMPLETE' else 2


if __name__=='__main__':raise SystemExit(main())

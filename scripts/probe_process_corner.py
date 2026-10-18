"""Q4 coupled orthogonal normal contacts, NOT a rotating AFSD tool solver."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from probe_process_contact import ContactBlock, reduced_correction
from probe_process_plastic import PlasticBlock, CandidateFailure, ROOT
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG=ROOT/'configs/benchmarks/process_corner_q4.json'


def resolve_active_set(branch, initial, limits):
    """Same coupled active-set logic for the independent spring and 3D body."""
    active=np.asarray(initial,dtype=bool).copy(); trace=[]
    for _ in range(limits['maximum_active_sets']):
        c=branch(active)
        trace.append(dict(active=active.tolist(),iterations=c.get('iteration_trace',[])))
        tensile=active & (c['normal_force'] < -limits['contact_force_tolerance_N'])
        penetrating=(~active) & (c['gap'] < -limits['gap_tolerance_m'])
        if not tensile.any() and not penetrating.any():
            c.update(active_set_trace=trace,event_kind='REQUESTED')
            return c
        active=(active & ~tensile) | penetrating
    raise CandidateFailure('bounded corner active sets exhausted',c,trace)


def check_work(ww, wc, accumulated_wall, accumulated_material, wr, ws, tolerance):
    """Never allow opposite signed defects on separate walls to cancel."""
    values=np.r_[ww,wc,accumulated_wall,accumulated_material,wr,ws]
    if not np.isfinite(values).all():
        raise ValueError('nonfinite corner work')
    if (np.any(np.abs(ww-wc)>tolerance) or
            np.any(np.abs(accumulated_wall+ww-accumulated_material-wc)>tolerance) or
            abs(float(np.sum(wc))+wr-ws)>tolerance):
        raise ValueError('per-wall contact work requires event split; no residual heat')


class CornerBlock(ContactBlock):
    def __init__(self,law='j2',schedule='staggered'):
        super().__init__(law,3)
        self.cornercfg=json.loads(CONFIG.read_text(encoding='utf-8'))
        p=self.cornercfg['parameters']
        actual=dict(length=self.mesh.p['length'],wall_minimum=self.cap_displacement(.1),
                    peak_time=.1,end_time=.2,staggered_y_return_start=.15)
        if (self.cornercfg['production_authorized'] is not False or
                self.cornercfg['wall_axes'] != [0,1] or
                self.cornercfg['wall_normals'] != [[-1,0,0],[0,-1,0]] or
                schedule not in ('symmetric','staggered')):
            raise ValueError('Q4 declared scope/normal/program mismatch')
        for key,value in actual.items():
            if p[key]['value']!=value:
                raise ValueError('Q4 declared/used parameter mismatch: '+key)
        self.schedule=schedule
        self.nodes=np.concatenate([np.flatnonzero(self.X[:,a]==.01) for a in (0,1)])
        self.faces=np.repeat(np.arange(2),9)
        self.axes=self.faces.copy()
        self.contact_dofs=3*self.nodes+self.axes
        if len(np.unique(self.contact_dofs))!=len(self.contact_dofs) or np.any(~self.free[self.contact_dofs]):
            raise ValueError('duplicate or permanently constrained contact DOF')
        self.active=np.zeros(len(self.nodes),dtype=bool)
        self.normal_force=np.zeros(len(self.nodes))
        self.face_force=np.zeros((2,len(self.X),3))
        self.wall_work_by_face=np.zeros(2); self.material_work_by_face=np.zeros(2)
        self.last_wall_by_face=np.zeros(2); self.last_material_by_face=np.zeros(2)

    def wall_offsets(self,t):
        x=self.cap_displacement(t)
        if self.schedule=='symmetric' or t<=.1:
            y=x
        elif t<=.15:
            y=-4e-5
        else:
            y=-4e-5*(.2-t)/(.2-.15)
        return np.array([x,y])

    def _equilibrate(self,t,active,maximum_iterations=None):
        active=np.asarray(active,dtype=bool)
        if active.shape!=self.active.shape:
            raise ValueError('bad corner contact inventory')
        unknown=self.free.copy(); unknown[self.contact_dofs[active]]=False
        u=self.displacement.copy()
        u.ravel()[self.contact_dofs[active]]=self.wall_offsets(t)[self.faces[active]]
        maximum=self.qcfg['numerical']['maximum_iterations'] if maximum_iterations is None else maximum_iterations
        if maximum<0: raise ValueError('negative iteration limit')
        trace=[]
        for iteration in range(maximum+1):
            c=self.trial(u,t); force=c['internal_force'].ravel()
            norm=float(np.linalg.norm(force[unknown])/self.force_scale)
            trace.append(dict(iteration=iteration,free_force_relative=norm))
            if norm<=self.ccfg['numerical']['force_relative']:
                c.update(active=active.copy(),iterations=iteration,iteration_trace=trace)
                return self.contact_state(c)
            if iteration==maximum: break
            try:
                step,condition=reduced_correction(self.tangent(c),force,unknown,
                                                  self.ccfg['numerical']['tangent_condition_limit'])
            except (ValueError,np.linalg.LinAlgError) as exc:
                raise CandidateFailure(str(exc),c,trace) from exc
            trace[-1]['tangent_condition']=condition
            for halving in range(self.qcfg['numerical']['line_search_halvings']+1):
                test_u=u+2.**(-halving)*step.reshape(u.shape)
                probe=self.trial(test_u,t)
                if np.linalg.norm(probe['internal_force'].ravel()[unknown])/self.force_scale < norm:
                    u=test_u; trace[-1]['line_search_factor']=2.**(-halving); break
            else: raise CandidateFailure('corner branch line search failed',c,trace)
        raise CandidateFailure('corner branch equilibrium did not converge',c,trace)

    def contact_state(self,c):
        active=np.asarray(c['active'],dtype=bool)
        if active.shape!=self.active.shape: raise ValueError('bad corner contact inventory')
        normal=np.zeros(len(self.nodes))
        normal[active]=c['internal_force'].ravel()[self.contact_dofs[active]]
        per_face=np.zeros_like(self.face_force)
        for face in range(2):
            mask=self.faces==face
            np.add.at(per_face[face].ravel(),self.contact_dofs[mask],-normal[mask])
        contact=per_face.sum(axis=0)
        support=np.zeros_like(self.X); support[self.fixed]=-(c['internal_force']+contact)[self.fixed]
        gap=self.wall_offsets(c['time_s'])[self.faces]-c['displacement'].ravel()[self.contact_dofs]
        c.update(normal_force=normal,face_force=per_face,contact_force=contact,support_force=support,gap=gap,
                 equilibrium_relative=float(np.linalg.norm(c['internal_force']+contact+support)/self.force_scale))
        return c

    def solve_candidate(self,t,maximum_iterations=None):
        lim=self.ccfg['numerical']
        guess=self.wall_offsets(t)[self.faces]-self.displacement.ravel()[self.contact_dofs]
        initial=(self.active & (self.normal_force>lim['contact_force_tolerance_N'])) | (guess < -lim['gap_tolerance_m'])
        return resolve_active_set(lambda a:self._equilibrate(t,a,maximum_iterations),initial,lim)

    # Q3's _release_event locates the first loss of positive normal force across
    # all currently loaded contact rows; our branch shares ONE material state.
    def advance(self,t,on_accept=None):
        knots=(.1,.15) if self.schedule=='staggered' else (.1,)
        if any(self.time<k-1e-14 and t>k+1e-14 for k in knots):
            raise ValueError('requested interval crosses a wall program corner')
        for _ in range(self.ccfg['numerical']['maximum_events_per_requested_interval']):
            c=self.solve_candidate(t)
            loaded=self.active & (self.normal_force>self.ccfg['numerical']['contact_force_tolerance_N'])
            if np.any(loaded & ~c['active']): c=self._release_event(t,loaded)
            self.commit(c)
            if on_accept is not None: on_accept()
            if self.time==t: return
        raise CandidateFailure('bounded corner release events exhausted')

    def commit(self,candidate):
        if candidate['base_step']!=self.count: raise ValueError('stale candidate')
        c=self.trial(candidate['displacement'],candidate['time_s'])
        c['active']=np.asarray(candidate['active'],dtype=bool).copy(); c=self.contact_state(c)
        lim=self.ccfg['numerical']
        if (c['equilibrium_relative']>lim['force_relative'] or
                np.any(c['normal_force'] < -lim['contact_force_tolerance_N']) or
                np.any(c['gap'] < -lim['gap_tolerance_m']) or
                np.any(np.abs(c['gap'][c['active']])>lim['gap_tolerance_m'])):
            raise CandidateFailure('corner equilibrium/complementarity rejected',c)
        if np.max(np.linalg.norm(c['displacement'],axis=1))/.01 >= self.qcfg['numerical']['displacement_over_length_limit']:
            raise CandidateFailure('left declared small-strain domain',c)
        du=c['displacement']-self.displacement
        dd=self.wall_offsets(c['time_s'])-self.wall_offsets(self.time)
        wc=np.sum(.5*(self.face_force+c['face_force'])*du[None],axis=(1,2))
        normal_mid=.5*(self.normal_force+c['normal_force'])
        ww=-np.array([normal_mid[self.faces==face].sum() for face in range(2)])*dd
        wr=float(np.sum(.5*(self.support_force+c['support_force'])*du))
        ws=float(np.einsum('i,iab,iab->',self.V,.5*(self.stress+c['stress']),self.mesh.gradient(du)))
        try: check_work(ww,wc,self.wall_work_by_face,self.material_work_by_face,wr,ws,lim['work_tolerance_J'])
        except ValueError as exc: raise CandidateFailure(str(exc),c) from exc
        kind=candidate.get('event_kind','DIRECT_COMMIT'); event=None
        if kind not in ('REQUESTED','RELEASE_LIMIT','DIRECT_COMMIT'): raise ValueError('unknown event metadata')
        if kind=='RELEASE_LIMIT':
            bracket=np.asarray(candidate.get('event_bracket_s',[]),dtype=float)
            if (bracket.shape!=(2,) or not np.isfinite(bracket).all() or bracket[0]!=c['time_s'] or
                    bracket[1]<bracket[0] or bracket[1]>.2 or bracket[1]-bracket[0]>lim['event_time_tolerance_s']):
                raise ValueError('invalid release event bracket metadata')
            at_limit=c['active'] & (c['normal_force']<=lim['contact_force_tolerance_N'])
            event=dict(step=self.count+1,physical_time_s=c['time_s'],bracket_s=bracket.tolist(),kind=kind,
                near_zero_contact_rows=np.flatnonzero(at_limit).tolist(),
                near_zero_feature_ids=(np.unique(self.faces[at_limit])+1).tolist())
        # One shared commit; no wall is permitted to advance material separately.
        self.displacement=c['displacement'].copy(); self.stress=c['stress']
        self.pa.x[:],self.pa.y[:],self.pa.z[:]=(self.X+self.displacement).T
        self.pa.S[:]=c['S'].ravel(); self.pa.ep[:]=c['ep']
        self.pa.plastic_work[:]+=c['plastic_increment_J']; self.pa.return_dissipation[:]+=c['return_increment_J']
        self.pa.strain_rate[:]=c['strain_rate']; self.rate=du/c['dt_s']
        self.time=c['time_s']; self.count+=1; self.active=c['active']; self.normal_force=c['normal_force']
        self.face_force=c['face_force']; self.contact_force=c['contact_force']; self.support_force=c['support_force']
        self.reaction=self.contact_force+self.support_force
        self.wall_work_by_face+=ww; self.material_work_by_face+=wc
        self.last_wall_by_face=ww; self.last_material_by_face=wc
        self.wall_work=float(self.wall_work_by_face.sum()); self.contact_work=float(self.material_work_by_face.sum())
        self.support_work+=wr; self.internal_work+=ws; self.constraint_work=self.contact_work+self.support_work
        self.last_residual=c['equilibrium_relative']; self.last_event=kind
        self.last_iterations=candidate.get('iterations',0)
        self.last_trace=candidate.get('active_set_trace',candidate.get('iteration_trace',[]))
        force=c['internal_force']
        self.last_force_defect=float(np.linalg.norm(force.sum(axis=0))/self.force_scale)
        self.last_torque_defect=float(np.linalg.norm(np.cross(self.X,force).sum(axis=0))/(self.force_scale*.01))
        if event is not None: self.events.append(event)

    def record(self):
        row=PlasticBlock.record(self)
        for name in ('actuator_work_J','cap_reaction_x_N','cap_displacement_m'): row.pop(name)
        gap=self.wall_offsets(self.time)[self.faces]-self.displacement.ravel()[self.contact_dofs]
        row.update(wall_displacement_m=self.wall_offsets(self.time).tolist(),
            wall_force_on_material_N=self.face_force.sum(axis=1).tolist(),
            wall_work_J=self.wall_work_by_face.tolist(),contact_material_work_J=self.material_work_by_face.tolist(),
            interval_wall_work_J=self.last_wall_by_face.tolist(),interval_contact_material_work_J=self.last_material_by_face.tolist(),
            wall_minus_material_work_J=(self.wall_work_by_face-self.material_work_by_face).tolist(),
            support_work_J=self.support_work,active_contacts_by_wall=[int(self.active[self.faces==f].sum()) for f in range(2)],
            minimum_gap_by_wall_m=[float(gap[self.faces==f].min()) for f in range(2)],
            maximum_gap_by_wall_m=[float(gap[self.faces==f].max()) for f in range(2)],
            minimum_normal_force_by_wall_N=[float(self.normal_force[self.faces==f].min()) for f in range(2)],
            maximum_abs_complementarity_J=float(np.max(np.abs(self.normal_force*gap))),event_kind=self.last_event)
        return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--law',choices=('elastic','j2'),required=True)
    p.add_argument('--schedule',choices=('symmetric','staggered'),required=True)
    p.add_argument('--increments-per-phase',type=int,choices=(8,16),required=True)
    p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False); start=time.perf_counter()
    b=CornerBlock(a.law,a.schedule); setup=time.perf_counter()-start
    frames=[]; records=[]
    def capture():
        records.append(b.record())
        frames.append(dict(position=b.X+b.displacement,u=b.displacement.copy(),stress=b.stress.copy(),
            ep=b.pa.ep.copy(),plastic=b.pa.plastic_work.copy(),returned=b.pa.return_dissipation.copy(),
            face_force=b.face_force.copy(),support=b.support_force.copy(),normal=b.normal_force.copy(),
            active=b.active.copy(),gap=b.wall_offsets(b.time)[b.faces]-b.displacement.ravel()[b.contact_dofs],rate=b.rate.copy()))
    capture(); status='CORNER_PATH_COMPLETE'; error=None; tick=time.perf_counter()
    try:
        for t in np.linspace(0.,.2,2*a.increments_per_phase+1)[1:]: b.advance(float(t),capture)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        status='BLOCKED'; error=str(exc); c=getattr(exc,'candidate',None)
        if c is not None:
            np.savez_compressed(a.output/'REJECTED_candidate.npz',**{k:v for k,v in c.items() if isinstance(v,np.ndarray)})
        rejection=dict(classification='REJECTED',restart_allowed=False,error=error,attempted_time_s=float(t),
            last_accepted_step=b.count,candidate_available=c is not None,
            candidate_time_s=None if c is None else c.get('time_s'),candidate_dt_s=None if c is None else c.get('dt_s'),
            candidate_base_step=None if c is None else c.get('base_step'),iterations=getattr(exc,'iterations',[]))
        (a.output/'REJECTED_manifest.json').write_text(json.dumps(rejection,indent=2)+'\n',encoding='utf-8')
    advance=time.perf_counter()-tick
    config=dict(**b.cornercfg,selected_law=a.law,schedule=a.schedule,increments_per_phase=a.increments_per_phase,
                geometry_elastic_inputs=b.mesh.cfg,plastic_parent=b.qcfg,contact_parent=b.ccfg)
    (a.output/'resolved_config.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(a.output/'history.npz',particle_id=b.pa.particle_id.copy(),tracer_type=b.pa.tracer_type.copy(),
        reference_position_m=b.X,reference_volume_m3=b.V,mass_kg=b.pa.m.copy(),support_components=b.fixed,
        contact_particle_id=b.pa.particle_id[b.nodes],contact_feature_id=b.faces+1,contact_axis=b.axes,
        contact_normal=-np.eye(3)[b.axes],physical_time_s=np.array([r['physical_time_s'] for r in records]),
        position_m=np.array([f['position'] for f in frames]),displacement_m=np.array([f['u'] for f in frames]),
        stress_Pa=np.array([f['stress'] for f in frames]),equivalent_plastic_strain=np.array([f['ep'] for f in frames]),
        plastic_work_J=np.array([f['plastic'] for f in frames]),return_dissipation_diagnostic_J=np.array([f['returned'] for f in frames]),
        contact_force_by_wall_N=np.array([f['face_force'] for f in frames]),support_force_N=np.array([f['support'] for f in frames]),
        normal_force_N=np.array([f['normal'] for f in frames]),active_contacts=np.array([f['active'] for f in frames]),
        gap_m=np.array([f['gap'] for f in frames]),kinematic_interval_rate_m_s=np.array([f['rate'] for f in frames]))
    checkpoint(a.output/'checkpoint.npz',b.pa,dict(classification=b.cornercfg['classification'],restart_allowed=False,
        physical_time_s=b.time,step=b.count,production_authorized=False,not_solved=b.cornercfg['not_solved']))
    scalar=dict(particle_id=b.pa.particle_id.tolist(),tracer_type=b.pa.tracer_type.tolist(),
                equivalent_plastic_strain=b.pa.ep.tolist(),plastic_work_J=b.pa.plastic_work.tolist())
    for i,j,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        scalar['stress_'+name+'_Pa']=b.stress[:,i,j].tolist()
    write_point_cloud(a.output/'final.vtk',(b.X+b.displacement).tolist(),scalar,
        dict(reference_position_m=b.X.tolist(),displacement_m=b.displacement.tolist()),
        title='Q4_ACCEPTED_TWO_WALL_BENCHMARK_NOT_AFSD_THERMAL_NOT_SOLVED')
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    source=['scripts/probe_process_corner.py','configs/benchmarks/process_corner_q4.json',
        'scripts/probe_process_contact.py','scripts/probe_process_plastic.py','scripts/probe_process_time.py',
        'scripts/probe_engineering_material.py','scripts/probe_engineering_block.py','scripts/probe_boundary_sbp.py',
        'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']
    report=dict(status=status,error=error,classification=b.cornercfg['classification'],production_authorized=False,
        selected_law=a.law,schedule=a.schedule,particles=len(b.X),contact_rows=len(b.nodes),
        permanent_free_dofs=int(b.free.sum()),increments_per_phase=a.increments_per_phase,
        accepted_steps=b.count,accepted_time_s=b.time,requested_end_s=.2,
        setup_wall_s=setup,advance_wall_s=advance,wall_before_report_s=time.perf_counter()-start,
        selected_array_bytes=sum(prop.get_npy_array().nbytes for prop in b.pa.properties.values())+b.mesh.D.nbytes+b.Lmap.nbytes+b.K.nbytes+b.face_force.nbytes,
        memory_note='selected retained arrays, NOT transient solver buffers or process RSS',
        condition_M_max=b.mesh.condition,events=b.events,records=records,not_solved=b.cornercfg['not_solved'],
        source_sha256={s:sha(ROOT/s) for s in source},artifact_sha256={f.name:sha(f) for f in a.output.iterdir() if f.is_file()})
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','source_sha256','artifact_sha256')},indent=2))
    print(json.dumps({k:v for k,v in records[-1].items() if k!='iteration_trace'},indent=2))
    return 0 if status=='CORNER_PATH_COMPLETE' else 2


if __name__=='__main__':raise SystemExit(main())

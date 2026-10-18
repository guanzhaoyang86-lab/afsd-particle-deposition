"""Q7 existing J2 + vector cone, fixed-normal small block; NOT AFSD/FKS."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from probe_process_friction_cone import ConeBlock, feasible_cone, vector_contact_work
from probe_process_plastic import CandidateFailure, ROOT
from probe_engineering_material import PerfectPlasticBenchmark
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import write_point_cloud

CONFIG=ROOT/'configs/benchmarks/process_plastic_friction_q7.json'
SOURCES=['scripts/probe_process_plastic_friction.py','configs/benchmarks/process_plastic_friction_q7.json',
    'scripts/probe_process_friction_cone.py','configs/benchmarks/process_friction_cone_q6.json','scripts/probe_process_friction.py',
    'configs/benchmarks/process_friction_q5.json','scripts/probe_process_contact.py','configs/benchmarks/process_contact_q3.json',
    'scripts/probe_process_plastic.py','configs/benchmarks/process_plastic_q2.json','scripts/probe_process_time.py',
    'scripts/probe_engineering_material.py','scripts/probe_engineering_block.py','scripts/probe_boundary_sbp.py',
    'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class PlasticFrictionBlock(ConeBlock):
    def __init__(self):
        super().__init__(0.)
        self.pfcfg=json.loads(CONFIG.read_text(encoding='utf-8'))
        p=self.pfcfg['parameters']
        expected=dict(yield_stress=90e6,temperature=293.15,young=70e9,poisson=0.,density=2700.,
            mu=.1,normal_displacement=-4e-5,tangent_amplitude=2e-5,phase_duration=.05)
        if (self.pfcfg['production_authorized'] is not False or
                any(p[k]['value']!=v for k,v in expected.items()) or
                p['tangent_direction_global_yz']['value']!=[.6,.8]):
            raise ValueError('Q7 declared manufactured physical parameters changed')
        self.law=PerfectPlasticBenchmark(p['yield_stress']['value'],p['temperature']['value'])
        self.law_name='j2'
        self.last_plastic=np.zeros(len(self.X)); self.last_return=np.zeros(len(self.X))

    def commit(self,candidate):
        """One transaction: never call elastic parent commit then add plasticity."""
        if candidate['base_step']!=self.count:raise ValueError('stale candidate')
        c=self.state(self.trial(candidate['displacement'],candidate['time_s']),candidate['released'])
        if (not feasible_cone(c,np.ones(len(self.nodes),dtype=bool),self.mu,self.limits) or
                np.max(np.linalg.norm(c['law_residual_N'],axis=1))>self.limits['force_tolerance_N']):
            raise CandidateFailure('cone contact law/root rejected',c)
        if (c['equilibrium_relative']>self.ccfg['numerical']['force_relative'] or
                np.any(np.abs(c['gap'])>self.limits['slip_tolerance_m'])):
            raise CandidateFailure('cone equilibrium/normal constraint rejected',c)
        if (np.any(c['ep']<self.pa.ep) or np.any(c['plastic_increment_J']<0) or
                np.any(c['return_increment_J']<0)):
            raise CandidateFailure('negative raw plastic increment/work',c)
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
        plastic=self.pa.plastic_work+c['plastic_increment_J']
        returned=self.pa.return_dissipation+c['return_increment_J']
        if not np.isfinite(np.r_[q,wc,ww,wr,ws,plastic,returned,c['ep']]).all():
            raise CandidateFailure('nonfinite coupled work/material',c)
        if np.any(q < -self.limits['negative_dissipation_tolerance_J']):
            raise CandidateFailure('negative raw cone friction dissipation',c)
        if max(abs(ww-wc-q.sum()),abs(wc+wr-ws),
               abs(self.wall_work+ww-self.contact_work-wc-(self.friction.sum()+q.sum())),
               abs(self.contact_work+wc+self.support_work+wr-self.internal_work-ws))>self.limits['work_tolerance_J']:
            raise CandidateFailure('coupled contact work identity failed',c)
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
        # All validation and ledger candidates computed before any mutation.
        self.displacement=c['displacement'].copy(); self.stress=c['stress']
        self.pa.S[:]=c['S'].ravel(); self.pa.ep[:]=c['ep']
        self.pa.plastic_work[:]=plastic; self.pa.return_dissipation[:]=returned
        self.last_plastic=c['plastic_increment_J']; self.last_return=c['return_increment_J']
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
        row.update(interval_plastic_work_J=float(self.last_plastic.sum()),
            interval_return_diagnostic_J=float(self.last_return.sum()),
            material_work_warning=abs(row['work_minus_elastic_minus_plastic_J'])>
                self.pfcfg['numerical']['energy_diagnostic_warning_J'])
        return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--increments-per-phase',type=int,choices=(4,8),required=True)
    p.add_argument('--output',type=Path,required=True); a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
    provenance=dict(classification='Q7_BOUNDED_ATTEMPT_NOT_PRODUCTION',increments_per_phase=a.increments_per_phase,
        source_sha256={s:sha(ROOT/s) for s in SOURCES},configuration_sha256=sha(CONFIG))
    (a.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    start=time.perf_counter(); b=PlasticFrictionBlock(); setup=time.perf_counter()-start
    frames=[]; records=[]
    def capture():
        records.append(b.record()); frames.append(dict(position=b.X+b.displacement,u=b.displacement.copy(),stress=b.stress.copy(),
            normal=b.normal_force.copy(),tau=b.tangent_force.copy(),modes=b.modes.copy(),force=b.contact_force.copy(),support=b.support_force.copy(),
            q=b.friction.copy(),dq=b.last_friction.copy(),slip=b.last_slip.copy(),impulse=b.last_impulse.copy(),relrate=b.last_relative_rate.copy(),
            rate=b.rate.copy(),lawres=b.last_law_residual.copy(),ep=b.pa.ep.copy(),plastic=b.pa.plastic_work.copy(),
            returned=b.pa.return_dissipation.copy(),dplastic=b.last_plastic.copy(),dreturn=b.last_return.copy(),strain_rate=b.pa.strain_rate.copy()))
        # Unique accepted-only journal survives an external timeout. A partially
        # written last file is invalid, never promoted; prior snapshots remain.
        np.savez_compressed(a.output/f'accepted_{b.count:04d}.npz',**frames[-1],
            particle_id=b.pa.particle_id.copy(),tracer_type=b.pa.tracer_type.copy(),
            physical_time_s=b.time,step=b.count,record_json=json.dumps(records[-1],allow_nan=False),
            classification='ACCEPTED_DIAGNOSTIC_NOT_RESTART',restart_allowed=False)
    capture(); tick=time.perf_counter(); status='J2_CONE_PATH_COMPLETE'; error=None
    try:
        for t in np.linspace(0.,.2,4*a.increments_per_phase+1)[1:]:b.advance(float(t),capture)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        status='BLOCKED'; error=str(exc); c=getattr(exc,'candidate',None)
        if c is not None:np.savez_compressed(a.output/'REJECTED_candidate.npz',**{k:v for k,v in c.items() if isinstance(v,np.ndarray)})
        detail=dict(classification='REJECTED',restart_allowed=False,error=error,attempted_time_s=float(t),last_accepted_step=b.count,
            last_accepted_time_s=b.time,candidate_available=c is not None,candidate_time_s=None if c is None else c.get('time_s'),
            candidate_base_step=None if c is None else c.get('base_step'),iterations=getattr(exc,'iterations',[]))
        (a.output/'REJECTED_manifest.json').write_text(json.dumps(detail,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    advance=time.perf_counter()-tick
    config=dict(**b.pfcfg,increments_per_phase=a.increments_per_phase,basis_angle_deg=0.,
        cone_parent=b.conecfg,friction_parent=b.fcfg,geometry_inputs=b.mesh.cfg,plastic_parent=b.qcfg,contact_parent=b.ccfg)
    (a.output/'resolved_config.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    data=dict(particle_id=b.pa.particle_id.copy(),tracer_type=b.pa.tracer_type.copy(),reference_position_m=b.X,reference_volume_m3=b.V,
        mass_kg=b.pa.m.copy(),support_components=b.fixed,contact_particle_id=b.pa.particle_id[b.nodes],contact_feature_id=np.ones(len(b.nodes),dtype=int),
        physical_time_s=np.array([r['physical_time_s'] for r in records]))
    mapping=dict(position_m='position',displacement_m='u',stress_Pa='stress',normal_force_N='normal',tangent_force_global_yz_N='tau',
        interval_modes='modes',contact_force_N='force',support_force_N='support',friction_dissipation_J='q',interval_friction_J='dq',
        interval_relative_slip_global_yz_m='slip',tangent_impulse_global_yz_N_s='impulse',relative_kinematic_rate_global_yz_m_s='relrate',
        kinematic_interval_rate_m_s='rate',contact_law_residual_global_yz_N='lawres',equivalent_plastic_strain='ep',
        plastic_work_J='plastic',return_dissipation_diagnostic_J='returned',interval_plastic_work_J='dplastic',
        interval_return_diagnostic_J='dreturn',equivalent_strain_rate_per_s='strain_rate')
    data.update({key:np.array([f[field] for f in frames]) for key,field in mapping.items()}); np.savez_compressed(a.output/'history.npz',**data)
    checkpoint(a.output/'checkpoint.npz',b.pa,dict(restart_allowed=False,classification=b.pfcfg['classification'],physical_time_s=b.time,step=b.count,production_authorized=False))
    scalar=dict(particle_id=b.pa.particle_id.tolist(),tracer_type=b.pa.tracer_type.tolist(),equivalent_plastic_strain=b.pa.ep.tolist(),
        plastic_work_J=b.pa.plastic_work.tolist(),return_dissipation_diagnostic_J=b.pa.return_dissipation.tolist())
    for i,j,name in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):scalar['stress_'+name+'_Pa']=b.stress[:,i,j].tolist()
    write_point_cloud(a.output/'final.vtk',(b.X+b.displacement).tolist(),scalar,dict(reference_position_m=b.X.tolist(),displacement_m=b.displacement.tolist()),title='Q7_ACCEPTED_J2_VECTOR_CONE_NOT_AFSD_THERMAL_NOT_SOLVED')
    report=dict(status=status,error=error,classification=b.pfcfg['classification'],production_authorized=False,particles=len(b.X),contact_rows=len(b.nodes),
        permanent_free_dofs=int(b.free.sum()),increments_per_phase=a.increments_per_phase,basis_angle_deg=0.,
        accepted_steps=b.count,accepted_time_s=b.time,requested_end_s=.2,setup_wall_s=setup,advance_wall_s=advance,
        wall_before_report_s=time.perf_counter()-start,events=b.events,records=records,not_solved=b.pfcfg['not_solved'],condition_M_max=b.mesh.condition,
        material_work_warning_steps=[r['step'] for r in records if r['material_work_warning']],
        selected_array_bytes=sum(prop.get_npy_array().nbytes for prop in b.pa.properties.values())+b.K.nbytes+b.Lmap.nbytes+b.mesh.D.nbytes,
        memory_note='selected arrays NOT process peak RSS',source_sha256=provenance['source_sha256'],
        artifact_sha256={f.name:sha(f) for f in a.output.iterdir() if f.is_file()})
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','source_sha256','artifact_sha256')},indent=2)); print(json.dumps(records[-1],indent=2))
    return 0 if status=='J2_CONE_PATH_COMPLETE' else 2


if __name__=='__main__':raise SystemExit(main())

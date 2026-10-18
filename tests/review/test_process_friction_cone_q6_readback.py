"""Q6 reader attacks on explicitly synthetic I/O fixtures, NOT mechanics runs."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from inspect_process_friction_cone import inspect, wall_reference
from afsd_scaffold.vtk_legacy import read_point_cloud, write_point_cloud
from afsd_backend.pysph.output import checkpoint
from pysph.base.utils import get_particle_array

SOURCES = ['scripts/probe_process_friction_cone.py','configs/benchmarks/process_friction_cone_q6.json',
    'scripts/probe_process_friction.py','configs/benchmarks/process_friction_q5.json',
    'scripts/probe_process_contact.py','scripts/probe_process_plastic.py','scripts/probe_process_time.py',
    'scripts/probe_engineering_material.py','scripts/probe_engineering_block.py','scripts/probe_boundary_sbp.py',
    'src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_schema_fixture(directory):
    """Three fake accepted frames exercise serialization/accounting only.

    Stress/internal force are NOT solved. A BLOCKED manifest and explicit
    evidence_kind prevent this artifact being represented as a physics run.
    Independent literal expected Q=9*.1*5e-6=4.5e-6 J is just vector work.
    """
    p = Path(directory)/'SCHEMA_FIXTURE_NOT_PHYSICAL'; p.mkdir()
    X = np.array([[x,y,z] for x in (0.,.005,.01) for y in (0.,.005,.01) for z in (0.,.005,.01)])
    nodes = np.flatnonzero(X[:,0]==.01); support_nodes = np.flatnonzero(X[:,0]==0)
    ids = 101+np.arange(27); tags = np.arange(27)%5
    t = np.array([0.,.05,.0625]); wall = np.array([wall_reference(a) for a in t])
    u = np.zeros((3,27,3)); u[:,:,0] = wall[:,0,None]*(X[:,0]/.01)[None]
    fixed = np.broadcast_to((X[:,0]==0)[:,None],X.shape).copy()
    V = .005**3*np.prod(np.where((X==0)|(X==.01),.5,1.),axis=1)
    normal = np.zeros((3,9)); normal[1:]=1.
    tau = np.zeros((3,9,2)); tau[1:] = [.06,.08]
    modes = np.zeros((3,9),dtype=int); modes[2]=1
    force = np.zeros((3,27,3)); force[:,nodes,0]=-normal; force[:,nodes,1:]=tau
    support = np.zeros_like(force); support[:,support_nodes] = -force[:,nodes]
    stress = np.zeros((3,27,3,3))
    slip = np.zeros((3,9,2)); impulse=np.zeros_like(slip); relrate=np.zeros_like(slip)
    rate=np.zeros_like(u); dq=np.zeros((3,9)); Q=np.zeros_like(dq)
    rows=[]; Ww=Wc=Wr=0.
    for i,time in enumerate(t):
        ww=wc=wr=0.
        if i:
            dt=t[i]-t[i-1]; du=u[i]-u[i-1]; dw=wall[i]-wall[i-1]
            mid=.5*(force[i]+force[i-1])
            slip[i]=du[nodes,1:]-dw[1:]; impulse[i]=mid[nodes,1:]*dt
            relrate[i]=slip[i]/dt; rate[i]=du/dt
            dq[i]=-np.einsum('ij,ij->i',impulse[i],relrate[i])
            Q[i]=Q[i-1]+dq[i]
            ww=float(mid.sum(axis=0)@dw); wc=float(np.sum(mid*du))
            wr=float(np.sum(.5*(support[i]+support[i-1])*du))
            Ww+=ww; Wc+=wc; Wr+=wr
        rows.append(dict(step=i,physical_time_s=float(time),mode_counts={str(m):int((modes[i]==m).sum()) for m in (0,1)},
            wall_displacement_m=wall[i].tolist(),wall_force_on_material_N=force[i].sum(axis=0).tolist(),
            interval_wall_work_J=ww,interval_contact_material_work_J=wc,wall_work_J=Ww,
            contact_material_work_J=Wc,support_work_J=Wr,stress_work_J=Wc+Wr,
            friction_dissipation_J=float(Q[i].sum()),interval_friction_J=float(dq[i].sum()),
            minimum_raw_contact_dissipation_J=float(dq[i].min()),
            friction_work_identity_J=Ww-Wc-float(Q[i].sum()),elastic_diagnostic_J=0.,plastic_work_J=0.,maximum_ep=0.,
            maximum_contact_law_residual_N=0.,minimum_cone_margin_N=0.,
            event_kind='CONE_LIMIT' if i==1 else 'REQUESTED' if i else 'INITIAL'))
    data=dict(particle_id=ids,tracer_type=tags,reference_position_m=X,reference_volume_m3=V,
        mass_kg=2700*V,support_components=fixed,contact_particle_id=ids[nodes],
        contact_feature_id=np.ones(9,dtype=int),physical_time_s=t,equivalent_plastic_strain=np.zeros(27),
        position_m=X[None]+u,displacement_m=u,stress_Pa=stress,normal_force_N=normal,tangent_force_global_yz_N=tau,
        interval_modes=modes,contact_force_N=force,support_force_N=support,friction_dissipation_J=Q,
        interval_friction_J=dq,interval_relative_slip_global_yz_m=slip,tangent_impulse_global_yz_N_s=impulse,
        relative_kinematic_rate_global_yz_m_s=relrate,kinematic_interval_rate_m_s=rate,
        contact_law_residual_global_yz_N=np.zeros_like(slip))
    np.savez_compressed(p/'history.npz',**data)
    cfg=json.loads((ROOT/'configs/benchmarks/process_friction_cone_q6.json').read_text(encoding='utf-8'))
    cfg.update(increments_per_phase=4,basis_angle_deg=0.,tangent_basis_global_yz=np.eye(2).tolist(),
               evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
    (p/'resolved_config.json').write_text(json.dumps(cfg),encoding='utf-8')
    pa=get_particle_array(name='solid',x=data['position_m'][-1,:,0],y=data['position_m'][-1,:,1],z=data['position_m'][-1,:,2])
    for key,values in (('particle_id',ids),('tracer_type',tags)):
        pa.add_property(key,type='long',data=values)
    pa.add_property('temperature',data=np.full(27,293.15));pa.add_property('ep',data=np.zeros(27))
    checkpoint(p/'checkpoint.npz',pa,dict(restart_allowed=False,production_authorized=False,physical_time_s=float(t[-1]),step=2,
        classification=cfg['classification'],evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL'))
    scalars=dict(particle_id=ids.tolist(),tracer_type=tags.tolist(),equivalent_plastic_strain=np.zeros(27).tolist())
    for name in ('xx','yy','zz','xy','xz','yz'):scalars['stress_'+name+'_Pa']=np.zeros(27).tolist()
    write_point_cloud(p/'final.vtk',data['position_m'][-1].tolist(),scalars,
        dict(reference_position_m=X.tolist(),displacement_m=u[-1].tolist()),title='SCHEMA_FIXTURE_NOT_PHYSICAL')
    rejected=dict(classification='REJECTED',restart_allowed=False,error='SYNTHETIC_SCHEMA_ONLY_NOT_PHYSICS',
        last_accepted_step=2,candidate_available=False,candidate_time_s=None,candidate_base_step=None,attempted_time_s=.075)
    (p/'REJECTED_manifest.json').write_text(json.dumps(rejected),encoding='utf-8')
    report=dict(status='BLOCKED',error=rejected['error'],production_authorized=False,classification=cfg['classification'],
        evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL',particles=27,contact_rows=9,permanent_free_dofs=54,
        increments_per_phase=4,basis_angle_deg=0.,accepted_steps=2,accepted_time_s=float(t[-1]),requested_end_s=.2,
        records=rows,not_solved=cfg['not_solved'],
        events=[dict(step=1,physical_time_s=.05,kind='CONE_LIMIT',bracket_s=[.05,.05+5e-13],
            ending_interval_modes=modes[1].tolist(),released_rows=[False]*9,cone_margin_N=[0.]*9)],
        source_sha256={name:digest(ROOT/name) for name in SOURCES},
        artifact_sha256={f.name:digest(f) for f in p.iterdir() if f.is_file()})
    (p/'report.json').write_text(json.dumps(report),encoding='utf-8')
    return p


class ConeReadbackQ6(unittest.TestCase):
    def report(self,p):return json.loads((p/'report.json').read_text(encoding='utf-8'))
    def save_report(self,p,r):(p/'report.json').write_text(json.dumps(r),encoding='utf-8')
    def rehash(self,p,name):
        r=self.report(p);r['artifact_sha256'][name]=digest(p/name);self.save_report(p,r)
    def change(self,p,func):
        with np.load(p/'history.npz',allow_pickle=False) as saved:h={k:saved[k].copy() for k in saved.files}
        func(h);np.savez_compressed(p/'history.npz',**h);self.rehash(p,'history.npz')

    def test_01_explicit_schema_fixture_is_only_partial_with_independent_vector_work(self):
        with tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
            p=make_schema_fixture(d);r=inspect(p)
            self.assertEqual(r['status'],'PASS_PARTIAL_READBACK')
            self.assertEqual(r['completion_status'],'BLOCKED')
            self.assertAlmostEqual(sum(r['independent_friction_dissipation_J']),4.5e-6,delta=1e-18)
            self.assertFalse(r['production_authorized'])

    def test_02_hash_and_rehashed_identity_inventory(self):
        for kind,message in [('artifact','artifact changed'),('source','source changed'),('id','identity')]:
            with self.subTest(kind=kind),tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
                p=make_schema_fixture(d)
                if kind=='artifact':(p/'final.vtk').write_bytes((p/'final.vtk').read_bytes()+b'\nTAMPER\n')
                elif kind=='source':
                    r=self.report(p);r['source_sha256'][SOURCES[0]]='0'*64;self.save_report(p,r)
                else:self.change(p,lambda h:h['particle_id'].__setitem__(1,h['particle_id'][0]))
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_03_circular_cone_not_separate_component_bounds(self):
        with tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
            p=make_schema_fixture(d)
            def modify(h):
                # Each component .08<.1, but Euclidean norm exceeds .1.
                h['tangent_force_global_yz_N'][-1,0]=[.08,.08]
                node=np.flatnonzero(h['reference_position_m'][:,0]==.01)[0]
                h['contact_force_N'][-1,node,1:]=[.08,.08]
            self.change(p,modify)
            with self.assertRaisesRegex(ValueError,'circular Coulomb'):inspect(p)

    def test_04_traction_norm_alone_cannot_hide_wrong_vector_direction(self):
        with tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
            p=make_schema_fixture(d)
            def modify(h):
                h['tangent_force_global_yz_N'][-1,0]=[-.06,-.08]
                node=np.flatnonzero(h['reference_position_m'][:,0]==.01)[0]
                h['contact_force_N'][-1,node,1:]=[-.06,-.08]
            self.change(p,modify)
            r=self.report(p);r['records'][-1]['wall_force_on_material_N'][1]-=.12;r['records'][-1]['wall_force_on_material_N'][2]-=.16
            self.save_report(p,r)
            with self.assertRaisesRegex(ValueError,'vector slip opposition'):inspect(p)

    def test_05_per_contact_errors_cannot_cancel_in_sum(self):
        for field,message in [('interval_friction_J','per-contact friction work'),('friction_dissipation_J','per-contact accumulated')]:
            with self.subTest(field=field),tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
                p=make_schema_fixture(d)
                def modify(h):
                    before=h[field][-1].sum();h[field][-1,0]-=.01;h[field][-1,1]+=.01
                    self.assertAlmostEqual(float(before),float(h[field][-1].sum()),delta=1e-17)
                self.change(p,modify)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_06_actual_mode_root_and_relative_motion_ledger(self):
        for field,message in [('interval_modes','physical slip mode'),('contact_law_residual_global_yz_N','root ledger'),
                              ('interval_relative_slip_global_yz_m','relative slip')]:
            with self.subTest(field=field),tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
                p=make_schema_fixture(d)
                self.change(p,lambda h:h[field].__setitem__((-1,0),0 if field=='interval_modes' else [1.,0.]))
                if field=='interval_modes':
                    r=self.report(p);r['records'][-1]['mode_counts']={'0':1,'1':8};self.save_report(p,r)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_07_event_inventory_and_released_mask(self):
        for kind in ('missing','duplicate','wrong_margin','wrong_released'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
                p=make_schema_fixture(d);r=self.report(p)
                if kind=='missing':r['events']=[]
                elif kind=='duplicate':r['events']*=2
                elif kind=='wrong_margin':r['events'][0]['cone_margin_N'][0]=1.
                else:r['events'][0]['released_rows']=[1]*9
                self.save_report(p,r)
                with self.assertRaisesRegex(ValueError,'event'):inspect(p)

    def test_08_unsolved_and_completion_claims_rejected(self):
        for kind,message in [('temperature','unsolved field'),('production','cannot authorize production'),('complete','terminal/requested')]:
            with self.subTest(kind=kind),tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
                p=make_schema_fixture(d)
                if kind=='temperature':
                    self.change(p,lambda h:h.update(temperature_K=np.full((3,27),293.15)))
                else:
                    r=self.report(p)
                    if kind=='production':r['production_authorized']=True
                    else:r['status']='TWO_TANGENT_CONE_PATH_COMPLETE'
                    self.save_report(p,r)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_09_basis_rotation_is_representation_only_and_checked(self):
        with tempfile.TemporaryDirectory(prefix='q6_reader_') as d:
            p=make_schema_fixture(d);r=self.report(p);r['basis_angle_deg']=37.;self.save_report(p,r)
            cfg=json.loads((p/'resolved_config.json').read_text(encoding='utf-8'));a=np.deg2rad(37.)
            cfg.update(basis_angle_deg=37.,tangent_basis_global_yz=[[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
            (p/'resolved_config.json').write_text(json.dumps(cfg),encoding='utf-8');self.rehash(p,'resolved_config.json')
            self.assertEqual(inspect(p)['status'],'PASS_PARTIAL_READBACK')
            cfg['tangent_basis_global_yz'][0][0]=2.
            (p/'resolved_config.json').write_text(json.dumps(cfg),encoding='utf-8');self.rehash(p,'resolved_config.json')
            with self.assertRaisesRegex(ValueError,'tangent basis'):inspect(p)


if __name__=='__main__':unittest.main(verbosity=2)


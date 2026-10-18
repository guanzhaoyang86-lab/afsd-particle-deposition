"""Q7 I/O adversarial fixtures. SCHEMA_FIXTURE_NOT_PHYSICAL, not solve evidence."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from inspect_process_plastic_friction import inspect, JOURNAL_FIELDS
from test_process_friction_cone_q6_readback import make_schema_fixture as q6_fixture, digest
from afsd_backend.pysph.output import checkpoint
from afsd_scaffold.vtk_legacy import read_point_cloud, write_point_cloud
from pysph.solver.utils import load

SOURCES=['scripts/probe_process_plastic_friction.py','configs/benchmarks/process_plastic_friction_q7.json',
    'scripts/probe_process_friction_cone.py','configs/benchmarks/process_friction_cone_q6.json',
    'scripts/probe_process_friction.py','configs/benchmarks/process_friction_q5.json',
    'scripts/probe_process_contact.py','configs/benchmarks/process_contact_q3.json',
    'scripts/probe_process_plastic.py','configs/benchmarks/process_plastic_q2.json',
    'scripts/probe_process_time.py','scripts/probe_engineering_material.py','scripts/probe_engineering_block.py',
    'scripts/probe_boundary_sbp.py','src/afsd_backend/pysph/material.py','src/afsd_backend/pysph/mechanics.py']


def read_history(p):
    with np.load(p/'history.npz',allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}


def make_schema_fixture(directory,initial_only=False):
    """Reuse serialized vector fixture, add synthetic independent P=Y V ep ledgers.

    Constitutive stresses/forces are NOT physically solved. No result from this
    fixture is allowed to certify J2/friction coupled mechanics.
    """
    p=q6_fixture(directory); h=read_history(p)
    cfg=json.loads((ROOT/'configs/benchmarks/process_plastic_friction_q7.json').read_text(encoding='utf-8'))
    cfg.update(increments_per_phase=4,basis_angle_deg=0.,evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
    ep=np.broadcast_to(np.array([0.,1e-4,2e-4])[:,None],(3,27)).copy()
    dP=90e6*h['reference_volume_m3'][None]*np.vstack((ep[0],np.diff(ep,axis=0)))
    dR=1.5*35e9*h['reference_volume_m3'][None]*np.vstack((ep[0],np.diff(ep,axis=0)))**2
    h.update(equivalent_plastic_strain=ep,plastic_work_J=np.cumsum(dP,axis=0),interval_plastic_work_J=dP,
        return_dissipation_diagnostic_J=np.cumsum(dR,axis=0),interval_return_diagnostic_J=dR,
        equivalent_strain_rate_per_s=np.broadcast_to(np.array([0.,.002,.008])[:,None],(3,27)).copy())
    r=json.loads((p/'report.json').read_text(encoding='utf-8'))
    for i,row in enumerate(r['records']):
        P=float(h['plastic_work_J'][i].sum()); R=float(h['return_dissipation_diagnostic_J'][i].sum())
        defect=row['stress_work_J']-row['elastic_diagnostic_J']-P
        row.update(plastic_work_J=P,return_dissipation_diagnostic_J=R,
            interval_plastic_work_J=float(dP[i].sum()),interval_return_diagnostic_J=float(dR[i].sum()),
            maximum_ep=float(ep[i].max()),volume_mean_ep=float(np.sum(h['reference_volume_m3']*ep[i])/h['reference_volume_m3'].sum()),
            maximum_strain_rate_per_s=float(h['equivalent_strain_rate_per_s'][i].max()),
            work_minus_elastic_minus_plastic_J=defect,material_work_warning=bool(abs(defect)>1e-6))
    if initial_only:
        for key in tuple(h):
            if h[key].ndim and h[key].shape[0]==3:h[key]=h[key][:1].copy()
        r['records']=r['records'][:1];r['events']=[];r['accepted_steps']=0;r['accepted_time_s']=0.
    r.update(classification=cfg['classification'],not_solved=cfg['not_solved'],
        material_work_warning_steps=[row['step'] for row in r['records'] if row['material_work_warning']])
    np.savez_compressed(p/'history.npz',**h)
    (p/'resolved_config.json').write_text(json.dumps(cfg),encoding='utf-8')
    saved=load(str(p/'checkpoint.npz'));pa=saved['arrays']['solid'];state=saved['solver_data']
    pa.x[:],pa.y[:],pa.z[:]=h['position_m'][-1].T;pa.ep[:]=h['equivalent_plastic_strain'][-1]
    for name,field in [('plastic_work','plastic_work_J'),('return_dissipation','return_dissipation_diagnostic_J')]:
        pa.add_property(name,data=h[field][-1])
    state.update(classification=cfg['classification'],step=r['accepted_steps'],physical_time_s=r['accepted_time_s'])
    checkpoint(p/'checkpoint.npz',pa,state)
    vtk=read_point_cloud(p/'final.vtk');scalars=vtk['scalars']
    for field in ('equivalent_plastic_strain','plastic_work_J','return_dissipation_diagnostic_J'):scalars[field]=h[field][-1].tolist()
    write_point_cloud(p/'final.vtk',h['position_m'][-1].tolist(),scalars,
        dict(reference_position_m=h['reference_position_m'].tolist(),displacement_m=h['displacement_m'][-1].tolist()),
        title='SCHEMA_FIXTURE_NOT_PHYSICAL')
    rejected=json.loads((p/'REJECTED_manifest.json').read_text(encoding='utf-8'))
    rejected.update(last_accepted_step=r['accepted_steps'],last_accepted_time_s=r['accepted_time_s'],
        attempted_time_s=.0125 if initial_only else .075)
    (p/'REJECTED_manifest.json').write_text(json.dumps(rejected),encoding='utf-8')
    provenance=dict(classification='Q7_BOUNDED_ATTEMPT_NOT_PRODUCTION',increments_per_phase=4,
        source_sha256={name:digest(ROOT/name) for name in SOURCES},
        configuration_sha256=digest(ROOT/'configs/benchmarks/process_plastic_friction_q7.json'))
    (p/'provenance.json').write_text(json.dumps(provenance),encoding='utf-8')
    for i,row in enumerate(r['records']):
        np.savez_compressed(p/f'accepted_{i:04d}.npz',**{name:h[field][i] for name,field in JOURNAL_FIELDS.items()},
            particle_id=h['particle_id'],tracer_type=h['tracer_type'],physical_time_s=h['physical_time_s'][i],step=i,
            record_json=json.dumps(row),classification='ACCEPTED_DIAGNOSTIC_NOT_RESTART',restart_allowed=False)
    r.update(source_sha256={name:digest(ROOT/name) for name in SOURCES},
        artifact_sha256={f.name:digest(f) for f in p.iterdir() if f.is_file() and f.name!='report.json'})
    (p/'report.json').write_text(json.dumps(r),encoding='utf-8')
    return p


class PlasticFrictionReadbackQ7(unittest.TestCase):
    def report(self,p):return json.loads((p/'report.json').read_text(encoding='utf-8'))
    def save_report(self,p,r):(p/'report.json').write_text(json.dumps(r),encoding='utf-8')
    def rehash(self,p,name):
        r=self.report(p);r['artifact_sha256'][name]=digest(p/name);self.save_report(p,r)
    def change(self,p,modify):
        h=read_history(p);modify(h);np.savez_compressed(p/'history.npz',**h);self.rehash(p,'history.npz')

    def test_01_schema_fixture_separates_material_and_friction_and_warns(self):
        with tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
            result=inspect(make_schema_fixture(d))
            self.assertEqual(result['status'],'PASS_PARTIAL_READBACK')
            self.assertEqual(result['evidence_kind'],'SCHEMA_FIXTURE_NOT_PHYSICAL')
            self.assertAlmostEqual(result['independent_plastic_work_J'],.018,delta=1e-15)
            self.assertAlmostEqual(result['independent_return_diagnostic_J'],.00105,delta=1e-15)
            self.assertAlmostEqual(sum(result['independent_friction_dissipation_J']),4.5e-6,delta=1e-18)
            self.assertEqual(result['material_work_warning_steps'],[1,2])
            self.assertFalse(result['material_energy_certified']);self.assertFalse(result['production_authorized'])

    def test_02_initial_only_blocked_is_honest_partial_not_activity(self):
        with tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
            result=inspect(make_schema_fixture(d,True))
            self.assertEqual(result['status'],'PASS_PARTIAL_READBACK');self.assertEqual(result['accepted_steps'],0)
            self.assertFalse(result['plastic_activity_observed']);self.assertFalse(result['friction_activity_observed'])

    def test_03_hash_provenance_and_rehashed_id_attacks(self):
        for attack,message in [('source','source changed'),('artifact','artifact changed'),('id','identity')]:
            with self.subTest(attack=attack),tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
                p=make_schema_fixture(d)
                if attack=='source':
                    r=self.report(p);r['source_sha256'][SOURCES[0]]='0'*64;self.save_report(p,r)
                elif attack=='artifact':(p/'final.vtk').write_bytes((p/'final.vtk').read_bytes()+b'\nTAMPER')
                else:self.change(p,lambda h:h['particle_id'].__setitem__(1,h['particle_id'][0]))
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_04_material_must_be_history_and_raw_nonnegative(self):
        for attack,message in [('shape','time history'),('negative','negative raw'),('decreasing','decreasing')]:
            with self.subTest(attack=attack),tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
                p=make_schema_fixture(d)
                def change(h):
                    if attack=='shape':h['equivalent_plastic_strain']=h['equivalent_plastic_strain'][-1]
                    elif attack=='negative':h['interval_plastic_work_J'][-1,0]=-1e-20
                    else:h['equivalent_plastic_strain'][-1,0]=0.
                self.change(p,change)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_05_plastic_dissipation_cannot_be_hidden_in_total(self):
        with tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
            p=make_schema_fixture(d)
            def change(h):h['plastic_work_J'][-1,0]+=1e-6;h['plastic_work_J'][-1,1]-=1e-6
            self.change(p,change)
            with self.assertRaisesRegex(ValueError,'Y V ep'):inspect(p)

    def test_06_return_increment_has_independent_quadratic_reference(self):
        with tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
            p=make_schema_fixture(d)
            def change(h):h['return_dissipation_diagnostic_J'][-1,0]+=1e-6;h['interval_return_diagnostic_J'][-1,0]+=1e-6
            self.change(p,change)
            with self.assertRaisesRegex(ValueError,'independent return'):inspect(p)

    def test_07_interval_ledger_and_energy_warning_cannot_be_relabelled(self):
        for attack,message in [('interval','plastic cumulative'),('warning','warning mismatch'),('defect','material work diagnostic')]:
            with self.subTest(attack=attack),tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
                p=make_schema_fixture(d)
                if attack=='interval':self.change(p,lambda h:h['interval_plastic_work_J'].__setitem__((-1,0),0.))
                else:
                    r=self.report(p)
                    r['records'][-1]['material_work_warning' if attack=='warning' else 'work_minus_elastic_minus_plastic_J']=False if attack=='warning' else 0.
                    self.save_report(p,r)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_08_rejected_candidate_cannot_be_accepted_and_partial_is_not_complete(self):
        for attack,message in [('metadata','BLOCKED/rejected'),('missing','presence'),('complete','terminal/requested')]:
            with self.subTest(attack=attack),tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
                p=make_schema_fixture(d)
                if attack=='complete':
                    r=self.report(p);r['status']='J2_CONE_PATH_COMPLETE';r['error']=None;self.save_report(p,r)
                else:
                    path=p/'REJECTED_manifest.json';rejected=json.loads(path.read_text(encoding='utf-8'))
                    if attack=='metadata':rejected['last_accepted_time_s']=rejected['attempted_time_s']
                    else:rejected['candidate_available']=True
                    path.write_text(json.dumps(rejected),encoding='utf-8');self.rehash(p,path.name)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_09_circle_direction_and_per_contact_work_are_independent_checks(self):
        for attack,message in [('circle','circular Coulomb'),('direction','vector slip opposition'),('work','per-contact friction work')]:
            with self.subTest(attack=attack),tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
                p=make_schema_fixture(d)
                def change(h):
                    if attack=='work':h['interval_friction_J'][-1,0]+=1e-6;h['interval_friction_J'][-1,1]-=1e-6
                    else:
                        value=[.08,.08] if attack=='circle' else [-.06,-.08]
                        h['tangent_force_global_yz_N'][-1,0]=value
                        node=np.flatnonzero(h['reference_position_m'][:,0]==.01)[0];h['contact_force_N'][-1,node,1:]=value
                self.change(p,change)
                if attack=='direction':
                    r=self.report(p);r['records'][-1]['wall_force_on_material_N'][1]-=.12;r['records'][-1]['wall_force_on_material_N'][2]-=.16;self.save_report(p,r)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_10_parameter_units_and_unsolved_claims_are_checked(self):
        for attack,message in [('unit','manufactured physical'),('thermal','unsolved field'),('production','cannot authorize')]:
            with self.subTest(attack=attack),tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
                p=make_schema_fixture(d)
                if attack=='unit':
                    cfg=json.loads((p/'resolved_config.json').read_text(encoding='utf-8'));cfg['parameters']['yield_stress']['unit']='MPa'
                    (p/'resolved_config.json').write_text(json.dumps(cfg),encoding='utf-8');self.rehash(p,'resolved_config.json')
                elif attack=='thermal':self.change(p,lambda h:h.update(temperature_K=np.full((3,27),293.15)))
                else:r=self.report(p);r['production_authorized']=True;self.save_report(p,r)
                with self.assertRaisesRegex(ValueError,message):inspect(p)

    def test_11_journal_preserves_last_accepted_not_rejected_state(self):
        for attack,message in [('state','journal/history'),('metadata','journal metadata'),('identity','journal ID')]:
            with self.subTest(attack=attack),tempfile.TemporaryDirectory(prefix='q7_readback_') as d:
                p=make_schema_fixture(d);path=p/'accepted_0002.npz'
                with np.load(path,allow_pickle=False) as saved:j={k:saved[k].copy() for k in saved.files}
                if attack=='state':j['ep'][0]+=1e-3
                elif attack=='metadata':j['classification']=np.array('REJECTED')
                else:j['particle_id'][0]+=100
                np.savez_compressed(path,**j);self.rehash(p,path.name)
                with self.assertRaisesRegex(ValueError,message):inspect(p)


if __name__=='__main__':unittest.main(verbosity=2)

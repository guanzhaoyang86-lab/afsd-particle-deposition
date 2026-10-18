"""TOP fixtures: literal references and complete scalar/batch implementation checks.

No actual process state is loaded and no NNPS/material/AFSD solve is performed.
"""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts'))
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact
from afsd_backend.pysph.quasistatic_contact import ContactGeometryError
from batched_contact import BatchedChainRuleContact
from contact_shell_finite_enclosure import audit_contact_shell_finite_enclosure
from contact_shell_audit import CERTIFIED,FOUND,UNRESOLVED


def pose(theta=0.,origin=(0.,0.,0.)):
    c,s=np.cos(theta),np.sin(theta)
    return dict(origin=np.asarray(origin,float),rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]))


def make(kind,scale=1.):return kind(2*scale,scale,scale,.1*scale,100.,200.,.3,base_z=None)


class VisibleTopContactTests(unittest.TestCase):
    def close(self,a,b,tol=1e-12):np.testing.assert_allclose(a,b,atol=tol,rtol=0.)

    def evaluate(self,contact,first,last=None,history=None,p0=None,p1=None):
        x0=np.asarray(first,float).reshape(-1,3)
        x1=x0.copy() if last is None else np.asarray(last,float).reshape(-1,3)
        return contact.evaluate(x0,x1,pose() if p0 is None else p0,pose() if p1 is None else p1,.1,history)

    def identities(self,r):
        self.assertTrue(np.all(np.abs(r['dg_identity_error_J'])<=r['dg_identity_tolerance_J']))
        self.assertTrue(np.all(np.abs(r['tangential_work_identity_error_J'])<=r['tangential_work_identity_tolerance_J']))
        normal_relative=r['normal_particle_work_J']-r['normal_wall_work_J']
        self.close(normal_relative+r['normal_potential1_J']-r['normal_potential0_J'],0.)
        total_relative=r['particle_work_J']-r['wall_work_J']
        storage=(r['normal_potential1_J']-r['normal_potential0_J']+
                 r['tangential_potential1_J']-r['tangential_potential0_J'])
        self.close(total_relative+storage+r['friction_dissipation_J']+
                   r['algorithmic_dissipation_J']+r['release_energy_J'],0.)
        self.assertTrue(np.all(r['friction_dissipation_J']>=0.))

    def equivalent(self,a,b):
        self.assertEqual(a.keys(),b.keys())
        for name,expected in b.items():
            with self.subTest(field=name):
                if isinstance(expected,np.ndarray):
                    if expected.dtype==bool:np.testing.assert_array_equal(a[name],expected)
                    else:np.testing.assert_allclose(a[name],expected,atol=3e-13,rtol=3e-13)
                else:self.assertEqual(a[name],expected)

    def test_01_top_planar_compression_literal_average_force_and_signed_work(self):
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind),[.8,.1,1.06],[.8,.1,1.04])
            self.assertEqual(set(np.flatnonzero(r['active'][0])),{1})
            self.close(r['normal_force_by_feature_N'][0,1],[0.,0.,5.])
            self.close(r['normal_potential0_J'].sum(),.08)
            self.close(r['normal_potential1_J'].sum(),.18)
            self.close(r['normal_particle_work_J'].sum(),-.1);self.identities(r)

    def test_02_top_creation_literal_secant_and_one_dimensional_sliding(self):
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind),[.8,.1,1.14],[.8,.16,1.08])
            self.close(r['normal_force_by_feature_N'][0,1],[0.,0.,1/3])
            self.close(r['tangential_force_by_feature_N'][0,1],[0.,-.1,0.])
            self.close(r['history'][0,1],[0.,.0005,0.])
            self.close(r['normal_particle_work_J'].sum(),-.02)
            self.close(r['friction_dissipation_J'].sum(),.00595)
            self.close(r['algorithmic_dissipation_J'].sum(),.000025)
            self.assertTrue(r['sliding'][0,1]);self.identities(r)

    def test_03_top_release_restores_normal_storage_and_separates_spring_release(self):
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            contact=make(kind);f=self.evaluate(contact,[.8,.1,1.14],[.8,.16,1.08])
            r=self.evaluate(contact,[.8,.16,1.08],[.8,.1,1.14],f['history'])
            self.close(r['normal_force_by_feature_N'],f['normal_force_by_feature_N'])
            self.close(r['normal_particle_work_J'].sum(),.02)
            self.close(r['history'],0.);self.close(r['friction_dissipation_J'].sum(),0.)
            self.close(r['release_energy_J'].sum(),.000025)
            self.assertFalse(r['active'].any());self.assertFalse(r['friction_active'].any())
            self.assertEqual(r['release_energy_policy'],'UNASSIGNED_ACCOUNT_NOT_AUTOMATICALLY_HEAT');self.identities(r)

    def test_04_top_stick_literal_storage_and_algorithmic_loss(self):
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind),[.8,.1,1.05],[.8,.101,1.05])
            self.close(r['force_N'],[[0.,-.2,5.]])
            self.close(r['tangential_potential1_J'].sum(),.0001)
            self.close(r['algorithmic_dissipation_J'].sum(),.0001)
            self.close(r['friction_dissipation_J'].sum(),0.)
            self.assertFalse(r['sliding'].any());self.identities(r)

    def test_05_top_sliding_literal_coulomb_force_and_heat(self):
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind),[.8,.1,1.05],[.8,.12,1.05])
            self.close(r['force_N'],[[0.,-1.5,5.]])
            self.close(r['history'][0,1],[0.,.0075,0.])
            self.close(r['tangential_potential1_J'].sum(),.005625)
            self.close(r['algorithmic_dissipation_J'].sum(),.005625)
            self.close(r['friction_dissipation_J'].sum(),.01875);self.identities(r)

    def test_06_top_bore_lip_literal_duplicate_potential_and_chain_work(self):
        first=np.array([.46,.1,1.03]);last=np.array([.48,.105,1.05])
        def U(q):return 100*(.1-np.hypot(q[0]-.5,q[2]-1.))**2
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind),first,last)
            self.assertEqual(set(np.flatnonzero(r['active'][0])),{1,3})
            self.close(r['normal_potential0_J'].sum(),U(first))
            self.close(r['normal_potential1_J'].sum(),U(last))
            self.close(r['normal_particle_work_J'].sum(),U(first)-U(last))
            self.close(r['normal_force_by_feature_N'][0,1],r['normal_force_by_feature_N'][0,3])
            self.assertIn('DUPLICATE_LIP_CONTRIBUTIONS_RETAINED',r['normal_model']);self.identities(r)

    def test_07_top_to_bore_entry_and_exit_literal_separate_patch_distances(self):
        first=np.array([.46,.1,1.02]);last=np.array([.46,.1,.98])
        def potentials(q):
            top=np.hypot(q[0]-.5,q[2]-1.)
            bore=np.hypot(q[0]-.5,max(q[2]-1.,0.))
            return 50*(.1-np.array([top,bore]))**2
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            contact=make(kind);f=self.evaluate(contact,first,last);r=self.evaluate(contact,last,first,f['history'])
            for value,a,b in ((f,first,last),(r,last,first)):
                self.close(value['normal_potential0_J'][0,[1,3]],potentials(a))
                self.close(value['normal_potential1_J'][0,[1,3]],potentials(b))
                self.close(value['normal_particle_work_J'].sum(),potentials(a).sum()-potentials(b).sum())
                self.assertEqual(set(np.flatnonzero(value['active'][0])),{1,3});self.identities(value)
            self.close(f['normal_force_by_feature_N'],r['normal_force_by_feature_N'])

    def test_08_top_common_rigid_motion_preserves_history_and_has_no_relative_work(self):
        first=np.array([[.8,.1,1.05]]);p1=pose(.25,[.03,-.04,.02])
        last=first@p1['rotation'].T+p1['origin'];history=np.zeros((1,8,3));history[0,1,1]=.003
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind),first,last,history,p1=p1)
            self.close(r['history'],history)
            self.close(r['normal_force_by_feature_N'][0,1],np.array([0.,0.,5.])@p1['rotation'].T)
            self.close(r['friction_dissipation_J'].sum(),0.)
            self.close(r['particle_work_J'],r['wall_work_J']);self.identities(r)

    def test_09_moving_top_lip_signed_wall_work_and_general_frame_covariance(self):
        first=np.array([[.46,.1,1.03]]);body_last=np.array([[.48,.105,1.05]])
        p0=pose();p1=pose(.05,[.003,-.004,.002]);last=body_last@p1['rotation'].T+p1['origin']
        Q=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]]);shift=np.array([.2,-.1,.3])
        def moved(p):return dict(origin=p['origin']@Q.T+shift,rotation=Q@p['rotation'])
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            contact=make(kind);plain=self.evaluate(contact,first,last,p0=p0,p1=p1)
            other=self.evaluate(contact,first@Q.T+shift,last@Q.T+shift,p0=moved(p0),p1=moved(p1))
            self.close(other['force_N'],plain['force_N']@Q.T,2e-12)
            self.close(other['history'],plain['history'])
            for name in ('particle_work_J','wall_work_J','normal_potential0_J','normal_potential1_J',
                         'friction_dissipation_J','algorithmic_dissipation_J','release_energy_J'):
                self.close(other[name],plain[name])
            self.identities(plain);self.identities(other)

    def test_10_rotating_translating_top_scalar_batch_complete_fields_and_input_purity(self):
        first=np.array([[.8,.1,1.06],[.8,.1,1.14],[.46,.1,1.03],[.46,.1,1.02]])
        end=np.array([[.8,.12,1.04],[.8,.16,1.08],[.48,.105,1.05],[.46,.1,.98]])
        p0=pose();p1=pose(.02,[.001,-.002,.0003]);last=end@p1['rotation'].T+p1['origin']
        history=np.zeros((4,8,3));history[0,1,1]=.002;history[2,1,1]=.001;history[2,3,1]=.001
        inputs=(first,last,history,p0['origin'],p0['rotation'],p1['origin'],p1['rotation'])
        before=[v.copy() for v in inputs]
        scalar=make(ChainRuleDiscreteContact);batch=make(BatchedChainRuleContact)
        for points,p in ((first,p0),(last,p1)):self.equivalent(batch.geometry(points,p),scalar.geometry(points,p))
        expected=self.evaluate(scalar,first,last,history,p0,p1);actual=self.evaluate(batch,first,last,history,p0,p1)
        self.equivalent(actual,expected);self.identities(actual)
        for value,original in zip(inputs,before):np.testing.assert_array_equal(value,original)
        self.assertFalse(np.shares_memory(actual['history'],history))

    def test_11_top_finite_clear_path_certificate_retains_opening(self):
        contact=make(BatchedChainRuleContact);contact.particle_radius=.25
        first=np.array([[.35,.1,1.21]]);last=np.array([[.35,.12,1.21]])
        r=audit_contact_shell_finite_enclosure(contact,first,last,pose(),pose(),maximum_queries=2)
        self.assertEqual(r['status'],CERTIFIED)
        self.assertTrue(any(row['feature']=='TOP' for row in r['finite_patch_certificates']))
        self.assertEqual(r['additional_finite_patch_queries'],0)

    def test_12_hidden_rotating_top_contact_not_certified_from_inactive_endpoints(self):
        x=np.array([[.41,0.,1.04]]);p0,p1=pose(-np.pi/4),pose(np.pi/4)
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            contact=make(kind)
            self.assertFalse(contact.geometry(x,p0)['active'].any());self.assertFalse(contact.geometry(x,p1)['active'].any())
            self.close(contact.geometry(x,pose())['gaps_m'][0,1],np.hypot(.09,.04))
            limited=audit_contact_shell_finite_enclosure(contact,x,x,p0,p1,maximum_queries=2,maximum_path_evaluations=8)
            self.assertEqual(limited['status'],UNRESOLVED)
            full=audit_contact_shell_finite_enclosure(contact,x,x,p0,p1)
            self.assertEqual(full['status'],FOUND)
            self.assertTrue(any(row['feature']=='TOP' for row in full['hidden_contacts']))

    def test_13_open_bore_top_entry_has_zero_contact(self):
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind),[0.,0.,1.2],[0.,0.,.8])
            self.assertTrue(r['sweep']['certified']);self.assertFalse(r['active'].any())
            self.close(r['force_N'],0.);self.close(r['normal_potential1_J'],0.);self.identities(r)

    def test_14_top_medial_axis_stays_rejected(self):
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            with self.assertRaisesRegex(ContactGeometryError,'multiple closest'):
                self.evaluate(make(kind),[.47,.47,1.04])

    def test_15_top_millimetre_fixture_matches_literal_dimensional_scaling(self):
        scale=.001
        for kind in (ChainRuleDiscreteContact,BatchedChainRuleContact):
            r=self.evaluate(make(kind,scale),scale*np.array([.8,.1,1.06]),scale*np.array([.8,.1,1.04]))
            self.close(r['normal_force_by_feature_N'][0,1],[0.,0.,5*scale])
            self.close(r['normal_potential1_J'].sum(),.18*scale**2)
            self.close(r['normal_particle_work_J'].sum(),-.1*scale**2);self.identities(r)


if __name__=='__main__':unittest.main(verbosity=2)

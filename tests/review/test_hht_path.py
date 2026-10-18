"""Literal polynomial and event references independent of tested path functions."""
import sys
from pathlib import Path
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import hht_path as hp


def audit(v0,v1,a0,a1,dt,f0=(-1.,0.,0.),f1=(-1.,0.,0.),wall=(0.,0.,0.)):
    return hp.contact_audit([0,0,1],f0,f1,wall,wall,v0,v1,a0,a1,dt,1.,1.)


class PathTests(unittest.TestCase):
    def test_constant_acceleration_exact(self):
        x,v=hp.endpoint([2.,0,0],[1.,0,0],[-3.,0,0],[-3.,0,0],1/3)
        np.testing.assert_allclose(x,[13/6,0,0],atol=1e-12,rtol=0)
        np.testing.assert_allclose(v,0,atol=1e-12,rtol=0)

    def test_dense_path_integrates_newmark_displacement(self):
        x0=np.array([2.,-1.,.3]);v0=np.array([.1,.2,-.5]);a0=np.array([4.,-3.,2.]);a1=np.array([-1.,2.,3.])
        x1,v1=hp.endpoint(x0,v0,a0,a1,.2)
        vals=[hp.sample(x0,v0,a0,a1,.2,s)[1] for s in [0,.5,1]]
        np.testing.assert_allclose(x1-x0,.2*(vals[0]+4*vals[1]+vals[2])/6,atol=1e-12,rtol=0)
        np.testing.assert_allclose(hp.sample(x0,v0,a0,a1,.2,1)[0],x1,atol=1e-12,rtol=0)

    def test_linear_force_work_hand_polynomial(self):
        # a0=0,a1=2,h=1 -> dv=1.2, path vx=1+1.23s-.03s².
        # F=2+3s: integral =2+5.46/2+3.63/3-.09/4.
        w=hp.linear_work([2.,0,0],[5.,0,0],[1.,0,0],[2.2,0,0],[0.,0,0],[2.,0,0],1.)
        self.assertAlmostEqual(float(w['work_J']),5.9175,delta=1e-12)
        np.testing.assert_allclose(w['impulse_Ns'],[3.5,0,0],atol=1e-12,rtol=0)

    def test_hidden_double_reversal_even_when_integrated_heat_zero(self):
        r=audit([1.,0,0],[1.,0,0],[240.,0,0],[-160.,0,0],1.)
        self.assertAlmostEqual(r['friction_heat_raw_J'],0.,delta=1e-12)
        self.assertAlmostEqual(r['maximum_friction_power_W'],.5,delta=1e-12)
        self.assertFalse(r['whole_path_passive']);self.assertTrue(r['requires_event_split'])
        np.testing.assert_allclose(r['interior_zero_fractions'],[(3-np.sqrt(3))/6,(3+np.sqrt(3))/6],atol=1e-10,rtol=0)

    def test_event_split_reversal_literal_heat(self):
        first=audit([1.,0,0],[0.,0,0],[-3.,0,0],[-3.,0,0],1/3)
        last=audit([0.,0,0],[-2/3,0,0],[-1.,0,0],[-1.,0,0],2/3,f0=(1.,0,0),f1=(1.,0,0))
        self.assertTrue(first['whole_path_passive']);self.assertTrue(last['whole_path_passive'])
        self.assertAlmostEqual(first['friction_heat_raw_J']+last['friction_heat_raw_J'],7/18,delta=1e-12)
        self.assertFalse(first['requires_event_split']);self.assertFalse(last['requires_event_split'])

    def test_relative_wall_boost_and_normal_velocity(self):
        original=audit([1.,0,.7],[0.,0,.7],[-3.,0,0],[-3.,0,0],1/3)
        boosted=audit([3.,-1.,1.4],[2.,-1.,1.4],[-3.,0,0],[-3.,0,0],1/3,wall=(2.,-1.,.7))
        self.assertAlmostEqual(original['friction_heat_raw_J'],boosted['friction_heat_raw_J'],delta=1e-12)
        self.assertTrue(boosted['whole_path_passive'])

    def test_signed_algorithmic_kinetic_identity(self):
        m=np.array([2.]);f0=np.array([[3.,1.,0]]);f1=np.array([[-2.,4.,1.]])
        a0=np.array([[1.,.7,.2]]);a1=(.9*f1+.1*f0)/m[:,None]
        v0=np.array([[.3,-.2,.1]]);_,v1=hp.endpoint(np.zeros((1,3)),v0,a0,a1,.1)
        r=hp.kinetic_identity(m,v0,v1,a0,a1,f0,f1,.1)
        np.testing.assert_allclose(r['kinetic_identity_residual_J'],0,atol=1e-12,rtol=0)
        np.testing.assert_allclose(r['momentum_identity_residual_Ns'],0,atol=1e-12,rtol=0)

    def test_non_tangent_or_excess_traction_rejected(self):
        for f in ((-2.,0,0),(-1.,0,.1)):
            with self.assertRaisesRegex(ValueError,'tangent cone'):
                audit([1.,0,0],[1.,0,0],[0.,0,0],[0.,0,0],1.,f0=f)

    def test_tangent_touch_zero_is_not_lost_by_cubic_root_conditioning(self):
        # Independently v(s)=(s-.5)^2; power never assists, but zero exists.
        r=audit([.25,0,0],[.25,0,0],[40.,0,0],[-80/3,0,0],1.)
        self.assertTrue(r['whole_path_passive']);self.assertTrue(r['requires_event_split'])
        self.assertLess(min(abs(s-.5) for s in r['interior_zero_fractions']),1e-10)

    def test_invalid_step_and_nonfinite_physics_rejected(self):
        for h in (0.,-1.,np.nan,np.inf):
            with self.assertRaisesRegex(ValueError,'dt'):
                audit([.25,0,0],[.25,0,0],[0.,0,0],[0.,0,0],h)
        with self.assertRaisesRegex(ValueError,'finite'):
            hp.linear_work([np.nan,0,0],[1.,0,0],[1.,0,0],[1.,0,0],[0.,0,0],[0.,0,0],1.)
        with self.assertRaisesRegex(ValueError,'mass'):
            hp.kinetic_identity(np.array([np.nan]),np.zeros((1,3)),np.zeros((1,3)),
                np.zeros((1,3)),np.zeros((1,3)),np.zeros((1,3)),np.zeros((1,3)),1.)


if __name__=='__main__':unittest.main(verbosity=2)

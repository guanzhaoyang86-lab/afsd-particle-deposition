"""Independent literal/algebra checks for the slow RK4 oracle and event rules."""
from pathlib import Path
import sys
import time
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import reference_hht_slow as ref
from hht_slow_inputs import inputs


class SlowReference(unittest.TestCase):
    def test_literal_boundary_and_no_old_signed_speed_floor(self):
        d=inputs(); p=d['parameters']; x=d['position_m']; v=d['velocity_m_s'].copy(); T=d['temperature_K']
        b=ref.boundary_fields(0,x,v,T,d)
        np.testing.assert_allclose(b['face_force'],[[126,0,25.2],[-126,0,25.2]],rtol=0,atol=5e-11)
        np.testing.assert_allclose(b['friction'].sum(axis=1),[4.536,4.536],rtol=0,atol=1e-11)
        v[:,2]=.19 # genuine .01 m/s slip: not the old .05 fixture domain
        small=ref.boundary_fields(0,x,v,T,d)
        self.assertAlmostEqual(small['minimum_speed_m_s'],.01,delta=1e-16)
        v[:,1]=.03;v[:,2]=.21 # positive z AND transverse slip is not a stop
        turning=ref.boundary_fields(0,x,v,T,d)
        self.assertGreater(float(turning['friction'].sum()),0)

    def test_exact_wall_history_no_microsecond_control_jump(self):
        p=inputs()['parameters']; plane0,w0=ref.wall_program(0,p)
        plane,w=ref.wall_program(.05,p) # quarter of a .2 s period
        np.testing.assert_allclose(plane-plane0,[1e-7,-1e-7],atol=1e-18,rtol=0)
        np.testing.assert_allclose(w[:,0],[np.pi*1e-6,-np.pi*1e-6],atol=1e-20,rtol=0)
        for t in (0,1e-6,2e-6,20e-6):
            np.testing.assert_array_equal(ref.wall_program(t,p)[1][:,2],[.2,.2])

    def test_independent_dense_constant_acceleration_and_rk4_endpoint(self):
        h=.4;y=np.array([2.,3.]);k1=np.array([3.,5.]);k2=np.array([4.,5.]);k3=k2.copy();k4=np.array([5.,5.])
        c=ref.dense_coefficients(y,h,[k1,k2,k3,k4])
        for s in (0.,.125,.5,1.):
            exact=np.array([2+3*h*s+2.5*(h*s)**2,3+5*h*s])
            np.testing.assert_allclose(ref.dense_value(c,s),exact,rtol=0,atol=2e-15)
        # RK4 polynomial y'=y, h=.2. Literal independently expanded polynomial.
        h=.2; ks=[np.ones(1),np.array([1+h/2]),np.array([1+h/2+h*h/4]),np.array([1+h+h*h/2+h**3/4])]
        c=ref.dense_coefficients(np.ones(1),h,ks)
        self.assertAlmostEqual(float(ref.dense_value(c,1)[0]),1+h+h*h/2+h**3/6+h**4/24,delta=3e-16)

    def test_z_turn_is_not_zero_vector_and_true_zero_is_detected(self):
        d=inputs(); y=ref.packed_initial(d); c=np.zeros((4,len(y)));c[0]=y
        vc=c[:,3*ref.N:6*ref.N].reshape(4,ref.N,3)
        vc[0,:,1]=.03;vc[0,:,2]=.1;vc[1,:,2]=.2
        curved=ref.tangent_events(c,d['reference_position_m'],d)
        self.assertFalse(curved['events']);self.assertEqual(len(curved['component_turns']),18)
        for row in curved['component_turns']:
            self.assertAlmostEqual(row['fraction'],.5,delta=1e-15)
            self.assertAlmostEqual(row['full_speed_m_s'],.03,delta=1e-15)
            self.assertFalse(row['stop_certified'])
        vc[:,:,1]=0
        stopped=ref.tangent_events(c,d['reference_position_m'],d)
        self.assertEqual({row['particle_id'] for row in stopped['events']},
            {int(d['particle_id'][j]) for ids in d['face_nodes'] for j in ids})
        self.assertTrue(all(abs(row['fraction']-.5)<1e-10 for row in stopped['events']))

    def test_touching_zero_without_sign_change(self):
        d=inputs(); y=ref.packed_initial(d); c=np.zeros((4,len(y)));c[0]=y
        vc=c[:,3*ref.N:6*ref.N].reshape(4,ref.N,3)
        vc[0,:,2]=.45;vc[1,:,2]=-1;vc[2,:,2]=1 # ut=(s-.5)^2
        audited=ref.tangent_events(c,d['reference_position_m'],d)
        self.assertTrue(audited['events'])
        self.assertTrue(any(abs(row['fraction']-.5)<1e-8 for row in audited['events']))

    def test_physical_opening_and_singular_zero_rejected_without_clipping(self):
        d=inputs(); y=ref.packed_initial(d); before=y.copy()
        x,v,S,rho,T,ell=ref.unpack(y);x[0,0]=3e-6
        with self.assertRaises(ref.DomainBoundary) as ctx:ref.fields(0,y,d)
        self.assertEqual(ctx.exception.kind,'OPENING')
        y=before.copy(); ref.unpack(y)[1][0,2]=.2
        with self.assertRaises(ref.DomainBoundary) as ctx:ref.fields(0,y,d)
        self.assertEqual(ctx.exception.kind,'ZERO_SLIP_SINGULAR')
        np.testing.assert_array_equal(d['position_m'][0],[0,0,0])

    def test_rhs_literal_heat_and_power_preserve_nonzero_startup(self):
        d=inputs(); y=ref.packed_initial(d); before=y.copy();dy=ref.rhs(0,y,d)
        np.testing.assert_array_equal(y,before)
        xdot,acc,sdot,rhodot,tdot,ledger=ref.unpack(dy)
        np.testing.assert_array_equal(xdot,d['velocity_m_s'])
        self.assertGreater(float(np.max(np.abs(acc))),1e4)
        self.assertAlmostEqual(float(np.sum(d['mass_kg']*900*tdot)),6.3+.5*9.072,delta=1e-10)
        b=ref.boundary_fields(0,d['position_m'],d['velocity_m_s'],d['temperature_K'],d)
        np.testing.assert_allclose(b['normal_material_work']+b['tangent_material_work']-b['wall_work']+b['friction'].sum(axis=1),0,atol=1e-13)

    def test_rejected_candidate_is_not_published(self):
        d=inputs(); initial=ref.packed_initial(d)
        def fail(t,y,h,data,deadline):
            candidate=y.copy();candidate[0]=9.
            raise ref.DomainBoundary('OPENING',t+h,candidate,face=0,particle_id=1)
        with patch.object(ref,'rk4_candidate',side_effect=fail):
            result=ref.run_reference(d,ref.REFERENCE_STEPS_S[0],time.perf_counter()+1)
        self.assertEqual(result['accepted_steps'],0)
        np.testing.assert_array_equal(result['states'][0],initial)
        self.assertEqual(len(result['rejected']),13)
        self.assertEqual(result['terminal']['actual_halvings'],12)
        self.assertEqual(result['terminal']['candidate_attempts'],13)
        self.assertEqual(result['terminal']['status'],'EVENT_APPROACH_BLOCKED')
        self.assertTrue(all(not r['accepted'] and not r['restart_allowed'] for r in result['rejected']))

    def test_reestablished_positive_right_endpoint_is_not_a_fake_event(self):
        d=inputs(); calls=0
        def first_stage_overshoot_then_positive(t,y,h,data,deadline):
            nonlocal calls
            calls+=1
            if calls==1:
                raise ref.DomainBoundary('OPENING',t+h,y,face=0,particle_id=1)
            return y.copy(),dict(component_turns=[])
        with patch.object(ref,'rk4_candidate',side_effect=first_stage_overshoot_then_positive):
            result=ref.run_reference(d,ref.REFERENCE_STEPS_S[0],time.perf_counter()+1)
        term=result['terminal']
        self.assertEqual(term['status'],'EVENT_LOCALIZATION_INCONSISTENT')
        self.assertIsNone(term['interval_s'])
        self.assertEqual(term['actual_halvings'],1)
        self.assertEqual(term['candidate_attempts'],3)
        self.assertEqual(result['accepted_steps'],1)
        self.assertEqual(result['times_s'][-1],ref.REFERENCE_STEPS_S[0]/2)
        self.assertEqual(result['rejected'][-1]['kind'],'EVENT_LOCALIZATION_INCONSISTENT')
        self.assertFalse(ref.event_agreement(result,result)['passed'])


if __name__=='__main__':unittest.main(verbosity=2)

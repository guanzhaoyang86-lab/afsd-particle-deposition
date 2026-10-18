"""Literal analytic event/reversal references, independent of tested backend."""
from pathlib import Path
import math
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import reference_hht_event as r


class IndependentEventReference(unittest.TestCase):
    def test_scalar_reversal_literal_continuous_solution(self):
        got=r.scalar_reference()
        self.assertAlmostEqual(got['event_time_s'],1/3,delta=1e-12)
        self.assertAlmostEqual(got['event_state'][0],1/6,delta=1e-12)
        self.assertLessEqual(abs(got['event_state'][1]),got['root']['full_zero_bound_m_s'])
        np.testing.assert_allclose(got['final'],[-1/18,-2/3,7/18,1/9],rtol=0,atol=1e-12)
        self.assertAlmostEqual(got['delta_kinetic_J'],-5/18,delta=1e-12)
        self.assertFalse(got['static']['feasible'][0]);self.assertAlmostEqual(got['static']['friction_force_N'][0,0],1.)
        self.assertFalse(got['root']['velocity_projection_used'])
        self.assertAlmostEqual(float(got['final'][3]-got['final'][2]),got['delta_kinetic_J'],delta=1e-12)
        np.testing.assert_array_equal(got['physical_jump'],0)

    def test_feasible_stick_reaction_not_automatic_reverse(self):
        got=r.scalar_reference(other=-.5)
        self.assertAlmostEqual(got['event_time_s'],2/3,delta=1e-12)
        np.testing.assert_allclose(got['final'],[1/3,0.,1/3,-1/6],rtol=0,atol=1e-12)
        self.assertTrue(got['static']['feasible'][0])
        self.assertAlmostEqual(got['static']['friction_force_N'][0,0],.5,delta=1e-15)
        np.testing.assert_array_equal(got['static']['resultant_tangent_force_N'],0)
        self.assertAlmostEqual(got['delta_kinetic_J'],-.5,delta=1e-12)
        np.testing.assert_array_equal(got['physical_jump'],0)

    def test_vector_component_turn_analytic_constant_other_force(self):
        # Polar ODE integration gives k=sqrt(sqrt2-1), t(ux=0)=(k+2-sqrt2)/3,
        # uy=k. These formulas do not call any contact evaluator for answers.
        root2=math.sqrt(2);k=math.sqrt(root2-1);target=(k+2-root2)/3
        def fun(t,u):return np.array([-2.,0.])-u/np.linalg.norm(u)
        u=np.array([1.,1.]);h=target/1024
        for step in range(1024):u,_=r.rk4_step(fun,step*h,u,h)
        np.testing.assert_allclose(u,[0.,k],rtol=0,atol=1e-10)
        self.assertGreater(np.linalg.norm(u),.6)
        self.assertGreater(np.linalg.norm(u),r.zero_bound(u))

    def test_scalar_root_cannot_certify_nonzero_transverse_slip(self):
        y=np.array([1.,.1]);before=y.copy()
        with self.assertRaises(r.ReferenceStop) as ctx:
            r.locate_zero(lambda h:y+np.array([-2*h,0.]),0.,y,1.,lambda z:z[None,:],np.array([[1.,0.]]))
        self.assertEqual(ctx.exception.kind,'SCALAR_ROOT_NOT_FULL_VECTOR_ZERO')
        np.testing.assert_array_equal(y,before)

    def test_joint_symmetric_root_no_velocity_assignment(self):
        y=np.array([1.,1.]);before=y.copy()
        got=r.locate_zero(lambda h:y-3*h,0.,y,1.,lambda z:z[:,None],np.ones((2,1)))
        self.assertAlmostEqual(got['time_s'],1/3,delta=1e-12)
        np.testing.assert_array_equal(got['state'],y-3*got['time_s'])
        np.testing.assert_array_equal(y,before)
        self.assertLessEqual(got['evaluations'],32)
        self.assertLessEqual(got['interval_s'][1]-got['interval_s'][0],got['time_width_bound_s'])
        self.assertTrue(all(not row['accepted'] and not row['restart_allowed'] for row in got['candidates']))

    def test_mean_zero_does_not_hide_individual_nonzero_slip(self):
        # Initial real coarse-reference failure: opposite root residuals can
        # have exactly zero mean but neither contact meets full-zero bound.
        y=np.array([1.-5e-13,1.+5e-13]);before=y.copy()
        with self.assertRaises(r.ReferenceStop) as ctx:
            r.locate_zero(lambda h:y-h,0.,y,2.,lambda z:z[:,None],np.ones((2,1)))
        self.assertEqual(ctx.exception.kind,'SCALAR_ROOT_NOT_FULL_VECTOR_ZERO')
        self.assertGreater(np.max(np.linalg.norm(ctx.exception.details['tangential_velocity_m_s'],axis=1)),r.zero_bound(y))
        np.testing.assert_array_equal(y,before)

    def test_measured_event_order_not_particle_identifier_order(self):
        first=dict(time_s=1.,interval_s=[.999,1.001])
        later=dict(time_s=2.,interval_s=[1.999,2.001])
        self.assertEqual(r.select_first_event([(4,later),(22,first)])[0],22)
        with self.assertRaises(r.ReferenceStop) as ctx:
            r.select_first_event([(4,first),(22,dict(time_s=1.0005,interval_s=[1.,1.001]))])
        self.assertEqual(ctx.exception.kind,'INDIVIDUAL_EVENT_ORDER_NOT_CERTIFIED')

    def test_positive_near_root_slip_has_own_time_no_sync(self):
        y=np.array([5e-14]);before=y.copy();start=6e-6
        got=r.locate_zero(lambda h:y-8e4*h,start,y,1e-11,lambda z:z[:,None],
            np.ones((1,1)),allow_near_incoming=True)
        self.assertGreater(got['time_s'],start)
        self.assertAlmostEqual(got['time_s']-start,6.25e-19,delta=1e-21)
        np.testing.assert_array_equal(y,before)
        self.assertFalse(got['velocity_projection_used'])

    def test_two_force_driven_events_keep_actual_intervening_force(self):
        y=np.array([1.,1.0001]);before=y.copy();probes=[]
        for j in range(2):
            probes.append((j,r.locate_zero(lambda h:y-3*h,0.,y,1.,lambda z,j=j:z[[j],None],np.ones((1,1)))))
        j,first=r.select_first_event(probes);self.assertEqual(j,0)
        state=first['state'].copy();time=first['time_s'];force=np.array([-1.,-3.])
        second=r.locate_zero(lambda h:state+force*h,time,state,.001,
            lambda z:z[[1],None],np.ones((1,1)),allow_near_incoming=True)
        self.assertAlmostEqual(first['time_s'],1/3,delta=1e-12)
        self.assertAlmostEqual(second['time_s'],1.0001/3,delta=1e-12)
        self.assertAlmostEqual(second['state'][0],-1e-4/3,delta=1e-12)
        np.testing.assert_array_equal(y,before)

    def test_nonbracketing_and_invalid_static_data_reject(self):
        with self.assertRaises(r.ReferenceStop) as ctx:
            r.locate_zero(lambda h:np.array([1.+h]),0.,np.array([1.]),1.,lambda z:z[:,None],np.ones((1,1)))
        self.assertEqual(ctx.exception.kind,'EVENT_NOT_BRACKETED')
        for force,cap in (([[float('nan'),0]],[1.]),([[1.,0]],[-1.])):
            with self.assertRaises(ValueError):r.static_disk(force,cap)

    def test_right_limit_only_at_located_epoch_and_raw_signed_power(self):
        data=r.inputs();x=data['position_m'];v=data['velocity_m_s'].copy();v[3]=[0.,0.,.2-1e-14]
        y=np.r_[x.ravel(),v.ravel(),data['deviatoric_stress_Pa'].ravel(),data['density_kg_m3'],data['temperature_K'],np.zeros(r.LEDGER_SIZE)]
        before=y.copy();mode=dict(kind='RIGHT_EVENT_LIMIT',time_s=0.,time_bound_s=5e-19,directions={3:np.array([0.,0.,1.])})
        _,b,_=r.fields(0.,y,data,mode)
        self.assertEqual(len(b['right_limit_records']),1)
        self.assertLess(b['friction'][0,3],0) # retained raw root-roundoff residual, no abs/max
        self.assertLessEqual(abs(b['friction'][0,3]),b['right_limit_records'][0]['root_roundoff_power_bound_W'])
        np.testing.assert_array_equal(y,before)
        with self.assertRaises(r.ReferenceStop) as ctx:r.fields(1e-9,y,data,mode)
        self.assertEqual(ctx.exception.kind,'NEXT_FULL_VECTOR_ZERO')

    def test_positive_small_outgoing_uses_actual_vector_not_direction_override(self):
        data=r.inputs();v=data['velocity_m_s'].copy();v[3]=[0.,1e-14,.2+2e-14]
        y=np.r_[data['position_m'].ravel(),v.ravel(),data['deviatoric_stress_Pa'].ravel(),data['density_kg_m3'],data['temperature_K'],np.zeros(r.LEDGER_SIZE)]
        before=y.copy();mode=dict(kind='SEQUENTIAL_RIGHT',directions={3:np.array([0.,0.,1.])},
            event_epochs={3:0.},time_bound_s=5e-19,allow_small_outgoing=True)
        _,b,_=r.fields(1e-17,y,data,mode);u=v[3]-[0.,0.,.2]
        cap=data['parameters']['friction']*data['parameters']['normal_stiffness']*b['penetration'][0,3]
        np.testing.assert_allclose(b['force'][3,1:],-cap*u[1:]/np.linalg.norm(u),rtol=1e-15,atol=1e-15)
        self.assertAlmostEqual(b['friction'][0,3],cap*np.linalg.norm(u),delta=1e-25)
        self.assertEqual(b['right_limit_records'],[])
        np.testing.assert_array_equal(y,before)

    def test_endpoint_alias_requires_positive_one_sided_polynomial_certificate(self):
        poly=np.array([[0.,1e-18,1e-14],[0.,-1e-8,1e-4],[0.,0.,0.],[0.,0.,0.]])
        before=poly.copy();start=6e-6;epoch=start-5e-19;bound=512*np.finfo(float).eps
        got=r.endpoint_alias_certificate(poly,4e-9,start,epoch,np.array([0.,0.,1.]),1e-10,bound)
        self.assertIsNotNone(got)
        self.assertAlmostEqual(got['projected_acceleration_min_m_s2'],25000.,delta=1e-10)
        self.assertAlmostEqual(got['projected_velocity_min_m_s'],1e-14,delta=1e-28)
        self.assertGreater(got['projected_velocity_at_candidate_m_s'],0)
        self.assertLessEqual(got['elapsed_since_event_s'],got['derivative_derived_roundoff_window_s'])
        np.testing.assert_array_equal(poly,before)
        self.assertIsNone(r.endpoint_alias_certificate(poly,4e-9,start,start-1e-12,np.array([0.,0.,1.]),1e-10,bound))
        crossing=poly.copy();crossing[0,2]=-1e-14
        self.assertIsNone(r.endpoint_alias_certificate(crossing,4e-9,start,epoch,np.array([0.,0.,1.]),1e-10,bound))
        returning=poly.copy();returning[0,2]=2e-14;returning[1,2]=-1e-4
        self.assertIsNone(r.endpoint_alias_certificate(returning,4e-9,start,epoch,np.array([0.,0.,1.]),1e-10,bound))
        large=poly.copy();large[0,1]=.1
        self.assertIsNone(r.endpoint_alias_certificate(large,4e-9,start,epoch,np.array([0.,0.,1.]),1e-10,bound))

    def test_endpoint_alias_checks_derivative_over_whole_interval(self):
        # End derivative positive, start derivative negative: cannot certify
        # monotone departure by inspecting only the candidate endpoint.
        poly=np.zeros((4,3));poly[:,2]=[1e-14,-1e-14,2e-14,0.]
        self.assertIsNone(r.endpoint_alias_certificate(poly,1.,0.,0.,np.array([0.,0.,1.]),.5,512*np.finfo(float).eps))


if __name__=='__main__':unittest.main(verbosity=2)

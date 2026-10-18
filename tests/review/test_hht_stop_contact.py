"""Independent one-sided/PRIVATE prescribed-path contact integral references.

This is NOT a force-driven particle trajectory. The paths isolate quadrature,
one-sided endpoint limits and rejection contracts with literal analytic SI
references, registered in PROCESS_STOP_BENCHMARK_CONTRACT before execution.
"""
import copy
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import hht_stop_contact as tested
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.updated_trial import particle_hash


def fixture(u0=(1.,0.,0.),a0=(-3.,0.,0.),a1=None,dt=1./3.,
            wall0=(0.,0.,0.),wall_change=(0.,0.,0.),kind='left_zero',
            incoming=(1.,0.,0.),outgoing=None):
    wall0=np.asarray(wall0,float);change=np.asarray(wall_change,float)
    a0=np.asarray(a0,float)+change/dt
    a1=a0.copy() if a1 is None else np.asarray(a1,float)+change/dt
    p=dict(dx=.005,h_ratio=1.3,rho0=2700.,width=.01,marker_low=-.1,marker_high=.1,
        bar_length=.01,actuator_depth=.001,temperature=400.,traverse=0.,feed=0.,
        shear=35e9,bulk=70e9/3,cp=900.,conductivity=167.,room=293.15,melt=855.15,
        friction=1./3.,plastic_heat_fraction=1.,heat_partition=.5,
        contact_coulomb_relative_tolerance=.01,hht_weight=.9)
    pa=make_array(np.array([[.03,-.02,0.]]),np.zeros(1,int),p)
    pa.u[:],pa.v[:],pa.w[:]=(wall0+u0).reshape(3,1)
    t0=1.25
    def boundary(private,at):
        s=(at-t0)/dt;wall=wall0+s*change
        v=np.array([private.u[0],private.v[0],private.w[0]])
        slip=(v-wall)*[1.,1.,0.];speed=np.linalg.norm(slip)
        ft=np.zeros(3) if speed==0 else -slip/speed
        # The callback explicitly declares the selected one-sided or PRIVATE
        # law; the helper must validate it, not silently replace an opaque
        # ordinary-Coulomb callback to obtain a desired passing result.
        if kind=='left_zero' and at==t0+dt:ft=-np.asarray(incoming,float)
        elif kind=='right_zero' and at==t0 and outgoing is not None:ft=-np.asarray(outgoing,float)
        elif kind=='private_continuation':
            if speed==0:ft=-np.asarray(incoming,float)
            elif slip@np.asarray(incoming,float)<0:ft=slip/speed
        return dict(force_definition='NONCONTACT_APPLIED_FORCE_ONLY',force_N=np.zeros((1,3)),contacts=[dict(
            particle_index=0,contact_id='fixed:1',normal=np.array([0.,0.,1.]),normal_force_N=3.,
            tangential_force_N=ft,wall_velocity_m_s=wall,wall_temperature_K=410.,conductance_W_K=.025,
            penetration_m=1e-6,normal_potential_J=1.5e-6,normal_stiffness_N_m=3e6)])
    def gap(*args):
        return dict(minimum_overlap_m=1e-6,maximum_overlap_m=1e-6,
                    fixed_normals=True,affine_tangential_wall_velocity=True)
    context={} if kind is None else {'fixed:1':dict(incoming_direction=np.asarray(incoming,float),endpoint=kind)}
    if outgoing is not None:context['fixed:1']['outgoing_direction']=np.asarray(outgoing,float)
    return dict(pa=pa,p=p,boundary=boundary,t0=t0,acc0=a0[None,:],acc1=a1[None,:],
                dt=dt,weight=.9,gap_certificate=gap,event_context=context)


def emit(name,**values):
    print('HHT_STOP_CONTACT_EVIDENCE '+json.dumps(dict(name=name,
        classification='PRESCRIBED_PATH_COMPONENT_NOT_PHYSICAL_TRAJECTORY',**values),sort_keys=True,allow_nan=False),flush=True)


class StopContactReferences(unittest.TestCase):
    def close(self,a,b):np.testing.assert_allclose(a,b,rtol=0.,atol=1e-12)

    def integrate_unchanged(self,case):
        before=particle_hash(case['pa']);p=copy.deepcopy(case['p'])
        result=tested.integrate(**case)
        self.assertEqual(particle_hash(case['pa']),before);self.assertEqual(case['p'],p)
        self.assertFalse(result['accepted'])
        # No contact integral helper may grant continuation/publication authority.
        self.assertFalse(result['restart_allowed'])
        self.assertFalse(result['can_publish'])
        return result

    def reject_unchanged(self,case):
        before=particle_hash(case['pa']);p=copy.deepcopy(case['p'])
        with self.assertRaises(ValueError):tested.integrate(**case)
        self.assertEqual(particle_hash(case['pa']),before);self.assertEqual(case['p'],p)

    def reference(self,result,impulse,heat,wall_work):
        self.close(result['actual_impulse_Ns'][0],impulse)
        self.close(result['friction_raw_J'],[heat])
        self.close(result['actual_material_work_J'],[wall_work-heat])
        self.close(result['actual_wall_work_J'],[wall_work])
        self.close(result['delta_impulse_Ns'],np.zeros((1,3)))
        self.close(result['delta_material_work_J'],np.zeros(1))

    def test_01_incoming_left_zero_analytic_constant_acceleration(self):
        h=1./3.;J=np.array([-1./3.,0.,1.]);Q=1./6.
        for wall,change in ((np.zeros(3),np.zeros(3)),
                            (np.array([.2,-.3,0.]),np.array([.1,.4,0.]))):
            with self.subTest(wall=wall.tolist()):
                r=self.integrate_unchanged(fixture(dt=h,wall0=wall,wall_change=change))
                self.reference(r,J,Q,float(J@(wall+.5*change)))
        emit('left_zero',expected_impulse_Ns=J.tolist(),expected_raw_heat_J=Q)

    def test_02_outgoing_right_zero_analytic_constant_acceleration(self):
        h=2./3.;J=np.array([2./3.,0.,2.]);Q=2./9.
        for wall,change in ((np.zeros(3),np.zeros(3)),
                            (np.array([.2,-.3,0.]),np.array([.1,.4,0.]))):
            with self.subTest(wall=wall.tolist()):
                r=self.integrate_unchanged(fixture(u0=(0.,0.,0.),a0=(-1.,0.,0.),dt=h,
                    wall0=wall,wall_change=change,kind='right_zero',outgoing=(-1.,0.,0.)))
                self.reference(r,J,Q,float(J@(wall+.5*change)))
        emit('right_zero',expected_impulse_Ns=J.tolist(),expected_raw_heat_J=Q)

    def test_03_private_crossing_keeps_negative_node_heat_without_certification(self):
        r=self.integrate_unchanged(fixture(dt=.5,kind='private_continuation'))
        self.reference(r,np.array([-.5,0.,1.5]),.125,0.)
        row=r['contacts'][0];nodes=row['quadrature']
        negative=[];rawsum=0.
        for node in nodes:
            # h*w*C*u(s), with h=.5 and u(s)=1-1.5s. No tested-function reference.
            expected=.5*node['weight']*(1.-1.5*node['fraction'])
            self.close(node['friction_heat_raw_J'],expected);rawsum+=node['friction_heat_raw_J']
            if expected<0:negative.append(node['friction_heat_raw_J'])
        self.assertGreater(len(negative),0);self.assertLess(min(negative),0.);self.close(rawsum,.125)
        self.assertFalse(row['whole_path_passive'])
        self.assertTrue(r['private_continuation'])
        self.assertFalse(r['physical_contact_path_certified'])
        self.assertFalse(r['whole_path_passive'])
        emit('private_collinear_crossing',raw_heat_J=float(r['friction_raw_J'][0]),
             negative_node_count=len(negative),minimum_node_heat_J=float(min(negative)),physical_contact_certified=False)

    def test_04_negative_total_private_heat_and_unmarked_crossing_reject(self):
        self.reject_unchanged(fixture(dt=1.,kind='private_continuation'))
        self.reject_unchanged(fixture(dt=.5,kind=None))

    def test_05_transverse_nonzero_and_noncollinear_private_reject(self):
        self.reject_unchanged(fixture(u0=(1.,.25,0.)))
        self.reject_unchanged(fixture(a0=(-3.,.2,0.),dt=.5,kind='private_continuation'))

    def test_06_left_limit_rejects_an_earlier_complete_zero(self):
        # Independent dense polynomial: u(s)=1-3s+2s², zeros at .5 and1.
        self.reject_unchanged(fixture(a0=(79.,0.,0.),a1=(-163./3.,0.,0.),dt=1.,incoming=(-1.,0.,0.)))

    def test_07_right_start_derivative_must_match_declared_outgoing(self):
        # u(s)=(s-2s²,s-s²): fullzero only at start, endpoint(-1,0).
        # Declared -x matches far endpoint but disagrees with actual initial
        # direction (1,1)/sqrt2. An endpoint-only direction check would miss it.
        self.reject_unchanged(fixture(u0=(0.,0.,0.),a0=(-81.,-40.,0.),
            a1=(157./3.,80./3.,0.),dt=1.,kind='right_zero',outgoing=(-1.,0.,0.)))
        self.reject_unchanged(fixture(incoming=(-1.,0.,0.)))

    def test_08_invalid_event_context_rejects_without_mutation(self):
        self.reject_unchanged(fixture(u0=(0.,0.,0.),a0=(-1.,0.,0.),dt=2./3.,kind='right_zero'))
        self.reject_unchanged(fixture(kind='invented_mode'))
        self.reject_unchanged(fixture(incoming=(0.,0.,0.)))
        case=fixture();case['gap_certificate']=None;self.reject_unchanged(case)


if __name__=='__main__':unittest.main(verbosity=2)

"""Independent references for the new local kick, not production authorization."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.contact_kick import kick, domain_kick
from afsd_backend.pysph.contact import WallContact
from test_contact_review import wall_parameters
from afsd_backend.pysph.kick_diagnostics import clean
from afsd_backend.pysph.output import restore, checkpoint, digest
from test_coupling_benchmarks import REVERSAL_REFERENCE
from review_artifacts import artifact_path

class EventContact(unittest.TestCase):
    def check_identity(self,r):
        for a in r['intervals']:
            scale=abs(a['delta_kinetic_J'])+abs(a['other_work_J'])+abs(a['constraint_work_J'])
            self.assertLess(abs(a['kinetic_work_residual_J']),1e-12+1e-10*scale)
            self.assertTrue(np.all(a['friction_heat_raw_J'] >= -1e-24))
        # Normal relative work must be included for compressive penalty kicks.
        contact_normal_work=sum(float(np.sum(a['normal_work_J'])) for a in r['intervals'])
        normal_wall_work=sum(float(np.sum(a['normal_impulse']*self.walls)) for a in r['intervals'])
        balance=(r['other_work_J']+r['constraint_work_J']+r['wall_work_J']
                 +contact_normal_work-normal_wall_work-float(np.sum(r['friction_heat_J'])))
        self.assertAlmostEqual(r['delta_kinetic_J'],balance,delta=1e-12+1e-10*abs(balance))

    def run_case(self,v=[1.,0.,0.],external=[-2.,0.,0.],walls=None,normals=None,fn=None,dt=1.,**kw):
        self.walls=np.array([[0.,0.,0.]]) if walls is None else np.array(walls,float)
        n=np.array([[0.,0.,1.]]) if normals is None else np.array(normals,float)
        fn=np.ones(len(n)) if fn is None else np.array(fn,float)
        # Counterbalance normal forces so the independent tangent references
        # have no normal kinetic-energy contribution, unless explicitly tested.
        external=np.array(external,float)-np.sum(fn[:,None]*n,axis=0)
        r=kick(v,1.,external,n,fn,self.walls,1.,dt,**kw)
        self.check_identity(r)
        return r

    def test_analytic_reversal_not_reflection_or_stall(self):
        r=self.run_case()
        self.assertAlmostEqual(r['events'][0]['time_s'],1/3,delta=1e-10)
        values=dict(u_final=r['velocity'][0],Q_friction=float(r['friction_heat_J'].sum()),
                    W_other=r['other_work_J'],delta_K=r['delta_kinetic_J'])
        for name,value in values.items():
            self.assertAlmostEqual(value,REVERSAL_REFERENCE[name],delta=1e-10,msg=name)
        self.assertEqual(len(r['intervals']),2)

    def test_static_wall_stops_without_rebound_and_then_sticks(self):
        r=self.run_case(external=[0.,0.,0.],dt=2.)
        np.testing.assert_allclose(r['velocity'],[0.,0.,0.],atol=1e-12)
        self.assertAlmostEqual(float(r['friction_heat_J'].sum()),.5,delta=1e-12)
        self.assertAlmostEqual(r['events'][0]['time_s'],1.,delta=1e-12)

    def test_moving_wall_galilean_work_identity(self):
        r=self.run_case(v=[1.5,0.,0.],external=[0.,0.,0.],walls=[[.5,0.,0.]],dt=2.)
        np.testing.assert_allclose(r['velocity'],[.5,0.,0.],atol=1e-12)
        self.assertAlmostEqual(float(r['friction_heat_J'].sum()),.5,delta=1e-12)
        self.assertAlmostEqual(r['wall_work_J'],-.5,delta=1e-12)
        self.assertAlmostEqual(r['delta_kinetic_J'],-1.,delta=1e-12)

    def test_corotating_zero_slip_at_collocated_point(self):
        x=np.array([.4,-.7,.1]); omega=np.array([0.,0.,2.])
        wall=np.array([.3,0.,0.])+np.cross(omega,x)
        r=self.run_case(v=wall,external=[0.,0.,0.],walls=[wall],dt=.1)
        np.testing.assert_allclose(r['velocity'],wall,atol=1e-12)
        np.testing.assert_allclose(r['friction_heat_J'],[0.],atol=1e-12)

    def test_two_contacts_actuator_coupling_exact_reversal(self):
        # Normals x,z; prescribed x=z=0, only y free. Two faces each muFn=1.
        # Fother_y=-4, u0_y=1: stop=1/6; reverse with a=-2 for 5/6.
        # Q=2*(1/12+25/36)=14/9, Wother=22/9, dK=8/9.
        r=self.run_case(v=[0.,1.,0.],external=[0.,-4.,0.],
            normals=[[1.,0.,0.],[0.,0.,1.]],walls=[[0.,0.,0.],[0.,0.,0.]],
            fn=[1.,1.],free_axes=[0.,1.,0.],prescribed_velocity=[0.,1.,0.])
        self.assertAlmostEqual(r['events'][0]['time_s'],1/6,delta=1e-10)
        np.testing.assert_allclose(r['velocity'],[0.,-5/3,0.],atol=1e-10)
        self.assertAlmostEqual(float(r['friction_heat_J'].sum()),14/9,delta=1e-10)
        self.assertAlmostEqual(r['other_work_J'],22/9,delta=1e-10)
        self.assertAlmostEqual(r['delta_kinetic_J'],8/9,delta=1e-10)
        # Shared y inverse mass =1 also in OFF-diagonal face blocks.
        self.assertEqual(r['contact_mobility'][0,1,1,1],1.)
        self.assertEqual(r['contact_mobility'][0,0,0,0],0.)

    def test_prescribed_motion_reaction_work_not_discarded(self):
        # x prescribed 1m/s; static wall friction -1N needs +1N reaction.
        r=self.run_case(v=[1.,0.,0.],external=[0.,0.,0.],dt=2.,
                        free_axes=[0.,1.,1.],prescribed_velocity=[1.,0.,0.])
        self.assertAlmostEqual(r['constraint_work_J'],2.,delta=1e-12)
        self.assertAlmostEqual(float(r['friction_heat_J'].sum()),2.,delta=1e-12)
        self.assertAlmostEqual(r['delta_kinetic_J'],0.,delta=1e-12)

    def test_curved_slip_against_independent_rk4(self):
        # Speed stays far from zero; plain independent ODE integration is valid.
        def rhs(v): return np.array([-.4,.7])-1.2*v/np.linalg.norm(v)
        u=np.array([1.,.4]); dt=.6/8192
        q=0.; other=0.
        # RK4 of velocity plus integral heat and external work, no UUT contact.
        state=np.array([*u,0.,0.])
        def derivative(z):
            v=z[:2]
            return np.r_[rhs(v),1.2*np.linalg.norm(v),np.array([-.4,.7])@v]
        for _ in range(8192):
            a=derivative(state); b=derivative(state+dt*a/2)
            c=derivative(state+dt*b/2); d=derivative(state+dt*c)
            state+=dt*(a+2*b+2*c+d)/6
        self.walls=np.zeros((1,3))
        r=kick([1.,.4,0.],1.,[-.4,.7,-1.2],[[0.,0.,1.]],[1.2],
               self.walls,1.,.6,rtol=1e-9)
        self.check_identity(r)
        np.testing.assert_allclose(r['velocity'][:2],state[:2],atol=2e-6,rtol=0.)
        self.assertAlmostEqual(float(r['friction_heat_J'].sum()),state[2],delta=2e-6)
        self.assertAlmostEqual(r['other_work_J'],state[3],delta=2e-6)

    def test_original_first_failed_particle_frozen_kick(self):
        path=ROOT/'outputs/coupling_20260930/phase1_verified'
        e=json.loads((path/'first_negative_event.json').read_text())
        pa,_=restore(path/'last_accepted.npz')
        idx=int(np.flatnonzero(pa.particle_id==e['particle_id'])[0])
        dt=e['dt_s']; normal=np.array(e['normal'])
        fn=float(np.array(e['normal_contact_impulse'])@normal)/dt
        self.walls=np.array([e['actual_wall_velocity_used']])
        r=kick(e['v_before_internal'],pa.m[idx],np.array(e['internal_impulse'])/dt,
            [normal],[fn],self.walls,.3,dt)
        self.check_identity(r)
        target=artifact_path('first_particle_repaired.json')
        target.write_text(json.dumps(clean(r),indent=2)+'\n')
        self.assertGreaterEqual(float(r['friction_heat_J'].sum()),0.)
        print('FIRST_FROZEN_KICK',json.dumps(clean({k:v for k,v in r.items() if k!='intervals'})))

    def test_first_particle_independent_vector_stop_and_heat(self):
        path=ROOT/'outputs/coupling_20260930/phase1_verified'
        e=json.loads((path/'first_negative_event.json').read_text())
        pa,_=restore(path/'last_accepted.npz')
        idx=int(np.flatnonzero(pa.particle_id==e['particle_id'])[0])
        mass=pa.m[idx]; dt=e['dt_s']; normal=np.array(e['normal'])
        tangent=np.eye(3)-np.outer(normal,normal)
        other=np.array(e['internal_impulse'])/dt
        fn=float(np.array(e['normal_contact_impulse'])@normal)/dt
        wall=np.array(e['actual_wall_velocity_used'])
        u0=tangent@(np.array(e['v_before_internal'])-wall)
        drive=tangent@other
        capacity=.3*fn
        self.assertLess(np.linalg.norm(drive),capacity) # sticking is feasible
        # Closed polar-coordinate integral of m*du/dt=drive-c*u/|u|.
        # With a=c/|drive|, r(theta) proportional to
        # tan(theta/2)**a / sin(theta). Integrate time and c*r*dt to
        # theta=0. This reference does not invoke the tested cone/events.
        speed=np.linalg.norm(u0); drive_norm=np.linalg.norm(drive)
        cosine=float(u0@drive)/(speed*drive_norm)
        positive=(1+cosine)/2; negative=(1-cosine)/2
        expected_stop=mass*speed*(positive/(capacity-drive_norm)+negative/(capacity+drive_norm))
        expected_q=mass*capacity*speed**2*(
            positive**2/(2*(capacity-drive_norm))+positive*negative/capacity
            +negative**2/(2*(capacity+drive_norm)))
        evidence=[]
        for accuracy,tolerance in ((2e-7,1e-3),(1e-10,1e-5)):
            r=kick(e['v_before_internal'],mass,other,[normal],[fn],[wall],.3,dt,rtol=accuracy)
            heat=float(r['friction_heat_J'].sum())
            self.assertLess(abs(heat/expected_q-1),tolerance)
            self.assertLess(np.linalg.norm(tangent@(r['velocity']-wall)),1e-10)
            self.assertTrue(r['events'])
            actual_stop=r['events'][0]['time_s']
            self.assertLess(abs(actual_stop/expected_stop-1),tolerance)
            evidence.append(dict(rtol=accuracy,heat_J=heat,stop_s=actual_stop,
                relative_heat_error=abs(heat/expected_q-1),relative_stop_error=abs(actual_stop/expected_stop-1)))
        target=artifact_path('first_particle_independent_reference.json')
        target.write_text(json.dumps(dict(reference='independent closed polar integral of vector Coulomb ODE',
            Q_J=expected_q,t_stop_s=expected_stop,
            comparisons=evidence),indent=2)+'\n')

    def test_rejected_checkpoint_cannot_restart(self):
        path=ROOT/'outputs/coupling_20260930/phase1_verified/REJECTED_partial_integrator_state.npz'
        with self.assertRaisesRegex(ValueError,'REJECTED'):
            restore(path)

    def test_domain_adapter_reversal_heat_partition_and_work(self):
        p=wall_parameters(feed=0.)
        velocity=np.array([[1.,0.,0.]])
        acceleration,heat,ledger,actuator_work,traces=domain_kick(
            WallContact(100.,1.),np.array([[0.,0.,-.91]]),velocity,
            np.array([1.]),np.array([350.]),0.,1.,p,
            np.array([[-2.,0.,-1.]]),np.array([0]),np.array([47]))
        np.testing.assert_allclose(velocity+acceleration,[[-2/3,0.,0.]],atol=1e-10)
        self.assertAlmostEqual(heat[0],.4*7/18,delta=1e-10)
        self.assertAlmostEqual(ledger['wall_heat_rate'],.6*7/18,delta=1e-10)
        self.assertAlmostEqual(ledger['frozen_internal_impulse_work_J'],1/9,delta=1e-10)
        self.assertAlmostEqual(actuator_work,0.,delta=1e-10)
        self.assertEqual(traces[0]['particle_id'],47)

    def test_interrupted_checkpoint_keeps_previous_accepted_file(self):
        import tempfile
        from unittest.mock import patch
        reference=ROOT/'outputs/coupling_20260930/phase1_verified/last_accepted.npz'
        pa,state=restore(reference)
        with tempfile.TemporaryDirectory() as directory:
            target=Path(directory)/'accepted.npz'
            checkpoint(target,pa,state)
            before=digest(target)
            with patch('afsd_backend.pysph.output.dump',side_effect=RuntimeError('simulated interrupted write')):
                with self.assertRaisesRegex(RuntimeError,'interrupted'):
                    checkpoint(target,pa,dict(state,count=999))
            self.assertEqual(digest(target),before)
            _,recovered=restore(target)
            self.assertEqual(recovered['count'],31)

if __name__=='__main__':
    unittest.main(verbosity=2)


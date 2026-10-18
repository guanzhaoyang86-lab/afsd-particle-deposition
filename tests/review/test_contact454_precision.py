"""Independent near-zero/coupled mobility regressions for the 454 repair."""
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.contact_kick import kick, cone_solve, tangent_bases, ContactConvergenceError
from afsd_backend.pysph.kick_diagnostics import clean
from review_artifacts import artifact_path

def rotation():
    axis=np.array([2.,-1.,3.])/math.sqrt(14.)
    a=.713;c=math.cos(a);s=math.sin(a)
    K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    return c*np.eye(3)+(1-c)*np.outer(axis,axis)+s*K

def polar_reference(u,drive,mass,capacity):
    speed=float(np.linalg.norm(u)); force=float(np.linalg.norm(drive))
    if force==0:return mass*speed/capacity,.5*mass*speed*speed
    cosine=float(u@drive)/(speed*force)
    b=(1+cosine)/2;d=(1-cosine)/2
    time=mass*speed*(b/(capacity-force)+d/(capacity+force))
    heat=mass*capacity*speed**2*(b*b/(2*(capacity-force))+b*d/capacity+d*d/(2*(capacity+force)))
    return time,heat

class Contact454Precision(unittest.TestCase):
    def test_frozen_cone_decimal_reference_and_no_normal_leakage(self):
        record=json.loads((ROOT/'outputs/coupling_20260930/new_candidate_RED.json').read_text())
        c=record['failed_cone']
        raw=json.loads((ROOT/'outputs/coupling_20260930/event_dx05_cfl01/REJECTED_contact_debug.json').read_text())['failed_local_context']
        n=np.array(raw['faces'][0]['normal']);base=np.array(c['base'][0])
        # Decimal constructs an independent basis for this normal, nz=0.
        with localcontext() as context:
            context.prec=70
            nx,ny=[Decimal(str(v)) for v in n[:2]]
            length=(nx*nx+ny*ny).sqrt();nx/=length;ny/=length
            x,y,z=[Decimal(str(v)) for v in base]
            q=[-ny*x+nx*y,z]
            lam=Decimal(str(c['coefficient']))/Decimal(str(raw['mass_kg']))
            f=[-v/lam for v in q];norm=sum(v*v for v in f).sqrt()
            limit=Decimal(str(c['limits_N'][0]))
            if norm>limit:f=[v*limit/norm for v in f]
            expected=np.array([float(-ny*f[0]),float(nx*f[0]),float(f[1])])
        T=tangent_bases([n])
        b=np.einsum('nki,k->ni',T,base)
        force=cone_solve(b,T,np.array(c['mobility']),c['coefficient'],np.array(c['limits_N']))[0]
        np.testing.assert_allclose(force,expected,atol=1e-11,rtol=1e-10)
        self.assertLess(abs(float(force@n)),1e-14*np.linalg.norm(force))
        artifact_path('cone_precision.json').write_text(json.dumps(dict(
            reference_N=expected.tolist(),actual_N=force.tolist(),
            error_N=float(np.linalg.norm(force-expected)),normal_leakage_N=float(force@n)),indent=2)+'\n')

    def test_single_contact_constrained_effective_mass(self):
        # n=(x+y)/sqrt(2), x constrained, y,z free, mass=2.
        # other force=(0,2,2), zero slip: equilibrium y,z fixes ft_y=ft_z=-2;
        # tangency fixes ft_x=+2. Isolated scalar mass gives the WRONG y force.
        n=np.array([1.,1.,0.])/math.sqrt(2.)
        T=tangent_bases([n]);M=np.diag([0.,.5,.5])
        base=np.einsum('nki,k->ni',T,np.array([0.,1.,1.]))
        result=cone_solve(base,T,M,1.,np.array([10.]))
        np.testing.assert_allclose(result,[[2.,-2.,-2.]],atol=1e-10,rtol=1e-10)
        np.testing.assert_allclose(M@(np.array([0.,2.,2.])+result[0]),0.,atol=1e-10)

    def test_two_faces_share_one_free_dof_in_stick(self):
        # x,z constrained; y acceleration -2.5 must be cancelled.
        # Capacities 1 and 2. Split is NOT unique: (1,1.5) and (.75,1.75)
        # both satisfy equilibrium. Test physics, not a solver-specific split.
        T=tangent_bases([[1.,0.,0.],[0.,0.,1.]])
        M=np.diag([0.,1.,0.]);base=np.einsum('nki,k->ni',T,np.array([0.,-2.5,0.]))
        f=cone_solve(base,T,M,1.,np.array([1.,2.]))
        self.assertTrue(np.all(np.linalg.norm(f,axis=1)<=np.array([1.,2.])+1e-10))
        self.assertAlmostEqual(float(f[:,1].sum()),2.5,delta=1e-10)
        self.assertAlmostEqual(float(f[0,0]),0.,delta=1e-10)
        self.assertAlmostEqual(float(f[1,2]),0.,delta=1e-10)
        self.assertLess(np.linalg.norm(M@(np.array([0.,-2.5,0.])+f.sum(axis=0))),1e-10)

    def test_captured_454_rotations_tighter_accuracy_and_work(self):
        c=json.loads((ROOT/'outputs/coupling_20260930/event_dx05_cfl01/REJECTED_contact_debug.json').read_text())['failed_local_context']
        n=np.array(c['faces'][0]['normal']);wall=np.array(c['faces'][0]['wall_velocity'])
        v=np.array(c['velocity_before']);other=np.array(c['internal_force'])
        fn=2.63e12*.0005**2*c['faces'][0]['penetration'];mass=c['mass_kg']
        expected_time=5.985479017657219e-9
        expected_heat=8.164870224344708e-12
        evidence=[]
        for Q in (np.eye(3),rotation(),np.diag([-1.,-1.,1.])):
            for accuracy,tolerance in ((2e-7,1e-3),(1e-10,1e-5)):
                result=kick(Q@v,mass,Q@other,[Q@n],[fn],[Q@wall],.3,c['dt_s'],rtol=accuracy)
                self.assertTrue(result['events'])
                stop=result['events'][0]['time_s'];heat=float(result['friction_heat_J'].sum())
                self.assertLess(abs(stop/expected_time-1),tolerance)
                self.assertLess(abs(heat/expected_heat-1),tolerance)
                final=Q.T@result['velocity']-wall
                self.assertLess(np.linalg.norm(np.cross(n,final)),1e-10)
                for interval in result['intervals']:
                    raw=-float(np.sum(interval['tangential_impulse']*interval['tangential_relative_midpoint']))
                    self.assertAlmostEqual(float(interval['friction_heat_raw_J'].sum()),raw,delta=1e-25)
                # Full impulse identity independent of work evaluator.
                impulse=(c['dt_s']*(Q@other)+result['normal_impulse'].sum(axis=0)
                         +result['tangential_impulse'].sum(axis=0)+result['constraint_impulse'])
                np.testing.assert_allclose(mass*(result['velocity']-Q@v),impulse,atol=1e-22,rtol=1e-10)
                evidence.append(dict(rotation=Q,rtol=accuracy,heat_J=heat,stop_s=stop,
                    time_error=abs(stop/expected_time-1),heat_error=abs(heat/expected_heat-1),
                    intervals=len(result['intervals']),events=result['events']))
        artifact_path('candidate454_covariance.json').write_text(json.dumps(clean(evidence),indent=2)+'\n')

    def test_manufactured_curved_stops_do_not_rebound(self):
        evidence=[];Q=rotation()
        for fraction in (0.,.02,.5,.9):
            mass=.5;capacity=.8
            u=np.array([.7,.2,0.]);drive=capacity*fraction*np.array([-.6,.8,0.])
            stop,heat=polar_reference(u,drive,mass,capacity)
            wall=Q@np.array([.2,-.1,.3]);normal=Q@np.array([0.,0.,1.])
            try:
                result=kick(wall+Q@u,mass,Q@drive-capacity*normal,[normal],[capacity],[wall],1.,1.3*stop)
            except ContactConvergenceError as error:
                context=dict(fraction=fraction,reference_stop_s=stop,
                    error=str(error),interval_context=getattr(error,'interval_context',None))
                artifact_path('manufactured_failure.json').write_text(json.dumps(clean(context),indent=2)+'\n')
                print('MANUFACTURED_FAILURE',json.dumps(clean(context)))
                raise
            self.assertTrue(result['events'])
            self.assertLess(abs(result['events'][0]['time_s']/stop-1),1e-3)
            self.assertLess(abs(float(result['friction_heat_J'].sum())/heat-1),1e-3)
            self.assertLess(np.linalg.norm(np.cross(normal,result['velocity']-wall)),1e-10)
            evidence.append(dict(drive_over_capacity=fraction,reference_stop_s=stop,reference_heat_J=heat,
                actual_stop_s=result['events'][0]['time_s'],actual_heat_J=float(result['friction_heat_J'].sum())))
        artifact_path('manufactured_stops.json').write_text(json.dumps(evidence,indent=2)+'\n')

if __name__=='__main__':unittest.main(verbosity=2)


"""Fixed class-level references; no solver-generated expected answer."""
import itertools
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.contact_kick import cone_solve,kick,ContactConvergenceError
from afsd_backend.pysph.kick_diagnostics import clean
from review_artifacts import artifact_path


def independent_basis(normal):
    n=np.asarray(normal,float);n/=np.linalg.norm(n)
    axis=np.eye(3)[np.argmin(np.abs(n))]
    a=axis-n*float(axis@n);a/=np.linalg.norm(a)
    return np.column_stack((a,np.cross(n,a)))


def rotation():
    n=np.array([2.,-1.,3.])/math.sqrt(14.)
    K=np.array([[0.,-n[2],n[1]],[n[2],0.,-n[0]],[-n[1],n[0],0.]])
    angle=.713
    return np.eye(3)*math.cos(angle)+(1-math.cos(angle))*np.outer(n,n)+math.sin(angle)*K


def oracle(base,T,M,c,limits,cartesian):
    """Independent direct algebra and disk support function, not UUT stopping."""
    force=np.array([T[a].T@cartesian[a] for a in range(len(T))])
    A=np.concatenate(list(T),axis=1)
    H=c*(A.T@M@A)
    fs=max(float(np.max(limits)),np.finfo(float).tiny)
    # Use the full physical mobility scale, not a tangent block which can be
    # exactly zero in the continuum but retain cancellation noise after rotation.
    vs=max(float(np.linalg.norm(base)),c*float(np.linalg.norm(M,2))*fs)
    if vs==0.:vs=1. # identically zero objective; no tiny-denominator overflow
    x=force/fs;b=base/vs;h=H*(fs/vs);r=limits/fs
    g=(b.ravel()+h@x.ravel()).reshape(-1,2)
    gap=float(np.sum(g*x)+np.sum(r*np.linalg.norm(g,axis=1)))
    projected=[]
    for radius,y in zip(r,x-g):
        length=np.linalg.norm(y)
        projected.append(y if length<=radius else y*(radius/length))
    residual=float(np.linalg.norm(x-np.array(projected)))
    feasibility=float(max(0.,np.max(np.linalg.norm(x,axis=1)-r)))
    objective=float(np.sum(b*x)+.5*x.ravel()@h@x.ravel())
    normal_error=float(np.linalg.norm(cartesian-np.einsum('nki,ni->nk',T,force))/fs)
    return dict(gap=gap,stationarity=residual,feasibility=feasibility,
                objective=objective,normal_error=normal_error,
                objective_scale=fs*vs)


def matrix():
    for N,rank,mode,geometry,scale,rot in itertools.product(
            (1,2,3),range(4),('interior','mixed','sliding'),
            ('orthogonal','near_parallel'),(1e-6,1.,1e6),(False,True)):
        Q=rotation() if rot else np.eye(3)
        normals=np.eye(3)[:N] if geometry=='orthogonal' else np.array(
            [[1.,k*1e-4,(k%2)*2e-4] for k in range(N)])
        T=np.array([Q@independent_basis(n) for n in normals])
        M=Q@np.diag([1. if k<rank else 0. for k in range(3)])@Q.T
        limits=np.geomspace(1e-3,1.,N) if N>1 else np.ones(1)
        target=[];multipliers=[]
        for k,cap in enumerate(limits):
            d=np.array([math.cos(.4+k),math.sin(.4+k)])
            active=mode=='sliding' or (mode=='mixed' and k%2==0)
            target.append(cap*d*(1. if active else .4))
            multipliers.append(1.+.3*k if active else 0.)
        target=np.array(target);A=np.concatenate(list(T),axis=1)
        H=A.T@M@A
        # KKT construction BEFORE calling the tested solver.
        base=-(H@target.ravel()).reshape(N,2)-np.array(multipliers)[:,None]*target
        order=np.arange(N)[::-1] if rot else np.arange(N)
        yield dict(name=f'n{N}_rank{rank}_{mode}_{geometry}_s{scale}_r{rot}',
            T=T[order],M=M,c=1/(scale*scale),base=base[order]/scale,
            limits=limits[order]*scale,reference=target[order]*scale)


class ContactClass(unittest.TestCase):
    def test_frozen_multiface_holdout_full_kick_identity(self):
        # Holdout439 includes the full frozen kick, not only its failed QP.
        # This checks completion/impulse/work, not an invented exact trajectory.
        path=ROOT/'outputs/contact454_20260930/fixed_dx05_cfl01/REJECTED_contact_debug.json'
        c=json.loads(path.read_text())['failed_local_context'];faces=c['faces']
        mass=c['mass_kg'];dt=c['dt_s'];mu=.3
        fn=np.array(c['failed_cone']['limits_N'])/mu
        try:
            result=kick(c['velocity_before'],mass,c['internal_force'],
                [a['normal'] for a in faces],fn,[a['wall_velocity'] for a in faces],
                mu,dt,free_axes=c['free_axes'],prescribed_velocity=c['prescribed_velocity'])
        except ContactConvergenceError as error:
            artifact_path('holdout439_full_kick.json').write_text(json.dumps(clean(dict(
                status='FAIL',error=str(error),cone=getattr(error,'cone_context',None),
                interval=getattr(error,'interval_context',None))),indent=2)+'\n')
            self.fail(str(error))
        v0=np.array(c['velocity_before']);v1=result['velocity']
        impulse=dt*np.array(c['internal_force'])+result['normal_impulse'].sum(axis=0)+\
            result['tangential_impulse'].sum(axis=0)+result['constraint_impulse']
        np.testing.assert_allclose(mass*(v1-v0),impulse,atol=1e-22,rtol=1e-9)
        self.assertLess(abs(result['kinetic_work_residual_J']),1e-20)
        self.assertAlmostEqual(sum(a['dt_s'] for a in result['intervals']),dt,delta=1e-22)
        for a in result['intervals']:
            actual=-np.sum(a['tangential_impulse']*a['tangential_relative_midpoint'],axis=1)
            np.testing.assert_allclose(a['friction_heat_raw_J'],actual,atol=1e-28,rtol=1e-10)
        artifact_path('holdout439_full_kick.json').write_text(json.dumps(clean(dict(
            status='PASS_LOCAL_IDENTITY_NOT_TRAJECTORY_REFERENCE',result=result)),indent=2)+'\n')

    def test_manufactured_convex_kkt_matrix(self):
        evidence=[]
        for case in matrix():
            name=case['name'];T=case['T'];reference=case['reference']
            expected=np.einsum('nki,ni->nk',T,reference)
            cert=oracle(case['base'],T,case['M'],case['c'],case['limits'],expected)
            self.assertLess(cert['stationarity'],2e-14,'invalid independent KKT fixture')
            record=dict(name=name,reference_certificate=cert)
            try:
                f=cone_solve(case['base'],T,case['M'],case['c'],case['limits'])
                measured=oracle(case['base'],T,case['M'],case['c'],case['limits'],f)
                vs=cert['objective_scale']/max(case['limits'])
                endpoint_error=np.linalg.norm(case['c']*case['M']@(f-expected).sum(axis=0))/vs
                record.update(measured=measured,scaled_endpoint_error=float(endpoint_error),status='PASS')
                checks=(measured['feasibility']<=2e-12,
                        abs(measured['gap'])<=2e-10,measured['stationarity']<=2e-10,
                        abs(measured['objective']-cert['objective'])<=2e-10,
                        measured['normal_error']<=2e-12,endpoint_error<=2e-9)
                if not all(checks):record['status']='FAIL_CERTIFICATE'
            except (ContactConvergenceError,ValueError) as error:
                record.update(status='FAIL_SOLVE',error=str(error))
            evidence.append(record)
        artifact_path('contact_class_matrix.json').write_text(json.dumps(evidence,indent=2)+'\n')
        failures=[a['name']+':'+a['status'] for a in evidence if a['status']!='PASS']
        print('CLASS_MATRIX',len(evidence),'cases',len(failures),'failures')
        self.assertEqual(failures,[],str(failures[:12]))

    def test_two_distinct_static_force_decompositions_are_valid(self):
        T=np.array([independent_basis([1.,0.,0.]),independent_basis([0.,0.,1.])])
        M=np.diag([0.,1.,0.]);b=np.array([t.T@np.array([0.,-2.5,0.]) for t in T])
        certs=[]
        for forces in ([[0.,1.,0.],[0.,1.5,0.]],[[0.,.75,0.],[0.,1.75,0.]]):
            cert=oracle(b,T,M,1.,np.array([1.,2.]),np.array(forces))
            self.assertLess(abs(cert['gap']),1e-14)
            self.assertLess(cert['stationarity'],1e-14)
            self.assertLess(cert['feasibility'],1e-14)
            certs.append(cert)
        self.assertAlmostEqual(certs[0]['objective'],certs[1]['objective'],delta=1e-14)

    def test_zero_capacity_and_zero_mobility(self):
        T=np.array([independent_basis([0.,0.,1.])]*3)
        b=np.array([[1.,2.],[0.,0.],[-3.,4.]])
        r=np.array([0.,1.,2.])
        f=cone_solve(b,T,np.zeros((3,3)),1.,r)
        cert=oracle(b,T,np.zeros((3,3)),1.,r,f)
        self.assertLess(abs(cert['gap']),1e-14)
        self.assertLess(cert['feasibility'],1e-14)

    def test_bad_inputs_fail_and_zero_coefficient_is_linear(self):
        T=np.array([independent_basis([0.,0.,1.])]);b=np.array([[3.,4.]])
        for M,c,caps in ((-np.eye(3),1.,[1.]),(np.eye(3),-1.,[1.]),
                         (np.eye(3),1.,[-1.]),(np.eye(3)*np.nan,1.,[1.])):
            with self.subTest(c=c,caps=caps):
                with self.assertRaises(ValueError):cone_solve(b,T,M,c,np.array(caps))
        f=cone_solve(b,T,np.eye(3),0.,np.array([2.]))
        np.testing.assert_allclose(T[0].T@f[0],[-1.2,-1.6],atol=1e-13)

    def test_roundoff_heat_remains_raw_and_bound_is_frame_invariant(self):
        # Same manufactured curved stop, now audit the recorded rounding budget.
        Q=rotation();mass=.5;cap=.8;u=np.array([.7,.2,0.])
        drive=cap*.5*np.array([-.6,.8,0.]);normal=Q@np.array([0.,0.,1.])
        expected_time=.49834249077337645
        records=[]
        for boost in (np.zeros(3),np.array([1.,-2.,3.])):
            wall=Q@np.array([.2,-.1,.3])+boost
            result=kick(wall+Q@u,mass,Q@drive-cap*normal,[normal],[cap],[wall],1.,1.3*expected_time)
            for a in result['intervals']:
                v=a['v_before']-a['reference_frame_velocity']
                w=a['v_after']-a['reference_frame_velocity']
                J=a['tangential_impulse']
                # Walls stationary in the chosen frame for this fixture.
                bound=256*np.finfo(float).eps*np.linalg.norm(J,axis=1)*(
                    a['relative_speed_scale_m_s']+np.linalg.norm(v)+np.linalg.norm(w))
                np.testing.assert_allclose(a['friction_rounding_bound_J'],bound,rtol=1e-12,atol=1e-28)
                raw=-np.sum(J*a['tangential_relative_midpoint'],axis=1)
                np.testing.assert_allclose(a['friction_heat_raw_J'],raw,rtol=1e-12,atol=1e-28)
                self.assertTrue(np.all(a['friction_heat_raw_J']>=-bound-1e-28))
            records.append(dict(boost=boost.tolist(),heat_J=float(sum(result['friction_heat_J'])),
                minimum_raw_J=min(float(min(a['friction_heat_raw_J'])) for a in result['intervals'])))
        self.assertAlmostEqual(records[0]['heat_J'],records[1]['heat_J'],delta=1e-10)
        artifact_path('raw_roundoff_heat.json').write_text(json.dumps(records,indent=2)+'\n')

    def test_kick_parallel_faces_reversal_rotations_constraints(self):
        for N,constrained,rot in itertools.product((1,2,3),(False,True),(False,True)):
            # Axis permutation preserves representable diagonal actuator axes.
            Q=np.array([[0.,1.,0.],[0.,0.,1.],[1.,0.,0.]]) if rot else np.eye(3)
            n=Q@np.array([0.,0.,1.]);d=Q@np.array([1.,0.,0.])
            wall=Q@np.array([.5,-.25,.125]);v=wall+d
            free=np.abs(d) if constrained else np.ones(3)
            with self.subTest(N=N,constrained=constrained,rotated=rot):
                result=kick(v,1.,-2*d-n,[n]*N,[1/N]*N,[wall]*N,1.,1.,
                            free_axes=free,prescribed_velocity=v)
                np.testing.assert_allclose(result['velocity'],wall-2*d/3,atol=1e-10,rtol=0.)
                self.assertAlmostEqual(sum(result['friction_heat_J']),7/18,delta=1e-10)
                self.assertAlmostEqual(result['events'][0]['time_s'],1/3,delta=1e-10)
                self.assertLess(abs(result['kinetic_work_residual_J']),1e-10)

    def test_different_wall_velocities_counterbalance(self):
        # +/-1m/s walls exert +/-1N: particle stays at rest, heat=2*dt.
        result=kick([0.,0.,0.],1.,[0.,0.,-2.],[[0.,0.,1.]]*2,[1.,1.],
                    [[1.,0.,0.],[-1.,0.,0.]],1.,.2)
        np.testing.assert_allclose(result['velocity'],0.,atol=1e-12)
        np.testing.assert_allclose(result['friction_heat_J'],[.2,.2],atol=1e-12)
        self.assertAlmostEqual(result['wall_work_J'],.4,delta=1e-12)

    def test_normal_common_boost_does_not_change_tangential_reversal(self):
        # Large boost NORMAL to plane: x=1 is still represented exactly.
        # Do not permit an O(eps*boost) allowance on unaffected x/heat fields.
        evidence=[]
        for boost in (0.,1e3,1e14):
            result=kick([1.,0.,boost],1.,[-2.,0.,-1.],[[0.,0.,1.]],[1.],
                        [[0.,0.,boost]],1.,1.)
            row=dict(boost=boost,final_x=float(result['velocity'][0]),
                     heat_J=float(sum(result['friction_heat_J'])),events=result['events'])
            evidence.append(row)
        artifact_path('common_boost.json').write_text(json.dumps(evidence,indent=2)+'\n')
        for row in evidence:
            self.assertAlmostEqual(row['final_x'],-2/3,delta=1e-10)
            self.assertAlmostEqual(row['heat_J'],7/18,delta=1e-10)


if __name__=='__main__':unittest.main(verbosity=2)

"""Independent one-way J2-normal/vector-friction assembly; NOT an AFSD oracle."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from probe_engineering_material import PerfectPlasticBenchmark
from afsd_backend.pysph.material import update
from probe_process_friction_cone import (cone_response,cone_rows,advance_cone,
    vector_contact_work,feasible_cone,solve_cone_modes)
from probe_process_plastic_friction import PlasticFrictionBlock

A=np.array([.6,.8])
LIMITS=dict(force_tolerance_N=1e-12,slip_tolerance_m=1e-12,
    maximum_active_sets=16,maximum_events_per_interval=12,event_bisections=48,
    event_time_tolerance_s=1e-14,mode_probe_fraction=1e-4)


class NormalJ2TangentOwner:
    """Test assembly, not a replacement material/contact implementation.

    Actual material.update returns the normal force for each pure candidate.
    Tested Q6 cone rows solve two unknown tangent DOFs and the tested event
    helper locates sliding. All analytic references live only in assertions.
    Normal motion is prescribed; this is one-way capacity coupling, not one
    material point with simultaneous normal and shear strain.
    """
    def __init__(self,stationary=False,basis_angle=0.,kappa=1000.):
        self.time=0.;self.count=0;self.nodes=np.array([0]);self.mu=.2
        self.kappa=kappa;self.limits=dict(LIMITS);self.stationary=stationary
        self.S=np.zeros((1,3,3));self.ep=np.zeros(1);self.u=np.zeros(2)
        self.normal_force=np.zeros(1);self.tangent_force=np.zeros((1,2))
        self.plastic=self.returned=self.heat=self.wall_work=self.material_work=self.normal_work=0.
        self.rows=[];self.law=PerfectPlasticBenchmark(3.,293.15)
        angle=np.deg2rad(basis_angle)
        self.basis=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        # The normal placeholder is prescribed and removed from the solve.
        # dN/du_t=0 for this declared assembly; no frozen 3D normal law implied.
        self.K=np.diag([0.,1000.,1000.]);self.nd=np.array([0]);self.td=np.array([[1,2]])

    @staticmethod
    def strain(t):
        if not np.isfinite(t) or t<0 or t>3:raise ValueError('outside fixture program')
        return .001*t if t<=1 else .001+.002*(t-1) if t<=2 else .003-.001*(t-2)

    def wall(self,t):
        scalar=0. if t<=1 or self.stationary else .004*(t-1) if t<=2 else .004
        return np.r_[-self.strain(t),scalar*A]

    def normal_trial(self,t):
        if t<=self.time:raise ValueError('candidate must advance')
        dt=t-self.time; gradient=np.zeros((1,3,3))
        gradient[0,0,0]=-(self.strain(t)-self.strain(self.time))/dt
        S,ep,plastic,returned,rate=update(self.S,self.ep,gradient,np.array([293.15]),dt,1500.,self.law)
        sigma=S-1000.*self.strain(t)*np.eye(3)
        return dict(S=S,ep=ep,plastic=plastic,returned=returned,rate=rate,
            stress=sigma,normal_force=-sigma[:,0,0])

    def branch(self,t,released):
        released=np.asarray(released,dtype=bool);dw=self.wall(t)-self.wall(self.time)
        u=self.u.copy() if released[0] else self.u+dw[1:]
        material=self.normal_trial(t)
        def state(v):
            tau=1000.*v[None,:];slip=(v-self.u-dw[1:])[None,:]
            response=cone_response(tau,material['normal_force'],slip,self.mu,self.kappa,self.basis)
            return dict(base_step=self.count,time_s=t,u=v.copy(),tangent_force=tau,
                relative_increment=slip,response=response,**material)
        for _ in range(20):
            c=state(u)
            if not released[0] or np.linalg.norm(c['response']['residual'])<1e-13:return c
            residual=np.r_[0.,c['response']['residual'].ravel()]
            K,r=cone_rows(self.K,residual,self.nd,self.td,c['response'],self.mu,self.kappa,released)
            step=np.linalg.solve(K[1:,1:],r[1:]);before=np.linalg.norm(residual[1:])
            for half in range(12):
                trial=u+step*2.**(-half)
                if np.linalg.norm(state(trial)['response']['residual'])<before:
                    u=trial;break
            else:raise RuntimeError('fixture bounded cone line search')
        raise RuntimeError('fixture bounded cone Newton')

    def commit(self,c):
        if c['base_step']!=self.count:raise ValueError('stale fixture candidate')
        if not feasible_cone(c,np.array([True]),self.mu,self.limits):
            raise ValueError('fixture infeasible contact')
        if np.linalg.norm(c['response']['residual'])>self.limits['force_tolerance_N']:
            raise ValueError('fixture unconverged cone root')
        oldforce=np.r_[-self.normal_force[0],self.tangent_force[0]][None,:]
        force=np.r_[-c['normal_force'][0],c['tangent_force'][0]][None,:]
        de=self.strain(c['time_s'])-self.strain(self.time)
        du=np.r_[-de,c['u']-self.u][None,:]
        work=vector_contact_work(oldforce,force,du,self.wall(c['time_s'])-self.wall(self.time),self.nodes,c['time_s']-self.time)
        q=float(work['friction_increment_J'][0]);P=float(c['plastic'][0])
        self.normal_work+=float(.5*(self.normal_force[0]+c['normal_force'][0])*de)
        self.wall_work+=work['wall_work_J'];self.material_work+=work['material_work_J']
        self.heat+=q;self.plastic+=P;self.returned+=float(c['returned'][0])
        self.rows.append(dict(t=c['time_s'],kind=c['event_kind'],released=c['released'].copy(),
            normal=c['normal_force'].copy(),u=c['u'].copy(),tau=c['tangent_force'].copy(),q=q,plastic=P))
        self.S=c['S'].copy();self.ep=c['ep'].copy();self.u=c['u'].copy()
        self.normal_force=c['normal_force'].copy();self.tangent_force=c['tangent_force'].copy()
        self.time=c['time_s'];self.count+=1

    def snapshot(self):
        return dict(S=self.S.copy(),ep=self.ep.copy(),u=self.u.copy(),N=self.normal_force.copy(),
            tau=self.tangent_force.copy(),ledger=np.array([self.plastic,self.returned,self.heat,
                self.wall_work,self.material_work,self.normal_work]),time=self.time,count=self.count,rows=len(self.rows))


class IndependentPlasticFrictionQ7(unittest.TestCase):
    def close(self,a,b,atol=1e-10):np.testing.assert_allclose(a,b,atol=atol,rtol=0)

    def test_01_actual_j2_current_reaction_and_independent_normal_energy(self):
        b=NormalJ2TangentOwner(stationary=True)
        advance_cone(b,1.);advance_cone(b,2.)
        self.close(b.S,[np.diag([-2.,1.,1.])]);self.close(b.ep,[1/750],atol=1e-12)
        self.close(b.normal_force,[5.]);self.close(b.tangent_force,[[0.,0.]])
        self.close([b.normal_work,b.plastic,b.heat,b.wall_work],[.0095,.004,0.,.0095])
        En=float(np.sum(b.S*b.S)/(4*1500.)+.5*1000.*.003**2)
        self.close(En,.0055);self.close(b.normal_work-En-b.plastic,0.)

    def test_02_simultaneous_plasticity_friction_exact_event_and_signed_work(self):
        b=NormalJ2TangentOwner();advance_cone(b,1.);advance_cone(b,2.)
        self.close([r['t'] for r in b.rows],[1.,7/6,2.])
        self.close(b.S,[np.diag([-2.,1.,1.])]);self.close(b.ep,[1/750],atol=1e-12)
        self.close(b.normal_force,[5.]);self.close(b.tangent_force,[[.6,.8]])
        self.close(b.u,[.0006,.0008]);self.close(b.wall(2.)[1:],[.0024,.0032])
        self.close([b.plastic,b.heat,b.normal_work,b.wall_work,b.material_work],
            [.004,.0025,.0095,.0125,.01])
        self.close(b.wall_work-b.normal_work,.003)
        self.close(b.returned,13/4500)
        self.close(.5*1000.*(b.u@b.u),.0005)
        self.close(b.wall_work-(.0055+.0005+b.plastic+b.heat),0.)
        self.assertGreater(b.rows[-1]['q'],0.);self.assertGreater(b.rows[-1]['plastic'],0.)
        self.assertEqual([r['kind'] for r in b.rows],['REQUESTED','CONE_LIMIT','REQUESTED'])
        self.assertFalse(b.rows[1]['released'][0]);self.assertTrue(b.rows[-1]['released'][0])

    def test_03_normal_unloading_drives_slip_at_stationary_tangent_wall(self):
        b=NormalJ2TangentOwner();advance_cone(b,1.);advance_cone(b,2.)
        old=np.array([b.plastic,b.heat,b.normal_work,b.wall_work,b.returned]);advance_cone(b,3.)
        self.close(b.S,np.zeros((1,3,3)));self.close(b.ep,[1/750],atol=1e-12)
        self.close(b.normal_force,[2.]);self.close(b.tangent_force,[[.24,.32]])
        self.close(b.u,[.00024,.00032])
        self.close(np.array([b.plastic,b.heat,b.normal_work,b.wall_work,b.returned])-old,
            [0.,.00042,-.0035,-.0035,0.])
        self.close([b.wall_work,b.plastic,b.heat],[.009,.004,.00292])
        self.close(.5*1000.*(b.u@b.u),.00008)
        self.close(b.wall_work-(.002+.00008+b.plastic+b.heat),0.)

    def test_04_pure_repeated_normal_and_contact_candidates_do_not_accumulate(self):
        b=NormalJ2TangentOwner();advance_cone(b,1.);before=b.snapshot()
        first=b.branch(2.,np.array([True]))
        self.assertGreater(first['plastic'][0],0.)
        for t in (1.25,1.9,1.125,2.):
            b.branch(t,np.array([True]));b.branch(t,np.array([False]));b.normal_trial(t)
        for k,v in b.snapshot().items():np.testing.assert_array_equal(v,before[k],err_msg=k)
        second=b.branch(2.,np.array([True]))
        for key in ('S','ep','plastic','normal_force','u','tangent_force'):
            np.testing.assert_array_equal(first[key],second[key],err_msg=key)
        advance_cone(b,2.);self.close([b.plastic,b.heat],[.004,.0025])

    def test_05_basis_representation_does_not_change_coupled_reference(self):
        for angle in (0.,37.):
            b=NormalJ2TangentOwner(basis_angle=angle);advance_cone(b,1.);advance_cone(b,2.);advance_cone(b,3.)
            self.close(b.u,[.00024,.00032]);self.close(b.tangent_force,[[.24,.32]])
            self.close([b.wall_work,b.plastic,b.heat],[.009,.004,.00292])
            self.close([r['t'] for r in b.rows],[1.,7/6,2.,3.])

    def test_06_numerical_root_scaling_not_physical_stiffness_or_modulus_change(self):
        for kappa in (100.,1000.,10000.):
            b=NormalJ2TangentOwner(kappa=kappa);advance_cone(b,1.);advance_cone(b,2.)
            self.close(b.u,[.0006,.0008]);self.close(b.normal_force,[5.])
            self.close([b.wall_work,b.plastic,b.heat],[.0125,.004,.0025])

    def test_07_return_diagnostic_changes_with_partition_not_physical_heat(self):
        coarse=NormalJ2TangentOwner();fine=NormalJ2TangentOwner()
        for b in (coarse,fine):advance_cone(b,1.)
        advance_cone(coarse,2.)
        for t in (7/6,1.5,2.):advance_cone(fine,t)
        self.assertLess(fine.returned,coarse.returned)
        self.close([fine.wall_work,fine.plastic,fine.heat],[.0125,.004,.0025])
        self.close([fine.u@A,fine.ep[0]],[.001,1/750])


def actual_snapshot(b):
    particles={name:(p.get_npy_array().dtype.str,p.get_npy_array().shape,p.get_npy_array().tobytes())
        for name,p in b.pa.properties.items()}
    arrays={k:getattr(b,k).copy() for k in ('displacement','stress','normal_force','tangent_force',
        'contact_force','support_force','friction','last_plastic','last_return','rate')}
    scalars={k:copy.deepcopy(getattr(b,k)) for k in ('time','count','wall_work','contact_work',
        'support_work','internal_work','events')}
    return particles,arrays,scalars


class ActualPlasticFrictionQ7(unittest.TestCase):
    def unchanged(self,b,before):
        after=actual_snapshot(b)
        self.assertEqual(after[0],before[0]);self.assertEqual(after[2],before[2])
        for k in after[1]:np.testing.assert_array_equal(after[1][k],before[1][k],err_msg=k)

    def test_08_actual_plastic_affine_trial_and_tangent_are_pure(self):
        b=PlasticFrictionBlock();before=actual_snapshot(b)
        self.assertEqual(len(b.X),27);self.assertEqual(int(b.free.sum()),54)
        self.assertEqual(b.law.value,90e6);self.assertEqual(b.young,70e9)
        u=np.zeros_like(b.X);u[:,0]=-.002*b.X[:,0]
        first=b.trial(u,.025);self.assertTrue(np.all(first['ep']>0.))
        # Independent uniaxial-STRAIN, not the Q2 uniaxial-stress oracle.
        np.testing.assert_allclose(first['S'],np.broadcast_to(np.diag([-60e6,30e6,30e6]),(27,3,3)),rtol=0,atol=2e-5)
        np.testing.assert_allclose(first['ep'],50e6/(3*35e9),rtol=0,atol=2e-15)
        b.tangent(first);second=b.trial(u,.025)
        for k in first:np.testing.assert_array_equal(first[k],second[k],err_msg=k)
        self.unchanged(b,before)

    def test_09_actual_plastic_commit_recomputes_untrusted_fields_once(self):
        b=PlasticFrictionBlock();commit=b.commit;count_positive=[];saved=[]
        def corrupt_then_commit(c):
            expected=b.state(b.trial(c['displacement'],c['time_s']),c['released'])
            P=b.pa.plastic_work+expected['plastic_increment_J']
            R=b.pa.return_dissipation+expected['return_increment_J']
            count_positive.append(float(expected['plastic_increment_J'].sum())>0.)
            for name in ('stress','S','ep','plastic_increment_J','return_increment_J',
                         'normal_force','tangent_force','law_residual_N'):
                c[name][:]=123456789.
            commit(c)
            for name,actual in [('stress',b.stress),('S',b.pa.S.reshape(-1,3,3)),('ep',b.pa.ep)]:
                np.testing.assert_array_equal(actual,expected[name],err_msg=name)
            np.testing.assert_array_equal(b.pa.plastic_work,P)
            np.testing.assert_array_equal(b.pa.return_dissipation,R)
            saved.append(c)
        b.commit=corrupt_then_commit;b.advance(.025)
        self.assertTrue(any(count_positive));self.assertGreater(float(b.pa.ep.max()),0.)
        self.assertGreater(float(b.pa.plastic_work.sum()),0.);self.assertGreater(float(b.friction.sum()),0.)
        before=actual_snapshot(b)
        with self.assertRaisesRegex(ValueError,'stale'):commit(saved[-1])
        self.unchanged(b,before)

    def test_10_tags_after_plastic_loading_leave_future_coupled_candidate_identical(self):
        b=PlasticFrictionBlock();b.advance(.025)
        self.assertGreater(float(b.pa.ep.max()),0.)
        first=solve_cone_modes(lambda m:b.branch(.026,m),len(b.nodes),b.mu,b.limits)
        ids=b.pa.particle_id.copy();b.pa.tracer_type[:]=71+np.arange(27)%5
        before=actual_snapshot(b)
        second=solve_cone_modes(lambda m:b.branch(.026,m),len(b.nodes),b.mu,b.limits)
        for name in ('displacement','stress','S','ep','plastic_increment_J','return_increment_J',
                     'normal_force','tangent_force','contact_force','law_residual_N'):
            np.testing.assert_array_equal(first[name],second[name],err_msg=name)
        np.testing.assert_array_equal(b.pa.particle_id,ids);self.unchanged(b,before)


if __name__=='__main__':unittest.main(verbosity=2)

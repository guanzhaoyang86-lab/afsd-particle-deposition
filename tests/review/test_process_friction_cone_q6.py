"""Literal, independently derived vector-disk references; no AFSD certification."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from probe_process_friction_cone import (ConeBlock,CandidateFailure,cone_response,
    cone_rows,solve_cone_modes,advance_cone,vector_contact_work,feasible_cone)

E=np.array([.6,.8])
LIMITS=dict(force_tolerance_N=1e-12,slip_tolerance_m=1e-12,
    maximum_active_sets=16,maximum_events_per_interval=12,event_bisections=48,
    event_time_tolerance_s=1e-12,mode_probe_fraction=1e-4)


class ElasticVectorOwner:
    """Shared linear equilibrium fixture; tested cone rows determine actual root.

    References in assertions are hand-derived literals, never generated here.
    The first 2*count body DOFs belong to contacts, remaining body DOFs are free.
    Dummy held normals give lambda=2; coupled=True uses an actual 3DOF stiffness.
    """
    def __init__(self,K,count,wall,base=None,kappa=1.,basis_angle=0.,coupled=False):
        K=np.asarray(K,dtype=float); self.body_size=len(K); self.coupled=coupled
        self.nodes=np.arange(count); self.wall=wall; self.mu=.5; self.kappa=kappa
        self.time=0.; self.count=0; self.limits=dict(LIMITS)
        self.base=base or (lambda t:np.zeros(self.body_size))
        if coupled:
            self.K=K.copy(); self.nd=np.array([0]); self.td=np.array([[1,2]])
            self.u=np.array([-1.,0.,0.])
        else:
            n=len(K); self.K=np.zeros((n+count,n+count)); self.K[:n,:n]=K
            self.nd=np.arange(n,n+count); self.K[self.nd,self.nd]=1.
            self.td=np.arange(2*count).reshape(count,2)
            self.u=np.zeros(n+count); self.u[self.nd]=-2.
        a=np.deg2rad(basis_angle); self.basis=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        force=self.K@(self.u-self.fullbase(0.))
        self.normal_force=-force[self.nd]; self.tangent_force=force[self.td]
        self.heat=self.wall_work=self.material_work=self.base_work=0.; self.rows=[]

    def fullbase(self,t):
        b=np.zeros(len(self.K)); b[:self.body_size]=self.base(t)
        return b

    def branch(self,t,released):
        released=np.asarray(released,dtype=bool); unknown=np.ones(len(self.K),dtype=bool)
        unknown[self.nd]=False; unknown[self.td[~released].ravel()]=False
        dw=np.asarray(self.wall(t))-self.wall(self.time); u=self.u.copy()
        u[self.td[~released]]=self.u[self.td[~released]]+dw[~released]
        def state(v):
            f=self.K@(v-self.fullbase(t)); tau=f[self.td]; normal=-f[self.nd]
            slip=v[self.td]-self.u[self.td]-dw
            response=cone_response(tau,normal,slip,self.mu,self.kappa,self.basis)
            residual=-f; residual[self.td[released].ravel()]=response['residual'][released].ravel()
            return dict(base_step=self.count,time_s=t,u=v.copy(),normal_force=normal,
                tangent_force=tau,relative_increment=slip,response=response),residual
        for _ in range(32):
            c,r=state(u); norm=np.linalg.norm(r[unknown])
            if norm<1e-13:return c
            M,r=cone_rows(self.K,r,self.nd,self.td,c['response'],self.mu,self.kappa,released)
            step=np.zeros(len(u)); step[unknown]=np.linalg.solve(M[np.ix_(unknown,unknown)],r[unknown])
            for half in range(12):
                trial=u+step*2.**(-half)
                if np.linalg.norm(state(trial)[1][unknown])<norm:
                    u=trial;break
            else:raise RuntimeError('independent fixture Newton line search')
        raise RuntimeError('independent fixture Newton budget')

    def commit(self,c):
        mid=.5*(self.tangent_force+c['tangent_force']); du=c['u']-self.u
        dw=np.asarray(self.wall(c['time_s']))-self.wall(self.time)
        q=-np.sum(mid*(du[self.td]-dw),axis=1)
        oldf=self.K@(self.u-self.fullbase(self.time))
        newf=self.K@(c['u']-self.fullbase(c['time_s']))
        self.base_work-=float(.5*(oldf+newf)@(self.fullbase(c['time_s'])-self.fullbase(self.time)))
        self.wall_work+=float(np.sum(mid*dw)); self.material_work+=float(np.sum(mid*du[self.td]))
        self.heat+=float(q.sum())
        self.rows.append(dict(t=c['time_s'],w=np.asarray(self.wall(c['time_s'])),u=c['u'].copy(),
            tau=c['tangent_force'].copy(),q=q,kind=c['event_kind'],released=c['released'].copy()))
        self.u=c['u'].copy();self.time=c['time_s'];self.count+=1
        self.normal_force=c['normal_force'].copy();self.tangent_force=c['tangent_force'].copy()


def snapshot(b):
    return dict(u=b.displacement.copy(),stress=b.stress.copy(),S=b.pa.S.copy(),ep=b.pa.ep.copy(),
        xyz=np.array([b.pa.x,b.pa.y,b.pa.z]),force=b.contact_force.copy(),q=b.friction.copy(),
        tau=b.tangent_force.copy(),normal=b.normal_force.copy(),time=b.time,count=b.count,
        works=np.array([b.wall_work,b.contact_work,b.support_work,b.internal_work]))


class IndependentConeReferences(unittest.TestCase):
    def close(self,a,b,atol=1e-10):np.testing.assert_allclose(a,b,atol=atol,rtol=0)
    def unchanged(self,b,before):
        for k,v in snapshot(b).items():np.testing.assert_array_equal(v,before[k],err_msg=k)

    def test_01_oblique_reversal_events_and_signed_work(self):
        b=ElasticVectorOwner(np.eye(2),1,lambda t:(E*(t if t<=3 else 6-t))[None,:])
        advance_cone(b,3.);advance_cone(b,6.)
        self.close([r['w'][0] for r in b.rows],np.array([1,3,1,0])[:,None]*E)
        self.close([r['u'][:2] for r in b.rows],np.array([1,1,-1,-1])[:,None]*E)
        self.close([r['tau'][0] for r in b.rows],np.array([1,1,-1,-1])[:,None]*E)
        self.close([r['q'][0] for r in b.rows],[0,2,0,1])
        np.testing.assert_array_equal([r['released'][0] for r in b.rows],[False,True,False,True])
        self.close([b.heat,b.wall_work,b.material_work],[3,3.5,.5])

    def test_02_two_reversals_and_intermediate_stick_not_wall_direction(self):
        def wall(t):return ((t if t<=3 else 6-t if t<=9 else t-12)*E)[None,:]
        b=ElasticVectorOwner(np.eye(2),1,wall)
        advance_cone(b,3.);advance_cone(b,4.)
        self.close(b.u[:2],[0,0]);self.close(b.tangent_force,[[0,0]])
        self.assertFalse(b.rows[-1]['released'][0])
        advance_cone(b,9.);advance_cone(b,12.)
        self.close(b.u[:2],[.6,.8]);self.close([b.heat,b.wall_work,b.material_work],[7,7.5,.5])

    def test_03_shared_free_two_vector_dofs_not_diagonal_contacts(self):
        K=np.kron([[2,0,-1],[0,2,-1],[-1,-1,3]],np.eye(2))
        wall=lambda t:np.tile(.75*(t if t<=3 else 6-t)*E,(2,1))
        b=ElasticVectorOwner(K,2,wall);advance_cone(b,3.);advance_cone(b,6.)
        self.close([r['u'][:6] for r in b.rows],
            [np.r_[.75*E,.75*E,.5*E],np.r_[.75*E,.75*E,.5*E],
             np.r_[-.75*E,-.75*E,-.5*E],np.r_[-.75*E,-.75*E,-.5*E]])
        self.close([b.heat,b.wall_work,b.material_work],[4.5,5.25,.75])
        self.close(.5*b.u[:6]@K@b.u[:6],.75)

    def test_04_capacity_uses_current_coupled_normal_reaction(self):
        K=[[2,.3,.4],[.3,1,0],[.4,0,1]]
        for sign,q,normal,tau in [(1,1.2,1.4,[.42,.56]),(-1,-2/3,7/3,[-.7,-14/15])]:
            b=ElasticVectorOwner(K,1,lambda t:(sign*t*E)[None,:],coupled=True)
            c=solve_cone_modes(lambda m:b.branch(2.,m),1,.5,LIMITS)
            self.close(c['u'][1:],q*E);self.close(c['normal_force'],[normal]);self.close(c['tangent_force'],[tau])

    def test_05_anisotropic_noncollinear_endpoint_reference(self):
        b=ElasticVectorOwner(np.diag([2.,1.]),1,lambda t:np.array([[1.5*t,2.4*t]]))
        c=solve_cone_modes(lambda m:b.branch(1.,m),1,.5,LIMITS)
        self.close(c['u'][:2],[.3,.8]);self.close(c['tangent_force'],[[.6,.8]])
        self.close(c['relative_increment'],[[-1.2,-1.6]])

    def test_06_coordinate_transform_and_internal_basis_leave_root_unchanged(self):
        for angle in (0.,37.,90.):
            b=ElasticVectorOwner(np.diag([2.,1.]),1,lambda t:np.array([[1.5*t,2.4*t]]),basis_angle=angle)
            c=solve_cone_modes(lambda m:b.branch(1.,m),1,.5,LIMITS);self.close(c['u'][:2],[.3,.8])
        b=ElasticVectorOwner(np.diag([1.,2.]),1,lambda t:np.array([[-2.4*t,1.5*t]]))
        c=solve_cone_modes(lambda m:b.branch(1.,m),1,.5,LIMITS)
        self.close(c['u'][:2],[-.8,.3]);self.close(c['tangent_force'],[[-.8,.6]])

    def test_07_numerical_kappa_does_not_change_physical_root(self):
        for kappa in (.1,1.,10.):
            b=ElasticVectorOwner(np.diag([2.,1.]),1,lambda t:np.array([[1.5*t,2.4*t]]),kappa=kappa)
            c=solve_cone_modes(lambda m:b.branch(1.,m),1,.5,LIMITS)
            self.close(c['u'][:2],[.3,.8]);self.close(c['tangent_force'],[[.6,.8]])

    def test_08_static_wall_and_common_translation_work(self):
        b=ElasticVectorOwner(np.eye(2),1,lambda t:np.zeros((1,2)),lambda t:t*E)
        advance_cone(b,3.);self.close(b.u[:2],[1.2,1.6])
        self.close([b.wall_work,b.material_work,b.heat,b.base_work],[0,-2,2,2.5])
        b=ElasticVectorOwner(np.eye(2),1,lambda t:(t*E)[None,:],lambda t:(t-.5)*E)
        advance_cone(b,2.)
        self.close([b.wall_work,b.material_work,b.heat,b.base_work],[1,1,0,-1])
        self.close(b.tangent_force,[[.3,.4]]);self.assertFalse(b.rows[-1]['released'][0])

    def test_09_circle_not_component_box_and_independent_projection_derivative(self):
        c=dict(normal_force=np.array([2.]),tangent_force=np.array([[.8,.8]]),relative_increment=np.zeros((1,2)))
        self.assertFalse(feasible_cone(c,np.array([False]),.5,LIMITS))
        c['tangent_force'][:]=[.6,.8];self.assertTrue(feasible_cone(c,np.array([False]),.5,LIMITS))
        r=cone_response([[3,4]],[2],[[0,0]],1.,1.)
        self.close(r['residual'],[[-1.8,-2.4]])
        self.close(r['Pq'],[[[.256,-.192],[-.192,.144]]]);self.close(r['Pc'],[[.6,.8]])
        inside=cone_response([[.1,.2]],[2],[[0,0]],1.,1.)
        self.close(inside['Pq'],[np.eye(2)]);self.close(inside['Pc'],[[0,0]])
        zero=cone_response([[0,0]],[0],[[0,0]],1.,1.)
        self.assertTrue(np.isfinite(zero['Pq']).all());self.close(zero['residual'],[[0,0]])

    def test_10_vector_work_retains_negative_raw_term_without_cancellation(self):
        force=np.array([[0,.6,.8],[0,-.6,-.8]])
        du=np.array([[0,.6,.8],[0,.3,.4]])
        w=vector_contact_work(force,force,du,np.zeros(3),np.array([0,1]),.25)
        self.close(w['friction_increment_J'],[-1,.5]);self.close(w['tangent_impulse_N_s'],[[.15,.2],[-.15,-.2]])
        self.close(w['relative_kinematic_rate_m_s'],[[2.4,3.2],[1.2,1.6]])
        w=vector_contact_work(force[:1],force[:1],du[:1],np.array([0,1.8,2.4]),np.array([0]),2.)
        self.close([w['wall_work_J'],w['material_work_J'],w['friction_increment_J'][0]],[3,1,2])

    def test_11_invalid_inputs_and_normal_detachment_do_not_clip(self):
        for tau,normal,slip in [([[np.nan,0]],[1],[[0,0]]),([[0,0]],[-1],[[0,0]]),([[0,0]],[1],[[0,np.inf]])]:
            with self.assertRaises(ValueError):cone_response(tau,normal,slip,.5,1.)
        with self.assertRaisesRegex(ValueError,'orthonormal'):cone_response([[0,0]],[1],[[0,0]],.5,1.,np.ones((2,2)))
        with self.assertRaisesRegex(ValueError,'positive'):vector_contact_work(np.zeros((1,3)),np.zeros((1,3)),np.zeros((1,3)),np.zeros(3),np.array([0]),0.)

    def test_12_near_start_event_and_opposite_boundary_probe_are_not_skipped(self):
        b=ElasticVectorOwner(np.eye(2),1,lambda t:(t*E)[None,:])
        advance_cone(b,1.-1e-6);advance_cone(b,3.)
        self.close([r['w'][0]@E for r in b.rows],[1.-1e-6,1.,3.])
        self.close([b.heat,b.wall_work],[2,2.5])
        b=ElasticVectorOwner(np.eye(2),1,lambda t:((t if t<=3 else 6-t)*E)[None,:])
        advance_cone(b,3.);b.limits['mode_probe_fraction']=.9
        advance_cone(b,6.)
        self.close([r['w'][0]@E for r in b.rows],[1,3,1,0]);self.close(b.heat,3.)

    def test_13_actual_body_inventory_trial_purity_and_passive_tags(self):
        b=ConeBlock();self.assertEqual(len(b.X),27);self.assertEqual(len(b.nodes),9)
        self.assertEqual(b.tangent_dofs.shape,(9,2));self.assertEqual(int(b.free.sum()),54)
        before=snapshot(b);c=b.branch(.005,np.zeros(9,dtype=bool));self.unchanged(b,before)
        b.pa.tracer_type[:]=71+np.arange(27)%5
        other=b.branch(.005,np.zeros(9,dtype=bool))
        for key in ('displacement','stress','normal_force','tangent_force','contact_force'):
            np.testing.assert_array_equal(c[key],other[key],err_msg=key)
        self.unchanged(b,before)

    def test_14_actual_failed_metadata_and_program_corner_are_atomic(self):
        b=ConeBlock();before=snapshot(b)
        c=b.branch(.005,np.zeros(9,dtype=bool));c['event_kind']='CONE_LIMIT'
        with self.assertRaisesRegex(ValueError,'bracket'):b.commit(c)
        self.unchanged(b,before)
        with self.assertRaisesRegex(ValueError,'program corner'):b.advance(.075)
        self.unchanged(b,before)

    def test_15_actual_commit_recomputes_fields_and_rejects_stale_candidate(self):
        b=ConeBlock();c=b.branch(.005,np.zeros(9,dtype=bool));expected=c['stress'].copy()
        c['stress'][:]=123.;c['normal_force'][:]=-100.;c['tangent_force'][:]=999.;c['ep'][:]=5.
        b.commit(c);self.close(b.stress,expected,atol=1e-8)
        np.testing.assert_array_equal(b.pa.ep,np.zeros(27));self.assertGreater(b.normal_force.min(),0.)
        with self.assertRaisesRegex(ValueError,'stale'):b.commit(c)


if __name__=='__main__':unittest.main(verbosity=2)

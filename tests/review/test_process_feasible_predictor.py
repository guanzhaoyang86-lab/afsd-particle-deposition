"""Independent SPD elimination reference; no solver-produced reference values."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from probe_process_friction_cone import ConeBlock, feasible_cone
from probe_process_plastic import CandidateFailure
from probe_process_feasible_predictor import solve_with_feasible_predictor


class LinearFixture:
    """SI two-node vector embedding of (x,s,q); remaining DOFs are uncoupled.

    s is prescribed rigid-wall tangent travel, NOT a material trajectory.
    q is solved, initial prestress retained. Pure linear reference, not AFSD.
    """
    state = ConeBlock.state
    branch = ConeBlock.branch

    def __init__(self):
        self.X = np.zeros((2,3))
        self.nodes = np.array([0])
        self.contact_dofs = np.array([0])
        self.tangent_dofs = np.array([[1,2]])
        self.fixed = np.array([[False,False,False],[False,True,True]])
        self.free = ~self.fixed.ravel()
        self.K = np.eye(6)
        indices = [0,1,3]
        self.K[np.ix_(indices,indices)] = [[2,1,-1],[1,3,-1],[-1,-1,2]]
        self.displacement = np.array([[-1.,0.,0.],[-.5,0.,0.]])
        self.time = 0.
        self.count = 0
        self.mu = 100.
        self.kappa = 1.
        self.basis = np.eye(2)
        self.force_scale = 1.
        self.ccfg = dict(numerical=dict(force_relative=1e-12,tangent_condition_limit=1e12))
        self.qcfg = dict(numerical=dict(maximum_iterations=14,line_search_halvings=8))
        self.limits = dict(force_tolerance_N=1e-11,slip_tolerance_m=1e-12)
        self.normal_force = np.array([1.5])
        self.tangent_force = np.array([[-.5,0.]])

    def wall(self,t):
        return np.array([-1.,t,0.])

    def trial(self,u,t):
        if np.shape(u)!=(2,3) or not np.isfinite(u).all():
            raise ValueError('invalid linear displacement')
        if np.any(u[self.fixed] != 0):
            raise ValueError('support changed')
        force = -(self.K@u.ravel()).reshape(2,3)
        return dict(displacement=u.copy(),internal_force=force,residual=force.copy(),
                    time_s=t,dt_s=t-self.time,base_step=self.count)

    def tangent(self,c):
        return self.K.copy()

    def commit(self,c):
        raise AssertionError('pure fixture cannot commit')


def snapshot(b):
    return (b.displacement.tobytes(), b.normal_force.tobytes(),
            b.tangent_force.tobytes(), b.time, b.count)


class LegacyNegativeSeedRed(unittest.TestCase):
    @unittest.expectedFailure
    def test_original_guard_cannot_reach_existing_positive_normal_equilibrium(self):
        b = LinearFixture()
        c = b.branch(2.5,np.array([False]))
        np.testing.assert_allclose(c['normal_force'],[.25],atol=1e-11,rtol=0)


class FeasiblePredictorTests(unittest.TestCase):
    def close(self,a,b):
        np.testing.assert_allclose(a,b,atol=1e-11,rtol=0)

    def test_false_opening_initial_guess_reaches_independent_positive_root(self):
        b = LinearFixture(); before = snapshot(b); raw=[]
        c = solve_with_feasible_predictor(b,2.5,np.array([False]),
                                         observe=lambda stage,c:raw.append((stage,c)))
        self.close(raw[0][1]['internal_force'][0,0],-1.)
        self.close(c['displacement'],[[-1.,2.5,0.],[.75,0.,0.]])
        self.close(c['normal_force'],[.25])
        self.close(c['tangent_force'],[[5.75,0.]])
        self.close(c['residual'],np.zeros((2,3)))
        self.assertTrue(feasible_cone(c,np.array([False]),b.mu,b.limits))
        self.assertEqual(c['predictor_corrections'],1)
        self.assertEqual(c['cone_corrections'],0)
        self.assertEqual(c['iterations'],1)
        self.assertEqual(snapshot(b),before)

    def test_signed_tensile_equilibrium_is_blocked_not_clipped(self):
        b = LinearFixture(); before = snapshot(b)
        with self.assertRaisesRegex(CandidateFailure,'negative normal at interior equilibrium') as caught:
            solve_with_feasible_predictor(b,4.,np.array([False]))
        c = caught.exception.candidate
        self.close(c['displacement'],[[-1.,4.,0.],[1.5,0.,0.]])
        self.close(c['internal_force'][0,0],-.5)
        self.close(c['internal_force'][1,0],0.)
        self.assertEqual(snapshot(b),before)

    def test_independent_linear_work_no_predictor_energy_turned_into_heat(self):
        b=LinearFixture(); c=solve_with_feasible_predictor(b,2.5,np.array([False]))
        u0=b.displacement.ravel(); u=c['displacement'].ravel()
        E0=.5*u0@b.K@u0; E=.5*u@b.K@u
        # External forces = K u. q has zero force at BOTH equilibria.
        work=.5*(b.K@u0+b.K@u)@(u-u0)
        self.close(E0,.75); self.close(E,117/16); self.close(work,105/16)
        self.close(E-E0,work)
        self.assertNotIn('heat',c)
        self.assertEqual(b.time,0.)

    def test_positive_seed_is_exact_old_branch_without_predictor(self):
        b=LinearFixture()
        old=b.branch(.5,np.array([False]))
        new=solve_with_feasible_predictor(b,.5,np.array([False]))
        for key in ('displacement','normal_force','tangent_force','residual','law_residual_N'):
            np.testing.assert_array_equal(old[key],new[key])
        self.assertEqual(new['predictor_corrections'],0)
        self.assertEqual(old['iterations'],new['iterations'])

    def test_zero_budget_and_invalid_budgets_fail_without_mutation(self):
        b=LinearFixture(); before=snapshot(b)
        with self.assertRaisesRegex(CandidateFailure,'shared correction budget'):
            solve_with_feasible_predictor(b,2.5,np.array([False]),0)
        for value in (-1,15,True,1.5):
            with self.subTest(value=value),self.assertRaises(ValueError):
                solve_with_feasible_predictor(b,2.5,np.array([False]),value)
        self.assertEqual(snapshot(b),before)

    def test_no_projection_to_zero_normal_accepts_free_residual(self):
        b=LinearFixture()
        raw=b.trial(np.array([[-1.,4.,0.],[2.,0.,0.]]),4.)
        self.close(raw['internal_force'][0,0],0.)
        self.close(raw['internal_force'][1,0],-1.)
        # Starting from this N=0 candidate must NOT be accepted without equilibrium.
        b.displacement=raw['displacement'].copy()
        before=snapshot(b)
        with self.assertRaises(CandidateFailure):
            solve_with_feasible_predictor(b,4.,np.array([False]))
        self.assertEqual(snapshot(b),before)

    def test_all_contacts_held_only_in_predictor_then_released_coupled(self):
        b=LinearFixture()
        # Force the initial raw seed negative even when tangent row is released.
        # x=-1,s=2.5,q=-.5 => N=-1. Accepted-like fixture solely for seeding.
        b.displacement[0,1]=2.5
        before=snapshot(b); raw=[]
        c=solve_with_feasible_predictor(b,2.5,np.array([True]),
                          observe=lambda stage,c:raw.append((stage,c)))
        stages=[u for stage,u in raw if stage=='PREDICTOR_LINE_SEARCH']
        self.assertTrue(stages)
        self.close(stages[0]['displacement'][0],[-1.,2.5,0.])
        self.assertGreater(c['cone_corrections'],0)
        self.assertLessEqual(c['iterations'],14)
        self.assertEqual(c['iterations'],c['predictor_corrections']+c['cone_corrections'])
        self.assertEqual(snapshot(b),before)

    def test_predictor_consumes_shared_budget_before_cone(self):
        b=LinearFixture(); b.displacement[0,1]=2.5
        before=snapshot(b)
        with self.assertRaisesRegex(CandidateFailure,'shared correction budget exhausted in cone'):
            solve_with_feasible_predictor(b,2.5,np.array([True]),1)
        self.assertEqual(snapshot(b),before)

    def test_singular_interior_rejected_without_pseudoinverse(self):
        b=LinearFixture(); before=snapshot(b)
        # Deliberately corrupted solver tangent; force oracle stays independent.
        original=b.tangent
        def singular(c):
            k=original(c); k[3,3]=0.; return k
        b.tangent=singular
        with self.assertRaisesRegex(CandidateFailure,'ill-conditioned'):
            solve_with_feasible_predictor(b,2.5,np.array([False]))
        self.assertEqual(snapshot(b),before)

    def test_original_negative_normal_guard_unchanged(self):
        b=LinearFixture()
        raw=b.trial(np.array([[-1.,2.5,0.],[-.5,0.,0.]]),2.5)
        with self.assertRaisesRegex(ValueError,'negative normal'):
            b.state(raw,np.array([False]))

    def test_invalid_time_and_release_inventory_rejected(self):
        b=LinearFixture(); before=snapshot(b)
        for t in (0.,-1.,float('nan'),float('inf')):
            with self.subTest(t=t),self.assertRaises(ValueError):
                solve_with_feasible_predictor(b,t,np.array([False]))
        for mask in (np.array([0]),np.array([False,True])):
            with self.assertRaises(ValueError):
                solve_with_feasible_predictor(b,2.5,mask)
        self.assertEqual(snapshot(b),before)

    def test_normal_boundary_sides_do_not_replace_cone_feasibility(self):
        b=LinearFixture()
        eps=2.**-20
        c=solve_with_feasible_predictor(b,3.-eps,np.array([False]))
        self.close(c['normal_force'],[2.**-21])
        # A fixed-stick branch root can be outside the friction disk; NOT accepted.
        self.assertFalse(feasible_cone(c,np.array([False]),b.mu,b.limits))
        with self.assertRaises(CandidateFailure) as caught:
            solve_with_feasible_predictor(b,3.+eps,np.array([False]))
        self.close(caught.exception.candidate['internal_force'][0,0],-2.**-21)

    def test_nonfinite_raw_force_explicitly_rejected(self):
        b=LinearFixture(); original=b.trial; before=snapshot(b)
        def nonfinite(u,t):
            c=original(u,t); c['internal_force'][1,0]=np.nan; return c
        b.trial=nonfinite
        with self.assertRaisesRegex(CandidateFailure,'nonfinite'):
            solve_with_feasible_predictor(b,2.5,np.array([False]))
        self.assertEqual(snapshot(b),before)


if __name__=='__main__':
    unittest.main(verbosity=2)

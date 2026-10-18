"""Independent literal quasi-static references; not the old inertial AFSD test."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from probe_process_friction import (FrictionBlock, CandidateFailure,
    advance_contact, solve_friction_modes, contact_work)


LIMITS = dict(force_tolerance_N=1e-12, slip_tolerance_m=1e-12,
    maximum_active_sets=16, maximum_events_per_interval=12,
    event_bisections=48, event_time_tolerance_s=1e-12, mode_probe_fraction=1e-4)


class ElasticSpringOwner:
    """Pure linear branch equilibrium; production helper selects contact modes.

    Contact DOFs are the first count components. All remaining DOFs are free.
    Literal expected answers in tests come from the separate piecewise analysis,
    not this fixture or the tested solver. No per-contact diagonal approximation.
    """
    def __init__(self, K, count, wall, base=None):
        self.K = np.asarray(K, dtype=float)
        self.nodes = np.arange(count)
        self.wall = wall
        self.base = base or (lambda t: np.zeros(len(self.K)))
        self.mu = .5
        self.limits = dict(LIMITS)
        self.time = 0.
        self.count = 0
        self.u = np.zeros(len(self.K))
        self.force = (self.K @ (self.u-self.base(0.)))[:count]
        self.normal_force = np.full(count,2.)
        self.tangent_force = self.force.copy()
        self.wall_work = self.material_work = self.base_work = self.heat = 0.
        self.rows = []

    def branch(self, t, modes):
        modes = np.asarray(modes)
        n = len(modes)
        stick = modes == 0
        fixed = np.zeros(len(self.K), dtype=bool)
        fixed[:n] = stick
        u = self.u.copy()
        dw = np.asarray(self.wall(t)) - self.wall(self.time)
        u[:n][stick] = (self.u[:n] + dw)[stick]
        free = ~fixed
        imposed = np.zeros(len(self.K))
        imposed[:n][~stick] = -modes[~stick]
        rhs = imposed + self.K @ self.base(t)
        u[free] = np.linalg.solve(self.K[np.ix_(free, free)],
            rhs[free] - self.K[np.ix_(free, fixed)] @ u[fixed])
        force = (self.K @ (u-self.base(t)))[:n]
        return dict(base_step=self.count, time_s=t, u=u,
            normal_force=np.full(n, 2.), tangent_force=force,
            relative_increment=u[:n]-self.u[:n]-dw)

    def commit(self, c):
        # This fixture only stores the algorithm's actual accepted endpoints.
        n = len(self.nodes)
        mid = .5*(self.force + c['tangent_force'])
        du = c['u']-self.u
        dw = np.asarray(self.wall(c['time_s'])) - self.wall(self.time)
        q = -mid*(du[:n]-dw)
        old_spring_force = self.K @ (self.u-self.base(self.time))
        new_spring_force = self.K @ (c['u']-self.base(c['time_s']))
        self.base_work += float(-.5*(old_spring_force+new_spring_force) @
            (self.base(c['time_s'])-self.base(self.time)))
        self.wall_work += float(mid @ dw)
        self.material_work += float(mid @ du[:n])
        self.heat += float(q.sum())
        self.rows.append(dict(t=c['time_s'], w=np.array(self.wall(c['time_s'])),
            u=c['u'].copy(), force=c['tangent_force'].copy(), q=q.copy(),
            kind=c['event_kind'], modes=c['modes'].copy()))
        self.time = c['time_s']
        self.count += 1
        self.u = c['u'].copy()
        self.force = c['tangent_force'].copy()
        self.normal_force = c['normal_force'].copy()
        self.tangent_force = self.force.copy()


def snapshot(b):
    return dict(u=b.displacement.copy(), stress=b.stress.copy(), S=b.pa.S.copy(),
        ep=b.pa.ep.copy(), x=b.pa.x.copy(), y=b.pa.y.copy(), z=b.pa.z.copy(),
        q=b.friction.copy(), force=b.contact_force.copy(), time=b.time, count=b.count,
        works=np.array([b.wall_work,b.contact_work,b.support_work,b.internal_work]))


class IndependentFrictionReferences(unittest.TestCase):
    def assert_accepted_snapshot_equal(self, before, b):
        for key, value in snapshot(b).items():
            np.testing.assert_array_equal(value, before[key], err_msg=key)

    def test_01_scalar_reversal_events_and_signed_work(self):
        b = ElasticSpringOwner([[1.]], 1,
            lambda t: np.array([t if t <= 3. else 6.-t]))
        advance_contact(b, 3.)
        advance_contact(b, 6.)
        # Constants are hand-derived, not produced by the contact implementation.
        np.testing.assert_allclose([r['w'][0] for r in b.rows], [1.,3.,1.,0.], atol=1e-10, rtol=0)
        np.testing.assert_allclose([r['u'][0] for r in b.rows], [1.,1.,-1.,-1.], atol=1e-10, rtol=0)
        np.testing.assert_allclose([r['force'][0] for r in b.rows], [1.,1.,-1.,-1.], atol=1e-10, rtol=0)
        np.testing.assert_array_equal([r['modes'][0] for r in b.rows], [0,-1,0,1])
        np.testing.assert_allclose([r['q'][0] for r in b.rows], [0.,2.,0.,1.], atol=1e-10, rtol=0)
        np.testing.assert_allclose([b.heat,b.wall_work,b.material_work], [3.,3.5,.5], atol=1e-10, rtol=0)
        self.assertAlmostEqual(.5*b.u@b.K@b.u, .5, delta=1e-10)

    def test_02_two_reversals_preserve_correct_retained_displacement(self):
        def wall(t):
            return np.array([t if t<=3 else 6-t if t<=9 else t-12])
        b=ElasticSpringOwner([[1.]],1,wall)
        for t in (3.,9.,12.): advance_contact(b,t)
        np.testing.assert_allclose([r['w'][0] for r in b.rows], [1,3,1,-3,-1,0], atol=1e-10, rtol=0)
        np.testing.assert_allclose([r['u'][0] for r in b.rows], [1,1,-1,-1,1,1], atol=1e-10, rtol=0)
        np.testing.assert_array_equal([r['modes'][0] for r in b.rows], [0,-1,0,1,0,-1])
        np.testing.assert_allclose([b.heat,b.wall_work,b.material_work], [7,7.5,.5], atol=1e-10, rtol=0)

    def test_03_shared_free_dof_coupling_not_diagonal_spring(self):
        K = [[2,0,-1],[0,2,-1],[-1,-1,3]]
        b = ElasticSpringOwner(K,2,lambda t: np.full(2,.75*(t if t<=3 else 6-t)))
        advance_contact(b,3.); advance_contact(b,6.)
        np.testing.assert_allclose([r['w'][0] for r in b.rows],[.75,2.25,.75,0],atol=1e-10,rtol=0)
        np.testing.assert_allclose([r['u'] for r in b.rows],
            [[.75,.75,.5],[.75,.75,.5],[-.75,-.75,-.5],[-.75,-.75,-.5]],atol=1e-10,rtol=0)
        np.testing.assert_allclose([b.heat,b.wall_work,b.material_work],[4.5,5.25,.75],atol=1e-10,rtol=0)
        self.assertAlmostEqual(.5*b.u@b.K@b.u,.75,delta=1e-10)

    def test_04_coulomb_capacity_changes_with_coupled_normal_reaction(self):
        def branch(w, modes):
            mode = int(modes[0])
            u = w if mode==0 else (.5-mode)/(1-.25*mode)
            return dict(normal_force=np.array([2-.5*u]),tangent_force=np.array([-.5+u]),
                relative_increment=np.array([u-w]),u=u)
        forward=solve_friction_modes(lambda m:branch(2.,m),1,.5,LIMITS)
        reverse=solve_friction_modes(lambda m:branch(-2.,m),1,.5,LIMITS)
        np.testing.assert_allclose([forward['u'],forward['normal_force'][0],forward['tangent_force'][0]],
            [1.2,1.4,.7],atol=1e-10,rtol=0)
        np.testing.assert_allclose([reverse['u'],reverse['normal_force'][0],reverse['tangent_force'][0]],
            [-2/3,7/3,-7/6],atol=1e-10,rtol=0)
        np.testing.assert_array_equal(forward['modes'],[-1])
        np.testing.assert_array_equal(reverse['modes'],[1])

    def test_05_static_wall_and_driven_support_work(self):
        b=ElasticSpringOwner([[1.]],1,lambda t:np.array([0.]),lambda t:np.array([t]))
        advance_contact(b,3.)
        np.testing.assert_allclose([r['t'] for r in b.rows],[1.,3.],atol=1e-10,rtol=0)
        np.testing.assert_allclose([b.u[0],b.force[0],b.wall_work,b.material_work,b.heat,b.base_work],
            [2.,-1.,0.,-2.,2.,2.5],atol=1e-10,rtol=0)

    def test_06_common_translation_zero_relative_slip_not_zero_material_work(self):
        b=ElasticSpringOwner([[1.]],1,lambda t:np.array([t]),lambda t:np.array([t-.5]))
        advance_contact(b,2.)
        np.testing.assert_allclose([b.u[0],b.force[0],b.wall_work,b.material_work,b.heat,b.base_work],
            [2.,.5,1.,1.,0.,-1.],atol=1e-10,rtol=0)
        np.testing.assert_array_equal(b.rows[-1]['modes'],[0])

    def test_07_actual_contact_work_retains_raw_sign_and_moving_wall_identity(self):
        old=np.array([[0.,1.,0.],[0.,-2.,0.]])
        du=np.array([[0.,1.,0.],[0.,.5,0.]])
        # Contact0 has physically negative Q; helper must expose, not rectify it.
        c=contact_work(old,old,du,np.zeros(3),np.array([0,1]),.25)
        np.testing.assert_allclose(c['friction_increment_J'],[-1.,1.],atol=1e-10,rtol=0)
        np.testing.assert_allclose(c['tangent_impulse_N_s'],[.25,-.5],atol=1e-10,rtol=0)
        np.testing.assert_allclose(c['relative_kinematic_rate_m_s'],[4.,2.],atol=1e-10,rtol=0)
        moving=contact_work(old[:1],old[:1],np.array([[0.,1.,0.]]),np.array([0.,3.,0.]),np.array([0]),2.)
        np.testing.assert_allclose([moving['wall_work_J'],moving['material_work_J'],moving['friction_increment_J'][0]],
            [3.,1.,2.],atol=1e-10,rtol=0)

    def test_08_normal_detachment_blocks_and_mode_budget_is_bounded(self):
        def detached(mode):
            return dict(normal_force=np.array([-1.]),tangent_force=np.array([0.]),relative_increment=np.array([0.]))
        with self.assertRaisesRegex(CandidateFailure,'normal detachment'):
            solve_friction_modes(detached,1,.5,LIMITS)
        def cycle(mode):
            return dict(normal_force=np.array([2.]),tangent_force=np.array([2.]),relative_increment=np.array([1.]))
        with self.assertRaisesRegex(CandidateFailure,'bounded friction modes'):
            solve_friction_modes(cycle,1,.5,{**LIMITS,'maximum_active_sets':3})

    def test_09_actual_block_inventory_trial_purity_and_passive_labels(self):
        b=FrictionBlock()
        self.assertEqual(len(b.X),27); self.assertEqual(len(b.nodes),9)
        self.assertEqual(int(b.free.sum()),54)
        np.testing.assert_array_equal(b.fixed,np.broadcast_to((b.X[:,0]==0)[:,None],b.X.shape))
        before=snapshot(b); c=b.branch(.005,np.zeros(9,dtype=int))
        self.assert_accepted_snapshot_equal(before,b)
        original_labels=b.pa.tracer_type.copy()
        b.pa.tracer_type[:]=71+np.arange(27)%5
        self.assertFalse(np.array_equal(original_labels,b.pa.tracer_type))
        changed=b.branch(.005,np.zeros(9,dtype=int))
        for name in ('displacement','stress','normal_force','tangent_force','contact_force'):
            np.testing.assert_array_equal(changed[name],c[name],err_msg=name)
        self.assert_accepted_snapshot_equal(before,b)

    def test_10_invalid_event_metadata_rolls_back_actual_material(self):
        b=FrictionBlock(); before=snapshot(b)
        c=b.branch(.005,np.zeros(9,dtype=int)); c['event_kind']='FRICTION_LIMIT'
        with self.assertRaisesRegex(ValueError,'invalid friction event bracket'): b.commit(c)
        self.assert_accepted_snapshot_equal(before,b)

    def test_11_actual_commit_recomputes_untrusted_derived_fields(self):
        b=FrictionBlock(); c=b.branch(.005,np.zeros(9,dtype=int)); expected=c['stress'].copy()
        c['stress'][:]=123456.; c['tangent_force'][:]=999.; c['normal_force'][:]=-111.
        c['ep'][:]=9.; c['plastic_increment_J'][:]=1e9
        b.commit(c)
        np.testing.assert_allclose(b.stress,expected,atol=1e-8,rtol=0)
        np.testing.assert_array_equal(b.pa.ep,np.zeros(27))
        self.assertGreater(float(b.normal_force.min()),0.)
        self.assertEqual(b.count,1)
        with self.assertRaisesRegex(ValueError,'stale candidate'): b.commit(c)

    def test_12_program_corner_cannot_be_hidden_in_one_increment(self):
        b=FrictionBlock(); before=snapshot(b)
        with self.assertRaisesRegex(ValueError,'program corner'): b.advance(.075)
        self.assert_accepted_snapshot_equal(before,b)

    def test_13_mode_probe_must_not_jump_a_near_start_yield_event(self):
        b=ElasticSpringOwner([[1.]],1,lambda t:np.array([t]))
        advance_contact(b,1.-1e-6)
        advance_contact(b,3.)
        np.testing.assert_allclose([r['w'][0] for r in b.rows],
            [1.-1e-6,1.,3.],atol=1e-10,rtol=0)
        np.testing.assert_allclose([b.heat,b.wall_work,b.material_work],
            [2.,2.5,.5],atol=1e-10,rtol=0)


if __name__=='__main__': unittest.main(verbosity=2)

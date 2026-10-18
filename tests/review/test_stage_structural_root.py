"""Independent manufactured roots; no claimed AFSD accuracy from these tests."""
from pathlib import Path
import sys
import unittest
from unittest import mock
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.stage_structural_root import solve_stage
from afsd_backend.pysph.stage_iteration import StageFailure


class Map:
    def __init__(self): self.calls=0; self.times=[]
    def candidate(self):
        return dict(position_m=np.full((2,3),.002),velocity_m_s=np.zeros((2,3)),
            density_kg_m3=np.full(2,2700.),temperature_K=np.full(2,723.),
            mean_velocity_m_s=np.zeros((2,3)))
    def evaluate(self,c,dt):
        self.calls+=1;self.times.append(dt)
        result={k:v.copy() for k,v in c.items()}
        result['velocity_m_s']=-2*c['velocity_m_s']+.002
        return dict(map=result)


class Proposal:
    def __init__(self,scale=1.):self.scale=scale
    def apply(self,c,g):
        return dict(proposal={k:c[k]+self.scale*(g[k]-c[k]) for k in c},
                    diagnostics={'scope':'MANUFACTURED_NUMERICAL_PROPOSAL_ONLY'})


class StructuralRootTests(unittest.TestCase):
    def test_literal_scalar_root_independent_of_nonsingular_scaling(self):
        for factor in (1.,.125,2.):
            obj=Map();r=solve_stage(obj,Proposal(factor),.125)
            np.testing.assert_allclose(r['candidate']['velocity_m_s'],1/1500,atol=1e-9,rtol=0)
            self.assertEqual(r['status'],'ROOT_CONVERGED_NOT_ACCEPTED')
            self.assertEqual(obj.times,[.125]*obj.calls)

    def test_zero_and_tiny_proposal_cannot_accept_nonzero_raw_residual(self):
        for factor in (0.,1e-20):
            with self.assertRaises(StageFailure) as ctx:
                solve_stage(Map(),Proposal(factor),.125,maximum_evaluations=3)
            self.assertGreater(ctx.exception.trace[-1]['maximum_scaled_residual'],1e6)
            self.assertLess(ctx.exception.trace[-1]['maximum_scaled_proposal'],1.)

    def test_five_fields_with_literal_independent_event_mean(self):
        class Coupled(Map):
            def evaluate(self,c,dt):
                g={k:v.copy() for k,v in c.items()}
                a=c['velocity_m_s']/1e-3;b=(c['position_m']-.002)/1e-5
                g['velocity_m_s']=1e-3*(-2*a+.5*b+3.5)
                g['position_m']=.002+1e-5*(.4*a-1.5*b+4.6)
                g['density_kg_m3'][:]=2700.27;g['temperature_K'][:]=723.02
                g['mean_velocity_m_s'][:]=-.004
                return dict(map=g)
        r=solve_stage(Coupled(),Proposal(.5),.001)['candidate']
        np.testing.assert_allclose(r['velocity_m_s'],1e-3*221/146,atol=1e-9,rtol=0)
        np.testing.assert_allclose(r['position_m'],.002+1e-5*152/73,atol=1e-12,rtol=0)
        np.testing.assert_allclose(r['mean_velocity_m_s'],-.004,atol=1e-9,rtol=0)
        np.testing.assert_allclose(r['temperature_K'],723.02,atol=5e-12,rtol=0)
        np.testing.assert_allclose(r['density_kg_m3'],2700.27,atol=2.7e-7,rtol=0)
        self.assertGreater(np.max(abs(r['mean_velocity_m_s']-.5*r['velocity_m_s'])),.004)

    def test_harmonic_midpoint_literal_root_actual_structural_proposal(self):
        from afsd_backend.pysph.structural_mobility import ConstrainedResponse
        from afsd_backend.pysph.stage_preconditioner import StagePreconditioner
        class Oscillator(Map):
            def candidate(self):
                c=super().candidate();c['velocity_m_s'][:]=.001;c['mean_velocity_m_s'][:]=.001
                return c
            def evaluate(self,c,dt):
                g={k:v.copy() for k,v in c.items()}
                # Literal force=-4*xmid; mass=1. This is NOT a tested-force call.
                v1=.001-4*dt*.5*(.002+c['position_m'])
                mean=.5*(.001+v1)
                g['velocity_m_s']=v1;g['mean_velocity_m_s']=mean;g['position_m']=.002+dt*mean
                return dict(map=g)
        response=ConstrainedResponse(np.ones(6),4*np.eye(6),1.,np.ones(6,dtype=bool))
        pre=StagePreconditioner.from_response(response,np.array([10,20]))
        r=solve_stage(Oscillator(),pre,1.)['candidate']
        np.testing.assert_allclose(r['position_m'],.0005,atol=1e-12,rtol=0)
        np.testing.assert_allclose(r['velocity_m_s'],-.004,atol=1e-9,rtol=0)
        np.testing.assert_allclose(r['mean_velocity_m_s'],-.0015,atol=1e-9,rtol=0)

    def test_raw_growth_rejection_preserves_best_anchor(self):
        class Growth(Map):
            def evaluate(self,c,dt):
                g={k:v.copy() for k,v in c.items()};u=c['velocity_m_s']
                g['velocity_m_s']=u+np.where(u<.0009,.002-2*u,.02+u)
                return dict(map=g)
        with self.assertRaises(StageFailure) as ctx:
            solve_stage(Growth(),Proposal(),.125,maximum_evaluations=7)
        self.assertTrue(any(r['candidate_disposition']=='REJECTED_ACCELERATED_CANDIDATE' for r in ctx.exception.trace))
        self.assertTrue(any(r['fallback']=='RAW_GROWTH_GT_2_BEST_REJECTED_RESET_HISTORY' for r in ctx.exception.trace))

    def test_failed_trial_keeps_actual_context_and_does_not_commit(self):
        obj=Map()
        def bad(c,dt):
            e=ValueError('contact failure');e.contact_context={'particle_id':20};raise e
        obj.evaluate=bad;rows=[]
        with self.assertRaises(StageFailure) as ctx:
            solve_stage(obj,Proposal(),.125,on_evaluation=lambda row,c,e:rows.append((row,e)))
        self.assertEqual(ctx.exception.original_context,{'particle_id':20})
        self.assertIsNone(ctx.exception.evaluation);self.assertEqual(len(rows),1)

    def test_nan_raw_map_never_replaced_by_preconditioner(self):
        obj=Map();actual=obj.evaluate
        def bad(c,dt):
            r=actual(c,dt);r['map']['temperature_K'][0]=np.nan;return r
        obj.evaluate=bad
        with self.assertRaises(StageFailure) as ctx:solve_stage(obj,Proposal(),.125)
        self.assertTrue(np.isnan(ctx.exception.evaluation['map']['temperature_K'][0]))

    def test_exact_evaluation_budget_and_input_purity(self):
        obj=Map();shared=obj.candidate();original={k:v.copy() for k,v in shared.items()}
        obj.candidate=lambda:shared
        with self.assertRaises(StageFailure) as ctx:solve_stage(obj,Proposal(),.125,maximum_evaluations=2)
        self.assertEqual(obj.calls,2);self.assertEqual(len(ctx.exception.trace),2)
        for k in shared:np.testing.assert_array_equal(shared[k],original[k])

    def test_wall_budget_before_any_evaluation(self):
        obj=Map()
        with mock.patch('afsd_backend.pysph.stage_structural_root.time.perf_counter',side_effect=[0.,1.]):
            with self.assertRaises(StageFailure):solve_stage(obj,Proposal(),.125,wall_seconds=.1)
        self.assertEqual(obj.calls,0)

    def test_invalid_settings_fail_early(self):
        for kwargs in ({'maximum_evaluations':0},{'maximum_evaluations':True},{'wall_seconds':0},
                       {'wall_seconds':float('inf')}):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):solve_stage(Map(),Proposal(),.125,**kwargs)


if __name__=='__main__':unittest.main(verbosity=2)

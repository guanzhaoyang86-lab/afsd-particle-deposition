"""Independent algebra references, not an AFSD or contact validation claim."""
from pathlib import Path
import json
import sys
import unittest
from unittest import mock
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.stage_anderson import solve_stage, _accelerate, _Coordinates
from afsd_backend.pysph.stage_iteration import StageFailure, solve_stage as picard


class Map:
    def __init__(self, slope=-2.): self.slope=slope; self.calls=0; self.times=[]
    def candidate(self):
        return dict(position_m=np.full((2,3), .002), velocity_m_s=np.zeros((2,3)),
                    density_kg_m3=np.full(2,2700.), temperature_K=np.full(2,723.),
                    mean_velocity_m_s=np.zeros((2,3)))
    def evaluate(self,c,dt):
        self.calls+=1; self.times.append(dt)
        result={k:v.copy() for k,v in c.items()}
        result['velocity_m_s']=self.slope*c['velocity_m_s']+.002
        return dict(map=result)


class AndersonTests(unittest.TestCase):
    def test_negative_slope_literal_root_where_picard_diverges(self):
        obj=Map(); result=solve_stage(obj,.125)
        np.testing.assert_allclose(result['candidate']['velocity_m_s'],1/1500,atol=1e-9,rtol=0)
        self.assertEqual(result['status'],'ROOT_CONVERGED_NOT_ACCEPTED')
        self.assertLessEqual(obj.calls,4)
        self.assertEqual(obj.times,[.125]*obj.calls)
        with self.assertRaises(StageFailure):picard(Map(),.125)

    def test_coupled_multidimensional_and_multiunit_literal_root(self):
        class Coupled(Map):
            def evaluate(self,c,dt):
                result={k:v.copy() for k,v in c.items()}
                # Hand equations: a'=-2a+0.5b+3.5, b'=0.4a-1.5b+4.6.
                # Independent elimination below yields the exact root.
                a=c['velocity_m_s']/1e-3
                b=(c['position_m']-.002)/1e-5
                result['velocity_m_s']=1e-3*(-2*a+.5*b+3.5)
                result['position_m']=.002+1e-5*(.4*a-1.5*b+4.6)
                result['temperature_K'][:]=723.+.02
                result['density_kg_m3'][:]=2700.+.27
                result['mean_velocity_m_s'][:]=-.004
                return dict(map=result)
        # 3 a - .5 b = 3.5; -.4 a + 2.5 b = 4.6.
        # Elimination gives a=221/146, b=152/73.
        r=solve_stage(Coupled(),.001)
        np.testing.assert_allclose(r['candidate']['velocity_m_s'],1e-3*221/146,atol=1e-9,rtol=0)
        np.testing.assert_allclose(r['candidate']['position_m'],.002+1e-5*152/73,atol=1e-12,rtol=0)
        np.testing.assert_allclose(r['candidate']['temperature_K'],723.02,atol=5e-12,rtol=0)
        np.testing.assert_allclose(r['candidate']['density_kg_m3'],2700.27,atol=2.7e-7,rtol=0)
        self.assertLessEqual(len(r['trace']),32)

    def test_piecewise_constraint_root_independent(self):
        class Piecewise(Map):
            def evaluate(self,c,dt):
                result={k:v.copy() for k,v in c.items()}
                u=c['velocity_m_s']
                # The two continuous branches meet at u=0.001.
                result['velocity_m_s']=np.where(u<.001, .5*u+.002, -.25*u+.00275)
                return dict(map=result)
        r=solve_stage(Piecewise(),.125)
        np.testing.assert_allclose(r['candidate']['velocity_m_s'],.0022,atol=1e-9,rtol=0)

    def test_rank_deficient_secants_explicitly_discarded_not_pinv(self):
        dx=np.array([[1.,2.],[2.,4.],[3.,6.]])
        df=-3*dx
        nxt,info=_accelerate(np.zeros(3),np.array([3.,6.,9.]),dx,df)
        np.testing.assert_allclose(nxt,[1.,2.,3.],atol=2e-15,rtol=0)
        self.assertEqual(info['secant_depth'],1)
        self.assertEqual(len(info['discarded_history']),1)
        self.assertGreater(info['discarded_history'][0]['conditioning'],1e10)

    def test_rate_kink_manufactured_map_literal_root(self):
        class Kink(Map):
            def evaluate(self,c,dt):
                result={k:v.copy() for k,v in c.items()}
                z=c['velocity_m_s']/1e-3
                # max is ONLY this manufactured test map's specified kink.
                # z=1.1-9(z-1) in z>1 gives exactly z=1.01.
                result['velocity_m_s']=1e-3*(1.1-9*np.maximum(z-1,0))
                return dict(map=result)
        r=solve_stage(Kink(),.125)
        np.testing.assert_allclose(r['candidate']['velocity_m_s'],.00101,atol=1e-9,rtol=0)

    def test_zero_secants_fallback_is_logged(self):
        nxt,info=_accelerate(np.zeros(2),np.ones(2),np.zeros((2,2)),np.zeros((2,2)))
        np.testing.assert_array_equal(nxt,[.75,.75])
        self.assertEqual(info['fallback'],'NO_WELL_CONDITIONED_SECANT')
        self.assertEqual(len(info['discarded_history']),2)

    def test_singular_history_diagnostic_is_explicit_strict_json(self):
        _,info=_accelerate(np.zeros(2),np.ones(2),np.zeros((2,2)),np.zeros((2,2)))
        trace=[dict(evaluation=1,**info)]
        reread=json.loads(json.dumps(trace,allow_nan=False))
        for row in reread[0]['discarded_history']:
            self.assertIsNone(row['conditioning'])
            self.assertEqual(row['conditioning_status'],'NONFINITE_SINGULAR')
        self.assertEqual(reread[0]['fallback'],'NO_WELL_CONDITIONED_SECANT')

    def test_large_absolute_temperature_has_zero_reference_coordinate(self):
        reference=Map().candidate(); coords=_Coordinates(reference)
        np.testing.assert_array_equal(coords.pack(reference),0.)
        candidate={k:v.copy() for k,v in reference.items()}; candidate['temperature_K']+=.02
        np.testing.assert_array_equal(coords.unpack(coords.pack(candidate))['temperature_K'],candidate['temperature_K'])
        reference['temperature_K'][:]=300.
        np.testing.assert_array_equal(coords.reference['temperature_K'],723.)

    def test_growth_rejection_callback_and_best_picard_restart(self):
        class Growth(Map):
            def evaluate(self,c,dt):
                self.calls+=1
                result={k:v.copy() for k,v in c.items()}
                # Deliberately noncontractive, deterministic piecewise scalar map.
                u=c['velocity_m_s']; f=np.where(u<.0009,.002-2*u,.02+u)
                result['velocity_m_s']=u+f
                return dict(map=result)
        events=[]
        with self.assertRaises(StageFailure) as ctx:
            solve_stage(Growth(),.125,maximum_evaluations=7,
                        on_evaluation=lambda row,c,e:events.append((dict(row),c['velocity_m_s'].copy())))
        rejected=[row for row,c in events if row['candidate_disposition']=='REJECTED_ACCELERATED_CANDIDATE']
        self.assertTrue(rejected)
        self.assertTrue(any(row['fallback']=='GROWTH_GT_2_BEST_REJECTED_RESET_HISTORY' for row,c in events))
        self.assertEqual(len(events),len(ctx.exception.trace))

    def test_failed_trial_preserves_actual_candidate_original_context(self):
        obj=Map()
        def reject(c,dt):
            exc=ValueError('known candidate failure');exc.contact_context={'particle_id':268};raise exc
        obj.evaluate=reject
        callbacks=[]
        with self.assertRaises(StageFailure) as ctx:
            solve_stage(obj,.125,on_evaluation=lambda row,c,e:callbacks.append((row,e)))
        self.assertEqual(ctx.exception.original_context,{'particle_id':268})
        self.assertIsNone(ctx.exception.evaluation)
        self.assertEqual(ctx.exception.trace[-1]['candidate_disposition'],'REJECTED_TRIAL_ERROR')
        self.assertEqual(callbacks[0][0]['evaluation_status'],'FAILED_TRIAL')
        self.assertIsNone(callbacks[0][1])
        np.testing.assert_array_equal(ctx.exception.candidate['temperature_K'],723.)

    def test_nan_result_rejected_not_clipped(self):
        obj=Map(); original=obj.evaluate
        def nan(c,dt):
            e=original(c,dt);e['map']['temperature_K'][0]=np.nan;return e
        obj.evaluate=nan
        with self.assertRaises(StageFailure) as ctx:solve_stage(obj,.125)
        self.assertTrue(np.isnan(ctx.exception.evaluation['map']['temperature_K'][0]))

    def test_evaluation_budget_and_trace_are_exact(self):
        obj=Map()
        with self.assertRaises(StageFailure) as ctx:solve_stage(obj,.125,maximum_evaluations=2)
        self.assertEqual(obj.calls,2);self.assertEqual(len(ctx.exception.trace),2)
        self.assertIn('evaluation limit',str(ctx.exception))

    def test_wall_budget_stops_before_any_physics_evaluation(self):
        obj=Map()
        with mock.patch('afsd_backend.pysph.stage_anderson.time.perf_counter',side_effect=[0.,1.]):
            with self.assertRaises(StageFailure):solve_stage(obj,.125,wall_seconds=.1)
        self.assertEqual(obj.calls,0)

    def test_readonly_initial_state_and_invalid_settings(self):
        shared=Map().candidate(); before={k:v.copy() for k,v in shared.items()}
        obj=Map();obj.candidate=lambda:shared
        solve_stage(obj,.125)
        for k in shared:np.testing.assert_array_equal(shared[k],before[k])
        for kwargs in ({'maximum_evaluations':0},{'maximum_evaluations':True},{'wall_seconds':0},
                       {'wall_seconds':float('inf')}):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):solve_stage(Map(),.125,**kwargs)


if __name__=='__main__':unittest.main(verbosity=2)

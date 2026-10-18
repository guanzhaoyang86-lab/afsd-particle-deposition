"""Independent fixed-point algebra; not an AFSD solution or contact reference."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.stage_iteration import initial_candidate, solve_stage, residual_norms, StageFailure


class AffineMap:
    def __init__(self, slope=0.): self.slope=slope; self.calls=0
    def candidate(self):
        return dict(position_m=np.zeros((2,3)), velocity_m_s=np.zeros((2,3)),
                    density_kg_m3=np.full(2, 2700.), temperature_K=np.full(2, 400.),
                    mean_velocity_m_s=np.zeros((2,3)))
    def evaluate(self,c,dt):
        self.calls+=1; result={k:v.copy() for k,v in c.items()}
        result['velocity_m_s']=self.slope*c['velocity_m_s']+.002
        return dict(map=result)


class StageIterationTests(unittest.TestCase):
    def test_constant_map_literal_root(self):
        obj=AffineMap(); r=solve_stage(obj,.1)
        np.testing.assert_allclose(r['candidate']['velocity_m_s'],.002,atol=1e-9,rtol=0)
        self.assertEqual(r['status'],'ROOT_CONVERGED_NOT_ACCEPTED')
        self.assertLessEqual(obj.calls,32)
    def test_negative_slope_literal_root_and_not_physical_damping(self):
        r=solve_stage(AffineMap(-.5),.1)
        np.testing.assert_allclose(r['candidate']['velocity_m_s'],1/750,atol=1e-9,rtol=0)
    def test_divergence_preserves_every_trace_not_accept(self):
        obj=AffineMap(2.)
        with self.assertRaises(StageFailure) as ctx: solve_stage(obj,.1,maximum_evaluations=3)
        self.assertEqual(obj.calls,3); self.assertEqual(len(ctx.exception.trace),3)
        self.assertGreater(ctx.exception.trace[-1]['maximum_scaled_residual'],1)
    def test_callback_keeps_physical_dt_no_auto_shrink(self):
        times=[]; obj=AffineMap(); original=obj.evaluate
        def wrapped(c,dt): times.append(dt); return original(c,dt)
        obj.evaluate=wrapped; record=[]
        solve_stage(obj,.125,on_evaluation=lambda row,c,e:record.append(row['evaluation']))
        self.assertEqual(times,[.125]*len(times)); self.assertEqual(record,list(range(1,len(times)+1)))
    def test_domain_failure_context_retained(self):
        obj=AffineMap()
        def reject(c,dt):
            exc=ValueError('known rejected contact'); exc.contact_context={'particle_id':123}; raise exc
        obj.evaluate=reject
        with self.assertRaises(StageFailure) as ctx: solve_stage(obj,.1)
        self.assertEqual(ctx.exception.original_context,{'particle_id':123})
    def test_residual_density_is_relative_temperature_absolute(self):
        c=AffineMap().candidate(); mapping={k:v.copy() for k,v in c.items()}
        mapping['density_kg_m3']+=.0027; mapping['temperature_K']+=.125
        norms,_=residual_norms(c,mapping,c['density_kg_m3'])
        self.assertAlmostEqual(norms['density_kg_m3'],1e-6,places=14)
        self.assertEqual(norms['temperature_K'],.125)
    def test_zero_budget_and_bad_relaxation_rejected(self):
        for kwargs in ({'maximum_evaluations':0},{'maximum_evaluations':True},{'relaxation':0},
                       {'relaxation':1.1},{'wall_seconds':0}):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):solve_stage(AffineMap(),.1,**kwargs)
    def test_initial_guess_does_not_mutate_input(self):
        class Shared(AffineMap):
            def candidate(self):return shared
        shared=AffineMap().candidate(); shared['mean_velocity_m_s'][:]=2.
        result=initial_candidate(Shared(),.1)
        np.testing.assert_array_equal(shared['position_m'],0.)
        np.testing.assert_allclose(result['position_m'],.2,atol=1e-15)
    def test_nan_map_is_rejected_not_clipped(self):
        obj=AffineMap(); original=obj.evaluate
        def bad(c,dt):
            r=original(c,dt);r['map']['temperature_K'][0]=np.nan;return r
        obj.evaluate=bad
        with self.assertRaises(StageFailure):solve_stage(obj,.1)


if __name__=='__main__':unittest.main(verbosity=2)

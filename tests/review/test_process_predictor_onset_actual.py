"""Actual integration guards; no full physical path inside this test module."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from probe_process_predictor_onset import PredictorOnsetBlock
from probe_process_plastic_friction_onset import OnsetBlock,boundary_eligibility
from probe_process_plastic_friction import PlasticFrictionBlock
from probe_process_plastic import CandidateFailure
from diagnose_process_plastic_friction_q7 import state_signature


class PredictorOnsetIntegration(unittest.TestCase):
    def test_01_only_branch_and_named_metadata_change(self):
        old,new=OnsetBlock(),PredictorOnsetBlock()
        a,b=state_signature(old),state_signature(new)
        b['metadata']['pfcfg'].pop('predictor_variant')
        self.assertEqual(a,b)
        self.assertIs(PredictorOnsetBlock.advance,OnsetBlock.advance)
        for name in ('commit','trial','tangent','state','wall'):
            self.assertIs(getattr(PredictorOnsetBlock,name),getattr(PlasticFrictionBlock,name))
        self.assertEqual(new.qcfg['numerical']['maximum_iterations'],14)
        self.assertEqual(new.pfcfg['predictor_variant']['run_timeout_s'],180)

    def test_02_preload_events_history_and_tags(self):
        b=PredictorOnsetBlock(); records=[]; trace=[]
        b.trace_sink=trace.append
        b.advance(.025,lambda:records.append(b.record()))
        self.assertEqual(b.time,.025)
        self.assertGreater(b.count,1)
        self.assertGreater(len(b.events),0)
        self.assertGreater(float(b.pa.ep.max()),0.)
        self.assertGreater(float(b.pa.plastic_work.sum()),0.)
        self.assertGreater(float(b.friction.sum()),0.)
        for e in b.events:
            self.assertLessEqual(e['bracket_s'][1]-e['bracket_s'][0],1e-12)
        self.assertEqual(len(records),b.count)
        self.assertIsNone(trace[-1]['error'])
        mask,_=boundary_eligibility(b)
        before=state_signature(b)
        first=b.branch(.0251,mask)
        self.assertEqual(before,state_signature(b))
        b.pa.tracer_type[:]=np.arange(len(b.X))%5+20
        b.pa.tag[:]=np.arange(len(b.X))%3
        ids=b.pa.particle_id.copy(); before=state_signature(b)
        second=b.branch(.0251,mask)
        for name in ('displacement','stress','S','ep','plastic_increment_J',
                     'normal_force','tangent_force','relative_increment'):
            np.testing.assert_array_equal(first[name],second[name])
        self.assertEqual(before,state_signature(b))
        np.testing.assert_array_equal(ids,b.pa.particle_id)

    def test_03_original_commit_rejects_bad_candidate_atomically(self):
        b=PredictorOnsetBlock()
        c=b.branch(.0125,np.zeros(len(b.nodes),dtype=bool))
        before=state_signature(b)
        c['displacement']=c['displacement'].copy()
        c['displacement'][b.nodes[0],0]+=1e-6
        with self.assertRaises((CandidateFailure,ValueError)):
            b.commit(c)
        self.assertEqual(before,state_signature(b))
        stale=dict(c,base_step=-1)
        with self.assertRaisesRegex(ValueError,'stale'):
            b.commit(stale)
        self.assertEqual(before,state_signature(b))

    def test_04_valid_first_commit_recomputes_not_trusted_force(self):
        b=PredictorOnsetBlock()
        c=b.branch(.0125,np.zeros(len(b.nodes),dtype=bool))
        expected=c['normal_force'].copy()
        c['normal_force'][:]=123456.
        c['contact_force'][:]=np.nan
        c['event_kind']='REQUESTED'
        b.commit(c)
        self.assertEqual(b.count,1); self.assertEqual(b.time,.0125)
        np.testing.assert_array_equal(b.normal_force,expected)
        self.assertTrue(np.isfinite(b.contact_force).all())
        self.assertEqual(b.events,[])

    def test_05_budget_corner_guards_and_positive_seed_compatibility(self):
        old,b=OnsetBlock(),PredictorOnsetBlock()
        mask=np.zeros(len(b.nodes),dtype=bool)
        a=old.branch(.0125,mask); c=b.branch(.0125,mask)
        for name in ('displacement','stress','S','ep','normal_force','tangent_force'):
            np.testing.assert_array_equal(a[name],c[name])
        self.assertEqual(c['predictor_corrections'],0)
        before=state_signature(b)
        for v in (-1,15,True,1.5):
            with self.subTest(v=v),self.assertRaises(ValueError):
                b.branch(.0125,mask,maximum_iterations=v)
        with self.assertRaisesRegex(ValueError,'crosses'):
            b.advance(.1)
        self.assertEqual(before,state_signature(b))


if __name__=='__main__':
    unittest.main(verbosity=2)

"""Actual small-block guards for the onset entry; no full .2s path claim."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from diagnose_process_plastic_friction_q7 import state_signature
from probe_process_plastic_friction import PlasticFrictionBlock
from probe_process_plastic_friction_recovery import RecoveryBlock
from probe_process_plastic_friction_onset import OnsetBlock,boundary_eligibility
from probe_process_plastic import CandidateFailure


class ActualOnset(unittest.TestCase):
    def test_01_only_entry_metadata_changes_and_zero_capacity_not_released(self):
        old,new=RecoveryBlock(),OnsetBlock()
        a,b=state_signature(old),state_signature(new)
        b['metadata']['pfcfg'].pop('onset_variant')
        b['metadata'].pop('trace_sink')
        self.assertEqual(a,b)
        mask,margin=boundary_eligibility(new)
        self.assertFalse(mask.any())
        np.testing.assert_array_equal(margin,np.zeros(9))
        for name in ('commit','trial','tangent','state','wall'):
            self.assertIs(getattr(OnsetBlock,name),getattr(PlasticFrictionBlock,name))
        self.assertEqual(new.qcfg['numerical']['maximum_iterations'],14)

    def test_02_boundary_classification_does_not_normalize_raw_forces(self):
        b=OnsetBlock()
        # Classification-only diagnostic fixture; not an equilibrium solution.
        b.normal_force[:]=10.
        b.tangent_force[:]=np.array([1.,0.])
        b.tangent_force[:3,0]=[1.+4e-5,1.-4e-5,1.-2e-4]
        before=state_signature(b)
        mask,margin=boundary_eligibility(b)
        np.testing.assert_array_equal(mask[:3],[True,True,False])
        np.testing.assert_allclose(margin[:3],[-4e-5,4e-5,2e-4],atol=1e-15,rtol=0)
        self.assertEqual(state_signature(b),before)
        b.tangent_force[0,0]=1.001
        before=state_signature(b)
        with self.assertRaisesRegex(CandidateFailure,'outside cone tolerance'):
            boundary_eligibility(b)
        self.assertEqual(state_signature(b),before)

    def test_03_actual_zero_load_preload_events_and_passive_tags(self):
        b=OnsetBlock(); ids=b.pa.particle_id.copy()
        saved=[]; diagnostics=[]
        b.trace_sink=diagnostics.append
        b.advance(.025,lambda:saved.append(b.record()))
        self.assertEqual(b.time,.025)
        self.assertGreater(float(b.pa.ep.max()),0.)
        self.assertGreater(float(b.pa.plastic_work.sum()),0.)
        self.assertGreater(float(b.friction.sum()),0.)
        self.assertTrue(b.events)
        self.assertEqual(diagnostics[0]['error'],None)
        self.assertEqual(diagnostics[0]['last_accepted_step'],len(saved))
        firstclass=next(x for x in diagnostics[0]['trace'] if x['stage']=='CLASSIFICATION')
        self.assertFalse(any(firstclass['eligible']))
        for event in b.events:
            self.assertLessEqual(event['bracket_s'][1]-event['bracket_s'][0],1e-12)
        mask,_=boundary_eligibility(b)
        first=b.branch(.0251,mask)
        b.pa.tracer_type[:]=np.arange(27)%5+44
        before=state_signature(b)
        second=b.branch(.0251,mask)
        for key in ('displacement','S','ep','plastic_increment_J','normal_force','tangent_force'):
            np.testing.assert_array_equal(first[key],second[key])
        self.assertEqual(state_signature(b),before)
        np.testing.assert_array_equal(b.pa.particle_id,ids)


if __name__=='__main__':
    unittest.main(verbosity=2)

"""New-model history isolation and candidate publication invariants."""
from pathlib import Path
import sys
import json
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.discrete_process import (
    independent_initial_state,DiscreteContactTrial,accepted_copy)
from afsd_backend.pysph import quasistatic_process as old


def parameters():
    p=json.loads((ROOT/'outputs/contact_material_search_20261004/root/resolved_config.json').read_text())['parameters']
    p['bore_width']=.002504
    return p


class DiscreteProcessTests(unittest.TestCase):
    def test_old_model_and_reset_interval_cannot_be_resumed(self):
        p=parameters();pa,state=old.independent_initial_state(p)
        with self.assertRaisesRegex(ValueError,'own independent history'):
            DiscreteContactTrial(pa,state,p)
        pa,state=independent_initial_state(p);state['step']=1
        with self.assertRaisesRegex(ValueError,'cannot be reset'):
            DiscreteContactTrial(pa,state,p)

    def test_rest_commit_keeps_identity_and_publishes_single_interval_record(self):
        p=parameters();p.update(feed=0.,omega=0.,traverse=0.,tool_temperature=723.15,base_temperature=723.15)
        pa,state=independent_initial_state(p);trial=DiscreteContactTrial(pa,state,p)
        r=trial.evaluate(np.zeros((len(pa.x),3)),.0005)
        other,new,fresh=accepted_copy(trial,r)
        self.assertEqual(state['step'],0);self.assertEqual(state['interval_energy']['steps'],0)
        self.assertEqual(new['step'],1);self.assertEqual(new['interval_energy']['steps'],1)
        self.assertEqual(new['crossing_ids'],[])
        np.testing.assert_array_equal(other.particle_id,pa.particle_id)
        np.testing.assert_array_equal(other.m,pa.m)
        self.assertLess(fresh['force_residual_N'],1e-10)

    def test_mutated_unconverged_candidate_is_recomputed_and_rejected(self):
        p=parameters();pa,state=independent_initial_state(p);trial=DiscreteContactTrial(pa,state,p)
        r=trial.evaluate(trial.prescribed(.0005),.0005)
        r['scaled_residual']=0.;r['force'][:]=0.;r['stress'][:]=0.
        before=pa.x.copy()
        with self.assertRaises(ValueError): accepted_copy(trial,r)
        np.testing.assert_array_equal(pa.x,before)
        self.assertEqual(state['step'],0)


if __name__=='__main__':unittest.main()

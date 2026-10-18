"""Q5 commit and inventory guards; no complete process rerun."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from probe_process_friction import FrictionBlock, CandidateFailure, solve_friction_modes


def accepted_snapshot(b):
    result={}
    # Include all directly retained arrays/scalars, not just positions/stress.
    for key,value in vars(b).items():
        if isinstance(value,np.ndarray): result['array:'+key]=value.copy()
        elif isinstance(value,(float,int,str,bool,type(None))): result['scalar:'+key]=value
    for key,prop in b.pa.properties.items():
        result['particle:'+key]=prop.get_npy_array().copy()
    result['events']=deepcopy(b.events)
    result['last_trace']=deepcopy(b.last_trace)
    return result


class FrictionGuardRegression(unittest.TestCase):
    def test_01_actual_unsplit_slip_candidate_is_rejected_atomically(self):
        b=FrictionBlock()
        b.advance(.05)
        before=accepted_snapshot(b)
        c=solve_friction_modes(lambda modes:b.branch(.1,modes),len(b.nodes),b.mu,b.limits)
        self.assertTrue(np.any(c['modes']!=0),'candidate must exercise actual sliding')
        self.assertTrue(np.any(np.abs(c['relative_increment'])>b.limits['slip_tolerance_m']))
        with self.assertRaisesRegex(CandidateFailure,'unsplit stick/slip transition'):
            b.commit(c)
        after=accepted_snapshot(b)
        self.assertEqual(set(before),set(after))
        for key in before:
            if isinstance(before[key],np.ndarray):
                np.testing.assert_array_equal(after[key],before[key],err_msg=key)
            else: self.assertEqual(after[key],before[key],key)

    def test_02_nonfinite_inventory_never_counts_as_feasible_contact(self):
        b=FrictionBlock()
        for field in ('normal_force','tangent_force','relative_increment'):
            for invalid in (np.nan,np.inf,-np.inf):
                with self.subTest(field=field,invalid=invalid):
                    data=dict(normal_force=np.array([2.]),tangent_force=np.array([0.]),relative_increment=np.array([0.]))
                    data[field][0]=invalid
                    with self.assertRaisesRegex(CandidateFailure,'invalid/nonfinite contact inventory'):
                        solve_friction_modes(lambda modes:data,1,.5,b.limits)

    def test_03_mismatched_inventory_shape_is_rejected(self):
        b=FrictionBlock()
        for field in ('normal_force','tangent_force','relative_increment'):
            with self.subTest(field=field):
                data=dict(normal_force=np.array([2.]),tangent_force=np.array([0.]),relative_increment=np.array([0.]))
                data[field]=np.zeros(2)
                with self.assertRaisesRegex(CandidateFailure,'invalid/nonfinite contact inventory'):
                    solve_friction_modes(lambda modes:data,1,.5,b.limits)


if __name__=='__main__': unittest.main(verbosity=2)

"""Fake metadata/coverage regression; no ProcessThermalTrial.evaluate call."""
import copy
from pathlib import Path
import sys,unittest
from types import SimpleNamespace
import numpy as np

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from probe_afsd_process_bridge import signature,mechanism_coverage


class Driver(unittest.TestCase):
    def example(self):
        return dict(boundary=dict(traces=[dict(particle_id=7,free_axes=[1,1,0],
            contact_rows=[dict(contact_id='7:A',feature_id=1),dict(contact_id='7:B',feature_id=2)],
            result=dict(intervals=[dict(kind='SLIDE')],events=[dict(event='ZERO_SLIP',contacts=[0],time_s=.1)]))]),
            material=dict(ep=np.array([.02]),strain_rate_per_s=np.array([2.])),
            normal_before=dict(active_contact_ids=['7:A','7:B']),normal_after=dict(active_contact_ids=['7:A','7:B']))

    def test_event_face_change_cannot_borrow_same_branch(self):
        a=self.example();b=copy.deepcopy(a);b['boundary']['traces'][0]['result']['events'][0]['contacts']=[1]
        pa=SimpleNamespace(ep=np.array([.01]));p=dict(rate_floor=1.)
        sa=signature(a,pa,p);sb=signature(b,pa,p)
        self.assertEqual(sa['traces'][0]['events'][0]['contacts'],[('7:A',1)])
        self.assertEqual(sb['traces'][0]['events'][0]['contacts'],[('7:B',2)])
        self.assertNotEqual(sa['traces'][0]['events'],sb['traces'][0]['events'])

    def test_event_time_move_not_itself_branch_switch(self):
        a=self.example();b=copy.deepcopy(a);b['boundary']['traces'][0]['result']['events'][0]['time_s']=.1001
        pa=SimpleNamespace(ep=np.array([.01]));p=dict(rate_floor=1.)
        self.assertEqual(signature(a,pa,p)['traces'][0]['events'],signature(b,pa,p)['traces'][0]['events'])
        b['boundary']['traces'][0]['result']['events'][0]['contacts']=[2]
        with self.assertRaises(ValueError):signature(b,pa,p)

    def test_coverage_all_three_mechanisms_required(self):
        s=dict(positive_plastic_particle_count=1,positive_plastic_work_J=1e-9,contact_count=2,
            multicontact_particle_count=1,actuator_active_particle_count=1,
            nonzero_constraint_impulse_particle_count=1,constraint_work_J=-1e-10)
        self.assertTrue(mechanism_coverage(s))
        for key in s:
            bad=dict(s);bad[key]=0
            with self.subTest(key=key):self.assertFalse(mechanism_coverage(bad))


if __name__=='__main__':unittest.main(verbosity=2)

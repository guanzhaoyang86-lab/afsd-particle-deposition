"""Manufactured algebra only: no NNPS, equilibrium, or thermal solver."""
import copy
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from test_inlet_source_work import fixture
import fixed_inventory_source_work as fixed


def case(**options):
    _,prepared,actual,_,args=fixture(**options)
    # All three particles already exist. This is NOT the fixture's birth event.
    parent=prepared.extract_particles(np.arange(3,dtype=np.uint32))
    args['parameters']['room']=300.;args['parent_contact_storage_J']=np.zeros((3,2))
    return parent,prepared,actual,args


class FixedInventoryWorkTests(unittest.TestCase):
    def test_real_reaction_work_and_heat_close_without_mass_source(self):
        parent,prepared,actual,args=case();r=fixed.evaluate(parent,prepared,actual,**args)
        self.assertEqual(r['source_increment']['mass_kg'],0.)
        self.assertEqual(r['source_increment']['sensible_heat_J'],0.)
        self.assertIsNone(r['birth_event'])
        self.assertGreater(r['computed_terms']['mechanical_reaction_work']['value'],0.)
        self.assertAlmostEqual(r['open_energy']['open_raw_residual_J'],0.,delta=1e-12)

    def test_negative_reaction_and_unrepaired_defect_are_retained(self):
        a,b,c,args=case(unload=True);r=fixed.evaluate(a,b,c,**args)
        self.assertLess(r['computed_terms']['mechanical_reaction_work']['value'],0.)
        a,b,c,args=case(raw_defect=.012);r=fixed.evaluate(a,b,c,**args)
        self.assertAlmostEqual(r['open_energy']['open_raw_residual_J'],.012,delta=1e-12)

    def test_actual_sorting_preserves_calculation(self):
        a,b,c,args=case();r=fixed.evaluate(a,b,c,**args)
        c=c.extract_particles(np.array([2,0,1],dtype=np.uint32));s=fixed.evaluate(a,b,c,**args)
        self.assertEqual(r,s)

    def test_birth_or_changed_old_state_is_rejected(self):
        a,b,c,args=case();a=a.extract_particles(np.array([0,1],dtype=np.uint32))
        with self.assertRaises(ValueError):fixed.evaluate(a,b,c,**args)
        a,b,c,args=case();b.temperature[0]+=1.
        with self.assertRaises(ValueError):fixed.evaluate(a,b,c,**args)

    def test_wrong_reaction_or_nonzero_exterior_traction_is_rejected(self):
        for mode in ('reaction','traction'):
            a,b,c,args=case()
            if mode=='reaction':args['mechanical_energy']['prescribed_boundary_work_J']+=.01
            else:args['external_port']['traction_Pa'][2]=1.
            with self.subTest(mode=mode),self.assertRaises(ValueError):fixed.evaluate(a,b,c,**args)

    def test_artifact_is_bound_and_tampering_rejected(self):
        a,b,c,args=case();r=fixed.evaluate(a,b,c,**args)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'work.json';path.write_text(json.dumps(r),encoding='utf-8')
            w=fixed.record(r,path);self.assertNotEqual(w['mechanical_reaction_work']['value'],0.)
            changed=copy.deepcopy(r);changed['computed_terms']['mechanical_reaction_work']['value']=0.
            path.write_text(json.dumps(changed),encoding='utf-8')
            with self.assertRaises(ValueError):fixed.record(changed,path)


if __name__=='__main__':unittest.main(verbosity=2)

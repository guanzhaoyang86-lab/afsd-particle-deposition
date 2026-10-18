"""Independent arithmetic/metadata fixtures, NOT a continuum simulation."""
from pathlib import Path
import copy
import json
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import inspect_axial_feed as reader


def interval(base=False,released=False):
    vb=[1.5,-3.,-.5] if released else [1.,-2.,-.5]
    va=[2.,-4.,.25] if released else [1.5,-3.,-.5]
    nc=int(base);mid=.5*(np.asarray(vb)+np.asarray(va));tangent=np.tile(mid,(nc,1));tangent[:,2]=0.
    return dict(dt_s=.25,start_s=0.,v_before=vb,v_after=va,normal_impulse=[[0.,0.,.5]] if base else [],
        tangential_impulse=np.zeros((nc,3)).tolist(),constraint_impulse=[0.,0.,0.] if released else [0.,0.,-2. if base else -1.5],
        tangential_relative_midpoint=tangent.tolist(),friction_heat_raw_J=[0.]*nc,friction_rounding_bound_J=[0.]*nc,
        other_work_J=8.5625 if released else 5.5,constraint_work_J=0. if released else 1. if base else .75,
        wall_work_J=0.,delta_kinetic_J=8.5625 if released else 6.25,normal_work_J=[-.25] if base else [],
        kinetic_work_residual_J=0.)


def audit(a,base=False,released=False):
    return reader.audit_interval(a,np.array([1.5,-3.,-.5] if released else [1.,-2.,-.5]),np.array([4.,-8.,6.]),
        2. if base else 0.,[1.,1.,1.] if released else [1.,1.,0.],[0.,0.,-.5],2.,.25,0.)


class AxialIndependentReadback(unittest.TestCase):
    def test_01_literal_reaction_and_actuator_work(self):
        r=audit(interval());np.testing.assert_array_equal(r['v_after'],[1.5,-3.,-.5])
        np.testing.assert_array_equal(r['actuator_impulse'],[0.,0.,-1.5])
        self.assertEqual(r['actuator'],.75);self.assertEqual(r['internal'],5.5);self.assertEqual(r['kinetic'],6.25)

    def test_02_normal_contact_reaction_and_negative_normal_work(self):
        r=audit(interval(base=True),base=True)
        np.testing.assert_array_equal(r['normal_impulse'],[0.,0.,.5])
        np.testing.assert_array_equal(r['actuator_impulse'],[0.,0.,-2.])
        self.assertEqual(r['actuator'],1.);self.assertEqual(r['normal'],-.25)
        self.assertEqual(r['kinetic'],6.25)

    def test_03_release_preserves_velocity_and_free_acceleration(self):
        r=audit(interval(released=True),released=True)
        np.testing.assert_array_equal(r['v_after'],[2.,-4.,.25])
        np.testing.assert_array_equal(r['actuator_impulse'],[0.,0.,0.])
        self.assertEqual(r['internal'],8.5625);self.assertEqual(r['kinetic'],8.5625);self.assertEqual(r['actuator'],0.)

    def test_04_omitted_reaction_and_locked_transverse_dofs_rejected(self):
        for wrong in ([0.,0.,0.],[1.,0.,-1.5]):
            a=interval();a['constraint_impulse']=wrong
            with self.assertRaisesRegex(ValueError,'actual impulse balance'):audit(a)

    def test_05_negative_normal_work_not_clipped_to_heat(self):
        for key,value in [('normal_work_J',[0.]),('friction_heat_raw_J',[.25]),('constraint_work_J',.75)]:
            a=interval(base=True);a[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):audit(a,base=True)

    def test_06_wall_work_cannot_be_confused_with_material_normal_work(self):
        a=interval(base=True);a['wall_work_J']=-.25
        with self.assertRaisesRegex(ValueError,'saved wall_work'):audit(a,base=True)

    def test_07_release_zeroing_and_bad_timing_rejected(self):
        a=interval(released=True);a['v_before']=[0.,0.,0.]
        with self.assertRaisesRegex(ValueError,'velocity continuity'):audit(a,released=True)
        for key,value in [('dt_s',-.25),('start_s',.25),('v_after',[1.5,float('nan'),-.5])]:
            a=interval();a[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):audit(a)

    def test_08_world_geometry_nonzero_initial_kinetic_and_momentum(self):
        x=reader.literal_geometry();self.assertEqual(x.shape,(64,3))
        self.assertEqual(np.linalg.matrix_rank(x-x.mean(axis=0)),3)
        z=np.unique(x[:,2]);np.testing.assert_allclose(z,[-.00975,-.00925,-.00875,-.00825],atol=2e-18,rtol=0.)
        active=x[:,2]>-.0085;v=np.zeros((64,3));v[active,2]=-.00212
        self.assertEqual(np.count_nonzero(active),16)
        self.assertAlmostEqual(.5*3.375e-7*np.sum(v*v),1.213488e-11,delta=1e-26)
        np.testing.assert_allclose(3.375e-7*v.sum(axis=0),[0.,0.,-1.1448e-8],rtol=0.,atol=1e-23)

    def test_09_fixed_material_selection_and_right_hand_release_phase(self):
        z0=reader.literal_geometry()[:,2]
        held=reader.actuator_reference(z0,3.999999e-7);released=reader.actuator_reference(z0,4e-7)
        self.assertEqual(np.count_nonzero(held['active']),16);self.assertEqual(np.count_nonzero(released['active']),0)
        np.testing.assert_array_equal(held['eligible'],released['eligible'])
        np.testing.assert_array_equal(released['free_axes'],np.ones((64,3)))
        self.assertAlmostEqual(released['commanded_displacement_m'],-8.48e-10,delta=1e-25)
        self.assertEqual(released['phase'],'RELEASED')

    def test_10_control_tamper_cannot_lock_x_or_drive_current_free_positions(self):
        z0=reader.literal_geometry()[:,2];good=reader.actuator_reference(z0,0.)
        good={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in good.items()}
        good['displacement_status']='ACTUATOR_COMMAND_ONLY_NOT_FREE_MATERIAL_PATH'
        reader.audit_control(good,z0,0.)
        bad=copy.deepcopy(good);bad['free_axes'][0][0]=0.
        with self.assertRaisesRegex(ValueError,'actuator free_axes'):reader.audit_control(bad,z0,0.)
        bad=copy.deepcopy(good);bad['eligible'][0]=True
        with self.assertRaisesRegex(ValueError,'original-material eligible'):reader.audit_control(bad,z0,0.)

    def test_11_base_location_and_unilateral_signed_gap(self):
        r=reader.normal_reference(np.array([[0.,0.,-.0098],[0.,0.,-.0097]]))
        np.testing.assert_allclose(r['gap'],[.0002,.0003],rtol=0.,atol=2e-18)
        np.testing.assert_allclose(r['force'],[32.875,0.],rtol=0.,atol=1e-12)
        self.assertLess(r['signed'][1],0.);self.assertFalse(r['active'][1]);self.assertEqual(r['potential'][1],0.)
        self.assertAlmostEqual(r['potential'][0],.000821875,delta=3e-17)

    def test_12_reader_does_not_import_model_as_reference(self):
        source=(ROOT/'scripts/inspect_axial_feed.py').read_text(encoding='utf-8')
        for forbidden in ('from afsd_backend','import pysph','import probe_axial','from probe_axial'):
            self.assertNotIn(forbidden,source)


if __name__=='__main__':unittest.main(verbosity=2)

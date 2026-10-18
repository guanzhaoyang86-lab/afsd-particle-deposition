"""Actual old kick plus new passive-selection boundary interface regressions."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from afsd_backend.pysph.axial_actuator import control,validate_interval
from afsd_backend.pysph.contact_kick import kick
from afsd_backend.pysph.geometry import surfaces
from afsd_backend.pysph.rigid_contact_kinematics import contact_rows
from probe_axial_actuator_reference import benchmark


class AxialActuatorTests(unittest.TestCase):
    def c(self,t=0,z=None):
        if z is None:z=-.01+(np.arange(4)+.5)*.0005
        return control(z,t,-.01,.002,.0005,.00212,4e-7)

    def test_01_original_single_layer_and_only_z_constrained(self):
        r=self.c();np.testing.assert_array_equal(r['eligible'],[False,False,False,True])
        np.testing.assert_array_equal(r['free_axes'],[[1,1,1],[1,1,1],[1,1,1],[1,1,0]])
        np.testing.assert_array_equal(r['prescribed_velocity'][-1],[0,0,-.00212])

    def test_02_release_right_hand_state_has_no_prescribed_motion(self):
        r=self.c(4e-7);self.assertEqual(r['phase'],'RELEASED')
        np.testing.assert_array_equal(r['active'],np.zeros(4,bool))
        np.testing.assert_array_equal(r['free_axes'],np.ones((4,3)))
        self.assertAlmostEqual(r['commanded_displacement_m'],-8.48e-10,delta=1e-24)

    def test_03_stable_original_subset_under_permutation(self):
        z=-.01+(np.arange(4)+.5)*.0005;order=[3,1,0,2]
        a=self.c(z=z);b=self.c(z=z[order])
        np.testing.assert_array_equal(a['eligible'][order],b['eligible'])

    def test_04_labels_absent_from_boundary_api(self):
        import inspect
        names=set(inspect.signature(control).parameters)
        self.assertFalse(names&{'tracer_type','as_rs_tag','material_id','particle_id'})

    def test_05_release_straddling_rejected_not_silently_averaged(self):
        with self.assertRaisesRegex(ValueError,'step boundary'):validate_interval(3.99e-7,5e-9,4e-7)
        validate_interval(3.95e-7,5e-9,4e-7);validate_interval(4e-7,5e-9,4e-7)

    def test_06_parameter_domain_guards(self):
        for t in (-1.,np.nan):
            with self.assertRaises(ValueError):self.c(t)
        with self.assertRaises(ValueError):control([0],0,0,1,1,.1,1)
        with self.assertRaises(ValueError):control([np.nan],0,0,1,.1,.1,1)

    def test_07_literal_driven_and_released_reference(self):
        r=benchmark(0.);self.assertEqual(r['status'],'PASS_LITERAL_AXIAL_ACTUATOR_REFERENCE_NOT_AFSD')
        np.testing.assert_allclose(r['records'][-1]['velocity_m_s'],[2,-4,.25],rtol=0,atol=1e-12)

    def test_08_literal_base_normal_and_actuator_work(self):
        r=benchmark(2.);first=r['traces'][0]['result']
        self.assertAlmostEqual(first['constraint_work_J'],1.)
        self.assertAlmostEqual(sum(sum(s['normal_work_J']) for s in first['intervals']),-.25)

    def test_09_wrong_initial_velocity_is_not_overwritten(self):
        with self.assertRaisesRegex(ValueError,'explicit actuator event'):
            kick(np.zeros(3),2.,np.array([4.,-8.,6.]),[],[],[],0.,.25,
                free_axes=[1,1,0],prescribed_velocity=[0,0,-.5])

    def test_10_free_motion_not_clamped_to_zero_prescription(self):
        r=kick(np.array([1.,2.,-.5]),2.,np.zeros(3),[],[],[],0.,.25,
            free_axes=[1,1,1],prescribed_velocity=np.zeros(3))
        np.testing.assert_array_equal(r['velocity'],[1,2,-.5])
        np.testing.assert_array_equal(r['constraint_impulse'],[0,0,0])

    def test_11_geometry_base_coordinate_and_open_bore(self):
        p=dict(dx=.0005,traverse=0.,omega=0.,layer=-.009,tool_top=.02,
            tool_radius=.02,bore_width=.008,substrate_depth=.01,base_temperature=293.15,tool_temperature=293.15)
        x=np.array([[0.,0.,-.00976],[0.,0.,-.00825]])
        r=contact_rows(x,np.array([11,12]),0.,p)
        np.testing.assert_array_equal(r['particle_id'],[11]);np.testing.assert_array_equal(r['feature_id'],[0])
        np.testing.assert_allclose(r['normal'],[[0,0,1]],rtol=0,atol=0)
        self.assertAlmostEqual(r['gap'][0],.00024,delta=1e-17)
        np.testing.assert_allclose(r['contact_point'],[[0,0,-.01]],rtol=0,atol=1e-17)

    def test_12_initial_energy_and_momentum_literal(self):
        mass=2700*.0005**3
        self.assertAlmostEqual(.5*16*mass*.00212**2,1.213488e-11,delta=1e-25)
        self.assertAlmostEqual(-16*mass*.00212,-1.1448e-8,delta=1e-23)


if __name__=='__main__':unittest.main(verbosity=2)

"""Independent reader references and hand-written impulse ledger; no solver imports."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from inspect_normal_contact import analytic_impact, normal_at, audit_normal_interval


def literal_interval():
    # m=1, dt=1/8, Fn=2.5, mu=.2, static wall, v0=(1,0,0).
    # Jt=(-1/16,0,0), Jn=(0,0,-5/16), v1=(15/16,0,-5/16).
    # K1-K0=-3/256, Q=31/512, Wn=25/512, exactly K+Q=Wn.
    # Binary-exact fixture respects the component-scale 1e-17 momentum audit;
    # the first decimal fixture's 4.16e-17 subtraction roundoff is archived.
    return dict(dt_s=.125,start_s=0.,v_before=[1.,0.,0.],v_after=[.9375,0.,-.3125],
        normal_impulse=[[0.,0.,-.3125]],tangential_impulse=[[-.0625,0.,0.]],
        constraint_impulse=[0.,0.,0.],tangential_relative_midpoint=[[.96875,0.,0.]],
        friction_heat_raw_J=[.060546875],reference_frame_velocity=[0.,0.,0.],
        relative_speed_scale_m_s=1.,
        friction_rounding_bound_J=[256*np.finfo(float).eps*.0625*(2+np.sqrt(.9765625))],
        other_work_J=0.,constraint_work_J=0.,wall_work_J=0.,delta_kinetic_J=-.01171875,
        normal_work_J=[.048828125])


def inspect(row):
    return audit_normal_interval(row,np.array([1.,0.,0.]),np.zeros(3),[[0.,0.,0.]],
        2.5,[1,1,1],1.,.125,0.,np.array([1.,0.,0.]))


class IndependentReadbackTests(unittest.TestCase):
    def test_01_literal_initial_peak_and_release(self):
        a=analytic_impact(0.,0.); self.assertEqual(float(a['z']),1.); self.assertEqual(float(a['v']),-1.)
        b=analytic_impact(np.pi/4,0.)
        self.assertEqual(float(b['delta']),.5); self.assertEqual(float(b['force']),2.)
        self.assertEqual(float(b['potential']),.5)
        c=analytic_impact(np.pi/2,0.)
        self.assertEqual(float(c['normal_impulse']),2.); self.assertEqual(float(c['v']),1.)
        self.assertEqual(float(c['potential']),0.)

    def test_02_literal_moving_wall_work_and_galilean_offset(self):
        a=analytic_impact(2.,0.); b=analytic_impact(2.,.2)
        self.assertAlmostEqual(float(a['z']),1.4292036732051034)
        self.assertAlmostEqual(float(b['z']),1.8292036732051034)
        self.assertAlmostEqual(float(b['v']),1.2)
        self.assertEqual(float(b['wall_work']),.4)
        self.assertAlmostEqual(.5*float(b['v'])**2-.32,float(b['wall_work']))

    def test_03_literal_spring_energy_and_lift_geometry(self):
        x=np.array([[.01,0.,.00175]])
        a=normal_at(x,0.); b=normal_at(x,6e-7)
        self.assertAlmostEqual(float(a['force'][0]),1.,delta=1e-12)
        self.assertAlmostEqual(float(a['potential'][0]),1/1315000.,delta=1e-18)
        self.assertAlmostEqual(b['layer']-a['layer'],8e-6,delta=1e-18)
        self.assertEqual(b['wall_vz'],20.); self.assertFalse(b['active'][0])
        self.assertEqual(float(b['force'][0]),0.); self.assertEqual(float(b['potential'][0]),0.)

    def test_04_handwritten_normal_tangential_kinetic_identity(self):
        a=inspect(literal_interval())
        self.assertEqual(a['kinetic'],-.01171875)
        self.assertEqual(a['normal'],.048828125)
        self.assertEqual(a['q'],.060546875)

    def test_05_catches_wrong_normal_impulse(self):
        row=literal_interval(); row['normal_impulse'][0][2]=-.24
        with self.assertRaisesRegex(ValueError,'normal impulse'): inspect(row)

    def test_06_catches_hidden_normal_clamp(self):
        row=literal_interval(); row['v_after'][2]=0.; row['constraint_impulse'][2]=.3125
        with self.assertRaisesRegex(ValueError,'free-DOF constraint'): inspect(row)

    def test_07_catches_friction_heat_absolute_or_false_value(self):
        row=literal_interval(); row['friction_heat_raw_J']=[-.060546875]
        with self.assertRaisesRegex(ValueError,'raw signed friction'): inspect(row)

    def test_08_catches_missing_normal_material_work(self):
        row=literal_interval(); row['normal_work_J']=[0.]
        with self.assertRaisesRegex(ValueError,'normal material work'): inspect(row)

    def test_09_negative_reference_time_rejected(self):
        with self.assertRaisesRegex(ValueError,'negative reference time'): analytic_impact(-.1,0.)


if __name__=='__main__': unittest.main(verbosity=2)

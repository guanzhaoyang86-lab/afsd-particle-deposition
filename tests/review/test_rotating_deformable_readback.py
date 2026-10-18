"""SCHEMA/ARITHMETIC_FIXTURE_NOT_PHYSICAL: no solver or trajectory run.

Independent literals: m=3.375e-7 kg, dt=5e-9 s, Ft=.2 N, wall_y=.2 m/s.
Jy=1e-9, vy=2/675, W=2e-10, K=1/675*1e-9, Q=W-K.
"""
import copy
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import inspect_rotating_deformable_block as audit


def interval():
    speed = 2/675
    rounding = 256*np.finfo(float).eps*1e-9*(.2+.2+(.2-speed))
    return dict(dt_s=5e-9,start_s=0.,v_before=[0.,0.,0.],v_after=[0.,speed,0.],
        normal_impulse=[[0.,0.,-5e-9]],tangential_impulse=[[0.,1e-9,0.]],
        constraint_impulse=[0.,0.,5e-9],tangential_relative_midpoint=[[0.,1/675-.2,0.]],
        friction_heat_raw_J=[2e-10-1e-9/675],friction_rounding_bound_J=[rounding],
        reference_frame_velocity=[0.,.2,0.],relative_speed_scale_m_s=.2,
        other_work_J=0.,wall_work_J=2e-10,constraint_work_J=0.,normal_work_J=[0.],
        delta_kinetic_J=1e-9/675,kinetic_work_residual_J=0.)


def inspect_interval(a):
    return audit.audit_interval(a,np.zeros(3),np.zeros(3),np.array([[0.,.2,0.]]),
        [1,1,0],3.375e-7,5e-9,0.,np.zeros(3))


def configuration():
    source = json.loads((ROOT/'configs/benchmarks/rotating_deformable_block.json').read_text(encoding='utf-8'))
    return dict(classification='ENGINEERING_APPROXIMATE_ROTATING_DEFORMABLE_BLOCK_NOT_AFSD',
        production_authorized=False,source_configuration=source,
        parameters={k:v['value'] for k,v in source['parameters'].items()},dt_s=5e-9,
        physical_scaling='NONE',tags_are_observation_only=True,spatial_operator='corrected_work_conjugate',
        normal_law='PRESCRIBED 1 N PER TOP PARTICLE; NORMAL DOF CLAMPED; NOT PENALTY',
        field_status=dict(temperature='PRESCRIBED_INPUT_NOT_SOLVED',Es='DIAGNOSTIC_ONLY_NOT_CERTIFIED_FREE_ENERGY',
            friction_dissipation='SOLVED_RAW_NOT_DEPOSITED_AS_HEAT'))


class IndependentDeformableReadback(unittest.TestCase):
    def test_01_literal_interval_impulse_and_signed_work(self):
        result = inspect_interval(interval())
        self.assertAlmostEqual(result['q'],1.985185185185185e-10,delta=1e-25)
        self.assertAlmostEqual(result['wall'],2e-10,delta=1e-25)
        self.assertAlmostEqual(result['dk'],1.4814814814814816e-12,delta=1e-26)
        np.testing.assert_allclose(result['v_after'],[0.,.002962962962962963,0.],rtol=0.,atol=1e-17)

    def test_02_impulse_and_reaction_tamper_rejected(self):
        for key,value in [('tangential_impulse',[[0.,2e-9,0.]]),('constraint_impulse',[0.,0.,0.])]:
            a=interval();a[key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'actual kick momentum'):
                inspect_interval(a)

    def test_03_saved_raw_friction_cannot_replace_signed_work(self):
        for value in (0.,-1.985185185185185e-10,2e-10):
            a=interval();a['friction_heat_raw_J']=[value]
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'signed friction work'):
                inspect_interval(a)

    def test_04_saved_roundoff_bound_cannot_be_enlarged(self):
        a=interval();a['friction_rounding_bound_J']=[1e-10]
        with self.assertRaisesRegex(ValueError,'original friction rounding bound'):inspect_interval(a)

    def test_05_internal_force_must_be_in_momentum_equation(self):
        with self.assertRaisesRegex(ValueError,'actual kick momentum'):
            audit.audit_interval(interval(),np.zeros(3),np.array([1.,0.,0.]),np.array([[0.,.2,0.]]),
                [1,1,0],3.375e-7,5e-9,0.,np.zeros(3))

    def test_06_time_midpoint_and_nonfinite_rejected(self):
        for key,value in [('dt_s',-5e-9),('start_s',5e-9),
                ('tangential_relative_midpoint',[[0.,-.2,0.]]),('v_after',[0.,float('nan'),0.])]:
            a=interval();a[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):inspect_interval(a)

    def test_07_independent_rest_step_has_free_inner_points(self):
        x=np.array([[.01,0.,.00175],[.01,0.,.00125],[.01,0.,.00025]])
        r=audit.independent_first_step(x,5e-9)
        np.testing.assert_allclose(r['velocity'],[[0.,2/675,0.],[0.,0.,0.],[0.,0.,0.]],rtol=0.,atol=1e-17)
        expected=x.copy();expected[0,1]=5e-9/675
        np.testing.assert_allclose(r['position'],expected,rtol=0.,atol=1e-22)
        self.assertAlmostEqual(r['wall_work'],2e-10,delta=1e-25)
        self.assertAlmostEqual(r['q'],2e-10-1e-9/675,delta=1e-25)

    def test_08_literal_pair_velocity_gradient_uses_neighbor_volumes(self):
        velocity=np.array([[0.,0.,0.],[1.,2.,3.]])
        i,j=np.array([0]),np.array([1]);g=np.array([[2.,0.,0.]])
        B=np.repeat(np.eye(3)[None],2,axis=0)
        L=audit.reconstruct_gradient(velocity,i,j,g,np.array([.25,.5]),B)
        expected=np.array([[[1.,0.,0.],[2.,0.,0.],[3.,0.,0.]],
                           [[.5,0.,0.],[1.,0.,0.],[1.5,0.,0.]]])
        np.testing.assert_array_equal(L,expected)

    def test_09_physical_parameter_changes_are_not_dt_sensitivity(self):
        a=configuration();self.assertEqual(audit.audit_configuration(a),5e-9)
        for key,value in [('rho0',27.),('shear',26.3e7),('normal_force',2.),('dx',.00025)]:
            wrong=copy.deepcopy(a);wrong['parameters'][key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'frozen physical parameters'):
                audit.audit_configuration(wrong)

    def test_10_fake_temperature_or_energy_certification_rejected(self):
        for key,value in [('temperature','SOLVED'),('Es','CERTIFIED_FREE_ENERGY'),
                          ('friction_dissipation','HEAT_RESIDUAL_CORRECTION')]:
            a=configuration();a['field_status'][key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'field honesty policy'):
                audit.audit_configuration(a)

    def test_11_failed_or_wrong_kernel_neighborhood_rejected(self):
        x=np.array([[0.,0.,0.],[.0005,0.,0.]])
        with self.assertRaisesRegex(ValueError,'independent radial kernel derivative'):
            audit.audit_operator(x,np.zeros((2,3)),np.ones(2),np.ones(2),[0],[1],
                np.zeros((1,3)),np.repeat(np.eye(3)[None],2,axis=0),np.zeros((2,3,3)),
                np.zeros((2,3)),np.zeros((2,3,3)))

    def test_12_reader_imports_no_dynamics_reference_generator(self):
        source=(ROOT/'scripts/inspect_rotating_deformable_block.py').read_text(encoding='utf-8')
        for forbidden in ('from afsd_backend','import pysph','import probe_rotating','from probe_rotating'):
            self.assertNotIn(forbidden,source)


if __name__=='__main__':unittest.main(verbosity=2)

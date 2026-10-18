"""Independent literal-input/schedule checks; no mechanical trajectory run.

The references below are literal SI data and elementary volume/count formulas,
not contact/material/integrator results. Array preservation is additionally
checked against the frozen fast fixture, without relabeling its output.
Roundoff allowances: mass2e-19kg/node, schedule2e-20s, volume5e-22m3.
These are wiring tests, NOT validation of slow coupled dynamics or events.
"""
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))

import hht_slow_inputs as tested
import reference_hht_block as old_reference


def evidence(name, **values):
    print('HHT_SLOW_INPUT_EVIDENCE '+json.dumps(dict(name=name,
          trajectory_advanced=False, **values), sort_keys=True, allow_nan=False), flush=True)


class SlowInputs(unittest.TestCase):
    def test_01_literal_three_dimensional_geometry_and_physical_mass(self):
        data = tested.inputs()
        ijk = np.array([(i,j,k) for i in range(3) for j in range(3) for k in range(3)])
        points = ijk*.005
        volume = np.array([(.005**3)*(.5 if i != 1 else 1.)*
            (.5 if j != 1 else 1.)*(.5 if k != 1 else 1.) for i,j,k in ijk])
        np.testing.assert_array_equal(data['position_m'], points)
        np.testing.assert_array_equal(data['reference_position_m'], points)
        np.testing.assert_allclose(data['reference_volume_m3'], volume, rtol=0., atol=5e-22)
        np.testing.assert_allclose(data['mass_kg'], 2700.*volume, rtol=0., atol=2e-19)
        np.testing.assert_array_equal(data['particle_id'], np.arange(1,28))
        np.testing.assert_array_equal(data['smoothing_length_m'], np.full(27,.005*1.3))
        self.assertAlmostEqual(float(data['mass_kg'].sum()), .0027, delta=2e-18)
        self.assertAlmostEqual(float(data['reference_volume_m3'].sum()), 1e-6, delta=2e-21)
        for axis in range(3): self.assertEqual(len(np.unique(points[:,axis])), 3)
        self.assertEqual(data['position_m'].shape, (27,3))
        evidence('literal_geometry', particle_count=27, material_coordinate_dofs=81,
                 mass_kg=float(data['mass_kg'].sum()), size_m=.01)

    def test_02_literal_nonzero_startup_and_unscaled_material(self):
        data = tested.inputs(); p = data['parameters']
        expected_v = np.zeros((27,3)); expected_v[:,2] = .02
        expected_S = np.tile(np.diag([20000.,-10000.,-10000.]), (27,1,1))
        expected_T = np.repeat([400.,401.,402.], 9)
        np.testing.assert_array_equal(data['velocity_m_s'], expected_v)
        np.testing.assert_array_equal(data['deviatoric_stress_Pa'], expected_S)
        np.testing.assert_array_equal(data['density_kg_m3'], np.full(27,2700.))
        np.testing.assert_array_equal(data['ep'], np.zeros(27))
        np.testing.assert_array_equal(data['temperature_K'], expected_T)
        expected = dict(rho0=2700.,young=70e9,poisson=0.,shear=35e9,bulk=70e9/3,
            a1=317.6e6,a2=1e6,a3=70.,b1=12.44,b2=5.65,c1=13.,c2=.01,
            room=293.15,melt=855.15,reference_rate=1.,rate_floor=1.)
        for key,value in expected.items(): self.assertEqual(p[key], value, key)
        self.assertGreater(float(np.linalg.norm(expected_S)), 0.)
        self.assertGreater(float(np.linalg.norm(expected_v)), 0.)
        self.assertNotIn('mass_scale',p); self.assertNotIn('time_scale',p)
        self.assertNotIn('modulus_scale',p)
        evidence('unscaled_startup', velocity_z_m_s=.02, stress_diagonal_Pa=[20000,-10000,-10000],
                 young_Pa=70e9, initial_kinetic_J=.5*.0027*.02**2)

    def test_03_literal_contact_cohorts_and_heat_inputs(self):
        data = tested.inputs(); p = data['parameters']
        np.testing.assert_array_equal(data['face_nodes'][0], np.arange(9))
        np.testing.assert_array_equal(data['face_nodes'][1], np.arange(18,27))
        self.assertEqual(set(data['face_nodes'][0]) & set(data['face_nodes'][1]), set())
        expected = dict(dx=.005,normal_preload=2e-6,normal_stroke=2e-7,
            normal_stiffness=7e6,friction=.2,wall_speed=.2,
            conductivity=167.,cp=900.,contact_conductance=1000.,
            heat_partition=.5,plastic_heat_fraction=1.)
        for key,value in expected.items(): self.assertEqual(p[key], value, key)
        np.testing.assert_array_equal(data['wall_temperature_K'], [410.,420.])
        self.assertAlmostEqual(p['normal_stiffness']*p['normal_preload'],14.,delta=2e-15)
        self.assertEqual(p['contact_conductance']*p['dx']**2,.025)
        evidence('contact_and_heat', contacts=18, initial_normal_force_N_each=14.,
                 reservoir_conductance_W_K_each=.025, original_physical_parameters=True)

    def test_04_all_old_physical_arrays_preserved_and_detached(self):
        slow, old = tested.inputs(), old_reference.fixture()
        for key,value in old.items():
            if isinstance(value,np.ndarray):
                np.testing.assert_array_equal(slow[key],value,err_msg=key)
                self.assertEqual(slow[key].dtype,value.dtype,key)
            elif key == 'face_nodes':
                for actual, expected in zip(slow[key],value):
                    np.testing.assert_array_equal(actual,expected)
            elif key not in ('classification','parameters'):
                self.assertEqual(slow[key],value,key)
        for key,value in old['parameters'].items():
            if key != 'end_time': self.assertEqual(slow['parameters'][key],value,key)
        allowed_new = {'desired_end_time_s','normal_frequency_Hz','wall_tangential_program',
                       'maximum_translation_m','maximum_pair_relative_change'}
        self.assertEqual(set(slow['parameters'])-set(old['parameters']),allowed_new)
        slow['position_m'][0,0] = 123.
        slow['velocity_m_s'][0,2] = 456.
        slow['face_nodes'][0][0] = 26
        slow['parameters']['young'] = 1.
        new = tested.inputs()
        np.testing.assert_array_equal(new['position_m'],old['position_m'])
        np.testing.assert_array_equal(new['velocity_m_s'],old['velocity_m_s'])
        np.testing.assert_array_equal(new['face_nodes'][0],np.arange(9))
        self.assertEqual(new['parameters']['young'],70e9)
        self.assertEqual(old_reference.fixture()['parameters']['end_time'],2e-6)
        evidence('frozen_physical_arrays', unchanged=True, independent_copy=True,
                 old_fixture_time_unchanged_s=2e-6)

    def test_05_distinct_rigid_program_not_prescribed_free_motion(self):
        data = tested.inputs(); p = data['parameters']
        self.assertEqual(data['classification'],'LOCAL_SLOW_LOAD_FULL_STARTUP_FEASIBILITY_NOT_AFSD')
        self.assertFalse(data['production_authorized'])
        self.assertEqual(p['normal_frequency_Hz'],5.)
        self.assertEqual(p['wall_tangential_program'],'CONSTANT_PLUS_Z')
        self.assertEqual(p['end_time'],20e-6)
        self.assertEqual(p['desired_end_time_s'],1e-3)
        self.assertEqual(p['maximum_translation_m'],250e-6)
        self.assertEqual(p['maximum_pair_relative_change'],.005)
        self.assertIn('MANUFACTURED_SLOW_LOADING_NOT_AFSD_TOOL',data['program_source'])
        for mapping in (data,p):
            for forbidden in ('prescribed_velocity','prescribed_displacement','counterforce',
                              'actuator','cancel_internal_force','tracer_force','AS_force','RS_force'):
                self.assertNotIn(forbidden,mapping)
        # Inputs alone cannot prove how the future driver uses free DOFs.
        evidence('rigid_program_input_only', normal_frequency_Hz=5., wall_tangent_m_s=.2,
                 free_material_motion_driver_validation='NOT_RUN',not_afsd_tool=True)

    def test_06_all_registered_schedules_and_every_interval(self):
        for level,total in ((0,28),(1,56),(2,112)):
            result = np.asarray(tested.schedule(level))
            self.assertEqual(result.shape,(total,2))
            self.assertEqual(result[0,0],0.)
            self.assertEqual(result[-1,1],20e-6)
            self.assertTrue(np.all(result[:,1]>result[:,0]))
            np.testing.assert_allclose(result[1:,0],result[:-1,1],atol=2e-20,rtol=0.)
            counts=np.array([8,4,16])*2**level
            expected_dt=np.repeat(np.array([250e-9,500e-9,1000e-9])/2**level,counts)
            np.testing.assert_allclose(result[:,1]-result[:,0],expected_dt,atol=2e-20,rtol=0.)
            self.assertEqual(result[counts[0]-1,1],2e-6)
            self.assertEqual(result[counts[:2].sum()-1,1],4e-6)
            self.assertAlmostEqual(float(np.sum(result[:,1]-result[:,0])),20e-6,delta=2e-20)
            evidence('schedule', level=level, intervals=total, final_time_s=float(result[-1,1]),
                     interval_counts=counts.tolist(), full_startup_included=True)

    def test_07_schedule_rejects_unregistered_refinements(self):
        for value in (-1,3,4,False,True,0.,1.,'0',None,np.int64(0)):
            with self.subTest(value=str(value)):
                with self.assertRaises(ValueError): tested.schedule(value)

    def test_08_reference_grids_and_snapshot_times_not_stretched(self):
        self.assertEqual(tested.REFERENCE_STEPS_S,(7.8125e-9,3.90625e-9))
        self.assertEqual(tested.SNAPSHOT_STEP_S,62.5e-9)
        self.assertEqual(tested.SCOUT_END_S,20e-6)
        self.assertEqual(tested.DESIRED_END_S,1e-3)
        for dt,count in zip(tested.REFERENCE_STEPS_S,(2560,5120)):
            self.assertAlmostEqual(tested.SCOUT_END_S/dt,count,delta=1e-12)
        self.assertAlmostEqual(tested.SCOUT_END_S/tested.SNAPSHOT_STEP_S,320.,delta=1e-12)
        evidence('reference_grids', reference_dt_s=list(tested.REFERENCE_STEPS_S),
                 reference_steps=[2560,5120],common_snapshots_including_initial=321,
                 one_ms_reference_status='NOT_RUN')


if __name__ == '__main__':
    unittest.main(verbosity=2)

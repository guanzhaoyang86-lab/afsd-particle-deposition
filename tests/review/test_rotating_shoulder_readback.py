"""SCHEMA/ARITHMETIC_FIXTURE_NOT_PHYSICAL: no solver or reference integration.

Literal first-step values are from m=.001, Ft=.2, dt=2e-5, wall_y=.2.
Mutations must be rejected by the independent saved-array reader.
"""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import inspect_rotating_shoulder_component as audit


def interval():
    return dict(dt_s=2e-5, start_s=0., v_before=[0., 0., 0.], v_after=[0., .004, 0.],
        normal_impulse=[[0., 0., -2e-5]], tangential_impulse=[[0., 4e-6, 0.]],
        constraint_impulse=[0., 0., 2e-5], tangential_relative_midpoint=[[0., -.198, 0.]],
        friction_heat_raw_J=[7.92e-7], friction_rounding_bound_J=[1.3551471056416633e-19],
        relative_speed_scale_m_s=.2, reference_frame_velocity=[0., .2, 0.],
        wall_work_J=8e-7, constraint_work_J=0., normal_work_J=[0.], other_work_J=0.,
        delta_kinetic_J=8e-9, kinetic_work_residual_J=0.)


def trace():
    a = interval()
    result = dict(velocity=[0., .004, 0.], constraint_projector=[[1, 0, 0], [0, 1, 0], [0, 0, 0]],
        intervals=[a], friction_heat_J=[7.92e-7], normal_impulse=[[0., 0., -2e-5]],
        tangential_impulse=[[0., 4e-6, 0.]], constraint_impulse=[0., 0., 2e-5],
        wall_work_J=8e-7, constraint_work_J=0., other_work_J=0., delta_kinetic_J=8e-9,
        kinetic_work_residual_J=0.)
    row = dict(particle_id=101, feature_id=1, contact_id='101:1',
        application_point=[.01, 0., .0013], contact_point=[.01, 0., .0015],
        normal=[0., 0., -1.], body_origin=[0., 0., .0015], translation_velocity=[0., 0., 0.],
        omega=[0., 0., 20.], gap=.0002, penetration=.00005,
        wall_velocity=[0., .2, 0.], wall_velocity_at_face=[0., .2, 0.])
    return dict(particle_id=101, feature_id=1, step=1, stage='MIDPOINT_FROZEN_GEOMETRY_KICK',
        time_mid_s=1e-5, position_m=[.01, 0., .0013], contact_row=row,
        material_velocity_before_m_s=[0., 0., 0.], result=result)


def audit_fixture(a):
    return audit.audit_trace(a, np.array([.01, 0., .0013]), np.zeros(3),
        np.array([.01, 4e-8, .0013]), np.array([0., .004, 0.]), 2e-5, 1, 101)


class IndependentReadbackFixtures(unittest.TestCase):
    def test_01_literal_subinterval_work_and_impulses(self):
        result = audit.audit_interval(interval(), np.zeros(3), np.array([0., .2, 0.]), 2e-5, 0.)
        self.assertAlmostEqual(result['q'], 7.92e-7, delta=1e-20)
        self.assertAlmostEqual(result['wall'], 8e-7, delta=1e-20)
        self.assertAlmostEqual(result['dk'], 8e-9, delta=1e-20)

    def test_02_tangent_impulse_tamper_rejected(self):
        a = interval(); a['tangential_impulse'][0][1] = 5e-6
        with self.assertRaisesRegex(ValueError, 'impulse momentum'):
            audit.audit_interval(a, np.zeros(3), np.array([0., .2, 0.]), 2e-5, 0.)

    def test_03_normal_reaction_tamper_rejected(self):
        a = interval(); a['constraint_impulse'][2] = 0.
        with self.assertRaisesRegex(ValueError, 'actuator reaction'):
            audit.audit_interval(a, np.zeros(3), np.array([0., .2, 0.]), 2e-5, 0.)

    def test_04_negative_signed_value_cannot_be_clipped(self):
        a = interval()
        a.update(v_before=[1e-3, 0., 0.], v_after=[1e-3, 0., 0.],
            tangential_impulse=[[1e-24, 0., 0.]], tangential_relative_midpoint=[[1e-3, 0., 0.]],
            friction_heat_raw_J=[-1e-27], friction_rounding_bound_J=[2e-27],
            wall_work_J=0., delta_kinetic_J=0.)
        result = audit.audit_interval(a, np.array([1e-3, 0., 0.]), np.zeros(3), 2e-5, 0.)
        self.assertEqual(result['negative_raw'], 1)
        for wrong in (0., 1e-27):
            a['friction_heat_raw_J'] = [wrong]
            with self.assertRaisesRegex(ValueError, 'clipped or made absolute'):
                audit.audit_interval(a, np.array([1e-3, 0., 0.]), np.zeros(3), 2e-5, 0.)

    def test_05_interval_time_and_nonfinite_rejected(self):
        for key, value in [('dt_s', -2e-5), ('start_s', 1e-5), ('v_after', [0., float('nan'), 0.])]:
            a = interval(); a[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                audit.audit_interval(a, np.zeros(3), np.array([0., .2, 0.]), 2e-5, 0.)

    def test_06_relative_velocity_and_friction_heat_checked_independently(self):
        for key, value in [('tangential_relative_midpoint', [[0., -.2, 0.]]), ('friction_heat_raw_J', [8e-7])]:
            a = interval(); a[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                audit.audit_interval(a, np.zeros(3), np.array([0., .2, 0.]), 2e-5, 0.)

    def test_07_literal_trace_angular_impulse(self):
        result = audit_fixture(trace())
        self.assertAlmostEqual(result['torque_impulse_z'], 4e-8, delta=1e-21)
        self.assertAlmostEqual(20*result['torque_impulse_z'], result['wall'], delta=1e-20)

    def test_08_feature_and_identity_swap_rejected(self):
        for field, value in [('feature_id', 3), ('particle_id', 102), ('step', 2)]:
            a = trace(); a[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError): audit_fixture(a)

    def test_09_uniform_translation_cannot_replace_rotation(self):
        a = trace(); a['contact_row']['wall_velocity'] = [.2, 0., 0.]
        with self.assertRaisesRegex(ValueError, 'rotating wall velocity'): audit_fixture(a)

    def test_10_wrong_contact_point_or_drift_rejected(self):
        a = trace(); a['contact_row']['contact_point'][2] = .0013
        with self.assertRaisesRegex(ValueError, 'projection'): audit_fixture(a)
        b = trace(); b['position_m'][1] = 4e-8
        with self.assertRaisesRegex(ValueError, 'first drift'): audit_fixture(b)

    def test_11_missing_interval_and_duration_mismatch_rejected(self):
        a = trace(); a['result']['intervals'] = []
        with self.assertRaisesRegex(ValueError, 'empty accepted'): audit_fixture(a)
        b = trace(); b['result']['intervals'].append(copy.deepcopy(interval()))
        with self.assertRaises(ValueError): audit_fixture(b)

    def test_12_reference_refinement_limits_are_not_fitted(self):
        a = dict(position_m=[[0., 0., 0.]]*4, velocity_m_s=[[0., 0., 0.]]*4,
                 friction_raw_J=[1e-5]*4, wall_work_J=[2e-5]*4)
        audit.reference_refinement(a, copy.deepcopy(a))
        b = copy.deepcopy(a); b['position_m'][0] = [2e-11, 0., 0.]
        with self.assertRaisesRegex(ValueError, 'reference refinement'): audit.reference_refinement(a, b)

    def test_13_component_reference_vector_error_and_energy_limit(self):
        reference = dict(position_m=[[0., 0., 0.]]*4, velocity_m_s=[[0., 0., 0.]]*4,
                         friction_raw_J=[1e-5]*4, wall_work_J=[2e-5]*4)
        run = dict(final_position_m=[[0., 0., 0.]]*4, final_velocity_m_s=[[0., 0., 0.]]*4,
                   friction_raw_J=4e-5, wall_work_J=8e-5)
        audit.compare_to_reference(run, reference)
        for key, value in [('final_position_m', [[2e-8, 0., 0.]]*4), ('wall_work_J', 9e-5)]:
            wrong = copy.deepcopy(run); wrong[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'component/reference'):
                audit.compare_to_reference(wrong, reference)

    def test_14_incomplete_physical_run_not_compared(self):
        with patch.object(audit, 'load', return_value=dict(status='BLOCKED', configuration={})):
            with self.assertRaisesRegex(ValueError, 'failed endpoint'): audit.audit_run('SCHEMA_ONLY', 2e-5)

    def test_15_reader_import_does_not_load_solver_or_reference(self):
        source = (ROOT/'scripts/inspect_rotating_shoulder_component.py').read_text(encoding='utf-8')
        for forbidden in ('import probe_rotating', 'from afsd_backend', 'import pysph', 'reference_rk4('):
            self.assertNotIn(forbidden, source)

    def test_16_sliding_law_checks_cap_and_direction_not_just_positive_heat(self):
        a = interval(); audit.audit_sliding_friction(a, np.zeros(3), np.array([0., .2, 0.]))
        for wrong in ([[0., 2e-6, 0.]], [[4e-6, 0., 0.]]):
            a['tangential_impulse'] = wrong
            with self.assertRaisesRegex(ValueError, 'Coulomb sliding impulse'):
                audit.audit_sliding_friction(a, np.zeros(3), np.array([0., .2, 0.]))

    def test_17_saved_roundoff_bound_cannot_be_arbitrarily_enlarged(self):
        a = interval(); a['friction_rounding_bound_J'] = [1e-10]
        with self.assertRaisesRegex(ValueError, 'original roundoff bound'):
            audit.audit_sliding_friction(a, np.zeros(3), np.array([0., .2, 0.]))


if __name__ == '__main__': unittest.main(verbosity=2)

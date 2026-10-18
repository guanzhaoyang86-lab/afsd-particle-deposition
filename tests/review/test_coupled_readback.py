"""Literal and tamper tests for independent saved-array checks, no physics runs."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import inspect_coupled_afsd as reader


def plastic_fixture():
    return dict(dt_s=1., mass_kg=[.0027], rho_before_kg_m3=[2700.],
        ep_before=[0.], ep_after=[.001], S_before_Pa=[np.zeros((3,3)).tolist()],
        S_after_Pa=[np.diag([60e6,-30e6,-30e6]).tolist()], drho_kg_m3_s=[0.],
        sigma_mid_Pa=[np.diag([30e6,-15e6,-15e6]).tolist()],
        plastic_work_J=[.09], return_dissipation_J=[.03945])


def thermal_fixture():
    # One particle, no pair conduction: P=.09, etaQ=.01, contact=-.005.
    # Independent sum .095 J; return energy intentionally not an input.
    temperature = 400+.095/2.43
    measured = 2.43*(temperature-400)
    source = dict(plastic_work_J=[.09], plastic_heat_J=[.09], plastic_unthermalized_J=[0.],
        friction_raw_J=[.02], friction_material_heat_J=[.01], friction_wall_heat_J=[.01],
        contact_in_J=[-.005], internal_conduction_J=[0.], total_heat_input_J=[.095],
        measured_thermal_increment_J=[measured], mc_delta_T_residual_J=[measured-.095])
    ledger = {k: v[0] for k, v in source.items()}; ledger['wall_heat_out_J'] = .015
    return dict(dt_s=1., mass_kg=[.0027], rho_before_kg_m3=[2700.], x_mid_m=[[0,0,0]],
        temperature_before_K=[400.], temperature_after_K=[temperature], i=[], j=[], g=[],
        thermal=dict(particle_energy_J=source, ledger=ledger, pair_power_W=[0.], conduction_row_W_K=[0.]))


def interval_fixture():
    # Unit mass, Fn=1, mu=.3, h=.1: Jt=-.03 ex, Jn=.1 ez.
    # v0=(1,0,0), v1=(.97,0,.1), raw Q=.02955, normal work=.005 J.
    vb = np.array([.97,0,.1]); bound = 256*np.finfo(float).eps*.03*(2+np.linalg.norm(vb))
    return dict(dt_s=.1, v_before=[1,0,0], v_after=vb.tolist(), normal_impulse=[[0,0,.1]],
        tangential_impulse=[[-.03,0,0]], constraint_impulse=[0,0,0],
        friction_heat_raw_J=[.02955], reference_frame_velocity=[0,0,0],
        relative_speed_scale_m_s=1., friction_rounding_bound_J=[bound],
        other_work_J=0., constraint_work_J=0., wall_work_J=0., delta_kinetic_J=-.02455,
        normal_work_J=[.005])


def short_window_fixture():
    # Literal measured fixture: dt=10ns, slowest=4s, 70 steps at 700ns.
    # 1.5*4*70+30=450s; the next 100ns quantum requires 80 steps/510s.
    config = dict(registered_end_times_s=[2e-6,4e-6], forecast_launch_limit_seconds=450,
        forecast_multiplier=1.5, forecast_overhead_seconds=30, internal_wall_seconds=540,
        external_wall_seconds=570)
    original = dict(parameters=dict(rho0=2700., shear=26.3e9), explicit_configuration=config)
    probe = dict(status='COMPLETE_COST_PROBE', accepted_steps=5, accepted_dt_s=[1e-8]*5,
        step_wall_seconds=[2.,3.,4.,3.,2.], configuration_sha256=reader.fingerprint(original))
    choice = dict(end_time_s=7e-7, estimated_steps=70, forecast_wall_seconds=450., permitted=True)
    override = dict(classification='AFFORDABLE_SHORT_WINDOW_NOT_ORIGINAL_2US_COMPLETION',
        contract='docs/COUPLED_AFSD_SHORT_WINDOW.md', choice=choice,
        origin_probe_manifest_sha256='literal_receipt_hash', physical_parameters_unchanged=True)
    resolved = dict(original, explicit_configuration=dict(config,registered_end_times_s=[7e-7]),
        original_explicit_configuration=config, short_window_budget_override=override)
    return resolved, original, probe


class IndependentCoupledReadback(unittest.TestCase):
    def thermal(self, journal):
        return reader.audit_thermal(journal, np.array([.09]), np.array([.02]), np.array([-.005]), np.array([.00025]))

    def interval(self, interval):
        return reader.audit_interval(interval, np.array([1.,0,0]), 1., np.zeros(3),
            np.array([[0.,0,1]]), np.array([1.]), np.zeros((1,3)), np.ones(3))

    def test_independent_j2_power_and_separate_return(self):
        dep, plastic, returned = reader.audit_plastic(plastic_fixture())
        self.assertAlmostEqual(dep[0], .001, places=15)
        self.assertAlmostEqual(plastic[0], .09, places=14)
        self.assertAlmostEqual(returned[0], .03945, places=14)

    def test_return_dissipation_cannot_be_folded_into_plastic_heat(self):
        case = plastic_fixture(); case['plastic_work_J'][0] += .03945
        with self.assertRaisesRegex(ValueError, 'qnew dep V'): reader.audit_plastic(case)

    def test_negative_ep_increment_rejected_not_clipped(self):
        case = plastic_fixture(); case['ep_after'] = [-.001]
        with self.assertRaisesRegex(ValueError, 'negative plastic'): reader.audit_plastic(case)

    def test_literal_thermal_partition_and_signed_contact(self):
        result = self.thermal(thermal_fixture())
        self.assertAlmostEqual(result['plastic_heat_J'], .09)
        self.assertAlmostEqual(result['friction_material_heat_J'], .01)
        self.assertAlmostEqual(result['contact_in_J'], -.005)

    def test_temperature_tampering_is_not_hidden_by_matching_saved_ledger(self):
        case = thermal_fixture(); case['temperature_after_K'][0] += 1.
        with self.assertRaisesRegex(ValueError, 'mcDeltaT'): self.thermal(case)

    def test_clipping_signed_contact_heat_caught(self):
        case = thermal_fixture(); case['thermal']['particle_energy_J']['contact_in_J'][0] = 0.
        with self.assertRaisesRegex(ValueError, 'contact_in_J'): self.thermal(case)

    def test_double_friction_heating_caught(self):
        case = thermal_fixture(); case['temperature_after_K'][0] += .01/2.43
        with self.assertRaisesRegex(ValueError, 'mcDeltaT'): self.thermal(case)

    def test_literal_friction_and_normal_work(self):
        result = self.interval(interval_fixture())
        self.assertAlmostEqual(result[-1]['friction_raw_J'], .02955, places=14)
        self.assertAlmostEqual(result[-1]['normal_material_work_J'], .005, places=14)
        self.assertAlmostEqual(result[-1]['delta_kinetic_J'], -.02455, places=14)

    def test_posthoc_velocity_edit_caught_by_momentum(self):
        case = interval_fixture(); case['v_after'][1] += .01
        with self.assertRaisesRegex(ValueError, 'kick momentum'): self.interval(case)

    def test_reaction_on_free_axis_caught(self):
        case = interval_fixture(); case['constraint_impulse'] = [0,.01,0]; case['v_after'][1] += .01
        with self.assertRaisesRegex(ValueError, 'reaction on free'): self.interval(case)

    def test_signed_friction_tampering_caught(self):
        case = interval_fixture(); case['friction_heat_raw_J'][0] *= -1
        with self.assertRaisesRegex(ValueError, 'signed actual'): self.interval(case)

    def test_geometry_is_342_distinct_3d_centres_and_25_original_controls(self):
        x = reader.initial_geometry()
        self.assertEqual(x.shape, (342,3)); self.assertEqual(len(np.unique(x, axis=0)), 342)
        self.assertEqual(np.linalg.matrix_rank(x-x.mean(axis=0)), 3)
        self.assertEqual(np.count_nonzero(x[:,2] < 0), 192)
        self.assertEqual(np.count_nonzero(x[:,2] > .0025), 25)

    def test_open_bore_and_finite_shoulder_not_solid_cap(self):
        x = np.array([[0,0,.00151], [.0014,0,.00151], [.0014,0,.003]])
        faces = reader.reference_faces(x, [1,2,3], 0.)
        self.assertFalse(any(key.startswith('1:') for key in faces))
        self.assertIn('2:1', faces)
        self.assertFalse(any(key.startswith('3:') for key in faces))

    def test_planar_neighbor_inventory_cannot_pass_3d_operator(self):
        # Correct literal radial derivative first, so failure is specifically
        # the rank-one moment, not an earlier malformed-gradient shortcut.
        gx = -15/(62*np.pi*.00065**4)*(3*(10/13)**2-6)
        case = dict(x_mid_m=[[0,0,0],[.0005,0,0]], i=[0], j=[1],
            g=[[gx,0,0]], mass_kg=[reader.MASS]*2, rho_before_kg_m3=[2700.]*2,
            force_density_kg_m3=[2700.]*2)
        with self.assertRaisesRegex(ValueError, 'non-3D/ill-conditioned'): reader.audit_spatial(case)

    def test_short_budget_independent_quantized_literal(self):
        _, _, probe = short_window_fixture()
        actual = reader.independent_short_choice(probe)
        self.assertEqual(actual, dict(end_time_s=7e-7, estimated_steps=70, forecast_wall_seconds=450., permitted=True))

    def test_short_budget_does_not_lower_safety_when_nothing_affordable(self):
        _, _, probe = short_window_fixture(); probe['step_wall_seconds'] = [2000.]*5
        with self.assertRaisesRegex(ValueError, 'no affordable'): reader.independent_short_choice(probe)

    def test_short_window_only_schedule_changes(self):
        resolved, original, probe = short_window_fixture()
        self.assertEqual(reader.audit_short_configuration(resolved,original,probe,'literal_receipt_hash')['estimated_steps'], 70)

    def test_short_window_hidden_material_change_is_rejected(self):
        resolved, original, probe = short_window_fixture(); resolved = copy.deepcopy(resolved)
        resolved['parameters']['shear'] *= .1
        with self.assertRaisesRegex(ValueError, 'any physical parameter'):
            reader.audit_short_configuration(resolved,original,probe,'literal_receipt_hash')

    def test_short_window_outcome_based_smaller_horizon_is_rejected(self):
        resolved, original, probe = short_window_fixture()
        resolved['short_window_budget_override']['choice'] = dict(end_time_s=6e-7, estimated_steps=60, forecast_wall_seconds=390., permitted=True)
        with self.assertRaisesRegex(ValueError, 'longest affordable'):
            reader.audit_short_configuration(resolved,original,probe,'literal_receipt_hash')


if __name__ == '__main__': unittest.main(verbosity=2)

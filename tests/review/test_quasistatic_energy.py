"""Independent hand-calculated work identities; not AFSD validation."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))

from afsd_backend.pysph.quasistatic_energy import material_energy_breakdown


def account(S0, S1, H, G=30., K=50., rho1=2., plastic=0., returned=0.):
    return material_energy_breakdown(
        np.asarray(S0)[None], np.asarray(S1)[None], np.array([2.]),
        np.array([rho1]), np.array([6.]), np.asarray(H)[None]/.2,
        .2, G, K, 2., np.array([plastic]), np.array([returned]))


def check_unloaded_elastic_step_keeps_raw_fifty_percent_work_defect():
    G, strain, volume = 30., .01, 3.
    H = np.diag([strain, -strain, 0.])
    result = account(np.zeros((3, 3)), 2*G*H, H)
    storage = 2*volume*G*strain**2
    np.testing.assert_allclose(result['deviatoric_storage_change_J'], storage)
    np.testing.assert_allclose(result['deviatoric_endpoint_work_J'], 2*storage)
    np.testing.assert_allclose(result['material_raw_energy_defect_J'], storage)
    np.testing.assert_allclose(result['be_elastic_deviatoric_dissipation_J'], storage)
    np.testing.assert_allclose(result['material_discrete_identity_error_J'], 0., atol=1e-14)
    assert result['algorithmic_terms_are_heat'] is False
    assert result['raw_acceptance_threshold_modified'] is False


def check_combined_spin_and_stretch_retains_negative_signed_spin_defect():
    G, strain, theta, volume = 30., .01, .3, 3.
    H = np.array([[strain, -theta, 0.], [theta, -strain, 0.], [0., 0., 0.]])
    # A half-angle rotation rotates this traceless tensor's components by theta.
    S1 = 2*G*strain*np.array([[np.cos(theta), np.sin(theta), 0.],
                            [np.sin(theta), -np.cos(theta), 0.], [0., 0., 0.]])
    result = account(np.zeros((3, 3)), S1, H)
    expected_spin = 4*volume*G*strain**2*(np.cos(theta)-1.)
    np.testing.assert_allclose(result['deviatoric_spin_defect_J'], expected_spin)
    assert result['deviatoric_spin_defect_J'] < 0.
    np.testing.assert_allclose(result['deviatoric_discrete_identity_error_J'], 0., atol=1e-14)
    assert abs(result['material_raw_energy_defect_J']) > 1e-14


def check_density_be_work_differs_from_true_volume_increment_by_signed_term():
    K, volume0, dilation = 50., 3., -.04
    result = account(np.zeros((3, 3)), np.zeros((3, 3)),
                     np.eye(3)*dilation/3., rho1=2./(1+dilation))
    storage = K*volume0*(dilation-np.log1p(dilation))
    defect = K*volume0*dilation**3/(1+dilation)
    BE = K*volume0*(np.log1p(dilation)-dilation/(1+dilation))
    np.testing.assert_allclose(result['volumetric_endpoint_work_J'], K*volume0*dilation**2)
    np.testing.assert_allclose(result['volumetric_storage_change_J'], storage)
    np.testing.assert_allclose(result['volumetric_density_be_defect_J'], defect)
    np.testing.assert_allclose(result['be_volumetric_dissipation_J'], BE)
    assert result['volumetric_density_be_defect_J'] < 0.
    np.testing.assert_allclose(result['density_closure_relative_max'], 0., atol=1e-14)
    np.testing.assert_allclose(result['density_be_work_identity_error_J'], 0., atol=1e-14)


def check_known_radial_return_diagnostic_is_not_extra_endpoint_dissipation():
    G, strain, retained, volume = 30., .2, .4, 3.
    H = np.diag([strain, -strain, 0.])
    trial = 2*G*H
    qtrial = 2*np.sqrt(3)*G*strain
    q1 = retained*qtrial
    dep = (qtrial-q1)/(3*G)
    result = account(np.zeros((3, 3)), retained*trial, H,
                     plastic=q1*dep, returned=1.5*G*dep**2)
    stored = 2*volume*G*strain**2*retained**2
    np.testing.assert_allclose(result['return_map_identity_error_J'], 0., atol=1e-13)
    np.testing.assert_allclose(result['material_raw_energy_defect_J'], stored)
    np.testing.assert_allclose(result['be_elastic_deviatoric_dissipation_J'], stored)
    np.testing.assert_allclose(result['return_map_diagnostic_J'], volume*1.5*G*dep**2)
    assert abs(result['return_map_diagnostic_J']-stored) > 1e-13
    np.testing.assert_allclose(result['material_discrete_identity_error_J'], 0., atol=1e-13)


class QuasistaticEnergyTests(unittest.TestCase):
    def test_unloaded_elastic_step_keeps_raw_fifty_percent_work_defect(self):
        check_unloaded_elastic_step_keeps_raw_fifty_percent_work_defect()

    def test_combined_spin_and_stretch_retains_negative_signed_spin_defect(self):
        check_combined_spin_and_stretch_retains_negative_signed_spin_defect()

    def test_density_be_work_differs_from_true_volume_increment_by_signed_term(self):
        check_density_be_work_differs_from_true_volume_increment_by_signed_term()

    def test_known_radial_return_diagnostic_is_not_extra_endpoint_dissipation(self):
        check_known_radial_return_diagnostic_is_not_extra_endpoint_dissipation()


if __name__ == '__main__':
    unittest.main()

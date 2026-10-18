"""Independent literal references for a metered existing material/thermal chain."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.coupled_thermal import advance_thermal
from afsd_backend.pysph.material import update


class ConstantYield:
    """Manufactured material-point oracle input, NOT AA6061/FKS replacement."""
    def yield_stress(self, ep, rate, temperature):
        return np.full(np.broadcast(ep, rate, temperature).shape, 90e6)


def inputs():
    return dict(mass=np.array([2., 3.]), rho=np.array([4., 5.]),
        T=np.array([400., 300.]), conductivity=np.array([6., 8.]),
        x=np.array([[0., 0., 0.], [.25, 0., 0.]]),
        i=np.array([0], dtype=int), j=np.array([1], dtype=int), g=np.array([[2., 0., 0.]]),
        dt=.01, cp=np.array([10., 20.]), plastic_work_J=np.array([10., 0.]),
        friction_raw_J=np.array([4., -1e-12]), contact_power_W=np.array([-5., 7.]),
        wall_row_W_K=np.array([1., 2.]), chi=.9, partition=.25, room=293.15, melt=855.15)


class CoupledThermalTests(unittest.TestCase):
    def test_01_known_j2_return_heats_only_chi_plastic_work(self):
        strain = np.array([np.diag([.002, -.001, -.001])])
        s, ep, plastic_density, returned, _ = update(np.zeros((1, 3, 3)), np.zeros(1),
            strain, np.array([400.]), 1., 30e9, ConstantYield())
        # q_trial=180MPa, dep=.001, q_new=90MPa, V=1e-6m3.
        np.testing.assert_allclose(s[0], np.diag([60e6, -30e6, -30e6]), rtol=0, atol=2e-6)
        self.assertAlmostEqual(ep[0], .001, delta=1e-16)
        self.assertAlmostEqual(plastic_density[0]*1e-6, .09, delta=1e-14)
        self.assertAlmostEqual(returned[0]*1e-6, .045, delta=1e-14)
        a = inputs(); a.update(mass=np.array([.0027]), rho=np.array([2700.]),
            T=np.array([400.]), conductivity=167., x=np.zeros((1, 3)),
            i=np.array([], dtype=int), j=np.array([], dtype=int), g=np.zeros((0, 3)),
            dt=1., cp=900., plastic_work_J=plastic_density*1e-6,
            friction_raw_J=np.zeros(1), contact_power_W=np.zeros(1), wall_row_W_K=np.zeros(1))
        r = advance_thermal(**a)
        self.assertAlmostEqual(r['temperature'][0], 400.+1./30., delta=1e-12)
        self.assertAlmostEqual(r['ledger']['plastic_heat_J'], .081, delta=1e-14)
        self.assertAlmostEqual(r['ledger']['plastic_unthermalized_J'], .009, delta=1e-14)
        self.assertNotIn('return_dissipation_J', r['particle_energy_J'])

    def test_02_independent_mixed_sources_signed_contact_and_conduction(self):
        r = advance_thermal(**inputs())
        # Literal old-operator coefficient C=1152/35 W/K, powers -/+100C.
        np.testing.assert_allclose(r['pair_power_W'], [-115200/35, 115200/35], rtol=1e-14)
        expected = np.array([9.+1.-.05-1152/35, -.25e-12+.07+1152/35])
        np.testing.assert_allclose(r['particle_energy_J']['total_heat_input_J'], expected, rtol=0, atol=1e-13)
        np.testing.assert_allclose(r['temperature'], [400.+expected[0]/20., 300.+expected[1]/60.], rtol=0, atol=1e-12)
        self.assertAlmostEqual(r['ledger']['contact_in_J'], .02, delta=1e-16)
        self.assertEqual(r['ledger']['internal_conduction_J'], 0.)
        self.assertAlmostEqual(r['ledger']['total_heat_input_J'], 10.02-.25e-12, delta=1e-13)
        self.assertAlmostEqual(r['ledger']['wall_heat_out_J'], 2.98-.75e-12, delta=1e-13)
        self.assertLess(abs(r['ledger']['mc_delta_T_residual_J']), 1e-11)

    def test_03_raw_negative_friction_remains_signed_without_clipping(self):
        a = inputs(); a['friction_raw_J'] = np.array([-4., -2.])
        r = advance_thermal(**a)
        np.testing.assert_array_equal(r['particle_energy_J']['friction_raw_J'], [-4., -2.])
        np.testing.assert_array_equal(r['particle_energy_J']['friction_material_heat_J'], [-1., -.5])
        self.assertEqual(r['ledger']['friction_raw_J'], -6.)

    def test_04_quarter_capacity_bound_includes_wall_row(self):
        a = inputs(); r = advance_thermal(**a)
        expected = .25*20./(1152/35+1.)
        self.assertAlmostEqual(r['thermal_dt_bound_s'], expected, delta=1e-16)
        a['dt'] = expected*(1.+1e-10)
        with self.assertRaisesRegex(ValueError, 'time-step bound'): advance_thermal(**a)

    def test_05_wall_row_alone_limits_step(self):
        a = inputs(); a.update(i=np.array([], int), j=np.array([], int), g=np.zeros((0, 3)))
        a['wall_row_W_K'] = np.array([100., 0.])
        self.assertAlmostEqual(advance_thermal(**a)['thermal_dt_bound_s'], .05, delta=1e-16)

    def test_06_success_and_failure_do_not_mutate_inputs(self):
        a = inputs(); before = {k: v.copy() for k, v in a.items() if isinstance(v, np.ndarray)}
        result = advance_thermal(**a)
        result['temperature'][:] = 0.
        for key, value in before.items(): np.testing.assert_array_equal(a[key], value)
        a['plastic_work_J'][:] = 1e20
        before = {k: v.copy() for k, v in a.items() if isinstance(v, np.ndarray)}
        with self.assertRaisesRegex(ValueError, 'NOT clipped'): advance_thermal(**a)
        for key, value in before.items(): np.testing.assert_array_equal(a[key], value)

    def test_07_invalid_positive_parameters_and_nonfinite_sources_rejected(self):
        for key in ('mass', 'rho', 'conductivity', 'cp'):
            for value in (0., -1., np.nan, np.inf):
                a = inputs(); a[key] = np.full(2, value)
                with self.subTest(key=key, value=value):
                    with self.assertRaises(ValueError): advance_thermal(**a)
        for key in ('T', 'plastic_work_J', 'friction_raw_J', 'contact_power_W', 'wall_row_W_K'):
            a = inputs(); a[key][0] = np.nan
            with self.subTest(key=key):
                with self.assertRaises(ValueError): advance_thermal(**a)

    def test_08_invalid_partition_and_chi_rejected_without_clamping(self):
        for key in ('chi', 'partition'):
            for value in (-.01, 1.01, np.inf, np.nan, [0.5]):
                a = inputs(); a[key] = value
                with self.subTest(key=key, value=value):
                    with self.assertRaises(ValueError): advance_thermal(**a)

    def test_09_negative_plastic_work_wall_row_and_temperature_domain_fail(self):
        for key, value in (('plastic_work_J', -1.), ('wall_row_W_K', -1.), ('T', 293.14), ('T', 855.15)):
            a = inputs(); a[key][0] = value
            with self.subTest(key=key, value=value):
                with self.assertRaises(ValueError): advance_thermal(**a)
        a = inputs(); a['contact_power_W'][:] = -1e20
        with self.assertRaisesRegex(ValueError, 'NOT clipped'): advance_thermal(**a)

    def test_10_return_and_residual_are_not_accepted_heat_interfaces(self):
        for key in ('return_dissipation_J', 'energy_residual_J'):
            a = inputs(); a[key] = np.ones(2)
            with self.assertRaises(TypeError): advance_thermal(**a)

    def test_11_pair_integrity_and_raw_radial_gradient_guards(self):
        changes = [dict(i=np.array([0]), j=np.array([0])),
            dict(i=np.array([0, 1]), j=np.array([1, 0]), g=np.array([[2., 0., 0.], [-2., 0., 0.]])),
            dict(i=np.array([2])), dict(i=np.array([0.0])),
            dict(g=np.array([[2., 1., 0.]])), dict(g=np.array([[-2., 0., 0.]])),
            dict(x=np.zeros((2, 3)))]
        for change in changes:
            a = inputs(); a.update(change)
            with self.subTest(change=change):
                with self.assertRaises(ValueError): advance_thermal(**a)

    def test_12_uniform_temperature_no_sources_and_no_conduction_rows(self):
        a = inputs(); a['T'][:] = 400.
        for key in ('plastic_work_J', 'friction_raw_J', 'contact_power_W'): a[key][:] = 0.
        r = advance_thermal(**a)
        np.testing.assert_array_equal(r['temperature'], [400., 400.])
        np.testing.assert_array_equal(r['pair_power_W'], [0., 0.])
        a.update(i=np.array([], int), j=np.array([], int), g=np.zeros((0, 3)))
        a['wall_row_W_K'][:] = 0.
        self.assertEqual(advance_thermal(**a)['thermal_dt_bound_s'], float('inf'))


if __name__ == '__main__': unittest.main(verbosity=2)

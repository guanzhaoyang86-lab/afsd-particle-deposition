"""Independent E1 references; no production or AFSD acceptance implied."""
from pathlib import Path
import copy
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from probe_engineering_material import PerfectPlasticBenchmark, assess, load_config, run_history


class EngineeringMaterialE1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config()
        cls.histories = {n: run_history(cls.cfg, n) for n in (16, 32, 64, 17, 33, 65)}

    def test_independent_peak_and_reverse_yield(self):
        # Explicit tensors, not generated from the updater or its invariant.
        a = np.array([np.diag([1., -.5, -.5]), np.diag([-.5, 1., -.5]),
                      np.diag([-.5, -.5, 1.]),
                      [[-1/3, 1/3, 1/3], [1/3, 1/6, 2/3], [1/3, 2/3, 1/6]]])
        for n, h in self.histories.items():
            np.testing.assert_allclose(h['stress_Pa'][n], 60e6*a, rtol=0, atol=.009)
            np.testing.assert_allclose(h['stress_Pa'][-1], -60e6*a, rtol=0, atol=.009)
            np.testing.assert_allclose(h['plastic_strain'][n], .003, rtol=0, atol=1e-12)
            np.testing.assert_allclose(h['plastic_strain'][-1], .005, rtol=0, atol=1e-12)
            np.testing.assert_array_equal(h['prescribed_F'][-1], np.broadcast_to(np.eye(3), (4, 3, 3)))

    def test_entire_piecewise_history_and_raw_work(self):
        for n, h in self.histories.items():
            r = assess(self.cfg, n, h)
            self.assertEqual(r['status'], 'PASS', r)
            np.testing.assert_allclose(h['equivalent_strain_rate_per_s'], .004, rtol=0, atol=1e-14)

    def test_independent_energy_oracles(self):
        for n, h in self.histories.items():
            np.testing.assert_allclose(h['cumulative_plastic_work_J'][-1], .45, rtol=0, atol=1e-10)
            np.testing.assert_allclose(h['elastic_shear_diagnostic_J'][-1], .045, rtol=0, atol=1e-10)
            ref = .495 if n % 4 == 0 else .495-.315/n**2
            np.testing.assert_allclose(h['cumulative_stress_work_J'][-1], ref, rtol=0, atol=1e-10)

    def test_unaligned_work_converges_without_adding_return_heat(self):
        errors = [float(np.max(np.abs(self.histories[n]['cumulative_stress_work_J'][-1]-.495))) for n in (17, 33, 65)]
        for old, new in zip(errors, errors[1:]):
            self.assertLess(new/old, .4)
        self.assertLess(errors[-1]/.495, 1e-3)
        returns = [self.histories[n]['cumulative_algorithmic_return_J'][-1] for n in (16, 32, 64)]
        self.assertTrue(np.all(returns[0] > returns[1]))
        self.assertTrue(np.all(returns[1] > returns[2]))
        self.assertTrue(np.all(returns[2] > 0))

    def test_benchmark_law_refuses_unprescribed_inputs(self):
        law = PerfectPlasticBenchmark(90e6, 293.15)
        for ep, rate, temp in [(-1., 1., 293.15), (0., -1., 293.15),
                               (0., 1., 300.), (np.nan, 1., 293.15)]:
            with self.assertRaises(ValueError):
                law.yield_stress(ep, rate, temp)
        self.assertIs(self.cfg['production_authorized'], False)
        self.assertEqual(float(law.yield_stress(0., 0., 293.15)), 90e6)

    def test_fixed_oracle_rejects_changed_problem_or_grid(self):
        with self.assertRaises(ValueError):
            run_history(self.cfg, 18)
        with self.assertRaises(ValueError):
            assess(self.cfg, 18, self.histories[16])
        for key, changed in [('shear', 29e9), ('peak_log_strain', .003), ('segment_duration', 2.)]:
            cfg = copy.deepcopy(self.cfg)
            cfg['parameters'][key]['value'] = changed
            with self.assertRaises(ValueError):
                run_history(cfg, 16)


if __name__ == '__main__':
    unittest.main(verbosity=2)

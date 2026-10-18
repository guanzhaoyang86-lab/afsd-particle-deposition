"""Read-only 2026-10-01 accepted-data snapshot, NOT general material validation.

The expected particle identities/times describe this archived run only. The
thermal reference is independent SI arithmetic on solved particle temperatures
and recorded sources; no mechanics/contact/thermal update function is called.

Before inspecting the thermal residual, the review declared the absolute budget
16 * eps64 * (accepted_steps + 1) * sum(m*cp) * max(abs(T0), abs(T)).
This conservatively allows accumulated binary64 temperature addition/reduction
roundoff at the large absolute-temperature scale. It is a snapshot accounting
check, not a rigorous error bound for the entire spatial thermal discretization
and does not modify any solver tolerance.
"""
import hashlib
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.output import restore
from afsd_scaffold.vtk_legacy import read_point_cloud

RUN = ROOT/'outputs/continuation_e3_20261001/to_2us'
ORIGIN = ROOT/'outputs/short_window_e3_20261001/dx05_cfl01_05us'


class ContinuationPhysicalEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads((RUN/'manifest.json').read_text())
        cls.config = json.loads((RUN/'resolved_config.json').read_text())
        cls.pa, cls.state = restore(RUN/'checkpoint.npz')
        cls.diagnostics = [json.loads(line) for line in
                           (RUN/'diagnostics.jsonl').read_text().splitlines() if line.strip()]
        debug = cls.config['parameters']['debug']
        cls.cp = float(debug['cp']['value'])
        cls.t0 = float(debug['temperature']['value'])
        cls.chi = float(debug['plastic_heat_fraction']['value'])
        cls.ledger = cls.state['ledger']
        cls.plastic_work = math.fsum(float(v) for v in cls.pa.plastic_work)
        cls.temperature_energy = math.fsum(
            float(m)*cls.cp*(float(t)-cls.t0)
            for m, t in zip(cls.pa.m, cls.pa.temperature))
        cls.source_energy = math.fsum((cls.chi*cls.plastic_work,
                                      cls.ledger['friction_heat_J'],
                                      cls.ledger['contact_heat_J']))
        thermal_capacity = math.fsum(float(m)*cls.cp for m in cls.pa.m)
        absolute_temperature = max(abs(cls.t0), float(np.abs(cls.pa.temperature).max()))
        cls.thermal_budget = (16*np.finfo(np.float64).eps*(cls.state['count']+1)
                              *thermal_capacity*absolute_temperature)
        print(json.dumps(dict(
            scope='ACCEPTED_DATA_SNAPSHOT_NOT_GENERAL_MATERIAL_VALIDATION',
            temperature_energy_J=cls.temperature_energy,
            independent_source_energy_J=cls.source_energy,
            thermal_residual_J=cls.temperature_energy-cls.source_energy,
            declared_roundoff_budget_J=cls.thermal_budget,
            plastic_work_J=cls.plastic_work,
            numerical_return_dissipation_J=cls.ledger['numerical_return_dissipation_J'],
            final_total_energy_residual_J=cls.diagnostics[-1]['energy_residual_J']),
            indent=2), flush=True)

    def test_accepted_checkpoint_identity_and_exact_snapshot_endpoint(self):
        self.assertEqual(self.manifest['status'], 'COMPLETE')
        self.assertEqual(self.state['count'], 191)
        self.assertEqual(self.state['t'], 2e-6)
        self.assertEqual(self.manifest['checkpoint']['step'], 191)
        self.assertEqual(hashlib.sha256((RUN/'checkpoint.npz').read_bytes()).hexdigest(),
                         self.manifest['checkpoint']['sha256'])
        self.assertEqual(hashlib.sha256((RUN/'resolved_config.json').read_bytes()).hexdigest(),
                         self.manifest['configuration']['sha256'])
        self.assertEqual(self.diagnostics[-1]['step'], 191)
        self.assertEqual(self.diagnostics[-1]['time_s'], 2e-6)
        self.assertEqual(self.diagnostics[-1]['alarms'], [])
        self.assertFalse(list(RUN.glob('REJECTED*.npz')))
        self.assertFalse((RUN/'HARD_TIMEOUT.json').exists())

    def test_three_plastic_ids_are_free_initial_bar_below_shoulder(self):
        pa = self.pa
        positive = pa.ep > 1e-10  # observation floor, not a constitutive limit
        self.assertEqual(set(int(v) for v in pa.particle_id[positive]), {541, 542, 543})
        self.assertTrue(np.all(pa.ep >= 0))
        layer = self.config['parameters']['process']['layer_height']['value']
        reference_z = {541: .00025, 542: .00075, 543: .00125}
        for i in np.flatnonzero(positive):
            self.assertEqual(int(pa.component_id[i]), 1)
            self.assertEqual(int(pa.actuator[i]), 0)
            self.assertAlmostEqual(float(pa.x0[i]), 0., delta=1e-18)
            self.assertAlmostEqual(float(pa.y0[i]), 0., delta=1e-18)
            self.assertAlmostEqual(float(pa.z0[i]), reference_z[int(pa.particle_id[i])], delta=1e-18)
            self.assertLess(float(pa.z0[i]), layer)
            self.assertGreater(float(pa.plastic_work[i]), 0.)
        self.assertEqual(float(pa.m[pa.feed_crossed == 1].sum()), 0.)

    def test_uniform_initial_temperature_and_si_capacity_are_actual_inputs(self):
        debug = self.config['parameters']['debug']
        self.assertEqual(debug['cp']['unit'], 'J/(kg*K)')
        self.assertEqual(debug['temperature']['unit'], 'K')
        self.assertEqual(debug['plastic_heat_fraction']['unit'], '1')
        self.assertGreater(self.cp, 0.)
        origin_series = json.loads((ORIGIN/'series.json').read_text())
        frame = origin_series['frames'][0]
        self.assertEqual((frame['step'], frame['time_s']), (0, 0.))
        path = ORIGIN/frame['file']
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), frame['sha256'])
        data = read_point_cloud(path)
        np.testing.assert_array_equal(data['scalars']['temperature_K'],
                                      np.full(len(self.pa.m), self.t0))

    def test_particle_plastic_work_matches_separate_ledger_without_return_heat(self):
        self.assertEqual(self.chi, 1.)  # specific captured configuration, not a general assumption
        self.assertTrue(np.all(self.pa.plastic_work >= 0))
        budget = (64*np.finfo(np.float64).eps*(self.state['count']+1)
                  *max(abs(self.plastic_work), abs(self.ledger['plastic_heat_J'])))
        self.assertGreater(self.plastic_work, 0.)
        self.assertAlmostEqual(self.plastic_work, self.ledger['plastic_heat_J'], delta=budget)
        self.assertGreater(self.ledger['numerical_return_dissipation_J'], 0.)

    def test_independent_global_thermal_source_identity(self):
        # Internal pair conduction redistributes heat; global net is zero.
        # Contact conduction may be negative and is added with its raw sign.
        self.assertLess(self.ledger['contact_heat_J'], 0.)
        self.assertGreater(self.ledger['friction_heat_J'], 0.)
        self.assertAlmostEqual(self.temperature_energy, self.source_energy,
                               delta=self.thermal_budget)

    def test_snapshot_excludes_return_dissipation_or_energy_residual_as_heat(self):
        # Counterfactual source additions must be resolvable above roundoff.
        # A PASS rejects those additions for this evidence; it does not certify
        # every future branch of the solver or its constitutive free energy.
        for name, extra in (
            ('numerical return dissipation', self.ledger['numerical_return_dissipation_J']),
            ('total energy residual', self.diagnostics[-1]['energy_residual_J'])):
            with self.subTest(forbidden_source=name):
                self.assertGreater(abs(extra), 10*self.thermal_budget)
                self.assertGreater(abs(self.temperature_energy-(self.source_energy+extra)),
                                   self.thermal_budget)


if __name__ == '__main__':
    unittest.main(verbosity=2)

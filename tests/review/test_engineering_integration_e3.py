"""E3 evidence checks preserve failed-run classification, never certify AFSD."""
import sys
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from inspect_engineering_integration import accepted_endpoint, inspect
from afsd_backend.pysph.output import restore


class IntegrationEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = ROOT/'outputs/engineering_e3_20261001/dx05_cfl01'
        cls.report = inspect(cls.directory)

    def test_accepted_not_rejected_state_is_reported(self):
        r = self.report
        self.assertEqual(r['accepted_step'], 21)
        self.assertEqual(r['accepted_frames'], 22)
        self.assertAlmostEqual(r['accepted_time_s'], 2.2029619823669545e-7, delta=1e-20)
        self.assertEqual(r['rejected']['candidate_step'], 22)
        self.assertEqual(r['rejected']['particle_id'], 476)

    def test_partial_run_never_claims_complete_or_deposition(self):
        r = self.report
        self.assertEqual(r['runner_status'], 'FAILED')
        self.assertFalse(r['reached_requested_end'])
        for field in ('paper_reproduction', 'deposition_validated', 'mixing_validated', 'production_authorized'):
            self.assertIs(r[field], False)
        self.assertEqual(r['observed']['plastified_particle_count'], 0)

    def test_rejected_checkpoint_cannot_be_restarted(self):
        with self.assertRaises(ValueError):
            restore(self.directory/'REJECTED_candidate.npz')
        pa, state = restore(self.directory/'checkpoint.npz')
        self.assertEqual(state['count'], 21)
        self.assertEqual(len(np.unique(pa.particle_id)), 657)

    def test_tool_motion_has_independent_units_and_accepted_time(self):
        r = self.report; t = r['accepted_time_s']
        self.assertAlmostEqual(r['observed']['tool_rotation_rad'], 10*np.pi*t, delta=1e-18)
        self.assertAlmostEqual(r['observed']['tool_travel_m'], .00212*t, delta=1e-22)
        self.assertAlmostEqual(r['observed']['prescribed_cap_feed_m'], .00106*t, delta=1e-22)

    def test_invalid_frame_order_rejected(self):
        with self.assertRaises(ValueError):
            accepted_endpoint({'frames': [{'step':1,'time_s':.1}, {'step':1,'time_s':.2}]})
        with self.assertRaises(ValueError):
            accepted_endpoint({'frames': []})


if __name__ == '__main__':
    unittest.main(verbosity=2)

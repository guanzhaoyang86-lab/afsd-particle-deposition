"""Pure acceptance/event regression; saved data replay is not an AFSD rerun."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.interval_acceptance import (
    initial_interval_energy_state, assess_interval_energy, require_interval_energy,
    validate_interval_energy_state, initial_bore_crossing_state,
    observe_square_bore_crossings, validate_bore_crossing_state)


def energy(raw, scale):
    return dict(global_raw_residual_J=raw, absolute_work_scale_J=scale)


def tool(t, *, speed=0., omega=0.):
    c, s = np.cos(omega*t), np.sin(omega*t)
    return dict(origin=np.array([speed*t, 0., 0.]),
                rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


class IntervalEnergyTests(unittest.TestCase):
    def test_absolute_allowance_is_not_reset_each_step_and_rejection_is_pure(self):
        old = initial_interval_energy_state()
        saved = deepcopy(old)
        first = require_interval_energy(old, energy(8e-8, 1e-7))
        self.assertEqual(old, saved)
        second = assess_interval_energy(first, energy(8e-8, 1e-7))
        self.assertTrue(second['gates']['step'])
        self.assertFalse(second['gates']['interval_signed'])
        self.assertFalse(second['gates']['interval_no_cancellation'])
        self.assertIsNone(second['next_state'])
        self.assertEqual(first['steps'], 1)
        self.assertAlmostEqual(second['candidate_interval']['interval_limit_J'], 1.1e-7)

    def test_signed_residual_cancellation_cannot_qualify_history(self):
        first = require_interval_energy(initial_interval_energy_state(), energy(8e-8, 1e-7))
        second = assess_interval_energy(first, energy(-8e-8, 1e-7))
        self.assertTrue(second['gates']['step'])
        self.assertTrue(second['gates']['interval_signed'])
        self.assertFalse(second['gates']['interval_no_cancellation'])
        self.assertEqual(second['candidate_interval']['cumulative_raw_residual_J'], 0.)
        with self.assertRaisesRegex(ValueError, 'no_cancellation'):
            require_interval_energy(first, energy(-8e-8, 1e-7))

    def test_nonfinite_and_negative_scale_rejected(self):
        original = initial_interval_energy_state()
        for raw, scale in [(np.nan, 1.), (np.inf, 1.), (-np.inf, 1.),
                           (0., np.nan), (0., np.inf), (0., -1.)]:
            with self.subTest(raw=raw, scale=scale), self.assertRaises(ValueError):
                assess_interval_energy(original, energy(raw, scale))
        self.assertEqual(original, initial_interval_energy_state())

    def test_full_checkpoint_validation_rejects_tampered_contract_or_summary(self):
        accepted = require_interval_energy(initial_interval_energy_state(), energy(2e-8, 1e-6))
        self.assertTrue(validate_interval_energy_state(json.loads(json.dumps(accepted))))
        for key, value in [('absolute_allowance_J', 2e-7), ('steps', 0),
                           ('cumulative_raw_residual_J', 0.), ('interval_limit_J', np.inf)]:
            bad = deepcopy(accepted); bad[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_interval_energy_state(bad)
        failed = assess_interval_energy(accepted, energy(2e-7, 1e-8))['candidate_interval']
        with self.assertRaisesRegex(ValueError, 'failed interval energy prefix'):
            validate_interval_energy_state(failed)

    def test_saved_old_refinement_replay_rejects_second_prefix(self):
        path = ROOT/'outputs/hot_quasistatic_probe_20261004/refinement'
        original_records = []
        for attempt in [1, 2, 3, 7, 8, 9]:
            data = json.loads((path/f'attempt_{attempt:02d}/energy.json').read_text(encoding='utf-8'))
            record = data['mechanical']
            self.assertLessEqual(abs(record['global_raw_residual_J']), record['registered_energy_limit_J'])
            original_records.append(record)
        state = require_interval_energy(initial_interval_energy_state(), original_records[0])
        second = assess_interval_energy(state, original_records[1])
        self.assertFalse(second['accepted'])
        self.assertTrue(second['gates']['step'])
        self.assertIsNone(second['next_state'])
        self.assertAlmostEqual(second['candidate_interval']['cumulative_raw_residual_J'],
                               1.710267547846782e-7, places=20)
        # Historical accepted files are read-only evidence; do not rewrite them
        # or continue from a prefix that this new contract rejects.


class BoreObservationTests(unittest.TestCase):
    def make(self, x, components=None, masses=None, pose_at_time=tool):
        x = np.asarray(x, float)
        ids = np.arange(10, 10+len(x))
        components = np.ones(len(x), int) if components is None else np.asarray(components)
        masses = np.ones(len(x)) if masses is None else np.asarray(masses, float)
        state = initial_bore_crossing_state(x, ids, components, masses, 0., pose_at_time(0.), 2.)
        return x, ids, components, masses, state

    def test_counts_only_initially_upper_bar_inside_open_square(self):
        x, ids, comp, mass, state = self.make([[0, 0, 1], [3, 0, 1], [0, 0, -1],
                                              [0, 0, 1], [1, 0, 1]],
                                             components=[1, 1, 1, 0, 1], masses=[2, 3, 4, 5, 6])
        end = x.copy(); end[:, 2] = -2
        before = deepcopy(state); xcopy = x.copy(); endcopy = end.copy()
        observed = observe_square_bore_crossings(state, x, end, ids, comp, mass, 0., 1., tool)
        self.assertEqual(observed['crossing_ids'], [10])
        self.assertEqual(observed['crossing_mass_kg'], 2.)
        self.assertAlmostEqual(observed['new_events'][0]['time_s'], 1/3)
        self.assertEqual(state, before)
        np.testing.assert_array_equal(x, xcopy); np.testing.assert_array_equal(end, endcopy)

    def test_uses_actual_rotated_translated_pose_at_intersection(self):
        pose_at_time = lambda t: tool(t, speed=2., omega=np.pi/2)
        # At t=.5, relative global xy=(1.1,0) rotates to (.778,-.778), inside.
        # At either endpoint the same point is outside or evaluated wrongly.
        x, ids, comp, mass, state = self.make([[2.1, 0, 1]], pose_at_time=pose_at_time)
        end = x.copy(); end[:, 2] = -1
        result = observe_square_bore_crossings(state, x, end, ids, comp, mass,
                                               0., 1., pose_at_time)
        self.assertEqual(result['crossing_ids'], [10])
        event = result['new_events'][0]
        np.testing.assert_allclose(event['tool_local_position_m'],
                                   [1.1/np.sqrt(2), -1.1/np.sqrt(2), 0.], atol=1e-15)

    def test_crossing_geometry_uses_intersection_not_final_particle_location(self):
        x, ids, comp, mass, state = self.make([[3, 0, 1]])
        end = np.array([[0., 0., -1.]])
        result = observe_square_bore_crossings(state, x, end, ids, comp, mass, 0., 1., tool)
        self.assertEqual(result['crossing_ids'], [])  # plane xy=(1.5,0), outside

    def test_returning_particle_is_not_counted_twice_and_initial_lower_stays_ineligible(self):
        x, ids, comp, mass, state = self.make([[0, 0, 1], [0, 0, -1]])
        down = np.array([[0., 0., -1.], [0., 0., -1.]])
        first = observe_square_bore_crossings(state, x, down, ids, comp, mass, 0., 1., tool)
        upper = np.array([[0., 0., 1.], [0., 0., 1.]])
        up = observe_square_bore_crossings(first['next_state'], down, upper, ids, comp,
                                          mass, 1., 2., tool)
        last = observe_square_bore_crossings(up['next_state'], upper, down, ids, comp,
                                            mass, 2., 3., tool)
        self.assertEqual(last['crossing_ids'], [10])
        self.assertEqual(last['new_events'], [])
        self.assertEqual(last['crossing_mass_kg'], 1.)

    def test_sorting_is_identity_safe_and_passive_tag_not_an_input(self):
        x, ids, comp, mass, state = self.make([[0, 0, 1], [3, 0, 1]], masses=[2, 3])
        end = x.copy(); end[:, 2] = -1
        order = np.array([1, 0])
        result = observe_square_bore_crossings(state, x[order], end[order], ids[order],
                                               comp[order], mass[order], 0., 1., tool)
        self.assertEqual(result['crossing_ids'], [10])
        self.assertEqual(result['crossing_mass_kg'], 2.)
        self.assertTrue(validate_bore_crossing_state(json.loads(json.dumps(result['next_state']))))

    def test_bad_geometry_time_mass_and_checkpoint_evidence_reject(self):
        x, ids, comp, mass, state = self.make([[0, 0, 1]])
        end = x.copy(); end[:, 2] = -1
        with self.assertRaisesRegex(ValueError, 'accepted time'):
            observe_square_bore_crossings(state, x, end, ids, comp, mass, .1, 1., tool)
        with self.assertRaisesRegex(ValueError, 'mass changed'):
            observe_square_bore_crossings(state, x, end, ids, comp, mass*2, 0., 1., tool)
        badend = end.copy(); badend[0, 0] = np.nan
        with self.assertRaises(ValueError):
            observe_square_bore_crossings(state, x, badend, ids, comp, mass, 0., 1., tool)
        nextstate = observe_square_bore_crossings(state, x, end, ids, comp, mass,
                                                  0., 1., tool)['next_state']
        nextstate['events'][0]['tool_local_position_m'][0] = 5.
        with self.assertRaisesRegex(ValueError, 'through-square'):
            validate_bore_crossing_state(nextstate)


if __name__ == '__main__':
    unittest.main()

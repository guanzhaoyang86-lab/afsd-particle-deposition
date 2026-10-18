"""Manufactured administrative scheduling tests; no physics or processes."""
from pathlib import Path
import copy
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import afsd_cpu_local_interval_scheduler_20261015 as schedule


class PendingIntervals(unittest.TestCase):
    def test_start_clips_the_nominal_interval_at_phase_end(self):
        state = schedule.start(.2199, .001)
        self.assertEqual(state['pending_targets'], (.22,))
        self.assertEqual(schedule.current_dt(state, .2199), .22 - .2199)

    def test_start_at_target_is_complete(self):
        state = schedule.start(.22, .001)
        self.assertEqual(state['pending_targets'], ())
        for operation, value in ((schedule.current_dt, .22), (schedule.split, .22),
                                 (schedule.accept, .22)):
            with self.subTest(operation=operation.__name__), self.assertRaises(ValueError):
                operation(state, value)

    def test_more_than_four_splits_have_no_counter_limit(self):
        state = schedule.start(.08, .00011865234375)
        original = state['pending_targets'][0]
        for _ in range(12):
            previous = copy.deepcopy(state)
            state = schedule.split(state, .08)
            self.assertEqual(state['pending_targets'][-1], original)
            self.assertGreater(schedule.current_dt(state, .08), 0.)
            self.assertEqual(state['pending_targets'][1:], previous['pending_targets'])
        self.assertEqual(len(state['pending_targets']), 13)

    def test_split_and_accept_do_not_mutate_input(self):
        first = schedule.start(.08, .001)
        first_copy = copy.deepcopy(first)
        second = schedule.split(first, .08)
        self.assertEqual(first, first_copy)
        second_copy = copy.deepcopy(second)
        third = schedule.accept(second, second['pending_targets'][0])
        self.assertEqual(second, second_copy)
        self.assertEqual(len(third['pending_targets']), 1)

    def test_success_restores_remaining_original_target_then_nominal_dt(self):
        original = schedule.start(.08, .00011865234375)
        smaller = schedule.split(original, .08)
        middle, end = smaller['pending_targets']
        remaining = schedule.accept(smaller, middle)
        self.assertEqual(remaining['pending_targets'], (end,))
        self.assertEqual(schedule.current_dt(remaining, middle), end - middle)
        complete = schedule.accept(remaining, end)
        self.assertEqual(complete['pending_targets'], ())
        next_macro = schedule.start(complete['accepted_time_s'], complete['nominal_dt'])
        self.assertEqual(next_macro['nominal_dt'], original['nominal_dt'])
        self.assertGreater(schedule.current_dt(next_macro, end), end - middle)

    def test_one_accept_consumes_only_one_of_many_boundaries(self):
        state = schedule.start(.08, .001)
        for _ in range(5):
            state = schedule.split(state, .08)
        pending = state['pending_targets']
        advanced = schedule.accept(state, pending[0])
        self.assertEqual(advanced['pending_targets'], pending[1:])

    def test_json_list_state_is_accepted_without_aliasing(self):
        state = schedule.start(.08, .001)
        state['pending_targets'] = list(state['pending_targets'])
        smaller = schedule.split(state, .08)
        self.assertIsInstance(smaller['pending_targets'], tuple)
        self.assertEqual(len(state['pending_targets']), 1)

    def test_stale_parent_time_is_rejected(self):
        state = schedule.start(.08, .001)
        for operation in (schedule.current_dt, schedule.split):
            with self.subTest(operation=operation.__name__), self.assertRaises(ValueError):
                operation(state, math.nextafter(.08, math.inf))

    def test_small_endpoint_roundoff_retains_actual_saved_time(self):
        state = schedule.start(.08, .001)
        actual = math.nextafter(state['pending_targets'][0], math.inf)
        accepted = schedule.accept(state, actual)
        self.assertEqual(accepted['accepted_time_s'], actual)
        self.assertEqual(accepted['pending_targets'], ())

    def test_unfinished_endpoint_is_not_accepted(self):
        state = schedule.start(.08, .001)
        with self.assertRaises(ValueError):
            schedule.accept(state, .0809)

    def test_accept_cannot_skip_a_pending_boundary(self):
        state = schedule.split(schedule.start(.08, .001), .08)
        with self.assertRaises(ValueError):
            schedule.accept(state, state['pending_targets'][1])

    def test_accept_cannot_go_beyond_phase_end(self):
        state = schedule.start(.2199, .001)
        with self.assertRaises(ValueError):
            schedule.accept(state, math.nextafter(.22, math.inf))

    def test_nonrepresentable_fine_midpoint_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'NOT_REPRESENTABLE'):
            schedule.start(1., math.ulp(1.), target=2.)

    def test_half_underflow_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'NOT_REPRESENTABLE'):
            schedule.start(0., math.ulp(0.), target=1.)

    def test_split_stops_only_when_fine_children_are_unrepresentable(self):
        state = schedule.start(1., 2. * math.ulp(1.), target=2.)
        self.assertEqual(schedule.current_dt(state, 1.), 2. * math.ulp(1.))
        with self.assertRaisesRegex(ValueError, 'NOT_REPRESENTABLE'):
            schedule.split(state, 1.)

    def test_nan_infinity_bool_and_negative_inputs_are_rejected(self):
        for invalid in (float('nan'), float('inf'), -1., True, '0.08'):
            for name, arguments in (
                ('time', (invalid, .001, .22)),
                ('dt', (.08, invalid, .22)),
                ('target', (.08, .001, invalid)),
            ):
                with self.subTest(field=name, value=invalid), self.assertRaises(ValueError):
                    schedule.start(*arguments)
            state = schedule.start(.08, .001)
            with self.subTest(field='end', value=invalid), self.assertRaises(ValueError):
                schedule.accept(state, invalid)

    def test_zero_dt_and_start_beyond_phase_are_rejected(self):
        with self.assertRaises(ValueError):
            schedule.start(.08, 0.)
        with self.assertRaises(ValueError):
            schedule.start(.23, .001)

    def test_corrupt_state_order_or_nonfinite_boundary_is_rejected(self):
        base = schedule.start(.08, .001)
        for pending in ((.081, .0805), (.081, .081), (float('nan'),), (.23,), (.08,), (True,)):
            state = dict(base, pending_targets=pending)
            with self.subTest(pending=pending), self.assertRaises(ValueError):
                schedule.current_dt(state, .08)

    def test_roundoff_cannot_swallow_an_adjacent_machine_boundary(self):
        state = schedule.start(1., 4. * math.ulp(1.), target=2.)
        state = schedule.split(state, 1.)
        with self.assertRaises(ValueError):
            schedule.accept(state, state['pending_targets'][1])
        with self.assertRaises(ValueError):
            schedule.accept(state, math.nextafter(state['pending_targets'][0], math.inf))


if __name__ == '__main__':
    unittest.main()

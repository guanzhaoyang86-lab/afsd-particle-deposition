"""Proof-binding regression fixtures; no physical solve or actual receipt."""
import copy
import unittest
from unittest.mock import patch
import readback_afsd_deposition_20261007 as reader


class PairBindingTests(unittest.TestCase):
    def setUp(self):
        self.screen = dict(passed=True, global_screen={'pass': True}, local_screen={'passed': True},
            crossing_ids_match=True, macro_index=1, fine_checkpoint_sha256=['fine1', 'fine2'], time_s=.001)
        self.rows = [dict(half_index=i, macro_index=1, sha256='fine'+str(i), time_s=i*.0005,
            screen_path='fixture/pair.json', screen_sha256='screen', observer_path='fixture/observer.json',
            observer_sha256='observer') for i in (1, 2)]

    def check(self, rows=None, screen=None):
        screen = self.screen if screen is None else screen
        with patch.object(reader.io, 'read', return_value=screen), patch.object(reader.io.core, 'sha',
                side_effect=lambda p: 'screen' if p.name == 'pair.json' else 'observer'):
            return reader.checked_screened_rows(dict(screened_history=self.rows if rows is None else rows,
                                                     pair_screens=[screen]))

    def test_exact_two_fine_states_are_bound_to_their_screen(self):
        self.assertEqual(self.check(), self.rows)

    def test_coarse_or_unrelated_candidate_cannot_borrow_passing_screen(self):
        rows = copy.deepcopy(self.rows)
        rows[0]['sha256'] = 'coarse'
        with self.assertRaises(ValueError):
            self.check(rows=rows)

    def test_partial_or_repeated_pair_is_rejected(self):
        with self.assertRaises(ValueError):
            self.check(rows=self.rows[:1])
        with self.assertRaises(ValueError):
            self.check(rows=[self.rows[1], self.rows[1]])

    def test_overall_pass_cannot_hide_failed_physical_comparison(self):
        screen = copy.deepcopy(self.screen)
        screen['local_screen']['passed'] = False
        with self.assertRaises(ValueError):
            self.check(screen=screen)


if __name__ == '__main__':
    unittest.main()

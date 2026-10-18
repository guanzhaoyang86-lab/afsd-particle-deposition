"""Pure geometry/cost-reference regression, no solver evaluation."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from assess_nonlinear_continuation_cost import geometric_gap

class GeometryCostTests(unittest.TestCase):
    def test_bar_straddling_plane_uses_positive_upstream_gap_not_negative_global_min(self):
        z=np.array([.25,.75,1.25,1.75,2.25,2.75])*1e-3
        r=geometric_gap(z,.0015,.00005)
        self.assertAlmostEqual(r['minimum_upstream_original_gap_m'],.00025)
        self.assertEqual(r['initially_below_plane_count'],3)
        self.assertEqual(r['initially_on_or_above_plane_count'],3)
        self.assertFalse(r['any_center_crossing_permitted_by_current_guard'])
        self.assertAlmostEqual(r['minimum_absolute_gap_under_guard_m'],.0002)
        self.assertAlmostEqual(r['minimum_upstream_original_gap_m']/.00212,25/212)

    def test_guard_reaching_either_side_is_not_claimed_to_prevent_crossing(self):
        self.assertTrue(geometric_gap([.0014,.0018],.0015,.0002)['any_center_crossing_permitted_by_current_guard'])

    def test_no_upstream_cohort_does_not_produce_negative_feed_time(self):
        with self.assertRaises(ValueError):geometric_gap([.00025,.00075],.0015,.00005)

if __name__=='__main__':unittest.main(verbosity=2)

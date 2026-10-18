"""Saved-data comparator must not compare altered physics or partial endpoints."""
from pathlib import Path
import copy
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import compare_process_friction_cone as comparison


class ComparisonGuards(unittest.TestCase):
    def test_changed_physics_and_event_substitution_cannot_masquerade_as_sensitivity(self):
        case=comparison.load_verified(ROOT/'outputs/process_friction_cone_q6_20261001/oblique_i4')
        changed=copy.deepcopy(case)
        changed['config']['parameters']['mu']['value']=.2
        with self.assertRaisesRegex(ValueError,'changed physical configuration'):
            comparison.equal_physics(case,changed)
        changed=copy.deepcopy(case)
        last=len(changed['report']['records'])-1
        changed['report']['records'][last]['event_kind']='CONE_LIMIT'
        with self.assertRaisesRegex(ValueError,'unique saved requested state'):
            comparison.requested_indices(changed,np.array([.2]))

    def test_partial_readback_cannot_be_promoted_into_same_end_comparison(self):
        # Isolate the comparator's gate; the reader's partial-state tests are separate.
        with patch.object(comparison,'inspect',return_value=dict(status='PASS_PARTIAL_READBACK',accepted_time_s=.199)):
            with self.assertRaisesRegex(ValueError,'complete common 0.2 s endpoint'):
                comparison.load_verified(Path('unread-partial-run'))


if __name__=='__main__':unittest.main(verbosity=2)

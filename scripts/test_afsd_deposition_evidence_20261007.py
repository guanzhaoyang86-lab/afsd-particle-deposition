"""Synthetic observation fixtures; these are not AFSD physics results."""
import unittest
import numpy as np
from afsd_deposition_evidence_20261007 import summarize_frames

P = dict(tool_radius=2.7, dx=.5, bore_width=2.504)


def frame(t=0.):
    return dict(time_s=t, ids=np.array([1, 2, 3, 4]),
        position=np.array([[2., 0., .5], [-3.1, 0., .5], [0., 3.1, .5], [-3.1, 0., .5]]),
        mass=np.ones(4), ep=np.zeros(4), eligible_ids=[1, 2, 3], crossing_ids=[1, 2, 3],
        tool_origin=[0., 0., 1.5], tool_rotation=np.eye(3))


class EvidenceTests(unittest.TestCase):
    def test_initial_lower_material_and_lateral_extrusion_are_excluded(self):
        r = summarize_frames([frame()], P)['final']
        self.assertEqual(r['supplied_under_shoulder']['particle_ids'], [1])
        self.assertEqual(r['supplied_trailing']['particle_ids'], [2])

    def test_persistence_tracks_ids_through_permutation(self):
        a, b = frame(), frame(.1)
        order = [2, 0, 3, 1]
        b['ep'][1] = .25
        b['position'][1, 0] -= .02
        for key in ('ids', 'position', 'mass', 'ep'):
            b[key] = b[key][order]
        r = summarize_frames([a, b], P)['final']['trailing_residence'][0]
        self.assertEqual(r['particle_id'], 2)
        self.assertEqual(r['consecutive_saved_frames'], 2)
        self.assertAlmostEqual(r['sampled_residence_s'], .1)
        self.assertAlmostEqual(r['displacement_since_first_trailing_sample_m'], .02)
        self.assertAlmostEqual(r['ep_increment_since_first_saved_crossing'], .25)

    def test_leaving_trailing_region_resets_residence(self):
        a, b, c = frame(), frame(.1), frame(.2)
        b['position'][1, 0] = 0.
        r = summarize_frames([a, b, c], P)['final']['trailing_residence'][0]
        self.assertEqual(r['consecutive_saved_frames'], 1)
        self.assertEqual(r['sampled_residence_s'], 0.)

    def test_changing_per_particle_mass_is_rejected(self):
        a, b = frame(), frame(.1)
        b['mass'][0] += .2
        b['mass'][1] -= .2
        with self.assertRaises(ValueError):
            summarize_frames([a, b], P)

    def test_relabelled_eligibility_is_rejected(self):
        a, b = frame(), frame(.1)
        b['eligible_ids'].append(4)
        with self.assertRaises(ValueError):
            summarize_frames([a, b], P)

    def test_lost_crossing_history_is_rejected(self):
        a, b = frame(), frame(.1)
        b['crossing_ids'].remove(2)
        with self.assertRaises(ValueError):
            summarize_frames([a, b], P)

    def test_duplicate_time_is_rejected(self):
        with self.assertRaises(ValueError):
            summarize_frames([frame(), frame()], P)

    def test_no_saved_frame_does_not_claim_deposition(self):
        with self.assertRaises(ValueError):
            summarize_frames([], P)
        r = summarize_frames([frame()], P)
        self.assertFalse(r['complete_deposition_proven'])
        self.assertFalse(r['metallurgical_bonding_proven'])


if __name__ == '__main__':
    unittest.main()

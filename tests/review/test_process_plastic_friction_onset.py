"""Independent onset references; old RED retained, no actual AFSD path run."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'tests/review'))
from test_process_friction_cone_q6 import ElasticVectorOwner
from probe_process_friction_cone import CandidateFailure, advance_cone
from probe_process_plastic_friction_onset import advance_onset

E = np.array([.6, .8])
LIMITS = dict(force_tolerance_N=9e-5, slip_tolerance_m=1e-11,
              work_tolerance_J=1e-6, negative_dissipation_tolerance_J=1e-10,
              maximum_active_sets=16, maximum_events_per_interval=12,
              event_bisections=40, event_time_tolerance_s=1e-12,
              mode_probe_fraction=1e-4)


class GuardedElasticOwner(ElasticVectorOwner):
    """Actual cone-root fixture plus independent production-style commit gates."""
    def __init__(self, K, count, wall, base, **kwargs):
        super().__init__(K, count, wall, base, **kwargs)
        self.limits = dict(LIMITS)

    def commit(self, c):
        if c['base_step'] != self.count:
            raise ValueError('stale fixture candidate')
        capacity = self.mu * c['normal_force']
        tau = c['tangent_force']
        slip = c['relative_increment']
        ft = self.limits['force_tolerance_N']
        if np.any(np.linalg.norm(tau, axis=1) > capacity + ft):
            raise CandidateFailure('fixture cone exceeded', c)
        if np.max(np.linalg.norm(c['response']['residual'], axis=1)) > ft:
            raise CandidateFailure('fixture root exceeded', c)
        moving = np.linalg.norm(slip, axis=1) > self.limits['slip_tolerance_m']
        off = np.abs(np.linalg.norm(self.tangent_force, axis=1)
                     - self.mu * self.normal_force) > ft
        opposite = np.sum(self.tangent_force * tau, axis=1) < 0
        if np.any(moving & (off | opposite)):
            raise CandidateFailure('unsplit cone stick/slip transition', c)
        raw_q = -np.sum(.5 * (self.tangent_force + tau) * slip, axis=1)
        if np.any(raw_q < -self.limits['negative_dissipation_tolerance_J']):
            raise CandidateFailure('negative raw fixture friction', c)
        super().commit(c)


def single(initial=1., direction=E, wall_rate=.2 * E):
    direction = np.asarray(direction, dtype=float)
    wall_rate = np.asarray(wall_rate, dtype=float)
    return GuardedElasticOwner(np.eye(2), 1,
        lambda t: (t * wall_rate)[None, :],
        lambda t: -initial * direction)


def coupled():
    return GuardedElasticOwner(np.kron([[2., -1.], [-1., 2.]], np.eye(2)), 2,
        lambda t: np.array([.4, .2])[:, None] * t * E,
        lambda t: -np.r_[14/15 * E, 13/15 * E])


class OldConeOnsetRed(unittest.TestCase):
    @unittest.expectedFailure
    def test_legacy_boundary_onset_should_advance(self):
        """RED: old tiny mode probe misses exact boundary onset; no lower exists."""
        b = single()
        advance_cone(b, 1.)
        np.testing.assert_allclose(b.u[:2], [0., 0.], atol=1e-10, rtol=0)
        self.assertAlmostEqual(b.heat, .2, places=10)


class IndependentBoundaryOnset(unittest.TestCase):
    def close(self, actual, expected, atol=1e-10):
        np.testing.assert_allclose(actual, expected, atol=atol, rtol=0)

    def test_01_boundary_outward_starts_sliding_without_false_interior_event(self):
        b = single(); calls = []
        advance_onset(b, 1., lambda: calls.append(b.time))
        self.close(b.u[:2], [0., 0.]); self.close(b.tangent_force, [E])
        self.close([b.heat, b.wall_work, b.material_work], [.2, .2, 0.])
        self.assertEqual(calls, [1.]); self.assertEqual(len(b.rows), 1)
        self.assertTrue(b.rows[0]['released'][0])

    def test_02_boundary_inward_eligibility_does_not_force_sliding(self):
        b = single(wall_rate=-.2 * E); advance_onset(b, 1.)
        self.close(b.u[:2], -.2 * E); self.close(b.tangent_force, [.8 * E])
        self.close([b.heat, b.wall_work, b.material_work], [0., -.18, -.18])
        self.close(b.rows[-1]['q'], [0.]); self.assertEqual(b.time, 1.)

    def test_03_genuinely_interior_contact_keeps_exact_event(self):
        b = single(initial=.5, wall_rate=E); advance_onset(b, 1.)
        self.close([r['t'] for r in b.rows], [.5, 1.])
        self.assertEqual(b.rows[0]['kind'], 'CONE_LIMIT')
        self.assertFalse(b.rows[0]['released'][0])
        self.close(b.u[:2], .5 * E); self.close(b.tangent_force, [E])
        self.close([b.heat, b.wall_work, b.material_work], [.5, .875, .375])

    def test_04_coupled_neighbor_activates_at_two_thirds_not_one_half(self):
        b = coupled(); advance_onset(b, 1.)
        self.close([r['t'] for r in b.rows], [2/3, 1.])
        self.close(b.u[:4], np.r_[E/15, 2*E/15])
        self.close(b.tangent_force, [E, E])
        self.close(np.sum([r['q'] for r in b.rows], axis=0), [1/3, 1/15])
        self.close([b.heat, b.wall_work, b.material_work], [2/5, 44/75, 14/75])
        self.assertTrue(b.rows[0]['released'][0])
        self.assertFalse(b.rows[0]['released'][1])
        self.assertTrue(b.rows[-1]['released'].all())

    def test_05_signed_initial_offset_is_not_projected_and_raw_ledger_preserved(self):
        for delta, q, ww, wm in [(4e-5, .2000440008, .200004, -.0000400008),
                                (-4e-5, .1999560008, .199996, .0000399992)]:
            with self.subTest(delta=delta):
                b = single(initial=1+delta)
                self.close(b.tangent_force, [(1+delta)*E])
                before = b.tangent_force.copy(); advance_onset(b, 1.)
                self.close(b.u[:2], -delta*E); self.close(b.tangent_force, [E])
                self.close([b.heat, b.wall_work, b.material_work], [q, ww, wm])
                self.close(b.wall_work-b.material_work-b.heat, 0.)
                self.assertGreater(abs(float(np.linalg.norm(before)-1)), 1e-5)

    def test_06_transverse_discrete_endpoint_rotates_tangent_force(self):
        b = single(direction=[1., 0.], wall_rate=[.2, 1.6])
        advance_onset(b, 1.)
        self.close(b.u[:2], [-.4, .8]); self.close(b.tangent_force, [[.6, .8]])
        self.close([b.heat, b.wall_work, b.material_work], [.8, .8, 0.])
        self.close(b.u[:2]-np.array([.2, 1.6]), [-.6, -.8])

    def test_07_straight_onset_and_coupled_work_partition_invariance(self):
        for factory, expected in [(single, [.2, .2, 0.]),
                                  (coupled, [2/5, 44/75, 14/75])]:
            b = factory()
            for t in np.linspace(0., 1., 11)[1:]:
                advance_onset(b, float(t))
            self.close([b.heat, b.wall_work, b.material_work], expected)

    def test_08_reversal_root_rejected_atomically_then_driver_splits_exact_event(self):
        b = single(wall_rate=-3*E)
        before = b.u.copy(), b.tangent_force.copy(), b.time, b.count
        c = b.branch(1., np.array([True]))
        c.update(released=np.array([True]), event_kind='REQUESTED')
        with self.assertRaisesRegex(CandidateFailure, 'unsplit|opposite|reversal'):
            b.commit(c)
        np.testing.assert_array_equal(b.u, before[0])
        np.testing.assert_array_equal(b.tangent_force, before[1])
        self.assertEqual((b.time, b.count), before[2:])
        advance_onset(b, 1.)
        self.close([r['t'] for r in b.rows], [2/3, 1.])
        self.close(b.u[:2], -2*E); self.close(b.tangent_force, [-E])
        self.close([b.heat, b.wall_work, b.material_work], [1., 1., 0.])
        self.assertFalse(b.rows[0]['released'][0])

    def test_09_invalid_initial_overcone_state_not_silently_projected(self):
        b = single(initial=1.001); old = b.tangent_force.copy()
        with self.assertRaises((CandidateFailure, ValueError)):
            advance_onset(b, 1.)
        np.testing.assert_array_equal(b.tangent_force, old)
        self.assertEqual(b.count, 0); self.assertEqual(b.time, 0.)

    def test_10_shared_free_vector_dofs_allow_simultaneous_sliding_and_sticking(self):
        K = np.kron([[2., 0., -1.], [0., 2., -1.], [-1., -1., 3.]], np.eye(2))
        b = GuardedElasticOwner(K, 2,
            lambda t: np.array([.3, -.3])[:, None] * t * E,
            lambda t: -np.r_[.75*E, .75*E, .5*E])
        self.close(b.tangent_force, [E, E])
        advance_onset(b, 1.)
        self.close(b.u[:6], np.r_[-.06*E, -.3*E, -.12*E])
        self.close(b.tangent_force, [E, .52*E])
        self.close(b.u[b.td]-b.wall(1.), [-.36*E, np.zeros(2)])
        self.close(b.rows[-1]['q'], [.36, 0.])
        self.close([b.heat, b.wall_work, b.material_work], [.36, .072, -.288])
        force = b.K @ (b.u-b.fullbase(1.))
        self.close(force[4:6], [0., 0.])


if __name__ == '__main__':
    unittest.main(verbosity=2)

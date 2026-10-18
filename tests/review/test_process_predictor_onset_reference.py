"""Independent SPD spring + original onset driver; NOT an AFSD simulation.

Literal references and tolerances frozen before execution in
docs/PROCESS_PREDICTOR_ONSET_REFERENCES.md. No reference uses contact solver.
"""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'tests/review'))
from test_process_feasible_predictor import LinearFixture
from test_process_plastic_friction_onset import LIMITS
from probe_process_feasible_predictor import solve_with_feasible_predictor
from probe_process_friction_cone import (ConeBlock, feasible_cone,
                                         vector_contact_work)
from probe_process_plastic import CandidateFailure
from probe_process_plastic_friction_onset import advance_onset


class PredictorEventFixture(LinearFixture):
    """Original coupled state/root plus independent atomic fixture ledger.

    q is FREE, released tangent s is FREE. Only wall motion is prescribed.
    Commit is a fixture contract, not a replacement for production commit.
    """
    def __init__(self, use_predictor=True):
        super().__init__()
        self.mu = 10.
        self.limits = dict(LIMITS)
        self.use_predictor = use_predictor
        self.wall_work = self.material_work = self.heat = 0.
        self.rows = []
        self.branch_calls = []

    def wall(self, t):
        return np.array([-1., 2.5*t, 0.])

    def branch(self, t, released, maximum_iterations=None):
        maximum = 14 if maximum_iterations is None else maximum_iterations
        if self.use_predictor:
            c = solve_with_feasible_predictor(self, t, released, maximum)
        else:
            c = ConeBlock.branch(self, t, released, maximum)
        self.branch_calls.append(dict(time_s=t,
            released=np.asarray(released).copy(),
            predictor_corrections=c.get('predictor_corrections', 0),
            iterations=c['iterations']))
        return c

    def energy(self):
        u = self.displacement.ravel()
        return float(.5*u@self.K@u)

    def commit(self, candidate):
        # No state changes until every gate has passed.
        if candidate['base_step'] != self.count:
            raise ValueError('stale fixture candidate')
        if not np.isfinite(candidate['time_s']) or candidate['time_s'] <= self.time:
            raise ValueError('fixture time must advance')
        c = self.state(self.trial(candidate['displacement'],
                                 candidate['time_s']), candidate['released'])
        lim = self.limits
        if (not feasible_cone(c, np.ones(1, dtype=bool), self.mu, lim)
                or np.max(np.linalg.norm(c['law_residual_N'], axis=1))
                > lim['force_tolerance_N']):
            raise CandidateFailure('fixture cone contact law/root rejected', c)
        if (c['equilibrium_relative'] > self.ccfg['numerical']['force_relative']
                or np.any(np.abs(c['gap']) > lim['slip_tolerance_m'])):
            raise CandidateFailure('fixture equilibrium/gap rejected', c)
        moving = c['modes'] != 0
        off = (np.abs(np.linalg.norm(self.tangent_force, axis=1)
                      - self.mu*self.normal_force) > lim['force_tolerance_N'])
        opposite = np.sum(self.tangent_force*c['tangent_force'], axis=1) < 0
        if np.any(moving & (off | opposite)):
            raise CandidateFailure('unsplit fixture stick/slip transition', c)
        oldforce = np.zeros_like(self.X)
        oldforce[self.nodes, 0] = -self.normal_force
        oldforce[self.nodes, 1:] = self.tangent_force
        du = c['displacement'] - self.displacement
        dw = self.wall(c['time_s']) - self.wall(self.time)
        work = vector_contact_work(oldforce, c['contact_force'], du, dw,
                                   self.nodes, c['dt_s'])
        q = work['friction_increment_J']
        ww, wm = work['wall_work_J'], work['material_work_J']
        if not np.isfinite(np.r_[q, ww, wm]).all():
            raise CandidateFailure('nonfinite fixture work', c)
        if np.any(q < -lim['negative_dissipation_tolerance_J']):
            raise CandidateFailure('negative raw fixture friction', c)
        new_u = c['displacement'].ravel()
        delta_es = float(.5*new_u@self.K@new_u-self.energy())
        if max(abs(ww-wm-q.sum()), abs(wm-delta_es)) > lim['work_tolerance_J']:
            raise CandidateFailure('fixture work identity rejected', c)
        kind = candidate.get('event_kind', 'DIRECT_COMMIT')
        if kind not in ('DIRECT_COMMIT', 'REQUESTED', 'CONE_LIMIT'):
            raise ValueError('invalid fixture event kind')
        bracket = candidate.get('event_bracket_s')
        if kind == 'CONE_LIMIT':
            a = np.asarray(bracket, dtype=float)
            if (a.shape != (2,) or not np.isfinite(a).all()
                    or a[0] != c['time_s']
                    or not 0 <= a[1]-a[0] <= lim['event_time_tolerance_s']):
                raise ValueError('invalid fixture event bracket')
        self.rows.append(dict(time_s=c['time_s'], kind=kind,
            bracket=bracket, released=c['released'].copy(),
            displacement=c['displacement'].copy(), normal=c['normal_force'].copy(),
            tau=c['tangent_force'].copy(), q=q.copy(),
            slip=c['relative_increment'].copy(), ww=ww, wm=wm))
        self.displacement = c['displacement'].copy()
        self.normal_force = c['normal_force'].copy()
        self.tangent_force = c['tangent_force'].copy()
        self.time = c['time_s']
        self.count += 1
        self.wall_work += ww
        self.material_work += wm
        self.heat += float(q.sum())


def accepted_signature(b):
    return (b.displacement.tobytes(), b.normal_force.tobytes(),
            b.tangent_force.tobytes(), b.time, b.count,
            b.wall_work, b.material_work, b.heat, len(b.rows))


class LegacyPredictorOnsetRed(unittest.TestCase):
    @unittest.expectedFailure
    def test_original_onset_branch_stops_at_infeasible_target_guess(self):
        b = PredictorEventFixture(use_predictor=False)
        advance_onset(b, 1.)
        np.testing.assert_allclose(b.displacement[0, 1], 31/15,
                                   atol=1e-10, rtol=0)


class IndependentPredictorOnset(unittest.TestCase):
    def close(self, actual, expected, atol=1e-10):
        np.testing.assert_allclose(actual, expected, atol=atol, rtol=0)

    def test_01_initial_guess_recovered_but_overcone_target_not_committed(self):
        b = PredictorEventFixture(); before = accepted_signature(b)
        c = b.branch(1., np.array([False]))
        self.close(c['displacement'], [[-1., 5/2, 0.], [3/4, 0., 0.]])
        self.close(c['normal_force'], [1/4])
        self.close(c['tangent_force'], [[23/4, 0.]])
        self.assertEqual(c['predictor_corrections'], 1)
        self.assertFalse(feasible_cone(c, np.array([False]), b.mu, b.limits))
        with self.assertRaisesRegex(CandidateFailure, 'cone contact law/root'):
            b.commit(c)
        self.assertEqual(accepted_signature(b), before)

    def test_02_original_onset_locates_analytic_event_and_keeps_free_dofs(self):
        b = PredictorEventFixture(); accepted=[]; trace=[]
        advance_onset(b, 1., lambda: accepted.append(b.time), trace.append)
        self.close(accepted, [62/75, 1.], atol=1e-12)
        self.assertEqual([r['kind'] for r in b.rows], ['CONE_LIMIT', 'REQUESTED'])
        self.close(b.rows[0]['displacement'], [[-1., 31/15, 0.], [8/15, 0., 0.]])
        self.close(b.displacement, [[-1., 31/15, 0.], [8/15, 0., 0.]])
        self.close(b.normal_force, [7/15]); self.close(b.tangent_force, [[14/3, 0.]])
        self.close(b.rows[-1]['slip'], [[-13/30, 0.]])
        self.assertFalse(b.rows[0]['released'][0])
        self.assertTrue(b.rows[-1]['released'][0])
        lo, hi = b.rows[0]['bracket']
        self.assertLessEqual(lo, 62/75)
        self.assertGreaterEqual(hi, 62/75)
        self.assertLessEqual(hi-lo, 1e-12)
        self.assertEqual(sum(r['stage']=='ACCEPTED' for r in trace), 2)
        self.assertGreaterEqual(sum(r['stage']=='INTERIOR_BISECTION' for r in trace), 1)
        self.assertTrue(any(c['predictor_corrections'] == 1 for c in b.branch_calls))
        self.assertTrue(all(c['iterations'] <= 14 for c in b.branch_calls))

    def test_03_independent_integrated_work_and_raw_signed_heat(self):
        b = PredictorEventFixture(); self.close(b.energy(), 3/4)
        advance_onset(b, 1.)
        self.close([b.rows[0]['q'][0], b.rows[-1]['q'][0]], [0., 91/45])
        self.close([b.heat, b.wall_work, b.material_work, b.energy()],
                   [91/45, 1139/180, 155/36, 91/18])
        self.close(b.wall_work-b.material_work-b.heat, 0.)
        self.close(b.energy()-3/4-b.material_work, 0.)

    def test_04_direct_fully_released_root_cannot_skip_interior_event(self):
        b = PredictorEventFixture(); before = accepted_signature(b)
        c = b.branch(1., np.array([True]))
        self.close(c['displacement'], [[-1., 31/15, 0.], [8/15, 0., 0.]])
        self.assertTrue(feasible_cone(c, np.array([True]), b.mu, b.limits))
        with self.assertRaisesRegex(CandidateFailure, 'unsplit fixture'):
            b.commit(c)
        self.assertEqual(accepted_signature(b), before)

    def test_05_independent_same_path_partition_invariance(self):
        b = PredictorEventFixture()
        for t in (.2, .5, .9, 1.):
            advance_onset(b, t)
        events = [r for r in b.rows if r['kind']=='CONE_LIMIT']
        self.assertEqual(len(events), 1)
        self.close(events[0]['time_s'], 62/75, atol=1e-12)
        self.close(b.displacement, [[-1., 31/15, 0.], [8/15, 0., 0.]])
        self.close([b.heat, b.wall_work, b.material_work],
                   [91/45, 1139/180, 155/36])

    def test_06_tensile_stick_equilibrium_remains_rejected_without_mutation(self):
        b = PredictorEventFixture(); before = accepted_signature(b)
        with self.assertRaisesRegex(CandidateFailure, 'negative normal at interior equilibrium') as caught:
            b.branch(1.6, np.array([False]))
        self.close(caught.exception.candidate['internal_force'][0, 0], -.5)
        self.close(caught.exception.candidate['displacement'][1, 0], 1.5)
        self.assertEqual(accepted_signature(b), before)

    def test_07_actual_work_utility_retains_negative_component(self):
        force = np.array([[0., 1., 0.], [0., -1., 0.]])
        du = np.array([[0., 1., 0.], [0., .5, 0.]])
        work = vector_contact_work(force, force, du, np.zeros(3),
                                   np.array([0, 1]), .25)
        self.close(work['friction_increment_J'], [-1., .5])
        self.close([work['wall_work_J'], work['material_work_J']], [0., .5])

    def test_08_event_metadata_and_stale_candidate_rejections_are_atomic(self):
        b = PredictorEventFixture(); before = accepted_signature(b)
        c = b.branch(.1, np.array([False]))
        c.update(event_kind='CONE_LIMIT', event_bracket_s=[.1, .11])
        with self.assertRaisesRegex(ValueError, 'event bracket'):
            b.commit(c)
        self.assertEqual(accepted_signature(b), before)
        c['event_kind'] = 'REQUESTED'; c['base_step'] = 1
        with self.assertRaisesRegex(ValueError, 'stale fixture'):
            b.commit(c)
        self.assertEqual(accepted_signature(b), before)


if __name__ == '__main__':
    unittest.main(verbosity=2)

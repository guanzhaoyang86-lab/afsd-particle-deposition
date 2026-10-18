"""Bounded numerical response refinement; no physical process acceptance."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import event_response as original
from afsd_backend.pysph import event_response_refined as refined
import test_event_response as fixtures


def manufactured(curvatures=(3.,), cusp=False):
    """Literal polynomial or absolute-value map, independently evaluated values.

    This fixture represents a completed original response, not a physical kick.
    Values vbar=.5F+aF^2, v1=F+4F^2 give exact central B=.5,C=1 at zero.
    """
    def evaluate(relative):
        h = relative
        zero = np.zeros(3)
        base = dict(mean_velocity_m_s=zero.copy(), velocity_m_s=zero.copy(), axis=None,
                    step_N=None, status='PROBE_COMPLETED', raw_heat_certificate_pass=True)
        r = dict(available=False, accepted=False, B=None, C=None, baseline=base,
            probes=[base], column_evidence=[], force_step_N=np.full(3, h),
            kick_calls=1, kick_seconds=0., response_seconds=0.,
            scope='MANUFACTURED_POLYNOMIAL_NOT_PHYSICAL_KICK')
        matrices = np.zeros((2, 3, 3))
        for axis, curvature in enumerate(curvatures):
            levels = []
            for factor in (1., .5):
                pair = []
                for sign in (1., -1.):
                    f = sign*factor*h
                    values = np.array([abs(f), abs(f)]) if cusp else np.array([.5*f+curvature*f*f, f+4*f*f])
                    vectors = np.zeros((2, 3)); vectors[:, axis] = values
                    probe = dict(mean_velocity_m_s=vectors[0], velocity_m_s=vectors[1],
                        axis=axis, step_N=f, status='PROBE_COMPLETED',
                        raw_heat_certificate_pass=True, physical_branch_change=None)
                    r['probes'].append(probe); r['kick_calls'] += 1
                    pair.append(vectors)
                forward = pair[0]/(factor*h); backward = -pair[1]/(factor*h)
                levels.append((.5*(forward+backward), forward, backward))
            coarse = levels[0][0]; fine, forward, backward = levels[1]
            central_bound = .002*np.maximum(abs(coarse), abs(fine))+2e-5
            bound = .002*np.maximum(abs(forward), abs(backward))+2e-5
            gap = abs(forward-backward)
            column = dict(axis=axis, refinement_pass=bool(np.all(abs(fine-coarse) <= central_bound)),
                forward_half_h_s_kg=forward, backward_half_h_s_kg=backward,
                one_sided_allowed_s_kg=bound, one_sided_difference_s_kg=gap,
                one_sided_pass=bool(np.all(gap <= bound)))
            r['column_evidence'].append(column)
            if not column['one_sided_pass']:
                r['reason'] = refined.RETRY_REASON
                return r
            matrices[:, :, axis] = fine
        r.update(available=True, B=matrices[0], C=matrices[1], reason='MANUFACTURED_ORIGINAL_BOUNDS_PASSED')
        return r
    return evaluate


class RefinedResponse(unittest.TestCase):
    def record(self, name, value):
        print('REFINED_EVENT_EVIDENCE '+json.dumps(dict(name=name, result=original._plain(value)),
            sort_keys=True, allow_nan=False), flush=True)
        self.assertFalse(value['accepted'])
        self.assertEqual(value['kick_calls'], sum(a['kick_calls'] for a in value['attempts']))
        return value

    def test_literal_smooth_curvature_two_halvings_original_bounds(self):
        r = self.record('literal_polynomial', refined._bounded(manufactured(), 1e-3))
        self.assertTrue(r['available']); self.assertEqual(r['refinement_halvings'], 2)
        self.assertAlmostEqual(r['B'][0, 0], .5, delta=1e-12)
        self.assertAlmostEqual(r['C'][0, 0], 1., delta=1e-12)
        self.assertEqual([a['force_step_N'][0] for a in r['attempts']], [1e-3, 5e-4, 2.5e-4])
        self.assertFalse(r['attempts'][0]['available'])

    def test_abs_cusp_rejects_without_retry(self):
        r = self.record('literal_abs_cusp', refined._bounded(manufactured(cusp=True), 1e-3))
        self.assertFalse(r['available']); self.assertEqual(r['refinement_halvings'], 0)
        self.assertEqual(r['reason'], 'REFINEMENT_SMOOTHNESS_CERTIFICATE_FAILED')
        self.assertIsNone(r['B']); self.assertIsNone(r['C'])

    def test_unrelated_original_failure_never_retries(self):
        calls = []
        def failure(relative):
            calls.append(relative)
            return dict(available=False, accepted=False, B=None, C=None, kick_calls=2,
                kick_seconds=.1, reason='PHYSICAL_ZERO_SLIP_SEQUENCE_CHANGED')
        r = self.record('unrelated_failure', refined._bounded(failure, 1e-3))
        self.assertFalse(r['available']); self.assertEqual(calls, [1e-3])
        self.assertEqual(r['reason'], 'PHYSICAL_ZERO_SLIP_SEQUENCE_CHANGED')

    def test_first_failing_axes_change_compared_by_identity(self):
        r = self.record('axis_switch', refined._bounded(manufactured((1.5, 3.)), 1e-3))
        self.assertTrue(r['available']); self.assertEqual(r['refinement_halvings'], 2)
        self.assertEqual([a['column_evidence'][-1]['axis'] for a in r['attempts'][:2]], [0, 1])
        compared = [c['axis'] for c in r['refinement_certificates'] if c['kind'].startswith('CROSS')]
        self.assertEqual(compared, [0, 1])

    def test_budget_exhaustion_is_unavailable(self):
        r = self.record('fixed_budget_exhausted', refined._bounded(manufactured((100.,)), 1e-3))
        self.assertFalse(r['available']); self.assertEqual(len(r['attempts']), 3)
        self.assertEqual(r['reason'], 'REFINEMENT_BUDGET_EXHAUSTED')
        self.assertIsNone(r['B']); self.assertIsNone(r['C'])

    def test_literal_reversal_unchanged_no_retry(self):
        r = self.record('literal_reversal', refined.response(**fixtures.reversal()))
        self.assertTrue(r['available']); self.assertEqual(r['refinement_halvings'], 0)
        self.assertAlmostEqual(r['B'][0, 0], 19/54, delta=2e-7)
        self.assertAlmostEqual(r['C'][0, 0], 7/9, delta=2e-7)
        self.assertAlmostEqual(r['baseline']['friction_raw_J'], 7/18, delta=2e-11)

    def test_real_pid1_old_unavailable_refined_available_and_pure(self):
        source = ROOT/'src/afsd_backend/pysph/event_response.py'
        path = ROOT/'outputs/process_global_20261003/root_attempt/latest_trial.json.gz'
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        with gzip.open(path, 'rt', encoding='utf-8') as stream:
            trace = json.load(stream)['evaluation']['boundary']['traces'][0]
        self.assertEqual(trace['particle_id'], 1)
        before = original._digest(trace)
        with (ROOT/'outputs/process_event_20261003/response_probe/resolved_config.json').open(encoding='utf-8') as stream:
            mu = json.load(stream)['parameters']['friction']
        r = self.record('real_pid1', refined.response_from_trace(trace, mu=mu))
        self.assertFalse(r['attempts'][0]['available'])
        self.assertEqual(r['attempts'][0]['reason'], refined.RETRY_REASON)
        self.assertTrue(r['available'], r['reason'])
        self.assertLessEqual(r['refinement_halvings'], 2)
        self.assertTrue(all(a['saved_replay_comparison']['passed'] for a in r['attempts']))
        self.assertTrue(all(p['raw_heat_certificate_pass'] for a in r['attempts'] for p in a['probes']))
        self.assertEqual(before, original._digest(trace))
        self.assertEqual(source_hash, hashlib.sha256(source.read_bytes()).hexdigest())


if __name__ == '__main__':
    unittest.main(verbosity=2)

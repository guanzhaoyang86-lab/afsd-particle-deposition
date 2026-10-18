"""Q7R bounded root recovery: archived failure is preserved, never accepted.

Only pure branch candidates are evaluated here. The checkpoint reconstruction
is the existing diagnostic-only helper, not a production restart. Independent
contact checks below use Euclidean norms and the literal Coulomb coefficient;
they do not ask the tested feasibility function to supply reference answers.
"""
from pathlib import Path
import copy
import json
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from diagnose_process_plastic_friction_q7 import restore_diagnostic, state_signature
from inspect_process_plastic_friction import inspect
from probe_process_plastic import CandidateFailure
from probe_process_plastic_friction import PlasticFrictionBlock
from probe_process_plastic_friction_recovery import RecoveryBlock

RUN = ROOT / 'outputs/process_plastic_friction_q7_20261001/n3_i4'


class BoundedPlasticFrictionRecovery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # This reader verifies archived source/config/artifact digests, accepted
        # journals, and the rejected/accepted separation before raw loading.
        cls.readback = inspect(RUN)
        if cls.readback['status'] != 'PASS_PARTIAL_READBACK':
            raise AssertionError('requires the original verified BLOCKED run')
        cls.report = json.loads((RUN / 'report.json').read_text(encoding='utf-8'))
        cls.rejected = json.loads((RUN / 'REJECTED_manifest.json').read_text(encoding='utf-8'))
        with np.load(RUN / 'history.npz', allow_pickle=False) as data:
            cls.history = {key: data[key].copy() for key in data.files}
        with np.load(RUN / 'REJECTED_candidate.npz', allow_pickle=False) as data:
            cls.rejected_arrays = {key: data[key].copy() for key in data.files}
        cls.released = cls.rejected_arrays['released'].copy()

    def restored(self):
        body = restore_diagnostic(RUN, self.report, self.history)
        self.assertEqual((body.count, body.time), (6, .05))
        return body

    def test_01_only_declared_numerical_budget_changes(self):
        baseline, recovery = PlasticFrictionBlock(), RecoveryBlock()
        before, after = state_signature(baseline), state_signature(recovery)
        self.assertEqual(before['metadata']['qcfg']['numerical']['maximum_iterations'], 12)
        self.assertEqual(after['metadata']['qcfg']['numerical']['maximum_iterations'], 14)
        after['metadata']['qcfg']['numerical']['maximum_iterations'] = 12
        variant = after['metadata']['pfcfg'].pop('recovery_variant')
        self.assertEqual(variant['variant'], 'Q7R_BOUNDED_NEWTON14')
        self.assertIs(variant['production_authorized'], False)
        self.assertEqual(after, before)
        # No overriding the event search, force/material law, or commit guards.
        for method in ('advance', 'commit', 'state', 'trial', 'tangent', 'wall'):
            self.assertIs(getattr(RecoveryBlock, method), getattr(PlasticFrictionBlock, method))

    def test_02_original_twelve_correction_failure_replays_exactly(self):
        body = self.restored()
        before = state_signature(body)
        with self.assertRaisesRegex(CandidateFailure, 'cone branch equilibrium/root did not converge') as caught:
            body.branch(.0625, self.released, maximum_iterations=12)
        failure = caught.exception
        self.assertEqual(failure.iterations, self.rejected['iterations'])
        for key, reference in self.rejected_arrays.items():
            np.testing.assert_array_equal(failure.candidate[key], reference, err_msg=key)
        self.assertGreater(failure.iterations[-1]['root_relative'], 1e-8)
        self.assertEqual(state_signature(body), before)

    def test_03_fourteen_budget_root_is_not_contact_feasibility(self):
        body = self.restored()
        before = state_signature(body)
        candidate = body.branch(.0625, self.released, maximum_iterations=14)
        self.assertEqual(candidate['iterations'], 13)
        # Assemble the actually free/released rows independently, not using
        # assess_candidate/feasible_cone or the branch's stored convergence flag.
        unknown = body.free.copy()
        unknown[body.contact_dofs] = False
        unknown[body.tangent_dofs[~self.released].ravel()] = False
        residual = candidate['residual'].ravel().copy()
        residual[body.tangent_dofs[self.released].ravel()] = candidate['law_residual_N'][self.released].ravel()
        root = np.linalg.norm(residual[unknown]) / body.force_scale
        self.assertLess(root, 1e-8)
        # Literal diagnostic reference, not regenerated with tested functions.
        self.assertAlmostEqual(root, 1.637588335220718e-10, delta=1e-12)
        row = int(np.flatnonzero(body.pa.particle_id[body.nodes] == 21)[0])
        self.assertFalse(self.released[row])
        excess = np.linalg.norm(candidate['tangent_force'][row]) - .1 * candidate['normal_force'][row]
        self.assertGreater(excess, 5.)
        self.assertAlmostEqual(excess, 5.691232925791425, delta=1e-8)
        self.assertEqual(state_signature(body), before)

    def test_04_real_commit_rejects_infeasible_root_without_partial_write(self):
        body = self.restored()
        candidate = body.branch(.0625, self.released, maximum_iterations=14)
        # An isolated test-only binding exercises the real guard. The diagnostic
        # helper and production restore policy remain unchanged on disk.
        body.commit = types.MethodType(PlasticFrictionBlock.commit, body)
        before = state_signature(body)
        with self.assertRaisesRegex(CandidateFailure, 'cone contact law/root rejected'):
            body.commit(candidate)
        self.assertEqual(state_signature(body), before)
        self.assertEqual((body.count, body.time), (6, .05))
        with self.assertRaisesRegex(RuntimeError, 'diagnostic cannot commit/advance'):
            body.advance(.0625)
        # A new diagnostic load still has both of its original safety guards.
        fresh = self.restored()
        with self.assertRaisesRegex(RuntimeError, 'diagnostic cannot commit/advance'):
            fresh.commit(candidate)

    def test_05_budget_override_cannot_expand_or_silently_coerce(self):
        body = RecoveryBlock()
        before = state_signature(body)
        for invalid in (-1, 15, 14., True, np.int64(14), float('nan'), '14'):
            with self.subTest(value=repr(invalid)):
                with self.assertRaisesRegex(ValueError, 'integer correction budget'):
                    body.branch(.0625, self.released, maximum_iterations=invalid)
        # Interface-only spy: proves the wrapper forwards the declared bound;
        # this is not counted as a mechanics solution or an accepted interval.
        sentinel = object()
        with patch.object(PlasticFrictionBlock, 'branch', return_value=sentinel) as inherited:
            self.assertIs(body.branch(.0625, self.released), sentinel)
            inherited.assert_called_once_with(.0625, self.released, maximum_iterations=14)
        self.assertEqual(state_signature(body), before)

    def test_06_short_actual_preload_propagates_budget_through_event_search(self):
        # Actual .025 s preload/event regression only, not the .2 s coupled path.
        body = RecoveryBlock()
        ids = body.pa.particle_id.copy()
        labels = body.pa.tracer_type.copy()
        original = PlasticFrictionBlock.branch
        calls, accepted = [], []
        previous_ep = body.pa.ep.copy()

        def recording_branch(owner, t, released, maximum_iterations=None):
            calls.append(dict(time=t, maximum_iterations=maximum_iterations,
                              base_time=owner.time, base_step=owner.count))
            return original(owner, t, released, maximum_iterations=maximum_iterations)

        def inspect_accepted():
            nonlocal previous_ep
            np.testing.assert_array_equal(body.pa.particle_id, ids)
            np.testing.assert_array_equal(body.pa.tracer_type, labels)
            self.assertTrue(np.all(body.pa.ep >= previous_ep))
            previous_ep = body.pa.ep.copy()
            self.assertTrue(np.all(body.normal_force >= 0.))
            excess = np.linalg.norm(body.tangent_force, axis=1) - .1 * body.normal_force
            self.assertLessEqual(float(excess.max()), body.limits['force_tolerance_N'])
            self.assertLessEqual(float(np.linalg.norm(body.last_law_residual, axis=1).max()),
                                 body.limits['force_tolerance_N'])
            self.assertGreaterEqual(float(body.last_friction.min()),
                                    -body.limits['negative_dissipation_tolerance_J'])
            self.assertTrue(np.all(body.last_plastic >= 0.))
            np.testing.assert_allclose(body.displacement[body.fixed], 0., atol=0., rtol=0.)
            np.testing.assert_allclose(body.displacement[body.nodes, 0], body.wall(body.time)[0],
                                       atol=body.limits['slip_tolerance_m'], rtol=0.)
            accepted.append(body.record())

        with patch.object(PlasticFrictionBlock, 'branch', new=recording_branch):
            body.advance(.025, inspect_accepted)
        self.assertEqual(body.time, .025)
        self.assertEqual(len(accepted), body.count)
        self.assertTrue(calls)
        self.assertTrue(all(call['maximum_iterations'] == 14 for call in calls))
        self.assertGreater(len(calls), len(accepted))  # probes/bisections, not just accepted endpoints
        self.assertGreater(float(body.pa.plastic_work.sum()), 0.)
        self.assertGreater(float(body.pa.ep.max()), 0.)
        self.assertGreater(float(body.friction.sum()), 0.)
        self.assertTrue(body.events)
        self.assertTrue(any(row['event_kind'] == 'CONE_LIMIT' for row in accepted))
        for event in body.events:
            self.assertEqual(event['kind'], 'CONE_LIMIT')
            lower, upper = event['bracket_s']
            self.assertEqual(lower, event['physical_time_s'])
            self.assertGreaterEqual(upper - lower, 0.)
            self.assertLessEqual(upper - lower, body.limits['event_time_tolerance_s'])
            self.assertLessEqual(upper, .025)


if __name__ == '__main__':
    unittest.main(verbosity=2)

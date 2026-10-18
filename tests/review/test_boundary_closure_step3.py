"""Independent polynomial-moment counterexample; PASS means correct rejection."""
from pathlib import Path
from contextlib import redirect_stdout
import hashlib
import io
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from probe_boundary_perturbation import evaluate, checked_project, IncompatibleMoments, main


class InteriorPerturbation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report, cls.control, cls.moved = evaluate()

    def test_only_interior_position_changes_not_weights_or_boundary(self):
        a, b = self.control, self.moved
        for key in ('volume_m3', 'mass_kg', 'rho_kg_m3', 'h_m', 'signed_face_area_m2'):
            np.testing.assert_array_equal(a[key], b[key])
        outside = np.any(np.abs(a['position_m']) == 1., axis=1)
        np.testing.assert_array_equal(a['position_m'][outside], b['position_m'][outside])
        delta = b['position_m']-a['position_m']
        center = np.flatnonzero(np.all(a['position_m'] == 0., axis=1))
        expected = np.zeros((125, 3))
        expected[center[0]] = [.025, -.015, .010]
        np.testing.assert_array_equal(delta, expected)
        self.assertAlmostEqual(b['volume_m3'].sum(), 8., delta=1e-12)
        for k in range(3):
            e = b['signed_face_area_m2'][k]
            self.assertAlmostEqual(e[e > 0].sum(), 4., delta=1e-12)
            self.assertAlmostEqual(e[e < 0].sum(), -4., delta=1e-12)

    def test_literal_independent_reference_moments(self):
        # Single center weight is 0.125 m^3; only its first moment changes.
        expected_mu = np.array([.003125, -.001875, .00125])
        np.testing.assert_allclose(self.report['perturbed']['first_volume_moment_m4'],
                                   expected_mu, atol=1e-12, rtol=0.)
        np.testing.assert_allclose(self.report['perturbed']['quadrature_centroid_m'],
                                   [.000390625, -.000234375, .00015625], atol=1e-12, rtol=0.)
        # Independently expanded mu_l delta_mk + delta_lk mu_m for ell=1m.
        reference = np.array([
            [[0, 0, 0, 0], [0, .00625, -.001875, .00125],
             [0, -.001875, 0, 0], [0, .00125, 0, 0]],
            [[0, 0, 0, 0], [0, 0, .003125, 0],
             [0, .003125, -.00375, .00125], [0, 0, .00125, 0]],
            [[0, 0, 0, 0], [0, 0, 0, .003125],
             [0, 0, 0, -.001875], [0, .003125, -.001875, .0025]]])
        np.testing.assert_allclose(self.moved['moment_defect_m2'], reference, atol=1e-12, rtol=0.)
        np.testing.assert_allclose(self.control['moment_defect_m2'], 0., atol=1e-12, rtol=0.)

    def test_no_projection_after_incompatible_moments(self):
        a = self.moved
        with patch('probe_boundary_perturbation.project', side_effect=AssertionError('must not run')) as spy:
            with self.assertRaisesRegex(IncompatibleMoments, 'projection NOT_RUN'):
                checked_project(a['position_m'], a['volume_m3'], a['signed_face_area_m2'],
                                a['original_D_per_m'], a['neighbor_i'], a['neighbor_j'])
            spy.assert_not_called()

    def test_regular_control_remains_acceptable(self):
        row = self.report['control']
        self.assertEqual(row['method_status'], 'FIXTURE_PASS')
        self.assertTrue(row['accepted_operator'])
        for value in row['checks'].values():
            self.assertLessEqual(value, 1e-10)

    def test_original_corrected_gradient_is_not_degenerate(self):
        a = self.moved
        measured = np.einsum('kij,ja->iak', a['original_D_per_m'], a['position_m'])
        np.testing.assert_allclose(measured, np.broadcast_to(np.eye(3), (125, 3, 3)),
                                   atol=1e-10, rtol=0.)
        self.assertLess(self.report['perturbed']['original_max_moment_condition'], 1e6)
        self.assertGreater(self.report['perturbed']['max_moment_defect_m2'], 1e-10)

    def test_rejection_is_not_presented_as_physical_state_or_method_pass(self):
        r = self.report['perturbed']
        self.assertEqual(r['method_status'], 'BLOCKED_BY_MOMENT_COMPATIBILITY')
        self.assertEqual(r['projection_status'], 'NOT_RUN')
        self.assertEqual(r['zero_mode_test'], 'NOT_RUN')
        self.assertFalse(r['accepted_operator'])
        self.assertNotIn('candidate_D_per_m', self.moved)
        self.assertEqual(str(self.moved['operator_status']), 'REJECTED_FOR_OPERATOR_CONSTRUCTION')
        self.assertFalse(self.report['production_authorized'])
        self.assertEqual(self.report['physically_accepted_terminal_state'], 'NONE')

    def test_saved_config_byte_hash_and_rejected_artifact(self):
        # Exercise real write/read, so Windows newline translation is visible.
        with tempfile.TemporaryDirectory(prefix='step3_artifact_',
                dir=ROOT/'logs/boundary_step3_20261001') as directory:
            out = Path(directory)/'probe'
            with patch.object(sys, 'argv', ['probe', '--output', str(out)]), redirect_stdout(io.StringIO()):
                self.assertEqual(main(), 1)  # Method remains BLOCKED.
            report = json.loads((out/'report.json').read_text(encoding='utf-8'))
            raw = (out/'fixture_config.json').read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(), report['configuration_sha256'])
            self.assertNotIn(b'\r\n', raw)
            with np.load(out/'REJECTED_sampling_fixture.npz', allow_pickle=False) as data:
                self.assertNotIn('candidate_D_per_m', data.files)
                self.assertEqual(str(data['operator_status']), 'REJECTED_FOR_OPERATOR_CONSTRUCTION')


if __name__ == '__main__':
    unittest.main(verbosity=2)


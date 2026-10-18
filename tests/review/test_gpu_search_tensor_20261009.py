"""CPU-only validation of an explicit search derivative domain extension.

Finite differences use the ORIGINAL physical material.update. Saved failure
inputs are hash-pinned and read only. No equilibrium or GPU solve is executed.
"""
from pathlib import Path
import ast
import copy
import hashlib
import inspect
import sys
import textwrap
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import afsd_search_tensor_20261009 as search
import afsd_gpu_search_trial_20261009 as model
from afsd_backend.pysph.material import update, FKS
from afsd_backend.pysph.material_response import _deviatoric
from afsd_backend.pysph.quasistatic_process import _geometry


def elastic_case(trace=True):
    s = np.array([[[4e5, 2e4, -1e4], [2e4, 2e5, 3e4], [-1e4, 3e4, 1e5]]])
    if not trace:
        s -= np.trace(s, axis1=1, axis2=2)[:, None, None]*np.eye(3)/3
    gradient = np.array([[[.07, -1.2, .2], [1.1, -.03, .3], [-.2, .1, -.04]]])
    return s, np.array([.002]), gradient, np.array([293.15]), 0.0002373046875, 26.3e9, FKS()


def assert_finite_difference(test, args, rtol=2e-6, atol=.02, expect_elastic=True):
    result = search.material_tangent(*args)
    if expect_elastic:
        test.assertFalse(result['plastic_active'].any(), 'elastic fixture must stay inactive')
    test.assertTrue(np.isfinite(result['tangent_Pa_s']).all())
    h = 1e-7
    finite = np.zeros_like(result['tangent_Pa_s'])
    for a in range(3):
        for b in range(3):
            plus = list(args); minus = list(args)
            plus[2] = args[2].copy(); minus[2] = args[2].copy()
            plus[2][:, a, b] += h; minus[2][:, a, b] -= h
            finite[:, :, 3*a+b] = ((update(*plus)[0]-update(*minus)[0])/(2*h)).reshape(-1, 9)
    np.testing.assert_allclose(result['tangent_Pa_s'], finite, rtol=rtol, atol=atol)
    np.testing.assert_array_equal(result['stress_Pa'], update(*args)[0])
    plus = list(args); minus = list(args)
    plus[3] = args[3]+1e-4; minus[3] = args[3]-1e-4
    # Room-temperature elastic fixtures need a forward difference to stay in
    # the original law's supported temperature domain.
    if np.any(minus[3] < args[-1].room):
        thermal = (update(*plus)[0]-update(*args)[0])/1e-4
    else:
        thermal = (update(*plus)[0]-update(*minus)[0])/2e-4
    np.testing.assert_allclose(result['temperature_tangent_Pa_K'], thermal.reshape(-1, 9),
                               rtol=rtol, atol=atol)
    return result


class SearchDomainTests(unittest.TestCase):
    def test_nonzero_trace_elastic_derivative_matches_original_update(self):
        args = elastic_case()
        with self.assertRaisesRegex(ValueError, 'deviatoric'):
            search.original_tangent.material_tangent(*args)
        assert_finite_difference(self, args)

    def test_hydrostatic_trace_is_retained_without_projection(self):
        args = list(elastic_case()); args[0] = np.eye(3)[None]*2e5
        result = assert_finite_difference(self, args)
        self.assertAlmostEqual(float(np.trace(result['stress_Pa'][0])), 6e5, delta=1e-7)

    def test_zero_rate_composite_tangent_matches_original_update(self):
        args = list(elastic_case()); args[2] = np.zeros_like(args[2])
        assert_finite_difference(self, args)

    def test_old_domain_is_numerically_identical(self):
        args = elastic_case(trace=False)
        old = search.original_tangent.material_tangent(*args)
        new = search.material_tangent(*args)
        for key in old:
            if key != 'classification':
                np.testing.assert_array_equal(new[key], old[key], err_msg=key)
        direction = np.arange(9.).reshape(1, 3, 3)/17
        old = search.original_response.trial_direction(args[0], args[2], direction, args[4], args[5])
        new = search.trial_direction(args[0], args[2], direction, args[4], args[5])
        for key in old:
            if key != 'classification':
                np.testing.assert_array_equal(new[key], old[key], err_msg=key)

    def test_input_and_original_providers_are_not_mutated(self):
        args = elastic_case(); before = copy.deepcopy(args)
        bindings = (search.original_response.trial_direction, search.original_tangent.material_tangent,
                    model.parent.ProcessTrial.tangent, model.parent.ProcessTrial.evaluate)
        search.material_tangent(*args)
        for original, saved in zip(args[:-1], before[:-1]):
            np.testing.assert_array_equal(original, saved)
        self.assertEqual(bindings, (search.original_response.trial_direction,
                         search.original_tangent.material_tangent, model.parent.ProcessTrial.tangent,
                         model.parent.ProcessTrial.evaluate))

    def test_old_plastic_domain_is_numerically_identical(self):
        args = list(elastic_case(trace=False))
        args[0] = np.diag([8e8, -4e8, -4e8])[None]
        old = search.original_tangent.material_tangent(*args)
        new = search.material_tangent(*args)
        self.assertTrue(new['plastic_active'].all())
        for key in old:
            if key != 'classification':
                np.testing.assert_array_equal(new[key], old[key], err_msg=key)

    def test_small_force_matrix_uses_unchanged_assembly_and_is_finite(self):
        args = elastic_case(); tangent = search.material_tangent(*args)['tangent_Pa_s']
        D = np.arange(27., dtype=float).reshape(9, 3)/30
        J = search.original_tangent.assemble_frozen_force_tangent(
            D, np.array([1e-10]), tangent, np.array([2700.]), 2700., 68e9,
            args[4], np.zeros((1, 3)))
        self.assertEqual(J.shape, (3, 3))
        self.assertTrue(np.isfinite(J).all())
        self.assertGreater(float(np.linalg.norm(J)), 0.)

    def test_nonfinite_and_asymmetric_stress_remain_rejected(self):
        for field in ('nan', 'inf', 'asymmetric'):
            with self.subTest(field=field):
                args = list(elastic_case()); args[0] = args[0].copy()
                args[0][0, 0, 1] = dict(nan=np.nan, inf=np.inf, asymmetric=2e4+1.)[field]
                with self.assertRaises(ValueError):
                    search.material_tangent(*args)

    def test_symmetry_threshold_is_original_128_epsilon(self):
        args = elastic_case(); s = args[0].copy()
        tol = 128*np.finfo(float).eps*np.max(abs(s))
        s[0, 0, 1] += 2*tol
        with self.assertRaisesRegex(ValueError, 'symmetric'):
            search.trial_direction(s, args[2], np.zeros_like(args[2]), args[4], args[5])

    def test_active_radial_guard_is_preserved(self):
        args = list(elastic_case()); args[0] = np.diag([8e8, -4e8, -3.9e8])[None]
        with self.assertRaisesRegex(ValueError, 'deviatoric'):
            search.material_tangent(*args)
        self.assertIs(search.radial_direction, search.original_response.radial_direction)


class SavedFailureTests(unittest.TestCase):
    def test_actual_failed_parent_matches_fd_without_projecting_trace(self):
        base = ROOT/'outputs/gpu_acceleration_20261007/campaign_10h_20261008_resume_v1/attempt_0005'
        pins = {
            'fine_1/NOT_RESTARTABLE_candidate.npz': 'fb0de26851b8d0000a36e4a296193d1900262302527f66b4bd6faf3d0cd5c948',
            'fine_2/broyden/last_UNCOMMITTED_candidate.npz': 'c7e09f1197893245415942db1fe77ec9934b82d4226a85a13285e6062dad8b76',
            'fine_2/broyden/last_UNCOMMITTED_candidate_graph.npz': 'e12c8f82037d526a7c7c034d07fe4a07a94c385637bd02861d4f2253275faea6',
        }
        for name, digest in pins.items():
            self.assertEqual(hashlib.sha256((base/name).read_bytes()).hexdigest(), digest)
        with np.load(base/'fine_1/NOT_RESTARTABLE_candidate.npz', allow_pickle=True) as stored:
            d = stored['particles'].item()['solid']['arrays']; state = stored['solver_data'].item()
        with np.load(base/'fine_2/broyden/last_UNCOMMITTED_candidate.npz') as candidate, np.load(
                base/'fine_2/broyden/last_UNCOMMITTED_candidate_graph.npz') as graph:
            H, _ = _geometry(candidate['position_m'], candidate['displacement_m'], d['m'],
                             candidate['density_kg_m3'], graph['i'], graph['j'], graph['g'])
            idx = int(np.flatnonzero(d['particle_id'] == 1216)[0]); dt = .00023730468750000005
            s = d['S'].reshape(-1, 3, 3)[idx:idx+1].copy(); p = state['parameters']
            args = (s, d['ep'][idx:idx+1].copy(), H[idx:idx+1]/dt,
                    (d['temperature'][idx:idx+1]+candidate['temperature_K'][idx:idx+1])/2,
                    dt, p['shear'], FKS(**{k:p[k] for k in FKS.__dataclass_fields__}))
        with self.assertRaisesRegex(ValueError, 'deviatoric'):
            _deviatoric(s, 'stress')
        original_bytes = s.tobytes()
        self.assertAlmostEqual(float(np.trace(s[0])), -3.601599019020796e-10, delta=1e-20)
        assert_finite_difference(self, args, expect_elastic=False)
        self.assertEqual(original_bytes, s.tobytes())


class StaticScopeTests(unittest.TestCase):
    def test_static_adaptations_pass_and_physical_bindings_unchanged(self):
        self.assertTrue(search.static_adaptation_audit()['physical_update_unchanged'])
        result = model.operator_reuse_audit()
        self.assertFalse(result['physical_gates_changed'])
        self.assertIs(model.ProcessTrial.evaluate, model.parent.ProcessTrial.evaluate)
        self.assertIs(model.advance_candidate, model.parent.advance_candidate)
        self.assertTrue(issubclass(model.ProcessTrial, model.parent.ProcessTrial))

    def test_removed_branch_gate_is_detected(self):
        tree = ast.parse(Path(search.__file__).read_text())
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'material_tangent')
        function.body = [n for n in function.body if not (isinstance(n, ast.If) and 'returned plastic branch' in ast.unparse(n))]
        with self.assertRaisesRegex(ValueError, 'outside explicit'):
            search.static_adaptation_audit(tree)

    def test_changed_symmetry_gate_is_detected(self):
        tree = ast.parse(Path(search.__file__).read_text())
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == '_symmetric_tensor')
        for node in ast.walk(function):
            if isinstance(node, ast.Constant) and node.value == 128:
                node.value = 1024
        with self.assertRaisesRegex(ValueError, 'symmetry check changed'):
            search.static_adaptation_audit(tree)

    def test_tangent_contact_step_change_is_detected(self):
        tree = ast.parse(textwrap.dedent(inspect.getsource(model.ProcessTrial.tangent))).body[0]
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == 1e-9:
                node.value = 1e-7
        with self.assertRaises(ValueError):
            model.operator_reuse_audit(tree)


if __name__ == '__main__':
    unittest.main(verbosity=2)

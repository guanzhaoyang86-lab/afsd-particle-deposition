"""Static adaptation, independent origin, and in-domain complete-map equality."""
import ast
import copy
import unittest
import numpy as np
import afsd_deposition_remote_setup_20261007 as model
import afsd_process_framework_20261007 as reference_model
import run_afsd_process_framework_20261007 as shared


class RemoteSetupTests(unittest.TestCase):
    def test_independent_origin_has_new_numerical_protocol(self):
        pa, state = model.independent_initial_state()
        self.assertEqual(len(pa.x), 1381)
        self.assertEqual(state['time_s'], 0.)
        self.assertEqual(state['step'], 0)
        self.assertEqual(state['numerical_domain_protocol'], model.NUMERICAL_PROTOCOL)
        self.assertFalse(state['accepted'])
        self.assertFalse(state['restart_allowed'])
        self.assertIsNone(state['parent_checkpoint'])
        model.ProcessTrial(pa, state, .0005)
        state.pop('numerical_domain_protocol')
        with self.assertRaises(ValueError):
            model.ProcessTrial(pa, state, .0005)

    def test_physical_equations_and_providers_are_unchanged(self):
        audit = model.operator_reuse_audit()
        self.assertTrue(audit['endpoint_equations_ast_identical_after_removing_only_remote_edge_assertion'])
        self.assertTrue(audit['endpoint_filter_unchanged'])
        self.assertFalse(audit['runtime_monkeypatch'])

    def test_old_and_new_complete_maps_agree_in_existing_domain(self):
        # Two unsolved residual maps on an already qualified existing geometry;
        # this is a component equivalence benchmark, not a new AFSD time step.
        pa, state = reference_model.independent_initial_state()
        trial = reference_model.ProcessTrial(pa, state, .0005)
        u = trial.prescribed(.0005)
        original = trial.evaluate(u, .0005)
        adapted = model.ProcessTrial.evaluate(trial, u, .0005)
        for key in ('position', 'displacement', 'i', 'j', 'g', 'density', 'gradient',
                    'force', 'stress', 'ep', 'temperature'):
            np.testing.assert_array_equal(original[key], adapted[key])
        differences = shared.batch.candidate_equivalence(
            dict(original, energy_comparison=shared.core.energy_diagnostics(trial, original)),
            dict(adapted, energy_comparison=shared.core.energy_diagnostics(trial, adapted)))
        self.assertTrue(differences)
        np.testing.assert_array_equal(original['thermal']['pair_conductance_W_K'],
                                      adapted['thermal']['pair_conductance_W_K'])

    def test_event_router_is_only_a_static_domain_adaptation(self):
        source = ast.parse((model.ROOT/'scripts/visible_support_events_20261007.py').read_text(encoding='utf-8'))
        actual = ast.parse((model.ROOT/'scripts/visible_support_remote_events_20261007.py').read_text(encoding='utf-8'))
        for node in actual.body:
            if isinstance(node, ast.ImportFrom) and node.module == 'visible_support_remote_domain_20261007':
                node.module = 'visible_support_numerics'
                node.names = [name for name in node.names if name.name != 'require_remote_filter_consistency']
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'CLASSIFICATION' for t in node.targets):
                node.value = ast.Constant('OBSERVED_VISIBLE_SUPPORT_EVENT_ROUTER_20261007')
        class Normalize(ast.NodeTransformer):
            def visit_Assign(self, node):
                if (len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and
                        node.targets[0].id == 'numerical_domain'):
                    return ast.Expr(value=node.value)
                return self.generic_visit(node)
            def visit_Expr(self, node):
                if isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == 'require_remote_filter_consistency':
                    return None
                return self.generic_visit(node)
        actual = Normalize().visit(actual)
        self.assertEqual(ast.dump(source, include_attributes=False), ast.dump(actual, include_attributes=False))


if __name__ == '__main__':
    unittest.main()

"""Read-only AST audit of the GPU search adapter; no physics or GPU execution.

The full module must consist only of the listed explicit imports and four frozen
CPU functions. The only normalized changes are the three dense linalg calls,
the Newton entry-point name, and its explicitly retained original journal.
"""
from pathlib import Path
import ast
import copy
import unittest

ROOT = Path(__file__).resolve().parents[1]
ADAPTER = ROOT/'scripts/afsd_gpu_equilibrium_20261007.py'
FUNCTIONS = {
    '_full_matrix': ('scripts/secant_equilibrium.py', '_full_matrix', 'inv'),
    'solve': ('scripts/visibility_equilibrium.py', 'solve', 'solve'),
    'newton_direction': ('src/afsd_backend/pysph/quasistatic_process.py', 'newton_direction', 'solve'),
    'newton_solve': ('scripts/run_discrete_contact_process.py', 'solve', None),
}
IMPORTS = '''from pathlib import Path
import time
import numpy as np
import afsd_gpu_linalg_20261007 as gpu_linalg
from secant_equilibrium import (MAX_REFRESHES, MAX_EVALUATIONS, SearchBudgetError,
    _condition, _broyden, summary, write)
from visibility_equilibrium import CONTROLLER, journal
from run_discrete_contact_process import journal as newton_journal
'''


def _dump(node):
    return ast.dump(node, include_attributes=False)


def _function(path, name):
    tree = ast.parse((ROOT/path).read_text(encoding='utf-8'))
    return next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name)


def audit_gpu_solver_adaptation(new_tree=None):
    """Fail if any unregistered execution, gate, budget, or binding changes."""
    tree = ast.parse(ADAPTER.read_text(encoding='utf-8')) if new_tree is None else copy.deepcopy(new_tree)
    if tree.body and isinstance(tree.body[0], ast.Expr) and isinstance(tree.body[0].value, ast.Constant) and isinstance(tree.body[0].value.value, str):
        tree.body.pop(0)
    imports = ast.parse(IMPORTS).body
    if len(tree.body) != len(imports)+len(FUNCTIONS):
        raise ValueError('unregistered adapter top-level execution or declaration')
    if [_dump(node) for node in tree.body[:len(imports)]] != [_dump(node) for node in imports]:
        raise ValueError('unregistered imported binding; CPU helpers and explicit GPU module must be retained')
    counts = {}
    for node, (new_name, (path, old_name, operation)) in zip(tree.body[len(imports):], FUNCTIONS.items()):
        if not isinstance(node, ast.FunctionDef) or node.name != new_name:
            raise ValueError('unexpected static GPU function set/order')
        original = _function(path, old_name)
        node.name = old_name
        replacements = {'dense': 0, 'journal': 0}

        class Normalize(ast.NodeTransformer):
            def visit_Call(self, call):
                function = call.func
                if isinstance(function, ast.Attribute) and isinstance(function.value, ast.Name) and function.value.id == 'gpu_linalg':
                    if operation is None or function.attr != operation:
                        raise ValueError('unregistered GPU operation')
                    function.value = ast.Attribute(value=ast.Name(id='np', ctx=ast.Load()), attr='linalg', ctx=ast.Load())
                    replacements['dense'] += 1
                if isinstance(function, ast.Name) and function.id == 'newton_journal':
                    if new_name != 'newton_solve':
                        raise ValueError('original journal binding changed')
                    function.id = 'journal'
                    replacements['journal'] += 1
                return self.generic_visit(call)

        normalized = Normalize().visit(node)
        expected_journals = sum(isinstance(item, ast.Call) and isinstance(item.func, ast.Name) and item.func.id == 'journal' for item in ast.walk(original)) if new_name == 'newton_solve' else 0
        expected = {'dense': int(operation is not None), 'journal': expected_journals}
        if replacements != expected:
            raise ValueError('each registered GPU/journal boundary must be replaced exactly once: '+new_name)
        if _dump(normalized) != _dump(original):
            raise ValueError('unregistered control/gate/budget change: '+new_name)
        counts[new_name] = replacements
    return dict(status='PASS_STATIC_GPU_DENSE_LINALG_ADAPTATION_ONLY', replacements=counts,
        dense_inverse_calls=1, dense_solve_calls=2, runtime_global_rebinding=False,
        unchanged='physical maps, tangent assembly, CPU Broyden updates, condition/backward checks, graph policy, search and budgets',
        scope='Static structure only; numerical backend and same-parent physical GPU comparisons require separate tests.')


class GPUEquilibriumAdaptationTests(unittest.TestCase):
    def tree(self):
        return ast.parse(ADAPTER.read_text(encoding='utf-8'))

    def function(self, tree, name):
        return next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name)

    def test_only_registered_static_boundaries_differ(self):
        result = audit_gpu_solver_adaptation()
        self.assertEqual(result['status'], 'PASS_STATIC_GPU_DENSE_LINALG_ADAPTATION_ONLY')
        self.assertEqual(result['dense_inverse_calls'], 1)
        self.assertEqual(result['dense_solve_calls'], 2)

    def test_changed_broyden_budget_is_detected(self):
        tree = self.tree()
        self.function(tree, 'solve').args.defaults[-1] = ast.Constant(value=25)
        with self.assertRaisesRegex(ValueError, 'control/gate/budget'):
            audit_gpu_solver_adaptation(tree)

    def test_changed_newton_backward_threshold_is_detected(self):
        tree = self.tree()
        for node in ast.walk(self.function(tree, 'newton_direction')):
            if isinstance(node, ast.Constant) and node.value == 1e-9:
                node.value = 1e-8
                break
        else:
            self.fail('original Newton backward threshold missing')
        with self.assertRaisesRegex(ValueError, 'control/gate/budget'):
            audit_gpu_solver_adaptation(tree)

    def test_hidden_cpu_solve_is_detected(self):
        tree = self.tree()
        for node in ast.walk(self.function(tree, 'solve')):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name) and node.func.value.id == 'gpu_linalg':
                node.func.value = ast.Attribute(value=ast.Name(id='np', ctx=ast.Load()), attr='linalg', ctx=ast.Load())
                break
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            audit_gpu_solver_adaptation(tree)

    def test_cpu_rank_one_binding_change_is_detected(self):
        tree = self.tree()
        for node in tree.body:
            if isinstance(node, ast.ImportFrom) and node.module == 'secant_equilibrium':
                for entry in node.names:
                    if entry.name == '_broyden':
                        entry.name = '_full_matrix'
        with self.assertRaisesRegex(ValueError, 'imported binding'):
            audit_gpu_solver_adaptation(tree)

    def test_top_level_rebinding_is_detected(self):
        tree = self.tree()
        tree.body.append(ast.parse('gpu_linalg.solve = np.linalg.solve').body[0])
        with self.assertRaisesRegex(ValueError, 'top-level execution'):
            audit_gpu_solver_adaptation(tree)

    def test_newton_original_journal_is_required(self):
        tree = self.tree()
        for node in ast.walk(self.function(tree, 'newton_solve')):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'newton_journal':
                node.func.id = 'journal'
                break
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            audit_gpu_solver_adaptation(tree)

    def test_graph_change_guard_is_preserved(self):
        tree = self.tree()
        for node in ast.walk(self.function(tree, 'solve')):
            if isinstance(node, ast.Constant) and node.value == 'ACTUAL_SUPPORT_GRAPH_CHANGED':
                node.value = 'IGNORE_GRAPH_CHANGE'
                break
        else:
            self.fail('original graph guard missing')
        with self.assertRaisesRegex(ValueError, 'control/gate/budget'):
            audit_gpu_solver_adaptation(tree)


if __name__ == '__main__':
    unittest.main()

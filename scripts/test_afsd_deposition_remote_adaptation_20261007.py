"""Exact AST boundary audit for the static remote-domain runner adaptation.

No solver is imported or executed. Only the model import, explicit path/event
audit calls/imports, default receipt and two added test requirements may differ.
All force/thermal gates, comparisons, controller, time/budget limits, history
selection, failure behavior, and file/source checks must remain identical.
"""
from pathlib import Path
import ast
import copy
import unittest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'scripts/run_afsd_deposition_stage_20261007.py'
REMOTE = ROOT/'scripts/run_afsd_deposition_stage_remote_20261007.py'
NEW_TESTS = ('scripts/test_visible_support_remote_domain_20261007.py',
             'scripts/test_afsd_deposition_remote_adaptation_20261007.py')


def audit_runner_adaptation(old_tree=None, new_tree=None):
    original = ast.parse(BASE.read_text(encoding='utf-8')) if old_tree is None else copy.deepcopy(old_tree)
    adapted = ast.parse(REMOTE.read_text(encoding='utf-8')) if new_tree is None else copy.deepcopy(new_tree)
    counts = dict(model_import=0, path_import=0, events_import=0, default_receipt=0,
                  required_tests=0, path_call=0, events_call=0)
    class Normalize(ast.NodeTransformer):
        def visit_Import(self, node):
            for entry in node.names:
                if entry.name == 'afsd_deposition_remote_setup_20261007':
                    if entry.asname != 'model':
                        raise ValueError('remote setup must retain the explicit model alias')
                    entry.name = 'afsd_deposition_setup_20261007'; counts['model_import'] += 1
            return node

        def visit_ImportFrom(self, node):
            expected = {
                'visible_support_remote_domain_20261007': ('audit_visible_support_path', 'remote_audit_visible_support_path', 'path_import'),
                'visible_support_remote_events_20261007': ('audit_visibility_events', 'remote_audit_visibility_events', 'events_import'),
            }
            if node.module in expected:
                name, alias, counter = expected[node.module]
                if node.level != 0 or len(node.names) != 1 or (node.names[0].name, node.names[0].asname) != (name, alias):
                    raise ValueError('unexpected remote audit import')
                counts[counter] += 1
                return None
            return node

        def visit_Assign(self, node):
            names = [v.id for v in node.targets if isinstance(v, ast.Name)]
            if names == ['DEFAULT_CHECKS']:
                if not (isinstance(node.value, ast.BinOp) and isinstance(node.value.op, ast.Div)
                        and isinstance(node.value.right, ast.Constant)
                        and node.value.right.value == 'outputs/deposition_campaign_20261007/components_remote_v1.json'):
                    raise ValueError('unexpected remote default receipt')
                node.value.right.value = 'outputs/deposition_stages_20261007/components_v1.json'
                counts['default_receipt'] += 1
            if names == ['required_tests']:
                if not isinstance(node.value, ast.Tuple) or tuple(ast.literal_eval(v) for v in node.value.elts[-2:]) != NEW_TESTS:
                    raise ValueError('remote domain and adaptation tests must be required')
                node.value.elts = node.value.elts[:-2]
                counts['required_tests'] += 1
            return self.generic_visit(node)

        def visit_Call(self, node):
            if isinstance(node.func, ast.Name):
                mapping = {'remote_audit_visible_support_path': ('audit_visible_support_path', 'path_call'),
                           'remote_audit_visibility_events': ('audit_visibility_events', 'events_call')}
                if node.func.id in mapping:
                    attr, counter = mapping[node.func.id]
                    node.func = ast.Attribute(value=ast.Name(id='shared', ctx=ast.Load()), attr=attr, ctx=ast.Load())
                    counts[counter] += 1
            return self.generic_visit(node)
    normalized = Normalize().visit(adapted)
    if counts != {key: 1 for key in counts}:
        raise ValueError('remote adaptation did not replace each registered boundary exactly once: '+str(counts))
    if ast.dump(original, include_attributes=False) != ast.dump(normalized, include_attributes=False):
        raise ValueError('unregistered runner control/gate/history/budget change')
    return dict(status='PASS_STATIC_REMOTE_RUNNER_ADAPTATION_ONLY', replacements=counts,
        old_runner_modified=False, runtime_global_rebinding=False,
        unchanged='all candidate gates, adaptive comparisons, selection/history, continuation and budgets',
        scope='Runner AST only; remote geometry proof and evaluate adaptation require separate tests.')


class RemoteRunnerAdaptationTests(unittest.TestCase):
    def test_only_registered_static_boundaries_differ(self):
        self.assertEqual(audit_runner_adaptation()['status'], 'PASS_STATIC_REMOTE_RUNNER_ADAPTATION_ONLY')

    def test_changed_budget_is_detected(self):
        tree = ast.parse(REMOTE.read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == 450.:
                node.value = 451.; break
        else:
            self.fail('registered forecast limit missing')
        with self.assertRaisesRegex(ValueError, 'unregistered'):
            audit_runner_adaptation(new_tree=tree)

    def test_missing_new_domain_test_is_detected(self):
        tree = ast.parse(REMOTE.read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == NEW_TESTS[0]:
                node.value = 'scripts/unrelated.py'; break
        with self.assertRaisesRegex(ValueError, 'must be required'):
            audit_runner_adaptation(new_tree=tree)

    def test_bypassing_new_path_audit_is_detected(self):
        tree = ast.parse(REMOTE.read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'remote_audit_visible_support_path':
                node.func.id = 'fake_pass'; break
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            audit_runner_adaptation(new_tree=tree)


if __name__ == '__main__':
    unittest.main()

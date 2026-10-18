"""Only the certificate provider, registrations and read-only cost journal differ."""
from pathlib import Path
import ast
import inspect
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_visibility_event_secant_20261010 as previous
import run_gpu_visibility_event_box_20261010 as current
import afsd_event_blocked_box_20261010 as certificate


class FrozenPhysicalDriver(unittest.TestCase):
    def test_entire_driver_diff_is_only_provider_registration_and_cost_journal(self):
        old = ast.parse(Path(previous.__file__).read_text(encoding='utf-8'))
        new = ast.parse(Path(current.__file__).read_text(encoding='utf-8'))
        changed_imports = 0
        for node in new.body:
            if isinstance(node, ast.ImportFrom) and node.module == 'afsd_event_blocked_box_20261010':
                node.module = 'afsd_event_blocked_path_20261010'
                changed_imports += 1
        self.assertEqual(changed_imports, 1)
        solve = next(n for n in new.body if isinstance(n, ast.FunctionDef) and n.name == 'solve_physical')
        journals = [n for n in solve.body if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
                    and any(isinstance(q, ast.Constant) and q.value == 'physical_solve_cost.json'
                            for q in ast.walk(n))]
        self.assertEqual(len(journals), 1)
        journal = journals[0].value
        self.assertEqual(ast.unparse(journal.func), 'core.write')
        self.assertEqual(ast.unparse(journal.args[0]), "directory / 'physical_solve_cost.json'")
        fields = {k.arg: k.value for k in journal.args[1].keywords}
        self.assertEqual(ast.unparse(fields['wall_seconds']), 'time.perf_counter() - began')
        self.assertEqual(ast.unparse(fields['checkpoint_sha256']), 'core.sha(target)')
        self.assertEqual(ast.unparse(fields['parent_fingerprint']), 'fingerprint')
        self.assertEqual(ast.unparse(fields['dt_s']), 'dt')
        self.assertEqual(ast.literal_eval(fields['status']), 'MEASURED_COMPLETED_PHYSICAL_SOLVE')
        self.assertFalse(ast.literal_eval(fields['accepted']))
        solve.body.remove(journals[0])
        registrations = {'CONTRACT', 'CHECKS', 'BLOCKED_PREFIX_PROTOCOL'}
        def strip(tree):
            tree.body[:] = [node for node in tree.body if not (
                isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in registrations)]
            return ast.dump(tree)
        self.assertEqual(strip(new), strip(old))

    def test_every_physical_provider_and_limit_is_frozen(self):
        for name in ('model','solver','gpu','cpu','core','post','event','batch','old','shared',
                     'locate_activation'):
            self.assertIs(getattr(current, name), getattr(previous, name))
        for name in ('PROTOCOL','LOCATOR_PROTOCOL','MAXIMUM_PHYSICAL_SOLVES',
                     'MAXIMUM_EVENTS_PER_BRANCH','EVENT_TIME_TOLERANCE_S',
                     'EVENT_RELATIVE_TIME_TOLERANCE','EVENT_POSITION_TOLERANCE_M'):
            self.assertEqual(getattr(current, name), getattr(previous, name))

    def test_registered_new_provider_has_same_api_and_default_budget(self):
        self.assertIs(current.audit_blocked_prefix, certificate.audit_blocked_prefix)
        self.assertEqual(inspect.signature(current.audit_blocked_prefix),
                         inspect.signature(previous.audit_blocked_prefix))
        self.assertEqual(current.CONTRACT, ROOT/'docs/GPU_VISIBILITY_EVENT_BOX_20261010.md')
        self.assertEqual(current.CHECKS, current.BASE/'event_box_components_20261010_v1.json')
        self.assertEqual(current.BLOCKED_PREFIX_PROTOCOL,
                         'FIXED_WITNESS_BOX_CONTINUOUS_INCOMING_BLOCKED_PREFIX_20261010_V1')


if __name__ == '__main__':
    unittest.main()

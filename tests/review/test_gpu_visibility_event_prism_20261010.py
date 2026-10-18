"""Only the explicitly registered visibility proof provider changes from v3."""
from pathlib import Path
import ast
import inspect
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_visibility_event_box_20261010 as previous
import run_gpu_visibility_event_prism_20261010 as current


class FrozenPhysics(unittest.TestCase):
    def test_entire_driver_is_unchanged_except_provider_and_registration(self):
        before = ast.parse(Path(previous.__file__).read_text(encoding='utf-8'))
        after = ast.parse(Path(current.__file__).read_text(encoding='utf-8'))
        imports = [n for n in after.body if isinstance(n, ast.Import)
                   and any(a.name == 'visible_support_prism_20261010' and a.asname == 'path_proof' for a in n.names)]
        self.assertEqual(len(imports), 1)
        self.assertEqual(len(imports[0].names), 1)
        after.body.remove(imports[0])
        qualify = next(n for n in after.body if isinstance(n, ast.FunctionDef) and n.name == 'qualify')
        changed = 0
        for node in ast.walk(qualify):
            if isinstance(node, ast.Call) and ast.unparse(node.func) == 'path_proof.audit_visible_support_path':
                node.func = ast.parse('original.remote_audit_visible_support_path', mode='eval').body
                changed += 1
        self.assertEqual(changed, 1)
        registrations = {'VISIBLE_PATH_PROTOCOL', 'CONTRACT', 'CHECKS'}
        def normalized(tree):
            tree.body[:] = [n for n in tree.body if not (isinstance(n, ast.Assign)
                and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)
                and n.targets[0].id in registrations)]
            return ast.dump(tree)
        self.assertEqual(normalized(before), normalized(after))

    def test_same_solver_physics_event_location_and_blocked_prefix(self):
        for name in ('model', 'solver', 'gpu', 'cpu', 'core', 'post', 'event', 'batch', 'old',
                     'shared', 'locate_activation', 'audit_blocked_prefix'):
            self.assertIs(getattr(current, name), getattr(previous, name))
        for name in ('PROTOCOL', 'LOCATOR_PROTOCOL', 'BLOCKED_PREFIX_PROTOCOL',
                     'MAXIMUM_PHYSICAL_SOLVES', 'MAXIMUM_EVENTS_PER_BRANCH',
                     'EVENT_TIME_TOLERANCE_S', 'EVENT_RELATIVE_TIME_TOLERANCE',
                     'EVENT_POSITION_TOLERANCE_M'):
            self.assertEqual(getattr(current, name), getattr(previous, name))

    def test_same_whole_path_api_and_registered_evidence_paths(self):
        self.assertEqual(inspect.signature(current.path_proof.audit_visible_support_path),
                         inspect.signature(previous.original.remote_audit_visible_support_path))
        self.assertEqual(current.CONTRACT, ROOT/'docs/GPU_VISIBILITY_EVENT_PRISM_20261010.md')
        self.assertEqual(current.CHECKS, current.BASE/'event_prism_components_20261010_v1.json')
        self.assertEqual(current.VISIBLE_PATH_PROTOCOL,
                         'SECOND_ORDER_SWEPT_SEGMENT_PRISM_SEPARATION_20261010_V1')


if __name__ == '__main__':
    unittest.main()

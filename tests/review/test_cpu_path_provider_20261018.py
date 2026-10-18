"""Freeze all CPU physics and receipt checks around an explicit proof provider."""
from pathlib import Path
import ast
import copy
import inspect
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))


class ProviderBoundary(unittest.TestCase):
    def test_driver_changes_only_proof_import_and_protocol(self):
        old = ast.parse((ROOT/'scripts/run_cpu_visibility_event_20261015.py').read_text(encoding='utf-8'))
        new = ast.parse((ROOT/'scripts/run_cpu_visibility_event_20261018.py').read_text(encoding='utf-8'))
        imports = registrations = 0
        for node in new.body:
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == 'visible_support_anchored_20261018':
                        self.assertEqual(alias.asname, 'path_proof')
                        alias.name = 'visible_support_prism_20261010'
                        imports += 1
            if isinstance(node, ast.Assign) and ast.unparse(node.targets[0]) == 'VISIBLE_PATH_PROTOCOL':
                self.assertEqual(node.value.value, 'ENDPOINT_ANCHORED_SWEPT_SEGMENT_CERTIFICATE_20261018_V1')
                node.value.value = 'SECOND_ORDER_SWEPT_SEGMENT_PRISM_SEPARATION_20261010_V1'
                registrations += 1
        self.assertEqual((imports, registrations), (1, 1))
        self.assertEqual(ast.dump(old), ast.dump(new))

    def test_receipt_checks_only_switch_explicit_driver(self):
        old = ast.parse((ROOT/'scripts/verify_cpu_event_20261015.py').read_text(encoding='utf-8'))
        new = ast.parse((ROOT/'scripts/verify_cpu_event_20261018.py').read_text(encoding='utf-8'))
        count = 0
        for node in ast.walk(new):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == 'run_cpu_visibility_event_20261018':
                        self.assertEqual(alias.asname, 'physical')
                        alias.name = 'run_cpu_visibility_event_20261015'
                        count += 1
        self.assertEqual(count, 1)
        self.assertEqual(ast.dump(old), ast.dump(new))

    def test_physics_and_event_parameters_identical(self):
        import run_cpu_visibility_event_20261015 as old
        import run_cpu_visibility_event_20261018 as new
        for name in ('model', 'solver', 'cpu', 'core', 'post', 'event', 'batch', 'old', 'shared', 'audit_blocked_prefix'):
            self.assertIs(getattr(new, name), getattr(old, name), name)
        for name in ('PROTOCOL', 'LOCATOR_PROTOCOL', 'BLOCKED_PREFIX_PROTOCOL', 'MAXIMUM_PHYSICAL_SOLVES',
                     'MAXIMUM_EVENTS_PER_BRANCH', 'EVENT_TIME_TOLERANCE_S', 'EVENT_RELATIVE_TIME_TOLERANCE',
                     'EVENT_POSITION_TOLERANCE_M'):
            self.assertEqual(getattr(new, name), getattr(old, name), name)
        self.assertEqual(inspect.signature(new.path_proof.audit_visible_support_path),
                         inspect.signature(old.path_proof.audit_visible_support_path))

    def test_existing_adaptation_audit_still_passes(self):
        import run_cpu_visibility_event_20261018 as new
        report = new.adaptation_audit()
        self.assertEqual(report['status'], 'PASS_CPU_EXPLICIT_COUPLED_EVENT_INTEGRATION')
        self.assertFalse(report['gpu_provider_used'])
        self.assertTrue(report['endpoint_policy_unchanged'])


if __name__ == '__main__':
    unittest.main()

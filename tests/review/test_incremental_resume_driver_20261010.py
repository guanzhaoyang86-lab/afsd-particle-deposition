"""Resume administration keeps the previously qualified controller unchanged."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_incremental_resume_20261010 as driver


class ResumeDriverTests(unittest.TestCase):
    def test_entire_driver_only_changes_registered_administrative_names(self):
        source = (ROOT/'scripts/run_gpu_incremental_20261010.py').read_text(encoding='utf-8')
        replacements = {
            'run_gpu_incremental_20261010.py':'run_gpu_incremental_resume_20261010.py',
            'afsd_incremental_worker_20261010':'afsd_incremental_worker_resume_20261010',
            'GPU_INCREMENTAL_CONTINUATION_20261010.md':'GPU_INCREMENTAL_RESUME_20261010.md',
            'incremental_components_20261010_v1.json':'incremental_resume_components_20261010_v1.json',
            'LIVE_WINDOWS_LOCKED_HISTORY_INCREMENTAL_VERIFICATION_20261010':
                'LIVE_WINDOWS_LOCKED_HISTORY_RESUME_BOUNDED_GRANT_20261010'}
        for old, new in replacements.items():
            self.assertIn(old, source)
            source = source.replace(old, new)
        actual = (ROOT/'scripts/run_gpu_incremental_resume_20261010.py').read_text(encoding='utf-8')
        self.assertEqual(ast.dump(ast.parse(source)), ast.dump(ast.parse(actual)))

    def test_original_preflight_and_physical_pair_static_audit(self):
        audit = driver.static_audit()
        self.assertEqual(audit['status'], 'PASS_UNCHANGED_INITIAL_RECEIPT_ASSERTIONS_AND_PHYSICAL_PAIR')
        self.assertTrue(audit['inherited_physical_audit']['pair_gates_ast_identical'])
        self.assertTrue(audit['inherited_physical_audit']['complete_candidate_ast_identical'])
        self.assertFalse(audit['runtime_global_rebinding'])

    def test_replay_uses_exact_complete_parent_with_unchanged_physical_source(self):
        base = driver.BASE/'incremental_deposition_20261010_v1'
        previous = json.loads((base/'attempt_0055/manifest.json').read_text(encoding='utf-8'))
        target = json.loads((base/'attempt_0056/manifest.json').read_text(encoding='utf-8'))
        self.assertEqual(previous['status'], driver.CONTINUE_PASS)
        self.assertEqual(target['status'], driver.CONTINUE_PASS)
        self.assertEqual(target['parent_sha256'], previous['screened_continuation_sha256'])
        self.assertEqual(target['dt_s'], previous['dt_s'])
        self.assertEqual(target['source_sha256'], previous['source_sha256'])
        self.assertAlmostEqual(target['screened_end_time_s'], .05500244140624995)
        self.assertEqual(target['screened_end_step'], 370)
        for folder, meta in ((base/'attempt_0055', previous), (base/'attempt_0056', target)):
            supervision = json.loads((folder/'supervision.json').read_text(encoding='utf-8'))
            self.assertEqual(supervision['exit_code'], 0)
            self.assertTrue(supervision['source_inputs_still_match'])
            self.assertEqual(hashlib.sha256((folder/'manifest.json').read_bytes()).hexdigest(),
                supervision['manifest_sha256'])
            self.assertEqual(hashlib.sha256((ROOT/meta['screened_continuation_checkpoint']).read_bytes()).hexdigest(),
                meta['screened_continuation_sha256'])

    def test_incomplete_57_is_not_mistaken_for_complete_parent(self):
        folder = driver.BASE/'incremental_deposition_20261010_v1/attempt_0057'
        meta = json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
        supervision = json.loads((folder/'supervision.json').read_text(encoding='utf-8'))
        self.assertNotEqual(meta['status'], driver.CONTINUE_PASS)
        self.assertFalse(meta['screened_continuation_authorized'])
        self.assertNotEqual(supervision['exit_code'], 0)
        self.assertIn('did not authorize worker', meta['error'])
        self.assertTrue(meta['pair_screen']['passed'])


if __name__ == '__main__':
    unittest.main(verbosity=2)

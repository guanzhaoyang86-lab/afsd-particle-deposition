"""Independent administrative preflight-cache tests, without physical solves.

Only synthetic files are changed. Frozen production modules, real receipts,
checkpoints and the running campaign are never modified or executed here.
"""
from pathlib import Path
import ast
import copy
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'scripts/run_gpu_cached_preflight_20261008.py'
sys.path.insert(0, str(ROOT/'scripts'))
SPEC = importlib.util.spec_from_file_location('cached_preflight_under_test', SOURCE)
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)


def digest(data):
    return hashlib.sha256(data).hexdigest()


class CacheIntegrityTests(unittest.TestCase):
    def setUp(self):
        parent = ROOT/'tmp'; parent.mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix='preflight_cache_test_', dir=parent)
        self.folder = Path(self.temp.name).resolve()
        self.assertTrue(self.folder.is_relative_to(parent.resolve()))
        self.path = self.folder/'evidence.bin'; self.path.write_bytes(b'unchanged evidence')
        self.cache = adapter.CachedPreflight(root=self.folder)

    def tearDown(self):
        self.temp.cleanup()

    def test_repeated_and_aliased_hash_access_agree(self):
        expected = digest(self.path.read_bytes())
        self.assertEqual(self.cache.sha(self.path), expected)
        self.assertEqual(self.cache.sha(Path('evidence.bin')), expected)
        self.assertEqual(self.cache.sha(self.folder/'unused'/ '..'/'evidence.bin'), expected)
        self.cache.verify({'evidence.bin': expected})
        stats = self.cache.finish()
        self.assertEqual(stats['unique_files'], 1)
        self.assertEqual(stats['byte_hash_reads'], 2)
        self.assertGreater(stats['hash_requests'], 2)

    def test_conflicting_digest_is_rejected_after_a_cache_hit(self):
        expected = self.cache.sha(self.path)
        self.cache.verify({'evidence.bin': expected})
        with self.assertRaises(ValueError):
            self.cache.verify({str(self.path): '0'*64})

    def test_conflicting_aliases_in_one_mapping_cannot_overwrite(self):
        expected = digest(self.path.read_bytes())
        with self.assertRaises(ValueError):
            self.cache.verify({'evidence.bin': expected, str(self.path): 'f'*64})

    def test_content_change_during_preflight_is_rejected(self):
        self.cache.sha(self.path)
        self.path.write_bytes(b'changed evidence with a different size')
        with self.assertRaises((ValueError, RuntimeError)):
            self.cache.sha(self.path)

    def test_final_hash_detects_same_size_restored_timestamp_tamper(self):
        self.cache.sha(self.path)
        previous = self.path.stat()
        self.path.write_bytes(b'X'+self.path.read_bytes()[1:])
        os.utime(self.path, ns=(previous.st_atime_ns, previous.st_mtime_ns))
        with self.assertRaises((ValueError, RuntimeError)):
            self.cache.finish()

    def test_file_replacement_cannot_reuse_an_old_hash(self):
        self.cache.sha(self.path)
        replacement = self.folder/'replacement.bin'; replacement.write_bytes(b'replacement bytes')
        replacement.replace(self.path)
        with self.assertRaises((ValueError, RuntimeError)):
            self.cache.finish()

    def test_missing_file_is_not_served_from_cache(self):
        self.cache.sha(self.path); self.path.unlink()
        with self.assertRaises((OSError, ValueError, RuntimeError)):
            self.cache.finish()

    def test_json_cache_returns_independent_mutable_values(self):
        path = self.folder/'manifest.json'
        original = {'status': 'PASS', 'nested': {'ids': [1, 2]}, 'source_sha256': {'a': 'b'}}
        path.write_text(json.dumps(original), encoding='utf-8')
        first = self.cache.read(path)
        first['nested']['ids'].append(999)
        first['source_sha256'].clear()
        second = self.cache.read(path)
        self.assertEqual(second, original)
        self.assertIsNot(first, second)
        self.cache.finish()

    def test_json_read_alone_still_participates_in_final_hash_audit(self):
        path = self.folder/'receipt.json'; path.write_text('{"passed": true}', encoding='utf-8')
        self.cache.read(path)
        path.write_text('{"passed":false}', encoding='utf-8')
        with self.assertRaises((ValueError, RuntimeError)):
            self.cache.finish()

    def test_cache_is_not_shared_across_preflight_instances(self):
        old_hash = self.cache.sha(self.path); self.cache.finish()
        self.path.write_bytes(b'new independent input')
        fresh = adapter.CachedPreflight(root=self.folder)
        self.assertNotEqual(fresh.sha(self.path), old_hash)
        self.assertEqual(fresh.sha(self.path), digest(self.path.read_bytes()))
        fresh.finish()


class ReceiptAndGateTests(unittest.TestCase):
    setUp = CacheIntegrityTests.setUp
    tearDown = CacheIntegrityTests.tearDown

    def receipt(self, status, *, gpu=False, extra=None):
        folder = self.folder/'receipt'; folder.mkdir()
        artifact = folder/'artifact.bin'; artifact.write_bytes(b'checked diagnostic evidence')
        log = folder/'run.log'; log.write_text('completed', encoding='utf-8')
        meta = dict(status=status, classification=adapter.model.CLASSIFICATION,
            accepted=False, restart_allowed=False, diagnostic_only=True, accepted_time_advanced_s=0.,
            source_sha256={'evidence.bin': digest(self.path.read_bytes())},
            artifacts_sha256={'receipt/artifact.bin': digest(artifact.read_bytes())})
        if gpu:
            meta['reference_sha256'] = {}
        if extra:
            meta.update(extra)
        path = folder/'manifest.json'; path.write_text(json.dumps(meta), encoding='utf-8')
        supervision = dict(exit_code=0, source_inputs_still_match=True,
            manifest_sha256=digest(path.read_bytes()), log_sha256=digest(log.read_bytes()))
        (folder/'supervision.json').write_text(json.dumps(supervision), encoding='utf-8')
        return folder, meta

    def test_cpu_cache_hit_retains_callers_status_and_source_constraints(self):
        folder, meta = self.receipt(adapter.COST_PASS)
        result = self.cache.checked_receipt(folder, {adapter.COST_PASS}, meta['source_sha256'])
        result['source_sha256'].clear()
        self.assertEqual(self.cache.checked_receipt(folder, {adapter.COST_PASS})['source_sha256'], meta['source_sha256'])
        with self.assertRaises(ValueError):
            self.cache.checked_receipt(folder, {adapter.STAGE_PASS})
        with self.assertRaises(ValueError):
            self.cache.checked_receipt(folder, {adapter.COST_PASS}, {})

    def test_gpu_cache_hit_cannot_promote_benchmark_to_physical_parent(self):
        folder, _ = self.receipt(adapter.BENCH_PASS, gpu=True)
        self.cache.gpu_receipt(folder, {adapter.BENCH_PASS})
        with self.assertRaises(ValueError):
            self.cache.gpu_receipt(folder, {adapter.CONTINUE_PASS})

    def test_cached_receipt_artifact_tampering_fails_final_seal(self):
        folder, _ = self.receipt(adapter.BENCH_PASS, gpu=True)
        self.cache.gpu_receipt(folder, {adapter.BENCH_PASS})
        (folder/'artifact.bin').write_bytes(b'changed diagnostic evidence')
        with self.assertRaises(ValueError):
            self.cache.finish()

    def test_gpu_receipt_still_rejects_failed_physical_pair_screen(self):
        screen = dict(passed=False, global_screen={'pass': False}, local_screen={'passed': True}, crossing_ids_match=True)
        folder, _ = self.receipt(adapter.CONTINUE_PASS, gpu=True,
            extra=dict(screened_continuation_authorized=True, screened_history=[{}, {}], pair_screen=screen))
        with self.assertRaisesRegex(ValueError, 'complete passing GPU fine pair'):
            self.cache.gpu_receipt(folder, {adapter.CONTINUE_PASS})


class StaticBoundaryTests(unittest.TestCase):
    def test_original_contract_worker_and_model_provider_binding_are_retained(self):
        self.assertEqual(adapter.ORIGINAL_RUNNER, Path(adapter.original.__file__).resolve())
        self.assertEqual(adapter.CONTRACT, adapter.original.CONTRACT)
        self.assertIs(adapter.model.validate_state, adapter.cpu.model.validate_state)
        self.assertEqual(adapter.model.CLASSIFICATION, adapter.cpu.model.CLASSIFICATION)
        adapter.static_adaptation_audit()

    def test_removing_physical_pair_gate_fails_static_audit(self):
        tree = ast.parse(SOURCE.read_text(encoding='utf-8'))
        cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'CachedPreflight')
        method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == 'validate_gpu_screened_receipt')
        gate = next(node for node in method.body if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Name) and node.value.func.id == 'require')
        gate.value.args[0] = ast.Constant(value=True)
        with self.assertRaises(ValueError):
            adapter.static_adaptation_audit(tree)

    def test_loosened_forecast_budget_fails_static_audit(self):
        tree = ast.parse(SOURCE.read_text(encoding='utf-8'))
        cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'CachedPreflight')
        method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == 'preflight')
        changed = False
        for node in ast.walk(method):
            if isinstance(node, ast.Constant) and node.value == 450.:
                node.value = 451.; changed = True
        self.assertTrue(changed)
        with self.assertRaises(ValueError):
            adapter.static_adaptation_audit(tree)


if __name__ == '__main__':
    unittest.main()

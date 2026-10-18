"""Independent GPU candidate boundary audit; no CUDA or physics execution.

The copied candidate must retain the entire original qualification sequence.
Only its explicit model/search providers and Newton entry point may differ.
Tests mutate those boundaries to establish that gates cannot be dropped while
the copied functions are still labelled equivalent.
"""
from pathlib import Path
import ast
import copy
from types import SimpleNamespace
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ADAPTER = ROOT/'scripts/run_afsd_gpu_20261007.py'
REFERENCE = ROOT/'scripts/run_afsd_deposition_stage_remote_20261007.py'
PROVIDER_IMPORTS = {
    'cpu': 'run_afsd_deposition_stage_remote_20261007',
    'model': 'afsd_gpu_trial_20261007',
    'solver': 'afsd_gpu_equilibrium_20261007',
    'gpu': 'afsd_gpu_linalg_20261007',
}
INHERITED = ('core', 'post', 'event', 'batch', 'old', 'shared')


def _dump(node):
    return ast.dump(node, include_attributes=False)


def _function(tree, name):
    rows = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name]
    if len(rows) != 1:
        raise ValueError('one explicit candidate function required: '+name)
    return rows[0]


def audit_gpu_candidate_adaptation(new_tree=None):
    original = ast.parse(REFERENCE.read_text(encoding='utf-8'))
    actual = ast.parse(ADAPTER.read_text(encoding='utf-8')) if new_tree is None else copy.deepcopy(new_tree)
    bound = {}
    inherited_nodes = []
    protected = set(PROVIDER_IMPORTS) | set(INHERITED)
    for node in actual.body:
        if isinstance(node, ast.Import):
            for entry in node.names:
                alias = entry.asname or entry.name
                if alias in PROVIDER_IMPORTS:
                    if alias in bound or entry.name != PROVIDER_IMPORTS[alias]:
                        raise ValueError('candidate provider import changed')
                    bound[alias] = entry.name
        elif isinstance(node, ast.ImportFrom):
            if any((entry.asname or entry.name) in protected for entry in node.names):
                raise ValueError('candidate provider shadowed by from-import')
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            if node.name in protected:
                raise ValueError('candidate provider shadowed by declaration')
        elif isinstance(node, ast.Assign):
            names = {n.id for target in node.targets for n in ast.walk(target) if isinstance(n, ast.Name)}
            if names & protected:
                inherited_nodes.append(node)
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
            names = {n.id for n in ast.walk(node.target) if isinstance(n, ast.Name)}
            if names & protected:
                raise ValueError('candidate provider shadowed by assignment')
    if bound != PROVIDER_IMPORTS:
        raise ValueError('missing explicit CPU/GPU candidate provider')
    expected = ast.parse('core, post, event, batch, old, shared = cpu.core, cpu.post, cpu.event, cpu.batch, cpu.old, cpu.shared').body[0]
    if len(inherited_nodes) != 1 or _dump(inherited_nodes[0]) != _dump(expected):
        raise ValueError('original qualification providers rebound')
    replacements = {}
    for name in ('scalar_equivalence', 'complete_candidate'):
        reference = _function(original, name)
        changed = _function(actual, name)
        count = 0
        for node in ast.walk(changed):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == 'solver' and node.func.attr == 'newton_solve'):
                node.func.value.id, node.func.attr = 'core', 'solve'
                count += 1
        if count != int(name == 'complete_candidate') or _dump(changed) != _dump(reference):
            raise ValueError('candidate physical gate or retry boundary changed: '+name)
        replacements[name] = count
    return dict(status='PASS_INDEPENDENT_GPU_CANDIDATE_BOUNDARIES', replacements=replacements,
        no_physics_executed=True, no_cuda_executed=True,
        unchanged='force/material/history advance, scalar equivalence, geometry/shell/path, energy, output and diagnostic flags')


class GPUCandidateBoundaryTests(unittest.TestCase):
    def tree(self):
        return ast.parse(ADAPTER.read_text(encoding='utf-8'))

    def test_original_gates_with_explicit_gpu_providers(self):
        self.assertEqual(audit_gpu_candidate_adaptation()['replacements']['complete_candidate'], 1)

    def test_removed_scalar_equivalence_is_detected(self):
        tree = self.tree(); function = _function(tree, 'complete_candidate')
        function.body = [node for node in function.body if not (
            isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name) and node.value.func.id == 'scalar_equivalence')]
        with self.assertRaisesRegex(ValueError, 'physical gate'):
            audit_gpu_candidate_adaptation(tree)

    def test_removed_path_failure_is_detected(self):
        tree = self.tree(); function = _function(tree, 'complete_candidate')
        for node in ast.walk(function):
            if isinstance(node, ast.If) and any(isinstance(v, ast.Name) and v.id == 'path' for v in ast.walk(node.test)):
                node.test = ast.Constant(value=False); break
        else:
            self.fail('path failure boundary missing')
        with self.assertRaisesRegex(ValueError, 'physical gate'):
            audit_gpu_candidate_adaptation(tree)

    def test_fresh_advance_bypass_is_detected(self):
        tree = self.tree(); function = _function(tree, 'complete_candidate')
        for node in ast.walk(function):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'advance_candidate':
                node.func.attr = 'reuse_without_fresh_evaluation'; break
        else:
            self.fail('candidate advance boundary missing')
        with self.assertRaisesRegex(ValueError, 'physical gate'):
            audit_gpu_candidate_adaptation(tree)

    def test_formal_acceptance_upgrade_is_detected(self):
        tree = self.tree(); function = _function(tree, 'complete_candidate')
        for node in ast.walk(function):
            if isinstance(node, ast.keyword) and node.arg == 'accepted':
                node.value = ast.Constant(value=True); break
        else:
            self.fail('diagnostic acceptance flag missing')
        with self.assertRaisesRegex(ValueError, 'physical gate'):
            audit_gpu_candidate_adaptation(tree)

    def test_cpu_fallback_provider_is_detected(self):
        tree = self.tree()
        for node in tree.body:
            if isinstance(node, ast.Import):
                for entry in node.names:
                    if entry.asname == 'solver':
                        entry.name = 'visibility_equilibrium'
        with self.assertRaisesRegex(ValueError, 'provider import'):
            audit_gpu_candidate_adaptation(tree)

    def test_extra_global_provider_rebinding_is_detected(self):
        tree = self.tree(); tree.body.append(ast.parse('model = cpu.model').body[0])
        with self.assertRaisesRegex(ValueError, 'providers rebound'):
            audit_gpu_candidate_adaptation(tree)

    def test_newton_retry_budget_change_is_detected(self):
        tree = self.tree(); function = _function(tree, 'complete_candidate')
        for node in ast.walk(function):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'newton_solve':
                node.args[-1] = ast.Constant(value=48); break
        with self.assertRaisesRegex(ValueError, 'physical gate'):
            audit_gpu_candidate_adaptation(tree)


class GPUReceiptBindingTests(unittest.TestCase):
    """Run the real receipt validator against isolated in-memory disk records."""
    def fixture(self):
        folder = ROOT/'outputs/gpu_receipt_test/stage'
        parent_path = ROOT/'outputs/gpu_receipt_test/parent.npz'
        screen = dict(passed=True, global_screen={'pass': True}, local_screen={'passed': True},
            crossing_ids_match=True, fine_checkpoint_sha256=['first-sha', 'last-sha'],
            coarse_checkpoint_sha256='coarse-sha')
        states, hashes = {}, {}
        rows = []
        for half, digest in ((1, 'first-sha'), (2, 'last-sha')):
            checkpoint = folder/('fine_'+str(half))/'NOT_RESTARTABLE_candidate.npz'
            relative = checkpoint.relative_to(ROOT).as_posix(); time_s = 1.+half*.125
            rows.append(dict(checkpoint=relative, sha256=digest, time_s=time_s))
            states[checkpoint.resolve()] = dict(diagnostic_only=True, accepted=False, restart_allowed=False,
                time_s=time_s, step=100+half, classification='test-model', accepted_time_advanced_s=0.,
                machine_history=['old']*100+['new-first', 'new-last'][:half],
                interval_energy={'records': ['old-energy']*100+['first-energy', 'last-energy'][:half]})
            hashes[checkpoint.resolve()] = digest
        hashes[(folder/'coarse/NOT_RESTARTABLE_candidate.npz').resolve()] = 'coarse-sha'
        hashes[parent_path.resolve()] = 'parent-sha'
        states[parent_path.resolve()] = dict(diagnostic_only=True, accepted=False, restart_allowed=False,
            time_s=1., step=100, classification='test-model', accepted_time_advanced_s=0.,
            machine_history=['old']*100, interval_energy={'records': ['old-energy']*100})
        meta = dict(screened_history=rows, pair_screen=screen, screened_continuation_authorized=True,
            screened_continuation_checkpoint=rows[-1]['checkpoint'], screened_continuation_sha256='last-sha',
            screened_end_time_s=1.25, screened_end_step=102, parent_time_s=1., dt_s=.25,
            parent_checkpoint=parent_path.relative_to(ROOT).as_posix(), parent_sha256='parent-sha',
            parent_fingerprint='parent-fingerprint')
        documents = {(folder/'pair_screen.json').resolve(): copy.deepcopy(screen)}
        def require(value, reason):
            if not value:
                raise ValueError(reason)
        namespace = dict(ROOT=ROOT, require=require, np=np,
            read=lambda path: copy.deepcopy(documents[Path(path).resolve()]),
            load=lambda path: (None, copy.deepcopy(states[Path(path).resolve()])),
            core=SimpleNamespace(sha=lambda path: hashes[Path(path).resolve()]),
            cpu=SimpleNamespace(parent_fingerprint=lambda pa, state: 'parent-fingerprint'),
            model=SimpleNamespace(validate_state=lambda state: None, CLASSIFICATION='test-model'))
        actual = ast.parse(ADAPTER.read_text(encoding='utf-8'))
        function = _function(actual, 'validate_gpu_screened_receipt')
        # Only the function definition is compiled. Imports, main and physics
        # from the runner never execute; its actual guard logic does execute.
        exec(compile(ast.Module(body=[function], type_ignores=[]), '<receipt-validator-test>', 'exec'), namespace)
        return namespace['validate_gpu_screened_receipt'], meta, folder, states, hashes, documents

    def test_exact_last_fine_receipt_passes(self):
        check, meta, folder, *_ = self.fixture()
        check(meta, folder)

    def test_coarse_checkpoint_cannot_seed_continuation(self):
        check, meta, folder, *_ = self.fixture()
        meta['screened_continuation_checkpoint'] = (folder/'coarse/NOT_RESTARTABLE_candidate.npz').relative_to(ROOT).as_posix()
        meta['screened_continuation_sha256'] = 'coarse-sha'
        with self.assertRaises(ValueError):
            check(meta, folder)

    def test_pair_screen_cannot_bind_different_fine_state(self):
        check, meta, folder, states, hashes, documents = self.fixture()
        meta['pair_screen']['fine_checkpoint_sha256'][1] = 'another-sha'
        documents[(folder/'pair_screen.json').resolve()] = copy.deepcopy(meta['pair_screen'])
        with self.assertRaises(ValueError):
            check(meta, folder)

    def test_failed_local_screen_cannot_seed_continuation(self):
        check, meta, folder, states, hashes, documents = self.fixture()
        meta['pair_screen']['local_screen']['passed'] = False
        documents[(folder/'pair_screen.json').resolve()] = copy.deepcopy(meta['pair_screen'])
        with self.assertRaises(ValueError):
            check(meta, folder)

    def test_checkpoint_file_hash_change_is_rejected(self):
        check, meta, folder, states, hashes, _ = self.fixture()
        hashes[(ROOT/meta['screened_continuation_checkpoint']).resolve()] = 'tampered-sha'
        with self.assertRaises(ValueError):
            check(meta, folder)

    def test_screened_metadata_cannot_misstate_endpoint(self):
        check, meta, folder, *_ = self.fixture()
        meta['screened_end_step'] += 1
        with self.assertRaises(ValueError):
            check(meta, folder)

    def test_half_step_timing_cannot_skip_midpoint(self):
        check, meta, folder, states, hashes, _ = self.fixture()
        meta['screened_history'][0]['time_s'] = 1.1
        states[(ROOT/meta['screened_history'][0]['checkpoint']).resolve()]['time_s'] = 1.1
        with self.assertRaises(ValueError):
            check(meta, folder)


if __name__ == '__main__':
    unittest.main(verbosity=2)

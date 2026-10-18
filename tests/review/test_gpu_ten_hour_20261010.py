"""The authorized ten-hour continuation changes administration, not physics."""
import argparse
import ast
import inspect
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch, Mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_ten_hour_20261010 as driver


class TenHourDriverTests(unittest.TestCase):
    def test_entire_driver_only_changes_registered_names_and_wall_budget(self):
        source = (ROOT/'scripts/run_gpu_incremental_resume_20261010.py').read_text(encoding='utf-8')
        replacements = {
            'run_gpu_incremental_resume_20261010.py':'run_gpu_ten_hour_20261010.py',
            'GPU_INCREMENTAL_RESUME_20261010.md':'GPU_TEN_HOUR_CONTINUATION_20261010.md',
            'incremental_resume_components_20261010_v1.json':'ten_hour_components_20261010_v1.json',
            'LIVE_WINDOWS_LOCKED_HISTORY_RESUME_BOUNDED_GRANT_20261010':
                'LIVE_WINDOWS_LOCKED_HISTORY_TEN_HOUR_20261010',
            '14400.':'36000.',
        }
        for old, new in replacements.items():
            self.assertIn(old, source)
            source = source.replace(old, new)
        actual = ast.parse((ROOT/'scripts/run_gpu_ten_hour_20261010.py').read_text(encoding='utf-8'))
        counts = dict(imports=0, helpers=0, bindings=0, recovery_calls=0, dispatch=0)
        class Normalize(ast.NodeTransformer):
            def visit_Import(self, node):
                if ast.dump(node) == ast.dump(ast.parse('import afsd_geometry_refinement_20261010 as geometry_refinement').body[0]):
                    counts['imports'] += 1
                    return None
                return node
            def visit_FunctionDef(self, node):
                if node.name == 'geometry_budget_recovery':
                    counts['helpers'] += 1
                    return None
                return self.generic_visit(node)
            def visit_List(self, node):
                for item in node.elts[:]:
                    if isinstance(item, ast.Constant) and item.value == 'scripts/afsd_geometry_refinement_20261010.py':
                        node.elts.remove(item); counts['bindings'] += 1
                return self.generic_visit(node)
            def visit_Assign(self, node):
                expected = ast.parse('recovery=geometry_budget_recovery(pool,child,failure,meta,policy)').body[0]
                if ast.dump(node) == ast.dump(expected):
                    counts['recovery_calls'] += 1
                    return None
                return self.generic_visit(node)
            def visit_If(self, node):
                if ast.dump(node.test) == ast.dump(ast.parse('recovery is None', mode='eval').body):
                    self_outer.assertEqual(ast.dump(node.orelse[0]), ast.dump(ast.parse('policy,event=recovery').body[0]))
                    self_outer.assertEqual(len(node.orelse), 1)
                    counts['dispatch'] += 1
                    return node.body
                return self.generic_visit(node)
        self_outer = self
        normalized = Normalize().visit(actual)
        self.assertEqual(counts, dict(imports=1, helpers=1, bindings=1, recovery_calls=1, dispatch=1))
        self.assertEqual(ast.dump(ast.parse(source)), ast.dump(normalized))

    def test_original_preflight_and_physical_pair_static_audit(self):
        audit = driver.static_audit()
        self.assertEqual(audit['status'], 'PASS_UNCHANGED_INITIAL_RECEIPT_ASSERTIONS_AND_PHYSICAL_PAIR')
        self.assertTrue(audit['inherited_physical_audit']['pair_gates_ast_identical'])
        self.assertTrue(audit['inherited_physical_audit']['complete_candidate_ast_identical'])
        self.assertFalse(audit['runtime_global_rebinding'])

    def test_cli_defaults_to_ten_hours(self):
        argv = ['run_gpu_ten_hour_20261010.py', '--output', 'unused-new-session', '--parent', 'unused-parent']
        with patch.object(sys, 'argv', argv), patch.object(driver, 'run', return_value=0) as run:
            self.assertEqual(driver.main(), 0)
        args = run.call_args.args[0]
        self.assertEqual(args.wall_seconds, 36000.)
        self.assertEqual(args.max_stages, 1000)
        self.assertEqual(args.maximum_dt, .000474609375)
        self.assertEqual(args.target_time, .22)
        self.assertEqual(args.shoulder_observe_ms, 5.)

    def test_cli_accepts_explicit_ten_hours(self):
        argv = ['run_gpu_ten_hour_20261010.py', '--output', 'unused-new-session',
                '--parent', 'unused-parent', '--wall-seconds', '36000']
        with patch.object(sys, 'argv', argv), patch.object(driver, 'run', return_value=0) as run:
            self.assertEqual(driver.main(), 0)
        self.assertEqual(run.call_args.args[0].wall_seconds, 36000.)

    def test_executed_wall_budget_guard_accepts_only_authorized_range(self):
        # Evaluate the actual leading run guard without starting history I/O or a GPU.
        tree = ast.parse(inspect.getsource(driver.run))
        guards = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                  and isinstance(n.func, ast.Name) and n.func.id == 'require'
                  and any(isinstance(v, ast.Attribute) and v.attr == 'wall_seconds'
                          for v in ast.walk(n.args[0]))]
        self.assertEqual(len(guards), 1)
        expression = compile(ast.Expression(guards[0].args[0]), '<actual-run-wall-guard>', 'eval')
        for seconds, allowed in ((.001, True), (14400., True), (36000., True),
                                 (36000.001, False), (0., False), (-1., False),
                                 (float('inf'), False), (float('nan'), False)):
            with self.subTest(wall_seconds=seconds):
                args = argparse.Namespace(wall_seconds=seconds, max_stages=1000)
                self.assertIs(eval(expression, {'args': args}), allowed)
        for count in (0, 2001):
            with self.subTest(max_stages=count):
                args = argparse.Namespace(wall_seconds=36000., max_stages=count)
                self.assertFalse(eval(expression, {'args': args}))

    def geometry_fixture(self):
        folder = driver.BASE/'incremental_deposition_resume_20261010_v1/attempt_0019'
        failure = driver.read(folder/'manifest.json')
        parent = driver.read(folder.parent/'attempt_0018/manifest.json')
        state = driver.adaptive.initial_state(parent['screened_continuation_sha256'], failure['dt_s'])
        state = driver.adaptive.begin_attempt(state, failure['dt_s'])
        pool = Mock()
        pool.sha.side_effect = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
        return folder, failure, parent, state, pool

    def test_real_sealed_budget_failure_retries_same_complete_parent(self):
        folder, failure, parent, state, pool = self.geometry_fixture()
        before = copy.deepcopy(state)
        smaller, event = driver.geometry_budget_recovery(pool, folder, failure, parent, state)
        self.assertEqual(event['action'], 'REFINE')
        self.assertEqual(event['candidate_role'], 'coarse')
        self.assertEqual(smaller['parent_sha256'], parent['screened_continuation_sha256'])
        self.assertEqual(smaller['nominal_dt_s'], failure['dt_s']/2)
        self.assertEqual(smaller['passes'], 0)
        self.assertEqual(state, before)
        self.assertEqual(pool.seal.call_count, 3)

    def test_changed_parent_source_interval_or_unbound_report_rejects(self):
        for kind in ('parent','source','interval','unbound'):
            folder, failure, parent, state, pool = self.geometry_fixture()
            with self.subTest(kind=kind):
                if kind == 'parent': parent['screened_continuation_sha256'] = 'a'*64
                elif kind == 'source': parent['source_sha256'] = {'not-qualified':'a'*64}
                elif kind == 'interval': state['pending']['dt_s'] /= 2
                else:
                    failure['artifacts_sha256'].pop(driver.rel(folder/'coarse/visible_path_audit.json'))
                if kind == 'unbound':
                    with patch.object(driver, 'read', side_effect=lambda path: failure if path.name=='manifest.json' else json.loads(path.read_text(encoding='utf-8'))):
                        with self.assertRaises(ValueError): driver.geometry_budget_recovery(pool,folder,failure,parent,state)
                else:
                    with self.assertRaises(ValueError): driver.geometry_budget_recovery(pool,folder,failure,parent,state)

    def test_non_geometry_failure_cannot_enter_new_retry_route(self):
        folder, failure, parent, state, pool = self.geometry_fixture()
        failure['error'] = "ValueError('CUDA failed')"
        self.assertIsNone(driver.geometry_budget_recovery(pool,folder,failure,parent,state))
        pool.sha.assert_not_called(); pool.seal.assert_not_called()


if __name__ == '__main__':
    unittest.main(verbosity=2)

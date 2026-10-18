"""CPU event-adapter invariants and real failed-interval classification.

These are source/provider, manufactured locator and rejection checks. They do
not execute a nonlinear AFSD solve or authorize a saved failed candidate.
"""
from pathlib import Path
import ast
import copy
import hashlib
import inspect
import json
import sys
import tempfile
import textwrap
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import afsd_cpu_event_trial_20261015 as model
import afsd_deposition_remote_setup_20261007 as cpu_model
import run_cpu_visibility_event_20261015 as current
import run_afsd_deposition_stage_remote_20261007 as cpu
import afsd_event_locator_delay_20261010 as locator
from test_event_locator_delay_20261010 import ForcedParticle


def function_ast(function):
    return ast.parse(textwrap.dedent(inspect.getsource(function))).body[0]


def inherited_ast(name):
    tree = ast.parse((ROOT/'scripts/run_gpu_visibility_event_delay_20261010.py').read_text(encoding='utf-8'))
    return next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name)


def ast_equal(test, actual, expected):
    test.assertEqual(ast.dump(actual, include_attributes=False),
                     ast.dump(expected, include_attributes=False))


class FrozenCpuPhysics(unittest.TestCase):
    def test_no_event_equations_are_exactly_the_original_cpu_equations(self):
        actual = function_ast(model.ProcessTrial.evaluate)
        extra = [node for node in actual.body if isinstance(node, ast.If)
                 and ast.unparse(node.test) == 'self.incoming_pid_pairs']
        self.assertEqual(len(extra), 1)
        old_event = ast.parse((ROOT/'scripts/afsd_gpu_event_trial_20261010.py').read_text(encoding='utf-8'))
        old_evaluate = next(node for node in ast.walk(old_event)
                            if isinstance(node, ast.FunctionDef) and node.name == 'evaluate')
        old_extra = next(node for node in old_evaluate.body if isinstance(node, ast.If)
                         and ast.unparse(node.test) == 'self.incoming_pid_pairs')
        ast_equal(self, extra[0], old_extra)
        actual.body.remove(extra[0])
        ast_equal(self, actual, function_ast(cpu_model.ProcessTrial.evaluate))

    def test_physical_providers_tangent_and_fresh_advance_are_cpu_objects(self):
        for name in ('operators', '_geometry', 'filter_visible_support', 'update',
                     'advance_implicit_thermal', 'pressure', 'internal_force',
                     'edge_fingerprint', 'validate_state', 'advance_candidate'):
            self.assertIs(getattr(model, name), getattr(cpu_model, name), name)
        self.assertIs(model.ProcessTrial.tangent, cpu_model.ProcessTrial.tangent)
        self.assertIs(model.ProcessTrial.prescribed, cpu_model.ProcessTrial.prescribed)
        for name in ('solver', 'core', 'post', 'event', 'batch', 'old', 'shared'):
            self.assertIs(getattr(current, name), getattr(cpu, name), name)
        self.assertTrue(model.operator_reuse_audit()['physical_equations_unchanged'])

    def test_one_selected_graph_feeds_density_material_force_and_heat(self):
        body = function_ast(model.ProcessTrial.evaluate)
        selected = [node for node in ast.walk(body) if isinstance(node, ast.Assign)
                    and ast.unparse(node.targets[0]) == '(i, j, g)']
        self.assertEqual(len(selected), 1)
        self.assertEqual(ast.unparse(selected[0].value), "(graph['i'], graph['j'], graph['g'])")
        calls = {ast.unparse(node.func): node for node in ast.walk(body) if isinstance(node, ast.Call)}
        for name, start in (('_geometry', 4), ('advance_implicit_thermal', 5), ('internal_force', 2)):
            self.assertEqual([ast.unparse(v) for v in calls[name].args[start:start+3]], ['i', 'j', 'g'])
        self.assertEqual(ast.unparse(calls['update'].args[2]), 'L')
        gradients = [node for node in ast.walk(body) if isinstance(node, ast.Assign)
                     and ast.unparse(node.targets[0]) == '(L, volume)']
        self.assertEqual(ast.unparse(gradients[0].value.elts[0]), 'H / dt')

    def test_inherited_geometry_and_event_qualification_are_exact(self):
        for name in ('resolved_routing', 'probe_prefix_is_ready', 'scalar_equivalence',
                     'event_probe_guess', 'routing', 'qualify', 'leaf_row',
                     'copy_macro_endpoint', 'classify_activation'):
            with self.subTest(name=name):
                actual = function_ast(getattr(current, name))
                if name == 'qualify':
                    additions = [node for node in actual.body if isinstance(node, ast.Assign)
                                 and isinstance(node.value, ast.Call)
                                 and ast.unparse(node.value.func) == 'bind_parent_metadata']
                    self.assertEqual(len(additions), 1)
                    wanted = ast.parse("path = bind_parent_metadata(path, item['parent_checkpoint'], "
                        "item['parent_sha256'], cpu.parent_fingerprint(trial.pa, trial.state))").body[0]
                    ast_equal(self, additions[0], wanted)
                    actual.body.remove(additions[0])
                ast_equal(self, actual, inherited_ast(name))

    def test_branch_all_original_event_gates_survive_metadata_and_budget_adaptation(self):
        actual, expected = function_ast(current.complete_branch), inherited_ast('complete_branch')
        self.assertEqual(ast.unparse(actual.body[0]),
            "maximum_events = bounded_budget_limit(budget, 'maximum_events_per_branch', MAXIMUM_EVENTS_PER_BRANCH)")
        self.assertIn('actual branch parent checkpoint SHA changed', ast.unparse(actual.body[1]))
        del actual.body[:2]
        removed = dict(item_binding=0, proof_binding=0, locator_binding=0, budget_name=0)
        class AdministrativeOnly(ast.NodeTransformer):
            def visit_Expr(self, node):
                if isinstance(node.value, ast.Call) and ast.unparse(node.value.func) in ('candidate.update', 'probe.update'):
                    keys = {keyword.arg: ast.unparse(keyword.value) for keyword in node.value.keywords}
                    if keys == dict(parent_checkpoint='parent_checkpoint', parent_sha256='parent_sha'):
                        removed['item_binding'] += 1
                        return None
                return self.generic_visit(node)
            def visit_Assign(self, node):
                if isinstance(node.value, ast.Call) and ast.unparse(node.value.func) == 'bind_parent_metadata':
                    removed['proof_binding'] += 1
                    return None
                return self.generic_visit(node)
            def visit_Call(self, node):
                if ast.unparse(node.func) == 'locate_activation':
                    binding = [keyword for keyword in node.keywords if keyword.arg == 'parent_sha256']
                    if len(binding) == 1 and ast.unparse(binding[0].value) == 'parent_sha':
                        node.keywords.remove(binding[0]); removed['locator_binding'] += 1
                return self.generic_visit(node)
            def visit_Name(self, node):
                if node.id == 'maximum_events' and isinstance(node.ctx, ast.Load):
                    node.id = 'MAXIMUM_EVENTS_PER_BRANCH'; removed['budget_name'] += 1
                return node
        actual = AdministrativeOnly().visit(actual)
        self.assertEqual(removed, dict(item_binding=2, proof_binding=2, locator_binding=1, budget_name=1))
        ast_equal(self, actual, expected)

    def test_physical_solve_only_tightens_budget_and_uses_original_cpu_newton(self):
        actual, expected = function_ast(current.solve_physical), inherited_ast('solve_physical')
        # Replace only the new bounded admission statements with the inherited
        # cap. All nonlinear, inventory, fresh-state and save/load gates remain.
        self.assertEqual(ast.unparse(actual.body[1].targets[0]), 'maximum_solves')
        self.assertIn('0 <=', ast.unparse(actual.body[2]))
        actual.body[1:3] = [copy.deepcopy(expected.body[1])]
        changed = 0
        for node in ast.walk(expected):
            if isinstance(node, ast.Call) and ast.unparse(node.func) == 'solver.newton_solve':
                node.func = ast.parse('core.solve', mode='eval').body
                changed += 1
        self.assertEqual(changed, 1)
        ast_equal(self, actual, expected)

    def test_pair_global_local_crossing_and_history_gates_are_exact(self):
        actual, expected = function_ast(current.run_pair), inherited_ast('run_pair')
        self.assertIn('actual pair parent checkpoint SHA changed', ast.unparse(actual.body[0]))
        actual.body.pop(0)
        # Provider labels are the only changes after the extra actual SHA gate.
        for node in ast.walk(actual):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                node.value = node.value.replace('CPU interval', 'GPU interval').replace(
                    'no CPU event continuation', 'no GPU continuation')
        ast_equal(self, actual, expected)


class SavedShortBarEvent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        folder = ROOT/'outputs/short_bar_restart_20261015_v1/run/short_4p0mm/stage_036/macro_0001/attempt_00/coarse'
        files = {
            'routing': ('visibility_event_routing.json', '32e920d7b2100c55fd6e39d21169e4d21447614410305635e0bb13f6e68c3dae'),
            'path': ('visible_path_audit.json', 'e8fae34678270c3ea2293118b3f6aec7904f0e2a91ecd1e2c15e9a3011b0e7ff'),
            'shell': ('continuous_shell_audit.json', '83c9e012a8194aee4a0af5211207d595bf7f3ff851829d592a0b4ab7a0e9e49a'),
            'geometry': ('geometry_audit.json', '2a89710b069e61760af304ca82c8afd37a3c4e6ab8aea20cd6afbbcf59df4019')}
        cls.fixture = {}
        for key, (name, expected) in files.items():
            raw = (folder/name).read_bytes()
            if hashlib.sha256(raw).hexdigest() != expected:
                raise AssertionError('Frozen short-bar failed-event fixture changed: '+name)
            cls.fixture[key] = json.loads(raw)

    def test_actual_short_bar_failure_is_single_resolved_incoming_pair(self):
        pairs, fraction = current.classify_activation(copy.deepcopy(self.fixture))
        self.assertEqual(pairs, ((1059, 1085),))
        self.assertGreater(fraction, 0.)
        self.assertLess(fraction, 1.)

    def test_outgoing_unknown_additional_and_unresolved_events_are_rejected(self):
        def outgoing(item): item['routing']['transition_brackets'][0]['direction'] = 'VISIBLE_TO_BLOCKED'
        def two_events(item): item['routing']['transition_brackets'] *= 2
        def unknown_pair(item): item['path']['found_edges'][0]['particle_id_i'] = 1
        def budget_failure(item): item['path']['failure_kind'] = 'VISIBLE_PATH_QUERY_BUDGET_EXHAUSTED'
        def unresolved(item): item['routing']['unresolved_samples'] = [dict(reason='UNKNOWN')]
        def shell_failure(item): item['shell']['status'] = 'UNRESOLVED'
        def missed_contact(item): item['geometry']['transient_contacts_missed_by_endpoints'] = [1]
        for mutate in (outgoing, two_events, unknown_pair, budget_failure, unresolved, shell_failure, missed_contact):
            item = copy.deepcopy(self.fixture)
            mutate(item)
            with self.subTest(case=mutate.__name__), self.assertRaises(ValueError):
                current.classify_activation(item)

    def test_tightened_one_event_budget_stops_before_a_second_locator(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory); parent = base/'parent.npz'; parent.write_bytes(b'fixed parent')
            with patch.object(current, 'solve_physical', return_value=copy.deepcopy(self.fixture)), \
                 patch.object(current, 'qualify', return_value=False), \
                 patch.object(current, 'locate_activation') as locate:
                with self.assertRaisesRegex(ValueError, 'event count limit'):
                    current.complete_branch(None, None, .001, base/'branch', float('inf'), 1., {},
                        dict(solves=0, maximum_events_per_branch=1), parent, current.core.sha(parent), depth=1)
                locate.assert_not_called()


class ParentIdentityAndBudgets(unittest.TestCase):
    def test_administrative_caps_cannot_enlarge_registered_limits(self):
        self.assertEqual(current.MAXIMUM_PHYSICAL_SOLVES, 25)
        self.assertEqual(current.MAXIMUM_EVENTS_PER_BRANCH, 2)
        for key, default in (('maximum_solves', 25), ('maximum_events_per_branch', 2)):
            self.assertEqual(current.bounded_budget_limit({}, key, default), default)
            self.assertEqual(current.bounded_budget_limit({key: 1}, key, default), 1)
            for value in (True, False, 0, -1, default+1, float(default), None, '1'):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    current.bounded_budget_limit({key: value}, key, default)

    def test_exhausted_tightened_budget_stops_before_trial_or_output(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)/'not_created'
            for budget in (dict(solves=12, maximum_solves=12), dict(solves=25),
                           dict(solves=-1), dict(solves=True), dict(solves=0, maximum_solves=26)):
                with self.subTest(budget=budget), patch.object(model, 'ProcessTrial') as trial:
                    with self.assertRaises(ValueError):
                        current.solve_physical(None, None, .001, target, float('inf'), 1., {}, budget)
                    trial.assert_not_called()
                    self.assertFalse(target.exists())

    def test_branch_mismatched_actual_checkpoint_sha_stops_before_solve(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory); parent = base/'parent.npz'; parent.write_bytes(b'actual parent')
            with patch.object(current, 'solve_physical') as solve, self.assertRaisesRegex(ValueError, 'SHA changed'):
                current.complete_branch(None, None, .001, base/'branch', float('inf'), 1., {},
                    dict(solves=0), parent, '0'*64)
            solve.assert_not_called()
            self.assertFalse((base/'branch').exists())

    def test_last_admitted_solve_consumes_the_tightened_budget(self):
        class ReachedTrialConstruction(RuntimeError): pass
        with tempfile.TemporaryDirectory() as directory:
            budget = dict(solves=11, maximum_solves=12)
            with patch.object(current.cpu, 'parent_fingerprint', return_value='fixed'), \
                 patch.object(model, 'ProcessTrial', side_effect=ReachedTrialConstruction), \
                 self.assertRaises(ReachedTrialConstruction):
                current.solve_physical(None, None, .001, Path(directory)/'probe', float('inf'), 1., {}, budget)
            self.assertEqual(budget['solves'], 12)

    def test_geometry_report_parent_metadata_preserves_all_proof_decisions(self):
        original = dict(status='UNRESOLVED', failure_kind='ACTUAL_FAILURE',
                        frozen_parent_sha256='b'*64, query_log=[dict(query=1)])
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)/'parent.npz'; parent.write_bytes(b'actual parent')
            sha = current.core.sha(parent)
            with patch.object(current, 'rel', side_effect=lambda path: str(path)):
                bound = current.bind_parent_metadata(original, parent, sha, 'actual fingerprint')
                with self.assertRaisesRegex(ValueError, 'SHA changed'):
                    current.bind_parent_metadata(original, parent, '0'*64, 'actual fingerprint')
            for key in ('status', 'failure_kind', 'query_log'):
                self.assertEqual(bound[key], original[key])
            self.assertEqual(bound['parent_sha256'], sha)
            self.assertEqual(bound['actual_parent_fingerprint'], 'actual fingerprint')
            self.assertEqual(bound['legacy_fixture_parent_sha256'], 'b'*64)
            self.assertFalse(bound['legacy_fixture_parent_is_current'])
            self.assertNotIn('frozen_parent_sha256', bound)
            self.assertEqual(original['frozen_parent_sha256'], 'b'*64)

    def test_pair_mismatched_actual_checkpoint_sha_stops_before_load(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)/'parent.npz'; parent.write_bytes(b'actual parent')
            with patch.object(current, 'load') as load, self.assertRaisesRegex(ValueError, 'SHA changed'):
                current.run_pair(Path(directory), dict(parent_checkpoint=str(parent), parent_sha256='0'*64), {}, float('inf'))
            load.assert_not_called()

    def test_locator_binds_actual_parent_without_changing_algorithm_or_payloads(self):
        before, after = ForcedParticle(), ForcedParticle()
        kwargs = dict(time_tolerance=1e-7, maximum_solves=10)
        reference = locator.locate_activation(before, .001, .6, **kwargs)
        result = current.locate_activation(after, .001, .6, parent_sha256='a'*64, **kwargs)
        self.assertEqual(before.calls, after.calls)
        for name in ('status', 'reason', 'solve_count', 'bracket_width_s',
                     'maximum_position_difference_m', 'high_activation_prefix_certified', 'accepted', 'restart_allowed'):
            self.assertEqual(result[name], reference[name], name)
        self.assertEqual(result['parent_sha256'], 'a'*64)
        self.assertEqual(result['legacy_fixture_parent_sha256'], reference['frozen_parent_sha256'])
        self.assertFalse(result['legacy_fixture_parent_is_current'])
        self.assertNotIn('frozen_parent_sha256', result)
        self.assertIs(result['hi']['payload'], after.payloads[after.calls.index(result['hi']['dt_s'])])
        self.assertLessEqual(len(after.calls), 10)

    def test_locator_rejects_invalid_parent_and_enlarged_probe_budget(self):
        for parent in (None, '', 'a'*63, 'A'*64, 'z'*64, 123):
            probe = ForcedParticle()
            with self.subTest(parent=parent), self.assertRaises(ValueError):
                current.locate_activation(probe, .001, .6, parent_sha256=parent, time_tolerance=1e-7)
            self.assertFalse(probe.calls)
        probe = ForcedParticle()
        with self.assertRaises(ValueError):
            current.locate_activation(probe, .001, .6, parent_sha256='a'*64,
                                      time_tolerance=1e-7, maximum_solves=11)
        self.assertFalse(probe.calls)

    def test_same_parent_probe_guess_rejects_changed_parent_or_event(self):
        free = np.array([[True, False, True], [False, True, False]])
        state = dict(free=free, identity='fixed-parent')
        prescribed = np.full((2, 3), 2e-6)
        trial = types.SimpleNamespace(pa=object(), state=state, incoming_pid_pairs=((1059,1085),),
            x0=np.zeros((2,3)), p=dict(dx=.001), prescribed=lambda dt: prescribed*dt)
        prior = dict(trial=types.SimpleNamespace(pa=object(), state=copy.deepcopy(state),
            incoming_pid_pairs=((1059,1085),), motion=types.SimpleNamespace(dt_s=.001)),
            fresh=dict(displacement=np.full((2,3), 1e-7)))
        with patch.object(current.cpu, 'parent_fingerprint', side_effect=lambda pa,s: s['identity']):
            guess = current.event_probe_guess(trial, .0005, prior)
            np.testing.assert_array_equal(guess[free], np.full(np.count_nonzero(free), 5e-8))
            np.testing.assert_array_equal(guess[~free], (.0005*prescribed)[~free])
            self.assertFalse(np.shares_memory(guess, prior['fresh']['displacement']))
            prior['trial'].state['identity'] = 'other-parent'
            with self.assertRaisesRegex(ValueError, 'immutable parent'):
                current.event_probe_guess(trial, .0005, prior)
            prior['trial'].state['identity'] = 'fixed-parent'
            prior['trial'].incoming_pid_pairs = ((1059, 1086),)
            with self.assertRaisesRegex(ValueError, 'event identity'):
                current.event_probe_guess(trial, .0005, prior)


if __name__ == '__main__':
    unittest.main(verbosity=2)

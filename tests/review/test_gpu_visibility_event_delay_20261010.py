"""Preserve frozen physics and provide the locator its existing delay gate."""
from pathlib import Path
import ast
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_visibility_event_prism_20261010 as previous
import run_gpu_visibility_event_delay_20261010 as current


class FrozenPhysics(unittest.TestCase):
    def test_entire_driver_changes_only_locator_registration_and_delay_evidence(self):
        before = ast.parse(Path(previous.__file__).read_text(encoding='utf-8'))
        after = ast.parse(Path(current.__file__).read_text(encoding='utf-8'))
        imports = [n for n in after.body if isinstance(n, ast.ImportFrom)
                   and n.module == 'afsd_event_locator_delay_20261010']
        self.assertEqual(len(imports), 1)
        imports[0].module = 'afsd_event_locator_secant_20261010'
        after.body[:] = [n for n in after.body if not
                         (isinstance(n, ast.FunctionDef) and n.name == 'probe_prefix_is_ready')]
        callback = next(n for n in ast.walk(after)
                        if isinstance(n, ast.FunctionDef) and n.name == 'solve_probe')
        ret = next(n for n in callback.body if isinstance(n, ast.Return))
        extras = [k for k in ret.value.keywords if k.arg in
                  ('activation_delay_s', 'activation_prefix_certified')]
        self.assertEqual([(k.arg, ast.unparse(k.value)) for k in extras],
                         [('activation_delay_s', 'activation_delay'),
                          ('activation_prefix_certified', 'prefix_ready')])
        for extra in extras:
            ret.value.keywords.remove(extra)
        start = next(i for i, n in enumerate(callback.body) if isinstance(n, ast.Assign)
                     and ast.unparse(n.targets[0]) == 'activation_delay')
        expected = ast.parse('''activation_delay = None if blocked else probe_dt*(1.-changes[0]['start_fraction'])
prefix_ready = False
if activation_delay is not None and activation_delay <= tolerance:
    prefix_end = 1.-tolerance/probe_dt
    require(0. < prefix_end < 1., 'bounded activation probe requires a positive blocked prefix')
    prefix = audit_blocked_prefix(trial.x0, fresh['position'], i, j,
        trial.motion.pose(trial.motion.t0), trial.motion.pose(trial.motion.t1),
        trial.p, pa.particle_id, prefix_end)
    core.write(probe['directory']/'probe_incoming_blocked_prefix.json', prefix)
    prefix_ready = probe_prefix_is_ready(prefix)
''')
        self.assertEqual([ast.dump(n) for n in callback.body[start:start+3]],
                         [ast.dump(n) for n in expected.body])
        del callback.body[start:start+3]
        registrations = {'LOCATOR_PROTOCOL', 'CONTRACT', 'CHECKS'}
        def normalized(tree):
            tree.body[:] = [n for n in tree.body if not (isinstance(n, ast.Assign)
                and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)
                and n.targets[0].id in registrations)]
            return ast.dump(tree)
        self.assertEqual(normalized(before), normalized(after))

    def test_solver_geometry_and_independent_final_checks_unchanged(self):
        for name in ('model', 'solver', 'gpu', 'cpu', 'core', 'post', 'event', 'batch', 'old',
                     'shared', 'audit_blocked_prefix', 'path_proof'):
            self.assertIs(getattr(current, name), getattr(previous, name))
        for name in ('PROTOCOL', 'BLOCKED_PREFIX_PROTOCOL', 'VISIBLE_PATH_PROTOCOL',
                     'MAXIMUM_PHYSICAL_SOLVES', 'MAXIMUM_EVENTS_PER_BRANCH',
                     'EVENT_TIME_TOLERANCE_S', 'EVENT_RELATIVE_TIME_TOLERANCE',
                     'EVENT_POSITION_TOLERANCE_M'):
            self.assertEqual(getattr(current, name), getattr(previous, name))
        self.assertIsNot(current.locate_activation, previous.locate_activation)

    def test_delay_uses_earliest_router_endpoint_and_registered_paths(self):
        tree = ast.parse(Path(current.__file__).read_text(encoding='utf-8'))
        callback = next(n for n in ast.walk(tree)
                        if isinstance(n, ast.FunctionDef) and n.name == 'solve_probe')
        value = next(n.value for n in callback.body if isinstance(n, ast.Assign)
                     and ast.unparse(n.targets[0]) == 'activation_delay')
        code = compile(ast.Expression(value), '<delay-evidence>', 'eval')
        self.assertIsNone(eval(code, {}, dict(blocked=True, probe_dt=1., changes=[])))
        dt = 0.0002494350930402596
        start, end = .9997591115534306, .9997591124847531
        delay = eval(code, {}, dict(blocked=False, probe_dt=dt,
                                   changes=[dict(start_fraction=start, end_fraction=end)]))
        self.assertAlmostEqual(delay, 60.08603208237255e-9, places=20)
        self.assertGreater(delay, dt*(1.-end))
        self.assertGreater(delay, 50e-9)
        self.assertEqual(current.CONTRACT, ROOT/'docs/GPU_VISIBILITY_EVENT_DELAY_20261010.md')
        self.assertEqual(current.CHECKS, current.BASE/'event_delay_components_20261010_v1.json')
        self.assertEqual(current.LOCATOR_PROTOCOL, 'REAL_ENDPOINT_AND_PATH_DELAY_SECANT_20261010_V1')

    def test_prefix_gate_only_retries_pure_original_budget_exhaustion(self):
        report = dict(status='UNRESOLVED', failure_kind='MAXIMUM_ANALYTIC_BATCH_QUERIES',
                      analytic_batch_queries=256, maximum_queries=256, visible_samples=[],
                      degenerate_edges=[], numerical_domain_failures=[],
                      unresolved_intervals=[dict(reason='MAXIMUM_ANALYTIC_BATCH_QUERIES')])
        self.assertFalse(current.probe_prefix_is_ready(report))
        self.assertTrue(current.probe_prefix_is_ready(dict(report,
            status='CERTIFIED_BLOCKED_PREFIX', failure_kind=None, analytic_batch_queries=200,
            unresolved_intervals=[])))
        for change in (dict(status='OTHER'), dict(failure_kind='MAXIMUM_SUBDIVISION_DEPTH'),
                       dict(analytic_batch_queries=255), dict(maximum_queries=512),
                       dict(visible_samples=[1]), dict(degenerate_edges=[1]),
                       dict(unresolved_intervals=[]),
                       dict(unresolved_intervals=[dict(reason='NONFINITE_BODY_PATH')]),
                       dict(numerical_domain_failures=[1]), dict(error='bad geometry')):
            with self.subTest(change=change), self.assertRaises(ValueError):
                current.probe_prefix_is_ready(dict(report, **change))
        with self.assertRaises(ValueError):
            current.probe_prefix_is_ready(dict(report, status='CERTIFIED_BLOCKED_PREFIX'))


if __name__ == '__main__':
    unittest.main()

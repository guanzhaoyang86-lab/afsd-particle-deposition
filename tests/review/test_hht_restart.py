"""Bounded restart receipt/metadata tests; never advances the HHT trajectory."""
import copy
from contextlib import redirect_stderr
import io
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import restart_hht_probe as tested


def receipt_fixture():
    sources = {'src/afsd_backend/pysph/hht_coupled.py': 'a'*64}
    sources.update({name: 'b'*64 for name in tested.EXTRA_SOURCES})
    sweep = dict(classification=tested.driver.CLASS, status=tested.PASS_SWEEP,
        mode='sweep', production_authorized=False, afsd_progress=False, new_afsd_steps=0,
        receipts=dict(preflight_sha256='c'*64), source_sha256=sources.copy(),
        comparison=dict(passed=True), completed_grids=[dict(steps=n, time_s=2e-6,
                                                          wall_seconds=.01*n) for n in (8, 16, 32)])
    preflight = dict(status=tested.PASS_PREFLIGHT, source_sha256=sources.copy())
    return sweep, preflight, sources


def published_history():
    result, artifacts = [], {}
    for n in range(1, 33):
        relative = f'stage_{n:03d}/accepted.npz'
        result.append(dict(checkpoint=relative, checkpoint_sha256='d'*64,
                           state_hash='e'*64, time_s=n*(2e-6/32.)))
        artifacts['grid_032/'+relative] = 'd'*64
        artifacts[f'grid_032/stage_{n:03d}/integrals.json'] = 'f'*64
    return result, artifacts


def integrals():
    return {name: np.zeros(shape) if shape else 0.
            for name, shape in tested.INTEGRAL_SHAPES.items()}


class Property:
    def __init__(self, value):
        self.value = np.asarray(value)
    def get_npy_array(self):
        return self.value


class RestartGates(unittest.TestCase):
    def test_01_parser_has_fixed_halfway_scope(self):
        args = tested.parser().parse_args(['--sweep', 'sweep.json', '--preflight', 'preflight.json',
                                          '--output', 'fresh'])
        self.assertFalse(args.worker)
        self.assertEqual(args.sweep, Path('sweep.json'))
        self.assertFalse(hasattr(args, 'checkpoint'))
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            tested.parser().parse_args(['--sweep', 'sweep.json'])
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            tested.parser().parse_args(['--sweep', 's', '--preflight', 'p', '--output', 'o',
                                       '--checkpoint', 'unpublished.npz'])

    def test_02_receipt_success_and_fixed_cost_ceiling(self):
        sweep, preflight, sources = receipt_fixture()
        unit, forecast = tested.validate_receipt_payloads(sweep, preflight, sources, 'c'*64)
        self.assertEqual(unit, .01)
        self.assertAlmostEqual(forecast, 5.32, delta=1e-14)
        sweep['completed_grids'][-1]['wall_seconds'] = 45.
        self.assertEqual(tested.validate_receipt_payloads(sweep, preflight, sources, 'c'*64)[1], 50.)
        sweep['completed_grids'][-1]['wall_seconds'] = 45.0001
        with self.assertRaisesRegex(ValueError, 'COST_EXCEEDS'):
            tested.validate_receipt_payloads(sweep, preflight, sources, 'c'*64)

    def test_03_failed_stale_partial_or_nonlocal_receipts_reject(self):
        changes = [lambda s, p: s.update(status='FAILED_LOCAL_BENCHMARK'),
                   lambda s, p: s.update(classification='AFSD'),
                   lambda s, p: s.update(mode='pilot'),
                   lambda s, p: s.update(production_authorized=True),
                   lambda s, p: s.update(afsd_progress=True),
                   lambda s, p: s.update(new_afsd_steps=1),
                   lambda s, p: p.update(status='FAIL_HHT_PREFLIGHT'),
                   lambda s, p: s['receipts'].update(preflight_sha256='wrong'),
                   lambda s, p: s['comparison'].update(passed=False),
                   lambda s, p: s.update(completed_grids=s['completed_grids'][:2]),
                   lambda s, p: s['completed_grids'][-1].update(time_s=1.5e-6),
                   lambda s, p: s['completed_grids'][-1].update(wall_seconds=float('nan')),
                   lambda s, p: s['source_sha256'].pop(tested.EXTRA_SOURCES[0]),
                   lambda s, p: p['source_sha256'].pop(tested.EXTRA_SOURCES[2])]
        for number, change in enumerate(changes):
            with self.subTest(case=number):
                sweep, preflight, sources = receipt_fixture()
                change(sweep, preflight)
                with self.assertRaises(ValueError):
                    tested.validate_receipt_payloads(sweep, preflight, sources, 'c'*64)

    def test_04_only_complete_published_checkpoint_inventory_is_allowed(self):
        history, artifacts = published_history()
        tested.validate_history(history, artifacts, 2e-6/32.)
        cases = [lambda h, a: h.pop(),
                 lambda h, a: h[15].update(checkpoint='private/accepted.npz'),
                 lambda h, a: h[15].update(checkpoint='../stage_016/accepted.npz'),
                 lambda h, a: h[15].update(checkpoint_sha256='changed'),
                 lambda h, a: h[15].update(state_hash=''),
                 lambda h, a: h[15].update(time_s=0.),
                 lambda h, a: a.pop('grid_032/stage_017/integrals.json')]
        for number, change in enumerate(cases):
            with self.subTest(case=number):
                modified, manifest = copy.deepcopy(history), dict(artifacts)
                change(modified, manifest)
                with self.assertRaises(ValueError):
                    tested.validate_history(modified, manifest, 2e-6/32.)

    def test_05_manifest_escape_missing_and_byte_tamper_reject(self):
        with tempfile.TemporaryDirectory(prefix='hht_restart_receipt_') as temporary:
            base = Path(temporary)
            artifact = base/'receipt.json'; artifact.write_text('{}', encoding='utf-8')
            digest = tested.driver.sha(artifact)
            tested.verify_manifest(base, {'receipt.json': digest}, 'test artifact')
            with self.assertRaisesRegex(ValueError, 'escaping'):
                tested.checked_path(base, '../external.npz')
            with self.assertRaisesRegex(ValueError, 'escaping'):
                tested.checked_path(base, str(artifact.resolve()))
            with self.assertRaisesRegex(ValueError, 'drift'):
                tested.verify_manifest(base, {'missing.npz': digest}, 'test artifact')
            artifact.write_text('{"changed":true}', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'drift'):
                tested.verify_manifest(base, {'receipt.json': digest}, 'test artifact')

    def test_06_complete_finite_integrals_and_bitwise_difference(self):
        ledger = integrals()
        self.assertTrue(tested.integral_comparison(ledger, tested.driver.clean(ledger))['bitwise_equal'])
        changed = copy.deepcopy(ledger); changed['friction_raw'][0] = 1e-30
        result = tested.integral_comparison(ledger, changed)
        self.assertFalse(result['bitwise_equal'])
        self.assertEqual(result['maximum_absolute_differences']['friction_raw'], 1e-30)
        for change in (lambda x: x.pop('root_residual_work_J'),
                       lambda x: x.update(impulse=np.zeros(6)),
                       lambda x: x.update(root_residual_work_J=float('nan'))):
            invalid = copy.deepcopy(ledger); change(invalid)
            with self.assertRaises(ValueError):
                tested.parse_integrals(invalid)

    def test_07_full_particle_semantic_and_ledger_comparison_is_conjunctive(self):
        first = SimpleNamespace(properties={'x': Property([0., .01]), 'm': Property([1., 2.])})
        second = copy.deepcopy(first)
        state = dict(state_hash='a'*64, acceleration_m_s2=np.zeros((2, 3)), force_N=np.zeros((2, 3)))
        copied = copy.deepcopy(state)
        ledger = integrals()
        result = tested.compare_candidate(first, state, ledger, second, copied, ledger)
        self.assertTrue(result['passed'])
        second.properties['x'].value[0] = 1e-30
        result = tested.compare_candidate(first, state, ledger, second, copied, ledger)
        self.assertFalse(result['passed'])
        self.assertEqual(result['maximum_property_differences']['x'], 1e-30)
        second = copy.deepcopy(first); copied['state_hash'] = 'b'*64
        self.assertFalse(tested.compare_candidate(first, state, ledger, second, copied, ledger)['passed'])
        copied = copy.deepcopy(state); wrong = copy.deepcopy(ledger); wrong['wall_work'][1] = 1e-30
        self.assertFalse(tested.compare_candidate(first, state, ledger, second, copied, wrong)['passed'])


if __name__ == '__main__':
    unittest.main(verbosity=2)

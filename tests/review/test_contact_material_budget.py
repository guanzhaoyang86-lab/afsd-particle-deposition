"""Cost gate arithmetic, before the single bounded root."""
from pathlib import Path
import json
import subprocess
import sys
import time
import unittest
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from run_contact_material_root import evaluation_budget
from probe_updated_stage import sha


class BudgetTests(unittest.TestCase):
    def test_measured_cost_reduces_evaluations_without_extending_time(self):
        b = evaluation_budget(9.064526899999919)
        self.assertEqual(b['maximum_evaluations'], 30)
        self.assertAlmostEqual(b['forecast_wall_seconds'], 437.9037104999964)
        self.assertLessEqual(b['forecast_wall_seconds'], b['launch_limit_s'])
        self.assertEqual(b['external_limit_s'], 600)

    def test_limits_and_nonfinite(self):
        self.assertEqual(evaluation_budget(1.)['maximum_evaluations'], 32)
        for bad in (0., -1., float('nan'), float('inf'), 150.):
            with self.assertRaises(ValueError): evaluation_budget(bad)


if __name__ == '__main__':
    if '--child' in sys.argv:
        sys.argv.remove('--child')
        unittest.main(verbosity=2)
    else:
        logs = ROOT/'logs/contact_material_search_20261004'
        log = logs/'budget_test.log'; receipt = logs/'budget_test.json'
        if log.exists() or receipt.exists(): raise ValueError('fresh budget test evidence required')
        start = time.perf_counter()
        command = [sys.executable, '-B', str(Path(__file__).resolve()), '--child']
        with log.open('wb') as f:
            r = subprocess.run(command, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, timeout=60)
        sources = ['scripts/run_contact_material_root.py', 'tests/review/test_contact_material_budget.py',
                   'docs/CONTACT_MATERIAL_BUDGET_ADDENDUM.md']
        row = dict(exit_code=r.returncode, command=command, wall_seconds=time.perf_counter()-start,
            source_sha256={s: sha(ROOT/s) for s in sources}, log_sha256=sha(log))
        receipt.write_text(json.dumps(row, indent=2), encoding='utf-8')
        print(json.dumps(row)); raise SystemExit(r.returncode)

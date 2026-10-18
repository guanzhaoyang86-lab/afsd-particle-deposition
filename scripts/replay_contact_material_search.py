"""Replay the frozen numerical package into fresh output directories.

Only output-directory bindings change. Each invocation runs one bounded gate;
it never starts a suite of physical experiments or alters physical settings.
The time-refinement branch has its own archived driver and is not dispatched.
"""
from pathlib import Path
import argparse
import importlib
import json
import shutil
import subprocess
import sys
import time
ROOT = Path(__file__).resolve().parents[1]
import run_contact_material_search as package


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage', choices=['tests', 'response', 'cost', 'root'])
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = ap.parse_args(); output = args.output.resolve()
    if (not output.is_relative_to(ROOT/'outputs') or output == ROOT/'outputs' or
            output == package.BASE or package.BASE.is_relative_to(output) or
            output.is_relative_to(package.BASE)):
        raise ValueError('fresh distinct output subdirectory required; archived package is protected')
    logs = output/'execution_logs'; label = 'tests_review1' if args.stage == 'tests' else args.stage
    if args.worker:
        # The frozen functions resolve these storage paths at invocation time.
        # Import root driver only after binding its imported storage constants.
        package.BASE = output; package.LOGS = logs
        if args.stage == 'root':
            module = importlib.import_module('run_contact_material_root')
            original_logs = ROOT/'logs/contact_material_search_20261004'
            for name in ('budget_test.json', 'budget_test.log'):
                source = original_logs/name; target = logs/name
                if target.exists():
                    if package.sha(target) != package.sha(source): raise ValueError('budget test evidence drift')
                else:
                    shutil.copyfile(source, target)
            return module.worker()
        return package.worker(args.stage, label)
    logs.mkdir(parents=True, exist_ok=True)
    log = logs/(label+'.log'); receipt = logs/(label+'_supervision.json')
    if log.exists() or receipt.exists() or (output/label).exists():
        raise ValueError('stage already exists; use another output directory')
    command = [sys.executable, '-B', str(Path(__file__).resolve()), args.stage,
               '--output', str(output), '--worker']
    timeout = dict(tests=240, response=240, cost=90, root=600)[args.stage]
    start = time.perf_counter()
    with log.open('wb') as stream:
        try:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = 124
    result = dict(command=command, exit_code=code, timeout_s=timeout,
        wall_seconds=time.perf_counter()-start, log_sha256=package.sha(log),
        storage_dispatcher_sha256=package.sha(Path(__file__).resolve()),
        scope='OUTPUT_PATH_BINDINGS_ONLY; FROZEN_PHYSICAL_AND_NUMERICAL_IMPLEMENTATION')
    receipt.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result)); return code


if __name__ == '__main__': raise SystemExit(main())

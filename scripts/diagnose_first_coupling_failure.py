"""Bounded legacy replay: diagnostics must not change one accepted bit."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import threading
import os
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.config_adapter import load_settings
from afsd_backend.pysph.particles import create_geometry
from afsd_backend.pysph.solver import SolidSolver
from afsd_backend.pysph.kick_diagnostics import KickRecorder, write_json, clean
from afsd_backend.pysph.output import digest

def signature(solver):
    h = hashlib.sha256()
    for name, value in sorted(solver.pa.properties.items()):
        h.update(name.encode()); h.update(value.get_npy_array().tobytes())
    h.update(json.dumps(solver.state(), sort_keys=True).encode())
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', default='configs/debug/starved.json')
    parser.add_argument('--output', default='outputs/coupling_20260930/phase1')
    parser.add_argument('--steps', type=int, default=40)
    parser.add_argument('--wall-seconds', type=float, default=600)
    args = parser.parse_args()
    if not 0 < args.wall_seconds <= 600 or not 0 < args.steps <= 40:
        raise ValueError('diagnostic cap: 40 kicks, 600 seconds')
    watchdog = threading.Timer(args.wall_seconds, lambda: os._exit(124))
    watchdog.daemon = True
    watchdog.start()
    config, p, consumed = load_settings(args.case)
    hashes = {str(f.relative_to(ROOT)): digest(f) for f in sorted((ROOT/'src/afsd_backend').rglob('*.py'))}
    provenance = dict(source_hashes=hashes, runner_sha256=digest(__file__),
        configuration_file_sha256=digest(args.case),
        resolved_configuration_sha256=hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
        backend_commit='1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5',
        configuration=config, consumed=consumed, physics='UNCHANGED_LEGACY_DIAGNOSTIC_ONLY')
    start = time.perf_counter()
    baseline = SolidSolver(create_geometry(p), p)
    expected = []
    legacy_error = None
    for _ in range(args.steps):
        try:
            baseline.advance()
            expected.append(signature(baseline))
        except (RuntimeError, ValueError) as error:
            legacy_error = str(error)
            break
    baseline_seconds = time.perf_counter()-start
    if legacy_error is None:
        raise RuntimeError('no legacy failure found within cap')
    recorder = KickRecorder(args.output, provenance)
    observer = SolidSolver(create_geometry(p), p, diagnostic_recorder=recorder)
    recorder.accept(observer)
    matched = 0
    try:
        for _ in range(args.steps):
            observer.advance()
            if signature(observer) != expected[observer.count-1]:
                raise AssertionError('observer changed accepted state at '+str(observer.count))
            matched += 1
            recorder.accept(observer)
        raise AssertionError('observed run unexpectedly succeeded')
    except (RuntimeError, ValueError) as error:
        recorder.reject(observer, error)
        if str(error) != legacy_error:
            raise AssertionError('different failure: '+str(error))
    elapsed = time.perf_counter()-start
    report = dict(status='PHYSICS_FAILED_DIAGNOSTIC_CAPTURED',
        observer_noninterference='PASS', bitwise_identical_accepted_steps=matched,
        accepted_time_s=observer.time, failed_candidate_step=observer.count+1,
        failure=legacy_error, baseline_wall_seconds=baseline_seconds,
        total_wall_seconds=elapsed, particle_count=observer.pa.get_number_of_particles(),
        initial_dt_s=SolidSolver(create_geometry(p), p).stable_dt(),
        cap_seconds=args.wall_seconds, source_hashes=hashes,
        first_negative=recorder.first_negative)
    write_json(Path(args.output)/'report.json', report)
    print(json.dumps(clean(report), indent=2))
    watchdog.cancel()
    return 0

if __name__ == '__main__':
    raise SystemExit(main())


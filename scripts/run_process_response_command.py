"""Hard60s supervision and fresh receipts for the response work package."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
ALLOWED_COMMANDS = ('scripts/check_process_response.py', 'scripts/probe_hht_response.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.label.replace('_', '').isalnum():
        raise ValueError('simple fresh label required')
    cmd = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not cmd or cmd[0] not in ALLOWED_COMMANDS:
        raise ValueError('only registered component gate or PRIVATE saved diagnostic command')
    directory = ROOT/'logs/process_response_20261003'
    directory.mkdir(exist_ok=True)
    receipt = directory/(args.label+'_supervision.json')
    log = directory/(args.label+'.log')
    if receipt.exists() or log.exists():
        raise FileExistsError('never overwrite command evidence')
    command = [sys.executable, '-B']+cmd
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    started = time.perf_counter()
    timed_out = False
    with log.open('x', encoding='utf-8') as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=env,
                                   stdout=stream, stderr=subprocess.STDOUT)
        try:
            code = process.wait(timeout=60.)
        except subprocess.TimeoutExpired:
            timed_out = True
            # The preflight has its own55s shared child-test deadline. Kill
            # any remaining Windows descendants as well as the named command
            # if an unexpected stall reaches the external deadline.
            if os.name == 'nt':
                try:
                    subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                   timeout=5., check=False)
                except (OSError, subprocess.TimeoutExpired):
                    pass
            if process.poll() is None:
                process.kill()
            process.wait()
            stream.write('\nEXTERNAL60S_TIMEOUT\n')
            code = 124
    elapsed = time.perf_counter()-started
    data = dict(command=command, return_code=code, wall_seconds=elapsed,
        external_limit_seconds=60., timed_out=timed_out,
        classification='BOUNDED_RESPONSE_PACKAGE_COMMAND_NOT_PROCESS_AUTHORIZATION',
        production_authorized=False,
        log=log.name, log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
        supervisor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with receipt.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(data, indent=2, allow_nan=False)+'\n')
    print(json.dumps(data), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())

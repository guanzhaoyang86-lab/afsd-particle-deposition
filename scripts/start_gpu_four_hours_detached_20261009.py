"""Hidden, explicit Windows Job breakaway for a qualified four-hour campaign.

The detached child records its real Job policy before invoking the fixed
supervisor. A remaining kill-on-close policy fails closed; no fallback launch
or change to system Job policy is attempted. This module imports no solver.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import runpy
import subprocess
import sys
import start_gpu_resume_detached_20261008 as original

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/gpu_acceleration_20261007'
SUPERVISOR = ROOT/'scripts/supervise_gpu_four_hours_20261009.py'
CREATION_FLAGS = 16777736
job_policy, write_new = original.job_policy, original.write_new


def inside(value):
    path = (ROOT/Path(value)).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError('project-local path required')
    return path


def arguments(args):
    output, parent, runner = (inside(getattr(args, name)) for name in ('output', 'parent', 'runner'))
    if not output.is_relative_to(BASE) or output == BASE or output.exists():
        raise ValueError('new project-local GPU campaign required')
    if not parent.is_relative_to(BASE) or not parent.is_dir() or not runner.is_file():
        raise ValueError('existing project-local parent and runner required')
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 14400:
        raise ValueError('at most four authorized hours required')
    if (not math.isfinite(args.minimum_dt) or not 0 < args.minimum_dt <= 1e-7
            or not 0 <= args.max_halvings <= 8):
        raise ValueError('registered temporal refinement bounds required')
    if args.initial_dt is not None and (not math.isfinite(args.initial_dt)
            or not args.minimum_dt <= args.initial_dt <= .005):
        raise ValueError('finite registered initial dt required')
    if not SUPERVISOR.is_file():
        raise ValueError('fixed four-hour supervisor missing')
    values = ['--output', output.relative_to(ROOT).as_posix(),
              '--parent', parent.relative_to(ROOT).as_posix(),
              '--runner', runner.relative_to(ROOT).as_posix(),
              '--wall-seconds', repr(args.wall_seconds),
              '--minimum-dt', repr(args.minimum_dt),
              '--max-halvings', str(args.max_halvings)]
    if args.runner_qualification is not None:
        qualification = inside(args.runner_qualification)
        if not qualification.is_file():
            raise ValueError('existing runner qualification required')
        values += ['--runner-qualification', qualification.relative_to(ROOT).as_posix()]
    if args.initial_dt is not None:
        values += ['--initial-dt', repr(args.initial_dt)]
    if args.preflight_only:
        values += ['--preflight-only']
    return output, values


def source_hashes():
    paths = (Path(__file__).resolve(), Path(original.__file__).resolve(), SUPERVISOR)
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    parser.add_argument('--parent')
    parser.add_argument('--runner')
    parser.add_argument('--runner-qualification')
    parser.add_argument('--wall-seconds', type=float, default=14400.)
    parser.add_argument('--initial-dt', type=float)
    parser.add_argument('--minimum-dt', type=float, default=1e-7)
    parser.add_argument('--max-halvings', type=int, default=8)
    parser.add_argument('--preflight-only', action='store_true')
    parser.add_argument('--detached-child', action='store_true')
    parser.add_argument('--inspect-job', action='store_true')
    args = parser.parse_args(argv)
    if args.inspect_job:
        print(json.dumps(job_policy())); return 0
    for name in ('output', 'parent', 'runner'):
        if getattr(args, name) is None:
            parser.error('--'+name+' required')
    output, forwarded = arguments(args)
    bootstrap = output.with_name(output.name+'.bootstrap.json')
    if args.detached_child:
        policy = job_policy()
        write_new(bootstrap, dict(status='PASS_DETACHED_JOB_POLICY' if not policy['kill_on_job_close']
            else 'FAIL_KILL_ON_JOB_CLOSE_STILL_PRESENT', **policy,
            supervisor=SUPERVISOR.relative_to(ROOT).as_posix(), supervisor_arguments=forwarded,
            wall_limit_s=args.wall_seconds, source_sha256=source_hashes()))
        if policy['kill_on_job_close']:
            raise RuntimeError('detached supervisor still inherits kill-on-job-close policy')
        sys.argv = [str(SUPERVISOR), *forwarded]
        sys.path.insert(0, str(ROOT/'scripts'))
        runpy.run_path(str(SUPERVISOR), run_name='__main__')
        return 0
    receipt = output.with_name(output.name+'.launch.json')
    log = output.with_name(output.name+'.stdout.log')
    if any(p.exists() for p in (receipt, log, bootstrap)):
        raise ValueError('new detached launch records required')
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--detached-child', *forwarded]
    startup = subprocess.STARTUPINFO()
    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW; startup.wShowWindow = 0
    flags = subprocess.CREATE_BREAKAWAY_FROM_JOB | subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
    if flags != CREATION_FLAGS:
        raise RuntimeError('registered Windows creation flags differ')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    parent_policy = job_policy(); sources = source_hashes()
    with log.open('xb', buffering=0) as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
            stdout=stream, stderr=subprocess.STDOUT, startupinfo=startup,
            creationflags=flags, close_fds=True)
    evidence = dict(status='STARTED_DETACHED_SUPERVISOR_PENDING_BOOTSTRAP_CHECK', pid=process.pid,
        command=command, creationflags=flags, parent_job_policy=parent_policy,
        wall_limit_s=args.wall_seconds, output=output.relative_to(ROOT).as_posix(),
        supervisor=SUPERVISOR.relative_to(ROOT).as_posix(), supervisor_arguments=forwarded,
        source_sha256=sources)
    write_new(receipt, evidence)
    print(json.dumps(evidence)); return 0


if __name__ == '__main__':
    raise SystemExit(main())

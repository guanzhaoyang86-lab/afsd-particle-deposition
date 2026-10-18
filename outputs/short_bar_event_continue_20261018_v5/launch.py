"""Hidden CPU campaign launch with a fresh, explicitly supplied wall-time budget.

Importing this module and --inspect-job never write launch records or start work.
An actual parent launch freezes one configuration and all descendants inherit it.
"""
from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
import hashlib
import json
import math
import os
import runpy
import subprocess
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from start_gpu_resume_detached_20261008 import job_policy, write_new
import start_gpu_resume_detached_20261008 as original

CAMPAIGN = BASE / 'controller.py'
CONFIG = BASE / 'run.configuration.json'
CONFIGURATION_SCHEMA = 'CPU_LOCAL_INTERVAL_RUN_CONFIGURATION_V1'
BOOTSTRAP = BASE / 'run.bootstrap.json'
RECEIPT = BASE / 'run.launch.json'
LOG = BASE / 'run.stdout.log'
EXIT = BASE / 'run.exit.json'


def now():
    return datetime.now(timezone.utc)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    return {p.relative_to(ROOT).as_posix(): digest(p) for p in
            (CONFIG, Path(__file__).resolve(), CAMPAIGN, Path(original.__file__).resolve())}


def duration(value):
    """No arbitrary long-run cap; the datetime range remains representable."""
    if isinstance(value, bool):
        raise ValueError('--wall-seconds must be finite and at least 650')
    try:
        seconds = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError('--wall-seconds must be finite and at least 650') from error
    if not math.isfinite(seconds) or seconds < 650.:
        raise ValueError('--wall-seconds must be finite and at least 650')
    return seconds


def validate_configuration(configuration):
    if not isinstance(configuration, dict) or set(configuration) != {
            'schema', 'started_utc', 'deadline_utc', 'wall_seconds'}:
        raise ValueError('exact run configuration schema required')
    if configuration['schema'] != CONFIGURATION_SCHEMA:
        raise ValueError('unsupported run configuration schema')
    if type(configuration['wall_seconds']) not in (int, float):
        raise ValueError('configuration wall_seconds must be numeric')
    seconds = duration(configuration['wall_seconds'])
    try:
        started = datetime.fromisoformat(configuration['started_utc'])
        deadline = datetime.fromisoformat(configuration['deadline_utc'])
    except (TypeError, ValueError) as error:
        raise ValueError('UTC ISO configuration times required') from error
    if started.tzinfo is None or deadline.tzinfo is None or \
            started.utcoffset() != timedelta(0) or deadline.utcoffset() != timedelta(0):
        raise ValueError('configuration timestamps must use UTC')
    if not math.isclose((deadline-started).total_seconds(), seconds, rel_tol=0., abs_tol=1e-6):
        raise ValueError('configuration deadline must equal its original start plus duration')
    return deadline


def load_configuration(expected_sha):
    if not isinstance(expected_sha, str) or expected_sha != digest(CONFIG):
        raise RuntimeError('Frozen run configuration SHA mismatch')
    configuration = json.loads(CONFIG.read_text(encoding='utf-8'))
    deadline = validate_configuration(configuration)
    return configuration, deadline


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--detached-child', action='store_true')
    parser.add_argument('--expected-campaign-sha')
    parser.add_argument('--expected-configuration-sha')
    # Parse as text so --inspect-job precedes even numeric duration validation.
    parser.add_argument('--wall-seconds')
    parser.add_argument('--inspect-job', action='store_true')
    args = parser.parse_args(argv)
    if args.inspect_job:
        print(json.dumps(job_policy()))
        return 0
    if args.detached_child:
        if args.wall_seconds is not None:
            parser.error('detached child must inherit its frozen configuration, not a new duration')
        if not args.expected_campaign_sha or not args.expected_configuration_sha:
            parser.error('detached child requires both expected source and configuration hashes')
    else:
        if args.wall_seconds is None:
            parser.error('--wall-seconds is required for an actual launch')
        if args.expected_campaign_sha is not None or args.expected_configuration_sha is not None:
            parser.error('expected hashes are internal detached-child arguments')
        try:
            wall_seconds = duration(args.wall_seconds)
        except ValueError as error:
            parser.error(str(error))
    if os.name != 'nt' or not CAMPAIGN.is_file():
        raise RuntimeError('Registered Windows campaign must exist')
    if Path(sys.executable).resolve() != (ROOT / '.envs/pysph/python.exe').resolve():
        raise RuntimeError('Use the registered PySPH Python interpreter')
    if (BASE / 'run').exists():
        raise RuntimeError('Fresh run directory required; refusing duplicate launch')
    if args.detached_child:
        if BOOTSTRAP.exists() or EXIT.exists():
            raise RuntimeError('Fresh detached-child records required')
        configuration, deadline = load_configuration(args.expected_configuration_sha)
        if deadline <= now():
            raise RuntimeError('Frozen run configuration deadline already exhausted')
        policy = job_policy()
        matched = args.expected_campaign_sha == digest(CAMPAIGN)
        qualified = matched and policy['kill_on_job_close'] is False
        write_new(BOOTSTRAP, dict(status='PASS_DETACHED_JOB_POLICY_AND_SOURCE' if qualified else
            'FAIL_DETACHED_JOB_POLICY_OR_SOURCE', created_utc=now().isoformat(),
            **policy, campaign_source_matches=matched, source_sha256=sources(),
            configuration_sha256=args.expected_configuration_sha,
            authorized_started_utc=configuration['started_utc'],
            authorized_deadline_utc=configuration['deadline_utc'], physics_started=False))
        if not qualified:
            raise RuntimeError('Detached policy or frozen campaign source mismatch')
        code = 2
        error = None
        previous_argv = sys.argv
        try:
            sys.argv = [str(CAMPAIGN), '--run', '--configuration', str(CONFIG),
                        '--expected-configuration-sha', args.expected_configuration_sha]
            runpy.run_path(str(CAMPAIGN), run_name='__main__')
            code = 0
        except SystemExit as exc:
            code = 0 if exc.code is None else exc.code if isinstance(exc.code, int) else 2
        except BaseException as exc:
            error = repr(exc)
            raise
        finally:
            sys.argv = previous_argv
            write_new(EXIT, dict(exit_code=code, error=error, finished_utc=now().isoformat(),
                                pid=os.getpid(), configuration_sha256=args.expected_configuration_sha,
                                authorized_started_utc=configuration['started_utc'],
                                authorized_deadline_utc=configuration['deadline_utc']))
        return code
    if any(path.exists() for path in (CONFIG, BOOTSTRAP, RECEIPT, LOG, EXIT)):
        raise RuntimeError('New launch records required; refusing duplicate launch')
    policy = job_policy()
    started = now()
    try:
        deadline = started + timedelta(seconds=wall_seconds)
    except (OverflowError, ValueError) as error:
        raise ValueError('requested duration is outside the representable datetime range') from error
    configuration = dict(schema=CONFIGURATION_SCHEMA, started_utc=started.isoformat(),
                         deadline_utc=deadline.isoformat(), wall_seconds=wall_seconds)
    validate_configuration(configuration)
    # Exclusive creation also arbitrates concurrent parent launch attempts.
    write_new(CONFIG, configuration)
    configuration_sha = digest(CONFIG)
    source_sha256 = sources()
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--detached-child',
               '--expected-campaign-sha', digest(CAMPAIGN),
               '--expected-configuration-sha', configuration_sha]
    startup = subprocess.STARTUPINFO()
    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startup.wShowWindow = 0
    flags = subprocess.CREATE_BREAKAWAY_FROM_JOB | subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    try:
        with LOG.open('xb', buffering=0) as stream:
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                stdout=stream, stderr=subprocess.STDOUT, startupinfo=startup,
                creationflags=flags, close_fds=True)
    except BaseException as error:
        write_new(RECEIPT, dict(status='FAIL_DETACHED_CAMPAIGN_LAUNCH', error=repr(error),
            created_utc=now().isoformat(), configuration_sha256=configuration_sha,
            authorized_started_utc=configuration['started_utc'],
            authorized_deadline_utc=configuration['deadline_utc'], source_sha256=source_sha256,
            physics_started=False, configuration_retained=True))
        raise
    evidence = dict(status='STARTED_DETACHED_CAMPAIGN_PENDING_BOOTSTRAP', pid=process.pid,
        created_utc=now().isoformat(), command=command, creationflags=flags,
        configuration_sha256=configuration_sha,
        authorized_started_utc=configuration['started_utc'],
        authorized_deadline_utc=configuration['deadline_utc'], parent_job_policy=policy,
        source_sha256=source_sha256, progress=(BASE/'run/progress.json').relative_to(ROOT).as_posix())
    write_new(RECEIPT, evidence)
    print(json.dumps(evidence))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

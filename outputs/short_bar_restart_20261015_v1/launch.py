"""Hidden detached launcher for the explicitly authorized short-bar CPU campaign."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import runpy
import subprocess
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from start_gpu_resume_detached_20261008 import job_policy, write_new
import start_gpu_resume_detached_20261008 as original

CAMPAIGN = BASE / 'campaign.py'
DEADLINE = datetime(2026, 10, 15, 10, 15, 41, tzinfo=timezone.utc)
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
            (Path(__file__).resolve(), CAMPAIGN, Path(original.__file__).resolve())}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--detached-child', action='store_true')
    parser.add_argument('--expected-campaign-sha')
    parser.add_argument('--inspect-job', action='store_true')
    args = parser.parse_args(argv)
    if args.inspect_job:
        print(json.dumps(job_policy()))
        return 0
    if os.name != 'nt' or not CAMPAIGN.is_file():
        raise RuntimeError('Registered Windows campaign must exist')
    if Path(sys.executable).resolve() != (ROOT / '.envs/pysph/python.exe').resolve():
        raise RuntimeError('Use the registered PySPH Python interpreter')
    if (BASE / 'run').exists() or (DEADLINE - now()).total_seconds() < 650:
        raise RuntimeError('New campaign and sufficient original authorized time required')
    if args.detached_child:
        policy = job_policy()
        matched = args.expected_campaign_sha == digest(CAMPAIGN)
        qualified = matched and not policy['kill_on_job_close']
        write_new(BOOTSTRAP, dict(status='PASS_DETACHED_JOB_POLICY_AND_SOURCE' if qualified else
            'FAIL_DETACHED_JOB_POLICY_OR_SOURCE', created_utc=now().isoformat(),
            **policy, campaign_source_matches=matched, source_sha256=sources(),
            authorized_deadline_utc=DEADLINE.isoformat(), physics_started=False))
        if not qualified:
            raise RuntimeError('Detached policy or frozen campaign source mismatch')
        code = 2
        error = None
        try:
            sys.argv = [str(CAMPAIGN), '--run']
            runpy.run_path(str(CAMPAIGN), run_name='__main__')
            code = 0
        except SystemExit as exc:
            code = 0 if exc.code is None else exc.code if isinstance(exc.code, int) else 2
        except BaseException as exc:
            error = repr(exc)
            raise
        finally:
            write_new(EXIT, dict(exit_code=code, error=error, finished_utc=now().isoformat(),
                                pid=os.getpid(), authorized_deadline_utc=DEADLINE.isoformat()))
        return code
    if any(path.exists() for path in (BOOTSTRAP, RECEIPT, LOG, EXIT)):
        raise RuntimeError('New launch records required; refusing duplicate launch')
    policy = job_policy()
    source_sha256 = sources()
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--detached-child',
               '--expected-campaign-sha', digest(CAMPAIGN)]
    startup = subprocess.STARTUPINFO()
    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startup.wShowWindow = 0
    flags = subprocess.CREATE_BREAKAWAY_FROM_JOB | subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    with LOG.open('xb', buffering=0) as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
            stdout=stream, stderr=subprocess.STDOUT, startupinfo=startup,
            creationflags=flags, close_fds=True)
    evidence = dict(status='STARTED_DETACHED_CAMPAIGN_PENDING_BOOTSTRAP', pid=process.pid,
        created_utc=now().isoformat(), command=command, creationflags=flags,
        authorized_deadline_utc=DEADLINE.isoformat(), parent_job_policy=policy,
        source_sha256=source_sha256, progress=(BASE/'run/progress.json').relative_to(ROOT).as_posix())
    write_new(RECEIPT, evidence)
    print(json.dumps(evidence))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

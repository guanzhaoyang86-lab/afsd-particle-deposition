"""Hidden, explicitly permitted Job breakaway for the bounded GPU resume.

The child records its actual Job policy before invoking the tested supervisor.
No elevation, system configuration change, or fallback to a kill-on-close Job.
"""
from pathlib import Path
from ctypes import wintypes
import argparse
import ctypes
import hashlib
import json
import os
import runpy
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/gpu_acceleration_20261007'
SUPERVISOR = ROOT/'scripts/resume_gpu_campaign_20261008.py'


def job_policy():
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    kernel.IsProcessInJob.argtypes = [wintypes.HANDLE, wintypes.HANDLE, ctypes.POINTER(wintypes.BOOL)]
    member = wintypes.BOOL()
    if not kernel.IsProcessInJob(kernel.GetCurrentProcess(), None, ctypes.byref(member)):
        raise ctypes.WinError(ctypes.get_last_error())

    class Basic(ctypes.Structure):
        _fields_ = [('process_time', ctypes.c_int64), ('job_time', ctypes.c_int64),
            ('flags', wintypes.DWORD), ('minimum_working_set', ctypes.c_size_t),
            ('maximum_working_set', ctypes.c_size_t), ('active_process_limit', wintypes.DWORD),
            ('affinity', ctypes.c_size_t), ('priority', wintypes.DWORD), ('scheduling', wintypes.DWORD)]

    class Counters(ctypes.Structure):
        _fields_ = [(name, ctypes.c_uint64) for name in
            ('read_ops', 'write_ops', 'other_ops', 'read_bytes', 'write_bytes', 'other_bytes')]

    class Extended(ctypes.Structure):
        _fields_ = [('basic', Basic), ('io', Counters), ('process_memory', ctypes.c_size_t),
            ('job_memory', ctypes.c_size_t), ('peak_process_memory', ctypes.c_size_t),
            ('peak_job_memory', ctypes.c_size_t)]

    flags = 0
    if member.value:
        info = Extended(); returned = wintypes.DWORD()
        kernel.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int,
            ctypes.c_void_p, wintypes.DWORD, ctypes.POINTER(wintypes.DWORD)]
        if not kernel.QueryInformationJobObject(None, 9, ctypes.byref(info), ctypes.sizeof(info), ctypes.byref(returned)):
            raise ctypes.WinError(ctypes.get_last_error())
        flags = int(info.basic.flags)
    return dict(pid=os.getpid(), in_job=bool(member.value), limit_flags=flags,
        breakaway_allowed=bool(flags & 0x800), kill_on_job_close=bool(flags & 0x2000))


def write_new(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    parser.add_argument('--detached-child', action='store_true')
    parser.add_argument('--inspect-job', action='store_true')
    args = parser.parse_args()
    if args.inspect_job:
        print(json.dumps(job_policy())); return 0
    if not args.output:
        parser.error('--output required')
    output = (ROOT/args.output).resolve()
    if not output.is_relative_to(BASE) or output == BASE or output.exists():
        raise ValueError('new project-local GPU campaign required')
    prefix = output.parent/output.name
    bootstrap = prefix.with_name(prefix.name+'.bootstrap.json')
    if args.detached_child:
        policy = job_policy()
        write_new(bootstrap, dict(status='PASS_DETACHED_JOB_POLICY' if not policy['kill_on_job_close']
            else 'FAIL_KILL_ON_JOB_CLOSE_STILL_PRESENT', **policy))
        if policy['kill_on_job_close']:
            raise RuntimeError('detached supervisor still inherits kill-on-job-close policy')
        sys.argv = [str(SUPERVISOR), '--output', output.relative_to(ROOT).as_posix()]
        sys.path.insert(0, str(ROOT/'scripts'))
        runpy.run_path(str(SUPERVISOR), run_name='__main__')
        return 0
    receipt = prefix.with_name(prefix.name+'.launch.json')
    log = prefix.with_name(prefix.name+'.stdout.log')
    if any(p.exists() for p in (receipt, log, bootstrap)):
        raise ValueError('new detached launch record required')
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--detached-child',
        '--output', output.relative_to(ROOT).as_posix()]
    startup = subprocess.STARTUPINFO()
    startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW; startup.wShowWindow = 0
    flags = subprocess.CREATE_BREAKAWAY_FROM_JOB | subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONPATH='')
    parent_policy = job_policy()
    with log.open('xb', buffering=0) as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
            stdout=stream, stderr=subprocess.STDOUT, startupinfo=startup,
            creationflags=flags, close_fds=True)
    evidence = dict(status='STARTED_DETACHED_SUPERVISOR_PENDING_BOOTSTRAP_CHECK', pid=process.pid,
        command=command, creationflags=flags, parent_job_policy=parent_policy,
        original_deadline='2026-10-09T02:19:45+08:00', output=output.relative_to(ROOT).as_posix(),
        source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__).resolve(), SUPERVISOR)})
    write_new(receipt, evidence)
    print(json.dumps(evidence)); return 0


if __name__ == '__main__':
    raise SystemExit(main())

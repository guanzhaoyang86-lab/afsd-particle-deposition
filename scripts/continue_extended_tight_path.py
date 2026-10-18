"""Bounded extension from the fully audited smaller-step recovery endpoint.

This script does not import a solver, change a step size, retry a failed stage,
or promote a partial accepted prefix. Physical evidence stays in the new
registered extension series. --output contains only this invocation's supervision files.
Create OUTPUT/STOP_REQUESTED to stop at the next stage boundary.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'outputs/extended_tight_path_20261007'
LEGACY_LOGS = ROOT / 'logs/extended_tight_path_20261007'
RUNNER = ROOT / 'scripts/run_extended_tight_path_20261007.py'
READER = ROOT / 'scripts/readback_extended_tight_path_20261007.py'
PASS_COARSE = 'QUALIFIED_CONTINUATION_COARSE_NOT_ACCEPTED'
PASS_FINE = 'PASS_VISIBLE_CONTINUATION_WINDOW_NOT_DEPOSITION_VALIDATION'
PASS_READBACK = 'PASS_EXTENDED_TIGHT_PATH_READBACK_NOT_PHYSICAL_VALIDATION'
FIRST_DIRECTORY = ROOT / 'outputs/tight_path_refinement_20261007/window_004'
FIRST_READBACK_STATUS = 'PASS_TIGHT_PATH_REFINEMENT_READBACK_NOT_PHYSICAL_VALIDATION'
FIRST_TIME_S = 0.0521875
FIRST_GLOBAL_STEP = 823
STAGE_TIMEOUT = 630.0  # Original worker 540 s / original supervisor 600 s remain unchanged.
READBACK_TIMEOUT = 600.0
FINALIZE_RESERVE = 30.0
WINDOW_S = 16 * 0.0000625


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def atomic_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def project_path(value):
    path = Path(value)
    path = (ROOT / path).resolve() if not path.is_absolute() else path.resolve()
    require(path.is_relative_to(ROOT), 'path must stay inside the project')
    return path


def audit_endpoint(path, expected_status=PASS_READBACK):
    """A passing readback can certify a prefix without completing its window."""
    report = read(path)
    require(report.get('status') == expected_status and report.get('exit_code') == 0,
            'independent readback did not pass: ' + str(path))
    fields = report['final_accepted_fields']
    history = report['accepted_history']
    checkpoint = history[-1]['checkpoint'] if history else report['parent_checkpoint']
    checkpoint_sha = history[-1]['sha256'] if history else report['parent_sha256']
    require(digest(project_path(checkpoint)) == checkpoint_sha, 'audited checkpoint changed')
    return report, dict(readback_path=str(path), readback_sha256=digest(path),
        checkpoint_path=str(project_path(checkpoint)), checkpoint_sha256=checkpoint_sha,
        time_s=report['accepted_time_s'], window_completed=report['window_completed'],
        crossing_count=len(fields['crossing_ids']), crossing_ids=fields['crossing_ids'],
        crossing_mass_kg=fields['crossing_mass_kg'])


def parent_directory(number):
    return FIRST_DIRECTORY if number == 0 else BASE / f'window_{number:03d}'


def parent_endpoint(number):
    return audit_endpoint(parent_directory(number) / 'readback/summary.json',
        FIRST_READBACK_STATUS if number == 0 else PASS_READBACK)


def complete_parent(number):
    directory = parent_directory(number)
    report, endpoint = parent_endpoint(number)
    require(report.get('window_completed') is True,
            f'window {number} is incomplete; a PASS prefix cannot start the next window')
    require(report.get('top_event_process_screen_passed') is True,
            'previous window has no passing TOP process screen')
    require(abs(report['accepted_time_s'] - (FIRST_TIME_S + number * WINDOW_S)) < 1e-14,
            'parent time differs from the registered continuation endpoint')
    require(report['accepted_global_steps'] == FIRST_GLOBAL_STEP + number * 32,
            'parent global step differs from the registered continuation endpoint')
    if number == 0:
        require(project_path(endpoint['checkpoint_path']) == FIRST_DIRECTORY / 'fine/step_032/accepted.npz',
                'first parent is not the registered recovery window 4 endpoint')
    costs = {}
    for stage, expected, count in (('coarse', PASS_COARSE, 16), ('fine', PASS_FINE, 32)):
        manifest_path = directory / stage / 'manifest.json'
        manifest = read(manifest_path)
        require(manifest['status'] == expected and manifest['candidate_steps'] == count,
                'previous complete stage is absent: ' + stage)
        history_path = directory / stage / 'candidate_history.json'
        history = read(history_path)
        require(len(history) == count, 'previous candidate history is incomplete: ' + stage)
        measured_dt = (0.000015625 if stage == 'coarse' else 0.0000078125) if number == 0 else (
            0.0000625 if stage == 'coarse' else 0.00003125)
        require(all(row['dt_s'] == measured_dt for row in history),
                'parent candidate history has an unexpected time step: ' + stage)
        cost = float(manifest['wall_seconds'])
        require(math.isfinite(cost) and cost > 0, 'invalid measured stage cost')
        costs[stage] = dict(manifest_path=str(manifest_path), manifest_sha256=digest(manifest_path),
            candidate_history_sha256=digest(history_path), measured_wall_s=cost,
            candidate_steps=count, measured_s_per_step=cost / count,
            maximum_candidate_wall_s=max(float(row['wall_seconds']) for row in history),
            prior_forecast_wall_s=manifest.get('forecast_wall_s'),
            peak_process_working_set_bytes=manifest.get('peak_process_working_set_bytes'),
            measured_parent_dt_s=measured_dt,
            requested_next_dt_s=0.0000625 if stage == 'coarse' else 0.00003125)
        if stage == 'fine':
            accepted = read(directory / 'fine/accepted_history.json')
            require(manifest['accepted_steps'] == len(accepted) == len(report['accepted_history']) == 32,
                    'previous fine accepted history is not a full window')
            require(accepted[-1]['sha256'] == endpoint['checkpoint_sha256'],
                    'previous fine endpoint does not match the readback')
    costs['readback_measured_wall_s'] = report['wall_seconds']
    costs['particle_count'] = report['final_accepted_fields']['particle_count']
    return endpoint, costs


def stop_tree(process):
    """Only terminate the process group/tree started by this supervisor."""
    if os.name == 'nt':
        result = subprocess.run(['taskkill.exe', '/PID', str(process.pid), '/T', '/F'],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=20,
            creationflags=subprocess.CREATE_NO_WINDOW)
        outcome = dict(command=['taskkill.exe', '/PID', str(process.pid), '/T', '/F'],
                       exit_code=result.returncode, output=result.stdout)
    else:
        try:
            os.killpg(process.pid, signal.SIGKILL)
            outcome = dict(signal='SIGKILL', process_group=process.pid)
        except ProcessLookupError:
            outcome = dict(process_group=process.pid, already_exited=True)
    if process.poll() is None:
        process.kill()
    process.wait(timeout=20)
    return outcome


class WindowsJob:
    """Kill descendant workers even if their supervisor exits unexpectedly."""
    def __init__(self):
        self.handle = None
        if os.name != 'nt':
            return
        import ctypes
        from ctypes import wintypes

        class BasicLimits(ctypes.Structure):
            _fields_ = [('process_time', ctypes.c_int64), ('job_time', ctypes.c_int64),
                        ('flags', wintypes.DWORD), ('min_working_set', ctypes.c_size_t),
                        ('max_working_set', ctypes.c_size_t), ('active_limit', wintypes.DWORD),
                        ('affinity', ctypes.c_size_t), ('priority', wintypes.DWORD),
                        ('scheduling', wintypes.DWORD)]

        class Counters(ctypes.Structure):
            _fields_ = [(name, ctypes.c_uint64) for name in
                        ('read_ops', 'write_ops', 'other_ops', 'read_bytes', 'write_bytes', 'other_bytes')]

        class ExtendedLimits(ctypes.Structure):
            _fields_ = [('basic', BasicLimits), ('io', Counters),
                        ('process_memory', ctypes.c_size_t), ('job_memory', ctypes.c_size_t),
                        ('peak_process_memory', ctypes.c_size_t), ('peak_job_memory', ctypes.c_size_t)]

        self.api = ctypes.WinDLL('kernel32', use_last_error=True)
        self.api.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
        self.api.CreateJobObjectW.restype = wintypes.HANDLE
        self.api.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                                     ctypes.c_void_p, wintypes.DWORD]
        self.api.SetInformationJobObject.restype = wintypes.BOOL
        self.api.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        self.api.AssignProcessToJobObject.restype = wintypes.BOOL
        self.api.CloseHandle.argtypes = [wintypes.HANDLE]
        self.api.CloseHandle.restype = wintypes.BOOL
        self.handle = self.api.CreateJobObjectW(None, None)
        if not self.handle:
            raise ctypes.WinError(ctypes.get_last_error())
        limits = ExtendedLimits()
        limits.basic.flags = 0x2000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not self.api.SetInformationJobObject(self.handle, 9, ctypes.byref(limits), ctypes.sizeof(limits)):
            error = ctypes.WinError(ctypes.get_last_error())
            self.close()
            raise error

    def assign(self, process):
        if self.handle and not self.api.AssignProcessToJobObject(self.handle, int(process._handle)):
            import ctypes
            raise ctypes.WinError(ctypes.get_last_error())

    def close(self):
        if self.handle:
            self.api.CloseHandle(self.handle)
            self.handle = None


class Supervisor:
    def __init__(self, args):
        self.args = args
        self.output = project_path(args.output)
        require(self.output.is_relative_to(ROOT / 'outputs') and
                not self.output.is_relative_to(BASE), 'supervision output must be a separate outputs directory')
        self.output.mkdir(parents=True, exist_ok=True)
        # Exclusive creation prevents two invocations from claiming the same receipt.
        with (self.output / 'supervisor.pid').open('x', encoding='utf-8') as stream:
            stream.write(str(os.getpid()) + '\n')
        require(not (self.output / 'progress.json').exists(), 'refuse supervision receipt overwrite')
        self.started = time.monotonic()
        self.deadline = self.started + args.wall_limit
        self.progress = dict(status='PREFLIGHT', command=[sys.executable, *sys.argv],
            supervisor_pid=os.getpid(), active_child_pid=None, started_utc=timestamp(),
            start_window=args.start_window, end_window=args.end_window, target_s=None,
            wall_limit_s=args.wall_limit, stage_timeout_s=STAGE_TIMEOUT,
            original_worker_limit_s=540, original_runner_limit_s=600,
            readback_timeout_s=READBACK_TIMEOUT, physical_series=str(BASE),
            latest_audited=None, stages=[], source_sha256={str(path): digest(path)
                for path in (Path(__file__), RUNNER, READER)},
            deposition_validated=False, original_preflight_and_physical_gates_unchanged=True)
        self.save()

    def save(self):
        self.progress.update(updated_utc=timestamp(), elapsed_wall_s=time.monotonic() - self.started,
                             remaining_wall_s=max(0.0, self.deadline - time.monotonic()))
        atomic_json(self.output / 'progress.json', self.progress)

    def stop_requested(self):
        return (self.output / 'STOP_REQUESTED').exists()

    def finish(self, status, code, reason):
        self.progress.update(status=status, exit_code=code, reason=reason,
                             active_child_pid=None, finished_utc=timestamp())
        self.save()
        print(json.dumps({key: self.progress[key] for key in
            ('status', 'exit_code', 'reason', 'target_s', 'latest_audited', 'elapsed_wall_s')}, indent=2))
        return code

    def stage(self, number, stage):
        require(all(digest(Path(path)) == saved for path, saved in self.progress['source_sha256'].items()),
                'supervisor, runner, or reader source changed after launch')
        directory = BASE / f'window_{number:03d}'
        if stage == 'readback':
            require(not (directory / 'readback').exists(), 'refuse existing readback output')
            command = [sys.executable, '-B', str(READER), '--window', str(directory),
                       '--output', str(directory / 'readback')]
            evidence = directory / 'readback/summary.json'
            timeout = READBACK_TIMEOUT
        else:
            require(not (directory / stage).exists(), 'refuse existing physical stage: ' + stage)
            for suffix in ('.log', '_supervision.json'):
                require(not (LEGACY_LOGS / directory.name / (stage + suffix)).exists(),
                        'refuse existing original stage supervision')
            command = [sys.executable, '-B', str(RUNNER), '--window', str(number), '--stage', stage]
            evidence = directory / stage / 'manifest.json'
            timeout = STAGE_TIMEOUT
        require(self.deadline - time.monotonic() >= timeout + FINALIZE_RESERVE,
                'remaining wall budget cannot contain the complete stage')
        log = self.output / f'{directory.name}_{stage}.log'
        entry = dict(window=number, stage=stage, command=command, cwd=str(ROOT),
            started_utc=timestamp(), status='STARTING', exit_code=None, pid=None,
            timeout_s=timeout, log=str(log), evidence_path=str(evidence))
        self.progress['stages'].append(entry)
        self.progress.update(status='RUNNING', active_window=number, active_stage=stage)
        self.save()
        start = time.monotonic()
        process = None
        job = WindowsJob()
        try:
            with log.open('x', encoding='utf-8') as stream:
                process = subprocess.Popen(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                    env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'),
                    creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0,
                    start_new_session=os.name != 'nt')
                job.assign(process)
                entry.update(status='RUNNING', pid=process.pid)
                self.progress['active_child_pid'] = process.pid
                self.save()
                try:
                    entry['exit_code'] = process.wait(timeout=timeout)
                    entry['status'] = 'EXITED' if entry['exit_code'] == 0 else 'FAILED'
                except subprocess.TimeoutExpired:
                    entry.update(status='TIMEOUT', exit_code=124, cleanup=stop_tree(process))
        except BaseException:
            entry.update(status='SUPERVISION_ERROR', error=traceback.format_exc())
            if process is not None:
                entry['cleanup'] = stop_tree(process)
                entry['exit_code'] = process.returncode
            raise
        finally:
            job.close()
            entry.update(finished_utc=timestamp(), elapsed_wall_s=time.monotonic() - start)
            if evidence.exists():
                try:
                    manifest = read(evidence)
                    entry.update(evidence_status=manifest.get('status'), evidence_sha256=digest(evidence))
                except (OSError, ValueError):
                    entry['evidence_status'] = 'UNREADABLE'
            self.progress['active_child_pid'] = None
            self.save()
        return entry

    def run(self):
        previous = self.args.start_window - 1
        self.progress['target_s'] = FIRST_TIME_S + self.args.end_window * WINDOW_S
        # Record an independently accepted partial endpoint even when continuation is blocked.
        _, self.progress['latest_audited'] = parent_endpoint(previous)
        self.save()
        for number in range(self.args.start_window, self.args.end_window + 1):
            endpoint, costs = complete_parent(number - 1)
            self.progress.update(latest_audited=endpoint, measured_parent_costs=costs)
            self.save()
            if self.stop_requested():
                return self.finish('STOP_REQUESTED', 0, 'Stopped between stages before a new window.')
            # Reserve all three complete bounded stages before starting new physics.
            if self.deadline - time.monotonic() < 2 * STAGE_TIMEOUT + READBACK_TIMEOUT + FINALIZE_RESERVE:
                return self.finish('WALL_BUDGET_STOP', 0, 'Insufficient budget for another full window and readback.')
            directory = BASE / f'window_{number:03d}'
            require(not directory.exists(), 'refuse an existing window; no implicit resume or retry')
            stopped = None
            for stage, expected in (('coarse', PASS_COARSE), ('fine', PASS_FINE)):
                if self.stop_requested():
                    stopped = 'STOP_REQUESTED'
                    break
                result = self.stage(number, stage)
                if result['exit_code'] != 0 or result.get('evidence_status') != expected:
                    stopped = 'PHYSICAL_STAGE_FAILED'
                    break
            if not directory.exists():
                return self.finish(stopped or 'STOP_REQUESTED',
                    2 if stopped == 'PHYSICAL_STAGE_FAILED' else 0,
                    'No physical evidence directory was created; continuation stopped.')
            # A failed/partial stage is audited once, then stops. It is never retried.
            result = self.stage(number, 'readback')
            if result['exit_code'] != 0 or result.get('evidence_status') != PASS_READBACK:
                return self.finish('READBACK_FAILED', 2, 'Independent readback failed; no later window is allowed.')
            report, endpoint = audit_endpoint(directory / 'readback/summary.json')
            self.progress['latest_audited'] = endpoint
            self.save()
            if stopped:
                return self.finish(stopped, 0 if stopped == 'STOP_REQUESTED' else 2,
                    'Saved evidence audited; continuation stopped without promoting an incomplete window.')
            if report.get('window_completed') is not True:
                return self.finish('INCOMPLETE_WINDOW', 2, 'Readback passed only a partial prefix; continuation blocked.')
        return self.finish('TARGET_WINDOWS_COMPLETED', 0, 'Every requested window passed independent full-window readback.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start-window', type=int, default=1)
    parser.add_argument('--end-window', type=int, default=48)
    parser.add_argument('--output', default='outputs/extended_tight_path_supervision_20261007/v1')
    parser.add_argument('--wall-limit', type=float, default=18000.0, help='total bounded wall seconds (default 5 hours)')
    args = parser.parse_args()
    require(1 <= args.start_window <= args.end_window <= 48, 'require 1 <= start <= end <= 48')
    require(math.isfinite(args.wall_limit) and args.wall_limit > 0, 'wall limit must be finite and positive')
    supervisor = Supervisor(args)
    try:
        return supervisor.run()
    except KeyboardInterrupt:
        return supervisor.finish('INTERRUPTED', 130, 'Interrupted; the active child process tree was terminated.')
    except Exception:
        return supervisor.finish('BLOCKED_OR_FAILED', 2, traceback.format_exc())


if __name__ == '__main__':
    raise SystemExit(main())

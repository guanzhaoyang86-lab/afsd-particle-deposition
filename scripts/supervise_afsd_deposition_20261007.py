"""Finite stdlib-only supervision of the frozen screened deposition runner.

No physical module is imported. Every child performs the existing independent
source, history, geometry, energy and local coarse/fine gates. Only a complete
verified PASS stage becomes the next diagnostic parent.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import re
import signal
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT/'scripts/run_afsd_deposition_stage_remote_20261007.py'
CONTRACT = ROOT/'docs/DEPOSITION_SUPERVISION_CONTRACT_20261007.md'
TEST = ROOT/'scripts/test_supervise_afsd_deposition_20261007.py'
PASS_STATUSES = frozenset((
    'PASS_SCREENED_DEPOSITION_STAGE_DIAGNOSTIC_CONTINUATION',
    'PASS_FIRST_NEW_MATERIAL_UNDER_SHOULDER_DIAGNOSTIC_CONTINUATION',
    'PASS_SCREENED_PREFIX_BUDGET_BOUND_NOT_TARGET'))
COST_PASS = 'PASS_DEPOSITION_COST_MAP_NOT_PHYSICAL_VALIDATION'
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def inside(path, directory=ROOT):
    resolved = (directory/path).resolve()
    require(resolved.is_relative_to(directory.resolve()), 'path outside registered directory: '+str(path))
    return resolved


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.'+str(os.getpid())+'.pending')
    with temporary.open('w', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, path)


def verify(mapping):
    for name, digest in mapping.items():
        require(sha(inside(Path(name))) == digest, 'source/evidence hash changed: '+name)


def candidate_budget(cost_forecast_s, parent_maximum_s, maximum_quota=18):
    values = (float(cost_forecast_s), float(parent_maximum_s))
    require(all(math.isfinite(v) and v >= 0 for v in values), 'finite nonnegative measured cost required')
    measured = max(values)
    require(measured > 0 and isinstance(maximum_quota, int) and 3 <= maximum_quota <= 18,
            'positive measured cost and quota cap3..18 required')
    quota = min(maximum_quota, math.floor((450.-30.)/(1.5*measured)))
    while quota > 0 and 30.+1.5*quota*measured > 450.:
        quota -= 1
    return dict(measured_candidate_s=measured, candidate_quota=quota,
                max_macrosteps=min(6, quota//3), forecast_s=30.+1.5*quota*measured,
                runnable=quota >= 3)


def stop_after_stage(previous_time, current_time, status, shoulder_mass, target_time):
    if status not in PASS_STATUSES:
        return 'CHILD_FAILED_OR_NOT_SCREENED'
    if not (math.isfinite(current_time) and current_time > previous_time):
        return 'NO_PHYSICAL_TIME_INCREASE'
    if shoulder_mass > 0:
        return 'FIRST_SCREENED_NEW_MATERIAL_UNDER_SHOULDER'
    if current_time >= target_time:
        return 'GLOBAL_DIAGNOSTIC_TARGET_REACHED'
    return None


def verify_receipt(folder, statuses, frozen=None):
    meta, supervision = read(folder/'manifest.json'), read(folder/'supervision.json')
    require(meta['status'] in statuses and supervision['exit_code'] == 0
        and supervision['source_inputs_still_match']
        and supervision['manifest_sha256'] == sha(folder/'manifest.json')
        and supervision['log_sha256'] == sha(folder/'run.log'), 'failed, incomplete or altered child receipt')
    require(meta['accepted'] is False and meta['restart_allowed'] is False
        and meta['accepted_time_advanced_s'] == 0., 'formal acceptance is outside this supervisor')
    verify(meta['source_sha256']); verify(meta['artifacts_sha256'])
    if frozen is not None:
        require(meta['source_sha256'] == frozen, 'physical implementation changed across supervised stages')
    return meta


def verify_parent(folder, cost):
    meta = verify_receipt(folder, PASS_STATUSES, cost['source_sha256'])
    require(meta['classification'] == cost['classification'] and meta['screened_continuation_authorized']
        and meta['initial_checkpoint'] == cost['initial_checkpoint'] and meta['initial_sha256'] == cost['initial_sha256'],
        'screened parent belongs to a different independent model/origin')
    rows = meta['screened_history']
    require(rows and len(rows) == 2*len(meta['pair_screens']), 'complete screened fine pairs required')
    for offset in range(0, len(rows), 2):
        first, last = rows[offset:offset+2]
        require(first['half_index'] == 1 and last['half_index'] == 2
            and first['macro_index'] == last['macro_index'] and first['screen_path'] == last['screen_path'],
            'screened pair row identity mismatch')
        screen = read(inside(Path(first['screen_path'])))
        require(screen == meta['pair_screens'][offset//2] and screen['passed']
            and screen['global_screen']['pass'] and screen['local_screen']['passed'] and screen['crossing_ids_match']
            and screen['fine_checkpoint_sha256'] == [first['sha256'], last['sha256']]
            and screen['time_s'] == last['time_s'], 'unrelated PASS screen cannot endorse a candidate')
        for row in (first, last):
            require(sha(inside(Path(row['checkpoint']))) == row['sha256']
                and sha(inside(Path(row['observer_path']))) == row['observer_sha256']
                and sha(inside(Path(row['screen_path']))) == row['screen_sha256'], 'screened evidence changed')
    require(meta['screened_continuation_checkpoint'] == rows[-1]['checkpoint']
        and meta['screened_continuation_sha256'] == rows[-1]['sha256']
        and meta['screened_end_time_s'] == rows[-1]['time_s'], 'continuation must be the last screened fine state')
    return meta


def snapshot(meta, folder):
    last = meta['screened_history'][-1]
    observer = read(inside(Path(last['observer_path'])))
    supplied = observer['supplied_material']
    return dict(parent_stage=rel(folder), parent_manifest_sha256=sha(folder/'manifest.json'),
        physical_time_s=meta['screened_end_time_s'], diagnostic_steps=meta['screened_end_step'],
        last_screened_checkpoint=meta['screened_continuation_checkpoint'],
        last_screened_checkpoint_sha256=meta['screened_continuation_sha256'],
        eligible_crossed_mass_kg=supplied['eligible_crossed']['mass_kg'],
        new_material_under_shoulder_mass_kg=supplied['crossed_under_annular_shoulder']['mass_kg'],
        maximum_stage_new_material_under_shoulder_mass_kg=max(float(r['under_shoulder_mass_kg']) for r in meta['screened_history']),
        next_proposed_macro_dt_s=meta.get('next_proposed_macro_dt_s'))


def launch_options():
    if os.name == 'nt':
        return dict(creationflags=subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP)
    return dict(start_new_session=True)


def terminate_owned_tree(process, stream):
    """Stop only the child tree created by this supervisor, retaining all files."""
    if os.name == 'nt':
        result = subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
            stdout=stream, stderr=subprocess.STDOUT, timeout=15.,
            creationflags=subprocess.CREATE_NO_WINDOW)
        if result.returncode != 0 and process.poll() is None:
            process.kill()
        return dict(method='TASKKILL_OWNED_CHILD_PID_TREE', exit_code=result.returncode)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    return dict(method='KILL_OWNED_CHILD_PROCESS_GROUP', exit_code=0)


def plan(args):
    output = inside(Path(args.output)); cost_folder = inside(args.cost); parent = inside(args.parent)
    require(output.is_relative_to(ROOT/'outputs') and output != ROOT/'outputs' and not output.exists(),
            'new supervisor output directory required')
    require(1 <= args.max_stages <= 128 and math.isfinite(args.wall_seconds) and 650. <= args.wall_seconds <= 43200.,
            'bounded1..128 stages and650..43200 seconds required')
    require(math.isfinite(args.target_time) and args.target_time > 0.
        and 0 < args.min_dt <= args.max_dt <= .02 and 0 <= args.max_halvings <= 12,
        'finite bounded physical target and step limits required')
    python = (ROOT/args.python).resolve()
    require(python.is_file(), 'existing Python interpreter required')
    cost = verify_receipt(cost_folder, {COST_PASS})
    meta = verify_parent(parent, cost)
    match = re.fullmatch(r'stage_(\d+)', parent.name)
    require(match is not None, 'parent stage directory must use stage_NNN naming')
    first_index = int(match.group(1))+1
    next_folder = parent.parent/f'stage_{first_index:03d}'
    require(not next_folder.exists(), 'next stage already exists; use the newest verified parent')
    source = {rel(p): sha(p) for p in (Path(__file__).resolve(), CONTRACT, TEST, RUNNER)}
    source.update(cost['source_sha256'])
    inputs = {rel(p): sha(p) for p in (cost_folder/'manifest.json', cost_folder/'supervision.json',
        parent/'manifest.json', parent/'supervision.json')}
    return dict(schema='FINITE_STDLIB_DEPOSITION_SUPERVISOR_20261007', source_sha256=source, input_sha256=inputs,
        cost_directory=rel(cost_folder), initial_parent=rel(parent), stages_directory=rel(parent.parent),
        supervisor_output=rel(output), first_stage_index=first_index, python=str(python),
        configuration=cost['configuration'], component_receipt=cost['component_receipt'],
        max_stages=args.max_stages, total_wall_limit_s=args.wall_seconds, target_time_s=args.target_time,
        maximum_child_external_s=650., runner_internal_s=540., runner_external_s=600.,
        max_dt_s=args.max_dt, min_dt_s=args.min_dt, max_halvings=args.max_halvings,
        stop_on_first_shoulder=True, no_physical_module_imported=True,
        initial_snapshot=snapshot(meta, parent),
        next_budget=candidate_budget(cost['candidate_forecast_s'], meta.get('maximum_candidate_wall_s', 0.))), cost, meta


def supervise(protocol, cost, meta):
    output = ROOT/protocol['supervisor_output']; output.mkdir(parents=True, exist_ok=False)
    atomic_json(output/'protocol.json', protocol)
    began = time.perf_counter(); parent = ROOT/protocol['initial_parent']; stages = []
    progress = dict(schema=protocol['schema'], status='STARTING', completed_stage_count=0,
        physical_acceptance=False, accepted_time_advanced_s=0., wall_seconds=0., stages=stages,
        **protocol['initial_snapshot'])
    def update(**values):
        progress.update(values); progress['wall_seconds'] = time.perf_counter()-began
        atomic_json(output/'progress.json', progress)
    update()
    code = 0
    try:
        for offset in range(protocol['max_stages']):
            verify(protocol['source_sha256']); verify(protocol['input_sha256'])
            if progress['maximum_stage_new_material_under_shoulder_mass_kg'] > 0:
                update(status='STOPPED_FIRST_SCREENED_NEW_MATERIAL_UNDER_SHOULDER'); break
            if progress['physical_time_s'] >= protocol['target_time_s']:
                update(status='STOPPED_GLOBAL_DIAGNOSTIC_TARGET_REACHED'); break
            remaining = protocol['total_wall_limit_s']-(time.perf_counter()-began)
            if remaining < protocol['maximum_child_external_s']:
                update(status='STOPPED_TOTAL_WALL_RESERVE_BOUND', remaining_wall_seconds=remaining); break
            budget = candidate_budget(cost['candidate_forecast_s'], meta.get('maximum_candidate_wall_s', 0.))
            if not budget['runnable']:
                update(status='STOPPED_COST_REASSESSMENT_REQUIRED', next_budget=budget); break
            index = protocol['first_stage_index']+offset
            child = ROOT/protocol['stages_directory']/f'stage_{index:03d}'
            require(not child.exists(), 'next stage path became occupied; refusing overwrite or stale-parent branch')
            command = [protocol['python'], '-B', str(RUNNER), '--stage', 'stage', '--output', rel(child),
                '--checks', protocol['component_receipt'], '--config', protocol['configuration'],
                '--cost', protocol['cost_directory'], '--parent', rel(parent),
                '--end-time', str(protocol['target_time_s']), '--max-dt', str(protocol['max_dt_s']),
                '--min-dt', str(protocol['min_dt_s']), '--max-halvings', str(protocol['max_halvings']),
                '--max-candidate-attempts', str(budget['candidate_quota']),
                '--max-macrosteps', str(budget['max_macrosteps']), '--stop-on-first-shoulder']
            launch_log = output/f'stage_{index:03d}_launcher.log'
            update(status='RUNNING_CHILD', current_stage=rel(child), next_stage_index=index+1,
                   next_budget=budget, current_command=command)
            launched = time.perf_counter(); timeout = False; termination = None
            env = dict(os.environ, **{key: '1' for key in THREADS}); env['PYTHONPATH'] = ''
            with launch_log.open('x', encoding='utf-8') as stream:
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream,
                                           stderr=subprocess.STDOUT, **launch_options())
                update(child_pid=process.pid)
                while process.poll() is None:
                    if time.perf_counter()-launched >= protocol['maximum_child_external_s']:
                        timeout = True; termination = terminate_owned_tree(process, stream)
                        process.wait(timeout=20.); break
                    try:
                        process.wait(timeout=5.)
                    except subprocess.TimeoutExpired:
                        live = {}
                        manifest = child/'manifest.json'
                        if manifest.is_file():
                            try:
                                m = read(manifest)
                                live = dict(child_manifest_status=m.get('status'),
                                    child_screened_time_s=m.get('screened_end_time_s'),
                                    child_candidate_attempt_count=m.get('candidate_attempt_count'))
                            except (OSError, ValueError):
                                pass
                        update(**live)
                exit_code = process.returncode
            row = dict(stage=rel(child), command=command, child_pid=process.pid, exit_code=exit_code,
                timeout=timeout, termination=termination, wall_seconds=time.perf_counter()-launched,
                launcher_log=rel(launch_log), launcher_log_sha256=sha(launch_log), parent=rel(parent))
            stages.append(row); atomic_json(output/f'stage_{index:03d}_launch_receipt.json', row)
            if timeout or exit_code != 0:
                update(status='STOPPED_CHILD_TIMEOUT' if timeout else 'STOPPED_CHILD_FAILURE',
                       child_exit_code=exit_code, completed_stage_count=len(stages)-1)
                code = 2; break
            next_meta = verify_parent(child, cost)
            require(next_meta['parent_stage_manifest'] == rel(parent/'manifest.json')
                and next_meta['parent_stage_manifest_sha256'] == sha(parent/'manifest.json')
                and next_meta['parent_checkpoint_sha256'] == meta['screened_continuation_sha256'],
                'child did not continue this exact verified parent')
            new_snapshot = snapshot(next_meta, child)
            stop = stop_after_stage(progress['physical_time_s'], new_snapshot['physical_time_s'],
                next_meta['status'], new_snapshot['maximum_stage_new_material_under_shoulder_mass_kg'], protocol['target_time_s'])
            require(new_snapshot['diagnostic_steps'] > progress['diagnostic_steps'], 'child has no screened step increase')
            row.update(verified=True, manifest_sha256=sha(child/'manifest.json'),
                supervision_sha256=sha(child/'supervision.json'), status=next_meta['status'], **new_snapshot)
            atomic_json(output/f'stage_{index:03d}_launch_receipt.json', row)
            if stop == 'NO_PHYSICAL_TIME_INCREASE':
                update(status='STOPPED_NO_PHYSICAL_TIME_INCREASE'); code = 2; break
            parent, meta = child, next_meta
            update(status='VERIFIED_PASS_READY_FOR_NEXT_STAGE', completed_stage_count=offset+1,
                   current_stage=rel(child), **new_snapshot)
            if stop:
                update(status='STOPPED_'+stop); break
        else:
            update(status='STOPPED_DECLARED_STAGE_COUNT_BOUND')
        verify(protocol['source_sha256']); verify(protocol['input_sha256'])
    except Exception as error:
        code = 2; update(status='STOPPED_SUPERVISION_VERIFICATION_FAILURE', error=repr(error), traceback=traceback.format_exc())
    update(supervisor_exit_code=code, finished=True)
    print(json.dumps(dict(status=progress['status'], completed_stages=progress['completed_stage_count'],
        physical_time_s=progress['physical_time_s'], progress=rel(output/'progress.json'), exit_code=code)), flush=True)
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True); parser.add_argument('--cost', type=Path, required=True)
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--python', type=Path, default=ROOT/'.envs/pysph/python.exe')
    parser.add_argument('--max-stages', type=int, default=12); parser.add_argument('--wall-seconds', type=float, default=3600.)
    parser.add_argument('--target-time', type=float, default=.22); parser.add_argument('--max-dt', type=float, default=.005)
    parser.add_argument('--min-dt', type=float, default=.00000390625); parser.add_argument('--max-halvings', type=int, default=8)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args(); protocol, cost, meta = plan(args)
    if args.preflight_only:
        print(json.dumps({key: protocol[key] for key in ('initial_parent', 'first_stage_index', 'max_stages',
            'total_wall_limit_s', 'target_time_s', 'initial_snapshot', 'next_budget')}, indent=2)); return 0
    return supervise(protocol, cost, meta)


if __name__ == '__main__':
    raise SystemExit(main())

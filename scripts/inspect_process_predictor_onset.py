"""Predictor/onset saved-data provenance, correction-budget and physics audit.

No solver/material update is imported or called. Chronology checks are derived
from the frozen onset reader; only the branch correction trace has two stages.
The unchanged original Q7 physics reader audits actual saved physical fields.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
import inspect_process_plastic_friction as parent
from inspect_process_plastic_friction_recovery import (
    BASELINE, BASELINE_SHA256, BASELINE_SOURCES, RECOVERY_CONFIG, digest, read_json)

ROOT = Path(__file__).resolve().parents[1]
ONSET_CONFIG = 'configs/benchmarks/process_plastic_friction_onset.json'
PREDICTOR_CONFIG = 'configs/benchmarks/process_predictor_onset.json'
ONSET_SOURCES = {'scripts/probe_process_plastic_friction_onset.py', ONSET_CONFIG,
    'scripts/probe_process_plastic_friction_recovery.py', RECOVERY_CONFIG}
SUPPLEMENTAL_SOURCES = {'scripts/probe_process_predictor_onset.py', PREDICTOR_CONFIG,
    'scripts/probe_process_feasible_predictor.py'}


def mask(value, size, label):
    if not isinstance(value, list) or len(value) != size or any(type(v) is not bool for v in value):
        raise ValueError(label)
    return np.asarray(value, dtype=bool)


def finite(value, shape, label):
    a = np.asarray(value, dtype=float)
    if a.shape != shape or not np.isfinite(a).all():
        raise ValueError(label)
    return a


def audit_iterations(iterations, nc, status):
    """One shared 14-correction budget, not 14 predictor plus 14 cone."""
    if not isinstance(iterations, list) or len(iterations) > 15:
        raise ValueError('undeclared shared branch correction trace/budget')
    cone_seen = False
    factors = {2.**(-k) for k in range(9)}
    for index, row in enumerate(iterations):
        stage = row.get('stage')
        normals = finite(row.get('normal_force_N'), (nc,), 'branch normal inventory')
        if stage == 'INTERIOR_PREDICTOR':
            if cone_seen or type(row.get('corrections_used')) is not int or row['corrections_used'] != index:
                raise ValueError('predictor order/shared correction counter mismatch')
            residual = row.get('interior_force_relative')
            if ('root_relative' in row or not isinstance(residual, (int,float))
                    or not np.isfinite(residual) or residual < 0 or not np.any(normals < 0)):
                raise ValueError('invalid predictor trace or predictor promoted to root')
        elif stage == 'CONE_ROOT':
            cone_seen = True
            residual = row.get('root_relative')
            if (type(row.get('iteration')) is not int or row['iteration'] != index
                    or not isinstance(residual, (int,float)) or not np.isfinite(residual)
                    or residual < 0 or np.any(normals < 0)):
                raise ValueError('invalid cone root/shared correction trace')
        else:
            raise ValueError('unknown shared branch correction stage')
        factor = row.get('line_search_factor')
        if factor is not None:
            condition = row.get('tangent_condition')
            if (factor not in factors or index >= 14
                    or not isinstance(condition, (int,float))
                    or not np.isfinite(condition) or not 0 < condition <= 1e12):
                raise ValueError('invalid correction line search or shared budget overflow')
        elif index < len(iterations)-1:
            raise ValueError('continued branch without successful correction')
    if status == 'ROOT_RETURNED_NOT_ACCEPTED':
        if (not iterations or iterations[-1].get('stage') != 'CONE_ROOT'
                or iterations[-1]['root_relative'] > 1e-8
                or 'line_search_factor' in iterations[-1]):
            raise ValueError('root-returned candidate lacks converged cone trace')
    return dict(predictor_rows=sum(row['stage']=='INTERIOR_PREDICTOR' for row in iterations),
                corrections=sum('line_search_factor' in row for row in iterations))


def audit_attempts(path, report):
    """Validate trace chronology without treating a converged root as accepted."""
    files = {p.name:p for p in path.glob('attempt_*.json')}
    inventory = {name for name in report['artifact_sha256'] if name.startswith('attempt_') and name.endswith('.json')}
    if not files or set(files) != inventory:
        raise ValueError('attempt trace inventory missing or inconsistent')
    attempts = []
    for name, file in files.items():
        if report['artifact_sha256'][name] != digest(file):
            raise ValueError('attempt trace artifact changed')
        attempts.append((read_json(file), name))
    attempts.sort(key=lambda pair:pair[0]['requested_time_s'])
    n = report['increments_per_phase']
    grid = np.linspace(0., .2, 4*n+1)[1:]
    if len(attempts) > len(grid):
        raise ValueError('too many requested attempts')
    count, current, accepted, roots = 0, 0., [], 0
    rows = report['records']
    nc = report['contact_rows']
    for index, (attempt, name) in enumerate(attempts):
        requested = attempt.get('requested_time_s')
        if (attempt.get('classification') != 'ATTEMPT_TRACE_NOT_RESTART'
                or attempt.get('restart_allowed') is not False
                or requested != grid[index] or not current < requested
                or attempt.get('base_time_s') != current or attempt.get('base_step') != count
                or type(attempt.get('base_step')) is not int
                or name != f'attempt_{count:04d}_{requested:.17g}.json'):
            raise ValueError('attempt classification/base/request identity mismatch')
        trace = attempt.get('trace')
        if not isinstance(trace, list) or not trace:
            raise ValueError('missing attempt records')
        probe = selected = previous = None
        candidates = []
        lo, hi, lower_seen = current, requested, False
        for item in trace:
            stage = item.get('stage')
            if stage in ('BOUNDARY_PROBE', 'TARGET', 'INTERIOR_BISECTION'):
                time = item.get('time_s')
                if not isinstance(time, (int,float)) or not np.isfinite(time) or not current < time <= requested:
                    raise ValueError('candidate trace time outside current accepted interval')
                released = mask(item.get('released'), nc, 'candidate release inventory')
                if stage == 'BOUNDARY_PROBE':
                    if previous is not None:
                        raise ValueError('unexpected boundary probe order')
                    probe = item
                elif selected is None or not np.array_equal(released, selected):
                    raise ValueError('candidate release mask changed after classification')
                elif stage == 'TARGET' and (previous is None or previous.get('stage') != 'CLASSIFICATION' or time != requested):
                    raise ValueError('target trace order/time mismatch')
                elif stage == 'INTERIOR_BISECTION' and (previous is None or previous.get('stage') not in ('TARGET','BRACKET') or time != .5*(lo+hi)):
                    raise ValueError('interior bisection order/time mismatch')
                status = item.get('status')
                if status not in ('ROOT_RETURNED_NOT_ACCEPTED','BLOCKED'):
                    raise ValueError('candidate promoted or unknown candidate classification')
                iterations = item.get('iterations')
                audit_iterations(iterations, nc, status)
                if status == 'BLOCKED':
                    if item is not trace[-1] or attempt.get('error') != item.get('error'):
                        raise ValueError('blocked root followed by promoted state')
                else:
                    finite(item.get('cone_margin_N'), (nc,), 'candidate cone margin inventory')
                    finite(item.get('relative_increment'), (nc,2), 'candidate slip inventory')
                    candidates.append(item)
                    roots += 1
            elif stage == 'CLASSIFICATION':
                if previous is not probe or probe is None or probe.get('status') != 'ROOT_RETURNED_NOT_ACCEPTED' or item.get('base_time_s') != current:
                    raise ValueError('classification not based on current coupled probe')
                eligible = mask(item.get('eligible'), nc, 'eligibility inventory')
                inward = mask(item.get('resolved_inward'), nc, 'inward inventory')
                selected = mask(item.get('released'), nc, 'selected release inventory')
                margin = finite(item.get('initial_margin_N'), (nc,), 'initial cone margin inventory')
                if (not np.array_equal(eligible, probe['released']) or np.any(margin < -9e-5)
                        or np.any(eligible & (np.abs(margin) > 9e-5))
                        or not np.array_equal(inward, eligible & (np.asarray(probe['cone_margin_N']) > 9e-5))
                        or not np.array_equal(selected, eligible & ~inward)):
                    raise ValueError('onset classification arithmetic mismatch')
            elif stage == 'BRACKET':
                if previous is None or previous.get('stage') != 'INTERIOR_BISECTION' or previous.get('status') != 'ROOT_RETURNED_NOT_ACCEPTED':
                    raise ValueError('bracket not based on returned bisection root')
                valid = item.get('strict_feasible')
                if type(valid) is not bool:
                    raise ValueError('invalid strict bracket classification')
                if valid:
                    if np.any(np.asarray(previous['cone_margin_N'])[~selected] < 0):
                        raise ValueError('strict bracket accepts outside held-stick cone')
                    lo, lower_seen = previous['time_s'], True
                else:
                    hi = previous['time_s']
                if (item.get('lower_s') != lo or item.get('upper_s') != hi
                        or item.get('feasible_lower_seen') is not lower_seen):
                    raise ValueError('bracket state update mismatch')
            elif stage == 'ACCEPTED':
                time = item.get('time_s')
                kind = item.get('event_kind')
                if (type(item.get('step')) is not int or item['step'] != count+1
                        or not current < time <= requested
                        or not any(c['time_s'] == time for c in candidates)):
                    raise ValueError('accepted trace lacks corresponding returned candidate')
                if kind == 'REQUESTED':
                    if time != requested or previous is None or previous.get('stage') != 'TARGET':
                        raise ValueError('requested acceptance is not requested target')
                elif kind == 'CONE_LIMIT':
                    if not lower_seen or time != lo or not 0 <= hi-lo <= 1e-12:
                        raise ValueError('cone-limit acceptance lacks feasible narrow bracket')
                else:
                    raise ValueError('unknown accepted event kind')
                count += 1
                if count >= len(rows) or rows[count].get('step') != count or rows[count].get('physical_time_s') != time or rows[count].get('event_kind') != kind:
                    raise ValueError('accepted trace/report mismatch')
                accepted.append((count,time))
                current = time
                probe = selected = previous = None
                candidates = []
                lo, hi, lower_seen = current, requested, False
                continue
            else:
                raise ValueError('unknown attempt stage')
            previous = item
        error = attempt.get('error')
        if error is not None:
            if index != len(attempts)-1 or report.get('status') != 'BLOCKED' or error != report.get('error'):
                raise ValueError('rejected attempt was continued or failure metadata changed')
        elif current != requested:
            raise ValueError('successful attempt did not reach requested endpoint')
        if attempt.get('last_accepted_step') != count or attempt.get('last_accepted_time_s') != current:
            raise ValueError('attempt endpoint promoted rejected state')
    if count != report['accepted_steps'] or current != report['accepted_time_s'] or len(rows) != count+1:
        raise ValueError('accepted trace inventory differs from report')
    if report['status'] == 'J2_CONE_PATH_COMPLETE':
        if len(attempts) != len(grid) or attempts[-1][0].get('error') is not None:
            raise ValueError('complete path missing successful attempts')
    elif report['status'] != 'BLOCKED' or attempts[-1][0].get('error') != report.get('error') or report.get('error') is None:
        raise ValueError('blocked path missing failed attempt')
    return dict(status='PASS_ATTEMPT_CHRONOLOGY', attempts=len(attempts),
        accepted_trace_rows=len(accepted), returned_candidate_roots=roots,
        candidates_are_not_accepted_states=True)


def inspect(directory):
    path = Path(directory)
    report = read_json(path/'report.json')
    cfg = read_json(path/'resolved_config.json')
    sidecar_path = path/'predictor_provenance.json'
    sidecar = read_json(sidecar_path)
    onset, recovery = read_json(ROOT/ONSET_CONFIG), read_json(ROOT/RECOVERY_CONFIG)
    predictor = read_json(ROOT/PREDICTOR_CONFIG)
    budget = recovery.get('maximum_iterations', {})
    if (onset.get('variant') != 'Q7_BOUNDARY_ELIGIBILITY_V1'
            or onset.get('parent') != RECOVERY_CONFIG or onset.get('production_authorized') is not False
            or type(onset.get('maximum_iterations')) is not int or onset['maximum_iterations'] != 14
            or onset.get('run_timeout_s') != 60 or onset.get('permitted_cases') != [4,8]
            or recovery.get('variant') != 'Q7R_BOUNDED_NEWTON14'
            or recovery.get('production_authorized') is not False
            or type(budget.get('value')) is not int or budget['value'] != 14 or budget.get('unit') != 'count'
            or recovery.get('baseline_maximum_iterations') != 12 or recovery.get('run_timeout_s') != 60
            or recovery.get('permitted_cases') != [4,8]
            or recovery.get('parent') != 'configs/benchmarks/process_plastic_friction_q7.json'):
        raise ValueError('undeclared onset/recovery variant budget/scope')
    if (predictor.get('variant') != 'Q7_PREDICTOR_ONSET_V1'
            or predictor.get('parent') != ONSET_CONFIG or predictor.get('production_authorized') is not False
            or type(predictor.get('maximum_total_corrections')) is not int
            or predictor['maximum_total_corrections'] != 14
            or predictor.get('permitted_cases') != [4,8] or predictor.get('run_timeout_s') != 180):
        raise ValueError('undeclared predictor variant/shared budget/scope')
    if (sidecar.get('variant') != predictor['variant'] or sidecar.get('production_authorized') is not False
            or sidecar.get('configuration') != predictor):
        raise ValueError('predictor provenance/configuration mismatch')
    if sidecar.get('source_sha256') != {s:digest(ROOT/s) for s in SUPPLEMENTAL_SOURCES}:
        raise ValueError('supplemental source inventory/hash mismatch')
    if sidecar.get('onset_source_sha256') != {s:digest(ROOT/s) for s in ONSET_SOURCES}:
        raise ValueError('onset source inventory/hash mismatch')
    sources = sidecar.get('baseline_source_sha256')
    if (not isinstance(sources, dict) or set(sources) != BASELINE_SOURCES
            or sources != report.get('source_sha256') or sources != {s:digest(ROOT/s) for s in BASELINE_SOURCES}):
        raise ValueError('baseline source inventory/hash mismatch')
    if report.get('artifact_sha256', {}).get('predictor_provenance.json') != digest(sidecar_path):
        raise ValueError('predictor provenance missing from artifact inventory or changed')
    if digest(ROOT/BASELINE) != BASELINE_SHA256:
        raise ValueError('frozen original resolved configuration changed')
    expected = copy.deepcopy(read_json(ROOT/BASELINE))
    n = report.get('increments_per_phase')
    if type(n) is not int or n not in (4,8):
        raise ValueError('undeclared requested increment grid')
    if expected['plastic_parent']['numerical']['maximum_iterations'] != 12:
        raise ValueError('original baseline budget is not 12')
    expected['increments_per_phase'] = n
    expected['plastic_parent']['numerical']['maximum_iterations'] = 14
    expected['recovery_variant'] = recovery
    expected['onset_variant'] = onset
    expected['predictor_variant'] = predictor
    if cfg != expected:
        raise ValueError('resolved configuration changed outside declared budget/grid/metadata')
    trace_result = audit_attempts(path, report)
    result = parent.inspect(path)
    result['predictor_onset_audit'] = dict(status='PASS_PREDICTOR_ONSET_METADATA', variant=predictor['variant'],
        shared_correction_budget=14, onset_source_sha256=sidecar['onset_source_sha256'],
        baseline_resolved_config_sha256=BASELINE_SHA256,
        supplemental_source_sha256=sidecar['source_sha256'],
        physics_tolerances_boundaries_unchanged=True, attempt_trace=trace_result,
        caveat='Metadata/chronology and saved-data checks, not an independent continuum certificate.')
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inspect(args.run)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result['status'] == 'PASS_COMPLETE_READBACK' else 2


if __name__ == '__main__':
    raise SystemExit(main())

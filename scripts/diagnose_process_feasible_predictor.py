"""One frozen-target old/new predictor comparison; NO commit or path advance.

Every saved trial/root is REJECTED_DIAGNOSTIC_NOT_RESTART, including a returned
root. Predictor work is never accepted heat. Execute only under the externally
preregistered timeout; this script contains no retry or alternate target search.
"""
import argparse
import copy
import inspect as pyinspect
import json
from pathlib import Path
import time

import numpy as np

from diagnose_process_plastic_friction_onset import restore_onset_diagnostic
from diagnose_process_plastic_friction_q7 import assess_candidate, sha, state_signature
from inspect_process_plastic_friction_onset import inspect as inspect_onset
from probe_process_feasible_predictor import solve_with_feasible_predictor
from probe_process_plastic_friction_onset import ROOT


PREREG = 'docs/PROCESS_FEASIBLE_PREDICTOR_CONTRACT.md'
FROZEN_RUN = ROOT/'outputs/process_plastic_friction_onset_20261002/n3_i4'
EXPECTED_ERROR = 'negative normal outside fixed-contact cone scope'
TARGET = .125
BUDGET = 14
CLASSIFICATION = 'REJECTED_DIAGNOSTIC_NOT_RESTART'
DIAGNOSTIC_SOURCES = (
    'scripts/diagnose_process_feasible_predictor.py',
    'scripts/probe_process_feasible_predictor.py',
    'scripts/diagnose_process_plastic_friction_onset.py',
    'scripts/diagnose_process_plastic_friction_q7.py',
    'scripts/inspect_process_plastic_friction_onset.py',
    'scripts/inspect_process_plastic_friction_recovery.py',
    'scripts/inspect_process_plastic_friction.py',
    PREREG,
)


def json_safe(value):
    """Retain nonfinite evidence explicitly, never silently replace by zero."""
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, dict):
        return {str(k):json_safe(v) for k,v in value.items()}
    if isinstance(value, (tuple,list)):
        return [json_safe(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return dict(nonfinite_raw_value=repr(value))
    return value


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(json_safe(value), indent=2, allow_nan=False)+'\n')


def calling_context():
    """Record numerical call sites, not frame globals or sensitive locals."""
    result = []
    frame = pyinspect.currentframe().f_back
    try:
        while frame is not None:
            name = Path(frame.f_code.co_filename).name
            if name in ('probe_process_friction_cone.py', 'probe_process_feasible_predictor.py'):
                result.append(dict(source=name, function=frame.f_code.co_name,
                    source_line=frame.f_lineno,
                    context={key:copy.deepcopy(frame.f_locals[key])
                             for key in ('iteration','used','half','stage')
                             if key in frame.f_locals}))
            frame = frame.f_back
    finally:
        del frame
    return result


def save_candidate(path, b, c, released, kind, extra=None):
    arrays = {k:np.array(v,copy=True,order='K') for k,v in c.items()
              if isinstance(v,np.ndarray)}
    arrays.update(particle_id=b.pa.particle_id.copy(),
        contact_particle_id=b.pa.particle_id[b.nodes].copy(),
        contact_feature_id=np.ones(len(b.nodes),dtype=int), released=released.copy())
    if extra is not None:
        arrays.update(extra)
    np.savez_compressed(path, classification=np.array(CLASSIFICATION),
        restart_allowed=np.array(False), candidate_kind=np.array(kind),
        accepted_time_s=np.array(b.time), accepted_step=np.array(b.count),
        candidate_time_s=np.array(c['time_s']), **arrays)


def compare_branch(label, directory, run, report, history, released):
    """Exactly one branch/helper invocation on a freshly reconstructed owner."""
    b = restore_onset_diagnostic(run, report, history)
    before = state_signature(b)
    record_before = copy.deepcopy(b.record())
    actual_trial = b.trial
    raw_trace, observer_trace = [], []
    latest = None
    latest_object = None
    interior = b.free.copy()
    interior[b.contact_dofs] = False
    interior[b.tangent_dofs.ravel()] = False
    unknown = b.free.copy()
    unknown[b.contact_dofs] = False
    unknown[b.tangent_dofs[~released].ravel()] = False

    def capture_trial(displacement, t):
        nonlocal latest, latest_object
        c = actual_trial(displacement, t)
        latest = {k:np.array(v,copy=True,order='K') if isinstance(v,np.ndarray)
                  else copy.deepcopy(v) for k,v in c.items()}
        latest_object = c
        F, u = c['internal_force'], c['displacement']
        normal = F[b.nodes,0].copy()
        tau = -F[b.nodes,1:].copy()
        slip = u[b.nodes,1:]-b.displacement[b.nodes,1:]-(b.wall(t)-b.wall(b.time))[1:]
        contact = np.zeros_like(b.X)
        contact[b.nodes,0] = -normal
        contact[b.nodes,1:] = tau
        support = np.zeros_like(b.X)
        support[b.fixed] = -(F+contact)[b.fixed]
        mechanical = F+contact+support
        gap = b.wall(t)[0]-u[b.nodes,0]
        filename = f'REJECTED_{label}_raw_{len(raw_trace):04d}.npz'
        save_candidate(directory/filename, b, latest, released, 'RAW_TRIAL_NOT_ACCEPTED',
            dict(raw_normal_force_N=normal, raw_tangent_force_global_yz_N=tau,
                 raw_relative_increment_global_yz_m=slip, enforced_normal_gap_m=gap,
                 assigned_contact_force_N=contact, assigned_support_force_N=support,
                 mechanical_balance_if_contact_forces_assigned_N=mechanical))
        raw_trace.append(dict(trial_call=len(raw_trace)+1, stage='RAW_BEFORE_STATE',
            artifact=filename, classification=CLASSIFICATION, restart_allowed=False,
            target_time_s=float(t), candidate_base_step=int(c['base_step']),
            normal_force_N=normal.tolist(), negative_particle_id=b.pa.particle_id[b.nodes][normal<0].astype(int).tolist(),
            minimum_normal_N=float(np.min(normal)), enforced_normal_gap_m=gap.tolist(),
            interior_force_relative=float(np.linalg.norm(F.ravel()[interior])/b.force_scale),
            raw_internal_force_on_branch_unknown_relative=float(np.linalg.norm(F.ravel()[unknown])/b.force_scale),
            raw_assigned_contact_mechanical_relative=float(np.linalg.norm(mechanical.ravel()[b.free])/b.force_scale),
            residual_note='Raw mechanical diagnostics are NOT the coupled cone root residual.',
            finite_internal_force=bool(np.isfinite(F).all()), call_context=calling_context(),
            scalar_candidate_fields={k:v for k,v in latest.items()
                                     if isinstance(v,(str,int,float,bool,np.generic)) or v is None}))
        if t != TARGET or c['base_step'] != 16 or c['time_s'] != TARGET:
            raise ValueError('diagnostic left the single frozen target/base state')
        return c

    def observe(stage, c):
        if c is not latest_object or not raw_trace:
            raise ValueError('predictor observation has no matching actual raw trial')
        raw_trace[-1]['stage'] = stage
        observer_trace.append(dict(stage=stage, trial_call=raw_trace[-1]['trial_call']))

    b.trial = capture_trial
    tick = time.perf_counter()
    error = None
    candidate = None
    iterations = []
    status = 'BRANCH_ROOT_BLOCKED'
    assessment = None
    assessment_error = None
    try:
        if label == 'old':
            candidate = b.branch(TARGET, released.copy(), maximum_iterations=BUDGET)
        elif label == 'predictor':
            candidate = solve_with_feasible_predictor(b, TARGET, released.copy(),
                maximum_iterations=BUDGET, observe=observe)
        else:
            raise ValueError('undeclared diagnostic comparison')
        iterations = copy.deepcopy(candidate.get('iteration_trace', []))
        status = 'BRANCH_ROOT_RETURNED_NOT_ACCEPTED'
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        error = str(exc)
        candidate = getattr(exc, 'candidate', None)
        iterations = copy.deepcopy(getattr(exc, 'iterations', []))
    elapsed = time.perf_counter()-tick
    # A raw candidate is retained even when state() throws before constructing
    # contact inventory. No second solve or re-evaluation is used to invent it.
    if candidate is None:
        candidate = latest
    candidate_artifact = None
    if candidate is not None:
        candidate_artifact = f'REJECTED_{label}_final_candidate.npz'
        kind = 'RETURNED_ROOT_NOT_ACCEPTED' if status == 'BRANCH_ROOT_RETURNED_NOT_ACCEPTED' else 'BLOCKED_CANDIDATE_NOT_ACCEPTED'
        save_candidate(directory/candidate_artifact, b, candidate, released, kind)
        assess_keys = ('normal_force','tangent_force','relative_increment','residual',
                       'law_residual_N','contact_force','equilibrium_relative','plastic_increment_J')
        if all(key in candidate for key in assess_keys):
            try:
                # The original diagnostic checks forced-stick overcone,
                # unsplit boundary reversal and signed raw contact work.
                assessment = assess_candidate(b, candidate, released)
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                assessment_error = str(exc)
        else:
            assessment_error = 'NOT_AVAILABLE_RAW_TRIAL_ONLY_NO_CONTACT_STATE'
    else:
        assessment_error = 'NOT_AVAILABLE_NO_RETURNED_RAW_TRIAL'
    branch_status = 'ROOT_RETURNED' if status == 'BRANCH_ROOT_RETURNED_NOT_ACCEPTED' else 'BLOCKED'
    if branch_status == 'ROOT_RETURNED' and assessment is not None:
        if assessment['event_split_required']:
            status = 'ROOT_RETURNED_EVENT_SPLIT_REQUIRED_NOT_ACCEPTED'
        elif (not all(assessment[key] for key in ('root_pass','all_contact_root_pass',
                'forced_stick_feasible','released_flow_feasible','branch_cone_feasible'))
                or assessment['minimum_raw_candidate_friction_J'] < -b.limits['negative_dissipation_tolerance_J']):
            status = 'ROOT_RETURNED_INFEASIBLE_NOT_ACCEPTED'
    after = state_signature(b)
    record_after = b.record()
    unchanged = before == after and record_before == record_after == report['records'][-1]
    result = dict(comparison=label, status=status, branch_status=branch_status,
        error=error, classification=CLASSIFICATION,
        restart_allowed=False, branch_invocations=1, target_time_s=TARGET,
        maximum_total_corrections=BUDGET, accepted_step=b.count, accepted_time_s=b.time,
        accepted_record_exact=record_before == report['records'][-1],
        accepted_record_before=record_before, accepted_record_after=record_after,
        accepted_state_before=before, accepted_state_after=after,
        accepted_state_unchanged=unchanged, raw_trial_calls=len(raw_trace),
        raw_trace=raw_trace, observer_trace=observer_trace, iterations=iterations,
        returned_corrections=None if candidate is None else candidate.get('iterations'),
        predictor_corrections=None if candidate is None else candidate.get('predictor_corrections'),
        cone_corrections=None if candidate is None else candidate.get('cone_corrections'),
        final_candidate_artifact=candidate_artifact, assessment=assessment,
        assessment_error=assessment_error, wall_seconds=elapsed,
        candidate_work_is_accepted_heat=False, accepted_progress=False,
        warning='All trials/roots remain diagnostic. Predictor energy differences are not heat; original event/commit gates were not bypassed or executed.')
    write_json(directory/f'{label}_comparison.json', result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(argv)
    if args.run.resolve() != FROZEN_RUN.resolve():
        raise ValueError('only the frozen Onset n3_i4 diagnostic is authorized')
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    readback = inspect_onset(args.run)
    if readback['status'] != 'PASS_PARTIAL_READBACK':
        raise ValueError('requires frozen original partial run')
    report = json.loads((args.run/'report.json').read_text(encoding='utf-8'))
    manifest = json.loads((args.run/'REJECTED_manifest.json').read_text(encoding='utf-8'))
    attempt = json.loads((args.run/'attempt_0016_0.125.json').read_text(encoding='utf-8'))
    provenance = json.loads((args.run/'onset_provenance.json').read_text(encoding='utf-8'))
    if (report['accepted_steps'] != 16 or report['accepted_time_s'] != .1125
            or report['error'] != EXPECTED_ERROR or manifest['candidate_available'] is not False
            or manifest['attempted_time_s'] != TARGET or attempt['error'] != EXPECTED_ERROR):
        raise ValueError('not the preregistered first negative-normal failure')
    targets = [row for row in attempt['trace']
               if row['stage'] == 'TARGET' and row['status'] == 'BLOCKED']
    if len(targets) != 1 or targets[0]['time_s'] != TARGET:
        raise ValueError('expected exactly one frozen failing target')
    values = targets[0]['released']
    if len(values) != 9 or any(type(v) is not bool for v in values):
        raise ValueError('invalid saved contact release inventory')
    released = np.asarray(values, dtype=bool)
    sources = dict(report['source_sha256'], **provenance['source_sha256'])
    diagnostic_sources = {name:sha(ROOT/name) for name in DIAGNOSTIC_SOURCES}
    evidence = dict(report['artifact_sha256'])
    evidence['report.json'] = sha(args.run/'report.json')
    for name, digest in sources.items():
        if sha(ROOT/name) != digest:
            raise ValueError('frozen source mismatch: '+name)
    for name, digest in evidence.items():
        if sha(args.run/name) != digest:
            raise ValueError('input artifact mismatch: '+name)
    with np.load(args.run/'history.npz', allow_pickle=False) as archive:
        history = {k:archive[k].copy() for k in archive.files}
    input_record = dict(classification='FROZEN_INPUT_ONLY_NOT_RESTART',
        restart_allowed=False, run=str(args.run.resolve()), target_time_s=TARGET,
        accepted_step=16, accepted_time_s=.1125, released_mask=released.tolist(),
        saved_target=targets[0], input_sha256=evidence, source_sha256=sources,
        diagnostic_source_sha256=diagnostic_sources,
        configuration_sha256=evidence['resolved_config.json'],
        preregistration_sha256=diagnostic_sources[PREREG], readback=readback)
    write_json(args.output/'inputs.json', input_record)
    comparisons = []
    for label in ('old', 'predictor'):
        comparison = compare_branch(label, args.output, args.run, report, history, released)
        comparisons.append(comparison)
        if not comparison['accepted_state_unchanged']:
            # Never run the other branch after evidence of accepted mutation.
            break
    input_unchanged = all(sha(args.run/name) == digest for name,digest in evidence.items())
    source_unchanged = all(sha(ROOT/name) == digest for name,digest in dict(sources, **diagnostic_sources).items())
    old = comparisons[0]
    old_reproduced = (old['status'] == 'BRANCH_ROOT_BLOCKED' and old['error'] == EXPECTED_ERROR
                      and old['raw_trial_calls'] == 1)
    new = comparisons[-1] if len(comparisons) == 2 else None
    intact = (input_unchanged and source_unchanged
              and all(c['accepted_state_unchanged'] for c in comparisons))
    new_root = new is not None and new['status'] == 'BRANCH_ROOT_RETURNED_NOT_ACCEPTED'
    assessed = new is not None and new['assessment'] is not None and new['assessment_error'] is None
    exit_code = 0 if intact and old_reproduced and new_root and assessed else 2
    result_status = 'DIAGNOSTIC_COMPLETE_NOT_PATH_COMPLETE' if exit_code == 0 else 'DIAGNOSTIC_BLOCKED'
    if intact and old_reproduced and new is not None and new['status'] == 'ROOT_RETURNED_EVENT_SPLIT_REQUIRED_NOT_ACCEPTED':
        result_status = new['status']
    result = dict(status=result_status,
        classification='PURE_FROZEN_TARGET_PREDICTOR_COMPARISON_NOT_RESTART',
        production_authorized=False, restart_allowed=False, accepted_progress=False,
        accepted_time_s=.1125, accepted_step=16, target_candidate_time_s=TARGET,
        comparison_invocations=len(comparisons), maximum_total_corrections_per_branch=BUDGET,
        expected_old_guard_reproduced=old_reproduced,
        expected_old_guard_status='EXPECTED_OLD_GUARD_RED' if old_reproduced else 'UNEXPECTED_OLD_REPLAY',
        input_artifacts_unchanged=input_unchanged, source_unchanged=source_unchanged,
        accepted_states_unchanged=all(c['accepted_state_unchanged'] for c in comparisons),
        input_sha256=evidence, source_sha256=sources, diagnostic_source_sha256=diagnostic_sources,
        configuration_sha256=evidence['resolved_config.json'], preregistration_sha256=diagnostic_sources[PREREG],
        released_mask=released.tolist(), comparisons=comparisons,
        exit_code=exit_code, wall_seconds=time.perf_counter()-started,
        warning='A returned and assessed branch root is still NOT_ACCEPTED. Check forced_stick_feasible and event_split_required; raw candidate friction/plastic work is not accepted heat. No event search, commit, advance, target substitution or full path was run.',
        artifact_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    write_json(args.output/'diagnostic.json', result)
    print(json.dumps(json_safe({key:result[key] for key in ('status','exit_code','accepted_time_s',
        'accepted_step','target_candidate_time_s','comparison_invocations','expected_old_guard_reproduced',
        'input_artifacts_unchanged','source_unchanged','accepted_states_unchanged','wall_seconds')}), indent=2))
    for item in comparisons:
        print(json.dumps(json_safe({key:item[key] for key in ('comparison','status','error','raw_trial_calls',
            'returned_corrections','predictor_corrections','cone_corrections','assessment','assessment_error')}), indent=2))
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())

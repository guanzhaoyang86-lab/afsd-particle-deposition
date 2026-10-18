"""Pure Q7R event replay: no accepted progress, no altered contact checks."""
import argparse
import json
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np

import probe_process_friction_cone as cone
from diagnose_process_plastic_friction_q7 import restore_diagnostic, state_signature, sha
from inspect_process_plastic_friction_recovery import inspect
from probe_process_plastic_friction import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    readback = inspect(args.run)
    if readback['status'] != 'PASS_PARTIAL_READBACK':
        raise ValueError('this replay requires the original partial Q7R evidence')
    report = json.loads((args.run/'report.json').read_text(encoding='utf-8'))
    rejected = json.loads((args.run/'REJECTED_manifest.json').read_text(encoding='utf-8'))
    if (report['error'] != 'cone event not bracketed' or report['accepted_steps'] != 6
            or report['accepted_time_s'] != .05 or rejected['candidate_time_s'] != .0625):
        raise ValueError('this diagnostic is restricted to the first Q7R event failure')
    with np.load(args.run/'history.npz', allow_pickle=False) as archive:
        history = {k:archive[k].copy() for k in archive.files}
    with np.load(args.run/'REJECTED_candidate.npz', allow_pickle=False) as archive:
        target_arrays = {k:archive[k].copy() for k in archive.files}
    evidence = {name:sha(args.run/name) for name in report['artifact_sha256']}
    evidence['report.json'] = sha(args.run/'report.json')
    body = restore_diagnostic(args.run, report, history)
    body.qcfg['numerical']['maximum_iterations'] = 14
    recovery = json.loads((ROOT/'configs/benchmarks/process_plastic_friction_q7_recovery.json').read_text(encoding='utf-8'))
    body.pfcfg['recovery_variant'] = recovery
    before = state_signature(body)
    branch_original = body.branch
    feasible_original = cone.feasible_cone
    ids = body.pa.particle_id[body.nodes].astype(int)
    branches, checks = [], []
    last_branch = None
    target_replay = None

    def metrics(c, released):
        # Independent saved-vector arithmetic, not the contact response function.
        normal = np.asarray(c['normal_force'])
        tau = np.asarray(c['tangent_force'])
        slip = np.asarray(c['relative_increment'])
        distance = np.sqrt(np.sum(slip*slip, axis=1))
        capacity = body.mu*normal
        margin = capacity-np.sqrt(np.sum(tau*tau, axis=1))
        flow = np.sqrt(np.sum((distance[:,None]*tau+capacity[:,None]*slip)**2, axis=1))
        allowance = body.limits['force_tolerance_N']*distance+capacity*body.limits['slip_tolerance_m']
        released = np.asarray(released, dtype=bool)
        moving = released & (distance>body.limits['slip_tolerance_m'])
        normal_block = (~released) & (margin < -body.limits['force_tolerance_N'])
        strict_block = (~released) & (margin < 0.)
        released_block = released & (margin < -body.limits['force_tolerance_N'])
        flow_block = moving & (flow>allowance)
        return dict(particle_id=ids.tolist(), normal_force_N=normal.tolist(),
            tangent_force_global_yz_N=tau.tolist(), slip_global_yz_m=slip.tolist(),
            cone_margin_N=margin.tolist(), flow_residual_N_m=flow.tolist(),
            flow_allowance_N_m=allowance.tolist(), released=released.tolist(),
            normal_forced_stick_blocker_ids=ids[normal_block].tolist(),
            strict_forced_stick_blocker_ids=ids[strict_block].tolist(),
            released_cone_blocker_ids=ids[released_block].tolist(),
            released_flow_blocker_ids=ids[flow_block].tolist(),
            negative_normal_ids=ids[normal<0].tolist(),
            independent_normal_feasible=bool(not (normal_block|released_block|flow_block).any() and np.all(normal>=0)),
            independent_strict_feasible=bool(not (strict_block|released_block|flow_block).any() and np.all(normal>=0)),
            contact_root_norm_N=np.sqrt(np.sum(c['law_residual_N']**2, axis=1)).tolist())

    def branch(t, released, maximum_iterations=None):
        nonlocal last_branch, target_replay
        if maximum_iterations is not None and maximum_iterations>14:
            raise ValueError('diagnostic may not exceed the declared budget')
        entry = dict(call=len(branches), time_s=float(t), released=np.asarray(released).tolist(),
            budget=14 if maximum_iterations is None else maximum_iterations)
        branches.append(entry)
        last_branch = entry['call']
        try:
            c = branch_original(t, released, maximum_iterations=maximum_iterations)
        except Exception as exc:
            entry.update(status='BLOCKED_BRANCH', error=str(exc), iterations=getattr(exc,'iterations',[]))
            c = getattr(exc,'candidate',None)
            if c is not None:entry['candidate_metrics']=metrics(c,released)
            raise
        entry.update(status='RETURNED_ROOT_NOT_ACCEPTED', iterations=c['iteration_trace'],
            candidate_metrics=metrics(c,released))
        if t == .0625 and np.array_equal(released,target_arrays['released']):
            target_replay = {k:bool(c[k].dtype==v.dtype and c[k].shape==v.shape and c[k].tobytes()==v.tobytes())
                             for k,v in target_arrays.items()}
        return c

    def feasible(c, released, mu, limits, strict=False):
        entry = dict(call=len(checks), branch_call=last_branch, time_s=float(c['time_s']), strict=bool(strict),
            metrics=metrics(c,released))
        checks.append(entry)
        try:
            result = feasible_original(c,released,mu,limits,strict=strict)
        except Exception as exc:
            entry.update(error=str(exc));raise
        entry['actual_return'] = bool(result)
        independent = entry['metrics']['independent_strict_feasible' if strict else 'independent_normal_feasible']
        if independent != bool(result):raise ValueError('independent feasibility arithmetic differs')
        return result

    body.branch = branch
    caught = None
    with patch.object(cone,'feasible_cone',feasible):
        try:
            # body.advance and commit remain the forbidden reconstruction stubs.
            cone.advance_cone(body,.0625)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            caught = str(exc)
    after = state_signature(body)
    if before != after:raise ValueError('pure event replay mutated accepted state')
    if target_replay is None or not all(target_replay.values()):raise ValueError('target candidate did not replay bytewise')
    if any(sha(args.run/name)!=digest for name,digest in evidence.items()):raise ValueError('input artifact changed')
    lower,upper=.05,.0625
    lower_found=False
    inferred=[]
    for check in checks:
        if not check['strict']:continue
        if check['actual_return']:lower=check['time_s'];lower_found=True
        else:upper=check['time_s']
        inferred.append(dict(check_call=check['call'],lower_s=lower,upper_s=upper,
            width_s=upper-lower,feasible_lower_candidate_seen=lower_found))
    if caught != 'cone event not bracketed':raise ValueError('replay did not reproduce the recorded event failure: '+str(caught))
    out = dict(classification='PURE_EVENT_REPLAY_REJECTED_DIAGNOSTIC_NOT_RESTART',
        production_authorized=False,restart_allowed=False,status='DIAGNOSTIC_COMPLETE_NOT_PATH_COMPLETE',
        error=caught,accepted_time_s=body.time,accepted_step=body.count,accepted_state_unchanged=True,
        accepted_state_before=before,accepted_state_after=after,readback_status=readback['status'],
        target_candidate_byte_replay=target_replay,branches=branches,feasibility_checks=checks,
        inferred_bisection=dict(note='Reconstructed from observed strict checks and unchanged advance_cone source, not saved original locals.',
            trace=inferred,lower_candidate_seen=lower_found),
        input_sha256=evidence,source_sha256={name:sha(ROOT/name) for name in
            ('scripts/diagnose_process_plastic_friction_recovery.py','scripts/diagnose_process_plastic_friction_q7.py',
             'scripts/inspect_process_plastic_friction_recovery.py','scripts/inspect_process_plastic_friction.py',
             'scripts/probe_process_friction_cone.py','scripts/probe_process_plastic_friction_recovery.py',
             'configs/benchmarks/process_plastic_friction_q7_recovery.json')},
        preregistration_sha256=sha(ROOT/'docs/PROCESS_PLASTIC_FRICTION_Q7_RECOVERY_EVENT_DIAGNOSTIC.md'),
        wall_seconds=time.perf_counter()-started)
    with (args.output/'diagnostic.json').open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=out['status'],error=caught,branches=len(branches),checks=len(checks),
        strict_checks=len(inferred),last_inferred_bounds=inferred[-1] if inferred else None,
        accepted_unchanged=True,target_exact=True,wall_seconds=out['wall_seconds']),indent=2))
    return 0


if __name__=='__main__':raise SystemExit(main())

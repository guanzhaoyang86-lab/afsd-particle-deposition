"""One unchanged branch replay to retain first negative-normal raw candidate."""
import argparse
import copy
import inspect
import json
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np
from pysph.solver.utils import load

import probe_process_friction_cone as cone
from probe_process_plastic_friction_onset import OnsetBlock, ROOT
from diagnose_process_plastic_friction_q7 import state_signature, sha

PREREG = 'docs/PROCESS_PLASTIC_FRICTION_ONSET_FAILURE_DIAGNOSTIC.md'


def restore_onset_diagnostic(path, report, history):
    """Raw read-only diagnostic reconstruction; NOT a production restore API."""
    b = OnsetBlock()
    saved = load(str(path/'checkpoint.npz'))
    pa, meta = saved['arrays']['solid'], saved['solver_data']
    if (meta['restart_allowed'] is not False or meta['step'] != 16
            or meta['physical_time_s'] != .1125):
        raise ValueError('restricted diagnostic checkpoint differs')
    if set(pa.properties) != set(b.pa.properties):
        raise ValueError('particle inventory mismatch')
    for name in pa.properties:
        actual, original = pa.properties[name].get_npy_array(), b.pa.properties[name].get_npy_array()
        if actual.shape != original.shape or actual.dtype != original.dtype:
            raise ValueError('particle schema mismatch: '+name)
    for expected, key in [(b.X, 'reference_position_m'), (b.V, 'reference_volume_m3'),
                          (pa.particle_id, 'particle_id')]:
        if not np.array_equal(expected, history[key]):
            raise ValueError('reference/ID mismatch: '+key)
    b.pa = pa
    mapping = dict(displacement='displacement_m', stress='stress_Pa', normal_force='normal_force_N',
        tangent_force='tangent_force_global_yz_N', modes='interval_modes', contact_force='contact_force_N',
        support_force='support_force_N', friction='friction_dissipation_J', last_friction='interval_friction_J',
        last_slip='interval_relative_slip_global_yz_m', last_impulse='tangent_impulse_global_yz_N_s',
        last_relative_rate='relative_kinematic_rate_global_yz_m_s', rate='kinematic_interval_rate_m_s',
        last_law_residual='contact_law_residual_global_yz_N', last_plastic='interval_plastic_work_J',
        last_return='interval_return_diagnostic_J')
    for attr, key in mapping.items():
        # Match zeros_like(X) F-layout used by actual state(), including reductions.
        order = 'F' if attr in ('contact_force', 'support_force') else 'C'
        setattr(b, attr, np.array(history[key][-1], copy=True, order=order))
    row = report['records'][-1]
    scalar = dict(time='physical_time_s', count='step', wall_work='wall_work_J',
        contact_work='contact_material_work_J', support_work='support_work_J', internal_work='stress_work_J',
        last_wall_increment='interval_wall_work_J', last_contact_increment='interval_contact_material_work_J',
        last_residual='free_force_relative', last_iterations='corrections', last_trace='iteration_trace',
        last_force_defect='internal_force_relative', last_torque_defect='internal_torque_relative',
        accepted_floor_count='accepted_rate_floor_evaluations', last_event='event_kind')
    for attr, key in scalar.items():
        setattr(b, attr, copy.deepcopy(row[key]))
    b.constraint_work = b.contact_work+b.support_work
    b.reaction = b.contact_force+b.support_force
    b.events = copy.deepcopy(report['events'])
    reconstructed = b.record()
    if reconstructed != row:
        differences = {k:dict(reconstructed=reconstructed.get(k), saved=row.get(k))
                       for k in set(reconstructed)|set(row) if reconstructed.get(k) != row.get(k)}
        raise ValueError('accepted record mismatch: '+json.dumps(differences))
    # Require checkpoint material state to coincide with final accepted history.
    for actual, key in [(b.pa.ep, 'equivalent_plastic_strain'), (b.pa.plastic_work, 'plastic_work_J'),
        (b.pa.return_dissipation, 'return_dissipation_diagnostic_J'),
        (np.array([b.pa.x,b.pa.y,b.pa.z]).T, 'position_m')]:
        if not np.array_equal(actual, history[key][-1]):
            raise ValueError('checkpoint/history mismatch: '+key)
    with np.load(path/'accepted_0016.npz', allow_pickle=False) as journal:
        if json.loads(str(journal['record_json'])) != row:
            raise ValueError('accepted journal record mismatch')
        for attr, key in [('displacement','u'), ('stress','stress'), ('normal_force','normal'),
                          ('tangent_force','tau'), ('friction','q'), ('last_plastic','dplastic')]:
            if not np.array_equal(getattr(b,attr),journal[key]):
                raise ValueError('accepted journal array mismatch: '+key)
    def forbidden(*args, **kwargs):
        raise RuntimeError('diagnostic cannot commit/advance accepted state')
    b.commit = forbidden
    b.advance = forbidden
    return b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    report = json.loads((args.run/'report.json').read_text(encoding='utf-8'))
    manifest = json.loads((args.run/'REJECTED_manifest.json').read_text(encoding='utf-8'))
    attempt = json.loads((args.run/'attempt_0016_0.125.json').read_text(encoding='utf-8'))
    provenance = json.loads((args.run/'onset_provenance.json').read_text(encoding='utf-8'))
    expected_error = 'negative normal outside fixed-contact cone scope'
    if (report['accepted_steps'] != 16 or report['accepted_time_s'] != .1125
            or report['error'] != expected_error or manifest['candidate_available'] is not False
            or manifest['attempted_time_s'] != .125 or attempt['error'] != expected_error):
        raise ValueError('not the preregistered first onset negative-normal failure')
    source = dict(report['source_sha256'], **provenance['source_sha256'])
    for name, digest in source.items():
        if sha(ROOT/name) != digest:
            raise ValueError('source mismatch: '+name)
    evidence = dict(report['artifact_sha256'])
    evidence['report.json'] = sha(args.run/'report.json')
    for name, digest in evidence.items():
        if sha(args.run/name) != digest:
            raise ValueError('input artifact mismatch: '+name)
    targets = [r for r in attempt['trace'] if r['stage'] == 'TARGET' and r['status'] == 'BLOCKED']
    if len(targets) != 1:
        raise ValueError('expected exactly one failing target')
    target = targets[0]
    released = np.array(target['released'], dtype=bool)
    with np.load(args.run/'history.npz', allow_pickle=False) as archive:
        history = {k:archive[k].copy() for k in archive.files}
    b = restore_onset_diagnostic(args.run, report, history)
    before = state_signature(b)
    actual_trial, actual_response = b.trial, cone.cone_response
    latest_trial = None
    first = None
    trial_calls = 0
    response_calls = []

    def capture_trial(displacement, t):
        nonlocal latest_trial, trial_calls
        c = actual_trial(displacement,t)
        trial_calls += 1
        latest_trial = {k:np.array(v,copy=True,order='K') if isinstance(v,np.ndarray) else copy.deepcopy(v)
                        for k,v in c.items()}
        return c

    def capture_response(tau, normal, slip, mu, kappa, basis=None):
        nonlocal first
        branch_frame = inspect.currentframe().f_back
        while branch_frame is not None and branch_frame.f_code is not cone.ConeBlock.branch.__code__:
            branch_frame = branch_frame.f_back
        if branch_frame is None:
            raise ValueError('cone branch calling context unavailable')
        context = dict(source_line=branch_frame.f_lineno,
            iteration=branch_frame.f_locals.get('iteration'),
            half_local_may_be_stale=branch_frame.f_locals.get('half'),
            iteration_trace=copy.deepcopy(branch_frame.f_locals.get('trace',[])),
            call=len(response_calls), trial_call=trial_calls,
            minimum_normal_N=float(np.min(normal)), negative_normal_count=int(np.count_nonzero(normal<0)))
        response_calls.append(context)
        if np.any(normal<0) and first is None:
            if latest_trial is None:
                raise ValueError('no actual raw trial before negative normal')
            first = dict(candidate=latest_trial, tau=np.array(tau,copy=True), normal=np.array(normal,copy=True),
                         slip=np.array(slip,copy=True), context=context)
        # Retain original domain guard, exceptions and any original line-search retry.
        return actual_response(tau,normal,slip,mu,kappa,basis)

    b.trial = capture_trial
    caught = None
    with patch.object(cone,'cone_response',capture_response):
        try:
            b.branch(target['time_s'], released, maximum_iterations=14)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            caught = str(exc)
    after = state_signature(b)
    if before != after:
        raise ValueError('diagnostic mutated accepted state')
    if caught != expected_error or first is None:
        raise ValueError('failed to reproduce first negative-normal exception: '+str(caught))
    c = first['candidate']
    F, u = c['internal_force'], c['displacement']
    if not np.array_equal(first['normal'],F[b.nodes,0]):
        raise ValueError('captured normal not actual trial internal force')
    if not np.array_equal(first['tau'],-F[b.nodes,1:]):
        raise ValueError('captured tangent not actual trial internal force')
    gap = b.wall(c['time_s'])[0]-u[b.nodes,0]
    contact = np.zeros_like(b.X)
    contact[b.nodes,0] = -first['normal']
    contact[b.nodes,1:] = first['tau']
    support = np.zeros_like(b.X)
    support[b.fixed] = -(F+contact)[b.fixed]
    mechanical = F+contact+support
    ids = b.pa.particle_id[b.nodes].astype(int)
    negative = first['normal']<0
    arrays = {k:v for k,v in c.items() if isinstance(v,np.ndarray)}
    arrays.update(particle_id=b.pa.particle_id.copy(), contact_particle_id=ids,
        contact_feature_id=np.ones(len(ids),dtype=int), released=released,
        captured_normal_force_N=first['normal'], captured_tangent_force_global_yz_N=first['tau'],
        relative_increment_global_yz_m=first['slip'], enforced_normal_gap_m=gap,
        assigned_contact_force_N=contact, assigned_support_force_N=support,
        mechanical_balance_if_contact_forces_assigned_N=mechanical)
    np.savez_compressed(args.output/'REJECTED_first_negative_normal.npz',
        classification=np.array('REJECTED_DIAGNOSTIC_NOT_RESTART'), restart_allowed=np.array(False),
        accepted_time_s=np.array(b.time), accepted_step=np.array(b.count),
        candidate_time_s=np.array(c['time_s']), **arrays)
    if any(sha(args.run/name) != digest for name,digest in evidence.items()):
        raise ValueError('input artifacts changed')
    if any(sha(ROOT/name) != digest for name,digest in source.items()):
        raise ValueError('run source changed during diagnostic')
    result = dict(status='DIAGNOSTIC_COMPLETE_NOT_PATH_COMPLETE',
        classification='PURE_NEGATIVE_NORMAL_REPLAY_REJECTED_DIAGNOSTIC_NOT_RESTART',
        production_authorized=False,restart_allowed=False,error=caught,
        accepted_time_s=b.time,accepted_step=b.count,accepted_record_exact=True,
        accepted_state_unchanged=True,accepted_state_before=before,accepted_state_after=after,
        branch_invocations=1,maximum_iterations=14,target_time_s=target['time_s'],released=released.tolist(),
        trial_calls=trial_calls,response_calls=response_calls,first_negative_context=first['context'],
        contact_particle_id=ids.tolist(),normal_force_N=first['normal'].tolist(),
        accepted_normal_force_N=b.normal_force.tolist(),negative_particle_id=ids[negative].tolist(),
        negative_normal_force_N=first['normal'][negative].tolist(),feature_id=1,
        enforced_normal_gap_m=gap.tolist(),tangent_force_global_yz_N=first['tau'].tolist(),
        relative_increment_global_yz_m=first['slip'].tolist(),
        raw_free_mechanical_residual_relative=float(np.linalg.norm(mechanical.ravel()[b.free])/b.force_scale),
        residual_note='Assigned-contact mechanical balance only, NOT a cone root or accepted equilibrium.',
        interpretation='Nonconverged raw trial violates fixed-positive-normal domain; not proof of converged physical separation. No opening root, event, or accepted candidate solved.',
        input_sha256=evidence,source_sha256=source,
        diagnostic_sources_sha256={name:sha(ROOT/name) for name in
            ('scripts/diagnose_process_plastic_friction_onset.py','scripts/diagnose_process_plastic_friction_q7.py')},
        preregistration_sha256=sha(ROOT/PREREG),wall_seconds=time.perf_counter()-started,
        artifact_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    with (args.output/'diagnostic.json').open('x',encoding='utf-8') as f:
        f.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','error','accepted_time_s','accepted_step',
        'branch_invocations','trial_calls','first_negative_context','negative_particle_id',
        'negative_normal_force_N','enforced_normal_gap_m','accepted_state_unchanged','wall_seconds')},indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

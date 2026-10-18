"""Bounded read-only hidden-shell audit over actual exported checkpoint pairs.

No equilibrium/material/thermal evaluation and no new trajectory samples.
Each segment retains the frozen 256 finite-patch query limit. Stop on the first
FOUND/UNRESOLVED result; unsurveyed segments never inherit a passing status.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
import traceback

import numpy as np

import readback_discrete_process as reader
from batched_contact import BatchedChainRuleContact
from contact_shell_audit import CERTIFIED
from contact_shell_enclosure import audit_contact_shell_enclosure

ROOT = reader.ROOT
DEFAULT_READBACK = ROOT/'outputs/coupled_cost_20261004/batched/readback_21_875ms'
DEFAULT_OUTPUT = ROOT/'outputs/coupled_cost_20261004/chain_shell_21_875ms'


def write(path, value):
    pending = path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(reader.clean(value), indent=2, allow_nan=False), encoding='utf-8')
    pending.replace(path)


def sources():
    paths = [Path(__file__).resolve(), ROOT/'scripts/readback_discrete_process.py',
        ROOT/'scripts/batched_contact.py', ROOT/'scripts/contact_shell_audit.py',
        ROOT/'scripts/contact_shell_enclosure.py', ROOT/'tests/review/test_contact_shell_enclosure.py',
        ROOT/'src/afsd_backend/pysph/quasistatic_contact.py',
        ROOT/'src/afsd_backend/pysph/discrete_contact.py']
    return {p.relative_to(ROOT).as_posix(): reader.sha(p) for p in paths}


def worker(readback, output):
    reader.require(not output.exists(), 'refuse to overwrite existing audit evidence')
    output.mkdir(parents=True)
    begin = time.perf_counter()
    rows = []
    report = dict(status='RUNNING_READ_ONLY_GEOMETRY', physical_evaluations=0,
        physical_time_advanced_s=0., interpolation_performed=False,
        audited_segment_count=0, certified_segment_count=0, stopped_at_segment=None,
        geometry_budget_per_segment=256, analytic_path_budget_per_segment=256,
        internal_wall_limit_s=110., external_wall_limit_s=120.,
        scope='ENDPOINT_INACTIVE_PATCH_PAIRS_ON_ACTUAL_DECLARED_CHECKPOINT_SEGMENTS',
        limitation='Does not validate active-endpoint force accuracy, physical trajectory, or experiments.')
    def stop():
        try:
            write(output/'HARD_TIMEOUT.json', dict(status='UNRESOLVED_TOTAL_READ_ONLY_BUDGET',
                elapsed_wall_s=time.perf_counter()-begin, completed_segments=rows,
                unsurveyed_segments_are_not_certified=True, physical_evaluations=0))
        finally:
            os._exit(124)
    timer = threading.Timer(110., stop); timer.daemon=True; timer.start()
    code = 2
    try:
        frozen = sources(); report['source_sha256'] = frozen
        summary_path, trajectory_path = readback/'summary.json', readback/'accepted_trajectory.npz'
        summary = reader.read(summary_path)
        reader.require(summary['status'] == 'PASS_ACTUAL_CHAIN_READBACK_NOT_PHYSICAL_OR_EXPERIMENTAL_VALIDATION',
                       'required actual-chain readback audit absent')
        reader.require(summary['physical_evaluations'] == 0, 'unexpected readback generation mode')
        reader.require(reader.sha(trajectory_path) == summary['trajectory_sha256'], 'trajectory hash differs')
        reader.require(summary['readback_script_sha256'] == reader.sha(ROOT/'scripts/readback_batched_process.py'),
                       'readback implementation changed after evidence creation')
        for name, digest in summary['source_and_input_sha256'].items():
            reader.require(reader.sha(ROOT/name) == digest, 'bound readback source/input changed: '+name)
        final_file = ROOT/summary['final']['checkpoint']
        reader.require(reader.sha(final_file) == summary['final']['sha256'], 'final checkpoint hash differs')
        final_pa, final_state = reader.restore(final_file)
        p = final_state['parameters']; provenance = final_state['parameter_provenance']
        parameter_file = ROOT/provenance['parameter_source']
        reader.require(reader.sha(parameter_file) == provenance['source_sha256'], 'parameter source hash differs')
        bound_inputs = {path.relative_to(ROOT).as_posix(): reader.sha(path)
                        for path in (summary_path, trajectory_path, final_file, parameter_file)}
        with np.load(trajectory_path, allow_pickle=False) as data:
            reader.require(bool(data['actual_checkpoint_samples_only']) and not bool(data['interpolation_performed']),
                           'trajectory contains fabricated samples')
            times = data['time_s'].copy(); ids = data['particle_id'].copy()
            positions = data['position_m'].copy(); origins = data['tool_origin_m'].copy()
            rotations = data['tool_rotation'].copy()
            stored_ids = data['particle_property__particle_id'][..., 0]
            reader.require(np.all(stored_ids == ids[None, :]), 'persistent trajectory IDs changed')
        reader.require(len(times) == summary['samples'] == summary['accepted_steps']+1 and
                       len(summary['rows']) == summary['accepted_steps'], 'trajectory/index count mismatch')
        reader.require(positions.shape == (len(times), len(ids), 3) and origins.shape == (len(times), 3) and
                       rotations.shape == (len(times), 3, 3), 'trajectory dimensions differ')
        reader.require(np.isfinite(positions).all() and np.isfinite(times).all() and np.all(np.diff(times)>0),
                       'nonfinite or discontinuous exported positions/time')
        reader.require(times[0] == 0. and times[-1] == final_state['time_s'] == summary['accepted_time_s'],
                       'trajectory endpoints differ from accepted chain')
        reader.require(np.array_equal(positions[-1], reader.xyz(final_pa)) and
                       np.array_equal(ids, np.sort(final_pa.particle_id)), 'trajectory final IDs/positions differ')
        for number, time_s in enumerate(times):
            expected = reader.expected_pose(p, time_s)
            reader.require(np.allclose(origins[number], expected['origin'], atol=1e-14, rtol=0.) and
                           np.allclose(rotations[number], expected['rotation'], atol=1e-14, rtol=0.),
                           'exported tool pose differs from declared program')
        report.update(input_sha256=bound_inputs, parameters=p, parameter_provenance=provenance,
            total_segment_count=len(times)-1, trajectory_final_time_s=float(times[-1]),
            particle_count=len(ids), starting_time_s=float(times[0]))
        write(output/'manifest.json', report); write(output/'segments.json', rows)
        contact = BatchedChainRuleContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
            tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2,
            normal_stiffness=p['normal_stiffness'], tangential_stiffness=p['normal_stiffness'],
            friction=p['friction'], base_z=-p['substrate_depth'])
        for number in range(1, len(times)):
            if time.perf_counter()-begin >= 105.:
                report.update(status='UNRESOLVED_TOTAL_READ_ONLY_BUDGET', stopped_at_segment=number,
                    unexamined_segment_count=len(times)-number)
                code = 124
                break
            source_row = summary['rows'][number-1]
            reader.require(source_row['step'] == number and source_row['time_s'] == times[number] and
                           abs(source_row['dt_s']-(times[number]-times[number-1])) <= 1e-15,
                           'readback row/time/step mismatch')
            start = time.perf_counter()
            audit = audit_contact_shell_enclosure(contact, positions[number-1], positions[number],
                dict(origin=origins[number-1], rotation=rotations[number-1]),
                dict(origin=origins[number], rotation=rotations[number]))
            # The scientific algorithm uses row indices; publish real IDs next
            # to them without substituting identity into the geometric model.
            for key in ('hidden_contacts', 'planar_certificates', 'planar_unresolved_intervals'):
                for entry in audit.get(key, []):
                    entry['particle_id'] = int(ids[entry['particle_index']])
            for interval in audit.get('unresolved_intervals', []):
                interval['particle_id_feature_pairs'] = [[int(ids[i]), int(j)] for i, j in interval['pairs']]
            audit.update(global_step=number, time0_s=float(times[number-1]), time1_s=float(times[number]),
                accepted_checkpoint=source_row['checkpoint'], accepted_checkpoint_sha256=source_row['sha256'],
                physical_evaluations=0)
            path = output/f'segment_{number:03d}.json'; write(path, audit)
            row = dict(global_step=number, time0_s=float(times[number-1]), time1_s=float(times[number]),
                status=audit['status'], failure_kind=audit['failure_kind'],
                original_status=audit['original_audit']['status'], original_failure_kind=audit['original_failure_kind'],
                geometry_queries=audit['geometry_queries'], additional_finite_patch_queries=audit['additional_finite_patch_queries'],
                plane_path_evaluations=audit['plane_path_evaluations'], planar_certificates=len(audit['planar_certificates']),
                hidden_contacts=audit['hidden_contacts'], planar_unresolved_intervals=audit['planar_unresolved_intervals'],
                result_file=path.relative_to(ROOT).as_posix(), result_sha256=reader.sha(path), wall_seconds=time.perf_counter()-start)
            rows.append(row); write(output/'segments.json', rows)
            report.update(audited_segment_count=len(rows), certified_segment_count=sum(r['status']==CERTIFIED for r in rows),
                last_audited_time_s=float(times[number]), unexamined_segment_count=len(times)-1-len(rows))
            write(output/'manifest.json', report)
            print(json.dumps({key: row[key] for key in ('global_step', 'time1_s', 'status', 'geometry_queries', 'plane_path_evaluations', 'wall_seconds')}), flush=True)
            if audit['status'] != CERTIFIED:
                report.update(status=audit['status'], stopped_at_segment=number, failure_kind=audit['failure_kind'],
                    first_failed_result=row['result_file'])
                code = 2
                break
        else:
            report['status'] = CERTIFIED
            code = 0
        for name, digest in dict(frozen, **bound_inputs).items():
            reader.require(reader.sha(ROOT/name) == digest, 'source/input drift during geometry audit: '+name)
        report['segments_index_sha256'] = reader.sha(output/'segments.json')
    except Exception as error:
        traceback.print_exc()
        report.update(status='FAILED_READ_ONLY_AUDIT', error=repr(error))
        code = 2
    finally:
        timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,
        unsurveyed_segments_are_not_certified=True,
        no_new_physical_time_or_acceptance=True)
    write(output/'manifest.json', report)
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--readback', type=Path, default=DEFAULT_READBACK)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.readback, args.output = args.readback.resolve(), args.output.resolve()
    reader.require(args.readback.is_relative_to(ROOT) and args.output.is_relative_to(ROOT), 'evidence path outside project')
    if args.worker:
        raise SystemExit(worker(args.readback, args.output))
    log_path = args.output.with_name(args.output.name+'_supervision.log')
    receipt_path = args.output.with_name(args.output.name+'_supervision.json')
    reader.require(not args.output.exists() and not log_path.exists() and not receipt_path.exists(), 'refuse evidence overwrite')
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--readback', str(args.readback), '--output', str(args.output), '--worker']
    begin = time.perf_counter()
    with log_path.open('w', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                timeout=120., env=dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:
            code = 124
    write(receipt_path, dict(command=command, exit_code=code, wall_seconds=time.perf_counter()-begin,
        external_wall_limit_s=120., physical_evaluations=0, physical_time_advanced_s=0.))
    print(json.dumps(dict(exit_code=code, output=str(args.output), log=str(log_path))))
    raise SystemExit(code)

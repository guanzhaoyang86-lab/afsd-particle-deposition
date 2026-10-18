"""Read-only PID1 differentiation diagnosis; never advances a physical state.

PREREGISTRATION (before execution): use saved failed PID1, force axis x only,
unchanged physical inputs and old kick defaults, h=h0/2**k for ALL k=0..5.
No stop-on-green and no source/acceptance changes. Record raw probes, force
steps, all B/C forward/backward/central estimates and convergence ratios.
Compare Richardson one-sided estimates solely as a mathematical diagnostic.
The surrounding command must enforce timeout=60s and save a fresh log.
"""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph import event_response as er
from afsd_backend.pysph.contact_kick import kick


def read(path):
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        return json.load(stream)


def emit(kind, value):
    print(kind+' '+json.dumps(er._plain(value), sort_keys=True, allow_nan=False), flush=True)


def main():
    start = time.perf_counter()
    failed_path = ROOT/'outputs/process_event_20261003/response_probe/particle_000001.json.gz'
    anchor_path = ROOT/'outputs/process_global_20261003/root_attempt/latest_trial.json.gz'
    source_path = ROOT/'src/afsd_backend/pysph/event_response.py'
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (failed_path, anchor_path, source_path)}
    failed = read(failed_path)
    journal = read(anchor_path)['evaluation']['boundary']['traces'][0]
    assert journal['particle_id'] == failed['particle_id'] == 1
    with (ROOT/'outputs/process_event_20261003/response_probe/resolved_config.json').open(encoding='utf-8') as stream:
        mu = json.load(stream)['parameters']['friction']
    faces = journal['contact_rows']; dt = journal['dt_s']; mass = journal['mass_kg']
    normals = np.array([f['normal'] for f in faces])
    fn = np.array([f['normal_force_N'] for f in faces])
    walls = np.array([f['wall_velocity'] for f in faces])
    labels = [f['contact_id'] for f in faces]
    force = np.array(journal['internal_force_N'])
    h0 = failed['force_step_N'][0]
    scale = dt/mass
    emit('PREREGISTERED', dict(particle_id=1, axis=0, force_steps_N=[h0/2**k for k in range(6)],
         all_steps_run_without_adaptive_selection=True, accepted=False, hashes=before,
         initial_failure=failed['reason'], initial_column=failed['column_evidence'][0]))
    calls = 0; kick_time = 0.

    def evaluate(delta):
        nonlocal calls, kick_time
        shifted = force.copy(); shifted[0] += delta
        tick = time.perf_counter()
        physical = kick(journal['material_velocity_before_m_s'], mass, shifted, normals,
            fn, walls, mu, dt, free_axes=journal['free_axes'],
            prescribed_velocity=journal['prescribed_velocity'])
        calls += 1; kick_time += time.perf_counter()-tick
        probe = er._summary(physical, dt, normals, mu*fn, walls, labels)
        endpoint_values = np.array([a[key] for a in physical['intervals'] for key in ('v_before', 'v_after')])
        slip = endpoint_values[:, None, :]-walls[None, :, :]
        slip -= np.sum(slip*normals, axis=2)[:, :, None]*normals
        probe.update(delta_force_N=delta, minimum_all_interval_endpoint_slip_m_s=float(
            np.min(np.linalg.norm(slip, axis=2))),
            interval_kind_sequence=[a['kind'] for a in physical['intervals']],
            interval_dt_s=[a['dt_s'] for a in physical['intervals']])
        emit('RAW_LOCAL_PROBE', probe)
        return probe, np.stack([probe['mean_velocity_m_s'], probe['velocity_m_s']])

    anchor, value0 = evaluate(0.)
    assert anchor['result_sha256'] == failed['baseline']['result_sha256']
    rows = []
    for k in range(6):
        h = h0/2**k
        plus, pv = evaluate(h); minus, mv = evaluate(-h)
        forward = (pv-value0)/h; backward = (value0-mv)/h; central = .5*(forward+backward)
        gap = np.abs(forward-backward)
        bound = .002*np.maximum(np.abs(forward), np.abs(backward))+2e-5*scale
        branch = [er._branch_change(anchor, p, labels) for p in (plus, minus)]
        partitions = np.array([plus['interval_dt_s'], minus['interval_dt_s']])
        row = dict(k=k, h_N=h, forward_s_kg=forward, backward_s_kg=backward,
            central_s_kg=central, one_sided_gap_s_kg=gap, original_bound_s_kg=bound,
            original_gap_test_pass=bool(np.all(gap <= bound)), physical_branch_changes=branch,
            maximum_normalized_partition_change=float(np.max(np.abs(
                partitions-np.array(anchor['interval_dt_s']))))/dt,
            minimum_slip_m_s=min(p['minimum_all_interval_endpoint_slip_m_s'] for p in (plus, minus)))
        if rows:
            previous = rows[-1]
            rich_forward = 2*forward-previous['forward_s_kg']
            rich_backward = 2*backward-previous['backward_s_kg']
            rich_bound = .002*np.maximum(np.abs(rich_forward), np.abs(rich_backward))+2e-5*scale
            row.update(central_change_s_kg=central-previous['central_s_kg'],
                gap_ratio_Cyx=float(previous['one_sided_gap_s_kg'][1, 1]/gap[1, 1]),
                gap_ratio_Byx=float(previous['one_sided_gap_s_kg'][0, 1]/gap[0, 1]),
                richardson_forward_s_kg=rich_forward, richardson_backward_s_kg=rich_backward,
                richardson_gap_s_kg=np.abs(rich_forward-rich_backward),
                richardson_gap_pass=bool(np.all(np.abs(rich_forward-rich_backward) <= rich_bound)))
        rows.append(row)
        emit('DERIVATIVE_LADDER', row)
    after = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in (failed_path, anchor_path, source_path)}
    assert before == after
    emit('SUMMARY', dict(status='READ_ONLY_DIAGNOSIS_NOT_A_TANGENT_RECEIPT', accepted=False,
        source_and_initial_failure_preserved=True, kick_calls=calls, kick_seconds=kick_time,
        wall_seconds=time.perf_counter()-start, all_force_levels_run=6,
        baseline_minimum_slip_m_s=anchor['minimum_all_interval_endpoint_slip_m_s'],
        original_response_still_unavailable=True))


if __name__ == '__main__':
    main()

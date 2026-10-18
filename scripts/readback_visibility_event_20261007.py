"""Read a rejected saved trial; locate visibility events without rerunning it."""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import readback_discrete_process as reader
from visible_support_events_20261007 import audit_visibility_events

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT/'outputs/extended_tight_path_20261007/window_011/coarse/step_011'
DEFAULT_OUTPUT = ROOT/'outputs/visibility_events_20261007/window011_step011_v1.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(source, output):
    if output.exists():
        raise ValueError('refuse to overwrite event evidence')
    source = source.resolve()
    inputs = [source/'visible_graph.npz', source/'visible_graph.json',
              source/'NOT_RESTARTABLE_candidate.npz', source/'visible_path_audit.json']
    before = {str(path): sha(path) for path in inputs}
    old = json.loads(inputs[3].read_text(encoding='utf-8'))
    meta = json.loads(inputs[1].read_text(encoding='utf-8'))
    state = reader.load(str(inputs[2]))['solver_data']
    if state.get('restart_allowed') is not False:
        raise ValueError('this readback is specifically for an unaccepted candidate')
    p = state['parameters']
    t0, t1 = float(meta['parent_time_s']), float(meta['endpoint_time_s'])
    started = time.perf_counter()
    with np.load(inputs[0], allow_pickle=False) as saved:
        # Diagnose the exact retained graph rejected by the frozen path audit.
        events = audit_visibility_events(saved['position0'], saved['position1'],
            saved['i'], saved['j'], reader.expected_pose(p, t0), reader.expected_pose(p, t1),
            p, saved['particle_id'])
        if events['input_sha256']['world_position0'] != old['input_sha256']['world_position0']:
            raise ValueError('saved rejected path start fingerprint mismatch')
        if events['input_sha256']['world_position1'] != old['input_sha256']['world_position1']:
            raise ValueError('saved rejected path end fingerprint mismatch')
        if events['input_sha256']['supplied_topology'] != old['input_sha256']['retained_topology']:
            raise ValueError('saved rejected support topology fingerprint mismatch')
    for row in events['transition_brackets']:
        row['start_time_s'] = t0+row['start_fraction']*(t1-t0)
        row['end_time_s'] = t0+row['end_fraction']*(t1-t0)
    after = {str(path): sha(path) for path in inputs}
    if before != after:
        raise ValueError('old evidence changed during readback')
    code = [Path(__file__), ROOT/'scripts/visible_support_events_20261007.py',
            ROOT/'scripts/visible_support.py', ROOT/'scripts/visible_support_numerics.py',
            ROOT/'scripts/visible_support_path_tight.py', ROOT/'scripts/audit_solid_occluded_neighbors.py',
            ROOT/'src/afsd_backend/pysph/quasistatic_contact.py']
    result = dict(status='READBACK_COMPLETE_ORIGINAL_CANDIDATE_STILL_REJECTED',
        classification='SAVED_CANDIDATE_GEOMETRY_DIAGNOSTIC_ONLY',
        original_status=old['status'], original_failure_kind=old['failure_kind'],
        old_evidence_unchanged=True, acceptance_changed=False, restart_allowed=False,
        physical_solver_evaluations=0, nnps_evaluations=0, new_accepted_physical_time_s=0.,
        start_time_s=t0, end_time_s=t1, old_found_edges=old['found_edges'], events=events,
        input_sha256=before, source_sha256={str(path): sha(path) for path in code},
        wall_seconds=time.perf_counter()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
    print(json.dumps(dict(output=str(output), status=result['status'], event_status=events['status'],
        physical_solver_evaluations=0, wall_seconds=result['wall_seconds'],
        endpoint_changed_edges=events['endpoint_changed_edges'],
        transition_brackets=events['transition_brackets']), ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.source, args.output)

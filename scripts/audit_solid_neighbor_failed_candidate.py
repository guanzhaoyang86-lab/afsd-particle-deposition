"""Read raw rejected candidate data for geometry only; never restore/restart it."""
from pathlib import Path
import json
import time
import numpy as np
import readback_discrete_process as reader
from audit_solid_occluded_neighbors import audit_solid_edge_visibility
from audit_solid_neighbor_snapshots import sorted_edges
from afsd_backend.pysph.mechanics import pairs

ROOT = reader.ROOT
FOLDER = ROOT/'outputs/coupled_cost_20261004/batched/shoulder_window_22_30ms/coarse/step_019'
OUT = ROOT/'outputs/coupled_cost_20261004/solid_occluded_neighbors_unaccepted26_75ms.json'


def main():
    begin = time.perf_counter()
    reader.require(not OUT.exists(), 'refuse to overwrite failed-candidate geometry audit')
    receipt_path = ROOT/'outputs/coupled_cost_20261004/solid_edge_visibility_receipt.json'
    receipt = reader.read(receipt_path)
    guard_path = ROOT/'scripts/audit_solid_occluded_neighbors.py'
    reader.require(reader.sha(guard_path) == receipt['source_sha256'][guard_path.relative_to(ROOT).as_posix()], 'frozen guard changed')
    meta = reader.read(FOLDER.parent/'manifest.json')
    reader.require(meta['status'] == 'FAILED_NO_NEW_ACCEPTED_STATE' and meta['accepted_steps'] == 0, 'unexpected old failed branch status')
    raw = reader.load(str(FOLDER/'NOT_RESTARTABLE_candidate.npz'))
    pa, state = raw['arrays']['solid'], raw['solver_data']
    reader.require(state['restart_allowed'] is False, 'candidate unexpectedly restartable')
    reader.require(abs(state['time_s']-.02675) < 1e-15, 'unexpected candidate time')
    with np.load(FOLDER/'FULL_UNCOMMITTED_candidate.npz', allow_pickle=False) as candidate:
        reader.require(not bool(candidate['accepted']) and not bool(candidate['restart_allowed']), 'candidate publication flag changed')
        x = candidate['position_m'].copy()
        reader.require(np.array_equal(candidate['particle_id'], pa.particle_id), 'raw candidate IDs differ')
        reader.require(np.array_equal(x, np.column_stack((pa.x, pa.y, pa.z))), 'raw candidate positions differ')
    with np.load(FOLDER/'candidate_evidence.npz', allow_pickle=False) as evidence:
        reader.require(np.array_equal(evidence['particle_id'], pa.particle_id), 'graph raw IDs differ')
        i, j = evidence['i'].copy(), evidence['j'].copy()
    rebuilt_i, rebuilt_j = pairs(pa)
    reader.require(np.array_equal(sorted_edges(i, j), sorted_edges(rebuilt_i, rebuilt_j)), 'saved actual candidate NNPS graph differs from read-only rebuild')
    pose = reader.expected_pose(state['parameters'], state['time_s'])
    result = audit_solid_edge_visibility(x, i, j, pose, state['parameters'])
    for row in result['occluded_edges']:
        row.update(particle_id_i=int(pa.particle_id[row['i']]), particle_id_j=int(pa.particle_id[row['j']]),
            component_id_i=int(pa.component_id[row['i']]), component_id_j=int(pa.component_id[row['j']]))
    paths = [Path(__file__), guard_path, receipt_path, FOLDER.parent/'manifest.json',
        FOLDER/'NOT_RESTARTABLE_candidate.npz', FOLDER/'FULL_UNCOMMITTED_candidate.npz', FOLDER/'candidate_evidence.npz']
    record = dict(status='DIAGNOSTIC_UNACCEPTED_CANDIDATE_ONLY',
        accepted=False, restart_allowed=False, official_accepted_time_s=meta['accepted_time_s'],
        candidate_time_s=state['time_s'], candidate_nominal_step=state['step'],
        saved_actual_graph_exactly_matches_rebuilt=True, raw_candidate_snapshot_matches_fields=True,
        physical_solver_evaluations=0, real_pysph_nnps_rebuilds=1,
        pose_from_saved_parameters_and_candidate_time=reader.clean(pose), endpoint=result,
        source_sha256={path.relative_to(ROOT).as_posix(): reader.sha(path) for path in paths},
        qualification='This geometry result does not overturn the old FAIL, publish the candidate, or extend the 22 ms accepted chain.',
        wall_seconds=time.perf_counter()-begin)
    OUT.write_text(json.dumps(reader.clean(record), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(output=str(OUT), candidate_time_s=state['time_s'], edge_count=result['edge_count'],
        occluded_edge_count=result['occluded_edge_count'], max_inside_length_m=result['max_occluded_interval_length_m'],
        wall_seconds=record['wall_seconds'])))


if __name__ == '__main__':
    main()

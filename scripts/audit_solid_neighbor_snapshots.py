"""Read-only three-snapshot NNPS visibility audit; never evaluates a trial."""
from pathlib import Path
import json
import time
import numpy as np
import readback_discrete_process as reader
from audit_solid_occluded_neighbors import audit_solid_edge_visibility
from afsd_backend.pysph.mechanics import pairs

ROOT = reader.ROOT
OUT = ROOT/'outputs/coupled_cost_20261004/solid_occluded_neighbors_audit.json'
OLD = ROOT/'outputs/discrete_contact_process_20261004/v3'
NEW = ROOT/'outputs/coupled_cost_20261004/batched'


def sorted_edges(i, j):
    value = np.column_stack((i, j))
    return value[np.lexsort((value[:, 1], value[:, 0]))]


def main():
    begin = time.perf_counter()
    reader.require(not OUT.exists(), 'refuse to overwrite independent audit')
    final_folder = NEW/'first_shoulder_enclosed/halves/step_002'
    final_pa, final_state = reader.restore(final_folder/'accepted.npz')
    p = final_state['parameters']
    physics_paths = [ROOT/'src/afsd_backend/pysph'/name for name in (
        'mechanics.py', 'kernel.py', 'quasistatic_process.py', 'implicit_thermal.py')]
    for path in physics_paths:
        name = path.relative_to(ROOT).as_posix()
        reader.require(final_state['source_sha256'][name] == reader.sha(path), 'frozen physics source mismatch: '+name)
    receipt_file = ROOT/'outputs/coupled_cost_20261004/solid_edge_visibility_receipt.json'
    receipt = reader.read(receipt_file)
    for name, digest in receipt['source_sha256'].items():
        reader.require(reader.sha(ROOT/name) == digest, 'guard receipt source changed: '+name)
    source_paths = [Path(__file__), ROOT/'scripts/audit_solid_occluded_neighbors.py',
        ROOT/'scripts/check_solid_edge_visibility.py', receipt_file,
        ROOT/'scripts/run_batched_process.py']+physics_paths
    rows = []
    for label, checkpoint, evidence_path in (
        ('independent_initial', OLD/'quarters/independent_initial.npz', None),
        ('qualified_5ms', OLD/'post_start/pilot_refinement/step_016/accepted.npz', None),
        ('qualified_22ms', final_folder/'accepted.npz', final_folder/'candidate_evidence.npz')):
        pa, state = reader.restore(checkpoint)
        x = np.column_stack((pa.x, pa.y, pa.z))
        source_paths.append(checkpoint)
        reader.require(np.array_equal(pa.particle_id, final_pa.particle_id), 'snapshot persistent raw ID ordering differs')
        reader.require(np.array_equal(pa.component_id, final_pa.component_id), 'snapshot components differ')
        reader.require(np.array_equal(pa.m, final_pa.m), 'snapshot masses differ')
        if 'parameters' in state:
            reader.require(reader.clean(state['parameters']) == reader.clean(p), 'physical parameters changed')
        expected_pose = reader.expected_pose(p, state['time_s'])
        pose_path = checkpoint.parent/'tool_pose.json'
        if label != 'independent_initial':
            saved_pose = reader.read(pose_path)
            for key in expected_pose:
                reader.require(np.allclose(expected_pose[key], saved_pose[key], rtol=0., atol=1e-15), 'saved tool pose mismatch')
            pose = saved_pose
            source_paths.append(pose_path)
        else:
            pose = expected_pose
        rebuild_begin = time.perf_counter()
        rebuilt_i, rebuilt_j = pairs(pa)
        rebuild_wall = time.perf_counter()-rebuild_begin
        saved_match = None
        if evidence_path is not None:
            with np.load(evidence_path, allow_pickle=False) as evidence:
                reader.require(np.array_equal(evidence['particle_id'], pa.particle_id), 'saved graph ID order differs')
                i, j = evidence['i'].copy(), evidence['j'].copy()
            saved_match = np.array_equal(sorted_edges(i, j), sorted_edges(rebuilt_i, rebuilt_j))
            reader.require(saved_match, 'real rebuilt NNPS differs from saved candidate graph')
            source_paths.append(evidence_path)
            candidate_path = checkpoint.parent/'FULL_UNCOMMITTED_candidate.npz'
            with np.load(candidate_path, allow_pickle=False) as candidate:
                reader.require(np.array_equal(candidate['particle_id'], pa.particle_id), 'candidate raw IDs differ')
                reader.require(np.array_equal(candidate['position_m'], x), 'candidate x differs from accepted endpoint')
            source_paths.append(candidate_path)
        else:
            i, j = rebuilt_i, rebuilt_j
        audit_begin = time.perf_counter()
        result = audit_solid_edge_visibility(x, i, j, pose, p)
        audit_wall = time.perf_counter()-audit_begin
        radius = np.linalg.norm(x[j]-x[i], axis=1)
        q = radius/(.5*(pa.h[i]+pa.h[j]))
        participant_ids, participant_components = set(), set()
        for edge in result['occluded_edges']:
            left, right = edge['i'], edge['j']
            edge.update(particle_id_i=int(pa.particle_id[left]), particle_id_j=int(pa.particle_id[right]),
                component_id_i=int(pa.component_id[left]), component_id_j=int(pa.component_id[right]),
                kernel_support_ratio_q=float(q[edge['edge_index']]),
                raw_kernel_gradient_nonzero=bool(0. < q[edge['edge_index']] < 2.))
            participant_ids.update((int(pa.particle_id[left]), int(pa.particle_id[right])))
            participant_components.update((int(pa.component_id[left]), int(pa.component_id[right])))
        components, counts = np.unique(np.sort(np.column_stack((pa.component_id[i], pa.component_id[j])), axis=1), axis=0, return_counts=True)
        result.update(snapshot=label, checkpoint=checkpoint.relative_to(ROOT).as_posix(),
            accepted_time_s=state['time_s'], accepted_step=state['step'], particle_count=len(x),
            edge_source='SAVED_FULL_CANDIDATE_EVIDENCE' if evidence_path else 'READ_ONLY_REAL_PYSPH_NNPS_REBUILD',
            saved_actual_graph_exactly_matches_rebuilt=saved_match,
            kernel_nonzero_edge_count=int(np.count_nonzero((q > 0.) & (q < 2.))),
            edge_counts_by_component_pair=[dict(component_pair=pair.tolist(), count=int(count)) for pair, count in zip(components, counts)],
            occluded_participant_particle_ids=sorted(participant_ids),
            occluded_participant_component_ids=sorted(participant_components),
            worst_occluded_edge=max(result['occluded_edges'], key=lambda row: row['maximum_contiguous_inside_length_m'], default=None),
            tool_pose=reader.clean(pose), nnps_rebuild_wall_seconds=rebuild_wall,
            visibility_guard_wall_seconds=audit_wall)
        rows.append(result)
    source_hashes = {path.relative_to(ROOT).as_posix(): reader.sha(path) for path in source_paths}
    result = dict(status='PASS_SNAPSHOT_ENDPOINT_VISIBILITY_ONLY' if all(row['occluded_edge_count'] == 0 for row in rows) else 'FAIL_SOLID_OCCLUDED_NNPS_EDGE',
        physical_solver_evaluations=0, constitutive_or_force_evaluations=0, thermal_solves=0,
        real_pysph_nnps_rebuilds=3, snapshots=rows, geometry_fixture_count=receipt['fixture_count'],
        source_sha256=source_hashes, physical_dimensions_m={name: p[name] for name in ('layer', 'tool_top', 'tool_radius', 'bore_width', 'dx')},
        particle_smoothing_length_range_m=[float(final_pa.h.min()), float(final_pa.h.max())],
        source_participation_audit=dict(visibility_filter_present=False,
            endpoint_graph='quasistatic_process.py:93 operators(cp) uses candidate x1.',
            correction_gradient='quasistatic_process.py:46-54 accumulates every supplied i,j into moment/raw displacement gradient.',
            internal_force='quasistatic_process.py:138 passes identical i,j,g to mechanics.py:81-98 paired force.',
            heat_conduction='quasistatic_process.py:122-128 passes identical i,j,g to implicit_thermal.py:15-50 conductance.',
            batched_variant='scripts/run_batched_process.py:33-40 replaces contact only; no graph visibility mask.',
            qualification='All current NNPS edges with nonzero raw kernel gradient participate in gradient, force, and heat. Actual net pair force may cancel; no solver was evaluated.'),
        limitations=['These three snapshots do not certify every past interval or any future geometry.',
            'No occlusion filtering is implemented in mechanics or heat; a future occluded edge is a stop condition, not silently removed.',
            'Endpoint NNPS graph visibility does not establish free-surface SPH consistency or converged material transport.',
            'This geometry audit is not AFSD physics validation or paper reproduction.'],
        wall_seconds=time.perf_counter()-begin)
    OUT.write_text(json.dumps(reader.clean(result), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(status=result['status'], output=str(OUT),
        snapshots=[dict(name=row['snapshot'], edges=row['edge_count'], occluded=row['occluded_edge_count'],
            saved_graph_match=row['saved_actual_graph_exactly_matches_rebuilt']) for row in rows], wall_seconds=result['wall_seconds'])))


if __name__ == '__main__':
    main()

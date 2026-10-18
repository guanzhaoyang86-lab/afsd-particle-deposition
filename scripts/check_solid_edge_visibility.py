"""Geometry-only fixtures and saved 22 ms endpoint guard cost; no physics."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from audit_solid_occluded_neighbors import audit_solid_edge_visibility

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/coupled_cost_20261004/solid_edge_visibility_receipt.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixtures():
    p = dict(tool_radius=3., bore_width=2., layer=0., tool_top=2.)
    pose = dict(origin=np.zeros(3), rotation=np.eye(3))
    cases = [
        ('inside_open_bore', [[.2, 0., .5], [.8, 0., 1.5]], 0.),
        ('cross_tool_wall', [[0., 0., 1.], [4., 0., 1.]], 2.),
        ('touch_shoulder_plane_only', [[.5, 0., 0.], [2.5, 0., 0.]], 0.),
        ('around_exact_aperture_lip', [[1.5, 0., -1.], [.5, 0., 1.]], 0.),
        ('clip_aperture_lip', [[1.5, 0., -1.], [.75, 0., 1.]], np.sqrt(.75**2+2.**2)/6.),
        ('touch_cylinder_only', [[-2., 3., 1.], [2., 3., 1.]], 0.),
        ('along_bore_wall_only', [[1., -.5, .5], [1., .5, 1.5]], 0.),
        ('through_two_walls', [[-4., 0., 1.], [4., 0., 1.]], 2.)]
    rows = []
    for name, xyz, expected in cases:
        result = audit_solid_edge_visibility(np.array(xyz), np.array([0]), np.array([1]), pose, p)
        assert result['occluded_edge_count'] == int(expected > 1e-9), name
        assert abs(result['max_occluded_interval_length_m']-expected) <= 2e-14, name
        rows.append(dict(name=name, expected_maximum_length_m=float(expected),
            actual_maximum_length_m=result['max_occluded_interval_length_m'], status='PASS_GEOMETRY_FIXTURE'))
    angle = .37
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.], [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
    translation = np.array([.71, -.22, .13])
    result = audit_solid_edge_visibility(np.array(cases[1][1])@rotation.T+translation,
        np.array([0]), np.array([1]), dict(origin=translation, rotation=rotation), p)
    assert result['occluded_edge_count'] == 1 and abs(result['max_occluded_interval_length_m']-2.) < 2e-14
    rows.append(dict(name='rigid_pose_covariance', status='PASS_GEOMETRY_FIXTURE'))
    # Threshold is physical contiguous length, not root-parameter fraction.
    for length, expected in ((.5e-9, 0), (1.5e-9, 1)):
        result = audit_solid_edge_visibility(np.array([[1.-length, 0., 1.], [1.+length, 0., 1.]]),
            np.array([0]), np.array([1]), pose, p)
        assert result['occluded_edge_count'] == expected
        rows.append(dict(name='length_threshold_'+str(length), status='PASS_GEOMETRY_FIXTURE'))
    return rows


def main():
    assert not OUT.exists(), 'refuse to overwrite geometry receipt'
    fixture_rows = fixtures()
    import readback_discrete_process as reader
    folder = ROOT/'outputs/coupled_cost_20261004/batched/first_shoulder_enclosed/halves/step_002'
    pa, state = reader.restore(folder/'accepted.npz')
    evidence = np.load(folder/'candidate_evidence.npz', allow_pickle=False)
    assert np.array_equal(pa.particle_id, evidence['particle_id'])
    x = np.column_stack((pa.x, pa.y, pa.z))
    candidate = np.load(folder/'FULL_UNCOMMITTED_candidate.npz', allow_pickle=False)
    assert np.array_equal(x, candidate['position_m'])
    pose = reader.read(folder/'tool_pose.json')
    times = []
    for _ in range(7):
        begin = time.perf_counter()
        result = audit_solid_edge_visibility(x, evidence['i'], evidence['j'], pose, state['parameters'])
        times.append(time.perf_counter()-begin)
    paths = [Path(__file__), ROOT/'scripts/audit_solid_occluded_neighbors.py', folder/'accepted.npz',
        folder/'candidate_evidence.npz', folder/'FULL_UNCOMMITTED_candidate.npz', folder/'tool_pose.json']
    receipt = dict(status='PASS_GEOMETRY_FIXTURES_AND_22MS_ENDPOINT_VISIBILITY',
        physical_solver_evaluations=0, fixtures=fixture_rows, fixture_count=len(fixture_rows),
        accepted_time_s=state['time_s'], accepted_steps=state['step'],
        raw_particle_id_order_verified=True, saved_candidate_position_matches_accepted=True,
        endpoint=result, guard_wall_seconds=times, maximum_guard_wall_seconds=max(times),
        median_guard_wall_seconds=float(np.median(times)),
        source_sha256={str(path.relative_to(ROOT)).replace('\\', '/'): sha(path) for path in paths})
    OUT.write_text(json.dumps(receipt, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(receipt=str(OUT), fixtures=len(fixture_rows), endpoint_status=result['status'],
        edges=result['edge_count'], occluded=result['occluded_edge_count'], guard_max_s=max(times),
        module_sha256=sha(ROOT/'scripts/audit_solid_occluded_neighbors.py'))))


if __name__ == '__main__':
    main()

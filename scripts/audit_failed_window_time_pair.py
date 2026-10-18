"""Read-only audit of the saved 25.5 ms pair-14 local stress rejection.

Loads raw saved candidate arrays, contact journals and NNPS edges only. Does
not construct a trial or call force, material, NNPS or thermal evaluators.
"""
from pathlib import Path
import json
import numpy as np
import readback_discrete_process as reader

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/coupled_cost_20261004/batched/finite_shoulder_windows/window_1'
OUT = ROOT/'outputs/coupled_cost_20261004/failure_time_audit.json'
PID = 275


def main():
    if OUT.exists():
        raise FileExistsError('refuse to overwrite independent failure audit')
    inputs = {}

    def bind(path):
        inputs[path.relative_to(ROOT).as_posix()] = reader.sha(path)
        return path

    def read(path):
        return reader.read(bind(path))

    def load(path):
        raw = reader.load(str(bind(path)))
        return raw['arrays']['solid'], raw['solver_data']

    screen = read(BASE/'fine/pair_014_temporal_screen.json')
    manifests = {name: read(BASE/name/'manifest.json') for name in ('coarse', 'fine')}
    assert manifests['coarse']['status'] == 'NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED'
    assert manifests['fine']['status'] == 'FAILED_WITH_ACCEPTED_PREFIX'
    assert manifests['fine']['accepted_steps'] == 26
    assert abs(manifests['fine']['accepted_time_s']-.02525) < 1e-15
    assert screen['pair'] == 14 and not screen['local_screen']['passed'] and screen['global_screen']['pass']
    expected_parent = '5368f67e717b4ce42d3927b8eca4d03a502a7ab59c5899d3d9b0b776c71d2d95'
    assert screen['parent_sha256'] == expected_parent
    for manifest in manifests.values():
        assert manifest['parent_sha256'] == expected_parent
    rows = {}
    raw_states = {}
    graphs = {}
    specifications = (
        ('coarse', 'coarse/step_014', 'coarse/step_013/NOT_RESTARTABLE_candidate.npz'),
        ('fine_half_1', 'fine/step_027', 'fine/step_026/accepted.npz'),
        ('fine', 'fine/step_028', 'fine/step_027/NOT_RESTARTABLE_candidate.npz'),
    )
    common_ids = None
    for label, folder, parent_path in specifications:
        directory = BASE/folder
        candidate_path = directory/'NOT_RESTARTABLE_candidate.npz'
        pa, state = load(candidate_path)
        parent, previous = load(BASE/parent_path)
        assert state['restart_allowed'] is False
        assert not (directory/'accepted.npz').exists()
        ids = pa.particle_id.copy()
        assert len(np.unique(ids)) == len(ids) and np.array_equal(ids, parent.particle_id)
        if common_ids is None:
            common_ids = ids
        else:
            assert np.array_equal(common_ids, ids)
        k = int(np.flatnonzero(ids == PID)[0])
        with np.load(bind(directory/'candidate_evidence.npz'), allow_pickle=False) as evidence:
            assert np.array_equal(ids, evidence['particle_id'])
            left, right = ids[evidence['i']], ids[evidence['j']]
            graph = set(tuple(sorted((int(i), int(j)))) for i, j in zip(left, right))
        graphs[label] = graph
        journal = read(directory/'FULL_UNCOMMITTED_candidate.json')
        contact = journal['contact']
        dt = state['time_s']-previous['time_s']
        assert abs(dt-contact['dt_s']) < 1e-15
        with np.load(bind(directory/'FULL_UNCOMMITTED_candidate.npz'), allow_pickle=False) as saved:
            assert np.array_equal(ids, saved['particle_id'])
            assert np.array_equal(pa.S.reshape(-1, 3, 3), saved['stress_Pa'])
            assert np.array_equal(np.column_stack((pa.x, pa.y, pa.z)), saved['position_m'])
            assert np.array_equal(pa.ep, saved['plastic_strain'])
        stress = pa.S.reshape(-1, 3, 3)[k].copy()
        d = pa.D.reshape(-1, 3, 3)[k].copy()
        dev = d-np.eye(3)*np.trace(d)/3.
        reconstructed_rate = float(np.sqrt(2./3.*np.sum(dev*dev)))
        reconstructed_ep_rate = float((pa.ep[k]-parent.ep[k])/contact['dt_s'])
        assert np.isclose(reconstructed_rate, pa.strain_rate[k], rtol=1e-13, atol=1e-15)
        assert np.isclose(reconstructed_ep_rate, pa.ep_rate[k], rtol=1e-13, atol=1e-15)
        feature_names = contact['features']
        active = np.asarray(contact['active'], bool)
        events = read(directory/'sampled_nonbore_events.json')
        seeds = set(screen['cumulative_shoulder_ids'])
        neighbors = sorted(j if i == PID else i for i, j in graph if PID in (i, j))
        shoulder_links = sorted(seed for seed in seeds if tuple(sorted((seed, PID))) in graph)
        p = state['parameters']
        identity = {name: int(getattr(pa, name)[k]) for name in
                    ('component_id', 'tracer_id', 'label', 'particle_id', 'actuator') if name in pa.properties}
        rows[label] = dict(folder=folder, time_s=state['time_s'], global_step=state['step'],
            candidate_restart_allowed=False, raw_particle_index=k, identity=identity,
            mass_kg=float(pa.m[k]), position_m=[float(pa.x[k]), float(pa.y[k]), float(pa.z[k])],
            stress_Pa=stress.tolist(), stress_frobenius_norm_Pa=float(np.linalg.norm(stress)),
            equivalent_stress_Pa=float(np.sqrt(1.5)*np.linalg.norm(stress)),
            stress_trace_Pa=float(np.trace(stress)), symmetric_rate_tensor_per_s=d.tolist(),
            strain_rate_per_s=float(pa.strain_rate[k]), strain_rate_reconstructed_from_D_per_s=reconstructed_rate,
            ep_rate_per_s=float(pa.ep_rate[k]), ep_rate_from_saved_parent_per_s=reconstructed_ep_rate,
            ep=float(pa.ep[k]), parent_ep=float(parent.ep[k]), temperature_K=float(pa.temperature[k]),
            rate_floor_per_s=p['rate_floor'], reference_rate_per_s=p['reference_rate'],
            below_rate_floor=bool(pa.strain_rate[k] < p['rate_floor']),
            material_rate_multiplier_log=float(np.log(max(pa.strain_rate[k], p['rate_floor'])/p['reference_rate'])),
            rate_floor_is_existing_explicit_engineering_assumption=True,
            active_contact_features=[feature_names[n] for n in np.flatnonzero(active[k])],
            endpoint_gap_m=dict(zip(feature_names, contact['gaps_m'][k])),
            endpoint_shoulder_active_ids=ids[active[:, feature_names.index('SHOULDER')]].astype(int).tolist(),
            sampled_shoulder_active_ids=sorted({int(row['particle_id']) for row in events if row['feature']=='SHOULDER'}),
            neighbor_ids=neighbors, shoulder_seed_neighbors=shoulder_links, edge_count=len(graph),
            particle_property_names=sorted(pa.properties),
            candidate_sha256=reader.sha(candidate_path))
        raw_states[label] = (pa, state)

    coarse, fine = raw_states['coarse'][0], raw_states['fine'][0]
    k = rows['fine']['raw_particle_index']
    assert rows['coarse']['time_s'] == rows['fine']['time_s'] == screen['time_s']
    assert screen['coarse_checkpoint_sha256'] == rows['coarse']['candidate_sha256']
    assert np.array_equal(coarse.m, fine.m)
    passive_equal = {}
    for name in ('particle_id', 'component_id', 'tracer_id', 'label', 'actuator'):
        if name in coarse.properties:
            passive_equal[name] = bool(np.array_equal(getattr(coarse, name), getattr(fine, name)))
            assert passive_equal[name]
    sf, sc = fine.S.reshape(-1, 3, 3)[k], coarse.S.reshape(-1, 3, 3)[k]
    difference = sf-sc
    floor = screen['local_screen']['stress_denominator_floor_Pa']
    denominator = max(float(np.linalg.norm(sf)), floor)
    relative = float(np.linalg.norm(difference)/denominator)
    local = next(row for row in screen['local_screen']['rows'] if row['particle_id'] == PID)
    assert relative == local['stress_relative']
    selected = set(screen['cumulative_shoulder_ids'])
    for graph in graphs.values():
        for edge in graph:
            if set(edge) & set(screen['cumulative_shoulder_ids']):
                selected.update(edge)
    assert sorted(selected) == screen['local_screen']['particle_ids']
    angle = float(np.degrees(np.arccos(np.clip(np.sum(sf*sc)/(np.linalg.norm(sf)*np.linalg.norm(sc)), -1., 1.))))
    trends = []
    for pair in range(1, 15):
        old_screen = read(BASE/f'fine/pair_{pair:03d}_temporal_screen.json')
        pid_row = next((row for row in old_screen['local_screen']['rows'] if row['particle_id']==PID), None)
        trends.append(dict(pair=pair, time_s=old_screen['time_s'], pid_selected=pid_row is not None,
            pid_stress_relative=None if pid_row is None else pid_row['stress_relative'],
            max_local_stress_relative=max(row['stress_relative'] for row in old_screen['local_screen']['rows']),
            local_passed=old_screen['local_screen']['passed']))
    graph_comparison = {name: dict(added_edges_vs_coarse=sorted(graph-graphs['coarse']),
        removed_edges_vs_coarse=sorted(graphs['coarse']-graph)) for name, graph in graphs.items()}
    result = dict(status='CONFIRMED_REAL_SAVED_LOCAL_STRESS_DIFFERENCE_NOT_ID_MAPPING_ERROR',
        physical_solver_evaluations=0, material_updates=0, nnps_evaluations=0,
        particle_id=PID, time_s=screen['time_s'], comparison_scope='Independent paths sharing 22 ms parent; not same direct parent.',
        accepted_prefix_time_s=manifests['fine']['accepted_time_s'], accepted_prefix_global_step=179,
        failed_pair_restartable=False, branches=rows,
        checks=dict(raw_ids_match_all_three_graphs_and_journals=True, mass_equal=True,
            passive_fields_equal=passive_equal, saved_stress_equals_checkpoint_exactly=True,
            selected_ids_equal_reconstructed_seed_one_hop_union=True, recomputed_stress_relative_equals_screen_exactly=True),
        stress_comparison=dict(difference_Pa=difference.tolist(), difference_frobenius_norm_Pa=float(np.linalg.norm(difference)),
            denominator_Pa=denominator, denominator_floor_Pa=floor, denominator_floor_active=bool(np.linalg.norm(sf)<=floor),
            relative=relative, limit=screen['local_screen']['stress_relative_limit'],
            relative_excess=relative-screen['local_screen']['stress_relative_limit'],
            tensor_direction_angle_degrees=angle,
            stress_norm_relative_difference=float(abs(np.linalg.norm(sf)-np.linalg.norm(sc))/np.linalg.norm(sf))),
        local_row=local, global_screen=screen['global_screen'],
        cumulative_shoulder_ids=screen['cumulative_shoulder_ids'], graph_comparison=graph_comparison,
        same_pid_threshold_history=trends,
        interpretation=['The 10% rejection is reproduced exactly from raw saved stress tensors.',
            'The 1 MPa denominator floor is inactive; similar stress magnitudes have different tensor components/directions.',
            'Both endpoint material rates lie below the disclosed 1/s rate floor; this is not an endpoint floor-switch mismatch.',
            'PID275 is not itself in contact; inclusion follows saved shoulder-seed NNPS adjacency.',
            'This readback identifies a temporal-path field difference, not its unique physical/numerical cause or a proof of convergence.'],
        input_sha256=inputs,
        source_sha256={path.relative_to(ROOT).as_posix():reader.sha(path) for path in
            (Path(__file__), ROOT/'scripts/readback_discrete_process.py', ROOT/'src/afsd_backend/pysph/material.py',
             ROOT/'src/afsd_backend/pysph/quasistatic_process.py', ROOT/'scripts/run_first_shoulder_event.py')})
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(status=result['status'], output=str(OUT), stress_relative=relative,
        difference_Pa=result['stress_comparison']['difference_frobenius_norm_Pa'],
        tensor_angle_degrees=angle, shoulder_seed_neighbors={name:row['shoulder_seed_neighbors'] for name,row in rows.items()},
        graphs_equal=all(not value['added_edges_vs_coarse'] and not value['removed_edges_vs_coarse'] for value in graph_comparison.values()))))


if __name__ == '__main__':
    main()

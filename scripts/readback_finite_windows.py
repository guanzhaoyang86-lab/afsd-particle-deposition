"""Read-only accepted finite-window chain, with independent pair qualification.

Coarse and unpublished fine data are read only to audit references/screens.
Only accepted fine checkpoints become trajectory samples. No trial is created.
"""
from pathlib import Path
import argparse
import json
import math
import time
import numpy as np
import readback_batched_process as base
import readback_shoulder_process as shoulder
import run_discrete_post_start as post
import run_first_shoulder_event as event
from audit_solid_occluded_neighbors import audit_solid_edge_visibility
from afsd_backend.pysph.mechanics import pairs

ROOT, BASE = base.ROOT, base.BASE/'finite_shoulder_windows'
require, read, sha, clean, reader = base.require, base.read, base.sha, base.clean, base.reader
_ORIGINAL_DEFAULT_STAGES = base._default_stages
WINDOWS = {1: (.022, .030), 2: (.030, .038), 3: (.038, .046)}
COARSE_DT, FINE_DT = .00025, .000125
INITIAL_PARENT = base.BASE/'first_shoulder_enclosed/halves/step_002/accepted.npz'
INITIAL_SHA = '5368f67e717b4ce42d3927b8eca4d03a502a7ab59c5899d3d9b0b776c71d2d95'


def same_numbers(first, second, context):
    """Exact structure/flags; allow only floating reduction roundoff."""
    if isinstance(first, dict):
        require(isinstance(second, dict) and set(first) == set(second), context+': keys differ')
        for key in first:
            same_numbers(first[key], second[key], context+'.'+key)
    elif isinstance(first, list):
        require(isinstance(second, list) and len(first) == len(second), context+': lengths differ')
        for number, (a, b) in enumerate(zip(first, second)):
            same_numbers(a, b, context+'.'+str(number))
    elif isinstance(first, float):
        require(math.isfinite(first) and isinstance(second, (float, int)) and
                math.isclose(first, second, rel_tol=5e-13, abs_tol=1e-15), context+': numerical value differs')
    else:
        require(first == second, context+': value differs')


def crossing(state):
    observation = state['bore_observation']
    return {key: observation[key] for key in ('crossing_ids', 'crossing_mass_kg', 'events', 'last_time_s')}


def ordered_edges(i, j):
    edges = np.column_stack((i, j))
    return edges[np.lexsort((edges[:, 1], edges[:, 0]))]


class Audit:
    def __init__(self):
        self.hashes = {}
        self.candidates = 0
        self.screens = 0

    def bind(self, path, expected=None):
        path = Path(path).resolve()
        require(path.is_relative_to(ROOT.resolve()), 'audit path escapes project')
        digest = sha(path)
        require(expected is None or digest == expected, 'hash differs: '+str(path))
        name = path.relative_to(ROOT).as_posix()
        require(name not in self.hashes or self.hashes[name] == digest, 'file changed during audit: '+name)
        self.hashes[name] = digest
        return digest

    def document(self, path, expected=None):
        self.bind(path, expected)
        return read(path)

    def sources(self, mapping):
        for name, digest in mapping.items():
            self.bind(ROOT/name, digest)

    def candidate(self, directory, row, previous_pa, previous_state, p, frozen):
        number = row['step']; folder = directory/f'step_{number:03d}'
        require(row['checkpoint'] == f'step_{number:03d}/NOT_RESTARTABLE_candidate.npz' and
                row['evidence'] == f'step_{number:03d}/candidate_evidence.npz' and
                row['events'] == f'step_{number:03d}/sampled_nonbore_events.json', 'candidate controlled path differs')
        path = directory/row['checkpoint']; self.bind(path, row['sha256'])
        raw = reader.load(str(path)); pa, state = raw['arrays']['solid'], raw['solver_data']
        require(state.get('restart_allowed') is False and state.get('event_candidate_accepted') is False and
                row['restart_allowed'] is False, 'numerical reference became accepted/restartable')
        require(state['parameters'] == p and state['implementation_source_sha256'] == frozen and
                state['source_sha256'] == previous_state['source_sha256'] and
                state['parameter_provenance'] == previous_state['parameter_provenance'] and
                state['step'] == previous_state['step']+1 == row['global_step'] and
                state['time_s'] == row['time_s'] and
                abs(state['time_s']-previous_state['time_s']-row['dt_s']) <= 1e-15,
                'candidate model/global step/time differs')
        reader.validate_interval_energy_state(state['interval_energy'])
        reader.validate_bore_crossing_state(state['bore_observation'])
        require(state['interval_energy']['records'][:-1] == previous_state['interval_energy']['records'] and
                state['interval_energy']['steps'] == state['step'], 'candidate cumulative energy prefix differs')
        require(crossing(state) == row['bore_observation'], 'candidate crossing record differs')
        before, fields = reader.fields(previous_pa), reader.fields(pa)
        require(set(fields) == set(before), 'candidate property inventory differs')
        for name in fields:
            require(np.isfinite(fields[name]).all(), 'nonfinite candidate field: '+name)
        for name in ('particle_id', 'component_id', 'material_id', 'tracer_type', 'wholebar_tag', 'actuator', 'm', 'x0', 'y0', 'z0'):
            require(np.array_equal(fields[name], before[name]), 'candidate persistent field differs: '+name)
        require(np.array_equal(pa.particle_id, previous_pa.particle_id), 'candidate raw ID ordering changed')
        observed = base.observe_square_bore_crossings(previous_state['bore_observation'], reader.xyz(previous_pa), reader.xyz(pa),
            fields['particle_id'].ravel(), fields['component_id'].ravel(), fields['m'].ravel(),
            previous_state['time_s'], state['time_s'], lambda t: reader.expected_pose(p, t))
        require(observed['next_state'] == state['bore_observation'], 'candidate crossing replay differs')
        evidence_path = directory/row['evidence']; self.bind(evidence_path, row['evidence_sha256'])
        with np.load(evidence_path, allow_pickle=False) as graph:
            ids, i, j = (graph[name].copy() for name in ('particle_id', 'i', 'j'))
        require(np.array_equal(ids, pa.particle_id), 'candidate graph raw IDs differ')
        x = np.column_stack((pa.x, pa.y, pa.z))
        journal_path = folder/'FULL_UNCOMMITTED_candidate.npz'; self.bind(journal_path)
        with np.load(journal_path, allow_pickle=False) as saved:
            require(not bool(saved['accepted']) and not bool(saved['restart_allowed']) and
                    np.array_equal(saved['particle_id'], ids) and np.array_equal(saved['position_m'], x),
                    'full candidate journal position or IDs differ')
            checks = dict(density_kg_m3=pa.rho, stress_Pa=pa.S.reshape(-1, 3, 3),
                temperature_K=pa.temperature, plastic_strain=pa.ep,
                velocity_m_s=np.column_stack((pa.u, pa.v, pa.w)), contact_history=state['contact_history'])
            for name, value in checks.items():
                require(np.array_equal(saved[name], value), 'full candidate field differs: '+name)
            require(np.allclose(saved['displacement_m'], x-np.column_stack((previous_pa.x, previous_pa.y, previous_pa.z)),
                                rtol=0., atol=1e-15), 'candidate displacement does not join previous state')
        journal = self.document(folder/'FULL_UNCOMMITTED_candidate.json')
        require(journal['metrics']['force_residual_N'] <= journal['metrics']['force_limit_N'] and
                journal['contact']['sweep']['certified'], 'candidate force/solid sweep qualification absent')
        events = self.document(directory/row['events'], row['events_sha256'])
        require(all(item['feature'] == 'SHOULDER' for item in events), 'candidate event outside shoulder scope')
        audits = {}
        require(set(row['audit_artifacts']) == {'continuous_shell_audit.json', 'scalar_batch_equivalence.json',
                'geometry_audit.json', 'solid_edge_visibility.json'}, 'candidate audit inventory differs')
        for name, artifact in row['audit_artifacts'].items():
            require(artifact['path'] == f'step_{number:03d}/'+name, 'candidate audit controlled path differs')
            audits[name] = self.document(directory/artifact['path'], artifact['sha256'])
        shell, geometry, visible = (audits[name] for name in ('continuous_shell_audit.json', 'geometry_audit.json', 'solid_edge_visibility.json'))
        require(shell['status'] == event.CERTIFIED and shell['additional_finite_patch_queries'] == 0 and
                audits['scalar_batch_equivalence.json']['status'] == 'PASS_SAME_CANDIDATE_COMPLETE_MAP' and
                not geometry['transient_contacts_missed_by_endpoints'] and
                set(geometry['new_non_bore_features']) <= {'SHOULDER'}, 'candidate finite geometry qualification absent')
        require(visible['candidate_sha256'] == row['sha256'] and visible['evidence_sha256'] == row['evidence_sha256'] and
                visible['endpoint_time_s'] == row['time_s'] and visible['physical_solver_evaluations'] == 0,
                'visibility report is not bound to this candidate')
        rebuilt_i, rebuilt_j = pairs(pa)
        require(np.array_equal(ordered_edges(i, j), ordered_edges(rebuilt_i, rebuilt_j)), 'candidate saved graph differs from real endpoint NNPS')
        fresh = audit_solid_edge_visibility(x, i, j, reader.expected_pose(p, state['time_s']), p)
        require(fresh['status'] == 'PASS_NO_SOLID_OCCLUDED_ENDPOINT_EDGE' and fresh['occluded_edge_count'] == 0,
                'candidate endpoint now fails pure visibility replay')
        for key, value in fresh.items():
            same_numbers(value, visible[key], 'visibility.'+key)
        self.candidates += 1
        return dict(pa=pa, state=state, graph=(ids, i, j), events=events, path=path, row=row)


def qualify_windows(audit, count):
    require(count in WINDOWS, 'only explicit first 1, 2, or 3 windows supported')
    audit.bind(INITIAL_PARENT, INITIAL_SHA)
    parent_pa, parent_state = reader.restore(INITIAL_PARENT)
    parent_path = INITIAL_PARENT
    p, provenance = parent_state['parameters'], parent_state['parameter_provenance']
    old_screen = audit.document(INITIAL_PARENT.parent.parent/'temporal_comparison.json')
    seeds = set(old_screen['event_particle_ids'])
    stages, reports = [], []
    for number in range(1, count+1):
        start, end = WINDOWS[number]
        parent_hash = sha(parent_path)
        require(abs(parent_state['time_s']-start) <= 1e-15 and not parent_state['bore_observation']['crossing_ids'],
                'window starts from wrong time or after an earlier crossing')
        coarse_dir, fine_dir = (BASE/f'window_{number}'/stage for stage in ('coarse', 'fine'))
        coarse_meta = audit.document(coarse_dir/'manifest.json')
        coarse_rows = audit.document(coarse_dir/'candidate_history.json')
        require('RUNNING' not in coarse_meta['status'] and coarse_meta['accepted_steps'] == 0 and
                not list(coarse_dir.glob('step_*/accepted.npz')), 'coarse still running or wrongly published')
        frozen = coarse_meta['source_sha256']; audit.sources(frozen)
        require(coarse_meta['parent_sha256'] == parent_hash and
                (ROOT/coarse_meta['parent_checkpoint']).resolve() == parent_path.resolve() and
                coarse_meta['provenance'] == provenance and coarse_meta['registered_window'] == number and
                abs(coarse_meta['starting_time_s']-start) <= 1e-15 and abs(coarse_meta['accepted_time_s']-start) <= 1e-15 and
                coarse_meta['registered_start_s'] == start and coarse_meta['requested_time_s'] == end and coarse_meta['dt_s'] == COARSE_DT and
                coarse_meta['candidate_steps'] == len(coarse_rows), 'coarse window origin or registration differs')
        coarse, previous_pa, previous_state, previous_hash = [], parent_pa, parent_state, parent_hash
        coarse_seeds = set(seeds)
        for index, row in enumerate(coarse_rows, 1):
            require(row['step'] == index and row['dt_s'] == COARSE_DT and row['parent_sha256'] == previous_hash and
                    abs(row['time_s']-(start+index*COARSE_DT)) <= 1e-15, 'coarse numerical reference chain differs')
            item = audit.candidate(coarse_dir, row, previous_pa, previous_state, p, frozen)
            coarse_seeds.update(event_row['particle_id'] for event_row in item['events'])
            require(sorted(coarse_seeds) == row['cumulative_shoulder_ids'], 'coarse cumulative event IDs differ')
            coarse.append(item); previous_pa, previous_state, previous_hash = item['pa'], item['state'], row['sha256']
        if coarse_meta['status'] != 'NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED':
            require(number == count and coarse_meta['status'] == 'FAILED_NO_NEW_ACCEPTED_STATE' and
                    not (fine_dir/'manifest.json').exists(), 'failed coarse can only terminate the requested series')
            reports.append(dict(window=number, coarse_status=coarse_meta['status'], fine_status='NOT_RUN', accepted_steps=0))
            break
        require(len(coarse) == 32 and abs(coarse_meta['candidate_time_s']-end) <= 1e-15, 'coarse reference incomplete')
        meta = audit.document(fine_dir/'manifest.json')
        history = audit.document(fine_dir/'accepted_history.json')
        candidate_rows = audit.document(fine_dir/'candidate_history.json')
        status = meta['status']; accepted = len(history)
        require(status in ('PASS_WINDOW_NOT_TRANSPORT', 'PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION',
                'FAILED_WITH_ACCEPTED_PREFIX', 'FAILED_NO_NEW_ACCEPTED_STATE'), 'fine stage is not finalized')
        require(meta['source_sha256'] == frozen and meta['provenance'] == provenance and
                meta['parent_sha256'] == parent_hash and (ROOT/meta['parent_checkpoint']).resolve() == parent_path.resolve() and
                meta['accepted_steps'] == accepted and meta['candidate_steps'] == len(candidate_rows) and
                meta['registered_window'] == number and meta['dt_s'] == FINE_DT and
                abs(meta['starting_time_s']-start) <= 1e-15 and meta['registered_start_s'] == start and meta['requested_time_s'] == end and
                accepted <= len(candidate_rows) <= min(64, accepted+2),
                'fine stage source/count/origin differs')
        require((status != 'FAILED_NO_NEW_ACCEPTED_STATE' or accepted == 0) and
                (status != 'FAILED_WITH_ACCEPTED_PREFIX' or accepted > 0), 'failure status misstates accepted prefix')
        accepted_files = {path.resolve() for path in fine_dir.glob('step_*/accepted.npz')}
        require(accepted_files == {(fine_dir/row['checkpoint']).resolve() for row in history}, 'unlisted or missing accepted checkpoint')
        fine, pair_seeds = [], {}
        previous_pa, previous_state = parent_pa, parent_state
        fine_seeds = set(seeds)
        for index, row in enumerate(candidate_rows, 1):
            expected_parent = (candidate_rows[index-2]['sha256'] if index % 2 == 0 else
                parent_hash if index == 1 else history[index-2]['sha256'])
            require(row['step'] == index and row['dt_s'] == FINE_DT and row['parent_sha256'] == expected_parent and
                    abs(row['time_s']-(start+index*FINE_DT)) <= 1e-15, 'fine numerical candidate parent/time differs')
            item = audit.candidate(fine_dir, row, previous_pa, previous_state, p, frozen)
            fine_seeds.update(event_row['particle_id'] for event_row in item['events'])
            require(sorted(fine_seeds) == row['cumulative_shoulder_ids'], 'fine cumulative event IDs differ')
            fine.append(item); previous_pa, previous_state = item['pa'], item['state']
            if index % 2 == 0:
                pair = index//2
                fine_seeds.update(coarse_rows[pair-1]['cumulative_shoulder_ids'])
                pair_seeds[pair] = sorted(fine_seeds)
        screens_path = fine_dir/'screen_history.json'
        screen_rows = audit.document(screens_path) if screens_path.exists() else []
        pair_count = (accepted+1)//2
        require(len(screen_rows) == accepted//2 or
                (status == 'FAILED_WITH_ACCEPTED_PREFIX' and accepted > 0 and len(screen_rows) == pair_count-1),
                'screen history does not cover the accepted prefix')
        screen_hashes = {}
        for pair in range(1, pair_count+1):
            require(len(fine) >= 2*pair, 'accepted half lacks both prepublication pair candidates')
            screen_file = fine_dir/f'pair_{pair:03d}_temporal_screen.json'
            screen = audit.document(screen_file); digest = sha(screen_file)
            if pair <= len(screen_rows):
                entry = screen_rows[pair-1]
                require(entry == dict(pair=pair, screen=screen_file.name, sha256=digest), 'screen history hash/path differs')
            else:
                require(status == 'FAILED_WITH_ACCEPTED_PREFIX' and pair == pair_count,
                        'missing screen index permitted only for interrupted final publication')
            reference, first, second = coarse[pair-1], fine[2*pair-2], fine[2*pair-1]
            selected = event.one_hop_ids(pair_seeds[pair], [reference['graph'], first['graph'], second['graph']])
            expected_global = post.compare(second['pa'], reference['pa'], p)
            expected_local = event.local_compare(second['pa'], reference['pa'], selected)
            require(screen['pair'] == pair and abs(screen['time_s']-(start+pair*COARSE_DT)) <= 1e-15 and
                    screen['parent_sha256'] == parent_hash and screen['coarse_checkpoint_sha256'] == reference['row']['sha256'] and
                    screen['cumulative_shoulder_ids'] == pair_seeds[pair], 'pair screen time/origin/coarse/seeds differ')
            same_numbers(expected_global, screen['global_screen'], 'recomputed_global_screen')
            same_numbers(expected_local, screen['local_screen'], 'recomputed_local_screen')
            require(expected_global['pass'] and expected_local['passed'] and screen['crossing_ids_match'] and
                    screen['coarse_bore_observation'] == crossing(reference['state']) and
                    screen['fine_bore_observation'] == crossing(second['state']) and
                    crossing(reference['state'])['crossing_ids'] == crossing(second['state'])['crossing_ids'],
                    'accepted pair does not pass both time and strict crossing screens')
            screen_hashes[pair] = digest; audit.screens += 1
        unaccepted_screen_rows = []
        for pair in range(pair_count+1, len(fine)//2+1):
            screen_file = fine_dir/f'pair_{pair:03d}_temporal_screen.json'
            if not screen_file.exists():
                continue
            screen = audit.document(screen_file)
            reference, first, second = coarse[pair-1], fine[2*pair-2], fine[2*pair-1]
            selected = event.one_hop_ids(pair_seeds[pair], [reference['graph'], first['graph'], second['graph']])
            expected_global = post.compare(second['pa'], reference['pa'], p)
            expected_local = event.local_compare(second['pa'], reference['pa'], selected)
            require(screen['pair'] == pair and screen['parent_sha256'] == parent_hash and
                    screen['coarse_checkpoint_sha256'] == reference['row']['sha256'] and
                    screen['cumulative_shoulder_ids'] == pair_seeds[pair] and
                    abs(screen['time_s']-(start+pair*COARSE_DT)) <= 1e-15, 'unaccepted pair reference differs')
            same_numbers(expected_global, screen['global_screen'], 'unaccepted_global_screen')
            same_numbers(expected_local, screen['local_screen'], 'unaccepted_local_screen')
            same_crossing_ids = crossing(reference['state'])['crossing_ids'] == crossing(second['state'])['crossing_ids']
            require(screen['coarse_bore_observation'] == crossing(reference['state']) and
                    screen['fine_bore_observation'] == crossing(second['state']) and
                    screen['crossing_ids_match'] == same_crossing_ids, 'unaccepted crossing comparison differs')
            if 'last_screen' in meta:
                require(meta['last_screen'] == screen, 'manifest last failed/unpublished screen differs')
            unaccepted_screen_rows.append(dict(pair=pair, time_s=screen['time_s'], global_screen=expected_global,
                local_screen=expected_local, crossing_ids_match=same_crossing_ids,
                screen_sha256=sha(screen_file), accepted=False, exported_samples=0))
        direct_hash = parent_hash
        for index, row in enumerate(history, 1):
            require(row['step'] == index and row['checkpoint'] == f'step_{index:03d}/accepted.npz' and
                    row['parent_sha256'] == direct_hash and row['dt_s'] == FINE_DT, 'accepted fine parent/path differs')
            path = fine_dir/row['checkpoint']; audit.bind(path, row['sha256'])
            pa, state = reader.restore(path); item = fine[index-1]; pair = (index+1)//2
            reader.identical_fields(reader.fields(pa), reader.fields(item['pa']), 'accepted/reference candidate fields')
            for key in ('time_s', 'step', 'energy', 'interval_energy', 'bore_observation', 'contact_history', 'free'):
                require(clean(state[key]) == clean(item['state'][key]), 'accepted/candidate metadata differs: '+key)
            require(state['window_pair_index'] == pair and state['window_temporal_screen_sha256'] == row['temporal_screen_sha256'] == screen_hashes[pair] and
                    state['window_coarse_reference_sha256'] == row['coarse_reference_sha256'] == coarse_rows[pair-1]['sha256'],
                    'accepted checkpoint is not bound to screened same-time reference')
            direct_hash = row['sha256']
        if accepted:
            parent_path = fine_dir/history[-1]['checkpoint']; parent_pa, parent_state = reader.restore(parent_path)
            require(meta['accepted_time_s'] == parent_state['time_s'] and meta['interval_energy'] == parent_state['interval_energy'] and
                    meta['crossing_mass_kg'] == parent_state['bore_observation']['crossing_mass_kg'], 'fine manifest final accepted state differs')
            stages.append((f'finite_window_{number}_fine', fine_dir, accepted, status, 'event'))
        else:
            require(abs(meta['accepted_time_s']-start) <= 1e-15, 'zero accepted prefix advanced official time')
        crossing_ids = parent_state['bore_observation']['crossing_ids']
        if status == 'PASS_WINDOW_NOT_TRANSPORT':
            require(accepted == 64 and len(screen_rows) == 32 and not crossing_ids and
                    abs(parent_state['time_s']-end) <= 1e-15, 'successful finite window incomplete or already crossed')
        elif status == 'PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION':
            require(accepted % 2 == 0 and accepted > 0 and crossing_ids and
                    meta['first_passage_pair'] == accepted//2 and meta['bore_observation'] == crossing(parent_state),
                    'first passage status lacks final screened pair')
            require(all(not crossing(fine[2*pair-1]['state'])['crossing_ids'] for pair in range(1, accepted//2)),
                    'first crossing occurred before the declared stopping pair')
        if number < count:
            require(status == 'PASS_WINDOW_NOT_TRANSPORT' and not crossing_ids,
                    'cannot append a later window after a partial/failed/crossing window')
        reports.append(dict(window=number, coarse_status=coarse_meta['status'], fine_status=status,
            accepted_steps=accepted, candidate_steps=len(candidate_rows), independently_recomputed_qualified_pairs=pair_count,
            screen_index_entries=len(screen_rows), final_accepted_time_s=parent_state['time_s'],
            interrupted_final_pair_publication=bool(accepted % 2 or len(screen_rows) < pair_count),
            coarse_samples_exported=0, first_passage_ids=list(crossing_ids),
            independently_recomputed_unaccepted_pair_screens=unaccepted_screen_rows))
        if pair_count:
            seeds.update(pair_seeds[pair_count])
    return stages, reports, parent_state


def package(output, windows):
    begin = time.perf_counter(); output = Path(output).resolve()
    require(not output.exists(), 'refuse to overwrite finite-window readback')
    audit = Audit()
    for path in (Path(__file__), Path(base.__file__), Path(shoulder.__file__), Path(reader.__file__),
                 Path(event.__file__), Path(post.__file__), ROOT/'scripts/audit_solid_occluded_neighbors.py'):
        audit.bind(path)
    for name in ('summary.json', 'accepted_trajectory.npz'):
        audit.bind(base.BASE/'readback_22ms'/name)
    extra_stages, reports, final_state = qualify_windows(audit, windows)
    stages = shoulder.shoulder_stages()+extra_stages
    require(base._default_stages is _ORIGINAL_DEFAULT_STAGES, 'another reader stage adapter is installed')
    base._default_stages = lambda: list(stages)
    try:
        summary = base.package(output)
    finally:
        base._default_stages = _ORIGINAL_DEFAULT_STAGES
    require(summary['accepted_steps'] == final_state['step'] and summary['samples'] == final_state['step']+1 and
            summary['accepted_time_s'] == final_state['time_s'], 'exported actual chain endpoint differs from qualified prefix')
    for name, digest in audit.hashes.items():
        require(sha(ROOT/name) == digest, 'source/evidence changed during finite readback: '+name)
    summary.update(finite_window_reader_sha256=sha(Path(__file__)), requested_window_count=windows,
        finite_window_qualification=reports, independently_recomputed_window_pair_count=audit.screens,
        endpoint_visibility_replayed_candidate_count=audit.candidates, real_nnps_rebuilds=audit.candidates,
        added_physical_solver_evaluations=0, coarse_frames_exported=0,
        finite_window_source_and_evidence_sha256=audit.hashes,
        controlled_stage_list=[dict(name=name, directory=directory.relative_to(ROOT).as_posix(),
            accepted_count=count, required_status=status, audit_mode=mode)
            for name, directory, count, status, mode in stages],
        original_default_stages_restored=True, append_argument_used=False,
        finite_window_audit_wall_seconds=time.perf_counter()-begin)
    (output/'summary.json').write_text(json.dumps(clean(summary), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(status='PASS_FINITE_WINDOW_ACCEPTED_CHAIN_READBACK', output=str(output),
        accepted_steps=summary['accepted_steps'], samples=summary['samples'], accepted_time_s=summary['accepted_time_s'],
        independent_pair_screens=audit.screens, candidate_visibility_replays=audit.candidates)))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--windows', type=int, choices=(1, 2, 3), required=True)
    parser.add_argument('--output', type=Path, required=True)
    arguments = parser.parse_args()
    package(arguments.output, arguments.windows)

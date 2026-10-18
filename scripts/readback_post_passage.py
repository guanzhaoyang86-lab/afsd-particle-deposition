"""Read-only post-passage response, with crossed-ID local scope and raw fields.

The accepted first-passage chain is retained intact. The single new window is
37.375->41.375 ms, with at most 32 coarse pairs / 64 accepted fine increments.
Extra crossed-particle observations are independently recomputed from fields.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import readback_refined_continuation as continued

refined, prior = continued.refined, continued.prior
base, shoulder, reader = prior.base, prior.shoulder, prior.reader
ROOT = prior.ROOT
BASE = base.BASE/'post_passage_short_response'
COARSE_DT, FINE_DT = .000125, .0000625
INITIAL_PARENT = base.BASE/'refined_shoulder_continuation/window_1/fine/step_118/accepted.npz'
INITIAL_SHA = '9f78a757bb9897e97598b09bee5842f29f6d6d7a2d217dfb39add0fb4a096b67'
same_numbers, crossing = prior.same_numbers, prior.crossing
require, read, sha, clean, event, post = prior.require, prior.read, prior.sha, prior.clean, prior.event, prior.post
_ORIGINAL_DEFAULT_STAGES = base._default_stages
FROZEN_CONTINUATION_READER_SHA256 = '584c5b3b5a5badbcfa05731d26e5ecbb6e77e2a3dd4277acce85e0d2c6b27e17'


class PostAudit(prior.Audit):
    def set_post_parent(self, pa, state, path, digest):
        require(state['step'] == 399 and abs(state['time_s']-.037375) <= 1e-15 and
                state['bore_observation']['crossing_ids'] == [262] and sha(path) == digest == INITIAL_SHA,
                'wrong post-passage parent state')
        self.post_parent_pa, self.post_parent_state, self.post_parent_sha = pa, state, digest

    def parent_response(self, pa, state):
        initial, parent = self.post_parent_pa, self.post_parent_state
        p = state['parameters']
        pose, pose0 = reader.expected_pose(p, state['time_s']), reader.expected_pose(p, parent['time_s'])
        body = (np.column_stack((pa.x, pa.y, pa.z))-pose['origin'])@pose['rotation']
        body0 = (np.column_stack((initial.x, initial.y, initial.z))-pose0['origin'])@pose0['rotation']
        rows = []
        for pid in parent['bore_observation']['crossing_ids']:
            index = np.flatnonzero(pa.particle_id == pid); index0 = np.flatnonzero(initial.particle_id == pid)
            require(len(index) == len(index0) == 1, 'crossed particle identity absent or repeated')
            i, j = int(index[0]), int(index0[0])
            radius, radius0 = float(np.linalg.norm(body[i, :2])), float(np.linalg.norm(body0[j, :2]))
            rows.append(dict(particle_id=pid, plastic_strain=float(pa.ep[i]),
                plastic_strain_increment_from_parent=float(pa.ep[i]-initial.ep[j]),
                temperature_K=float(pa.temperature[i]), temperature_increment_from_parent_K=float(pa.temperature[i]-initial.temperature[j]),
                tool_body_position_m=body[i].tolist(), tool_radial_m=radius, tool_radial_increment_from_parent_m=radius-radius0,
                tool_z_m=float(body[i, 2]), tool_z_increment_from_parent_m=float(body[i, 2]-body0[j, 2])))
        return clean(dict(classification='OBSERVATION_ONLY_NOT_DEPOSITION_OR_FORCE_INPUT', time_s=state['time_s'],
            parent_time_s=parent['time_s'], parent_sha256=self.post_parent_sha, tool_pose=pose,
            parent_crossed_ids=parent['bore_observation']['crossing_ids'], current_bore_observation=crossing(state), rows=rows))

    def candidate(self, directory, row, previous_pa, previous_state, p, frozen):
        item = super().candidate(directory, row, previous_pa, previous_state, p, frozen)
        if directory.resolve().parent == BASE.resolve():
            expected_path = f"step_{row['step']:03d}/parent_crossed_response.json"
            require(row['parent_crossed_response'] == expected_path, 'post-passage observation path differs')
            response = self.document(directory/expected_path, row['parent_crossed_response_sha256'])
            expected = self.parent_response(item['pa'], item['state'])
            same_numbers(expected, response, 'candidate_parent_crossed_response')
            item['parent_response'] = expected
        return item


def post_local_scope(audit, screen, reference, first, second, seeds):
    crossed = set(reference['state']['bore_observation']['crossing_ids'])
    crossed.update(first['state']['bore_observation']['crossing_ids'])
    crossed.update(second['state']['bore_observation']['crossing_ids'])
    require(screen['cumulative_crossing_ids'] == sorted(crossed) and
            screen['local_scope'] == 'CUMULATIVE_SHOULDER_AND_CROSSING_IDS_ONE_HOP_UNION',
            'post-passage local observation scope omitted crossed IDs')
    same_numbers(reference['parent_response'], screen['coarse_parent_crossed_response'], 'coarse_screen_parent_response')
    same_numbers(second['parent_response'], screen['fine_parent_crossed_response'], 'fine_screen_parent_response')
    return event.one_hop_ids(set(seeds)|crossed, [reference['graph'], first['graph'], second['graph']])


def qualify_post_passage(audit):
    count = 1  # Only the registered four-millisecond post-passage window.
    audit.bind(INITIAL_PARENT, INITIAL_SHA)
    parent_pa, parent_state = reader.restore(INITIAL_PARENT)
    parent_path = INITIAL_PARENT
    p, provenance = parent_state['parameters'], parent_state['parameter_provenance']
    old_screen = audit.document(INITIAL_PARENT.parent.parent/'pair_059_temporal_screen.json')
    seeds = set(old_screen['cumulative_shoulder_ids'])
    audit.set_post_parent(parent_pa, parent_state, parent_path, INITIAL_SHA)
    stages, reports = [], []
    for number in (1,):
        start = parent_state['time_s']; end = start+.004
        parent_hash = sha(parent_path)
        require(abs(parent_state['time_s']-start) <= 1e-15 and parent_state['bore_observation']['crossing_ids'] == [262],
                'post-passage registered parent time or crossed identity differs')
        coarse_dir, fine_dir = (BASE/stage for stage in ('coarse', 'fine'))
        coarse_meta = audit.document(coarse_dir/'manifest.json')
        coarse_rows = audit.document(coarse_dir/'candidate_history.json')
        require('RUNNING' not in coarse_meta['status'] and coarse_meta['accepted_steps'] == 0 and
                not list(coarse_dir.glob('step_*/accepted.npz')), 'coarse still running or wrongly published')
        frozen = coarse_meta['source_sha256']; audit.sources(frozen)
        require(coarse_meta['parent_crossed_ids'] == [262] and coarse_meta['passage_window'] == 1, 'coarse post-passage origin differs')
        require(coarse_meta['parent_sha256'] == parent_hash and
                (ROOT/coarse_meta['parent_checkpoint']).resolve() == parent_path.resolve() and
                coarse_meta['provenance'] == provenance and coarse_meta['registered_refinement_level'] == number and
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
            reports.append(dict(post_passage_window=1, coarse_status=coarse_meta['status'], fine_status='NOT_RUN', accepted_steps=0))
            break
        require(len(coarse) == 32 and abs(coarse_meta['candidate_time_s']-end) <= 1e-15, 'coarse reference incomplete')
        meta = audit.document(fine_dir/'manifest.json')
        history = audit.document(fine_dir/'accepted_history.json')
        candidate_rows = audit.document(fine_dir/'candidate_history.json')
        status = meta['status']; accepted = len(history)
        require(meta['parent_crossed_ids'] == [262] and meta['passage_window'] == 1, 'fine post-passage origin differs')
        require(status in ('POST_PASSAGE_SHORT_RESPONSE_NOT_DEPOSITION',
                'FAILED_WITH_ACCEPTED_PREFIX', 'FAILED_NO_NEW_ACCEPTED_STATE'), 'fine stage is not finalized')
        require(meta['source_sha256'] == frozen and meta['provenance'] == provenance and
                meta['parent_sha256'] == parent_hash and (ROOT/meta['parent_checkpoint']).resolve() == parent_path.resolve() and
                meta['accepted_steps'] == accepted and meta['candidate_steps'] == len(candidate_rows) and
                meta['registered_refinement_level'] == number and meta['dt_s'] == FINE_DT and
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
            selected = post_local_scope(audit, screen, reference, first, second, pair_seeds[pair])
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
            selected = post_local_scope(audit, screen, reference, first, second, pair_seeds[pair])
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
            stages.append(('post_passage_fine', fine_dir, accepted, status, 'event'))
        else:
            require(abs(meta['accepted_time_s']-start) <= 1e-15, 'zero accepted prefix advanced official time')
        crossing_ids = parent_state['bore_observation']['crossing_ids']
        if status == 'POST_PASSAGE_SHORT_RESPONSE_NOT_DEPOSITION':
            require(accepted == 64 and len(screen_rows) == 32 and
                    abs(parent_state['time_s']-end) <= 1e-15, 'successful finite window incomplete or already crossed')
        if status == 'POST_PASSAGE_SHORT_RESPONSE_NOT_DEPOSITION':
            same_numbers(audit.parent_response(parent_pa, parent_state), meta['final_parent_crossed_response'], 'final_parent_crossed_response')
            require(meta['bore_observation'] == crossing(parent_state), 'post-passage final crossing inventory differs')
        if number < count:
            require(status == 'POST_PASSAGE_SHORT_RESPONSE_NOT_DEPOSITION' and not crossing_ids,
                    'cannot append a later window after a partial/failed/crossing window')
        reports.append(dict(post_passage_window=1, coarse_status=coarse_meta['status'], fine_status=status,
            accepted_steps=accepted, candidate_steps=len(candidate_rows), independently_recomputed_qualified_pairs=pair_count,
            screen_index_entries=len(screen_rows), final_accepted_time_s=parent_state['time_s'],
            interrupted_final_pair_publication=bool(accepted % 2 or len(screen_rows) < pair_count),
            coarse_samples_exported=0, current_crossing_ids=list(crossing_ids),
            parent_crossed_response_at_last_accepted=audit.parent_response(parent_pa, parent_state),
            independently_recomputed_unaccepted_pair_screens=unaccepted_screen_rows))
        if pair_count:
            seeds.update(pair_seeds[pair_count])
    return stages, reports, parent_state


def package(output):
    begin = time.perf_counter(); output = Path(output).resolve()
    require(not output.exists(), 'refuse to overwrite post-passage readback')
    require(sha(Path(prior.__file__)) == refined.FROZEN_PRIOR_READER_SHA256 and
            sha(Path(refined.__file__)) == continued.FROZEN_REFINEMENT_READER_SHA256 and
            sha(Path(continued.__file__)) == FROZEN_CONTINUATION_READER_SHA256, 'frozen reader helper changed')
    for stage in ('coarse', 'fine'):
        path = BASE/stage/'manifest.json'
        if path.exists():
            require('RUNNING' not in read(path)['status'], 'post-passage stage is still running')
    audit = PostAudit()
    for path in (Path(__file__), Path(prior.__file__), Path(refined.__file__), Path(continued.__file__),
                 Path(base.__file__), Path(shoulder.__file__), Path(reader.__file__), Path(event.__file__), Path(post.__file__),
                 ROOT/'scripts/audit_solid_occluded_neighbors.py', ROOT/'scripts/run_post_passage_short_response.py',
                 ROOT/'docs/POST_PASSAGE_SHORT_RESPONSE_CONTRACT_20261004.md'):
        audit.bind(path)
    for old_directory in ('readback_22ms', 'readback_25_25ms', 'readback_30ms_refined', 'readback_first_passage_37_375ms'):
        for name in ('summary.json', 'accepted_trajectory.npz'):
            audit.bind(base.BASE/old_directory/name)
    refined_stages, refined_reports, refined_state = refined.qualify_refinement(audit)
    require(len(refined_stages) == 1 and refined_reports[0]['fine_status'] == 'PASS_WINDOW_NOT_TRANSPORT' and
            refined_state['step'] == 281, 'complete 30 ms refined chain absent')
    continued_stages, continued_reports, passage_state = continued.qualify_continuation(audit, 1)
    require(len(continued_stages) == 1 and continued_reports[0]['fine_status'] ==
            'PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION' and passage_state['step'] == 399 and
            passage_state['bore_observation']['crossing_ids'] == [262], 'complete qualified first-passage chain absent')
    inherited_candidates, inherited_screens = audit.candidates, audit.screens
    extra_stages, reports, final_state = qualify_post_passage(audit)
    stages = shoulder.shoulder_stages()+refined_stages+continued_stages+extra_stages
    require(base._default_stages is _ORIGINAL_DEFAULT_STAGES, 'another stage adapter is installed')
    base._default_stages = lambda: list(stages)
    try:
        summary = base.package(output)
    finally:
        base._default_stages = _ORIGINAL_DEFAULT_STAGES
    require(summary['accepted_steps'] == final_state['step'] and summary['samples'] == final_state['step']+1 and
            summary['accepted_time_s'] == final_state['time_s'], 'post-passage export endpoint differs')
    require(all('finite_shoulder_windows/' not in row['checkpoint'] and '/coarse/' not in row['checkpoint']
                for row in summary['rows']), 'old failed prefix or coarse path was mixed into actual samples')
    passage = continued.accepted_passage_evidence(summary, output, audit)
    parent_response_rows = [audit.parent_response(audit.post_parent_pa, audit.post_parent_state)]
    for row in summary['rows']:
        checkpoint = ROOT/row['checkpoint']
        if checkpoint.resolve().is_relative_to(BASE.resolve()):
            parent_response_rows.append(audit.document(checkpoint.parent/'parent_crossed_response.json'))
    require(len(parent_response_rows) == 1+sum(row['accepted_steps'] for row in reports),
            'post-passage response history includes unpublished or omits accepted samples')
    for name, digest in audit.hashes.items():
        require(sha(ROOT/name) == digest, 'source/evidence changed during post-passage readback: '+name)
    summary.update(post_passage_reader_sha256=sha(Path(__file__)),
        first_passage_parent_checkpoint=INITIAL_PARENT.relative_to(ROOT).as_posix(), first_passage_parent_sha256=INITIAL_SHA,
        inherited_refinement_qualification=refined_reports, inherited_continuation_qualification=continued_reports,
        post_passage_qualification=reports, independently_recomputed_inherited_pair_count=inherited_screens,
        independently_recomputed_post_passage_pair_count=audit.screens-inherited_screens,
        inherited_candidate_visibility_replays=inherited_candidates,
        post_passage_candidate_visibility_replays=audit.candidates-inherited_candidates,
        real_nnps_rebuilds=audit.candidates, added_physical_solver_evaluations=0,
        coarse_frames_exported=0, old_25_25ms_branch_frames_exported=0,
        accepted_centre_passage_evidence=passage, parent_crossed_response_actual_samples=parent_response_rows,
        post_passage_source_and_evidence_sha256=audit.hashes,
        controlled_stage_list=[dict(name=name, directory=directory.relative_to(ROOT).as_posix(),
            accepted_count=count, required_status=status, audit_mode=mode)
            for name, directory, count, status, mode in stages],
        original_default_stages_restored=True, append_argument_used=False,
        post_passage_audit_wall_seconds=time.perf_counter()-begin,
        post_passage_interpretation='Raw signed ep/T/tool-radial/z changes of previously crossed particles are observations; nonzero increments alone do not validate spreading, bonding, deposition, or experiment.')
    (output/'summary.json').write_text(json.dumps(clean(summary), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(status='PASS_POST_PASSAGE_ACCEPTED_CHAIN_READBACK', output=str(output),
        accepted_steps=summary['accepted_steps'], samples=summary['samples'], accepted_time_s=summary['accepted_time_s'],
        post_passage_pair_screens=audit.screens-inherited_screens, crossing_ids=summary['crossing_ids'],
        crossing_mass_kg=summary['final']['crossing_mass_kg'])))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    arguments = parser.parse_args()
    package(arguments.output)

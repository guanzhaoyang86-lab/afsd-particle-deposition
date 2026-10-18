"""Verify real event substeps without changing any frozen receipt validator.

Two screened fine macro endpoints remain required. Their internal physical
steps retain every original machine and energy record instead of pretending
that an event-resolved macro interval was one physical state transition.
"""
import math
import numbers
from pathlib import Path
import numpy as np

import run_gpu_cached_preflight_20261008 as legacy
from visible_support import filter_visible_support, _hash_array, _topology_hash

SCHEMA = 'EXPLICIT_COUPLED_EVENT_SUBSTEPS_V1'
QUALIFIED = 'PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN'
ROLES = ('coarse', 'fine_1', 'fine_2')
require = legacy.require


def _finite(value):
    return isinstance(value, numbers.Real) and not isinstance(value, (bool, np.bool_)) and math.isfinite(value)


def _positions(pa):
    return np.column_stack((pa.x, pa.y, pa.z))


def _validate_event_resolution(transaction, meta, context, following, artifact, path, state_validator):
    """Read back bounds, chosen physical states, blocked proof and restored graph."""
    row, pa, state, parent_pa, parent_state, qualification, role_dir = context
    require(following is not None, 'activation has no selected post-event physical substep')
    next_row, next_pa, next_state, *_ = following
    leaf = path(row['checkpoint']).parent
    resolution_path = artifact(leaf/'event_resolution.json', leaf)
    require('event_resolution_path' in qualification and 'event_resolution_sha256' in qualification
            and path(qualification['event_resolution_path']) == resolution_path
            and transaction.sha(resolution_path) == qualification['event_resolution_sha256'],
            'activation qualification must bind exact resolution')
    resolution = transaction.read(resolution_path)
    require(resolution['status'] == 'PASS_BOUNDED_PHYSICALLY_RESOLVED_ACTIVATION'
            and resolution['protocol'] == meta['event_protocol']
            and resolution['incoming_restored'] is True and resolution['original_active_path_certified'] is True
            and resolution['interpolation_accepted'] is False
            and resolution['accepted'] is False and resolution['restart_allowed'] is False,
            'invalid activation resolution flags or protocol')
    pairs = resolution['incoming_pid_pairs']
    require(isinstance(pairs, list) and len(pairs) == 1 and len(pairs[0]) == 2
            and all(isinstance(v, numbers.Integral) and not isinstance(v, bool) for v in pairs[0])
            and 0 <= pairs[0][0] < pairs[0][1]
            and qualification.get('delayed_activation_pairs') == pairs,
            'exact registered single incoming pair required')
    require(path(resolution['high_checkpoint']) == path(row['checkpoint'])
            and path(resolution['first_post_event_checkpoint']) == path(next_row['checkpoint'])
            and path(next_row['parent_checkpoint']) == path(row['checkpoint'])
            and next_row['parent_sha256'] == row['sha256'], 'activation chosen high/next leaf differs')
    location_path = artifact(resolution['location'], role_dir)
    require(transaction.sha(location_path) == resolution['location_sha256'], 'activation location hash differs')
    location = transaction.read(location_path)
    require(location['status'] == 'BRACKETED_ACTIVATION' and location['accepted'] is False
            and location['restart_allowed'] is False and location['diagnostic_only'] is True
            and location['interpolated_states'] == 0 and location['assumed_zero_endpoint_used_as_result'] is False
            and location['lo']['blocked'] is True and location['hi']['blocked'] is False,
            'activation lacks two real blocked/visible bracket endpoints')
    tolerance = resolution['tolerance_s']; scale = meta.get('event_tolerance_scale', 1.)
    require(scale in (1., .5, .25) and _finite(tolerance) and tolerance > 0.
            and _finite(location['full_dt_s']) and location['full_dt_s'] >= row['dt_s']
            and abs(tolerance-min(1e-7, location['full_dt_s']/1024.)*scale) <= 1e-19
            and location['time_tolerance_s'] == tolerance
            and location['position_difference_limit_m'] == 5e-8,
            'activation tolerance differs from registered bounds')
    high = location['hi']; low = location['lo']
    low_path = artifact(resolution['low_checkpoint'], role_dir)
    require(path(low['checkpoint']) == low_path and transaction.sha(low_path) == low['sha256']
            and path(high['checkpoint']) == path(row['checkpoint']) and high['sha256'] == row['sha256']
            and high['dt_s'] == row['dt_s'] and _finite(low['dt_s']) and 0. < low['dt_s'] < high['dt_s'],
            'activation location must bind selected high and real earlier low')
    width = high['dt_s']-low['dt_s']
    require(width <= tolerance and location['bracket_width_s'] == width
            and resolution['bracket_width_s'] == width
            and resolution['real_solution_bracket_width_s'] == width,
            'activation physical bracket exceeds tolerance or reports differ')
    low_pa, low_state = transaction.load(low_path); state_validator(low_state)
    require(low_state['step'] == parent_state['step']+1
            and low_state['implementation_source_sha256'] == meta['source_sha256']
            and abs(low_state['time_s']-parent_state['time_s']-low['dt_s']) <= 1e-13
            and low_state['machine_history'][:-1] == parent_state['machine_history']
            and low_state['interval_energy']['records'][:-1] == parent_state['interval_energy']['records']
            and len(low_state['machine_history']) == len(parent_state['machine_history'])+1
            and len(low_state['interval_energy']['records']) == len(parent_state['interval_energy']['records'])+1,
            'activation low was not independently advanced from same full parent')
    low_qualification = transaction.read(artifact(low_path.parent/'candidate_qualification.json', low_path.parent))
    require(low_qualification['status'] == QUALIFIED and low_qualification['accepted'] is False
            and low_qualification['restart_allowed'] is False
            and low_qualification.get('delayed_activation_pairs') == pairs,
            'real low endpoint lacks original physical qualification')
    ids = np.asarray(pa.particle_id)
    require(np.array_equal(ids, low_pa.particle_id) and np.array_equal(ids, parent_pa.particle_id)
            and np.array_equal(ids, next_pa.particle_id), 'activation material identity order differs')
    where = {int(pid): index for index, pid in enumerate(ids)}
    require(all(pid in where for pair in pairs for pid in pair), 'activation incoming pair absent from actual particles')
    i = np.array([min(where[a], where[b]) for a,b in pairs], dtype=int)
    j = np.array([max(where[a], where[b]) for a,b in pairs], dtype=int)
    lo_x, hi_x, parent_x = _positions(low_pa), _positions(pa), _positions(parent_pa)
    require(np.isfinite(lo_x).all() and np.isfinite(hi_x).all(), 'nonfinite activation positions')
    delta = float(np.max(np.linalg.norm(hi_x-lo_x, axis=1)))
    require(delta <= 5e-8 and np.isclose(delta, location['maximum_position_difference_m'], rtol=1e-12, atol=1e-18)
            and np.isclose(delta, resolution['maximum_position_difference_m'], rtol=1e-12, atol=1e-18),
            'actual activation endpoint displacement exceeds bound or saved value')
    motion = state['machine_history'][-1]; parameters = state['parameters']
    pose0 = dict(origin=motion['origin0_m'], rotation=motion['rotation0'])
    pose1 = dict(origin=motion['origin1_m'], rotation=motion['rotation1'])
    low_motion = low_state['machine_history'][-1]
    low_pose = dict(origin=low_motion['origin1_m'], rotation=low_motion['rotation1'])
    low_visibility = filter_visible_support(lo_x,i,j,np.zeros((len(i),3)),low_pose,parameters,ids)
    high_visibility = filter_visible_support(hi_x,i,j,np.zeros((len(i),3)),pose1,parameters,ids)
    require(len(low_visibility['blocked_raw_edge_indices']) == len(i)
            and len(high_visibility['blocked_raw_edge_indices']) == 0,
            'actual low/high endpoint visibility contradicts location')
    prefix_path = artifact(resolution['blocked_prefix_path'], leaf)
    require(transaction.sha(prefix_path) == resolution['blocked_prefix_sha256'], 'blocked-prefix certificate hash differs')
    certificate = transaction.read(prefix_path)
    endpoint = certificate['prefix_end_fraction']
    require(certificate['status'] == 'CERTIFIED_BLOCKED_PREFIX'
            and certificate['prefix_start_fraction'] == 0. and _finite(endpoint) and 0. < endpoint <= 1.
            and endpoint == resolution['blocked_prefix_end_fraction']
            and endpoint >= 1.-tolerance/row['dt_s']-1e-12
            and certificate['edge_count'] == len(i) and certificate['maximum_queries'] == 256
            and certificate['failure_kind'] is None and not certificate['visible_samples']
            and not certificate['degenerate_edges'] and not certificate['numerical_domain_failures']
            and not certificate['unresolved_intervals'], 'continuous blocked prefix is incomplete or outside delay bound')
    edge_ids = np.column_stack((ids[i],ids[j])).astype('<i8',copy=False)
    expected_inputs = dict(world_position0=_hash_array(parent_x),world_position1=_hash_array(hi_x),
        particle_id_order=_hash_array(ids.astype('<i8')),incoming_edge_pid_order=_hash_array(edge_ids),
        incoming_topology=_topology_hash(edge_ids),origin0=_hash_array(np.asarray(pose0['origin'],float)),
        origin1=_hash_array(np.asarray(pose1['origin'],float)),rotation0=_hash_array(np.asarray(pose0['rotation'],float)),
        rotation1=_hash_array(np.asarray(pose1['rotation'],float)),
        geometry=_hash_array(np.array([parameters['tool_radius'],.5*parameters['bore_width'],
                                      parameters['tool_top']-parameters['layer']])),
        prefix_endpoint=_hash_array(np.array([endpoint])))
    require(certificate['input_sha256'] == expected_inputs, 'blocked-prefix proof does not bind actual high chord/geometry')
    intervals = [[] for _ in i]
    for proof in certificate['certificates']:
        a,b = proof['start_fraction'],proof['end_fraction']
        require(_finite(a) and _finite(b) and 0. <= a < b <= endpoint,
                'invalid blocked-prefix certified interval')
        witnesses={w['edge_index']:w for w in proof['witnesses']}
        require(set(witnesses) == set(proof['edge_indices']), 'blocked-prefix witnesses do not bind certified edges')
        for index in proof['edge_indices']:
            require(isinstance(index,int) and 0 <= index < len(i), 'invalid blocked-prefix incoming edge index')
            witness=witnesses[index];sa,sb=witness['segment_parameter_interval'];value=witness['witness_segment_parameter']
            require(0. <= sa < value < sb <= 1., 'blocked-prefix lacks a strictly interior segment witness')
            intervals[index].append((a,b))
    for coverage in intervals:
        end=0.
        for a,b in sorted(coverage):
            require(a <= end+1e-15, 'gap in continuous blocked-prefix proof')
            end=max(end,b)
        require(end == endpoint, 'continuous blocked-prefix proof does not reach declared endpoint')
    delay = resolution['observed_high_path_activation_delay_s']
    require(resolution['activation_delay_bound_s'] == tolerance and _finite(delay) and 0. <= delay <= tolerance,
            'activation delay exceeds registered bound')
    routing = transaction.read(artifact(leaf/'union_visibility_events.json',leaf))
    require(routing['status'] == 'VISIBILITY_EVENT_REQUIRES_COUPLED_RESOLVE'
            and not routing['unresolved_samples'] and len(routing['transition_brackets']) == 1,
            'high leaf lacks its own resolved activation route')
    bracket = routing['transition_brackets'][0]
    require(bracket['direction'] == 'BLOCKED_TO_VISIBLE'
            and sorted([bracket['particle_id_i'],bracket['particle_id_j']]) == pairs[0]
            and 0. <= bracket['start_fraction'] <= bracket['end_fraction'] <= 1.
            and np.isclose(row['dt_s']*(1.-bracket['start_fraction']),delay,rtol=1e-12,atol=1e-18),
            'reported high-path activation delay or incoming identity differs')
    active = transaction.read(artifact(leaf/'visible_path_audit.json',leaf))
    require(active['status'] == 'CERTIFIED_VISIBLE_SUPPORT_PATH', 'high active support path was not certified')
    graph_path = artifact(path(next_row['checkpoint']).parent/'visible_graph.npz',path(next_row['checkpoint']).parent)
    with np.load(graph_path, allow_pickle=False) as graph:
        gi,gj,gids = graph['i'],graph['j'],graph['particle_id']
        require(np.array_equal(gids, next_pa.particle_id)
                and np.array_equal(graph['position0'],hi_x) and np.array_equal(graph['position1'],_positions(next_pa))
                and gi.dtype.kind in 'iu' and gj.dtype.kind in 'iu' and gi.shape == gj.shape
                and np.all(gi>=0) and np.all(gj<len(gids)) and np.all(gi<gj),
                'post-event graph does not bind actual selected physical states')
        restored={tuple(sorted((int(gids[a]),int(gids[b])))) for a,b in zip(gi,gj)}
        require({tuple(pair) for pair in pairs} <= restored, 'incoming pair absent from actual post-event graph')
    return resolution


def validate_gpu_screened_receipt(transaction, meta, folder, *, legacy_validator,
                                  state_validator=None, parent_fingerprint=None):
    """Dispatch old receipts unchanged; check explicit event chains in new ones.

    ``legacy_validator`` is the unbound frozen method, taking transaction,
    manifest and folder. The caller still executes the original outer receipt
    checks (status, supervision, source/reference/artifact SHA and diagnostics).
    This function does not authorize a source transition or relax a gate.
    """
    if 'event_substeps' not in meta:
        return legacy_validator(transaction, meta, folder)
    state_validator = state_validator or legacy.model.validate_state
    parent_fingerprint = parent_fingerprint or legacy.cpu.parent_fingerprint
    root = Path(transaction.root).resolve()
    folder = Path(folder).resolve()
    require(folder.is_relative_to(root), 'event receipt outside project root')
    details = meta['event_substeps']
    require(isinstance(details, dict) and details.get('schema') == SCHEMA
            and set(details) == {'schema', *ROLES}, 'registered event substep schema required')
    artifacts = meta['artifacts_sha256']

    def path(value):
        result = (root/str(value)).resolve()
        require(result.is_relative_to(root), 'event evidence outside project root')
        return result

    def artifact(value, within=None):
        result = path(value)
        if within is not None:
            require(result.is_relative_to(within), 'event evidence outside candidate role')
        key = result.relative_to(root).as_posix()
        require(key in artifacts and transaction.sha(result) == artifacts[key],
                'event evidence is absent or changed: '+key)
        return result

    def finite(value):
        return isinstance(value, numbers.Real) and not isinstance(value, bool) and math.isfinite(value)

    def valid_step(value):
        return isinstance(value, numbers.Integral) and not isinstance(value, bool) and value >= 0

    history, screen = meta['screened_history'], meta['pair_screen']
    require(meta['screened_continuation_authorized'] and len(history) == 2
            and screen['passed'] and screen['global_screen']['pass']
            and screen['local_screen']['passed'] and screen['crossing_ids_match'],
            'complete passing event fine pair required')
    require(transaction.read(artifact(folder/'pair_screen.json')) == screen,
            'event manifest and pair screen differ')
    require(screen['fine_checkpoint_sha256'] == [row['sha256'] for row in history]
            and meta['screened_continuation_checkpoint'] == history[-1]['checkpoint']
            and meta['screened_continuation_sha256'] == history[-1]['sha256'],
            'event continuation must bind exact last fine endpoint')
    require(finite(meta['dt_s']) and meta['dt_s'] > 0., 'finite positive event macro interval required')
    parent_path = path(meta['parent_checkpoint'])
    parent_pa, parent_state = transaction.load(parent_path)
    require(transaction.sha(parent_path) == meta['parent_sha256']
            and parent_fingerprint(parent_pa, parent_state) == meta['parent_fingerprint'],
            'event full parent identity changed')
    require(parent_state['time_s'] == meta['parent_time_s'], 'event parent time differs')
    state_validator(parent_state)
    endpoints = {}
    counts = {}
    seen = set()
    require(isinstance(meta.get('event_records'),dict) and set(meta['event_records']) == set(ROLES),
            'event receipt needs exact per-role resolution inventory')
    for role in ROLES:
        rows = details[role]
        require(isinstance(rows, list) and rows, 'event role needs a nonempty substep chain')
        role_dir = folder/role
        if role == 'fine_2':
            previous_path, previous_pa, previous_state, previous_sha = endpoints['fine_1']
        else:
            previous_path, previous_pa, previous_state, previous_sha = (
                parent_path, parent_pa, parent_state, meta['parent_sha256'])
        role_start = previous_state['time_s']
        contexts = []
        for row in rows:
            required = {'checkpoint', 'sha256', 'parent_checkpoint', 'parent_sha256',
                        'parent_time_s', 'parent_step', 'time_s', 'step', 'dt_s',
                        'observer_path', 'qualification_path'}
            require(isinstance(row, dict) and required.issubset(row), 'incomplete event substep row')
            current_path = artifact(row['checkpoint'], role_dir)
            require(current_path not in seen, 'reused event substep checkpoint')
            seen.add(current_path)
            require(path(row['parent_checkpoint']) == previous_path
                    and row['parent_sha256'] == previous_sha
                    and transaction.sha(path(row['parent_checkpoint'])) == previous_sha,
                    'event substep exact parent chain differs')
            require(row['parent_time_s'] == previous_state['time_s']
                    and row['parent_step'] == previous_state['step']
                    and finite(row['dt_s']) and row['dt_s'] > 0.
                    and finite(row['time_s']) and row['time_s'] > row['parent_time_s']
                    and abs(row['time_s']-row['parent_time_s']-row['dt_s']) <= 1e-13
                    and valid_step(row['step']) and row['step'] == row['parent_step']+1,
                    'event substep time or physical step is discontinuous')
            require(transaction.sha(current_path) == row['sha256'], 'event substep checkpoint hash differs')
            current_pa, current_state = transaction.load(current_path)
            state_validator(current_state)
            require(current_state['diagnostic_only'] and not current_state['accepted']
                    and not current_state['restart_allowed']
                    and current_state['implementation_source_sha256'] == meta['source_sha256']
                    and current_state['time_s'] == row['time_s']
                    and current_state['step'] == row['step'], 'event substep metadata differs')
            old_machine = previous_state['machine_history']
            new_machine = current_state['machine_history']
            old_energy = previous_state['interval_energy']['records']
            new_energy = current_state['interval_energy']['records']
            require(len(new_machine) == len(old_machine)+1 and new_machine[:-1] == old_machine
                    and len(new_energy) == len(old_energy)+1 and new_energy[:-1] == old_energy,
                    'event substep must append exactly one untouched machine/energy record')
            motion = new_machine[-1]
            require(motion['t0_s'] == row['parent_time_s'] and motion['t1_s'] == row['time_s']
                    and motion['dt_s'] == row['dt_s'], 'event substep motion differs from physical interval')
            artifact(row['observer_path'], current_path.parent)
            qualification_path = artifact(row['qualification_path'], current_path.parent)
            require(qualification_path == current_path.parent/'candidate_qualification.json',
                    'event substep qualification path differs')
            qualification = transaction.read(qualification_path)
            require(qualification['status'] == QUALIFIED and not qualification['accepted']
                    and not qualification['restart_allowed'], 'event substep is not physically qualified')
            for source_key in ('source_sha256', 'implementation_source_sha256'):
                if source_key in qualification:
                    require(qualification[source_key] == meta['source_sha256'],
                            'event substep qualification source differs')
            resolution_path = current_path.parent/'event_resolution.json'
            resolution_key = resolution_path.relative_to(root).as_posix()
            if resolution_key in artifacts or resolution_path.exists():
                artifact(resolution_path, current_path.parent)
            if 'event_resolution_path' in qualification:
                require(path(qualification['event_resolution_path']) == resolution_path,
                        'event resolution qualification path differs')
                artifact(resolution_path, current_path.parent)
            if 'event_resolution_sha256' in qualification:
                artifact(resolution_path, current_path.parent)
                require(transaction.sha(resolution_path) == qualification['event_resolution_sha256'],
                        'event resolution qualification hash differs')
            contexts.append((row,current_pa,current_state,previous_pa,previous_state,qualification,role_dir))
            previous_path, previous_pa, previous_state, previous_sha = (
                current_path, current_pa, current_state, row['sha256'])
        resolutions = []
        for index, context in enumerate(contexts):
            candidate_row, _, _, _, _, qualification, _ = context
            resolution_path = path(candidate_row['checkpoint']).parent/'event_resolution.json'
            declared = bool(qualification.get('delayed_activation_pairs'))
            present = (resolution_path.relative_to(root).as_posix() in artifacts
                       or resolution_path.exists() or 'event_resolution_path' in qualification)
            require(not declared or present, 'delayed incoming support lacks event resolution')
            if present:
                following = contexts[index+1] if index+1 < len(contexts) else None
                resolutions.append(_validate_event_resolution(transaction,meta,context,following,artifact,path,state_validator))
        require(meta['event_records'][role] == resolutions, 'manifest event inventory differs from selected leaf resolutions')
        expected_dt = meta['dt_s'] if role == 'coarse' else meta['dt_s']/2
        require(abs(previous_state['time_s']-role_start-expected_dt) <= 1e-13,
                'event substeps do not cover complete macro interval')
        alias = artifact(role_dir/'NOT_RESTARTABLE_candidate.npz', role_dir)
        require(transaction.sha(alias) == previous_sha, 'event macro endpoint is not exact last substep copy')
        # Fine_2 starts from the published fine_1 macro alias, not a hidden branch.
        endpoints[role] = (alias, previous_pa, previous_state, previous_sha)
        counts[role] = len(rows)
        if role != 'coarse':
            half = 1 if role == 'fine_1' else 2
            row = history[half-1]
            require(path(row['checkpoint']) == alias and row['sha256'] == previous_sha
                    and row['time_s'] == previous_state['time_s']
                    and abs(row['time_s']-meta['parent_time_s']-half*meta['dt_s']/2) <= 1e-13,
                    'event screened fine macro endpoint differs')
            macro_observer = artifact(row['observer_path'], role_dir)
            require(transaction.sha(macro_observer) == transaction.sha(path(rows[-1]['observer_path'])),
                    'event macro observer differs from final physical substep')
    _, _, final_state, final_sha = endpoints['fine_2']
    require(final_state['time_s'] == meta['screened_end_time_s']
            and final_state['step'] == meta['screened_end_step']
            and final_state['step'] == parent_state['step']+counts['fine_1']+counts['fine_2']
            and final_sha == meta['screened_continuation_sha256']
            and abs(final_state['time_s']-meta['parent_time_s']-meta['dt_s']) <= 1e-13,
            'event screened endpoint time/step differs')
    require(endpoints['coarse'][3] == screen['coarse_checkpoint_sha256'], 'event coarse screen binding differs')
    return dict(status='PASS_EXPLICIT_EVENT_SUBSTEP_RECEIPT', substep_counts=counts,
                fine_physical_steps=counts['fine_1']+counts['fine_2'], numerical_gates_changed=False)

"""Explicit left-side graph for a bounded visibility-activation event solve.

Only declared incoming edges blocked by the original start filter may have
their endpoint activation delayed. All raw endpoints pass the unchanged
numerical-domain checks and filter. This module does not solve, advance time,
certify a path, authorize acceptance, or persist a disabled-edge list. The
caller must use this one returned graph for density, correction, material
gradient, force and heat, and independently solve/localize/certify the event.
"""
import copy
import hashlib
import json
import numbers

import numpy as np

from visible_support import filter_visible_support, _hash_array, _topology_hash
from visible_support_remote_domain_20261007 import (
    require_visible_support_numerical_domain, require_remote_filter_consistency)


CLASSIFICATION = 'EXPLICIT_LEFT_SIDE_DELAYED_VISIBILITY_ACTIVATION_20261010_V1'
POLICY = 'ORIGINAL_ENDPOINT_FILTER_MINUS_DECLARED_START_BLOCKED_INCOMING_PAIRS'


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _incoming_pairs(value, ids, raw_pairs):
    _require(isinstance(value, (list, tuple)), 'explicit incoming PID-pair list required')
    known = set(int(pid) for pid in ids)
    raw_lookup = {tuple(sorted((int(a), int(b)))): index
        for index, (a, b) in enumerate(raw_pairs)}
    seen = set()
    result = []
    for pair in value:
        _require(isinstance(pair, (list, tuple)) and len(pair) == 2
            and all(isinstance(pid, numbers.Integral) and not isinstance(pid, (bool, np.bool_)) for pid in pair),
            'incoming pair requires two integer persistent IDs')
        first, second = (int(pid) for pid in pair)
        _require(first < second, 'incoming PID pairs must be canonical and non-self')
        key = (first, second)
        _require(first in known and second in known, 'unknown incoming persistent particle ID')
        _require(key not in seen, 'duplicate incoming PID pair')
        _require(key in raw_lookup, 'incoming PID pair absent from supplied raw support')
        seen.add(key)
        result.append((key, raw_lookup[key]))
    return sorted(result)


def restrict_activation_graph(x0, x1, raw_i, raw_j, raw_g, pose0, pose1, p,
        particle_ids, incoming_pid_pairs):
    """Return one common graph with only explicitly declared activations delayed.

    The incoming list uses canonical persistent PID pairs (PID order, not array
    indices). Every requested pair must be present in raw support and blocked
    at x0 by the original strict/roundoff filter. A requested edge that remains
    blocked at x1 is unchanged. With an empty list, every original physical
    array and original fingerprint is unchanged. Inputs are never modified.

    A fresh call with an empty mask restores the entire original endpoint graph;
    no graph state is retained in this module. Metadata distinguishes endpoint
    occlusion from an explicit, unaccepted event-side activation delay.
    """
    first, last = np.asarray(x0), np.asarray(x1)
    _require(first.shape == last.shape, 'matching start and endpoint positions required')
    start_domain = require_visible_support_numerical_domain(first, raw_i, raw_j,
        pose0, p, particle_ids)
    start_graph = filter_visible_support(first, raw_i, raw_j, raw_g, pose0, p, particle_ids)
    require_remote_filter_consistency(start_domain, start_graph)
    endpoint_domain = require_visible_support_numerical_domain(last, raw_i, raw_j,
        pose1, p, particle_ids)
    endpoint_graph = filter_visible_support(last, raw_i, raw_j, raw_g, pose1, p, particle_ids)
    require_remote_filter_consistency(endpoint_domain, endpoint_graph)
    incoming = _incoming_pairs(incoming_pid_pairs, np.asarray(particle_ids),
        endpoint_graph['raw_edge_pid_pairs'])
    start_blocked = set(int(index) for index in start_graph['blocked_raw_edge_indices'])
    _require(all(index in start_blocked for _, index in incoming),
        'incoming activation delay requires originally blocked start edge')

    endpoint_kept = endpoint_graph['kept_raw_edge_indices']
    endpoint_kept_set = set(int(index) for index in endpoint_kept)
    delayed = np.array(sorted(index for _, index in incoming if index in endpoint_kept_set), dtype=np.int64)
    keep = endpoint_kept[~np.isin(endpoint_kept, delayed)]
    blocked = np.sort(np.r_[endpoint_graph['blocked_raw_edge_indices'], delayed])
    i, j, gradients = np.asarray(raw_i), np.asarray(raw_j), np.asarray(raw_g)
    pairs = endpoint_graph['raw_edge_pid_pairs']
    selected = copy.deepcopy(endpoint_graph)
    selected.update(i=i[keep].copy(), j=j[keep].copy(), g=gradients[keep].copy(),
        kept_raw_edge_indices=keep.copy(), blocked_raw_edge_indices=blocked.copy(),
        visible_edge_pid_pairs=pairs[keep].copy(), blocked_edge_pid_pairs=pairs[blocked].copy(),
        activation_delayed_raw_edge_indices=delayed.copy())
    metadata = dict(classification=CLASSIFICATION, policy=POLICY,
        incoming_pid_pairs=[list(pair) for pair, _ in incoming],
        requested_raw_edge_indices=[int(index) for _, index in incoming],
        delayed_raw_edge_indices=delayed.tolist(),
        delayed_pid_pairs=[sorted((int(a), int(b))) for a, b in pairs[delayed]],
        delayed_edge_count=len(delayed),
        start_blocked_proven_by_original_filter=True,
        endpoint_visible_before_activation_delay=True,
        endpoint_graph_fingerprints=copy.deepcopy(endpoint_graph['fingerprints']),
        start_graph_fingerprints=copy.deepcopy(start_graph['fingerprints']),
        gradient_values_changed=False, physical_solver_evaluations=0,
        material_updates=0, time_advanced_s=0., path_certified=False,
        event_time_resolved=False, accepted=False, restart_allowed=False,
        removal_after_solve=False, persistent_disabled_edges=False)
    start_records = {row['raw_edge_index']: row for row in start_graph['diagnostics']['blocked_edges']}
    records = copy.deepcopy(endpoint_graph['diagnostics']['blocked_edges'])
    for index in delayed:
        index = int(index)
        records.append(dict(raw_edge_index=index, i=int(i[index]), j=int(j[index]),
            particle_id_i=int(pairs[index, 0]), particle_id_j=int(pairs[index, 1]),
            reason='EXPLICIT_EVENT_SIDE_ACTIVATION_DELAY_NOT_ENDPOINT_OCCLUSION',
            strict_maximum_contiguous_length_m=0., strict_total_length_m=0.,
            roundoff_expanded_maximum_contiguous_length_m=0., strict_intervals=[],
            start_filter_evidence=copy.deepcopy(start_records[index])))
    records.sort(key=lambda row: row['raw_edge_index'])
    selected['diagnostics'].update(classification=CLASSIFICATION, policy=POLICY,
        scope='DECLARED_EVENT_LEFT_SIDE_COMMON_GRAPH_REQUIRES_COUPLED_RESOLVE_AND_PATH_CERTIFICATE',
        visible_edge_count=len(keep), blocked_edge_count=len(blocked),
        endpoint_blocked_edge_count=len(endpoint_graph['blocked_raw_edge_indices']),
        activation_delayed_edge_count=len(delayed), blocked_edges=records,
        activation_delay=copy.deepcopy(metadata),
        endpoint_filter_classification=endpoint_graph['diagnostics']['classification'])
    selected['diagnostics']['limitations'] += [
        'This named event-side graph delays only declared incoming activation; it is not the original endpoint-only graph.',
        'A root must be solved again using this same graph in density, gradient, force and heat; removing an edge from a finished result is invalid.',
        'The event time and old-side/new-side path certificates are not established by this selector.',
        'Original endpoint geometric occlusion counts exclude the separately reported activation-delayed edges.']
    fingerprints = selected['fingerprints']
    fingerprints.update(visible_topology_sha256=_topology_hash(pairs[keep]),
        blocked_topology_sha256=_topology_hash(pairs[blocked]),
        visible_gradient_sha256=_hash_array(gradients[keep]))
    policy_payload = dict(classification=CLASSIFICATION, policy=POLICY,
        incoming_pid_pairs=metadata['incoming_pid_pairs'],
        start_world_positions_sha256=start_graph['fingerprints']['world_positions_sha256'],
        start_tool_origin_sha256=start_graph['fingerprints']['tool_origin_sha256'],
        start_tool_rotation_sha256=start_graph['fingerprints']['tool_rotation_sha256'],
        endpoint_geometry_policy_sha256=endpoint_graph['fingerprints']['geometry_policy_sha256'])
    fingerprints['activation_policy_sha256'] = hashlib.sha256(
        json.dumps(policy_payload, sort_keys=True, separators=(',', ':')).encode('utf-8')).hexdigest()
    selected.update(endpoint_graph=copy.deepcopy(endpoint_graph), activation_delay=metadata,
        start_numerical_domain=copy.deepcopy(start_domain),
        endpoint_numerical_domain=copy.deepcopy(endpoint_domain))
    _require(np.array_equal(np.sort(np.r_[keep, blocked]), np.arange(len(i)))
        and len(set(blocked.tolist())) == len(blocked), 'common graph must partition all raw support')
    _require(set(endpoint_kept.tolist())-set(keep.tolist()) == set(delayed.tolist()),
        'graph difference must be exactly the declared endpoint-visible incoming set')
    return selected

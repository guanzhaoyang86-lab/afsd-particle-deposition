"""Independent readback of CPU event branches and same-parent sensitivity.

No physical solve is performed here. The frozen explicit-event receipt reader
and the original three-endpoint sensitivity comparisons are reused unchanged.
Administrative publication and complete ancestor authorization belong to the
caller; partial branch results never authorize a continuation.
"""
from pathlib import Path
import itertools
import math
import numpy as np

import afsd_deposition_remote_setup_20261007 as cpu_model
import run_afsd_deposition_stage_remote_20261007 as cpu
import afsd_event_receipt_20261010 as event_receipt
import run_gpu_event_delay_campaign_20261010 as sensitivity

ROOT = cpu.ROOT
ROLES = ('coarse', 'fine_1', 'fine_2')
SCALES = (1., .5, .25)
require, read, load, rel = cpu.require, cpu.read, cpu.load, cpu.rel
sha = cpu.core.sha


def path(value):
    result = (ROOT / str(value)).resolve()
    require(result.is_relative_to(ROOT), 'CPU event evidence outside project')
    return result


def positive(value):
    return type(value) in (int, float) and math.isfinite(value) and value > 0.


def verify_map(mapping, within=None):
    require(isinstance(mapping, dict) and bool(mapping), 'nonempty SHA evidence map required')
    for name, digest in mapping.items():
        target = path(name)
        require(rel(target) == name and isinstance(digest, str) and len(digest) == 64,
                'noncanonical evidence path or digest')
        require(within is None or target.is_relative_to(within), 'artifact outside its branch')
        require(sha(target) == digest, 'changed CPU event evidence: ' + name)


class Transaction:
    root = ROOT

    @staticmethod
    def sha(value):
        return sha(path(value))

    @staticmethod
    def read(value):
        return read(path(value))

    @staticmethod
    def load(value):
        return load(path(value))


def _legacy_forbidden(*args):
    raise ValueError('CPU event publication requires explicit physical substeps')


def _supervision(folder, *, protocol_required):
    result = read(folder / 'supervision.json')
    require(result['exit_code'] == 0 and result['timed_out'] is False
            and result['source_inputs_still_match'] is True,
            'CPU event child or aggregate did not complete unchanged')
    require(result['manifest_sha256'] == sha(folder / 'manifest.json'),
            'CPU event supervision does not bind manifest')
    if protocol_required:
        require(result['protocol_sha256'] == sha(folder / 'protocol.json'),
                'CPU event supervision does not bind protocol')
    require(positive(result['wall_seconds']), 'positive supervised wall time required')
    return result


def _branch(folder, role, meta, expected_parent, expected_sha):
    branch = folder / role
    protocol = read(branch / 'protocol.json')
    record = read(branch / 'manifest.json')
    supervised = _supervision(branch, protocol_required=True)
    expected_dt = meta['dt_s'] if role == 'coarse' else meta['dt_s'] / 2.
    identities = dict(role=role, event_tolerance_scale=meta['event_tolerance_scale'],
        dt_s=meta['dt_s'], role_dt_s=expected_dt,
        parent_checkpoint=meta['parent_checkpoint'], parent_sha256=meta['parent_sha256'],
        role_parent_checkpoint=rel(expected_parent), role_parent_sha256=expected_sha,
        source_sha256=meta['source_sha256'])
    for name, expected in identities.items():
        require(protocol.get(name) == expected and record.get(name) == expected,
                'CPU event branch identity mismatch: ' + role + '/' + name)
    require(record['status'] == 'PASS_CPU_EVENT_BRANCH_DIAGNOSTIC_ONLY'
            and record['accepted'] is False and record['restart_allowed'] is False,
            'CPU event branch must remain qualified diagnostic evidence')
    require(record['leaves'] == meta['event_substeps'][role]
            and record['event_records'] == meta['event_records'][role],
            'CPU event branch leaf/resolution inventory differs')
    require(len(record['event_records']) <= 1,
            'CPU bounded branch cannot exceed its registered single-event scope')
    verify_map(record['artifacts_sha256'], branch)
    for name, digest in record['artifacts_sha256'].items():
        require(meta['artifacts_sha256'].get(name) == digest,
                'aggregate omitted or changed a branch artifact')
    for name in ('protocol.json', 'manifest.json', 'supervision.json'):
        target = branch / name
        require(meta['artifacts_sha256'].get(rel(target)) == sha(target),
                'aggregate does not bind complete supervised branch')
    endpoint = branch / 'NOT_RESTARTABLE_candidate.npz'
    require(path(record['checkpoint']) == endpoint and record['sha256'] == sha(endpoint),
            'branch endpoint does not bind exact published alias')
    pa, state = load(endpoint)
    cpu_model.validate_state(state)
    require(state['time_s'] == record['state_time_s'] and state['step'] == record['state_step'],
            'branch saved endpoint time/step differs')
    require(positive(record['wall_seconds']) and record['wall_seconds'] <= supervised['wall_seconds'] + 1.,
            'branch execution cost exceeds supervised wall time')
    solves = record['solves']
    journals = [name for name in record['artifacts_sha256'] if name.endswith('/physical_solve_cost.json')]
    require(type(solves) is int and 1 <= solves <= 12 and len(journals) == solves,
            'bounded branch must bind every completed physical solve cost')
    seen = set()
    for name in journals:
        journal = read(path(name)); target = path(journal['checkpoint'])
        require(target == path(name).parent / 'NOT_RESTARTABLE_candidate.npz'
                and target not in seen and record['artifacts_sha256'].get(rel(target)) == sha(target)
                and journal['checkpoint_sha256'] == sha(target),
                'physical solve journal does not bind unique saved state')
        seen.add(target)
        require(journal['status'] == 'MEASURED_COMPLETED_PHYSICAL_SOLVE'
                and journal['accepted'] is False and positive(journal['wall_seconds'])
                and positive(journal['dt_s']) and journal['dt_s'] <= expected_dt * (1. + 1e-12)
                and journal['scope'] == 'NONLINEAR_SOLVE_FRESH_ADVANCE_AND_SAVED_STATE_EXCLUDES_LATER_GEOMETRY_QUALIFICATION',
                'invalid complete physical solve cost journal')
    return dict(manifest=record, endpoint=(pa, state), supervision=supervised)


def verify_pair(folder, expected_source=None):
    """Fully read one complete aggregate; does not authorize a partial group."""
    import run_cpu_visibility_event_20261018 as physical
    folder = path(folder)
    meta = read(folder / 'manifest.json')
    _supervision(folder, protocol_required=False)
    require(meta['status'] == physical.CONTINUE_PASS
            and meta['screened_continuation_authorized'] is True
            and meta['accepted'] is False and meta['restart_allowed'] is False,
            'complete diagnostic CPU event pair required')
    require(positive(meta['dt_s']) and meta['event_tolerance_scale'] in SCALES,
            'registered positive macro interval and tolerance scale required')
    require(meta['event_protocol'] == physical.PROTOCOL
            and meta['event_locator_protocol'] == physical.LOCATOR_PROTOCOL
            and meta['blocked_prefix_protocol'] == physical.BLOCKED_PREFIX_PROTOCOL
            and meta['visible_path_protocol'] == physical.VISIBLE_PATH_PROTOCOL,
            'CPU event physical protocols differ')
    if expected_source is not None:
        require(meta['source_sha256'] == expected_source, 'unexpected CPU event implementation')
    verify_map(meta['source_sha256'])
    require(meta['source_sha256'].get(rel(Path(__file__).resolve())) == sha(Path(__file__)),
            'CPU event source map must bind this exact receipt reader')
    verify_map(meta['artifacts_sha256'], folder)
    parent = path(meta['parent_checkpoint'])
    require(sha(parent) == meta['parent_sha256'], 'CPU event actual common parent changed')
    pa, state = load(parent)
    cpu_model.validate_state(state)
    inherited = state['implementation_source_sha256']
    require(inherited == meta['parent_implementation_source_sha256']
            and all(meta['source_sha256'].get(k) == v for k, v in inherited.items()),
            'CPU source transition must retain the exact original implementation map')
    branches = {}
    for role in ROLES:
        previous = folder / 'fine_1/NOT_RESTARTABLE_candidate.npz' if role == 'fine_2' else parent
        branches[role] = _branch(folder, role, meta, previous, sha(previous))
    receipt = event_receipt.validate_gpu_screened_receipt(Transaction(), meta, folder,
        legacy_validator=_legacy_forbidden, state_validator=cpu_model.validate_state,
        parent_fingerprint=cpu.parent_fingerprint)
    _actual_parent_bindings(meta)
    # Reloaded particle states, not manifest booleans, determine the pair gates.
    coarse, _ = branches['coarse']['endpoint']
    fine, fine_state = branches['fine_2']['endpoint']
    coarse_state = branches['coarse']['endpoint'][1]
    glob = cpu.post.compare(fine, coarse, state['parameters'])
    local = cpu.event.local_compare(fine, coarse, sorted(int(v) for v in pa.particle_id))
    matching = fine_state['bore_observation']['crossing_ids'] == coarse_state['bore_observation']['crossing_ids']
    screen = meta['pair_screen']
    require(glob == screen['global_screen'] and local == screen['local_screen']
            and matching == screen['crossing_ids_match']
            and glob['pass'] and local['passed'] and matching,
            'saved CPU pair screen differs from recomputed original gates')
    for branch in branches.values():
        other, other_state = branch['endpoint']
        require(np.array_equal(other.particle_id, pa.particle_id) and np.array_equal(other.m, pa.m)
                and other_state['parameters'] == state['parameters'],
                'CPU event endpoint changed inventory or physical parameters')
    require(meta['physical_solve_count'] == sum(b['manifest']['solves'] for b in branches.values()),
            'aggregate physical solve inventory differs from bounded branches')
    return dict(status='PASS_CPU_EVENT_PAIR_READBACK', folder=rel(folder),
        manifest_sha256=sha(folder / 'manifest.json'), receipt=receipt,
        original_pair_gates_recomputed=True, source_transition_verified=True,
        selected_checkpoint=meta['screened_continuation_checkpoint'],
        selected_sha256=meta['screened_continuation_sha256'], meta=meta)


def _actual_parent_bindings(meta):
    """Check the explicit CPU stamps in addition to the old geometric proof."""
    artifacts = meta['artifacts_sha256']
    for role in ROLES:
        for row in meta['event_substeps'][role]:
            previous = path(row['parent_checkpoint'])
            previous_sha = row['parent_sha256']
            parent_pa, parent_state = load(previous)
            fingerprint = cpu.parent_fingerprint(parent_pa, parent_state)
            leaf = path(row['checkpoint']).parent
            reports = [leaf / 'visible_path_audit.json']
            resolution = leaf / 'event_resolution.json'
            if resolution.exists():
                event = read(resolution)
                reports.append(path(event['blocked_prefix_path']))
                reports.append(path(event['low_checkpoint']).parent / 'visible_path_audit.json')
                location = read(path(event['location']))
                require(location.get('parent_sha256') == previous_sha
                        and location.get('legacy_fixture_parent_is_current') is False
                        and location.get('parent_binding') == 'CALLER_VERIFIED_ACTUAL_CHECKPOINT_SHA256'
                        and location.get('inherited_locator_algorithm_unchanged') is True,
                        'CPU locator does not bind the selected actual parent')
            for target in reports:
                require(artifacts.get(rel(target)) == sha(target),
                        'CPU actual-parent report is not bound by the pair')
                report = read(target)
                require(path(report['parent_checkpoint']) == previous
                        and report['parent_sha256'] == previous_sha
                        and report['actual_parent_fingerprint'] == fingerprint
                        and report['legacy_fixture_parent_is_current'] is False
                        and report['parent_binding'] == 'CALLER_VERIFIED_ACTUAL_CHECKPOINT_SHA256'
                        and report['inherited_geometry_algorithm_unchanged'] is True,
                        'CPU geometric proof actual-parent stamp differs')


def compare_group(folders, expected_source=None):
    """Read all three siblings then reuse original endpoint/timing thresholds."""
    values = list(folders.values()) if isinstance(folders, dict) else list(folders)
    require(len(values) == 3, 'exactly three registered tolerance siblings required')
    verified = [verify_pair(value, expected_source) for value in values]
    stages = {}
    for result in verified:
        meta = result['meta']; scale = meta['event_tolerance_scale']
        require(scale not in stages, 'duplicate CPU tolerance scale')
        stages[scale] = (path(result['folder']), meta)
    require(set(stages) == set(SCALES), 'missing registered CPU event tolerance sibling')
    baseline = stages[1.][1]
    parent_pa, parent_state = load(path(baseline['parent_checkpoint']))
    identities = ('parent_checkpoint', 'parent_sha256', 'parent_fingerprint',
        'parent_time_s', 'parent_implementation_source_sha256', 'source_sha256', 'dt_s')
    endpoints, timing = {}, {}
    for scale in SCALES:
        folder, meta = stages[scale]
        require(all(meta[key] == baseline[key] for key in identities),
                'CPU tolerance variants must use identical complete parent/source/macro interval')
        require(any(meta['event_records'][role] for role in ROLES),
                'CPU tolerance qualification must contain an actual resolved event')
        endpoints[scale] = {role: load(folder / role / 'NOT_RESTARTABLE_candidate.npz') for role in ROLES}
        timing[scale] = sensitivity.event_timing(meta, reader=read, loader=load)
    comparisons = []
    for a, b in itertools.combinations(SCALES, 2):
        screens = {role: sensitivity.compare_endpoint(endpoints[b][role], endpoints[a][role],
            parent_pa.particle_id, parent_pa.m, parent_state['parameters']) for role in ROLES}
        events = sensitivity.compare_event_timing(timing[a], timing[b])
        comparisons.append(dict(scales=[a, b], macro_endpoint_screens=screens, event_timing=events))
    return dict(status='PASS_SAME_PARENT_CPU_EVENT_TOLERANCE_SENSITIVITY',
        common_parent_checkpoint=baseline['parent_checkpoint'], common_parent_sha256=baseline['parent_sha256'],
        dt_s=baseline['dt_s'], stages=[dict(scale=scale, directory=rel(stages[scale][0]),
            manifest_sha256=sha(stages[scale][0] / 'manifest.json'),
            checkpoint_sha256=stages[scale][1]['screened_continuation_sha256']) for scale in SCALES],
        selected_scale=.25, selected_stage=rel(stages[.25][0]),
        selected_checkpoint=stages[.25][1]['screened_continuation_checkpoint'],
        selected_sha256=stages[.25][1]['screened_continuation_sha256'],
        comparisons=comparisons, event_timing=timing, accepted=False,
        original_endpoint_and_event_timing_thresholds_reused=True)

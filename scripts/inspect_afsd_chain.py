"""Read accepted AFSD restart segments without advancing or modifying physics.

Trusted project checkpoints only: upstream NPZ loading uses pickle.  This audits
saved evidence, not constitutive correctness or restart/discretization accuracy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from afsd_backend.pysph.output import restore
from afsd_scaffold.vtk_legacy import read_point_cloud

IDENTITY = ('particle_id', 'material_id', 'component_id', 'tracer_type',
            'wholebar_tag', 'actuator')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    def invalid(value):
        raise ValueError('nonfinite JSON number: ' + value)
    return json.loads(Path(path).read_text(encoding='utf-8'), parse_constant=invalid)


def finite(value, label):
    if isinstance(value, dict):
        for name, child in value.items():
            finite(child, label + '.' + name)
    elif isinstance(value, (list, tuple)):
        for child in value:
            finite(child, label)
    elif isinstance(value, (int, float)) and not np.isfinite(value):
        raise ValueError('nonfinite ' + label)


def ordered_frames(frames):
    if not frames:
        raise ValueError('empty accepted series')
    for frame in frames:
        finite(frame, 'frame')
        if frame['step'] < 0 or int(frame['step']) != frame['step'] or frame['time_s'] < 0:
            raise ValueError('invalid accepted step/time')
        if 'REJECTED' in frame['file'].upper():
            raise ValueError('REJECTED frame in accepted series')
    if any(b['step'] <= a['step'] or b['time_s'] <= a['time_s']
           for a, b in zip(frames, frames[1:])):
        raise ValueError('duplicate or nonmonotone accepted step/time')
    return frames


def id_keyed(data):
    """Canonicalize *all* saved point fields, with exact integral ID checks."""
    if 'REJECTED' in data.get('title', '').upper():
        raise ValueError('REJECTED data in accepted series')
    ids = np.asarray(data['scalars']['particle_id'])
    if not np.isfinite(ids).all() or not np.equal(ids, np.floor(ids)).all():
        raise ValueError('noninteger/nonfinite persistent ID')
    if len(np.unique(ids)) != len(ids):
        raise ValueError('duplicate persistent ID')
    order = np.argsort(ids)
    fields = {'position_m': np.asarray(data['points'])}
    for group in ('scalars', 'vectors'):
        for name, values in data[group].items():
            if name in fields:
                raise ValueError('duplicate output field ' + name)
            fields[name] = np.asarray(values)
    for name, values in fields.items():
        if len(values) != len(ids) or not np.isfinite(values).all():
            raise ValueError('nonfinite or inconsistent output field ' + name)
        fields[name] = values[order]
    if fields['position_m'].shape != (len(ids), 3):
        raise ValueError('position must be genuinely three-component data')
    return fields


def exact_fields(a, b, names=None):
    """No tolerance/interpolation at a checkpoint restart boundary."""
    all_fields = names is None
    names = set(a) if all_fields else set(names)
    if not names <= set(a) or not names <= set(b) or (all_fields and set(a) != set(b)):
        raise ValueError('restart field inventory mismatch')
    for name in sorted(names):
        if not np.array_equal(a[name], b[name]):
            raise ValueError('ID-keyed field mismatch: ' + name)


def validate_metadata(sources, config, state, root=ROOT):
    for path, expected in sources.items():
        actual_path = root / Path(path.replace('\\', '/'))
        if sha(actual_path) != expected:
            raise ValueError('source digest mismatch: ' + path)
    fingerprint = hashlib.sha256(json.dumps(config['parameters'], sort_keys=True).encode()).hexdigest()
    if state.get('source_hashes') != sources:
        raise ValueError('checkpoint source metadata mismatch')
    if state.get('configuration_fingerprint') != fingerprint:
        raise ValueError('checkpoint configuration fingerprint mismatch')
    return fingerprint


def checkpoint_fields(pa, data):
    """Check saved primary state, rates and stress against final VTK by ID."""
    order = np.argsort(pa.particle_id)
    properties = {name: getattr(pa, name)[order] for name in IDENTITY}
    for prop, name in (('rho', 'density_kg_per_m3'), ('temperature', 'temperature_K'),
                       ('ep', 'equivalent_plastic_strain'), ('m', 'mass_kg'),
                       ('h', 'smoothing_length_m')):
        properties[name] = getattr(pa, prop)[order]
    for props, name in ((('x', 'y', 'z'), 'position_m'),
                         (('u', 'v', 'w'), 'velocity_m_per_s'),
                         (('x0', 'y0', 'z0'), 'initial_position_m')):
        properties[name] = np.column_stack([getattr(pa, p) for p in props])[order]
    exact_fields(properties, data, list(properties))
    stress = pa.S.reshape(-1, 3, 3)[order] - data['pressure_Pa'][:, None, None] * np.eye(3)
    for a, b, suffix in ((0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'),
                         (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')):
        if not np.array_equal(stress[:, a, b], data['stress_' + suffix + '_Pa']):
            raise ValueError('checkpoint stress mismatch: ' + suffix)
    if 'equivalent_deviatoric_strain_rate_per_s' in data:
        if not np.array_equal(pa.strain_rate[order], data['equivalent_deviatoric_strain_rate_per_s']):
            raise ValueError('checkpoint strain rate mismatch')
        if not np.array_equal(pa.ep_rate[order], data['equivalent_plastic_strain_rate_per_s']):
            raise ValueError('checkpoint plastic rate mismatch')
        strain_rate = pa.D.reshape(-1, 3, 3)[order]
        for a, b, suffix in ((0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'),
                             (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')):
            if not np.array_equal(strain_rate[:, a, b], data['strain_rate_' + suffix + '_per_s']):
                raise ValueError('checkpoint rate tensor mismatch: ' + suffix)


def load_segment(directory, root=ROOT):
    directory = Path(directory).resolve()
    manifest = read_json(directory / 'manifest.json')
    performance = read_json(directory / 'performance.json')
    sources = read_json(directory / 'source_hashes.json')
    config = read_json(directory / 'resolved_config.json')
    frames = ordered_frames(read_json(directory / 'series.json')['frames'])
    if sha(directory / 'resolved_config.json') != manifest['configuration']['sha256']:
        raise ValueError('configuration artifact hash mismatch')
    if sha(directory / 'checkpoint.npz') != manifest['checkpoint']['sha256']:
        raise ValueError('checkpoint artifact hash mismatch; incomplete manifest needs explicit recovery')
    pa, state = restore(directory / 'checkpoint.npz')
    finite(state, 'checkpoint state')
    for name, prop in pa.properties.items():
        if not np.isfinite(prop.get_npy_array()).all():
            raise ValueError('nonfinite checkpoint property ' + name)
    fingerprint = validate_metadata(sources, config, state, root)
    if state['count'] != frames[-1]['step'] or state['t'] != frames[-1]['time_s']:
        raise ValueError('checkpoint/accepted endpoint mismatch')
    if manifest['checkpoint']['step'] != frames[-1]['step']:
        raise ValueError('manifest/accepted endpoint mismatch')
    if state['contact_variant'] != manifest['provenance']['contact_variant']:
        raise ValueError('contact variant metadata mismatch')
    first = last = None
    for frame in frames:
        path = (directory / frame['file']).resolve()
        if path.parent != directory or sha(path) != frame['sha256']:
            raise ValueError('frame location/hash mismatch')
        data = id_keyed(read_point_cloud(path))
        if first is None:
            first = data
        exact_fields(first, data, IDENTITY + ('initial_position_m', 'mass_kg'))
        last = data
    checkpoint_fields(pa, last)
    diagnostics = []
    for row in (directory / 'diagnostics.jsonl').read_text(encoding='utf-8').splitlines():
        if row.strip():
            value = json.loads(row)
            finite(value, 'diagnostics')
            diagnostics.append(value)
    keys = [d['step'] for d in diagnostics]
    if keys != sorted(set(keys)):
        raise ValueError('duplicate or nonmonotone diagnostic step')
    by_step = {d['step']: d for d in diagnostics}
    if any(row['alarms'] for row in diagnostics if row['step'] <= frames[-1]['step']):
        raise ValueError('alarmed diagnostic inside accepted segment')
    # Restart initial frame has no newly evaluated diagnostic row; parent covers it.
    for frame in frames[1:]:
        row = by_step.get(frame['step'])
        if row is None or row['time_s'] != frame['time_s'] or row['alarms']:
            raise ValueError('missing/misaligned/alarmed accepted diagnostics')
    rejected = []
    for path in sorted(directory.glob('REJECTED*.npz')):
        try:
            restore(path)
        except ValueError as error:
            if 'REJECTED candidate' not in str(error):
                raise
        else:
            raise ValueError('REJECTED candidate can be restored')
        rejected.append({'file': str(path), 'sha256': sha(path), 'restart_allowed': False})
    return dict(directory=directory, manifest=manifest, performance=performance,
                sources=sources, config=config, fingerprint=fingerprint, state=state,
                frames=frames, first=first, last=last, diagnostics=diagnostics,
                rejected=rejected, checkpoint_sha256=sha(directory / 'checkpoint.npz'))


def join_segments(segments):
    """Require exact lineage, then omit only the verified shared boundary frame."""
    combined = []
    for index, segment in enumerate(segments):
        if index:
            parent = segments[index - 1]
            checkpoint_parent = segment['manifest']['provenance'].get('checkpoint_parent')
            if checkpoint_parent is None or Path(checkpoint_parent).resolve() != parent['directory'] / 'checkpoint.npz':
                raise ValueError('broken checkpoint parent lineage')
            if segment['sources'] != parent['sources'] or segment['config'] != parent['config']:
                raise ValueError('restart source/config mismatch')
            if segment['state']['contact_variant'] != parent['state']['contact_variant']:
                raise ValueError('restart contact variant mismatch')
            a, b = parent['frames'][-1], segment['frames'][0]
            if (a['step'], a['time_s']) != (b['step'], b['time_s']):
                raise ValueError('restart boundary step/time mismatch')
            for field in ('tool_translation_m', 'tool_rotation_rad', 'tool_axis',
                          'rate_evaluation_time_s', 'unavailable_fields_reason'):
                if a.get(field) != b.get(field):
                    raise ValueError('restart boundary metadata mismatch: ' + field)
            exact_fields(parent['last'], segment['first'])
        for frame in segment['frames'][1 if index else 0:]:
            combined.append({**frame, 'run_directory': str(segment['directory'])})
    ordered_frames(combined)
    return combined


def observation(start, end):
    free = end['actuator'] == 0
    delta = end['position_m'] - start['position_m']
    result = dict(displacement_max_m=float(np.linalg.norm(delta, axis=1).max()),
                  free_displacement_max_m=float(np.linalg.norm(delta[free], axis=1).max()) if free.any() else None,
                  free_speed_max_m_s=float(np.linalg.norm(end['velocity_m_per_s'][free], axis=1).max()) if free.any() else None,
                  maximum_equivalent_plastic_strain=float(end['equivalent_plastic_strain'].max()),
                  plastified_particle_count=int(np.count_nonzero(end['equivalent_plastic_strain'] > 0)),
                  temperature_min_K=float(end['temperature_K'].min()),
                  temperature_max_K=float(end['temperature_K'].max()))
    return result


def validate_parent_receipt(path, parent, target_end_time_s):
    """Check a separately saved pre-run receipt without inventing missing metadata."""
    receipt = read_json(path)
    if receipt.get('scope') != 'BEFORE_RUN_PARENT_RECEIPT':
        raise ValueError('receipt is not identified as a before-run parent receipt')
    if Path(receipt['checkpoint_path']).resolve() != parent['directory']/'checkpoint.npz':
        raise ValueError('receipt parent path mismatch')
    for key, file in (('checkpoint_sha256', 'checkpoint.npz'), ('manifest_sha256', 'manifest.json'),
                       ('series_sha256', 'series.json'), ('source_hashes_sha256', 'source_hashes.json'),
                       ('configuration_sha256', 'resolved_config.json')):
        if receipt[key] != sha(parent['directory']/file):
            raise ValueError('parent receipt digest mismatch: ' + key)
    final = parent['frames'][-1]
    if receipt['accepted_step'] != final['step'] or receipt['accepted_time_s'] != final['time_s']:
        raise ValueError('parent receipt accepted state mismatch')
    if receipt['requested_end_time_s'] != target_end_time_s:
        raise ValueError('parent receipt target changed')
    return dict(path=str(Path(path).resolve()), sha256=sha(path),
                recorded_utc=receipt['recorded_utc'], status='PASS',
                parent_checkpoint_sha256=receipt['checkpoint_sha256'])


def accepted_target_complete(endpoint_time_s, target_time_s, runner_status, manifest_status):
    """An accepted timestamp alone cannot reclassify FAILED/INCOMPLETE output."""
    return (endpoint_time_s >= target_time_s and runner_status == 'COMPLETE'
            and manifest_status == 'COMPLETE')


def inspect_chain(directories, target_end_time_s, root=ROOT, parent_receipts=None):
    if not directories or not np.isfinite(target_end_time_s) or target_end_time_s <= 0:
        raise ValueError('segments and explicit positive finite physical target required')
    segments = [load_segment(path, root) for path in directories]
    combined = join_segments(segments)
    if parent_receipts is not None and len(parent_receipts) != len(segments)-1:
        raise ValueError('provide one before-run receipt for each parent-child link')
    receipts = [validate_parent_receipt(path, parent, target_end_time_s)
                for path, parent in zip(parent_receipts or [], segments[:-1])]
    first, last = segments[0], segments[-1]
    summaries = []
    for index, segment in enumerate(segments):
        start, end = segment['frames'][0], segment['frames'][-1]
        summaries.append(dict(directory=str(segment['directory']),
            runner_status=segment['performance']['status'], start_step=start['step'], end_step=end['step'],
            start_time_s=start['time_s'], end_time_s=end['time_s'], accepted_duration_s=end['time_s']-start['time_s'],
            accepted_additional_steps=end['step']-start['step'],
            checkpoint_sha256=segment['checkpoint_sha256'], configuration_fingerprint=segment['fingerprint'],
            observed_segment=observation(segment['first'], segment['last']),
            cumulative_ledger_at_endpoint=segment['state']['ledger'],
            segment_ledger_delta={k: v-segments[index-1]['state']['ledger'][k]
                                  for k, v in segment['state']['ledger'].items()} if index else None,
            rejected_candidates_excluded=segment['rejected'],
            evidence_sha256={name: sha(segment['directory']/name) for name in
                ('series.json', 'manifest.json', 'performance.json', 'source_hashes.json', 'resolved_config.json')},
            performance=segment['performance']))
    reached = accepted_target_complete(combined[-1]['time_s'], target_end_time_s,
                                       last['performance']['status'], last['manifest']['status'])
    origin_available = combined[0]['step'] == 0 and combined[0]['time_s'] == 0
    return dict(readback_status='PASS', status='ACCEPTED_TARGET_REACHED' if reached else 'VALID_PARTIAL_CHAIN',
        target_end_time_s=target_end_time_s, reached_target=reached,
        accepted_start_time_s=combined[0]['time_s'], accepted_end_time_s=combined[-1]['time_s'],
        accepted_end_step=combined[-1]['step'], accepted_unique_frames=len(combined),
        duplicated_boundary_frames_removed=len(segments)-1, origin_frame_available=origin_available,
        observed_chain_start_to_end=observation(first['first'], last['last']),
        observed_origin_to_end=observation(first['first'], last['last']) if origin_available else None,
        before_run_parent_receipts=receipts,
        segments=summaries, composed_series=combined,
        paper_reproduction=False, deposition_validated=False, mixing_validated=False, production_authorized=False,
        limitations=[
            'Saved accepted evidence only; no solver advancement, interpolation, or restart-equivalence proof.',
            'The existing runner does not embed an input-time parent digest. Any supplied separate before-run receipt is checked; otherwise only readback parent hash, manifest, lineage, source/config and boundary fields are checked.',
            'Child initial checkpoint internal/scratch state was not separately archived; exact boundary comparison covers all exported particle fields, not unsaved internal state.',
            'No state here establishes macroscopic deposition or quantitative paper reproduction.',
        ])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='append', required=True, type=Path,
                        help='accepted segments in parent-to-child order; repeat per segment')
    parser.add_argument('--target-end-time', required=True, type=float,
                        help='absolute SI physical target; not silently changed for partial results')
    parser.add_argument('--parent-receipt', action='append', type=Path,
                        help='optional separately archived before-run receipt, one per parent-child link')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = inspect_chain(args.run, args.target_end_time, parent_receipts=args.parent_receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'composed_series'}, indent=2, allow_nan=False))
    raise SystemExit(0 if report['reached_target'] else 2)

"""Read-only process observables and cost scales; not a deposition certificate."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud
from afsd_backend.pysph.output import restore
from inspect_afsd_chain import inspect_chain

EP_OBSERVATION_FLOOR = 1e-10  # observation only; never used by mechanics


def observe(data, parameters, time_s):
    scalar = {k: np.asarray(v) for k, v in data['scalars'].items()}
    x = np.asarray(data['points'], float)
    origin = np.asarray(data['vectors']['initial_position_m'], float)
    velocity = np.asarray(data['vectors']['velocity_m_per_s'], float)
    required = ('particle_id', 'component_id', 'actuator', 'mass_kg',
                'density_kg_per_m3', 'equivalent_plastic_strain')
    if x.ndim != 2 or x.shape[1] != 3 or origin.shape != x.shape or velocity.shape != x.shape:
        raise ValueError('3D positions, original positions and velocity required')
    if any(k not in scalar or scalar[k].shape != (len(x),) for k in required):
        raise ValueError('missing/malformed solved particle fields')
    if not all(np.isfinite(v).all() for v in (x, origin, velocity, *(scalar[k] for k in required))):
        raise ValueError('nonfinite observed fields')
    if np.any(scalar['mass_kg'] <= 0) or np.any(scalar['density_kg_per_m3'] <= 0):
        raise ValueError('positive mass/density required, not clipped')
    if len(np.unique(scalar['particle_id'])) != len(x):
        raise ValueError('duplicate persistent IDs')
    free = scalar['actuator'] == 0
    bar = scalar['component_id'] == 1
    ep = scalar['equivalent_plastic_strain']
    theta = parameters['omega']*time_s
    xc, yc = x[:, 0]-parameters['traverse']*time_s, x[:, 1]
    # LAB -> rotating bore coordinates, for region query ONLY.
    local_x = xc*math.cos(theta)+yc*math.sin(theta)
    local_y = -xc*math.sin(theta)+yc*math.cos(theta)
    outside_bore = np.maximum(np.abs(local_x), np.abs(local_y)) > parameters['bore_width']/2
    candidate = bar & free & outside_bore & (x[:, 2] > 0) & (x[:, 2] <= parameters['layer'])
    net_crossed = bar & (origin[:, 2] >= parameters['layer']) & (x[:, 2] < parameters['layer'])
    plastic = free & (ep > EP_OBSERVATION_FLOOR)
    return dict(time_s=float(time_s), free_plastic_raw_positive_count=int(np.count_nonzero(free & (ep > 0))),
        free_plastic_above_floor_count=int(plastic.sum()),
        free_bar_plastic_count=int(np.count_nonzero(plastic & bar)),
        free_substrate_plastic_count=int(np.count_nonzero(plastic & ~bar)),
        free_plastic_particle_ids=sorted(scalar['particle_id'][plastic].astype(int).tolist()),
        maximum_ep=float(ep.max()), negative_ep_count=int(np.count_nonzero(ep < 0)),
        free_max_displacement_from_origin_m=float(np.linalg.norm(x[free]-origin[free], axis=1).max()) if free.any() else None,
        net_initially_above_layer_currently_below_mass_kg=float(scalar['mass_kg'][net_crossed].sum()),
        outside_bore_layer_candidate_particles=int(candidate.sum()),
        outside_bore_layer_candidate_mass_kg=float(scalar['mass_kg'][candidate].sum()),
        outside_bore_layer_candidate_volume_m3=float(np.sum(scalar['mass_kg'][candidate]/scalar['density_kg_per_m3'][candidate])),
        prescribed_feed_m=float(parameters['feed']*time_s),
        tool_rotation_rad=float(theta),
        sustained_deposition_validated=False)


def process_costs(p, dt_s, wall_per_step_s):
    if dt_s <= 0 or wall_per_step_s <= 0 or not np.isfinite([dt_s, wall_per_step_s]).all():
        raise ValueError('positive finite measured time scales required')
    if p['dx'] <= 0 or p['layer'] <= 0:
        raise ValueError('positive geometry lengths required')
    events = [('tool_1degree', math.pi/180, abs(p['omega']), 'rad'),
              ('tool_1revolution', 2*math.pi, abs(p['omega']), 'rad'),
              ('feed_0p1dx', .1*p['dx'], abs(p['feed']), 'm'),
              ('feed_1dx', p['dx'], abs(p['feed']), 'm'),
              ('feed_1layer', p['layer'], abs(p['feed']), 'm'),
              ('traverse_1dx', p['dx'], abs(p['traverse']), 'm')]
    rows = []
    for name, extent, rate, unit in events:
        if rate == 0:
            rows.append(dict(event=name, status='NOT_AVAILABLE_ZERO_PRESCRIBED_RATE'))
            continue
        duration = extent/rate
        steps = math.ceil(duration/dt_s)
        rows.append(dict(event=name, status='CONSTANT_CURRENT_COST_EXTRAPOLATION',
                         extent=extent, extent_unit=unit, physical_seconds=duration,
                         estimated_steps=steps, estimated_wall_seconds=steps*wall_per_step_s))
    return rows


def measure(runs, target_end_time_s, parent_receipts=None):
    chain = inspect_chain(runs, target_end_time_s, parent_receipts=parent_receipts)
    # Each exported VTK is hashed here too; no solver advance or rejection migration.
    records = []
    regular_dt = []
    total_wall = 0.
    total_steps = 0
    prior_identity = None
    for directory in runs:
        directory = Path(directory)
        series = json.loads((directory/'series.json').read_text())
        config = json.loads((directory/'resolved_config.json').read_text())
        debug = config['parameters']['debug']; process = config['parameters']['process']
        p = dict(dx=debug['dx']['value'], layer=process['layer_height']['value'],
                 bore_width=debug['bore_width']['value'], omega=process['tool_angular_speed']['value'],
                 traverse=process['tool_traverse_speed']['value'], feed=process['feedstock_axial_speed']['value'])
        performance = json.loads((directory/'performance.json').read_text())
        total_wall += performance['wall_seconds']
        frames = series['frames']
        total_steps += frames[-1]['step']-frames[0]['step']
        for frame in frames:
            path = directory/frame['file']
            if hashlib.sha256(path.read_bytes()).hexdigest() != frame['sha256']:
                raise ValueError('frame hash mismatch')
            if records and frame['step'] == records[-1]['step'] and frame['time_s'] == records[-1]['time_s']:
                continue  # repeated accepted seam; separate chain checker validates equality
            data = read_point_cloud(path)
            ids = np.sort(np.asarray(data['scalars']['particle_id'], dtype=np.int64))
            if prior_identity is not None and not np.array_equal(ids, prior_identity):
                raise ValueError('particle IDs changed; unsupported inlet observation')
            prior_identity = ids
            if records and (frame['step'] <= records[-1]['step'] or frame['time_s'] <= records[-1]['time_s']):
                raise ValueError('non-monotone accepted chain')
            row = observe(data, p, frame['time_s'])
            row.update(step=frame['step'], file=str(path))
            records.append(row)
        accepted = {f['step'] for f in frames}
        diagnostics = [json.loads(s) for s in (directory/'diagnostics.jsonl').read_text().splitlines() if s.strip()]
        regular_dt += [d['dt_s'] for d in diagnostics if d['step'] in accepted]
    pa, state = restore(Path(runs[-1])/'checkpoint.npz')
    # Accepted state only. Units of plastic_work are cumulative per-particle J.
    free = pa.actuator == 0
    bar = pa.component_id == 1
    first = next((r for r in records if r['free_plastic_above_floor_count']), None)
    dt = float(np.median(regular_dt)) if regular_dt else None
    cost = total_wall/total_steps if total_steps else None
    return dict(schema_version='1.1', classification='ACCEPTED_STATE_MECHANISM_OBSERVATION', production_authorized=False,
        chain_status=chain['status'], reached_requested_target=chain['reached_target'],
        paper_comparison=False,
        geometry_roi='original bar identity(component_id=1),free(actuator=0),particle centres:0<z<=layer,outside rotating square bore; NOT necessarily newly fed material or all deposited material',
        volume_operator='sum(m/rho),3Dmaterial volume proxy; NOT area ratio',
        ep_observation_floor=EP_OBSERVATION_FLOOR, first_free_plastic_observation=first,
        accepted_time_s=state['t'], accepted_step=state['count'], final_observation=records[-1],
        final_cumulative_plastic_work_J=float(pa.plastic_work.sum()),
        free_bar_cumulative_plastic_work_J=float(pa.plastic_work[free & bar].sum()),
        free_substrate_cumulative_plastic_work_J=float(pa.plastic_work[free & ~bar].sum()),
        negative_raw_plastic_work_count=int(np.count_nonzero(pa.plastic_work < 0)),
        final_ever_feed_plane_crossed_mass_kg=float(pa.m[pa.feed_crossed == 1].sum()),
        final_ledger=state['ledger'], frames=records,
        measured_cost_basis=dict(accepted_new_steps=total_steps, total_wall_seconds=total_wall,
                                median_dt_s=dt, wall_per_accepted_step_s=cost),
        cost_extrapolation=process_costs(p, dt, cost) if dt is not None and cost is not None else [],
        caution='Initial interference can cause startup plasticity; nonzero ep is NOT sustained feeding/deposition. Finite bar/prebonded interface; no new physics or long-run authorization.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='append', required=True, type=Path)
    parser.add_argument('--target-end-time', required=True, type=float)
    parser.add_argument('--parent-receipt', action='append', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = measure(args.run, args.target_end_time, args.parent_receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'frames'}, indent=2, allow_nan=False))
    raise SystemExit(0 if report['reached_requested_target'] else 2)

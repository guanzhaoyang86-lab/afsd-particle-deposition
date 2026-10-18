"""Plot complete saved diagnostic states only; no solver or GPU is imported.

The input session may be running. progress.json identifies its last complete
PASS; any newer, incomplete attempt is ignored. The changing JSON inputs are
copied to immutable per-plot snapshots. All positions are actual saved particle
centres. The first geometry is the session's complete initial parent, not a
fabricated time-zero state. Run in an existing NumPy/Matplotlib environment.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from afsd_deposition_progress_20261010 import frame_from_saved, summarize_progress


ROOT = Path(__file__).resolve().parents[1]


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inside(value):
    path = Path(value)
    path = (ROOT/path).resolve() if not path.is_absolute() else path.resolve()
    require(path.is_relative_to(ROOT/'outputs'), 'saved evidence must be inside project outputs')
    return path


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def complete_stage(folder):
    """Check this complete receipt, without rereading the historical closure."""
    manifest, supervision = folder/'manifest.json', folder/'supervision.json'
    m, s = read(manifest), read(supervision)
    require(m['status'] == 'PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION'
        and s['exit_code'] == 0 and s['source_inputs_still_match'] is True
        and s['manifest_sha256'] == sha(manifest) and s['log_sha256'] == sha(folder/'run.log'),
        'plot needs a complete, intact, successful stage receipt')
    require(m['accepted'] is False and m['restart_allowed'] is False and m['diagnostic_only'] is True
        and m['screened_continuation_authorized'] is True, 'diagnostic screened evidence required')
    screen = read(folder/'pair_screen.json')
    require(screen == m['pair_screen'] and screen['passed'] is True
        and screen['global_screen']['pass'] is True and screen['local_screen']['passed'] is True
        and screen['crossing_ids_match'] is True, 'complete passing coarse/two-half screen required')
    checkpoint = inside(m['screened_continuation_checkpoint'])
    require(checkpoint == (folder/'fine_2/NOT_RESTARTABLE_candidate.npz').resolve()
        and m['screened_continuation_checkpoint'] == m['screened_history'][-1]['checkpoint']
        and m['screened_continuation_sha256'] == m['screened_history'][-1]['sha256']
        and screen['fine_checkpoint_sha256'][-1] == m['screened_continuation_sha256']
        and sha(checkpoint) == m['screened_continuation_sha256'], 'last complete fine-two checkpoint binding differs')
    return m


def saved_geometry(meta):
    checkpoint = inside(meta['screened_continuation_checkpoint'])
    # These are trusted local solver checkpoints; their existing NumPy object
    # payload requires allow_pickle. No external arbitrary NPZ is accepted.
    with np.load(checkpoint, allow_pickle=True) as saved:
        state = saved['solver_data'].item()
        arrays = saved['particles'].item()['solid']['arrays']
        fields = {key: np.asarray(arrays[key]).copy() for key in
            ('particle_id', 'component_id', 'x', 'y', 'z', 'm', 'ep')}
    require(state['time_s'] == meta['screened_end_time_s'] and state['step'] == meta['screened_end_step'],
        'saved checkpoint endpoint differs from complete receipt')
    frame = frame_from_saved(SimpleNamespace(**fields), state, None, meta['observer'])
    observation = summarize_progress([frame], state['parameters'])
    require(observation['final']['saved_observer_matches'] is True, 'saved endpoint observer differs')
    return frame, state, observation['final']


def render(history, initial, latest, parameters, target):
    t = np.array([row['time_s'] for row in history])*1000
    gap = np.array([np.nan if row['nearest_eligible_uncrossed_gap_m'] is None
        else row['nearest_eligible_uncrossed_gap_m']*1e6 for row in history])
    crossed = np.array([row['supplied_crossed']['mass_kg'] for row in history])*1e6
    under = np.array([row['supplied_under_shoulder']['mass_kg'] for row in history])*1e6
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.7), layout='constrained')
    ax = axes[0, 0]
    ax.plot(t, gap, 'o-', color='#0072b2', ms=3, lw=1.5, label='Nearest eligible, not yet crossed')
    ax.axhline(0., color='#555555', ls='--', lw=1., label='Outlet plane')
    ax.set(xlabel='Simulated time (ms)', ylabel='Particle-centre gap (µm)', title='A. Approach to the bore outlet')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    ax = axes[0, 1]
    ax.step(t, crossed, where='post', color='#009e73', lw=1.8, label='Historically crossed through bore')
    ax.plot(t, under, 'o-', color='#d55e00', ms=3, lw=1.4, label='Currently under annular shoulder')
    ax.set(xlabel='Simulated time (ms)', ylabel='Eligible supplied-material mass (mg)',
        title='B. Passage and shoulder occupancy')
    if not np.any(crossed) and not np.any(under):
        ax.set_ylim(-.02, .1); ax.text(.5, .6, 'No supplied-material passage or shoulder occupancy yet',
            ha='center', transform=ax.transAxes, fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=.2)

    ids = np.asarray(latest['ids']); component = np.asarray(latest['component_ids'])
    eligible = np.isin(ids, latest['eligible_ids']); supplied = np.isin(ids, latest['crossing_ids'])
    x = np.asarray(latest['position'])*1000; x0 = np.asarray(initial['position'])*1000
    relative = (np.asarray(latest['position'])-latest['tool_origin'])@latest['tool_rotation']
    radius = float(parameters['tool_radius']); half = float(parameters['bore_width'])/2
    top = float(parameters.get('initial_substrate_top_z_m', 0.))
    under_mask = supplied & (latest['position'][:, 2] > top) & (relative[:, 2] < 0.)
    under_mask &= np.max(np.abs(relative[:, :2]), axis=1) >= half
    under_mask &= np.linalg.norm(relative[:, :2], axis=1) <= radius
    categories = [(component != 1, '#b5b5b5', 'Substrate', 9),
        ((component == 1) & ~eligible, '#8c6d31', 'Initially lower bar (ineligible)', 12),
        (eligible & ~supplied, '#0072b2', 'Eligible, not yet crossed', 14),
        (supplied & ~under_mask, '#009e73', 'Crossed supplied material', 20),
        (under_mask, '#d55e00', 'Supplied material under shoulder', 24)]
    for ax, coordinates, title, ylabel in ((axes[1, 0], (0, 2), 'C. Actual particle centres: XZ projection', 'World z (mm)'),
            (axes[1, 1], (0, 1), 'D. Actual particle centres: XY projection', 'World y (mm)')):
        a, b = coordinates
        ax.scatter(x0[:, a], x0[:, b], s=15, facecolors='none', edgecolors='#747474', lw=.45,
            alpha=.35, label=f'Initial session parent ({initial["time_s"]*1000:.3f} ms)')
        for mask, color, label, size in categories:
            if np.any(mask):
                ax.scatter(x[mask, a], x[mask, b], s=size, c=color, alpha=.65, edgecolors='none', label=label)
        ax.set(xlabel='World x (mm)', ylabel=ylabel, title=title)
        ax.set_aspect('equal', adjustable='box'); ax.grid(alpha=.15)
    origin = np.asarray(latest['tool_origin'])*1000
    axes[1, 0].plot([origin[0]-radius*1000, origin[0]+radius*1000], [origin[2]]*2,
        ls='--', color='#252525', lw=1., label='Outlet / shoulder plane')
    axes[1, 0].legend(loc='upper left', fontsize=6.8)
    square = np.array([[-half, -half, 0.], [half, -half, 0.], [half, half, 0.], [-half, half, 0.]])
    square = square @ np.asarray(latest['tool_rotation']).T + latest['tool_origin']
    axes[1, 1].add_patch(Polygon(square[:, :2]*1000, closed=True, fill=False, lw=1.1,
        edgecolor='#252525', label='Rotated open square bore'))
    axes[1, 1].add_patch(Circle(origin[:2], radius*1000, fill=False, lw=1., ls='--',
        edgecolor='#555555', label='Shoulder outer radius'))
    axes[1, 1].legend(loc='lower right', fontsize=6.8)
    fig.suptitle(f'Actual saved AFSD diagnostic states | latest complete pair: {latest["time_s"]*1000:.3f} ms\n'
        '450 °C hot start · finite bar · geometric observations, no metallurgical-bonding claim', fontsize=12)
    fig.supxlabel('Time plots connect saved samples only. Gap follows the nearest eligible centre; it is not one particle trajectory.\n'
        'Initial outlines and latest coloured points are direct projections, with no inferred material path.', fontsize=8)
    fig.savefig(target, dpi=180, bbox_inches='tight'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', required=True)
    parser.add_argument('--output', required=True, help='New output directory, never overwritten')
    args = parser.parse_args()
    session, output = inside(args.session), inside(args.output)
    require(session.is_dir() and not output.exists(), 'existing session and new plot output required')
    progress_bytes = (session/'progress.json').read_bytes()
    observation_bytes = (session/'deposition_progress.json').read_bytes()
    progress, observation = json.loads(progress_bytes), json.loads(observation_bytes)
    protocol = read(session/'protocol.json')
    require('latest_pass_stage' in progress, 'session has not published a complete parent yet')
    latest_folder, initial_folder = inside(progress['latest_pass_stage']), inside(protocol['initial_parent'])
    latest_meta, initial_meta = complete_stage(latest_folder), complete_stage(initial_folder)
    require(latest_meta['screened_continuation_sha256'] == progress['checkpoint_sha256']
        and latest_meta['screened_end_time_s'] == progress['physical_time_s'], 'live progress must name an exact complete PASS')
    latest, state, fresh = saved_geometry(latest_meta)
    initial, initial_state, _ = saved_geometry(initial_meta)
    require(state['parameters'] == initial_state['parameters'], 'initial/latest physical parameters differ')
    for key in ('ids', 'component_ids', 'mass'):
        left, right = np.argsort(initial['ids']), np.argsort(latest['ids'])
        require(np.array_equal(np.asarray(initial[key])[left], np.asarray(latest[key])[right]),
            'initial/latest immutable particle inventory differs')
    history = [row for row in observation['history'] if row['time_s'] <= latest['time_s']]
    require(history and history[-1]['time_s'] == latest['time_s']
        and all(a['time_s'] < b['time_s'] for a, b in zip(history, history[1:])),
        'observer publication is incomplete; retry after the next complete stage update')
    for key in ('supplied_crossed', 'supplied_under_shoulder', 'nearest_eligible_uncrossed_gap_m',
            'nearest_eligible_uncrossed_particle_id'):
        require(history[-1][key] == fresh[key], 'passive history endpoint disagrees with saved particles: '+key)
    output.mkdir(parents=True)
    (output/'progress_snapshot.json').write_bytes(progress_bytes)
    (output/'deposition_progress_snapshot.json').write_bytes(observation_bytes)
    target = output/'deposition_progress.png'
    render(history, initial, latest, state['parameters'], target)
    sources = {p.relative_to(ROOT).as_posix():sha(p) for p in
        (Path(__file__).resolve(), ROOT/'scripts/afsd_deposition_progress_20261010.py')}
    evidence = {}
    for folder, meta in ((initial_folder, initial_meta), (latest_folder, latest_meta)):
        for path in (folder/'manifest.json', folder/'supervision.json', folder/'pair_screen.json',
                folder/'run.log', inside(meta['screened_continuation_checkpoint'])):
            evidence[path.relative_to(ROOT).as_posix()] = sha(path)
        require(sha(inside(meta['screened_continuation_checkpoint'])) == meta['screened_continuation_sha256'],
            'saved checkpoint changed during plotting')
    receipt = dict(status='PASS_COMPLETE_SAVED_STATE_VISUALIZATION_NOT_PHYSICAL_VALIDATION',
        source_sha256=sources, evidence_sha256=evidence,
        session=session.relative_to(ROOT).as_posix(), session_status_at_capture=progress['status'],
        initial_parent=initial_folder.relative_to(ROOT).as_posix(), latest_complete_stage=latest_folder.relative_to(ROOT).as_posix(),
        initial_checkpoint=initial_meta['screened_continuation_checkpoint'],
        initial_checkpoint_sha256=initial_meta['screened_continuation_sha256'],
        latest_checkpoint=latest_meta['screened_continuation_checkpoint'],
        latest_checkpoint_sha256=latest_meta['screened_continuation_sha256'],
        initial_time_s=initial['time_s'], latest_complete_time_s=latest['time_s'], plotted_frame_count=len(history),
        first_plotted_time_s=history[0]['time_s'], final_gap_m=fresh['nearest_eligible_uncrossed_gap_m'],
        final_crossed=fresh['supplied_crossed'], final_under_shoulder=fresh['supplied_under_shoulder'],
        captured_progress_sha256=hashlib.sha256(progress_bytes).hexdigest(),
        captured_observation_sha256=hashlib.sha256(observation_bytes).hexdigest(),
        output_png=target.relative_to(ROOT).as_posix(), output_png_sha256=sha(target),
        physical_solver_evaluations=0, nnps_evaluations=0, material_trajectory_prescribed=False,
        accepted=False, restart_allowed=False, metallurgical_bonding_proven=False,
        initial_means='FIRST_COMPLETE_PARENT_OF_THIS_SESSION_NOT_TIME_ZERO',
        source_verification='CURRENT_RECEIPTS_AND_CHECKPOINTS_ONLY_NOT_NEW_FULL_HISTORY_PHYSICAL_QUALIFICATION',
        plotted_positions='ACTUAL_SAVED_CENTRES_PROJECTION_NO_INTERPOLATED_POSITIONS',
        runtime=dict(python=sys.executable, numpy=np.__version__, matplotlib=matplotlib.__version__), exit_code=0)
    with (output/'plot_receipt.json').open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(status=receipt['status'], plot=str(target), receipt=str(output/'plot_receipt.json'),
        latest_complete_time_s=latest['time_s'], final_gap_m=receipt['final_gap_m'])))


if __name__ == '__main__':
    main()

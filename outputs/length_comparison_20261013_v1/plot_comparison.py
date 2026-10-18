"""Plot saved equal-time length comparisons, without importing a physics solver.

Usage (a matplotlib environment is required):
  python plot_comparison.py [--output-dir outputs/length_comparison_20261013_v1/plots_NAME]
The default creates a new timestamped output folder. Existing outputs are never
overwritten. New comparison_*.json files are included on each invocation.
"""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
BASE = Path(__file__).resolve().parent
COLORS = {'long': '#2166AC', 'short': '#D97928'}
LABELS = {'long': '6.5 mm bar', 'short': '4.0 mm bar'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def source_comparisons():
    records = []
    for path in BASE.glob('comparison_*.json'):
        item = read(path)
        require(item.get('classification') in {
            'PASSIVE_QUALIFIED_EQUAL_TIME_FREE_BAR_LENGTH_COMPARISON',
            'PASSIVE_EQUAL_TIME_LENGTH_COMPARISON_WITH_ARCHIVED_CPU_REMOTE_REFERENCE'},
            'unrecognized comparison classification: ' + str(path))
        require(item.get('physical_solver_evaluations') == 0, 'comparison must be passive')
        time_s = float(item['time_s'])
        require(math.isfinite(time_s), 'finite physical time required')
        for case in ('long', 'short'):
            receipt = item['receipt_checks'][case]
            stage = (ROOT / receipt['stage']).resolve()
            require(stage.is_relative_to(ROOT / 'outputs'), 'stage outside project outputs')
            manifest_path = stage / 'manifest.json'
            require(sha(manifest_path) == receipt['manifest_sha256'], 'source manifest changed')
            manifest = read(manifest_path)
            supervision = read(stage / 'supervision.json')
            require(supervision['exit_code'] == 0 and supervision['source_inputs_still_match']
                and supervision['manifest_sha256'] == receipt['manifest_sha256']
                and supervision['log_sha256'] == sha(stage / 'run.log'), 'source supervision mismatch')
            checkpoint = ROOT / manifest['screened_continuation_checkpoint']
            require(sha(checkpoint) == receipt['final_checkpoint_sha256']
                == manifest['screened_continuation_sha256'], 'qualified endpoint changed')
            require(receipt['latest_pair_passed'], 'unscreened endpoint')
            for branch in ('fine', 'coarse'):
                endpoint = item['endpoint_metrics'][case][branch]
                require(abs(float(endpoint['time_s']) - time_s) <= 1e-13,
                        'comparison mixes physical times')
                require(endpoint['groups']['initial_eligible_near_outlet']['particle_count'] == 25
                    and abs(endpoint['groups']['initial_eligible_near_outlet']['initial_z_min_m']-.00175)<1e-12
                    and abs(endpoint['groups']['initial_eligible_near_outlet']['initial_z_max_m']-.00175)<1e-12,
                    'plot requires the same fixed initial outlet layer')
            params = item['endpoint_metrics'][case]['fine']['actuator']
            require(abs(params['actual_mean_downward_feed_m']-params['commanded_downward_feed_m']) < 1e-12,
                    'actuator command differs from observed imposed displacement')
        feeds = [item['endpoint_metrics'][case]['fine']['actuator']['commanded_downward_feed_m']
                 for case in ('long', 'short')]
        require(abs(feeds[0]-feeds[1]) < 1e-12, 'commands differ at matched time')
        records.append(dict(path=path, sha256=sha(path), data=item))
    records.sort(key=lambda record: record['data']['time_s'])
    require(records, 'no qualified comparison files found')
    require(len({round(record['data']['time_s'], 12) for record in records}) == len(records),
            'duplicate physical times; explicitly resolve duplicate evidence instead of selecting silently')
    return records


def layer(endpoint, height=.00275):
    selected = [entry for entry in endpoint['initial_z_layers']
                if abs(entry['initial_z_m']-height) < 1e-12]
    require(len(selected) == 1 and selected[0]['prescribed_z_particle_count'] == 0
            and selected[0]['metrics']['particle_count'] == 25,
            'common free layer must have 25 particles and no prescribed-z particles')
    return selected[0]['metrics']


def make_plot(records, output):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.titlesize': 12, 'axes.labelsize': 10, 'axes.linewidth': .8,
        'xtick.labelsize': 9, 'ytick.labelsize': 9, 'pdf.fonttype': 42,
        'savefig.facecolor': 'white'})
    times = np.array([record['data']['time_s']*1e3 for record in records])
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 8.5))
    figure.subplots_adjust(left=.085, right=.975, bottom=.17, top=.80, wspace=.28, hspace=.47)
    a, b, c, d = axes.ravel()
    feed = [record['data']['endpoint_metrics']['long']['fine']['actuator']['commanded_downward_feed_m']*1e6
            for record in records]
    a.plot(times, feed, '--', color='#777777', linewidth=1.6, label='Imposed feed (both)')
    points = []
    for case in ('long', 'short'):
        fine = [record['data']['endpoint_metrics'][case]['fine'] for record in records]
        coarse = [record['data']['endpoint_metrics'][case]['coarse'] for record in records]
        outlet = [-entry['groups']['initial_eligible_near_outlet']['mean_dz_m']*1e6 for entry in fine]
        upper = [entry['free_upper_segment']['shortening_m']/entry['free_upper_segment']['initial_length_m']*1e6
                 for entry in fine]
        radius = [layer(entry)['centroid_rms_radius_change_m']*1e6 for entry in fine]
        reaction = [entry['actuator']['reaction_on_actuator_z_N'] for entry in fine]
        coarse_reaction = [entry['actuator']['reaction_on_actuator_z_N'] for entry in coarse]
        for axis, values in ((a, outlet), (b, upper), (c, radius), (d, reaction)):
            axis.plot(times, values, '-o', color=COLORS[case], markersize=4.8,
                      linewidth=1.8, label=LABELS[case], zorder=3)
        d.plot(times, coarse_reaction, linestyle='none', marker='x', markersize=6.5,
               markeredgewidth=1.15, color=COLORS[case], alpha=.85, zorder=4)
        for i, entry in enumerate(fine):
            points.append(dict(case=case, time_s=entry['time_s'],
                imposed_feed_m=entry['actuator']['commanded_downward_feed_m'],
                outlet_mean_downward_displacement_m=-entry['groups']['initial_eligible_near_outlet']['mean_dz_m'],
                upper_free_segment_shortening_strain=upper[i]*1e-6,
                upper_free_segment=entry['free_upper_segment'],
                common_z0_2p75mm_layer_rms_radius_change_m=radius[i]*1e-6,
                actuator_reaction_z_N=reaction[i], latest_macro_coarse_reaction_z_N=coarse_reaction[i],
                crossed_count=entry['deposition']['supplied_crossed']['particle_count'],
                under_shoulder_count=entry['deposition']['supplied_under_shoulder']['particle_count']))
    a.set_title('(a) Feed and fixed outlet-layer motion', loc='left', pad=10)
    a.set_ylabel('Downward displacement (µm)')
    a.legend(frameon=False, fontsize=8.5, loc='upper left')
    a.text(.98, .045, r'Outlet cohort: $z_0=1.75$ mm, 25 particles', transform=a.transAxes,
           ha='right', va='bottom', fontsize=8, color='#555555')
    b.set_title('(b) Upper free-segment shortening', loc='left', pad=10)
    b.set_ylabel('Shortening / initial span (µstrain)')
    b.text(.98, .965, 'Different free spans (initial z):\n'
           '6.5 mm bar: 2.75–5.75 mm; 4.0 mm bar: 2.75–3.25 mm',
           transform=b.transAxes, va='top', ha='right', fontsize=7.8, color='#555555',
           bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=2.))
    c.set_title('(c) Common upper layer: lateral change', loc='left', pad=10)
    c.set_ylabel('Centroid RMS radius change (µm)')
    c.text(.98, .045, r'$z_0=2.75$ mm, 25 free particles in both cases',
           transform=c.transAxes, ha='right', va='bottom', fontsize=8, color='#555555')
    d.set_title('(d) Reaction on the axial actuator', loc='left', pad=10)
    d.set_ylabel('Upward reaction force (N)')
    for axis in axes.ravel():
        axis.set_xlabel('Physical time (ms)')
        axis.spines[['top', 'right']].set_visible(False)
        axis.grid(axis='y', color='#e3e7eb', linewidth=.65, zorder=0)
        axis.ticklabel_format(axis='y', style='plain', useOffset=False)
        axis.margins(x=.075, y=.16)
    handles = [Line2D([0], [0], color=COLORS[case], marker='o', linewidth=1.8,
                      markersize=5, label=LABELS[case]) for case in ('long', 'short')]
    handles.append(Line2D([0], [0], color='#555555', marker='x', linestyle='none',
                          markersize=6, label='Last-macro coarse endpoint (panel d)'))
    figure.legend(handles=handles, loc='upper center', bbox_to_anchor=(.53, .884),
                  ncol=3, frameon=False, fontsize=10)
    figure.suptitle('Hot-bar length sensitivity', x=.085, ha='left', y=.97, fontsize=20, fontweight='medium')
    figure.text(.085, .914, 'Independent hot-start controls  ·  Particle spacing 0.5 mm  ·  Early ramp (0–20 ms)',
                fontsize=11, color='#555555')
    last = records[-1]['data']
    counts = {case: last['endpoint_metrics'][case]['fine']['deposition'] for case in ('long', 'short')}
    count_text = ' / '.join(f"{LABELS[case]}: {counts[case]['supplied_crossed']['particle_count']} crossed, "
                  f"{counts[case]['supplied_under_shoulder']['particle_count']} under shoulder" for case in ('long','short'))
    figure.text(.085, .108, 'Matched saved times; lines connect observations. × shows a final-macro local difference, not a confidence interval.',
                fontsize=8.6, color='#555555')
    figure.text(.085, .077, count_text + '. No deposited layer established by this comparison.', fontsize=8.6, color='#555555')
    figure.text(.085, .045, 'Long-bar points may use the verified archived CPU reference with identical used commands; no GPU continuation state is substituted.',
                fontsize=8.6, color='#555555')
    figure.savefig(output/'length_comparison.png', dpi=180)
    figure.savefig(output/'length_comparison.pdf')
    plt.close(figure)
    summary = dict(classification='PASSIVE_EQUAL_TIME_LENGTH_COMPARISON_FIGURE',
        physical_solver_evaluations=0, sources=[dict(path=str(record['path'].relative_to(ROOT)).replace('\\','/'),
            sha256=record['sha256'], time_s=record['data']['time_s']) for record in records],
        plotting_script_sha256=sha(__file__), figure_files={name:sha(output/name)
            for name in ('length_comparison.png','length_comparison.pdf')}, plotted_points=points,
        latest_time_s=last['time_s'], latest_effects=last['equal_time_effect_vs_local_temporal_difference'],
        limitations=['Only equal-time qualified saved endpoints are plotted; curves between samples are visual connectors.',
            'Final-macro coarse/fine differences are local indicators, not confidence intervals or accumulated error bounds.',
            'The two upper free segments have different lengths and are not identical material cohorts.',
            'An early-ramp length response does not establish steady feeding, deposition, or the cause of the 120 ms GPU behavior.'])
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    output = (ROOT/args.output_dir).resolve() if args.output_dir else BASE/('plots_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S'))
    require(output.is_relative_to(BASE) and output != BASE and not output.exists(), 'a new subdirectory is required')
    records = source_comparisons()
    output.mkdir(parents=True, exist_ok=False)
    result = make_plot(records, output)
    print(json.dumps(dict(output=str(output), matched_time_points=len(records),
                         latest_time_ms=result['latest_time_s']*1e3), indent=2))


if __name__ == '__main__':
    main()

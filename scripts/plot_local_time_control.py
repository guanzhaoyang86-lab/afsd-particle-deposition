"""Plot saved local-time screens only; no solver, PySPH or particle import.

Run only after the requested study has stopped and its saved-only stage audits
exist. A saved readback PASS may retain a temporal FAIL; this renderer never
changes the physical or process qualification attached to that evidence.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PASS = 'PASS_SAVED_LOCAL_TIME_READBACK_NOT_ACCEPTED'
PRIMARY_PASS = 'PASS_LOCAL_PAIR_DIAGNOSTIC_NOT_ACCEPTED'
DIAGNOSTIC = dict(accepted=False, restart_allowed=False,
                  accepted_time_advanced_s=0., production_authorized=False)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def finite(value, name):
    value = float(value)
    require(math.isfinite(value) and value >= 0., 'invalid '+name)
    return value


def diagnostic(record, name):
    require(record.get('accepted') is False and record.get('restart_allowed') is False,
            name+' is not an unaccepted diagnostic')


class SavedJSON:
    def __init__(self):
        self.hashes = {}

    def read(self, path):
        path = Path(path).resolve()
        require(path.is_relative_to(ROOT/'outputs') and path.suffix == '.json',
                'only project saved JSON inputs are allowed')
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        name = rel(path)
        require(name not in self.hashes or self.hashes[name] == digest,
                'input changed during plotting: '+name)
        self.hashes[name] = digest
        return json.loads(data.decode('utf-8'))

    def unchanged(self):
        require(all(sha(ROOT/name) == digest for name, digest in self.hashes.items()),
                'saved inputs changed during rendering')


def screen_metrics(screen):
    diagnostic(screen, 'screen')
    local = screen['local_screen']; glob = screen['global_screen']
    rows = local['rows']; ids = local['particle_ids']
    require(len(rows) == len(ids) == len(set(ids)) == 918 and
            {row['particle_id'] for row in rows} == set(ids),
            'screen does not contain all 918 unique IDs')
    require(local['stress_relative_limit'] == .1 and local['plastic_relative_limit'] == .2,
            'original stress/plastic comparison limits differ')
    local_stress = max(finite(row['stress_relative'], 'local stress error') for row in rows)
    local_plastic = max(finite(row['plastic_relative'], 'local plastic error') for row in rows)
    global_stress = finite(glob['stress_relative'], 'global stress error')
    global_plastic = finite(glob['plastic_relative'], 'global plastic error')
    return dict(time_s=finite(screen['time_s'], 'screen time'), passed=screen['passed'],
        local_stress_relative=local_stress, local_plastic_relative=local_plastic,
        global_stress_relative=global_stress, global_plastic_relative=global_plastic,
        local_stress_normalized=local_stress/.1, local_plastic_normalized=local_plastic/.2,
        global_stress_normalized=global_stress/.1, global_plastic_normalized=global_plastic/.2,
        all_gate_maximum_normalized_error=finite(screen['normalized_error'], 'overall screen error'),
        crossing_ids_match=screen['crossing_ids_match'], local_particle_count=918,
        comparison_scope=screen['comparison_scope'])


def collect(base, reader):
    primary_folders = sorted(base.glob('attempt_[0-9][0-9][0-9]'))
    replay_folders = sorted(base.glob('replay_[0-9][0-9][0-9]'))
    require(primary_folders, 'no saved primary attempts')
    primary = []; replay = []; stages = []; elapsed = 0.; candidates = 0
    last_controller = None
    for phase, folders, prefix in (('primary', primary_folders, 'attempt'),
                                    ('replay', replay_folders, 'replay')):
        require([p.name for p in folders] == [f'{prefix}_{i:03d}' for i in range(1, len(folders)+1)],
                'nonconsecutive saved '+phase+' stages')
        for index, folder in enumerate(folders, 1):
            manifest_path = folder/'manifest.json'; audit_path = folder/'readback/summary.json'
            manifest = reader.read(manifest_path)
            supervision = reader.read(folder/'supervision.json')
            audit = reader.read(audit_path)
            require('RUNNING' not in manifest['status'] and supervision['exit_code'] in (0, 2),
                    'stage is incomplete or timed out; do not render an active study')
            require(audit['status'] == AUDIT_PASS and audit['exit_code'] == 0 and
                    audit['phase'] == phase and audit['index'] == index,
                    'saved-only audit absent or inconsistent')
            require(supervision['manifest_sha256'] == reader.hashes[rel(manifest_path)] and
                    audit['manifest_sha256'] == reader.hashes[rel(manifest_path)] and
                    audit['stage_status'] == manifest['status'], 'saved manifest binding differs')
            for name, record in (('manifest', manifest), ('supervision', supervision), ('readback', audit)):
                diagnostic(record, name)
            wall = finite(supervision['wall_seconds'], 'supervised wall time')
            require(wall == audit['supervised_wall_seconds'], 'supervised cost differs from readback')
            elapsed += wall; candidates += audit['candidate_count']
            require(abs(elapsed-audit['total_supervised_wall_seconds']) <= 1e-9 and
                    candidates == audit['total_candidate_count'], 'cumulative saved cost ledger differs')
            stage = dict(phase=phase, index=index, status=manifest['status'],
                readback_status=audit['status'], supervised_wall_seconds=wall,
                candidate_count=audit['candidate_count'], study_qualified=audit['study_qualified'],
                temporal_screens_all_pass=audit['temporal_screens_all_pass'])
            stages.append(stage)
            if phase == 'primary':
                last_controller = audit['controller_after']
                item = dict(attempt=index, stage_status=manifest['status'],
                    outcome='PASS' if manifest['status'] == PRIMARY_PASS else
                            'REJECT' if manifest['status'].startswith('TEMPORAL_REJECT') else 'PHYSICAL FAIL',
                    t0_s=manifest['proposal']['t0_s'], t1_s=manifest['proposal']['t1_s'],
                    coarse_dt_s=manifest['proposal']['coarse_dt_s'],
                    fine_dt_s=manifest['proposal']['fine_dt_s'], metrics=None)
                if manifest.get('screen') is not None:
                    path = folder/'screen.json'; screen = reader.read(path)
                    require(screen == manifest['screen'] and manifest['screens'] == [screen],
                            'primary screen file differs from manifest')
                    require(audit['evidence_sha256'].get(rel(path)) == reader.hashes[rel(path)],
                            'primary screen lacks audited hash binding')
                    item['metrics'] = screen_metrics(screen)
                primary.append(item)
            else:
                file_screens = []
                for path in sorted(folder.glob('screen_[0-9][0-9][0-9].json')):
                    screen = reader.read(path)
                    require(audit['evidence_sha256'].get(rel(path)) == reader.hashes[rel(path)],
                            'replay screen lacks audited hash binding')
                    file_screens.append(screen)
                    replay.append(dict(stage=index, screen_file=rel(path), **screen_metrics(screen)))
                require(file_screens == manifest['screens'], 'replay screen files differ from manifest')
    return dict(primary_attempts=primary, replay_screens=replay, stages=stages,
        qualified_primary_intervals=last_controller['qualified_intervals'],
        qualified_primary_endpoint_s=last_controller['time_s'],
        terminal_primary_controller_status=last_controller['status'],
        total_supervised_wall_seconds=elapsed, total_candidate_count=candidates,
        independent_replay_status='NOT_RUN' if not replay_folders else
            'PASS_REGISTERED_INTERVAL' if stages[-1]['study_qualified'] else 'NOT_QUALIFIED',
        study_qualified=stages[-1]['study_qualified'])


STYLES = (
    ('local_stress_normalized', '#0072B2', 'o', '-', 'All-ID maximum stress / 0.1'),
    ('local_plastic_normalized', '#D55E00', 's', '-', 'All-ID maximum plastic / 0.2'),
    ('global_stress_normalized', '#56B4E9', '^', '--', 'Global stress / 0.1'),
    ('global_plastic_normalized', '#CC79A7', 'D', '--', 'Global plastic / 0.2'),
)


def panel(ax, rows, coordinates):
    ax.set_facecolor('#FBFCFD'); ax.grid(axis='y', color='#DFE5E9', linewidth=.65)
    ax.axhline(1., color='#8F3830', linestyle=(0, (5, 3)), linewidth=1.2, zorder=1)
    for key, color, marker, linestyle, label in STYLES:
        ax.plot(coordinates, [row[key] for row in rows], color=color, marker=marker,
                linestyle=linestyle, linewidth=1.7, markersize=5.5, label=label,
                markeredgecolor='white', markeredgewidth=.55, zorder=3)
    maximum = max([1.]+[row[key] for row in rows for key, *_ in STYLES])
    ax.set_ylim(0., 1.18*maximum)
    ax.text(.015, 1., ' gate = 1', transform=ax.get_yaxis_transform(), color='#8F3830',
            ha='left', va='bottom', fontsize=9)
    ax.set_ylabel('Error / unchanged comparison limit')


def render(data, png):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
        'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 200})
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.6), gridspec_kw={'width_ratios': [1.15, 1.]})
    fig.subplots_adjust(left=.075, right=.975, bottom=.235, top=.705, wspace=.25)
    fig.text(.075, .952, 'Local timestep screens and independent replay',
             fontsize=20, weight='bold', color='#1C3040')
    fig.text(.075, .912, 'Independent hot-QS diagnostic  |  918 particles  |  accepted=false  |  restart_allowed=false',
             fontsize=11.5, color='#435663')
    fig.text(.075, .876,
        f"Qualified primary endpoint: {1000*data['qualified_primary_endpoint_s']:.6f} ms   |   "
        f"{data['total_candidate_count']} saved candidates   |   "
        f"{data['total_supervised_wall_seconds']:.3f} s supervised wall time   |   "
        f"Replay: {data['independent_replay_status']}", fontsize=10.4, color='#435663')
    handles = [Line2D([0], [0], color=color, marker=marker, linestyle=style,
                      linewidth=1.7, label=label) for _, color, marker, style, label in STYLES]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.068, .838),
               ncol=2, frameon=False, fontsize=10.5, columnspacing=3., handlelength=3.2)
    primary = data['primary_attempts']; actual = [item for item in primary if item['metrics'] is not None]
    panel(axes[0], [item['metrics'] for item in actual], [item['attempt'] for item in actual])
    axes[0].set_title('A  Same-immediate-parent comparisons', loc='left', pad=20, fontsize=12)
    axes[0].set_xlabel('Primary attempt number; second row = shared endpoint [ms]')
    axes[0].set_xticks([item['attempt'] for item in primary],
        [f"{item['attempt']}\n{1000*item['t1_s']:.4f}" for item in primary], fontsize=8.5)
    axes[0].set_xlim(.5, len(primary)+.5)
    for item in primary:
        axes[0].text(item['attempt'], 1.025, item['outcome'],
            transform=axes[0].get_xaxis_transform(), ha='center', va='bottom',
            fontsize=8.4, weight='bold', color='#24725D' if item['outcome'] == 'PASS' else '#A13C2D')
    replay = data['replay_screens']
    if replay:
        panel(axes[1], replay, [1000*row['time_s'] for row in replay])
        axes[1].set_xlabel('Independent replay shared endpoint time [ms]')
        if all(row['passed'] for row in replay):
            axes[1].text(.5, 1.025, f'All {len(replay)} shared endpoints PASS',
                transform=axes[1].transAxes, ha='center', va='bottom', fontsize=9,
                weight='bold', color='#24725D')
        else:
            for row in replay:
                axes[1].text(1000*row['time_s'], 1.025, 'PASS' if row['passed'] else 'REJECT',
                    transform=axes[1].get_xaxis_transform(), ha='center', va='bottom', fontsize=8,
                    color='#24725D' if row['passed'] else '#A13C2D')
    else:
        panel(axes[1], [], [])
        axes[1].set_xlim(3., 3.5); axes[1].set_xlabel('Independent replay shared endpoint time [ms]')
        message = 'NOT_RUN\nNo independent replay results exist.' if data['independent_replay_status'] == 'NOT_RUN' else \
                  'NO COMPLETED SHARED SCREEN\nReplay remains NOT_QUALIFIED.'
        axes[1].text(.5, .5, message, transform=axes[1].transAxes, ha='center', va='center',
                     fontsize=13, weight='bold', color='#596C78', linespacing=1.8)
    axes[1].set_title('B  Half-primary-fine-grid replay from 3 ms', loc='left', pad=20, fontsize=12)
    fig.text(.075, .145, 'Attempt endpoints and parent states differ. The left panel is a controller history, not a convergence-order plot.',
             fontsize=10.3, color='#435663')
    fig.text(.075, .113, 'Lines connect saved error metrics only. PASS also requires position, temperature, crossing IDs and all physical gates.',
             fontsize=10.3, color='#435663')
    fig.text(.075, .079, 'All states remain unaccepted diagnostics. Local qualification is not deposition or accepted process advancement.',
             fontsize=10.5, weight='bold', color='#A13C2D')
    fig.text(.075, .046, 'The earlier common-initial 500/250 us failure remains unchanged. No synthetic trajectory or unrun replay is displayed.',
             fontsize=10., color='#435663')
    fig.savefig(png, facecolor='white'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--terminal-record', help='Optional saved terminal JSON; bound and reported without inferred semantics')
    args = parser.parse_args(); started = time.perf_counter()
    base = (ROOT/args.base).resolve(); output = (ROOT/args.output).resolve()
    require(base.is_relative_to(ROOT/'outputs') and base.is_dir(), 'existing project output study required')
    require(output.is_relative_to(ROOT/'outputs') and output != ROOT/'outputs' and not output.exists(),
            'fresh project output directory required')
    reader = SavedJSON(); data = collect(base, reader)
    terminal = reader.read((ROOT/args.terminal_record).resolve()) if args.terminal_record else None
    source_hash = sha(Path(__file__))
    reader.unchanged()
    output.mkdir(parents=True, exist_ok=False)
    png = output/'local_time_control_screens.png'; render(data, png)
    reader.unchanged()
    require(sha(Path(__file__)) == source_hash, 'renderer source changed')
    summary = dict(status='RENDERED_SAVED_LOCAL_TIME_DIAGNOSTIC_NOT_PROCESS_ACCEPTANCE',
        input_study=rel(base), command=[sys.executable, *sys.argv],
        source_sha256={rel(Path(__file__)): source_hash}, input_json_sha256=reader.hashes,
        artifacts_sha256={rel(png): sha(png)}, matplotlib_version=matplotlib.__version__,
        terminal_record=terminal, physical_evaluations=0, material_evaluations=0,
        nnps_evaluations=0, plotted_fields_interpolated=False, convergence_order_claimed=False,
        wall_seconds=time.perf_counter()-started, **data, **DIAGNOSTIC)
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=summary['status'], png=rel(png),
        qualified_primary_endpoint_s=data['qualified_primary_endpoint_s'],
        total_supervised_wall_seconds=data['total_supervised_wall_seconds'],
        independent_replay_status=data['independent_replay_status'])))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

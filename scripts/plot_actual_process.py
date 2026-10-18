"""Plot saved, independently audited checkpoint samples; no simulation."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(directory):
    directory = directory.resolve()
    source = directory/'accepted_trajectory.npz'
    summary_path = directory/'summary.json'
    summary = json.loads(summary_path.read_text(encoding='utf-8'))
    if digest(source) != summary['trajectory_sha256']:
        raise ValueError('trajectory hash mismatch')
    target = directory/'process_progress.png'
    if target.exists():
        raise ValueError('refuse plot overwrite')
    with np.load(source, allow_pickle=False) as data:
        x, t, ids = data['position_m'], data['time_s'], data['particle_id']
        ep, temp = data['plastic_strain'], data['temperature_K']
        origin, rotation = data['tool_origin_m'], data['tool_rotation']
        crossing = data['eligible_crossing_count']
        eligible = np.isin(ids, summary['eligible_ids'])
        if not bool(data['actual_checkpoint_samples_only']) or bool(data['interpolation_performed']):
            raise ValueError('only actual samples may be plotted')
    fig = plt.figure(figsize=(12.5, 8.3), layout='constrained')
    a = fig.add_subplot(2, 2, 1, projection='3d')
    cloud = a.scatter(*(x[-1]*1000).T, c=ep[-1], s=10, cmap='viridis', depthshade=False)
    a.set(xlabel='x (mm)', ylabel='y (mm)', zlabel='z (mm)',
          title=f'Actual final particles at {t[-1]*1000:.3f} ms')
    a.view_init(elev=22, azim=-53)
    fig.colorbar(cloud, ax=a, shrink=.7, label='Equivalent plastic strain')
    # This is only the shoulder reference height; it is not a filled tool face.
    a.text2D(.02, .98, '342 particles; shoulder contact is not centre passage',
             transform=a.transAxes, fontsize=8, va='top')
    b = fig.add_subplot(2, 2, 2)
    gap = np.min(x[:, eligible, 2]-origin[:, None, 2], axis=1)*1e6
    b.plot(t*1000, gap, color='#176b9b', marker='.', markersize=2)
    b.axhline(0., color='black', linewidth=.8)
    b.set(xlabel='Physical time (ms)', ylabel='Centre above shoulder (um)',
          title='Nearest initially eligible bar centre')
    b.text(.03, .08, f'Eligible square-bore crossings: {int(crossing[-1])}',
           transform=b.transAxes, fontsize=11, weight='bold')
    c = fig.add_subplot(2, 2, 3)
    records = summary['interval_energy']['records']
    raw = np.array([row['global_raw_residual_J'] for row in records])
    work = np.array([row['absolute_work_scale_J'] for row in records])
    c.semilogy(t[1:]*1000, np.cumsum(np.abs(raw))*1e6, label='Sum of absolute raw defects')
    c.semilogy(t[1:]*1000, (1e-7+.05*np.cumsum(work))*1e6, label='Registered interval limit')
    c.set(xlabel='Physical time (ms)', ylabel='Energy (uJ)', title='Uncorrected mechanical energy ledger')
    c.legend(fontsize=8)
    d = fig.add_subplot(2, 2, 4)
    d.plot(t*1000, temp.min(axis=1)-273.15, label='Minimum')
    d.plot(t*1000, temp.max(axis=1)-273.15, label='Maximum')
    d.set(xlabel='Physical time (ms)', ylabel='Temperature (deg C)', title='Solved temperature range')
    d.legend(fontsize=8)
    for ax in (b, c, d):
        ax.grid(alpha=.2)
    fig.suptitle('ENGINEERING_APPROXIMATE: independent hot quasi-static process\n'
                 'Actual checkpoint samples; no experimental or spatial validation', fontsize=12)
    fig.savefig(target, dpi=160)
    plt.close(fig)
    record = dict(status='PLOT_OF_SAVED_FIELDS_ONLY', physical_evaluations=0,
        input_sha256={str(source):digest(source), str(summary_path):digest(summary_path)},
        script_sha256=digest(Path(__file__)), output_sha256=digest(target),
        accepted_time_s=float(t[-1]), samples=len(t), crossing_count=int(crossing[-1]),
        lines_join_actual_samples=True, additional_trajectory_samples_created=False)
    (directory/'process_progress_receipt.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    print(json.dumps(dict(image=str(target), **{key:record[key] for key in ('samples','accepted_time_s','crossing_count')})))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('directory', type=Path)
    main(parser.parse_args().directory)

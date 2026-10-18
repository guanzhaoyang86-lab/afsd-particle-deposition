"""Plot actual accepted samples of the first strictly eligible crossing ID."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(directory):
    directory = Path(directory).resolve()
    source, report = directory/'accepted_trajectory.npz', directory/'summary.json'
    summary = json.loads(report.read_text(encoding='utf-8'))
    if sha(source) != summary['trajectory_sha256']:
        raise ValueError('Saved trajectory hash differs')
    event = summary['accepted_centre_passage_evidence']['first_event']
    if event is None:
        raise ValueError('No qualified first passage')
    observed = event['accepted_event']
    pid = observed['particle_id']
    target = directory/'first_passage_actual_samples.png'
    if target.exists():
        raise ValueError('Refuse plot overwrite')
    with np.load(source, allow_pickle=False) as data:
        if not bool(data['actual_checkpoint_samples_only']) or bool(data['interpolation_performed']):
            raise ValueError('Only original accepted samples may be plotted')
        index = int(np.flatnonzero(data['particle_id'] == pid)[0])
        t = data['time_s']*1000
        points = data['position_m'][:, index]
        local = np.einsum('ni,nij->nj', points-data['tool_origin_m'], data['tool_rotation'])
        ep = data['plastic_strain'][:, index]
        temp = data['temperature_K'][:, index]-273.15
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8), layout='constrained')
    a, b, c, d = axes.ravel()
    a.plot(t, local[:, 2]*1e6, '.-', markersize=2)
    a.axhline(0., color='black', linewidth=.8)
    a.set(xlabel='Physical time (ms)', ylabel='Centre above shoulder (um)',
          title=f'Initially eligible bar particle {pid}')
    zoom = t >= 36.
    b.plot(t[zoom], local[zoom, 2]*1e6, '.-', markersize=4, label='Actual accepted samples')
    b.axhline(0., color='black', linewidth=.8)
    b.scatter([observed['time_s']*1000], [0.], marker='x', s=60, color='#b43c23',
              label='Passive segment intersection')
    lo, hi = event['actual_bracket_time_s']
    b.axvspan(lo*1000, hi*1000, color='#8cb9d7', alpha=.2)
    b.set(xlabel='Physical time (ms)', ylabel='Centre above shoulder (um)',
          title='Saved endpoint bracket around first passage')
    b.legend(fontsize=8)
    c.plot(t, ep, '.-', markersize=2)
    c.set(xlabel='Physical time (ms)', ylabel='Equivalent plastic strain',
          title='Actual material response of the same particle')
    d.plot(t, temp, '.-', markersize=2)
    d.set(xlabel='Physical time (ms)', ylabel='Temperature (deg C)',
          title='Solved temperature of the same particle')
    for ax in axes.ravel():
        ax.grid(alpha=.2)
    fig.suptitle('ENGINEERING_APPROXIMATE: first accepted square-bore centre passage\n'
                 'Independent hot quasi-static model; centre passage is not deposition validation', fontsize=12)
    fig.savefig(target, dpi=150)
    plt.close(fig)
    receipt = dict(status='PLOT_OF_ACTUAL_ACCEPTED_FIELDS_ONLY', physical_evaluations=0,
        source_sha256=sha(Path(__file__)), input_sha256={str(source):sha(source), str(report):sha(report)},
        output_sha256=sha(target), particle_id=pid, samples=len(t),
        event_time_s=observed['time_s'], actual_endpoint_bracket_s=[lo, hi],
        mass_kg=observed['mass_kg'], additional_trajectory_samples_created=False,
        lines_join_actual_samples=True, crossing_marker_is_observer_intersection_not_additional_sample=True)
    (directory/'first_passage_plot_receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps(dict(image=str(target), particle_id=pid, samples=len(t))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    main(parser.parse_args().directory)

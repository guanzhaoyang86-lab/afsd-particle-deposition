"""Plot hash-bound passive observations, with no physical map evaluation."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle


ROOT = Path(__file__).resolve().parents[1]


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def require(value, message):
    if not value: raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--observation', required=True)
    options = parser.parse_args(); directory = Path(options.observation).resolve()
    require(directory.is_relative_to(ROOT/'outputs'), 'observation outside outputs')
    target = directory/'deposition_progress.png'; receipt = directory/'plot_receipt.json'
    require(not target.exists() and not receipt.exists(), 'refuse to overwrite plot evidence')
    summary_path = directory/'summary.json'; input_path = directory/'observation_inputs.npz'
    summary = json.loads(summary_path.read_text(encoding='utf-8'))
    require(summary['status']=='PASS_PASSIVE_ACCEPTED_ARRAY_OBSERVATION' and summary['exit_code']==0, 'observation failed')
    require(sha(input_path)==summary['observation_inputs_sha256'], 'observation inputs changed')
    checkpoint = ROOT/summary['accepted_checkpoint']
    require(sha(checkpoint)==summary['accepted_checkpoint_sha256'], 'accepted checkpoint changed')
    sources = {p.relative_to(ROOT).as_posix():sha(p) for p in (Path(__file__).resolve(),summary_path,input_path,checkpoint)}
    with np.load(input_path,allow_pickle=False) as saved: data = {k:saved[k] for k in saved.files}
    with np.load(checkpoint,allow_pickle=True) as saved:
        state = saved['solver_data'].item(); arrays = saved['particles'].item()['solid']['arrays']
    p = state['parameters']; ids = data['particle_id']; times = data['time_s']
    local = np.einsum('fni,fij->fnj',data['position_m']-data['tool_origin_m'][:,None,:],data['tool_rotation'])
    require(np.array_equal(np.column_stack([arrays[k] for k in ('x','y','z')]),data['position_m'][-1]), 'final particle fields differ')
    require(np.array_equal(arrays['particle_id'],ids), 'final ID ordering differs')
    eligible = np.isin(ids,state['bore_observation']['eligible_ids']); bar = data['component_id']==1
    initial_lower = bar & ~eligible; remaining = eligible & ~np.isin(ids,state['bore_observation']['crossing_ids'])
    colors = ['#d45500','#009e73','#b34db2','#0072b2','#d6a300']
    fig, axes = plt.subplots(2,2,figsize=(13,9),layout='constrained')
    ax, height_ax, ep_ax, thermal_ax = axes.flat
    ax.scatter(local[-1,initial_lower,0]*1e3,local[-1,initial_lower,1]*1e3,s=15,c='#a4a4a4',alpha=.7,label='Initially lower bar (ineligible)')
    ax.scatter(local[-1,remaining,0]*1e3,local[-1,remaining,1]*1e3,s=13,c='#c2d9ec',label='Eligible, not yet crossed')
    half = p['bore_width']*500; radius = p['tool_radius']*1000
    ax.add_patch(Rectangle((-half,-half),2*half,2*half,fill=False,color='#303030',lw=1.5,label='Open square bore'))
    ax.add_patch(Circle((0,0),radius,fill=False,color='#666666',ls='--',lw=1.1,label='Shoulder outer edge'))
    duration = []
    for index,row in enumerate(summary['crossed_particle_rows']):
        k = int(np.flatnonzero(ids==row['particle_id'])[0]); color = colors[index%len(colors)]
        ax.plot(local[:,k,0]*1e3,local[:,k,1]*1e3,color=color,lw=1.2)
        ax.scatter(local[-1,k,0]*1e3,local[-1,k,1]*1e3,color=color,s=40,zorder=5)
        ax.annotate(str(row['particle_id']),(local[-1,k,0]*1e3,local[-1,k,1]*1e3),xytext=(5,6),textcoords='offset points',color=color,fontsize=9)
        observed = times>=.030
        height_ax.plot(times[observed]*1e3,local[observed,k,2]*1e6,color=color,lw=1.4,label=f'PID {row["particle_id"]}')
        et = row['crossing_time_s']; post = np.flatnonzero(times>=et); first = int(post[0])
        ep_ax.plot(times[post]*1e3,data['plastic_strain'][post,k]-data['plastic_strain'][first,k],color=color,lw=1.5,label=f'PID {row["particle_id"]}')
        height_ax.scatter(times[first]*1e3,local[first,k,2]*1e6,color=color,s=24,zorder=5)
        delta = np.diff(data['plastic_strain'][post,k]); dt = np.diff(times[post])
        tolerance = 128*np.finfo(float).eps*np.maximum(1.,np.maximum(np.abs(data['plastic_strain'][post[:-1],k]),np.abs(data['plastic_strain'][post[1:],k])))
        positive = delta>tolerance; run = longest = 0.
        for present,increment in zip(positive,dt):
            run = run+float(increment) if present else 0.; longest = max(longest,run)
        duration.append(dict(particle_id=row['particle_id'],reported_crossing_time_s=et,
            actual_accepted_bracket_s=[float(times[first-1]),float(times[first])],
            first_post_crossing_actual_step=int(first),post_crossing_baseline_time_s=float(times[first]),
            post_crossing_observed_time_span_s=float(times[-1]-times[first]),
            longest_consecutive_positive_plastic_increment_sampled_interval_s=longest,
            positive_plastic_interval_count=int(positive.sum()),negative_plastic_interval_count=int(np.count_nonzero(delta < -tolerance)),
            total_post_crossing_sampled_intervals=len(delta),post_crossing_plastic_increment=row['post_crossing_plastic_increment'],
            current_square_projection_exit_distance_m=row['distance_to_square_projection_exit_m']))
    ax.set(xlabel='Tool-frame x (mm)',ylabel='Tool-frame y (mm)',title='Actual trajectories; three crossed centres remain inside bore',xlim=(-3,3),ylim=(-3,3))
    ax.set_aspect('equal');ax.legend(loc='lower left',fontsize=7)
    height_ax.axhline(0,color='#444444',lw=1,ls='--');height_ax.set(xlabel='Accepted time (ms)',ylabel='Centre height above outlet (um)',title='Outlet passage from accepted positions')
    height_ax.legend(fontsize=8);height_ax.grid(alpha=.18)
    ep_ax.set(xlabel='Accepted time (ms)',ylabel='Plastic strain increment',title='Plasticity after first accepted post-crossing sample')
    ep_ax.legend(fontsize=8);ep_ax.grid(alpha=.18)
    scatter = thermal_ax.scatter(local[-1,:,0]*1e3,local[-1,:,2]*1e3,c=arrays['temperature']-273.15,cmap='inferno',s=22,edgecolors='none')
    height = (p['tool_top']-p['layer'])*1000
    for left,width in ((-radius,radius-half),(half,radius-half)):
        thermal_ax.add_patch(Rectangle((left,0),width,height,facecolor='#dedede',edgecolor='#777777',alpha=.5,zorder=0))
    thermal_ax.axhline(0,color='#666666',lw=.8,ls='--')
    thermal_ax.set(xlabel='Tool-frame x (mm)',ylabel='Height above outlet (mm)',title='Final actual temperature; tool central section shaded',xlim=(-3,3))
    fig.colorbar(scatter,ax=thermal_ax,label='Temperature (C)',shrink=.82)
    fig.suptitle(f'Independent 450 C hot quasi-static branch | {times[-1]*1000:.3f} ms | NOT FULL DEPOSITION',fontsize=15)
    fig.text(.5,-.012,f'{len(times)} actual accepted frames. Lines connect stored samples; no new states or physical solutions. Newly crossed material outside the bore: 0 mg.',ha='center',fontsize=9)
    fig.savefig(target,dpi=170,bbox_inches='tight');plt.close(fig)
    for name,expected in sources.items(): require(sha(ROOT/name)==expected,'source changed while plotting')
    result = dict(status='PASS_SAVED_ARRAY_VISUALIZATION_NOT_PHYSICAL_VALIDATION',command=[sys.executable,'-B',str(Path(__file__).resolve()),*sys.argv[1:]],
        source_and_input_sha256=sources,output_png_sha256=sha(target),physical_solver_evaluations=0,nnps_evaluations=0,
        actual_frame_count=len(times),accepted_time_s=float(times[-1]),crossed_particle_sampled_plastic_duration=duration,
        duration_interpretation='Positive increments across consecutive actual sampled intervals; no claim about a continuous substep path or continuous plastic rate.',
        lines_connect_actual_samples_only=True,exit_code=0)
    receipt.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],plot=str(target),duration=duration),indent=2))


if __name__=='__main__':main()

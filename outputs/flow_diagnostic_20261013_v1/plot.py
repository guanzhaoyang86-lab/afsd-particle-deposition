"""Render measured data only; no solver or interpolated particle states."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
result = json.loads((OUT/'diagnostic.json').read_text(encoding='utf-8'))
records = result['records']
tracked = [v['particle_id'] for v in records[0]['tracked_crossed']]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig, ax = plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
last = records[-1]['layers']
z = [v['initial_z_mm'] for v in last]
ax[0,0].plot([v['mean_dz_um'] for v in last],z,'o-',color='#245b86')
ax[0,0].axvline(0,color='0.7',lw=.7)
ax[0,0].set(xlabel='Mean axial displacement over interval (um)',ylabel='Initial layer height z0 (mm)',
    title='A  Axial motion of persistent layers')
ax[0,1].plot([v['centered_rms_radius_change_um'] for v in last],z,'o-',color='#a95f24')
ax[0,1].axvline(0,color='0.7',lw=.7)
ax[0,1].set(xlabel='Change in centroid-based RMS radius (um)',ylabel='Initial layer height z0 (mm)',
    title='B  Radial size change of the same layers')
times = [v['time_s']*1e3 for v in records]
for k in [3,4,5,8,12]:
    ax[1,0].plot(times,[v['layers'][k]['mean_dz_um'] for v in records],label=f'z0={z[k]:.2f} mm')
ax[1,0].set(xlabel='Physical time (ms)',ylabel='Mean axial displacement (um)',
    title='C  Axial displacement at saved endpoints')
ax[1,0].legend(fontsize=8,ncol=2)
for j,pid in enumerate(tracked):
    ax[1,1].plot(times,[v['tracked_crossed'][j]['radius_mm'] for v in records],label=str(pid))
ax[1,1].set(xlabel='Physical time (ms)',ylabel='Radius of particle centre (mm)',
    title='D  The four previously crossed particles')
ax[1,1].legend(fontsize=8,ncol=2,title='Persistent ID',title_fontsize=8)
fig.suptitle('AFSD saved-state diagnostic | 89.590-120.439 ms | no new simulation',fontsize=13)
fig.savefig(OUT/'flow_diagnostic.png',dpi=180)
fig.savefig(OUT/'flow_diagnostic.pdf')
print('Wrote flow_diagnostic.png and flow_diagnostic.pdf')

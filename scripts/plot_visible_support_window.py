"""Plot audited saved arrays only; no interpolation or physics evaluation."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/visible_support_20261004/readback_v1'
BASE=ROOT/'outputs/visible_support_20261004/process_v1/fine'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))


def main():
    target=OUT/'accepted_window.png';receipt=OUT/'plot_receipt.json'
    if target.exists() or receipt.exists():raise ValueError('refuse plot evidence overwrite')
    summary=read(OUT/'summary.json');trajectory=ROOT/summary['trajectory']['path']
    if summary['status']!='PASS_SAVED_EVIDENCE_READBACK_NOT_FRESH_PHYSICS' or sha(trajectory)!=summary['trajectory']['sha256']:
        raise ValueError('audited actual trajectory absent')
    rows=read(BASE/'candidate_history.json')
    screens=[read(BASE/f'pair_{n:03d}_temporal_screen.json') for n in range(1,9)]
    inputs=[OUT/'summary.json',trajectory,BASE/'candidate_history.json']+[BASE/f'pair_{n:03d}_temporal_screen.json' for n in range(1,9)]
    bound={p.relative_to(ROOT).as_posix():sha(p) for p in inputs}
    with np.load(trajectory,allow_pickle=False) as data:
        t=1000*data['time_s'];ids=data['particle_id'];x=data['position_m'];ep=data['plastic_strain'];origin=data['tool_origin_m']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    a,b,c,d=axes.ravel()
    a.step(t,np.r_[0,[r['blocked_edge_count'] for r in rows]],where='post',color='#ba4d2b',lw=2)
    a.scatter(t[1:],[r['blocked_edge_count'] for r in rows],s=15,color='#ba4d2b')
    a.axvline(39.875,color='0.45',ls=':',label='Old model failed here')
    a.set(title='Tool-occluded raw edges excluded from ALL operators',ylabel='Excluded edge count',yticks=[0,1,2]);a.legend(fontsize=8)
    for pid,color in [(238,'#247ba0'),(262,'#8164a3')]:
        k=int(np.flatnonzero(ids==pid)[0]);label=f'PID {pid}'
        b.plot(t,1e6*(x[:,k,2]-origin[:,2]),'o-',color=color,label=label,ms=3)
        c.plot(t,ep[:,k]-ep[0,k],'o-',color=color,label=label,ms=3)
    b.axhline(0,color='0.6',ls=':');b.set(title='Accepted after-passage motion',ylabel='Centre height from outlet plane (um)');b.legend()
    c.set(title='Plastic strain continues after passage',ylabel='Plastic strain increase over this window');c.legend()
    errors=[100*max(row['stress_relative'] for row in screen['local_screen']['rows']) for screen in screens]
    d.plot([1000*s['time_s'] for s in screens],errors,'o-',color='#338361',ms=4,label='Recomputed local maximum')
    d.axhline(10,color='0.5',ls='--',label='Registered numerical screen: 10%')
    d.set(title='Same-time coarse / fine stress comparison',ylabel='Local stress difference (%)',ylim=(0,11));d.legend(fontsize=8)
    for ax in axes.ravel():ax.set_xlabel('Independent hot quasi-static time (ms)');ax.set_xlim(t[0],t[-1]);ax.grid(alpha=.2)
    fig.suptitle('Visible-support branch: 16 accepted steps to 40.375 ms\nTwo crossed particles remain inside the bore projection; deposition is not established',fontsize=13)
    fig.savefig(target,dpi=170);plt.close(fig)
    for name,digest in bound.items():
        if sha(ROOT/name)!=digest:raise ValueError('plot input changed')
    receipt.write_text(json.dumps(dict(status='PLOT_OF_AUDITED_SAVED_ARRAYS_ONLY',physical_evaluations=0,
        input_sha256=bound,script_sha256=sha(Path(__file__)),image_sha256=sha(target)),indent=2),encoding='utf-8')
    print(target)


if __name__=='__main__':main()

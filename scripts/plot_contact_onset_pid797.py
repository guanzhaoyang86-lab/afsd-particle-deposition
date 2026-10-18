"""Matplotlib rendering of already-exported scientific evidence; no PySPH import."""
from pathlib import Path
import hashlib
import json
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/contact_onset_time_study_20261005/v1/visualization_new'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source=OUT/'data.json';data=json.loads(source.read_text(encoding='utf-8'));rows=data['rows']
    assert sha(OUT/'pid797_saved_states.csv')==data['csv_sha256']
    parent=next(row for row in rows if row['branch']=='shared_parent')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.titlesize':11.5,
        'axes.labelsize':11.5,'xtick.labelsize':10,'ytick.labelsize':10,'savefig.dpi':200,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(3,1,figsize=(10.8,9.2),sharex=True)
    fig.subplots_adjust(left=.115,right=.965,top=.825,bottom=.19,hspace=.27)
    fig.text(.115,.963,'Contact-onset time-step study · PID 797',fontsize=18,weight='bold',color='#152838')
    fig.text(.115,.925,'Independent hot quasi-static model  |  saved diagnostic states only',fontsize=11.5,color='#435663')
    fig.text(.115,.893,'Accepted process remains at 62.5 μs / 14 steps; diagnostic endpoints reach 70.3125 μs.',fontsize=10.5,color='#435663')
    styles=[('h','#D55E00','s',r'$h=3.90625$ μs · 2 steps'),
        ('h_over_2','#0072B2','o',r'$h/2=1.953125$ μs · 4 steps'),
        ('h_over_4','#009E73','^',r'$h/4=0.9765625$ μs · 8 steps')]
    panels=[('q_Pa',1e-6,r'Equivalent stress $q$ [MPa]','a'),
        ('ep',1e5,r'Plastic strain $\bar\epsilon^p$ [$10^{-5}$]','b'),
        ('strain_rate_per_s',1.,r'Equivalent strain rate [s$^{-1}$]','c')]
    for ax,(key,scale,label,letter) in zip(axes,panels):
        ax.set_facecolor('#FCFDFE');ax.grid(axis='y',color='#DCE3E8',linewidth=.6)
        ax.axvline(parent['time_s']*1e6,color='#606B73',linewidth=.9,linestyle=':',zorder=0)
        for branch,color,marker,_ in styles:
            selected=[row for row in rows if row['branch']==branch]
            ax.plot([row['time_s']*1e6 for row in selected],[row[key]*scale for row in selected],
                color=color,marker=marker,markersize=5.5,linewidth=1.55,markeredgecolor='white',markeredgewidth=.7)
        ax.scatter([parent['time_s']*1e6],[parent[key]*scale],c='black',s=42,zorder=8,edgecolors='white',linewidths=.7)
        ax.set_ylabel(label,labelpad=9);ax.text(-.085,1.02,letter,transform=ax.transAxes,fontsize=12,weight='bold')
        ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    axes[1].axhline(1.,color='#6B6574',linestyle=(0,(5,3)),linewidth=1.)
    axes[1].text(70.2,1.005,r'comparison denominator floor $10^{-5}$',ha='right',va='bottom',fontsize=9,color='#5D5664')
    axes[2].axhline(1.,color='#6B6574',linestyle=(0,(5,3)),linewidth=1.)
    axes[2].text(70.2,1.007,r'material low-rate floor $1$ s$^{-1}$',ha='right',va='bottom',fontsize=9,color='#5D5664')
    handles=[Line2D([0],[0],color=color,marker=marker,lw=1.5,label=label) for _,color,marker,label in styles]
    handles.append(Line2D([0],[0],color='black',marker='o',lw=0,label='Common accepted parent'))
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.11,.875),ncol=2,frameon=False,
        columnspacing=2.4,handlelength=2.5,fontsize=10)
    axes[-1].set_xlabel('Physical time [μs]');axes[-1].set_xlim(62.15,70.65)
    axes[-1].set_xticks([62.5,64.453125,66.40625,68.359375,70.3125],['62.5','64.453','66.406','68.359','70.3125'])
    e0,e1=(v*1e6 for v in data['terminal_global_adjacent_temperature_errors_K'])
    fig.text(.115,.126,f'Temperature error is not monotone: terminal adjacent max |ΔT| = {e0:.3f} → {e1:.3f} μK.',
        fontsize=10.5,color='#9B3E29',weight='bold')
    fig.text(.115,.096,'Straight segments join stored endpoints only. Black dots show the common saved initial state at 62.5 μs.',fontsize=9.6,color='#435663')
    fig.text(.115,.070,'Rate values are unmodified outputs; dashed floors are model/comparison conventions. No extrapolated solution is plotted.',fontsize=9.6,color='#435663')
    fig.text(.115,.044,'Not an experimental validation, accepted continuation, or material-flow / deposition trajectory.',fontsize=9.6,color='#435663')
    png=OUT/'pid797_contact_onset_time_study.png';pdf=OUT/'pid797_contact_onset_time_study.pdf'
    fig.savefig(png,facecolor='white');fig.savefig(pdf,facecolor='white');plt.close(fig)
    receipt=dict(status='RENDERED_SAVED_DIAGNOSTIC_SCIENTIFIC_PLOT',command=[sys.executable,*sys.argv],
        matplotlib_version=matplotlib.__version__,input_data_sha256=sha(source),renderer_sha256=sha(Path(__file__)),
        artifacts_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in (png,pdf)},
        physical_evaluations=0,interpolated_physical_fields=False,plot_segments='CONNECT_SAVED_POINTS_ONLY')
    (OUT/'render_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'],png=str(png),pdf=str(pdf))))


if __name__=='__main__':main()

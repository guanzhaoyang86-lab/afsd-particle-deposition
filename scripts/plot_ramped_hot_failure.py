"""Render actual saved startup diagnostics; no solver/PySPH import."""
from pathlib import Path
import hashlib
import json
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/ramped_hot_startup_20261005/v1/visualization'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source=OUT/'data.json';data=json.loads(source.read_text(encoding='utf-8'))
    assert data['status']=='PASS_SAVED_RAMP_FAILURE_PLOT_EXPORT_NOT_PHYSICAL_VALIDATION'
    assert sha(OUT/'pid791_saved_states.csv')==data['csv_sha256']
    forensic_path=OUT.parent/'method_review/summary.json'
    forensic=json.loads(forensic_path.read_text(encoding='utf-8'))
    assert forensic['status']=='SAVED_ONLY_RAMP_POST_CONTACT_FAILURE_OBSERVATION_NOT_PHYSICAL_VALIDATION'
    assert forensic['maximum_stress_row']==data['terminal_local_screen']
    for branch,history in forensic['histories'].items():
        for step in history:
            if step['time_s']>.0035+2e-18:continue
            particle=next(p for p in step['particles'] if p['particle_id']==791)
            assert not any(v['active'] for v in particle['contact'].values())
            row=next(r for r in data['rows'] if r['branch']==branch and r['time_s']==step['time_s'])
            assert abs(particle['equivalent_stress_Pa']-row['q_Pa'])<1e-9
            assert particle['ep']==row['ep'] and particle['strain_rate_per_s']==row['strain_rate_per_s']
    png=OUT/'pid791_ramp_contact_failure.png';pdf=OUT/'pid791_ramp_contact_failure.pdf'
    assert not png.exists() and not pdf.exists() and not (OUT/'render_receipt.json').exists()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.labelsize':11,
        'axes.titlesize':11,'xtick.labelsize':10,'ytick.labelsize':10,'savefig.dpi':200,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(3,2,figsize=(12.5,10.5),sharex='col',sharey='row',
                          gridspec_kw={'width_ratios':[1.55,1.]})
    fig.subplots_adjust(left=.105,right=.97,top=.80,bottom=.22,hspace=.25,wspace=.12)
    fig.text(.105,.96,'Tool-load onset and time-step failure · PID 791',fontsize=19,weight='bold',color='#152838')
    fig.text(.105,.926,'ENGINEERING_APPROXIMATE · independent 20 ms ramp startup diagnostic',fontsize=11.8,color='#435663')
    stress=data['terminal_local_screen']['stress_relative']*100
    plastic=data['terminal_local_screen']['plastic_relative']*100
    fig.text(.105,.894,f'3.5 ms comparison FAIL: tensor-stress metric {stress:.2f}% (limit 10%); plastic metric {plastic:.2f}% (limit 20%).',
             fontsize=10.6,color='#A13C2D',weight='bold')
    styles=[('coarse','#D55E00','s','Coarse · 500 μs'),('fine','#0072B2','o','Fine · 250 μs')]
    handles=[Line2D([0],[0],color=color,marker=marker,lw=1.6,label=label) for _,color,marker,label in styles]
    handles += [Line2D([0],[0],color='black',marker='o',lw=0,label='Common initial stress / plastic strain'),
                Line2D([0],[0],color='#606B73',lw=1.1,ls=':',label='3.0 ms: first saved system tool load')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.10,.873),ncol=2,frameon=False,
               columnspacing=2.3,handlelength=2.7,fontsize=10)
    initial=next(r for r in data['rows'] if r['branch']=='independent_initial')
    panels=[('q_Pa',1e-6,r'Equivalent stress $q$ [MPa]'),
            ('ep',1e6,r'Plastic strain $\bar\epsilon^p$ [$10^{-6}$]'),
            ('strain_rate_per_s',1.,r'Equivalent strain rate [s$^{-1}$]')]
    for index,(key,scale,label) in enumerate(panels):
        for col,ax in enumerate(axes[index]):
            ax.set_facecolor('#FCFDFE');ax.grid(axis='y',color='#DCE3E8',linewidth=.6)
            ax.axvline(3.,color='#606B73',linewidth=1.,linestyle=':',zorder=0)
            ax.axvline(3.5,color='#AE6155',linewidth=.85,linestyle=(0,(4,3)),zorder=0)
            for branch,color,marker,_ in styles:
                selected=[r for r in data['rows'] if r['branch']==branch]
                if key!='strain_rate_per_s':selected=[initial,*selected]
                ax.plot([r['time_s']*1000 for r in selected],[r[key]*scale for r in selected],
                    color=color,marker=marker,markersize=5.0,linewidth=1.5,
                    markeredgecolor='white',markeredgewidth=.7,zorder=3)
            if col==0 and key!='strain_rate_per_s':
                ax.scatter([0.],[0.],color='black',s=34,edgecolors='white',linewidths=.65,zorder=8)
            if col==0:ax.set_ylabel(label,labelpad=8)
            if key=='strain_rate_per_s':
                ax.set_yscale('log');ax.set_ylim(4e-5,2.)
                ax.axhline(1.,color='#6B6574',linestyle=(0,(5,3)),linewidth=1.)
                ax.text(.97,1.09,'material low-rate floor 1 s$^{-1}$',transform=ax.get_yaxis_transform(),
                        ha='right',va='bottom',fontsize=9,color='#5D5664')
            elif key=='ep':
                ax.set_ylim(-.45,11.1);ax.axhline(10.,color='#9A939F',linestyle=(0,(5,3)),linewidth=.8)
                if col==1:ax.text(.97,10.18,r'comparison floor $10^{-5}$',transform=ax.get_yaxis_transform(),
                                 ha='right',va='bottom',fontsize=9,color='#5D5664')
            else:ax.set_ylim(-.07,1.66)
            ax.text(-.095 if col==0 else -.02,1.025,'abcdef'[2*index+col],transform=ax.transAxes,
                    fontsize=12,weight='bold')
    axes[0,0].set_title('Complete saved interval',loc='left',pad=10)
    axes[0,1].set_title('Zoom near contact onset',loc='left',pad=10)
    for row in axes:
        row[0].set_xlim(-.05,3.62);row[1].set_xlim(2.45,3.54)
    axes[-1,0].set_xticks([0,.5,1.,1.5,2.,2.5,3.,3.5])
    axes[-1,1].set_xticks([2.5,2.75,3.,3.25,3.5])
    for ax in axes[-1]:ax.set_xlabel('Independent ramp model time [ms]')
    fig.text(.105,.154,'All displayed states are unaccepted diagnostics. No deposition or accepted continuation is shown.',
             fontsize=10.8,weight='bold',color='#A13C2D')
    fig.text(.105,.123,'Lines join saved points only; both paths stop at 3.5 ms here. The 3.0 ms guide is a saved-load endpoint, not an exact onset time.',
             fontsize=9.9,color='#435663')
    fig.text(.105,.099,'PID 791 itself has no active contact. Its saved material response differs after neighboring particles engage the tool.',
             fontsize=9.9,color='#435663')
    fig.text(.105,.075,r'Terminal comparison floors (1 MPa / $10^{-5}$) are active. Rates are unmodified on a log axis; unevaluated initial scratch rate is omitted.',
             fontsize=9.9,color='#435663')
    fig.text(.105,.048,'This independent study adds 0 accepted time to the pre-existing 68.359375 μs process branch.',
             fontsize=9.9,color='#435663')
    fig.savefig(png,facecolor='white');fig.savefig(pdf,facecolor='white');plt.close(fig)
    receipt=dict(status='RENDERED_SAVED_RAMP_DIAGNOSTIC_FAILURE_SCIENTIFIC_PLOT',command=[sys.executable,*sys.argv],
        matplotlib_version=matplotlib.__version__,input_data_sha256=sha(source),renderer_sha256=sha(Path(__file__)),
        independent_forensic_path=str(forensic_path.relative_to(ROOT)),independent_forensic_sha256=sha(forensic_path),
        artifacts_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in (png,pdf)},
        physical_evaluations=0,nnps_evaluations=0,interpolated_physical_fields=False,
        plotted_time_range_s=[0.,.0035],unpaired_coarse_4ms_excluded=True)
    (OUT/'render_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'],png=str(png),pdf=str(pdf))))


if __name__=='__main__':main()

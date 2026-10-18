"""Independent saved-only readback of the bounded contact-onset time study.

Reuse frozen reader arithmetic, never the physical controller or its globals.
All diagnostic levels stay unaccepted even if their temporal screens pass.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import readback_expanded_hot_campaign_window as prior

ROOT=prior.ROOT
old=prior.old
Reader,require,clean,digest=prior.Reader,prior.require,prior.clean,prior.digest
path=prior.project_path
BASELINE=ROOT/'outputs/expanded_hot_campaign_20261005/window_001'
BASELINE_AUDIT=BASELINE/'readback/summary.json'
BASELINE_AUDIT_SHA='89f999f398a9f3de70535500121bbba3a359ab3b39cbd351a0e8b00960e88658'
PARENT_SHA='12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b'
CONTROLLER=ROOT/'scripts/run_contact_onset_time_study.py'
CONTRACT=ROOT/'docs/CONTACT_ONSET_TIME_STUDY_CONTRACT_20261005.md'
COMPLETE='COMPLETED_CONTACT_ONSET_RESOLUTION_DIAGNOSTIC_NOT_ACCEPTED'
COUNTS={'four':4,'eight':8}
SPAN=7.8125e-6


def bind_tree(reader,folder):
    return {f.relative_to(ROOT).as_posix():digest(reader.bind(f))
        for f in sorted(folder.rglob('*')) if f.is_file()}


def baseline(reader):
    audit=reader.read(reader.bind(BASELINE_AUDIT,BASELINE_AUDIT_SHA))
    require(audit['status']==prior.PASS_READBACK and audit['exit_code']==0 and
        audit['added_accepted_steps']==0 and audit['accepted_global_steps']==14 and
        audit['window_completed'] is False and audit['stages']['fine']['status']=='FAILED_NO_NEW_ACCEPTED_STATE',
        'original failed comparison was changed/promoted')
    reader.sources(audit['evidence_sha256'])
    context=prior.parent_context(reader,BASELINE,reader.read(BASELINE/'window_protocol.json'))
    require(context['parent_sha256']==PARENT_SHA and context['parent']['step']==14 and
        abs(context['parent']['time_s']-62.5e-6)<1e-15,'incorrect real accepted origin')
    previous,state,direct=context['initial'],context['parent'],PARENT_SHA
    rows=reader.read(BASELINE/'fine/candidate_history.json');cases=[]
    require(len(rows)==2,'saved initial resolution is incomplete')
    for number,row in enumerate(rows,1):
        require(row['step']==number and row['parent_sha256']==direct,'baseline uncommitted reference chain differs')
        item=prior.candidate_audit(reader,'fine',BASELINE/'fine',row,previous,state,context,context['frozen'])
        previous,state,direct=item['pa'],item['state'],row['sha256'];cases.append(item)
    measurements=list(context['protocol']['budget']['measurements'])
    for stage in ('coarse','fine'):
        folder=BASELINE/stage;meta=reader.read(folder/'manifest.json');rows=reader.read(folder/'candidate_history.json')
        require(meta['candidate_steps']==len(rows)==2,'baseline cost scope changed')
        measurements.append(dict(stage='baseline_'+stage,
            maximum_complete_candidate_s=max(r['wall_seconds'] for r in rows),
            stage_wall_per_complete_candidate_s=meta['wall_seconds']/len(rows)))
    return context,audit,cases,measurements


def compare_screens(reader,folder,cases,reference,context):
    seeds=set(context['seeds']);screens=[]
    for pair in range(1,len(cases)//2+1):
        a,b=cases[2*pair-2:2*pair];c=reference[pair-1]
        for item in (a,b,c):
            seeds.update(v['particle_id'] for v in item['row']['events'])
            seeds.update(item['row']['blocked_particle_ids'])
            seeds.update(item['state']['bore_observation']['crossing_ids'])
        selected=set(seeds)
        for item in (a,b,c):
            graph=item['graph'];ids=graph['particle_id']
            for left,right in (('i','j'),('raw_i','raw_j')):
                x,y=ids[graph[left]],ids[graph[right]]
                mask=np.isin(x,list(seeds))|np.isin(y,list(seeds))
                selected.update(np.r_[x[mask],y[mask]].tolist())
        glob,local=old.compare_arrays(b['pa'],c['pa'],sorted(selected),context['p'])
        absolute=[]
        for k in np.flatnonzero(np.isin(b['pa'].particle_id,sorted(selected))):
            fine,coarse=b['pa'],c['pa']
            absolute.append(dict(particle_id=int(fine.particle_id[k]),
                stress_norm_difference_Pa=float(np.linalg.norm(fine.S.reshape(-1,9)[k]-coarse.S.reshape(-1,9)[k])),
                fine_stress_norm_Pa=float(np.linalg.norm(fine.S.reshape(-1,9)[k])),
                plastic_strain_difference=float(abs(fine.ep[k]-coarse.ep[k])),fine_plastic_strain=float(fine.ep[k])))
        crossing=b['state']['bore_observation']['crossing_ids']==c['state']['bore_observation']['crossing_ids']
        require(abs(b['state']['time_s']-c['state']['time_s'])<1e-15,'reference/shared endpoint differs')
        expected=dict(pair=pair,time_s=b['state']['time_s'],parent_sha256=PARENT_SHA,
            coarse_reference_sha256=c['row']['sha256'],fine_candidate_sha256=[a['row']['sha256'],b['row']['sha256']],
            cumulative_event_and_blocked_ids=sorted(seeds),global_screen=glob,local_screen=local,
            local_absolute_differences=absolute,crossing_ids_match=crossing,
            passed=glob['pass'] and local['passed'] and crossing,accepted=False,restart_allowed=False)
        screen=reader.read(folder/f'pair_{pair:03d}_screen.json')
        require(clean(expected)==screen,'shared endpoint screen does not independently replay')
        screens.append(screen)
    require(len(list(folder.glob('pair_*_screen.json')))==len(screens),'orphan/unbound temporal screen')
    return screens


def reduction(reader,stage,base,current):
    pp=BASELINE/'fine/pair_001_temporal_screen.json' if stage=='four' else base/'four/pair_002_screen.json'
    previous=reader.read(pp)
    def metrics(screen):
        result={k:screen['global_screen'][k] for k in ('position_max_m','temperature_max_K','stress_relative','plastic_relative')}
        result.update({'local_'+key:max(r[key] for r in screen['local_screen']['rows']) for key in ('stress_relative','plastic_relative')})
        return result
    before,after=metrics(previous),metrics(current)
    return dict(previous_screen=pp.relative_to(ROOT).as_posix(),previous_screen_sha256=digest(pp),
        previous_errors=before,current_errors=after,
        current_over_previous={k:after[k]/before[k] if before[k]>0 else None for k in before},
        identical_local_id_set=previous['local_screen']['particle_ids']==current['local_screen']['particle_ids'],
        denominator_floors=dict(stress_Pa=1e6,plastic_strain=1e-5),
        scope='EMPIRICAL_ADJACENT_RESOLUTION_DIFFERENCES_NOT_ASSUMED_ORDER_OR_EXTRAPOLATED_SOLUTION')


def audit_stage(reader,stage,base,context,audit,reference,measurements,previous_result):
    folder=base/stage;protocol=reader.read(folder/'protocol.json');meta=reader.read(folder/'manifest.json')
    supervision=reader.read(folder/'supervision.json')
    reader.bind(folder/'run.log',supervision['log_sha256']);reader.bind(folder/'manifest.json',supervision['manifest_sha256'])
    require(supervision['worker_limit_s']==540 and supervision['external_limit_s']==600 and
        supervision['wall_seconds']<=605 and supervision['source_inputs_still_match'] is True and
        supervision['source_error'] is None and supervision['accepted'] is False and supervision['restart_allowed'] is False,
        'diagnostic supervision/limits/source integrity differ')
    expected_command=[supervision['command'][0],'-B',str(CONTROLLER),'--stage',stage,'--output',base.relative_to(ROOT).as_posix(),'--worker']
    require(supervision['command']==expected_command,'diagnostic worker command differs')
    refs=dict(audit['evidence_sha256']);refs[BASELINE_AUDIT.relative_to(ROOT).as_posix()]=BASELINE_AUDIT_SHA
    reference_folder=BASELINE/'fine'
    if stage=='eight':
        require(previous_result['physical_trajectory_completed'],'eight started before four physically completed')
        reference_folder=base/'four';refs.update(bind_tree(reader,reference_folder))
        rows=reader.read(reference_folder/'candidate_history.json');previous_meta=reader.read(reference_folder/'manifest.json')
        measurements=measurements+[dict(stage='four',maximum_complete_candidate_s=max(r['wall_seconds'] for r in rows),
            stage_wall_per_complete_candidate_s=previous_meta['wall_seconds']/len(rows))]
    frozen=dict(context['frozen'])
    for p in (CONTROLLER,CONTRACT,BASELINE_AUDIT):
        reader.bind(p);name=p.relative_to(ROOT).as_posix();value=digest(p)
        require(name not in frozen or frozen[name]==value,'source union conflicts');frozen[name]=value
    measured=max(max(r['maximum_complete_candidate_s'],r['stage_wall_per_complete_candidate_s']) for r in measurements)
    count=COUNTS[stage];cost=measured+1.;forecast=30+1.5*count*cost
    expected=dict(schema='BOUNDED_CONTACT_ONSET_TIME_STUDY_V1',stage=stage,
        source_sha256=frozen,reference_sha256=refs,parent_checkpoint=context['parent_path'].relative_to(ROOT).as_posix(),parent_sha256=PARENT_SHA,
        independent_origin_sha256=prior.INDEPENDENT_SHA,parent_time_s=context['parent']['time_s'],parent_global_step=14,
        baseline_readback=BASELINE_AUDIT.relative_to(ROOT).as_posix(),baseline_readback_sha256=BASELINE_AUDIT_SHA,
        reference_directory=reference_folder.relative_to(ROOT).as_posix(),reference_count=len(reference),step_count=count,dt_s=SPAN/count,
        endpoint_s=context['parent']['time_s']+SPAN,measurements=measurements,measured_cost_s=measured,per_candidate_cost_s=cost,
        forecast_s=forecast,forecast_limit_s=450.,internal_limit_s=540.,external_limit_s=600.,
        status='REGISTERED_DIAGNOSTIC' if forecast<=450 else 'BLOCKED_COST_NO_RUN',
        accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.,
        physical_evaluations_in_preflight=0,nnps_evaluations_in_preflight=0)
    require(protocol==expected and protocol['status']=='REGISTERED_DIAGNOSTIC','immutable protocol/cost/ancestry does not replay')
    reader.sources(frozen);reader.sources(refs)
    for key,value in protocol.items():
        if key!='status':require(meta[key]==value,'manifest changed protocol field: '+key)
    complete=meta['status']==COMPLETE
    require((supervision['exit_code']==0)==complete and meta['status'] in
        (COMPLETE,'FAIL_CONTACT_ONSET_DIAGNOSTIC_NOT_ACCEPTED','RUNNING_DIAGNOSTIC'),
        'status/exit code differs')
    if meta['status']=='RUNNING_DIAGNOSTIC':
        require(supervision['exit_code']==124,'nonterminal diagnostic cannot be audited')
        reader.read(folder/'HARD_TIMEOUT.json')
    else:require(meta['wall_seconds']<=540.5,'diagnostic exceeded internal budget')
    rows=reader.read(folder/'candidate_history.json');cases=[]
    require(meta['candidate_steps']==len(rows)<=count,'diagnostic candidate count differs')
    require(not list(folder.rglob('accepted.npz')),'diagnostic published an accepted checkpoint')
    dynamic=dict(context,fine_dt=protocol['dt_s'])
    previous,state,direct=context['initial'],context['parent'],PARENT_SHA
    for number,row in enumerate(rows,1):
        require(row['step']==number and row['parent_sha256']==direct,'diagnostic immediate-parent hash chain differs')
        item=prior.candidate_audit(reader,'fine',folder,row,previous,state,dynamic,frozen)
        require(item['state']['branch_origin_sha256']==prior.INDEPENDENT_SHA and
            item['state']['window_origin_sha256']==context['parent']['window_origin_sha256'],
            'unpublished diagnostic changed inherited origin metadata')
        previous,state,direct=item['pa'],item['state'],row['sha256'];cases.append(item)
    screens=compare_screens(reader,folder,cases,reference,context)
    require(meta['completed_screen_count']==len(screens),'manifest screen count differs')
    artifacts={f.relative_to(ROOT).as_posix():digest(reader.bind(f)) for f in sorted(folder.rglob('*'))
        if f.is_file() and f.name not in ('manifest.json','run.log','supervision.json')}
    require(meta['artifacts_sha256']==artifacts,'complete/failed stage evidence ledger differs')
    if complete:
        require(len(cases)==count and meta['physical_trajectory_completed'] is True and
            meta['candidate_time_s']==state['time_s'] and abs(state['time_s']-protocol['endpoint_s'])<1e-15 and
            meta['pair_screens']==screens and meta['temporal_screens_all_pass']==all(s['passed'] for s in screens) and
            meta['interval_energy']==clean(state['interval_energy']) and
            meta['crossing_ids']==state['bore_observation']['crossing_ids'] and
            path(meta['candidate_end_checkpoint'])==folder/rows[-1]['checkpoint'] and meta['candidate_end_sha256']==direct,
            'completed trajectory/energy/screen endpoint differs')
        require(meta['terminal_error_reduction']==reduction(reader,stage,base,screens[-1]),'empirical adjacent errors do not replay')
    else:require(meta['physical_trajectory_completed'] is False,'failed stage claims completion')
    partial=[f.relative_to(ROOT).as_posix() for f in sorted(folder.glob('step_*'))
        if f.is_dir() and f.name not in {f"step_{r['step']:03d}" for r in rows}]
    particle=[]
    for item in cases:
        pa=item['pa'];i=int(np.flatnonzero(pa.particle_id==797)[0])
        j=reader.read(path(item['summary']['checkpoint']).parent/'FULL_UNCOMMITTED_candidate.json')
        feature=j['contact']['features'].index('BORE_Y_PLUS')
        particle.append(dict(time_s=item['state']['time_s'],ep=float(pa.ep[i]),ep_rate_per_s=float(pa.ep_rate[i]),
            strain_rate_per_s=float(pa.strain_rate[i]),stress_norm_Pa=float(np.linalg.norm(pa.S.reshape(-1,9)[i])),
            temperature_K=float(pa.temperature[i]),bore_y_plus_active=j['contact']['active'][i][feature],
            bore_y_plus_sliding=j['contact']['sliding'][i][feature],
            bore_y_plus_normal_load_N=j['contact']['normal_load_N'][i][feature]))
    result=dict(stage=stage,status=meta['status'],exit_code=supervision['exit_code'],
        physical_trajectory_completed=complete,temporal_screens_all_pass=all(s['passed'] for s in screens) if complete else False,
        forecast_s=forecast,wall_seconds=meta['wall_seconds'],external_wall_seconds=supervision['wall_seconds'],
        candidate_count=len(cases),accepted_steps=0,accepted_time_advanced_s=0.,
        diagnostic_candidate_end_time_s=state['time_s'],diagnostic_global_step=state['step'],
        full_prefix_preserved=True,complete_interval_energy=clean(state['interval_energy']),
        screens=screens,terminal_error_reduction=meta.get('terminal_error_reduction'),
        candidate_summaries=[v['summary'] for v in cases],PID797_saved_path=particle,
        failed_partial_directories=partial,manifest_sha256=digest(folder/'manifest.json'))
    return cases,result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',required=True);parser.add_argument('--through',choices=tuple(COUNTS),required=True)
    parser.add_argument('--output',required=True);args=parser.parse_args()
    base,output=path(args.base),path(args.output)
    require(base.is_relative_to(ROOT/'outputs') and base.is_dir() and output.is_relative_to(ROOT/'outputs') and
        not output.exists() and not any(output.is_relative_to(base/s) for s in COUNTS),
        'readback requires a fresh output outside physical stage folders')
    output.mkdir(parents=True,exist_ok=False);reader=Reader();started=time.perf_counter()
    report=dict(status='READBACK_RUNNING',physical_solver_evaluations=0,nnps_evaluations=0,
        material_evaluations=0,contact_evaluations=0,new_accepted_steps=0,accepted_time_advanced_s=0.,
        accepted=False,restart_allowed=False,production_authorized=False,command=[sys.executable,*sys.argv])
    try:
        reader.bind(Path(__file__));reader.bind(Path(prior.__file__));reader.bind(Path(old.__file__))
        context,audit,reference,measurements=baseline(reader);results=[];previous_result=None
        for stage in ('four','eight'):
            reference,result=audit_stage(reader,stage,base,context,audit,reference,measurements,previous_result)
            results.append(result);previous_result=result
            if stage==args.through:break
        reader.finish()
        report.update(status='PASS_CONTACT_ONSET_TIME_STUDY_READBACK_NOT_PHYSICAL_VALIDATION',
            source_entry_checks_including_repeats=reader.source_entries,evidence_sha256=reader.evidence,
            stages=results,actual_accepted_parent=context['parent_path'].relative_to(ROOT).as_posix(),
            actual_accepted_parent_sha256=PARENT_SHA,actual_accepted_global_step=14,
            actual_accepted_time_s=context['parent']['time_s'],
            temporal_screen_passes_are_not_publication=True,
            limitations=['Readback PASS verifies saved evidence consistency, including faithfully retained failures.',
                'No diagnostic state is accepted or restartable, even when adjacent temporal screens pass.',
                'Empirical adjacent differences are not experimental error, an assumed convergence order, or an extrapolated solution.',
                'Contact-history, whole-history and spatial convergence remain unestablished.']);code=0
    except Exception:
        report.update(status='FAIL_CONTACT_ONSET_TIME_STUDY_READBACK',error=traceback.format_exc(),evidence_sha256=reader.evidence);code=2
    report.update(exit_code=code,wall_seconds=time.perf_counter()-started)
    (output/'summary.json').write_text(json.dumps(clean(report),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report.get(k) for k in ('status','exit_code','actual_accepted_global_step','actual_accepted_time_s','wall_seconds','error')},indent=2))
    return code


if __name__=='__main__':raise SystemExit(main())

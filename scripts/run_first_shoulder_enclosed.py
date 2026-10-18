"""Explicit adapter: retry one first-shoulder interval with a stronger audit.

All physical trials, scalar equivalence, local screens, and two-half commit
ordering are inherited from the unchanged first-shoulder controller. The
only numerical acceptance extension is a separately tested conservative
plane-curvature proof for the old shell audit's unresolved planar intervals.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time
import run_first_shoulder_event as original
from contact_shell_enclosure import audit_contact_shell_enclosure, CLASSIFICATION as ENCLOSURE

batch,core,post=original.batch,original.core,original.post
BASE=batch.BASE/'first_shoulder_enclosed'
OLD_BASE=batch.BASE/'first_shoulder_event'
CONTRACT=core.ROOT/'docs/FIRST_SHOULDER_ENCLOSED_CONTRACT_20261004.md'
ENCLOSURE_RECEIPT=core.ROOT/'logs/coupled_cost_20261004/shell_enclosure/receipt.json'
ENCLOSURE_READBACK=core.ROOT/'outputs/coupled_cost_20261004/shoulder_shell_enclosure.json'
REVIEW=core.ROOT/'outputs/coupled_cost_20261004/shell_curvature_review.json'
PARENT_SHA='6d5284d82600fe89dac405d2ca26545b448397ec230d08fa43bb3a12c836699a'
FAILED_CANDIDATE_SHA='904714c47961a5325e259a69b61468d9d6213e289c914ab00ab1d8144526fd94'
LEGACY_PREFLIGHT=original.preflight
ORIGINAL_WRITE=core.write
require,checked_file=original.require,original.checked_file
ADAPTER=dict(classification='EXPLICIT_FIRST_SHOULDER_CONTROLLER_ADAPTER_V1',
    reused_controller='scripts/run_first_shoulder_event.py',
    runtime_replacements=['BASE','preflight','audit_contact_shell','MAX_APPROACH_STEPS=1'],
    manifest_annotation='core.write wrapper adds this metadata to new-BASE manifests only',
    acceptance_extension=ENCLOSURE,physical_equations_changed=False,
    old_failed_run_reclassified=False,origin_time_s=.021875,origin_global_step=151,
    maximum_approach_steps=1,approach_dt_s=.000125)


def preflight(stage):
    """Bind the old accepted prefix and old failed candidate; evaluate no force."""
    legacy=LEGACY_PREFLIGHT('approach')
    p,provenance,baseline=legacy['p'],legacy['provenance'],core.frozen_sources()
    old_manifest=OLD_BASE/'approach/manifest.json';old_history=OLD_BASE/'approach/accepted_history.json'
    meta,rows=post.read(old_manifest),post.read(old_history)
    require(meta['status']=='FAILED_WITH_ACCEPTED_PREFIX' and meta['accepted_steps']==len(rows)==15 and
            abs(meta['accepted_time_s']-.021875)<1e-15 and meta['source_sha256']==legacy['frozen'] and
            meta['provenance']==provenance,'old failed approach/prefix/source differs')
    post.check_sources(meta['source_sha256'])
    require(not (OLD_BASE/'approach/event.json').exists(),'old failed run unexpectedly contains an event qualification')
    previous_sha=legacy['parent_hash'];prefix=list(legacy['state']['interval_energy']['records'])
    parent=None;pa=state=None
    references={old_manifest:core.sha(old_manifest),old_history:core.sha(old_history)}
    for number,row in enumerate(rows,1):
        parent=OLD_BASE/'approach'/row['checkpoint']
        require(core.sha(parent)==row['sha256'] and row['parent_sha256']==previous_sha,
                'old accepted direct-parent chain/hash differs')
        pa,state=core.restore(parent)
        require(state['step']==row['global_step']==136+number and state['time_s']==row['time_s'] and
                abs(state['time_s']-(.020+number*original.DT))<1e-15 and row['dt_s']==original.DT and
                state['source_sha256']==baseline and state['implementation_source_sha256']==legacy['frozen'] and
                state['parameters']==p and state['parameter_provenance']==provenance and
                state['continuation_parent_sha256']==previous_sha and state['branch_origin_sha256']==legacy['parent_hash'],
                'old accepted checkpoint metadata differs')
        batch.make_trial(pa,state,p)  # serialized histories only, no evaluate
        require(state['interval_energy']['records'][:-1]==prefix,'old accepted cumulative prefix differs')
        prefix=list(state['interval_energy']['records']);previous_sha=row['sha256'];references[parent]=previous_sha
    require(previous_sha==PARENT_SHA and state['step']==151,'registered 21.875ms parent differs')
    origin,origin_hash=parent,PARENT_SHA
    failure_dir=OLD_BASE/'approach/step_016';failed=failure_dir/'NOT_RESTARTABLE_candidate.npz'
    require(core.sha(failed)==FAILED_CANDIDATE_SHA,'old failed candidate hash differs')
    failed_data=core.load(str(failed));failed_state=failed_data['solver_data']
    require(failed_state['restart_allowed'] is False and failed_state['event_candidate_accepted'] is False and
            failed_state['step']==152 and abs(failed_state['time_s']-.022)<1e-15 and
            failed_state['interval_energy']['records'][:-1]==prefix,
            'old candidate is not the same-parent nonrestartable interval')
    require(not (failure_dir/'accepted.npz').exists(),'old rejected candidate was published')
    old_shell_file=failure_dir/'continuous_shell_audit.json';old_shell=post.read(old_shell_file)
    require(old_shell['status']=='UNRESOLVED' and old_shell['failure_kind']=='MAXIMUM_GEOMETRY_QUERIES' and
            not old_shell['hidden_contacts'],'old failure is not the registered unresolved shell case')
    references.update({failed:FAILED_CANDIDATE_SHA,old_shell_file:core.sha(old_shell_file)})
    receipt,certificate,review=post.read(ENCLOSURE_RECEIPT),post.read(ENCLOSURE_READBACK),post.read(REVIEW)
    require(receipt['exit_code']==0,'new enclosure fixtures failed')
    post.check_sources(receipt['source_sha256'])
    require(certificate['status']==original.CERTIFIED and certificate['physical_evaluations']==0 and
            certificate['additional_finite_patch_queries']==0 and
            certificate['original_audit']==old_shell,'new read-only certificate does not preserve the actual failed audit')
    post.check_sources(certificate['source_sha256'])
    expected_inputs={origin.resolve():PARENT_SHA,failed.resolve():FAILED_CANDIDATE_SHA}
    for evidence in (certificate,review):
        supplied={checked_file(name,digest).resolve():digest for name,digest in evidence['input_sha256'].items()}
        require(all(supplied.get(path)==digest for path,digest in expected_inputs.items()),
                'readback/review does not bind both original checkpoint inputs')
    require(review['scope']=='READ_ONLY_SAVED_POSITIONS_AND_JOURNALS_NO_PHYSICAL_EVALUATION' and
            review['old_status']=='UNRESOLVED' and all(row['plane_lower_clearance_nm']>0 for row in review['rows']),
            'independent read-only planar review absent')
    extra=[Path(__file__).resolve(),CONTRACT,core.ROOT/'scripts/contact_shell_enclosure.py',
           core.ROOT/'tests/review/test_contact_shell_enclosure.py',ENCLOSURE_RECEIPT,ENCLOSURE_READBACK,REVIEW,
           old_manifest,old_history,old_shell_file]
    frozen=dict(legacy['frozen'])
    frozen.update(receipt['source_sha256']);frozen.update(certificate['source_sha256'])
    frozen.update({path.relative_to(core.ROOT).as_posix():core.sha(path) for path in extra})
    references.update({ENCLOSURE_RECEIPT:core.sha(ENCLOSURE_RECEIPT),ENCLOSURE_READBACK:core.sha(ENCLOSURE_READBACK),REVIEW:core.sha(REVIEW)})
    per_root=max(legacy['per_root'],max(row['wall_seconds'] for row in rows),
                 meta['wall_seconds']-sum(row['wall_seconds'] for row in rows))
    event=None
    if stage!='approach':
        approach=post.read(BASE/'approach/manifest.json');event_file=BASE/'approach/event.json';event=post.read(event_file)
        require(approach['status']=='EVENT_CANDIDATE_NOT_ACCEPTED' and approach['source_sha256']==frozen and
                approach['provenance']==provenance and approach['controlled_adapter']==ADAPTER,
                'matching enclosed approach event absent')
        require(core.sha(event_file)==approach['event_record_sha256'] and event['source_sha256']==frozen and
                event['dt_s']==original.DT and event['new_features']==['SHOULDER'] and event['event_particle_ids'],
                'enclosed event provenance/feature differs')
        parent=checked_file(event['parent_checkpoint'],event['parent_sha256'])
        require(parent.resolve()==origin.resolve() and event['parent_sha256']==PARENT_SHA and
                event['parent_time_s']==state['time_s'] and event['parent_global_step']==151,
                'enclosed event must retain exactly the 21.875ms parent')
        references[event_file]=core.sha(event_file)
        references[checked_file(event['candidate_checkpoint'],event['candidate_sha256'])]=event['candidate_sha256']
        per_root=max(per_root,event['wall_seconds'])
        if stage=='halves':
            macro_file=BASE/'macro/manifest.json';macro=post.read(macro_file)
            require(macro['status']=='MACRO_NUMERICALLY_QUALIFIED_NOT_ACCEPTED' and macro['source_sha256']==frozen and
                    macro['parent_sha256']==PARENT_SHA and macro['event_record_sha256']==core.sha(event_file) and
                    macro['controlled_adapter']==ADAPTER and abs(macro['candidate_time_s']-.022)<1e-15,
                    'same-parent enclosed macro qualification absent')
            references[macro_file]=core.sha(macro_file)
            references[checked_file(macro['candidate_checkpoint'],macro['candidate_sha256'])]=macro['candidate_sha256']
            references[checked_file(macro['candidate_evidence'],macro['candidate_evidence_sha256'])]=macro['candidate_evidence_sha256']
            per_root=max(per_root,macro['candidate_wall_seconds'])
    forecast=30+1.5*(2 if stage=='halves' else 1)*per_root
    require(forecast<=450,'enclosed event forecast exceeds 450s')
    return dict(p=p,provenance=provenance,pa=pa,state=state,parent=parent,parent_hash=PARENT_SHA,
                origin=origin,origin_hash=origin_hash,frozen=frozen,per_root=per_root,
                forecast=forecast,event=event,references=references)


def write_with_adapter(path,value):
    destination=Path(path).resolve()
    if destination.name=='manifest.json' and destination.is_relative_to(BASE.resolve()) and isinstance(value,dict):
        value=dict(value,controlled_adapter=ADAPTER)
    return ORIGINAL_WRITE(path,value)


def configure():
    """Process-local adapter only; no old source/evidence file is modified."""
    original.BASE=BASE
    original.preflight=preflight
    original.audit_contact_shell=audit_contact_shell_enclosure
    original.MAX_APPROACH_STEPS=1
    original.core.write=write_with_adapter


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['approach','macro','halves'],required=True)
    parser.add_argument('--worker',action='store_true');parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    configure()
    if args.preflight_only:
        context=preflight(args.stage)
        print(json.dumps(dict(stage=args.stage,forecast_wall_s=context['forecast'],parent_sha256=context['parent_hash'],
            parent_time_s=context['state']['time_s'],parent_global_step=context['state']['step'],
            controlled_adapter=ADAPTER,physical_evaluations=0)))
        raise SystemExit(0)
    if args.worker:raise SystemExit(original.worker(args.stage))
    logs=core.ROOT/'logs/coupled_cost_20261004/batched/first_shoulder_enclosed';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/args.stage).exists(),'refuse enclosed evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    ORIGINAL_WRITE(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,
        wall_seconds=time.perf_counter()-started,limit_s=600,controlled_adapter=ADAPTER))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)

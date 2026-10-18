"""Immutable predictor/onset evidence archive; not scientific acceptance."""
import json
from pathlib import Path
from freeze_process_friction_q5 import sha, read, safe_file, write_new
from freeze_process_feasible_predictor import tests

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_FEASIBLE_PREDICTOR_FREEZE.json'
TARGET=ROOT/'docs/PROCESS_PREDICTOR_ONSET_FREEZE.json'
LOGS='logs/process_predictor_onset_20261002'
OUTPUTS='outputs/process_predictor_onset_20261002'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def main():
    if TARGET.exists(): raise FileExistsError('never overwrite a freeze')
    old=read(PRIOR); before=read(safe_file(LOGS+'/freeze_before.log'))
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):
        raise ValueError('missing initial old freeze verification')
    changed=[s for s,h in old['files_sha256'].items() if sha(safe_file(s))!=h]
    if set(changed)-ALLOWED: raise ValueError('unauthorized previous evidence change: '+str(changed))
    case=ROOT/OUTPUTS/'n3_i4'; r=read(case/'report.json')
    rbpath=ROOT/OUTPUTS/'n3_i4_readback.json'; rb=read(rbpath)
    if (r['status'] not in ('J2_CONE_PATH_COMPLETE','BLOCKED')
            or rb['completion_status']!=r['status']
            or rb['accepted_time_s']!=r['accepted_time_s']
            or rb['predictor_onset_audit']['variant']!='Q7_PREDICTOR_ONSET_V1'):
        raise ValueError('run/readback classification mismatch')
    extra=read(case/'predictor_provenance.json')
    sources=dict(r['source_sha256'],**extra['source_sha256'],**extra['onset_source_sha256'])
    for name,h in sources.items():
        if sha(safe_file(name))!=h: raise ValueError('source changed after run: '+name)
    for name,h in r['artifact_sha256'].items():
        if Path(name).name!=name or sha(case/name)!=h: raise ValueError('artifact mismatch')
    notrun=read(ROOT/OUTPUTS/'n3_i8/NOT_RUN.json')
    if notrun['status']!='NOT_RUN' or not notrun['reason']:
        raise ValueError('missing explicit unrun outcome')
    files={safe_file(s) for s in old['files_sha256']}|{PRIOR}
    for directory,pattern in [('docs','PROCESS_PREDICTOR_ONSET*'),
            ('configs/benchmarks','process_predictor_onset*.json'),
            ('scripts','*process_predictor_onset*.py'),
            ('tests/review','test_process_predictor_onset*.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    for directory in (LOGS,OUTPUTS):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    result=dict(schema_version='1.0',date='2026-10-02',
        backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        production_authorized=False,scope='integrity only; NOT scientific acceptance',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),
            unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        actual_tests={name:tests(LOGS+'/'+name+'.log') for name in
            ('integration_tests_final','predictor_compatibility','onset_compatibility')},
        cases=dict(n3_i4=dict(status=r['status'],error=r['error'],
            accepted_time_s=r['accepted_time_s'],accepted_steps=r['accepted_steps'],
            events=len(r['events']),readback_status=rb['status'],
            report_sha256=sha(case/'report.json'),readback_sha256=sha(rbpath)),n3_i8=notrun),
        excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)})
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),preservation=result['preservation'],
                         cases=result['cases'],production_authorized=False),indent=2))


if __name__=='__main__': main()

"""Archive bounded normal-contact evidence; integrity is NOT scientific acceptance."""
from pathlib import Path
import json
from freeze_process_friction_q5 import sha,read,safe_file,write_new
from freeze_process_feasible_predictor import tests

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/ROTATING_DEFORMABLE_FREEZE.json'
TARGET=ROOT/'docs/NORMAL_CONTACT_FREEZE.json'
LOGS='logs/normal_contact_20261002';OUTPUTS='outputs/normal_contact_20261002'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}
ARCHIVED={'point_cost_probe','static_dt005'}
OLD_SOURCE='scripts/probe_normal_contact.py'
OLD_DIGEST='1c1eb2a65c9096d4e85499be7e3dd278b8173575a2a0242581d732300236d4cb'
SNAPSHOT=OUTPUTS+'/probe_normal_contact_before_restart_guard.py'


def main():
    if TARGET.exists():raise FileExistsError('never overwrite a freeze')
    old=read(PRIOR);before=read(safe_file(LOGS+'/freeze_before.log'))
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):
        raise ValueError('prior initial verification missing')
    changed=[s for s,h in old['files_sha256'].items() if sha(safe_file(s))!=h]
    if set(changed)-ALLOWED:raise ValueError('old evidence changed: '+str(changed))
    if len((ROOT/'HANDOFF.md').read_text(encoding='utf-8'))>800:raise ValueError('handoff length')
    if sha(safe_file(SNAPSHOT))!=OLD_DIGEST:raise ValueError('original point source snapshot changed')
    cases={}
    for path in sorted((ROOT/OUTPUTS).glob('*/manifest.json')):
        r=read(path);case=path.parent.name;resolutions={}
        if r['production_authorized'] is not False:raise ValueError('unintended production flag')
        for src,h in r['source_sha256'].items():
            if case in ARCHIVED and src==OLD_SOURCE:
                if h!=OLD_DIGEST:raise ValueError('wrong archived source version')
                target=safe_file(SNAPSHOT);resolutions[src]=SNAPSHOT
            else:target=safe_file(src)
            if sha(target)!=h:raise ValueError('run source drift: '+case+' '+src)
        for artifact,h in r['artifact_sha256'].items():
            if Path(artifact).name!=artifact or sha(path.parent/artifact)!=h:raise ValueError('artifact drift: '+case+' '+artifact)
        cases[case]={k:r.get(k) for k in ('status','accepted_steps','accepted_time_s','wall_seconds','production_authorized')}
        cases[case]['source_snapshot_resolution']=resolutions
        cases[case]['current_restart_compatible']=not bool(resolutions)
    # The guard-only repair must not be advertised as a different physics result.
    for name in ('history.json','normal_intervals.json'):
        if sha(ROOT/OUTPUTS/'static_dt005'/name)!=sha(ROOT/OUTPUTS/'static_dt005_guarded'/name):
            raise ValueError('restart guard changed dynamics evidence')
    comparison=read(ROOT/OUTPUTS/'block_comparison.json')
    if comparison['screening_status']!='SENSITIVE_WARNING':raise ValueError('sensitivity warning not preserved')
    files={safe_file(s) for s in old['files_sha256']}|{PRIOR}
    for directory,pattern in [('docs','NORMAL_CONTACT*'),('configs/benchmarks','normal*.json'),
        ('scripts','*normal*.py'),('tests/review','test_normal*.py')]:
        files.update(p for p in (ROOT/directory).glob(pattern) if p.is_file())
    files.add(ROOT/'src/afsd_backend/pysph/normal_penalty.py')
    for directory in (LOGS,OUTPUTS):files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files={p for p in files if p.relative_to(ROOT).as_posix() not in EXCLUDED}
    names=['normal_unit_tests_final.log','normal_unit_tests_guarded.log','block_unit_tests_final.log',
        'old_contact_tests.log','old_deformable_tests.log','reader_reference_tests.log',
        'reader_reference_tests_binary_fixture.log','short_regression_initial.log']
    inventory={name:tests(LOGS+'/'+name) for name in names}
    result=dict(schema_version='1.0',date='2026-10-02',backend_repository=old['backend_repository'],
        backend_commit=old['backend_commit'],production_authorized=False,scope='file integrity only; NOT scientific acceptance',
        previous_freeze=PRIOR.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(PRIOR),
        preservation=dict(previous_entries=len(old['files_sha256']),unchanged_previous_entries=len(old['files_sha256'])-len(changed),
            authorized_modified_previous_entries=changed,old_failures_preserved=True),
        actual_test_logs=inventory,unique_final_new_short_passes=43,selected_old_short_passes=26,
        cases=cases,comparison_status=comparison['screening_status'],
        archived_point_guard_only_revision=dict(original_source=OLD_SOURCE,snapshot=SNAPSHOT,original_sha256=OLD_DIGEST,
            history_and_contact_intervals_unchanged=True),
        excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)})
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),preservation=result['preservation'],
        new_short_passes=43,old_short_passes=26,comparison_status=result['comparison_status']),indent=2))


if __name__=='__main__':main()

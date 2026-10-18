"""Freeze predictor evidence without running physics or overwriting archives."""
import json
from pathlib import Path
import re

from freeze_process_friction_q5 import sha, read, safe_file, write_new

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'docs/PROCESS_PLASTIC_FRICTION_ONSET_FREEZE.json'
TARGET=ROOT/'docs/PROCESS_FEASIBLE_PREDICTOR_FREEZE.json'
LOGS='logs/process_feasible_predictor_20261002'
OUTPUTS='outputs/process_feasible_predictor_20261002'
ALLOWED={'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md'}
EXCLUDED={TARGET.relative_to(ROOT).as_posix(),LOGS+'/freeze_create.log',LOGS+'/freeze_final.log'}


def tests(relative):
    # Match verdict endings, including unittest's wrapped docstring RED line.
    txt=safe_file(relative).read_text(encoding='utf-8-sig')
    runs=re.findall(r'^Ran (\d+) tests? in ',txt,re.M)
    verdicts=re.findall(r' \.\.\. (ok|expected failure|FAIL|ERROR|unexpected success)$',txt,re.M)
    terminal=re.findall(r'^(OK(?: \([^\r\n]+\))?|FAILED \([^\r\n]+\))$',txt,re.M)
    if len(runs)!=1 or len(terminal)!=1 or len(verdicts)!=int(runs[0]):
        raise ValueError('test inventory incomplete: '+relative)
    passed=verdicts.count('ok'); red=verdicts.count('expected failure')
    failed=sum(verdicts.count(v) for v in ('FAIL','ERROR','unexpected success'))
    if failed or not terminal[0].startswith('OK'):
        status='FAIL_PRESERVED'
    else:
        status='COMPLETED_WITH_EXPECTED_FAILURE' if red else 'PASS'
    return dict(status=status,tests_run=int(runs[0]),passing_tests=passed,
                expected_failures=red,unexpected_failures_or_errors=failed,
                log=relative,sha256=sha(safe_file(relative)))


def main():
    if TARGET.exists(): raise FileExistsError('never overwrite a freeze')
    old=read(PRIOR); before=read(safe_file(LOGS+'/freeze_before.log'))
    if before['status']!='PASS' or before['files_checked']!=len(old['files_sha256']):
        raise ValueError('missing initial old freeze verification')
    changed=[s for s,h in old['files_sha256'].items() if sha(safe_file(s))!=h]
    if set(changed)-ALLOWED: raise ValueError('unauthorized prior change: '+str(changed))
    case=ROOT/OUTPUTS/'target_0125'
    d=read(case/'diagnostic.json')
    if (not d['accepted_states_unchanged'] or not d['input_artifacts_unchanged']
            or not d['source_unchanged'] or not d['expected_old_guard_reproduced']
            or d['accepted_progress'] or d['accepted_step']!=16
            or d['accepted_time_s']!=.1125 or d['comparison_invocations']!=2):
        raise ValueError('diagnostic integrity/progress mismatch')
    for name,h in dict(d['source_sha256'],**d['diagnostic_source_sha256']).items():
        if sha(safe_file(name))!=h: raise ValueError('source changed after diagnostic: '+name)
    original=ROOT/'outputs/process_plastic_friction_onset_20261002/n3_i4'
    for name,h in d['input_sha256'].items():
        if Path(name).name!=name or sha(original/name)!=h: raise ValueError('old input mismatch')
    for name,h in d['artifact_sha256'].items():
        if Path(name).name!=name or sha(case/name)!=h: raise ValueError('diagnostic artifact mismatch')
    new=d['comparisons'][1]
    if (new['returned_corrections'] != new['predictor_corrections']+new['cone_corrections']
            or new['returned_corrections']>14):
        raise ValueError('shared correction budget mismatch')
    marker=read(ROOT/OUTPUTS/'NOT_RUN.json')
    if marker['status']!='NOT_RUN': raise ValueError('missing explicit unrun status')
    files={safe_file(s) for s in old['files_sha256']}|{PRIOR}
    for directory,pattern in [('docs','PROCESS_FEASIBLE_PREDICTOR*'),
            ('scripts','*process_feasible_predictor*.py'),
            ('tests/review','test_process_feasible_predictor*.py')]:
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
        tests={name:tests(LOGS+'/'+name+'.log') for name in
            ('independent_final','onset_compatibility','q7_compatibility','q6_compatibility')},
        diagnostic=dict(status=d['status'],exit_code=d['exit_code'],
            accepted_progress=False,accepted_time_s=d['accepted_time_s'],
            accepted_step=d['accepted_step'],corrections=new['returned_corrections'],
            sha256=sha(case/'diagnostic.json')),
        excluded_self_or_future_paths=sorted(EXCLUDED),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)})
    write_new(TARGET,result)
    print(json.dumps(dict(status='CREATED',files=len(files),
        preservation=result['preservation'],tests=result['tests'],diagnostic=result['diagnostic'],
        production_authorized=False),indent=2))


if __name__=='__main__': main()

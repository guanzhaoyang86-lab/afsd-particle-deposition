"""Create a NEW authorized evidence freeze; never rewrite the previous freeze."""
import argparse
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
ALLOWED={
'src/afsd_backend/pysph/application.py','src/afsd_backend/pysph/output.py',
'src/afsd_backend/pysph/solver.py','docs/SCIENTIFIC_REVIEW.md',
'docs/VALIDATION_CONTRACT.md','docs/IMPLEMENTATION_PLAN.md',
'docs/EQUATION_CODE_TEST_MAP.md','docs/OUTPUT_SCHEMA.md','HANDOFF.md'}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='docs/COUPLING_REVIEW_FREEZE.json')
    args=parser.parse_args();target=ROOT/args.output
    if target.exists():raise FileExistsError('never overwrite a freeze')
    oldpath=ROOT/'docs/SCIENTIFIC_FREEZE.json'
    old=json.loads(oldpath.read_text(encoding='utf-8'))
    changes=[k for k,v in old['files_sha256'].items() if sha(ROOT/k)!=v]
    unexpected=set(changes)-ALLOWED
    if unexpected:raise RuntimeError('unexpected old-evidence/physics edits: '+str(unexpected))
    unchanged_core={}
    for name in ('contact.py','material.py','mechanics.py','geometry.py','integrator.py'):
        before=ROOT/'logs/coupling_20260930/source_before'/name
        now=ROOT/'src/afsd_backend/pysph'/name
        unchanged_core[name]=sha(before)==sha(now)
    if not all(unchanged_core.values()):raise RuntimeError('out-of-scope core change')
    preservation=dict(status='PASS',previous_freeze_sha256=sha(oldpath),
        previous_entries=len(old['files_sha256']),authorized_modified_entries=changes,
        unmodified_previous_entries=len(old['files_sha256'])-len(changes),
        unchanged_core=unchanged_core,all_original_failure_logs_retained=True,
        scope='byte integrity, NOT scientific acceptance')
    report=ROOT/'logs/coupling_20260930/preservation_report.json'
    if report.exists():raise FileExistsError(report)
    report.write_text(json.dumps(preservation,indent=2)+'\n',encoding='utf-8')
    files={ROOT/k for k in old['files_sha256']}
    files.add(oldpath)
    for directory in ('src/afsd_backend','tests/review','logs/coupling_20260930',
                      'outputs/coupling_20260930'):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts and not p.name.endswith('.pending.npz'))
    for name in ('COUPLING_FAILURE_REVIEW.md','COUPLING_DIAGNOSTIC_CONTRACT.md','FINITE_STRAIN_VARIANT_DESIGN.md'):
        files.add(ROOT/'docs'/name)
    for name in ('diagnose_first_coupling_failure.py','run_coupling_review.py','freeze_coupling_evidence.py'):
        files.add(ROOT/'scripts'/name)
    data=dict(schema_version='1.0',date='2026-09-30',production_authorized=False,
        backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        scientific_status='APPROXIMATE/BLOCKED; 101 tests: 95 PASS, 6 FAIL',
        scope='authorized diagnostic/contact revision; integrity is NOT scientific acceptance',
        previous_freeze='docs/SCIENTIFIC_FREEZE.json',previous_freeze_sha256=sha(oldpath),
        post_trial_code_changes='iteration diagnostics, atomic checkpoint and conservative rejected trace metadata only; no second 2us run',
        files_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sorted(files)})
    target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='CREATED',files=len(files),preservation=preservation),indent=2))
if __name__=='__main__':main()


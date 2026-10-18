"""New integrity snapshot; old evidence and non-contact physics must be intact."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ALLOWED={
    'HANDOFF.md','docs/IMPLEMENTATION_PLAN.md','docs/OUTPUT_SCHEMA.md',
    'docs/SCIENTIFIC_REVIEW.md','docs/VALIDATION_CONTRACT.md',
    'scripts/run_coupling_review.py','src/afsd_backend/pysph/contact_kick.py',
    'tests/review/test_contact454_precision.py','tests/review/test_contact439_open_risk.py',
}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target=ROOT/'docs/METHOD_CLASS_FREEZE.json'
    report=ROOT/'logs/contact_class_20260930/preservation_report.json'
    if target.exists() or report.exists():raise FileExistsError('never overwrite a freeze/report')
    prior=ROOT/'docs/CONTACT454_REPAIR_FREEZE.json'
    old=json.loads(prior.read_text(encoding='utf-8'))
    changes=[key for key,value in old['files_sha256'].items() if sha(ROOT/key)!=value]
    unexpected=set(changes)-ALLOWED
    if unexpected:raise RuntimeError('out-of-scope edits: '+str(unexpected))
    for key,name in (
        ('src/afsd_backend/pysph/contact_kick.py','contact_kick.py'),
        ('tests/review/test_contact454_precision.py','test_contact454_precision.py')):
        if sha(ROOT/'logs/contact_class_20260930/source_before'/name)!=old['files_sha256'][key]:
            raise RuntimeError('historical snapshot mismatch: '+key)
    core={k:sha(ROOT/k)==v for k,v in old['files_sha256'].items()
          if k.startswith(('configs/','src/')) and k not in ALLOWED}
    if not all(core.values()):raise RuntimeError('non-contact physics/config modified')
    summary=json.loads((ROOT/'logs/contact_class_20260930/final_suite_with_holdout/summary.json').read_text())
    for key,value in summary['source_hashes'].items():
        if sha(ROOT/key)!=value:raise RuntimeError('source edited after final test: '+key)
    matrix=json.loads((ROOT/'logs/contact_class_20260930/final_suite_with_holdout/artifacts/contact_class_matrix.json').read_text())
    preservation=dict(status='PASS',scope='byte integrity only, NOT scientific acceptance',
        previous_entries=len(old['files_sha256']),authorized_modified_entries=changes,
        unmodified_previous_entries=len(old['files_sha256'])-len(changes),
        preserved_old_evidence_entries=sum(k.startswith(('outputs/','logs/')) for k in old['files_sha256']),
        unchanged_non_contact_sources_and_parameters=core,old_failure_records_preserved=True,
        tested_source_unchanged=True,real_afsd_run='NOT_RUN')
    report.write_text(json.dumps(preservation,indent=2)+'\n',encoding='utf-8')
    files={ROOT/k for k in old['files_sha256']};files.add(prior)
    for directory in ('src/afsd_backend','tests/review','logs/contact_class_20260930',
                      'outputs/contact_class_20260930'):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts and not p.name.endswith('.pending.npz'))
    for name in ('docs/CONTACT_CLASS_CONTRACT.md','docs/METHOD_CLASS_REVIEW.md',
                 'scripts/run_contact_class.py','scripts/probe_method_boundaries.py',
                 'scripts/freeze_method_class_evidence.py'):
        files.add(ROOT/name)
    result=dict(schema_version='1.0',date='2026-09-30',production_authorized=False,
        backend_repository=old['backend_repository'],backend_commit=old['backend_commit'],
        scientific_status='APPROXIMATE/BLOCKED; local contact fixed-matrix PASS; material/boundary unresolved',
        tests_run=summary['tests_run'],failures=summary['failures'],errors=summary['errors'],
        cone_cases=len(matrix),cone_failures=sum(a['status']!='PASS' for a in matrix),
        previous_freeze=prior.relative_to(ROOT).as_posix(),previous_freeze_sha256=sha(prior),
        files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)})
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='CREATED',files=len(files),preservation=preservation),indent=2))


if __name__=='__main__':main()

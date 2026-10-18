"""Freeze this bounded repair, checking every prior evidence entry first."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ALLOWED = {
    'src/afsd_backend/pysph/contact_kick.py',
    'tests/review/test_contact_events.py',
    'tests/review/test_event_solver_smoke.py',
    'tests/review/test_new_candidate_red.py',
    'scripts/run_coupling_review.py',
    'docs/COUPLING_FAILURE_REVIEW.md', 'docs/IMPLEMENTATION_PLAN.md',
    'docs/SCIENTIFIC_REVIEW.md', 'docs/OUTPUT_SCHEMA.md', 'HANDOFF.md',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = ROOT/'docs/CONTACT454_REPAIR_FREEZE.json'
    report = ROOT/'logs/contact454_20260930/preservation_report.json'
    if target.exists() or report.exists():
        raise FileExistsError('never overwrite an evidence freeze/report')
    oldpath = ROOT/'docs/COUPLING_REVIEW_FREEZE.json'
    old = json.loads(oldpath.read_text(encoding='utf-8'))
    changes = [k for k, v in old['files_sha256'].items() if sha(ROOT/k) != v]
    unexpected = set(changes)-ALLOWED
    if unexpected:
        raise RuntimeError('out-of-scope changes: '+str(unexpected))
    core = {k: sha(ROOT/k) == v for k, v in old['files_sha256'].items()
            if (k.startswith('src/') and k not in ALLOWED) or k.startswith('configs/')}
    evidence = [k for k in old['files_sha256']
                if k.startswith(('outputs/', 'logs/'))]
    if not all(core.values()):
        raise RuntimeError('unexpected core or parameter changes')
    before = ROOT/'logs/contact454_20260930/source_before/contact_kick.py'
    if sha(before) != old['files_sha256']['src/afsd_backend/pysph/contact_kick.py']:
        raise RuntimeError('pre-repair source backup mismatch')
    run = ROOT/'outputs/contact454_20260930/fixed_dx05_cfl01'
    run_hashes = json.loads((run/'source_hashes.json').read_text(encoding='utf-8'))
    # No physics-source edit after the real trial is allowed in this freeze.
    if not isinstance(run_hashes, dict):
        raise TypeError('unexpected run hash schema')
    for key, expected in run_hashes.items():
        if sha(ROOT/key) != expected:
            raise RuntimeError('post-trial source change: '+key)
    preservation = dict(
        status='PASS', previous_entries=len(old['files_sha256']),
        authorized_modified_entries=changes,
        unmodified_previous_entries=len(old['files_sha256'])-len(changes),
        unchanged_core=core, preserved_old_evidence_entries=len(evidence),
        all_original_failure_logs_retained=True,
        physics_source_unchanged_since_trial=True,
        scope='byte integrity only, NOT scientific acceptance')
    report.write_text(json.dumps(preservation, indent=2)+'\n', encoding='utf-8')
    files = {ROOT/k for k in old['files_sha256']}
    files.add(oldpath)
    for directory in ('tests/review', 'logs/contact454_20260930',
                      'outputs/contact454_20260930'):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts and not p.name.endswith('.pending.npz'))
    for path in ('docs/CONTACT454_REPAIR_CONTRACT.md',
                 'docs/CONTACT454_REPAIR_REPORT.md',
                 'scripts/probe_contact454.py', 'scripts/freeze_contact454_evidence.py'):
        files.add(ROOT/path)
    data = dict(
        schema_version='1.0', date='2026-09-30', production_authorized=False,
        backend_repository=old['backend_repository'], backend_commit=old['backend_commit'],
        scientific_status='APPROXIMATE/BLOCKED; 107 tests: 101 PASS, 6 FAIL',
        scope='454 bounded numerical repair; new 439 multi-contact blocker',
        previous_freeze=str(oldpath.relative_to(ROOT)).replace('\\', '/'),
        previous_freeze_sha256=sha(oldpath),
        files_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files)})
    target.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status='CREATED', files=len(files), preservation=preservation), indent=2))


if __name__ == '__main__':
    main()

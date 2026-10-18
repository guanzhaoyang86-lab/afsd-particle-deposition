"""Read-only diagnostic of a rejected saved chord; never solve or accept state."""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'src')]
from visible_support_remote_domain_20261007 import audit_visible_support_path
from afsd_backend.pysph.quasistatic_contact import _rotation_path
from visible_support import _hash_array, _topology_hash


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def save(path, value):
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


base = ROOT / 'outputs/gpu_acceleration_20261007'
source = base / 'incremental_deposition_resume_20261010_v1/attempt_0019'
coarse = source / 'coarse'
target = Path(__file__).parent
manifest = json.loads((source / 'manifest.json').read_text(encoding='utf-8'))
original = json.loads((coarse / 'visible_path_audit.json').read_text(encoding='utf-8'))
motion = json.loads((coarse / 'motion.json').read_text(encoding='utf-8'))
shell = json.loads((coarse / 'continuous_shell_audit.json').read_text(encoding='utf-8'))
geometry = json.loads((coarse / 'geometry_audit.json').read_text(encoding='utf-8'))
router = json.loads((coarse / 'visibility_event_routing.json').read_text(encoding='utf-8'))
with np.load(coarse / 'NOT_RESTARTABLE_candidate.npz', allow_pickle=True) as saved:
    parameters = saved['solver_data'].item()['parameters']
with np.load(coarse / 'visible_graph.npz') as saved:
    graph = {key: saved[key] for key in saved.files}
x0, x1 = graph['position0'], graph['position1']
i, j, ids = graph['i'], graph['j'], graph['particle_id']
o0, o1 = np.asarray(motion['origin0_m']), np.asarray(motion['origin1_m'])
R0, R1 = np.asarray(motion['rotation0']), np.asarray(motion['rotation1'])
angle, rotation = _rotation_path(R0, R1)
pairs = np.column_stack((ids[i], ids[j])).astype('<i8', copy=False)
input_arrays = dict(world_position0=x0, world_position1=x1,
    particle_id_order=ids.astype('<i8'), retained_edge_pid_order=pairs,
    origin0=o0, origin1=o1, rotation0=R0, rotation1=R1,
    geometry=np.array([parameters['tool_radius'], parameters['bore_width']/2,
        parameters['tool_top']-parameters['layer']]))
input_hashes = {key: _hash_array(value) for key, value in input_arrays.items()}
input_hashes['retained_topology'] = _topology_hash(pairs)
assert input_hashes == original['input_sha256']
modules = [
    'scripts/visible_support_remote_domain_20261007.py',
    'scripts/visible_support_path_tight.py',
    'scripts/visible_support.py',
    'scripts/visible_support_numerics.py',
    'scripts/audit_solid_occluded_neighbors.py',
    'src/afsd_backend/pysph/quasistatic_contact.py',
]
source_hashes = {name: sha(ROOT/name) for name in modules}
assert all(manifest['source_sha256'][name] == digest for name, digest in source_hashes.items())
evidence_files = [source/'manifest.json', source/'protocol.json',
    coarse/'visible_path_audit.json', coarse/'motion.json', coarse/'visible_graph.npz',
    coarse/'NOT_RESTARTABLE_candidate.npz', coarse/'continuous_shell_audit.json',
    coarse/'geometry_audit.json', coarse/'visibility_event_routing.json']
reports = []
for number, (a, b) in enumerate(((0., .5), (.5, 1.)), start=1):
    xa = x0 if a == 0. else x0+a*(x1-x0)
    xb = x1 if b == 1. else x0+b*(x1-x0)
    qa = dict(origin=o0 if a == 0. else o0+a*(o1-o0),
        rotation=R0 if a == 0. else rotation(a))
    qb = dict(origin=o1 if b == 1. else o0+b*(o1-o0),
        rotation=R1 if b == 1. else rotation(b))
    began = time.perf_counter()
    report = audit_visible_support_path(xa, xb, i, j, qa, qb, parameters, ids)
    elapsed = time.perf_counter()-began
    path = target/f'half_{number}_visible_path_audit.json'
    save(path, report)
    reports.append(dict(saved_chord_start_fraction=a, saved_chord_end_fraction=b,
        wall_seconds=elapsed, report=rel(path), report_sha256=sha(path),
        status=report['status'], failure_kind=report['failure_kind'],
        analytic_batch_queries=report['analytic_batch_queries'],
        maximum_queries=report['maximum_queries'], deepest_subdivision=report['deepest_subdivision'],
        input_sha256=report['input_sha256'],
        found_edge_count=len(report['found_edges']),
        degenerate_edge_count=len(report['degenerate_edges']),
        numerical_domain_failure_count=len(report['numerical_domain_failures']),
        unresolved_interval_count=len(report['unresolved_intervals'])))
summary = dict(
    classification='SAVED_REJECTED_CHORD_READ_ONLY_GEOMETRY_BUDGET_DIAGNOSTIC',
    status='BOTH_DIAGNOSTIC_HALF_CHORD_CERTIFICATES_PASS' if all(
        x['status']=='CERTIFIED_VISIBLE_SUPPORT_PATH' for x in reports) else 'DIAGNOSTIC_HALF_CHORD_NOT_CERTIFIED',
    accepted=False, restart_allowed=False, diagnostic_only=True,
    physical_solver_evaluations=0, nnps_evaluations=0, material_updates=0,
    force_evaluations=0, gpu_evaluations=0, accepted_time_advanced_s=0.,
    original_parent_checkpoint=manifest['parent_checkpoint'], original_parent_sha256=manifest['parent_sha256'],
    original_time_interval_s=[motion['t0_s'],motion['t1_s']], original_dt_s=manifest['dt_s'],
    original_status=original['status'], original_failure_kind=original['failure_kind'],
    original_queries=original['analytic_batch_queries'], original_unresolved_intervals=original['unresolved_intervals'],
    original_found_edges=original['found_edges'], original_degenerate_edges=original['degenerate_edges'],
    original_numerical_domain_failures=original['numerical_domain_failures'],
    original_shell_status=shell['status'], original_geometry_audit=geometry,
    original_router_status=router['status'], original_router_requires_resolve=router['requires_resolve'],
    original_router_changed_edges=router['changed_edges'],
    original_router_geometry_evaluations=router['geometry_evaluations'],
    original_input_sha256=input_hashes, exact_saved_input_hash_match=True,
    original_evidence_sha256={rel(path):sha(path) for path in evidence_files},
    diagnostic_source_sha256={rel(Path(__file__).resolve()):sha(Path(__file__).resolve()), **source_hashes},
    source_matches_original_physical_manifest=True, halves=reports,
    interpretation='The original whole-chord certificate exhausted its per-call query budget; both diagnostic half chords satisfy the unchanged certificate. This supports trying a bounded smaller physical step but does not accept the rejected candidate or certify a newly solved step.',
    required_next_action='From the last complete qualified parent, physically solve a new coarse and two-half pair at dt/2, retaining all original candidate and pair gates. Preserve the old failure. Stop on observed event, real intersection, numerical ambiguity, other failed gates, or bounded retry exhaustion.',
    limitations=[
        'Intermediate diagnostic positions are interpolated saved chord values, never accepted states or checkpoints.',
        'Separate half-path SO(3) evaluations use the inherited rotation routine and its original engineering roundoff envelope; they are not a claim of bitwise equality of every full-path fractional rotation.',
        'No guarantee that the newly solved smaller physical step will pass; its actual trajectory and fields must be recomputed.',
        'No support edge is removed, no physical tolerance is relaxed, and no topology event integrator is implemented.',
    ])
destination = base/'geometry_budget_review_20261010_v1.json'
save(destination, summary)
print(json.dumps(dict(summary=rel(destination),status=summary['status'],halves=[
    {k:value[k] for k in ['status','analytic_batch_queries','wall_seconds','report']} for value in reports]), indent=2))

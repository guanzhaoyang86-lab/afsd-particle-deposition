"""Independent Q7R provenance/config gate, then frozen Q7 saved-data audit.

No solver or material update is imported. Configuration is compared with the
hash-pinned original Q7 output, permitting only the declared 12 -> 14 budget,
the explicit recovery metadata, and the preregistered 4/8 requested grid.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import inspect_process_plastic_friction as parent

ROOT = Path(__file__).resolve().parents[1]
BASELINE = 'outputs/process_plastic_friction_q7_20261001/n3_i4/resolved_config.json'
BASELINE_SHA256 = '651f7173006cd7123d8ccdf5b43d1bde81796841814b22666c468b0b2495a847'
RECOVERY_CONFIG = 'configs/benchmarks/process_plastic_friction_q7_recovery.json'
SUPPLEMENTAL_SOURCES = {'scripts/probe_process_plastic_friction_recovery.py', RECOVERY_CONFIG}
BASELINE_SOURCES = {'scripts/probe_process_plastic_friction.py',
    'configs/benchmarks/process_plastic_friction_q7.json',
    'configs/benchmarks/process_contact_q3.json', 'configs/benchmarks/process_plastic_q2.json',
    'scripts/probe_process_friction_cone.py', 'configs/benchmarks/process_friction_cone_q6.json',
    'scripts/probe_process_friction.py', 'configs/benchmarks/process_friction_q5.json',
    'scripts/probe_process_contact.py', 'scripts/probe_process_plastic.py',
    'scripts/probe_process_time.py', 'scripts/probe_engineering_material.py',
    'scripts/probe_engineering_block.py', 'scripts/probe_boundary_sbp.py',
    'src/afsd_backend/pysph/material.py', 'src/afsd_backend/pysph/mechanics.py'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def inspect(directory):
    path = Path(directory)
    r = read_json(path/'report.json')
    cfg = read_json(path/'resolved_config.json')
    sidecar_path = path/'recovery_provenance.json'
    sidecar = read_json(sidecar_path)
    recovery = read_json(ROOT/RECOVERY_CONFIG)
    budget = recovery.get('maximum_iterations', {})
    if (recovery.get('variant') != 'Q7R_BOUNDED_NEWTON14'
            or recovery.get('production_authorized') is not False
            or type(budget.get('value')) is not int or budget['value'] != 14
            or budget.get('unit') != 'count'
            or recovery.get('baseline_maximum_iterations') != 12
            or recovery.get('run_timeout_s') != 60
            or recovery.get('permitted_cases') != [4, 8]
            or recovery.get('parent') != 'configs/benchmarks/process_plastic_friction_q7.json'):
        raise ValueError('undeclared recovery variant/budget/scope')
    if (sidecar.get('variant') != recovery['variant']
            or sidecar.get('production_authorized') is not False
            or sidecar.get('configuration') != recovery):
        raise ValueError('recovery provenance/configuration mismatch')
    if sidecar.get('source_sha256') != {name:digest(ROOT/name) for name in SUPPLEMENTAL_SOURCES}:
        raise ValueError('supplemental source inventory/hash mismatch')
    sources = sidecar.get('baseline_source_sha256')
    if (not isinstance(sources, dict) or set(sources) != BASELINE_SOURCES
            or sources != r.get('source_sha256')
            or sources != {name:digest(ROOT/name) for name in BASELINE_SOURCES}):
        raise ValueError('baseline source inventory/hash mismatch')
    if r.get('artifact_sha256', {}).get('recovery_provenance.json') != digest(sidecar_path):
        raise ValueError('recovery provenance missing from artifact inventory or changed')
    if digest(ROOT/BASELINE) != BASELINE_SHA256:
        raise ValueError('frozen original resolved configuration changed')
    expected = read_json(ROOT/BASELINE)
    n = r.get('increments_per_phase')
    if type(n) is not int or n not in (4, 8):
        raise ValueError('undeclared requested increment grid')
    if expected['plastic_parent']['numerical']['maximum_iterations'] != 12:
        raise ValueError('original baseline budget is not 12')
    expected = copy.deepcopy(expected)
    expected['increments_per_phase'] = n
    expected['plastic_parent']['numerical']['maximum_iterations'] = 14
    expected['recovery_variant'] = recovery
    if cfg != expected:
        raise ValueError('resolved configuration changed outside declared budget/grid/metadata')
    result = parent.inspect(path)
    result['recovery_audit'] = dict(status='PASS_RECOVERY_METADATA',
        variant=recovery['variant'], baseline_maximum_iterations=12,
        maximum_iterations=14, baseline_resolved_config_sha256=BASELINE_SHA256,
        supplemental_source_sha256=sidecar['source_sha256'],
        physics_tolerances_boundaries_unchanged=True,
        caveat='Metadata and saved-data verification only; not independent continuum validation.')
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inspect(args.run)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result['status'] == 'PASS_COMPLETE_READBACK' else 2


if __name__ == '__main__':
    raise SystemExit(main())

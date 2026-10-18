"""Read-only qualification of the saved failed chord; no physical solve."""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import visible_support_anchored_20261018 as proof


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    began=time.perf_counter()
    saved=ROOT/'outputs/short_bar_event_continue_20261015_v4/run/piece_000002_attempt_000001/fine_1/after_event/proposal'
    names=('visible_graph.npz','motion.json','visible_path_audit.json','NOT_RESTARTABLE_candidate.npz')
    before={str((saved/name).relative_to(ROOT)).replace('\\','/'):sha(saved/name) for name in names}
    original=json.loads((saved/'visible_path_audit.json').read_text())
    assert original['status']==proof.UNRESOLVED and original['analytic_batch_queries']==256
    assert not original['found_edges'] and not original['degenerate_edges'] and not original['numerical_domain_failures']
    cfg=json.loads((ROOT/'configs/deposition_fill_then_track_20261007.json').read_text())
    p=json.loads((ROOT/cfg['base_parameters']).read_text())['parameters']
    p.update({k:v['value'] for k,v in cfg['parameter_overrides'].items()})
    motion=json.loads((saved/'motion.json').read_text())
    with np.load(saved/'visible_graph.npz',allow_pickle=False) as z:
        report=proof.audit_visible_support_path(z['position0'],z['position1'],z['i'],z['j'],
            dict(origin=motion['origin0_m'],rotation=motion['rotation0']),
            dict(origin=motion['origin1_m'],rotation=motion['rotation1']),p,z['particle_id'])
    assert report['input_sha256']==original['input_sha256']
    assert report['status']==proof.CERTIFIED and report['analytic_batch_queries']<=256
    assert all(sha(ROOT/path)==value for path,value in before.items())
    complete=BASE/'saved_geometry_replay_full.json'
    complete.write_text(json.dumps(report,indent=2)+'\n')
    supplementary=[w for row in report['certificates'] for w in row.get('witnesses',[])
        if w['endpoint_enclosure_method']=='ANGLE_AWARE_ENDPOINT_ANCHORED_CURVATURE_PRISM']
    assert supplementary
    source_paths=('scripts/visible_support_anchored_20261018.py',
        'scripts/visible_support_prism_20261010.py',
        'tests/review/test_visible_support_anchored_20261018.py',
        'docs/VISIBLE_SUPPORT_ANCHORED_CONTRACT_20261018.md',
        str(Path(__file__).resolve().relative_to(ROOT)).replace('\\','/'))
    result=dict(status='PASS_READONLY_SAVED_CHORD_ENDPOINT_ANCHORED_CERTIFICATE',
        physical_solves=0,nnps_evaluations=0,physical_state_promoted=False,restart_allowed=False,
        old_status=original['status'],old_query_count=original['analytic_batch_queries'],
        new_status=report['status'],new_query_count=report['analytic_batch_queries'],
        edge_count=report['edge_count'],new_deepest_subdivision=report['deepest_subdivision'],
        input_sha256=report['input_sha256'],old_input_sha256=original['input_sha256'],
        endpoint_anchored_enclosure=report['endpoint_anchored_enclosure'],
        supplemental_witnesses=len(supplementary),minimum_supplemental_projection_gap_m=min(
            wall['strict_gap_m'] for witness in supplementary for wall in witness['wall_separations']),
        saved_artifact_sha256=before,full_report_sha256=sha(complete),
        source_sha256={path:sha(ROOT/path) for path in source_paths},
        wall_seconds=time.perf_counter()-began)
    (BASE/'saved_geometry_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('status','edge_count','new_query_count',
        'new_deepest_subdivision','supplemental_witnesses','minimum_supplemental_projection_gap_m','wall_seconds')}))


if __name__=='__main__':main()

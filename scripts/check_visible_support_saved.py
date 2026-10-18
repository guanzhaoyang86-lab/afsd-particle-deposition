"""Read-only saved-graph filter audit; no NNPS or physical solve."""
from pathlib import Path
import json
import time
import numpy as np
import readback_discrete_process as reader
from afsd_backend.pysph.kernel import HyperbolicSpline
from visible_support import filter_visible_support

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/coupled_cost_20261004/visible_support_geometry_saved.json'


def main():
    reader.require(not OUT.exists(),'refuse saved-graph audit overwrite')
    cases=[('accepted_39_375ms',ROOT/'outputs/coupled_cost_20261004/batched/pre_occlusion_response/fine/step_032/accepted.npz'),
        ('old_rejected_39_875ms',ROOT/'outputs/coupled_cost_20261004/batched/post_passage_short_response/coarse/step_020/NOT_RESTARTABLE_candidate.npz')]
    rows=[];inputs={}
    for name,path in cases:
        raw=reader.load(str(path));pa,state=raw['arrays']['solid'],raw['solver_data']
        evidence_path=path.parent/'candidate_evidence.npz'
        with np.load(evidence_path,allow_pickle=False) as saved:
            reader.require(np.array_equal(pa.particle_id,saved['particle_id']),'saved graph particle order differs')
            i,j=saved['i'].copy(),saved['j'].copy()
        x=np.column_stack((pa.x,pa.y,pa.z))
        # A pure radial kernel geometry calculation on already-saved pairs,
        # exactly the original operators formula; no neighborhood is rebuilt.
        g=HyperbolicSpline().gradient(x[i]-x[j],(pa.h[i]+pa.h[j])/2)
        p=state['parameters'];pose=reader.expected_pose(p,state['time_s'])
        started=time.perf_counter();result=filter_visible_support(x,i,j,g,pose,p,pa.particle_id)
        elapsed=time.perf_counter()-started
        if name=='accepted_39_375ms':
            reader.require(not len(result['blocked_raw_edge_indices']) and
                result['i'].tobytes()==i.tobytes() and result['j'].tobytes()==j.tobytes() and result['g'].tobytes()==g.tobytes(),
                'accepted saved graph is not an exact no-filter identity')
        else:
            reader.require([sorted(pair) for pair in result['blocked_edge_pid_pairs'].tolist()]==[[219,251]],
                'old failed edge was not independently excluded')
        rows.append(dict(name=name,time_s=state['time_s'],global_step=state['step'],
            originally_accepted=path.name=='accepted.npz',original_restart_allowed=state.get('restart_allowed',True),
            graph_i_j_g_bitwise_identical=bool(np.array_equal(result['i'],i) and np.array_equal(result['j'],j) and np.array_equal(result['g'],g)),
            diagnostics=result['diagnostics'],fingerprints=result['fingerprints'],wall_seconds=elapsed))
        inputs.update({v.relative_to(ROOT).as_posix():reader.sha(v) for v in (path,evidence_path)})
    source_paths=[Path(__file__),ROOT/'scripts/visible_support.py',ROOT/'scripts/test_visible_support.py',
        ROOT/'scripts/audit_solid_occluded_neighbors.py',ROOT/'src/afsd_backend/pysph/kernel.py',Path(reader.__file__)]
    report=dict(status='PASS_SAVED_ENDPOINT_GRAPH_IDENTITY_AND_OLD_BLOCKED_EDGE_DIAGNOSTIC',
        physical_solver_evaluations=0,nnps_evaluations=0,material_updates=0,force_evaluations=0,
        original_failed_candidate_status_unchanged=True,rows=rows,input_sha256=inputs,
        source_sha256={p.relative_to(ROOT).as_posix():reader.sha(p) for p in source_paths},
        limitations=['This proves an endpoint graph selector identity, not a recomputed coupled-field or full-residual equivalence.',
            'Raw radial gradients were recomputed geometrically on the saved pairs; no NNPS or physics was run.',
            'No retained-edge path certificate is implied; the old rejected candidate remains unaccepted.'])
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(output=str(OUT),rows=[dict(name=r['name'],raw=r['diagnostics']['raw_edge_count'],
        visible=r['diagnostics']['visible_edge_count'],blocked=r['diagnostics']['blocked_edge_count'],
        strict=r['diagnostics']['strict_blocked_edge_count'],roundoff=r['diagnostics']['roundoff_only_blocked_edge_count'],
        identity=r['graph_i_j_g_bitwise_identical'],wall_seconds=r['wall_seconds']) for r in rows])))


if __name__=='__main__':
    main()

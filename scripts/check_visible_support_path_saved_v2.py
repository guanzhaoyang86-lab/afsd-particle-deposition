"""Pure saved-state path audit, with fixed endpoint-retained support edges."""
from pathlib import Path
import json
import time
import numpy as np
import readback_discrete_process as reader
from afsd_backend.pysph.kernel import HyperbolicSpline
from visible_support import filter_visible_support
from visible_support_path import audit_visible_support_path
from visible_support_numerics import require_visible_support_numerical_domain

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/coupled_cost_20261004/visible_support_last_accepted_path_v2.json'


def main():
    reader.require(not OUT.exists(),'refuse saved-path audit overwrite')
    base=ROOT/'outputs/coupled_cost_20261004/batched/pre_occlusion_response/fine'
    first=base/'step_031/accepted.npz';last=base/'step_032/accepted.npz'
    raw0=reader.load(str(first));raw1=reader.load(str(last))
    pa0,state0=raw0['arrays']['solid'],raw0['solver_data'];pa1,state1=raw1['arrays']['solid'],raw1['solver_data']
    reader.require(np.array_equal(pa0.particle_id,pa1.particle_id),'particle order changed')
    reader.require(state0['parameters']==state1['parameters'],'parameters changed')
    evidence=last.parent/'candidate_evidence.npz'
    with np.load(evidence,allow_pickle=False) as saved:
        reader.require(np.array_equal(pa1.particle_id,saved['particle_id']),'saved graph particle order differs')
        i,j=saved['i'].copy(),saved['j'].copy()
    x0=np.column_stack((pa0.x,pa0.y,pa0.z));x1=np.column_stack((pa1.x,pa1.y,pa1.z))
    p=state1['parameters'];pose0=reader.expected_pose(p,state0['time_s']);pose1=reader.expected_pose(p,state1['time_s'])
    g=HyperbolicSpline().gradient(x1[i]-x1[j],(pa1.h[i]+pa1.h[j])/2)
    domain=require_visible_support_numerical_domain(x1,i,j,pose1,p,pa1.particle_id)
    endpoint=filter_visible_support(x1,i,j,g,pose1,p,pa1.particle_id)
    reader.require(not len(endpoint['blocked_raw_edge_indices']),'accepted endpoint graph unexpectedly filtered')
    started=time.perf_counter()
    audit=audit_visible_support_path(x0,x1,endpoint['i'],endpoint['j'],pose0,pose1,p,pa1.particle_id)
    elapsed=time.perf_counter()-started
    source_paths=[Path(__file__),ROOT/'scripts/visible_support.py',ROOT/'scripts/visible_support_path.py',
        ROOT/'scripts/test_visible_support.py',ROOT/'scripts/test_visible_support_path.py',ROOT/'scripts/visible_support_numerics.py',ROOT/'scripts/test_visible_support_numerics.py',
        ROOT/'scripts/audit_solid_occluded_neighbors.py',ROOT/'src/afsd_backend/pysph/quasistatic_contact.py',
        ROOT/'src/afsd_backend/pysph/kernel.py',Path(reader.__file__)]
    report=dict(status=audit['status'],physical_solver_evaluations=0,nnps_evaluations=0,
        force_evaluations=0,material_updates=0,old_acceptance_unchanged=True,
        start_time_s=state0['time_s'],end_time_s=state1['time_s'],start_global_step=state0['step'],end_global_step=state1['step'],
        endpoint_numerical_domain=domain,endpoint=endpoint['diagnostics'],endpoint_fingerprints=endpoint['fingerprints'],audit=audit,wall_seconds=elapsed,
        input_sha256={path.relative_to(ROOT).as_posix():reader.sha(path) for path in (first,last,evidence)},
        source_sha256={path.relative_to(ROOT).as_posix():reader.sha(path) for path in source_paths},
        limitations=['Only this saved final segment and its endpoint-retained graph are covered.',
            'No existing acceptance or failed-run classification is changed by this read-only diagnostic.',
            'The declared path certificate does not establish inertia accuracy or unknown exact material trajectories.'])
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(dict(output=str(OUT),status=audit['status'],wall_seconds=elapsed,edge_count=audit['edge_count'],
        queries=audit['analytic_batch_queries'],depth=audit['deepest_subdivision'],found=len(audit['found_edges']),
        unresolved=len(audit['unresolved_intervals']),endpoint_blocked=endpoint['diagnostics']['blocked_edge_count'])))


if __name__=='__main__':main()

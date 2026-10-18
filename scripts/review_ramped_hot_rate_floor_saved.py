"""All-particle saved rate-floor and contact-neighborhood inventory; no physics."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import review_expanded_pid797_saved as pure

ROOT=pure.ROOT
BASE=ROOT/'outputs/ramped_hot_startup_20261005/v1'
METHOD_SHA='d27a4a84dfe5c7fd716cfc64d48020ce705afdbf92e9d2320378f80d4aedef54'


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True);args=parser.parse_args()
    out=(ROOT/args.output).resolve();assert out.is_relative_to(ROOT/'outputs') and not out.exists()
    t=time.perf_counter();saved=pure.Saved();saved.bind(Path(__file__));saved.bind(Path(pure.__file__))
    method_path=BASE/'method_review/summary.json';method=saved.read(method_path);assert pure.sha(method_path)==METHOD_SHA
    ids=sorted({method['maximum_stress_row']['particle_id'],method['maximum_plastic_row']['particle_id']})
    stages={}
    for stage in ('coarse','fine'):
        initial=pure.load(str(saved.bind(BASE/'initial/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz')))
        previous=initial['arrays']['solid'];p=initial['solver_data']['parameters'];rows=[]
        for chunk in range(1,5):
            folder=BASE/f'chunk_{chunk:03d}'/stage
            for row in saved.read(folder/'candidate_history.json'):
                path=folder/row['checkpoint'];data=pure.load(str(saved.bind(path)));pa,s=data['arrays']['solid'],data['solver_data']
                assert pure.sha(path)==row['sha256'] and not s['accepted'] and not s['restart_allowed']
                sd=path.parent;j=saved.read(sd/'FULL_UNCOMMITTED_candidate.json');c=j['contact']
                force=np.asarray(c['force_by_feature_N']);active=np.asarray(c['active'])
                loaded=pa.particle_id[(np.linalg.norm(force[:,:7],axis=2).max(axis=1)>j['metrics']['force_limit_N']) & active[:,:7].any(axis=1)]
                neighborhoods=[]
                with np.load(saved.bind(sd/'visible_graph.npz'),allow_pickle=False) as g:
                    for pid in ids:
                        near={}
                        for prefix in ('','raw_'):
                            a,b=g[prefix+'i'],g[prefix+'j'];mask=(pa.particle_id[a]==pid)|(pa.particle_id[b]==pid)
                            values=set(pa.particle_id[np.r_[a[mask],b[mask]]].astype(int).tolist());values.discard(pid)
                            near['visible' if prefix=='' else 'raw']=sorted(values)
                        nearby=sorted(set(near['visible']) & set(loaded.astype(int).tolist()))
                        neighborhoods.append(dict(particle_id=pid,one_hop=near,loaded_one_hop_ids=nearby,
                            loaded_neighbor_contacts={str(v):pure.contact_row(j,int(np.flatnonzero(pa.particle_id==v)[0])) for v in nearby}))
                crossings=pa.particle_id[(previous.strain_rate<=p['rate_floor'])!=(pa.strain_rate<=p['rate_floor'])].astype(int).tolist()
                rows.append(dict(time_s=s['time_s'],checkpoint=path.relative_to(ROOT).as_posix(),checkpoint_sha256=row['sha256'],
                    accepted=False,restart_allowed=False,particle_count=len(pa.x),rate_floor_per_s=p['rate_floor'],
                    minimum_total_equivalent_rate_per_s=float(pa.strain_rate.min()),maximum_total_equivalent_rate_per_s=float(pa.strain_rate.max()),
                    above_floor_ids=pa.particle_id[pa.strain_rate>p['rate_floor']].astype(int).tolist(),
                    floor_branch_changed_ids=crossings,loaded_nonbase_ids=loaded.astype(int).tolist(),neighborhoods=neighborhoods,
                    diagnostic_unused_after_failed_shared_endpoint=bool(s['time_s']>.0035+1e-15)))
                previous=pa
        stages[stage]=rows
    through=[r for rows in stages.values() for r in rows if not r['diagnostic_unused_after_failed_shared_endpoint']]
    no_crossing=all(not r['above_floor_ids'] and not r['floor_branch_changed_ids'] for r in through)
    for name,digest in saved.files.items():assert pure.sha(ROOT/name)==digest
    report=dict(status='SAVED_ONLY_RAMP_ALL_PARTICLE_RATE_FLOOR_AND_ONE_HOP_INVENTORY',
        physical_solver_evaluations=0,nnps_evaluations=0,contact_evaluations=0,constitutive_update_evaluations=0,
        accepted=False,restart_allowed=False,new_accepted_steps=0,accepted_time_advanced_s=0.,particle_ids=ids,stages=stages,
        all_918_points_below_rate_floor_through_failed_3_5ms=no_crossing,
        interpretation=[
            'Through the failed 3.5 ms comparison every saved particle rate on both trajectories is below 1/s; the original FKS rate multiplier remains one.',
            'There is no saved rate-floor branch crossing at this failure. The unused coarse 4 ms state must not be used to explain a prior 3.5 ms comparison.',
            'PID791 has no direct active contact. Actual loaded particles 796/797 enter its retained one-hop neighborhood by fine 3.25 ms and coarse 3.5 ms.',
            'The coarse interval has no extra PID791 plastic strain, while the fine first half plastifies and the second half unloads; this is evidence of different resolved material histories.',
            'These observations locate sensitivity in the coupled loading/plastic response; they do not prove a unique causal mechanism or a remedy.'
        ],input_and_source_sha256=saved.files,wall_seconds=time.perf_counter()-t)
    out.mkdir(parents=True,exist_ok=False);(out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('status','all_918_points_below_rate_floor_through_failed_3_5ms','wall_seconds')},indent=2));return 0


if __name__=='__main__':raise SystemExit(main())

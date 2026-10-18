"""Bounded saved-geometry diagnosis; never retroactively accept a failed path."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time
import numpy as np
from visible_support_path import audit_visible_support_path, _rotation_path
from visible_support import _hash_array
from audit_solid_occluded_neighbors import _partitions
from visible_support_numerics import check_partition_numerical_domain, PASS

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'outputs/visible_local_refinement_20261004/window_004/coarse/step_008'
OUTPUT=ROOT/'outputs/visible_path_diagnostic_20261004/quarter_turn_v1'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def relative(path):return path.relative_to(ROOT).as_posix()


def worker():
    started=time.perf_counter();protocol=read(OUTPUT/'protocol.json');bound=protocol['source_and_input_sha256']
    for name,value in bound.items():
        if sha(ROOT/name)!=value:raise ValueError('bound input changed')
    old=read(INPUT/'visible_path_audit.json');meta=read(INPUT/'visible_graph.json')
    if sha(INPUT/'visible_graph.npz')!=old['graph_sha256'] or sha(INPUT/'NOT_RESTARTABLE_candidate.npz')!=old['candidate_sha256']:
        raise ValueError('original failed-path input binding differs')
    with np.load(INPUT/'visible_graph.npz',allow_pickle=False) as saved:data={k:saved[k] for k in saved.files}
    with np.load(INPUT/'NOT_RESTARTABLE_candidate.npz',allow_pickle=True) as saved:state=saved['solver_data'].item()
    if state['restart_allowed'] is not False:raise ValueError('failed reference became restartable')
    p=state['parameters'];t0=meta['parent_time_s'];t1=meta['endpoint_time_s']
    def pose(t):
        a=p['omega']*t;c,s=np.cos(a),np.sin(a)
        return dict(origin=np.array([p['traverse']*t,0.,p['layer']]),rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]))
    first,last=data['position0'],data['position1'];po0,po1=pose(t0),pose(t1)
    for key,value in [('world_position0',first),('world_position1',last),('origin0',po0['origin']),('origin1',po1['origin']),('rotation0',po0['rotation']),('rotation1',po1['rotation'])]:
        if _hash_array(value)!=old['input_sha256'][key]:raise ValueError('original array hash differs '+key)
    indices=sorted({edge for interval in old['unresolved_intervals'] for edge in interval['edge_indices']})
    ids=data['particle_id'];i,j=data['i'][indices],data['j'][indices]
    if np.column_stack((ids[i],ids[j])).tolist()!=[[339,341]]:raise ValueError('registered unresolved edge changed')
    angle,rotation_at=_rotation_path(po0['rotation'],po1['rotation'])
    def subpose(f):return dict(origin=po0['origin']+f*(po1['origin']-po0['origin']),rotation=rotation_at(f))
    audits=[]
    for label,lo,hi in [('same_edge_full_step',0.,1.),('old_pending_interval_1',.875,.9375),('old_pending_interval_2',.9375,1.)]:
        if time.perf_counter()-started>17:raise TimeoutError('geometry diagnostic soft budget')
        a=first+lo*(last-first);b=first+hi*(last-first)
        result=audit_visible_support_path(a,b,i,j,subpose(lo),subpose(hi),p,ids,maximum_depth=32,maximum_queries=256)
        write(OUTPUT/(label+'.json'),result)
        audits.append(dict(label=label,original_fraction_interval=[lo,hi],status=result['status'],failure_kind=result['failure_kind'],
            analytic_batch_queries=result['analytic_batch_queries'],deepest_subdivision=result['deepest_subdivision'],
            found_edges=result['found_edges'],unresolved_intervals=result['unresolved_intervals'],
            certificate_file=relative(OUTPUT/(label+'.json')),certificate_sha256=sha(OUTPUT/(label+'.json'))))
    a=p['bore_width']/2;R=p['tool_radius'];H=p['tool_top']-p['layer'];samples=[]
    for f in (0.,.5,.875,.9375,1.):
        po=subpose(f);q=(first+f*(last-first)-po['origin'])@po['rotation'];start=q[i];direction=q[j]-start
        def hit(d):
            radius=np.nextafter(R+d,np.inf);domain=check_partition_numerical_domain(start,direction,radius,min(R,a,H))
            if domain['status']!=PASS:raise ValueError('static partition numerical domain unresolved')
            return bool(np.any(_partitions(start,direction,radius,np.nextafter(a-d,-np.inf),np.nextafter(-d,-np.inf),np.nextafter(H+d,np.inf))[3]>0.))
        strict=_partitions(start,direction,R,a,0.,H)[3]
        lower,upper=0.,1e-5
        if not hit(upper):raise ValueError('registered clearance bracket insufficient')
        for _ in range(48):
            mid=.5*(lower+upper)
            if hit(mid):upper=mid
            else:lower=mid
        samples.append(dict(fraction=f,time_s=t0+f*(t1-t0),body_endpoints_m=[start[0].tolist(),(start+direction)[0].tolist()],
            strict_inside_length_m=float(strict[0]),axis_aligned_expanded_solid_hit_threshold_bracket_m=[lower,upper],
            clearance_interpretation='Inflated-cylinder/shrunken-square/z-expanded partition threshold; NOT Euclidean distance or a continuous-path certificate.'))
    parent_graph=ROOT/'outputs/visible_local_refinement_20261004/window_004/coarse/step_007/visible_graph.npz'
    with np.load(parent_graph,allow_pickle=False) as saved:
        pairs=np.sort(np.column_stack((saved['particle_id'][saved['i']],saved['particle_id'][saved['j']])),axis=1)
        present=bool(np.any(np.all(pairs==[339,341],axis=1)))
    report=dict(status='COMPLETED_SAVED_GEOMETRY_DIAGNOSTIC_NOT_ORIGINAL_PATH_ACCEPTANCE',physical_solver_evaluations=0,nnps_evaluations=0,
        original_status=old['status'],original_failure_kind=old['failure_kind'],original_failure_preserved=True,original_reclassified_pass=False,
        source_and_input_sha256=bound,original_start_time_s=t0,original_endpoint_time_s=t1,
        original_found_edges=old['found_edges'],original_numerical_domain_failures=old['numerical_domain_failures'],
        original_maximum_queries=old['maximum_queries'],original_deepest_subdivision=old['deepest_subdivision'],
        original_maximum_depth=old['maximum_depth'],unresolved_edge_ids=[[339,341]],retained_in_original_parent_graph=present,
        independent_subset_audits=audits,static_geometry_samples=samples,wall_seconds=time.perf_counter()-started,
        conclusion_scope='Subset/subinterval diagnostics cannot retroactively change the old all-edge 256-query acceptance contract; no force/NNPS or extra particle path was generated.',
        quarter_turn_causality='NOT_ESTABLISHED_FROM_ANGLE_COINCIDENCE',exit_code=0)
    for name,value in bound.items():
        if sha(ROOT/name)!=value:raise ValueError('input changed during diagnosis')
    write(OUTPUT/'summary.json',report)
    print(json.dumps({k:report[k] for k in ('status','original_failure_kind','retained_in_original_parent_graph','independent_subset_audits','static_geometry_samples','wall_seconds')},indent=2))


if __name__=='__main__':
    if '--worker' in sys.argv:worker();raise SystemExit(0)
    if OUTPUT.exists():raise ValueError('refuse to overwrite independent diagnostic')
    OUTPUT.mkdir(parents=True)
    watched=[Path(__file__).resolve(),ROOT/'scripts/visible_support_path.py',ROOT/'scripts/visible_support.py',ROOT/'scripts/visible_support_numerics.py',
        ROOT/'scripts/audit_solid_occluded_neighbors.py',ROOT/'src/afsd_backend/pysph/quasistatic_contact.py',
        *(INPUT/name for name in ('visible_path_audit.json','visible_graph.json','visible_graph.npz','NOT_RESTARTABLE_candidate.npz')),
        INPUT.parent/'step_007/visible_graph.npz']
    protocol=dict(classification='INDEPENDENT_SAVED_GEOMETRY_SUBSET_DIAGNOSTIC_ONLY',registered_before_worker=True,
        maximum_external_wall_seconds=20,soft_geometry_budget_seconds=17,source_and_input_sha256={relative(p):sha(p) for p in watched},
        original_query_limit_unchanged=256,registered_calls=[dict(edge_ids=[339,341],interval=span,maximum_queries=256,maximum_depth=32)
            for span in ([0.,1.],[.875,.9375],[.9375,1.])],physical_solver_evaluations=0,nnps_evaluations=0,
        static_clearance_sample_fractions=[0.,.5,.875,.9375,1.],static_expansion_bisections_per_sample=48,
        no_old_acceptance_or_failures_changed=True)
    write(OUTPUT/'protocol.json',protocol)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--worker'];begin=time.perf_counter()
    with (OUTPUT/'run.log').open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=20,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
        except subprocess.TimeoutExpired:code=124
    write(OUTPUT/'supervision.json',dict(command=command,exit_code=code,external_wall_seconds=time.perf_counter()-begin,limit_seconds=20))
    print((OUTPUT/'run.log').read_text(encoding='utf-8'));raise SystemExit(code)

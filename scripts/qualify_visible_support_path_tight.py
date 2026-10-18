"""Supervised components and entire saved failed-path replay; no physics."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/visible_path_tight_20261004/qualification_v1'
INPUT=ROOT/'outputs/visible_local_refinement_20261004/window_004/coarse/step_008'
CONTRACT=ROOT/'docs/VISIBLE_SUPPORT_TIGHT_BOUND_CONTRACT_20261004.md'
LIMIT=120


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(p):return json.loads(p.read_text(encoding='utf-8'))


def verify_bound():
    protocol=read(OUT/'protocol.json')
    for name,digest in protocol['source_and_input_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('source/input hash changed: '+name)
    return protocol


def components():
    verify_bound()
    import test_visible_support_path_tight as module
    start=time.perf_counter()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(module))
    verify_bound()
    report=dict(status='PASS_COMPONENTS' if result.wasSuccessful() else 'FAIL_COMPONENTS',
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),
        wall_seconds=time.perf_counter()-start,physical_solver_evaluations=0,nnps_evaluations=0)
    write(OUT/'components.json',report)
    return 0 if result.wasSuccessful() else 1


def replay():
    verify_bound()
    from visible_support_path_tight import audit_visible_support_path,CERTIFIED
    from visible_support import _hash_array,_topology_hash
    old=read(INPUT/'visible_path_audit.json');meta=read(INPUT/'visible_graph.json')
    if old['status']!='UNRESOLVED_VISIBLE_SUPPORT_PATH' or old['failure_kind']!='MAXIMUM_ANALYTIC_BATCH_QUERIES':
        raise ValueError('registered failure differs')
    if sha(INPUT/'visible_graph.npz')!=old['graph_sha256'] or sha(INPUT/'NOT_RESTARTABLE_candidate.npz')!=old['candidate_sha256']:
        raise ValueError('original failed-path binding mismatch')
    with np.load(INPUT/'visible_graph.npz',allow_pickle=False) as f:data={k:f[k] for k in f.files}
    with np.load(INPUT/'NOT_RESTARTABLE_candidate.npz',allow_pickle=True) as f:state=f['solver_data'].item()
    if state['restart_allowed'] is not False:raise ValueError('original failed candidate restart flag changed')
    p=state['parameters'];t0=meta['parent_time_s'];t1=meta['endpoint_time_s']
    def pose(t):
        angle=p['omega']*t;c,s=np.cos(angle),np.sin(angle)
        return dict(origin=np.array([p['traverse']*t,0.,p['layer']]),rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]))
    first,last=data['position0'],data['position1'];a,b=pose(t0),pose(t1)
    ids=data['particle_id'];i,j=data['i'],data['j'];pairs=np.column_stack((ids[i],ids[j])).astype('<i8',copy=False)
    if len(i)!=7789:raise ValueError('entire registered graph must contain 7789 retained edges')
    arrays=dict(world_position0=first,world_position1=last,origin0=a['origin'],origin1=b['origin'],
        rotation0=a['rotation'],rotation1=b['rotation'],particle_id_order=ids.astype('<i8'),retained_edge_pid_order=pairs,
        geometry=np.array([p['tool_radius'],p['bore_width']/2,p['tool_top']-p['layer']]))
    for key,value in arrays.items():
        if _hash_array(value)!=old['input_sha256'][key]:raise ValueError('original array hash mismatch: '+key)
    if _topology_hash(pairs)!=old['input_sha256']['retained_topology']:raise ValueError('old topology differs')
    start=time.perf_counter()
    result=audit_visible_support_path(first,last,i,j,a,b,p,ids,maximum_depth=32,maximum_queries=256)
    result.update(diagnostic_only=True,original_failure_preserved=True,original_candidate_restart_allowed=False,
        original_result=rel(INPUT/'visible_path_audit.json'),original_result_sha256=sha(INPUT/'visible_path_audit.json'),
        graph_sha256=sha(INPUT/'visible_graph.npz'),candidate_sha256=sha(INPUT/'NOT_RESTARTABLE_candidate.npz'),
        parent_time_s=t0,endpoint_time_s=t1,accepted_time_added_s=0.,wall_seconds=time.perf_counter()-start)
    verify_bound();write(OUT/'replay.json',result)
    print(json.dumps({k:result[k] for k in ('status','edge_count','analytic_batch_queries','deepest_subdivision','failure_kind','wall_seconds')},indent=2))
    return 0 if result['status']==CERTIFIED else 2


def main():
    if '--components-worker' in sys.argv:return components()
    if '--replay-worker' in sys.argv:return replay()
    if OUT.exists():raise ValueError('refuse to overwrite a qualification')
    OUT.mkdir(parents=True)
    source_names=['scripts/qualify_visible_support_path_tight.py','scripts/visible_support_path_tight.py',
        'scripts/test_visible_support_path_tight.py','scripts/test_visible_support_path.py','scripts/test_visible_support_numerics.py',
        'scripts/visible_support_path.py','scripts/visible_support.py','scripts/visible_support_numerics.py',
        'scripts/audit_solid_occluded_neighbors.py','src/afsd_backend/pysph/quasistatic_contact.py',rel(CONTRACT)]
    watched=[ROOT/name for name in source_names]+[INPUT/name for name in
        ('visible_path_audit.json','visible_graph.json','visible_graph.npz','NOT_RESTARTABLE_candidate.npz')]
    protocol=dict(classification='INDEPENDENT_TIGHT_PATH_COMPONENTS_AND_ENTIRE_FAILED_GEOMETRY_REPLAY',
        registered_before_workers=True,contract=rel(CONTRACT),external_limit_seconds_per_worker=LIMIT,
        maximum_queries=256,maximum_depth=32,expected_replay_edge_count=7789,
        source_and_input_sha256={rel(p):sha(p) for p in watched},no_force_nnps_or_real_process_evaluations=True,
        no_retroactive_acceptance_or_restart=True,worker_order=['components','replay_only_if_components_pass'])
    write(OUT/'protocol.json',protocol);receipts=[];passed=True
    for worker in ('components','replay'):
        if not passed:break
        verify_bound();command=[sys.executable,'-B',str(Path(__file__).resolve()),'--'+worker+'-worker']
        begin=time.perf_counter()
        with (OUT/(worker+'.log')).open('w',encoding='utf-8') as stream:
            try:code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=LIMIT,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')).returncode
            except subprocess.TimeoutExpired:code=124
        receipt=dict(worker=worker,command=command,exit_code=code,external_wall_seconds=time.perf_counter()-begin,
            external_limit_seconds=LIMIT,log=rel(OUT/(worker+'.log')),log_sha256=sha(OUT/(worker+'.log')))
        receipts.append(receipt);write(OUT/(worker+'_supervision.json'),receipt)
        verify_bound();passed=code==0
    results={name:read(OUT/(name+'.json')) if (OUT/(name+'.json')).exists() else None for name in ('components','replay')}
    passed=passed and results['components'] is not None and results['components']['status']=='PASS_COMPONENTS' and \
        results['replay'] is not None and results['replay']['status']=='CERTIFIED_VISIBLE_SUPPORT_PATH'
    summary=dict(status='PASS_VISIBLE_PATH_TIGHT_COMPONENTS_AND_SAVED_REPLAY_NOT_PHYSICAL_VALIDATION' if passed else
        'FAIL_VISIBLE_PATH_TIGHT_QUALIFICATION',qualification_passed=passed,contract=rel(CONTRACT),
        protocol_sha256=sha(OUT/'protocol.json'),source_and_input_sha256=protocol['source_and_input_sha256'],
        components=results['components'],replay_status=None if results['replay'] is None else results['replay']['status'],
        replay_file=rel(OUT/'replay.json') if results['replay'] is not None else None,
        replay_sha256=sha(OUT/'replay.json') if results['replay'] is not None else None,workers=receipts,
        physical_solver_evaluations=0,nnps_evaluations=0,accepted_time_added_s=0.,original_failure_preserved=True,
        original_candidate_restart_allowed=False,before_after_source_input_hashes_match=True,
        limitation='Independent geometry qualification and saved-path diagnostic; no material/thermal/physical or timestep convergence validation.')
    write(OUT/'summary.json',summary);print(json.dumps(summary,indent=2));return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())

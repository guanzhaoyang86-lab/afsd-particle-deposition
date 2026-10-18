"""Read saved visible-support branches only: no NNPS, residual or force solve.

Arithmetic energy, observer, fingerprints and time screens are independently
replayed from arrays. Qualified coarse/cost are never appended as accepted.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import traceback
import numpy as np
from pysph.solver.utils import load

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.interval_acceptance import (
    validate_interval_energy_state, validate_bore_crossing_state,
    observe_square_bore_crossings)

BASE = ROOT/'outputs/visible_support_20261004/process_v1'
PARENT = ROOT/'outputs/coupled_cost_20261004/batched/pre_occlusion_response/fine/step_032/accepted.npz'
PARENT_SHA = '553f57b1dd5f0555d5a3726f7942705bd485ffd688ef3bf34f16ff342fd9a804'
PRIOR = ROOT/'outputs/coupled_cost_20261004/batched/readback_pre_occlusion_39_375ms/summary.json'
PRIOR_SHA = '4632c3cfc7430ad00942aeea58edb227607333c5b078610ca66af1838595b600'
CLASS = 'HOT_QS_VISIBLE_SUPPORT_CHAIN_RULE_CONTACT_APPROXIMATE_V1'
START, END, OLD_STEPS = .039375, .040375, 431
PASS = dict(cost='PASS_COST_ROOT_NOT_ACCEPTED', coarse='QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED',
            fine='PASS_VISIBLE_SUPPORT_WINDOW_TIME_SCREEN_NOT_PHYSICAL_VALIDATION')
INVARIANTS = ('particle_id','component_id','material_id','tracer_type','wholebar_tag','actuator','m','h','x0','y0','z0')


def require(value, message):
    if not bool(value): raise ValueError(message)


def clean(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k:clean(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [clean(v) for v in value]
    return value


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def array_hash(value):
    a=np.ascontiguousarray(value);h=hashlib.sha256()
    h.update(str(a.dtype).encode('ascii'));h.update(str(a.shape).encode('ascii'));h.update(a.tobytes())
    return h.hexdigest()


def ordered_pairs(ids,i,j):
    return np.column_stack((ids[i],ids[j])).astype('<i8')


def sorted_pairs(a):
    a=np.sort(np.asarray(a,dtype='<i8').reshape(-1,2),axis=1)
    return a[np.lexsort((a[:,1],a[:,0]))] if len(a) else a


def topology(a): return array_hash(sorted_pairs(a))
def positions(pa): return np.column_stack((pa.x,pa.y,pa.z))


def pose(p,t):
    a=p['omega']*t;c,s=np.cos(a),np.sin(a)
    return dict(origin=np.array([p['traverse']*t,0.,p['layer']]),
                rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]))


class Reader:
    def __init__(self): self.evidence={};self.cache={};self.source_entries=0
    def bind(self,path,expected=None):
        path=Path(path).resolve();require(path.is_relative_to(ROOT),'evidence escaped project')
        key=path.relative_to(ROOT).as_posix()
        if key not in self.evidence: self.evidence[key]=digest(path)
        if expected is not None: require(self.evidence[key]==expected,'SHA mismatch: '+key)
        return path
    def read(self,path):
        path=self.bind(path)
        return json.loads(path.read_text(encoding='utf-8'))
    def sources(self,mapping):
        for path,value in mapping.items(): self.bind(ROOT/path,value)
        self.source_entries+=len(mapping)
    def checkpoint(self,path):
        path=self.bind(path);key=str(path)
        if key not in self.cache:
            saved=load(key);self.cache[key]=(saved['arrays']['solid'],saved['solver_data'])
        return self.cache[key]
    def npz(self,path):
        with np.load(self.bind(path),allow_pickle=False) as value:
            return {key:value[key].copy() for key in value.files}
    def finish(self):
        for path,expected in self.evidence.items():require(digest(ROOT/path)==expected,'evidence changed during readback: '+path)


def properties_equal(a,b):
    require(set(a.properties)==set(b.properties),'particle property names changed')
    for key in a.properties:
        require(np.array_equal(a.properties[key].get_npy_array(),b.properties[key].get_npy_array()),'particle field differs: '+key)


def state_stats(pa,state):
    q=np.sqrt(1.5*np.sum(pa.S.reshape(-1,3,3)**2,axis=(1,2)))
    return dict(time_s=state['time_s'],global_step=state['step'],particle_count=len(pa.x),
        mass_kg=math.fsum(pa.m),temperature_min_K=float(pa.temperature.min()),temperature_max_K=float(pa.temperature.max()),
        plastic_strain_max=float(pa.ep.max()),equivalent_stress_max_Pa=float(q.max()),
        net_displacement_from_independent_initial_max_m=float(np.linalg.norm(positions(pa)-np.column_stack((pa.x0,pa.y0,pa.z0)),axis=1).max()),
        crossing_ids=state['bore_observation']['crossing_ids'],crossing_mass_kg=state['bore_observation']['crossing_mass_kg'])


def validate_state(pa,state,initial,prefix,p,number,dt,restart):
    require(state['classification']==CLASS and state['restart_allowed'] is restart,'model/restart flag mismatch')
    require(state['step']==OLD_STEPS+number and abs(state['time_s']-(START+number*dt))<1e-14,'candidate time/step mismatch')
    tr=state['spatial_model_transition']
    require(tr['parent_checkpoint_sha256']==PARENT_SHA and tr['diagnostic_only'] is False and
        tr['fields_reset'] is False and tr['interval_energy_reset'] is False and tr['observer_reset'] is False and
        tr['qualification']=='PASS_PARENT_MAP_TANGENT_AND_COMPONENT_GATES','invalid explicit transition')
    require(state['parameters']==p,'candidate parameter mismatch')
    for key in INVARIANTS:require(np.array_equal(getattr(pa,key),getattr(initial,key)),'passive identity/mass changed: '+key)
    for key,prop in pa.properties.items():require(np.isfinite(prop.get_npy_array()).all(),'nonfinite particle field: '+key)
    require(np.all(pa.rho>0) and np.all(pa.ep>=0) and np.all(pa.temperature>=p['room']) and np.all(pa.temperature<p['melt']),'stored physical-domain failure')
    validate_interval_energy_state(state['interval_energy']);validate_bore_crossing_state(state['bore_observation'])
    records=state['interval_energy']['records']
    require(len(records)==state['step'] and records[:OLD_STEPS]==prefix,'old complete interval prefix lost')
    require(state['crossing_ids']==state['bore_observation']['crossing_ids'],'observer crossing IDs differ')


def graph_audit(reader,directory,pa,state,previous,p,row):
    saved=reader.npz(directory/'visible_graph.npz');meta=reader.read(directory/'visible_graph.json')
    ids=pa.particle_id.astype('<i8');i,j,g=saved['i'],saved['j'],saved['g'];ri,rj,rg=saved['raw_i'],saved['raw_j'],saved['raw_g']
    kept,blocked=saved['kept_raw_edge_indices'],saved['blocked_raw_edge_indices']
    require(np.array_equal(saved['particle_id'],pa.particle_id),'graph particle ordering differs')
    require(np.array_equal(np.sort(np.r_[kept,blocked]),np.arange(len(ri))) and len(np.unique(np.r_[kept,blocked]))==len(ri),'raw graph partition invalid')
    for key,value in (('i',ri),('j',rj),('g',rg)):require(np.array_equal(saved[key],value[kept]),'graph is not common raw subset: '+key)
    require(np.all(ri>=0) and np.all(ri<rj) and np.all(rj<len(ids)) and len(np.unique(np.column_stack((ri,rj)),axis=0))==len(ri),'invalid raw undirected graph')
    require(np.array_equal(saved['position0'],positions(previous)) and np.array_equal(saved['position1'],positions(pa)),'graph positions differ from actual path')
    require(np.array_equal(saved['volume'],pa.m/pa.rho),'graph current volumes differ')
    ev=reader.npz(directory/'candidate_evidence.npz')
    require(all(np.array_equal(ev[key],saved[key]) for key in ('particle_id','i','j')),'journal active graph differs')
    rawpairs=ordered_pairs(ids,ri,rj);vpairs=ordered_pairs(ids,i,j)
    h=dict(particle_id_order_sha256=array_hash(ids),raw_edge_pid_order_sha256=array_hash(rawpairs),
        raw_topology_sha256=topology(rawpairs),visible_topology_sha256=topology(vpairs),
        blocked_topology_sha256=topology(rawpairs[blocked]),raw_gradient_sha256=array_hash(rg),
        visible_gradient_sha256=array_hash(g),world_positions_sha256=array_hash(saved['position1']))
    o1=pose(p,state['time_s']);t0=state['time_s']-row['dt_s'];o0=pose(p,t0)
    # Producer formed time endpoint by actual parent + dt; use the saved parent
    # time to retain identical floating-point pose evaluation.
    t0=meta['parent_time_s'];o0=pose(p,t0);o1=pose(p,t0+row['dt_s'])
    h.update(tool_origin_sha256=array_hash(o1['origin']),tool_rotation_sha256=array_hash(o1['rotation']))
    for key,value in h.items():require(meta['fingerprints'][key]==value,'graph fingerprint mismatch: '+key)
    fp=hashlib.sha256(sorted_pairs(vpairs).tobytes()).hexdigest()
    require(fp==row['graph_fingerprint']==meta['graph_fingerprint'],'graph topology fingerprint mismatch')
    require(meta['numerical_domain']['status']=='PASS_VISIBLE_SUPPORT_NUMERICAL_DOMAIN','numerical geometry domain uncertified')
    diag=meta['diagnostics'];blockedids=sorted(set(rawpairs[blocked].ravel().tolist()))
    require(diag['blocked_edge_count']==len(blocked)==row['blocked_edge_count'] and row['blocked_particle_ids']==blockedids,'blocked graph statistics differ')
    require(diag['visible_edge_count']==len(i) and diag['raw_edge_count']==len(ri),'graph counts differ')
    require(sorted(e['raw_edge_index'] for e in diag['blocked_edges'])==sorted(blocked.tolist()),'blocked edge records differ')
    path=reader.read(directory/'visible_path_audit.json')
    require(path['status']=='CERTIFIED_VISIBLE_SUPPORT_PATH' and not path['found_edges'] and not path['unresolved_intervals'] and
        not path['degenerate_edges'] and not path['numerical_domain_failures'] and path['edge_count']==len(i),'whole path not certified')
    require(path['candidate_sha256']==row['sha256'] and path['graph_sha256']==digest(directory/'visible_graph.npz'),'path certificate candidate/graph binding differs')
    ph=dict(world_position0=array_hash(saved['position0']),world_position1=array_hash(saved['position1']),
        particle_id_order=array_hash(ids),retained_edge_pid_order=array_hash(vpairs),retained_topology=topology(vpairs),
        origin0=array_hash(o0['origin']),origin1=array_hash(o1['origin']),rotation0=array_hash(o0['rotation']),rotation1=array_hash(o1['rotation']),
        geometry=array_hash(np.array([p['tool_radius'],.5*p['bore_width'],p['tool_top']-p['layer']])))
    require(path['input_sha256']==ph,'path certificate actual inputs differ')
    intervals=[[] for _ in range(len(i))]
    for cert in path['certificates']:
        lo,hi=cert['start_fraction'],cert['end_fraction'];require(0<=lo<hi<=1,'invalid certified interval')
        for k in cert['edge_indices']:
            require(0<=k<len(i),'certificate edge index outside graph');intervals[k].append((lo,hi))
    for edges in intervals:
        covered=0.
        for lo,hi in sorted(edges):require(lo<=covered+1e-15,'path certificate coverage gap');covered=max(covered,hi)
        require(covered>=1.-1e-15,'path certificate lacks complete [0,1] coverage')
    return saved,dict(raw_edge_count=len(ri),visible_edge_count=len(i),blocked_edge_count=len(blocked),blocked_particle_ids=blockedids,
        maximum_strict_intersection_length_m=diag['maximum_strict_intersection_length_m'],path_queries=path['analytic_batch_queries'],
        path_seconds=path['wall_seconds'],graph_fingerprint=fp)


def candidate_audit(reader,stage,directory,row,previous,previous_state,initial,prefix,p,frozen):
    path=reader.bind(directory/row['checkpoint'],row['sha256']);pa,state=reader.checkpoint(path)
    for relative,expected in row['artifacts'].items():reader.bind(directory/relative,expected)
    dt=.0000625 if stage=='fine' else .000125
    require(row['dt_s']==dt and row['global_step']==OLD_STEPS+row['step'] and row['time_s']==state['time_s'],'candidate history state mismatch')
    validate_state(pa,state,initial,prefix,p,row['step'],dt,False)
    require(state['implementation_source_sha256']==frozen,'candidate implementation provenance differs')
    require(state['interval_energy']['records'][:-1]==previous_state['interval_energy']['records'],'candidate interval direct parent differs')
    sd=path.parent;saved,graph=graph_audit(reader,sd,pa,state,previous,p,row)
    journal=reader.npz(sd/'FULL_UNCOMMITTED_candidate.npz');j=reader.read(sd/'FULL_UNCOMMITTED_candidate.json');energy=reader.read(sd/'energy.json')
    require(journal['accepted'].item() is False and journal['restart_allowed'].item() is False and j['accepted'] is False and j['restart_allowed'] is False,'uncommitted journal became accepted')
    values=dict(particle_id=pa.particle_id,position_m=positions(pa),velocity_m_s=np.column_stack((pa.u,pa.v,pa.w)),
        density_kg_m3=pa.rho,stress_Pa=pa.S.reshape(-1,3,3),temperature_K=pa.temperature,plastic_strain=pa.ep,
        contact_history=state['contact_history'])
    for key,value in values.items():require(np.array_equal(journal[key],value),'journal field differs: '+key)
    require(np.allclose(journal['displacement_m'],positions(pa)-positions(previous),rtol=0,atol=1e-18),'journal displacement differs')
    require(clean(energy['mechanical'])==clean(state['energy'])==clean(j['energy']) and
        clean(energy['thermal'])==clean(j['thermal']) and clean(energy['interval'])==clean(state['interval_energy']),'saved energy/thermal journals differ')
    actual=float(np.max(np.abs(journal['force_N'][np.asarray(state['free'],bool)])))
    require(actual<=j['metrics']['force_limit_N'] and j['metrics']['scaled_residual']<=1,'saved complete force not converged')
    e=energy['mechanical'];th=energy['thermal'];last=state['interval_energy']['records'][-1]
    require(last=={k:e[k] for k in last},'interval last row differs from actual raw energy')
    thermal=math.fsum(pa.m*p['cp']*(pa.temperature-previous.temperature))
    plastic=math.fsum(pa.plastic_work-previous.plastic_work)
    limit=1e-12+1e-10*abs(th['total_heat_input_J'])
    require(abs(thermal-th['measured_thermal_increment_J'])<=limit and abs(thermal-th['total_heat_input_J'])<=limit,'global saved heat capacity balance failed')
    require(abs(plastic-th['plastic_work_J'])<=limit and abs(th['plastic_heat_J']-p['plastic_heat_fraction']*plastic)<=limit,'actual plastic work/heat ledger differs')
    require(abs(th['internal_conduction_J'])<=limit,'internal pair heat not conservative')
    observation=observe_square_bore_crossings(previous_state['bore_observation'],positions(previous),positions(pa),pa.particle_id,
        pa.component_id,pa.m,previous_state['time_s'],state['time_s'],lambda t:pose(p,t))
    require(clean(observation['next_state'])==clean(state['bore_observation']),'saved strict-cross observer does not replay')
    shell=reader.read(sd/'continuous_shell_audit.json');geometry=reader.read(sd/'geometry_audit.json')
    require(shell['status']=='CERTIFIED_NO_HIDDEN_SHELL_CONTACT' and not shell['hidden_contacts'],'particle shell certificate failed')
    require(not geometry['transient_contacts_missed_by_endpoints'] and set(geometry['new_non_bore_features'])<={'SHOULDER'},'unqualified particle contact event')
    require(reader.read(sd/'scalar_batch_equivalence.json')['status']=='PASS_SAME_FILTERED_COMPLETE_MAP','scalar/batch gate absent')
    require(reader.read(sd/'sampled_nonbore_events.json')==row['events'],'contact event history differs')
    return dict(pa=pa,state=state,graph=saved,row=row,summary=dict(**state_stats(pa,state),**graph,
        checkpoint=path.relative_to(ROOT).as_posix(),checkpoint_sha256=row['sha256'],accepted=False,wall_seconds=row['wall_seconds'],
        force_residual_N=actual,global_raw_residual_J=e['global_raw_residual_J'],absolute_work_scale_J=e['absolute_work_scale_J'],
        external_work_J=e['prescribed_boundary_work_J']+e['tool_work_J'],plastic_work_J=plastic,
        thermal_capacity_increment_J=thermal,heat_balance_residual_J=thermal-th['total_heat_input_J']))


def compare_arrays(fine,coarse,selected,p):
    require(np.array_equal(fine.particle_id,coarse.particle_id),'comparison particle order differs')
    norm=np.linalg.norm;dx=norm(positions(fine)-positions(coarse),axis=1);dT=np.abs(fine.temperature-coarse.temperature)
    global_values=dict(position_max_m=float(dx.max()),temperature_max_K=float(dT.max()),
        stress_relative=float(norm(fine.S-coarse.S)/max(norm(fine.S),1e6*np.sqrt(len(fine.x)))),
        plastic_relative=float(norm(fine.ep-coarse.ep)/max(norm(fine.ep),1e-5*np.sqrt(len(fine.x)))))
    global_values['pass']=bool(global_values['position_max_m']<=.01*p['dx'] and global_values['temperature_max_K']<=.5 and
        global_values['stress_relative']<=.1 and global_values['plastic_relative']<=.2)
    index=np.flatnonzero(np.isin(fine.particle_id,selected));require(len(index)==len(selected)>0,'local screen identity incomplete')
    sr=norm(fine.S.reshape(-1,9)[index]-coarse.S.reshape(-1,9)[index],axis=1)/np.maximum(norm(fine.S.reshape(-1,9)[index],axis=1),1e6)
    ep=np.abs(fine.ep[index]-coarse.ep[index])/np.maximum(np.abs(fine.ep[index]),1e-5)
    rows=[dict(particle_id=int(fine.particle_id[k]),position_error_m=float(dx[k]),temperature_error_K=float(dT[k]),
               stress_relative=float(sr[n]),plastic_relative=float(ep[n])) for n,k in enumerate(index)]
    return global_values,dict(particle_ids=selected,rows=rows,position_limit_m=5e-6,temperature_limit_K=.5,
        stress_relative_limit=.1,stress_denominator_floor_Pa=1e6,plastic_relative_limit=.2,plastic_denominator_floor=1e-5,
        passed=bool(np.all(dx[index]<=5e-6) and np.all(dT[index]<=.5) and np.all(sr<=.1) and np.all(ep<=.2)))


def replay_screens(reader,cases,p,parent_state):
    parent_screen=reader.read(PARENT.parent.parent/'pair_016_temporal_screen.json')
    seeds=set(parent_screen['cumulative_shoulder_ids'])|set(parent_state['bore_observation']['crossing_ids']);reports=[]
    fine=cases.get('fine',[]);coarse=cases.get('coarse',[])
    for path in sorted((BASE/'fine').glob('pair_*_temporal_screen.json')):
        screen=reader.read(path);number=screen['pair'];require(number==len(reports)+1 and 2*number<=len(fine) and number<=len(coarse),'unbound time screen')
        a,b=fine[2*number-2:2*number];c=coarse[number-1]
        require(screen['parent_sha256']==PARENT_SHA and screen['coarse_checkpoint_sha256']==c['row']['sha256'] and
            screen['fine_candidate_sha256']==[a['row']['sha256'],b['row']['sha256']] and
            screen['time_s']==b['state']['time_s']==c['state']['time_s'],'pair endpoint/hash binding differs')
        for item in (a,b,c):
            seeds.update(e['particle_id'] for e in item['row']['events']);seeds.update(item['row']['blocked_particle_ids']);seeds.update(item['state']['bore_observation']['crossing_ids'])
        require(screen['cumulative_event_and_blocked_ids']==sorted(seeds),'pair cumulative locality seeds differ')
        selected=set(seeds)
        for item in (a,b,c):
            g=item['graph'];ids=g['particle_id']
            for left,right in (('i','j'),('raw_i','raw_j')):
                x,y=ids[g[left]],ids[g[right]];mask=np.isin(x,list(seeds))|np.isin(y,list(seeds));selected.update(np.r_[x[mask],y[mask]].tolist())
        glob,local=compare_arrays(b['pa'],c['pa'],sorted(selected),p)
        require(clean(glob)==screen['global_screen'] and clean(local)==screen['local_screen'],'numeric global/local screen does not replay')
        match=b['state']['bore_observation']['crossing_ids']==c['state']['bore_observation']['crossing_ids']
        require(screen['crossing_ids_match']==match and screen['fine_crossing_ids']==b['state']['bore_observation']['crossing_ids'] and
            screen['coarse_crossing_ids']==c['state']['bore_observation']['crossing_ids'],'strict crossing comparison differs')
        same=b['row']['graph_fingerprint']==c['row']['graph_fingerprint'];require(screen['endpoint_graphs_match']==same,'topology comparison differs')
        reports.append(dict(pair=number,passed=glob['pass'] and local['passed'] and match,sha256=digest(path),
            global_screen=glob,local_particle_count=len(selected),local_max_stress_relative=max(x['stress_relative'] for x in local['rows']),
            local_max_plastic_relative=max(x['plastic_relative'] for x in local['rows']),endpoint_graphs_match=same))
    return reports


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='outputs/visible_support_20261004/readback_v1');args=parser.parse_args()
    output=(ROOT/args.output).resolve();require(output.is_relative_to(ROOT),'output outside project');output.mkdir(parents=True,exist_ok=False)
    begin=time.perf_counter();reader=Reader();report=dict(status='READBACK_RUNNING',classification='SAVED_VISIBLE_PROCESS_EVIDENCE_AUDIT_NOT_FRESH_PHYSICS',
        command=[sys.executable,*sys.argv],physical_solver_evaluations=0,nnps_evaluations=0,production_authorized=False)
    try:
        reader.bind(Path(__file__));reader.bind(PARENT,PARENT_SHA);reader.bind(PRIOR,PRIOR_SHA)
        initial,parent=reader.checkpoint(PARENT);p=parent['parameters'];prefix=parent['interval_energy']['records']
        require(len(prefix)==parent['step']==OLD_STEPS and abs(parent['time_s']-START)<1e-14,'old true parent differs')
        validate_interval_energy_state(parent['interval_energy']);validate_bore_crossing_state(parent['bore_observation'])
        reader.sources(reader.read(PRIOR)['pre_occlusion_source_and_evidence_sha256'])
        cases={};stages={};common=None
        for stage in ('cost','coarse','fine'):
            directory=BASE/stage
            if not directory.exists():stages[stage]=dict(status='NOT_RUN',accepted_steps=0);cases[stage]=[];continue
            meta=reader.read(directory/'manifest.json');require(meta['status']!='RUNNING','stage is not terminal: '+stage)
            history=reader.read(directory/'accepted_history.json');rows=reader.read(directory/'candidate_history.json')
            source=meta['source_sha256'];reader.sources(source)
            if common is None:common=source
            else:require(common==source,'stage source collections differ')
            require(meta['classification']==CLASS and meta['parent_sha256']==PARENT_SHA,'stage parent/model differs')
            supervision=reader.read(ROOT/f'logs/visible_support_20261004/{stage}_supervision.json')
            require(supervision['limit_s']==600 and supervision['external_wall_s']<=605,'external budget evidence differs')
            require((supervision['exit_code']==0)==(meta['status']==PASS[stage]),'manifest/exit terminal status differs')
            require(meta['forecast_wall_s']<=450 and meta['wall_seconds']<=540.5,'stage budget exceeded without rejection')
            transfer,ts=reader.checkpoint(directory/'MODEL_TRANSITION_NOT_RESTARTABLE.npz');properties_equal(transfer,initial)
            validate_state(transfer,ts,initial,prefix,p,0,.000125,False)
            require(ts['interval_energy']==parent['interval_energy'] and clean(ts['bore_observation'])==clean(parent['bore_observation']),'transition reset an observer/energy history')
            previous,previous_state=initial,parent;records=[];direct=PARENT_SHA
            for index,row in enumerate(rows,1):
                require(row['step']==index and row['parent_sha256']==direct,'candidate direct parent mismatch')
                item=candidate_audit(reader,stage,directory,row,previous,previous_state,initial,prefix,p,source);records.append(item)
                previous,previous_state=item['pa'],item['state']
                direct=row['sha256']
                if stage=='fine' and index%2==0:
                    published=directory/f'step_{index:03d}/accepted.npz'
                    if published.exists():direct=digest(published)
            require(meta.get('candidate_steps',0)==len(rows),'manifest candidate count differs')
            if rows:require(meta['candidate_time_s']==rows[-1]['time_s'],'manifest candidate endpoint differs')
            if meta['status']==PASS[stage]:require(len(rows)=={'cost':1,'coarse':8,'fine':16}[stage],'complete fixed stage length differs')
            accepted_files=sorted(directory.glob('step_*/accepted.npz'))
            if stage!='fine':require(not history and not accepted_files and meta['accepted_steps']==0,'cost/coarse became accepted')
            extra=[]
            known={Path(r['checkpoint']).parent.as_posix() for r in rows}
            for sd in sorted(directory.glob('step_*')):
                if sd.name not in known:
                    files=list(sd.rglob('*'));bound=[reader.bind(f).relative_to(ROOT).as_posix() for f in files if f.is_file()]
                    extra.append(dict(directory=sd.relative_to(ROOT).as_posix(),qualification='UNREGISTERED_FAILED_OR_PARTIAL_CANDIDATE_NOT_ACCEPTED',evidence_files=bound))
            cases[stage]=records;stages[stage]=dict(status=meta['status'],wall_seconds=meta['wall_seconds'],external_wall_seconds=supervision['external_wall_s'],
                exit_code=supervision['exit_code'],forecast_wall_s=meta['forecast_wall_s'],candidate_count=len(rows),accepted_steps=meta['accepted_steps'],
                accepted_time_s=meta['accepted_time_s'],candidates=[x['summary'] for x in records],unregistered_partial_evidence=extra,
                manifest_sha256=digest(directory/'manifest.json'),history=history)
        if stages['coarse']['status']!=PASS['coarse']:require(stages['fine']['status']=='NOT_RUN','fine ran after unqualified coarse')
        if stages['cost']['status']!=PASS['cost']:require(stages['coarse']['status']=='NOT_RUN','coarse ran after unqualified cost')
        screens=replay_screens(reader,cases,p,parent);accepted=[];direct=PARENT_SHA
        for index,path in enumerate(sorted((BASE/'fine').glob('step_*/accepted.npz')),1):
            require(path.parent.name==f'step_{index:03d}' and index<=len(cases['fine']),'accepted prefix is not contiguous')
            candidate=cases['fine'][index-1];pa,state=reader.checkpoint(path);properties_equal(pa,candidate['pa'])
            validate_state(pa,state,initial,prefix,p,index,.0000625,True)
            require(state['interval_energy']==candidate['state']['interval_energy'] and clean(state['bore_observation'])==clean(candidate['state']['bore_observation']),'published physical histories differ from candidate')
            pair=(index+1)//2;require(pair<=len(screens) and screens[pair-1]['passed'],'accepted half has no qualified pair')
            require(state['continuation_parent_sha256']==direct and state['branch_origin_sha256']==PARENT_SHA and
                state['implementation_source_sha256']==common and state['controller_source_sha256']==common and
                state['window_temporal_screen_sha256']==screens[pair-1]['sha256'] and
                state['window_coarse_reference_sha256']==cases['coarse'][pair-1]['row']['sha256'],'published parent/source/screen binding differs')
            direct=digest(path);entry=dict(step=index,global_step=state['step'],time_s=state['time_s'],sha256=direct,
                checkpoint=path.relative_to(ROOT).as_posix());accepted.append(entry)
            history=stages['fine']['history']
            if index<=len(history):
                row=history[index-1]
                require(row['step']==index and row['global_step']==state['step'] and row['time_s']==state['time_s'] and
                    row['sha256']==direct and row['parent_sha256']==state['continuation_parent_sha256'] and
                    row['temporal_screen_sha256']==state['window_temporal_screen_sha256'] and row['coarse_reference_sha256']==state['window_coarse_reference_sha256'],'published history entry differs')
        history=stages['fine'].get('history',[]);require(len(history)<=len(accepted),'history lists missing accepted file')
        if stages['fine']['status']==PASS['fine']:require(len(accepted)==len(history)==16 and all(x['passed'] for x in screens) and len(screens)==8,'successful fine stage incomplete')
        endpa,endstate=reader.checkpoint(ROOT/accepted[-1]['checkpoint']) if accepted else (initial,parent)
        new=cases['fine'][:len(accepted)];energy=endstate['interval_energy'];newrows=energy['records'][OLD_STEPS:]
        frames=[(initial,parent)]+[reader.checkpoint(ROOT/x['checkpoint']) for x in accepted]
        times=np.array([state['time_s'] for _,state in frames]);poses=[pose(p,t) for t in times]
        trajectory=output/'trajectory.npz'
        np.savez_compressed(trajectory,time_s=times,particle_id=initial.particle_id,
            position_m=np.array([positions(pa) for pa,_ in frames]),temperature_K=np.array([pa.temperature for pa,_ in frames]),
            stress_Pa=np.array([pa.S.reshape(-1,3,3) for pa,_ in frames]),plastic_strain=np.array([pa.ep for pa,_ in frames]),
            tool_origin_m=np.array([v['origin'] for v in poses]),tool_rotation=np.array([v['rotation'] for v in poses]),
            mass_kg=initial.m,component_id=initial.component_id,tracer_type=initial.tracer_type,wholebar_tag=initial.wholebar_tag)
        response=[];before_pose=poses[0];after_pose=poses[-1]
        before_local=(positions(initial)-before_pose['origin'])@before_pose['rotation']
        after_local=(positions(endpa)-after_pose['origin'])@after_pose['rotation']
        for pid in endstate['bore_observation']['crossing_ids']:
            k=int(np.flatnonzero(initial.particle_id==pid)[0]);r0=float(np.linalg.norm(before_local[k,:2]));r1=float(np.linalg.norm(after_local[k,:2]))
            response.append(dict(particle_id=pid,mass_kg=float(endpa.m[k]),
                plastic_strain_start=float(initial.ep[k]),plastic_strain_end=float(endpa.ep[k]),plastic_strain_increment=float(endpa.ep[k]-initial.ep[k]),
                temperature_start_K=float(initial.temperature[k]),temperature_end_K=float(endpa.temperature[k]),temperature_increment_K=float(endpa.temperature[k]-initial.temperature[k]),
                tool_radius_start_m=r0,tool_radius_end_m=r1,tool_radius_increment_m=r1-r0,
                tool_z_start_m=float(before_local[k,2]),tool_z_end_m=float(after_local[k,2]),tool_z_increment_m=float(after_local[k,2]-before_local[k,2]),
                strict_square_bore_projection_start=bool(np.all(np.abs(before_local[k,:2])<p['bore_width']/2)),
                strict_square_bore_projection_end=bool(np.all(np.abs(after_local[k,:2])<p['bore_width']/2)),
                square_bore_projection_margin_end_m=float(p['bore_width']/2-np.max(np.abs(after_local[k,:2])))))
        report.update(status='PASS_SAVED_EVIDENCE_READBACK_NOT_FRESH_PHYSICS',stages=stages,pair_screens=screens,
            initial_accepted_global_steps=OLD_STEPS,added_accepted_steps=len(accepted),accepted_global_steps=endstate['step'],
            accepted_time_s=endstate['time_s'],accepted_time_advanced_s=endstate['time_s']-parent['time_s'],
            actual_accepted_frame_count_including_independent_initial=endstate['step']+1,
            inherited_interval_never_reset=True,accepted_history=accepted,unrecorded_published_halves_reconciled=len(accepted)-len(history),
            final_accepted_fields=state_stats(endpa,endstate),crossing_events=endstate['bore_observation']['events'],
            crossed_particle_response=response,trajectory=dict(path=trajectory.relative_to(ROOT).as_posix(),sha256=digest(trajectory),
                frame_count=len(frames),scope='ACTUAL_TRUE_PARENT_PLUS_NEW_FINE_ACCEPTED_ONLY_NO_INTERPOLATION_NO_COARSE',
                inherited_432_frame_readback_summary=PRIOR.relative_to(ROOT).as_posix(),inherited_summary_sha256=PRIOR_SHA),
            complete_prefix_energy={k:v for k,v in energy.items() if k!='records'},
            new_accepted_energy=dict(signed_raw_residual_J=math.fsum(r['global_raw_residual_J'] for r in newrows),
                absolute_step_raw_residual_J=math.fsum(abs(r['global_raw_residual_J']) for r in newrows),
                absolute_work_scale_J=math.fsum(r['absolute_work_scale_J'] for r in newrows),
                external_work_J=math.fsum(x['summary']['external_work_J'] for x in new),
                plastic_work_J=math.fsum(x['summary']['plastic_work_J'] for x in new),
                measured_thermal_increment_J=math.fsum(x['summary']['thermal_capacity_increment_J'] for x in new)),
            limitations=['Saved graph and whole-path certificate inputs/coverage replayed; no new NNPS or geometry certificate evaluation.',
                'Full force and heat journals are bound to checkpoints; no fresh force/material/thermal solve was performed.',
                'Time screens replayed from actual saved arrays; no experimental or spatial-process validation is established.',
                'Cost and coarse candidates, and all failed/partial uncommitted candidates, are excluded from accepted states.'])
        reader.finish();report['evidence_sha256']=reader.evidence;report['source_entry_checks_including_repeats']=reader.source_entries
        code=0
    except Exception:
        report.update(status='FAIL_SAVED_EVIDENCE_READBACK',error=traceback.format_exc(),evidence_sha256=reader.evidence);code=2
    report.update(exit_code=code,wall_seconds=time.perf_counter()-begin)
    (output/'summary.json').write_text(json.dumps(clean(report),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report.get(k) for k in ('status','exit_code','added_accepted_steps','accepted_global_steps','accepted_time_s','wall_seconds','error')},indent=2))
    return code


if __name__=='__main__':raise SystemExit(main())

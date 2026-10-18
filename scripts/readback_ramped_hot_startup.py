"""Independent saved-only ramp-startup readback; diagnostic states never publish."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import readback_expanded_hot_campaign_window_v2 as prior
old=prior.old
ROOT=prior.ROOT
Reader,require,clean,digest=prior.Reader,prior.require,prior.clean,prior.digest
positions=old.positions
array_hash,ordered_pairs,sorted_pairs,topology=old.array_hash,old.ordered_pairs,old.sorted_pairs,old.topology
CLASS='INDEPENDENT_RAMPED_20MS_HOT_QS_VISIBLE_DIAGNOSTIC_V1'
SCOPE='NEW_918_PARTICLE_HOT_ZERO_TIME_SMOOTH_VELOCITY_STARTUP_ASSUMPTION_V1'
MOTION_CLASS='SMOOTHSTEP_VELOCITY_EXACT_ENDPOINTS_PIECEWISE_LINEAR_GEODESIC_STEP_V1'
PASS_STAGE='PASS_RAMP_PHYSICAL_STAGE_NOT_ACCEPTED'
PASS_READBACK='PASS_RAMP_CHUNK_SAVED_READBACK_NOT_ACCEPTED'
CONTRACT=ROOT/'docs/RAMPED_HOT_STARTUP_DIAGNOSTIC_CONTRACT_20261005.md'
CONTROLLER=ROOT/'scripts/run_ramped_hot_startup.py'
CHECK=ROOT/'outputs/ramped_hot_motion_20261005/v1/receipt.json'
OLD_COST=ROOT/'outputs/expanded_hot_domain_cost_20261004/v3/dt_15_625us/cost/manifest.json'
OLD_RANK=ROOT/'outputs/expanded_hot_qualification_20261005/v1/summary.json'
OLD_ROOT=ROOT/'outputs/expanded_hot_campaign_v2_20261005/window_001/fine/candidate_history.json'
INITIAL_REFERENCE=ROOT/'outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
INITIAL_REFERENCE_SHA='31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3'
DTS={'coarse':.0005,'fine':.00025}
COUNTS={'coarse':2,'fine':4}
ALLOWED_FEATURES=['SHOULDER','TOP']
project_path=prior.project_path
finite_properties=prior.finite_properties
correction_audit=prior.correction_audit
tight_path_arithmetic_audit=prior.tight_path_arithmetic_audit


class Motion:
    """Independent scalar arithmetic for the registered interval path."""
    def __init__(self,p,t0,dt):
        self.p=p;self.t0=t0;self.dt=dt;self.t1=t0+dt
        require(t0>=0 and dt>0 and np.isfinite([t0,dt,self.t1]).all(),'invalid ramp interval')
        self.e0=self.exposure(t0);self.e1=self.exposure(self.t1)
        self.de=self.increment(t0,self.t1)
        require(abs(p['omega']*(self.e1-self.e0))<math.pi,'ambiguous geodesic')
    @staticmethod
    def exposure(t):
        if t>=.02:return t-.01
        s=t/.02;return .02*s*s*s*(1.-.5*s)
    @staticmethod
    def factor(t):
        if t>=.02:return 1.
        s=t/.02;return s*s*(3.-2.*s)
    @staticmethod
    def increment(a,b):
        if a>=.02:return b-a
        if b>.02:return Motion.increment(a,.02)+b-.02
        x,y=a/.02,b/.02
        return (b-a)*((y*y+x*y+x*x)-.5*(y*y*y+x*y*y+x*x*y+x*x*x))
    def pose(self,t):
        require(self.t0<=t<=self.t1,'observer extrapolated ramp interval')
        if t==self.t0:e=self.e0
        elif t==self.t1:e=self.e1
        else:
            f=(t-self.t0)/(self.t1-self.t0);e=(1-f)*self.e0+f*self.e1
        angle=self.p['omega']*e;c,s=math.cos(angle),math.sin(angle)
        return dict(origin=np.array([self.p['traverse']*e,0.,self.p['layer']]),rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]]))
    def record(self):
        a,b=self.pose(self.t0),self.pose(self.t1)
        return dict(classification=MOTION_CLASS,t0_s=self.t0,t1_s=self.t1,dt_s=self.dt,
            exposure0_s=self.e0,exposure1_s=self.e1,exposure_increment_s=self.de,
            origin0_m=a['origin'].tolist(),origin1_m=b['origin'].tolist(),rotation0=a['rotation'].tolist(),rotation1=b['rotation'].tolist(),
            unwrapped_angle0_rad=self.p['omega']*self.e0,unwrapped_angle1_rad=self.p['omega']*self.e1,
            feed_displacement_m=-self.p['feed']*self.de,nominal_velocity_factor0=self.factor(self.t0),
            nominal_velocity_factor1=self.factor(self.t1),actual_step_average_velocity_factor=self.de/self.dt,
            path='LINEAR_WORLD_ORIGIN_Z_GEODESIC_ROTATION_LINEAR_ACTUATOR',global_ramp_not_used_inside_step=True)


def validate_state(pa,state,context,number,dt,restart):
    initial,parent,p=context['initial'],context['parent'],context['p']
    require(restart is False and state['classification']==CLASS and state['model_scope']==SCOPE and
        state['accepted'] is False and state['restart_allowed'] is False and state['diagnostic_only'] is True and
        state['accepted_time_advanced_s']==0.,'ramp diagnostic was relabelled or published')
    require(state['step']==parent['step']+number and abs(state['time_s']-parent['time_s']-number*dt)<1e-15,
        'branch/chunk time or step differs')
    for key in ('spatial_model_transition','branch_origin_sha256','expanded_initial_qualification'):
        require(key not in state,'old process origin or qualification was fabricated')
    for key in ('operator_variant','independent_initialization','parameters','parameter_provenance','nominal_parameters',
        'nominal_parameter_provenance','independent_domain_initialization','ramp_motion','source_sha256','ramp_component_source_sha256'):
        require(clean(state[key])==clean(parent[key]),'ramp invariant metadata changed: '+key)
    require(set(pa.properties)==set(initial.properties),'particle property schema changed')
    for key in old.INVARIANTS:require(np.array_equal(getattr(pa,key),getattr(initial,key)),'identity/mass changed: '+key)
    require(np.array_equal(state['free'],parent['free']),'physical axes changed');finite_properties(pa)
    require(np.all(pa.rho>0) and np.all(pa.ep>=0) and np.all(pa.temperature>=p['room']) and np.all(pa.temperature<p['melt']),
        'saved material domain failed')
    old.validate_interval_energy_state(state['interval_energy']);old.validate_bore_crossing_state(state['bore_observation'])
    rows=state['interval_energy']['records'];history=state['ramped_boundary_history']
    require(len(rows)==len(history)==state['step'] and rows[:parent['step']]==context['prefix'] and
        history[:parent['step']]==parent['ramped_boundary_history'],'energy or motion prefix reset')
    t=0.
    for r in history:
        m=Motion(p,t,r['dt_s']);require(r==m.record(),'ramp interval arithmetic does not replay');t=m.t1
    require(t==state['time_s']==state['bore_observation']['last_time_s'] and
        state['last_ramped_interval']==history[-1],'motion/observer terminal time differs')
    require(state['crossing_ids']==state['bore_observation']['crossing_ids'] and
        state['crossing_mass_kg']==state['bore_observation']['crossing_mass_kg'],'passive observer differs')


def graph_audit(reader,directory,pa,state,previous,p,row,motion):
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
    o0,o1=motion.pose(motion.t0),motion.pose(motion.t1)
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


def candidate_audit(reader, stage, directory, row, previous, previous_state, context, frozen):
    path = reader.bind(directory/row['checkpoint'], row['sha256'])
    pa, state = reader.checkpoint(path)
    require(path.parent == directory/f"step_{row['step']:03d}" and
            path.name == 'NOT_RESTARTABLE_candidate.npz', 'candidate path is not the registered uncommitted file')
    required_artifacts = ('NOT_RESTARTABLE_candidate.npz', 'FULL_UNCOMMITTED_candidate.npz',
        'FULL_UNCOMMITTED_candidate.json', 'candidate_evidence.npz', 'energy.json',
        'visible_graph.npz', 'visible_graph.json', 'visible_path_audit.json',
        'continuous_shell_audit.json', 'geometry_audit.json', 'scalar_batch_equivalence.json',
        'sampled_events.json', 'motion.json', 'tool_pose.json', 'response.json', 'transport.json', 'fields.vtk')
    for name in required_artifacts:
        require((path.parent/name).relative_to(directory).as_posix() in row['artifacts'],
                'candidate lacks required artifact hash: '+name)
    for relative, expected in row['artifacts'].items():
        actual = (directory/relative).resolve()
        require(actual.is_relative_to(path.parent.resolve()), 'candidate artifact escaped its step')
        reader.bind(actual, expected)
    dt = context['dt'];motion=Motion(context['p'],previous_state['time_s'],dt)
    require(row['dt_s'] == dt and row['global_step'] == context['parent']['step']+row['step'] and
            row['time_s'] == state['time_s'], 'candidate history time/step differs')
    require(abs(state['time_s']-previous_state['time_s']-dt) < 1e-15,
            'candidate actual immediate-parent time differs')
    validate_state(pa, state, context, row['step'], dt, False)
    require(state['implementation_source_sha256'] == frozen, 'candidate implementation source differs')
    require(state['interval_energy']['records'][:-1] == previous_state['interval_energy']['records'],
            'candidate interval direct parent differs')
    sd = path.parent
    graph_meta = reader.read(sd/'visible_graph.json')
    require(graph_meta['parent_time_s'] == previous_state['time_s'] and
            graph_meta['endpoint_time_s'] == state['time_s'], 'graph path time differs')
    saved, graph = graph_audit(reader, sd, pa, state, previous, context['p'], row,motion)
    graph.update(correction_audit(saved))
    path_report = reader.read(sd/'visible_path_audit.json')
    tight_path_arithmetic_audit(path_report, saved)
    require(path_report['analytic_batch_queries'] == len(path_report['query_log']) and
            path_report['analytic_batch_queries'] == path_report['strict_geometry_batch_queries']+
                path_report['enclosure_batch_queries'] and
            path_report['analytic_batch_queries'] <= path_report['maximum_queries'] <= 256 and
            path_report['deepest_subdivision'] <= path_report['maximum_depth'] <= 32 and
            path_report['edge_deletions'] == 0 and path_report['physical_solver_evaluations'] == 0,
            'visible path certificate scope/budget differs')
    require(path_report['certified_edge_interval_count'] ==
            sum(len(c['edge_indices']) for c in path_report['certificates']),
            'certificate interval count differs')
    journal = reader.npz(sd/'FULL_UNCOMMITTED_candidate.npz')
    j = reader.read(sd/'FULL_UNCOMMITTED_candidate.json')
    energy = reader.read(sd/'energy.json')
    require(journal['accepted'].item() is False and journal['restart_allowed'].item() is False and
            j['accepted'] is False and j['restart_allowed'] is False, 'journal became accepted')
    values = dict(particle_id=pa.particle_id, position_m=positions(pa),
        velocity_m_s=np.column_stack((pa.u, pa.v, pa.w)), density_kg_m3=pa.rho,
        stress_Pa=pa.S.reshape(-1, 3, 3), temperature_K=pa.temperature, plastic_strain=pa.ep,
        contact_history=state['contact_history'])
    for key, value in values.items():
        require(np.array_equal(journal[key], value), 'full journal field differs: '+key)
    for key, value in journal.items():
        if value.dtype.kind in 'fc': require(np.isfinite(value).all(), 'nonfinite journal: '+key)
    require(np.allclose(journal['displacement_m'], positions(pa)-positions(previous), rtol=0, atol=1e-18),
            'journal displacement differs')
    require(np.array_equal(journal['force_N'], journal['internal_force_N']+journal['contact_force_N']),
            'complete force journal components do not sum')
    require(np.array_equal(journal['contact_force_N'], np.asarray(j['contact']['force_N'])) and
            np.array_equal(journal['contact_history'], np.asarray(j['contact']['history'])),
            'contact arrays differ between complete journals')
    require(clean(energy['mechanical']) == clean(state['energy']) == clean(j['energy']) and
            clean(energy['thermal']) == clean(j['thermal']) and
            clean(energy['interval']) == clean(state['interval_energy']), 'energy/thermal journal differs')
    actual = float(np.max(np.abs(journal['force_N'][np.asarray(state['free'], bool)])))
    require(actual <= j['metrics']['force_limit_N'] and j['metrics']['scaled_residual'] <= 1,
            'saved complete residual did not converge')
    e, th, p = energy['mechanical'], energy['thermal'], context['p']
    require(j['contact']['sweep']['certified'] is True and
        float(np.asarray(j['contact']['penetration_m']).max())<=.1*p['dx'] and
        float(np.max(np.abs(pa.rho/p['rho0']-1)))<=p['density_limit'], 'original contact/density gates failed')
    require(e['registered_energy_limit_J']==1e-7+.05*e['absolute_work_scale_J'] and
        abs(e['global_raw_residual_J'])<=e['registered_energy_limit_J'], 'raw mechanical gate failed')
    require(np.all(pa.ep>=previous.ep) and np.max(np.linalg.norm(journal['displacement_m'],axis=1))<=.1*p['dx'],
        'plastic strain decreased or incremental geometry domain exceeded')
    last = state['interval_energy']['records'][-1]
    require(last == {k:e[k] for k in last}, 'interval last record differs from raw energy')
    thermal = math.fsum(pa.m*p['cp']*(pa.temperature-previous.temperature))
    plastic = math.fsum(pa.plastic_work-previous.plastic_work)
    limit = 1e-12+1e-10*abs(th['total_heat_input_J'])
    require(abs(thermal-th['measured_thermal_increment_J']) <= limit and
            abs(thermal-th['total_heat_input_J']) <= limit, 'saved heat capacity balance failed')
    require(abs(plastic-th['plastic_work_J']) <= limit and
            abs(th['plastic_heat_J']-p['plastic_heat_fraction']*plastic) <= limit,
            'actual plastic work/heat differ')
    require(abs(th['internal_conduction_J']) <= limit, 'internal pair heat is not conservative')
    observation = old.observe_square_bore_crossings(previous_state['bore_observation'],
        positions(previous), positions(pa), pa.particle_id, pa.component_id, pa.m,
        previous_state['time_s'], state['time_s'], motion.pose)
    require(clean(observation['next_state']) == clean(state['bore_observation']),
            'strict crossing observer does not replay')
    shell = reader.read(sd/'continuous_shell_audit.json')
    geometry = reader.read(sd/'geometry_audit.json')
    require(shell['status'] == 'CERTIFIED_NO_HIDDEN_SHELL_CONTACT' and not shell['hidden_contacts'],
            'continuous particle shell certificate failed')
    require(not geometry['transient_contacts_missed_by_endpoints'] and
            set(geometry['new_non_bore_features']) <= set(ALLOWED_FEATURES), 'unqualified geometry event')
    require(reader.read(sd/'scalar_batch_equivalence.json')['status'] == 'PASS_SAME_RAMP_FILTERED_COMPLETE_MAP',
            'same-filtered scalar/batch qualification absent')
    require(reader.read(sd/'sampled_events.json') == row['events'], 'event history differs')
    require(reader.read(sd/'motion.json')==motion.record()==state['last_ramped_interval'], 'motion journal mismatch')
    require(reader.read(sd/'tool_pose.json')==clean(motion.pose(motion.t1)), 'tool pose journal mismatch')
    require(reader.read(sd/'transport.json')==clean(state['bore_observation']), 'observer journal mismatch')
    prescribed=np.zeros_like(positions(pa));prescribed[previous.actuator.astype(bool),2]=-p['feed']*motion.de
    free=np.asarray(state['free'],bool)
    require(np.allclose(journal['displacement_m'][~free],prescribed[~free],rtol=0,atol=1e-16),'ramped physical axes violated')
    loaded=np.linalg.norm(np.asarray(j['contact']['force_by_feature_N'])[:,:7],axis=2).max(axis=1)>j['metrics']['force_limit_N']
    nonbase=np.asarray(j['contact']['active'])[:,:7].any(axis=1)
    response=dict(nonbase_loaded_ids=pa.particle_id[loaded & nonbase].astype(int).tolist(),
        new_plastic_ids=pa.particle_id[pa.ep>previous.ep].astype(int).tolist(),
        max_plastic_increment=float(np.max(pa.ep-previous.ep)),max_plastic_strain=float(pa.ep.max()),
        temperature_min_K=float(pa.temperature.min()),temperature_max_K=float(pa.temperature.max()),
        crossing_ids=state['bore_observation']['crossing_ids'])
    require(response==row['response']==reader.read(sd/'response.json'),'saved loaded/plastic response mismatch')
    # Bind sampled endpoint events to the saved complete contact response;
    # this arithmetic check never recalculates geometry or contact forces.
    features = j['contact']['features']
    endpoint_load = np.asarray(j['contact']['endpoint_normal_load_N'])
    require(endpoint_load.shape == (len(pa.particle_id), len(features)),
            'contact endpoint load shape differs')
    observed_endpoint = {(int(e['particle_id']), e['feature']) for e in row['events'] if e['fraction'] == 1. and e['feature'] in ALLOWED_FEATURES}
    actual_endpoint = {(int(pa.particle_id[k]), features[f])
        for k, f in np.argwhere(endpoint_load > 0) if features[f] in ALLOWED_FEATURES}
    require(observed_endpoint == actual_endpoint, 'sampled endpoint contact events do not match raw contact load')
    top_index = features.index('TOP')
    top_contact = {}
    for key in ('normal_potential0_J', 'normal_potential1_J', 'tangential_potential0_J',
            'tangential_potential1_J', 'friction_dissipation_J', 'algorithmic_dissipation_J',
            'release_energy_J', 'normal_particle_work_J', 'normal_wall_work_J',
            'particle_work_J', 'wall_work_J'):
        values = np.asarray(j['contact'][key])
        require(values.shape == endpoint_load.shape and np.isfinite(values).all(),
                'invalid saved contact feature field: '+key)
        top_contact[key] = math.fsum(values[:, top_index])
    return dict(pa=pa, state=state, graph=saved, row=row, response=response,
        summary=dict(**old.state_stats(pa, state), **graph,
            checkpoint=path.relative_to(ROOT).as_posix(), checkpoint_sha256=row['sha256'],
            top_contact=top_contact,
            accepted=False, wall_seconds=row['wall_seconds'], force_residual_N=actual,
            global_raw_residual_J=e['global_raw_residual_J'], absolute_work_scale_J=e['absolute_work_scale_J'],
            external_work_J=e['prescribed_boundary_work_J']+e['tool_work_J'], plastic_work_J=plastic,
            thermal_capacity_increment_J=thermal, heat_balance_residual_J=thermal-th['total_heat_input_J']))


SOURCE_MAP_SHA='7f8c68c55925bb1909eddd675257adffa95fe4d39cd83de6c1e6edef144b485d'
CHECK_SHA='b0110de0d92d0ed364862bbe8a5acb8b2eaafc3b1e2e16131d7956cd5e6e41f8'


def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def mapping_add(mapping,values):
    for name,value in values.items():
        require(name not in mapping or mapping[name]==value,'conflicting evidence map: '+name)
        mapping[name]=value


def tree_hashes(folder):
    return {rel(p):digest(p) for p in folder.rglob('*') if p.is_file()}


def protocol_audit(reader,base,folder,stage,chunk,frozen,parent_path,completed):
    protocol=reader.read(folder/'protocol.json');meta=reader.read(folder/'manifest.json')
    supervision=reader.read(folder/'supervision.json');reader.bind(folder/'run.log',supervision['log_sha256'])
    reader.bind(folder/'manifest.json',supervision['manifest_sha256'])
    require(supervision['source_inputs_still_match'] and supervision['source_error'] is None,'source drift during physical stage')
    require(protocol['source_sha256']==frozen and all(meta[k]==v for k,v in protocol.items()),'stage protocol changed')
    reader.sources(protocol['source_sha256']);reader.sources(protocol['reference_sha256']);reader.sources(meta['artifacts_sha256'])
    refs={};elapsed=0.;measured=max(v['wall_seconds'] for v in reader.read(OLD_ROOT))
    for previous in completed:
        pd,pm,ps=previous
        elapsed+=ps['wall_seconds']
        for name in ('supervision.json','manifest.json','candidate_history.json'):
            q=pd/name
            if q.exists():refs[rel(q)]=digest(q)
        rows=reader.read(pd/'candidate_history.json')
        if rows:measured=max(measured,max(r['wall_seconds'] for r in rows),pm['wall_seconds']/len(rows))
    if stage=='cost':
        c,r=reader.read(OLD_COST),reader.read(OLD_RANK)
        forecast=30+1.5*(2*c['full_map_wall_s']+c['tangent_linear_wall_s']+r['wall_seconds'])+120
    else:
        if chunk>1:
            ap=base/f'chunk_{chunk-1:03d}/readback/summary.json';a=reader.read(ap)
            require(a['status']==PASS_READBACK and a['temporal_screens_all_pass'] and a['chunk']==chunk-1 and
                a['base']==rel(base) and a['accepted_time_advanced_s']==0.,'previous chunk lacks matching independent gate')
            reader.sources(a['evidence_sha256']);mapping_add(refs,a['evidence_sha256']);refs[rel(ap)]=digest(ap)
        if stage=='fine':mapping_add(refs,tree_hashes(base/f'chunk_{chunk:03d}/coarse'))
        mapping_add(refs,tree_hashes(base/'cost'));refs[rel(parent_path)]=digest(parent_path)
        forecast=30+1.5*COUNTS[stage]*(measured+1.)
    require(protocol['reference_sha256']==refs,'stage reference set does not exactly replay')
    internal=min(540.,1200.-elapsed-5.);external=min(600.,1200.-elapsed)
    expected=dict(schema='RAMPED_HOT_STARTUP_DIAGNOSTIC_V1',stage=stage,chunk=chunk,
        parent_checkpoint=rel(parent_path) if parent_path else None,parent_sha256=digest(parent_path) if parent_path else None,
        step_count=0 if stage=='cost' else COUNTS[stage],dt_s=.0005 if stage=='cost' else DTS[stage],
        endpoint_s=0. if stage=='cost' else chunk*.001,comparison_scope='COMMON_INITIAL_CONTINUOUS_TWO_RESOLUTION',
        measured_candidate_cost_s=measured,forecast_s=forecast,previous_supervised_wall_s=elapsed,total_wall_limit_s=1200.,
        internal_limit_s=internal,external_limit_s=external,accepted=False,restart_allowed=False,accepted_time_advanced_s=0.)
    for k,v in expected.items():require(protocol[k]==v,'registered stage/budget differs: '+k)
    require(forecast<=450 and elapsed+forecast<=1200 and forecast<=internal-15,'forecast budget failed')
    require(supervision['accepted'] is False and supervision['restart_allowed'] is False and
        supervision['internal_limit_s']==internal and supervision['external_limit_s']==external and
        meta['accepted_steps']==0 and meta['accepted_time_advanced_s']==0.,'diagnostic publication or cap change')
    cmd=supervision['command']
    require(cmd[1:]==['-B',str(CONTROLLER),'--stage',stage,'--chunk',str(chunk),'--output',rel(base),'--worker'],
        'supervised worker command differs')
    if supervision['exit_code']==0:require(meta['wall_seconds']<internal and supervision['wall_seconds']<external,'successful stage exceeded cap')
    actual={rel(p):digest(p) for p in folder.rglob('*') if p.is_file() and p.name not in ('manifest.json','run.log','supervision.json')}
    if stage=='cost':mapping_add(actual,tree_hashes(base/'initial'))
    require(meta['artifacts_sha256']==actual,'terminal artifact ledger is incomplete')
    completed.append((folder,meta,supervision))
    return protocol,meta,supervision


def initial_audit(reader,base,completed):
    folder=base/'cost';frozen=reader.read(folder/'protocol.json')['source_sha256']
    require(hashlib.sha256(json.dumps(frozen,sort_keys=True,separators=(',',':')).encode()).hexdigest()==SOURCE_MAP_SHA,
        'registered frozen ramp source set differs')
    reader.bind(CHECK,CHECK_SHA);q=reader.read(CHECK)
    require(q['status']=='PASS_RAMPED_MOTION_COMPONENT_QUALIFICATION_NOT_PROCESS_ACCEPTANCE' and
        len(q['checks'])==8 and all(v['status']=='PASS' for v in q['checks']), 'motion component qualification absent')
    reader.sources(q['source_and_input_sha256'])
    for path in (CONTROLLER,CONTRACT,CHECK,OLD_COST,OLD_RANK,OLD_ROOT):
        require(frozen[rel(path)]==digest(path),'required implementation/input absent')
    protocol,meta,supervision=protocol_audit(reader,base,folder,'cost',1,frozen,None,completed)
    require(meta['status']=='PASS_RAMP_COST_INITIAL_RANK_SAME_MAP_NOT_ACCEPTED' and supervision['exit_code']==0 and
        meta['initial_state_unchanged'] and meta['particle_count']==918 and meta['free_dof']==1961 and
        meta['physical_full_map_evaluations']==2 and meta['stable_timestep_s']=='NOT_ESTABLISHED','cost/rank stage did not qualify')
    path=base/'initial/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
    require(meta['initial_checkpoint']==rel(path),'initial path mismatch');reader.bind(path,meta['initial_sha256'])
    pa,state=reader.checkpoint(path);require(clean(state)==reader.read(base/'initial/state.json'),'initial state JSON mismatch')
    reader.bind(INITIAL_REFERENCE,INITIAL_REFERENCE_SHA);reference,rs=reader.checkpoint(INITIAL_REFERENCE)
    old.properties_equal(pa,reference)
    p=state['parameters']
    require(p==state['nominal_parameters']==rs['parameters'] and state['nominal_parameter_provenance']==rs['parameter_provenance'],
        'nominal physical parameters/provenance differ')
    motion=dict(classification=MOTION_CLASS,ramp_duration_s=.02,
        velocity_factor='3*s^2-2*s^3 for 0<=s=t/T<=1; 1 thereafter',
        integrated_factor='T*(s^3-s^4/2) until T; t-T/2 thereafter',nominal_omega_rad_s=p['omega'],
        nominal_traverse_m_s=p['traverse'],nominal_feed_m_s=p['feed'],shoulder_height_m=p['layer'],
        step_path='LINEAR_WORLD_ORIGIN_Z_GEODESIC_ROTATION_LINEAR_ACTUATOR',
        interpretation='NEW_STARTUP_BOUNDARY_HISTORY_ASSUMPTION_NOT_PHYSICAL_VALIDATION')
    provenance=dict(scope=SCOPE,nominal_parameter_reference=rs['parameter_provenance'],boundary_condition_history=motion,
        changed_boundary_history_fields=['omega','traverse','feed'],nominal_values_are_not_instantaneous_startup_velocities=True,
        physical_qualification='NOT_RUN_NEW_STARTUP_BOUNDARY_ASSUMPTION',old_state_read=False,parent_checkpoint=None,
        mass_scale=1.,time_scale=1.,modulus_scale=1.)
    initialization=dict(rs['independent_domain_initialization'],kind=SCOPE,id_namespace='independent_ramped_hot_8mm',
        physical_qualification='NOT_RUN_NEW_STARTUP_BOUNDARY_ASSUMPTION',factory_scope=rs['model_scope'],motion_protocol=motion)
    require(state['ramp_motion']==motion and state['parameter_provenance']==provenance and
        state['independent_domain_initialization']==initialization,'new boundary provenance mismatch')
    require(state['time_s']==0. and state['step']==0 and state['classification']==CLASS and state['model_scope']==SCOPE and
        state['independent_initialization'] is True and state['accepted'] is False and state['restart_allowed'] is False and
        state['diagnostic_only'] is True and state['accepted_time_advanced_s']==0. and state['ramped_boundary_history']==[] and
        state['ramped_physical_qualification']=='NOT_RUN','initial diagnostic scope mismatch')
    for k in ('spatial_model_transition','branch_origin_sha256','expanded_initial_qualification'):
        require(k not in state,'old process ancestry injected')
    for k in ('contact_history','free','interval_energy','bore_observation','crossing_ids','crossing_mass_kg','operator_variant'):
        require(clean(state[k])==clean(rs[k]),'unloaded initial invariant differs: '+k)
    expected_sources=dict(rs['source_sha256'])
    for name in ('ramped_hot_motion.py','ramped_hot_trial.py'):expected_sources['scripts/'+name]=frozen['scripts/'+name]
    require(state['source_sha256']==state['ramp_component_source_sha256']==expected_sources,'initial model source mismatch')
    reader.sources(expected_sources);old.validate_interval_energy_state(state['interval_energy']);finite_properties(pa)
    rank=reader.read(folder/'initial_elastic_rank.json');require(rank==meta['initial_elastic_rank'],'initial rank journal mismatch')
    require(rank['status']=='PASS_INITIAL_FREE_ELASTIC_RANK_NOT_FUTURE_FLOW_STABILITY' and
        rank['raw_visible_graph_identical'] and rank['free_displacement_dof']==1961 and
        rank['relative_min_eigenvalue']>rank['rank_relative_threshold']==1e-12 and
        not any(rank[k] for k in ('mass_shift','pseudoinverse','artificial_axis_constraints','matrix_symmetrized')),
        'actual t0 elastic-rank qualification failed')
    arrays=reader.npz(folder/'initial_elastic_rank_arrays.npz')
    old_arrays=reader.npz(OLD_RANK.parent/'initial_elastic_rank_arrays.npz')
    require(set(arrays)==set(old_arrays) and all(np.array_equal(arrays[k],old_arrays[k]) for k in arrays),
        'same t0 physical geometry/stiffness/eigen evidence differs from frozen full-rank reference')
    require(arrays['accepted'].item() is False and arrays['restart_allowed'].item() is False,'rank arrays were published')
    require(meta['static_adaptation']['status']=='PASS_STATIC_RAMPED_BOUNDARY_ONLY_EVALUATE_AND_TANGENT_ADAPTATION' and
        meta['static_adaptation']['runtime_monkeypatch'] is False,'registered boundary-only adaptation audit failed')
    require(all(np.isfinite(v) and v<=1e-10 for v in meta['scalar_batch_differences'].values()),'initial scalar/batch map differs')
    return pa,state,path,frozen,dict(initial_checkpoint=rel(path),initial_sha256=digest(path),
        source_count=len(frozen),cost_wall_seconds=meta['wall_seconds'],full_map_wall_s=meta['full_map_wall_s'],
        tangent_wall_s=meta['tangent_wall_s'],initial_elastic_rank=rank,initial_reference_is_geometry_only_not_process_parent=True)


def audit_stage(reader,base,chunk,stage,initial,parent,frozen,completed):
    folder=base/f'chunk_{chunk:03d}'/stage;pa,state,parent_path=parent
    protocol,meta,supervision=protocol_audit(reader,base,folder,stage,chunk,frozen,parent_path,completed)
    require(abs(state['time_s']-(chunk-1)*.001)<1e-15 and state['step']==(chunk-1)*COUNTS[stage],
        'path was restarted from wrong branch/time')
    context=dict(initial=initial,parent=state,p=state['parameters'],dt=DTS[stage],prefix=state['interval_energy']['records'])
    rows=reader.read(folder/'candidate_history.json');require(len(rows)==meta['candidate_steps']<=COUNTS[stage],'candidate count mismatch')
    direct=digest(parent_path);cases=[]
    for number,row in enumerate(rows,1):
        require(row['step']==number and row['parent_sha256']==direct,'candidate immediate ancestry broken')
        item=candidate_audit(reader,stage,folder,row,pa,state,context,frozen);cases.append(item)
        pa,state=item['pa'],item['state'];direct=row['sha256']
    full=meta['status']==PASS_STAGE
    require(full==(supervision['exit_code']==0),'status/exit inconsistency')
    if full:
        require(len(rows)==COUNTS[stage] and abs(state['time_s']-chunk*.001)<1e-15 and
            meta['candidate_end_checkpoint']==rel(folder/rows[-1]['checkpoint']) and meta['candidate_end_sha256']==direct and
            clean(state['interval_energy'])==meta['interval_energy'] and meta['crossing_ids']==state['crossing_ids'],
            'complete branch endpoint/energy mismatch')
    else:require(meta['status']=='FAIL_RAMP_STARTUP_DIAGNOSTIC_NOT_ACCEPTED' and 'error' in meta,'unknown nonterminal physical status')
    partial=[rel(p) for p in sorted(folder.glob('step_*')) if p.is_dir() and int(p.name.split('_')[-1])>len(rows)]
    require(len(partial)<=1,'unexpected extra rejected candidate attempts')
    result=dict(stage=stage,status=meta['status'],exit_code=supervision['exit_code'],physical_stage_completed=full,
        parent_checkpoint=rel(parent_path),parent_sha256=digest(parent_path),candidate_count=len(rows),
        diagnostic_endpoint=old.state_stats(pa,state),complete_interval_energy=clean(state['interval_energy']),
        full_prefix_preserved=True,wall_seconds=meta['wall_seconds'],supervised_wall_seconds=supervision['wall_seconds'],
        forecast_s=protocol['forecast_s'],candidate_summaries=[v['summary'] for v in cases],
        failed_partial_directories=partial,error=meta.get('error'),accepted_steps=0,accepted_time_advanced_s=0.)
    terminal=(pa,state,folder/rows[-1]['checkpoint']) if cases else parent
    return cases,result,terminal,meta


def replay_screens(reader,base,chunk,coarse,fine,p,fine_meta):
    folder=base/f'chunk_{chunk:03d}/fine';paths=sorted(folder.glob('pair_*_screen.json'));reports=[]
    require(len(paths)==len(fine)//2,'complete fine pair lacks its registered screen')
    saved=[]
    for number,path in enumerate(paths,1):
        s=reader.read(path);saved.append(s);a,b=fine[2*number-2:2*number];c=coarse[number-1]
        require(s['pair']==number and s['time_s']==b['state']['time_s'] and abs(s['time_s']-c['state']['time_s'])<1e-15 and
            s['coarse_checkpoint']==rel(base/f'chunk_{chunk:03d}/coarse'/c['row']['checkpoint']) and
            s['coarse_sha256']==c['row']['sha256'] and s['fine_checkpoint_sha256']==[a['row']['sha256'],b['row']['sha256']] and
            s['comparison_scope']=='COMMON_INITIAL_CONTINUOUS_TWO_RESOLUTION' and s['accepted'] is False and
            s['restart_allowed'] is False,'shared endpoint binding mismatch')
        glob,local=old.compare_arrays(b['pa'],c['pa'],sorted(int(v) for v in b['pa'].particle_id),p)
        same=b['state']['crossing_ids']==c['state']['crossing_ids'];passed=bool(glob['pass'] and local['passed'] and same)
        require(clean(glob)==s['global_screen'] and clean(local)==s['local_screen'] and
            same==s['crossing_ids_match'] and passed==s['passed'],'all-particle temporal comparison does not replay')
        reports.append(dict(pair=number,time_s=s['time_s'],passed=passed,global_screen=glob,
            local_particle_count=len(local['rows']),local_max_stress_relative=max(v['stress_relative'] for v in local['rows']),
            local_max_plastic_relative=max(v['plastic_relative'] for v in local['rows']),
            crossing_ids_match=same,screen_sha256=digest(path)))
        require(number==len(paths) or passed,'continued past a rejected temporal screen')
    require(fine_meta['pair_screens']==saved,'manifest screen ledger differs')
    if fine_meta['status']==PASS_STAGE:require(len(reports)==2 and all(v['passed'] for v in reports) and fine_meta['temporal_screens_all_pass'],
        'complete stage lacks both passing comparisons')
    return reports


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--base',required=True)
    parser.add_argument('--chunk',type=int,required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    base=(ROOT/args.base).resolve();output=(ROOT/args.output).resolve()
    require(1<=args.chunk<=4 and base.is_relative_to(ROOT/'outputs') and base.is_dir() and
        output.is_relative_to(ROOT/'outputs') and not output.exists() and
        not any(output.is_relative_to(base/f'chunk_{k:03d}'/s) for k in range(1,5) for s in COUNTS),
        'fresh readback folder outside physical stage directories required')
    output.mkdir(parents=True,exist_ok=False);reader=Reader();started=time.perf_counter()
    report=dict(status='READBACK_RUNNING',base=rel(base),chunk=args.chunk,physical_solver_evaluations=0,
        nnps_evaluations=0,material_evaluations=0,contact_evaluations=0,accepted_steps=0,accepted_time_advanced_s=0.,
        accepted=False,restart_allowed=False,production_authorized=False,command=[sys.executable,*sys.argv])
    try:
        reader.bind(Path(__file__));reader.bind(Path(prior.__file__));reader.bind(Path(old.__file__))
        completed=[];initial,state,path,frozen,cost=initial_audit(reader,base,completed)
        parents={s:(initial,state,path) for s in COUNTS};history={s:[] for s in COUNTS};chunks=[];allpass=True
        first_loaded=None;loaded_pass=False
        for chunk in range(1,args.chunk+1):
            cases={};stages={};metas={}
            for stage in COUNTS:
                folder=base/f'chunk_{chunk:03d}'/stage
                if not folder.exists():
                    require(stage=='fine' and not stages['coarse']['physical_stage_completed'],'missing scheduled stage')
                    cases[stage]=[];stages[stage]=dict(status='NOT_RUN_AFTER_PHYSICAL_FAILURE',physical_stage_completed=False);continue
                cases[stage],stages[stage],parents[stage],metas[stage]=audit_stage(reader,base,chunk,stage,initial,parents[stage],frozen,completed)
                history[stage].extend(cases[stage])
            screens=replay_screens(reader,base,chunk,cases['coarse'],cases['fine'],state['parameters'],metas['fine']) if 'fine' in metas else []
            passed=all(v['physical_stage_completed'] for v in stages.values()) and len(screens)==2 and all(v['passed'] for v in screens)
            allpass=allpass and passed
            for screen in screens:
                loaded={}
                for stage in COUNTS:
                    eligible=[v for v in history[stage] if v['state']['time_s']<=screen['time_s']+1e-15]
                    loaded[stage]=any(v['response']['nonbase_loaded_ids'] for v in eligible) and any(v['response']['new_plastic_ids'] for v in eligible)
                if screen['passed'] and all(loaded.values()) and first_loaded is None:first_loaded=screen['time_s']
                elif screen['passed'] and first_loaded is not None and screen['time_s']>=first_loaded+.0005-1e-15:loaded_pass=True
            chunks.append(dict(chunk=chunk,stages=stages,screens=screens,temporal_screens_all_pass=passed))
            require(chunk==args.chunk or passed,'attempted later chunk after failure')
        reader.finish()
        report.update(status=PASS_READBACK,initial_cost=cost,chunks=chunks,temporal_screens_all_pass=allpass,
            continuation_within_registered_experiment_allowed=allpass and args.chunk<4,
            loaded_feasibility_screen_passed=bool(allpass and loaded_pass),first_both_paths_loaded_and_plastic_shared_time_s=first_loaded,
            physical_interpretation=('LOADED_POST_CONTACT_TEMPORAL_SCREEN_PASSED_NOT_PROCESS_ACCEPTANCE' if allpass and loaded_pass else
                'PRE_TOOL_CONTACT_ELASTIC_STARTUP' if not any(v['response']['nonbase_loaded_ids'] for h in history.values() for v in h) else
                'TOOL_CONTACT_OBSERVED_LOADED_POST_CONTACT_SCREEN_NOT_ESTABLISHED'),
            total_supervised_wall_seconds=sum(s['wall_seconds'] for _,_,s in completed),
            source_entry_checks_including_repeats=reader.source_entries,evidence_sha256=reader.evidence,
            limitations=['Readback PASS verifies saved evidence, including faithfully retained failures; it is not process acceptance.',
                'Both branches remain unaccepted and cannot be used as production or accepted-process restarts.',
                'Later chunks have separate immediate parents and a common independent t=0 origin.',
                'The 20 ms boundary history is a new engineering assumption; low tool load does not mean zero actuator work.',
                'Contact-history, whole-history, spatial and experimental convergence are not established.']);code=0
    except Exception:
        report.update(status='FAIL_RAMP_CHUNK_SAVED_READBACK',error=traceback.format_exc(),evidence_sha256=reader.evidence,
            temporal_screens_all_pass=False,continuation_within_registered_experiment_allowed=False);code=2
    report.update(exit_code=code,wall_seconds=time.perf_counter()-started)
    (output/'summary.json').write_text(json.dumps(clean(report),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report.get(k) for k in ('status','chunk','temporal_screens_all_pass','loaded_feasibility_screen_passed','physical_interpretation','wall_seconds','error')},indent=2))
    return code


if __name__=='__main__':raise SystemExit(main())


"""Static bounded adaptation of the two-birth diagnostic; no source is patched."""
import copy
import numpy as np
import two_insertion_diagnostic as old
import fixed_inventory_source_work as fixed_work
from afsd_backend.pysph.interval_acceptance import initial_bore_crossing_state

ROOT,require,write,read,sha,rel,verify=old.ROOT,old.require,old.write,old.read,old.sha,old.rel,old.verify
H=old.DT/2
STUDY='TWO_BIRTH_TIME_AND_PRELOADED_REFERENCE_V1'
CONTRACT=ROOT/'docs/INLET_REFERENCE_BENCHMARK_CONTRACT_20261005.md'


def new_layer(p,height):
    xy=np.arange(-2,3)*p['dx'];points=np.array([[x,y,height] for x in xy for y in xy])
    pa=old.make_array(points,np.ones(25,dtype=int),p);pa.particle_id[:]=0
    pa.u[:]=0.;pa.v[:]=0.;pa.w[:]=-p['feed'];pa.temperature[:]=723.15
    return pa


def initialize(branch,frozen):
    require(branch in ('inlet_fine','preloaded_fine'),'registered branch required')
    pa,state=old.initial_state(frozen);p=state['parameters']
    state.update(study_scope=STUDY,study_branch=branch,study_contract_sha256=sha(CONTRACT))
    if branch=='preloaded_fine':
        for k,height in enumerate((.00325,.00375)):
            layer=new_layer(p,height);layer.particle_id[:]=np.arange(919+25*k,944+25*k)
            pa.append_parray(layer)
        pa.set_output_arrays(list(pa.properties));free=np.ones((968,3),bool)
        free[:918]=old.inlet.process_state_for(pa.extract_particles(np.arange(918,dtype=np.uint32)),state['inlet_acceptance'])['free']
        active=(pa.component_id==1)&(pa.z>=.0025);pa.actuator[:]=active.astype(int);free[active,2]=False
        physical=dict(time_s=0.,free=free,contact_history=np.zeros((968,8,3)),
            bore_observation=initial_bore_crossing_state(old.positions(pa),pa.particle_id,pa.component_id,pa.m,0.,old.core.pose(p,0.),p['bore_width']))
        registry=old.inlet.inventory.initialize_registry(pa,rel(CONTRACT),acknowledgement=old.inlet.inventory.ACKNOWLEDGEMENT)
        v2=old.inlet.initialize_acceptance_state(pa,registry,physical,
            dict(frame='WORLD_Z',plane_z_m=.0025,unit='m',source=rel(CONTRACT)),rel(CONTRACT))
        state.update(inlet_acceptance=v2,bore_observation=v2['bore_observation'],initialization_kind='PRELOADED_50_EXTRA_INITIAL_PARTICLES_NO_BIRTH',
            preloaded_initial_ids=list(range(919,969)),original918_array_audit=state.pop('initial_array_audit'))
        position=old.positions(pa);extra=pa.particle_id>=919
        distances=np.sum((position[extra,None]-position[None,:])**2,axis=2)
        for row,col in enumerate(np.flatnonzero(extra)):distances[row,col]=np.inf
        minimum=float(np.sqrt(distances.min()));margin=128*np.finfo(float).eps*p['substrate_width']
        require(minimum>=p['dx']-margin,'preloaded initial nominal cells overlap')
        geometry=old.contact_factory(old.BatchedChainRuleContact,p).geometry(position,old.core.pose(p,0.))
        require(not np.any(geometry['potential_J'][extra]) and not np.any(pa.S.reshape(-1,3,3)[extra]) and
            np.all(pa.rho[extra]==p['rho0']) and free.sum()==2061,'preloaded initial contact/material/axes differ')
        state['preloaded_initial_geometry']=dict(minimum_extra_center_distance_m=minimum,roundoff_margin_m=margin,
            extra_material_and_contact_storage_J=0.,particle_count=968,free_dof=2061,eligible_ids=v2['bore_observation']['eligible_ids'])
    else:state['initialization_kind']='ORIGINAL_918_WITH_REGISTERED_BIRTHS_AT_0_AND_2H'
    state['study_initial_pa_sha256']=old.inlet.inventory._array_digest(pa)
    state['study_initial_mass_kg']=float(pa.m.sum())
    state['study_initial_sensible_heat_J']=float(np.sum(pa.m*p['cp']*(pa.temperature-p['room'])))
    return pa,state


def prepare(parent,state,step):
    require(1<=step<=4 and state['step']==step-1 and state['time_s']==(step-1)*H,'unregistered step/time')
    p=state['parameters'];v2=state['inlet_acceptance'];event=None
    if state['study_branch']=='inlet_fine' and step in (1,3):
        index=1 if step==1 else 2;incoming=new_layer(p,(.00325,.00375)[index-1])
        fields={key:dict(unit=old.inlet.inventory.REQUIRED.get(key,'1'),source=rel(CONTRACT),
            status='INITIAL_CONDITION' if key in old.inlet.inventory.REQUIRED else 'UNRESOLVED') for key in incoming.properties}
        spec=dict(event_id='diagnostic-inventory-birth-'+str(index),time_s=state['time_s'],source=rel(CONTRACT),fields=fields,
            heat_capacity=dict(value=p['cp'],unit='J/(kg K)',source=rel(CONTRACT),status='DECLARED_ASSUMPTION'),
            reference_temperature=dict(value=p['room'],unit='K',source=rel(CONTRACT),status='DECLARED_ASSUMPTION'))
        pa,pending,event=old.inlet.propose_insertion(parent,v2,incoming,spec,old.core.pose(p,state['time_s']))
    else:pa,pending=old.inlet.propose_mechanical_step(parent,v2)
    aligned=old.inlet.process_state_for(pa,pending);prior=old.inlet.process_state_for(parent,v2)
    inner=copy.deepcopy(state);inner.update(aligned,inlet_acceptance=pending,
        bore_observation=old.inlet._legacy_observer(pending['bore_observation']),observer_adapter='INTERNAL_FIXED_N_LEGACY_ONLY')
    require(not pending['control_events'][-1]['released_ids'],'release outside this narrow registered scope')
    expected=(943 if step<=2 else 968) if state['study_branch']=='inlet_fine' else 968
    require(len(pa.x)==expected and aligned['free'].sum()==1961+2*(expected-918),'unexpected inventory/control degrees')
    contact=old.contact_factory(old.BatchedChainRuleContact,p);t=state['time_s']
    start_geometry=contact.geometry(old.positions(parent),old.core.pose(p,t))
    storage=np.column_stack((start_geometry['potential_J'].sum(axis=1),.5*contact.kt*np.sum(prior['contact_history']**2,axis=(1,2))))
    geometric=dict(status='PASS',kind='UNCHANGED_PREVIOUSLY_SEALED_INVENTORY',newborn_ids=[],overlap_test='NO_NEW_PARTICLES')
    if event is not None:
        born=np.isin(pa.particle_id,event['allocated_particle_ids']);pos=old.positions(pa)
        squared=np.sum((pos[born,None]-pos[None,:])**2,axis=2)
        for k,j in enumerate(np.flatnonzero(born)):squared[k,j]=np.inf
        distance=float(np.sqrt(squared.min()));margin=128*np.finfo(float).eps*p['substrate_width']
        require(distance>=p['dx']-margin,'newborn initial cell overlap')
        geometry=contact.geometry(pos,old.core.pose(p,t));require(not np.any(geometry['potential_J'][born]),'newborn contact storage nonzero')
        geometric.update(kind='REAL_BIRTH_INITIAL_CELL_AND_CONTACT_CHECK',newborn_ids=event['allocated_particle_ids'],
            minimum_center_distance_m=distance,roundoff_margin_m=margin,initial_contact_storage_J=0.)
    return pa,inner,event,storage,geometric


def qualified_step(parent,trial,result,event,parent_storage,geometric,directory,frozen):
    pa,state,fresh=old.fixed_inventory_copy(trial,result);p=trial.p;t=trial.state['time_s'];dt=result['dt_s']
    require(dt==H and state['interval_energy']['records'][:-1]==trial.state['interval_energy']['records'],'step/prefix changed')
    require(result['graph_fingerprint']==fresh['graph_fingerprint'] and np.array_equal(result['i'],fresh['i']) and
        np.array_equal(result['j'],fresh['j']),'fresh actual graph differs')
    old.save_map(directory,trial,fresh);energy=old.core.energy_diagnostics(trial,fresh)
    scalar=old.InletTrial(trial.pa,trial.state,p);scalar.contact=old.contact_factory(old.ChainRuleDiscreteContact,p)
    sr=scalar.evaluate(fresh['displacement'],dt)
    difference=old.batch.candidate_equivalence(dict(sr,energy_comparison=old.core.energy_diagnostics(scalar,sr)),
        dict(fresh,energy_comparison=energy))
    require(sr['graph_fingerprint']==fresh['graph_fingerprint'] and np.array_equal(sr['g'],fresh['g']) and
        np.array_equal(sr['thermal']['pair_conductance_W_K'],fresh['thermal']['pair_conductance_W_K']),'scalar graph/heat differs')
    write(directory/'scalar_equivalence.json',dict(status='PASS',differences=difference))
    geometry=old.post.geometry_audit(trial,fresh);write(directory/'sampled_geometry.json',geometry)
    require(not geometry['transient_contacts_missed_by_endpoints'] and set(geometry['new_non_bore_features'])<=set(('SHOULDER','TOP')),
        'unqualified finite contact event')
    pose0,pose1=old.core.pose(p,t),old.core.pose(p,t+dt)
    shell=old.audit_contact_shell_finite_enclosure(trial.contact,trial.x0,fresh['position'],pose0,pose1)
    write(directory/'continuous_shell.json',shell);require(shell['status']==old.SHELL_CERTIFIED,'continuous shell unresolved')
    path=old.audit_visible_support_path(trial.x0,fresh['position'],fresh['i'],fresh['j'],pose0,pose1,p,trial.pa.particle_id,
        maximum_depth=32,maximum_queries=256)
    write(directory/'visible_path.json',path);require(path['status']==old.PATH_CERTIFIED,'retained support path unresolved')
    rank=old.reference_rank(trial,fresh,directory)
    port=dict(contract=old.source.PORT_CONTRACT,surface='EXTERIOR_BIRTH_BOUNDARY_NOT_INTERNAL_CONTROL_PLANE',
        internal_interface_work='NOT_ADDED_INTERNAL_PAIR_WORK',traction_Pa=[0.,0.,0.],source=rel(CONTRACT))
    arguments=dict(parameters=p,free=trial.state['free'],result=fresh,mechanical_energy=energy,
        parent_contact_storage_J=parent_storage,external_port=port)
    if event is None:
        calculated=fixed_work.evaluate(parent,trial.pa,pa,**arguments);recorder=fixed_work.record
    else:
        calculated=old.source.evaluate_source_work(parent,trial.pa,pa,event,**arguments);recorder=old.source.source_work_record
    write(directory/'source_work.json',calculated);work=recorder(calculated,directory/'source_work.json');write(directory/'source_work_v2.json',work)
    pending=trial.state['inlet_acceptance'];history=fresh['contact']['history'];free=trial.state['free']
    binding=old.inlet.mechanical_binding(pa,pending,history,free,state['time_s'])
    require(float(np.max(np.abs(fresh['force'][free])))<=fresh['force_limit_N'] and fresh['scaled_residual']<=1,'actual force failed')
    details=dict(force=dict(actual_free_force_N=float(np.max(np.abs(fresh['force'][free]))),limit_N=fresh['force_limit_N']),
        heat=dict(ledger=fresh['thermal']['ledger'],feedback_residual_K=fresh['thermal_residual_K']),
        energy=dict(actual=energy,interval=state['interval_energy']),geometry_overlap=geometric,
        NNPS=dict(graph_fingerprint=fresh['graph_fingerprint'],edge_count=len(fresh['i']),visibility=fresh['visibility_graph']['diagnostics']),
        visible_support=dict(status=path['status'],sha256=sha(directory/'visible_path.json')),gradient_rank=rank,
        contact=dict(status=shell['status'],sha256=sha(directory/'continuous_shell.json'),sweep=fresh['contact']['sweep'],
            scalar_sha256=sha(directory/'scalar_equivalence.json')))
    gates={}
    for key,value in details.items():
        path=directory/('gate_'+key+'.json');write(path,dict(status='PASS',binding=binding,scope=old.SCOPE,
            source_sha256=frozen,source_work_report_sha256=sha(directory/'source_work.json'),evidence=value))
        gates[key]=dict(status='PASS',artifact=rel(path),evidence_sha256=sha(path))
    receipt=dict(scope=old.SCOPE,step_id=state['study_branch']+'-step-'+str(state['step']),binding=binding,gates=gates,
        source_work_sha256=old.inlet._sha(work),source_kind='FIXED_INVENTORY_NO_BIRTH' if event is None else 'REAL_BIRTH',
        joint_inlet_benchmark='NOT_RUN',accepted=False,restart_allowed=False,
        audit_artifacts_sha256={rel(path):sha(path) for path in directory.rglob('*') if path.is_file()})
    write(directory/'mechanical_receipt.json',receipt)
    def verifier(given):
        require(given==read(directory/'mechanical_receipt.json'),'actual receipt changed');verify(frozen);verify(given['audit_artifacts_sha256'])
        for key in old.inlet.GATES:
            item=given['gates'][key];path=ROOT/item['artifact'];row=read(path)
            require(sha(path)==item['evidence_sha256'] and row['status']=='PASS' and row['binding']==binding and
                row['source_sha256']==frozen and row['source_work_report_sha256']==sha(directory/'source_work.json'),'gate binding changed')
        require(recorder(calculated,directory/'source_work.json')==work and read(directory/'source_work_v2.json')==work,'source work artifact changed')
        require(binding==old.inlet.mechanical_binding(pa,pending,history,free,state['time_s']),'actual PA/history changed')
        return True
    sealed=old.inlet.seal_mechanical_acceptance(pa,pending,history,free,state['time_s'],receipt,work,verifier,lambda time:old.core.pose(p,time))
    require(old.inlet._legacy_observer(sealed['bore_observation'])==old.clean(state['bore_observation']),'V2 and actual fixed-step observers differ')
    state.update(inlet_acceptance=sealed,bore_observation=sealed['bore_observation'],accepted=False,restart_allowed=False,process_acceptance=False,
        observer_adapter='BIRTH_AWARE_V2_PUBLISHED_INTERNAL_LEGACY_FIXED_N_ONLY',mechanical_receipt_sha256=sha(directory/'mechanical_receipt.json'),
        source_work_report_sha256=sha(directory/'source_work.json'))
    return pa,state,fresh,calculated

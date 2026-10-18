"""Saved-only audits and registered comparisons; never call a physical map."""
import copy
import json
import math
import time
import traceback
import numpy as np
import inlet_reference_benchmark as model
from afsd_backend.pysph.interval_acceptance import require_interval_energy
from readback_two_insertion_diagnostic import near

ROOT,require,read,write,sha,rel,verify=model.ROOT,model.require,model.read,model.write,model.sha,model.rel,model.verify
OLD_BASE=ROOT/'outputs/two_insertion_diagnostic_20261005/v1/diagnostic'


def checkpoint(directory,step):
    folder=directory/('step_'+str(step).zfill(3));path=folder/'DIAGNOSTIC_NOT_RESTARTABLE.npz'
    saved=model.old.domain.load(str(path));pa=saved['arrays']['solid'];state=saved['solver_data']
    require(state['accepted'] is False and state['restart_allowed'] is False,'comparison input became process accepted')
    return pa,state,folder


def matching_rows(pa,ids):
    index={int(pid):k for k,pid in enumerate(pa.particle_id)}
    require(set(ids)<=set(index),'comparison loses registered IDs')
    return pa.extract_particles(np.array([index[pid] for pid in ids],dtype=np.uint32))


def field_screen(fine,reference,graph_fine,graph_reference,*,original_only):
    ids=sorted(set(fine.particle_id.astype(int)) & set(reference.particle_id.astype(int)))
    if original_only:ids=[pid for pid in ids if pid<=918]
    else:require(set(ids)==set(fine.particle_id)==set(reference.particle_id),'temporal inventories differ')
    f,r=matching_rows(fine,ids),matching_rows(reference,ids)
    global_screen=model.old.post.compare(f,r,dict(dx=.0005))
    seeds=set(range(769,919));selected=set(seeds)
    for graph in (graph_fine,graph_reference):
        particle=graph['particle_id'];left,right=particle[graph['i']],particle[graph['j']]
        adjacent=np.isin(left,list(seeds))|np.isin(right,list(seeds))
        selected.update(int(value) for value in np.r_[left[adjacent],right[adjacent]])
    selected=sorted(selected&set(ids));index=np.array([ids.index(pid) for pid in selected])
    dx=np.linalg.norm(model.old.positions(f)-model.old.positions(r),axis=1)
    T=np.abs(f.temperature-r.temperature);sf=f.S.reshape(-1,9);sr=r.S.reshape(-1,9)
    Sn=np.linalg.norm(sf,axis=1);Rn=np.linalg.norm(sr,axis=1);Se=np.linalg.norm(sf-sr,axis=1)
    er=np.abs(f.ep-r.ep)/np.maximum(np.abs(f.ep),1e-5);srel=Se/np.maximum(Sn,1e6)
    passed=bool(global_screen['pass'] and np.all(dx[index]<=5e-6) and np.all(T[index]<=.5) and
        np.all(srel[index]<=.1) and np.all(er[index]<=.2))
    row=dict(passed=passed,global_screen=global_screen,comparison_particle_ids=ids,local_particle_ids=selected,
        local_position_max_m=float(dx[index].max()),local_temperature_max_K=float(T[index].max()),
        local_stress_relative_max=float(srel[index].max()),local_plastic_relative_max=float(er[index].max()),
        actual_stress_absolute_max_difference_Pa=float(Se.max()),actual_fine_stress_frobenius_norm_Pa=float(np.linalg.norm(sf)),
        actual_reference_stress_frobenius_norm_Pa=float(np.linalg.norm(sr)),
        stress_denominator_floor_Pa=1e6,plastic_denominator_floor=1e-5,
        stress_floor_active_particle_count=int(np.sum(Sn<1e6)),local_stress_floor_active_count=int(np.sum(Sn[index]<1e6)),
        local_rows=[dict(particle_id=pid,position_error_m=float(dx[k]),temperature_error_K=float(T[k]),
            absolute_stress_difference_Pa=float(Se[k]),fine_stress_norm_Pa=float(Sn[k]),reference_stress_norm_Pa=float(Rn[k]),
            stress_relative_with_floor=float(srel[k])) for pid,k in zip(selected,index)])
    return row


def temporal_screen(base):
    rows=[];last_fine=last_coarse=None
    for half,coarse in ((2,1),(4,2)):
        f,fs,fd=checkpoint(base/'inlet_fine',half);c,cs,cd=checkpoint(OLD_BASE,coarse)
        require(fs['time_s']==cs['time_s'] and fs['time_s']==half*model.H,'time screen endpoint mismatch')
        require([(e['time_s'],e['allocated_particle_ids']) for e in fs['inlet_acceptance']['source_events']]==
            [(e['time_s'],e['allocated_particle_ids']) for e in cs['inlet_acceptance']['source_events']], 'birth schedule/IDs changed with dt')
        gf,gc=dict(np.load(fd/'visible_graph.npz')),dict(np.load(cd/'visible_graph.npz'))
        screen=field_screen(f,c,gf,gc,original_only=False);screen['time_s']=fs['time_s']
        fine_prefix,coarse_prefix=fs['interval_energy'],cs['interval_energy']
        allowance=512*np.finfo(float).eps*(fine_prefix['cumulative_absolute_work_scale_J']+coarse_prefix['cumulative_absolute_work_scale_J'])
        prefix_pass=fine_prefix['cumulative_absolute_step_residual_J']<=coarse_prefix['cumulative_absolute_step_residual_J']+allowance
        screen.update(fine_prefix_absolute_raw_J=fine_prefix['cumulative_absolute_step_residual_J'],
            coarse_prefix_absolute_raw_J=coarse_prefix['cumulative_absolute_step_residual_J'],
            no_worse_raw_roundoff_J=allowance,no_worse_raw_pass=prefix_pass,passed=screen['passed'] and prefix_pass)
        rows.append(screen)
        last_fine,last_coarse=fs,cs
    fine=last_fine['interval_energy'];coarse=last_coarse['interval_energy']
    margin=512*np.finfo(float).eps*(fine['cumulative_absolute_work_scale_J']+coarse['cumulative_absolute_work_scale_J'])
    raw_pass=fine['cumulative_absolute_step_residual_J']<=coarse['cumulative_absolute_step_residual_J']+margin
    return dict(passed=all(row['passed'] for row in rows) and raw_pass,rows=rows,
        fine_absolute_raw_J=fine['cumulative_absolute_step_residual_J'],coarse_absolute_raw_J=coarse['cumulative_absolute_step_residual_J'],
        no_worse_raw_roundoff_J=margin,no_worse_raw_pass=raw_pass,
        original_force_heat_energy_gates_not_relaxed=True,physical_evaluations=0)


def audit_branch(base,branch):
    directory=base/branch;meta=read(directory/'manifest.json');external=read(directory/'supervision.json')
    require(meta['status']=='PASS_FOUR_STEP_REFERENCE_BRANCH_NOT_PROCESS_ACCEPTANCE' and meta['diagnostic_steps']==4 and
        external['exit_code']==0 and external['before_after_hashes_match'] and external['manifest_sha256']==sha(directory/'manifest.json') and
        external['log_sha256']==sha(directory/'run.log'),'completed actual branch receipt absent')
    verify(meta['source_sha256']);verify(meta['reference_sha256']);verify(meta['artifacts_sha256'])
    initial=ROOT/meta['initial_checkpoint'];require(sha(initial)==meta['initial_checkpoint_sha256'],'study initial changed')
    saved=model.old.domain.load(str(initial));parent=saved['arrays']['solid'];initialstate=saved['solver_data'];p=initialstate['parameters']
    v2=model.old.inlet.restore_acceptance_state(json.dumps(initialstate['inlet_acceptance']),parent)
    require(not v2['source_events'] and not v2['accepted_bindings'] and v2['mechanical_step_count']==0,'independent study initial has previous inlet history')
    require(len(parent.x)==(918 if branch=='inlet_fine' else 968),'initial inventory differs')
    prefix=copy.deepcopy(initialstate['interval_energy']);previous_sha=sha(initial);previous_contact=None;reports=[];rows=[]
    for step in range(1,5):
        actual,state,sd=checkpoint(directory,step);event=read(sd/'source_event.json');receipt=read(sd/'mechanical_receipt.json')
        require(state['study_scope']==model.STUDY and state['study_branch']==branch and state['time_s']==step*model.H and
            state['step']==step and state['parent_checkpoint_sha256']==previous_sha and state['source_sha256']==meta['source_sha256'],
            'branch provenance/time/parent differs')
        require(state['classification']==model.old.CLASSIFICATION and state['model_scope']==model.old.SCOPE and
            state['process_acceptance'] is False and state['parameters']==p and
            state['parameter_provenance']==initialstate['parameter_provenance'],'diagnostic model/parameters changed')
        if event is None:prepared,pending=model.old.inlet.propose_mechanical_step(parent,v2)
        else:
            height=.00325 if step==1 else .00375;incoming=model.new_layer(p,height)
            spec=dict(event_id=event['event_id'],time_s=event['time_s'],source=event['source'],fields=event['field_declarations'],
                heat_capacity=event['heat_capacity'],reference_temperature=event['reference_temperature'])
            prepared,pending,again=model.old.inlet.propose_insertion(parent,v2,incoming,spec,model.old.core.pose(p,(step-1)*model.H))
            require(again==event,'readback birth event differs')
        aligned=model.old.inlet.process_state_for(prepared,pending)
        f=dict(np.load(sd/'actual_fields.npz'));g=dict(np.load(sd/'visible_graph.npz'))
        require(np.array_equal(prepared.particle_id,f['particle_id']) and np.array_equal(actual.particle_id,prepared.particle_id) and
            np.array_equal(g['particle_id'],prepared.particle_id),'actual field inventory/row order differs')
        expected=prepared.extract_particles(np.arange(len(prepared.x),dtype=np.uint32))
        expected.x[:],expected.y[:],expected.z[:]=f['position'].T;expected.u[:],expected.v[:],expected.w[:]=(f['displacement']/model.H).T
        expected.rho[:]=f['density'];expected.S[:]=f['stress'].ravel();expected.ep[:]=f['ep'];expected.temperature[:]=f['temperature']
        expected.strain_rate[:]=f['rate'];expected.ep_rate[:]=(f['ep']-prepared.ep)/model.H
        expected.D[:]=(.5*(f['gradient']+f['gradient'].swapaxes(1,2))).ravel()
        expected.plastic_work[:]+=f['volume']*f['plastic_work_density'];expected.return_dissipation[:]+=f['volume']*f['return_dissipation_density']
        for name,prop in expected.properties.items():require(np.array_equal(prop.get_npy_array(),actual.properties[name].get_npy_array()),'raw checkpoint property differs: '+name)
        require(np.array_equal(f['position'],model.old.positions(prepared)+f['displacement']) and np.array_equal(f['free'],aligned['free']),
            'actual displacement/control differs')
        require(np.array_equal(state['contact_history'],f['contact_history']) and np.array_equal(state['free'],aligned['free']),
            'raw checkpoint physical history/free differs from actual map')
        journal=read(sd/'FULL_MAP_DIAGNOSTIC_NOT_ACCEPTED.json');energy_file=read(sd/'actual_energy.json');energy=energy_file['mechanical']
        contact={key:np.asarray(value) if isinstance(value,list) else value for key,value in journal['contact'].items()}
        require(journal['energy']==energy and journal['thermal']==energy_file['thermal']['ledger'],'actual energy ledgers differ')
        lookup={int(pid):i for i,pid in enumerate(prepared.particle_id)};oldrows=[lookup[int(pid)] for pid in parent.particle_id]
        initial_contact=np.column_stack((contact['normal_potential0_J'].sum(axis=1),contact['tangential_potential0_J'].sum(axis=1)))
        parent_storage=initial_contact[oldrows] if previous_contact is None else np.array([previous_contact[int(pid)] for pid in parent.particle_id])
        arguments=dict(parameters=p,free=aligned['free'],result=dict(displacement=f['displacement'],position=f['position'],
            force=f['force'],density=f['density'],stress=f['stress'],temperature=f['temperature'],volume=f['volume'],dt_s=model.H,
            material=(f['stress'],f['ep'],f['plastic_work_density'],f['return_dissipation_density'],f['rate']),contact=contact,thermal=energy_file['thermal']),
            mechanical_energy=energy,parent_contact_storage_J=parent_storage,external_port=read(sd/'source_work.json')['external_port'])
        if event is None:calculated=model.fixed_work.evaluate(parent,prepared,actual,**arguments);recorder=model.fixed_work.record
        else:calculated=model.old.source.evaluate_source_work(parent,prepared,actual,event,**arguments);recorder=model.old.source.source_work_record
        require(calculated==read(sd/'source_work.json'),'saved source-work recalculation differs')
        work=recorder(calculated,sd/'source_work.json');require(work==read(sd/'source_work_v2.json'),'V2 source-work artifact differs')
        prefix=require_interval_energy(prefix,energy);require(model.old.clean(state['interval_energy'])==prefix,'full energy prefix changed')
        binding=model.old.inlet.mechanical_binding(actual,pending,f['contact_history'],aligned['free'],step*model.H)
        require(receipt['binding']==binding and receipt['scope']==model.old.SCOPE and receipt['source_work_sha256']==model.old.inlet._sha(work),
            'actual mechanical/source receipt binding differs')
        verify(receipt['audit_artifacts_sha256'])
        for key in model.old.inlet.GATES:
            item=receipt['gates'][key];row=read(ROOT/item['artifact'])
            require(item['status']=='PASS' and sha(ROOT/item['artifact'])==item['evidence_sha256'] and row['binding']==binding and
                row['status']=='PASS' and row['source_sha256']==meta['source_sha256'],'actual gate not qualified')
        def verifier(given):
            require(given==read(sd/'mechanical_receipt.json'),'receipt changed');verify(given['audit_artifacts_sha256']);return True
        resealed=model.old.inlet.seal_mechanical_acceptance(actual,pending,f['contact_history'],aligned['free'],step*model.H,
            receipt,work,verifier,lambda time:model.old.core.pose(p,time))
        restored=model.old.inlet.restore_acceptance_state(json.dumps(state['inlet_acceptance']),actual)
        require(resealed==restored and state['bore_observation']==restored['bore_observation'],'V2 actual history/observer differs')
        require(restored['accepted_bindings'][:-1]==v2['accepted_bindings'] and restored['source_work_events'][:-1]==v2['source_work_events'],
            'prior accepted diagnostic history reset')
        require(float(np.max(np.abs(f['force'][aligned['free']])))<=journal['metrics']['force_limit_N'],'actual saved force fails gate')
        measured=prepared.m*p['cp']*(actual.temperature-prepared.temperature);heat=np.asarray(energy_file['thermal']['particle_energy_J']['total_heat_input_J'])
        require(np.max(np.abs(measured-heat))<=1e-12+1e-10*np.abs(heat).sum(),'actual saved per-particle heat fails')
        storage1=np.column_stack((contact['normal_potential1_J'].sum(axis=1),contact['tangential_potential1_J'].sum(axis=1)))
        previous_contact={int(pid):value for pid,value in zip(actual.particle_id,storage1)}
        reload=read(sd/'reload_and_sort.json');path=sd/'DIAGNOSTIC_NOT_RESTARTABLE.npz'
        require(reload['checkpoint_sha256']==sha(path) and reload['particle_ids_after']==actual.particle_id[::-1].tolist(),'registered reload/sort changed')
        parent=actual.extract_particles(np.arange(len(actual.x)-1,-1,-1,dtype=np.uint32))
        require(model.old.inlet.restore_acceptance_state(json.dumps(restored),parent)==restored,'sorted V2 registry differs')
        v2=restored;previous_sha=sha(path);reports.append(calculated)
        W=calculated['computed_terms']['mechanical_reaction_work']['value'];raw=energy['global_raw_residual_J']
        rows.append(dict(step=step,time_s=step*model.H,raw_mechanical_J=raw,reaction_work_J=W,
            signed_raw_over_reaction=None if W==0 else raw/W,ep_max=float(actual.ep.max()),
            temperature_max_change_K=float(np.max(np.abs(actual.temperature-prepared.temperature))),
            source_increment=calculated['source_increment'],particle_count=len(actual.x),
            source_kind='FIXED_INVENTORY' if event is None else 'BIRTH',released_ids=v2['control_events'][-1]['released_ids']))
    require(len(v2['source_events'])==(2 if branch=='inlet_fine' else 0) and v2['mechanical_step_count']==4,'final inlet histories differ')
    from run_two_insertion_diagnostic import cumulative_source_ledger
    cumulative=cumulative_source_ledger(reports,state)
    require(cumulative==meta['cumulative_source_ledger'],'independent cumulative source ledger differs')
    verify(meta['source_sha256']);verify(meta['reference_sha256']);verify(meta['artifacts_sha256'])
    W=math.fsum(row['reaction_work_J'] for row in rows);raw=prefix['cumulative_absolute_step_residual_J']
    return dict(status='PASS_SAVED_BRANCH_RECOMPUTATION',branch=branch,rows=rows,physical_evaluations=0,
        cumulative_source_ledger=cumulative,initial_sensible_heat_J=meta['initial_sensible_heat_J'],
        cumulative_absolute_raw_over_signed_reaction=None if W==0 else raw/W,
        final_checkpoint_sha256=previous_sha,worker_wall_seconds=meta['wall_seconds'])


def compare(base):
    directory=base/'saved_comparison';require(not directory.exists(),'comparison output exists');directory.mkdir()
    started=time.perf_counter()
    try:
        fine=audit_branch(base,'inlet_fine');reference=audit_branch(base,'preloaded_fine');temporal=temporal_screen(base)
        require(temporal['passed'],'temporal field/raw screen failed')
        fields=[]
        for step in (2,4):
            f,fs,fd=checkpoint(base/'inlet_fine',step);r,rs,rd=checkpoint(base/'preloaded_fine',step)
            require(fs['time_s']==rs['time_s'],'reference endpoint differs')
            row=field_screen(f,r,dict(np.load(fd/'visible_graph.npz')),dict(np.load(rd/'visible_graph.npz')),original_only=True)
            row['time_s']=fs['time_s'];fields.append(row)
        require(all(row['passed'] for row in fields),'preloaded finite-boundary field screen failed')
        source_heat=fine['cumulative_source_ledger']['source_sensible_heat_J']
        near(reference['initial_sensible_heat_J']-fine['initial_sensible_heat_J'],source_heat,
            'preloaded heat stock and real birth source differ',reference['initial_sensible_heat_J']+fine['initial_sensible_heat_J'])
        require(reference['cumulative_source_ledger']['source_sensible_heat_J']==0.,'preloaded reference invented a heat source')
        allrows=fine['rows']+reference['rows'];elastic=all(row['ep_max']==0. for row in allrows);isothermal=all(row['temperature_max_change_K']==0. for row in allrows)
        require(elastic and isothermal and all(not row['released_ids'] for row in allrows),'outside registered elastic/isothermal/no-release comparison scope')
        inputs={rel(path):sha(path) for branch in ('inlet_fine','preloaded_fine') for path in (base/branch).rglob('*') if path.is_file()}
        inputs[rel(__import__('pathlib').Path(__file__).resolve())]=sha(__import__('pathlib').Path(__file__).resolve())
        report=dict(status='PASS_ELASTIC_ISOTHERMAL_TWO_BIRTH_TIME_AND_PRELOADED_REFERENCE_NOT_PROCESS_ACCEPTANCE',
            accepted=False,restart_allowed=False,accepted_time_advanced_s=0.,source_and_input_sha256=inputs,
            temporal_screen=temporal,preloaded_field_screens=fields,inlet_fine=fine,preloaded_fine=reference,
            physical_evaluations=0,nnps_evaluations=0,root_evaluations=0,diagnostic_time_s=4*model.H,
            control_release='NOT_EXERCISED',plastic_response='NOT_EXERCISED',temperature_change='NOT_EXERCISED',
            continuous_inlet_process_qualification='NOT_GRANTED',wall_seconds=time.perf_counter()-started)
        verify(inputs);write(directory/'summary.json',report);print(json.dumps(dict(status=report['status'],output=rel(directory),physical_evaluations=0)))
    except Exception as error:
        write(directory/'FAILURE.json',dict(status='FAIL_SAVED_REFERENCE_COMPARISON',error=repr(error),traceback=traceback.format_exc(),physical_evaluations=0))
        raise

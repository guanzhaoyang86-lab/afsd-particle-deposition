"""Saved-only two-birth audit; no trial, NNPS, root, force, contact or heat evaluation."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import math
import sys
import time
import traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:sys.path.insert(0,str(ROOT/'src'))
from pysph.solver.utils import load
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.mechanics import positions
from afsd_backend.pysph.quasistatic_process import pose
from afsd_backend.pysph.interval_acceptance import require_interval_energy
import inlet_acceptance_v2 as inlet
import inlet_source_work as source

DT=3.90625e-6
CLASSIFICATION='HOT_QS_TWO_INVENTORY_BIRTH_DIAGNOSTIC_V1'
INITIAL=ROOT/'outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz'
INITIAL_SHA='31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3'
CONTRACT='docs/TWO_INSERTION_DIAGNOSTIC_CONTRACT_20261005.md'


def require(value,message):
    if not bool(value):raise ValueError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()
def plain(value):
    if isinstance(value,(np.ndarray,np.generic)):return value.tolist()
    if isinstance(value,dict):return {key:plain(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [plain(item) for item in value]
    return value
def write(path,value):
    with Path(path).open('x',encoding='utf-8') as handle:json.dump(plain(value),handle,indent=2,allow_nan=False);handle.write('\n')
def verify(mapping):
    for name,digest in mapping.items():
        path=(ROOT/name).resolve();require(path.is_relative_to(ROOT) and sha(path)==digest,'bound file changed: '+name)
def near(a,b,label,scale=0.):
    a,b=np.asarray(a,float),np.asarray(b,float)
    bound=1024*np.finfo(float).eps*max(float(np.max(np.abs(a),initial=0.)),
        float(np.max(np.abs(b),initial=0.)),scale,np.finfo(float).tiny)
    require(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all() and np.all(np.abs(a-b)<=bound),label)
    return bound


def prepared_from_saved_event(parent,v2,event,p,index):
    # Reconstruct declared birth data only. No geometry/map method is imported or called.
    xy=np.arange(-2,3)*p['dx'];height=(.00325,.00375)[index-1]
    incoming=make_array(np.array([[x,y,height] for x in xy for y in xy]),np.ones(25,dtype=int),p)
    incoming.particle_id[:]=0;incoming.u[:]=0.;incoming.v[:]=0.;incoming.w[:]=-p['feed'];incoming.temperature[:]=723.15
    spec=dict(event_id=event['event_id'],time_s=event['time_s'],source=event['source'],fields=event['field_declarations'],
        heat_capacity=event['heat_capacity'],reference_temperature=event['reference_temperature'])
    prepared,pending,recreated=inlet.propose_insertion(parent,v2,incoming,spec,pose(p,(index-1)*DT))
    require(recreated==event,'reconstructed source event differs from actual event')
    require(event['allocated_particle_ids']==list(range(919+25*(index-1),944+25*(index-1))), 'allocator reused/skipped IDs')
    old=inlet.process_state_for(parent,v2);now=inlet.process_state_for(prepared,pending)
    lookup={int(pid):k for k,pid in enumerate(prepared.particle_id)}
    retained=np.array([lookup[int(pid)] for pid in parent.particle_id])
    require(np.array_equal(old['contact_history'],now['contact_history'][retained]),'birth changed old contact history')
    require(np.array_equal(old['free'],now['free'][retained]),'unexpected old control change in short two-birth scope')
    born=np.isin(prepared.particle_id,event['allocated_particle_ids'])
    require(not np.any(now['contact_history'][born]) and not np.any(prepared.S.reshape(-1,3,3)[born]),'nonzero virgin history')
    require(pending['source_events'][:-1]==v2['source_events'] and
        pending['accepted_bindings']==v2['accepted_bindings'] and pending['source_work_events']==v2['source_work_events'],
        'birth reset historical source or mechanical ledger')
    return prepared,pending,now,retained,born


def fields_and_graph(stepdir,prepared,actual,state,aligned,p):
    f=dict(np.load(stepdir/'actual_fields.npz'));g=dict(np.load(stepdir/'visible_graph.npz'))
    require(np.array_equal(f['particle_id'],prepared.particle_id) and np.array_equal(g['particle_id'],prepared.particle_id),
        'saved field/graph row IDs differ from prepared PA')
    require(np.array_equal(actual.particle_id,prepared.particle_id),'raw checkpoint row order differs')
    require(np.array_equal(f['free'],aligned['free']),'actual free axes differ from V2 pending control')
    expected=prepared.extract_particles(np.arange(len(prepared.x),dtype=np.uint32))
    expected.x[:],expected.y[:],expected.z[:]=f['position'].T
    expected.u[:],expected.v[:],expected.w[:]=(f['displacement']/DT).T
    expected.rho[:]=f['density'];expected.S[:]=f['stress'].ravel();expected.ep[:]=f['ep']
    expected.temperature[:]=f['temperature'];expected.strain_rate[:]=f['rate'];expected.ep_rate[:]=(f['ep']-prepared.ep)/DT
    expected.D[:]=(.5*(f['gradient']+f['gradient'].swapaxes(1,2))).ravel()
    expected.plastic_work[:]+=f['volume']*f['plastic_work_density']
    expected.return_dissipation[:]+=f['volume']*f['return_dissipation_density']
    require(set(expected.properties)==set(actual.properties),'checkpoint property schema changed')
    for key,prop in expected.properties.items():
        require(np.array_equal(prop.get_npy_array(),actual.properties[key].get_npy_array()),'actual checkpoint field differs: '+key)
    require(np.array_equal(positions(prepared)+f['displacement'],f['position']),'saved displacement does not reach actual endpoint')
    require(np.array_equal(state['contact_history'],f['contact_history']) and np.array_equal(state['free'],aligned['free']),
        'checkpoint physical history/free not actual field values')
    for key,field in (('position0',positions(prepared)),('position1',f['position']),('volume',f['volume'])):
        require(np.array_equal(g[key],field),'graph field differs: '+key)
    kept=g['kept_raw_edge_indices'];blocked=g['blocked_raw_edge_indices']
    require(np.array_equal(np.sort(np.r_[kept,blocked]),np.arange(len(g['raw_i']))),'raw edge partition incomplete')
    for key in ('i','j','g'):require(np.array_equal(g[key],g['raw_'+key][kept]),'not exact common filtered graph: '+key)
    require(np.all(g['i']<g['j']) and np.all(g['i']>=0) and np.all(g['j']<len(prepared.x)), 'noncanonical graph indices')
    edges=np.sort(np.column_stack((prepared.particle_id[g['i']],prepared.particle_id[g['j']])).astype('<i8'),axis=1)
    edges=edges[np.lexsort((edges[:,1],edges[:,0]))]
    fingerprint=hashlib.sha256(edges.tobytes()).hexdigest()
    require(fingerprint==read(stepdir/'visible_graph.json')['graph_fingerprint'],'persistent-ID graph fingerprint differs')
    require(f['pair_conductance'].shape==(len(g['i']),) and np.all(f['pair_conductance']>=0),'invalid saved heat graph')
    return f,g,fingerprint


def audit(base,output):
    started=time.perf_counter();directory=base/'diagnostic';manifest=read(directory/'manifest.json');sup=read(directory/'supervision.json')
    require(manifest['status']=='PASS_TWO_INSERTION_DIAGNOSTIC_NOT_JOINT_OR_PROCESS_ACCEPTANCE' and
        manifest['diagnostic_steps']==2 and manifest['diagnostic_time_s']==2*DT,'complete two-step diagnostic absent')
    require(sup['exit_code']==0 and sup['before_after_hashes_match'] and sup['manifest_sha256']==sha(directory/'manifest.json') and
        sup['log_sha256']==sha(directory/'run.log'),'external execution receipt mismatch')
    bindings=dict(manifest['source_sha256']);bindings.update(manifest['artifacts_sha256']);bindings.update(manifest['reference_sha256'])
    bindings.update({rel(path):sha(path) for path in (directory/'manifest.json',directory/'supervision.json',directory/'run.log',Path(__file__))})
    verify(bindings);require(sha(INITIAL)==INITIAL_SHA,'independent origin differs')
    saved=load(str(INITIAL));parent=saved['arrays']['solid'];origin=saved['solver_data'];p=origin['parameters']
    registry=inlet.inventory.initialize_registry(parent,CONTRACT,acknowledgement=inlet.inventory.ACKNOWLEDGEMENT)
    plane=dict(frame='WORLD_Z',plane_z_m=.0025,unit='m',source=CONTRACT)
    v2=inlet.initialize_acceptance_state(parent,registry,origin,plane,CONTRACT)
    interval=copy.deepcopy(origin['interval_energy']);previous_sha=INITIAL_SHA;previous_contact=None;source_reports=[];rows=[]
    initial_mass=math.fsum(parent.m);initial_T=parent.temperature.copy();initial_position={int(pid):point.copy() for pid,point in zip(parent.particle_id,positions(parent))}
    frames=[dict(time_s=0.,particle_id=parent.particle_id.copy(),position=positions(parent).copy(),temperature=parent.temperature.copy(),
        stress=parent.S.reshape(-1,3,3).copy(),ep=parent.ep.copy())]
    for index in (1,2):
        stepdir=directory/('step_'+str(index).zfill(3));event=read(stepdir/'source_event.json');receipt=read(stepdir/'mechanical_receipt.json')
        checkpoint=stepdir/'DIAGNOSTIC_NOT_RESTARTABLE.npz';actual_saved=load(str(checkpoint));actual=actual_saved['arrays']['solid'];state=actual_saved['solver_data']
        require(state['classification']==CLASSIFICATION and state['model_scope']==inlet.DIAGNOSTIC_SCOPE and
            state['accepted'] is False and state['restart_allowed'] is False and state['process_acceptance'] is False,
            'diagnostic became a restartable process')
        require(state['parent_checkpoint_sha256']==previous_sha and state['time_s']==index*DT and state['step']==index,
            'checkpoint parent/time chain differs')
        require(state['source_sha256']==manifest['source_sha256'] and state['parameters']==p and
            state['parameter_provenance']==origin['parameter_provenance'],'checkpoint source/parameters differ')
        require('spatial_model_transition' not in state and state['bore_observation']['contract']==inlet.OBSERVER_CONTRACT,
            'published fake transition or legacy observer')
        prepared,pending,aligned,retained,born=prepared_from_saved_event(parent,v2,event,p,index)
        f,g,fingerprint=fields_and_graph(stepdir,prepared,actual,state,aligned,p)
        journal=read(stepdir/'FULL_MAP_DIAGNOSTIC_NOT_ACCEPTED.json');contact={key:np.asarray(value) if isinstance(value,list) else value
            for key,value in journal['contact'].items()}
        energy_file=read(stepdir/'actual_energy.json');energy=energy_file['mechanical'];thermal=energy_file['thermal']
        require(journal['energy']==energy and journal['thermal']==thermal['ledger'],'saved complete-map ledgers differ')
        require(np.array_equal(contact['history'],f['contact_history']),'saved actual contact history differs')
        require(not np.any(contact['normal_potential0_J'][born]) and not np.any(contact['tangential_potential0_J'][born]),
            'nonzero newborn contact storage')
        contact0=np.column_stack((contact['normal_potential0_J'].sum(axis=1),contact['tangential_potential0_J'].sum(axis=1)))
        if previous_contact is None:
            # Initial per-particle spring storage is read from the hash-bound map, not a re-evaluated contact law.
            parent_contact=contact0[retained]
            require(not np.any(parent_contact[:,1]),'initial tangential spring storage is nonzero')
        else:
            parent_contact=np.array([previous_contact[int(pid)] for pid in parent.particle_id])
            near(contact0[retained],parent_contact,'contact storage does not telescope across real restart')
        result=dict(displacement=f['displacement'],position=f['position'],force=f['force'],density=f['density'],stress=f['stress'],
            temperature=f['temperature'],volume=f['volume'],dt_s=DT,
            material=(f['stress'],f['ep'],f['plastic_work_density'],f['return_dissipation_density'],f['rate']),
            contact=contact,thermal=thermal)
        old_source=read(stepdir/'source_work.json')
        calculated=source.evaluate_source_work(parent,prepared,actual,event,parameters=p,free=aligned['free'],result=result,
            mechanical_energy=energy,parent_contact_storage_J=parent_contact,external_port=old_source['external_port'])
        require(calculated==old_source,'saved-field source-work recomputation differs')
        work=source.source_work_record(calculated,stepdir/'source_work.json')
        require(work==read(stepdir/'source_work_v2.json') and receipt['source_work_sha256']==inlet._sha(work),'V2 source work unbound')
        interval=require_interval_energy(interval,energy)
        require(plain(state['interval_energy'])==interval,'whole mechanical energy prefix reset/changed')
        binding=inlet.mechanical_binding(actual,pending,f['contact_history'],aligned['free'],index*DT)
        require(receipt['binding']==binding and receipt['scope']==inlet.DIAGNOSTIC_SCOPE and
            receipt['accepted'] is False and receipt['restart_allowed'] is False and receipt['joint_inlet_benchmark']=='NOT_RUN',
            'mechanical receipt class or actual field binding differs')
        verify(receipt['audit_artifacts_sha256'])
        for key in inlet.GATES:
            item=receipt['gates'][key];path=ROOT/item['artifact'];gate=read(path)
            require(item['status']=='PASS' and sha(path)==item['evidence_sha256'] and gate['status']=='PASS' and
                gate['binding']==binding and gate['source_sha256']==manifest['source_sha256'] and
                gate['source_work_report_sha256']==sha(stepdir/'source_work.json'),'unbound actual gate '+key)
        def verifier(given):
            require(given==read(stepdir/'mechanical_receipt.json'),'receipt changed during reseal')
            verify(given['audit_artifacts_sha256']);return True
        resealed=inlet.seal_mechanical_acceptance(actual,pending,f['contact_history'],aligned['free'],index*DT,
            receipt,work,verifier,lambda t:pose(p,t))
        restored=inlet.restore_acceptance_state(json.dumps(state['inlet_acceptance']),actual)
        require(resealed==restored and state['bore_observation']==restored['bore_observation'],'actual V2 history/observer reseal differs')
        require(restored['source_events'][:-1]==v2['source_events'] and restored['accepted_bindings'][:-1]==v2['accepted_bindings'] and
            restored['source_work_events'][:-1]==v2['source_work_events'],'mechanical seal discarded prior history')
        actual_force=float(np.max(np.abs(f['force'][aligned['free']])))
        require(actual_force<=journal['metrics']['force_limit_N'],'saved actual free force fails original limit')
        dT=actual.temperature-prepared.temperature
        measured=prepared.m*p['cp']*dT;particle_thermal=thermal['particle_energy_J'];heat=np.asarray(particle_thermal['total_heat_input_J'])
        require(np.max(np.abs(measured-heat))<=1e-12+1e-10*np.abs(heat).sum(),'saved per-particle heat balance failed')
        rank=read(stepdir/'reference_rank.json');eig=np.load(stepdir/'reference_rank_arrays.npz')['stiffness_eigenvalues_N_m']
        require(rank['status']=='PASS' and eig[0]>1e-12*eig[-1]>0 and rank['graph_fingerprint']==fingerprint,'actual graph rank receipt differs')
        require(read(stepdir/'continuous_shell.json')['status']=='CERTIFIED_NO_HIDDEN_SHELL_CONTACT' and
            read(stepdir/'visible_path.json')['status']=='CERTIFIED_VISIBLE_SUPPORT_PATH','continuous geometry certificate absent')
        require(restored['bore_observation']['birth_events'][-1]['crossing_added_by_birth']==[],'birth counted as crossing')
        reload=read(stepdir/'reload_and_sort.json');require(reload['checkpoint_sha256']==sha(checkpoint) and
            reload['particle_ids_before']==actual.particle_id.tolist() and reload['particle_ids_after']==actual.particle_id[::-1].tolist(),
            'actual reload permutation mismatch')
        reordered=actual.extract_particles(np.arange(len(actual.x)-1,-1,-1,dtype=np.uint32))
        require(inlet.restore_acceptance_state(json.dumps(restored),reordered)==restored and
            inlet.inventory._array_digest(reordered)==reload['actual_pa_sha256'],'sorted checkpoint loses fields/registry')
        transaction=read(stepdir/'transaction.json');require(transaction['checkpoint_sha256']==sha(checkpoint) and
            transaction['parent_checkpoint_sha256']==previous_sha and transaction['graph_fingerprint']==fingerprint,'transaction provenance differs')
        source_reports.append(calculated)
        storage1=np.column_stack((contact['normal_potential1_J'].sum(axis=1),contact['tangential_potential1_J'].sum(axis=1)))
        previous_contact={int(pid):value for pid,value in zip(actual.particle_id,storage1)}
        frames.append(dict(time_s=index*DT,particle_id=actual.particle_id.copy(),position=positions(actual).copy(),
            temperature=actual.temperature.copy(),stress=actual.S.reshape(-1,3,3).copy(),ep=actual.ep.copy()))
        rows.append(dict(step=index,particle_count=len(actual.x),birth_ids=event['allocated_particle_ids'],
            maximum_free_force_N=actual_force,maximum_temperature_change_K=float(np.max(np.abs(dT))),
            temperature_min_K=float(actual.temperature.min()),temperature_max_K=float(actual.temperature.max()),
            stress_frobenius_max_Pa=float(np.linalg.norm(f['stress'],axis=(1,2)).max()),ep_max=float(actual.ep.max()),
            plastic_work_J=energy['plastic_work_J'],raw_mechanical_J=energy['global_raw_residual_J'],
            raw_open_energy_J=calculated['open_energy']['open_raw_residual_J'],
            reaction_work_J=calculated['computed_terms']['mechanical_reaction_work']['value'],
            activated_ids=restored['control_events'][-1]['activated_ids'],released_ids=restored['control_events'][-1]['released_ids'],
            checkpoint_sha256=sha(checkpoint),source_increment=calculated['source_increment']))
        parent,v2,previous_sha=reordered,restored,sha(checkpoint)
    source_mass=math.fsum(row['source_increment']['mass_kg'] for row in rows)
    source_heat=math.fsum(row['source_increment']['sensible_heat_J'] for row in rows)
    near(math.fsum(parent.m)-initial_mass,source_mass,'mass does not telescope',initial_mass)
    first,last=source_reports[0]['open_energy'],source_reports[-1]['open_energy']
    storage_change=math.fsum(last[key+'_storage_after_J']-first[key+'_storage_before_J'] for key in ('material','sensible','contact'))
    near(storage_change,math.fsum(row['open_energy']['total_storage_increment_J'] for row in source_reports),
        'total material/thermal/contact storage does not telescope',abs(first['sensible_storage_before_J'])+abs(last['sensible_storage_after_J']))
    raw_open=math.fsum(row['open_energy']['open_raw_residual_J'] for row in source_reports)
    raw_abs=math.fsum(abs(row['open_energy']['open_raw_residual_J']) for row in source_reports)
    manifest_sum=manifest['cumulative_source_ledger']
    near(raw_open,manifest_sum['sum']['open_raw_residual_J'],'manifest source raw differs')
    near(source_heat,manifest_sum['source_sensible_heat_J'],'manifest source heat differs')
    require(v2['protocol_registry']['issued_particle_count']==50 and len(v2['source_events'])==2 and
        v2['mechanical_step_count']==2 and v2['scope_lock']==inlet.DIAGNOSTIC_SCOPE,'source count/history scope differs')
    verify(bindings)
    arrays={}
    for index,frame in enumerate(frames):
        for name,value in frame.items():arrays['frame_'+str(index)+'_'+name]=value
    # Ragged real frames: unborn particles are absent, not zero-filled or interpolated.
    np.savez_compressed(output/'actual_ragged_frames.npz',**arrays)
    report=dict(status='PASS_SAVED_TWO_INSERTION_DIAGNOSTIC_READBACK_NOT_JOINT_OR_PROCESS_ACCEPTANCE',
        source_and_input_sha256=bindings,accepted=False,restart_allowed=False,accepted_time_advanced_s=0.,
        physical_evaluations=0,nnps_evaluations=0,force_evaluations=0,thermal_evaluations=0,contact_evaluations=0,
        diagnostic_time_s=2*DT,diagnostic_steps=2,initial_particle_count=918,final_particle_count=len(parent.x),
        initial_mass_kg=initial_mass,final_mass_kg=math.fsum(parent.m),added_mass_kg=source_mass,source_heat_J=source_heat,
        total_storage_increment_J=storage_change,source_registry=v2['protocol_registry']['cumulative'],
        mechanical_interval_energy=interval,raw_open_energy_J=raw_open,absolute_raw_open_sum_J=raw_abs,
        final_bore_observation=v2['bore_observation'],rows=rows,control_release='NOT_EXERCISED',
        simulated_plastic_response='NOT_EXERCISED' if all(row['ep_max']==0. for row in rows) else 'OBSERVED',
        thermal_change='NOT_EXERCISED' if all(row['maximum_temperature_change_K']==0. for row in rows) else 'OBSERVED',
        worker_wall_seconds=manifest['wall_seconds'],external_wall_seconds=sup['wall_seconds'],readback_wall_seconds=time.perf_counter()-started,
        initial_normal_contact_storage_basis='HASH_BOUND_SAVED_FIRST_MAP; CONTACT_LAW_NOT_REEVALUATED',
        limitations=['Two coupled elastic startup steps with explicit hot inventory births; not a continuous mass-flow inlet.',
            'Preloaded-long-bar comparison and time refinement were not run.',
            'Saved gate artifacts are rehashed; continuous geometry certificates and rank eigensolve are not rerun.',
            'Source-work algebra, PA update, thermodynamic stocks, per-particle heat, receipt seal, IDs and complete history are independently read back.'],
        ragged_frames_sha256=sha(output/'actual_ragged_frames.npz'))
    write(output/'summary.json',report);return report


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    base=(ROOT/args.input).resolve() if not Path(args.input).is_absolute() else Path(args.input).resolve()
    output=(ROOT/args.output).resolve() if not Path(args.output).is_absolute() else Path(args.output).resolve()
    require(base.is_relative_to(ROOT/'outputs') and output.is_relative_to(ROOT/'outputs') and not output.exists(),'new project output directory required')
    output.mkdir(parents=True)
    try:
        result=audit(base,output);print(json.dumps({key:result[key] for key in
            ('status','diagnostic_steps','diagnostic_time_s','added_mass_kg','source_heat_J','raw_open_energy_J','readback_wall_seconds')}))
    except Exception as error:
        write(output/'FAILURE.json',dict(status='FAIL_SAVED_ONLY_READBACK',error=repr(error),traceback=traceback.format_exc(),physical_evaluations=0))
        raise


if __name__=='__main__':main()

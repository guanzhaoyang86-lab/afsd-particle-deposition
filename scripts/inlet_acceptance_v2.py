"""COMPONENT_ONLY acceptance transactions; no inlet physics is implemented.

See docs/INLET_ACCEPTANCE_V2_COMPONENT_CONTRACT_20261004.md. The external
receipt verifier is a trust boundary, not a solver or a self-issued PASS.
"""
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import inlet_protocol as inventory
from afsd_backend.pysph import interval_acceptance as crossing

CLASSIFICATION='INLET_ACCEPTANCE_V2_COMPONENT_ONLY'
OBSERVER_CONTRACT='INITIAL_OR_BIRTH_ABOVE_SHOULDER_THROUGH_OPEN_SQUARE_V2'
FIXTURE_SCOPE='MANUFACTURED_COMPONENT_FIXTURE_NOT_AFSD'
DIAGNOSTIC_SCOPE='DIAGNOSTIC_COUPLED_INLET_BENCHMARK_NOT_PROCESS_ACCEPTANCE'
PROCESS_SCOPE='QUALIFIED_COUPLED_PROCESS_STEP'
WORK_CONVENTION='PRESSURE_FLOW_SEPARATE_FROM_ACTUATOR_REACTION_NO_DUPLICATION'
GATES=('force','heat','energy','geometry_overlap','NNPS','visible_support','gradient_rank','contact')
WORK_TERMS=('pressure_flow_work','elastic_convected_energy','mechanical_reaction_work')
IMMUTABLE=('m','x0','y0','z0','material_id','component_id','tracer_type','wholebar_tag')


def _json(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
def _sha(value):return hashlib.sha256(_json(value).encode()).hexdigest()
def _positions(pa):return np.column_stack((pa.x,pa.y,pa.z))
def _order(pa):return np.argsort(pa.particle_id)


def _seal(value):
    result=copy.deepcopy(value);result.pop('state_sha256',None)
    result['state_sha256']=_sha(result)
    return result


def _hex(value):
    return isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value)


def _identity(pa):
    order=_order(pa);keys=list(IMMUTABLE)
    if 'as_rs_tag' in pa.properties:keys.append('as_rs_tag')
    return {key:np.asarray(getattr(pa,key))[order].tolist() for key in keys}


def _schema(pa):
    return {key:[prop.get_npy_array().dtype.str,int(pa.stride.get(key,1))] for key,prop in sorted(pa.properties.items())}


def _history(pa,contact_history,free):
    n=pa.get_number_of_particles();h=np.asarray(contact_history);f=np.asarray(free)
    if h.shape!=(n,8,3) or not np.isfinite(h).all() or f.shape!=(n,3) or f.dtype.kind!='b':
        raise ValueError('finite N by 8 by 3 contact history and boolean N by 3 free required')
    order=_order(pa)
    return dict(contact_history=h[order].tolist(),free=f[order].tolist())


def _legacy_observer(observer):
    result=copy.deepcopy(observer)
    for key in ('initial_ids','initial_eligible_ids','birth_events'):result.pop(key,None)
    result['contract']=crossing.CROSSING_CONTRACT
    return result


def _validate_observer(observer):
    if observer.get('contract')!=OBSERVER_CONTRACT:raise ValueError('V2 birth observer contract required')
    crossing.validate_bore_crossing_state(_legacy_observer(observer))
    initial=observer['initial_ids'];born=[];eligible=set(observer['initial_eligible_ids'])
    if initial!=sorted(set(initial)) or not eligible.issubset(initial):raise ValueError('invalid initial observer identities')
    for event in observer['birth_events']:
        ids=event['particle_ids'];qualified=event['eligible_ids'];positions=np.asarray(event['positions_m'],float)
        if (ids!=sorted(set(ids)) or not set(qualified).issubset(ids) or positions.shape!=(len(ids),3)
                or not np.isfinite(positions).all() or event['crossing_added_by_birth']!=[]
                or not observer['initial_time_s']<=event['time_s']<=observer['last_time_s']):
            raise ValueError('invalid birth observer event')
        expected=[pid for pid,z in zip(ids,positions[:,2]) if z>event['shoulder_z_m']]
        if qualified!=expected:raise ValueError('birth above-shoulder qualification mismatch')
        born.extend(ids);eligible.update(qualified)
    if (len(born)!=len(set(born)) or set(initial)&set(born)
            or sorted(initial+born)!=observer['particle_ids'] or sorted(eligible)!=observer['eligible_ids']):
        raise ValueError('observer birth inventory or eligibility mismatch')


def _validate_state(state,pa):
    if not isinstance(state,dict) or _seal(state)!=state:raise ValueError('acceptance state checksum mismatch')
    if (state.get('version')!=2 or state.get('classification')!=CLASSIFICATION
            or state.get('integration_status')!='COMPONENT_ONLY' or state.get('phase') not in ('SEALED','PENDING')):
        raise ValueError('unsupported inlet acceptance state')
    inventory.restore_registry(_json(state['protocol_registry']),pa)
    ids=np.asarray(pa.particle_id)[_order(pa)].astype(int).tolist()
    if ids!=state['particle_ids'] or _identity(pa)!=state['immutable_inventory'] or _schema(pa)!=state['property_schema']:
        raise ValueError('particle identity, mass, or passive history changed')
    n=len(ids)
    h=np.asarray(state['contact_history']);f=np.asarray(state['free']);base=np.asarray(state['uncontrolled_free'])
    if h.shape!=(n,8,3) or not np.isfinite(h).all() or f.shape!=(n,3) or f.dtype.kind!='b' or base.shape!=(n,3) or base.dtype.kind!='b':
        raise ValueError('invalid sealed contact/free layout')
    _validate_observer(state['bore_observation'])
    observer=state['bore_observation']
    if (observer['particle_ids']!=ids or observer['masses_kg']!=state['immutable_inventory']['m']
            or observer['component_ids']!=state['immutable_inventory']['component_id']
            or observer['last_time_s']!=state['time_s']):raise ValueError('observer and accepted inventory/time differ')
    events=state['source_events'];registry=state['protocol_registry']
    if [row['event_id'] for row in events]!=registry['event_ids']:raise ValueError('source event history mismatch')
    cumulative=dict(mass_kg=0.,momentum_kg_m_s=[0.,0.,0.],sensible_heat_J=0.,kinetic_energy_J=0.,volume_m3=0.)
    born=[]
    for row in events:
        born.extend(row['allocated_particle_ids'])
        for key,value in row['source_increment'].items():
            cumulative[key]=(np.asarray(cumulative[key])+np.asarray(value)).tolist()
    if cumulative!=registry['cumulative'] or len(born)!=registry['issued_particle_count']:
        raise ValueError('cumulative inlet source ledger mismatch')
    if born!=list(range(registry['initial_max_id']+1,registry['high_watermark']+1)):
        raise ValueError('allocated persistent ID sequence mismatch')
    if [b['event_id'] for b in observer['birth_events']]!=registry['event_ids']:
        raise ValueError('birth and source event IDs differ')
    if [pid for b in observer['birth_events'] for pid in b['particle_ids']]!=born:
        raise ValueError('birth and source allocated IDs differ')
    bindings=state['accepted_bindings'];works=state['source_work_events']
    if (len(bindings)!=state['mechanical_step_count'] or len(works)!=len(bindings)
            or len({row['step_id'] for row in bindings})!=len(bindings)
            or [(row['step_id'],row['receipt_sha256']) for row in bindings]!=[(row['step_id'],row['receipt_sha256']) for row in works]):
        raise ValueError('mechanical receipt/source-work sequence mismatch')
    if ((not bindings and state['scope_lock'] is not None) or
            (bindings and (state['scope_lock'] not in (FIXTURE_SCOPE,DIAGNOSTIC_SCOPE,PROCESS_SCOPE)
                or any(row['scope']!=state['scope_lock'] for row in bindings)))):
        raise ValueError('mechanical receipt scope history mismatch')
    if state['process_acceptance_by_this_module'] is not False:
        raise ValueError('component module cannot self-certify process acceptance')
    actual_active=sorted(pa.particle_id[np.asarray(pa.actuator).astype(bool)].astype(int).tolist())
    if actual_active!=state['control_active_ids']:raise ValueError('sealed control identities differ from PA')
    if state['phase']=='PENDING' and state['pending']['prepared_pa_sha256']!=inventory._array_digest(pa):
        raise ValueError('pending prepared PA mismatch')
    return state


def restore_acceptance_state(serialized_json,pa):
    """Restore PA and V2 state together; never infer an allocator from max(ID)."""
    return copy.deepcopy(_validate_state(json.loads(serialized_json),pa))


def initialize_acceptance_state(pa,registry,process_state,control_plane,provenance):
    old=inventory.restore_registry(_json(registry),pa)
    if old['transaction_count']!=0:raise ValueError('initialization cannot discard an existing inlet source history')
    plane=copy.deepcopy(control_plane)
    if (set(plane)!={'frame','plane_z_m','unit','source'} or plane['frame']!='WORLD_Z' or plane['unit']!='m'
            or not isinstance(plane['source'],str) or not plane['source'].strip()
            or not math.isfinite(float(plane['plane_z_m'])) or not isinstance(provenance,str) or not provenance.strip()):
        raise ValueError('declared world-z control plane and provenance required')
    hist=_history(pa,process_state['contact_history'],process_state['free']);order=_order(pa)
    active=(np.asarray(pa.component_id)==1)&(np.asarray(pa.z)>=float(plane['plane_z_m']))
    if not np.array_equal(np.asarray(pa.actuator),active.astype(int)):
        raise ValueError('initial actuator identity differs from explicit control plane')
    free=np.asarray(hist['free'],bool);active=active[order]
    if np.any(free[active,2]) or np.any(~free[active,:2]):raise ValueError('initial inlet control must constrain only z')
    base=free.copy();base[active,2]=True
    observation=copy.deepcopy(process_state['bore_observation'])
    crossing.validate_bore_crossing_state(observation)
    observation.update(contract=OBSERVER_CONTRACT,initial_ids=observation['particle_ids'].copy(),
        initial_eligible_ids=observation['eligible_ids'].copy(),birth_events=[])
    state=dict(version=2,classification=CLASSIFICATION,integration_status='COMPONENT_ONLY',phase='SEALED',
        provenance=provenance,time_s=float(process_state['time_s']),mechanical_step_count=0,protocol_registry=old,
        particle_ids=np.asarray(pa.particle_id)[order].astype(int).tolist(),immutable_inventory=_identity(pa),
        property_schema=_schema(pa),control_active_ids=np.asarray(pa.particle_id)[order][active].astype(int).tolist(),
        **hist,uncontrolled_free=base.tolist(),control_plane=plane,bore_observation=observation,
        source_events=[],control_events=[],accepted_bindings=[],source_work_events=[],pending=None,
        scope_lock=None,process_acceptance_by_this_module=False)
    state=_seal(state);_validate_state(state,pa)
    return state


def process_state_for(pa,state):
    """Recover ID-aligned history for a possibly sorted PA without mutating it."""
    _validate_state(state,pa)
    index={pid:k for k,pid in enumerate(state['particle_ids'])}
    rows=np.array([index[int(pid)] for pid in pa.particle_id])
    return dict(time_s=state['time_s'],contact_history=np.asarray(state['contact_history'])[rows].copy(),
        free=np.asarray(state['free'],bool)[rows].copy(),bore_observation=copy.deepcopy(state['bore_observation']))


def _prepare_control(candidate,state,parent_sha):
    result=copy.deepcopy(state);ids=np.asarray(candidate.particle_id);order=_order(candidate)
    previous=set(result['control_active_ids'])
    active=(candidate.component_id==1)&(candidate.z>=float(result['control_plane']['plane_z_m']))
    selected=set(ids[active].astype(int).tolist());candidate.actuator[:]=active.astype(int)
    base=np.asarray(result['uncontrolled_free'],bool);free=base.copy();free[active[order],2]=False
    result['free']=free.tolist()
    result['control_active_ids']=sorted(selected)
    control=dict(time_s=result['time_s'],active_ids=sorted(selected),activated_ids=sorted(selected-previous),
        released_ids=sorted(previous-selected),selection='COMPONENT_1_WORLD_Z_GREATER_OR_EQUAL_DECLARED_PLANE',
        mechanical_reaction_work='NOT_COMPUTED',status='CANDIDATE_CONTROL_TRANSITION_NOT_ACCEPTED')
    result['control_events'].append(control)
    reg=copy.deepcopy(result['protocol_registry']);reg['particle_state_sha256']=inventory._array_digest(candidate)
    result['protocol_registry']=inventory._sealed(reg)
    result.update(phase='PENDING',pending=dict(parent_state_sha256=parent_sha,
        prepared_pa_sha256=reg['particle_state_sha256'],prepared_positions_m=_positions(candidate)[order].tolist(),
        control_event_index=len(result['control_events'])-1))
    result=_seal(result);_validate_state(result,candidate)
    return result


def propose_mechanical_step(pa,state):
    _validate_state(state,pa)
    if state['phase']!='SEALED':raise ValueError('cannot prepare from an unaccepted pending transaction')
    candidate=pa.extract_particles(np.arange(pa.get_number_of_particles(),dtype=np.uint32))
    return candidate,_prepare_control(candidate,state,state['state_sha256'])


def propose_insertion(pa,state,incoming,specification,tool_pose):
    _validate_state(state,pa)
    if state['phase']!='SEALED':raise ValueError('previous pending transaction is not mechanically accepted')
    if float(specification['time_s'])!=state['time_s']:raise ValueError('birth time must equal accepted boundary time')
    candidate,registry,event=inventory.propose_insertion(pa,state['protocol_registry'],incoming,specification)
    result=copy.deepcopy(state);result['protocol_registry']=registry;result['source_events'].append(event)
    old_ids=state['particle_ids'];new_ids=event['allocated_particle_ids'];ids=np.asarray(candidate.particle_id)
    mapping={pid:k for k,pid in enumerate(old_ids)};rows=[];base=[]
    for pid in sorted(ids.astype(int).tolist()):
        rows.append(state['contact_history'][mapping[pid]] if pid in mapping else np.zeros((8,3)).tolist())
        base.append(state['uncontrolled_free'][mapping[pid]] if pid in mapping else [True,True,True])
    result.update(particle_ids=sorted(ids.astype(int).tolist()),immutable_inventory=_identity(candidate),
        contact_history=rows,uncontrolled_free=base,free=copy.deepcopy(base))
    obs=copy.deepcopy(state['bore_observation']);origin,_=crossing._pose(tool_pose)
    by_id={int(pid):k for k,pid in enumerate(ids)};born_positions=[_positions(candidate)[by_id[pid]].tolist() for pid in new_ids]
    eligible=[pid for pid,pos in zip(new_ids,born_positions) if pos[2]>origin[2]]
    birth=dict(event_id=event['event_id'],time_s=state['time_s'],particle_ids=new_ids.copy(),positions_m=born_positions,
        shoulder_z_m=float(origin[2]),eligible_ids=eligible,crossing_added_by_birth=[])
    obs['birth_events'].append(birth);obs['eligible_ids']=sorted(set(obs['eligible_ids'])|set(eligible))
    order=_order(candidate);obs.update(particle_ids=ids[order].astype(int).tolist(),
        component_ids=candidate.component_id[order].astype(int).tolist(),masses_kg=candidate.m[order].tolist())
    result['bore_observation']=obs
    result=_prepare_control(candidate,result,state['state_sha256'])
    return candidate,result,copy.deepcopy(event)


def mechanical_binding(actual_pa,pending,contact_history,free,time_s):
    """Pure digest template for an external solver's receipt, never a PASS."""
    return dict(pending_state_sha256=pending['state_sha256'],prepared_pa_sha256=pending['pending']['prepared_pa_sha256'],
        actual_pa_sha256=inventory._array_digest(actual_pa),actual_history_sha256=_sha(_history(actual_pa,contact_history,free)),
        start_time_s=pending['time_s'],end_time_s=float(time_s))


def seal_mechanical_acceptance(actual_pa,pending,contact_history,free,time_s,receipt,source_work,
        verify_receipt,pose_at_time):
    """Seal only after externally verified coupled receipts; no default zeros.

    The actual PA may have changed coordinates/stress/temperature/history and
    ordering. Existing masses, identity, passive labels, and intended free axes
    must be unchanged. The caller owns atomic PA/state checkpoint publication.
    """
    if _seal(pending)!=pending or pending.get('phase')!='PENDING':raise ValueError('sealed pending transaction required')
    ids=np.asarray(actual_pa.particle_id)[_order(actual_pa)].astype(int).tolist()
    if ids!=pending['particle_ids'] or _identity(actual_pa)!=pending['immutable_inventory'] or _schema(actual_pa)!=pending['property_schema']:
        raise ValueError('actual accepted PA lost or changed persistent inventory')
    hist=_history(actual_pa,contact_history,free)
    if hist['free']!=pending['free']:raise ValueError('actual free axes differ from prepared control')
    expected_active=pending['control_events'][-1]['active_ids']
    if sorted(actual_pa.particle_id[np.asarray(actual_pa.actuator).astype(bool)].astype(int).tolist())!=expected_active:
        raise ValueError('actual actuator identities differ from prepared control')
    binding=mechanical_binding(actual_pa,pending,contact_history,free,time_s)
    if not math.isfinite(binding['end_time_s']) or binding['end_time_s']<=binding['start_time_s']:
        raise ValueError('mechanical acceptance must advance finite time')
    if receipt.get('binding')!=binding:raise ValueError('qualified receipt does not bind actual PA/history/time')
    scope=receipt.get('scope')
    if scope not in (FIXTURE_SCOPE,DIAGNOSTIC_SCOPE,PROCESS_SCOPE):raise ValueError('explicit qualified receipt scope required')
    if pending['scope_lock'] is not None and scope!=pending['scope_lock']:
        raise ValueError('receipt scope history is locked; no fixture/diagnostic/process conversion')
    names=GATES+(('joint_inlet_benchmark',) if scope==PROCESS_SCOPE else ())
    for key in names:
        row=receipt.get('gates',{}).get(key,{})
        if row.get('status')!='PASS' or not _hex(row.get('evidence_sha256')):
            raise ValueError('missing qualified receipt gate: '+key)
    step_id=receipt.get('step_id')
    if not isinstance(step_id,str) or not step_id.strip() or step_id in [r['step_id'] for r in pending['accepted_bindings']]:
        raise ValueError('unique mechanical step receipt ID required')
    if not callable(verify_receipt) or verify_receipt(copy.deepcopy(receipt)) is not True:
        raise ValueError('external qualified receipt verification required')
    if source_work.get('convention')!=WORK_CONVENTION:raise ValueError('explicit source-work nonduplication convention required')
    for key in WORK_TERMS:
        row=source_work.get(key,{})
        if (row.get('status')!='COMPUTED' or row.get('unit')!='J' or not _hex(row.get('evidence_sha256'))
                or not isinstance(row.get('value'),(float,int)) or isinstance(row.get('value'),bool)
                or not math.isfinite(float(row['value']))):raise ValueError('source work NOT_COMPUTED or invalid: '+key)
    if receipt.get('source_work_sha256')!=_sha(source_work):raise ValueError('source work is not bound by qualified receipt')
    result=copy.deepcopy(pending);order=_order(actual_pa)
    x0=np.asarray(pending['pending']['prepared_positions_m']);x1=_positions(actual_pa)[order]
    observation=crossing.observe_square_bore_crossings(_legacy_observer(pending['bore_observation']),x0,x1,
        np.asarray(ids),actual_pa.component_id[order],actual_pa.m[order],pending['time_s'],float(time_s),pose_at_time)['next_state']
    observation.update(contract=OBSERVER_CONTRACT,initial_ids=pending['bore_observation']['initial_ids'],
        initial_eligible_ids=pending['bore_observation']['initial_eligible_ids'],birth_events=pending['bore_observation']['birth_events'])
    registry=copy.deepcopy(pending['protocol_registry']);registry['particle_state_sha256']=inventory._array_digest(actual_pa)
    result.update(protocol_registry=inventory._sealed(registry),**hist,bore_observation=observation,
        time_s=float(time_s),mechanical_step_count=pending['mechanical_step_count']+1,phase='SEALED',pending=None,scope_lock=scope)
    result['accepted_bindings'].append(dict(step_id=step_id,receipt_sha256=_sha(receipt),scope=scope,binding=binding,
        joint_process_certified_by_this_module=False,diagnostic_only=scope!=PROCESS_SCOPE,
        control_event_index=pending['pending']['control_event_index']))
    result['source_work_events'].append(dict(step_id=step_id,receipt_sha256=_sha(receipt),work=copy.deepcopy(source_work)))
    result=_seal(result);_validate_state(result,actual_pa)
    return result

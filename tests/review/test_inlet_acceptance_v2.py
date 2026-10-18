"""Manufactured transaction components; no AFSD/NNPS/force/thermal solve."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
from pysph.solver.utils import load

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts'))
import inlet_acceptance_v2 as v2
from afsd_backend.pysph import inlet_protocol as v1
from afsd_backend.pysph import interval_acceptance as crossing
from afsd_backend.pysph.output import checkpoint
from test_inlet_protocol import fixture,properties


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def pose(_time=0.):return dict(origin=np.zeros(3),rotation=np.eye(3))


def setup():
    old,registry,incoming,spec=fixture()
    h=np.arange(48,dtype=float).reshape(2,8,3)*1e-5
    free=np.array([[False,False,False],[True,True,True]])
    obs=crossing.initial_bore_crossing_state(np.column_stack((old.x,old.y,old.z)),old.particle_id,
        old.component_id,old.m,1.,pose(),1.)
    state=v2.initialize_acceptance_state(old,registry,dict(time_s=1.,contact_history=h,free=free,bore_observation=obs),
        dict(frame='WORLD_Z',plane_z_m=.4,unit='m',source='MANUFACTURED_GEOMETRIC_INLET_CONTROL_PLANE'),
        'MANUFACTURED_COMPONENT_FIXTURE_NOT_AFSD')
    return old,state,incoming,spec


def work():
    result=dict(convention=v2.WORK_CONVENTION)
    for term,value in zip(v2.WORK_TERMS,(-.15,.2,.31)):
        result[term]=dict(status='COMPUTED',unit='J',value=value,evidence_sha256=digest(['MANUFACTURED',term,value]))
    return result


def receipt(pa,state,h,free,time,step_id='manufactured-step-1',terms=None):
    terms=work() if terms is None else terms
    return dict(scope=v2.FIXTURE_SCOPE,step_id=step_id,binding=v2.mechanical_binding(pa,state,h,free,time),
        source_work_sha256=digest(terms),
        gates={gate:dict(status='PASS',evidence_sha256=digest(['MANUFACTURED_COMPONENT',gate])) for gate in v2.GATES})


def verifier(row):
    # Deliberately accepts only this fixture's evidence; NEVER a real process.
    return row.get('scope')==v2.FIXTURE_SCOPE and all(
        row['gates'][gate]['evidence_sha256']==digest(['MANUFACTURED_COMPONENT',gate]) for gate in v2.GATES)


def manufactured_update(candidate,pending):
    state=v2.process_state_for(candidate,pending)
    actual=candidate.extract_particles(np.arange(candidate.get_number_of_particles(),dtype=np.uint32))
    actual.x[:]+=.001;actual.temperature[:]+=2.;actual.S[:]+=.7;actual.ep[:]+=.002
    actual.plastic_work[:]+=.01;actual.return_dissipation[:]+=.003
    rows=np.flatnonzero(actual.particle_id==8)
    if len(rows):actual.z[rows[0]]=-.1  # manufactured accepted path, not an imposed AFSD trajectory
    return actual,state['contact_history']+1e-6,state['free']


class InletAcceptanceV2Tests(unittest.TestCase):
    def test_two_insertions_actual_update_restart_sort_and_no_duplicate_sources(self):
        old,state,incoming,spec=setup();parent=copy.deepcopy(state);old_fields=properties(old)
        a,pending,event=v2.propose_insertion(old,state,incoming,spec,pose())
        self.assertEqual(event['allocated_particle_ids'],[8,9])
        self.assertEqual(pending['bore_observation']['crossing_ids'],[])
        self.assertEqual(pending['bore_observation']['birth_events'][0]['eligible_ids'],[8,9])
        self.assertEqual(pending['control_events'][-1]['activated_ids'],[8,9])
        h=v2.process_state_for(a,pending)['contact_history']
        np.testing.assert_array_equal(h[:2],np.asarray(state['contact_history'])[[1,0]])
        np.testing.assert_array_equal(h[2:],np.zeros((2,8,3)))
        actual,h,free=manufactured_update(a,pending)
        terms=work();rec=receipt(actual,pending,h,free,2.,terms=terms)
        accepted=v2.seal_mechanical_acceptance(actual,pending,h,free,2.,rec,terms,verifier,pose)
        self.assertEqual(accepted['protocol_registry']['cumulative'],pending['protocol_registry']['cumulative'])
        self.assertEqual(accepted['bore_observation']['crossing_ids'],[8])
        self.assertEqual(accepted['bore_observation']['crossing_mass_kg'],2.)
        self.assertEqual(accepted['integration_status'],'COMPONENT_ONLY')
        with self.assertRaises(ValueError):v1.restore_registry(json.dumps(pending['protocol_registry']),actual)
        v1.restore_registry(json.dumps(accepted['protocol_registry']),actual)
        permutation=np.array([3,1,0,2],np.uint32)
        actual=actual.extract_particles(permutation)
        with tempfile.TemporaryDirectory(prefix='inlet_acceptance_v2_component_') as directory:
            path=Path(directory)/'component_only_not_afsd.npz'
            checkpoint(path,actual,dict(classification=v2.CLASSIFICATION,inlet_state=accepted,
                restart_allowed=False,physical_process_accepted=False))
            loaded=load(str(path));actual=loaded['arrays']['solid']
            accepted=v2.restore_acceptance_state(json.dumps(loaded['solver_data']['inlet_state']),actual)
        aligned=v2.process_state_for(actual,accepted)
        np.testing.assert_array_equal(aligned['contact_history'],h[permutation])
        second=copy.deepcopy(spec);second.update(event_id='feed-2',time_s=2.)
        b,p2,e2=v2.propose_insertion(actual,accepted,incoming,second,pose())
        self.assertEqual(e2['allocated_particle_ids'],[10,11])
        self.assertEqual(p2['protocol_registry']['cumulative']['mass_kg'],10.)
        self.assertEqual(p2['protocol_registry']['cumulative']['sensible_heat_J'],-800.)
        self.assertEqual(p2['control_events'][-1]['released_ids'],[8])
        self.assertEqual(p2['control_events'][-1]['activated_ids'],[10,11])
        self.assertEqual(p2['control_events'][-1]['mechanical_reaction_work'],'NOT_COMPUTED')
        bstate=v2.process_state_for(b,p2)
        np.testing.assert_array_equal(bstate['contact_history'][:4],aligned['contact_history'])
        index=int(np.flatnonzero(b.particle_id==8)[0]);self.assertTrue(bstate['free'][index,2])
        b,h2,f2=manufactured_update(b,p2)
        accepted2=v2.seal_mechanical_acceptance(b,p2,h2,f2,3.,receipt(b,p2,h2,f2,3.,'manufactured-step-2'),work(),verifier,pose)
        self.assertEqual(accepted2['protocol_registry']['event_ids'],['feed-1','feed-2'])
        self.assertEqual(accepted2['protocol_registry']['high_watermark'],11)
        self.assertEqual(accepted2['protocol_registry']['cumulative']['mass_kg'],10.)
        self.assertEqual(accepted2['bore_observation']['crossing_ids'],[8])
        self.assertEqual(len(accepted2['source_work_events']),2)
        self.assertEqual(state,parent)
        for key,value in old_fields.items():np.testing.assert_array_equal(old.properties[key].get_npy_array(),value)

    def test_birth_below_shoulder_has_no_eligibility_or_crossing(self):
        old,state,inc,spec=setup();inc.z[:]=[-.2,.0];inc.z0[:]=inc.z
        _,pending,_=v2.propose_insertion(old,state,inc,spec,pose())
        self.assertEqual(pending['bore_observation']['eligible_ids'],[3])
        self.assertEqual(pending['bore_observation']['birth_events'][0]['eligible_ids'],[])
        self.assertEqual(pending['bore_observation']['crossing_ids'],[])

    def test_original_crossing_and_initial_eligibility_survive_initialization(self):
        old,registry,inc,spec=fixture();x0=np.column_stack((old.x,old.y,old.z));x1=x0.copy();x1[1,2]=-.1
        obs=crossing.initial_bore_crossing_state(x0,old.particle_id,old.component_id,old.m,0.,pose(),1.)
        obs=crossing.observe_square_bore_crossings(obs,x0,x1,old.particle_id,old.component_id,old.m,0.,1.,pose)['next_state']
        old.z[:]=x1[:,2];registry=v1.initialize_registry(old,'MANUFACTURED_INDEPENDENT_SEED',acknowledgement=v1.ACKNOWLEDGEMENT)
        state=v2.initialize_acceptance_state(old,registry,dict(time_s=1.,contact_history=np.zeros((2,8,3)),
            free=np.ones((2,3),bool),bore_observation=obs),dict(frame='WORLD_Z',plane_z_m=.4,unit='m',source='FIXTURE'),'FIXTURE')
        _,pending,_=v2.propose_insertion(old,state,inc,spec,pose())
        self.assertEqual(pending['bore_observation']['crossing_ids'],[3])
        self.assertEqual(pending['bore_observation']['events'],obs['events'])
        self.assertEqual(pending['bore_observation']['initial_eligible_ids'],[3])

    def test_missing_work_or_zero_placeholder_not_computed_fails_closed(self):
        old,state,inc,spec=setup();pa,pending,_=v2.propose_insertion(old,state,inc,spec,pose())
        actual,h,f=manufactured_update(pa,pending)
        for term in v2.WORK_TERMS:
            for mode in ('missing','uncomputed','wrongunit','unbound'):
                terms=work()
                if mode=='missing':del terms[term]
                elif mode=='uncomputed':terms[term].update(value=0.,status='NOT_COMPUTED')
                elif mode=='wrongunit':terms[term]['unit']='kJ'
                rec=receipt(actual,pending,h,f,2.,terms=terms)
                if mode=='unbound':terms[term]['value']+=1.
                with self.subTest(term=term,mode=mode),self.assertRaises(ValueError):
                    v2.seal_mechanical_acceptance(actual,pending,h,f,2.,rec,terms,verifier,pose)

    def test_missing_or_false_qualified_gate_and_external_verifier_reject(self):
        old,state,inc,spec=setup();pa,pending,_=v2.propose_insertion(old,state,inc,spec,pose())
        actual,h,f=manufactured_update(pa,pending)
        for gate in v2.GATES:
            rec=receipt(actual,pending,h,f,2.);rec['gates'][gate]['status']='NOT_COMPUTED'
            with self.subTest(gate=gate),self.assertRaises(ValueError):
                v2.seal_mechanical_acceptance(actual,pending,h,f,2.,rec,work(),verifier,pose)
        with self.assertRaises(ValueError):
            v2.seal_mechanical_acceptance(actual,pending,h,f,2.,receipt(actual,pending,h,f,2.),work(),None,pose)
        rec=receipt(actual,pending,h,f,2.);rec['scope']=v2.PROCESS_SCOPE
        with self.assertRaisesRegex(ValueError,'joint_inlet_benchmark'):
            v2.seal_mechanical_acceptance(actual,pending,h,f,2.,rec,work(),lambda _:True,pose)

    def test_actual_binding_cannot_reuse_candidate_or_stale_history(self):
        old,state,inc,spec=setup();pa,pending,_=v2.propose_insertion(old,state,inc,spec,pose())
        state_rows=v2.process_state_for(pa,pending)
        stale=receipt(pa,pending,state_rows['contact_history'],state_rows['free'],2.)
        actual,h,f=manufactured_update(pa,pending)
        with self.assertRaisesRegex(ValueError,'bind actual'):
            v2.seal_mechanical_acceptance(actual,pending,h,f,2.,stale,work(),verifier,pose)
        rec=receipt(actual,pending,h,f,2.);h[0,0,0]+=.1
        with self.assertRaisesRegex(ValueError,'bind actual'):
            v2.seal_mechanical_acceptance(actual,pending,h,f,2.,rec,work(),verifier,pose)

    def test_lost_particle_mass_tag_schema_and_control_changes_reject(self):
        for mode in ('lost','mass','tag','schema','free','actuator'):
            old,state,inc,spec=setup();pa,pending,_=v2.propose_insertion(old,state,inc,spec,pose())
            actual,h,f=manufactured_update(pa,pending)
            if mode=='lost':actual=actual.extract_particles(np.array([0,1,2],np.uint32));h=h[:3];f=f[:3]
            elif mode=='mass':actual.m[0]*=2
            elif mode=='tag':actual.tracer_type[0]=999
            elif mode=='schema':actual.remove_property('S')
            elif mode=='free':f[0,0]=not f[0,0]
            else:actual.actuator[0]=1
            rec=receipt(actual,pending,h,f,2.)
            with self.subTest(mode=mode),self.assertRaises(ValueError):
                v2.seal_mechanical_acceptance(actual,pending,h,f,2.,rec,work(),verifier,pose)

    def test_pending_reinsertion_duplicate_event_wrong_time_and_tamper_reject(self):
        old,state,inc,spec=setup();pa,pending,_=v2.propose_insertion(old,state,inc,spec,pose())
        with self.assertRaises(ValueError):v2.propose_insertion(pa,pending,inc,spec,pose())
        with self.assertRaises(ValueError):v2.propose_mechanical_step(pa,pending)
        damaged=copy.deepcopy(pending);damaged['protocol_registry']['cumulative']['mass_kg']=0.
        with self.assertRaises(ValueError):v2.restore_acceptance_state(json.dumps(damaged),pa)
        actual,h,f=manufactured_update(pa,pending)
        accepted=v2.seal_mechanical_acceptance(actual,pending,h,f,2.,receipt(actual,pending,h,f,2.),work(),verifier,pose)
        wrong=copy.deepcopy(spec);wrong['time_s']=2.
        with self.assertRaises(ValueError):v2.propose_insertion(actual,accepted,inc,wrong,pose())
        wrong.update(event_id='feed-2',time_s=2.1)
        with self.assertRaises(ValueError):v2.propose_insertion(actual,accepted,inc,wrong,pose())
        actual.ep[0]+=.001
        with self.assertRaises(ValueError):v2.restore_acceptance_state(json.dumps(accepted),actual)

    def test_mechanical_only_reseal_preserves_inventory_and_signed_work(self):
        old,state,_,_=setup();pa,pending=v2.propose_mechanical_step(old,state)
        actual,h,f=manufactured_update(pa,pending)
        accepted=v2.seal_mechanical_acceptance(actual,pending,h,f,2.,receipt(actual,pending,h,f,2.),work(),verifier,pose)
        self.assertEqual(accepted['protocol_registry']['transaction_count'],0)
        self.assertEqual(accepted['protocol_registry']['cumulative']['mass_kg'],0.)
        self.assertEqual(accepted['source_work_events'][0]['work']['pressure_flow_work']['value'],-.15)
        self.assertEqual(accepted['source_events'],[])

    def test_diagnostic_two_insertion_chain_needs_no_prior_joint_but_is_not_process(self):
        # Receipt payloads here remain manufactured; this tests the API path
        # a future real benchmark can use, not a real coupled benchmark.
        pa,state,inc,spec=setup()
        for step in (1,2):
            spec=copy.deepcopy(spec);spec.update(event_id='diagnostic-feed-'+str(step),time_s=float(step))
            candidate,pending,_=v2.propose_insertion(pa,state,inc,spec,pose())
            pa,h,f=manufactured_update(candidate,pending)
            rec=receipt(pa,pending,h,f,float(step+1),'diagnostic-fixture-step-'+str(step))
            rec['scope']=v2.DIAGNOSTIC_SCOPE
            self.assertNotIn('joint_inlet_benchmark',rec['gates'])
            state=v2.seal_mechanical_acceptance(pa,pending,h,f,float(step+1),rec,work(),
                lambda row:row['scope']==v2.DIAGNOSTIC_SCOPE,pose)
            state=v2.restore_acceptance_state(json.dumps(state),pa)
            self.assertEqual(state['scope_lock'],v2.DIAGNOSTIC_SCOPE)
            self.assertFalse(state['process_acceptance_by_this_module'])
            self.assertTrue(state['accepted_bindings'][-1]['diagnostic_only'])
        self.assertEqual(state['protocol_registry']['transaction_count'],2)
        self.assertEqual(state['protocol_registry']['cumulative']['mass_kg'],10.)
        self.assertEqual(state['integration_status'],'COMPONENT_ONLY')

    def test_scope_history_never_upgrades_or_mixes(self):
        scopes=(v2.FIXTURE_SCOPE,v2.DIAGNOSTIC_SCOPE,v2.PROCESS_SCOPE)
        for first in scopes:
            pa,state,_,_=setup();candidate,pending=v2.propose_mechanical_step(pa,state)
            pa,h,f=manufactured_update(candidate,pending)
            rec=receipt(pa,pending,h,f,2.);rec['scope']=first
            if first==v2.PROCESS_SCOPE:
                rec['gates']['joint_inlet_benchmark']=dict(status='PASS',evidence_sha256=digest('MANUFACTURED_API_GATE'))
            state=v2.seal_mechanical_acceptance(pa,pending,h,f,2.,rec,work(),lambda _:True,pose)
            candidate,pending=v2.propose_mechanical_step(pa,state);pa,h,f=manufactured_update(candidate,pending)
            for next_scope in scopes:
                if next_scope==first:continue
                changed=receipt(pa,pending,h,f,3.,'manufactured-step-2');changed['scope']=next_scope
                changed['gates']['joint_inlet_benchmark']=dict(status='PASS',evidence_sha256=digest('MANUFACTURED_API_GATE'))
                with self.subTest(first=first,next_scope=next_scope),self.assertRaisesRegex(ValueError,'scope history is locked'):
                    v2.seal_mechanical_acceptance(pa,pending,h,f,3.,changed,work(),lambda _:True,pose)


if __name__=='__main__':unittest.main(verbosity=2)

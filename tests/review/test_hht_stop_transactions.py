"""Independent state/receipt guards; manufactured components, no 3D scout.

Candidate tests use one small existing unit-fixture stage. Restart tests mock
only the instantaneous internal-force evaluator to isolate the transaction;
their constructed state is explicitly NOT a physical event trajectory.
"""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import hht_stop_stage as stage,hht_stop_flow as flow
from afsd_backend.pysph.hht_coupled import _hash,_boundary
from afsd_backend.pysph.updated_trial import particle_hash
from test_hht_coupled import fixture,body


def rehash(value):
    value.pop('state_hash',None);value['state_hash']=_hash(value);return value


def copies(pa):
    return {name:prop.get_npy_array().copy() for name,prop in pa.properties.items()}


class StopTransactions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old,cls.p=fixture();cls.p['hht_weight']=.9
        cls.old.u[:]=.02;cls.bc=staticmethod(body([2.,-3.,5.]))
        cls.parent=stage.initialize(cls.old,cls.p,cls.bc)
        cls.candidate,cls.trial,cls.diag=stage.advance_private(cls.old,cls.parent,cls.p,cls.bc,1e-8,None)

    def certify(self,trial=None,diag=None,candidate=None,parent=None):
        return flow.certify_interval(self.old,self.parent if parent is None else parent,
            self.candidate if candidate is None else candidate,self.trial if trial is None else trial,
            self.diag if diag is None else diag,self.p,self.bc,None)

    def assert_fields(self,pa,before):
        self.assertEqual(set(pa.properties),set(before))
        for key,value in before.items():np.testing.assert_array_equal(pa.properties[key].get_npy_array(),value,err_msg=key)

    def test_01_private_checkpoint_cannot_restore_or_advance(self):
        before=copies(self.candidate)
        for accepted,restart in ((False,False),(True,False),(False,True)):
            state=copy.deepcopy(self.trial);state.update(accepted=accepted,restart_allowed=restart);rehash(state)
            with self.subTest(accepted=accepted,restart=restart):
                with self.assertRaises(ValueError):stage.restore_checkpoint(self.candidate,state,self.p)
                with self.assertRaises(ValueError):stage.advance_private(self.candidate,state,self.p,self.bc,1e-8,None)
        self.assert_fields(self.candidate,before)

    def test_02_valid_certification_is_detached_and_preserves_all_fields(self):
        old_before=copies(self.old);candidate_before=copies(self.candidate)
        trial_hash=self.trial['state_hash'];diag_hash=_hash(self.diag)
        accepted,report=self.certify()
        self.assertTrue(accepted['accepted']);self.assertTrue(accepted['restart_allowed'])
        self.assertEqual(accepted['parent_state_hash'],self.parent['state_hash'])
        self.assertEqual(report['accepted_state_hash'],accepted['state_hash'])
        self.assertEqual(self.trial['state_hash'],trial_hash);self.assertEqual(_hash(self.diag),diag_hash)
        self.assertFalse(self.trial['accepted']);self.assert_fields(self.old,old_before);self.assert_fields(self.candidate,candidate_before)

    def test_03_particle_parent_and_diagnostic_candidate_bindings(self):
        changed=self.candidate.extract_particles(np.arange(27,dtype=np.uint32));changed.temperature[0]+=1.
        with self.assertRaises(ValueError):self.certify(candidate=changed)
        parent=copy.deepcopy(self.parent);parent['time_s']+=1.;rehash(parent)
        with self.assertRaises(ValueError):self.certify(parent=parent)
        diag=copy.deepcopy(self.diag);diag['candidate_state_hash']='not-this-candidate'
        with self.assertRaises(ValueError):self.certify(diag=diag)
        diag=copy.deepcopy(self.diag);diag['accepted_parent_hash']='not-this-parent'
        with self.assertRaises(ValueError):self.certify(diag=diag)

    def test_04_time_before_must_be_the_accepted_parent_time(self):
        diag=copy.deepcopy(self.diag);diag['time_before_s']+=1e-5
        with self.assertRaises(ValueError):self.certify(diag=diag)

    def test_05_candidate_end_must_equal_parent_time_plus_dt(self):
        # Rehash a semantically valid but chronologically inconsistent receipt;
        # the diagnostic agrees with the wrong endpoint time, so merely checking
        # the two against each other is insufficient.
        trial=copy.deepcopy(self.trial);trial['time_s']+=1e-5;rehash(trial)
        diag=copy.deepcopy(self.diag);diag['time_after_s']=trial['time_s'];diag['candidate_state_hash']=trial['state_hash']
        with self.assertRaises(ValueError):self.certify(trial=trial,diag=diag)

    def test_06_state_interval_integrals_must_match_the_strict_recomputed_path(self):
        for key in ('actual_contact_impulse_Ns','actual_contact_material_work_J','raw_friction_J',
                    'delta_impulse_Ns','delta_material_work_J'):
            with self.subTest(field=key):
                trial=copy.deepcopy(self.trial);trial['last_interval'][key].flat[0]+=1e-6;rehash(trial)
                diag=copy.deepcopy(self.diag);diag['candidate_state_hash']=trial['state_hash']
                with self.assertRaises(ValueError):self.certify(trial=trial,diag=diag)

    def test_07_diagnostic_integral_and_raw_residual_guard(self):
        diag=copy.deepcopy(self.diag);diag['contact']['actual_impulse_Ns'][0,0]+=1e-6
        with self.assertRaises(ValueError):self.certify(diag=diag)
        for key,value in [('velocity_residual_m_s',1e-8),('position_residual_m',1e-11),
                          ('density_relative_residual',1e-9),('temperature_residual_K',5e-11)]:
            diag=copy.deepcopy(self.diag);diag['residuals'][key]=value
            with self.subTest(field=key):
                with self.assertRaises(ValueError):self.certify(diag=diag)

    def restart_fixture(self):
        pa,p=fixture()  # Intentionally retain the supported default weight.
        pa.S[:]=np.tile(np.diag([100.,-50.,-50.]).ravel(),27)
        pa.temperature[:]=np.linspace(400.,401.,27);pa.tracer_type[:]=np.arange(27)%3
        pa.u[1]=.1
        other=np.zeros((27,3));other[0]=pa.m[0]*np.array([-.4,0.,-1.]);other[1,2]=-pa.m[1]
        normal=np.array([0.,0.,1.])
        def boundary_for(outgoing=None,alter=None):
            def boundary(private,at):
                rows=[]
                for idx in (0,1):
                    fn=private.m[idx];cap=.2*fn;cid='contact:'+str(idx)
                    direction=np.array([1.,0.,0.])
                    if outgoing and cid in outgoing:direction=outgoing[cid]['outgoing_direction']
                    row=dict(particle_index=idx,contact_id=cid,normal=normal.copy(),normal_force_N=fn,
                        tangential_force_N=-cap*direction,wall_velocity_m_s=np.zeros(3),
                        wall_temperature_K=410.,conductance_W_K=.01,penetration_m=1e-6,
                        normal_potential_J=.5*fn*1e-6,normal_stiffness_N_m=fn/1e-6)
                    if alter=='nontarget_heat' and idx==1:row['wall_temperature_K']+=1.
                    if alter=='nontarget_tangent' and idx==1:row['tangential_force_N']=-row['tangential_force_N']
                    if alter=='target_heat' and idx==0:row['conductance_W_K']*=2.
                    rows.append(row)
                force=other.copy()
                if alter=='noncontact_force':force[2,0]+=1e-7
                return dict(force_definition='NONCONTACT_APPLIED_FORCE_ONLY',force_N=force,contacts=rows)
            return boundary
        left_boundary=_boundary(pa,boundary_for(),.5,p)
        canonical=np.tile([-3.,.1,.2],(27,1))
        left=stage._state(pa,p,.5,canonical,left_boundary['total_boundary_force_N'],left_boundary,
                          parent='manufactured-parent',dt=.125,interval=dict(transaction_fixture_only=True))
        context={'contact:0':dict(endpoint='left_zero',incoming_direction=np.array([1.,0.,0.]))}
        return pa,p,left,context,boundary_for

    def test_08_empty_restart_cannot_reset_acceleration(self):
        pa,p,left,context,bc=self.restart_fixture();before=copies(pa)
        with patch.object(flow,'_force',side_effect=AssertionError('empty context reached force calculation')):
            with self.assertRaises(ValueError):flow.restart_sliding(pa,left,p,{},bc,np.ones((27,3)),{})
        self.assert_fields(pa,before)

    def test_09_restart_preserves_every_field_and_uses_old_canonical_history(self):
        pa,p,left,context,bc=self.restart_fixture();before=copies(pa);before_left=_hash(left)
        internal=np.zeros((27,3));spatial=dict(classification='MOCK_INTERNAL_FORCE_TRANSACTION_ONLY')
        with patch.object(flow,'_force',return_value=(internal,spatial)):
            state,outgoing,event=flow.restart_sliding(pa,left,p,context,bc,np.ones((27,3)),{'contact:0':np.zeros(3)})
        self.assert_fields(pa,before);self.assertEqual(_hash(left),before_left)
        self.assertEqual(state['time_s'],left['time_s']);self.assertEqual(state['parent_state_hash'],left['state_hash'])
        np.testing.assert_array_equal(event['old_canonical_acceleration_m_s2'],left['acceleration_m_s2'])
        self.assertFalse(np.array_equal(left['acceleration_m_s2'],left['force_N']/pa.m[:,None]))
        expected=np.zeros((27,3));expected[0,0]=-.2;expected[1,0]=-.2
        np.testing.assert_allclose(state['acceleration_m_s2'],expected,rtol=0.,atol=1e-14)
        old_acc=left['acceleration_m_s2'];d=(1.-.9)**2/4.
        jump=.5*.125**2*d*float(np.sum(pa.m[:,None]*(expected*expected-old_acc*old_acc)))
        self.assertAlmostEqual(event['signed_acceleration_storage']['signed_jump_J'],jump,delta=1e-18)
        np.testing.assert_array_equal(event['physical_impulse_Ns'],np.zeros((27,3)))
        self.assertEqual(event['physical_work_J'],0.);self.assertEqual(event['physical_heat_J'],0.)
        self.assertEqual(event['physical_state_hash_before'],event['physical_state_hash_after'])
        self.assertEqual(state['last_interval'],left['last_interval'])
        np.testing.assert_array_equal(outgoing['contact:0']['outgoing_direction'],[-1.,0.,0.])

    def test_10_restart_rejects_hidden_nontarget_or_other_boundary_changes(self):
        for alteration in ('nontarget_heat','nontarget_tangent','target_heat','noncontact_force'):
            with self.subTest(alteration=alteration):
                pa,p,left,context,bc=self.restart_fixture();before=copies(pa)
                factory=lambda outgoing:bc(outgoing,alteration)
                with patch.object(flow,'_force',return_value=(np.zeros((27,3)),{})):
                    with self.assertRaises(ValueError):flow.restart_sliding(pa,left,p,context,factory,np.ones((27,3)),{'contact:0':np.zeros(3)})
                self.assert_fields(pa,before)


if __name__=='__main__':unittest.main(verbosity=2)

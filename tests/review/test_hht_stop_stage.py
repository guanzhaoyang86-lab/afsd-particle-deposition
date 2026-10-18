"""Event transaction invariants, no process-physics certification."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from afsd_backend.pysph import hht_curved_coupled as old,hht_stop_stage as new
from afsd_backend.pysph.hht_stop_flow import certify_interval,locate_duration
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.hht_coupled import _hash,_boundary
from afsd_backend.pysph.mechanics import positions,velocities
from test_hht_coupled import fixture,body
import probe_hht_curved as driver
from hht_stop_boundary import boundary_for


class StopStage(unittest.TestCase):
    def test_private_trial_not_restartable_and_strict_eventfree_reduction(self):
        pa,p=fixture();pa.u[:]=.02;bc=body([2.,-3.,5.]);h=1e-8
        prior=old.initialize(pa,p,bc);state=new.migrate_published(pa,prior,p)
        before=particle_hash(pa)
        a,sa,d=old.advance(pa,prior,p,bc,h,None)
        b,sb,db=new.advance_private(pa,state,p,bc,h,None)
        self.assertFalse(sb['accepted']);self.assertFalse(sb['restart_allowed'])
        with self.assertRaises(ValueError):new.restore_checkpoint(b,sb,p)
        with self.assertRaises(ValueError):new.advance_private(b,sb,p,bc,h,None)
        for key in pa.properties:np.testing.assert_array_equal(getattr(a,key),getattr(b,key))
        certified,detail=certify_interval(pa,state,b,sb,db,p,bc,None)
        restored,verified=new.restore_checkpoint(b,certified,p)
        self.assertEqual(verified['state_hash'],certified['state_hash'])
        self.assertEqual(particle_hash(restored),particle_hash(b));self.assertEqual(before,particle_hash(pa))
        self.assertTrue(detail['accepted'])
        np.testing.assert_allclose(velocities(b)-velocities(pa),h*np.array([2.,-3.,5.])[None,:]*np.ones((27,1)),atol=2e-14,rtol=0)

    def test_private_promotion_and_changed_integrals_rejected(self):
        pa,p=fixture();bc=body([0.,0.,0.]);state=new.initialize(pa,p,bc)
        candidate,cs,d=new.advance_private(pa,state,p,bc,1e-8,None)
        with self.assertRaises(ValueError):certify_interval(pa,state,candidate,cs,d,p,bc,None,
            {'bad':dict(endpoint='private_continuation')})
        changed=copy.deepcopy(d);changed['contact']['friction_raw_J'][0]+=1e-9
        with self.assertRaises(ValueError):certify_interval(pa,state,candidate,cs,changed,p,bc,None)
        damaged=copy.deepcopy(cs);damaged['acceleration_m_s2'][0,0]+=1
        with self.assertRaises(ValueError):certify_interval(pa,state,candidate,damaged,d,p,bc,None)

    def test_boundary_eventfree_matches_frozen_physics_and_labels_observational(self):
        pa,p=driver.setup();before=particle_hash(pa)
        a=_boundary(pa,driver.boundary_for(p),1e-6,p)
        b=_boundary(pa,boundary_for(p),1e-6,p)
        self.assertEqual(_hash(a),_hash(b));self.assertEqual(before,particle_hash(pa))
        pa.tracer_type[:]=7;pa.wholebar_tag[:]=99
        self.assertEqual(_hash(b),_hash(_boundary(pa,boundary_for(p),1e-6,p)))

    def test_duration_solver_locates_without_velocity_projection(self):
        # Independent force law du/dt=-3, actual stop1/3. The duration solver
        # receives only incoming state and callable candidates, not that root.
        pa,p=fixture();pa.u[:]=1.;pa.v[:]=0.;pa.w[:]=0.
        row=dict(contact_id='face:0',particle_index=0,normal=np.array([0.,1.,0.]),wall_velocity_m_s=np.zeros(3))
        parent=dict(boundary_snapshot=dict(contacts=[row]),acceleration_m_s2=np.tile([-3.,0.,0.],(27,1)))
        calls=[]
        def evaluate(h):
            calls.append(h);result=pa.extract_particles(np.arange(27,dtype=np.uint32));result.u[:]=1.-3.*h
            state=dict(boundary_snapshot=dict(contacts=[row]),acceleration_m_s2=parent['acceleration_m_s2'],
                state_hash=str(h),time_s=h)
            diag=dict(residuals=dict(velocity_residual_m_s=0.),evaluations=1)
            return result,state,diag
        context={'face:0':dict(endpoint='private_continuation',incoming_direction=np.array([1.,0.,0.]))}
        result,state,d,receipt=locate_duration(evaluate,pa,parent,context,.5)
        self.assertAlmostEqual(state['time_s'],1./3.,places=15)
        self.assertEqual(calls[0],.5);self.assertEqual(len(calls),2)
        self.assertFalse(receipt['physical_state_was_projected']);self.assertFalse(receipt['reference_time_was_used'])
        np.testing.assert_array_equal(result.u,1.-3.*state['time_s'])
        np.testing.assert_array_equal(pa.u,1.)


if __name__=='__main__':unittest.main(verbosity=2)

"""New-schema physical transaction tests, not process/temporal acceptance."""
from pathlib import Path
import copy
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from test_hht_coupled import fixture,empty,body,heat_reference
from afsd_backend.pysph import hht_coupled as old,hht_curved_coupled as new
from afsd_backend.pysph.mechanics import positions,velocities
from afsd_backend.pysph.updated_trial import particle_hash
import probe_hht_curved as driver
from afsd_backend.pysph.hht_curved_contact import CurvedContactError


class CurvedCoupled(unittest.TestCase):
    def test_constant_acceleration_and_new_checkpoint(self):
        pa,p=fixture();pa.u[:]=.03;pa.v[:]=-.02;pa.w[:]=.01
        acceleration=np.array([2.,-3.,5.]);bc=body(acceleration);h=2e-8
        state=new.initialize(pa,p,bc);before=particle_hash(pa)
        after,ns,d=new.advance(pa,state,p,bc,h,None)
        np.testing.assert_allclose(positions(after),positions(pa)+h*velocities(pa)+.5*h*h*acceleration,atol=2e-18,rtol=0)
        np.testing.assert_allclose(velocities(after),velocities(pa)+h*acceleration,atol=2e-14,rtol=0)
        self.assertEqual(particle_hash(pa),before)
        restored,rs=new.restore_checkpoint(pa,ns,p)
        self.assertEqual(particle_hash(restored),particle_hash(after));self.assertEqual(rs['state_hash'],ns['state_hash'])
        with self.assertRaises(ValueError):old.restore_checkpoint(pa,ns,p)
        bad=copy.deepcopy(ns);bad['force_N'][0,0]+=1
        with self.assertRaises(ValueError):new.restore_checkpoint(pa,bad,p)
        self.assertLess(float(np.max(np.abs(d['kinetic']['kinetic_identity_residual_J']))),1e-12)

    def test_no_contact_exact_baseline_reduction(self):
        pa,p=fixture();pa.u[:]=.025;bc=body([2.,0.,-3.]);h=1e-8
        a,sa,_=old.advance(pa,old.initialize(pa,p,bc),p,bc,h)
        b,sb,diag=new.advance(pa,new.initialize(pa,p,bc),p,bc,h,None)
        for name in pa.properties:np.testing.assert_array_equal(a.properties[name].get_npy_array(),b.properties[name].get_npy_array())
        np.testing.assert_array_equal(sa['force_N'],sb['force_N'])
        np.testing.assert_array_equal(diag['contact']['delta_impulse_Ns'],0.)

    def test_heat_reference_and_no_algorithmic_heat(self):
        pa,p=fixture();pa.temperature[:]=np.linspace(400.,430.,27);h=1e-4;z=np.zeros(27)
        reference=heat_reference(pa,p,h,z,z,z)
        after,_,diag=new.advance(pa,new.initialize(pa,p,empty),p,empty,h,None)
        np.testing.assert_allclose(after.temperature,reference,atol=2e-12,rtol=0)
        self.assertFalse(diag['algorithmic_terms_are_material_heat'])

    def test_explicit_migration_preserves_every_physical_array(self):
        pa,p=driver.setup();bc=driver.boundary_for(p);state=old.initialize(pa,p,bc)
        before=particle_hash(pa);before_state=copy.deepcopy(state)
        migrated,event=new.migrate_diagnostic(pa,state,p)
        self.assertEqual(before,particle_hash(pa));self.assertEqual(state['state_hash'],before_state['state_hash'])
        self.assertEqual(migrated['parent_state_hash'],state['state_hash'])
        self.assertEqual(event['physical_energy_jump_J'],0.)
        for key in pa.properties:np.testing.assert_array_equal(migrated['particle_snapshot'][key],state['particle_snapshot'][key])
        np.testing.assert_array_equal(migrated['acceleration_m_s2'],state['acceleration_m_s2'])
        np.testing.assert_array_equal(migrated['force_N'],state['force_N'])
        with self.assertRaises(ValueError):new.advance(pa,state,p,bc,1e-8,driver.slow_gap_certificate)

    def test_actual_small_coupled_step_and_passive_tags(self):
        pa,p=driver.setup();bc=driver.boundary_for(p);h=1e-8
        state=new.initialize(pa,p,bc);before=particle_hash(pa)
        after,ns,diag=new.advance(pa,state,p,bc,h,driver.slow_gap_certificate)
        driver.checks(after,pa,pa,ns,diag,p)
        tag=pa.extract_particles(np.arange(27,dtype=np.uint32));tag.tracer_type[:]=2;tag.wholebar_tag[:]=7
        result,_,_=new.advance(tag,new.initialize(tag,p,bc),p,bc,h,driver.slow_gap_certificate)
        for key in ('x','y','z','u','v','w','rho','S','temperature','ep'):
            np.testing.assert_array_equal(getattr(result,key),getattr(after,key))
        self.assertEqual(before,particle_hash(pa))
        self.assertGreater(float(diag['contact']['friction_raw_J'].sum()),0.)
        self.assertFalse(diag['classical_hht_stability_certified'])
        altered=copy.deepcopy(ns);altered['last_interval']['raw_friction_J'][0]+=1e-9
        with self.assertRaises(ValueError):new.restore_checkpoint(pa,altered,p)

    def test_early_contact_failure_retains_identity_and_private_candidate(self):
        pa,p=driver.setup();bc=driver.boundary_for(p);state=new.initialize(pa,p,bc)
        before=particle_hash(pa)
        error=CurvedContactError('synthetic exact event audit',dict(contact_id='0:4',zero_slip_audit={'synthetic':True}))
        with patch.object(new,'integrate_contact',side_effect=error):
            with self.assertRaises(new.CoupledStageError) as caught:
                new.advance(pa,state,p,bc,1e-8,driver.slow_gap_certificate)
        diagnostic=caught.exception.diagnostics
        self.assertEqual(diagnostic['contact_id'],'0:4')
        self.assertFalse(diagnostic['pending_unequilibrated_candidate']['accepted'])
        self.assertFalse(diagnostic['restart_allowed'])
        self.assertEqual(particle_hash(pa),before)


if __name__=='__main__':unittest.main(verbosity=2)

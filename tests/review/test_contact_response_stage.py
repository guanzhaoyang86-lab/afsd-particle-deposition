"""Private numerical preconditioner plumbing; NOT a process certificate."""
from pathlib import Path
import ast
import copy
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from afsd_backend.pysph import hht_stop_stage as old,hht_response_stage as new
from afsd_backend.pysph.contact_response_preconditioner import assemble,STRATEGY
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.mechanics import velocities
from test_hht_coupled import fixture,body
import probe_hht_curved as driver


def synthetic(private=False):
    n=np.array([0.,0.,1.]);u=np.array([2.,0.,0.]);sign=-1 if private else 1
    end=dict(contact_id='plane:1',particle_index=0,particle_id=1,normal=n,
        normal_force_N=10.,normal_stiffness_N_m=20.,wall_velocity_m_s=np.zeros(3),
        tangential_force_N=np.array([-sign*3.,0.,0.]))
    node=dict(fraction=.5,weight=1.,normal_force_N=10.,tangential_relative_velocity_m_s=u,
        tangential_force_N=end['tangential_force_N'].copy())
    cert=(dict(mode='private_continuation',incoming_direction=np.array([-1.,0.,0.]),zero_slip_bound_m_s=1e-13)
        if private else None)
    contact=dict(contacts=[dict(end,quadrature=[node],event_path_certificate=cert)])
    return contact,dict(contacts=[end]),u[None,:],dict(friction=.3,hht_weight=.9)


class ResponseStage(unittest.TestCase):
    def test_literal_discrete_block_and_private_sign(self):
        # gamma=.6,beta=.3025, h=.2 and one weight1 node s=.5:
        # full FORCE root derivative = -.00589*Kx -.10925*Kv.
        for private in (False,True):
            c,b,v,p=synthetic(private);before=copy.deepcopy(c)
            matrix,detail=assemble(c,b,v,.2,p);sigma=-1 if private else 1
            kx=np.array([[0.,0.,sigma*6.],[0.,0.,0.],[0.,0.,-20.]])
            kv=np.diag([0.,-sigma*1.5,0.])
            np.testing.assert_allclose(matrix,-.00589*kx-.10925*kv,atol=2e-14,rtol=0)
            self.assertFalse(detail['continuous_cross_zero_jacobian_certified'])
            self.assertFalse(detail['can_publish']);self.assertEqual(detail['strategy'],STRATEGY)
            np.testing.assert_array_equal(c['contacts'][0]['quadrature'][0]['tangential_force_N'],
                before['contacts'][0]['quadrature'][0]['tangential_force_N'])

    def test_response_rejects_zero_override_or_mismatched_force(self):
        c,b,v,p=synthetic(True);c['contacts'][0]['quadrature'][0]['tangential_relative_velocity_m_s']=np.zeros(3)
        with self.assertRaises(ValueError):assemble(c,b,v,.2,p)
        c,b,v,p=synthetic();c['contacts'][0]['quadrature'][0]['tangential_force_N']*=.5
        with self.assertRaises(ValueError):assemble(c,b,v,.2,p)
        c,b,v,p=synthetic(True);c['contacts'][0]['event_path_certificate']['mode']='right_zero'
        with self.assertRaises(ValueError):assemble(c,b,v,.2,p)

    def test_actual_nonlinear_evaluator_ast_unchanged(self):
        def evaluation(path,name):
            tree=ast.parse(path.read_text(encoding='utf-8'))
            advance=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name==name)
            inner=next(node for node in advance.body if isinstance(node,ast.FunctionDef) and node.name=='evaluate')
            return ast.dump(inner,include_attributes=False)
        folder=ROOT/'src/afsd_backend/pysph'
        self.assertEqual(evaluation(folder/'hht_stop_stage.py','advance_private'),
            evaluation(folder/'hht_response_stage.py','advance_private_response'))
        self.assertEqual(new.MAX_EVALUATIONS,24)

    def test_empty_contact_exact_baseline_reduction_and_nonrestartable(self):
        pa,p=fixture();bc=body([2.,-3.,5.]);state=old.initialize(pa,p,bc);before=particle_hash(pa)
        a,sa,da=old.advance_private(pa,state,p,bc,1e-8,None)
        b,sb,db=new.advance_private_response(pa,state,p,bc,1e-8,None)
        for name in pa.properties:np.testing.assert_array_equal(getattr(a,name),getattr(b,name))
        self.assertEqual(sa['state_hash'],sb['state_hash']);self.assertEqual(before,particle_hash(pa))
        self.assertFalse(sb['accepted']);self.assertFalse(sb['restart_allowed']);self.assertFalse(db['can_publish'])
        with self.assertRaises(ValueError):old.restore_checkpoint(b,sb,p)
        np.testing.assert_array_equal(velocities(a),velocities(b))

    def test_actual_small_stage_keeps_old_raw_gates_and_parent(self):
        pa,p=driver.setup();bc=driver.boundary_for(p);state=old.initialize(pa,p,bc);before=particle_hash(pa)
        candidate,ns,diag=new.advance_private_response(pa,state,p,bc,1e-8,driver.slow_gap_certificate)
        driver.checks(candidate,pa,pa,ns,diag,p)
        self.assertEqual(before,particle_hash(pa));self.assertFalse(ns['restart_allowed'])
        self.assertEqual(diag['numerical_preconditioner'],STRATEGY)
        self.assertGreater(len(diag['contact_response_history']),0)
        self.assertGreater(sum(diag['contact']['friction_raw_J']),0)
        self.assertLessEqual(diag['residuals']['velocity_residual_m_s'],1e-9)
        self.assertLessEqual(diag['residuals']['temperature_residual_K'],5e-12)
        self.assertFalse(diag['algorithmic_terms_are_material_heat'])
        changed=pa.extract_particles(np.arange(27,dtype=np.uint32));changed.actuator[0]=1
        with self.assertRaises(ValueError):new.advance_private_response(changed,state,p,bc,1e-8,driver.slow_gap_certificate)


if __name__=='__main__':unittest.main(verbosity=2)

"""Root orchestration regression, not process transport evidence."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
from unittest import mock
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import run_process_event as runner
from test_updated_trial import cloud
from afsd_backend.pysph.process_thermal_trial import ProcessThermalTrial,acceptance
from afsd_backend.pysph.updated_trial import particle_hash
from afsd_backend.pysph.stage_event_preconditioner import StageEventPreconditioner
from afsd_backend.pysph.stage_structural_root import solve_stage
from afsd_backend.pysph.event_response import response_from_trace


class ProcessEventRunnerTests(unittest.TestCase):
    def test_coupled_spring_and_actual_reversal_kick_literal_root_and_energy(self):
        from afsd_backend.pysph.event_response import response
        from afsd_backend.pysph.contact_kick import kick
        from afsd_backend.pysph.structural_mobility import FrozenElastic
        initial_x=19/36
        class Spring:
            def candidate(self):
                return dict(position_m=np.array([[initial_x,0.,0.]]),velocity_m_s=np.array([[1.,0.,0.]]),
                    mean_velocity_m_s=np.array([[1.,0.,0.]]),density_kg_m3=np.ones(1),temperature_K=np.full(1,300.))
            def evaluate(self,c,dt):
                f=-4*initial_x-2*dt*c['mean_velocity_m_s'][0,0]
                actual=kick([1.,0.,0.],1.,[f,0.,-1.],[[0.,0.,1.]],[1.],[[0.,0.,0.]],1.,dt,
                    free_axes=[1.,0.,0.],prescribed_velocity=[1.,0.,0.])
                mean=sum(a['dt_s']*.5*(a['v_before']+a['v_after']) for a in actual['intervals'])/dt
                mapped={k:v.copy() for k,v in c.items()};mapped['velocity_m_s']=actual['velocity'][None,:]
                mapped['mean_velocity_m_s']=mean[None,:];mapped['position_m']=np.array([[initial_x,0.,0.]])+dt*mean
                return dict(map=mapped,contact=actual)
        local=response([1.,0.,0.],1.,[-2.,0.,-1.],[[0.,0.,1.]],[1.],[[0.,0.,0.]],1.,1.,
            free_axes=[1.,0.,0.],prescribed_velocity=[1.,0.,0.])
        self.assertTrue(local['available'])
        elastic=FrozenElastic(np.zeros((9,3)),np.diag([4.,0.,0.]),[1.],[1.],
            np.array([[initial_x,0.,0.]]),np.array([7]),1.,1.,0)
        pre=StageEventPreconditioner.from_response_blocks(elastic,local['B'][None,:,:],local['C'][None,:,:],
            1.,np.array([[1,0,0]]),[1.],{'fixture':'INDEPENDENT_SPRING_KICK_NOT_AFSD'})
        r=solve_stage(Spring(),pre,1.);c=r['candidate'];actual=r['evaluation']['contact']
        np.testing.assert_allclose(c['position_m'][0,0],17/36,atol=1e-12,rtol=0)
        np.testing.assert_allclose(c['velocity_m_s'][0,0],-2/3,atol=1e-9,rtol=0)
        np.testing.assert_allclose(c['mean_velocity_m_s'][0,0],-1/18,atol=1e-9,rtol=0)
        q=float(np.sum(actual['friction_heat_J']))
        self.assertAlmostEqual(q,7/18,delta=1e-9)
        change=.5*(c['velocity_m_s'][0,0]**2-1.)+2*(c['position_m'][0,0]**2-initial_x**2)
        self.assertAlmostEqual(change,-7/18,delta=1e-9)
        self.assertAlmostEqual(change+q,0.,delta=1e-9)

    def test_literal_budget_and_no_new_physical_steps_from_components(self):
        r=runner.cost_budget(dict(preconditioner_build_wall_seconds=2.,trial_proposal_journal_wall_seconds=4.))
        self.assertTrue(r['allowed']);self.assertEqual(r['forecast_wall_seconds'],225.)
        r=runner.cost_budget(dict(preconditioner_build_wall_seconds=2.,trial_proposal_journal_wall_seconds=10.))
        self.assertFalse(r['allowed']);self.assertEqual(r['forecast_wall_seconds'],513.)
        for value in (0.,-1.,float('nan'),float('inf')):
            with self.assertRaises(ValueError):runner.cost_budget(dict(preconditioner_build_wall_seconds=2.,trial_proposal_journal_wall_seconds=value))

    def test_manifest_binds_rejected_anchor_and_accepted_parent_separately(self):
        with tempfile.TemporaryDirectory(prefix='afsd_event_receipt_') as folder:
            root=Path(folder);path=root/'manifest.json';sources={};watched={};parent={'sha256':'accepted'}
            anchor={'sha256':'rejected-only'};preflight={}
            report=dict(status='COMPLETE_FROZEN_EVENT_RESPONSE_NOT_PHYSICAL_STATE',settings=runner.SETTINGS,
                classification=runner.CLASSIFICATION,accepted_steps=0,production_authorized=False,
                source_sha256=sources,configuration_sha256='cfg',input_sha256=watched,parent_reference=parent,
                numerical_anchor=anchor,preflight=preflight,artifact_sha256={})
            path.write_text(json.dumps(report),encoding='utf-8')
            with mock.patch.object(runner,'ROOT',root):
                runner.checked_report(path,sources,'cfg',watched,parent,anchor,preflight,report['status'])
                for badparent,badanchor in ((anchor,anchor),(parent,parent)):
                    with self.assertRaisesRegex(ValueError,'parent/anchor'):
                        runner.checked_report(path,sources,'cfg',watched,badparent,badanchor,preflight,report['status'])
                with self.assertRaises(ValueError):runner.checked_report(path,sources,'different',watched,parent,anchor,preflight,report['status'])

    def test_partial_response_cannot_be_used_for_root(self):
        with tempfile.TemporaryDirectory(prefix='afsd_partial_response_') as folder:
            root=Path(folder);path=root/'manifest.json';pa,p=cloud()
            r=dict(status='BLOCKED_NO_NEW_ACCEPTED_STATE',settings=runner.SETTINGS,classification=runner.CLASSIFICATION,
                accepted_steps=0,production_authorized=False)
            path.write_text(json.dumps(r),encoding='utf-8')
            with mock.patch.object(runner,'ROOT',root):
                with self.assertRaises(ValueError):runner.response_arrays(path,pa,{},'cfg',{}, {},{}, {})

    def test_real_zero_equilibrium_joint_chain_preserves_accepted_history(self):
        pa,p=cloud();before=particle_hash(pa);trial=ProcessThermalTrial(pa,p,0.)
        c=trial.candidate();evaluation=trial.evaluate(c,1e-9);b=[];end=[]
        for row in evaluation['boundary']['traces']:
            r=response_from_trace(row,dt=1e-9,mu=p['friction'])
            self.assertTrue(r['available']);b.append(r['B']);end.append(r['C'])
        pre=StageEventPreconditioner(pa,p,1e-9,0.,np.asarray(b),np.asarray(end),{'fixture':'ZERO_NOT_AFSD'})
        result=solve_stage(trial,pre,1e-9)
        checks=acceptance(pa,p,trial,result['candidate'],result['evaluation'])
        self.assertEqual(checks['status'],'ACCEPTED_NUMERICAL_STAGE_APPROXIMATE_NOT_ENERGY_CERTIFIED')
        self.assertEqual(particle_hash(pa),before)
        self.assertEqual(len(result['trace']),1)
        for key in c:np.testing.assert_array_equal(result['candidate'][key],c[key])

    def test_driver_will_not_overwrite_evidence(self):
        with tempfile.TemporaryDirectory(prefix='afsd_event_existing_') as folder:
            root=Path(folder)
            with mock.patch.object(runner,'ROOT',root):
                with self.assertRaises(ValueError):runner.run(root,root/'none.json','response')


if __name__=='__main__':unittest.main(verbosity=2)

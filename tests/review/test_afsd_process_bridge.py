"""Independent bridge interface tests, fake evaluator/clock only.

No checkpoint, PySPH evaluator, solver, contact/material or trajectory call.
Preregistration: docs/PROCESS_BRIDGE_TEST_CONTRACT.md.
"""
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import afsd_process_bridge as tested

FIELDS=('position_m','velocity_m_s','density_kg_m3','temperature_K','mean_velocity_m_s')
LABELS=('base','plus','minus','base_repeat')
SIGNS=(0,1,-1,0)
DT=5e-8


def geometry():
    return (np.array([[0.,0.,0.],[2.,0.,0.],[0.,2.,0.],[0.,0.,2.]]),
            np.array([1.,2.,3.,4.]),np.array([40,10,30,20]),
            np.array([[1.,1.,1.],[1.,1.,0.],[1.,1.,0.],[1.,1.,1.]]))


def candidate():
    return dict(position_m=np.array([[.002,-.001,.003],[-.004,.003,.002]]),
                velocity_m_s=np.array([[.1,-.2,.3],[.4,.5,-.6]]),
                density_kg_m3=np.array([2701.,2699.]),
                temperature_K=np.array([400.,405.]),
                mean_velocity_m_s=np.array([[.04,-.07,.1],[.2,.15,-.21]]))


def trial_direction():
    return np.array([[.6,.8,0.],[-.2,.3,-.4]])


def literal_perturb(base,direction,sign):
    result={key:np.asarray(base[key],float).copy() for key in FIELDS}
    increment=sign*1e-7*direction
    result['mean_velocity_m_s']=base['mean_velocity_m_s']+increment
    result['velocity_m_s']=base['velocity_m_s']+2*increment
    result['position_m']=base['position_m']+DT*increment
    return result


class FakeClock:
    def __init__(self,start=0.):self.value=float(start)
    def __call__(self):return self.value
    def advance(self,seconds):self.value+=seconds


class BridgeInterfaces(unittest.TestCase):
    def fields_equal(self,left,right):
        self.assertEqual(set(left),set(right))
        for key in FIELDS:np.testing.assert_array_equal(left[key],right[key])

    def test_01_literal_nonuniform_mass_direction_and_constraints(self):
        args=geometry();before=copy.deepcopy(args)
        result=tested.direction(*args)
        raw=np.array([[-.325,.075,.15625],[.675,.325,0.],
                      [-.075,-.425,0.],[-.2,.1375,-.34375]])
        expected=raw/np.sqrt(.675**2+.325**2)
        np.testing.assert_allclose(result,expected,atol=2e-14,rtol=0.)
        self.assertAlmostEqual(float(np.linalg.norm(result,axis=1).max()),1.,places=14)
        np.testing.assert_array_equal(result[args[3]==0],0.)
        for actual,original in zip(args,before):np.testing.assert_array_equal(actual,original)
        self.assertFalse(np.shares_memory(result,args[0]))

    def test_02_id_reordering_and_geometry_translation_scale(self):
        args=geometry();original=tested.direction(*args)
        for permutation in (np.array([3,1,0,2]),np.array([2,0,3,1])):
            reordered=tested.direction(*(value[permutation] for value in args))
            np.testing.assert_array_equal(reordered,original[permutation])
        x,m,ids,free=args
        np.testing.assert_allclose(tested.direction(x+[.5,-.25,1.],m,ids,free),
                                   original,atol=2e-14,rtol=0.)
        np.testing.assert_allclose(tested.direction(8*x,4*m,ids,free),
                                   original,atol=2e-14,rtol=0.)

    def test_03_direction_rejects_invalid_identity_geometry_or_mobility(self):
        x,m,ids,free=geometry()
        changes=[(0,np.zeros((4,2))),(0,np.full((4,3),np.nan)),(0,np.ones((4,3))),
                 (1,np.array(1.)),(1,np.array([1.,2.,3.])),(1,np.array([1.,0.,3.,4.])),
                 (1,np.array([1.,-2.,3.,4.])),(1,np.array([1.,np.inf,3.,4.])),
                 (2,np.array([1,1,3,4])),(2,np.array([1.,2.5,3.,4.])),
                 (2,np.array([1.,np.nan,3.,4.])),(3,np.full((4,3),.5)),
                 (3,np.zeros((4,3))),(3,np.ones((3,3))),(3,np.full((4,3),np.nan))]
        for index,value in changes:
            args=[x,m,ids,free];args[index]=value
            with self.subTest(index=index,value=str(value)):
                with self.assertRaises(ValueError):tested.direction(*args)
        with self.assertRaises(ValueError):
            tested.direction(np.empty((0,3)),np.empty(0),np.empty(0,int),np.empty((0,3)))

    def test_04_literal_five_field_perturbations_detached_same_base(self):
        base=candidate();before=copy.deepcopy(base);direction=trial_direction();saved=[]
        for sign in (0,1,-1):
            result=tested.perturb(base,direction,DT,sign);saved.append(result)
            self.fields_equal(result,literal_perturb(before,direction,sign))
            self.fields_equal(base,before)
            for key in FIELDS:self.assertFalse(np.shares_memory(result[key],base[key]))
        for first in range(3):
            for second in range(first):
                for key in FIELDS:self.assertFalse(np.shares_memory(saved[first][key],saved[second][key]))
        saved[1]['density_kg_m3'][0]=999.
        self.fields_equal(base,before)
        np.testing.assert_array_equal(saved[2]['density_kg_m3'],before['density_kg_m3'])

    def test_05_perturb_rejects_invalid_five_field_or_scalar_domain(self):
        base=candidate();direction=trial_direction();before=copy.deepcopy(base)
        for dt in (0.,-1.,np.inf,np.nan,np.array([DT])):
            with self.subTest(dt=str(dt)):
                with self.assertRaises(ValueError):tested.perturb(base,direction,dt,1)
        for sign in (True,np.bool_(True),2,.5,np.nan,np.array([1])):
            with self.subTest(sign=str(sign)):
                with self.assertRaises(ValueError):tested.perturb(base,direction,DT,sign)
        for wrong in (np.ones((2,2)),np.full((2,3),np.nan)):
            with self.assertRaises(ValueError):tested.perturb(base,wrong,DT,1)
        bads=[]
        missing=copy.deepcopy(base);missing.pop('mean_velocity_m_s');bads.append(missing)
        extra=copy.deepcopy(base);extra['passive_label']=np.zeros(2);bads.append(extra)
        for key,value in (('density_kg_m3',np.array(2700.)),('temperature_K',np.zeros((2,1))),
                          ('position_m',np.ones((3,3))),('velocity_m_s',np.full((2,3),np.nan))):
            wrong=copy.deepcopy(base);wrong[key]=value;bads.append(wrong)
        for wrong in bads:
            with self.assertRaises(ValueError):tested.perturb(wrong,direction,DT,1)
        self.fields_equal(base,before)

    def test_06_budget_old_floor_maximum_and_inclusive_boundary(self):
        reference_cost=4.04492670000036
        cases=[(0.,4,[]),(20.,3,[1.,2.]),(40.,1,[1.]),(9.,3,[8.]),
               (9.+1e-10,3,[8.]),(12.,3,[8.,1.,1.]),(45.,0,[]),(45.+1e-10,0,[])]
        for elapsed,remaining,measured in cases:
            expected=elapsed+1.5*remaining*max([reference_cost]+measured)+5<=50
            with self.subTest(elapsed=elapsed,remaining=remaining,measured=measured):
                self.assertEqual(bool(tested.may_continue(elapsed,remaining,measured)),expected)
        self.assertTrue(tested.may_continue(9.,3,[8.]))
        self.assertFalse(tested.may_continue(12.,3,[8.,1.,1.]))

    def test_07_budget_invalid_values_reject(self):
        for elapsed,remaining,measured in ((-1.,1,[]),(np.inf,1,[]),(np.nan,1,[]),
            (0.,-1,[]),(0.,5,[]),(0.,1.,[]),(0.,True,[]),
            (0.,1,[-1.]),(0.,1,[np.inf]),(0.,1,[np.nan])):
            with self.subTest(elapsed=elapsed,remaining=remaining,measured=measured):
                with self.assertRaises(ValueError):tested.may_continue(elapsed,remaining,measured)

    def test_08_four_calls_same_base_archive_before_next_and_no_map_iteration(self):
        base=candidate();before=copy.deepcopy(base);direction=trial_direction();clock=FakeClock()
        calls=[];events=[];evaluations=[]
        def evaluator(c,dt):
            index=len(calls)+1;self.assertEqual(dt,DT)
            self.assertEqual(len(events),2*(index-1))
            self.fields_equal(c,literal_perturb(before,direction,SIGNS[index-1]))
            calls.append(copy.deepcopy(c));events.append(('evaluate',index));clock.advance(1.)
            result=dict(map={key:value+123. for key,value in c.items()},
                        accepted=False,restart_allowed=False,can_publish=False)
            evaluations.append(copy.deepcopy(result));return result
        def on_result(index,label,c,evaluation,seconds):
            self.assertEqual(label,LABELS[index-1]);self.assertEqual(seconds,1.)
            self.fields_equal(c,calls[index-1]);self.assertFalse(evaluation['accepted'])
            events.append(('archive',index));clock.advance(.25)
            return dict(index=index,label=label,detached=True)
        result=tested.run_four(evaluator,base,direction,DT,clock,on_result)
        self.assertEqual(result['status'],'FOUR_EVALUATIONS_COMPLETED_NOT_ACCEPTANCE')
        self.assertEqual(result['completed'],4);self.assertEqual(result['not_run'],[])
        self.assertEqual(len(calls),4);self.assertEqual(len(result['records']),4)
        self.assertEqual(events,[(kind,index) for index in range(1,5) for kind in ('evaluate','archive')])
        self.fields_equal(calls[0],calls[3]);self.fields_equal(evaluations[0]['map'],evaluations[3]['map'])
        self.fields_equal(base,before);self.assertEqual(result['elapsed_seconds'],5.)

    def test_09_budget_partial_includes_loading_measured_cost_and_archive_time(self):
        for loading,duration,archive,expected in ((30.,1.,0.,0),(0.,18.,0.,1),(0.,1.,30.,1)):
            with self.subTest(loading=loading,duration=duration,archive=archive):
                clock=FakeClock(loading);calls=[];saved=[]
                def evaluator(c,dt):calls.append(copy.deepcopy(c));clock.advance(duration);return {}
                def on_result(index,label,c,evaluation,seconds):
                    saved.append(index);clock.advance(archive);return dict(index=index)
                result=tested.run_four(evaluator,candidate(),trial_direction(),DT,clock,on_result,started=0.)
                self.assertEqual(result['status'],'PARTIAL_BUDGET_NOT_RUN_REMAINDER')
                self.assertEqual(result['completed'],expected);self.assertEqual(len(calls),expected)
                self.assertEqual(saved,list(range(1,expected+1)))
                self.assertEqual(result['not_run'],list(LABELS[expected:]))
                self.assertEqual(len(result['records']),expected)

    def test_10_evaluator_exception_keeps_prior_receipts_no_retry(self):
        base=candidate();before=copy.deepcopy(base);clock=FakeClock();calls=[];saved=[]
        failure=ValueError('synthetic rejected domain')
        def evaluator(c,dt):
            calls.append(copy.deepcopy(c));clock.advance(1.)
            if len(calls)==2:raise failure
            return {}
        def on_result(index,label,c,evaluation,seconds):saved.append(index);return dict(index=index)
        with self.assertRaises(tested.BridgeFailure) as caught:
            tested.run_four(evaluator,base,trial_direction(),DT,clock,on_result)
        error=caught.exception
        self.assertEqual(error.index,2);self.assertEqual(error.label,'plus')
        self.assertIs(error.original,failure);self.assertEqual(len(calls),2);self.assertEqual(saved,[1])
        self.assertEqual(len(error.records),1);self.fields_equal(error.candidate,calls[1])
        self.fields_equal(base,before)

    def test_11_callback_exception_and_input_mutation_are_terminal(self):
        for mode in ('callback_exception','evaluator_mutation','callback_mutation'):
            with self.subTest(mode=mode):
                base=candidate();before=copy.deepcopy(base);clock=FakeClock();counts=[0,0]
                def evaluator(c,dt):
                    counts[0]+=1;clock.advance(1.)
                    if mode=='evaluator_mutation':c['density_kg_m3'][0]+=1.
                    return {}
                def on_result(index,label,c,evaluation,seconds):
                    counts[1]+=1
                    if mode=='callback_exception':raise OSError('synthetic archive failure')
                    c['temperature_K'][0]+=1.;return {}
                with self.assertRaises(tested.BridgeFailure) as caught:
                    tested.run_four(evaluator,base,trial_direction(),DT,clock,on_result)
                self.assertEqual(caught.exception.index,1);self.assertEqual(caught.exception.records,[])
                self.assertEqual(counts,[1,0 if mode=='evaluator_mutation' else 1])
                self.fields_equal(base,before)

    def test_12_nonmonotonic_measurement_rejects_without_archive(self):
        clock=FakeClock();counts=[0,0]
        def evaluator(c,dt):counts[0]+=1;clock.advance(-1.);return {}
        def on_result(*args):counts[1]+=1;return {}
        with self.assertRaises(tested.BridgeFailure):
            tested.run_four(evaluator,candidate(),trial_direction(),DT,clock,on_result)
        self.assertEqual(counts,[1,0])

    def test_13_forecast_uses_full_evaluation_archive_and_audit_cost(self):
        base=candidate();before=copy.deepcopy(base);clock=FakeClock();calls=[];saved=[]
        def evaluator(c,dt):
            calls.append(copy.deepcopy(c));clock.advance(1.)
            return dict(accepted=False,restart_allowed=False,can_publish=False)
        def on_result(index,label,c,evaluation,seconds):
            self.assertEqual(index,1);self.assertEqual(label,'base')
            self.assertEqual(seconds,1.);self.fields_equal(c,before)
            # Nine seconds represent saving AND independently auditing this
            # result, which future calls will also need. They cannot count
            # only as sunk elapsed time while disappearing from the forecast.
            clock.advance(9.);saved.append(index);return dict(saved_and_audited=True)
        result=tested.run_four(evaluator,base,trial_direction(),DT,clock,on_result)
        self.assertEqual(result['status'],'PARTIAL_BUDGET_NOT_RUN_REMAINDER')
        self.assertEqual(result['completed'],1);self.assertEqual(len(calls),1)
        self.assertEqual(saved,[1]);self.assertEqual(len(result['records']),1)
        self.assertEqual(result['not_run'],['plus','minus','base_repeat'])
        self.assertEqual(result['records'][0]['evaluation_seconds'],1.)
        self.assertEqual(result['records'][0]['complete_measurement_seconds'],10.)
        self.assertEqual(result['next_forecast_seconds'],60.)
        self.assertEqual(clock(),10.);self.fields_equal(base,before)

    def test_14_explicit_budget_110_preserves_default50_and_rejects_invalid_limits(self):
        self.assertFalse(tested.may_continue(10.,3,[10.]))
        self.assertFalse(tested.may_continue(10.,3,[10.],limit_seconds=50.))
        self.assertTrue(tested.may_continue(10.,3,[10.],limit_seconds=110.))
        for limit,expected in ((50.,1),(110.,4)):
            with self.subTest(limit_seconds=limit):
                base=candidate();before=copy.deepcopy(base);clock=FakeClock();calls=[];saved=[]
                def evaluator(c,dt):
                    index=len(calls);self.fields_equal(c,literal_perturb(before,trial_direction(),SIGNS[index]))
                    calls.append(copy.deepcopy(c));clock.advance(1.);return {}
                def on_result(index,label,c,evaluation,seconds):
                    self.assertEqual(seconds,1.);self.assertEqual(label,LABELS[index-1])
                    saved.append(index);clock.advance(9.);return dict(index=index)
                result=tested.run_four(evaluator,base,trial_direction(),DT,clock,on_result,
                                       limit_seconds=limit)
                self.assertEqual(result['completed'],expected);self.assertEqual(len(calls),expected)
                self.assertEqual(saved,list(range(1,expected+1)))
                self.assertEqual(result['not_run'],list(LABELS[expected:]))
                self.assertEqual(clock(),10.*expected);self.fields_equal(base,before)
                for row in result['records']:
                    self.assertEqual(row['evaluation_seconds'],1.)
                    self.assertEqual(row['complete_measurement_seconds'],10.)
                self.assertEqual(result['status'],
                    'FOUR_EVALUATIONS_COMPLETED_NOT_ACCEPTANCE' if expected==4
                    else 'PARTIAL_BUDGET_NOT_RUN_REMAINDER')
        for invalid in (601.,0.,np.nan,True,np.bool_(True)):
            with self.subTest(invalid_limit=str(invalid)):
                with self.assertRaises(ValueError):
                    tested.may_continue(0.,4,[],limit_seconds=invalid)
                calls=[]
                def evaluator(c,dt):calls.append(c);return {}
                with self.assertRaises(ValueError):
                    tested.run_four(evaluator,candidate(),trial_direction(),DT,FakeClock(),
                                    lambda *args:{},limit_seconds=invalid)
                self.assertEqual(calls,[])


if __name__=='__main__':
    paths=['scripts/afsd_process_bridge.py','tests/review/test_afsd_process_bridge.py',
           'docs/PROCESS_BRIDGE_TEST_CONTRACT.md']
    print('BRIDGE_FAKE_TEST_SOURCE_RECEIPT '+json.dumps(dict(
        classification='FAKE_INTERFACE_TESTS_NO_PHYSICAL_EVALUATION',
        source_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}),
        sort_keys=True),flush=True)
    unittest.main(verbosity=2)

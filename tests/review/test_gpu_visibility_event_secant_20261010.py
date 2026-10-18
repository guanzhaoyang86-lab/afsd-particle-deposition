"""Verify unchanged ordinary physics and unaccepted same-parent probe guesses."""
from pathlib import Path
import ast
import copy
import inspect
import sys
import types
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import run_gpu_visibility_event_20261010 as previous
import run_gpu_visibility_event_secant_20261010 as current


class FrozenOrdinaryPath(unittest.TestCase):
    def test_global_call_budget_covers_two_independent_bounded_event_locations(self):
        self.assertEqual(current.MAXIMUM_PHYSICAL_SOLVES,3+2*10+2)
        self.assertEqual(current.MAXIMUM_EVENTS_PER_BRANCH,previous.MAXIMUM_EVENTS_PER_BRANCH)
        self.assertEqual(current.EVENT_TIME_TOLERANCE_S,previous.EVENT_TIME_TOLERANCE_S)
        self.assertEqual(current.EVENT_POSITION_TOLERANCE_M,previous.EVENT_POSITION_TOLERANCE_M)

    def test_trial_solver_and_physical_providers_are_same_objects(self):
        for name in ('model','solver','gpu','cpu','core','post','event','batch','old','shared'):
            self.assertIs(getattr(current,name),getattr(previous,name))

    def test_all_original_gate_and_geometry_functions_remain_exact_ast(self):
        for name in ('resolved_routing','scalar_equivalence','routing','qualify','leaf_row',
                     'copy_macro_endpoint','classify_activation'):
            self.assertEqual(ast.dump(ast.parse(inspect.getsource(getattr(current,name)))),
                             ast.dump(ast.parse(inspect.getsource(getattr(previous,name)))))

    def test_ordinary_physical_solve_only_adds_optional_probe_guess(self):
        actual=ast.parse(inspect.getsource(current.solve_physical))
        function=actual.body[0]
        self.assertEqual(function.args.args[-1].arg,'initial_guess')
        function.args.args.pop();function.args.defaults.pop()
        extra=[n for n in function.body if isinstance(n,ast.If)
               and isinstance(n.test,ast.Compare) and isinstance(n.test.left,ast.Name)
               and n.test.left.id=='initial_guess']
        self.assertEqual(len(extra),1);function.body.remove(extra[0])
        self.assertEqual(ast.dump(actual),ast.dump(ast.parse(inspect.getsource(previous.solve_physical))))

    def test_original_pair_gate_and_history_code_only_add_locator_manifest_tag(self):
        actual=ast.parse(inspect.getsource(current.run_pair))
        found=0
        for node in ast.walk(actual):
            if isinstance(node,ast.Call):
                extra=[k for k in node.keywords if k.arg=='event_locator_protocol']
                found+=len(extra)
                node.keywords[:]=[k for k in node.keywords if k.arg!='event_locator_protocol']
        self.assertEqual(found,1)
        self.assertEqual(ast.dump(actual),ast.dump(ast.parse(inspect.getsource(previous.run_pair))))


class ProbeInitialGuess(unittest.TestCase):
    def setUp(self):
        self.free=np.array([[True,False,True],[False,True,False]])
        self.state=dict(free=self.free,identity='original-parent')
        self.fixed=np.array([[7.,8.,9.],[4.,5.,6.]])*1e-5
        self.trial=types.SimpleNamespace(pa=object(),state=self.state,
            incoming_pid_pairs=((11,13),),x0=np.zeros((2,3)),p=dict(dx=.01),
            prescribed=lambda dt:self.fixed*dt)
        self.before=np.array([[1.,2.,3.],[4.,5.,6.]])*1e-7
        prior_trial=types.SimpleNamespace(pa=object(),state=copy.deepcopy(self.state),
            incoming_pid_pairs=((11,13),),motion=types.SimpleNamespace(dt_s=.01))
        self.prior=dict(trial=prior_trial,fresh=dict(displacement=self.before.copy()))
        self.fingerprint=patch.object(current.cpu,'parent_fingerprint',side_effect=lambda pa,s:s['identity'])
        self.fingerprint.start();self.addCleanup(self.fingerprint.stop)

    def test_only_free_guess_scales_prescribed_axes_are_recomputed(self):
        result=current.event_probe_guess(self.trial,.02,self.prior)
        np.testing.assert_array_equal(result[self.free],2*self.before[self.free])
        np.testing.assert_array_equal(result[~self.free],(.02*self.fixed)[~self.free])
        np.testing.assert_array_equal(self.prior['fresh']['displacement'],self.before)
        self.assertFalse(np.shares_memory(result,self.prior['fresh']['displacement']))
        self.assertNotIn('time_s',self.state)

    def test_different_parent_is_rejected(self):
        self.prior['trial'].state['identity']='other-parent'
        with self.assertRaises(ValueError):current.event_probe_guess(self.trial,.02,self.prior)

    def test_different_or_unmasked_event_is_rejected(self):
        for pairs in ((),((11,15),)):
            self.prior['trial'].incoming_pid_pairs=pairs
            with self.assertRaises(ValueError):current.event_probe_guess(self.trial,.02,self.prior)

    def test_nonfinite_and_out_of_domain_displacements_are_rejected(self):
        for value in (np.nan,1.):
            self.prior['fresh']['displacement'][0,0]=value
            with self.assertRaises(ValueError):current.event_probe_guess(self.trial,.02,self.prior)

    def test_invalid_time_and_shape_are_rejected(self):
        for dt in (0.,-1.,np.nan,np.inf):
            with self.assertRaises(ValueError):current.event_probe_guess(self.trial,dt,self.prior)
        self.prior['fresh']['displacement']=np.zeros((1,3))
        with self.assertRaises(ValueError):current.event_probe_guess(self.trial,.02,self.prior)


if __name__=='__main__':
    unittest.main(verbosity=2)

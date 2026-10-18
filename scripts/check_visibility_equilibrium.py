"""Manufactured graph-switch search checks, not AFSD physical validation."""
from pathlib import Path
import copy
import json
import tempfile
import time
import unittest
from unittest.mock import patch
import numpy as np
import visibility_equilibrium as search


class PiecewiseForce:
    def __init__(self, oscillating=False, forged_zero=False):
        self.state={'free':np.array([[True,False,False]]),'time_s':0.,'history':np.array([3.])}
        self.oscillating=oscillating;self.forged_zero=forged_zero;self.calls=0;self.builds=[]
    def evaluate(self,u,dt):
        u=np.array(u,copy=True);self.calls+=1;x=float(u[0,0])
        f=x+2*max(x-.5,0)-1.4
        if self.oscillating:f=np.exp(x)-3
        graph=str(self.calls) if self.oscillating else ('left' if x<=.5 else 'right')
        force=np.array([[f,0.,0.]])
        return dict(displacement=u,force=force,force_limit_N=1e-11,
            scaled_residual=0. if self.forged_zero else abs(f)/1e-11,
            graph_fingerprint=graph,dt_s=dt)
    def tangent(self,r):
        x=float(r['displacement'][0,0]);self.builds.append(r['graph_fingerprint'])
        return np.diag([np.exp(x) if self.oscillating else 1. if x<=.5 else 3.,1.,1.])


class GraphSwitchTests(unittest.TestCase):
    def run_case(self,trial,maximum=24):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        directory=Path(self.temp.name)
        with patch.object(search,'journal'),patch.object(search,'summary',return_value={}):
            result=search.solve(trial,.1,np.array([[0.,.2,.3]]),directory,time.perf_counter()+20,maximum)
        return result,json.loads((directory/'iterations.json').read_text())
    def test_piecewise_force_rebuilds_instead_of_secant_across_graph_change(self):
        trial=PiecewiseForce();before=copy.deepcopy(trial.state)
        (result,metrics),rows=self.run_case(trial)
        self.assertAlmostEqual(result['displacement'][0,0],.8,places=11)
        self.assertEqual(trial.builds,['left','right'])
        self.assertEqual(metrics['jacobian_refreshes'],1)
        changed=[r for r in rows if r.get('secant',{}).get('status')=='DISCARDED_ACROSS_GRAPH_CHANGE']
        self.assertEqual(len(changed),1)
        np.testing.assert_array_equal(result['displacement'][0,1:],[.2,.3])
        np.testing.assert_array_equal(trial.state['history'],before['history'])
    def test_actual_force_overrides_forged_convergence(self):
        (result,metrics),_=self.run_case(PiecewiseForce(forged_zero=True))
        self.assertGreater(metrics['iterations'],0)
        self.assertLessEqual(abs(result['force'][0,0]),1e-11)
    def test_repeated_graph_change_exhausts_original_refresh_budget(self):
        trial=PiecewiseForce(oscillating=True)
        with self.assertRaisesRegex(ValueError,'refresh budget exhausted'):
            self.run_case(trial)
        self.assertEqual(len(trial.builds),3)
        report=json.loads((Path(self.temp.name)/'secant_failure.json').read_text())
        self.assertFalse(report['accepted']);self.assertFalse(report['restart_allowed'])
        self.assertEqual(report['jacobian_refreshes'],2)
        self.assertLessEqual(report['evaluations'],100)


if __name__=='__main__':unittest.main(verbosity=2)

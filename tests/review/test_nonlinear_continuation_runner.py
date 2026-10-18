"""Manufactured stationary workflow and budget/transaction regressions, not AFSD."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(Path(__file__).parent))
import run_nonlinear_continuation as r
from test_updated_trial import cloud


def fixture():
    pa,p=cloud();pa.add_property('as_rs_tag',type='long',data=np.zeros(len(pa.x),dtype=int))
    return pa,{'t':0.},p,{}, {},'MANUFACTURED_ZERO_EQUILIBRIUM',dict(path='MANUFACTURED_NOT_AFSD',source_time_s=0.,source_step=0)


class RunnerTests(unittest.TestCase):
    def test_01_budget_arithmetic_scope_and_hash(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as directory:
            d=Path(directory)
            names=['resolved_config.json','last_accepted_reference.json','probe_trial.json.gz']
            for name in names:(d/name).write_text('METADATA_FIXTURE_NOT_PHYSICS',encoding='utf-8')
            report=dict(status='COMPLETE_UNCOMMITTED_COST_PROBE',mode='COST_PROBE',classification=r.CLASSIFICATION,
                production_authorized=False,accepted_steps=0,settings=r.SETTINGS,source_sha256={},configuration_sha256='cfg',
                complete_trial_wall_seconds=3.,artifact_sha256={name:r.sha(d/name) for name in names})
            p=d/'manifest.json';p.write_text(json.dumps(report),encoding='utf-8')
            self.assertEqual(r.budget(p,{},'cfg')['forecast_wall_seconds'],354.)
            report['complete_trial_wall_seconds']=4.;p.write_text(json.dumps(report),encoding='utf-8')
            self.assertEqual(r.budget(p,{},'cfg')['evaluation_budget'],70)
            report['complete_trial_wall_seconds']=10.;p.write_text(json.dumps(report),encoding='utf-8')
            self.assertFalse(r.budget(p,{},'cfg')['allowed'])
            with self.assertRaises(ValueError):r.budget(p,{'a':'changed'},'cfg')
            (d/names[0]).write_text('changed',encoding='utf-8')
            with self.assertRaises(ValueError):r.budget(p,{},'cfg')

    def test_02_state_quantities_are_not_summed_as_work(self):
        ledger=dict(boundary={},thermal={},energy_diagnostic_residual_J=0.)
        e=dict(boundary=dict(ledger=dict(wall_work_J=2.,normal_potential_J=3.,max_penetration_m=.01,
            contacts=2,minimum_raw_friction_J=-1e-25,tool_impulse_Ns=[1.,2.,3.],contact_subintervals=4)),
            thermal=dict(ledger=dict(plastic_heat_J=5.)))
        check=dict(total_energy_diagnostic_residual_J=.5)
        r.accumulate(ledger,e,check);e['boundary']['ledger']['normal_potential_J']=7.
        r.accumulate(ledger,e,check)
        self.assertEqual(ledger['boundary']['wall_work_J'],4.)
        self.assertNotIn('normal_potential_J',ledger['boundary'])
        self.assertEqual(ledger['boundary_state_latest']['normal_potential_J'],7.)
        self.assertEqual(ledger['boundary_state_latest']['contacts'],2)
        self.assertEqual(ledger['thermal']['plastic_heat_J'],10.)
        self.assertEqual(ledger['boundary']['contact_subintervals'],8)
        self.assertEqual(ledger['minimum_raw_friction_J'],-1e-25)

    def test_03_unknown_ledger_field_is_not_silently_accumulated(self):
        e=dict(boundary=dict(ledger=dict(unknown_temperature_K=300.)),thermal=dict(ledger={}))
        with self.assertRaises(ValueError):r.accumulate(dict(boundary={},thermal={}),e,{})

    def test_04_existing_directory_is_preserved(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as d,patch.object(r,'setup') as setup:
            with self.assertRaises(ValueError):r.run(Path(d),True)
            setup.assert_not_called()

    def test_05_manufactured_four_stage_transaction_real_roundtrip(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp,patch.object(r,'setup',side_effect=fixture),patch.object(r,'budget',return_value={'allowed':True,'evaluation_budget':72}):
            d=Path(temp)/'manufactured_zero';self.assertEqual(r.run(d,False,Path('FIXTURE')),0)
            m=json.loads((d/'manifest.json').read_text());history=json.loads((d/'accepted_history.json').read_text())
            self.assertEqual(m['accepted_steps'],4);self.assertEqual(m['total_evaluations'],4)
            self.assertEqual(len(history),4);self.assertEqual(m['accepted_time_s'],2e-7)
            for idx in range(1,5):
                x=json.loads((d/f'stage_{idx:03d}'/'transport.json').read_text())
                self.assertEqual(x['crossing']['status'],'NO_CENTER_CROSSING_NOT_DEPOSITION')
                self.assertEqual(x['groups']['substrate']['max_lab_displacement_m'],0.)

    def test_06_failed_publication_does_not_count_orphan_checkpoint(self):
        writer=r.write_json
        def failing(path,value,*args,**kwargs):
            if Path(path).name=='accepted_history.json' and value:raise OSError('INJECTED_PUBLICATION_FAILURE')
            return writer(path,value,*args,**kwargs)
        with tempfile.TemporaryDirectory(dir=ROOT/'outputs') as temp,patch.object(r,'setup',side_effect=fixture),patch.object(r,'budget',return_value={'allowed':True,'evaluation_budget':72}),patch.object(r,'write_json',side_effect=failing):
            d=Path(temp)/'manufactured_failure';self.assertEqual(r.run(d,False,Path('FIXTURE')),2)
            m=json.loads((d/'manifest.json').read_text());self.assertEqual(m['accepted_steps'],0)
            self.assertEqual(json.loads((d/'accepted_history.json').read_text()),[])
            self.assertTrue((d/'stage_001'/'accepted.npz').exists())
            self.assertTrue((d/'stage_001'/'REJECTED_candidate.json.gz').exists())


if __name__=='__main__':unittest.main(verbosity=2)

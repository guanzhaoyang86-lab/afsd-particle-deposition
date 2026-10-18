"""Storage and budget regressions; manufactured zero states are not AFSD progress."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
from unittest import mock
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from test_updated_trial import cloud
from afsd_backend.pysph.process_thermal_trial import ProcessThermalTrial,acceptance,CLASSIFICATION as TRIAL_CLASSIFICATION
from afsd_backend.pysph import process_global_state as storage
from afsd_backend.pysph.continuation_state import restore as old_restore
from afsd_backend.pysph.updated_trial import particle_hash
import run_process_global as runner

SOURCES={'MANUFACTURED_STORAGE_ONLY':'a'*64};CONFIG='b'*64


def zero():
    pa,p=cloud();pa.add_property('as_rs_tag',type='long',data=np.arange(len(pa.x))%3-1)
    pa.set_output_arrays(list(pa.properties))
    trial=ProcessThermalTrial(pa,p,0.);c=trial.candidate();e=trial.evaluate(c,1e-9)
    checks=acceptance(pa,p,trial,c,e)
    return pa,p,c,e,checks


class ProcessGlobalRunnerTests(unittest.TestCase):
    def test_named_actual_zero_commit_roundtrip_vtk_and_old_loader_rejection(self):
        pa,p,c,e,checks=zero();before=particle_hash(pa);after=storage.commit(pa,p,c,e,checks)
        self.assertEqual(particle_hash(pa),before)
        s=storage.state(after,p,1e-9,1,SOURCES,CONFIG,{'fixture':'ZERO_NOT_AFSD'},{},{})
        with tempfile.TemporaryDirectory(prefix='afsd_global_snapshot_') as folder:
            out=Path(folder);storage.save(out,after,p,s);digest=runner.sha(out/'accepted.npz')
            loaded,re=storage.restore(out/'accepted.npz',p,SOURCES,CONFIG,digest)
            self.assertEqual(re,s);self.assertEqual(particle_hash(loaded),particle_hash(after))
            vtk=(out/'accepted.vtk').read_text()
            self.assertIn(storage.CLASSIFICATION,vtk);self.assertNotIn('SCALARS e ',vtk)
            self.assertIn('particle_id',vtk);self.assertIn('stress_xx_Pa',vtk)
            with self.assertRaises(ValueError):old_restore(out/'accepted.npz',p,SOURCES,CONFIG,digest)
            with self.assertRaises(ValueError):storage.restore(out/'accepted.npz',p,SOURCES,CONFIG,'0'*64)

    def test_commit_requires_new_acceptance_and_parent_identity(self):
        pa,p,c,e,checks=zero()
        for bad in ({},{'status':checks['status']},{'variant':TRIAL_CLASSIFICATION,'status':'REJECTED'}):
            with self.assertRaises(ValueError):storage.commit(pa,p,c,e,bad)
        bad=copy.deepcopy(e);bad['accepted_state_hash']='0'*64
        with self.assertRaises(ValueError):storage.commit(pa,p,c,bad,checks)
        bad=copy.deepcopy(e);bad['particle_id']=bad['particle_id'][::-1]
        with self.assertRaises(ValueError):storage.commit(pa,p,c,bad,checks)

    def test_budget_counts_actual_build_all32_evaluations_and_safety_margin(self):
        with tempfile.TemporaryDirectory(prefix='afsd_global_budget_') as folder:
            root=Path(folder);file=root/'manifest.json';preflight={'fixture':'PREFLIGHT_ONLY'}
            report=dict(status='COMPLETE_UNCOMMITTED_COST_PROBE',classification=runner.CLASSIFICATION,
                settings=runner.SETTINGS,accepted_steps=0,production_authorized=False,source_sha256=SOURCES,
                configuration_sha256=CONFIG,preflight=preflight,artifact_sha256={},
                input_sha256={'input':'hash'},parent_reference={'sha256':'parent'},
                preconditioner_build_wall_seconds=2.,trial_and_proposal_wall_seconds=4.)
            file.write_text(json.dumps(report),encoding='utf-8')
            with mock.patch.object(runner,'ROOT',root):r=runner.budget(file,SOURCES,CONFIG,preflight,{'input':'hash'},{'sha256':'parent'})
            self.assertEqual(r['forecast_wall_seconds'],225.);self.assertTrue(r['allowed'])
            report['trial_and_proposal_wall_seconds']=10.;file.write_text(json.dumps(report),encoding='utf-8')
            with mock.patch.object(runner,'ROOT',root):r=runner.budget(file,SOURCES,CONFIG,preflight,{'input':'hash'},{'sha256':'parent'})
            self.assertEqual(r['forecast_wall_seconds'],513.);self.assertFalse(r['allowed'])
            report['configuration_sha256']='stale';file.write_text(json.dumps(report),encoding='utf-8')
            with mock.patch.object(runner,'ROOT',root):
                with self.assertRaises(ValueError):runner.budget(file,SOURCES,CONFIG,preflight,{'input':'hash'},{'sha256':'parent'})

    def test_probe_cost_cannot_be_reused_with_different_parent(self):
        with tempfile.TemporaryDirectory(prefix='afsd_global_wrong_parent_') as folder:
            root=Path(folder);file=root/'manifest.json'
            report=dict(status='COMPLETE_UNCOMMITTED_COST_PROBE',classification=runner.CLASSIFICATION,
                settings=runner.SETTINGS,accepted_steps=0,production_authorized=False,source_sha256=SOURCES,
                configuration_sha256=CONFIG,preflight={},input_sha256={},parent_reference={'sha256':'old'})
            file.write_text(json.dumps(report),encoding='utf-8')
            with mock.patch.object(runner,'ROOT',root):
                with self.assertRaisesRegex(ValueError,'initial history'):
                    runner.budget(file,SOURCES,CONFIG,{}, {},{'sha256':'different'})

    def test_run_does_not_reuse_existing_directory(self):
        with tempfile.TemporaryDirectory(prefix='afsd_global_existing_') as folder:
            root=Path(folder)
            with mock.patch.object(runner,'ROOT',root):
                with self.assertRaises(ValueError):runner.run(root,root/'missing.json',probe=True)


if __name__=='__main__':unittest.main(verbosity=2)

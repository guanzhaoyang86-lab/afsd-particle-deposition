"""Local event screening and provisional-state fixtures; no AFSD solve."""
from pathlib import Path
import copy
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import run_first_shoulder_event as event


def particles():
    return SimpleNamespace(particle_id=np.array([10,20,30,40]),x=np.zeros(4),y=np.zeros(4),z=np.zeros(4),
        temperature=np.ones(4)*723.,S=np.zeros(36),ep=np.zeros(4))


class FirstShoulderControllerTests(unittest.TestCase):
    def test_one_hop_unions_real_graphs_without_recursive_expansion(self):
        ids=np.array([10,20,30,40,50])
        graphs=[(ids,np.array([0,1,2]),np.array([1,2,3])),(ids,np.array([1]),np.array([4]))]
        self.assertEqual(event.one_hop_ids({20},graphs),[10,20,30,50])

    def test_zero_difference_passes_and_reports_all_selected_ids(self):
        result=event.local_compare(particles(),particles(),[10,20])
        self.assertTrue(result['passed'])
        self.assertEqual([r['particle_id'] for r in result['rows']],[10,20])

    def test_one_particle_stress_failure_cannot_be_diluted_by_others(self):
        fine,coarse=particles(),particles();coarse.S.reshape(-1,9)[1,0]=100001.
        self.assertFalse(event.local_compare(fine,coarse,[10,20,30,40])['passed'])
        self.assertTrue(event.local_compare(fine,coarse,[10,30,40])['passed'])

    def test_plastic_strain_floor_applies_particlewise(self):
        fine,coarse=particles(),particles();coarse.ep[1]=2.0001e-6
        self.assertFalse(event.local_compare(fine,coarse,[20])['passed'])

    def test_position_limit_uses_vector_norm(self):
        fine,coarse=particles(),particles();coarse.x[1]=4e-6;coarse.y[1]=4e-6
        self.assertFalse(event.local_compare(fine,coarse,[20])['passed'])

    def test_temperature_failure_is_local(self):
        fine,coarse=particles(),particles();coarse.temperature[1]+=.5001
        self.assertFalse(event.local_compare(fine,coarse,[20])['passed'])

    def test_missing_or_reordered_particle_ids_are_rejected(self):
        fine,coarse=particles(),particles()
        with self.assertRaisesRegex(ValueError,'IDs incomplete'):
            event.local_compare(fine,coarse,[999])
        coarse.particle_id=coarse.particle_id[::-1]
        with self.assertRaisesRegex(ValueError,'ID order'):
            event.local_compare(fine,coarse,[20])

    def test_uncommitted_save_does_not_poison_live_half_step_state(self):
        # Half2 uses the detached numerical Half1 state in memory. Writing a
        # nonrestartable candidate must not set restart_allowed=False on it or
        # replace its baseline model source/provenance for the next trial.
        state=dict(source_sha256={'model':'original'},implementation_source_sha256={'parent':'frozen'},
            interval_energy={'records':[{'raw':1.}]},bore_observation={'eligible_ids':[10]},
            energy={'raw':1.},time_s=.0200625,step=137)
        before=copy.deepcopy(state);saved={};pa=particles();frozen={'controller':'new'}
        result=dict(i=np.array([0]),j=np.array([1]),thermal={'ledger':{'balance':0.}})
        def checkpoint(path,array,data):saved.update(copy.deepcopy(data))
        with tempfile.TemporaryDirectory() as tmp, patch.object(event.core,'journal'), \
             patch.object(event.core,'checkpoint',side_effect=checkpoint), \
             patch.object(event.core,'verify_snapshot'), \
             patch.object(event.core,'load',side_effect=lambda p:{'solver_data':saved}):
            target,evidence=event.save_candidate(Path(tmp),None,pa,state,result,frozen)
            self.assertEqual(target.name,'NOT_RESTARTABLE_candidate.npz')
            self.assertTrue(evidence.exists())
        self.assertEqual(state,before)
        self.assertNotIn('restart_allowed',state)
        self.assertFalse(saved['restart_allowed'])
        self.assertEqual(saved['source_sha256'],before['source_sha256'])
        self.assertEqual(saved['implementation_source_sha256'],frozen)

    def test_scalar_metadata_and_solid_certificate_are_compared(self):
        for changed_key,changed_value in (('normal_model','OTHER'),('sweep',{'certified':False})):
            contact=dict(normal_model='CHAIN_RULE',sweep={'certified':True})
            old=dict(contact=copy.deepcopy(contact));new=dict(displacement=np.zeros((1,3)),dt_s=.001,
                contact=copy.deepcopy(contact));new['contact'][changed_key]=changed_value
            trial=SimpleNamespace(pa=None,state={},p={})
            with tempfile.TemporaryDirectory() as tmp, \
                 patch.object(event.core,'DiscreteContactTrial',return_value=SimpleNamespace(evaluate=lambda *a:old)), \
                 patch.object(event.core,'energy_diagnostics',return_value={}), \
                 patch.object(event.batch,'candidate_equivalence',return_value={}):
                with self.assertRaisesRegex(ValueError,'contact metadata differs'):
                    event.scalar_equivalence(trial,new,Path(tmp))


if __name__=='__main__':unittest.main(verbosity=2)

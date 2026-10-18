"""Strict administrative recovery checks; no GPU solve or live startup."""
import ast
import copy
import inspect
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import run_gpu_nonlinear_resume_20261011 as subject


class RecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.meta=subject.read(subject.PARENT/'manifest.json')
        cls.pa,cls.state=subject.load(subject.ROOT/cls.meta['screened_continuation_checkpoint'])
        base=subject.FAILED/'coarse/proposal/broyden'
        cls.report=subject.read(base/'last_UNCOMMITTED_candidate.json')
        with np.load(base/'last_UNCOMMITTED_candidate.npz',allow_pickle=False) as npz:cls.arrays={k:npz[k].copy() for k in npz.files}

    def pool(self):
        return SimpleNamespace(sha=subject.bridge.core.sha,seal=lambda mapping:None,verify=lambda mapping:None)

    def validate_arrays(self,arrays=None,report=None,state=None):
        return subject.validate_candidate_arrays(arrays if arrays is not None else self.arrays,report if report is not None else self.report,
            self.pa,state if state is not None else self.state,.000474609375)

    def test_worker_supervisor_and_loop_have_only_registered_change(self):
        self.assertFalse(subject.static_audit()['physical_source_changed'])
        self.assertIs(subject.cost_decision,subject.previous.cost_decision)
        self.assertIs(subject.measured_envelope,subject.previous.measured_envelope)
        self.assertIs(subject.authorized_wall_seconds,subject.previous.authorized_wall_seconds)

    def test_authorization_not_renewed_on_resume(self):
        self.assertEqual(subject.AUTHORIZED_AT.isoformat(),'2026-10-11T15:41:47+08:00')
        self.assertEqual(subject.ABSOLUTE_STOP.isoformat(),'2026-10-11T23:41:47+08:00')
        self.assertIn("**admin.snapshot(read(PARENT/'manifest.json'),PARENT)",inspect.getsource(subject.run))

    def test_parent_sha_time_step_exact(self):
        subject.verify_parent(self.meta)
        for name,value in (('screened_continuation_sha256','f'*64),('screened_end_time_s',.089589599609),('screened_end_step',598)):
            row=copy.deepcopy(self.meta);row[name]=value
            with self.assertRaises(ValueError):subject.verify_parent(row)

    def test_real_audit_explains_exact_additional_failed_cost_file(self):
        evidence=subject.verify_old_audit()
        self.assertEqual(evidence['final_locked_files'],30284);self.assertEqual(evidence['final_expected_files'],30283)
        self.assertEqual(evidence['expected_final_hash_reads'],60568)
        self.assertEqual(evidence['expected_final_hash_bytes'],140877341052)
        self.assertTrue(evidence['final_audit_covers_all_held_files'])

    def test_other_audit_count_or_bytes_gap_cannot_be_accepted(self):
        original=subject.read;path=subject.OLD_CAMPAIGN/'final_history_audit.json'
        for key in ('expected_files','locked_files','hash_reads','hash_bytes'):
            row=original(path);row[key]+=1
            with self.subTest(key=key),patch.object(subject,'read',side_effect=lambda p:row if Path(p)==path else original(p)):
                with self.assertRaises(ValueError):subject.verify_old_audit()

    def test_failed_cost_ledger_cannot_claim_complete_pair(self):
        original=subject.read;path=subject.OLD_CAMPAIGN/'attempt_0040_failed_cost.json';row=original(path);row['complete_pair']=True
        with patch.object(subject,'read',side_effect=lambda p:row if Path(p)==path else original(p)):
            with self.assertRaises(ValueError):subject.verify_old_audit()

    def test_force_threshold_uses_actual_parent_free_mask(self):
        result=self.validate_arrays()
        self.assertGreater(result['scaled_residual'],60.)
        self.assertEqual(result['force_limit_N'],self.report['metrics']['force_limit_N'])

    def test_nonfinite_candidate_is_rejected(self):
        arrays=copy.deepcopy(self.arrays);arrays['force_N'][0,0]=float('nan')
        with self.assertRaises(ValueError):self.validate_arrays(arrays=arrays)

    def test_force_limit_cannot_be_relaxed(self):
        report=copy.deepcopy(self.report);report['metrics']['force_limit_N']*=100
        with self.assertRaises(ValueError):self.validate_arrays(report=report)

    def test_recorded_force_residual_cannot_replace_actual_array_value(self):
        report=copy.deepcopy(self.report);report['metrics']['force_residual_N']/=2
        with self.assertRaises(ValueError):self.validate_arrays(report=report)

    def test_free_mask_shape_and_type_must_be_real_parent_domain(self):
        for value in (np.ones((1381,3)),np.zeros((1381,3),dtype=bool),np.ones(1381,dtype=bool)):
            state=copy.deepcopy(self.state);state['free']=value
            with self.assertRaises(ValueError):self.validate_arrays(state=state)

    def test_candidate_cannot_change_particle_inventory(self):
        arrays=copy.deepcopy(self.arrays);arrays['particle_id'][0]=-999
        with self.assertRaises(ValueError):self.validate_arrays(arrays=arrays)

    def test_failed_candidate_cannot_be_marked_accepted(self):
        arrays=copy.deepcopy(self.arrays);arrays['accepted']=np.array(True)
        with self.assertRaises(ValueError):self.validate_arrays(arrays=arrays)

    def test_candidate_position_and_velocity_bound_to_parent_and_dt(self):
        for name in ('position_m','velocity_m_s'):
            arrays=copy.deepcopy(self.arrays);arrays[name][0,0]+=1e-6
            with self.assertRaises(ValueError):self.validate_arrays(arrays=arrays)

    def test_force_components_must_match_saved_contact_report(self):
        report=copy.deepcopy(self.report);report['contact']['force_N'][0][0]+=1
        with self.assertRaises(ValueError):self.validate_arrays(report=report)

    def test_negative_material_domain_rejected(self):
        for name in ('density_kg_m3','temperature_K','plastic_strain'):
            arrays=copy.deepcopy(self.arrays);arrays[name][0]=-1
            with self.assertRaises(ValueError):self.validate_arrays(arrays=arrays)

    def actual_pending(self):
        state=copy.deepcopy(subject.read(subject.OLD_CAMPAIGN/'final_summary.json')['adaptive_state'])
        state.update(status='READY',pending=None,attempts=state['attempts']-1)
        return subject.adaptive.begin_attempt(state,.000474609375)

    def test_actual_full_nonlinear_evidence_refines_same_parent(self):
        failure=subject.read(subject.FAILED/'manifest.json')
        policy,event=subject.nonlinear_recovery(self.pool(),subject.FAILED,failure,self.meta,self.actual_pending(),self.meta['source_sha256'])
        self.assertEqual(policy['nominal_dt_s'],.0002373046875);self.assertEqual(policy['parent_sha256'],subject.PARENT_SHA256)
        self.assertEqual(policy['passes'],56);self.assertEqual(policy['attempts'],62)
        self.assertEqual(policy['status'],'READY');self.assertIsNone(policy['pending']);self.assertEqual(event['action'],'REFINE')

    def test_unrelated_nonlinear_error_is_not_recoverable(self):
        failure=subject.read(subject.FAILED/'manifest.json');failure['error']="ValueError('thermal/material fixed point did not converge')"
        self.assertIsNone(subject.nonlinear_recovery(self.pool(),subject.FAILED,failure,self.meta,self.actual_pending(),self.meta['source_sha256']))

    def test_missing_registered_newton_evidence_does_not_qualify(self):
        failure=subject.read(subject.FAILED/'manifest.json')
        failure['artifacts_sha256'].pop(subject.rel(subject.FAILED/'coarse/proposal/newton_retry/iterations.json'))
        self.assertIsNone(subject.nonlinear_recovery(self.pool(),subject.FAILED,failure,self.meta,self.actual_pending(),self.meta['source_sha256']))

    def test_wrong_full_parent_or_pending_interval_rejected(self):
        failure=subject.read(subject.FAILED/'manifest.json');pending=self.actual_pending();pending['pending']['dt_s']/=2
        with self.assertRaises(ValueError):subject.nonlinear_recovery(self.pool(),subject.FAILED,failure,self.meta,pending,self.meta['source_sha256'])

    def test_modified_candidate_sha_rejected(self):
        original=subject.read;failure=original(subject.FAILED/'manifest.json')
        path=subject.FAILED/'coarse/proposal/broyden/last_UNCOMMITTED_candidate.npz';failure['artifacts_sha256'][subject.rel(path)]='0'*64
        with patch.object(subject,'read',side_effect=lambda p:failure if Path(p)==subject.FAILED/'manifest.json' else original(p)):
            with self.assertRaises(ValueError):subject.nonlinear_recovery(self.pool(),subject.FAILED,failure,self.meta,self.actual_pending(),self.meta['source_sha256'])

    def test_actual_history_replays_35_passes_four_multi_failures_and_final_stop(self):
        # Reuse actual saved states; only prior/full receipt dispatch is mocked.
        v=subject.v5;oldv5=subject.previous.previous.OLD_CAMPAIGN;old90=subject.previous.OLD_CAMPAIGN
        session=subject.read(subject.OLD_CAMPAIGN/'protocol.json');frames=[]
        observation,state=subject.frozen.observe(subject.read(v.COMMON_PARENT/'manifest.json'),frames)
        folders=[oldv5/'sensitivity_quarter']+[oldv5/f'attempt_{i:04d}' for i in range(1,8)]+[old90/f'attempt_{i:04d}' for i in range(1,14)]
        for folder in folders:observation,state=subject.frozen.observe(subject.read(folder/'manifest.json'),frames)
        seed=dict(template=subject.read(subject.OLD_CAMPAIGN/'attempt_0001/protocol.json'),refs={},parent=subject.previous.PARENT,
            meta=subject.read(subject.previous.PARENT/'manifest.json'),costs=session['inherited_cost_evidence'],
            policy=session['restored_adaptive_state'],state=state,frames=frames,observation=observation,
            previous_failed_costs=session['inherited_failed_costs'])
        with patch.object(subject.previous,'restore_history',return_value=seed),patch.object(subject.frozen,'validated_stage',
                side_effect=lambda pool,folder,*args:subject.read(folder/'manifest.json')) as validated:
            result=subject.restore_history(self.pool(),SimpleNamespace(maximum_dt=.000474609375,target_time=.22,shoulder_observe_ms=5.))
        self.assertEqual(validated.call_count,35);self.assertEqual(len(result['frames']),120);self.assertEqual(len(result['costs']),59)
        self.assertEqual(len(result['previous_failed_costs']),6);self.assertEqual(result['policy']['parent_sha256'],subject.PARENT_SHA256)
        self.assertEqual(result['policy']['nominal_dt_s'],.0002373046875);self.assertEqual(result['policy']['passes'],56)
        self.assertEqual(result['policy']['attempts'],62)


if __name__=='__main__':unittest.main()

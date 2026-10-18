"""SCHEMA_FIXTURE_NOT_PHYSICAL: provenance and attempt chronology only.

The frozen parent physics reader is mocked. These synthetic records are not
stress/contact solutions and are never counted as physical benchmark results.
"""
import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import inspect_process_predictor_onset as reader


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')


class PredictorOnsetMetadataTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name)
        self.onset = reader.read_json(ROOT/reader.ONSET_CONFIG)
        self.predictor = reader.read_json(ROOT/reader.PREDICTOR_CONFIG)
        self.recovery = reader.read_json(ROOT/reader.RECOVERY_CONFIG)
        self.cfg = reader.read_json(ROOT/reader.BASELINE)
        self.cfg['plastic_parent']['numerical']['maximum_iterations'] = 14
        self.cfg['recovery_variant'] = copy.deepcopy(self.recovery)
        self.cfg['onset_variant'] = copy.deepcopy(self.onset)
        self.cfg['predictor_variant'] = copy.deepcopy(self.predictor)
        self.sources = {s:reader.digest(ROOT/s) for s in reader.BASELINE_SOURCES}
        self.sidecar = dict(variant=self.predictor['variant'], production_authorized=False,
            configuration=copy.deepcopy(self.predictor), baseline_source_sha256=self.sources.copy(),
            source_sha256={s:reader.digest(ROOT/s) for s in reader.SUPPLEMENTAL_SOURCES},
            onset_source_sha256={s:reader.digest(ROOT/s) for s in reader.ONSET_SOURCES},
            evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
        self.make_path(4)
        delegate = patch.object(reader.parent, 'inspect', return_value={
            'status':'PASS_COMPLETE_READBACK', 'evidence_kind':'SCHEMA_FIXTURE_NOT_PHYSICAL'})
        self.mock_parent = delegate.start()
        self.addCleanup(delegate.stop)

    def make_path(self, n):
        self.cfg['increments_per_phase'] = n
        self.report = dict(status='J2_CONE_PATH_COMPLETE', error=None,
            increments_per_phase=n, accepted_steps=4*n, accepted_time_s=.2, contact_rows=9,
            source_sha256=self.sources.copy(), artifact_sha256={},
            records=[dict(step=0, physical_time_s=0., event_kind='INITIAL')],
            evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
        self.attempts = []
        base = 0.
        for step, t in enumerate(np.linspace(0.,.2,4*n+1)[1:], 1):
            def branch(stage, time):
                return dict(stage=stage,time_s=time,released=[False]*9,
                    status='ROOT_RETURNED_NOT_ACCEPTED',cone_margin_N=[1.]*9,
                    relative_increment=[[0.,0.]]*9,iterations=[dict(stage='CONE_ROOT',iteration=0,root_relative=0.,normal_force_N=[1.]*9)])
            trace = [branch('BOUNDARY_PROBE',base+(t-base)*1e-4),
                dict(stage='CLASSIFICATION',base_time_s=base,eligible=[False]*9,
                     initial_margin_N=[1.]*9,resolved_inward=[False]*9,released=[False]*9),
                branch('TARGET',float(t)),dict(stage='ACCEPTED',time_s=float(t),step=step,event_kind='REQUESTED')]
            self.attempts.append(dict(classification='ATTEMPT_TRACE_NOT_RESTART',restart_allowed=False,
                requested_time_s=float(t),base_time_s=base,base_step=step-1,
                last_accepted_time_s=float(t),last_accepted_step=step,error=None,trace=trace))
            self.report['records'].append(dict(step=step,physical_time_s=float(t),event_kind='REQUESTED'))
            base=float(t)
        self.write()

    def write(self):
        for old in self.path.glob('attempt_*.json'):
            old.unlink()
        save(self.path/'resolved_config.json',self.cfg)
        save(self.path/'predictor_provenance.json',self.sidecar)
        self.report['artifact_sha256'] = {'predictor_provenance.json':reader.digest(self.path/'predictor_provenance.json')}
        for attempt in self.attempts:
            name=f"attempt_{attempt['base_step']:04d}_{attempt['requested_time_s']:.17g}.json"
            save(self.path/name,attempt)
            self.report['artifact_sha256'][name]=reader.digest(self.path/name)
        save(self.path/'report.json',self.report)

    def rejected(self, pattern):
        with self.assertRaisesRegex(ValueError,pattern):
            reader.inspect(self.path)
        self.mock_parent.assert_not_called()

    def test_complete_metadata_trace_delegates_to_frozen_physics_reader(self):
        result=reader.inspect(self.path)
        self.assertEqual(result['status'],'PASS_COMPLETE_READBACK')
        self.assertEqual(result['evidence_kind'],'SCHEMA_FIXTURE_NOT_PHYSICAL')
        audit=result['predictor_onset_audit']['attempt_trace']
        self.assertEqual((audit['attempts'],audit['accepted_trace_rows'],audit['returned_candidate_roots']),(16,16,32))
        self.mock_parent.assert_called_once_with(self.path)

    def test_eight_grid_only_other_permitted_case(self):
        self.make_path(8)
        self.assertEqual(reader.inspect(self.path)['predictor_onset_audit']['attempt_trace']['attempts'],32)

    def test_source_omission_and_baseline_digest_are_rejected(self):
        self.sidecar['source_sha256'].pop(reader.PREDICTOR_CONFIG)
        self.write();self.rejected('supplemental source inventory')
        self.sidecar['source_sha256']={s:reader.digest(ROOT/s) for s in reader.SUPPLEMENTAL_SOURCES}
        self.sidecar['baseline_source_sha256'][next(iter(self.sources))]='0'*64
        self.write();self.rejected('baseline source inventory')

    def test_changed_physics_tolerance_budget_and_variant_are_rejected(self):
        original=copy.deepcopy(self.cfg)
        for section,key in [('parameters','young'),('friction_parent','numerical'),('plastic_parent','numerical')]:
            self.cfg=copy.deepcopy(original)
            if section=='parameters':self.cfg[section][key]['value']*=.9
            elif section=='friction_parent':self.cfg[section][key]['force_tolerance_N']*=10
            else:self.cfg[section][key]['maximum_iterations']=15
            self.write();self.rejected('configuration changed outside')
        self.cfg=copy.deepcopy(original);del self.cfg['onset_variant']
        self.write();self.rejected('configuration changed outside')

    def test_sidecar_hash_and_inventory_are_required(self):
        self.sidecar['note']='after recorded digest'
        save(self.path/'predictor_provenance.json',self.sidecar)
        self.rejected('missing from artifact inventory or changed')
        self.write();del self.report['artifact_sha256']['predictor_provenance.json']
        save(self.path/'report.json',self.report)
        self.rejected('missing from artifact inventory')

    def test_no_attempt_or_missing_file_is_rejected(self):
        one=next(self.path.glob('attempt_*.json'));one.unlink()
        self.rejected('attempt trace inventory')
        self.attempts=[];self.write()
        self.rejected('attempt trace inventory')

    def test_trace_tamper_rejected_before_parent(self):
        one=next(self.path.glob('attempt_*.json'))
        data=reader.read_json(one);data['note']='post-manifest'
        save(one,data)
        self.rejected('attempt trace artifact changed')

    def test_attempt_cannot_claim_restart_or_accepted_classification(self):
        self.attempts[0]['restart_allowed']=True
        self.write();self.rejected('attempt classification')
        self.attempts[0]['restart_allowed']=False
        self.attempts[0]['classification']='ACCEPTED'
        self.write();self.rejected('attempt classification')

    def test_root_return_cannot_be_promoted_or_nonconverged(self):
        root=self.attempts[0]['trace'][2]
        root['status']='ACCEPTED'
        self.write();self.rejected('candidate promoted')
        root['status']='ROOT_RETURNED_NOT_ACCEPTED';root['iterations'][0]['root_relative']=.1
        self.write();self.rejected('lacks converged cone trace')

    def test_classification_must_follow_coupled_probe_and_arithmetic(self):
        classification=self.attempts[0]['trace'][1]
        classification['released'][0]=True
        self.write();self.rejected('classification arithmetic')

    def test_accepted_trace_step_time_must_match_saved_report(self):
        self.report['records'][1]['physical_time_s']+=1e-8
        self.write();self.rejected('accepted trace/report mismatch')

    def test_blocked_root_cannot_be_followed_by_acceptance(self):
        self.attempts[0]['trace'][2].update(status='BLOCKED',error='synthetic root failure')
        self.write();self.rejected('blocked root followed by promoted')

    def test_rejected_candidate_cannot_be_claimed_as_attempt_endpoint(self):
        self.attempts[0]['last_accepted_time_s']+=1e-8
        self.write();self.rejected('endpoint promoted rejected')

    def test_blocked_partial_trace_exit_two_and_output_no_overwrite(self):
        self.attempts=self.attempts[:1]
        attempt=self.attempts[0]
        attempt.update(last_accepted_step=0,last_accepted_time_s=0.,error='synthetic commit failure')
        attempt['trace']=attempt['trace'][:-1]
        self.report.update(status='BLOCKED',error=attempt['error'],accepted_steps=0,accepted_time_s=0.,
                           records=self.report['records'][:1])
        self.write()
        self.mock_parent.return_value=dict(status='PASS_PARTIAL_READBACK',evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
        output=self.path/'partial_readback.json'
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(reader.main(['--run',str(self.path),'--output',str(output)]),2)
        self.assertEqual(reader.read_json(output)['status'],'PASS_PARTIAL_READBACK')
        with self.assertRaises(FileExistsError):
            reader.main(['--run',str(self.path),'--output',str(output)])

    @staticmethod
    def predictor_trace():
        return [
            dict(stage='INTERIOR_PREDICTOR',corrections_used=0,
                 interior_force_relative=.2,normal_force_N=[-1.]*9,
                 tangent_condition=2.,line_search_factor=1.),
            dict(stage='CONE_ROOT',iteration=1,root_relative=0.,normal_force_N=[1.]*9)]

    def test_predictor_then_converged_cone_preserves_shared_counter(self):
        self.attempts[0]['trace'][2]['iterations']=self.predictor_trace()
        self.write()
        result=reader.inspect(self.path)
        self.assertEqual(result['predictor_onset_audit']['shared_correction_budget'],14)
        self.mock_parent.assert_called_once_with(self.path)

    def test_cone_only_without_predictor_is_allowed(self):
        trace=self.attempts[0]['trace'][2]['iterations']
        trace.insert(0,dict(stage='CONE_ROOT',iteration=0,root_relative=.1,
            normal_force_N=[1.]*9,tangent_condition=2.,line_search_factor=.5))
        trace[1]['iteration']=1
        self.write()
        self.assertEqual(reader.inspect(self.path)['status'],'PASS_COMPLETE_READBACK')

    def test_predictor_cannot_reset_cone_budget(self):
        trace=self.predictor_trace();trace[-1]['iteration']=0
        self.attempts[0]['trace'][2]['iterations']=trace
        self.write();self.rejected('shared correction trace')

    def test_predictor_cannot_reappear_after_cone(self):
        trace=self.predictor_trace()
        cone=dict(stage='CONE_ROOT',iteration=0,root_relative=.1,
                  normal_force_N=[1.]*9,tangent_condition=2.,line_search_factor=1.)
        trace[0]['corrections_used']=1;trace[1]['iteration']=2
        self.attempts[0]['trace'][2]['iterations']=[cone]+trace
        self.write();self.rejected('predictor order')

    def test_predictor_cannot_be_reported_as_converged_root(self):
        trace=self.predictor_trace()[:1]
        trace[0].pop('line_search_factor')
        self.attempts[0]['trace'][2]['iterations']=trace
        self.write();self.rejected('lacks converged cone trace')
        trace[0]['root_relative']=0.
        self.write();self.rejected('predictor promoted to root')

    def test_shared_budget_cannot_exceed_fourteen(self):
        trace=[dict(stage='INTERIOR_PREDICTOR',corrections_used=i,
                    interior_force_relative=.2,normal_force_N=[-1.]*9,
                    tangent_condition=2.,line_search_factor=1.) for i in range(15)]
        self.attempts[0]['trace'][2]['iterations']=trace
        self.write();self.rejected('shared budget overflow')
        trace.append(dict(stage='CONE_ROOT',iteration=15,root_relative=0.,normal_force_N=[1.]*9))
        self.write();self.rejected('trace/budget')

    def test_normals_and_correction_linesearch_are_checked(self):
        root=self.attempts[0]['trace'][2]
        for change,pattern in (
            (lambda tr:tr[1].update(normal_force_N=[-1.]*9),'invalid cone root'),
            (lambda tr:tr[0].update(normal_force_N=[1.]*9),'invalid predictor'),
            (lambda tr:tr[0].update(line_search_factor=.3),'line search'),
            (lambda tr:tr[0].update(tangent_condition=1e13),'line search'),
            (lambda tr:tr[0].pop('line_search_factor'),'without successful correction')):
            trace=self.predictor_trace();change(trace);root['iterations']=trace
            self.write();self.rejected(pattern)

    def test_missing_predictor_helper_or_onset_source_is_rejected(self):
        original=copy.deepcopy(self.sidecar)
        self.sidecar['source_sha256'].pop('scripts/probe_process_feasible_predictor.py')
        self.write();self.rejected('supplemental source inventory')
        self.sidecar=copy.deepcopy(original)
        self.sidecar['onset_source_sha256'].pop(reader.ONSET_CONFIG)
        self.write();self.rejected('onset source inventory')

    def test_predictor_metadata_cannot_change_budget_or_timeout(self):
        for key,value in [('maximum_total_corrections',15),('run_timeout_s',600)]:
            self.cfg['predictor_variant']=copy.deepcopy(self.predictor)
            self.cfg['predictor_variant'][key]=value
            self.write();self.rejected('configuration changed outside')

    def test_physics_reader_failure_is_not_swallowed_or_promoted(self):
        self.mock_parent.side_effect=ValueError('physical power identity failed')
        with self.assertRaisesRegex(ValueError,'physical power identity failed'):
            reader.inspect(self.path)
        self.mock_parent.assert_called_once_with(self.path)

    def test_mid_predictor_failure_remains_blocked_partial(self):
        self.attempts=self.attempts[:1]
        attempt=self.attempts[0]
        attempt.update(last_accepted_step=0,last_accepted_time_s=0.,
                       error='synthetic predictor blocked')
        root=attempt['trace'][2]
        root.update(status='BLOCKED',error=attempt['error'],
                    iterations=self.predictor_trace()[:1])
        root['iterations'][0].pop('line_search_factor')
        attempt['trace']=attempt['trace'][:-1]
        self.report.update(status='BLOCKED',error=attempt['error'],accepted_steps=0,
                           accepted_time_s=0.,records=self.report['records'][:1])
        self.write()
        self.mock_parent.return_value=dict(status='PASS_PARTIAL_READBACK',
            evidence_kind='SCHEMA_FIXTURE_NOT_PHYSICAL')
        self.assertEqual(reader.inspect(self.path)['status'],'PASS_PARTIAL_READBACK')


if __name__=='__main__':
    unittest.main(verbosity=2)

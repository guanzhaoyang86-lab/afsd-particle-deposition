"""Bounded administration checks; all physical solvers remain uncalled."""
import ast
import copy
from datetime import timedelta
import inspect
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import run_gpu_event_resume_8h_20261011 as subject
import afsd_event_resume_8h_worker_20261011 as worker


class AdministrationTests(unittest.TestCase):
    def cost(self,seconds=590.5056227,single=113.8501983):
        return dict(maximum_candidate_wall_s=seconds/2,includes_event_probes=True,
            candidate_wall_s=[seconds/2,seconds/4,seconds/4],candidate_total_wall_s=seconds,
            worker_wall_s=seconds,supervised_wall_s=seconds,measured_pair_wall_s=seconds,event_count=1,
            maximum_single_solve_wall_s=single,single_solve_cost_is_proxy=False)

    def meta(self):
        return dict(source_sha256={'physical':'a'*64},event_protocol=subject.frozen.PROTOCOL,
            event_locator_protocol=subject.LOCATOR_PROTOCOL,blocked_prefix_protocol=subject.BLOCKED_PREFIX_PROTOCOL,
            visible_path_protocol=subject.VISIBLE_PATH_PROTOCOL,event_tolerance_scale=.25,
            screened_continuation_sha256=subject.PARENT_SHA256,screened_end_time_s=subject.PARENT_TIME_S,
            screened_end_step=527,screened_continuation_checkpoint='fixture/parent.npz')

    def protocol(self,costs=None,remaining=4000.):
        meta=self.meta();template=dict(source_sha256=meta['source_sha256'],session_deadline_epoch_s=subject.ABSOLUTE_STOP.timestamp(),
            qualification_wall_reserve_s=1100.)
        return subject.stage_protocol(template,meta,subject.PARENT,.0002373046875,{'new-admin':'b'*64},{'session':'c'*64},
            {'history':'d'*64},costs or [self.cost()],remaining)

    def test_fixed_eight_hours_includes_preparation(self):
        self.assertEqual((subject.ABSOLUTE_STOP-subject.AUTHORIZED_AT).total_seconds(),28800.)
        self.assertEqual(subject.authorized_wall_seconds(28800.,subject.AUTHORIZED_AT+timedelta(seconds=600)),28200.)
        self.assertEqual(subject.authorized_wall_seconds(1200.,subject.AUTHORIZED_AT+timedelta(seconds=600)),1200.)

    def test_expired_early_or_invalid_authorization_rejected(self):
        for now in (subject.ABSOLUTE_STOP,subject.AUTHORIZED_AT-timedelta(seconds=1),subject.AUTHORIZED_AT.replace(tzinfo=None)):
            with self.assertRaises(ValueError):subject.authorized_wall_seconds(28800.,now)
        for value in (True,0.,-1.,28800.01,float('inf'),float('nan')):
            with self.assertRaises(ValueError):subject.authorized_wall_seconds(value,subject.AUTHORIZED_AT)

    def test_no_fixed_900_or_1200_forecast_stop(self):
        for measured in (590.5056227,800.,3000.):
            budget=subject.dynamic_budget([self.cost(measured)],10000.)
            self.assertEqual(budget['action'],'RUN');self.assertEqual(budget['forecast_s'],30+1.5*measured)
            self.assertFalse(budget['forecast_is_admission_limit'])

    def test_dynamic_resource_formula_and_caps(self):
        for measured,soft in ((590.,1260.),(1000.,1590.),(3000.,3600.)):
            b=subject.dynamic_budget([self.cost(measured)],10000.)
            self.assertEqual((b['physical_soft_deadline_s'],b['internal_limit_s'],b['external_limit_s']),(soft,soft+30,soft+120))

    def test_all_stage_limits_cap_at_actual_remaining(self):
        b=subject.dynamic_budget([self.cost(3000.)],600.)
        self.assertEqual(b['action'],'RUN');self.assertTrue(b['best_effort'])
        self.assertEqual((b['physical_soft_deadline_s'],b['internal_limit_s'],b['external_limit_s']),(600.,600.,600.))

    def test_measured_single_reserve_controls_last_useful_stage(self):
        limit=1.5*113.8501983+20.
        self.assertEqual(subject.dynamic_budget([self.cost()],limit)['action'],'RUN')
        self.assertEqual(subject.dynamic_budget([self.cost()],limit-.01)['action'],'STOP_WALL_RESERVE')

    def test_large_old_cost_cannot_be_hidden_by_new_cheap_stage(self):
        b=subject.dynamic_budget([self.cost(3000.),self.cost(100.)],5000.)
        self.assertEqual(b['forecast_s'],4530.)
        self.assertEqual(len(b['event_cost_envelope']['evidence']),2)
        self.assertIs(subject.measured_envelope,subject.previous.measured_envelope)
        self.assertIs(subject.cost_evidence,subject.previous.cost_evidence)

    def test_nonfinite_or_inconsistent_cost_rejected(self):
        for remaining in (True,-1.,float('nan'),float('inf')):
            with self.assertRaises(ValueError):subject.dynamic_budget([self.cost()],remaining)
        row=self.cost();row['supervised_wall_s']+=1.
        with self.assertRaises(ValueError):subject.dynamic_budget([row],4000.)

    def test_phase_and_dt_gates_preserved(self):
        phases=[dict(start=0.,end=.22,stationary=True,name='fill')]
        self.assertEqual(subject.cost_decision(self.meta(),.0004,4000.,phases,[self.cost()])['action'],'RUN')
        meta=self.meta();meta['screened_end_time_s']=.22
        self.assertEqual(subject.cost_decision(meta,.0004,4000.,phases,[self.cost()])['action'],'STOP_REGISTERED_STATIONARY_FILL_END')
        with self.assertRaises(ValueError):subject.cost_decision(self.meta(),1e-8,4000.,phases,[self.cost()])

    def test_protocol_separates_administrative_and_physical_sources(self):
        p=self.protocol();self.assertEqual(p['source_sha256'],self.meta()['source_sha256'])
        self.assertIn('new-admin',p['reference_sha256']);self.assertNotIn('new-admin',p['source_sha256'])
        self.assertIsNone(p['implementation_transition']);self.assertEqual(p['event_tolerance_scale'],.25)
        self.assertEqual(p['resume_contract'],subject.rel(subject.CONTRACT))
        self.assertEqual(p['sensitivity_contract'],subject.rel(subject.v5.CONTRACT))

    def test_protocol_metadata_matches_assigned_dynamic_timers(self):
        p=self.protocol(costs=[self.cost(3000.)],remaining=600.)
        self.assertEqual(p['forecast_s'],4530.);self.assertIsNone(p['forecast_limit_s'])
        self.assertFalse(p['forecast_is_admission_limit']);self.assertNotIn('qualification_wall_reserve_s',p)
        self.assertEqual(p['candidate_forecast_s'],113.8501983)
        for key in ('internal_limit_s','external_limit_s','physical_soft_deadline_s','worker_deadline_s'):self.assertEqual(p[key],600.)

    def test_new_provider_or_old_locator_rejected(self):
        meta=self.meta();meta['event_locator_protocol']='wrong'
        with self.assertRaises(ValueError):subject.stage_protocol(dict(source_sha256=meta['source_sha256']),meta,subject.PARENT,
            .0004,{}, {}, {}, [self.cost()],4000.)

    def test_worker_absolute_epoch_caps_shorter_session(self):
        p=self.protocol();epoch=subject.ABSOLUTE_STOP.timestamp();p['absolute_deadline_epoch_s']=epoch-1000
        soft,hard=worker.worker_limits(p,10.,epoch-1010,11.)
        self.assertEqual((soft,hard),(21.,10.))

    def test_worker_rejects_extended_epoch_invalid_timers_and_metadata(self):
        base=self.protocol();epoch=subject.ABSOLUTE_STOP.timestamp()
        for key,value in (('absolute_deadline_epoch_s',epoch+1),('internal_limit_s',3631.),('physical_soft_deadline_s',True),
                          ('external_limit_s',100.),('worker_deadline_s',1259.)):
            p=copy.deepcopy(base);p[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):worker.worker_limits(p,0.,epoch-4000,0.)
        with self.assertRaises(ValueError):worker.worker_limits(base,0.,epoch,0.)

    def test_worker_physics_and_live_grants_ast_unchanged(self):
        self.assertFalse(subject.static_audit()['physical_source_changed'])
        self.assertIn('run_gpu_visibility_event_delay_20261010',Path(worker.__file__).read_text())

    def test_source_and_rounded_parent_identity_rejected(self):
        subject.verify_parent(self.meta())
        for key,value in (('screened_continuation_sha256','f'*64),('screened_end_step',526),('screened_end_time_s',.0813432617187497)):
            meta=self.meta();meta[key]=value
            with self.assertRaises(ValueError):subject.verify_parent(meta)

    def test_restore_calls_old_full_preflight_once_without_new_gpu_stages(self):
        calls=[ast.unparse(n.func) for n in ast.walk(ast.parse(inspect.getsource(subject.restore_history))) if isinstance(n,ast.Call)]
        self.assertEqual(calls.count('previous.restore_history'),1);self.assertNotIn('run_stage',calls)
        self.assertNotIn('frozen.LockedHistory',calls);self.assertIn('range(1,14)',inspect.getsource(subject.restore_history))
        self.assertIn("policy==summary['adaptive_state']",inspect.getsource(subject.restore_history))
        self.assertIn("canonical(old_event)==old_receipt['adaptive_event']",inspect.getsource(subject.restore_history))

    def test_shoulder_observation_is_milestone_and_never_early_stop(self):
        code=inspect.getsource(subject.continue_selected)
        self.assertIn('shoulder_milestone=True',code);self.assertNotIn('STOPPED_SHOULDER_OBSERVATION_WINDOW_COMPLETE',code)
        self.assertNotIn('shoulder_observe_ms/1000',code)
        self.assertIn('STOPPED_STATIONARY_TARGET_REACHED',code)

    def test_failure_dispatch_keeps_original_geometry_and_temporal_checks(self):
        code=inspect.getsource(subject.continue_selected)
        self.assertLess(code.index('multi_event_recovery('),code.index('frozen.geometry_budget_recovery('))
        self.assertLess(code.index('frozen.geometry_budget_recovery('),code.index('adaptive.resolve_failure('))
        self.assertLess(code.index('new=frozen.validated_stage('),code.index('parent,meta=child,new'))

    def pending(self):
        return subject.adaptive.begin_attempt(subject.adaptive.initial_state(subject.PARENT_SHA256,.0004),.0004)

    def timeout_fixture(self,folder,reason='WORKER_HARD_TIMEOUT'):
        p=self.protocol();p['dt_s']=.0004
        supervision=dict(termination_reason=reason,source_error=None,exit_code=124,interrupted=True,wall_seconds=p['external_limit_s'])
        (folder/'HARD_TIMEOUT.json').write_text(json.dumps(dict(status='TIMEOUT_NOT_ACCEPTED')))
        (folder/'supervision.json').write_text(json.dumps(supervision))
        return Mock(),p

    def test_hard_timeout_halves_same_full_parent_only(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);pool,p=self.timeout_fixture(folder)
            policy,event=subject.timeout_recovery(pool,folder,{},self.meta(),self.pending(),p)
            self.assertEqual(policy['nominal_dt_s'],.0002);self.assertEqual(policy['parent_sha256'],subject.PARENT_SHA256)
            self.assertEqual(policy['status'],'READY');self.assertEqual(event['action'],'REFINE');self.assertEqual(policy['passes'],0)

    def test_recorded_physics_failure_blocks_timeout_recovery(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);pool,p=self.timeout_fixture(folder)
            self.assertIsNone(subject.timeout_recovery(pool,folder,{'error':'physics failed'},self.meta(),self.pending(),p))
            (folder/'WORKER_PRE_END_FAILURE.json').write_text('{}')
            self.assertIsNone(subject.timeout_recovery(pool,folder,{},self.meta(),self.pending(),p))

    def test_unknown_or_authorization_timeout_cannot_refine(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            for reason in ('SESSION_DEADLINE','USER_STOP',None):
                pool,p=self.timeout_fixture(folder,reason)
                self.assertIsNone(subject.timeout_recovery(pool,folder,{},self.meta(),self.pending(),p))

    def test_timeout_marker_and_parent_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);pool,p=self.timeout_fixture(folder);(folder/'HARD_TIMEOUT.json').unlink()
            with self.assertRaises(ValueError):subject.timeout_recovery(pool,folder,{},self.meta(),self.pending(),p)
            pool,p=self.timeout_fixture(folder);p['parent_sha256']='f'*64
            with self.assertRaises(ValueError):subject.timeout_recovery(pool,folder,{},self.meta(),self.pending(),p)

    def test_timeout_cannot_exceed_shared_eight_halvings(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);pool,p=self.timeout_fixture(folder);pending=self.pending();pending['consecutive_halvings']=8
            policy,event=subject.timeout_recovery(pool,folder,{},self.meta(),pending,p)
            self.assertEqual(event['action'],'STOP');self.assertEqual(policy['status'],'STOPPED')

    def test_eof_after_hard_exit_keeps_timeout_recognizable(self):
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)/'attempt';server=Mock();channel=Mock();server.poll.return_value=True;server.recv.side_effect=EOFError
            child=Mock();child.pid=123;child.exitcode=124;child.is_alive.side_effect=[True,False,False]
            child.start.side_effect=lambda:(output/'HARD_TIMEOUT.json').write_text(json.dumps(dict(status='TIMEOUT_NOT_ACCEPTED')))
            context=Mock();context.Pipe.return_value=(server,channel);context.Process.return_value=child
            pool=Mock();pool.sha.return_value='a'*64;pool.stats.return_value={}
            p=dict(source_sha256={},reference_sha256={},external_limit_s=1380.)
            with patch.object(subject.multiprocessing,'get_context',return_value=context):
                code,_=subject.run_stage(pool,output,p,subject.time.perf_counter()+4000,Mock(),Path(temp)/'STOP')
            result=json.loads((output/'supervision.json').read_text());self.assertEqual(code,124)
            self.assertIsNone(result['source_error']);self.assertEqual(result['termination_reason'],'WORKER_HARD_TIMEOUT')
            server.recv.assert_called_once()

    def test_detach_and_singleton_remain(self):
        code=inspect.getsource(subject.main)
        for part in ('CREATE_BREAKAWAY_FROM_JOB','DETACHED_PROCESS','startup.wShowWindow = 0',"require(not policy['kill_on_job_close']"):
            self.assertIn(part,code)
        self.assertIn('admin.singleton_and_awake()',inspect.getsource(subject.run))
        self.assertIn('pool.final_audit()',inspect.getsource(subject.run))

    def test_no_rebinding_of_frozen_implementation(self):
        tree=ast.parse(Path(subject.__file__).read_text())
        self.assertFalse(any(isinstance(n,(ast.Global,ast.Nonlocal)) for n in ast.walk(tree)))
        for n in ast.walk(tree):
            if isinstance(n,(ast.Assign,ast.AnnAssign,ast.AugAssign)):
                for target in n.targets if isinstance(n,ast.Assign) else [n.target]:
                    self.assertFalse(isinstance(target,ast.Attribute) and ast.unparse(target).split('.')[0] in
                        ('previous','v5','frozen','bridge','admin','adaptive'))

    def actual_pool(self):
        return SimpleNamespace(sha=subject.bridge.core.sha,seal=lambda mapping:None,verify=lambda mapping:None)

    def test_actual_saved_history_reconstructs_old_pending_then_new_refinement(self):
        # Expensive original full receipt verification is exercised by startup;
        # this fixture independently exercises restoration using actual saved states.
        v=subject.v5;old=subject.previous.OLD_CAMPAIGN;old_session=subject.read(subject.OLD_CAMPAIGN/'protocol.json')
        frames=[];observation,state=subject.frozen.observe(subject.read(v.COMMON_PARENT/'manifest.json'),frames)
        for folder in [old/'sensitivity_quarter']+[old/f'attempt_{i:04d}' for i in range(1,8)]:
            observation,state=subject.frozen.observe(subject.read(folder/'manifest.json'),frames)
        seed=dict(template=subject.read(subject.OLD_CAMPAIGN/'attempt_0001/protocol.json'),refs={},
            parent=subject.previous.PARENT,meta=subject.read(subject.previous.PARENT/'manifest.json'),
            costs=old_session['inherited_cost_evidence'],policy=old_session['restored_adaptive_state'],
            state=state,frames=frames,observation=observation)
        args=SimpleNamespace(maximum_dt=.000474609375,target_time=.22,shoulder_observe_ms=5.)
        with patch.object(subject.previous,'restore_history',return_value=seed),patch.object(subject.frozen,'validated_stage',
                side_effect=lambda pool,folder,*args:subject.read(folder/'manifest.json')) as validation:
            restored=subject.restore_history(self.actual_pool(),args)
        self.assertEqual(validation.call_count,13);self.assertEqual(len(restored['costs']),24)
        self.assertEqual(len(restored['frames']),48);self.assertEqual(restored['policy']['passes'],21)
        self.assertEqual(restored['policy']['attempts'],22);self.assertEqual(restored['policy']['nominal_dt_s'],.0002373046875)
        self.assertEqual(restored['policy']['parent_sha256'],subject.PARENT_SHA256)
        self.assertIsNone(restored['policy']['pending']);self.assertEqual(restored['initial_recovery']['action'],'REFINE')
        self.assertEqual(len(restored['previous_failed_costs']),1)
        self.assertFalse(restored['previous_failed_costs'][0]['complete_pair'])

    def actual_failure(self):
        return (subject.read(subject.FAILED/'manifest.json'),subject.read(subject.PARENT/'manifest.json'),
            subject.read(subject.OLD_CAMPAIGN/'final_summary.json')['adaptive_state'])

    def test_unregistered_nested_report_layout_is_not_recoverable(self):
        failure,parent,policy=self.actual_failure()
        name=subject.rel(subject.FAILED/'coarse/proposal/union_visibility_events.json')
        failure['artifacts_sha256'].pop(name)
        original=subject.read
        with patch.object(subject,'read',side_effect=lambda path:failure if Path(path)==subject.FAILED/'manifest.json' else original(path)):
            self.assertIsNone(subject.multi_event_recovery(self.actual_pool(),subject.FAILED,failure,parent,policy,parent['source_sha256']))

    def test_registered_but_tampered_report_still_raises(self):
        failure,parent,policy=self.actual_failure()
        name=subject.rel(subject.FAILED/'coarse/proposal/union_visibility_events.json')
        failure['artifacts_sha256'][name]='0'*64
        original=subject.read
        with patch.object(subject,'read',side_effect=lambda path:failure if Path(path)==subject.FAILED/'manifest.json' else original(path)):
            with self.assertRaisesRegex(ValueError,'evidence missing or changed'):
                subject.multi_event_recovery(self.actual_pool(),subject.FAILED,failure,parent,policy,parent['source_sha256'])


if __name__=='__main__':unittest.main()

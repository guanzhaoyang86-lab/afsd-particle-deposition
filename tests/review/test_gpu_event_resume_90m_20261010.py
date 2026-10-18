"""Administrative-only continuation boundary checks; no GPU or physical solve."""
import ast
import copy
from datetime import timedelta
import inspect
from pathlib import Path
import sys
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_event_resume_90m_20261010 as subject


class ResumeTests(unittest.TestCase):
    def cost(self, measured=590.5056227):
        return dict(maximum_candidate_wall_s=400., includes_event_probes=True,
            candidate_wall_s=[400., 100., measured-500.], candidate_total_wall_s=measured,
            worker_wall_s=measured, supervised_wall_s=measured, measured_pair_wall_s=measured,
            event_count=2, maximum_single_solve_wall_s=113.85019830000238,
            single_solve_cost_is_proxy=False)

    def meta(self):
        return dict(source_sha256={'physical':'a'*64}, event_protocol=subject.frozen.PROTOCOL,
            event_locator_protocol=subject.LOCATOR_PROTOCOL, blocked_prefix_protocol=subject.BLOCKED_PREFIX_PROTOCOL,
            visible_path_protocol=subject.VISIBLE_PATH_PROTOCOL, event_tolerance_scale=.25,
            screened_continuation_sha256=subject.PARENT_SHA256, screened_end_time_s=subject.PARENT_TIME_S,
            screened_end_step=499, screened_continuation_checkpoint='fixture/parent.npz')

    def policy(self):
        return dict(parent_sha256=subject.PARENT_SHA256, status='READY', pending=None, passes=8, attempts=8)

    def phases(self):
        return [dict(start=0., end=.22, stationary=True, name='fill')]

    def protocol(self, meta=None, template=None, costs=None, scale=.25):
        meta = self.meta() if meta is None else meta
        template = dict(source_sha256=meta['source_sha256'], forecast_limit_s=900., internal_limit_s=990.,
            external_limit_s=1080., physical_soft_deadline_s=960., worker_deadline_s=960., qualification_wall_reserve_s=1100.) if template is None else template
        return subject.stage_protocol(template, meta, subject.PARENT, .000474609375,
            {'new-admin':'c'*64}, {'session':'d'*64}, {'history':'e'*64}, scale, [self.cost()] if costs is None else costs)

    def test_authorization_includes_preparation_and_startup(self):
        now = subject.AUTHORIZED_AT+timedelta(seconds=600)
        self.assertEqual(subject.authorized_wall_seconds(5400., now), 4800.)
        self.assertEqual(subject.authorized_wall_seconds(3600., now), 3600.)
        self.assertEqual((subject.ABSOLUTE_STOP-subject.AUTHORIZED_AT).total_seconds(), 5400.)
        self.assertLess(subject.ABSOLUTE_STOP, subject.frozen.ABSOLUTE_STOP)

    def test_no_new_ninety_minutes_from_launch(self):
        now = subject.ABSOLUTE_STOP-timedelta(seconds=30)
        self.assertEqual(subject.authorized_wall_seconds(5400., now), 30.)

    def test_expired_or_pre_request_time_rejected(self):
        for now in (subject.ABSOLUTE_STOP, subject.ABSOLUTE_STOP+timedelta(seconds=1),
                    subject.AUTHORIZED_AT-timedelta(seconds=1)):
            with self.subTest(now=now), self.assertRaises(ValueError): subject.authorized_wall_seconds(5400., now)

    def test_bad_wall_limits_rejected(self):
        for value in (True, 0., -1., 5400.01, float('nan'), float('inf')):
            with self.subTest(value=value), self.assertRaises(ValueError):
                subject.authorized_wall_seconds(value, subject.AUTHORIZED_AT)

    def test_timezone_required(self):
        with self.assertRaises(ValueError): subject.authorized_wall_seconds(5400., subject.AUTHORIZED_AT.replace(tzinfo=None))

    def test_old_real_cost_stop_now_fits_explicit_new_budget(self):
        row = self.cost()
        result = subject.cost_decision(self.meta(), .000474609375, 4000., self.phases(), [row])
        self.assertEqual(result['action'], 'RUN')
        self.assertGreater(result['forecast_s'], 900.)
        self.assertAlmostEqual(result['forecast_s'], 30.+1.5*row['measured_pair_wall_s'])
        self.assertEqual(result['child_timeout_s'], 1380.)
        self.assertFalse(result['event_cost_envelope']['guarantees_unknown_events'])

    def test_new_forecast_boundary_and_original_twenty_second_reserve(self):
        for measured, remaining, expected in ((780., 1300., 'RUN'), (780.001, 2000., 'STOP_COST_REASSESSMENT'),
            (780., 1219.9, 'STOP_WALL_RESERVE'), (780., 1220., 'RUN')):
            with self.subTest(measured=measured, remaining=remaining):
                self.assertEqual(subject.cost_decision(self.meta(), .0004, remaining, self.phases(),
                    [self.cost(measured)])['action'], expected)

    def test_outer_limit_capped_by_remaining_authorization(self):
        result = subject.cost_decision(self.meta(), .0004, 1000., self.phases(), [self.cost()])
        self.assertEqual(result['child_timeout_s'], 1000.)

    def test_cost_model_includes_full_historical_maximum(self):
        expensive = self.cost(780.); cheaper = self.cost(550.)
        self.assertEqual(subject.measured_envelope([expensive, cheaper])['forecast_s'], 1200.)
        self.assertIs(subject.measured_envelope, subject.previous.measured_envelope)
        self.assertIs(subject.cost_evidence, subject.previous.cost_evidence)

    def test_inconsistent_hidden_supervision_cost_rejected(self):
        row = self.cost(); row['supervised_wall_s'] += 1
        with self.assertRaises(ValueError): subject.cost_decision(self.meta(), .0004, 4000., self.phases(), [row])

    def test_phase_end_and_minimum_step_gates_remain(self):
        meta = self.meta(); meta['screened_end_time_s'] = .22
        self.assertEqual(subject.cost_decision(meta, .0004, 4000., self.phases(), [self.cost()])['action'],
            'STOP_REGISTERED_STATIONARY_FILL_END')
        with self.assertRaises(ValueError): subject.cost_decision(self.meta(), 1e-8, 4000., self.phases(), [self.cost()])

    def test_protocol_physical_source_does_not_gain_admin_files(self):
        meta = self.meta(); before = copy.deepcopy(meta)
        protocol = self.protocol(meta=meta)
        self.assertEqual(protocol['source_sha256'], before['source_sha256'])
        self.assertEqual(meta, before)
        self.assertIn('new-admin', protocol['reference_sha256'])
        self.assertNotIn('new-admin', protocol['source_sha256'])
        self.assertEqual(protocol['parent_implementation_source_sha256'], before['source_sha256'])
        self.assertIsNone(protocol['implementation_transition'])

    def test_protocol_limit_metadata_and_contract_match_actual_timers(self):
        protocol = self.protocol()
        for key, value in (('forecast_limit_s',1200.),('internal_limit_s',1290.),('external_limit_s',1380.),
                           ('physical_soft_deadline_s',1260.),('worker_deadline_s',1260.)):
            self.assertEqual(protocol[key], value)
        self.assertEqual(protocol['administrative_protocol'], subject.PROTOCOL)
        self.assertEqual(protocol['resume_contract'], subject.rel(subject.CONTRACT))
        self.assertEqual(protocol['sensitivity_contract'], subject.rel(subject.previous.CONTRACT))
        self.assertEqual(protocol['resume_absolute_stop'], subject.ABSOLUTE_STOP.isoformat())
        self.assertNotIn('qualification_wall_reserve_s', protocol)
        self.assertEqual(protocol['candidate_forecast_s'], 113.85019830000238)

    def test_no_new_baseline_or_half_continuation_allowed(self):
        for scale in (1., .5):
            with self.assertRaises(ValueError): self.protocol(scale=scale)

    def test_new_physical_source_generation_is_rejected(self):
        with self.assertRaises(ValueError): self.protocol(template=dict(source_sha256={'other':'b'*64}))

    def test_parent_protocol_identity_cannot_change(self):
        meta = self.meta(); meta['event_locator_protocol'] = 'old-locator'
        with self.assertRaises(ValueError): self.protocol(meta=meta)

    def test_over_budget_stage_protocol_rejected(self):
        with self.assertRaises(ValueError): self.protocol(costs=[self.cost(781.)])

    def test_wrong_checkpoint_time_or_step_rejected(self):
        subject.verify_parent(self.meta(), self.policy())
        for key, value in (('screened_continuation_sha256','f'*64), ('screened_end_time_s',.076), ('screened_end_step',498)):
            meta = self.meta(); meta[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): subject.verify_parent(meta, self.policy())

    def test_pending_or_detached_policy_rejected(self):
        for key, value in (('parent_sha256','f'*64),('status','PENDING'),('pending',dict(dt_s=.0004))):
            policy = self.policy(); policy[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): subject.verify_parent(self.meta(), policy)

    def restoration(self):
        policy = self.policy(); observation = dict(frame_count=20, history=[dict(time_s=.075)])
        costs = [self.cost()]; summary = dict(adaptive_state=copy.deepcopy(policy),
            next_decision=dict(event_cost_envelope=subject.canonical(subject.measured_envelope(costs))))
        return policy, observation, costs, summary, copy.deepcopy(observation)

    def test_exact_rebuilt_policy_observation_cost_accepted(self):
        subject.verify_restoration(*self.restoration())

    def old_session(self):
        summary = dict(finished=True, exit_code=0, status='STOP_COST_REASSESSMENT', completed_stage_count=7,
            attempt_count=7, latest_pass_stage=subject.rel(subject.PARENT), checkpoint_sha256=subject.PARENT_SHA256,
            physical_time_s=subject.PARENT_TIME_S, screened_steps=subject.PARENT_STEP)
        audit = dict(failed=False, final_audits=1, unique_files=24355, expected_files=24355, locked_files=24355)
        return summary, audit

    def test_complete_old_session_audit_accepted(self):
        subject.verify_old_session(*self.old_session())

    def test_incomplete_original_final_audit_rejected(self):
        for key, value in (('failed',True),('final_audits',0),('unique_files',24354),('expected_files',24354),('locked_files',0)):
            summary, audit = self.old_session(); audit[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): subject.verify_old_session(summary,audit)

    def test_original_summary_must_identify_exact_last_parent(self):
        for key,value in (('checkpoint_sha256','f'*64),('screened_steps',498),('physical_time_s',.076),('attempt_count',8),
                          ('finished',False),('status','STOPPED_CHILD_FAILURE')):
            summary,audit = self.old_session(); summary[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): subject.verify_old_session(summary,audit)

    def test_policy_counters_cannot_reset_on_resume(self):
        values = self.restoration(); values[0]['passes'] = 0
        with self.assertRaisesRegex(ValueError,'policy'): subject.verify_restoration(*values)

    def test_observer_cannot_drop_old_fine_leaves(self):
        values = self.restoration(); values[1]['frame_count'] = 2
        with self.assertRaisesRegex(ValueError,'fine-leaf'): subject.verify_restoration(*values)

    def test_old_expensive_cost_cannot_be_discarded(self):
        values = self.restoration(); values[2][0] = self.cost(550.)
        with self.assertRaisesRegex(ValueError,'cost evidence'): subject.verify_restoration(*values)

    def test_worker_and_supervision_have_only_named_timer_changes(self):
        report = subject.static_audit()
        self.assertFalse(report['physical_source_changed'])
        self.assertEqual(report['grant_timeout_unchanged_s'],180.)

    def test_ordinary_loop_exactly_frozen_v5(self):
        self.assertEqual(ast.dump(ast.parse(inspect.getsource(subject.continue_selected))),
            ast.dump(ast.parse(inspect.getsource(subject.previous.continue_selected))))

    def test_one_preflight_and_no_gpu_during_restoration(self):
        tree=ast.parse(inspect.getsource(subject.restore_history))
        calls=[ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)]
        self.assertEqual(calls.count('frozen._initial_preflight'),1)
        self.assertNotIn('run_stage',calls)
        self.assertNotIn('frozen.LockedHistory',calls)
        self.assertIn('previous.compare_group',calls)
        self.assertEqual(calls.count('frozen.observe'),3)
        self.assertIn('range(1, 8)',inspect.getsource(subject.restore_history))
        run_calls=[ast.unparse(n.func) for n in ast.walk(ast.parse(inspect.getsource(subject.run))) if isinstance(n,ast.Call)]
        self.assertEqual(run_calls.count('frozen.LockedHistory'),1)

    def test_no_monkeypatch_of_frozen_modules(self):
        tree=ast.parse(Path(subject.__file__).read_text(encoding='utf-8'))
        self.assertFalse(any(isinstance(n,(ast.Global,ast.Nonlocal)) for n in ast.walk(tree)))
        for n in ast.walk(tree):
            if isinstance(n,(ast.Assign,ast.AnnAssign,ast.AugAssign)):
                for target in n.targets if isinstance(n,ast.Assign) else [n.target]:
                    self.assertFalse(isinstance(target,ast.Attribute) and ast.unparse(target).split('.')[0] in
                        ('previous','frozen','admin','adaptive','bridge'))

    def test_last_valid_parent_snapshot_survives_failed_child(self):
        code=inspect.getsource(subject.run)
        self.assertIn("admin.snapshot(read(PARENT/'manifest.json'), PARENT)",code)
        self.assertIn('initial_parent_verified_in_this_session=False',code)
        self.assertIn('initial_parent_verified_in_this_session=True',code)
        loop=inspect.getsource(subject.continue_selected)
        self.assertLess(loop.index('next_meta = frozen.validated_stage'),loop.index('parent, meta = child, next_meta'))

    def test_hidden_detach_preserves_original_job_policy(self):
        code=inspect.getsource(subject.main)
        for text in ('startup.wShowWindow = 0','CREATE_BREAKAWAY_FROM_JOB','DETACHED_PROCESS',
                     "require(not policy['kill_on_job_close']",'close_fds=True'):
            self.assertIn(text,code)


if __name__=='__main__': unittest.main()

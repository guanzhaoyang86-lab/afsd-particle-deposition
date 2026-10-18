"""One interrupted-campaign authorization; synthetic receipts, no physics."""
from pathlib import Path
from datetime import datetime, timedelta
import importlib.util
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import resume_gpu_campaign_20261008 as resume
import supervise_gpu_campaign_fast_20261008 as fast

spec = importlib.util.spec_from_file_location('frozen_fast_components',
    ROOT/'tests/review/test_gpu_fast_supervisor_20261008.py')
prior_tests = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior_tests)


class ResumeIdentityTests(unittest.TestCase):
    def test_entire_supervisor_loop_is_ast_identical(self):
        self.assertTrue(resume.static_loop_audit()['loop_ast_identical'])

    def test_original_deadline_and_physical_helpers_are_reused(self):
        self.assertEqual(resume.AUTHORIZED_STOP_AT, '2026-10-09T02:19:45+08:00')
        for name in ('remaining_seconds', 'budget_remaining', 'decision', 'verify_parent',
                     'snapshot', 'classify_failure', 'singleton_and_awake', 'RUNNER'):
            self.assertIs(getattr(resume, name), getattr(fast, name))

    def test_handoff_rejection_precedes_any_child_launch(self):
        with mock.patch.object(resume, 'validate_handoff', side_effect=ValueError('rejected origin')), \
             mock.patch.object(resume.subprocess, 'Popen') as launch:
            with self.assertRaisesRegex(ValueError, 'rejected origin'):
                resume.supervise(dict(source_sha256={}), {})
        launch.assert_not_called()


class InterruptionGateTests(unittest.TestCase):
    def setUp(self):
        tmp = ROOT/'tmp'; tmp.mkdir(exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(prefix='gpu_resume_test_', dir=tmp)
        self.root = Path(self.tmp.name).resolve()
        self.assertTrue(self.root.is_relative_to(tmp.resolve()))
        self.addCleanup(self.tmp.cleanup)
        self.origin = self.root/'fast'; self.origin.mkdir()
        self.original = self.root/'original'; self.original.mkdir()
        self.parent = self.origin/'attempt_0003'; self.parent.mkdir()
        self.unfinished = self.origin/'attempt_0004'; self.unfinished.mkdir()
        for name, value in (('INTERRUPTED', self.origin), ('ORIGINAL_ORIGIN', self.original),
                            ('PARENT', self.parent), ('UNFINISHED', self.unfinished)):
            self.enterContext(mock.patch.object(resume, name, value))
        self.process = self.enterContext(mock.patch.object(resume, 'process_is_running', return_value=False))
        source = self.root/'frozen_source.txt'; source.write_text('frozen physics', encoding='utf-8')
        self.source = {resume.rel(source): resume.sha(source)}
        stop = datetime.fromisoformat(resume.AUTHORIZED_STOP_AT)
        original = dict(schema='GPU_TEN_HOUR_SUPERVISOR_20261008', output=resume.rel(self.original),
            source_sha256=self.source, child_source_sha256=self.source, wall_limit_s=36000.,
            started_at=(stop-timedelta(hours=10)).isoformat(), latest_stop_at=resume.AUTHORIZED_STOP_AT)
        prior = dict(schema='GPU_FIXED_DEADLINE_FAST_SUPERVISOR_20261008', output=resume.rel(self.origin),
            source_sha256=self.source, child_source_sha256=self.source, origin_campaign=resume.rel(self.original),
            latest_stop_at=resume.AUTHORIZED_STOP_AT, supervisor_pid=11888)
        supplied = dict(eligible_crossed=dict(particle_count=0, mass_kg=0.),
            crossed_under_annular_shoulder=dict(particle_count=0, mass_kg=0.))
        self.meta = dict(source_sha256=self.source, screened_continuation_sha256=resume.PARENT_SHA256,
            screened_continuation_checkpoint=resume.rel(self.parent/'fine_2/NOT_RESTARTABLE_candidate.npz'),
            screened_end_step=156, screened_end_time_s=.0256953125, observer=dict(supplied_material=supplied))
        self.parent_verifier = self.enterContext(mock.patch.object(resume, 'verify_parent', return_value=self.meta))
        expected = resume.snapshot(self.meta, self.parent)
        progress = dict(expected, finished=False, status='RUNNING_GPU_CHILD', completed_stage_count=3,
            attempt_count=4, current_parent=resume.rel(self.parent), current_attempt=resume.rel(self.unfinished),
            latest_stop_at=resume.AUTHORIZED_STOP_AT, supervisor_pid=11888, child_pid=51152)
        interrupted = dict(source_sha256=self.source, parent_gpu_directory=resume.rel(self.parent),
            parent_sha256=resume.PARENT_SHA256)
        for path, value in ((self.original/'protocol.json', original), (self.origin/'protocol.json', prior),
                (self.origin/'progress.json', progress), (self.unfinished/'protocol.json', interrupted),
                (self.parent/'manifest.json', dict(synthetic=True))):
            resume.atomic_json(path, value)
        for path in (self.origin/'attempt_0003_launcher.log', self.origin/'attempt_0004_launcher.log', self.unfinished/'run.log'):
            path.write_text('', encoding='utf-8')
        receipt = dict(expected, verified=True, exit_code=0, timeout=False, termination=None,
            stage=resume.rel(self.parent), manifest_sha256=resume.sha(self.parent/'manifest.json'),
            log_sha256=resume.sha(self.origin/'attempt_0003_launcher.log'), wall_seconds=196.4)
        resume.atomic_json(self.origin/'attempt_0003_receipt.json', receipt)
        directory = self.unfinished/'coarse/broyden'; directory.mkdir(parents=True)
        resume.atomic_json(directory/'secant_failure.json', dict(reason='ACTUAL_SUPPORT_GRAPH_CHANGED'))

    def audit(self):
        return resume.read_interruption(self.origin, self.parent)

    def change(self, path, **fields):
        data = resume.read(path); data.update(fields); resume.atomic_json(path, data)

    def test_only_latest_complete_pass_is_used_and_fallback_evidence_is_preserved(self):
        _, _, meta, _, evidence = self.audit()
        self.assertIs(meta, self.meta)
        self.assertIn(resume.rel(self.unfinished/'coarse/broyden/secant_failure.json'), evidence)
        self.parent_verifier.assert_called_once_with(self.parent, self.source, self.origin/'attempt_0002')
        self.assertEqual(self.process.call_args_list, [mock.call(11888), mock.call(51152)])

    def test_either_live_process_blocks_resume(self):
        self.process.side_effect = lambda pid: pid == 51152
        with self.assertRaisesRegex(ValueError, 'still running'):
            self.audit()

    def test_wrong_or_stale_parent_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'exact stage-3 PASS'):
            resume.read_interruption(self.origin, self.origin/'attempt_0002')
        self.change(self.origin/'progress.json', latest_pass_stage=resume.rel(self.unfinished))
        with self.assertRaisesRegex(ValueError, 'missing-completion'):
            self.audit()

    def test_source_drift_is_rejected(self):
        (self.root/'frozen_source.txt').write_text('changed', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'hash changed'):
            self.audit()

    def test_parent_checkpoint_hash_cannot_change(self):
        self.change(self.unfinished/'protocol.json', parent_sha256='different')
        with self.assertRaisesRegex(ValueError, 'PASS parent changed'):
            self.audit()

    def test_timeout_json_stop_and_completed_receipts_are_not_interruptions(self):
        paths = (self.unfinished/'HARD_TIMEOUT.json', self.origin/'HARD_TIMEOUT.json',
            self.origin/'STOP_REQUESTED', self.unfinished/'manifest.json', self.unfinished/'supervision.json',
            self.unfinished/'pair_screen.json', self.origin/'attempt_0004_receipt.json', self.origin/'final_summary.json')
        for path in paths:
            with self.subTest(path=path.name):
                path.write_text('{}', encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'separate review'):
                    self.audit()
                path.unlink()

    def test_terminal_failure_or_external_timeout_in_logs_is_rejected(self):
        for value in ('FAIL_GPU_STAGE_NOT_ACCEPTED', 'EXTERNAL_TIMEOUT_NOT_ACCEPTED', 'HARD_TIMEOUT'):
            (self.unfinished/'run.log').write_text(value, encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'terminal failure'):
                self.audit()

    def test_outer_pass_receipt_must_be_verified(self):
        self.change(self.origin/'attempt_0003_receipt.json', verified=False)
        with self.assertRaisesRegex(ValueError, 'outer PASS receipt'):
            self.audit()

    def test_original_deadline_cannot_be_extended(self):
        stop = datetime.fromisoformat(resume.AUTHORIZED_STOP_AT)+timedelta(seconds=1.)
        self.change(self.original/'protocol.json', latest_stop_at=stop.isoformat())
        with self.assertRaisesRegex(ValueError, 'absolute deadline'):
            self.audit()


def load_tests(loader, standard_tests, pattern):
    standard_tests.addTests(loader.loadTestsFromTestCase(prior_tests.DeadlineTests))
    return standard_tests


if __name__ == '__main__':
    unittest.main(verbosity=2)

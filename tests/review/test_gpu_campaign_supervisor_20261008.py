"""Bounded GPU campaign controller tests; no physical runner is imported.

Synthetic receipts and child processes exercise supervision only. They do not
constitute GPU physics evidence or authorize a failed candidate as a parent.
"""
from pathlib import Path
import ast
import copy
import importlib.util
import io
import json
import math
import argparse
import contextlib
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'scripts/supervise_gpu_campaign_20261008.py'
sys.path.insert(0, str(ROOT/'scripts'))
SPEC = importlib.util.spec_from_file_location('campaign_supervisor_under_test', SOURCE)
campaign = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(campaign)


def metadata(time_s=.012, maximum_s=50.):
    return dict(status='PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION',
        screened_end_time_s=time_s, maximum_candidate_wall_s=maximum_s,
        dt_s=.000474609375, accepted=False, restart_allowed=False, diagnostic_only=True)


PHASES = [dict(start=0., end=.02, name='entry', stationary=True),
          dict(start=.02, end=.22, name='stationary_feed', stationary=True),
          dict(start=.22, end=.24, name='traverse', stationary=False)]
SCREEN_ERROR = "ValueError('original coarse/fine screen failed; no GPU continuation')"


def failed_time_screen():
    return dict(passed=False, accepted=False, global_screen={'pass': False},
        local_screen={'passed': True}, crossing_ids_match=True)


class DecisionTests(unittest.TestCase):
    def test_measured_three_candidate_forecast_and_650_second_cap(self):
        result = campaign.decision(metadata(), .000474609375, 36000., PHASES)
        self.assertEqual(result['action'], 'RUN')
        self.assertEqual(result['forecast_s'], 30.+4.5*50.)
        self.assertEqual(result['child_timeout_s'], 650.)
        self.assertEqual(result['dt_s'], .000474609375)

    def test_child_timeout_cannot_exceed_campaign_remaining_time(self):
        result = campaign.decision(metadata(), .0004, 300., PHASES)
        self.assertEqual(result['action'], 'RUN')
        self.assertLessEqual(result['child_timeout_s'], 300.)
        self.assertLessEqual(result['child_timeout_s'], 650.)

    def test_total_wall_reserve_stops_before_another_child(self):
        for remaining in (0., 1., 274.999999):
            result = campaign.decision(metadata(), .0004, remaining, PHASES)
            self.assertTrue(result['action'].startswith('STOP'))
        result = campaign.decision(metadata(), .0004, 275., PHASES)
        self.assertEqual(result['action'], 'RUN')

    def test_expensive_candidate_does_not_bypass_450_second_forecast(self):
        result = campaign.decision(metadata(maximum_s=100.), .0004, 36000., PHASES)
        self.assertTrue(result['action'].startswith('STOP'))

    def test_interval_is_clipped_to_current_phase(self):
        current = .02-2.**-14
        result = campaign.decision(metadata(time_s=current), .0004, 36000., PHASES)
        self.assertEqual(result['action'], 'RUN')
        self.assertEqual(result['dt_s'], 2.**-14)
        self.assertEqual(current+result['dt_s'], .02)
        self.assertEqual((current+result['dt_s']/2)+result['dt_s']/2, .02)

    def test_unrepresentable_exact_two_half_boundary_stops(self):
        result = campaign.decision(metadata(time_s=.0199), .0004, 36000., PHASES)
        self.assertEqual(result['action'], 'STOP_UNREPRESENTABLE_PHASE_BOUNDARY')

    def test_exact_boundary_uses_next_phase_without_growing_dt(self):
        result = campaign.decision(metadata(time_s=.02), .0004, 36000., PHASES)
        self.assertEqual(result['action'], 'RUN')
        self.assertEqual(result['dt_s'], .0004)

    def test_program_end_and_subminimum_interval_stop(self):
        for now in (.24, .25, .02-5e-8):
            result = campaign.decision(metadata(time_s=now), .0004, 36000., PHASES, minimum_dt=1e-7)
            self.assertTrue(result['action'].startswith('STOP'))

    def test_ten_hour_authorization_cannot_be_extended(self):
        args = argparse.Namespace(output=str(campaign.BASE/'uncreated_test_campaign'),
            parent=str(campaign.INITIAL), minimum_dt=1e-7, max_halvings=8)
        for seconds in (36000.001, float('inf'), float('nan')):
            args.wall_seconds = seconds
            with self.assertRaisesRegex(ValueError, 'ten hours'):
                campaign.plan(args)

    def test_automatic_stationary_run_cannot_cross_into_traverse(self):
        result = campaign.decision(metadata(time_s=.22), .0004, 36000., PHASES)
        self.assertTrue(result['action'].startswith('STOP'))


class FailureClassificationTests(unittest.TestCase):
    def classify(self, error=SCREEN_ERROR, screen=None, status='FAIL_GPU_STAGE_NOT_ACCEPTED'):
        return campaign.classify_failure({'status': status, 'error': error}, screen)

    def test_only_a_completed_failed_temporal_screen_can_refine(self):
        screen = failed_time_screen()
        self.assertEqual(self.classify(screen=screen), 'REFINE')
        screen['global_screen']['pass'] = True
        screen['local_screen']['passed'] = False
        self.assertEqual(self.classify(screen=screen), 'REFINE')

    def test_missing_or_inconsistent_screen_never_refines(self):
        screens = [None, {}, failed_time_screen(), failed_time_screen(), failed_time_screen()]
        screens[2]['passed'] = True
        screens[3]['crossing_ids_match'] = False
        screens[4]['global_screen']['pass'] = True
        for screen in screens:
            with self.subTest(screen=screen):
                self.assertNotEqual(self.classify(screen=screen), 'REFINE')

    def test_geometry_energy_cuda_source_and_root_failures_stop(self):
        for error in ('geometry/shell/visibility path failed; separate event resolution required',
                      'mechanical energy balance failed', 'Explicit CUDA solve failed',
                      'source/input drift before commit', 'registered secant iteration budget exhausted',
                      'unrecognized failure', ''):
            with self.subTest(error=error):
                self.assertNotEqual(self.classify(error=error, screen=failed_time_screen()), 'REFINE')

    def test_success_or_unknown_status_cannot_be_a_failed_screen_retry(self):
        for status in ('PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION', 'UNKNOWN', ''):
            self.assertNotEqual(self.classify(screen=failed_time_screen(), status=status), 'REFINE')

    def test_malformed_nested_screen_is_a_stop(self):
        for value in (None, [], 'PASS'):
            screen = failed_time_screen(); screen['global_screen'] = value
            self.assertEqual(self.classify(screen=screen), 'STOP')


class TemporaryWorkspace(unittest.TestCase):
    def setUp(self):
        parent = ROOT/'tmp'; parent.mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix='gpu_campaign_test_', dir=parent)
        self.folder = Path(self.temporary.name).resolve()
        self.assertTrue(self.folder.is_relative_to(parent.resolve()))

    def tearDown(self):
        self.temporary.cleanup()


class AtomicAndPathTests(TemporaryWorkspace):
    def test_path_cannot_escape_workspace(self):
        with self.assertRaises(ValueError):
            campaign.inside(Path('../outside_workspace'))

    def test_incomplete_atomic_write_preserves_previous_checkpoint_pointer(self):
        target = self.folder/'progress.json'
        previous = dict(checkpoint='last_complete_PASS.npz', finished=False)
        campaign.atomic_json(target, previous)
        with self.assertRaises(ValueError):
            campaign.atomic_json(target, dict(checkpoint='failed.npz', time=float('nan')))
        self.assertEqual(campaign.read(target), previous)

    def test_existing_campaign_cannot_be_overwritten_or_launch_child(self):
        folder = self.folder/'already_recorded'; folder.mkdir()
        marker = folder/'progress.json'; marker.write_bytes(b'original evidence')
        protocol = dict(output=campaign.rel(folder), wall_limit_s=36000.)
        with mock.patch.object(campaign, 'singleton_and_awake') as mutex, \
             mock.patch.object(campaign.subprocess, 'Popen') as child:
            with self.assertRaises(FileExistsError):
                campaign.supervise(protocol, {})
        self.assertEqual(marker.read_bytes(), b'original evidence')
        mutex.assert_not_called(); child.assert_not_called()


class ParentReceiptTests(TemporaryWorkspace):
    def receipt(self, name, parent=None):
        folder = self.folder/name; folder.mkdir()
        source = self.folder/'frozen_source.txt'
        if not source.exists():
            source.write_text('frozen synthetic source', encoding='utf-8')
        if parent is None:
            origin = self.folder/'origin.npz'; origin.write_bytes(b'independent origin')
            parent_time, parent_step = .01, 100
        else:
            origin = ROOT/parent[1]['screened_continuation_checkpoint']
            parent_time, parent_step = parent[1]['screened_end_time_s'], parent[1]['screened_end_step']
        paths = {}
        for role in ('coarse', 'fine_1', 'fine_2'):
            directory = folder/role; directory.mkdir()
            paths[role] = directory/'NOT_RESTARTABLE_candidate.npz'
            paths[role].write_bytes((name+role).encode())
        history = [dict(checkpoint=campaign.rel(paths[role]), sha256=campaign.sha(paths[role]),
            time_s=parent_time+half*.0004/2) for half, role in enumerate(('fine_1', 'fine_2'), 1)]
        screen = dict(passed=True, global_screen={'pass': True}, local_screen={'passed': True},
            crossing_ids_match=True, coarse_checkpoint_sha256=campaign.sha(paths['coarse']),
            fine_checkpoint_sha256=[row['sha256'] for row in history])
        campaign.atomic_json(folder/'pair_screen.json', screen)
        (folder/'run.log').write_text('synthetic completed child', encoding='utf-8')
        meta = dict(status=campaign.PASS, accepted=False, restart_allowed=False, diagnostic_only=True,
            accepted_time_advanced_s=0., screened_continuation_authorized=True,
            source_sha256={campaign.rel(source): campaign.sha(source)}, reference_sha256={},
            artifacts_sha256={campaign.rel(path): campaign.sha(path) for path in paths.values()},
            screened_history=history, pair_screen=screen, dt_s=.0004,
            screened_continuation_checkpoint=history[-1]['checkpoint'],
            screened_continuation_sha256=history[-1]['sha256'],
            screened_end_time_s=history[-1]['time_s'], screened_end_step=parent_step+2,
            parent_checkpoint=campaign.rel(origin), parent_sha256=campaign.sha(origin),
            parent_time_s=parent_time, parent_gpu_directory=campaign.rel(parent[0]) if parent else 'origin')
        self.publish(folder, meta)
        return folder, meta

    def publish(self, folder, meta):
        campaign.atomic_json(folder/'manifest.json', meta)
        campaign.atomic_json(folder/'supervision.json', dict(exit_code=0, source_inputs_still_match=True,
            manifest_sha256=campaign.sha(folder/'manifest.json'), log_sha256=campaign.sha(folder/'run.log')))

    def test_latest_complete_parent_and_child_lineage_are_verified(self):
        parent = self.receipt('parent'); child = self.receipt('child', parent)
        actual = campaign.verify_parent(child[0], parent[1]['source_sha256'], parent[0])
        self.assertEqual(actual['screened_continuation_checkpoint'], child[1]['screened_history'][-1]['checkpoint'])

    def test_fine_checkpoint_tampering_is_rejected(self):
        folder, meta = self.receipt('parent')
        (ROOT/meta['screened_history'][0]['checkpoint']).write_bytes(b'changed failed candidate')
        with self.assertRaises(ValueError):
            campaign.verify_parent(folder)

    def test_unrelated_passing_screen_cannot_endorse_fine_candidate(self):
        folder, meta = self.receipt('parent')
        meta['pair_screen']['fine_checkpoint_sha256'][0] = meta['pair_screen']['coarse_checkpoint_sha256']
        campaign.atomic_json(folder/'pair_screen.json', meta['pair_screen']); self.publish(folder, meta)
        with self.assertRaisesRegex(ValueError, 'last screened fine'):
            campaign.verify_parent(folder)

    def test_incomplete_or_failed_receipt_is_not_a_parent(self):
        folder, meta = self.receipt('parent')
        meta['status'] = 'FAIL_GPU_STAGE_NOT_ACCEPTED'; self.publish(folder, meta)
        with self.assertRaisesRegex(ValueError, 'complete GPU receipt'):
            campaign.verify_parent(folder)

    def test_wrong_prior_parent_is_rejected(self):
        first = self.receipt('first'); second = self.receipt('second', first)
        wrong_child = self.receipt('wrong_child', first)
        with self.assertRaisesRegex(ValueError, 'exact prior screened parent'):
            campaign.verify_parent(wrong_child[0], first[1]['source_sha256'], second[0])

    def test_source_drift_is_rejected(self):
        folder, meta = self.receipt('parent')
        (self.folder/'frozen_source.txt').write_text('altered source', encoding='utf-8')
        with self.assertRaises(ValueError):
            campaign.verify_parent(folder)

    def test_inconsistent_half_time_is_rejected(self):
        folder, meta = self.receipt('parent')
        meta['screened_history'][0]['time_s'] += 1e-7; self.publish(folder, meta)
        with self.assertRaisesRegex(ValueError, 'half-step'):
            campaign.verify_parent(folder)


class FakeClock:
    def __init__(self):
        self.value = 0.

    def __call__(self):
        return self.value


class ChildControllerTests(TemporaryWorkspace):
    def scenario(self, *, wall_limit=300., kind='timeout'):
        clock = FakeClock(); kernel = mock.Mock(); waits = []; processes = []; commands = []
        parent = self.folder/'initial_parent'; parent.mkdir()
        checkpoint = parent/'last_complete_PASS.npz'; checkpoint.write_bytes(b'unchanged physical parent')
        before = campaign.sha(checkpoint)
        output = self.folder/'campaign'
        meta = metadata()
        supplied = dict(eligible_crossed=dict(particle_count=0, mass_kg=0.),
            crossed_under_annular_shoulder=dict(particle_count=0, mass_kg=0.))
        meta['observer'] = dict(supplied_material=supplied)
        meta['screened_continuation_sha256'] = before
        snapshot = dict(latest_pass_stage=campaign.rel(parent), physical_time_s=.012,
            checkpoint=campaign.rel(checkpoint), checkpoint_sha256=before, screened_steps=100)
        protocol = dict(output=campaign.rel(output), initial_parent=campaign.rel(parent),
            initial_snapshot=snapshot, requested_macro_dt_s=.0004, wall_limit_s=wall_limit,
            source_sha256={}, child_source_sha256={}, minimum_dt_s=1e-7, max_halvings=8,
            phases=copy.deepcopy(PHASES), python=sys.executable)
        if kind == 'temporal':
            protocol['max_halvings'] = 1

        class FakeProcess:
            pid = 987654321

            def __init__(self, command, **kwargs):
                self.returncode = None
                processes.append(self)
                commands.append(command)
                if kind == 'timeout':
                    clock.value = min(wall_limit, 650.)-.25
                elif kind == 'failure':
                    child = ROOT/command[command.index('--output')+1]
                    child.mkdir()
                    campaign.atomic_json(child/'manifest.json',
                        dict(status='FAIL_GPU_STAGE_NOT_ACCEPTED', error="ValueError('energy failed')"))
                    self.returncode = 1
                elif kind == 'stop':
                    (output/'STOP_REQUESTED').write_text('explicit stop', encoding='utf-8')
                elif kind == 'temporal':
                    child = ROOT/command[command.index('--output')+1]; child.mkdir()
                    campaign.atomic_json(child/'manifest.json',
                        dict(status='FAIL_GPU_STAGE_NOT_ACCEPTED', error=SCREEN_ERROR,
                             source_sha256={}, reference_sha256={}, artifacts_sha256={}, parent_sha256=before))
                    campaign.atomic_json(child/'pair_screen.json', failed_time_screen())
                    for role in ('coarse', 'fine_1', 'fine_2'):
                        (child/role).mkdir()
                        campaign.atomic_json(child/role/'candidate_qualification.json',
                            dict(status='PASS_CANDIDATE_GATES_REQUIRES_PAIR_SCREEN'))
                    self.returncode = 1

            def poll(self):
                return self.returncode

            def wait(self, timeout=None):
                if self.returncode is not None:
                    return self.returncode
                waits.append(timeout); clock.value += timeout
                raise subprocess.TimeoutExpired('synthetic child', timeout)

        def terminate(process, stream):
            process.returncode = -9
            return dict(method='TEST_OWNED_CHILD_ONLY', exit_code=0)

        with mock.patch.object(campaign, 'singleton_and_awake', return_value=(kernel, 123)), \
             mock.patch.object(campaign.time, 'perf_counter', side_effect=clock), \
             mock.patch.object(campaign, 'verify'), \
             mock.patch.object(campaign.subprocess, 'Popen', side_effect=FakeProcess), \
             mock.patch.object(campaign, 'terminate_owned_tree', side_effect=terminate) as termination, \
             contextlib.redirect_stdout(io.StringIO()):
            code = campaign.supervise(protocol, meta)
        progress = campaign.read(output/'progress.json')
        self.assertEqual(progress['checkpoint'], snapshot['checkpoint'])
        self.assertEqual(progress['checkpoint_sha256'], before)
        self.assertEqual(progress['latest_pass_stage'], snapshot['latest_pass_stage'])
        self.assertEqual(progress['completed_stage_count'], 0)
        self.assertEqual(campaign.sha(checkpoint), before)
        self.assertTrue(progress['finished'])
        self.assertIsNone(progress['child_pid'])
        kernel.CloseHandle.assert_called_once_with(123)
        self.assertEqual(len(processes), 2 if kind == 'temporal' else 1)
        self.commands = commands
        return code, progress, clock.value, waits, termination.call_count

    def test_global_deadline_caps_active_wait_and_retains_checkpoint(self):
        code, progress, elapsed, waits, kills = self.scenario(wall_limit=300.)
        self.assertEqual(code, 2)
        self.assertLessEqual(elapsed, 300.)
        self.assertLessEqual(max(waits), .25)
        self.assertEqual(kills, 1)

    def test_child_650_second_deadline_caps_active_wait(self):
        code, progress, elapsed, waits, kills = self.scenario(wall_limit=36000.)
        self.assertEqual(code, 2)
        self.assertLessEqual(elapsed, 650.)
        self.assertLessEqual(max(waits), .25)
        self.assertEqual(kills, 1)

    def test_failed_fake_child_never_becomes_parent(self):
        code, progress, elapsed, waits, kills = self.scenario(kind='failure')
        self.assertEqual(code, 2)
        self.assertEqual(progress['status'], 'STOPPED_CHILD_FAILURE')
        self.assertEqual(kills, 0)

    def test_explicit_stop_terminates_owned_child_and_retains_checkpoint(self):
        code, progress, elapsed, waits, kills = self.scenario(kind='stop')
        self.assertEqual(code, 0)
        self.assertEqual(progress['status'], 'STOPPED_USER_REQUEST')
        self.assertEqual(kills, 1)

    def test_temporal_failure_halves_once_from_exact_same_parent(self):
        code, progress, elapsed, waits, kills = self.scenario(kind='temporal')
        self.assertEqual(code, 2)
        first, second = self.commands
        self.assertEqual(first[first.index('--parent')+1], second[second.index('--parent')+1])
        self.assertEqual(float(second[second.index('--dt')+1]), float(first[first.index('--dt')+1])/2)
        self.assertEqual(progress['halvings_from_same_parent'], 1)
        self.assertEqual(progress['status'], 'STOPPED_CHILD_FAILURE')


if __name__ == '__main__':
    unittest.main()

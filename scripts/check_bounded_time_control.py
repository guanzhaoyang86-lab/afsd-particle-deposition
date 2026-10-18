"""Pure controller transaction/cost checks; performs zero physical evaluations."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import sys
import time
import traceback
import bounded_time_control as control


ROOT = Path(__file__).resolve().parents[1]
STATUS = 'PASS_BOUNDED_TIME_CONTROLLER_LOGIC_NOT_PROCESS_ACCEPTANCE'


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def rejected(function):
    try:
        function()
    except (ValueError, TypeError, KeyError):
        return
    raise AssertionError('invalid transaction was accepted')


def digest(number):
    return hashlib.sha256(str(number).encode()).hexdigest()


def advance(controller, error=.5):
    attempt = controller.propose()
    return controller.resolve(attempt, physical_passed=True, temporal_passed=True,
                              normalized_error=error, fine2_sha256=digest(attempt['attempt']))


def pass_transaction_and_stale_receipts():
    c = control.LocalTimeController(digest(0)); before = c.snapshot()
    attempt = c.propose()
    require(c.parent_sha256 == before['parent_sha256'] and c.time_s == .003,
            'proposal changed diagnostic parent')
    rejected(c.propose)
    tampered = dict(attempt, parent_sha256=digest(99))
    rejected(lambda: c.resolve(tampered, physical_passed=True, temporal_passed=True,
                               normalized_error=.1, fine2_sha256=digest(1)))
    c.resolve(attempt, physical_passed=True, temporal_passed=True,
              normalized_error=.1, fine2_sha256=digest(1))
    require(c.parent_sha256 == digest(1) and abs(c.time_s-.00325) < 1e-15,
            'PASS did not select exactly the fine-two endpoint')
    rejected(lambda: c.resolve(attempt, physical_passed=True, temporal_passed=True,
                               normalized_error=.1, fine2_sha256=digest(2)))
    require(c.snapshot()['accepted'] is False and c.snapshot()['restart_allowed'] is False,
            'diagnostic transaction became accepted')
    return dict(stale_modified_duplicate_receipts_rejected=True, exact_parent_committed=True)


def temporal_rejection_preserves_parent():
    c = control.LocalTimeController(digest(0)); expected = (.00025, .000125, .0000625)
    for dt in expected:
        attempt = c.propose(); require(attempt['coarse_dt_s'] == dt, 'wrong bounded level')
        rejected(lambda: c.resolve(attempt, physical_passed=True, temporal_passed=False,
                                   normalized_error=2., fine2_sha256=digest(1)))
        c.resolve(attempt, physical_passed=True, temporal_passed=False, normalized_error=2.)
        require(c.parent_sha256 == digest(0) and c.time_s == .003 and
                not c.qualified_intervals, 'temporal rejection leaked a candidate into parent')
    rejected(c.propose)
    require(c.status == 'STOP_MINIMUM_TIMESTEP_TEMPORAL_FAILURE_NOT_QUALIFIED',
            'minimum level failure did not stop')
    return dict(levels_us=[250., 125., 62.5], rejected_fine_parents=6, no_unbounded_halving=True)


def physical_failure_is_terminal():
    c = control.LocalTimeController(digest(0)); attempt = c.propose()
    c.resolve(attempt, physical_passed=False, temporal_passed=False)
    rejected(c.propose)
    require(c.level == 0 and c.time_s == .003 and c.parent_sha256 == digest(0),
            'physical failure triggered refinement or advancement')
    return dict(physical_failure_terminal=True, parent_unchanged=True)


def conservative_growth_and_reset():
    c = control.LocalTimeController(digest(0)); attempt = c.propose()
    c.resolve(attempt, physical_passed=True, temporal_passed=False, normalized_error=2.)
    c.resolve(c.propose(), physical_passed=True, temporal_passed=False, normalized_error=2.)
    advance(c, .25); require(c.level == 2, 'grew after only one small-error pair')
    advance(c, .251); require(c.level == 2 and c.small_error_streak == 0,
                              'error above .25 did not reset streak')
    advance(c, .2); require(c.level == 2, 'grew after reset and one pair')
    advance(c, .25); require(c.level == 1, 'two qualifying pairs did not grow by two')
    advance(c, .1); advance(c, .1)
    require(c.level == 0 and abs(c.time_s-.0035) < 1e-15, 'growth exceeded 250 us or target')
    rejected(c.propose)
    return dict(two_consecutive_errors_required=True, maximum_coarse_us=250.)


def snapshot_reconstruction_and_tamper():
    c = control.LocalTimeController(digest(0)); advance(c)
    c.propose(); saved = c.snapshot()
    require(control.LocalTimeController.from_snapshot(saved).snapshot() == saved,
            'pending transaction did not survive exact reconstruction')
    for field, value in (('accepted', True), ('time_s', .004), ('level', 2),
                         ('primary_candidate_reservation', 0)):
        bad = copy.deepcopy(saved); bad[field] = value
        rejected(lambda bad=bad: control.LocalTimeController.from_snapshot(bad))
    bad = copy.deepcopy(saved); bad['events'][0]['next_parent_sha256'] = digest(99)
    rejected(lambda: control.LocalTimeController.from_snapshot(bad))
    return dict(history_replayed=True, altered_counters_parent_and_acceptance_rejected=True)


def attempts_and_endpoint_are_bounded():
    c = control.LocalTimeController(digest(0))
    for unused in range(2):
        c.resolve(c.propose(), physical_passed=True, temporal_passed=False, normalized_error=2.)
    advance(c, .2); advance(c, .2)
    c.resolve(c.propose(), physical_passed=True, temporal_passed=False, normalized_error=2.)
    for unused in range(5):
        advance(c, .5)
    require(abs(c.time_s-.0034375) < 1e-15 and c.attempts_started == 10,
            'attempt budget setup differs')
    rejected(c.propose)
    require(c.snapshot()['primary_candidate_reservation'] == 30,
            'primary candidate cap differs')
    require(control.LocalTimeController.from_snapshot(c.snapshot()).status == c.status,
            'attempt-limit stop did not round trip')
    rejected(c.replay_grid)
    rejected(lambda: control.LocalTimeController(digest(0), target_s=.005))
    return dict(maximum_pair_attempts=10, maximum_primary_candidates=30,
                unfinished_budget_stop_not_qualified=True)


def independent_replay_grid():
    origin = .003000000000000001
    c = control.LocalTimeController(digest(0), t0_s=origin)
    # A rejected 250 us comparison creates a mixed qualified schedule.
    c.resolve(c.propose(), physical_passed=True, temporal_passed=False, normalized_error=2.)
    while c.status == 'READY_DIAGNOSTIC_ONLY':
        advance(c, .2)
    grid = c.replay_grid()
    require(grid[0]['t0_s'] == origin and grid[-1]['t1_s'] == origin+.0005,
            'independent replay lost the full study interval')
    require(all(a['t1_s'] == b['t0_s'] for a, b in zip(grid, grid[1:])),
            'independent replay grid has gaps')
    require(all(row['origin_sha256'] == digest(0) for row in grid),
            'replay inherited a later primary parent')
    for i, interval in enumerate(c.qualified_intervals, 1):
        rows = [r for r in grid if r['primary_interval'] == i]
        require(len(rows) == 4 and all(r['dt_s'] == interval['fine_dt_s']/2 for r in rows),
                'replay did not halve each qualified fine step')
        require([r['shared_primary_fine_endpoint'] for r in rows] == [False, True, False, True],
                'shared comparison endpoints differ')
    require(all(r['accepted'] is False and r['restart_allowed'] is False for r in grid),
            'replay schedule published acceptance')
    return dict(replay_steps=len(grid), all_primary_fine_endpoints_compared=True,
                shared_origin_sha256=digest(0))


def shared_budget_and_failed_cost_floor():
    budget = control.stage_budget(3, 0., 0.)
    require(abs(budget['forecast_s']-220.4652) < 1e-10, 'registered cost floor lost')
    slower = control.stage_budget(3, 60., 100.)
    require(slower['forecast_s'] == 304.5, 'failed/slow measured cost was ignored')
    rejected(lambda: control.stage_budget(3, 200., 0.))
    rejected(lambda: control.stage_budget(3, 60., 900.))
    rejected(lambda: control.stage_budget(5, 0., 0.))
    rejected(lambda: control.stage_budget(3, float('nan'), 0.))
    # Primary time consumes the same ledger later used by replay stages.
    rejected(lambda: control.stage_budget(2, 60., 1000.))
    rejected(lambda: control.stage_budget(4, 0., 0., previous_candidate_count=59))
    replay = control.stage_budget(4, 0., 0., previous_candidate_count=58)
    require(abs(replay['forecast_s']-283.9536) < 1e-10 and replay['total_candidate_limit'] == 62,
            'replay stage or shared candidate bound differs')
    return dict(initial_three_root_forecast_s=budget['forecast_s'], shared_total_wall_limit_s=1200.,
                replay_has_no_fresh_wall_budget=True)


def all_id_error_and_original_limits():
    global_screen = dict(position_max_m=0., temperature_max_K=0., stress_relative=.01,
                         plastic_relative=0., **{'pass': True})
    local = dict(particle_ids=list(range(918)),
        rows=[dict(particle_id=i, position_error_m=0., temperature_error_K=0.,
                   stress_relative=0., plastic_relative=0.) for i in range(918)],
        position_limit_m=5e-6, temperature_limit_K=.5, stress_relative_limit=.1,
        stress_denominator_floor_Pa=1e6, plastic_relative_limit=.2,
        plastic_denominator_floor=1e-5, passed=True)
    local['rows'][791]['plastic_relative'] = .15
    require(abs(control.normalized_screen_error(global_screen, local, .0005)-.75) < 1e-15,
            'local worst particle was hidden by global norm')
    bad = copy.deepcopy(local); bad['rows'].pop()
    rejected(lambda: control.normalized_screen_error(global_screen, bad, .0005))
    bad = copy.deepcopy(local); bad['plastic_relative_limit'] = .3
    rejected(lambda: control.normalized_screen_error(global_screen, bad, .0005))
    bad = copy.deepcopy(local); bad['rows'][791]['plastic_relative'] = float('nan')
    rejected(lambda: control.normalized_screen_error(global_screen, bad, .0005))
    bad = copy.deepcopy(local); bad['rows'][791]['plastic_relative'] = .3
    rejected(lambda: control.normalized_screen_error(global_screen, bad, .0005))
    return dict(particle_count=918, original_limits_unchanged=True,
                maximum_normalized_error=.75, incomplete_nonfinite_and_false_pass_rejected=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='outputs/bounded_time_control_20261005/v1')
    args = parser.parse_args(); out = (ROOT/args.output).resolve()
    require(out.is_relative_to(ROOT/'outputs') and out != ROOT/'outputs' and not out.exists(),
            'new project output directory required')
    out.mkdir(parents=True, exist_ok=False)
    before = time.perf_counter()
    paths = (Path(__file__).resolve(), Path(control.__file__).resolve())
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    report = dict(status='RUNNING_PURE_CONTROLLER_CHECKS', source_sha256=hashes,
        command=[sys.executable, *sys.argv], checks=[], root_evaluations=0,
        physical_evaluations=0, material_evaluations=0, thermal_evaluations=0,
        nnps_evaluations=0, scope='TRANSACTION_GRID_BUDGET_POLICY_ONLY', **control.DIAGNOSTIC)
    checks = (pass_transaction_and_stale_receipts, temporal_rejection_preserves_parent,
        physical_failure_is_terminal, conservative_growth_and_reset,
        snapshot_reconstruction_and_tamper, attempts_and_endpoint_are_bounded,
        independent_replay_grid, shared_budget_and_failed_cost_floor,
        all_id_error_and_original_limits)
    try:
        for check in checks:
            report['checks'].append(dict(name=check.__name__, status='PASS', metrics=check()))
        require(all(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == value
                    for name, value in hashes.items()), 'controller sources changed during checks')
        report.update(status=STATUS, exit_code=0)
    except Exception:
        report.update(status='FAIL_BOUNDED_TIME_CONTROLLER_LOGIC', exit_code=1,
                      error=traceback.format_exc())
    report['wall_seconds'] = time.perf_counter()-before
    (out/'receipt.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=report['status'], checks=len(report['checks']),
                          wall_seconds=report['wall_seconds'])))
    if report['exit_code']:
        print(report['error'], file=sys.stderr)
    return report['exit_code']


if __name__ == '__main__':
    raise SystemExit(main())

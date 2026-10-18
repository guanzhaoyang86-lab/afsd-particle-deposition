"""Pure transaction/budget policy for a bounded, unaccepted local time study.

No particle, force, thermal, material, motion, solver or publication code is
imported. A supervising runner supplies fully audited candidate receipts.
Only a passing coarse-versus-two-halves comparison changes the diagnostic
parent. This is a timestep controller, not a different integration formula.
"""
import copy
import math


COMPARISON_SCOPE = 'SAME_IMMEDIATE_PARENT_LOCAL_TWO_RESOLUTION'
TICK_S = 1. / 64000.  # 15.625 us; includes the mandatory independent replay.
COARSE_TICKS = (16, 8, 4)  # 250, 125, 62.5 us, with exactly two fine steps.
MAX_PAIR_ATTEMPTS = 10
MAX_PRIMARY_CANDIDATES = 30
MAX_TOTAL_CANDIDATES = 62
COST_FLOOR_S = 41.3256
TOTAL_WALL_LIMIT_S = 1200.
DIAGNOSTIC = dict(accepted=False, restart_allowed=False,
                  accepted_time_advanced_s=0., diagnostic_only=True)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _finite(value, name):
    require(not isinstance(value, bool), name+' must be a finite number')
    result = float(value)
    require(math.isfinite(result), name+' must be finite')
    return result


def _tick(value, name):
    value = _finite(value, name)
    result = round(value/TICK_S)
    require(abs(result*TICK_S-value) <= 1e-15, name+' is outside the registered grid')
    return result


def _hash(value):
    require(isinstance(value, str) and len(value) == 64 and
            all(c in '0123456789abcdef' for c in value), 'checkpoint SHA256 required')
    return value


def stage_budget(candidate_count, measured_candidate_cost_s, previous_supervised_wall_s,
                 previous_candidate_count=0):
    """One shared supervised-time ledger covers primary, rejects and replay.

    measured_candidate_cost_s must retain every completed candidate's cost,
    including candidates in failed attempts, and the stage-per-candidate cost.
    The runner must also account for partially completed/failed supervision.
    """
    require(type(candidate_count) is int and 1 <= candidate_count <= 4,
            'a supervised stage contains at most four candidate roots')
    require(type(previous_candidate_count) is int and previous_candidate_count >= 0 and
            previous_candidate_count+candidate_count <= MAX_TOTAL_CANDIDATES,
            'shared primary/replay candidate cap exceeded')
    measured = _finite(measured_candidate_cost_s, 'measured candidate cost')
    elapsed = _finite(previous_supervised_wall_s, 'previous supervised wall time')
    require(measured >= 0. and elapsed >= 0., 'negative cost is invalid')
    cost = max(COST_FLOOR_S, measured)
    forecast = 30.+1.5*candidate_count*(cost+1.)
    require(forecast <= 450., 'stage forecast exceeds 450 seconds')
    require(elapsed+forecast <= TOTAL_WALL_LIMIT_S, 'shared study budget exhausted')
    external = min(600., TOTAL_WALL_LIMIT_S-elapsed)
    internal = min(540., external-5.)
    require(forecast <= internal-15., 'insufficient supervised deadline margin')
    return dict(candidate_count=candidate_count, measured_candidate_cost_s=cost,
                previous_candidate_count=previous_candidate_count,
                total_candidate_limit=MAX_TOTAL_CANDIDATES,
                forecast_s=forecast, previous_supervised_wall_s=elapsed,
                total_wall_limit_s=TOTAL_WALL_LIMIT_S, internal_limit_s=internal,
                external_limit_s=external, **DIAGNOSTIC)


def normalized_screen_error(global_screen, local_screen, dx):
    """Return the maximum original-limit ratio across global and all 918 IDs.

    Crossing-ID equality is a separate exact gate and must be included by the
    caller in temporal_passed; it has no meaningful continuous error ratio.
    """
    dx = _finite(dx, 'particle spacing')
    require(dx > 0., 'positive particle spacing required')
    global_fields = (('position_max_m', .01*dx), ('temperature_max_K', .5),
                     ('stress_relative', .1), ('plastic_relative', .2))
    local_fields = (('position_error_m', 5e-6), ('temperature_error_K', .5),
                    ('stress_relative', .1), ('plastic_relative', .2))
    fixed = dict(position_limit_m=5e-6, temperature_limit_K=.5,
                 stress_relative_limit=.1, stress_denominator_floor_Pa=1e6,
                 plastic_relative_limit=.2, plastic_denominator_floor=1e-5)
    require(all(local_screen[k] == v for k, v in fixed.items()),
            'local temporal thresholds changed')
    ids = local_screen['particle_ids']; rows = local_screen['rows']
    require(len(ids) == len(rows) == len(set(ids)) == 918 and
            [row['particle_id'] for row in rows] == ids,
            'exactly all 918 unique particle IDs must be compared')
    def ratios(row, fields):
        result = [_finite(row[name], name)/limit for name, limit in fields]
        require(all(value >= 0. for value in result), 'negative error is invalid')
        return result
    global_error = max(ratios(global_screen, global_fields))
    local_error = max(value for row in rows for value in ratios(row, local_fields))
    require(global_screen['pass'] == (global_error <= 1.) and
            local_screen['passed'] == (local_error <= 1.),
            'declared temporal pass disagrees with original limits')
    return max(global_error, local_error)


class LocalTimeController:
    """Bounded diagnostic parent selection; checkpoints always stay unaccepted.

    propose() reserves one of at most ten pair attempts. resolve() must
    follow all physical gates and the original global/local/crossing screen.
    After a temporal rejection both fine candidates are discarded as parents;
    at most two smaller levels exist. A physical failure ends the study.
    """
    def __init__(self, parent_sha256, t0_s=.003, target_s=None):
        self.origin_sha256 = _hash(parent_sha256)
        self.parent_sha256 = self.origin_sha256
        self.origin_time_s = _finite(t0_s, 'start time')
        target_s = self.origin_time_s+.0005 if target_s is None else target_s
        self.target_time_s = _finite(target_s, 'target time')
        self.start_tick = _tick(t0_s, 'start time')
        self.target_tick = _tick(target_s, 'target time')
        require(self.start_tick == 192 and self.target_tick == 224 and
                abs(self.target_time_s-self.origin_time_s-.0005) <= 1e-15,
                'only the registered 3 to 3.5 ms study is allowed')
        self.current_tick = self.start_tick
        self.level = 0
        self.small_error_streak = 0
        self.attempts_started = 0
        self.pending = None
        self.events = []
        self.qualified_intervals = []
        self.status = 'READY_DIAGNOSTIC_ONLY'

    @property
    def time_s(self):
        return self.origin_time_s+(self.current_tick-self.start_tick)*TICK_S

    def _time_at(self, tick):
        return self.origin_time_s+(tick-self.start_tick)*TICK_S

    def propose(self):
        require(self.status == 'READY_DIAGNOSTIC_ONLY' and self.pending is None,
                'no new attempt after a pending, failed or completed study')
        if self.attempts_started >= MAX_PAIR_ATTEMPTS:
            self.status = 'STOP_PAIR_ATTEMPT_LIMIT_NOT_QUALIFIED'
            raise ValueError('ten-pair attempt cap reached')
        remaining = self.target_tick-self.current_tick
        while COARSE_TICKS[self.level] > remaining:
            self.level += 1
            require(self.level < len(COARSE_TICKS), 'endpoint cannot fit registered levels')
        ticks = COARSE_TICKS[self.level]
        self.attempts_started += 1
        self.pending = dict(attempt=self.attempts_started,
            parent_sha256=self.parent_sha256, t0_s=self.time_s,
            t1_s=self._time_at(self.current_tick+ticks),
            coarse_dt_s=ticks*TICK_S, fine_dt_s=ticks*TICK_S/2,
            start_tick=self.current_tick, end_tick=self.current_tick+ticks,
            level=self.level, comparison_scope=COMPARISON_SCOPE, **DIAGNOSTIC)
        return copy.deepcopy(self.pending)

    def resolve(self, attempt, *, physical_passed, temporal_passed,
                normalized_error=None, fine2_sha256=None):
        require(self.pending is not None and attempt == self.pending,
                'stale, modified or unrelated attempt receipt')
        require(type(physical_passed) is bool and type(temporal_passed) is bool,
                'explicit physical and temporal boolean gates required')
        error = None if normalized_error is None else _finite(normalized_error, 'normalized error')
        require(error is None or error >= 0., 'negative normalized error')
        if physical_passed and temporal_passed:
            require(error is not None and error <= 1., 'passing temporal error must be <= 1')
            _hash(fine2_sha256)
            require(fine2_sha256 != self.parent_sha256, 'new fine endpoint hash required')
        else:
            require(fine2_sha256 is None, 'failed candidates cannot be proposed as parents')
        require(physical_passed or not temporal_passed, 'physical failure cannot have temporal PASS')
        before_parent = self.parent_sha256
        if not physical_passed:
            outcome = 'STOP_PHYSICAL_FAILURE_NOT_QUALIFIED'
            self.status = outcome
            self.small_error_streak = 0
        elif not temporal_passed:
            self.small_error_streak = 0
            if self.level == len(COARSE_TICKS)-1:
                outcome = 'STOP_MINIMUM_TIMESTEP_TEMPORAL_FAILURE_NOT_QUALIFIED'
                self.status = outcome
            else:
                self.level += 1
                outcome = 'REJECT_BOTH_FINE_CANDIDATES_RETRY_SAME_PARENT'
        else:
            self.parent_sha256 = fine2_sha256
            self.current_tick = attempt['end_tick']
            interval = dict(attempt, fine2_sha256=fine2_sha256,
                            normalized_error=error)
            self.qualified_intervals.append(interval)
            self.small_error_streak = self.small_error_streak+1 if error <= .25 else 0
            if self.small_error_streak == 2:
                self.level = max(0, self.level-1)
                self.small_error_streak = 0
            outcome = 'PASS_LOCAL_PAIR_DIAGNOSTIC_PARENT_ONLY'
            if self.current_tick == self.target_tick:
                self.status = 'PRIMARY_COMPLETE_INDEPENDENT_REPLAY_REQUIRED'
        event = dict(attempt=copy.deepcopy(attempt), physical_passed=physical_passed,
            temporal_passed=temporal_passed, normalized_error=error,
            fine2_sha256=fine2_sha256, previous_parent_sha256=before_parent,
            next_parent_sha256=self.parent_sha256, outcome=outcome,
            next_level=self.level, next_time_s=self.time_s,
            small_error_streak=self.small_error_streak, **DIAGNOSTIC)
        self.events.append(event)
        self.pending = None
        return copy.deepcopy(event)

    def snapshot(self):
        return dict(schema='BOUNDED_LOCAL_TIME_CONTROLLER_V1',
            origin_sha256=self.origin_sha256, parent_sha256=self.parent_sha256,
            t0_s=self.origin_time_s, target_s=self.target_time_s,
            time_s=self.time_s, level=self.level,
            small_error_streak=self.small_error_streak,
            attempts_started=self.attempts_started,
            primary_candidate_reservation=3*self.attempts_started,
            pending=copy.deepcopy(self.pending), events=copy.deepcopy(self.events),
            qualified_intervals=copy.deepcopy(self.qualified_intervals),
            status=self.status, **DIAGNOSTIC)

    @classmethod
    def from_snapshot(cls, saved):
        """Reconstruct decisions from their ledger, never trust saved counters."""
        result = cls(saved['origin_sha256'], saved['t0_s'], saved['target_s'])
        for event in saved['events']:
            attempt = result.propose()
            require(attempt == event['attempt'], 'saved attempt ancestry or grid changed')
            actual = result.resolve(attempt, physical_passed=event['physical_passed'],
                temporal_passed=event['temporal_passed'],
                normalized_error=event['normalized_error'], fine2_sha256=event['fine2_sha256'])
            require(actual == event, 'saved controller decision changed')
        if saved['pending'] is not None:
            require(result.propose() == saved['pending'], 'saved pending attempt changed')
        elif saved['status'] == 'STOP_PAIR_ATTEMPT_LIMIT_NOT_QUALIFIED':
            try:
                result.propose()
            except ValueError:
                pass
        require(result.snapshot() == saved, 'saved controller metadata or counters changed')
        return result

    def replay_grid(self):
        """Halve every qualified fine step, independently from the 3 ms origin.

        This only returns a grid. It does not qualify the replay, bound its
        affordability or allow any checkpoint to be accepted/restarted.
        """
        require(self.status == 'PRIMARY_COMPLETE_INDEPENDENT_REPLAY_REQUIRED',
                'replay grid requires the complete primary diagnostic path')
        result = []
        for index, interval in enumerate(self.qualified_intervals, 1):
            ticks = (interval['end_tick']-interval['start_tick'])//4
            for substep in range(4):
                begin = interval['start_tick']+substep*ticks
                result.append(dict(t0_s=self._time_at(begin), t1_s=self._time_at(begin+ticks),
                    dt_s=ticks*TICK_S, primary_interval=index,
                    primary_fine_step=substep//2+1,
                    shared_primary_fine_endpoint=(substep % 2 == 1),
                    origin_sha256=self.origin_sha256, **DIAGNOSTIC))
        require(len(result) <= 32 and result[0]['t0_s'] == self.origin_time_s and
                abs(result[-1]['t1_s']-self.target_time_s) <= 1e-15,
                'replay grid exceeds the registered interval bound')
        return result

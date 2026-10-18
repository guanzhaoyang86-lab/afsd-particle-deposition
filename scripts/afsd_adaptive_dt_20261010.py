"""Pure administrative recovery of the original screened macro timestep.

No material, force, thermal, geometry, particle or GPU operator is imported.
The driver must independently verify all original candidate gates and receipts.
Every proposed larger interval still needs one coarse and two fine candidates.
The policy never turns a failed candidate into a continuation parent.
"""
import copy
import math


PROTOCOL = 'ORIGINAL_LIMITS_TWO_PASS_DT_RECOVERY_20261010_V1'
DEFAULT_MAXIMUM_DT_S = .000474609375
TEMPORAL_ERROR = "ValueError('original coarse/fine screen failed; no GPU continuation')"
LOCAL_LIMITS = dict(position_limit_m=5e-6, temperature_limit_K=.5,
    stress_relative_limit=.1, stress_denominator_floor_Pa=1e6,
    plastic_relative_limit=.2, plastic_denominator_floor=1e-5)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value, name):
    require(type(value) in (int, float) and math.isfinite(value),
        name+' must be a finite number, not a boolean')
    return float(value)


def digest(value):
    require(isinstance(value, str) and len(value) == 64
        and all(c in '0123456789abcdef' for c in value), 'checkpoint SHA256 required')
    return value


def normalized_screen(screen, dx, particle_ids):
    """Normalize every original global/local limit; exact ID gates stay exact.

    particle_ids must come from the independently verified complete parent,
    not from the screen itself. Local thresholds and denominator floors are
    checked against the frozen comparison implementation, not user input.
    """
    dx = finite(dx, 'particle spacing')
    require(dx > 0., 'positive particle spacing required')
    require(isinstance(particle_ids, (list, tuple)) and particle_ids
        and all(type(v) is int for v in particle_ids)
        and len(set(particle_ids)) == len(particle_ids), 'complete unique parent IDs required')
    expected = sorted(particle_ids)
    require(isinstance(screen, dict), 'complete pair screen required')
    glob = screen['global_screen']; local = screen['local_screen']
    require(isinstance(glob, dict) and isinstance(local, dict), 'global/local screens required')
    require(all(type(local[k]) in (int, float) and local[k] == value
        for k, value in LOCAL_LIMITS.items()), 'original local thresholds or floors changed')
    ids = local['particle_ids']; rows = local['rows']
    require(isinstance(ids, list) and isinstance(rows, list)
        and all(type(v) is int for v in ids)
        and ids == expected and len(rows) == len(expected)
        and [row['particle_id'] for row in rows] == expected
        and all(type(row['particle_id']) is int for row in rows),
        'local screen must cover every complete-parent particle ID exactly once')
    global_fields = dict(position_max_m=.01*dx, temperature_max_K=.5,
        stress_relative=.1, plastic_relative=.2)
    local_fields = dict(position_error_m=LOCAL_LIMITS['position_limit_m'],
        temperature_error_K=.5, stress_relative=.1, plastic_relative=.2)

    def ratios(row, fields):
        answer = {name: finite(row[name], name)/limit for name, limit in fields.items()}
        require(all(math.isfinite(value) and value >= 0. for value in answer.values()),
            'negative or nonfinite normalized comparison error')
        return answer

    global_ratios = ratios(glob, global_fields)
    local_ratios = {name: -1. for name in local_fields}
    local_worst_ids = {}
    for row in rows:
        for name, value in ratios(row, local_fields).items():
            if value > local_ratios[name]:
                local_ratios[name] = value
                local_worst_ids[name] = row['particle_id']
    global_ratio = max(global_ratios.values()); local_ratio = max(local_ratios.values())
    require(type(glob['pass']) is bool and glob['pass'] == (global_ratio <= 1.)
        and type(local['passed']) is bool and local['passed'] == (local_ratio <= 1.),
        'declared comparison pass disagrees with original limits')
    require(type(screen['crossing_ids_match']) is bool and type(screen['passed']) is bool
        and screen['passed'] == (glob['pass'] and local['passed'] and screen['crossing_ids_match']),
        'pair pass/crossing flags are inconsistent')
    return dict(ratio=max(global_ratio, local_ratio), global_ratio=global_ratio,
        local_ratio=local_ratio, global_ratios=global_ratios,
        local_ratios=local_ratios, local_worst_particle_ids=local_worst_ids,
        compared_particle_count=len(expected), crossing_ids_match=screen['crossing_ids_match'],
        passed=screen['passed'])


def validate_state(state):
    require(isinstance(state, dict) and state.get('protocol') == PROTOCOL,
        'registered adaptive policy state required')
    config = state['config']
    minimum = finite(config['minimum_dt_s'], 'minimum dt')
    maximum = finite(config['maximum_dt_s'], 'maximum dt')
    nominal = finite(state['nominal_dt_s'], 'nominal dt')
    require(0. < minimum <= 1e-7 and minimum <= nominal <= maximum <= .005,
        'registered timestep bounds required')
    require(1. < finite(config['growth_factor'], 'growth factor') <= 2.
        and 0. < finite(config['small_error_ratio'], 'small error ratio') <= .25,
        'growth must be at most two and require original-limit ratio at most 0.25')
    for key, lower, upper in (('required_small_passes', 2, 100), ('cooldown_passes', 1, 32),
            ('max_halvings', 0, 8)):
        require(type(config[key]) is int and lower <= config[key] <= upper,
            'invalid adaptive policy '+key)
    for key in ('small_error_streak', 'cooldown_remaining', 'growth_failures',
            'consecutive_halvings', 'attempts', 'passes'):
        require(type(state[key]) is int and state[key] >= 0, 'invalid adaptive state '+key)
    require(state['passes'] <= state['attempts'] and state['consecutive_halvings'] <= config['max_halvings']
        and state['cooldown_remaining'] <= 32 and type(state['growth_pending']) is bool,
        'inconsistent adaptive counters')
    digest(state['parent_sha256'])
    require(state['status'] in ('READY', 'RUNNING', 'STOPPED'), 'invalid adaptive status')
    pending = state['pending']
    require((state['status'] == 'RUNNING') == (pending is not None), 'inconsistent pending attempt')
    if pending is not None:
        actual = finite(pending['dt_s'], 'actual dt')
        require(pending['parent_sha256'] == state['parent_sha256']
            and pending['nominal_dt_s'] == nominal and minimum <= actual
            and actual <= nominal*(1.+1e-14)
            and type(pending['phase_clipped']) is bool and type(pending['growth_trial']) is bool,
            'pending attempt differs from registered parent or timestep')
        require(pending['phase_clipped'] or actual == nominal,
            'only an explicit phase boundary may shorten a proposal')
    return state


def initial_state(parent_sha256, nominal_dt_s, *, maximum_dt_s=DEFAULT_MAXIMUM_DT_S,
        minimum_dt_s=1e-7, required_small_passes=2, small_error_ratio=.25,
        growth_factor=2., cooldown_passes=4, max_halvings=8):
    """Fresh state: old PASS history does not silently grant growth eligibility."""
    result = dict(protocol=PROTOCOL, status='READY', parent_sha256=digest(parent_sha256),
        nominal_dt_s=finite(nominal_dt_s, 'nominal dt'), last_passed_dt_s=None,
        small_error_streak=0, cooldown_remaining=0, growth_failures=0,
        growth_pending=False, consecutive_halvings=0, attempts=0, passes=0, pending=None,
        config=dict(maximum_dt_s=maximum_dt_s, minimum_dt_s=minimum_dt_s,
            required_small_passes=required_small_passes, small_error_ratio=small_error_ratio,
            growth_factor=growth_factor, cooldown_passes=cooldown_passes, max_halvings=max_halvings))
    validate_state(result)
    return result


def begin_attempt(state, actual_dt_s, *, phase_clipped=False):
    """Bind the driver's phase/cost-approved actual interval to this exact parent.

    The driver computes the exact two-half phase endpoint before calling this.
    Passing a clipped endpoint does not replace the next nominal timestep.
    """
    validate_state(state)
    require(state['status'] == 'READY', 'policy must be ready without a pending attempt')
    require(type(phase_clipped) is bool, 'explicit phase clipping flag required')
    result = copy.deepcopy(state)
    result['pending'] = dict(parent_sha256=state['parent_sha256'],
        nominal_dt_s=state['nominal_dt_s'], dt_s=finite(actual_dt_s, 'actual dt'),
        phase_clipped=phase_clipped, growth_trial=bool(state['growth_pending'] and not phase_clipped))
    result['status'] = 'RUNNING'; result['attempts'] += 1
    validate_state(result)
    return result


def _event(state, result, action, reason, normalized=None):
    return dict(action=action, reason=reason, attempt=state['attempts'],
        parent_before_sha256=state['parent_sha256'], parent_after_sha256=result['parent_sha256'],
        attempted_dt_s=state['pending']['dt_s'], nominal_before_s=state['nominal_dt_s'],
        next_nominal_dt_s=result['nominal_dt_s'], phase_clipped=state['pending']['phase_clipped'],
        growth_trial=state['pending']['growth_trial'], small_error_streak=result['small_error_streak'],
        cooldown_remaining=result['cooldown_remaining'], normalized_screen=normalized,
        diagnostic_only=True, accepted=False, restart_allowed=False, accepted_time_advanced_s=0.)


def resolve_pass(state, screen, dx, particle_ids, fine2_sha256, *, original_candidate_gates_passed):
    """Advance only a fully verified fine-two endpoint, then propose grow/hold."""
    validate_state(state)
    require(state['status'] == 'RUNNING', 'a matching pending attempt is required')
    require(original_candidate_gates_passed is True, 'all original candidate gates must pass')
    normalized = normalized_screen(screen, dx, particle_ids)
    require(normalized['passed'] and normalized['crossing_ids_match'], 'complete original pair PASS required')
    digest(fine2_sha256)
    require(fine2_sha256 != state['parent_sha256']
        and len(screen['fine_checkpoint_sha256']) == 2
        and screen['fine_checkpoint_sha256'][-1] == fine2_sha256,
        'must advance exact new fine-two checkpoint')
    result = copy.deepcopy(state); config = state['config']; pending = state['pending']
    result.update(status='READY', parent_sha256=fine2_sha256, pending=None,
        last_passed_dt_s=pending['dt_s'], consecutive_halvings=0, passes=state['passes']+1)
    action = 'HOLD'
    if pending['phase_clipped']:
        result['small_error_streak'] = 0
        reason = 'PHASE_ENDPOINT_PRESERVES_NOMINAL_DT'
    else:
        result['growth_pending'] = False
        if pending['growth_trial']:
            result['growth_failures'] = 0
        if state['cooldown_remaining']:
            result['cooldown_remaining'] -= 1; result['small_error_streak'] = 0
            reason = 'GROWTH_RETRY_COOLDOWN'
        elif normalized['ratio'] <= config['small_error_ratio']:
            result['small_error_streak'] += 1
            if state['nominal_dt_s'] >= config['maximum_dt_s']:
                result['small_error_streak'] = 0; reason = 'REGISTERED_MAXIMUM_DT'
            elif result['small_error_streak'] >= config['required_small_passes']:
                result['nominal_dt_s'] = min(config['maximum_dt_s'],
                    state['nominal_dt_s']*config['growth_factor'])
                result['small_error_streak'] = 0; result['growth_pending'] = True
                action = 'GROW'; reason = 'CONSECUTIVE_SMALL_COMPLETE_PAIR_ERRORS'
            else:
                reason = 'ACCUMULATING_SMALL_COMPLETE_PAIR_ERRORS'
        else:
            result['small_error_streak'] = 0; reason = 'ERROR_ABOVE_GROWTH_MARGIN'
    validate_state(result)
    return result, _event(state, result, action, reason, normalized)


def resolve_failure(state, manifest, screen, dx, particle_ids, *, original_candidate_gates_passed):
    """Only verified temporal disagreement may halve; every other failure stops.

    A temporal rejection leaves the parent hash unchanged. Failed growth incurs
    four full-step cooldown passes by default; repeated failed growth doubles
    that cooldown up to 32, preventing alternating expensive large/small trials.
    """
    validate_state(state)
    require(state['status'] == 'RUNNING', 'a matching pending attempt is required')
    result = copy.deepcopy(state); result.update(status='STOPPED', pending=None)
    result['small_error_streak'] = 0
    reason = 'NON_TEMPORAL_FAILURE_REQUIRES_DIAGNOSIS'; normalized = None
    temporal = (original_candidate_gates_passed is True and isinstance(manifest, dict)
        and manifest.get('status') == 'FAIL_GPU_STAGE_NOT_ACCEPTED'
        and manifest.get('error') == TEMPORAL_ERROR and isinstance(screen, dict)
        and screen.get('passed') is False and screen.get('crossing_ids_match') is True)
    if temporal:
        try:
            normalized = normalized_screen(screen, dx, particle_ids)
            temporal = normalized['ratio'] > 1. and not normalized['passed']
        except (ValueError, KeyError, TypeError, OverflowError):
            temporal = False; reason = 'INVALID_TEMPORAL_SCREEN_REQUIRES_DIAGNOSIS'
    action = 'STOP'
    if temporal:
        config = state['config']; pending = state['pending']; smaller = pending['dt_s']/2.
        if state['consecutive_halvings'] >= config['max_halvings']:
            reason = 'TEMPORAL_REFINEMENT_LIMIT'
        elif smaller < config['minimum_dt_s']:
            reason = 'MINIMUM_DT_REACHED'
        else:
            if pending['growth_trial']:
                result['growth_failures'] += 1
            cooldown = min(32, config['cooldown_passes'] * (2**min(5, max(0, result['growth_failures']-1))))
            result.update(status='READY', nominal_dt_s=smaller,
                consecutive_halvings=state['consecutive_halvings']+1,
                cooldown_remaining=max(state['cooldown_remaining'], cooldown), growth_pending=False)
            action = 'REFINE'; reason = 'TEMPORAL_SCREEN_FAIL_RETRY_UNCHANGED_PARENT'
    validate_state(result)
    return result, _event(state, result, action, reason, normalized)

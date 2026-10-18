"""Bounded O(h)-curvature refinement around the frozen original response.

Candidate-only wrapper: no physical inputs/settings or final derivative bounds
are changed. Original failures and all attempt evidence remain in the result.
"""
import time

import numpy as np

from . import event_response as original

MAX_HALVINGS = 2
GAP_CONTRACTION_MAX = .65
RETRY_REASON = 'ONE_SIDED_DIFFERENCES_INCONSISTENT'


def _signed_gap(column):
    return np.asarray(column['forward_half_h_s_kg'])-np.asarray(column['backward_half_h_s_kg'])


def _eligibility(attempt):
    """Prove O(h) contraction only on components that actually failed."""
    certificate = dict(eligible=False, reason=attempt.get('reason'))
    if attempt.get('available') or attempt.get('reason') != RETRY_REASON:
        return certificate
    columns = attempt.get('column_evidence', [])
    if not columns or not all(c['refinement_pass'] for c in columns):
        certificate['reason'] = 'CENTRAL_REFINEMENT_NOT_PASSED'
        return certificate
    for probe in attempt['probes']:
        if (probe.get('status') != 'PROBE_COMPLETED' or
                not probe.get('raw_heat_certificate_pass') or probe.get('physical_branch_change')):
            certificate['reason'] = 'PHYSICAL_PROBE_CERTIFICATE_NOT_PASSED'
            return certificate
    column = columns[-1]; axis = int(column['axis'])
    h = float(np.asarray(attempt['force_step_N'])[axis])
    candidates = [p for p in attempt['probes'] if p.get('axis') == axis]
    plus = next((p for p in candidates if p['step_N'] == h), None)
    minus = next((p for p in candidates if p['step_N'] == -h), None)
    if plus is None or minus is None or h <= 0:
        certificate['reason'] = 'COARSE_ONE_SIDED_EVIDENCE_MISSING'
        return certificate
    keys = ('mean_velocity_m_s', 'velocity_m_s')
    coarse = np.stack([(np.asarray(plus[key])-2*np.asarray(attempt['baseline'][key])+
                        np.asarray(minus[key]))/h for key in keys])
    fine = _signed_gap(column)
    offending = np.abs(fine) > np.asarray(column['one_sided_allowed_s_kg'])
    contraction = bool(np.any(offending) and np.all(
        np.abs(fine[offending]) <= GAP_CONTRACTION_MAX*np.abs(coarse[offending]))
        and np.all(fine[offending]*coarse[offending] >= 0))
    certificate.update(axis=axis, offending_components=offending,
        coarse_signed_gap_s_kg=coarse, fine_signed_gap_s_kg=fine,
        contraction_upper=GAP_CONTRACTION_MAX, eligible=contraction,
        reason='SMOOTH_O_H_SIDE_GAP' if contraction else 'SIDE_GAP_NOT_O_H_CONTRACTING')
    return certificate


def _cross_attempt(previous_certificate, current):
    axis = previous_certificate['axis']
    column = next((c for c in current.get('column_evidence', []) if c['axis'] == axis), None)
    if column is None:
        return dict(passed=False, axis=axis, reason='PREVIOUS_FAILING_AXIS_EVIDENCE_MISSING')
    previous = np.asarray(previous_certificate['fine_signed_gap_s_kg'])
    following = _signed_gap(column)
    mask = np.asarray(previous_certificate['offending_components'], bool)
    passed = bool(np.all(np.abs(following[mask]) <= GAP_CONTRACTION_MAX*np.abs(previous[mask]))
                  and np.all(following[mask]*previous[mask] >= 0))
    return dict(passed=passed, axis=axis, previous_offending_components=mask,
        previous_signed_gap_s_kg=previous, following_signed_gap_s_kg=following,
        contraction_upper=GAP_CONTRACTION_MAX,
        reason='PREVIOUS_AXIS_GAP_CONTRACTED' if passed else 'PREVIOUS_AXIS_GAP_DID_NOT_CONTRACT')


def _bounded(evaluate, initial_relative):
    if not np.isfinite(initial_relative) or initial_relative <= 0:
        raise ValueError('positive finite initial relative force step required')
    started = time.perf_counter(); attempts = []; certificates = []
    pending = None

    def finish(status, force_unavailable=False):
        answer = dict(attempts[-1])
        answer.update(accepted=False, attempts=attempts, refinement_certificates=certificates,
            refinement_status=status, refinement_halvings=len(attempts)-1,
            maximum_refinement_halvings=MAX_HALVINGS,
            original_response_available=bool(attempts[0]['available']),
            kick_calls=sum(int(a['kick_calls']) for a in attempts),
            kick_seconds=sum(float(a['kick_seconds']) for a in attempts),
            response_seconds=time.perf_counter()-started)
        if force_unavailable:
            answer.update(available=False, B=None, C=None, reason=status)
        return answer

    for index in range(MAX_HALVINGS+1):
        attempt = evaluate(initial_relative/2**index)
        attempts.append(attempt)
        if not attempt['available'] and attempt.get('reason') != RETRY_REASON:
            return finish('STOPPED_ON_NONELIGIBLE_ORIGINAL_FAILURE')
        if pending is not None:
            cross = _cross_attempt(pending, attempt)
            certificates.append(dict(kind='CROSS_ATTEMPT_PREVIOUS_AXIS', attempt=index, **cross))
            if not cross['passed']:
                return finish('REFINEMENT_SMOOTHNESS_CERTIFICATE_FAILED', True)
        if attempt['available']:
            return finish('ORIGINAL_RESPONSE_AVAILABLE' if index == 0 else 'BOUNDED_REFINED_RESPONSE_AVAILABLE')
        eligible = _eligibility(attempt)
        certificates.append(dict(kind='WITHIN_ATTEMPT_H_AND_HALF_H', attempt=index, **eligible))
        if not eligible['eligible']:
            return finish('REFINEMENT_SMOOTHNESS_CERTIFICATE_FAILED', True)
        if index == MAX_HALVINGS:
            return finish('REFINEMENT_BUDGET_EXHAUSTED', True)
        pending = eligible
    raise AssertionError('bounded refinement loop exceeded fixed budget')


def response(*args, **options):
    """Original primitive API, with at most two certified force-step halvings."""
    settings = dict(options)
    initial = settings.pop('force_step_relative', 1e-3)
    return _bounded(lambda relative: original.response(*args, force_step_relative=relative, **settings), initial)


def response_from_trace(trace, dt=None, mu=None, **options):
    """Original journal API; sums cost and retains all failed attempts."""
    settings = dict(options)
    initial = settings.pop('force_step_relative', 1e-3)
    return _bounded(lambda relative: original.response_from_trace(trace, dt=dt, mu=mu,
                    force_step_relative=relative, **settings), initial)

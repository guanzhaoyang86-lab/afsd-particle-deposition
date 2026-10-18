# Bounded candidate differentiation refinement, registered before its tests

2026-10-03. One numerical repair, implemented only in new
`event_response_refined.py`; original response source, physical kick, contracts,
test failures and cloud failure remain unchanged. No physical state advance,
dt/material/contact change, residual tolerance relaxation or extrapolated B/C.

Evidence motivating this change: PID1's actual old kick remains strictly sliding
(minimum sampled slip >0.0038 m/s), while forward/backward finite differences
have their expected opposite O(h) curvature biases. Across a preregistered six
level force-step ladder, the Cyx side-gap decreases by factors near two and
central differences by factors near four. The original side-gap failure at
h/2 is 1.13702e-4 versus 6.04937e-5 s/kg; the h/h2 central difference is only
2.57e-7 s/kg. Initial cloud failure is preserved, not converted to an acceptance.

## New wrapper policy

Public `response` and `response_from_trace` retain original calling conventions.
The first call uses the original force-step setting (default relative 1e-3).
At most TWO further original-response calls halve only the force perturbation.
All actual kick settings and input physical dt remain fixed. Available B/C are
the final original h/2 central estimates, never Richardson-extrapolated values.

A retry is allowed only for `ONE_SIDED_DIFFERENCES_INCONSISTENT`, after the
original central-refinement check passed and every recorded probe kept its
physical branch and raw signed-Q certificate. Reconstruct the coarse forward
and backward differences from the original probes. For each offending fine
component, require same-signed coarse/fine side-gap and
`abs(fine_gap) <= 0.65*abs(coarse_gap)`. A true absolute-value cusp has a constant
side-gap and must fail this certificate. The 0.65 upper bound is registered
before tests, allowing modest higher-order effects around the smooth 0.5 ratio.

After a retry, locate the PREVIOUS failing force axis by identity in the new
column evidence (never compare different first-failure axes). Its previously
offending components must contract by the same factor and sign criterion.
If the next failure is on another axis, its own coarse/fine certificate is
required before another retry. An unrelated failure, failed smoothness
certificate or exhausted two-halving budget returns unavailable with B/C None.
There is no free-response fallback and no extra caller-controlled retry budget.

The original final componentwise derivative bounds remain exactly
`2e-3 max(abs(Dh),abs(Dh/2)) + 2e-5 dt/m`, including one-sided checks.
All original calls/results are retained under `attempts`; total `kick_calls`
and `kick_seconds` include every attempt. The wrapper records elapsed time,
halving count and both within-attempt and cross-attempt contraction certificates.
Every result remains `accepted=false` and candidate-only.

## Preregistered tests and limits

- Literal manufactured smooth functions vbar(F)=F/2+3F^2,
  endpoint(F)=F+4F^2 at F=0, derivative references B=.5,C=1: exactly two
  halvings under the unchanged registered bounds; error <1e-12.
- Literal abs(F) cusp: coarse/fine side-gap is constant; reject without retry.
  An unrelated physical-branch/kick failure must reject without retry.
- Real saved PID1: original stays unavailable; wrapper must be available
  within the two-halving bound and preserve source/journal hashes, v0, inputs,
  physical dt, raw Q and original force-step diagnostics in every attempt.
- Literal reversal B=19/54,C=7/9 remains within 2e-7 s/kg, using no retry.
- Switching first failing force axes must be compared by identity. A budget
  exhaustion remains explicit unavailable. Tests run under external 60 s cap.

This repair establishes only a numerical candidate derivative. The original
continuous-ODE derivative benchmark remains FAIL and process accuracy is not
validated. No actual 342-particle cloud or global root is authorized by this
component test alone; parent driver owns those separate bounded operations.

## Initial result: real full response remains blocked

Initial seven-test suite: six PASS, one FAIL, 1.276 s, exit 1; complete evidence
`logs/process_event_20261003/event_response_refined_initial.log`. Literal smooth
curvature, abs cusp, unrelated-failure rejection, changed first failing axis,
fixed-budget rejection and unchanged reversal pass. The real PID1 availability
requirement fails: `REFINEMENT_BUDGET_EXHAUSTED`, accepted=false, B=C=None.

The first failed Fx column becomes consistent after one halving, but the next
Fy column then fails. Its Cxy side-gap is 1.06909238e-4 s/kg at relative force
step 0.0005; at step 0.00025 it falls to 5.34546071e-5, still above the unchanged
4.89052338e-5 bound. Both contractions meet the registered smoothness check,
but the two-halving budget cannot produce a complete available response.
All three attempts (5+9+9=23 actual kicks) are retained. No third halving,
threshold change, partial derivative fallback, full cloud or global root is run.
The expected-availability regression remains FAIL; no assertion is rewritten.

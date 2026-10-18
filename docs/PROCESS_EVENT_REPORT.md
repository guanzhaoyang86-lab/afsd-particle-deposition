# Event-response candidate root report — 2026-10-03

Scope: numerical proposal only; PySPH remains pinned at
1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5. No mass/time/modulus scaling,
physical clipping, label forces, imposed free-material trajectories or GPU work.
Project remains ACTIVE/INCOMPLETE. Original physical trials and failures stay.

## Implemented versus validated

`event_response.py` replays the original contact kick and computes separate
B=d(actual event mean velocity)/dF and C=d(endpoint velocity)/dF. The independent
reversal references are B=19/54,C=7/9, not B=C/2. The new nonsymmetric
`stage_event_preconditioner.py` couples these responses with frozen corrected
elastic stiffness and continuity. It changes root proposals only; physical
trial, materials, contact, heat, actual constraints and raw acceptance stay.
`run_process_event.py` separates accepted physical parent from REJECTED numerical
derivative anchor and has separate cloud, cost and solve modes with hard limits.

Initial independent suites: response10 PASS/1 FAIL, discrete derivative4 PASS,
Schur15 PASS, runner6 PASS. Original103 predecessor tests remain PASS.
The combined initial preflight correctly FAILED on checker-source drift during
the run. After all authors stopped editing, stable preflight repeated139 tests:
138 PASS/1 retained accuracy FAIL, no errors/skips. Exact logs/report are in
`logs/process_event_20261003/preflight_initial` and `preflight_stable`.

The retained FAIL is
`ActualEventResponse.test_two_free_axis_cross_coupling_independent_ode`:
B error4.06763264e-5 s/kg >2e-5. Original test, physical kick, tolerances and
failure log remain. Independent discrete-midpoint reference agrees with B/C
to4.57e-10/2.10e-9 s/kg; this permits only a numerical-map derivative probe.
This gate distinction was documented AFTER the failure, BEFORE actual probes;
it is not an all-tests-PASS or process-accuracy certificate.

## Actual first cloud failure and bounded numerical repair

`outputs/process_event_20261003/response_probe`: exit2, outer1.1028801 s,
particle1/feature1:0,5 local kicks, zero full physical evaluations/accepted steps.
Reason ONE_SIDED_DIFFERENCES_INCONSISTENT. Cyx side gap1.137017704e-4 s/kg
exceeds6.049372979e-5, while central refinement difference2.56954e-7 passes.
All5 probes have no zero-slip event, endpoint SLIDING,122 accepted subintervals
and5 rejected subintervals. Parent t=850.393632910904 ns remains unchanged.
This is a differentiation failure, not new negative friction or AFSD failure.

The fixed six-level PID1 diagnostic ran ALL levels rather than stopping at green.
13 old kicks, kick wall.510892 s/total.8422274 s. Cyx side gaps are
2.27405e-4,1.13702e-4,5.68508e-5,2.84254e-5,1.42127e-5,7.10641e-6 s/kg;
successive ratios~2 while central errors decrease~4. Minimum sampled slip
>=.0038149 m/s, same physical branch. Evidence supports ordinary smooth O(h)
one-sided truncation, not a visible cusp. Initial source/journal hashes preserved.

New `event_response_refined.py` is a bounded wrapper, not an edit to the failed
implementation: at most2 force-step halvings, only certified contracting
ONE_SIDED failures, same final derivative bounds, all attempts saved. Cusp,
physical branch/heat failure or exhausted budget stays unavailable. The named
refined driver retains the original physical trial and500ns/32-evaluation limits.
See PROCESS_EVENT_REFINEMENT_CONTRACT and PROCESS_EVENT_REFINED_ROOT_CONTRACT.

## Independent temporal-accuracy diagnosis

`event_schur_smooth_state_diagnosis.log` checks the smooth two-face fixture against
independent RK4(4096/8192) and fixed-midpoint equations. It ran exit0 in.995 s.
The old kick DOES return the two-half-step solution: independent agreement
1.57e-16, not the coarse step. Local endpoint estimates1.70806e-7/1.71974e-7
match independently measured local errors and are below2.23608e-7 m/s.
But complete-interval endpoint error3.37611e-7 and mean error4.10210e-6 m/s are
1.51x/18.35x this LOCAL bound. These ratios are a scope diagnosis, not a newly
invented global tolerance. The current controller does not control the actual
mean trajectory integral or its derivative. Continuous accuracy FAIL remains.

## Final bounded-run outcome

Refined component initial:7 tests6 PASS/1 FAIL,exit1,1.276 s. Real PID1 uses23
kicks and exhausts the two-halving budget. At hrel=.00025 the Cxy side gap is
5.345460712391156e-5 versus4.890523379248915e-5 s/kg. The earlier Fx issue contracts,
but the Fy column still fails. Returned B/C remain None. No third halving,
relaxed bound or test rewrite followed; initial RED assertion remains.

Final unified refined preflight:146 executed,144 PASS/2 FAIL,no errors;exit2,
internal7.9012598 s. The two failures are the original continuous reference and
real PID1 refinement budget. The final wiring suite was skipped by fail-fast;
it then ran separately1/1 PASS,external.3783896 s. Distinct total147 tests:
145 PASS/2 FAIL. Separate wiring PASS does NOT override the failed launch gate.
The new refined cloud, cost probe, fullroot and time/space sensitivity are all
NOT_RUN. Source-independent status/hash readback is saved separately.

Effective state remains old stage004 at850.393632910904 ns, real shoulder
crossings0. New physical evaluations0, accepted steps0. No incomplete response
or REJECTED candidate is a valid physical state or restart point. This bounded
numerical work package stops here; it will not consume another turn tuning h.

## What this does not solve

A derivative receipt or converged500ns root is not process accuracy, actual
transport or deposition. Highest elastic omega*dt~10 remains under-resolved.
The next ONE work package is PROCESS_EVENT_NEXT_WORK_PACKAGE: retain inertia,
establish coupled time accuracy and separately auditable algorithmic energy loss,
with a same-end-time benchmark and cost intersection. No accumulation of
arbitrary microsteps as a surrogate for completion. Finite open-tool transitions,
accepted inlet transactions/IDs/source balances and actual transport remain
subsequent required integration tasks. All old material/free-energy/parameter
risks remain; no production or complete paper-reproduction label is granted.

# Candidate event response: preregistered component contract

2026-10-03. Scope: new numerical candidate response only. Existing physical kick,
cone law, material, geometry, time step, acceptance and frozen evidence remain
unchanged. Every return has `accepted=false`; response availability is not state
acceptance. No AFSD advance or actual 342-particle cloud is run by this component.

`response_from_trace(trace, dt=None, mu=...)` consumes the existing
`coupled_boundary.evaluate_boundary` journal. `response(...)` is the primitive.
Outputs `B=d(vbar_actual)/d(F_internal)` and `C=d(v1_actual)/d(F_internal)` are
detached 3 by 3 arrays in s/kg. Geometry, normals, normal forces, wall velocity,
v0, mass, dt, free-axis mask and prescribed velocity stay fixed. vbar is the
interval-duration-weighted mean of the actual event trajectory. It is never
replaced by the endpoint midpoint. No tangent is installed into the physical
kick. Constrained force columns are identically zero by the physical projector.

## Preregistered numerical choices, before initial tests

- Zero contacts: exact `B=dt/(2m) P`, `C=dt/m P`, zero kick calls. With contacts,
  replay the anchor kick once, then central force differences at h and h/2 for
  every free force axis: at most 13 kick calls. Old kick tolerances are unchanged.
- Force scale (N): max of Euclidean internal force, total friction capacity,
  mass times the maximum of initial velocity relative to the first wall and
  inter-wall relative speeds divided by dt, and smallest normal
  floating-point positive number. Relative speed uses the first wall frame, so
  a common velocity boost does not set the differentiation scale. h=1e-3 scale
  N by default, separately recorded for every axis; no implicit unit-sized floor.
- Both derivative estimates must satisfy componentwise
  `abs(D_h-D_h/2) <= 2e-3*max(abs(D_h),abs(D_h/2)) + 2e-5*dt/m`.
  Return the h/2 estimate, without extrapolation, pseudo-inverse or fallback.
  Also require agreement of forward/backward h/2 one-sided estimates under the
  same bound; central refinement alone can conceal an absolute-value cusp.
- Detect changed physical stop-event contact sequence and endpoint cone status.
  Normalize the two old zero-slip event labels; coalesce adjacent events within
  1e-9 dt. Numerical interval counts and rejection counts are recorded but never
  define a physical branch. Endpoint zero-slip threshold is 512 eps times the
  old kick relative speed scale. Zero-slip cone saturation uses 2e-8 capacity
  tolerance. Sliding direction must have positive dot product with baseline.
  Any event-sequence, zero/sliding, interior/boundary change, nonfinite result,
  raw heat certificate failure, or h refinement failure returns unavailable.
- Replay uses the same inputs and old kick settings. Saved-journal versus
  replay differences for endpoint, actual mean and raw Q are recorded. Endpoint
  and mean bounds are `2e-12 + 1e-9*relative_speed_scale` m/s; Q bound is
  `1e-24 + 1e-9*sum(abs(saved face Q))` J. These only test reproducibility,
  never authorize a physical step. Initial source/journal arrays remain intact.
- Each probe records force, h, endpoint, actual mean, event sequence, endpoint
  cone status, interval count, signed raw-Q sum/min, original rounding bounds,
  reconstructed `-J_t dot u_mid` maximum discrepancy, and SHA256 of the entire
  unmodified kick result. Raw Q is never absolute-valued, clipped or repaired.
  Reconstruction agreement is tested per interval/face against the original
  rounding bound plus `1024 eps (abs(raw Q)+norm(Jt)*norm(u_mid))`; this is a
  coordinate-conversion certificate, separate from the unchanged signed raw-Q
  rejection. Tiny Q cannot alone set a cancellation-error tolerance.
  Failure contexts, including old cone/interval details, are retained. Total
  kick calls and elapsed kick/response times are reported, including failures.

## Independent acceptance references

Literal 1D reversal: m=1, u0=1, F=-2, muFn=1, dt=1, only x free.
Reference vbar=-1/18, endpoint=-2/3, Q=7/18, Bxx=19/54, Cxx=7/9.
The reference constants do not invoke the tested response or kick. Physical
values tolerance 2e-11; response derivatives absolute 2e-7 s/kg.

Free response exact within 1e-14; static stick tangent within 1e-10 s/kg;
common wall/material boost response within 1e-8 s/kg. Two-face reversal with
only y free has m=1, u0y=1, Fy=-4, two capacities 1 each, dt=1:
vbar_y=-11/18, v1y=-5/3, Q=14/9, Byy=11/27, Cyy=8/9.
Response derivatives tolerance 3e-7 s/kg. Off-diagonal contact mobility is
checked explicitly, and contact permutation must leave B,C within 1e-8 s/kg.
A two-free-axis, two-face case with a prescribed z velocity exercises nonzero
cross-axis response; compare against independent smooth ODE sensitivity
integration at absolute 2e-5 s/kg, without invoking the tested kick/response in
the reference integrator. Contact order and common boosts retain that response.
An endpoint stick/slip threshold fixture must return unavailable explicitly.

Every test command has an external 60-second timeout. Initial evidence is kept;
at most two evidence-driven repair rounds, each with separately named output.
Full process, production and any physical acceptance remain NOT_RUN here.

Pre-test audit corrected a transcribed two-face Byy constant before any test
execution: differentiation gives 1/72+25/72+5/108=11/27. No observed failure or
tolerance was changed. Prescribed output derivative rows are analytically zero;
the raw finite-difference rows are retained and must be within 128 eps dt/m
before the identity is imposed on the returned candidate matrices.

## Initial evidence and explicit reference-scope distinction

Initial `test_event_response.py`: 11 tests, 10 PASS and 1 FAIL, 0.311 s.
`logs/process_event_20261003/event_response_initial.log` preserves all raw probe
metrics and the failure. The two-free-axis fixture's B differs from the
independent continuous-ODE sensitivity by 4.06763264e-5 s/kg, above the original
2e-5 bound. This remains an accuracy FAIL; dt, old kick tolerances, that test and
its bound remain unchanged. No physical/process accuracy claim follows.

The old numerical kick uses four 0.05-s midpoint intervals for this dt=0.2-s
fixture. Candidate B is intentionally the derivative of that actual discrete
map. An independent implicit midpoint variational solve gives a diagnostic
discrepancy of 4.57e-10 for B and 2.10e-9 for C. This explains but does not erase
the continuous-ODE accuracy FAIL.

New, separate `test_event_response_discrete.py` tests that discrete reference
with absolute derivative bound 1e-7 s/kg, original four-step partition checked
explicitly, independent implicit residual <=1e-14 m/s, both off-diagonal entries
nonzero, prescribed rows and columns exact zero, h/h2 and one-sided consistency,
contact permutation plus common boost, replay and input purity. This is a
separate candidate-derivative benchmark, registered before its initial suite,
not a replacement reference or a repair of the original failed test. No source
repair was made in response to the original continuum accuracy failure.

Separate discrete suite initial result: 4/4 PASS, 0.227 s, exit 0. Its full
evidence is `logs/process_event_20261003/event_response_discrete_initial.log`.
The original continuous benchmark remains 10 PASS / 1 FAIL, exit 1. Commands
are each `.envs/pysph/python.exe -B tests/review/<test_file>.py`, launched through
`subprocess.run(..., timeout=60)` and captured to distinct fresh log names.
`kick_seconds` includes the probe's summary and complete-result hash overhead;
it is a conservative kick-plus-evidence wall-time cost. `response_seconds`
records total component elapsed time. No physical/source repair round used.

# Independent slow-startup reference and first branch boundary

2026-10-03. Registered physical contract: `PROCESS_SLOW_CONTRACT.md`.
Classification `LOCAL_SLOW_LOAD_FULL_STARTUP_FEASIBILITY_NOT_AFSD`.
This oracle does not advance the actual AFSD checkpoint and does not implement
the contact branch after a zero-slip event.

## Independence and event meaning

`scripts/reference_hht_slow.py` imports only literal `hht_slow_inputs.inputs`
and the previous independent reference's all-pair spatial algebra, packing and
ledger indices. It imports no tested backend, HHT interpolation, material
return, contact evaluator or BE heat routine. It independently evaluates the
new normal 5 Hz programme, constant +z wall speed, full-vector Coulomb sliding,
Jaumann elastic ODE, density continuity and thermal ODE. Original nonzero
stress/velocity and all physical masses/moduli are retained. Accepted states
have no prescribed material velocity and no temperature/yield/heat clipping.

Both references use fixed classical RK4 grids, 7.8125 ns and 3.90625 ns.
Ordinary shared snapshots are exactly 62.5 ns apart. Classical RK4's cubic
dense extension uses coefficients

    c0=y0
    c1=h*k1
    c2=h*(-1.5*k1+k2+k3-.5*k4)
    c3=(2*h/3)*(k1-k2-k3+k4).

The oracle checks endpoint and intermediate RK stages, zero-slip roots/minima
of the full tangential dense polynomial, and normal-gap extrema. The analytic
slow wall is represented in the gap check by its cubic Taylor polynomial plus
an explicit fourth-derivative remainder bound. Nonlinear geometry, thermal
domain and FKS elastic scope are additionally checked at dense fractions .25,
.5 and .75; this is NOT a whole-path nonlinear FKS interval certificate.

A signed-z crossing with nonzero transverse slip is logged, never called a
stop. A full dense-vector root must meet 512*eps*max(1,relevant speed).
The observed first event is on the natural symmetry plane, but that symmetry
is measured, not prescribed. The old .05 m/s fast-fixture margin is absent.
RK candidates crossing an event are retained as REJECTED with state, time,
event margin and identity; they are not accepted post-event reference states.

The localization budget remains the original 13 candidate attempts, at most
12 actual halvings. An accepted positive subinterval may narrow the approach
interval. If positive candidates reach a previously rejected right endpoint,
the purported event is not established: return
`EVENT_LOCALIZATION_INCONSISTENT`, not a fabricated/inverted event bracket.
An approach interval is deliberately reported `EVENT_APPROACH_BLOCKED`;
there is no exact event state or post-event integration permission.

## Execution and retained repair

Initial independent suite: 8 PASS, .047 s test / .236310 s external.
Initial reference: 13.152015 s internal / 13.347775 s external, exit 0,
`PASS_PREFIX_REFERENCE_EVENT_BLOCKED`. Initial source and output remain intact.
An audit found `event_halvings=12` mislabeled the final candidate-attempt index;
it was not the actual number of halvings. Repair 1 separates both counts and
adds the inconsistent-bracket guard with a synthetic regression. It changes
no physical formula, input, timestep, error bound or evaluation budget.

Repair-1 suite: **9 PASS**, .059 s test / .301526 s external, exit 0.
Repair-1 reference: **13.125873 s internal / 13.324709 s external**, exit 0.
Each reference used 13 event candidate attempts and 9 actual halvings.
Accepted and rejected physical NPZ hashes match the initial run exactly;
only corrected metadata and source receipts differ.

Authoritative result: `logs/process_slow_20261003/reference_repair1/report.json`.
Source SHA256:
`de98a19d05a4c796969a9fa16b159b88771f374dea1be771959de7c667ce3d66`.
Configuration SHA256:
`82e024ab392788610c2c2131adea653905110eb43ed7a27c8133cbeaab4f93be`.

## What was actually established

Both grids identify the first full-vector zero approach at left face 0,
particle ID 4, initially y=.005 m. Coarse interval is
[5.974655151367281, 5.974670410156344] us; fine interval is
[5.974655151367073, 5.974662780761604] us. They overlap. The dense-root
transverse speeds are 1.02e-17 and 1.21e-18 m/s, well below the registered
1.14e-13 m/s roundoff bound. This is not a mislabeled z-only turn.
The opposite-face ID 22 has nearly identical approach slip; do not interpret
the selected minimum's ID as proof that its symmetry partner is a separate
macroscopic event.

The last accepted pre-event fine state is at 5.974655151367073 us. Its minimum
full slip remains positive, 1.26233e-7 m/s; minimum normal overlap is 1.61882 um
and minimum yield margin is 313.363 MPa. No normal opening or plastic yield
was observed before this event. Neither reference reached 20 us.

Only the **96 shared accepted snapshots through 5.9375 us** constitute the
same-time reference accuracy comparison. All nine oracle limits pass:

|Quantity|Maximum shared-time difference|Registered bound|
|---|---:|---:|
|weighted RMS position|5.34e-17 m|1e-11 m|
|weighted RMS velocity|9.72e-11 m/s|2e-6 m/s|
|weighted RMS full stress|.0011474 Pa|70 Pa|
|temperature|5.97e-12 K|1e-9 K|
|face force component|2.99e-9 N|.007 N|
|face impulse component|2.13e-15 N s|1.4e-9 N s|
|face wall work|3.46e-16 J|2e-9 J|
|face raw friction heat|1.65e-16 J|2e-9 J|
|face reservoir heat|4.78e-19 J|2e-10 J|

This verifies the semidiscrete independent prefix, not continuum spatial
accuracy, a full 20 us HHT result, AFSD transport, or production economics.
The desired 1 ms reference remains **NOT_RUN**, because its forecast exceeds
the reference's 60 s external cap. Positive plastic accuracy also remains
NOT_RUN. No post-event solution or stationary contact is supplied.

## Read-only indication of the required next branch

`reference_pre_event_force_diagnostic.log` evaluates the unchanged final
accepted state, without advancing it. At IDs 4 and 22, internal +z force is
4.67484334 N while Coulomb capacity is 2.73923931 N. To hold zero relative
acceleration against this instantaneous internal force would require 1.706621
times the available friction. This SUPPORTS continued reverse sliding rather
than feasible sticking. It is not an exact event-state solve or authorization
to prescribe a velocity. A future event method must locate the event, test
actual coupled feasibility and restart with the correct force/acceleration
history while retaining physical state and separate numerical work.

## Receipt interface

`run_reference(data,step_s,deadline)` returns every accepted state/time,
shared-grid sample state/time, rejected candidates, component-turn diagnostics,
terminal metadata and last accepted physical record. CLI report fields used
by the driver include `source_sha256`, `configuration_sha256`,
`artifact_sha256`, `prefix_comparison`, `certified_common_end_time_s`,
`event_status`, `coarse_terminal`, `fine_terminal` and the explicit status.

`reference_history.npz` contains identical shared `time_s`, fine `states`,
`coarse_states`, persistent IDs, literal masses and reference positions.
`coarse/fine_history.npz` retain every accepted ordinary and approach state.
Separate rejected JSON/NPZ pairs retain attempted candidates with
`accepted=false,restart_allowed=false`. A passed prefix receipt does not
authorize publication of a tested solver state beyond its certified prefix.

## Read-only S0 first-rejection audit: predictor versus converged traction

The actual S0 driver accepted eight steps through 2 us and rejected its first
500 ns interval [2,2.5] us at evaluation 1, PID 3/face 0. Saved evidence is
`outputs/process_slow_20261003/S0_initial/trajectory/REJECTED.json`.
Its raw Q=5.58823e-8 J is positive; whole interpolated path is passive, no zero
event was detected and minimum slip=.0460024 m/s. The sampled traction-direction
defect .0187148921 exceeds the unchanged .01 bound. This is a different event
from the independent reference's later genuine zero-slip approach at 5.97466 us.

Actual source order in `hht_coupled.advance.evaluate` is: candidate Newmark
positions/velocities -> current midpoint NNPS/material -> current endpoint
internal/boundary forces -> `_contact_work` -> BE heat -> mapped acceleration
and residual/merit. The initial candidate copies old acceleration/density/T.
Therefore this rejection is conclusively a **predictor-stage veto before any
nonlinear residual correction**, not evidence that a converged candidate fails.
It also does not prove that the final coupled root could pass the same guard.

`logs/process_slow_20261003/audit_reference_traction.py` independently audits
already-saved coarse/fine reference states on precisely [2,2.5] us, with no time
integration, candidate publication or source changes. It compares endpoint-
linear traction to the actual intermediate reference Coulomb direction, using
the same linearly interpolated capacity for the direction-only diagnostic.
It separately compares to actual physical force/capacity; those are different
quantities because the actual normal force is not temporally linear.

Results (`reference_traction_audit.json`, .086117 s internal):

- Fine 129 snapshots: maximum direction-only defect .00938707914, below .01.
  Coarse 65 snapshots: .00938490453. PID 3 attains the same fine value within
  roundoff; its midpoint defect is .00938490453.
- Fine maximum physical-force-relative defect is .00693857795.
- Minimum actual slip is .0489022178 m/s. All sampled interpolated friction
  powers remain negative (largest -.0956066567 W).

These values SUPPORT a premature predictor veto as a cause, but do not prove a
valid HHT root. The reference path is not the HHT compatible quadratic velocity
path, and the tested solver's old endpoint contains integration error. Nor is
the finite saved-grid passivity sample a whole-path certificate.

One suitable next bounded diagnostic is the same saved stage-8/500 ns problem
in a private process: temporarily defer ONLY the .01 direction-accuracy veto
during nonlinear candidate exploration; preserve all passivity/zero-event,
thermal, material, geometry and finite-state guards; then reapply the original
.01 accuracy audit to the converged candidate. Retain the result solely as
unpublished diagnostic data, whether it passes or fails, and preserve the
original rejected run. This isolates gate ordering without relaxing the
acceptance criterion, editing frozen production code or advancing physics.

### Follow-up diagnostic executed: gate deferral alone is insufficient

The main thread registered `PROCESS_SLOW_PREDICTOR_DIAGNOSTIC.md` and ran exactly
that private 500 ns diagnostic once. Its authoritative report is
`logs/process_slow_20261003/predictor_diagnostic_initial/report.json` (inspected
after execution). Result: **FINAL_ORIGINAL_GUARD_FAIL**. Twelve evaluations
converged the unchanged raw coupled equations, with velocity residual
6.17552e-10 m/s, but the final maximum direction defect is
.01252632750448 > .01 at PID 3/face 0. The original final guard rejected it.
Raw friction Q remains +5.78006e-8 J, minimum slip .0492045 m/s, and no interior
zero event is indicated. Runtime was .193194 s. The original source/parent
hashes stayed unchanged and the private NPZ is explicitly non-restartable.

Thus an early predictor veto is CONFIRMED, but **not the sole blocking cause**.
The stronger hypothesis that simply delaying the veto would make this interval
admissible is RULED_OUT by the measured final root. The saved independent
reference's .009387 direction defect does not overrule the tested HHT path's
.012526 defect: they are different trajectories/interpolations.

Do not implement gate deferral as the advertised fix or relax .01. The next
bounded dependency is a contact-time representation capable of representing
curved nonzero slip with consistent impulse and raw work, independently tested
before insertion in the actual joint map. Changing only reported friction heat
would not suffice. Later located full-zero reversal is a separate dependency;
the 5.97466 us reference approach remains useful evidence but does not explain
away this earlier nonzero-turning failure. No additional run or physical state
was published by this follow-up.

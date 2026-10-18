# Independent review of the proposed event-response root

2026-10-03. This review reads the old source and saved failed journal; it does
not call an AFSD trial, contact kick, material update, acceptance or commit.
The new preflight entry point is separately implemented but has not been run
by this reviewer. Proposed component and runner tests remain NOT_RUN here.

## Scope and evidence

The proposed coupling is a consistent approximate Newton correction if both
event-average and endpoint force responses are retained separately. Its use
changes numerical proposals only. It does not establish that the actual
342-particle anchor admits a sufficiently smooth force derivative, that the
new root converges, or that 500 ns resolves mechanical dynamics.

Read sources: AGENTS.md, HANDOFF.md, ENGINEERING_APPROX_CONTRACT,
PROCESS_GLOBAL_REPORT, PROCESS_GLOBAL_NONLINEAR_DIAGNOSIS,
PROCESS_GLOBAL_STEP_AUDIT, PROCESS_COMPLETION_CONTRACT; actual
contact_kick.py, coupled_boundary.py, process_thermal_trial.py,
stage_iteration.py, stage_structural_root.py, structural_mobility.py and
both old stage preconditioners. No external method is being imported or
claimed as an independently validated AFSD solver.

One read-only Python extraction of the saved JSON/gzip completed with exit 0
(tool wall approximately 1.1 s), without writing files:

- PROCESS_GLOBAL_FREEZE contains 3565 entries.
- The rejected artifact SHA256 is
  `6b25260780a9fffe76db235259006278ed36b8c074fcf0f7cfb4cc2f82cd6ed2`.
- Its classification is REJECTED, accepted=false, restart_allowed=false.
- Its evaluation starts at 8.50393632910904e-7 s, dt=5e-7 s and binds the
  accepted particle hash
  `8053dbda1773fb057115b322c92ed1d789029459207ee5f545c0d64b31feeeab`.
- The saved boundary has 342 particle traces, 14008 contact subintervals and
  40 zero-slip events. These counts are saved-candidate diagnostics only.

The previous 32-evaluation failure remains a failure, with no new accepted
state. The recorded 99.5357% residual concentration near stop events supports
testing this missing derivative; it does not prove a unique cause.

## Independent algebra and units

Let d=map(candidate)-candidate, h=dt, E=trace(D), V=m/rho at the explicitly
chosen frozen geometry, and rho_ref the old reference density. K is the full
unprestressed frozen elastic stiffness, including its volumetric part. The
old density/pressure sign follows from sigma=-p I and the adjoint internal
force; compression gives the positive density-force term below.

    H = 0.5 E^T diag(V bulk/rho_ref)
    delta_rho = d_rho - h diag(rho_ref) E delta_bar
    delta_F = H d_rho - (h/2) K delta_bar
    B = d(vbar_actual)/d(F_internal)
    C = d(v1_actual)/d(F_internal)
    A = I + (h/2) B K
    A delta_bar = d_bar + B H d_rho
    delta_x1 = d_x1 + h B delta_F
    delta_v1 = d_v1 + C delta_F
    delta_T1 = d_T1

The bulk part of K already results from the displayed density elimination.
Adding a second bulk feedback to delta_F would double count it. B and C are
block diagonal over particles for a frozen local kick, but each 3 by 3 block
can couple force axes through simultaneous contacts.

Units: B,C are s/kg; E is 1/m; K is N/m; H d_rho and delta_F are N;
A is dimensionless; its right-hand side and delta_bar are m/s. The position
correction is m and the density correction kg/m3. A residual from this solve
is a velocity residual, not an impulse or physical reaction.

Without contact B=h P/(2m), C=h P/m. Multiplying the unconstrained Schur
equation by M recovers M+h^2 K/4 and the old pressure right-hand side
(h/4) E^T diag(V bulk/rho_ref) d_rho. This is a useful reduction check.
For actual events C need not equal 2B. Therefore replacing the endpoint
correction by twice the mean correction changes the approximate Jacobian.

The exact kinematic map is x0+h*vbar_actual, so the formula retains d_x1
and also satisfies delta_x1=d_x1+h*(delta_bar-d_bar), up to the linear
solve error. It does not impose vbar=(v0+v1)/2. Keeping the direct defects is
essential: zero proposals cannot serve as proof of a root.

For prescribed axes, the actual projector gives zero B/C rows and force
columns. Mean corrections on prescribed axes equal their original d_bar,
which is generally nonzero during iteration. Solve
A_ff delta_bar_f = rhs_f-A_fc delta_bar_c. Dropping constrained columns of K
before this shift loses the effect of the prescribed-axis correction on
neighboring free axes. The old control still derives axes from original z0
and the old midpoint time; no new stick labels become permanent constraints.
Do not assume A is symmetric positive definite: B need not commute with K,
and the event derivative is not the original symmetric inertial response.
Use a bounded general solve with residual/conditioning diagnostics, no pinv.

## Local derivative references and branch limits

The literal single-face reversal is independently consistent:
m=1, u0=1, F=-2, capacity=1, h=1; tau=1/3, tau'=1/9;
vbar=-1/18, v1=-2/3, B=19/54, C=7/9 (all in the contract's SI units).
These constants must be literal references, not obtained by calling the
kick, the difference routine, or the Schur code under test.

An independent check caught a reference arithmetic error in the draft
two-face contract and reported it to the component author before this review
was written. With only y free, m=1, u0y=1, Fy=-4, two capacities of 1 and h=1,
tau=1/(2-F)=1/6, tau'=1/36, vbar=-11/18 and v1=-5/3. Consequently

    B = tau'/2 + (1-tau)^2/2 - (F+2)*(1-tau)*tau'
      = 1/72 + 25/72 + 5/108 = 11/27
    C = (1-tau) - (F+2)*tau' = 8/9.

The draft 91/216 for B is incorrect; 11/27=88/216 is the independent
reference. Fixing a reference before testing is not relaxation of tolerance;
if any earlier execution used it, retain that initial result and explanation.

The old kick recomputes stopping times, stick feasibility, curved-slip
subdivisions and cone solves. Derivatives must rerun these actual events with
the same old tolerances and max_intervals=4096; freezing event times misses
terms already exposed by 19/54. No-contact exact derivatives are affordable
and avoid subtractive noise, but do not justify an analytic fallback for a
contact whose derivative is unavailable.

Near a stop tie, event order can change under an arbitrarily small force
perturbation. The old selection checks a newly smaller tau before appending
a within-tolerance tie; event labels alone are not a differentiability proof.
Compare physical contact/stop groups, endpoint slip/stick and cone status;
record numerical subdivisions separately. Coalescing near ties is a declared
diagnostic convention, not proof that the original map is differentiable.
Central differences at hF and hF/2 can agree at a kink (the elementary |F|
counterexample). Existing plus/minus probes permit one-sided derivative
comparison without extra kick calls. If sides or physical branches disagree,
report unavailable or explicitly characterized generalized derivatives;
never silently average them into a claimed smooth tangent. This applies to
B and C independently, with finite values and raw event heat certificates.
The component author confirmed the two-face reference correction preceded
all tests. Read-back of the new response source shows forward/backward hF/2
checks for both response matrices and explicit one-sided failure evidence.

## Rejected numerical anchor is not a restart

Using the last rejected candidate as a frozen numerical linearization point
is consistent with an approximate Newton/chord method and the existing
same-history candidate solver. It is permissible only with explicit separate
provenance. It does not change the physical initial-value problem.

Bind the anchor artifact, old run manifest/configuration/source receipts,
dt, particle IDs, accepted_state_hash, old time, and its candidate/evaluation
pair. Derivative inputs may use the saved candidate's internal force and
frozen midpoint faces, normal force and wall velocity. Every old v0, mass,
original z0, stress, ep, temperature history and actual trial initialization
must still come from the accepted stage_004. Check those correspondences;
never deserialize the rejected endpoint as an accepted ParticleArray.

If elastic K/E/V uses accepted geometry while B/C uses the rejected midpoint,
record both choices. That mixed frozen approximate Jacobian is allowable,
but is not the exact derivative at the rejected candidate. A new run still
starts its real trial history at stage_004; it neither resumes evaluation33
of the old attempt nor adds the old failed increment to physical time.
Store anchor metadata as accepted=false, restart_allowed=false and
physical_map=false. Candidate logs, acceptance history and the final physical
ledger must never adopt anchor impulses, heat, stress or plastic history.

## Criteria for exactly one new 500 ns attempt

1. Preserve the 3565-file freeze, all old failures and the published stage_004.
   Pass independent B/C references, shared contact/axis/constraint tests,
   branch rejection, all five-field Schur back-substitutions, nonsymmetric
   linear solve checks, false-root guards, and the original nine suites.
   Bind complete source hashes and the three new suites in a fresh receipt.
2. Construct the actual cloud response in a separately bounded uncommitted
   probe. Every required contact derivative must be available under the
   preregistered branch/refinement rules. Report particle IDs, B/C norms,
   force scales/steps, failures and total calls/intervals/time. One base plus
   central differences at two scales is up to 13 old kick calls per fully
   free contacting particle; the 14008-interval saved candidate warns that
   total construction cost cannot be estimated from only one easy particle.
3. Measure complete new construction, physical trial, proposal, journal I/O,
   replay/validation and cleanup costs. Bind costs to the same source,
   parameters, parent, rejected anchor, dt and settings. Use the registered
   conservative forecast including construction and 32 full evaluations;
   launch only if 1.5*(build+32*evaluation_and_journal)+30 <=450 s, or a
   stricter already-registered launch gate. Inner540/outer570 s remains below
   the per-real-run 600 s ceiling. No larger contact budget or omitted FD
   cost may make the forecast pass. A probe timeout means NOT_RUN full root.
4. Keep dt=500 ns, maximum_steps=1, maximum_evaluations=32, original initial
   accepted-history candidate and the old Anderson/raw-growth policy.
   Freeze B/C once if that is what the cost probe measured; no unbudgeted
   derivative rebuild inside the root. A derivative exception, physical
   guard, cone failure, raw heat failure, NaN, cost overrun or exhausted
   evaluations ends the attempt with a fresh rejected journal. No automatic
   change of dt, tolerance, perturbation size or iteration budget and rerun.
5. Acceptance uses the actual reevaluated map: x<=1e-12 m, v1/vbar<=1e-9 m/s,
   rho-relative<=1e-10 and T<=5e-12 K, then the unchanged complete physical
   acceptance and transactional commit. Record all energy/normal/angular
   diagnostics with their old NOT_ENERGY_CERTIFIED limitations. If rejected,
   new accepted_steps=0 and old stage_004 remains the only continuation.

These are launch/interpretation criteria, not a claim that the probe or root
has run. The root driver owns the actual preregistration and receipts.

## New source integration read-back

Read-back of StageEventPreconditioner confirms the displayed signs and h/2
factor, separate B/C back-substitutions, unchanged direct defects,
prescribed-axis shift, general free-block solve, and full/free backward-error
checks. No algebraic blocker was found. The response producer owns checking
its contact branch and zero prescribed-force columns; the Schur consumer
rejects nonzero prescribed output rows and does not silently modify inputs.

The first run_process_event.py draft correctly creates every real trial from
the old accepted parent, using the rejected journal only in response mode.
The root predictor is the old accepted-history predictor. An unavailable
response aborts, a complete response receipt binds the cached B/C arrays,
and cost/solve modes check the same source, parent, settings and anchor.
The cloud construction is separately bounded at inner180/outer210 s; its
cached arrays are loaded rather than rebuilt in the full root. Record this
separate cost explicitly rather than pretending it is free or multiplying
it by every root evaluation.

One concrete first-draft regression was sent to the driver author: restore
the old pre-save source/input hash check before emitting accepted.npz, in
addition to the check before publishing accepted_history. Otherwise a
detected source change can leave a nominal accepted checkpoint before the
attempt reports rejection. Explicit anchor time/restart checks and trace
v0/mass/actuator correspondence were also recommended to make the numerical
anchor boundary directly auditable. Later fixes/tests belong to the driver's
separate evidence; this paragraph records the reviewed first draft.

New scripts/check_process_event.py verifies the 3565-file old freeze before
and after testing, permits only the three named handoff/index changes, hashes
all backend modules plus the twelve registered suites and driver's SOURCES,
and runs each suite with an external60 s limit. The old nine-suite total must
remain103. Missing registered files yields explicit
NOT_RUN_PREFLIGHT_INCOMPLETE_SOURCES before any suite starts. Passing status
is PASS_EVENT_PREFLIGHT_NOT_PROCESS_VALIDATION. It has not been executed by
this reviewer, per the integration coordinator's instruction.

## Connection to the process completion goal

Even a PASS would establish one larger-step numerical root, not process
transport. The existing frozen-elastic estimate omega_max*dt approximately10
already rules out calling this a resolved-wave time step without further
accuracy evidence. Preserve inertia and give any later high-frequency
algorithmic dissipation a named integrator and independent energy/accuracy
account; do not hide it as material heat or clear checkpoint velocities.

The next process integration must also handle finite open-bore/shoulder-edge
events with consistent work, an accepted-history-bound inlet/ID/source
transaction, real through-hole material and subsequent plastic transport,
and rerun/restart equivalence. The 50 um displacement domain and unintegrated
edge/inlet components remain practical limits. Repeated 500 ns or smaller
accepted windows do not discharge those limits. Keep zero observed crossing
and incomplete process status until actual accepted transport is measured.

## Post-failure review and qualified numerical gate

This section is added AFTER the original response suite's initial execution,
and before any actual full-cloud response/root experiment. It revises the
preflight interpretation above; it is not described as a pre-test criterion.

The original event_response suite produced 10 PASS and 1 FAIL among 11 tests.
The retained failure is
`ActualEventResponse.test_two_free_axis_cross_coupling_independent_ode`:
B[1,1] differed from the independent continuous ODE reference by
4.06763264e-5 s/kg, above its registered 2e-5 bound. The original method,
fixture, physical kick, dt and accuracy bound are retained. Its B assertion
fails before its later C/permutation assertions execute, so those later
checks cannot be counted as passed in that original method.

The separate test_event_response_discrete.py differentiates four explicit
implicit-midpoint equations independently. For each fixed h=.05 interval,
with A=d(acceleration)/dv, the variational update is

    (I-h A/2) J_next = (I+h A/2) J_previous + h P/m.

The average derivative sums h*(J_previous+J_next)/2 and divides by total
dt=.2. Source review confirms that neither the tested kick nor its finite
differences generates this reference; all actual probes separately verify
the assumed four-.05 partition. Reported initial discrepancies are
B=4.56891e-10 and C=2.09703e-9 s/kg, within the separate registered discrete
bound1e-7. This suite independently checks both cross-axis entries,
prescribed rows/columns, contact permutation plus common boost, replay/input
purity and explicit refinement failure. Its saved initial result is4/4 PASS.

This supports using the derivative to accelerate a root of the unchanged
DISCRETE map. It does not repair the failed continuous accuracy comparison.
Therefore a bounded numerical-only root experiment remains scientifically
interpretable under the engineering exploration scope, with accuracy
FAIL_RETAINED/process_accuracy NOT_VALIDATED. It cannot justify a physical
macrostep accuracy claim, transport validation, relaxed raw residuals or a
changed friction/time-integration model. The unavailable-cloud-response,
cost, physical-domain and complete acceptance gates remain compulsory.

The revised checker runs all original11 tests, all4 discrete tests, old103
tests, Schur tests and driver tests: thirteen suites total. It permits only
the original named assertion failure with exactly11 tests, exit1,
FAILED(failures=1), no errors and the original B-versus-ODE tolerance text.
A different method, different assertion in that method, extra failure,
changed test count or timeout blocks the numerical gate. The original
failure log is hashed and must match that same parsed exception. Every other
suite, including the discrete suite, must fully PASS. The report separately
counts executed, passing, failed and errored tests, always retains
all_tests_pass=false, and publishes only
`NUMERICAL_EVENT_PREFLIGHT_PASS_WITH_CONTINUOUS_REFERENCE_FAIL` on success.
It never relabels the original failed test as PASS or xfail.

Read-only parser checks reproduced10 passing/1 failed from the initial log
and rejected modified method, extra failure, different assertion, changed
count and timeout cases. The ordinary discrete PASS was accepted; an
unexpectedly passing original continuous suite requires renewed review.
An initial shell-quoting error prevented that parser command from starting;
the corrected command exited0. This was parser validation only, not a
preflight execution or test-suite rerun. The coordinator runs the unified
preflight after all source files are stable.

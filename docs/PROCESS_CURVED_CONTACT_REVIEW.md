# Curved contact/code/checkpoint review

2026-10-03. Independent integration review, not paper/process certification.
Scope: new hht_curved_contact, hht_curved_coupled, hht_curved_identity and
probe_hht_curved. The reviewer authored the standalone quadrature helper;
its8 reference tests were independently authored/executed by another agent.
The coupled wrapper, identities and driver were authored separately. This
discloses ownership rather than calling self-review fully independent.

This source review runs no extra trajectory. The separately preregistered
one-interval replay is reported in PROCESS_CURVED_REPLAY_REPORT, not counted
as a new accepted step. Old frozen source/FAIL evidence is not edited.

## Actual call and force/heat accounting — PASS within declared scope

`advance` validates the complete new state and builds a frozen numerical
preconditioner. Each candidate then evaluates:

1. Newmark endpoint and compatible dense midpoint x/v; current midpoint
   geometry/density NNPS and corrected gradient.
2. Existing Jaumann/FKS material update; any plastic increment rejects this
   elastic benchmark. Endpoint current-geometry physical internal force and
   endpoint boundary force are recomputed.
3. Curved contact force over the actual dense path, then the existing BE heat
   solve using signed integrated friction Q and specified reservoir fields.
4. Joint corrected acceleration, density and temperature residuals, bounded
   numerical corrections, and final acceptance/publication checks.

There is no hidden free-material velocity assignment or actuator step: all81
fixture material DOFs are free. Particle position changes through the coupled
integrator. The true endpoint total force F and numerical HHT acceleration A
remain separate; the variant does not relabel A as physical F/M.

Actual contact impulse correction DeltaJ enters the acceleration residual as
DeltaJ/(gamma*h). Actual material-work correction DeltaW enters the signed
filter-work identity as vbar·DeltaJ−DeltaW. The new ledger accumulates provided
actual contact impulses and wall work; it does NOT reuse the old two-Gauss-node
`.5*h` normal impulse assumption. Shared-node reservoir attribution explicitly
rejects in this driver; the disjoint-face fixture does not certify a future
multi-face thermal ledger.

Thermal inputs are only the original conduction/reservoir terms plus raw
friction heat and the actual elastic branch's zero plastic increment. No
DeltaJ/DeltaW, root residual, force-filter work or energy residual is assigned
to heat. Absolute values in error scales/quadrature convergence checks are
not absolute-value replacements of signed friction heat. No temperature,
yield, force or negative-Q clipping is used by the new path.

## Pointwise passivity and certificate limits

The full three-component relative velocity is projected into the contact
tangent plane. No z-only slip direction is used. The actually integrated law
is `Ft=-mu*Fn*u_t/||u_t||`; hence `Ft·u_t=-mu*Fn*||u_t||<=0` on a certified
positive-normal-force/nonzero-slip path. The supplied callback must first pass
the unchanged1% law gate and then agree with this direct formula within512eps
force-scale roundoff. Inconsistent callbacks reject, rather than being silently
replaced with a convenient force. Both callback and actual force are retained
in node diagnostics.

Whole-path gap and kinematics are **source-bound preconditions**, not facts
inferred from a handful of points. For the registered slow planes, the driver
explicitly encloses cubic particle-x extrema and the analytic normal wall
stroke range0..2e-7m, enforcing positive overlap and the original5um ceiling.
Normals are fixed and tangential wall velocity is constant. The helper checks
sampled endpoint/node identity, normal, penalty force/potential, gap and affine
tangential-wall consistency against that contract.

Old full-vector polynomial zero detection, including tangential double touch,
is reused with dummy zero tractions **only to inspect slip roots**. Those dummy
tractions never enter mechanics or heat. The obsolete endpoint-linear force
passivity condition is not reused for a curved force path. Endpoint/full-vector
zeros reject; a single component crossing zero while another remains nonzero
does not automatically reject.

The certificate interface is not an adversarial callback verifier. A caller
that lies about an arbitrary wall programme or hidden feature change is outside
the verified source contract. Source/config/artifact receipts must accompany
the checkpoint and selected callback. Rotating normals, finite-tool features,
opening, sticking, shared-particle contact constraints and arbitrary nonlinear
tangential wall histories are NOT certified by this helper. Finite8/16/32
Gauss differences are an error estimator, not a rigorous arbitrary-integrand
theorem. Failure to converge by32 nodes rejects explicitly.

## New state, hashing and restart — PASS for recorded interval

The schema/integration definition differs from legacy HHT. Legacy state cannot
silently be advanced by the new solver. Diagnostic migration preserves every
physical property, canonical acceleration and physical endpoint force, records
the old parent hash, and is explicitly not a published trajectory transition.

The complete semantic hash includes `last_interval`, not only physical arrays:
actual impulse/work, DeltaJ/DeltaW, raw Q, raw contact quadrature and kinetic
accounting are therefore tamper-detectable. Validation checks the new complete
hash first, then uses a temporary legacy-shaped view only to reuse physical
layout/parameter/history checks. That temporary view is never saved. Restart
retains canonical acceleration/history; it does not reinitialize A=F/M and
thereby lose numerical history. Runtime source and boundary-certificate
provenance remain bound through the report/receipt package, not through a
Python callable serialized inside the checkpoint.

The actual stage012→013 disk replay matches all40 PA properties, complete
semantic state/last_interval and cumulative integrals bitwise. It used the exact
stored dt, one actual advance and original source/config/boundary; no original
file changed. See PROCESS_CURVED_REPLAY_REPORT. Restart across events remains
NOT_RUN, as does any cross-platform determinism claim.

## Rejection and current stopping point

Early contact errors now preserve their diagnostics, contact identity and a
clearly rejected unequilibrated candidate. Rejected line-search attempts and
their contact diagnostics are retained. Candidate data is not a restartable
physical state simply because its arrays are finite.

The S0 driver publishes13 local stages to5us. Its next5–6us candidate is stopped
by the **independent-reference-prefix gate**, not by a declared successful stop
event. REJECTED_UNREFERENCED has an outer REJECTED/non-restartable envelope;
accepted_history remains at5us. A nested candidate-state dump is forensic
data, not authorization to extract it and bypass the outer rejection.

The path helper certifies absence of zero slip on its *numerical candidate
path*. It cannot prove absence of an event in the more accurate continuous
reference when temporal error shifts event timing. Thus the reference's
PID4/PID22 event near5.97466us and a nonzero-slip coarse candidate at6us are
not evidence that the physical event has been solved. The prefix gate is
materially necessary. Step completion, passivity and replay do not establish
temporal accuracy, classical-HHT stability or process feasibility.

## Required next event work

Extend the independent reference through located full-vector zero slip using
the actual static-disk feasibility. At the saved approach, required traction
exceeds capacity, supporting reverse sliding rather than wall-velocity clamping.
Locate joint/near-simultaneous events by bracket evidence, not iteration order.
The solver then needs an accepted event state and right-limit force/numerical
acceleration/history treatment preserving x/v/S/rho/T and separating algorithmic
storage change from heat. A post-event independent benchmark and a new bounded
replay must precede any claim that the branch is verified. Coupled effective
mass/reaction is still needed for shared faces/actuators, followed by positive
plasticity, actual open-tool/inlet integration and observable material transport.

No blocking accounting/hash defect was found within this fixed-normal,
event-free, elastic scope. The scope gaps above remain BLOCKED/NOT_RUN, not PASS.

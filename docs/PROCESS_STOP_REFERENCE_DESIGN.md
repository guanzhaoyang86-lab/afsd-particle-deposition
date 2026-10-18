# Independent first-stop/reverse reference — preregistered

2026-10-03, before unit tests or a new trajectory. New source
`scripts/reference_hht_event.py`; the frozen slow oracle is not edited. This
is the same27node current-geometry elastic/FKS-domain/thermal physical problem,
not AFSD, new friction physics, a velocity projection or time/mass scaling.

## Scope and parents

Read and hash-bind the previous independent `reference_repair1/report.json`,
its actual coarse/fine accepted histories and shared snapshot history. Verify
the old source/configuration/artifact receipts. Each grid starts from its OWN
last accepted state near5.974655us, including all physical fields and accumulated
work/heat/impulse ledger. The initial0..5.9375us shared history is retained
exactly, not regenerated or stretched. Copy only independent spatial algebra
and literal state layout; no tested backend/contact/material/thermal function
or HHT path supplies a reference answer.

First locate the previously observed full-vector zero. Discover the partner
from the initial x-reflection geometry and contact topology, then independently
check both full velocities. Mechanical mirror symmetry applies while both
remain on the common temperature-independent elastic branch; the initial
temperature gradient is retained, and its effect on the FKS scope check is
not removed. An ID is only used to associate a discovered physical contact,
never as a constitutive parameter or tracer-dependent force.

## Event localization and tolerances

Extend only the smooth incoming sliding branch for PRIVATE root candidates.
On the nearly symmetric event contacts use its incoming limit direction;
other contacts retain full-vector Coulomb. Candidate states beyond zero remain
REJECTED and are never post-event physical states. The final event state is
the full RK4 state evaluated at the located time: do not assign any component
of x/v/S/rho/T to a desired value. Store every root candidate and its margins.

Use a bracketed scalar root of the mean incoming-direction relative velocity,
but require BOTH contacts' full tangential norms <=512eps*max(1,relevant speed).
A signed-component root with nonzero transverse slip is not a full stop.
Use at most32 candidate evaluations for safeguarded secant/bisection plus
local sign-bracket certification. A supplied parent bracket may be expanded
at most4times, never past one old reference step. Desired event-time bracket
width is max(5e-19s,128eps*max(abs(start_time),initial_bracket_duration)).
The scalar SI unit tests have a larger floating-time floor by this same rule.
Require the mean root residual <=1/16 of the full-vector roundoff bound. If
group symmetry/full-zero/bracket certification fails, report BLOCKED, not a
forced simultaneous event.

At the computed event, independently recompute all internal/contact forces.
For this disjoint all-free single-face fixture, required tangent reaction is
`-F_other,t`; the static disk is feasible iff its norm <=mu*Fn. If infeasible,
the right-limit sliding direction is along F_other,t. This is not the general
shared-face/actuator mobility calculation. Feasible3D sticking remains an
explicit unsupported branch; the independent scalar feasible-stick unit is
implemented to check force balance and to prevent automatic reversal.

The right-limit direction is used only at the located event epoch and only
when the full speed is inside the registered roundoff bound. Otherwise use
ordinary full-vector Coulomb. Record and count these limit evaluations. Raw
`-ft dot u` is retained signed, including any tiny root-roundoff endpoint
residual; it is never clipped, made absolute, or converted into compensating
heat. Physical state and accumulated ledger have no jump at the event.

## New bounded independent continuation

After an infeasible-static event, use the same unscaled RK4 steps7.8125ns and
3.90625ns through8us or the first new unsupported opening/zero/yield/domain
boundary. Align new shared snapshots to62.5ns, joining old0..5.9375us and new
6us onward. Verify RK stages, endpoints and independent cubic dense paths;
the new event's roundoff-sized start alias is logged, not mistaken for a
second physical event. Preserve each rejected candidate and last accepted
state. Do not silently continue through a later event.

Coarse/fine event-time agreement <=1e-12s, same event identities and actual
static classification required. Shared-time comparison retains the existing
oracle limits: weighted RMSx1e-11m,v2e-6m/s,stress70Pa,T1e-9K,faceforce.007N,
impulse1.4e-9Ns,wallwork/rawQ2e-9J,reservoirheat2e-10J. These are numerical
reference uncertainty checks, not acceptance of the coarse HHT trajectory.
Report independently whether the event was located, the common prefix passed,
and the full8us horizon was reached. Do not use different terminal states as
a convergence comparison.

Estimated extension<20s; both references combined internal55s/external60s.
At most two evidence-based repairs, retaining initial sources/failures. If a
new branch blocks8us, stop and report it.1ms, AFSD and production are NOT_RUN.

## Independent unit references, before3D execution

- m1kg,u0+1m/s,Fother-2N,C1N,total1s: tstop1/3s, xstop1/6m,
  ufinal-2/3m/s, xfinal-1/18m, Q7/18J, Wother1/9J, DeltaK-5/18J.
- Feasible stick Fother-.5N, otherwise same: tstop2/3s, xstop1/3m,
  static friction+.5N, finalu0, Q1/3J,Wother-1/6J,DeltaK-.5J.
- Nonzero transverse constant-force case m1,u0=(1,1),Fother=(-2,0),C1:
  first ux=0 at t=(sqrt(sqrt(2)-1)+2-sqrt(2))/3, with
  uy=sqrt(sqrt(2)-1)>0; no full-vector stop is permitted.
- Malformed/nonbracketing/group-with-transverse-slip inputs reject; the input
  physical state cannot be mutated; no event ledger impulse/work/heat jump.

Scalar/algebra unit bounds1e-12 in these unit-scale SI fixtures; constant-force
2D RK4 comparison1e-10 on its fixed independent1024step integration. Literal
answers are in tests and do not call a tested contact function to generate
their expected values. Each suite has an external60s bound.

## Initial unit execution and repair 1

The first unit run retained six PASS and one FAIL in
`logs/process_stop_20261003/reference_tests_initial.log`. The literal scalar
reversal endpoint, heat, work, event and energy references all met the
preregistered 1e-12 SI bound. A redundant exact floating-point equality on the
same energy identity compared -0.27777777777777773 with -0.2777777777777778
(5.55e-17 J difference). Repair 1 changes only that redundant assertion to
the ALREADY preregistered 1e-12 bound, not any numerical model or reference.
The original test source is preserved as
`logs/process_stop_20261003/test_hht_event_reference_initial.py`. No trajectory
was run before this repair; the independent oracle source is unchanged.

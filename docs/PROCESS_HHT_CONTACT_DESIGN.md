# HHT contact/event compatibility audit

2026-10-03. Status: DESIGN / ALGEBRA_CHECKED; coupled implementation and
physical benchmark NOT_RUN by this audit. This note changes no existing solver,
accepted checkpoint, material, friction coefficient, tolerance, or backend.
It supports `PROCESS_EVENT_NEXT_WORK_PACKAGE.md` and is not a substitute for its
joint small-3D common-end-time accuracy and cost gate.

## Decision

The old event kick cannot be inserted into Newmark/HHT by relabeling its
interval-average acceleration as endpoint acceleration. A viable bounded
extension can reuse the existing simultaneous-face cone solver, but needs its
own HHT path, work certificate, event handling, and acceleration schema.

The preferred first extension is **endpoint-linear physical contact traction
with the accepted left endpoint force preserved**, exact path work, and full-path
event/passivity guards. It is exact in its Coulomb direction law for fixed
sliding direction with affine normal load, but not for general turning slip.
Turning slip requires a declared contact-quadrature approximation and an
independent direction/time-error gate. Stop if event or accuracy costs fail the
preregistered budget; positivity alone is not a Coulomb accuracy certificate.

The constant-traction cone calculation below is an audited alternative,
NOT_SELECTED for general intervals: solving a new constant traction can alter
the already accepted left endpoint force. It is permissible only when that
left force is identical or an actual registered force jump supplies its new
right limit. This note does not propose a second numerical backend or campaign.

There is one exact obstruction: continuous, unchanged acceleration history
cannot in general support transition from sliding to exact sticking. A
right-limit acceleration initialization at that physical switch is necessary.
It preserves position, velocity, stress, temperature, and plastic history; the
algorithmic-state change must be journaled separately from physical heat.

## What the existing implementation actually certifies

- `contact_kick.py::kick` uses frozen forces, normal directions, and wall
  velocities. Each accepted interval has a linear velocity path. Its impulse
  work is therefore exactly impulse dotted with endpoint-average velocity.
  Exact constant-direction stopping and bounded curved-slip midpoint intervals
  are properties of that map, not of HHT.
- `contact_cone.py::solve` solves the simultaneous-face disk QP using one
  projected mobility. Its certificate concerns the supplied affine velocity
  response; it does not choose the time discretization or certify its path.
- `coupled_boundary.py::evaluate_boundary` obtains displacement velocity by
  summing the old intervals' actual means. It labels reservoir conduction
  separately from friction heat and includes actuator reaction work.
- `process_thermal_trial.py::certify_boundary_heat` explicitly reconstructs
  endpoint-average slip on each old linear interval. Its acceptance also uses
  `M(v1-v0)=h F_frozen+J_contact+J_actuator` without an algorithmic impulse.
  These assertions cannot be inherited unchanged by an HHT variant.
- `ProcessThermalTrial.evaluate` uses one material L/T interval and conservative
  BE heat. Every new HHT/event candidate must still start from the same accepted
  material/thermal state; nonlinear trial history is not accepted history.

## Compatible kinematics and an exact missed-reversal counterexample

Use the registered convention, with acceleration history denoted `A`:

    a = .9; gamma = .6; beta = .3025
    delta = beta-gamma/2 = .0025
    M A1 = a F1+(1-a) F0
    v1-v0 = h[(1-gamma)A0+gamma A1]
    (x1-x0)/h = (v0+v1)/2+h delta (A1-A0).

A continuous path satisfying both endpoint velocities and the Newmark
displacement is, for `0 <= s <= 1`,

    v(s) = v0+s dv+6 h delta dA s(1-s)
    x(s) = x0+h[v0 s+dv s^2/2+h delta dA(3s^2-2s^3)].

Here `dv=v1-v0` and `dA=A1-A0`. `A0,A1` are HHT acceleration-history variables;
the derivative of this quadratic dense velocity need not equal them at the
endpoints. This reconstruction is a named numerical path, not proof of exact
continuous dynamics. Its mean is the actual displacement divided by `h`.

Endpoint slip signs and nonnegative total Q do not exclude an interior reversal.
In one tangent direction choose `h=1`, `v0=v1=1`, `A0=192`, `A1=-128`.
These satisfy Newmark's velocity relation, and

    u(s)=1-4.8 s(1-s), mean(u)=.2, u(.5)=-.2.

A constant friction force `f=-1` gives positive net `Q=.2 J`, but assists
the material between the two roots near `.295876` and `.704124`. Thus a check
of endpoint signs, midpoint/mean Q, or cumulative Q alone is insufficient.
For constant f and constant/linear wall speed, `f dot u(s)` is quadratic;
for endpoint-linear f it is cubic. Check its maximum at both endpoints and
every interior stationary point.
Keep the raw signed values and the independently derived roundoff budget.
Also locate zero-slip roots/minima of `|u(s)|^2`; that polynomial is quartic.

## Physical impulse, work, and the separate algorithmic impulse

For a declared linear endpoint-force reconstruction
`F(s)=F0+s(F1-F0)`, define `Fbar=(F0+F1)/2`, `dF=F1-F0`, and
`r0=M A0-F0`. Direct substitution gives the exact identity

    Jphysical = h Fbar
    Ialgorithmic = M dv-Jphysical
                 = h[(1-gamma)r0+(.9*.6-.5)dF]
                 = h[.4 r0+.04 dF].

When the previous endpoint obeys the same HHT force history,
`r0=-.1(F0-Fprevious)`, so

    Ialgorithmic = .04 h[dF-(F0-Fprevious)].

Do not report `M dv` as physical wall/actuator impulse. No artificial mass,
physical damping, or friction-capacity scaling is introduced by this identity.
At a history initialization or force jump use the general `r0` expression;
the inherited-history simplification need not hold.

The exact work along the quadratic velocity path is

    Wphysical = h Fbar dot mean(v)+(h/12) dF dot dv.

Two-point Gauss quadrature integrates this cubic integrand exactly. If used,
report the two physical quadrature impulses `Jq=h wq F(sq)` and their actual
path velocities; then `W=sum Jq dot v(sq)`. Do not compress this into a single
total impulse dotted with endpoint-average velocity. For constant contact
traction, `dFcontact=0`, and its work is exactly `h f dot mean(v)`.

An independently evaluable signed kinetic identity is

    dK-Wphysical = Ialgorithmic dot (v0+v1)/2
                  -h^2 delta Fbar dot dA-(h/12)dF dot dv.

This is a derived kinematic/force-quadrature account, not a claim of positive
algorithmic dissipation. For a conservative linear elastic force `-Kx`, use its
known endpoint potential: `dU=K(x0+x1)/2 dot dx`. Interpolated endpoint force
along a curved path is not automatically the actual `-Kx(s)`. Do not count its
quadrature mismatch twice or rename the total balance residual dissipation.
The independent smooth elastic HHT energy proof remains a compulsory separate
gate. For the existing hypoelastic material retain the existing diagnostic
classification and separate plastic work/return diagnostics.

## Complete polynomial guards for endpoint-linear contact

For a fixed tangent projector `Pt` and affine wall velocity
`w(s)=w0+s dw`, store ascending vector coefficients

    u0 = Pt(v0-w0)
    u1 = Pt(dv+6h delta dA-dw)
    u2 = Pt(-6h delta dA).

With `f(s)=f0+s df`, friction-power coefficients are

    p0=f0.u0
    p1=f0.u1+df.u0
    p2=f0.u2+df.u1
    p3=df.u2.

Evaluate power at `s=0,1` and every real root in `(0,1)` of
`p1+2p2 s+3p3 s^2=0`. Use original vector evaluations at those locations
for the reported power and an independently derived floating-point bound;
do not enlarge a physical passivity tolerance to cover a positive maximum.
If root accuracy cannot certify the maximum, reject as uncertified.

The squared-slip coefficients are

    q0=u0.u0; q1=2u0.u1; q2=u1.u1+2u0.u2
    q3=2u1.u2; q4=u2.u2.

Check endpoints and all real interior roots of
`q1+2q2 s+3q3 s^2+4q4 s^3=0` for minima. Confirm an apparent zero using
the full vector `u(s)` and the original scale-based velocity-roundoff bound.
An identically zero vector polynomial is the static branch. A near-zero
minimum is not authorization to overwrite velocity. Component roots may
provide additional brackets; a scalar projection zero alone is not a vector
zero-slip event. Retain tangent projector and wall-layer provenance.

If normal load `N(s)` is affine and nonnegative, the endpoint inequalities
`|f0|<=mu N0`, `|f1|<=mu N1` imply the same cone-capacity inequality inside the
interval by convexity. Endpoint tangency is likewise preserved for fixed
normal. Endpoint Coulomb alignment does **not** imply interior alignment or
saturation. Fixed slip direction plus `f(s)=-mu N(s) direction` does satisfy
the sliding law throughout. Any wider scope requires the independent error
gate mentioned above. A rotating normal needs its own geometry quadrature;
these polynomial statements do not silently cover rotating projectors.

## Audited constant-traction response, not selected generally

Let `g0,g1` contain the smooth material/applied forces; freeze a registered
normal geometry and normal traction `N` over one event interval. Evaluate their
approximation errors independently; use the unchanged physical penalty law
and coefficient. `P` is the free-axis projector, `M` the real physical mass,
and `f` the simultaneous tangential face forces. Put the same `N+sum(f)` into
both endpoint force layers. For constant prescribed actuator velocities,

    Rj = -(I-P)[gj+N+sum(f)]  (j=0,1)
    A1 = P M^-1[a g1+(1-a)g0+N+sum(f)].

The constrained components of carried `A0` must be consistent with the
prescribed motion. `Rj` are physical actuator endpoint forces; their impulse
and work use the declared physical interpolation. The constraint does not
silently absorb `Ialgorithmic` on free axes.

If the constant traction is genuinely the registered physical left/right
force, the actual mean relative tangent velocity is affine in f:

    mean(u_c) = T_c^T[v0+h(.5-beta)A0
                    +beta h P M^-1(a g1+(1-a)g0+N)-wall_v_c]
                + beta h T_c^T P M^-1 sum(f).

Under that condition the existing disk solver may be called with coefficient **`beta*h`**,
shared mobility `P/M`, bases `T_c`, and unchanged limits `mu*N_c`. The old
midpoint coefficient `h/2` is wrong for this response. This statement is local
to a candidate smooth endpoint force; the complete material/spatial/thermal
nonlinear problem still needs the same-history outer solve.

For each face store separately

    Jt_c = h f_c
    Qraw_c = -Jt_c dot mean(actual relative tangent velocity)
    Wmaterial_c = Jt_c dot mean(actual material velocity)
    Wwall_c = Jt_c dot wall_velocity_c
    Wmaterial_c = Wwall_c-Qraw_c.

Normal work/potential and actuator work have separate records. Static walls
do no wall work. A prescribed moving wall or actuator may supply signed
external work. Their algorithmic impulse contributions are not physical heat.

An unconstrained re-solve of this cone may change the accepted f0. Do not do
that at an ordinary smooth endpoint while claiming the same HHT force history.
For the preferred endpoint-linear extension, f0 is stored and preserved;
only f1 is a candidate unknown. Its physical quadrature impulses/work use
the preceding linear-force formulas, not this constant-force formula.

The constant-traction cone chooses a force opposite the **mean** slip. For turning slip it need
not oppose the instantaneous slip direction exactly, even if instantaneous
power is nonpositive. Therefore accept only under a registered turning-angle
and independent path/impulse/work error check; otherwise split/reject within
the fixed event budget. This is a new contact discretization and requires its
own independent reversal/two-face benchmark. No generic continuous-Coulomb
claim follows from the old cone KKT certificate.

## Stop/stick event and acceleration-history initialization

For a constant-velocity sticking wall, exact stick requires `v1=v0=wall_v`
and `x1-x0=h wall_v` in its constrained tangent direction. Newmark then requires

    A1 = -(1-gamma)/gamma A0 = -(2/3)A0
    A1 = -(.5-beta)/beta A0 = -(79/121)A0.

These differ unless `A0=0`. Consequently carrying nonzero pre-stop A0 and
claiming exact sticking is impossible at the selected parameters.

At a genuine stop or control discontinuity: locate the event by re-solving
the candidate interval from the same accepted left state (or jointly solving
event time and endpoint). Simply cutting an old dense interpolant does not
guarantee that the shorter Newmark solve has that endpoint. Preserve the event
position and velocity without snapping. Evaluate right-limit geometry,
loading, material forces, and simultaneous static feasibility. Initialize the
right-side acceleration from real-mass force balance and the active contact/
actuator constraints; when static feasibility fails, start the physically
admissible sliding branch. A reset at every ordinary smooth interval would
change the registered HHT method and is not proposed.

Journal pre/post acceleration and force histories, reason, affected DOFs,
event time, and the change in the separately derived HHT algorithmic storage.
Physical state energy has no jump when x,v,stress,T are continuous. The
algorithmic history change supplies no friction/plastic heat. Changes in h
also need the appropriate storage-account convention if storage depends on h.

## Coupled transaction and bounded completion gate

Every event subdivision must update its geometry/material time layer and BE
thermal interval explicitly. Compute the material velocity gradient from the
registered actual path/mean quadrature and paired spatial operators; do not
retain an old endpoint-average velocity field under a different label. The
chosen material/geometry quadrature remains an approximation requiring the
common-final-time reference. Friction enters BE only through the sum of
individually certified raw Q; algorithmic terms and residuals never enter heat.

The whole requested interval commits atomically only after every segment and
the joint nonlinear/thermal acceptance pass. Persist physical state, accepted
HHT acceleration, force history, active-event state, and enough journal data
to reproduce initialization/step-size changes after restart. The old average
acceleration field and a rejected candidate are not valid history sources.

Proceed directly to the one preregistered loaded 3D block with two faces,
nonzero stress/velocity, loading/unloading, actual slip reversal, BE heat, and
an independent common-end-time reference. The algebra below is only a fast
supporting check. If a proposed implementation requires unchanged acceleration
through a stop, uses only endpoint/total-Q passivity, or cannot distinguish
physical and algorithmic impulse/work, stop with that specific incompatibility.
If the valid event subdivisions erase the cost benefit, report the no-overlap
accuracy/cost result instead of accumulating equivalent microsteps.

## Evidence and provenance

The path, impulse/work identities, counterexample, and stop incompatibility are
project derivations. An independent algebra check and its command/result are
saved in `logs/process_hht_contact_design_20261003/`; they do not execute the new
coupled solver. All original FAILs and physical accepted states remain intact.

The official [OpenSees HHT documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/integrator/HHT.html)
was inspected on 2026-10-03 for the gamma/beta convention and its explicit
restriction of unconditional-stability evidence to linear problems. The
registered endpoint-force convention above is stated directly; nonlinear
interpolated-state force evaluation in another package is not substituted.
[MOOSE Dynamics](https://mooseframework.inl.gov/modules/solid_mechanics/Dynamics.html)
distinguishes inertia, numerical HHT weighting, and physical damping. Neither
source supplies a contact, thermal, or SPH validation for this proposal.

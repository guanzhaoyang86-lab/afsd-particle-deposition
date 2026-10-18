# Stop, feasible stick and force-driven reversal benchmark contract

2026-10-03. Independent derivation and test preregistration, before execution.
Scope: new event helper and an SI scalar force-driven component benchmark in
the existing backend. No 3D scout, new physical material, AFSD step, mass/time
scaling, regularized friction, velocity projection, or frozen-source edit.
Execution status at registration: **NOT_RUN**.

The existing `hht_path` endpoint map and `hht_curved_identity` were inspected.
Their force-history and numerical-acceleration histories are distinct; neither
permits carrying an unexplained old history through a force discontinuity.
The current curved-contact interval rejects full-vector zero and cannot itself
provide a post-event static/sliding mode. Tests must exercise the new helper,
not reinterpret a rejected interval as an accepted physical event.

## Independent exact scalar reference

A free point mass slides tangentially against a stationary wall with fixed
Coulomb capacity C. There is no stiffness or other stored material energy:

    m du/dt = Fother + Ft,
    Ft = -C sign(u) when u != 0,
    |Ft| <= C when u = 0.

Use m=1 kg, x0=0 m, u0=+1 m/s, Fother=-2 N, C=1 N, total duration 1 s.
Fother is an applied force, not a prescribed trajectory. On the incoming
positive-slip branch, Ft=-1 N and acceleration is -3 m/s2. Direct integration
gives u=1-3t and x=t-3t2/2. Thus the first full stop is

    t*=1/3 s, x*=1/6 m, u*=0 m/s.

The static required reaction is +2 N, outside [-1,+1] N. Sticking is
infeasible. The right-limit direction is negative, so Ft+=+1 N and total
force/acceleration are -1 N and -1 m/s2. With tau=t-t*:

    u=-tau, x=1/6-tau2/2, 0<=tau<=2/3.

The benchmark's independently derived literal targets are:

| Quantity | Exact SI value |
|---|---:|
| Stop time | 1/3 s |
| Final position | -1/18 m |
| Final velocity | -2/3 m/s |
| Contact impulse | +1/3 N s |
| Applied-force impulse | -2 N s |
| Total momentum change | -5/3 kg m/s |
| Incoming raw friction heat | 1/6 J |
| Outgoing raw friction heat | 2/9 J |
| Total raw friction heat Q | 7/18 J |
| Applied work Wother | 1/9 J |
| Contact material work | -7/18 J |
| Wall work | 0 J |
| Kinetic-energy change | -5/18 J |

In particular Wother-Q=-5/18 J. The applied force first removes energy and
then adds it; taking absolute values of applied work would be wrong. Reversal
is required physics, not a failure, and raw Q is never abs/max clipped.

## Event discovery and HHT composition

Use future-force weight a=.9, gamma=.6, beta=.3025; preserve the explicit
canonical acceleration. On a constant-force branch with A0=A1=F/m the
existing Newmark/HHT map reduces exactly to constant acceleration, for any
positive interval length. Both curved corrections DeltaJ and DeltaW vanish.

The scalar benchmark first forms a private continuing-incoming-branch trial.
It discovers the first root of its actual dense velocity (linear in this
constant-force case), recomputes the shortened incoming step, resolves static
feasibility, initializes the right history, and advances the remaining time.
The expected 1/3 s must not be fed into the event locator. A trial crossing
the stop cannot be published as one smooth Coulomb interval. Test requested
outer steps 1, 1/2 and 1/4 s: the event is not generally on an outer grid.

Only incoming/outgoing force laws prescribe dynamics. At the event, the
located physical velocity is retained; setting an approximately small
velocity to zero is not an allowed repair. A numerical root must satisfy a
separate roundoff-scale full-vector condition, not merely the coupled
stage's 1e-9 m/s force-residual allowance. Material x/v/S/rho/T/m and persistent
IDs, if present, are continuous; old accepted state and failed trials remain
unchanged. Scalar tests do not prove a 3D coupled event root or reaction solve.

## Right-limit histories and signed numerical reset ledger

For this free-mass benchmark, once the physical branch is selected, set

    Fhistory+ = Fother + Ft+    (the true right physical force),
    A+ = Fhistory+/m.

Archive Fhistory-, A-, the left parent, located event time/residual, selected
mode, static required reaction/capacity and both right-limit arrays. This is
an explicit zero-duration event restart, not permission to overwrite previous
history or reinitialize acceleration during ordinary smooth steps. A moving
wall or constrained/shared contact requires its actual coupled reaction and
wall acceleration; this scalar free-mass formula is not a general solver.

The event has dx=dv=0 and no impact law: physical impulse, physical work,
raw friction heat and kinetic-energy jump are all exactly zero. Never report
m(A+-A-) as an event impulse; its units are force. The finite force jump has
zero time integral. The kinetic identity alone does not uniquely define an
instantaneous algorithmic work; an explicit storage convention is required.

For an auditable diagnostic choose the acceleration-storage term used by the
smooth linear HHT derivation, with a *single declared reference step length*
h_ref held fixed on both sides of the event:

    d = beta-gamma/2 = (1-a)^2/4,
    Z(A;h_ref) = h_ref2*d*(A^T M A)/2,
    DeltaZ_reset = h_ref2*d*(A+^T M A+ - A-^T M A-)/2.

The signed reset ledger entry is `numerical_reset_storage_jump_J=DeltaZ_reset`.
If an augmented diagnostic K+Z is balanced, the corresponding signed reset
work is defined as this same DeltaZ_reset, so Delta(K+Z)=DeltaZ_reset while
DeltaK=0. This is numerical storage/work, not friction heat or guaranteed
dissipation. It is not a nonlinear Lyapunov or unconditional-stability proof.
The general curved/contact scheme can have other history effects, so this
partial diagnostic must not be relabeled its complete algorithmic energy.

With the single outer 1 s reversal trial, the accepted left interval has
h_ref=1/3 s, A-=-3 and A+=-1 m/s2. Therefore

    DeltaZ_reset = -1/900 J.

For other outer grids h_ref is the actual last accepted incoming subinterval;
the diagnostic jump scales with h_ref2, whereas all physical targets above
are grid-independent. If the subsequent interval uses a different h_next,
changing the storage scale contributes a *separate* bookkeeping change
`d*(h_next2-h_ref2)*A+^T M A+/2`; do not hide it in the event jump. With
h_ref=1/3 and h_next=2/3 in the reversal benchmark this is +1/2400 J.

Why restart matters: if Fhistory is changed to F+ but A- is silently carried
into the right step, its initial force-history defect is r0=M A--F+.
The existing kinetic identity then exposes extra filter impulse
`h*(1-gamma)*r0`. Here r0=-2 N; with h=2/3 s that is -8/15 N s, yielding
u_final=-6/5 m/s rather than -2/3. A small root residual would not make that
different history choice the intended physical benchmark. Conversely,
discarding both histories without an event receipt would erase this evidence.

## Feasible stick and zero-start cases

Keep m=C=u0=1 in SI units, but set Fother=-1/2 N. Incoming acceleration is
-3/2 m/s2, so t*=2/3 s and x*=1/3 m. At the event the static reaction +1/2 N
is feasible, total right force is zero and A+=0. Through final time 1 s:

    xfinal=1/3 m, ufinal=0,
    Q=1/3 J, Wother=-1/6 J, DeltaK=-1/2 J,
    Jcontact=-1/2 N s, Jother=-1/2 N s.

The static interval has nonzero reaction impulse but zero displacement,
relative work and raw heat. No kinetic sign or epsilon-sliding surrogate is
needed. For h_ref=2/3 s, DeltaZ_reset=-1/800 J. Also register initial-rest
u0=0 with Fother=-1/2 N: remain at rest, Ft=+1/2 N, no heat/work. The exact
capacity boundary |Fother|=C is statically feasible; a clearly out-of-cone
load, for example -1-1e-8 N, is not made feasible by clipping.

## Full-vector guard and one-sided direction

Use the norm of the complete tangential vector. A path with components
u(s)=(s-1/2,1) crosses a scalar zero but has norm at least 1 m/s and is never
a zero event. This is a prescribed-polynomial event-classifier unit case,
not a force-driven trajectory result. Rotate tangent coordinates and reorder
IDs where the helper's interface supports these transformations.

At a certified zero, obtain a one-sided direction from the first nonzero
coefficient of the actual local velocity expansion, not from division by
zero, a remembered scalar sign, or a velocity perturbation. For coefficients
c0=0,c1=0,c2=(1,2), the direction on either side is (1,2)/sqrt(5). For
c0=0,c1=0,c2=0,c3=(1,2), the right direction is +(1,2)/sqrt(5) and the left
direction is its negative. An identically zero expansion has no sliding
direction; it needs a static feasibility branch, not an arbitrary unit vector.

A stronger optional independent force-driven transverse reference is also
derived, but NOT_RUN: m=C=1, u0=(1,1), constant Fother=(-2,0). Polar integration
of r_dot=-2 cos(theta)-1 and theta_dot=2 sin(theta)/r gives, with
k=sqrt(sqrt(2)-1) and z=tan(theta/2),

    uy = k/sqrt(z),
    t = k/2 * [sqrt(z)-z^(-3/2)/3-k+1/(3*k^3)].

At theta=pi/2 the first component is zero at
`t_turn=(k+2-sqrt(2))/3`, while uy=k>0. This is an algebraically derived
reference, not a completed numerical integration or 3D event validation.

## Preregistered independent tests and execution limits

Own new tests: `tests/review/test_hht_zero_event.py`. API details will follow
the new helper, but literal reference results above must not be generated by
calling that helper. Planned assertions:

1. Force-driven scalar reversal with a kinematically discovered stop; exact
   independent event/final/impulse/work/raw-Q targets on the three outer grids.
2. Located feasible stick, initial rest, exact capacity boundary and a clearly
   outside-capacity case; static reaction impulse is distinct from heat.
3. Explicit right-history initialization, archived left values, zero physical
   event jump, and signed DeltaZ at a declared unchanged h_ref.
4. Reject a nonzero transverse vector passed as a zero event, unsupported
   shared-contact/actuator scope, invalid normal/mass/capacity/nonfinite input,
   and an all-zero polynomial incorrectly requested as a sliding direction.
5. One-sided linear/quadratic/higher-order direction signs, including a
   double-touch zero; no normalization of tiny nonzero data into a fake event.
6. Inputs remain unchanged on success and failure. No accepted state is
   mutated by event discovery or mode selection.

Predeclared scalar unit-scale absolute bounds: event/final time 1e-12 s;
position 1e-12 m; velocity 1e-12 m/s; impulse 1e-12 N s; work/heat/storage
1e-12 J; static equilibrium 1e-12 N; unit-direction difference 1e-12. Zero
eligibility must use the helper's separately declared roundoff-scale complete
vector bound, tested with values unambiguously above and below it. These are
component algebra bounds, not tolerances for the 3D scout.

API resolved before execution: `hht_zero_event.zero_event_limits` excludes
this contact's normal and tangential forces from `other_force_N`, then adds
both contact components internally. Scalar tests therefore use normal +z,
normal force 1 N and balancing other normal force -1 N; the physical free
tangential dynamics above are unchanged. `one_sided_slip_direction` accepts
quadratic coefficients only, so its higher-order tests cover exact second
order; the cubic parity observation above is a derivation, not a claim that
the present API accepts cubics. `acceleration_storage_jump` exposes the
declared partial numerical-storage diagnostic without mutating histories.

Each test execution has an external ceiling of 60 s. Preserve the initial
command, stdout/stderr, source hashes, exit code and every FAIL; at most two
evidence-driven repairs. No physical trajectory is authorized by this contract.

## Initial component execution

Command: `.envs/pysph/python.exe -B tests/review/test_hht_zero_event.py`.
Result: **8 PASS**, 0.079 s unittest, 0.3015393 s external, exit 0 under the
60 s external ceiling. Initial stdout/stderr and supervision metadata are
`logs/process_stop_20261003/zero_event_tests_initial.log`,
`zero_event_tests_initial_stderr.log` and `zero_event_tests_supervision.json`.
No evidence repair was needed. The expected NumPy overflow warning for the
deliberately extreme but finite 1e308 invalid-force input is preserved; the
helper rejects the nonfinite derived norm and that rejection test passes.

The force-driven scalar reversal and feasible-stick targets pass at outer
steps 1, 1/2 and 1/4 s. Stop time is discovered from the private kinematic
trial; the independent literal target is used only for comparison. In these
executed cases the located velocity is exactly zero in stored floating-point
output, with no projection. The reversal's single-outer-step raw Q is
0.38888888888888895 J; applied work is 0.11111111111111122 J; final velocity
is -0.6666666666666667 m/s. The signed acceleration-storage reset is
-0.0011111111111111367 J at h_ref=1/3 s, separate from physical work and heat.

The tests also cover feasible/marginal/infeasible static reactions, an
accelerating-wall force balance, complete-vector rejection, shared-contact
and actuator-scope rejection, unchanged inputs, resolved first-order and
exact second-order one-sided limits, unresolved derivative cancellation and
invalid input. The optional 2D force-driven polar reference remains NOT_RUN
in this suite. Passing these component tests does not locate, solve or
authorize the subsequent 3D coupled stop event.

## One-sided contact-integral test preregistration

New helper `hht_stop_contact.py`, owned separately, is to be tested by
`tests/review/test_hht_stop_contact.py` before any 3D event-driver use. These
are single-particle prescribed-path component tests, not a second claim of
force-driven dynamics. Use fixed normal +z, normal force 3 N, mu=1/3 and
capacity C=1 N. All analytic references below are calculated independently
of the tested integrator. The same 1e-12 N s/J absolute bounds apply.

- Incoming left-zero path: u(t)=1-3t along +x, h=1/3 s. Its one-sided force
  is (-1,0,3) N through the left limit, J=(-1/3,0,1) N s, raw Q=1/6 J and
  stationary-wall material work -1/6 J. The endpoint is zero, not an interior
  continuation through zero.
- Outgoing right-zero path: u(t)=-t along x, h=2/3 s. Outgoing direction -x
  selects force (+1,0,3) N, J=(2/3,0,2) N s, raw Q=2/9 J and stationary-wall
  material work -2/9 J. The declared direction must match the first nonzero
  derivative of the actual dense path, not only its far endpoint.
- For either path, an affine common wall velocity w(s)=w0+s*dw is added to
  material velocity, with matching acceleration. Then raw Q/J are unchanged,
  Wwall=J dot (w0+dw/2), and Wmaterial=Wwall-Q. Constant total force means
  DeltaJ and DeltaW against its endpoint chord must vanish to roundoff.
- PRIVATE incoming-branch continuation: u(t)=1-3t over h=1/2 s crosses zero
  at t=1/3 but deliberately retains incoming tangential force -1 N on the
  naturally collinear extension. Its signed dx=1/8 m, J=(-1/2,0,3/2) N s,
  Q=+1/8 J and material work -1/8 J. Some quadrature nodes have **negative**
  raw Q contributions; they must be retained, not clipped. The return must
  explicitly remain private, unaccepted/non-restartable and not physically
  passive/zero-event certified. Positive integrated Q does not authorize it.
- The same PRIVATE continuation to h=1 s has Q=-1/2 J and must reject the
  negative integrated heat instead of applying abs/max. Ordinary strict
  contact with an interior zero remains rejected.
- A nonzero transverse endpoint component is not left-zero. A PRIVATE path
  with transverse direction change is not a naturally collinear continuation.
  Reject both. Reject a left-zero polynomial with an earlier full-vector zero
  and a right-zero path whose initial dense derivative opposes its declared
  outgoing direction, even if a later velocity points in the declared sense.
- Success and failure leave the input ParticleArray and parameters unchanged;
  missing/bad event contexts and invalid direction data cannot manufacture an
  accepted or restartable state.

The helper's generic return interface and event metadata will be inspected
before tests execute. Component test ceiling is external 60 s. Preserve the
initial source hash, command and logs; every FAIL stays evidence, with at
most two evidence-based repairs. Registration status: **NOT_RUN**.

API inspection before execution: the caller's boundary callback must itself
declare the matching one-sided endpoint force or PRIVATE incoming-branch
force, and the helper checks consistency; it does not silently replace an
ordinary callback. The manufactured test callback explicitly provides this
declared law. PRIVATE return flags are `private_continuation=true`,
`physical_contact_path_certified=false`, `whole_path_passive=false`,
`accepted=false`, `restart_allowed=false`, `can_publish=false`.

## Initial one-sided contact execution

Command: `.envs/pysph/python.exe -B tests/review/test_hht_stop_contact.py`.
Result: **8 PASS**, 0.082 s unittest, 0.5209038 s external, exit 0 within the
60 s external ceiling. No repair was needed. Initial stdout/stderr and
source-bound supervision are preserved under `logs/process_stop_20261003/`
as `stop_contact_tests_initial.log`, `stop_contact_tests_initial_stderr.log`
and `stop_contact_tests_supervision.json`.

The exact incoming/outgoing impulse, raw heat and material/wall-work values
pass both stationary- and affine-wall cases. The private half-second
continuation returns raw Q=0.12499999999999999 J while retaining six negative
node contributions; the minimum is -0.009862928226875255 J. Its physical,
passivity, acceptance, restart and publication statuses remain explicitly
false. Negative aggregate private Q, an unmarked interior crossing, nonzero
transverse endpoint slip, noncollinear private extension, an earlier complete
zero, and wrong initial dense derivative all reject without changing input
particles/parameters. In particular, the derivative test has a far endpoint
aligned with the declared direction but a different actual initial direction,
with no other full-vector zero; endpoint-only direction checking cannot pass
that test accidentally. Missing outgoing direction and malformed contexts
also reject.

These are prescribed-path integral/validation component results. They are
not a 3D event location, constrained static solve, accepted HHT event-driver
trajectory, or physical process-accuracy certification.

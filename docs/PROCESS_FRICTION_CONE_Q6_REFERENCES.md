# Q6: independent two-tangent Coulomb cone references

Status at definition: **NOT_RUN**. Execution results belong in the Q6 report and
logs, not in these reference values. The following are independent static
force-balance and piecewise-integral derivations, not evaluations of the contact
implementation. SI units apply throughout. Constants are manufactured benchmark
data, not material or friction parameters inferred for AA6061 or from the paper.

## Scope and signs

For an orthonormal tangent basis, let u and w be two-component material and wall
displacements, r=u-w, and f the tangential force on material. Normal compression
is lambda>=0; Coulomb capacity is C=mu*lambda. The isotropic law is

    stick: delta r = (0,0), ||f||_2 <= C
    slip:  f = -C * delta r / ||delta r||_2.

The latter is the incremental endpoint condition. Its use over a finite interval
does not by itself prove a continuously turning slip path was resolved. Event
splitting and the actual acceptance contract must state their applicable paths.
Independent component bounds |f1|<=C and |f2|<=C define a square, not this disk.
For example f=(.8,.8) N violates a unit disk (norm=sqrt(1.28) N), although both
component tests pass. Conversely f=(.6,.8) N is exactly on the unit disk.

For each accepted interval, work is evaluated from actual signed vector forces:

    Wwall = f_mid dot delta w
    Wmaterial = f_mid dot delta u
    Q = -f_mid dot (delta u-delta w) = Wwall-Wmaterial.

There is no abs(Q), positive-part operator, or residual-energy heat source. Q is
mechanical dissipation, not solved temperature. Normal work is accounted for
separately; its vanishing in a fixed-normal example must not be generalized to
moving normal boundaries. Zero slip does not imply zero wall/material work.

## A. Oblique proportional sliding, reversal and return

Let k=1 N/m in both tangents, C=1 N and e=(3/5,4/5). Wall w=s*e follows scalar
program s:0 -> 3 -> 0 m. Static balance gives f=k*u. Write u=q*e. The circular
cone then gives exactly the scalar piecewise law, independently derived here:

1. s:0->1: stick; q=s, f=s*e.
2. s:1->3: slide; q=1, f=e, relative increment opposite e.
3. s:3->1: stick; q=s-2, f=(s-2)*e.
4. s:1->0: reverse slide; q=-1, f=-e, relative increment along e.

| endpoint | wall s (m) | u (m) | f (N) | interval Q (J) | interval Wwall (J) |
| --- | ---: | --- | --- | ---: | ---: |
| first threshold | 1 | (.6,.8) | (.6,.8) | 0 | .5 |
| wall reversal | 3 | (.6,.8) | (.6,.8) | 2 | 2 |
| reverse threshold | 1 | (-.6,-.8) | (-.6,-.8) | 0 | 0 |
| return | 0 | (-.6,-.8) | (-.6,-.8) | 1 | 1 |

Thus final u=(-.6,-.8) m, Q=3 J, Wwall=3.5 J, Wmaterial=delta U=.5 J.
An extra sample at wall s=2 during the return has u=(0,0), f=(0,0): it is still
sticking, not reverse sliding. Flipping f as soon as wall velocity reverses is
wrong. Checking positive Q alone would not detect that error.

For two reversals s:0->3->-3->0, endpoint scalar arrays are

    s = [0,1,3,1,-3,-1,0]
    q = [0,1,1,-1,-1,1,1]

with u=q*e, f=q*e at these endpoints. Total Q=7 J, Wwall=7.5 J, final U=.5 J,
and final u=(.6,.8) m. These terminal values cannot be interchanged with the
one-reversal path. Algebraic static PASS does not supersede the existing dynamic
m=1 kg friction-reversal test or the old real AFSD failure evidence.

## B. Tangent-coordinate change, not physical tool rotation

Let R=[[0,-1],[1,0]] be an orthogonal tangent-coordinate transformation. Under
u'=R*u, w'=R*w, f'=R*f and K'=R*K*R^T, all dot products and norms are unchanged.
The e of A becomes e'=(-.8,.6). All event scalar s, Q=3 J, Wwall=3.5 J and U=.5 J
must remain unchanged; final u'=(.8,-.6). A second useful rotation is
R=[[.6,.8],[-.8,.6]], for which e'=(1,0).

Transform stiffness, loads and wall path together: rotating only wall direction
in an anisotropic discrete body is a different physical experiment. A basis
covariance test is neither a rotating tool simulation nor a Jaumann objectivity
test with nonzero initial stress. Those physical tests remain NOT_RUN in Q6.

## C. Static wall and zero-relative-slip translating preload

Static wall w=(0,0), moving spring base b=s*e, k=C=1, s:0->3. Force balance is
f=u-b. For s<=1, u=0 and f=-s*e; subsequently u=(s-1)*e and f=-e. Terminal
u=(1.2,1.6) m, f=(-.6,-.8) N, Q=2 J, Wwall=0, Wmaterial=-2 J, Wbase=2.5 J,
U=.5 J. Negative material work at a static wall is not negative dissipation.

For common translation, choose initial spring extension .5*e and force .5*e.
Translate wall, material and spring base together by 2*e. Then delta r=0,
Wwall=Wmaterial=1 J, Wbase=-1 J, Q=0, and delta U=0. This is a work/algorithm
fixture, not a claim that the clamped 3D block admits rigid translation.

## D. Two contacts coupled through shared free tangent displacement

Each coordinate below is a two-vector. Take the symmetric positive-definite
6-by-6 stiffness, in N/m,

    K = [[2I, 0, -I],
         [0, 2I, -I],
         [-I,-I, 3I]],    I=identity(2).

Two contact displacements u1,u2 act against unit circular friction disks; z is
an internal free displacement. Equilibrium gives z=(u1+u2)/3. Eliminating z gives
effective blocks [[5I/3,-I/3],[-I/3,5I/3]]. Equal contact motion along e therefore
has stiffness 4/3 per contact, not isolated diagonal stiffness 2.

With both walls w=.75*s*e and s:0->3->0, the literal endpoint scalars are

    wall/e = [.75, 2.25, .75, 0]
    u1/e = u2/e = [.75, .75, -.75, -.75]
    z/e = [.5, .5, -.5, -.5]
    f1/e = f2/e = [1,1,-1,-1].

Each contact has Q=2.25 J, Wwall=2.625 J, Wmaterial=.375 J. Totals Q=4.5 J,
Wwall=5.25 J and terminal U=.75 J. The wrong isolated onset would be .5 m rather
than .75 m. This benchmark tests common DOF coupling, not merely two independent
copies of a one-contact solver.

## E. Normal reaction coupled to oblique tangential load

Use coordinates (normal displacement d, two tangential u components), mu=.5,
and the positive-definite stiffness

    K = [[2, .3, .4],
         [.3, 1, 0],
         [.4, 0, 1]] N/m.

Fix d=-1 m. From force-on-material (-lambda,f)=K*(d,u), along u=q*e:

    lambda = 2-.5*q,     f=(q-.5)*e.

Starting q=0 is feasible with lambda=2 N and f=-.5*e N. Increasing wall s with
stick q=s reaches q=6/5, lambda=7/5 N, f=(.42,.56) N, normf=.7 N. The opposite
threshold is q=-2/3, lambda=7/3 N, f=(-.7,-14/15) N, normf=7/6 N.
These values follow from the current coupled normal reaction, not a frozen
capacity. Normal wall work remains zero because d is fixed while lambda varies.
This is a prestressed initial fixture and must not be silently confused with A.

## F. Non-collinear endpoint: slip direction is not wall direction

This is an incremental endpoint KKT oracle, **not** a continuous turning-path or
exact work-integration reference. Set K=diag(2,1) N/m, C=1 N, accepted u=0 and
w=0, and request w=(1.5,2.4) m. The literal solution is

    u=(.3,.8) m, f=(.6,.8) N,
    delta r=(-1.2,-1.6) m, ||delta r||=2 m.

Indeed K*u=f and f=-delta r/||delta r||, with normf=1. The actual slip direction
differs from normalized wall direction: using wall velocity in place of solved
relative motion fails this reference. The associated endpoint solution is
unique because K is positive definite and the disk dissipation is convex.
Under the 90-degree transformation of B, K'=diag(1,2), w'=(-2.4,1.5),
u'=(-.8,.3), f'=(-.8,.6), delta r'=(1.6,-1.2). Both coordinates must give the
same physical solution. Do not infer an exact cumulative Q for the continuous
path from this single large incremental endpoint.

## Preregistered tolerances and limits of interpretation

- Algebraic unit-scale displacement/force/work references: absolute 1e-10 m,
  1e-10 N, 1e-10 J, relative tolerance zero. Event position: 1e-10 m. These are
  inherited Q5 roundoff-level allowances against O(1) benchmark signals.
- The linear algebra fixture may use force/slip solve tolerances 1e-12 N/m
  respectively and 48 event bisections for <=6 s requests, targeting brackets
  <=1e-12 s. This is a search budget, not mass or physical-time scaling.
- The actual 3D body retains the separately frozen Q6 contract thresholds and
  bounded run budget. An X-minus-clamped finite cube is not an exact uniform
  two-component spring: no K=G*A/L exact continuum equality is asserted there.
- Check each contact's raw signed work, disk KKT, direction, sticking intervals,
  event order, total equilibrium, accepted-state rollback and passive labels.
  A positive total Q does not permit cancellation of wrong individual contacts.
- Nonlinear 3D material/work, finite rotations, changing normal/feature, thermal
  conversion, large slip and new collision/recontact are not certified by these
  algebraic references. The old material-cycle and spatial FAILs remain intact.

No tests or real run are asserted PASS by creation of this document.

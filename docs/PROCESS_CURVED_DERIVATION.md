# Impulse-corrected HHT contact: independent derivation

2026-10-03. Design before a new trajectory. This derives one explicitly named
NUMERICAL_VARIANT inside the existing backend. It is not a new material law,
heat-only patch, proof of AFSD completion, or inherited classical HHT stability
certificate. The old endpoint-linear implementation and failed evidence remain
unchanged. No integration was performed to produce this derivation.

## Definitions and what changes

Let M be the unchanged diagonal physical mass matrix; h is the actual step in
seconds. Positions, velocities, forces and work use global SI units. Define

    eta=1-a; gamma=.5+eta; beta=(1+eta)^2/4
    d=beta-gamma/2=eta^2/4; q=gamma*a-.5=eta*(.5-eta).

Use the same Newmark endpoint kinematics and compatible dense velocity:

    v1-v0=h*((1-gamma)*A0+gamma*A1)
    x1-x0=h*v0+h^2*((.5-beta)*A0+beta*A1)
    v(s)=v0+s*(v1-v0)+6*h*d*(A1-A0)*s*(1-s), 0<=s<=1.

Thus `integral v dt=x1-x0` and

    dx=h*vbar+h^2*d*DeltaA, vbar=(v0+v1)/2.

F0 and F1 are the TRUE total physical endpoint forces, including actual contact,
internal forces and noncontact applied forces. Define the baseline linear
force `Fl(s)=F0+s*(F1-F0)`. For each corrected contact component, evaluate its
actual chosen interval force f(s), using the same candidate path and actual
wall kinematics, and subtract its own endpoint chord fl(s). Sum over contacts:

    deltaf(s)=sum_contacts(f(s)-fl(s))
    DeltaJ=h*integral_0^1 deltaf(s) ds                   [N s]
    DeltaW=h*integral_0^1 deltaf(s) dot v(s) ds          [J]
    Fphysical(s)=Fl(s)+deltaf(s).

The unmodified internal-force endpoint chord is still a temporal approximation;
this proposal does not secretly claim it is the true nonlinear internal force
everywhere along the path. If normal contact is evaluated from actual moving
geometry rather than its chord, its correction belongs in DeltaJ AND DeltaW,
not only the tangential correction. Every chosen corrected component must use
the same quadrature contributions in force, impulse and mechanical work.

For fixed-normal active Coulomb sliding, the intended contact force is
`fn(s)*n - mu*fn(s)*ut(s)/norm(ut(s))`. Open contact or full-vector zero slip
is outside this interval branch and must be located/rejected, not regularized.
The new evaluator must not mistake a nonzero vector's z sign change for zero.

## Proposed coupled residual

Define the force-valued residual

    R=M*A1 - a*F1 - (1-a)*F0 - DeltaJ/(gamma*h).

The changed term therefore enters the mechanical residual, not merely thermal
bookkeeping. DeltaJ depends on the same x/v/A candidate that determines actual
contact, material and thermal state and must be recomputed during the root.
Using an old/predictor DeltaJ after the root changes the equations and is not
this variant. The existing bounded nonlinear solver may use an approximate
numerical preconditioner, but acceptance must check the new raw residual.

Let `r0=M*A0-F0`. Direct substitution into Newmark gives

    M*(v1-v0) = Iphysical + Ifilter + Iroot
    Iphysical = h*(F0+F1)/2 + DeltaJ
    Ifilter = h*((1-gamma)*r0 + q*(F1-F0))
    Iroot = gamma*h*R.

The original endpoint-linear force-filter impulse formula is unchanged. The
actual contact impulse is now the force integral included in Iphysical. This
is the required algebraic reason the proposal is not a heat-only substitution.

If a code reports acceleration residual `ra=A1-Amap`, then `R=M*ra`; do not
confuse its units with force. Store and test per-particle momentum balance,
not just a total that can conceal cancellation. Multiple contacts contribute
vectorially before the joint residual; this does not supply independent scalar
stopping caps or solve multi-face static feasibility by itself.

## Exact kinetic/work identity

For the baseline chord and this dense path, polynomial integration gives

    Wlinear=Fbar dot dx + h/12*(F1-F0) dot (v1-v0)
    Wphysical=Wlinear+DeltaW
    Wfilter_linear=Ifilter dot vbar
                 -h^2*d*Fbar dot (A1-A0)
                 -h/12*(F1-F0) dot (v1-v0)
    Wfilter_curve=DeltaJ dot vbar-DeltaW
    Wroot=Iroot dot vbar.

Taking the dot product of momentum balance with vbar proves exactly

    DeltaK = Wphysical + Wfilter_linear + Wfilter_curve + Wroot.

`Wfilter_curve` is indispensable: keeping the old kinetic diagnostic after
adding DeltaJ would leave an unaccounted defect. All filter terms are signed
NUMERICAL work; they are not assumed nonnegative and never become heat.
The equivalent compact formula is
`Wfilter=(Ifilter+Iphysical) dot vbar-Wphysical`.
For implementation, retaining the expanded terms avoids hiding physical work
and makes the new term independently testable. Quadrature consistency makes
this algebra exact for the quadrature representation; independently estimated
quadrature error is a separate approximation, not absorbed in Wroot or heat.

For each actual tangential contact, independently accumulate

    Qraw=-integral ft(s) dot (v(s)-w(s)) dt
    Wwall=integral ft(s) dot w(s) dt
    Wmaterial=integral ft(s) dot v(s) dt
    Wmaterial-Wwall+Qraw=0.

The material receives only its declared partition of Qraw; wall reservoir
exchange is separate. Normal-contact work and penalty-potential changes stay
separate from Qraw. For actual fixed-normal penalty geometry, compatible x/v
and wall motion imply `DeltaUn + integral fn*n dot (v-w) dt=0` up to quadrature
and root error. A temporal normal-potential residual is diagnostic, not heat.

Under a common constant velocity boost c, DeltaJ is unchanged while
DeltaW changes by `c dot DeltaJ`; vbar changes by c. Therefore the new
`Wfilter_curve` is boost-invariant, and raw relative friction heat is invariant.
The other algorithmic work transforms with its corresponding filter impulse,
as required by the momentum balance.

## History and checkpoint semantics

F1 remains the physical endpoint force. NEVER rename
`F1+DeltaJ/(a*gamma*h)` as the physical force to make an old root appear to fit.
A1 is the numerical HHT acceleration and generally is NOT `F1/M`. On accepting
a new-variant step, retain both arrays unchanged together with DeltaJ, DeltaW,
quadrature definition/tolerance, raw interval force/velocity samples, all signed
filter/root work, parameter/source hashes, and parent state hash. Give the
checkpoint a new integration definition/schema; old files must not silently
acquire the new meaning.

An explicitly named import of an old accepted endpoint may preserve its x/v/
S/rho/T, true F0 and numerical A0, recording old schema and hash. The equations
above remain valid for arbitrary retained r0. This is not permission to accept
a previously rejected state or reset old A0 to improve a result.

At a subsequently located physical/control event, any right-limit numerical
history reset must preserve physical state and record both force/acceleration
histories and signed numerical-storage change. This derivation does not solve
zero-slip/static feasibility or certify the full nonlinear algorithmic-energy
jump. No mass, modulus, physical time or initial velocity is rescaled.

## Reduction, independent checks and limitations

If every corrected actual f(s) equals its endpoint chord, DeltaJ=DeltaW=0.
The residual, impulse, work and previous HHT formulas reduce exactly to the
old scheme. This includes constant-normal straight Coulomb slip and any
straight slip with linearly varying capacity. This is an algebraic identity,
not a visually similar result.

An independent analytic curved contact check that is compatible with a
linear-velocity HHT path is `ut(s)=(1,s,0)`, capacity C and duration h:

    Jx=-C*h*asinh(1), Jy=-C*h*(sqrt(2)-1), Jz=0
    Qraw=C*h/2*(sqrt(2)+asinh(1)).

With constant wall speed w, `Wwall=w dot J` and `Wmaterial=Wwall-Qraw`.
This is a prescribed-path contact integration benchmark, not permission to
drive free AFSD material along that path. A circular full-vector rotation
benchmark can supplement it. Tests should independently verify these values,
constant/linear reduction, momentum/work identities, boosts, negative raw
heat rejection, complete-zero detection and candidate nonpublication.

For a smooth force, chord error is O(h^2), DeltaJ is O(h^3), and the added
acceleration term is O(h^2). It is therefore consistent with second-order
kinematics, assuming the underlying smooth coupled method/root and adequate
quadrature accuracy. This is not a proof of order at events or of stability.

Crucially, the next step's r0 includes the previous curved correction. Even
the nominal a=1 undamped control can then have a nonzero filter impulse on the
next curved-force step. Classical HHT linear unconditional-stability and
augmented-energy proofs do NOT transfer automatically. Energy and momentum
terms must be reported and the named variant compared against independent
same-time references. Positive Q alone is insufficient. Near-zero direction
singularities, nonlinear root conditioning, temporal force curvature and
quadrature error may still block affordable large steps.

Replacing a chord by a pointwise Coulomb evaluator makes a same-function
pointwise cone check nearly tautological. Preserve the physical-law bound,
but add independent integration references and a controlled quadrature-error
estimate rather than declaring fidelity from that check alone. Do not change
the old 1% result or use this derivation to claim an unrun long process passes.

## Identity helper preregistration, before unit execution

New module `hht_curved_identity.py` reuses the old endpoint-linear decomposition
without changing it. Its API accepts per-particle `delta_impulse_Ns` with the
same shape as velocity and `delta_material_work_J` with the mass shape. It
returns corrected physical work as `physical_force_work_J`, retaining the
chord as `baseline_linear_force_work_J`, explicit correction impulse/material/
filter work, corrected root impulse/work, and actual returned-velocity
kinematic/momentum/kinetic residuals. It does not infer zero residual from the
nominal root formula or relabel work as heat.

Independent tests registered now: one hand-calculated loaded node with nonzero
old history and DeltaW != vbar dot DeltaJ; arbitrary root and returned-velocity
defects; exact zero-correction reduction against the old map; a=1 with nonzero
history; Galilean transformation; analytic ut=(1,s,0) Coulomb integrals; and
nonfinite/shape rejection. Unit-scale SI algebra uses absolute1e-12 tolerance
(about4500 machine eps, allowing several vector products/subtractions), with
explicit1e-13 on the small hand-literal terms. No physical trajectory or
acceptance threshold is changed. A separate external60s bound applies to the
suite; all first results and any repairs are retained.

Identity execution: **7 PASS**, .050 s unittest / .193986 s external, exit0,
under subprocess timeout60s. No repair was needed. Command:
`.envs/pysph/python.exe -B tests/review/test_hht_curved_identity.py`.
Evidence `logs/process_curved_20261003/identity_initial.log`.
Source SHA256 `121695a7f97f6d2ec6771b6a54eb79ac2de51d43b7f089912fff61eb03c46059`;
test SHA256 `211d979d1fdde1239a204ce72643bea75d9ad2ae2c2842e4cc92926de988bca0`.
This certifies the bookkeeping helper's stated algebra/validation tests only,
not a coupled root, contact quadrature accuracy, process trajectory or stability.

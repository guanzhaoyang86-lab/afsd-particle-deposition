# Saved342 action defect: independent read-only causal review

2026-10-04. No new physical evaluation, kick, material update, NNPS search,
matrix construction, preconditioner application or timestep ran. This review
reads the existing action NPZ and real342 base/plus/minus/parent archives.
The new NumPy-only reader and full signed results are under
`logs/process_bridge_action_review_20261004/`. One invocation: internal
0.7240623s, external0.9621477s, exit0. All172 bound source/evidence hashes
matched before and after. Old FAILs and the source
`NONSMOOTH_NOT_CERTIFIED` classification are retained.

## Readback definition and numerical limitations

The reader first checks source/artifact hashes, persistent-ID correspondence,
the stored `delta_candidate`, `delta_residual`, and exact action identity
`E=delta_C+delta_correction=delta_C-P*delta_R`. It uses no tested contact,
material, spatial or preconditioner function to construct a reference.
Its source declares signed finite-difference decompositions, not a derivative
or causal proof. Algebraic reconstruction defects below are reported raw as
floating-point consistency diagnostics; no new physical PASS threshold,
epsilon refinement or outcome-dependent filter is introduced.

The plus/minus pair inventories and midpoint densities are exactly equal.
Therefore volumes are equal and pressure difference is zero in this specific
direction. With saved `Gi=B_i^T*g`, average denoted by an overbar and difference
by delta, the stored paired force gives the symmetric product identity

```text
delta Fij = Vi*Vj*(delta sigma_i*Gibar + sigmabar_i*delta Gi
                  +delta sigma_j*Gjbar + sigmabar_j*delta Gj)
delta Gi  = Bibar^T*delta g + delta Bi^T*gbar
delta sigma = .5*delta S1                 (S0 fixed, delta p=0)
```

Scatter each pair with opposite signs to its two saved particle indices.
The factor **.5** is retained: the actual force uses `(S0+S1)/2`, not S1.
Stored velocities and Gi likewise separate the gradient difference into
velocity-change and operator-change terms. This is array algebra on saved
operators, not a reconstruction of K or P.

For an additional material decomposition, infer the original trial stress
from saved output stress and plastic increment:
`qtrial=qout+3G*dep`, `Strial=Sout*qtrial/qout`. No return solve is repeated.
Then split

```text
delta Sout = 2G*h*dev(delta L)
           +[delta Strial-2G*h*dev(delta L)]
           +[delta Sout-delta Strial].
```

The middle bracket is a corotation/prestress remainder, not an exact isolated
Jaumann tangent. The last is the finite difference of the radial-return
correction. The inversion and subtraction have floating-point error; the
reported small remainder must not be treated as an exact constitutive test.
The first-order prestress-spin term is separately diagnostic only.

Raw reconstruction maxima: original force3.40e-16N; symmetric stress/geometry
force difference7.93e-16N; full decomposed force sum8.62e-16N; gradient
split3.36e-13/s; midpoint/returned stress split9.31e-10Pa. At fully free,
no-contact particles, saved mapped velocity differences equal
`h*delta F/m` to2.36e-17m/s and mean-velocity differences equal
`.5*h*delta F/m` to3.08e-17m/s.

## CONFIRMED: location and source omissions

Both existing P variants place the maximum mechanical action defect at
**PID268,z**. It has no contact and no event in base, plus **or** minus,
zero parent ep, and nonzero parent stress norm3.26662MPa. It is a direct saved
pair neighbor of the sole newly plastic PID269, at distance0.499972mm.

|Field|Structural max|Continuity max|Particle/component|
|---|---:|---:|---|
|Position E|3.49100e-15m|3.57738e-15m|268,z|
|Endpoint velocity E|1.3963441e-7m/s|1.4308975e-7m/s|268,z|
|Mean velocity E|6.9817206e-8m/s|7.1544873e-8m/s|268,z|
|Density E|2.81120e-8kg/m3|6.53039e-9kg/m3|198 /270|
|Temperature E|3.26054e-10K|same|269|

PID269 plus its71 direct saved-pair neighbors contain96.9656% of structural
velocity E's squared L2 norm and99.4832% for the continuity variant. These
are support fractions of one signed direction, not physical error reductions.
The top eight structural mechanical-error particles are269 or its direct
neighbors; every one has no contact in all three saved evaluations.

Source evidence is unambiguous:

- `structural_mobility.py:172-209` builds `K=Bgrad^T*V*Celastic*Bgrad` from
  frozen accepted geometry, masses/densities and shear/bulk. It does not read
  S0, ep, rate or temperature; its own scope explicitly excludes prestress,
  plasticity and evolving geometry.
- `stage_preconditioner.py:47-51,96-115` uses that elastic matrix for the
  mean-velocity correction and lifts it into x/v. Density and T remain direct
  defects. `stage_continuity_preconditioner.py:25-34,73-97` adds reference-
  density pressure/continuity coupling, not a material plastic tangent.
- `process_thermal_trial.py:144-156` actually recomputes corrected operators
  at candidate midpoint geometry, updates from the same S0/ep0 with midpoint
  temperature and candidate gradient, and uses `.5*(S0+S1)-pI` in its force.
- `material.py:29,53-73` actually includes corotational stress transport,
  rate-sensitive FKS yield and48 bisection return iterations. The existing
  structural model omits all of these constitutive derivatives.

These are confirmed *scope omissions*, not proof each omission is important
to this direction. An exact linear solve of the approximate P does not make
P an exact inverse of the full nonlinear residual derivative.

## SUPPORTED: missing same-branch rate-sensitive plastic response is first

At PID268 the archived force difference is almost entirely the stress term:

|Signed z-force difference contribution|N|
|---|---:|
|Actual delta F|+9.8764834e-7|
|Radial-return correction difference|+9.7637182e-7|
|Linear elastic response to mean-velocity difference|+1.1276458e-8|
|Changed force-assembly geometry|+6.77387e-14|
|Linear-gradient geometry contribution|−8.05230e-15|
|Corotation/prestress remainder|+4.43174e-15|

Across all particles, return-force L2 is2.27323e-6N, compared with actual
delta-force L2 2.40236e-6N, elastic velocity-response5.54880e-7N,
force-geometry1.85653e-11N, gradient-geometry1.05797e-12N, and corotation
remainder2.01160e-12N. These signed components may cancel; their norms are
not an additive partition or a causal percentage.

The return difference originates only at PID269 and its force contribution
is confined to its saved-pair support. At269:

```text
rate+ =1.0478804395827024/s; rate- =1.0478573073321584/s
dep+  =3.0306131330378787e-5; dep- =3.0307293896752590e-5
delta S1,zz           =61.2063046808Pa
delta inferred trial,zz =.0553165199Pa
delta return,zz       =61.1509881609Pa
```

Both rates are strictly above the original floor1/s and both outcomes have
positive plastic increment. A larger rate raises the FKS yield stress and
reduces the returned plastic increment here; the returned-stress response is
far larger than the purely elastic trial-stress response. This supports
prioritizing the **same-branch FKS rate-aware plastic tangent**, not first
patching friction, Jaumann transport or geometry on the strength of a name.
The no-contact map-force identities also rule out a local friction impulse as
the direct source of PID268's saved map difference.

## POSSIBLE / not established

No separate return/geometry/spin force contribution was passed through P's
inverse. The return contribution's dominance in local delta F is therefore
not a unique quantitative attribution of `delta C-P delta R`. P also transmits
coupling from other particles. Frozen accepted versus midpoint geometry,
prestress transport, density Schur approximation and contact/event response
remain real omissions and may matter in other directions or later states.
The present rate/plastic branch evidence does not make the whole342 response
smooth: its original contact signature is NONSMOOTH_NOT_CERTIFIED.

## One smallest next patch and independent benchmark — NOT_RUN

Limit the next work item to a **named optional numerical material-response
block for a fixed active/plastic branch**, preserving the actual48-iteration
`material.update`, all physical parameters/forces/heat, same-history input,
and existing acceptance gates. Do not hardcode PID269 or import the unrelated
27-point fixed-normal HHT contact correction. No new backend, full Jacobian,
matrix build or342 evaluation is authorized by this review.

For the scalar return closure with qtr>0, rate r and fixed history ep0,

```text
qtr-3G*dep-Y(ep0+dep,r,T)=0
ddep=(dqtr-Y_r*dr-Y_T*dT)/(3G+Y_ep)
dq=dqtr-3G*ddep
dS=(q/qtr)*dStrial + Strial*(dq/qtr-q*dqtr/qtr^2).
```

For the FKS branch above its rate floor,
`Y_r=(a1+a2*atan(a3*ep))*theta*b1*Tstar^b2/r`;
`Y_ep=a2*a3*lambda*theta/(1+(a3*ep)^2)`.
Below the floor, the fixed-branch Y_r is zero. At **r=rate_floor**, a classical
two-sided derivative is not certified: do not difference across the switch
and call it a Jacobian. Active/inactive yield switching requires the same
restriction. A rate-independent radial tangent alone would miss the specific
supported term. Keep the denominator **3G+Y_ep**, not3G.

One independent material-point benchmark package should isolate three
orthogonal scalar-return input directions: qtr only, rate only, T only.
Use a simple affine-in-ep/r/T yield function as a closed-form oracle for all
three derivatives, followed by the literal FKS expression and a separately
implemented scalar reference root/implicit derivative. Reference generation
must not call the tested tangent. Fix all amplitudes and numerical comparison
bounds before execution, including the48-bisection root interval uncertainty.
The physical `update` remains the same48-iteration algorithm: no replacement
root, altered history, reduced iterations, yield clipping or residual heat.

Then compose the derivative with a zero-spin deviatoric material-point input
at nonzero S0, rate safely above1/s and positive plastic return, always reusing
the exact same S0/ep0 for each directional evaluation. Compare qtr/r/T
derivatives independently before their composed strain-rate response. Add
fixed-branch below-floor and elastic controls plus an explicit at-floor
NOT_DIFFERENTIABLE/rejection expectation; no cross-branch central difference.
The force-coupling design must later retain the midpoint stress factor1/2,
actual physical mass and pressure/continuity convention. Until this local
benchmark is independently passed, do not wire the proposed block into a
new342 solve or claim action improvement. All next-patch work is NOT_RUN.

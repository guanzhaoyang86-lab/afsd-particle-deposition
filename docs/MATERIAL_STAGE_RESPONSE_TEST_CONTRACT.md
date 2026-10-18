# Independent material-stage search-block tests

2026-10-04, registered before the first test run. Scope is the pure-array
`MaterialStageResponse` constructor and `apply` only. Do not construct the
342-particle matrix, call `from_journal`, NNPS, `material.update`, a physical
trial, root, or checkpoint. The separate material-response component tests
belong to another suite. CapturingTrial/refresh integration is not implemented
or tested here and remains NOT_RUN.

Reference answers below are literal fractions and/or the direct solve of an
independently assembled **uneliminated** material/pressure/density/temperature
block. Never use the tested Jv/Jr/JT/A or correction to construct that reference.
Comparison bounds: rtol=2e-12, atol=2e-13 for O(1) manufactured SI quantities,
as in the old structural/continuity component contract. Returned backward
error must meet the original1e-12 bound. Input copies, zero corrections,
authority flags and detached arrays are exact. Unit transformations are
manufactured dimensional tests, not scaling an actual process problem.

Scalar reference1: m4,V2,K3,rho_ref5,rho_parent5,h2,D(Lxx,vx)=1,
C=0,ST=0,divL=0. dbar=.5,drho2 -> delta_bar=.32,delta_rho=-1.2.
Pressure-only mean equation is `4 delta_bar-.6 delta_rho=2` and continuity
is `10 delta_bar+delta_rho=2`. This checks its physical signs and midpoint1/2.

Scalar reference2: m2,V1,K2,rho_ref3,rho_parent3,h=.5,divL=.8,
D(Lxx,vx)=1,C_xxxx=4,ST_xx=-4. Literal Jv=-2,Jr=1/3,JT=1,
kappa=-25/24. dbar=2/5,drho=3/5,dT=1/5 gives delta_bar=1296/3725,
delta_rho=177/745,delta_T=1/5. The independent full block is
`[[5/2,-1/12,-1/4],[25/24,1,0],[0,0,1]]` times (bar,rho,T)
= (4/5,3/5,1/5). It catches the exact rational continuity coefficient and
the stress-midpoint1/2 times temperature-midpoint1/2 factor.

Two-particle fully coupled reference: m=(2,3),V=(1,1),K2,rho_ref1,
rho_parent=(.5,.5),h2,div=(0,0). D at particle0 has Lxx=v0x and
Lyy=-v1x, all other entries zero. Its material tangent rows/columns xx,yy
are `[[2,2],[0,4]]`; zz output is minus their row sum, preserving deviatoric
stress output. ST=(Sxx,Syy,Szz)=(-4,-4,8); particle1 C/ST are zero.
Selected physical velocity DOFs therefore have independent literal
Jv=[[-1,1],[0,-2]], Jr=JT=(1,-1)^T, kappa*E=(-1,1).
For dbar=(1,2),drho=(.5,0),dT=(.25,0), expected selected delta_bar
=(27/22,95/88), delta_rho0=31/88; the other fields are direct/uncoupled.
The uneliminated selected full block is

```text
[[ 3,-1,-1,-1],
 [ 0, 5, 1, 1],
 [ 1,-1, 1, 0],
 [ 0, 0, 0, 1]] * [bar0,bar1,rho0,T0] = [2,6,.5,.25]
```

Prescribing the second selected increment to2 yields bar0=27/16,
rho0=13/16, second algebraic reaction81/16. All constrained columns must
remain in the density back-substitution. This is an explicit manufactured
gradient/constitutive block, not a standalone physically supported SPH cloud.

Tests also cover: nonsymmetry retained; empty free block; zero residual and
independent event-average velocity; detached immutable inputs and output
metadata; density/temperature proposals not clipped; direct five-field lifts;
permutation/ID-only changes; mass/length/time/temperature unit covariance;
invalid shapes/values/domains and an exactly singular free matrix. A patched
`numpy.linalg.pinv` must never be called. Algebraic reactions are never heat.
No test asserts that this incomplete response is a physical contact tangent.

One test invocation, external60 s maximum (no actual trajectory). Preserve
the initial source hash/stdout/stderr/exit. If it fails, report to the source
owner; no unrecorded source or tolerance repair by this reviewer. Status before
first execution: NOT_RUN.

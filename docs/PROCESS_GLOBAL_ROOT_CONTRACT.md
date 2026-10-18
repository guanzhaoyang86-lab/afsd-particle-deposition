# Joint process-root experiment — registered before execution, 2026-10-03

This is a named ENGINEERING_APPROXIMATE numerical integration experiment, not
production or deposition validation. Parent is the published stage_004 at
850.393632910904 ns. No physical parameters, material history, contact kick,
geometry, displacement/temperature/density guards or raw five-field tolerances
are changed. The new thermal map is explicitly backward Euler; it is not the
old explicit thermal discretization. Its acceptance and limitations are in
PROCESS_GLOBAL_THERMAL_VARIANT.md.

## Numerical proposal only

ContinuityStagePreconditioner transforms map-minus-candidate defects through the frozen
elastic constrained response. Type-II Anderson operates on these *proposals*,
with beta=.75, history<=6 and explicit QR condition<=1e10 as in the old solver.
The equations and accepted root are still the actual five-field trial map.
Independent event-average velocity remains a field. Preconditioner reactions
are algebraic diagnostics and never substituted for physical contact reactions.
The first direct-rho StagePreconditioner is retained as a tested limited variant;
the actual experiment uses the continuity Schur RHS and density back-substitution
from PROCESS_GLOBAL_PRECONDITIONER.md / the continuity supplement. It remains an
approximate tangent; no claim of an exact FKS/contact Jacobian is made.

Convergence and factor-two best-residual growth safeguard use **original raw
map residuals**, never the preconditioned correction norm. x:1e-12 m,
v/vbar:1e-9 m/s, rho relative:1e-10, T:5e-12 K. These are inherited algebraic
tolerances, not claims of process accuracy. The full mechanical/thermal
acceptance gate still follows root convergence. Every evaluated rejected
candidate remains on disk; root solver commits no state. Domain/contact
exceptions stop this registered experiment; no hidden line search or dt change.

## Independent root tests (run before AFSD)

Literal scalar root u=(-2u+.002) => 1/1500 m/s; coupled five-field literal
roots; a frozen harmonic midpoint problem m=1 kg,k=4 N/m,dt=1 s,
x0=.002 m,v0=.001 m/s => x1=.0005 m,v1=-.004 m/s,vbar=-.0015 m/s.
An intentionally zero or tiny proposal cannot make a nonzero raw residual PASS.
Retain event-average != endpoint-average manufactured reference, raw growth
rejection, domain/NaN failures, exact evaluation budget and immutable input.
All comparisons use the raw tolerances above (oscillator x <=1e-12 m).

## One actual registered step, conditional on cost

Choose dt=500 ns solely from independently read-back scale audit
PROCESS_GLOBAL_STEP_AUDIT.md, not from trial-and-error matching a desired result.
Maximum physical steps=1; maximum complete trial evaluations=32. First measure
one uncommitted same-source trial and preconditioner build/apply, inner80 s /
outer90 s. Launch only when 1.5*(build+32*measured evaluation-and-proposal)+30 s
<=450 s. Actual solve inner540 s / outer570 s, always below600 s.
Do not reduce dt and repeat after failure within this experiment. Publish
REJECTED candidate and last accepted reference if any gate fails.

500 ns is NOT a resolved elastic wave step: the highest frozen elastic omega*dt
is approximately10 and midpoint phase error is large. A converged step would
prove larger-step numerical solvability only, not transport, process accuracy,
or permission to run a long process. Keep inertia and all moduli unchanged.

Actual results and remaining project-completion work must be reported separately.

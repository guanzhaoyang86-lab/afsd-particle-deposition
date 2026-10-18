# Cross-zero contact response: continuous limits versus finite quadrature

2026-10-03. Independent mathematical investigation; no particle trajectory,
material update, solver advance, checkpoint publication or production run.
This report supplements, without editing, the frozen
`PROCESS_STOP_CAUSAL_REVIEW.md`: its proposed complete contact-response
investigation must distinguish a finite-node numerical preconditioner from
a classical derivative of the continuous cross-zero impulse.

## Decision and evidence

**CONFIRMED / NOT_DIFFERENTIABLE:** the PRIVATE folded cross-zero contact
impulse generally has a transverse absolute-value cusp. Smooth force
Jacobians at individual nonzero Gauss nodes cannot be integrated and called
the classical Jacobian of the complete continuous impulse.

**PASS, narrowly scoped:** all 130 preregistered pure-math checks passed in
one run, internal 2.626212 s / externally supervised 2.763430 s, exit 0.
All 4516 frozen entries matched before and after. No repair or second run.
The limits, independent direct force quadrature and finite-rule FD evidence
are in `logs/process_response_20261003/limits_initial/report.json`.

**CONDITIONALLY_ELIGIBLE:** a separately named fixed-Gauss approximate
preconditioner can still be tested as a search-direction diagnostic. That
does not require a continuous Jacobian to exist. Its final residual and all
old domain/heat/zero/collinearity/publication guards must remain unchanged.
A saved-state experiment can test numerical utility, not certify a physical
stop, continuous differentiability, or production readiness.

Actual invocation from the project directory:

```powershell
.envs\pysph\python.exe -B logs/process_response_20261003/run_limits_probe.py
```

The new supervisor ran only:

```text
.envs/pysph/python.exe -B scripts/probe_contact_response_limits.py --output logs/process_response_20261003/limits_initial
```

Preregistration: `docs/PROCESS_RESPONSE_LIMITS_CONTRACT.md`. Full execution
record: `logs/process_response_20261003/limits_initial_command.log` and
`limits_initial_supervision.json`. The script imports only the Python
standard library and NumPy; no tested contact function generated a reference.
The old supervisor and frozen solver sources were not changed.

## Closed-form counterexample

Let s be dimensionless in [0,1], physical dt = h ds, U > 0 a velocity scale,
C > 0 the constant Coulomb cap, and s0 an interior zero of the unperturbed
path. Use orthogonal tangential coordinates and a fixed normal:

```text
q = s - s0,  a = s0,  b = 1 - s0,  rho = epsilon/U
u(s) = (U*q, epsilon)
Ft(s) = -C*sign(q)*u/|u|.
```

This is an explicitly declared mathematical PRIVATE extension, not an
accepted physical crossing. Choosing the opposite incoming direction/sign
convention multiplies every force and impulse below by -1; the existence or
nonexistence of a derivative is unchanged. In particular, the actual helper
continuation with incoming slip along negative x uses that opposite overall
sign. The proof does not silently identify the folded law as ordinary sliding.

Direct integration on each side of q = 0 gives, for a,b > 0:

```text
Jx = -C*h*[sqrt(a*a+rho*rho)+sqrt(b*b+rho*rho)-2*|rho|]
Jy = -C*h*rho*[asinh(b/|rho|)-asinh(a/|rho|)]
J(0) = (-C*h, 0).
```

The two one-sided expansions are:

```text
Jx-Jx(0) = 2*C*h*|rho| - C*h*rho^2*(1/a+1/b)/2 + O(rho^4)
Jy       = -C*h*rho*log(b/a) + O(rho^3).
```

Consequently dJx/depsilon has right limit +2Ch/U and left limit -2Ch/U.
The full vector impulse has no classical transverse derivative at epsilon
zero. Jy alone has derivative -(Ch/U)log(b/a), but differentiating its
integrand first produces -C/[U*q]. The integrals on either side separately
diverge. Only a paired Cauchy principal value supplies log(b/a); this is not
an ordinary integrable derivative justified by dominated convergence.

For a = b the Jy cancellation is exact and can hide the issue, while the Jx
cusp remains. A centered finite difference of Jx is identically zero because
Jx is even. Agreement of that centered difference with zero is therefore
not a differentiability test. Both signed one-sided quotients are necessary.

At s0 = 0:

```text
Jx = -C*h*[sqrt(1+rho*rho)-|rho|]
Jy = -C*h*rho*asinh(1/|rho|).
```

The parallel one-sided slopes are +Ch/U and -Ch/U. Jy/epsilon tends to
-(Ch/U)log(2U/|epsilon|), which diverges; s0 = 1 reverses this transverse
logarithm's sign. An endpoint-zero continuation is not generically rescued
by moving the interior singularity to an endpoint.

For comparison only, the ordinary, nonfolded Coulomb law has

```text
Jx = -C*h*[sqrt(b*b+rho*rho)-sqrt(a*a+rho*rho)]
Jy = -C*h*rho*[asinh(b/|rho|)+asinh(a/|rho|)].
```

Its interior transverse difference quotient diverges logarithmically. This
mathematical comparison is not authorization to cross a physical stop
without locating the event and determining its admissible right limit.

## The obstruction concerns the complete impulse correction

The frozen contact helper returns
`DeltaJ = Jactual - h*(Fcontact0+Fcontact1)/2`, not endpoint force alone.
At an interior zero both endpoints of this counterexample remain nonzero.
Their force dependence on epsilon is smooth, so subtracting their linear
endpoint impulse cannot cancel the continuous |epsilon| cusp. The same
obstruction survives in a residual containing the complete DeltaJ.

The HHT dense path has, with the left state held fixed,

```text
d = beta - gamma/2
delta_v(s) = h*[gamma*s+6*d*s*(1-s)]*delta_A1.
```

Thus a generic transverse acceleration perturbation gives nonzero
delta_u_perp at an interior zero. The constant-epsilon counterexample
identifies the same local singular response, not merely an unrelated endpoint
force defect. It does not prove every restricted event-preserving direction
is nondifferentiable, nor diagnose the old trial005 failure uniquely.

## What actually passed numerically

The contract fixed C = 1 N, h = 1 s, U = 1 m/s before execution. Direct
quadrature integrates the actual force independently using
q = +/-|rho|sinh(z), separately on each side, including its Jacobian.
Orders 96/192 agree with the closed form and each other within the registered
2e-12 Ns per component. The analytic stable difference quotient uses
`rho^2/(hypot(a,rho)+a)` to avoid subtracting near-equal unit impulses.
The absolute value here describes the proven cusp, not clipping of heat:
this probe computes no friction heat at all.

For s0 = 0.4, the independent target transverse slope is
-log(1.5) = -0.4054651081081644 kg. Selected recorded values:

|epsilon (m/s)|Right parallel quotient (kg)|Left parallel quotient (kg)|Jy/epsilon (kg)|
|---:|---:|---:|---:|
|1e-4|1.9997916666691984|-1.9997916666691984|-0.4054650994276088|
|1e-8|1.9999999791666665|-1.9999999791666665|-0.4054651081081638|
|1e-10|1.9999999997916669|-1.9999999997916669|-0.4054651081081637|
|1e-14|1.9999999999999791|-1.9999999999999791|-0.4054651081081602|

All centered parallel quotients are zero. At the endpoint s0 = 0, Jy/epsilon
is -19.11382792451231, -23.7189981105004, -28.324168296488494 and
-32.92933848247659 kg for epsilon 1e-8, 1e-10, 1e-12 and 1e-14 m/s,
respectively, confirming the registered logarithmic asymptote.

## Fixed finite quadrature is a different map

For a fixed finite set of nodes avoiding q = 0, the sum is smooth near
epsilon zero. Its derivative in this example is

```text
dJ_N/depsilon = (0, -(C*h/U)*sum_k weights_k/(s_k-s0)).
```

An independent FD of that finite sum correctly agrees with this expression.
It need not agree with a continuous derivative; the latter vector derivative
does not exist. For s0 = 0.4 all registered finite-map FD tests passed:

|Gauss order|Exact discrete transverse derivative (kg)|Centered FD (kg)|
|---:|---:|---:|
|8|-22.11322040948376|-22.113220300017087|
|16|0.1738594671213921|0.17385946128418164|
|32|0.4360944584643667|0.4360944507173508|
|64|1.0035580103462038|1.0035579991738317|
|128|2.4743448206127727|2.474344802417295|

These are not ordinary quadrature convergence to the PV value -0.405465.
At epsilon = 1e-10 m/s the order-32 parallel one-sided quotient is
1.0961529243191383e-8 kg, versus the continuous quotient
1.9999999997916669 kg. Refinement and epsilon-to-zero limits do not commute.
Symmetric s0 = 0.5 gives near-zero discrete transverse derivatives, again
without eliminating the continuous parallel cusp.

The actual solver additionally chooses among 8/16/32-node rules. A derivative
of a named fixed-32 map is not automatically the exact derivative of that
adaptive map. Rule changes and roundoff endpoint overrides need their own
scope statements even where an individual finite sum is smooth.

## Legal finite differences and permitted next experiment

The frozen `hht_stop_contact.py` first restricts PRIVATE target paths to
naturally collinear polynomial coefficients, with bound
`512*eps_machine*max(1, velocity_scale, coefficient_scale)`. It checks the
full vector; no velocity projection is allowed. Its roundoff-zero force
override is also part of the implemented map.

- Longitudinal perturbations that retain that collinearity, gap and branch
  domain may legitimately test the existing map. An exactly collinear
  folded constant-cap impulse is smooth/constant in many such directions.
- Generic transverse perturbations exceeding that bound are outside the old
  helper's admitted domain. The independent counterexample investigates the
  mathematical extension, not permission to bypass this guard.
- Perturbations within the roundoff tube can test finite-node/override
  arithmetic, subject to cancellation and quadrature error. They do not
  establish a continuous transverse derivative or a physical crossing.
- An event-preserving perturbation with transverse velocity zero at the
  same event can have a regular restricted derivative. For example
  `u = U*q*ex + epsilon*s*q*ey` makes the folded integrand
  `-C*(ex+rho*s*ey)/sqrt(1+rho^2*s^2)` and has derivative
  `(0,-Ch/(2U))` at zero. This is a distinct event-manifold formulation;
  finite such paths are generally outside the current fixed-direction
  collinearity domain. It is a design observation, NOT_RUN in this probe.

A safe next numerical diagnostic can be separately registered as
**fixed-Gauss approximate contact-response preconditioning**, not continuous
Newton linearization. Eligibility requires independent finite-rule response
tests including the full endpoint-linear correction, this retained cusp
counterexample, source/artifact gates and an unchanged actual residual.
Use the exact saved trial005 parent and duration, at most its old 24 residual
evaluations and a single supervised sub-60-second attempt. Do not change
adaptive quadrature acceptance, collinearity/full-vector zero guards, raw
heat, friction law, scaling, state acceptance or publication. Outputs stay
PRIVATE even if the root improves. Preserve the old failure.

Those finite-rule FD tests can establish algebraic consistency of the new
search direction. They cannot establish that the old adaptive map has that
derivative, the new search stays within the old domain, a continuous zero
event exists, two almost-simultaneous events share a joint root, or accepted
post-event physics converges. Failure on any of these remains evidence.
No conclusion here implies the proposed preconditioner must work.

## Status and provenance

|Item|Status|
|---|---|
|Closed-form / independent direct quadrature / limits / finite-rule FD|PASS: 130 checks|
|Classical full-vector continuous cross-zero Jacobian|NOT_DIFFERENTIABLE|
|Claim of a classical continuous cross-zero Newton fix|BLOCKED|
|Separately named finite-rule numerical preconditioner|CONDITIONALLY_ELIGIBLE, subject to its own gates|
|Saved trial005 diagnostic in this subtask|NOT_RUN|
|Event-manifold restricted derivative implementation|NOT_RUN|
|Particle advancement / new accepted state / actual AFSD progress|NOT_RUN / none / none|

The limits script SHA256 is
`853c351837654ede99489b699bcdbcb10f04a5e44c2c235847f78887b32509ae`;
preregistered contract SHA256 is
`d92980fee17d0e291d70e0d6a9bed7c75afd39c1c4914f886750b8156a53abfb`.
The frozen manifest SHA256 checked both times is
`71f37691671eab2afb3a35d198a688df9d1b2780d93f417df2b4b31aaa77f13c`.
The supervisor receipt binds its own source and command-log hashes. This
report does not upgrade the project to validated AFSD reproduction.

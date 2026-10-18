# Contact-path acceleration response: independent component contract

Registered 2026-10-03 before execution of the new tests. This package changes
no physical law, event tolerance, frozen source, or accepted state. No scout,
saved-trial solve, ODE trajectory, or AFSD process is run by these tests.

## Scope and independent reference

The new `contact_path_response` helper is a pure local matrix calculation.
Its point inputs have SI units: unit normal n, tangential nonzero slip u in
m/s, Fn>0 in N, kn>0 in N/m, mu>=0, and fixed branch sign sigma=+1 or -1.
The physical force is `F=Fn*n-sigma*mu*Fn*u/|u|`; sigma=-1 is only a PRIVATE
smooth branch extension, never physical passivity or publication authority.
Position-independent prescribed wall velocity and fixed normal are assumed.

Tests construct force values independently from that literal expression and
`Fn(x)=Fn0-kn*n.dot(x-x0)`. No tested helper, contact integrator, path sampler,
kinetic identity, or force routine supplies the expected value. The literal
point case n=(0,0,1), u=(3,4,0), Fn=5, kn=7, mu=.2 gives

```text
Kx = [[0,0,.84*sigma],[0,0,1.12*sigma],[0,0,-7]] N/m
Kv = sigma*[[-.128,.096,0],[.096,-.072,0],[0,0,0]] N s/m.
```

Independent central differences include normal-gap and transverse-direction
perturbations, both branch signs, a non-axis-aligned normal, and a proper 3D
rotation. The zero-friction limit retains the normal penalty response.

For the complete path test, literal dense HHT polynomials are evaluated from
fixed x0,v0,A0 and trial A1. With force weight a, gamma=1.5-a,
beta=(2-a)^2/4, d=beta-gamma/2, h>0:

```text
dv=h*((1-gamma)*A0+gamma*A1)
v(s)=v0+s*dv+6*h*d*(A1-A0)*s*(1-s)
x(s)=x0+h*(s*v0+.5*s*s*dv+h*d*(A1-A0)*(3*s*s-2*s**3))
Jactual=h*integral F(x(s),v(s),s) ds
DeltaJ=Jactual-h*(F0+F1)/2
Rcontact=-a*F1-(1-a)*F0-DeltaJ/(gamma*h)
```

The independent reference central-differences these complete force values,
not just the analytic derivative formula. Gauss64 computes the reference
impulse, Gauss128 checks its quadrature stability, and the tested matrix uses
independently generated Gauss16 rows. The three returned matrices must match
`D_A1 Jactual` (kg s), `D_A1 F1` (kg), and `D_A1 Rcontact` (kg), respectively.
The initial API proposal's impulse-response unit was corrected before any
execution: its key is `actual_impulse_acceleration_response_kg_s`, not kg.
The matrix identity is

`D_A1 Rcontact=(1/(2*gamma)-a)*D_A1 F1-D_A1 Jactual/(gamma*h)`.

Manufactured cases use h=.2 and .7 s; a=.9,1,2/3; both signs; and a rotating
nonzero tangential vector with a changing positive normal gap. All finite-
difference perturbations keep gap and full slip strictly nonzero and hold
the sign fixed. PRIVATE transverse derivatives here describe only a frozen
smooth extension; they do not relax the existing natural-collinearity guard
or prove a derivative across its admissibility boundary. No path crosses zero.
The helper must keep `continuous_path_differentiability_certified=False`.

## Predeclared tolerances and falsification checks

- Literal matrices, proper-rotation covariance, and exact matrix identities:
  absolute 2e-12 in the compared SI unit plus relative 2e-12.
- Point and complete-residual central differences: absolute 2e-8 in the
  compared SI unit plus relative 2e-6, componentwise. Fixed perturbation
  magnitudes are 2e-5 and 1e-5 in m, m/s, or m/s^2, as appropriate. Both
  estimates must satisfy the same bound; no step-size tuning after a failure.
- Gauss64 versus Gauss128 force impulse: absolute 1e-12 N s plus relative
  1e-10. Raw maximum errors are printed to the preserved initial log.
- A deliberately missing DeltaJ endpoint correction must differ from the
  correct residual derivative by more than 1e-3 kg in at least one tested
  component. A deliberately physical sign substituted for PRIVATE must
  differ by more than 1e-3 in at least one matrix component. These are
  sensitivity controls, not modifications to the production implementation.
- Invalid zero/full-normal slip, bad unit normal, nonfinite data, nonpositive
  Fn/kn/h, negative mu, invalid/bool sign, empty/bad-weight nodes, and changed
  normal/kn/mu along a declared fixed path must raise ValueError. A strictly
  positive small slip is not clipped or regularized.
- Caller arrays/rows remain unchanged, and returned matrices are detached.

One initial supervised component invocation has an external 60 s cap. Keep
its stdout/stderr and source/test/contract hashes even if FAIL. Any repair is
evidence-driven, limited to two, and must preserve initial evidence. Until
execution, all checks above are NOT_RUN. Success would establish only this
bounded numerical component, not solve trial005, locate a stop, or validate
the coupled trajectory.

## Optional two-component force-driven stop reference — NOT_RUN

For m=1 kg, Coulomb capacity C=2 N, other force (1,0) N, and initial slip
(1,1) m/s, write u=r(cos(theta),sin(theta)), q=tan(theta/2),
q0=sqrt(2)-1 and B=(3+2*sqrt(2))/2 in the corresponding SI units. Independent
polar equations `r_dot=cos(theta)-2`, `theta_dot=-sin(theta)/r` give
`r=B*q*(1+q*q)` and `dt/dq=-B*(1+q*q)`. Thus

```text
tstop=B*(q0+q0**3/3)=(2*sqrt(2)+1)/3 s
Q=2*integral r dt=B**2*(q0**2+q0**4+q0**6/3)
 =(5/4+sqrt(2)/3) J
DeltaK=-1 J; Wother=Q-1 J.
```

The required static traction has magnitude 1 N<C, so sticking is feasible
after the full-vector stop. These are derived literal references, not an
executed dynamic test, and no PASS is assigned to the optional benchmark.

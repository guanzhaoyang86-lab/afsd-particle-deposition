# Independent folded-contact response limits: preregistration

2026-10-03. Before numerical execution. Pure mathematical path integrals;
no particle trajectory, backend contact/material call, accepted state,
checkpoint or VTK. Preserve every4516entry of PROCESS_STOP_FREEZE exactly.
Output only under new logs/process_response_20261003 subdirectories.

## Fixed mathematical problem and independent limits

s is dimensionless in[0,1], physical dt=h ds. Fixed normal is perpendicular
to the displayed two-component tangent plane. C=1N, h=1s, U=1m/s,
u(s)=(U*(s-s0),epsilon), rho=epsilon/U, a=s0,b=1-s0.
The declared PRIVATE folded mathematical extension is

    Ft=-C*sign(s-s0)*u/|u|.

This overall sign convention is explicit; using the opposite incoming
convention multiplies all impulses/derivatives by-1 without restoring
differentiability. It is not a physical accepted cross-zero trajectory.

For a,b>0, independent integration gives

    Jx=-C*h*(hypot(a,rho)+hypot(b,rho)-2*abs(rho))
    Jy=-C*h*rho*(asinh(b/abs(rho))-asinh(a/abs(rho))).

At epsilon=0: J=(-C*h,0). At s0=.4, right/left epsilon derivatives
of Jx are +2/-2 in SI; derivative of Jy is -log(1.5). The complete vector
impulse has NO classical transverse derivative. Central Jx difference is
exactly0 and is not evidence against this cusp. Subtracting the smooth
endpoint-linear contact impulse cannot remove an interior cusp.

At s0=0, Jx has one-sided slopes+1/-1 and Jy/epsilon=-asinh(1/abs(rho))
diverges logarithmically. At s0=1, the transverse quotient has opposite sign.
Also inspect s0=.5 to expose accidental symmetric principal-value cancellation.

## Numerical checks fixed before execution

The new standalone scripts/probe_contact_response_limits.py imports only
Python standard library and NumPy. A separate direct force quadrature uses
q=+/-abs(rho)*sinh(z) on each side and the actual transformation Jacobian;
it does not call the closed-form routine to construct reference forces.
Orders96/192 and epsilon=+/-(.2,.02,.002,.0002,1e-8), s0=.4,.5,0,1.
Absolute closed-form/quadrature agreement and quadrature96/192 difference:
2e-12Ns per component. These are unit-scale math checks, not AFSD gates.

Stable one-sided difference quotients use
rho^2/(hypot(a,rho)+a), avoiding subtracting two nearly equal unit impulses.
Use |epsilon|=1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14.
For |epsilon|<=1e-8, interior parallel one-sided quotients must be within
3e-8 of+2/-2 and transverse within1e-11 of-log(1.5). Endpoint parallel
quotient within1e-8 of+1/-1; logarithmic asymptote within1e-11.

Separately inspect the FINITE Gauss map for orders8,16,32,64,128. Its
node derivative at epsilon0 is(0,-sum(weights/(s-s0))). A centered FD
with epsilon=1e-4*minimum_node_distance must match that finite-map derivative
within3e-7*max(1,abs(component)), not the continuous-integral derivative.
At order32 and epsilon1e-10, its stable one-sided parallel quotient must be
below1e-6 while the continuous quotient exceeds1.99: this is an intentionally
demonstrated distinction, not a failing physics tolerance to relax.

Record the ordinary, nonfolded law's analytic transverse logarithmic response
as a separate mathematical comparison, not permission to traverse a physical
stop without a located event. No reference CSV is invented.

## Integrity, cost and result semantics

One numerical probe, internal55s/external60s. Expected pure NumPy math and
4516-file integrity checks under10s; abort on any frozen-file mismatch.
No installs, GPU, parallel trajectory, production, or saved-trial solve.
Record command, source/contract hashes, every check and raw observed value.
PASS means the mathematical limit evidence agrees; the relevant continuous
vector Jacobian is still NOT_DIFFERENTIABLE, and using it as a classical
cross-zero Newton correction is BLOCKED. A finite-node derivative is not
silently promoted to a continuous one. Unperformed checks remain NOT_RUN.

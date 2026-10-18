# Fixed-branch FKS numerical response and bounded integration

2026-10-04. This is the next minimal patch after PROCESS_BRIDGE_ACTION_CAUSAL_REVIEW,
not a replacement material model or acceptance relaxation. The user permits a
smaller explicitly approximate deliverable if the complete route is impractical;
this does not permit fabricated transport, trajectories, clipping or scaling.

## Change and stop boundaries

Add a pure numerical response module; keep material.py, contact, spatial operators,
accepted histories, masses/moduli/times and thermal sources unchanged. Derive the
same-branch return derivative, including strain hardening, FKS rate and temperature
partials. Composition may include the exact Rodrigues trial-stress directional
derivative. It returns a derivative/search ingredient, never a new physical state.
No stored particle ID or tracer controls the derivative; all points use the same
law. Rate floor and elastic/plastic switching are not smoothed to obtain a Jacobian.

Independent component tests have their own preregistered contract. Preserve first
FAILs and source snapshots; at most two bounded implementation/test-tool repair
rounds before recording the reason and switching to delivery of the verified part.
Do not keep adjusting tolerances or test amplitudes until green. Old failure
counterexamples remain. A component PASS does not certify AFSD material physics.

## Analytic definitions and provenance

The source of the derivative is the actual unchanged material.py equations, not
another constitutive assumption. At an already solved positive-plastic base,
qtr-3G*dep-Y(ep0+dep,r,T)=0. Thus ddep=(dqtr-Yr*dr-YT*dT)/(3G+Yep),
dq=dqtr-3G*ddep, dS=(q/qtr)dStr+Str*(dq/qtr-q*dqtr/qtr^2).
Preserve the actual 48-bisection physical return; do not replace it with the
implicit derivative. Its finite root bracket uncertainty enters reference bounds.

The temperature factor remains 1-(1+exp(-c1*Tstar))^(-1/c2), and b1 retains the
current explicitly approximate dimensionless interpretation. No PDF parameter
or disputed paper unit is changed in this package. FKS derivatives use the end ep;
below the existing rate floor Yr=0; at the floor reject classical differentiation.
Zero qtr, the yield switch, invalid/nonfinite inputs and inconsistent declared
active branches fail explicitly. The base return residual is reported as a
diagnostic, not clipped, solved away or converted to heat.

The Rodrigues derivative uses analytic coefficient derivatives, with Taylor
evaluation for angle squared<1e-6 to avoid cancellation. Terms are retained through
cubic order in coefficient derivatives, so the first omitted term is far below
double-precision error at that switch. This changes no actual rotation/trajectory.
At zero deviatoric rate, a nonzero deviatoric direction lacks a two-sided norm
derivative and rejects. Pure spin direction with zero deviatoric part is allowed.

## Important counterevidence, not a promised universal fix

The prior 50ns/650ns bridge has one plastic point, so it is the first integration
target. The old 500ns/850ns last failed candidate instead has no new plasticity;
its stagnation is concentrated at stopping contact particles. Its failure is not
explained merely by the plastic omission at a different parent/time. Preserve this
distinction in all reports. Do not launch another265s root run on that assumption.
Read MATERIAL_RESPONSE_INTEGRATION_DESIGN for journal binding and matrix factors.

After component PASS, a stored-data response/action check may reuse the real
50ns base journal, before any new physical call. A future actual diagnostic needs
its own preregistered candidate count/cost and<=600s watchdog. No three-condition
production, GPU work, long paid task or silent quasistatic switch is included.

## Status on registration

Component source newly implemented; independent tests NOT_RUN. Integration into
the real five-field root, macrostep improvement and transport NOT_RUN. The last
real AFSD accepted state remains850.393632910904ns with zero crossings. The goal's
automatic runner currently reports usageLimited; ordinary work performed in this
turn is saved, but unattended overnight continuation is not guaranteed.

# Independent fixed-branch material-response component tests

2026-10-04, registered before execution. This package tests a new numerical
response helper, not a replacement constitutive update. The old48-iteration
`material.update`, all histories, FKS parameters and physical gates remain
unchanged. No342-particle evaluation, matrix assembly or trajectory is run.

## API and mathematical reference

`fks_partials(law,ep,rate,T)` receives the **returned** ep=ep0+plastic_increment.
`radial_direction(Strial,plastic_increment,dStrial,dr,dT,G,partials,active)`
receives the already solved base increment, not its directional derivative.
It returns dS in Pa, dplastic_strain dimensionless, and dq in Pa for a
dimensionless directional parameter. Scalar and batched cases must agree.

The independent closure is

```text
qtrial-3G*plastic_increment-Y(ep0+plastic_increment,r,T)=0
dplastic_increment=(dqtrial-Y_r*dr-Y_T*dT)/(3G+Y_ep)
dq=dqtrial-3G*dplastic_increment
dS=(q/qtrial)*dStrial+Strial*(dq/qtrial-q*dqtrial/qtrial^2).
```

No tested helper generates expected values. The affine oracle has
G=7Pa, ep0=.2, qtrial=14Pa, r=1.5/s, T=300K and
`Y=5+3*ep+2*r-.25*(T-300)` Pa. Its closed return has dep=.225 and q=9.275Pa.
Separate unit qtrial/r/T directions have literal `(ddep,dq)` values
`(1/24,1/8)`, `(-1/12,7/4)`, `(1/96,-7/32)`, respectively. This catches
dropping Y_ep from the denominator or losing rate/temperature signs. A
nonproportional shear direction checks the transverse stress response, not
only the equivalent norm.

For literal FKS tests, the reference directly evaluates the documented
atan/log/exp expression using standard arithmetic, and independently solves
the scalar closure by100 bisection iterations. It does not call the tested
partials, tangent, or `FKS.yield_stress` to construct reference outputs.
Independent qtrial-only, rate-only, temperature-only and nonproportional
directions use the same ep0 and hold all other independent inputs fixed.

Above the floor, rate differentiation is local to that branch. Below the
floor the rate partial is exactly0. At rate==floor the helper must reject
with NonDifferentiableMaterial; no centered difference straddles the switch.
An inactive radial branch requires zero base dep and strict qtrial<Y; exact
yield equality and the equivalent-norm origin qtrial=0 reject. Positive
active base dep, positive G/qtrial/returned q and finite data are required.

## Fixed numerical amplitudes and bounds

- Literal affine directions, scalar/batch equivalence and elastic controls:
  absolute2e-12, relative2e-12, in each returned SI unit.
- Literal FKS values: absolute1e-6Pa plus relative2e-12.
- FKS partial finite differences use fixed steps ep1e-7, rate1e-5/s,
  temperature1e-3K, and half each step. Compare with absolute.2 in the
  corresponding derivative SI unit plus relative2e-6. No step-size fitting.
- Radial-map and unchanged48-iteration-update directional differences use
  fixed dimensionless steps1e-5 and5e-6. Stress/equivalent derivative bounds
  are absolute.05Pa plus relative2e-5; plastic-increment derivative bound is
  absolute2e-12 plus relative2e-5. Direction magnitudes are recorded in the
  tests: trial-stress1e5Pa, rate.2/s, temperature10K, or a declared mixed tensor.
- Independent root versus original base update: stress absolute1e-6Pa plus
  relative2e-12; returned ep absolute2e-15, relative0.

For qtrial approximately5MPa and G=26.3GPa, the physical48-step bisection
increment interval is at most `qtrial/(3G)*2^-48`; its induced stress scale is
about2e-8Pa before division by the difference step. The fixed derivative
absolute allowance includes this uncertainty and subtracting the unchanged
ep0, rather than silently changing the original root iterations.

## Same-history composition and optional trial-response checks

The original `update` comparison uses fixed nonzero
`S0=diag(-5e6/3,-5e6/3,10e6/3)` Pa, ep0=.01, G=26.3GPa, h=5e-8s,
T=723.15K and zero-spin `L=diag(-1.05/2,-1.05/2,1.05)/s`.
Each perturbed call starts from exactly the same copied S0/ep0. It compares
separate symmetric-gradient and temperature directions, not accumulated
material steps. The trial-stress and rate chain rules are independently
formed before calling the new radial helper. An inactive control remains
strictly inside its branch.

The new optional `trial_direction` helper is checked with nonzero prestress
and pure planar rotation against a literal trigonometric derivative. A
noncommuting spin/strain direction uses an independent Hermitian eigensystem
matrix exponential to evaluate the trial-stress reference, with centered
steps1e-5 and5e-6. Bounds: stress derivative absolute2e-7Pa plus relative2e-6;
rate derivative absolute2e-9/s plus relative2e-6. A zero equivalent rate with
nonzero deviatoric direction is nondifferentiable and rejects; a pure-rotation
direction at zero rate can return exactly zero rate derivative.

Inputs must remain unchanged and outputs detached. No response helper may
claim physical acceptance, restart or a full coupled Jacobian. An initial
component invocation is externally capped at60s. Preserve initial source
hashes and every failure; at most two evidence-driven repair rounds are
allowed. This contract authorizes no coupled run. All checks are NOT_RUN
until separately logged.

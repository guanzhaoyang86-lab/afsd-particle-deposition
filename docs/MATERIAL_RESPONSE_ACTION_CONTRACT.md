# Saved real342 action with the new material search block

2026-10-04, before any actual new342 matrix construction or action. This is a
pure saved-journal experiment, not ProcessThermalTrial.evaluate, a root or a
physical step. It follows independent material and Schur component tests.

Inputs are the already hashed real342 base/plus/minus, full650ns parent,
configuration and old structural/continuity action arrays from PROCESS_BRIDGE.
Use exactly the same50ns interval, physical parameters and epsilon1e-7m/s.
No new NNPS, contact/material update, temperature solve, direction or epsilon.

Construct exactly one new MaterialStageResponse at the saved base journal:
full analytic same-branch corotational/J2/FKS response, saved midpoint gradient
and physical volume, exact candidate continuity derivative, material-Tmid
coupling, actual axial free-axis mask. This is a named nonsymmetric numerical
search model. Do not symmetrize, truncate singular values or change mass/modulus.
It omits evolving geometry/volume, contact derivative and BE thermal feedback.
Compared with old P, midpoint geometry and continuity linearization also change;
any improvement must not be attributed uniquely to the FKS rate term.

Apply this fixed P to the already saved(C+,map+) and(C-,map-) once each. Save
dC,dR,dcorr,E=dC+dcorr=dC-P*dR for all five fields and norms using original
root scales. Compare with both old saved P results, including all worsened or
cross-response fields. Inputs with zero dC have no reduction ratio. Report
plastic/no-contact/multi-contact/zero-event groups and full persistent-ID arrays.
Group labels are observational only; none affects assembly, force or material.

Also compute one algebraic base proposal for future cost/scope assessment,
but do not evaluate, publish, restart or declare it accepted. Mark all outputs
NUMERICAL_ONLY, accepted/restart_allowed/can_publish=false. No physical endpoint
is advanced. A small action error alone does not prove actual residual decrease.

Checks and bounds are fixed before the run: parent/source/artifact/ID identity
exact; saved D reproduction1e-10/s+1e-12 relative; saved rate1e-12/s+1e-12 relative;
reconstructed48-return stress1e-6Pa+1e-12 relative (bisection/subtraction roundoff);
linear solve backward error unchanged1e-12. No threshold is imposed to make
the action ratio pass: record it, including deterioration. Matrix force derivative
unitskg/s, candidate pressure derivativeN/(kg/m3), material temperature derivative
N/K, final matrixkg and rhsNs. Algebraic reaction is not physical external work.

Budget: one invocation, internal55s/external60s, forecast<15s from the previous
two-matrix saved action2.639s. Stop on failed input/branch/solve/budget. No second
physical experiment or root retry is implied. At most one obvious pre-action
tooling repair may be separately registered with original FAIL kept; a numerical
or model failure ends this action experiment. Component repairs are separately
bounded by their contracts. No parameter/tolerance tuning or continuation.

Before execution: material11 tests PASS; Schur tests pending/NOT_RUN; saved new
matrix/action NOT_RUN. The old500ns tail is elastic/contact-dominated and will
not be rerun on the presumption that this different plastic-state patch fixes it.

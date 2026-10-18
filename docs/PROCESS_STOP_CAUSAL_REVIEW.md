# Located-stop scout: read-only causal review

2026-10-03. This review does not run a trajectory, re-evaluate a material
update, alter a tolerance, or promote a private state. It reads the saved
JSON/gzip/NPZ arrays and current source. The reviewer authored the two small
contact helpers; their independent tests were written/executed by another
agent. This is an independently derived diagnostic/Jacobian review, not a
claim of an entirely independent second implementation of the solver.

## Evidence and publication boundary

Primary evidence is `outputs/process_stop_20261003/scout_initial/`, especially
`report.json`, `event_predictor.json`, `private/trial_00*.json.gz`,
`private/trial_002_REJECTED.json`, `private/trial_005_REJECTED.json`, and
`trajectory/accepted_history.json`. Command and external budget evidence:
`logs/process_stop_20261003/scout_initial_supervision.json` and
`scout_initial.log`. The command exited2; external4.576397s, internal4.111982s,
against the registered60s/55s limits. Forecast was45.852379s.

Readback recomputed SHA256 for all90 recorded sources and19 artifacts: all
matched. It used NumPy for projections/norms and the PySPH loader only to
inspect stored arrays; no tested force, contact, material or advance function
generated the diagnostic reference quantities. Zero additional steps ran.

There are exactly two published records:

1. Same-fields numerical-schema transition at4.9999999999999996us, bound to
   the previously published5us checkpoint, state hash
   `0876486dea1c4f26e615fb003ff8939209b6d13fe72655d952a3eb0825efa176`.
2. A newly solved5→6us interval, state hash
   `cbe85e8f8f77718bbbf393502efb709fd119f62c66946438d8106501fe2c8bb5`.

The new6us interval was actually evaluated14times and strictly audited before
publication. It was not loaded from the old REJECTED6us candidate. No event-left
state, right-limit restart, or post-event interval was published. Private
trials001/003/004 have accepted=false and restart_allowed=false in their NPZ
states; only trial001 also has a separate strictly certified published copy.
The valid continuation state remains6us. This is a27particle local elastic
benchmark, not new AFSD transport; accepted AFSD steps remain0 for this run.

## Actual invocation order

The source path is `probe_hht_stop.trial` →
`hht_stop_stage.advance_private` → bounded residual evaluations:

```text
HHT endpoint and dense midpoint from candidate acceleration
  → midpoint current-geometry corrected SPH operators
  → unchanged elastic FKS/material update and density map
  → actual endpoint internal force
  → mode-aware endpoint penalty/Coulomb boundary
  → save pending UNEQUILIBRATED candidate
  → hht_stop_contact.integrate:
       endpoint/hash/gap checks → full polynomial event/collinearity check
       → actual contact quadrature (only if these checks pass)
  → BE thermal solve with signed raw integrated Q
  → acceleration/density/temperature residual
  → elastic/density Schur correction and bounded line search
```

Only a converged private result can reach `certify_interval`, which repeats
the strict physical contact audit and compares solved impulses/work/Q before
publication. `locate_duration` is an OUTER duration search: each new h calls
the complete INNER nonlinear stage again. A predictor polynomial zero is
not the accepted event state. No actuator is present in this fixture.

## What happened, in chronological order

|Trial|Interval duration|Evidence/result|
|---|---:|---|
|001|1.0000000000000006us|14 evaluations; velocity residual3.43364e-10m/s; strictly certified and published at6us|
|002|0.4999999999999994us|Evaluation1 rejects a full-vector zero on the predictor, contact0:4; no quadrature or new state accepted|
|003|0.30167152850758784us|9 evaluations; private extended-branch root converges, velocity residual6.15749e-10m/s; endpoint already beyond stop|
|004|0.2755694698357375us|9 evaluations; private root converges, velocity residual2.52202e-10m/s; endpoint still beyond stop|
|005|0.27254762326048973us|18 evaluations; seven final backtracks rejected; inner root fails before outer duration root can finish|

Trial002's actual predictor gives a first zero at fraction0.5027858808459803,
not an oracle-imposed event time. The separately discovered partner predicts
durations2.513929404229899e-7 and2.5139294042988897e-7s for0:4 and1:22.
The two are grouped only for an attempted search. Their mean root is not a
certificate that both complete tangential vectors vanish.

Trial003 endpoint z-slips are +0.0021456269037147813 and
+0.002145626903076514m/s. Trial004 gives +0.0002511542648083731 and
+0.00025115426416713604m/s. They are PRIVATE incoming-branch continuations,
not physical reverse-sliding solutions. Their approximately6.4e-13m/s
inter-contact difference already exceeds twice the1.13687e-13m/s zero bound.
This supports caution about a joint root; finite inner residuals mean it does
not prove that the fully converged joint equations have no common root.

Trial005's line-search evidence is specific:

|Evaluation(s)|Outcome|
|---|---|
|2,3,4|Full updates reduce scaled merit653627.16 →78648.93 →9361.16|
|5|First natural-collinearity rejection,1:22; transverse coefficient1.22747e-13m/s|
|6|Half update accepted as private search base, merit4122.13|
|7–10|Full,half,quarter,eighth updates reject for0:4 or1:22 collinearity|
|11|1/16 update becomes new search base, merit3833.75|
|12–18|All seven updates,1 through1/64, reject; line search fails|

Final rejected contact1:22 has polynomial coefficient c1_y=
-1.1943047577456654e-13m/s. Independently projecting each coefficient onto
the stored incoming direction gives

`max_k ||c_k-(c_k dot d_in)*d_in|| = 1.194664339435182e-13m/s`,

versus1.1368683772161603e-13m/s: ratio1.050837866. These are distinct numbers;
the raw y coefficient is not precisely the transverse projection. The final
checked minimum full slip is1.1641316643478888e-13m/s. The scope guard rejects
before Gauss integration for this candidate.

The last evaluable base is ALSO not an equilibrated solution:

|Raw residual|Measured|Unchanged bound|
|---|---:|---:|
|Velocity|3.833754755341113e-6m/s|1e-9m/s|
|Position|3.1607642588637583e-13m|1e-12m|
|Density relative|2.357949226818703e-15|1e-10|
|Temperature|2.0463630789890885e-12K|5e-12K|

The velocity maximum is at PID19,y; PIDs25,7,1 have almost the same corner
residual magnitude. It is NOT at event PID22. The event particles' own
`h*max(abs(acceleration residual))` are about4.28e-11m/s. Thus the global
mechanical root remains3833.75times outside its velocity gate even though
other residual gates pass. Relaxing a5percent collinearity excess cannot by
itself turn this state into an accepted coupled solution.

## Why this is not the original negative-friction-heat error

No recorded error in this run is `negative midpoint friction dissipation`,
negative integrated raw heat, or a temperature/yield clamp. Trial002 fails
the zero-event check; trial005 fails the registered private-continuation
scope and then the nonlinear line search.

Actual signed private quadrature is retained: trial003 has three negative
nodes per selected contact, and trial004 has one. Each selected contact's
net raw Q remains positive, approximately8.54689e-9J and8.52620e-9J,
respectively; total all-contact Q is1.41803833e-6J and1.29314209e-6J.
Negative local work here is a documented unphysical branch extension used
only to search for a shorter event interval. It was not clipped or published.

The zero detector's nested `zero_slip_audit` uses dummy ZERO tractions solely
to inspect a velocity polynomial. Its displayed heat=-0 is not actual
friction heat evidence. No actual quadrature result exists for the final
collinearity-rejected evaluation; it must not be reported as Q=0 or PASS.

## Root-cause classification

- **CONFIRMED:** the private naturally-collinear extension guard excludes
  successive Newton/search directions, eventually every allowed backtrack.
  This prevents completion of the INNER root, so the OUTER h search stops.
- **CONFIRMED:** the source preconditioner contains elastic and density/
  pressure response, but no derivative of the actual normalized tangential
  Coulomb force or the dense contact impulse correction. The actual residual
  does recompute these quantities; this is a linearization limitation, not
  omitted physical force in the residual.
- **SUPPORTED:** small transverse perturbations become increasingly important
  near zero slip. At selected Gauss nodes, max(C/|u|) rises from about120Ns/m
  in trial001 to5758 and21159Ns/m in trials003/004. Transverse coefficients
  rise from order1e-17 to order1e-13 during the failed search. These facts and
  the derivative below support testing missing tangential response.
- **POSSIBLE, NOT PROVEN:** a normalized-force response absent from the
  preconditioner amplifies transverse floating-point/asymmetry perturbations
  enough to cause the observed scope rejection. No Jacobian spectrum, complete
  directional derivative test or fixed-state corrected solve was executed.
  The corner-dominated global residual prevents attributing everything to
  the final selected contact alone.
- **POSSIBLE:** nested h solves with finite inner tolerance, changing line
  searches and a narrow collinearity domain are poorly suited to locating a
  roundoff-scale full-vector event. The run demonstrates inner failure, not
  that an augmented method necessarily succeeds.
- **BLOCKED:** no coupled event-left state, right-limit restart, common-end
  post-event accuracy or process feasibility has been established. Neither
  scalar tests nor an accurately located independent-reference event remove
  this solver block. The accepted6us trajectory still has negative slip
  where the independent reference has already stopped near5.974656588us;
  existing temporal error must remain visible.

## Smallest next numerical investigation — NOT_RUN

Do not change friction physics, tolerances, quadrature acceptance, mass/moduli,
the old material, or collinearity/zero guards. Do not rerun the actual scout
first. Register and test the COMPLETE contact contribution to the coupled
residual derivative; an endpoint-only Coulomb Jacobian is insufficient.

For fixed normal n, P=I-n*n^T, u=P(v-w), r=|u|>0, t=u/r, C=mu*Fn, let
sigma=+1 on ordinary sliding and sigma=-1 on the PRIVATE extended side:

```text
Ft = -sigma*C*t
Kv = dFt/dv = -sigma*(C/r)*(P-t*t^T)
Kx = d(Fn*n+Ft)/dx = -kn*n*n^T + sigma*mu*kn*t*n^T
```

The PRIVATE negative side has the OPPOSITE transverse Jacobian sign. It must
not be treated as ordinary dissipative physical sliding. These smooth
formulas hold away from zero and branch switches, not at the override itself.
They assume this fixture's fixed normal and position-independent prescribed
wall velocity; they are not a finite square-tool/actuator Jacobian.

With a=force weight, gamma,beta as registered and d=beta-gamma/2:

```text
V_A(s) = d v(s)/d A1 = h*[gamma*s + 6*d*s*(1-s)]*I
X_A(s) = d x(s)/d A1 = h^2*[(gamma/2+3*d)*s^2-2*d*s^3]*I
D_A Jactual = h*integral_0^1 [Kx(s)*X_A(s)+Kv(s)*V_A(s)] ds
```

The source's force-unit root is
`R=M*A1-a*F1-(1-a)*F0-DeltaJ/(gamma*h)`.
Holding the accepted left state fixed, its contact derivative is

```text
D_A R_contact = (1/(2*gamma)-a)*D_A F1 - D_A Jactual/(gamma*h)
```

Equivalently, for the impulse-unit root gamma*h*R, it is
`h*(1/2-a*gamma)*D_A F1-D_A Jactual`. This explicitly differentiates
DeltaJ=Jactual-h*(Fcontact0+Fcontact1)/2; dropping that endpoint correction
or using the wrong PRIVATE sign tests a different residual.

Minimum component checks before any saved-state experiment:

1. Independent analytic directional derivatives of physical and PRIVATE
   nonzero full-vector Coulomb, including transverse perturbations, normal
   gap perturbations and coordinate rotations. Use finite differences only
   where the complete perturbed path remains in the same declared domain.
2. Independent finite differences of the COMPLETE contact residual above on
   a manufactured turning HHT path with positive gap. Compare full DeltaJ,
   endpoint force and h dependence; retain quadrature/refinement uncertainty.
   A private cross-zero integral needs its own branch-aware limiting analysis;
   do not infer differentiability from smooth formulas at individual nodes.
3. Retain the existing literal1D reversal/stick checks. Add/compose a genuine
   two-component stop and a component-only sign turn so that a scalar root
   cannot create a false event. An available independent polar reference is
   m=1,C=2,Fother=(1,0),u0=(1,1):
   `t_stop=(2*sqrt(2)+1)/3`, `Q=5/4+sqrt(2)/3`, with feasible sticking.
   This is an independently integrated constant-force Coulomb ODE reference,
   not a prescribed free-material trajectory or a new AFSD parameter set.
   The already documented m=C=1,Fother=(-2,0),u0=(1,1) case crosses ux=0
   with uy=sqrt(sqrt(2)-1)>0 and must NOT become a full-zero event.
4. Only after these checks, preregister at most one private diagnostic on
   the SAME saved trial005 parent/duration. Keep24 inner evaluations and
   existing domain gates, no publication or physical advance. Determine
   whether the corrected linearization lowers the global residual without
   leaving the allowed branch. Failure remains evidence, not a reason to
   increase tolerances or claim the fix works.

If a separate augmented event root is justified, add duration h to the
unknowns while retaining all physical DOFs and the actual material/thermal
residuals. A scalar event row is permitted only with an independently checked
complete transverse-zero condition; otherwise use chronological separate
events or a genuinely coupled vector event formulation. For a collinear row
g=d_in dot P(v1-w1), its direct derivatives are
`g_A=gamma*h*d_in^T*P` and
`g_h=d_in dot P[(1-gamma)*A0+gamma*A1-a_wall]`.
The full R_h also differentiates material, density, thermal, contact path and
duration; it cannot be replaced by the kinematic derivative alone. Count
all derivative/line-search evaluations. No velocity projection, discarded
private failure or unrecorded right-history reset is allowed.

All four next checks and the augmented solve are **NOT_RUN** in this review.
The available evidence motivates a small linearization/event benchmark,
not a new backend, regularized friction, global constitutive rewrite, or a
production trajectory. There is no basis here to promise that adding one
Jacobian term will finish AFSD reproduction.

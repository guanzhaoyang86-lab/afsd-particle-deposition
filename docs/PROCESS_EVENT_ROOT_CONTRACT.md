# Event-response joint root — preregistered 2026-10-03

Goal remains affordable true process advancement and material transport, not
one more microstep or all-green unit tests. This experiment targets the measured
missing contact response in PROCESS_GLOBAL_NONLINEAR_DIAGNOSIS. It changes only
the numerical candidate preconditioner, not the physical map, heat, constitutive
law, contact kick, geometry, original actuator axes, or raw acceptance tolerances.

## Sources and states must not be conflated

The physical parent is published stage_004 of nonlinear_continuation_20261003,
t=850.393632910904 ns. The last REJECTED process_global_20261003 candidate is
ONLY a frozen numerical linearization anchor. Its source/artifact/parent/config
hashes, dt, material history hash and persistent ID order must match.
No rejected x,v,rho,T,S,ep becomes the initial physical state. The new root starts
with the same original predictor from the accepted parent, not the rejected root.

Use actual old-kick event response B=dvbar/dF and C=dv1/dF, holding contact geometry,
normals/normal forces/wall velocities/v0/dt/constraints fixed. Units s/kg.
For the independent reversal m=1,u0=1,F=-2,muFn=1,dt=1: B=19/54,C=7/9,
vbar=-1/18,v1=-2/3,Q=7/18. B is neither free dt/(2m) nor C/2.
Response thresholds and branch diagnostics are preregistered separately in
PROCESS_EVENT_RESPONSE_CONTRACT. Any unavailable cloud response stops this
experiment; no silent local-mass fallback or contact-DOF deletion.

K is frozen accepted-geometry elastic stiffness; E=trace(corrected D),
H=.5 E^T diag(V bulk/rho_ref). Eliminate approximate continuity:

    (I + dt/2 B K) delta_vbar = dbar + B H drho
    delta_F = H drho - dt/2 K delta_vbar
    delta_x = dx + dt B delta_F
    delta_v1 = dv1 + C delta_F
    delta_rho = drho - dt diag(rho_ref) E delta_vbar

All d are true map-minus-candidate defects. Singular B is legal; never invert B.
The resulting matrix is generally nonsymmetric; use a direct full free-block
solve and explicit prescribed-axis corrections, no SPD assumption/pinv.
Linear backward error<=1e-12; full raw five-field gates remain x1e-12 m,
v/vbar1e-9 m/s,rho relative1e-10,T5e-12 K, followed by unchanged full acceptance.
Frozen event/elastic tangents omit geometric/prestress/FKS tangent changes and
are numerical approximations, not physical replacement models.

## Execution gates, fixed before actual runs

1. Old evidence freeze check and old/new independent tests; each suite<=60 s.
   A coupled spring/actual-kick manufactured root additionally uses m=1 kg,
   k=4 N/m,x0=19/36 m,u0=1 m/s,dt=1 s,muFn=1 N. Literal root
   x1=17/36 m,vbar=-1/18 m/s,v1=-2/3 m/s,Fmid=-2 N,Q=7/18 J.
   Delta(K+Es)=-7/18 J; compare state to original raw tolerances and energy
   at absolute1e-9 J. This invokes the real contact kick in the map but uses
   literal independently derived references, not a manufactured AFSD claim.
2. One full-cloud response probe with no material/physical advance. Prior actual
   complete-map calculation costs about6 s; central h,h/2 uses up to13 local
   kicks per contacted particle, suggesting roughly78 s plus IO. Preregister
   conservative180 s internal /210 s external bound, always below600 s. This is
   a cost estimate, not measured derivative-cloud performance. Save each completed
   particle response and all failures; if budget/branch evidence fails, stop.
3. Only with a complete response artifact, measure one unchanged physical trial
   plus preconditioner build/apply and equivalent journal IO; inner80/outer90 s.
4. One new fixed500 ns root, <=32 evaluations, using cached detached numerical
   B/C. Launch only if 1.5*(build+32*trial_and_proposal_and_IO)+30 <=450 s.
   Inner540/outer570 s, no dt/iteration changes, no physical-history warm start.
   This is a new registered algorithm experiment, not continuation of the prior
   rejected run or extension of its consumed32-evaluation budget.
5. Root acceptance remains the original full ProcessThermalTrial acceptance.
   State class is the same BE physical map; settings/source hashes explicitly
   name the changed numerical preconditioner and derivative anchor.

Accepted/rejected journals remain separate. The constant response matrix is
saved once with a hash, never logged as a force or external impulse. Results
must be independently read back; no locally passing response certifies fullroot.

## Not claimed by success here

500 ns still under-resolves the highest elastic modes. A converged root alone
does not certify macrostep accuracy, stable process flow, finite-shoulder/bore
passage, inlet commit, free energy, recomputed restart equivalence or deposition.
After this bounded experiment, use its evidence to choose the next process
integration work package; do not simply accumulate dozens of tiny steps.

## Post-test scope correction, BEFORE any actual cloud/root attempt

The initial11 response tests include one genuine retained FAIL: the existing
smooth two-face kick at dt=.2 s uses four .05 s midpoint intervals, and its
discrete mean-response differs from the continuous ODE derivative by4.0676e-5
s/kg, beyond the preregistered2e-5. Do not change dt/tolerance, mark it PASS,
delete that test, or hide the initial log. The original suite continues to run.

A separate independently differentiated fixed-midpoint reference agrees with
the new response to B4.57e-10/C2.10e-9 s/kg in the read-only diagnostic. New
discrete-reference regressions must pass before probing. This distinguishes
the intended derivative of the *actual discrete map* from continuous-time
accuracy; it does not erase the failed accuracy claim or validate a macrostep.

Accordingly, the numerical-only gate explicitly records the one known continuous
reference FAIL and requires every other old/new check plus the new discrete
reference tests to pass. Any other failure blocks. This change is after the
initial test and before cloud/fullroot execution, not retroactively preregistered.
Receipt status is NUMERICAL_EVENT_PREFLIGHT_PASS_WITH_CONTINUOUS_REFERENCE_FAIL,
never an all-tests-PASS or process-accuracy certificate. The physical map and
full physical-state acceptance requirements remain unchanged.

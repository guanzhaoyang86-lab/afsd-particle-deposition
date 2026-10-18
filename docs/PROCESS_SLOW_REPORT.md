# Slow-load startup: measured limits, not process completion

2026-10-03. ENGINEERING_APPROXIMATE development remains active. Upstream PySPH
commit1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5. New AFSD accepted steps:0;
the actual stage004 remains850.393632910904ns, shoulder crossings0. No claim
of completed deposition, transport, paper reproduction or75% completion.

## Implemented and actually exercised

`hht_slow_inputs.py`, `reference_hht_slow.py`, `probe_hht_slow.py` and
`check_process_slow.py` implement the registered slow-startup experiment.
All81 DOFs of the27particle3D elastic block remain free; physical mass/moduli,
nonzero initial stress/velocity and the existing paired corrected operators,
FKS elastic update and implicit heat are retained. The5Hz normal programme
and constant+.2m/s tangential wall motion are explicitly NEW manufactured
loading, not a rescaling of AFSD or the earlier fast fixture. No geometry,
material, GPU/parallel or old AFSD source was changed.

The desired1ms horizon was not launched: forecast exceeded its60s reference
budget. The planned20us scout was bounded by its first unsupported condition.
All actual invocations were below60s; the longest reference was13.35s external.

|Evidence|Result|Actual scope|
|---|---|---|
|Previous4067-entry freeze before work|PASS|File integrity, not physics|
|New preflight|168 PASS,10.2445s|19 suites;145 prior +23 new tests|
|Independent reference grids|PASS_PREFIX_REFERENCE_EVENT_BLOCKED|96 same-time snapshots through5.9375us|
|Cost pilot|PASS,6 evaluations|One250ns local step;external.5718s|
|S0 planned28-step startup|FAIL/STOPPED|8 published steps through2us;external1.4902s|
|One private500ns root diagnosis|FINAL_ORIGINAL_GUARD_FAIL|12 evaluations;no published state;external.5783s|
|S1/S2,20us convergence,1ms|NOT_RUN|No completed common horizon|
|AFSD continuation/transport/production|NOT_RUN|No new AFSD physical state|

New tests cover literal physical inputs/schedules, independent dense event
identity, counter/bracket integrity, actual3D boundary/labels, and mocked
publication rollback. Mock publication files are explicitly MOCK and are not
physical trajectories. The old fast-HHT velocity accuracy FAIL and every
earlier raw/contact/energy FAIL remain unchanged, not reclassified.

Independent saved-output readback PASS(integrity only): all8S0 checkpoints,
40particle properties, canonical state and VTK fields match; pilot's first
state is bitwise the same as S0's first state, not an extra segment. Nine
common-time frames through2us have maximum weightedRMS position error
1.46649e-9m, velocity.00297608739m/s and full stress26104.565Pa against the
independent reference. These are coarse-prefix raw errors, NOT an accuracy
PASS; no completed three-grid sensitivity is available. Readback external
3.04976s and no integration. Full values are in PROCESS_SLOW_OUTPUT_REVIEW.

## Actual call order and first tested-solver failure

```
accepted checkpoint
  -> initial guess(A1=A0,rho1=rho0,T1=T0)
  -> dense midpoint/endpoint x,v
  -> current midpoint NNPS + paired gradient
  -> unchanged material update + continuity map
  -> actual endpoint internal force + rigid contact force
  -> whole-path friction work/passivity + sampled direction gate
  -> BE thermal map
  -> joint acceleration/density/temperature residual
  -> bounded Schur correction/line search
  -> guards + checkpoint readback + atomic history publication
```

There is no actuator/velocity overwrite in this all-free manufactured case.
`hht_coupled.advance.evaluate` calls `_contact_work` before forming the first
joint residual. S0's first500ns interval[2,2.5]us therefore stops at evaluation1:
face0/PID3, maximum sampled Coulomb interpolation defect.0187148920923>.01.
Its minimum checked vector slip.0460024m/s is far above roundoff4.98e-15m/s;
friction Q=+5.5882274e-8J; whole-path power remains negative; no slip-zero split.
This is neither the old negative-dissipation event nor the later true stop.
Independent5point polynomial work/defect readback reproduces this rejection.

## Distinguishing a predictor veto from a genuinely inadmissible root

Registered `PROCESS_SLOW_PREDICTOR_DIAGNOSTIC.md` before one private500ns trial.
Only the exploratory sampled-direction veto was deferred; raw force/work,
all passivity/event/domain/thermal checks and24-evaluation budget remained.
The original1% gate was reapplied to the converged candidate. No source file,
parent physical array or published history changed; the result is explicitly
PRIVATE_NOT_RESTARTABLE, even though its algebraic residual converged.

After12 evaluations, velocity residual6.1755e-10m/s passes its unchanged1e-9
bound. But converged direction defect.0125263275045 still fails.01, with
positive Q5.78006e-8J and minimum vector slip.0492045m/s. This disproves the
proposal that simply postponing the check fixes the registered interval.
The independent reference's endpoint-linear traction on its actual stored
path has maximum direction defect.00938708; that different trajectory cannot
be used to declare this HHT candidate valid. No additional dt attempt was made.

## Independent first physical event (a separate later obstruction)

Both fixedRK4 grids independently approach full-vector zero at face0/PID4
near5.97466us; opposite-facePID22 is nearly simultaneous. Coarse bracket
[5.974655151367281,5.974670410156344]us; fine
[5.974655151367073,5.974662780761604]us. The transverse dense-root component
is roundoff-small, not merely signed-z reversal. All9 prefix uncertainty
bounds pass; no post-event state is accepted.

At the last accepted near-event state, internal tangential force4.67484334N
exceeds Coulomb capacity2.73923931N. Static traction is infeasible by1.935604N;
reversed friction would still leave+22940.49m/s2 tangential acceleration.
This SUPPORTS genuine reverse sliding, not forced sticking. Opening and yield
do not occur first (minimum overlap1.6188um; elastic margin313.36MPa).
An exact event state/restart and post-event solution remain NOT_RUN.

The oracle repair corrected attempt-vs-halving metadata and prevented false
inverted event brackets. Actual13 attempts/9 halvings; every physical NPZ hash
is identical before/after repair. It did not change physics, dt or the budget.

## Evidence-ranked conclusions

- CONFIRMED: S0 predictor veto occurs before residual formation; the bounded
  converged private root also fails the original contact-direction gate.
- CONFIRMED: current independent prefix is true3D evolving mechanics/heat,
  not a prescribed velocity field; saved near-stop static disk is infeasible.
- SUPPORTED: interval traction representation/time path needs a consistent
  improvement; the later physical event requires genuine reverse sliding.
- RULED_OUT for these observations: negative heat, floating-point-only defect,
  z-only stop mislabel, whole-block kinematic forcing or a mere iteration-cap hit.
- BLOCKED: complete20us coupled scout, exact stop/restart, shared-face actuator
  coupling, positive plasticity, real finite-tool/entry/transport integration.
- NOT_RUN: S1/S2, spatial sensitivity,1ms reference, actual slow restart and
  all new AFSD production. Different failed times are not a convergence study.

## Evidence and next action

Commands: `logs/process_slow_20261003/COMMANDS.md` and REFERENCE_COMMANDS.md.
Read `PROCESS_SLOW_REFERENCE_DESIGN`, `EVENT_AUDIT`, `OUTPUT_REVIEW` and
`NEXT_WORK_PACKAGE`. Outputs are under `outputs/process_slow_20261003`;
private/rejected records under the matching logs directory never authorize restart.

Next work is a consistent curved-slip contact interval method, with independent
impulse/work/heat identities and unchanged final gates, followed by located
zero-slip reversal. Merely deferring the guard, raising its tolerance, adding
pre-event microsteps, or relabeling the private root cannot complete this task.
`PROCESS_SLOW_FREEZE` certifies preservation only, not scientific acceptance.

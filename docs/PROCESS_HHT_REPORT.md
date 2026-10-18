# HHT process-integration work package — local coupled evidence

2026-10-03. ENGINEERING_APPROXIMATE development, NOT_AFSD_PROCESS_VALIDATED.
PySPH pin1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5 remains the only backend.
The existing real AFSD stage004 at850.393632910904ns remains unchanged; this
work uses a separately named27-particle/all81-free-DOF manufactured benchmark.
No crossing/deposition or positive plastic-work claim follows from it.

## What is implemented

- `hht_linear`: inertia-retaining HHT with a=.9 and smooth a=1 control;
  analytic/discrete independent oscillator and augmented-energy checks.
- `hht_path`: compatible cubic position/quadratic velocity, exact linear-force
  impulse/work, explicit signed force-filter work, whole-path friction power
  audit and zero-slip events, including independently found double-touch zero.
- `hht_coupled`: detached same-history current3D NNPS/corrected paired force,
  unchanged Jaumann/FKS elastic branch, midpoint continuity, raw contact heat
  and conservative BE temperature in one residual solve. Plastic increments
  reject in this local benchmark scope. No frozen-force replacement.
- `probe_hht_coupled`: actual force-driven two-face loading/unloading with
  explicit rigid-wall speed reversal at1us; fixed grids to2us, reference
  comparison, source/config receipts, accepted-history publication, canonical
  checkpoints and solved3D VTK. Tool/AFSD inlet integration is NOT_RUN.

The independent reference implements current all-pair SPH/Jaumann/heat ODE
without the tested backend calls. It establishes temporal accuracy for the
selected semidiscretization, NOT continuum stress accuracy or the paper.
The original E2 continuum error and Jaumann closed-cycle risk remain.

## Initial evidence retained, not erased by later repairs

- Reference initial import: FAIL because SciPy is absent; no ODE evaluated.
  Before actual ODE execution selected independent RK4(128/256steps per half).
  Reference round1PASS: maxx difference1.99e-17m,v3.62e-11m/s,S2.99e-4Pa,
  T4.72e-12K. Whole reference1.333s, external1.660s. No library installed.
- Initial path suite8PASS did not find double-touch zero slip. Independent
  audit did; original source and failed probe remain. Audit repair1 adds
  component roots/vertices and finite-input validation; path10PASS.
- Initial combined preflight: old103PASS, then new linear test import FAIL
  because its standalone command lacked src path. A test-entry import fix
  gives preflight_round1:132PASS in8.110s. No physics changed.
- Pilot initial250ns step: PASS,11evaluations, stage.1854301s including output;
  external.6030927s; processpeak60,502,016bytes. Forecast56steps30.768s under
  preregistered50s limit. Pilot is not common-end-time evidence.
- Sweep initial: FAIL at step6 after5accepted local stages to1.25us;
  external1.6298052s. Velocity rootresidual1.0500136e-6m/s versus1e-9;
  densityresidual/rho0=1.66661e-9. No common-final-time accuracy comparison.
  All candidate/accepted evidence is preserved under sweep_initial.

## Proven numerical repair mechanism

At the failed candidate, the original acceleration matrix contained bulk
stiffness as though density followed the acceleration correction. The separate
density update omitted that coupling and the acceleration RHS omitted the
pressure-force response to the density residual. This is a numerical block
linearization error, not evidence to change friction, heat or material physics.

Independent block diagnosis: omitted pressure term norm.003485N versus force
residual.001816N; old block equation residual.003554, corrected2.52e-18.
On the actual candidate, a1% old direction raises scaled residual1050.0136
to1050.4333; corrected1% lowers it to1039.1241 and full correction to60.2470.
Evidence: `logs/process_hht_coupled_20261003/repair1/diagnosis.json`.
The source-before snapshot and initial failure remain; no max-evaluation,
physical parameter, error tolerance or grid is increased or relaxed.

## Final-source execution and independent readback

`preflight_final/report.json`: **145PASS**,9.051155s internal9.646681s shell.
This is103existing component tests plus42new HHT tests. It is NOT a claim
that every project test passes: previous event147/145PASS+2FAIL and all older
scientific failures remain unchanged. The restart script's last pre-execution
status-label correction invalidated the earlier preflight_repair1 receipt;
it is preserved, and the complete final-source suite was rerun before the pilot.

Repair1 pilot:6evaluations,.1387959s complete stage,.6045132s external.
Preregistered forecast25.54514s. Actual repaired sweep: **all8/16/32steps reach
the same2us**, internal5.1502901s/external6.3140948s, peak process65.0MiB.
There are56published local benchmark stages,0new accepted AFSD stages.

| dt | Steps | Whole trajectory wall time | Mean whole-step cost |
|---|---:|---:|---:|
|250ns|8|.9220313s|.115254s|
|125ns|16|1.5766014s|.098538s|
|62.5ns|32|2.6171180s|.081785s|

The finest solve takes3–4evaluations per step (127total). No physical mass,
time, modulus, temperature or yield clipping/scaling is used. These tiny
manufactured-block costs cannot be substituted for actual AFSD cost estimates.

### Registered common-time errors, including both wall-force limits

Each number is the maximum over that grid's common sample times, not just
its final snapshot. Reference uncertainty is independently below1% of limits.

| Quantity |250ns|125ns|62.5ns|Fine limit|Result|
|---|---:|---:|---:|---:|---|
|Volume RMS position,m|1.79344e-9|4.62587e-10|1.17207e-10|1e-9|PASS|
|Volume RMS velocity,m/s|3.02356e-3|8.02231e-4|2.05200e-4|2e-4|**FAIL**|
|Volume RMS full stress,Pa|31456.5|8166.92|2074.77|7000|PASS|
|Max temperature,K|1.87226e-8|4.87108e-9|1.23634e-9|1e-7|PASS|
|Face force component,N|.175738|.0451417|.0114732|.7|PASS|
|Face impulse component,Ns|1.33553e-7|3.51573e-8|8.97282e-9|1.4e-7|PASS|
|Face wall work,J|1.36012e-6|3.37513e-7|8.42188e-8|2e-7|PASS|
|Face raw friction Q,J|9.44263e-9|2.44803e-9|6.21242e-10|2e-7|PASS|
|Face reservoir heat,J|6.88761e-13|3.45278e-13|1.72874e-13|2e-8|PASS|

All9decrease criteria pass. Velocity refinement ratios are.265327/.255787;
finest velocity exceeds its fixed allowance by2.60024%, at2us. Therefore the
overall report stays **FAIL_LOCAL_TEMPORAL_ACCURACY**, not a rounded PASS.
No added finer grid, changed tolerance or extra numerical repair is attempted.
Observed convergence is meaningful but is not the failed absolute-accuracy gate.

The independent reader audits all56published snapshots, semantic arrays,
source/config/reference/artifact hashes and computes work with independent
five-point quadrature. Endpoint-force discrepancy<=2.842e-14N;
cumulative face work<=4.405e-20J,Q<=1.186e-20J,impulse<=1.085e-19Ns;
saved-temperature BE equation<=8.338e-15J/node. It reproduces the velocity
FAIL exactly. All3control switches preserve physical fields and initialize
right F/m correctly. Details in PROCESS_HHT_OUTPUT_REVIEW and
`logs/process_hht_20261003/reference_output_*repair1.json`.

Finest minimum certified signed slip=.111875m/s; conservative complete-path
overlap bounds1.50730–2.20595um; maximum Coulomb sampled defect8.86425e-5.
Minimum low-rate FKS margin312.096MPa, maximum temperature change3.65644e-5K.
Plastic work is zero because this registered path is elastic, not because a
plastic source was disabled or replaced with invented data. Positive plastic
work accuracy is NOT_RUN. Signed force-filter work can be positive or negative;
it is not always-positive dissipation or heat. Fine max normal temporal
work/potential residual1.65523e-8J is separately retained.

## Current decision and remaining work

- Coupled **code implemented**, independent component tests PASS, same-final
  time refinement observed; absolute temporal accuracy gate **FAIL**.
- Strict stage016→032 recomputation infrastructure and7unit guards PASS;
  actual recomputation **NOT_RUN** because its registered sweep gate failed.
- Coupled a=1 control, spatial refinement, shared-particle multi-face/actuator
  HHT contact, actual stop/stick/restart, positive plastic work, and AFSD finite
  tool/inlet/process transport are **NOT_RUN** under the new HHT variant.
- Full nonlinear stored free energy/algorithmic jump are NOT_CERTIFIED. Old
  Jaumann and spatial free-surface errors have not been erased.
- No permission for three-case production or long/paid/GPU tasks. Actual AFSD
  time remains850.393632910904ns, shoulder crossings0. Project NOT_COMPLETE.

Do not make the next work package another tolerance chase or more same-scale
microsteps. The present forcing itself completes a cosine stroke in2us, so
this test does not establish feasible timesteps for300rpm slow process motion.
Next must establish a **process-timescale accuracy/cost route retaining real
inertia and startup**, then include positive plasticity and located boundary
events before finite-tool/inlet integration. Its registration must separate
high-frequency filtering from low-frequency material motion and keep this
near-limit FAIL. A hot/quasistatic reinitialization is a separate physical
variant requiring user direction; do not clear the old checkpoint velocity.

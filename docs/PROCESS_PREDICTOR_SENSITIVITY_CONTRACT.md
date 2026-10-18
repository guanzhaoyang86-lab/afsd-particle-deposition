# Frozen 4/8 increment comparison (2026-10-02, before the new run)

Scope: ONE new eight-increment PredictorOnset run. No changes to solver,
material, contact, events, source/configuration, physics or previous outputs.
Baseline: outputs/process_predictor_onset_20261002/n3_i4 (complete .2 s).
New: outputs/process_predictor_sensitivity_20261002/n3_i8, from initial state.
Prior archive 1658 files verified PASS; only HANDOFF/IMPLEMENTATION_PLAN/
OUTPUT_SCHEMA may change among archived files. Old i8/NOT_RUN remains historical.

## Preconditions and cost

Use the same frozen Q7_PREDICTOR_ONSET_V1, 27 points, 1 cm cube, nine fixed-normal
contacts, J2 Y=90 MPa, E=70 GPa, nu=0, rho=2700 kg/m3, mu=.1, prescribed293.15 K.
Normal40 um and tangent20 um, direction(.6,.8), four .05 s phases unchanged.
Only increments_per_phase4 ->8: requested intervals16 ->32 and .0125 ->.00625 s.
Internal events produce nonuniform accepted intervals; this is NOT a CFL test.
Same14 shared corrections, line search8, event limit12/request, bisection40,
force9e-5 N, root1e-8, slip1e-11 m, signed-Q threshold-1e-10 J unchanged.
Normal detachment is still unsupported. No state/heat clipping, mass/time/modulus
scaling or source tuning. New postprocessor tests first, then ONE run only.
Measured baseline76.738 s/40 accepted/24 events/886 returned candidate roots.
Estimate90-150 s; external timeout180 s (<600 s); singlethread; no automatic
retry. Step count unknown until event detection;32 requested intervals is not32
accepted steps. Timeout/failure retains accepted/rejected data and stops attempts.
Full old physics suites are NOT_RUN this turn, preserved by hash; solver unchanged.

## Measurement contract (not calibration, no experimental reference)

Require BOTH runs complete to .2 s and actual frozen reader PASS_COMPLETE_READBACK.
Check artifact/source hashes and same complete source inventory; resolved configs
must be identical except increments_per_phase. Static ID, tracer, reference
position, volume, mass, constraints and contact ID/feature arrays identical.
All accepted frames finite, increasing times, unambiguous persistent IDs.
Compare17 common requested times k*.0125, k=0..16, including all phase ends.
Match actual accepted frames within1e-12 s, uniquely; NO interpolation,
extrapolation, nearest partial endpoint or matching by raw step number.
Except initial, matched frames must be REQUESTED, not an adjacent event frame.
Verify wall displacements on these times agree to1e-14 m. Use contact force SUM
(wall force on material); do not use record.reaction_N, which includes supports
and is almost zero at equilibrium. Material displacement and ep matched by ID.

For each field, peak norm over matched frames/particles gives denominator
D=max(peak4,peak8,.01*nominal_scale). Difference=peak norm(field8-field4)/D.
Vector Euclidean norms; tensor Frobenius norm; scalar abs only for DIFFERENCES
and diagnostic norms, NEVER replace signed physical work with its magnitude.
Report absolute, normalized, final and per-common-time differences plus raw final
values. Nominal scales: u40e-6 m; contact force9000 N=Y*L^2; stress90e6 Pa;
ep.004=normal displacement/L; cumulative energy.36 J=Y*L^3*.004.
Fields: u, contact resultant, contact tangential resultant (nominal900 N=mu*Y*L^2),
stress, ep, volume-mean ep, cumulative plastic workP, signed friction workQ,
wall work and stress work. P/Q independently sum saved particle/contact arrays.
Additional ep max and material Ws-Es-P are raw diagnostics, not extra heat.

Predeclared sensitivity flag: <=5% WITHIN_SCREEN, >5% SENSITIVE_WARNING.
This is an engineering screening choice (one part in20), NOT a solver tolerance,
accuracy bound, experimental error, convergence certificate or physical acceptance.
The1% nominal floor avoids near-zero division; always expose actual denominator.
Do not change flags after seeing the new result. Warnings do not mandate repeated
tuning; retain them as limitations and use bounded component work next.

## Events and signed work

Extract newly eligible contact IDs from each event: not previously released AND
abs(cone_margin)<=original9e-5 N, feature from saved contact inventory.
These are explicit CONE_LIMIT arrivals at the boundary, not proof of subsequent
sliding and not an inventory of every mode change at phase corners.
Compare by (phase,particle_id,feature_id,occurrence), not global event index:
simultaneous contacts can split/group. Phase is ending phase, boundaries within
1e-12 assigned to previous phase. Report complete inventory differences, grouping,
order changes and time shifts; missing counterparts must NOT be silently dropped.
Matching inventories: |delta t|<=.00125 s WITHIN_SCREEN (10% coarse macro interval),
otherwise SENSITIVE_WARNING; not the much tighter1e-12 root-bracketing tolerance.
Mismatched inventories: EVENT_INVENTORY_MISMATCH, no forced time pairing.

Independent saved-array audit, no solver/material imports:
Q=-tau_mid dot (du_contact-dwall)_yz, preserve raw signed entries;
Wwall=sum(Fcontact_mid) dot dwall; Wmaterial=sum(Fcontact_mid dot du).
Reconstruction errors<=1e-14 J, interval contact identity<=1e-12 J (float sums on
sub-J benchmark, consistent with prior independent audit); signed Q>=-1e-10 J
is the unchanged solver gate. Report minima/counts even if below0 by roundoff.
Ws-Es-P remains signed and material_energy_certified=false regardless of flags.
No inference that refinement fixes stored energy or old closed-cycle/spatial FAILs.

## Outcome boundaries

Partial or invalid runs: BLOCKED_NOT_COMPARABLE, no endpoint sensitivity claim.
Completed runs: COMPARISON_COMPLETE plus WITHIN_SCREEN/SENSITIVE_WARNING per field.
Two resolutions cannot establish convergence order, continuum truth or paper fit.
No new spatial/thermal/rotating-tool/AFSD/production runs; these remain NOT_RUN.
If bounded comparison completes, next item is rotating-wall/local-feed interface
design and independent component benchmarks, not tuning this completed small block.

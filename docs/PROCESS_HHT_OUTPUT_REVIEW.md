# Independent saved-output review of the joint HHT benchmark

2026-10-03. Read-only audit; no trajectory, solver correction, material update
or thermal advance was executed. Scope is the selected 27-particle local
benchmark, not an AFSD checkpoint or deposition run.

Latest outcome: the completed repair-1 sweep is independently confirmed
`FAIL_LOCAL_TEMPORAL_ACCURACY`. All three grids reached 2 microseconds; the
finest velocity error is 2.60024% above its unchanged limit. The full review
below preserves the original incomplete-sweep findings and then records the
completed sweep. Restart execution remains NOT_RUN under the accuracy gate.

## Initial pilot and incomplete sweep

Reviewed `outputs/process_hht_20261003/pilot_initial` and `sweep_initial`.
The pilot publishes one 250 ns stage. The initial coarse sweep publishes five
250 ns stages, ending at 1.25 microseconds, then rejects step 6 for bounded
line-search failure. No full common-end-time accuracy claim follows from it.
The original rejected record and all five accepted snapshots are preserved.

Reader: `logs/process_hht_20261003/reference_output_readback.py`.
Output: `logs/process_hht_20261003/reference_output_readback_initial.json`.
Command used the pinned environment and an external subprocess timeout of
60 s. Exit 0; reader wall time 0.031456 s, complete external process 0.542194 s.
It imports only the independent reference's data/pair formulas and PySPH's
checkpoint loader, never the tested HHT step/path, driver, material update,
or heat-advance implementation. It checks the saved artifact hashes and
reconstructs endpoint forces, interval polynomial paths, face impulses/work,
signed Q and the BE heat equation directly from saved arrays.

All recorded artifact hashes and published checkpoint hashes match. The
driver's active source is undergoing a separately recorded preconditioner
repair; this review does not claim old receipts certify a changed source.
No missing final `frames.json` was substituted with rejected output: the
reader uses only `accepted_history.json` and the checkpoints it publishes.

| Independent reconstruction check | Maximum discrepancy |
|---|---:|
| Actual endpoint total force | 1.24345e-14 N |
| Cumulative face impulse | 1.35526e-20 N s |
| Cumulative face wall work | 1.69407e-20 J |
| Cumulative face raw friction heat | 1.69407e-21 J |
| Cumulative face reservoir heat | 8.47033e-21 J |
| Saved force-interpolation work | 2.11759e-22 J/node |
| Signed force-filter work journal | 1.32349e-22 J |
| Root-residual work journal | 4.52941e-22 J |
| Per-particle BE heat equation | 5.26014e-15 J |

These checks use independently generated five-point Gauss integration of the
saved linear force / quadratic velocity polynomials. They do not reuse the
driver's two-point quadrature values as reference. Heat is reconstructed with
the raw radial pair formula at the actual cubic path midpoint, midpoint
density, saved end temperature and per-face reservoir data. Algorithmic work,
normal temporal residual and return diagnostics are absent from these heat
sources. Raw friction heat is positive on every reviewed contact interval;
no absolute-value or clipping operation is needed to reproduce it.

The face inventory is disjoint in this fixture. The driver's per-node wall
heat is therefore counted exactly once per node when accumulated by face.
This bookkeeping alone would not certify a future shared-particle corner
contact: multiple reservoirs on one particle would require separate face
heat values to avoid assigning the summed per-node heat to every face.

## The one-microsecond control switch

Stage 4 integrates to the switch with the left wall-velocity limit, including
its cumulative wall work, impulse and Q. It then records the left force and
acceleration and initializes the right physical force/acceleration. Stage 5
uses that right history. The switch journal preserves identical physical
particle hashes before/after, including x,v,S,rho,T and labels.

The independently recomputed right force agrees to 5.32908e-15 N, and saved
right acceleration equals right F/m exactly. The largest force-component
jump is 6.08148 N. The partial acceleration-storage jump independently reads
`+5.758893315640344e-11 J`, versus saved
`+5.758893315640472e-11 J`. This is a numerical storage term with full
nonlinear algorithmic-energy jump explicitly uncertified; it is not heat or
instantaneous material work. Physical energy/velocity are not reset.

The driver computes candidate integrals before checkpoint/VTK writing and
before the control restart. `frame` copies accepted particle state and
integrals for output. Neither VTK values nor visualization/diagram state are
inputs to the solve or the work/heat calculation. The independent checkpoint
and polynomial reconstruction above confirms the numerical ledger without
reading the VTK files.

## Meaningful errors beyond the nonlinear root failure

No force-counting or ledger sign defect was found in the reviewed saved
stages. Nevertheless, the coarse step has material temporal error relative
to the independently integrated continuous semidiscrete reference:

| Saved time | Velocity RMS error | Full stress RMS error | Max face wall-work error |
|---|---:|---:|---:|
| 0.25 microseconds | 2.99851e-4 m/s | 1.45685 kPa | 1.84094e-7 J |
| 0.75 microseconds | 7.25808e-4 m/s | 12.3723 kPa | 1.14787e-6 J |
| 1.00 microseconds | 6.68649e-4 m/s | 19.8544 kPa | 1.36012e-6 J |
| 1.25 microseconds | 7.96034e-4 m/s | 27.0592 kPa | 1.15952e-6 J |

At 1.25 microseconds, displacement RMS error is 1.41140e-9 m, temperature
error 5.65320e-9 K, face-force error .131327 N, raw-Q error 1.82773e-9 J,
and face impulse error 8.12265e-8 N s. These are partial-time observations,
not a final failed/passed accuracy decision against the registered finest
grid. They must be retained when assessing actual completed refinement.

The coarse normal potential/work temporal residual reaches 9.80029e-7 J in
stage 3. The component explicitly approximates wall velocity and traction
linearly within each step, whereas the continuous reference follows the exact
cosine normal wall programme. Thus exact interpolation-work bookkeeping does
not imply exact continuous wall work or zero penalty-energy mismatch. The
residual is recorded separately and is not assigned to algorithmic loss/heat.

Force-filter work is signed: the first four increments are positive
`8.04703e-9`, `1.43404e-8`, `1.45403e-8`, `3.82190e-9 J`; stage 5 contributes
`-1.31612517e-7 J`, leaving cumulative `-9.08629277e-8 J`.
Calling this per-step quantity nonnegative dissipation would be false. The
current field names and metadata correctly call it signed algorithmic filter
work. Root-residual work is separately about `-4.12358e-14 J` cumulatively.

Status: saved pilot/partial-sweep integrity and signed-ledger reconstruction
PASS; full common-end-time accuracy/cost and resumed-trajectory equivalence
remain NOT_RUN in this review. The registered independent reference and all
original failures are unchanged.

## Completed repair-1 sweep: independent final decision

Reviewed all 56 published checkpoint/diagnostic pairs in
`outputs/process_hht_20261003/sweep_repair1`: 8, 16 and 32 steps, each to
2 microseconds. This review ran no new trajectory and did not change any
input, tolerance, grid, solver or independent-reference file.

New raw readback: `logs/process_hht_20261003/reference_output_readback_repair1.json`.
Independent error/decision summary:
`logs/process_hht_20261003/reference_output_summary_repair1.json`.
Reader runtime 0.230876 s internal / 0.744512 s external; exit 0 under the
60 s external timeout. The follow-on scalar summary took 0.322472 s external.
The completed simulation's own recorded runtime is separate from reader cost.

Every source hash in the completed report, reference source hash, resolved
configuration hash, artifact hash and journal checkpoint hash matches.
Additionally, the reader checks every actual saved particle property against
the complete semantic state snapshot: no inventory or value discrepancies.
The reference remains `reference_round1`, unchanged from its preregistered run.

The independently reconstructed errors match `comparison.json` to roundoff;
velocity error agrees exactly. The decisive maximum occurs at the final time,
not an unmatched intermediate sample:

| Metric, maximum over each grid's saved common-reference times | Finest error | Registered limit | Result |
|---|---:|---:|---|
| Displacement RMS | 1.172067386e-10 m | 1e-9 m | PASS |
| Velocity RMS | 2.052004835e-4 m/s | 2e-4 m/s | FAIL |
| Full stress RMS | 2074.766827 Pa | 7000 Pa | PASS |
| Temperature | 1.236344360e-9 K | 1e-7 K | PASS |
| Face reaction component | .01147321958 N | .7 N | PASS |
| Face impulse component | 8.972824223e-9 N s | 1.4e-7 N s | PASS |
| Face wall work | 8.421880351e-8 J | 2e-7 J | PASS |
| Face raw friction heat | 6.212417437e-10 J | 2e-7 J | PASS |
| Face reservoir heat | 1.728737511e-13 J | 2e-8 J | PASS |

Velocity maxima are `3.023555866e-3`, `8.022311318e-4`, and
`2.052004835e-4 m/s` on the three grids. Successive ratios are .265327 and
.255787. All nine metrics satisfy the registered decreasing-error ratio
<=.7, but refinement trend does not override the absolute velocity limit.
The exceedance is `5.200483516e-6 m/s`, or `2.600241758%` of the allowance;
it is vastly larger than the independent-reference velocity uncertainty
`3.61540e-11 m/s`. It cannot be attributed to oracle noise or rounding.

The finest full-stress maximum occurs at 1.5625 microseconds; final stress
error is 1465.546238 Pa. Maximum wall-work error occurs at the 1-microsecond
switch; final wall-work error is 9.272025975e-9 J. Reporting only final stress
or final wall work would conceal the larger registered history errors.

Independent ledger reconstruction over all 56 steps gives maxima:
endpoint total force 2.84218e-14 N; face impulse 1.08421e-19 N s;
wall work 4.40458e-20 J; raw Q 1.18585e-20 J; reservoir heat 1.60937e-20 J;
signed filter-work journal 1.58819e-22 J; node BE heat equation 8.33706e-15 J.
No force-counting, heat-sign, output-side-effect or event-limit discrepancy was
found. The three switch records preserve the physical particle hash and use
the correct right-limit F/m. Their independently reconstructed partial
acceleration-storage jumps are +5.75889e-11, +2.38172e-11 and +7.19260e-12 J;
none is passed to the heat equation.

The coarse, middle and fine cumulative signed filter work is respectively
`-1.51282676e-7`, `-3.71229446e-8`, `-8.96307106e-9 J`. Normal temporal
residual sums are `2.05122006e-7`, `5.16570040e-8`, `1.29404113e-8 J`;
maximum absolute individual residuals are `9.80028649e-7`, `1.30742478e-7`,
`1.65522629e-8 J`. These separate diagnostics decrease under refinement;
they still do not establish a certified nonlinear hypoelastic free energy or
nonnegative per-step HHT dissipation.

Final independent decision: all three local trajectories completed within
the measured run budget and eight finest absolute metrics passed, but the
registered joint accuracy gate FAILED. No finer grid or threshold retuning
was performed. Resumed-trajectory execution remains
`NOT_RUN_ACCURACY_GATE_FAILED`. This evidence does not authorize or imply
physical AFSD continuation, deposition, plastic-flow accuracy, or completion
of the full process work package.

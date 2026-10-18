# Contact454 repair — pre-test contract (2026-09-30)

Scope: contact_kick.py only for numerical physics changes. Preserve all material,
spatial, thermal, geometry, parameters, original failures and prior freezes.
Prior COUPLING_REVIEW_FREEZE: 245 entries verified PASS before edits.
No new dependencies; Git safe.directory is process-scoped to the pinned checkout
because the desktop sandbox user differs from the workspace owner.

First replay particle454 and its failed cone in new outputs/contact454_20260930.
Root-cause hypothesis: 3D projector cancellation of a large normal velocity,
repeated projection and near-zero curved slip produce artificial normal leakage
and cone two-cycles. Confirm with independent Decimal arithmetic and frozen inputs.
If supported, use orthonormal 3x2 tangent bases, preproject the base once, and
preserve all face coupling T_a^T P T_b / m. No isolated scalar mass cap.
Keep force tolerance2e-12, max2000 iterations, max4096 accepted intervals,
rtol2e-7/atol1e-12, zero tolerance128*eps*original speed scale unchanged.
Record roundoff-resolved stop events without setting particle velocity to zero.

Acceptance BEFORE new fixed-case replay:
- Original454: independent polar integral tstop=5.985479017657219e-9s,
  Q=8.164870224344708e-12J; time/heat relative1e-3 (default),1e-5 (tight rtol1e-10).
  Final tangent speed abs1e-10m/s; no rebound/artificial reflection.
- Existing analytic reversal, static/moving wall, corotation, multi-contact and
  actuator, kinetic work, corrected-gradient power, label/ID/restart unchanged.
- Independent cone oracle/rotation: force abs1e-11+rel1e-10, tangent normal
  leakage abs1e-14 times force scale. Include constrained/off-diagonal coupling.
- New manufactured oblique near-stop case: independently rotated/reversed,
  relative heat/time1e-3; no permissive heat clipping or fake residual heat.
- Raw work for each accepted interval remains -Jt dot u_mid.
- Original RED observations/files remain intact. Changed historical454 test may
  become PASS only on its unchanged physical reference, not expectedFailure.

Write refreshed test artifacts only under a new run directory, never overwrite
outputs/coupling_20260930 (frozen evidence). Keep each failed repair attempt log.
After local/full-chain gates: measure initial dt/step cost, then ONLY fixed
dx=.5mm/CFL=.1 to2us, hard600s. Stop at the next physical/numerical alarm.
No smaller CFL/grid in this narrow turn, no production/GPU/parallel/material rewrite.

Pre-run addition after a new manufactured failure (record retained):
drive/capacity=.02: h=6.883136429e-15s at t=.451937680s, time ulp
5.551115123e-17s. The old 64*eps*full_dt guard rejects a 124-ulp step.
Replace that heuristic only with actual representability checks for t+h/2 and
t+h. Keep ALL physical/work/force/zero-velocity thresholds above unchanged.
No time or velocity snapping is introduced; genuine non-progress still fails.

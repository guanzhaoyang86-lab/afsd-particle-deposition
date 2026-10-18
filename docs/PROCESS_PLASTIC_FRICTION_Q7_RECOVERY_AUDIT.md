# Q7R independent saved-state and event audit (2026-10-02)

Scope: read original/Q7R saved output with NumPy only; no contact/material
function generated these arithmetic reference values. Then execute exactly one
preregistered pure event replay, with accepted-state mutation forbidden.
No physical rerun or repair was performed by this audit.

## Accepted prefix and arithmetic checks: PASS, limited scope

Original: `outputs/process_plastic_friction_q7_20261001/n3_i4`.
Q7R: `outputs/process_plastic_friction_q7r_20261002/n3_i4`.

All 30 `history.npz` fields have equal dtype, shape and `tobytes()` content.
All seven accepted records and the two event records are exactly equal.
There are no accepted events after .05 s; both end at accepted step6.

From saved vectors, independently compute capacity=.1*N,
margin=capacity-sqrt(sum(tau**2)), distance=sqrt(sum(slip**2)), and
flow=sqrt(sum((distance*tau+capacity*slip)**2)). Compare with the unchanged
allowance `9e-5*distance+capacity*1e-11` on moving contacts.

- Minimum accepted margin: -6.399288849934237e-7 N, within 9e-5 N tolerance.
- Maximum moving flow minus allowance: -3.857986577685776e-10 N m.
- Independently recomputed raw `Q=-.5*(tau_old+tau_new) dot slip` differs
  from recorded interval Q by at most 5.421010862427522e-20 J; minimum raw Q=-0.
- Sum of midpoint contact force dotted with particle displacement increments,
  versus wall-force midpoint dotted with wall increment plus Q, differs by at
  most 1.6046192152785466e-17 J.

These establish accepted-prefix preservation and contact ledger consistency,
not continuum accuracy or material energy closure. The old material diagnostic
residual remains -0.004285349977920017 J and is not converted to heat.

## Actual Q7R failure: CONFIRMED, still BLOCKED

New failure is `cone event not bracketed`, not Newton root failure. Manifest
retains rejected target .0625 s / base_step6 separately from accepted .05 s.
The target has ID21 forced-stick margin -5.691232925791425 N, slip
9.470116246213047e-22 m and raw Q=-3.012674265576716e-20 J. This is a genuine
cone violation with merely rounding-scale Q; a negative-heat mechanism is not
the trigger. The old manifest has no bisection trace, so that alone could not
identify the local failure mechanism.

## One pure replay: diagnostic PASS, no physical progress

Source: `scripts/diagnose_process_plastic_friction_recovery.py`.
Preregistration: `PROCESS_PLASTIC_FRICTION_Q7_RECOVERY_EVENT_DIAGNOSTIC.md`.
Full evidence: `outputs/process_plastic_friction_q7r_20261002/event_diagnostic/diagnostic.json`.
Log: `logs/process_plastic_friction_q7r_20261002/event_diagnostic.log`.

Actual command inside single-thread external `subprocess.run(timeout=60)`:

```
.envs/pysph/python.exe -B scripts/diagnose_process_plastic_friction_recovery.py --run outputs/process_plastic_friction_q7r_20261002/n3_i4 --output outputs/process_plastic_friction_q7r_20261002/event_diagnostic
```

Exit0, outer process 2.5708648 s, inner 2.1832167 s. All 37 actual branches
returned roots, all 37 actual feasibility outcomes matched independent vector
inequalities, and all 19 target arrays replayed byte-for-byte. Full accepted
state signatures before and after are equal. No commit was reached/performed.

The trace confirms a disagreement between the mode probe's tolerance and the
strict bisection's required lower state:

1. At probe .05000125 s, the first all-stick branch releases IDs19/20/22/25.
2. The coupled next probe passes ordinary feasibility. ID21 remains forced
   stick with margin -6.0954860572337566e-5 N: within the declared 9e-5 N, but
   negative. Its previous accepted margin was -6.167756083641507e-8 N.
3. At target .0625 s, ID21 violates the cone by 5.691232925791425 N. A strict
   event bracket is required by the unchanged implementation.
4. All 34 tested strict midpoints fail; ID21 is a strict blocker at every one.
   The first midpoint .05625 has margin -1.3768351255410423 N. The first
   midpoint satisfying ordinary tolerance is .05000152587890625, margin
   -7.443486602909388e-5 N; it still fails strict margin>=0.
5. Last midpoint .0500000000007276 still has ID21 margin
   -6.224623660955331e-8 N. IDs24/26/27 additionally violate only the exact
   zero bound near the accepted endpoint. None violates the ordinary force
   tolerance there, and none has a released flow blocker.

Inferred from observed strict results and unchanged source: lower remains .05,
upper becomes .0500000000007276, width=7.275985369759042e-13 s, but no feasible
lower candidate was ever found. Thus `lower is None` triggers the original
failure; it is not exhaustion of the 40-bisection width budget. These bounds
are explicitly a replay-trace inference, not locals stored by the original run.

## Evidence classification and next scope

- CONFIRMED: branch roots recover under budget14; mixed ordinary/strict
  feasibility leaves a tolerated over-cone forced-stick probe without any
  strict feasible bisection lower sample. No accepted progress is fabricated.
- SUPPORTED: this is an onset-at-the-existing-boundary active-set issue rather
  than an interior event that can be reached by more bisections. The trace
  samples support this, not a proof for every real time in the interval.
- RULED_OUT for this failure: new Newton nonconvergence, released-flow failure,
  contact-work acceptance, negative-dissipation rejection, or mutation of the
  prior accepted prefix. The replay never reached commit.
- NOT_RUN: a boundary-aware active-set repair and its independent regression;
  finer requested grid, dynamic AFSD, thermal and production.

Next work should distinguish an already-on-boundary onset from a true interior
event, with an independent coupled stick/slip benchmark and unchanged force/
slip/work acceptance. Merely raising Newton/bisection budgets cannot supply the
missing strict lower candidate. Do not globally relax the cone or accept the
5.691 N over-limit target. This audit makes no repair and authorizes no rerun.

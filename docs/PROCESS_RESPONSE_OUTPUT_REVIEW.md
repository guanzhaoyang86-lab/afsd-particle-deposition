# Independent readback of the saved-trial005 response experiment

2026-10-03. Numerical search improvement confirmed; physical event completion
is still BLOCKED. This audit loads existing output only. It does not call
advance, contact, material, root search, scout or production functions.

## Result and provenance

The one authorized experiment is
`outputs/process_response_20261003/saved_trial005_initial/`:
`PASS_PRIVATE_NUMERICAL_ROOT_NOT_EVENT_OR_PROCESS`, four evaluations,
three full search updates, core 0.2641362 s, external 3.5592058 s.
It starts from the published 6 us checkpoint, not a REJECTED candidate, and
uses precisely the old trial005 duration 2.7254762326048973e-7 s and the
stored two-contact incoming context. The attempted endpoint is
6.27254762326049 us, but it is NOT a newly accepted physical time.

Independent readback preregistration:
`docs/PROCESS_RESPONSE_OUTPUT_REVIEW_CONTRACT.md`. Actual command:

```powershell
.envs\pysph\python.exe -B logs/process_response_20261003/run_output_readback.py
```

One run, no repairs: 1218 readback assertions PASS, internal 1.662077 s,
external 2.0568106 s, exit 0. Results and all raw check errors/bounds are
`logs/process_response_20261003/output_readback_initial/report.json`;
the command log and supervisor receipt are alongside it. Counts include
individual fields, hashes and nodes, not 1218 independent physics benchmarks.

The audit verified 105 bound sources and four output artifacts. Parent and
candidate each contain exactly 40 ParticleArray properties. Every property
matches its canonical snapshot exactly in name, dtype, shape and bytes;
all values are finite. Independently implemented byte/semantic hash grammar
matches both states and the published report. Parent state:
`cbe85e8f8f77718bbbf393502efb709fd119f62c66946438d8106501fe2c8bb5`.
Candidate:
`a19de1e4eb6bbbb9943bd8912f8d4f462e960d939a4f6eb7e018a5a279ab876b`.

Mass, smoothing length, original positions, persistent/material/component IDs,
passive tracer labels, actuator flags, plastic strain and accumulated plastic/
return work are unchanged. Full stored last_interval impulses, work, raw Q,
quadrature, contacts and kinetic histories agree with the saved diagnostic.
All 4516 freeze entries were checked before and after: the only differences
were the already-authorized IMPLEMENTATION_PLAN and OUTPUT_SCHEMA index
updates being prepared by the main agent; frozen code/failure evidence did
not change. The actual experiment's own pre/post receipts reported no changes.

## What improved, and what did not

|Quantity|Old saved trial005|New PRIVATE experiment|Unchanged root gate|
|---|---:|---:|---:|
|Complete evaluations|18, failed|4, converged|24 maximum|
|Velocity residual (m/s)|3.833754755341113e-6|6.753457897691434e-14|1e-9|
|Position residual (m)|3.1607642588637583e-13|5.567932668886892e-21|1e-12|
|Relative density residual|2.3579492268187028e-15|3.3684988954552896e-16|1e-10|
|Temperature residual (K)|2.0463630789890885e-12|0|5e-12|

The independently reconstructed force-unit balance gives velocity residual
6.753009793940466e-14 m/s; its tiny difference from the saved value is due to
algebraic operation order, within the preregistered readback bound. No tested
solver function generated that reference.

The new full-step merits are 22206.1507255, 7.0713213063 and
6.7534578977e-5. The old failed run remains on disk and is not reclassified.
This is strong evidence that adding contact response to the search matrix
helps this particular finite-quadrature root. It is not proof that omitted
response was the unique original cause, nor that the continuous cross-zero
Jacobian exists. The counterexample in PROCESS_RESPONSE_LIMITS remains valid.

The preflight AST test compares the actual nonlinear evaluator against the
old evaluator; forces, material update, thermal solve and original gates are
unchanged. The only intended numerical change is the search matrix. The
recorded three contact-response matrix norms are 1.2604780759e-4,
1.1899360913e-3 and 1.1920675241e-3 kg. Their asymmetry norms are about
1.177e-7 kg; this asymmetry was not silently removed. Combined numerical
matrix condition estimates are 120.7775, 17.7514 and 17.7868.
Readback independently reassembled each complete response block as
`(1/(2*gamma)-alpha)*D_A F1 - D_A Jactual/(gamma*h)` and verified the
recorded matrix norms/asymmetry. This checks the full endpoint-plus-DeltaJ
algebra; it is not a new continuous-derivative certificate.

## Negative tail: no saved negative nodes does not mean passivity

All 18 contacts use 16 saved Gauss nodes. None of these saved nodes has
negative raw Q in this particular candidate. Their last fraction is
0.994700467495825, while both selected full-vector zeros lie near
0.998760642. Consequently the short post-zero tail contains no saved node.
Do not import the earlier trial003/004 negative-node counts into this run.

Readback independently forms the stored slip polynomial, locates its
longitudinal zero, checks the whole vector against the unchanged roundoff
bound, and integrates the actual declared PRIVATE signed power on each side
with 96/192 Gauss rules. Normal gap is reconstructed from the saved dense
particle x path and the physical 5 Hz wall program; cap is mu*kn*gap.
No backend callback or force function is called. Results:

|Contact|Independent zero fraction|Post-zero duration (s)|Signed post-zero Q (J)|
|---|---:|---:|---:|
|0:4|0.9987606420108932|3.3778407429995633e-10|-1.3038460258178836e-14|
|1:22|0.9987606420390307|3.3778406663114967e-10|-1.3038459665758984e-14|

The complete vector norms at these independently computed roots are
2.27e-16 and 6.34e-16 m/s, below the unchanged 1.13687e-13 m/s bound.
The final z slips are approximately +2.8167221e-5 m/s, plainly not zero.
The two durations/root fractions are not made identical by projection.

The signed tail is real within this PRIVATE mathematical path, not a
floating-point heat sign error. It is about 2.7e7 times each contact's saved
4.84e-22 J integrated rounding bound. Both tails remain nonphysical incoming-
direction extensions beyond the zero and cannot be published as sliding.

Importantly, the saved NET contact Q agrees with the split high-order total:
saved minus independent total is -6.45e-23 J and +1.04e-21 J. Thus the
stored nodes do not explicitly sample the negative tail, but the net
quadrature does not demonstrably omit or clip its signed contribution.
This happens because the nearly collinear continued traction/power is smooth
across the very short tail. It does not convert node-only passivity into a
whole-path proof. Existing `whole_path_passive=False` and
`physical_contact_path_certified=False` are necessary and correct.

## Work and thermal bookkeeping

Readback recomputed each saved node's `Q=-Jt dot ut`, material work and wall
work, reduced all contact rows, and checked normal-relative work explicitly.
Total raw contact Q is 1.2787169041232336e-6 J. This is signed PRIVATE
quadrature, not heat generated by a physically certified reversed interval.

Independent polynomial force/velocity integration and endpoint momentum
algebra confirm the kinetic identity. Recorded sums are:

|Term|J|
|---|---:|
|Delta kinetic energy|-2.908223428586501e-7|
|Physical endpoint/chord work plus actual contact correction|-4.5832456896442387e-7|
|Signed algorithmic filter work|+1.6750222610350537e-7|
|Root residual work|+2.2678570160070966e-18|
|Kinetic identity roundoff residual|+5.31906626858275e-22|

Algorithmic work is recorded separately, not fed into temperature or claimed
as thermodynamic dissipation. This bookkeeping identity is not a proof of
finite-strain material free-energy consistency or temporal accuracy.

The stored thermal ledger retains a 50/50 friction split of
6.393584520616169e-7 J each, contact heat into material 1.7170430420455705e-6 J,
internal pair-conduction sum 1.58705e-23 J, and zero plastic heat in this
explicitly elastic fixture. Its measured thermal increment is
2.356401501444339e-6 J versus input 2.3564014941071875e-6 J; the signed
7.337151632450852e-15 J discrepancy remains recorded. The thermal routine was
not rerun by this audit; these ledger values are not promoted into a new
independent material/thermal benchmark or erased by residual heat.

## Publication and final classification

- PASS: source/artifact binding, 40-field roundtrip, parent history, signed
  raw-node reductions, complete response algebra and kinetic bookkeeping.
- PASS, numerical only: the single unchanged-duration PRIVATE root now
  converges, with the accepted parent unchanged.
- FAIL if interpreted as a physically admissible entire interval: it contains
  a post-zero incoming-branch tail; the source correctly refuses this claim.
- BLOCKED: an accepted event-left state, chronological event handling,
  admissible right-limit restart and validated post-event continuation.
- NOT_RUN: another scout, actual AFSD continuation, process accuracy or
  production. This audit adds no physical step.

`PRIVATE_NONRESTARTABLE.npz` remains accepted=false/restart_allowed=false;
the experiment actually tested and logged rejection by the normal restore
entry point. This audit used raw serializer readback only, never bypassed that
guard to resume. No accepted_history or physical checkpoint was published.
The valid local continuation parent is still 6 us; the actual AFSD transport
state is unchanged. A 5.68e7-fold residual reduction is useful diagnostic
progress, not completion of the physical stop or AFSD reconstruction.

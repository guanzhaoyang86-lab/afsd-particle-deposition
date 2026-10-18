# Four-call bridge helper: independent interface contract

Registered 2026-10-03 before this package's interface tests. These tests use
NumPy arrays, fake evaluators and fake clocks only. They do not instantiate
ProcessThermalTrial, load the342-particle checkpoint, solve a material/contact
problem, run a scout or publish any state. The governing work package is
`PROCESS_RESPONSE_NEXT_WORK_PACKAGE.md`; no physical tolerance is changed.

## Independent fixed direction and perturbation references

The direction test uses the four positions (0,0,0), (2,0,0), (0,2,0), (0,0,2)
in m, masses (1,2,3,4)kg and persistent IDs (40,10,30,20). Their independently
computed mass center is (.4,.6,.8)m and global length is2m. Applying the
registered fixed matrix and masking the second/third particle z axes gives
the following literal, unnormalised reference, not produced by the helper:

```text
[[-.325,  .075,  .15625],
 [ .675,  .325,  0.     ],
 [-.075, -.425,  0.     ],
 [-.2,    .1375,-.34375]]
```

The normalisation divisor is `sqrt(.675**2+.325**2)`. Direction/reference and
translation/uniform-length-scale checks use absolute2e-14, relative0, in
dimensionless D. Reordering rows with IDs/mass/masks must reorder D exactly,
because the mass reduction uses persistent-ID order. Constrained axes are
exactly zero; input arrays stay unchanged. The interface has no passive-tag,
force, stress, or desired-outcome argument.

The candidate has exactly the original five fields. For literal epsilon
1e-7m/s and dt=5e-8s, expected changes are `sign*epsilon*D` in mean velocity,
twice that in endpoint velocity, and dt times that in position. Density and
temperature are copied exactly. All references start from the same original
candidate, including C- and the final C; no tested perturbation result supplies
another expected input. Comparisons of the explicitly constructed arrays,
unchanged inputs, scalar fields and zero-sign copies are exact. Returned
arrays must not alias the candidate or one another's results.

Invalid shapes, nonfinite arrays, nonpositive mass, duplicate/noninteger IDs,
nonbinary free axes, zero geometry extent and zero free direction reject.
Perturb rejects nonpositive/nonfinite/non-scalar dt, non-scalar/bool/invalid
sign, missing/extra fields, shape mismatch or nonfinite candidate/direction.
This pure helper does not replace the old physical density/temperature/domain
checks of the actual evaluator.

## Fake-clock budget and transaction references

The initial complete-call estimate is the recorded
4.04492670000036s, not a measurement from these fake tests. The independent
budget is `elapsed+1.5*remaining*max(old_cost,measured)+5 <= limit_seconds`.
The helper default remains50s and the original13 tests retain that limit.
Here each measured duration includes evaluation **and** the full result
callback (archive plus audit), not evaluation alone. The separate
`evaluation_seconds` journal remains available; the budget uses
`complete_measurement_seconds` for its estimate of future complete calls.
Tests cover the old-cost floor, the maximum rather than average measured
duration, the inclusive50s boundary, overhead/loaded-time inclusion, invalid
elapsed/durations/counts and a budget-blocked zero-call result. No physical
threshold or epsilon is adjusted after observing a result.

Fake `run_four` calls must be exactly C,C+,C-,C, at most four, using the same
base even when each evaluator's returned map is deliberately very different.
Every result callback must finish before the next call. Fake archive time
counts both in elapsed-budget checks and in the future-call forecast. A
measured slow call may stop the sequence
with the exact completed records and remaining NOT_RUN labels. Evaluator or
callback failure is terminal, stores the failing index/label/candidate and
prior records, and never retries. Mutating an input during either stage must
be detected; failures do not authorize another call. Nonmonotonic evaluation
timing rejects. The callback stores only detached test receipts.

The helper's completion status counts calls, not physical acceptance or a
certificate that the real full run met every source/history/budget gate. The real
driver owns source/config/archive/parent purity and its final budget/watchdog.
This fake suite does not claim those unexecuted real342 checks have passed.

## Evidence and boundaries

One initial interface invocation is externally limited to60s and writes fresh
`logs/process_bridge_20261003/fake_tests_initial*` evidence with source/test/
contract hashes. Preserve any initial failure; this reviewer does not repair
the helper or silently rerun it. No old frozen file is edited.

Source-review follow-up, registered before its execution: the original12-test
PASS and original three sources remain under `fake_tests_initial*` and
`fake_initial_sources`. They did not detect that the original helper forecast
only evaluation time rather than the complete callback-inclusive cost. The
source owner corrected that budget omission; the reviewer adds test13 with
fake evaluation1s plus archive/audit9s. After the first call the independent
forecast is `10+1.5*3*10+5=60s>50s`: exactly one result must be saved, the next
three must be NOT_RUN, `evaluation_seconds=1` and
`complete_measurement_seconds=10`. The old evaluation-only estimate would
incorrectly permit the next call. All these fake times and counts are checked
exactly. One fresh `fake_full_cost*` invocation, external60s maximum, checks
all13 tests; it does not replace or reclassify the original12-test evidence.

Explicit-budget follow-up, registered before its execution and before the
first real evaluation: the source owner has declared driver internal110s /
external120s, below the existing user600s cap, to include per-interval audit
cost. This is a compute-budget change, not a physical tolerance change, new
direction or additional evaluation. The pure helper default stays50s;
`may_continue` and `run_four` accept optional `limit_seconds` in(0,600].
Test14 repeats the literal1s-evaluation+9s-callback sequence: explicit50s
must stop after one call, while explicit110s must complete exactly four,
with elapsed40s and all complete-cost journals10s. Invalid601,0,NaN or bool
limits must reject before any fake evaluation. A fresh
`fake_budget_explicit*` invocation checks14 tests under its own external60s
test cap. The previous13-test logs and complete-cost50 source/test/contract
snapshots remain preserved and are not reclassified.

Branch-signature and independent mechanical/thermal ledger tests belong to
separate review tasks. This suite does not certify a smooth full residual
derivative, material accuracy, an accepted timestep or transport. Before
execution all interface results are NOT_RUN; fake-only PASS cannot substitute
for the independently audited real four-call response measurement.

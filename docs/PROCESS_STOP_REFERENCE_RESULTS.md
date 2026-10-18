# Independent stop-reference initial result

2026-10-03. The preregistration is `PROCESS_STOP_REFERENCE_DESIGN.md`. This
result file is separate so the tested contract/source receipt remains valid.
No AFSD physical step or production run was added.

## Actually run

- Initial independent unit suite: six PASS, one FAIL from a redundant exact
  floating-point energy equality (5.55e-17 J). Original test and log retained.
  Repair 1 applies the already registered 1e-12 SI bound; no model changed.
- Repaired independent suite: seven PASS.
- One actual paired-reference invocation: FAIL_EVENT_REFERENCE before any
  new state was accepted, `SCALAR_ROOT_NOT_FULL_VECTOR_ZERO`, 0.0252 s internal.
  The coarse branch failed its root certification, so fine extension was
  NOT_RUN. Neither branch reached 8 us. The old accepted prefix is unchanged.
- Additional regression for zero mean hiding opposite individual residuals:
  eight PASS. No tolerance or numerical-physics change was made.

## Exact first coarse-root evidence

The mean root candidate occurs at 5.974656588380895 us. Its two contact
particle IDs are 4 and 22. Their tangential relative velocities (m/s) are:

| ID | uy | uz |
|---|---:|---:|
| 4 | 1.0394320176046215e-18 | +5.247746681646959e-13 |
| 22 | 2.0715097423218773e-18 | -5.247746681646959e-13 |

Mean signed slip is exactly zero in floating point, but EACH full vector
exceeds the fixed bound 1.1368683772161603e-13 m/s. It is therefore incorrect
to certify simultaneous full-vector zero, even though the residual is tiny.
No velocity was assigned/projected to zero. All three root candidates remain
REJECTED_PRIVATE_INCOMING_BRANCH_CANDIDATE and are not restartable. Their
complete 531-entry states and metadata are retained in the diagnostic of
`logs/process_stop_20261003/reference_initial/report.json`.

The mismatch is consistent with accumulated mirror-asymmetric floating-point
error in the frozen coarse history. Using the local incoming acceleration
about 8.787e4 m/s2 gives an ESTIMATED individual-event separation of 1.19e-17 s,
larger than the registered 5e-19 s localization width. It does not justify
increasing the velocity bound, forcing joint timing, or changing the old
accepted history. Sequential individual roots could be a separately
preregistered extension, with real intervening force histories and no snap;
that extension is NOT_RUN here.

## Evidence and scope

Commands (each external60s bound):

```text
.envs/pysph/python.exe -B tests/review/test_hht_event_reference.py
.envs/pysph/python.exe -B scripts/reference_hht_event.py --parent logs/process_slow_20261003/reference_repair1 --output logs/process_stop_20261003/reference_initial
```

Logs: `reference_tests_initial.log`, `reference_tests_repair1.log`,
`reference_tests_mean_zero_regression.log`, `reference_initial.log`, all under
`logs/process_stop_20261003/`. Original test:
`logs/process_stop_20261003/test_hht_event_reference_initial.py`.

The independent scalar reversal, feasible-stick and nonzero-transverse
component-turn tests PASS against analytic literals. The 3D simultaneous
event reference is BLOCKED, not PASS. No new post-event reference interval is
certified and the failed attempt cannot authorize coarse HHT publication.
The numerical root/reference implementation is not general multiface contact
or actuator-constrained friction. General3D feasible sticking and a second
full-zero event remain unsupported. The RK dense full-slip polynomial checks
are exhaustive for that cubic; nonlinear gap/yield/domain checks are at RK
stages, endpoint and dense quarter-points, not a mathematical global bound.
1 ms, AFSD production, dt/space convergence and complete reproduction are
NOT_RUN / not established.

## Subsequent bounded repairs and final result

The initial FAIL above was not overwritten. Parent-authorized numerical
repair 1 is preregistered in `PROCESS_STOP_SEQUENTIAL_ADDENDUM.md` and selects
individual events by measured time, switches only their respective forces,
and preserves every physical state/ledger. The original oracle source is
`logs/process_stop_20261003/reference_hht_event_initial.py`. Twelve independent
unit tests passed. `reference_sequential1/report.json` retains the next FAIL:
both stops were located, but fine tripped a near-endpoint component-root
alias and had no post-event shared sample.

Read-only evidence in `sequential_alias_audit.json` demonstrated that this
fine candidate had strictly positive and increasing outgoing z-slip, not a
second full zero. Final numerical repair 2 was preregistered in
`PROCESS_STOP_ENDPOINT_ALIAS_ADDENDUM.md`: certify the cubic's outgoing
projection over the complete tiny interval using polynomial extrema and
the original velocity-bound divided by its minimum outgoing derivative.
No state, raw heat, physical law or tolerance changed. Repair1 source is
`logs/process_stop_20261003/reference_hht_event_sequential1.py`.

Final independent tests:14 PASS in0.078s. Final paired oracle:
`logs/process_stop_20261003/reference_endpoint2/report.json`,
PASS_EVENT_REFERENCE_PREFIX_NEXT_BRANCH_BLOCKED in1.8063s. The independently
certified common endpoint is **6.5us**,105 shared snapshots,9 new post-event
snapshots. Old0..5.9375us time/fine/coarse arrays remain bitwise identical.
Both actual first stops have full-vector norms below the original
1.1368683772161603e-13m/s bound, without assigning particle velocities.

| Event | Coarse time (s) | Fine time (s) | Absolute difference (s) |
|---|---:|---:|---:|
| PID4 | 5.974656588374923e-6 | 5.974656587951323e-6 | 4.2360e-16 |
| PID22 | 5.974656588386867e-6 | 5.974656587951857e-6 | 4.3501e-16 |

Each contact's static capacity is insufficient by about1.935596N; outgoing
sliding is force-driven. Event physical state and ledger jumps are bitwise
zero. Across the complete common prefix, maximum comparison errors are:

| Quantity | Measured difference | Original bound |
|---|---:|---:|
| Weighted RMS position | 5.342e-17m | 1e-11m |
| Weighted RMS velocity | 1.030e-10m/s | 2e-6m/s |
| Weighted RMS full stress | .001148Pa | 70Pa |
| Maximum temperature | 5.969e-12K | 1e-9K |
| Face force | 8.028e-7N | .007N |
| Contact impulse | 1.383e-14Ns | 1.4e-9Ns |
| Wall work | 1.006e-15J | 2e-9J |
| Raw friction heat | 4.022e-16J | 2e-9J |
| Reservoir heat | 9.504e-19J | 2e-10J |

The next rejected full-zero candidates differ: coarse PID4 at6.5657616us,
fine PID22 at6.5564031us. They are **not converged event-time results** and
must not extend certification past6.5us. Pointwise independent static checks
give requirements.49424N/.56798N versus capacities2.70526N/2.70558N: feasible
sticking is now relevant and remains unsupported in this3D oracle. Neither
candidate is accepted or restartable.8us remains BLOCKED;1ms is NOT_RUN.

Final readback `logs/process_stop_20261003/reference_final_readback.json`
PASS verifies all source/artifact/parent hashes, bitwise old prefix,
per-event full-vector zero, and no state/ledger jump. It does not rerun a
trajectory. Final source SHA256:
`34c8a4350785b2b09e231cfca054eae2c155abb0c391a0455e8f57d0e09bd7d2`.
Shared history SHA256:
`d461a79ab8a9197c37d627fa873092dd0a7b1e9ca28991ba131b23fa30836a8c`.

Re-run to a NEW output directory only:

```text
.envs/pysph/python.exe -B tests/review/test_hht_event_reference.py
.envs/pysph/python.exe -B scripts/reference_hht_event.py --parent logs/process_slow_20261003/reference_repair1 --output logs/process_stop_20261003/reference_recheck_new
```

The declared scope is an independent manufactured elastic/thermal/contact
reference through its common6.5us prefix, not AFSD advance, general sticking,
or complete paper reproduction. No further repair was made after repair2.

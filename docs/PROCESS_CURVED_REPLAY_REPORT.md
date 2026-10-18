# Curved-contact published interval replay

2026-10-03. **PASS_BITWISE_PUBLISHED_INTERVAL_REPLAY_NOT_PROCESS**.
This is the one invocation registered in PROCESS_CURVED_REPLAY_CONTRACT.
No new accepted state, AFSD step or transport result is published.

## Actual execution

```text
.envs/pysph/python.exe -B logs/process_curved_20261003/replay_published.py
```

Private report: `logs/process_curved_20261003/replay_initial/report.json`.
External command/stdout/stderr: `logs/process_curved_20261003/replay_command/`.
Exit0, one actual call to the new coupled `advance`,15 residual evaluations.
Actual solve0.878522s; complete internal audit1.000124s; external1.476402s,
within internal55s/external60s. No retry or numerical repair was made.

Loaded only published `outputs/process_curved_20261003/S0_initial/trajectory/
stage_012/accepted.npz` at4us and recomputed the already published stage013.
The exact stored dt is9.999999999999997e-7s; end time is
4.9999999999999996e-6s. It was not replaced with a rounded literal1e-6s.
Both endpoints lie inside the independent reference's certified prefix.

The original recorded stage013 solver-only time0.769345s excludes output;
it was not mislabeled as a measured whole-step cost. Conservatively, the
entire measured S0 trajectory8.671028s bounds one complete stage including
publication. The preregistered formula gives
`2*8.671028+5 = 22.342056s`, below50s, before execution.

## Preconditions and immutable evidence

- Current execution-source inventory equals the original S0 source hashes.
- Every original output artifact, checkpoint/journal hash and complete
  semantic state hash passes before replay.
- Preflight, independent reference and private saved-RED receipt hashes match
  the S0 report; reference physical inputs and reference artifacts are verified.
- Resolved parameters equal the original parameters; no physical array,
  mass, modulus, time scale, wall programme or certificate is changed.
- Stage013's parent hash is exactly stage012's state hash.
- After replay, all source hashes and all original S0 artifact hashes remain
  unchanged, including accepted_history, REJECTED files and report.

The replay writes only a private report/log. It does not call the publisher or
write another checkpoint/VTK, and does not overwrite either original state.

## Exact comparison result

All40 ParticleArray properties have matching dtype, shape and raw array bytes.
Every complete semantic-state field matches recursively, including numerical
acceleration, true endpoint force, physical snapshot, parameter/parent hashes,
new integration definition, last_dt and last_interval. Array comparison is
byte-level; scalar floats compare binary representations, including signed
zero. This is stronger than a tolerance or ordinary numerical-equality check.

The complete last_interval matches bitwise: actual curved contact impulse,
material work, DeltaJ/DeltaW, signed raw friction heat, all raw quadrature nodes,
and kinetic/filter/root accounting. Recomputed cumulative face impulses,
wall work, friction heat, reservoir heat and numerical-work ledger, starting
from the saved stage012 ledger, match stage013 bitwise.

Expected and recomputed semantic state hash:
`64c572c02f61a3706102245d0db6f0459ef73dbe4e90d2788ae47ab277ddb5aa`.

Root residuals reproduce exactly: velocity8.3257071e-10m/s,
position2.5185264e-16m, relative density3.3684989e-16 and temperature0K.
The initial in-memory accepted PA/state is also unchanged.

## What this does not prove

This validates deterministic restart of **one already published event-free
interval** in this environment and named numerical variant. It is not an
independent physical or temporal-accuracy oracle, a cross-platform bitwise
guarantee,20us convergence, positive plastic-work validation or restart through
stop/stick/slip/feature events. The old strict HHT replay NOT_RUN status is not
overridden. No rejected/unreferenced state was resumed. Project completion,
real AFSD transport and production remain outside this result.

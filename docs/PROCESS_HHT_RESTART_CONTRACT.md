# Bounded local HHT recomputation from a published halfway checkpoint

2026-10-03. Preregistered before execution. Status: NOT_RUN until the successful
same-source three-grid sweep and component preflight gates below are satisfied.
This probes the existing local 27-particle PySPH fixture only. It neither reads
nor resumes an AFSD checkpoint and does not perform an independent oracle solve.

`scripts/restart_hht_probe.py --sweep <successful-report.json> --preflight
<matching-preflight-report.json> --output <fresh-project-directory>` is fixed
to the 32-step grid's published `stage_016/accepted.npz`. It recomputes exactly
16 steps of 62.5 ns to 2 microseconds with the original post-switch boundary.
No step, physical parameter, nonlinear budget, tolerance, or control changes.

Before any advancement, require the complete same-source successful sweep,
its exact bound preflight receipt, the preflight test/source/artifact hashes,
current driver/backend hashes, configuration hash, every sweep artifact hash,
and the 32-entry published fine-grid history. Every checkpoint/integral path
must stay inside the sweep; fixed expected stage names exclude arbitrary or
unpublished checkpoints. A failed sweep, rejected/timeout marker, missing
publication, altered source/configuration, malformed ledger, or mismatched
NPZ hash rejects. The restart script, test and this contract must themselves
be included in the matching preflight's source inventory.

Use fine-grid complete trajectory wall time divided by 32 as the measured
whole-step cost. Launch only when `2*16*whole_step_seconds+5 <= 50 s`.
The worker has a hard 55 s deadline; its normal CLI supervisor imposes 60 s
externally. Do not increase either limit or silently continue a partial run.

Restore full physical properties and semantic HHT history at stage 16, including
the explicit right-limit control restart. Load cumulative integrals from that
same published stage. Preserve field dtypes and all acceleration/force/control
histories. Recompute each of stages 17..32 using the unchanged coupled stage,
whole-path/physical guards and integral accumulation. Compare its entire
ParticleArray hash, complete semantic HHT state hash, and every cumulative
integral with the corresponding original published stage BEFORE publication.
Acceptance is bitwise equality, not a new physical tolerance. Report maximum
absolute field/acceleration/force/ledger differences even if equality fails.

Each new stage is private until comparison, source verification, saved NPZ
readback and atomic `accepted_history.json` publication complete. Saved files
not present in that history are not restart points. Failure records the last
published stage/time/hash; it cannot promote a failed candidate. Source input
artifacts are checked again after completion. All original evidence is retained.

Unit tests exercise parser/receipt/history/ledger rejection and exact comparison
logic using isolated synthetic metadata only; they never call `advance` or run
a trajectory. Such tests are not successful recomputation evidence. A passing
restart establishes deterministic recomputation for this local fixture only;
positive-plastic-work accuracy and AFSD process continuation remain untested.

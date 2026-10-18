# Q7 fixed-budget first rejected branch diagnostic

Preregistered before execution: **NOT_RUN**. This diagnostic is not another
physical run and does not grant production restart. The original BLOCKED run,
its REJECTED candidate, and all old numerical limits remain unchanged.

Input: `outputs/process_plastic_friction_q7_20261001/n3_i4`. Verify every source
and artifact digest and the independent partial-history reader before using
its last ACCEPTED checkpoint. Reconstruct the accepted diagnostic state at
t=.05 s, step 6; never initialize material history from the REJECTED candidate.
The rejected file supplies only the observed branch's boolean released mask
and a reproduction comparison, not initial material state or Newton guess.

Call the exact existing `ConeBlock.branch` for the same target t=.0625 s and
same released mask, from the same accepted displacement initial guess. Compare
the original maximum_iterations=12 with a fixed maximum_iterations=14. This
allows at most two additional Newton corrections, chosen before this test;
do not raise this budget again to chase a root. All residual, condition,
line-search, friction, material and physical parameters remain unchanged.
The two calls are pure: no commit, no advance, no event search or timestep
subdivision, no additional accepted frame. Byte-fingerprint every raw accepted
particle property, reference property, displacement/stress/force/work and
metadata before and after each call. Reconstruct and compare the saved record.

External process timeout=60 s, expected <20 s from existing setup/branch costs.
The output must report each iteration, root criterion, per-contact cone margin,
forced-stick/released feasibility, flow opposition, all-contact natural-map
residual, previous-boundary/direction transition risks and raw candidate Q.
An equilibrium root is not a legal accepted contact interval and is not a
completed path. The accepted state must stay exactly at .05 s/step 6.

Use `scripts/diagnose_process_plastic_friction_q7.py --run <input> --output
<new-directory>` inside the same single-thread/60 s outer-process wrapper used
for Q7. All newly saved candidates are classified `REJECTED_DIAGNOSTIC` with
restart_allowed=false, even if their branch root converges.

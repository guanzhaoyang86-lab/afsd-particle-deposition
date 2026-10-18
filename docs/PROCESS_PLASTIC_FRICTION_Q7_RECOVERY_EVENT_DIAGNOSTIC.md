# Q7R pure event replay (2026-10-02, preregistered)

The first Q7R attempt is BLOCKED: `cone event not bracketed`; last ACCEPTED
state remains .05 s / step 6. This is one diagnostics-only replay, not another
physical path, not a repair and not a restart. No iteration budget exceeds 14.

Before reconstruction, run the independent Q7R reader and require its honest
PASS_PARTIAL_READBACK status. Reconstruct only the accepted checkpoint using
the existing restricted `restore_diagnostic`, then apply the declared Q7R
budget14 and metadata. Keep both `commit` and `advance` forbidden. Invoke the
unchanged `advance_cone(body,.0625)` directly only to observe its branch and
feasibility calls. A reached commit must stop at the forbidden stub; no
candidate may become accepted, and no physics, inequality or tolerance changes.

Wrap every actual `branch` invocation, recording time, mask, correction trace,
normal/tangent forces, relative slip and independently calculated cone margin,
flow-vector residual and allowances. Wrap the actual module `feasible_cone`
temporarily to log strict/default calls, return values and blocker IDs. Do not
replace its return value. Preserve root convergence versus contact feasibility.
Record full accepted-state signatures before/after and verify identical; compare
the target arrays bytewise with the saved REJECTED target. Reconstruct bisection
bounds from observed strict calls, explicitly labeling this as trace inference;
the old exception did not save those bounds.

Single process, single-thread, external timeout60 s; prior complete attempt was
10.198 s, so estimate <20 s. One replay only. Save full evidence under
`outputs/process_plastic_friction_q7r_20261002/event_diagnostic` with diagnostic
classification and restart_allowed=false. Source/input hashes are required.
All candidates remain rejected diagnostics; the last accepted state remains
.05 s / step6. Classify strict-start/finite-event hypotheses only after trace.

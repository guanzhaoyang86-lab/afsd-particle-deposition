# Onset first negative-normal diagnostic preregistration (2026-10-02)

The one actual n3_i4 run accepted step16 at manufactured time .1125 s, then
stopped at target .125 s with `negative normal outside fixed-contact cone scope`.
Its saved attempt trace has a successful BOUNDARY_PROBE at .11250125 and a
failing TARGET with only the last of nine tangent rows eligible. The original
exception has no candidate arrays. This document preregisters one diagnostics-
only reconstruction, not another physical path run or a solver/contact change.

## Frozen experiment

- Input: outputs/process_plastic_friction_onset_20261002/n3_i4, its report,
  history, raw diagnostic checkpoint, accepted0016 journal and attempt0016 trace.
- Require every source/artifact hash saved by the run to match. Construct the
  same OnsetBlock; read raw checkpoint explicitly for diagnostic material state.
  Do not call/weaken production restore, whose restart_allowed=false guard stays.
  Restore every recorded accepted field and require record()==saved final record;
  particle IDs/reference geometry/volume and journal must agree with history.
- Bind commit and advance to unconditional exceptions. Fingerprint every particle
  and reference property, block array (including shape/dtype/strides), material
  law and serializable metadata before and after. They must remain identical.
- Invoke only branch(.125, saved TARGET mask, maximum_iterations=14) once.
  No mode probe rerun, event search, added Newton budget or accepted step is allowed.
- Scoped wrappers call unchanged trial and cone_response. Capture raw trial
  arrays immediately before the first cone_response receiving normal<0; then
  call the original function, retaining its error and normal handling. Record
  branch source line plus iteration/line-search context rather than guessing
  from the exception string. No abs, force clipping, root change or heat repair.
- Save candidate arrays/IDs as REJECTED_DIAGNOSTIC_NOT_RESTART. Reconstruct the
  contact force assignment and gap arithmetically only for inspection; label
  these as invalid-candidate bookkeeping, not accepted equilibrated quantities.
  Do not evaluate a cone root outside its admitted normal-force domain.
- Compare captured N with actual internal_force[node,0], retain all signed forces,
  gap, displacement, stress, ep, plastic/return increments, and raw free-body
  mechanical balance. Negative normal at a nonconverged iterate does not prove
  physical separation. Fixed-contact geometry still prescribes zero normal gap.
- Single thread, external timeout60 s (expected<3 s); fresh diagnostic directory;
  no run retry or finer-grid path. Restore mismatch means diagnostic BLOCKED.

Expected outputs: diagnostic.json and REJECTED_first_negative_normal.npz in
outputs/process_plastic_friction_onset_20261002/negative_normal_diagnostic;
log in logs/process_plastic_friction_onset_20261002/negative_normal_diagnostic.log.
Any measured result is added separately after execution, not to this preregistration.

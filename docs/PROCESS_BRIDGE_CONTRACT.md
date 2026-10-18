# Real342 process bridge, registered before execution

2026-10-03. Implements the frozen PROCESS_RESPONSE_NEXT_WORK_PACKAGE, not a
new material/contact/geometry model. Preserve all4600 previous frozen files
except HANDOFF, IMPLEMENTATION_PLAN and OUTPUT_SCHEMA continuation indices.

Use the public650.393633ns nonzero-plasticity parent and the bound stage001
50ns accepted candidate. The latest850.393633ns history remains read-only.
The existing named ProcessThermalTrial BE map is used unchanged. Its thermal
map may differ from the old accepted explicit-thermal map; do not require
zero residual, old T-map equality or a new accepted state.

Exactly one supervised measurement process, at most C,C+,C-,C evaluations.
Direction, epsilon=1e-7m/s, dt=5e-8s, geometry/parameters and source history
are fixed as specified in NEXT_WORK_PACKAGE. Independent checks use mapped
velocities/temperature for map bookkeeping, not pretend C is equilibrated.
The pure helper cannot import or advance physical solvers. The real driver
reuses ProcessThermalTrial; no HHT derivative is transplanted into it.

Before launch: old245 plus new fake-interface regressions must pass; synthetic
audit checks must exercise signs/source-once and fail on altered records.
Retain all initial failures. Source/parent/archive hashes before and after;
normal output is PRIVATE, accepted=false/restart_allowed=false/can_publish=false.
No commit, checkpoint migration, transport event or physical publication.

Numerical accounting bounds are those in the prior NEXT_WORK_PACKAGE and
PROCESS_BRIDGE_AUDIT_CONTRACT, not new physical accuracy tolerances. Reentry
C must reproduce all numerical output exactly, except explicitly excluded
timing metadata. Physics/ID/history purity is exact. Record actual impulse,
constraint work, signed friction per subinterval, positive plastic work and
BE source ledger. Do not abs/clip Q or use numerical residual heat.

Report the complete five-field signed directional quotient plus contact,
constraint/material/thermal response, with SI units. The actual floating-point
C+-C- perturbations are retained. A changed contact/event/plastic/floor branch
signature is NONSMOOTH_NOT_CERTIFIED, not proof of a smooth Jacobian. One
epsilon is not an epsilon-convergence test; this is a response measurement.
Do not adjust epsilon or physical parameters after reading results. Event
signatures bind the actual contact_id/feature_id participants, not just event
names/counts. Exact fake metadata regressions must catch a face-A to face-B
switch while allowing event time movement. Both base calls must actually
cover positive plasticity, shared multi-contact and nonzero axial reaction/work;
the perturbed calls may change branches, which must be reported rather than
silently excluded. Zero-response records cannot grant mechanism coverage.

Prior full map cost4.0449267s predicts29.26956s for four calls including margin.
Before any real evaluation, the complete-cost source review showed that the
old cost sample excluded the new per-interval audit and BE certificate. The
measurement budget is therefore explicitly revised to internal110/external120s
(still below the user's600s cap), not inherited from the old design's50/60s.
The helper default remains50s for its already-tested generic budget behavior;
the driver passes110s explicitly. This is a pre-run compute-budget decision,
not a changed physical tolerance, extra evaluation or second attempt.
Update the remaining forecast after every complete evaluation and archive/check;
skip remaining calls when it exceeds110s. Initial29.27s is a lower-confidence
estimate; no claim is made that it includes all new audit cost.
Stop after domain/audit/purity failure, retain completed evidence, no retries
or second process to complete the four-call allowance. Whole-project limits
of600s, no GPU/paid/long/production run still apply.

Completion of this package means an actually measured real342 interface and
independent bookkeeping, or clear partial/failure evidence. It does not mean
new physical time, transport, process-scale accuracy/cost feasibility or
completion of AFSD. Those remain required by the active overall goal.

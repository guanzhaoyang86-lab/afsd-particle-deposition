# Next bounded package: process-timescale inertial accuracy/cost feasibility

2026-10-03. DESIGN / NOT_IMPLEMENTED / NOT_RUN. Read PROCESS_HHT_REPORT and
OUTPUT_REVIEW first. No permission for AFSD continuation, production, a new
backend, GPU/parallel optimization, quasistatic replacement or long paid work.

## Do not repeat the finished experiment

The registered8/16/32 HHT grids all reached2us. Velocity converges with observed
orders1.914/1.967 but finest2.052004835e-4m/s exceeds2e-4. All other8absolute
bounds and9refinement ratios pass. Root residual<=4.53e-10m/s and reference
disagreement3.62e-11m/s cannot explain the5.20e-6m/s excess. No temporal-map
bug was independently confirmed; no further patch, h refinement or relaxed
tolerance is justified solely to turn that FAIL green. The gated half-run
restart stays NOT_RUN. Preserve every source-before/rejected/log artifact.

## The actual missing decision

The current manufactured normal wall completes a cosine cycle in2us, i.e.
500kHz. Its required temporal resolution does not establish the economics
of300rpm (5Hz) physical forcing. Nor does a stable implicit solve prove that
unresolved physical startup/contact transients are harmless. The old AFSD
state is genuinely unbalanced; its velocity/stress must not be cleared.

Next develop ONE slow-forcing, real-inertia3D contact/thermal feasibility
benchmark, including its full physical preload/startup and a substantially
longer slow-load window. Reuse the same PySPH/corrected/FKS/BE/HHT components,
physical masses/moduli and force-driven free DOFs. Give this distinct
manufactured load history its own name/provenance, not a relabeling or time
scaling of the completed fixture or of the actual AFSD checkpoint.

Before implementation/execution, register exact geometry/material/initial
state, applied load/rigid program, observation interval and fixed candidate
step schedules, independent reference construction and raw accuracy limits.
Resolve full early transients in the independent reference and preserve
velocity/reaction/work/Q/T checks; do not compare only a slow average or
drop an inconvenient startup interval. Any algebraic modal/linear reference
must specify its validity against the actual current-geometry material map;
calling a finer copy of the tested HHT step is not independent validation.
Check conditioning/no pseudoinverse and physical energy terms separately.

Measure entire build/NNPS/material/contact/thermal/root/output cost. Establish
the largest step/schedule meeting the predeclared physical error bounds, then
report whether an accuracy/cost overlap exists. This is not another campaign
of equivalent62.5ns windows or a promise that all elastic acoustic modes need
to be temporally resolved over the complete AFSD process.

Stop and report if opening, stop/stick, shared-face/actuator constraints,
yielding or geometry exits require an unimplemented branch. Do not keep the
contact active by clipping gap, assigning a free velocity, faking sticking,
changing physical loads, or hiding that new dependency. That failure is an
actionable feature requirement, not permission to remove inertia/constraints.

Each independent test/reference <=60s external; any later actual AFSD trial
needs its own measured estimate and remains<=600s. At most two evidence-based
implementation repairs in the named experiment; keep initial failures and
stop rather than increase budgets/tune thresholds. No run starts if its
registered wall-time forecast exceeds its cap. NOT_RUN is a valid outcome.

## Dependencies after this cost decision, not completed by it

1. Positive-plastic-work FKS/material temporal accuracy, thermal feedback and
   legacy Jaumann energy limitations under a specified applicable strain range.
2. Located contact opening/stop/slip/feature events and compatible right-limit
   acceleration with shared-particle faces/actuator reactions and work.
3. Existing finite_tool_geometry/normal_penalty and inlet_protocol integration
   into accepted transactions, real shoulder crossing and persistent IDs.
4. Independently recomputed AFSD restart, same-end-time and same-geometry
   resolution tests, observable transport and appropriately scoped comparisons.

A separate equilibrated hot-start/quasistatic process variant may change the
cost question, but it changes the physical start/model. It requires explicit
user direction and its own initialized result, never reinterpretation of the
current850.393633ns AFSD checkpoint. Project remains NOT_COMPLETE.

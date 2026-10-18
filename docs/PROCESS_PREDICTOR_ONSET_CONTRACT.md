# Predictor + Onset integration preregistration (2026-10-02)

Scope: only integrate the already tested pure feasible predictor with original
Onset events and Q7 atomic commit. Preserve old source/config/tests/FAIL records.
Initial archive PROCESS_FEASIBLE_PREDICTOR_FREEZE.json verified1567 files.
Only previous HANDOFF/IMPLEMENTATION_PLAN/OUTPUT_SCHEMA may be updated.

New PredictorOnsetBlock inherits OnsetBlock; ONLY branch changes to the previously
verified solve_with_feasible_predictor. No new contact/material model or opening.
Onset's boundary eligibility, strict interior bracketing, normal guard, work
formula, signed dissipation and final commit unchanged.
Predictor+cone still share14 corrections; line search8, bisection40, eventlimit12,
force9e-5N/root1e-8/slip1e-11m/energy1e-6J all unchanged. No force/temperature/
yield/negative work clipping, residual heat, scaling or changes to tags.

## Gate sequence

1. Independent linear spring predictor + exact event reference and guarded commit.
   Reference/tolerance preregistered separately in PROCESS_PREDICTOR_ONSET_REFERENCES.
   Old original branch's negative initial normal counterexample is retained RED.
2. Actual integration tests: class identity for unchanged commit/advance/trial/
   tangent/state/wall, same initial state/config apart from named metadata;
   actual small .025 preload with events and physical history, passive tags;
   rejected untrusted force/displacement candidate cannot change accepted state;
   budget/corner guards. No full path hidden inside unit tests.
3. Reader: require parent16+Onset4+new3 source inventory, exact config difference,
   and stage-aware shared correction budget (predictor rows aren't contact roots).
   Accepted-only physics data still audited by original Q7 reader.
4. Only after independent+integration tests pass: ONE full4-increment attempt,
   original physical1cm/27particle/9fixednormal J2+friction block, four phases,
  16 requested intervals to common.2s. This is a manufactured quasi-static program,
   not AFSD deposition time. No restart from a rejected/diagnostic checkpoint.
   A failure is preserved and stops this turn's path attempts; no repeated tuning.
5. Read accepted/rejected outputs independently. 8 increments/space/AFSD remain
   NOT_RUN this turn, even if4 completes. Full historical failure suite NOT_RUN,
   old failures remain hash-preserved; selected compatibility tests may rerun.

## Cost and limits fixed before run

Original Onset got .1125/16 accepted intervals and265 branch roots in45.687s.
Linear physical-program extrapolation .2/.1125 gives81.2s, but extra stick/slip
events make it uncertain. Previous predictor correction target took<.7s including
two diagnostic branches; it does not estimate all events. Plan90–150s for one full
attempt, external timeout180s (<user600s), singlethread. No performance tuning.
This is an explicit wall-clock harness limit for NEW variant, not a numerical/
physical time change; frozen parent config60s is historical and unchanged.
Tests/readers each external60s. Outer processes must be stoppable and report
status at <=60s intervals while running. On timeout prior accepted journals
survive; incomplete outputs must not be labelled complete/readbackPASS.

## Output and outcome separation

Reuse original Q7 journal/history/VTK/checkpoint and Onset attempt trace.
New predictor_provenance.json stores new wrapper/config/helper hashes +Onset4+
parent16; resolved_config.predictor_variant names Q7_PREDICTOR_ONSET_V1.
Trace stages INTERIOR_PREDICTOR and CONE_ROOT separate numerical initialization
from the contact root; total correction count stays<=14. No predictor state is
accepted or counted as a physical step. Full paths keep original status strings
J2_CONE_PATH_COMPLETE/BLOCKED, reader PASS_COMPLETE/PASS_PARTIAL_READBACK.
Partial never equals convergence; accepted data must satisfy unchanged cone,
work/material/ID checks. Candidate fields remain REJECTED, no fabricated heat.
Material Ws-Es-P warning and prior Jaumann/space issues remain unresolved.

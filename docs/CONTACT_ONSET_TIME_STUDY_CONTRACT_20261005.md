# One bounded contact-onset temporal study

The independent 918-particle branch remains accepted only through 62.5 us,
global step 14, checkpoint SHA
12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b.
The repeated 7.8125 versus 2 x 3.90625 us comparison failed at PID797 even
after shortening the window and using the exact same parent. Its local stress
and plastic-strain differences were 16.4276% and 30.5838%. No tolerance is
changed and no failed state is promoted. This study is a finite method review,
not an automatic halving campaign or permission to continue indefinitely.

All trajectories start from that exact accepted parent, with its complete
14-step energy prefix, contact and material histories. The existing two fully
qualified but unaccepted fine candidates in expanded_hot_campaign window_001
are the 3.90625 us reference. They are never used as restart parents. The only
new trajectories allowed are four steps of 1.953125 us and, if the cost and
physical gates permit, eight steps of 0.9765625 us. All end at 70.3125 us.
The eight-step branch starts again from the accepted 62.5 us state. There is
no fourth refinement level. A failed physical gate stops that branch; a failed
temporal comparison is retained while the registered trajectory may finish,
because it is diagnostic throughout.

Use the unchanged frozen expanded_hot_continuation_candidate implementation:
full nonlinear equilibrium, fresh coupled force, material/thermal feedback,
scalar/batched equivalence, raw and complete-prefix energy, penetration,
continuous contact shell and retained-edge visibility path. Preserve its
existing bounded Broyden and one same-parent Newton fallback. No inverse-solve
optimization is introduced in this study. No physical constants, prescribed
axes, labels, contact laws, low-rate material assumptions or observers change.

At every matching endpoint compare against the previous resolution using the
original global and one-hop event-local field tests and exact crossing IDs.
Also store local absolute stress/strain differences, denominator floors, the
original failing comparison and ratios of terminal errors between successive
levels. Report empirical error reduction only; no assumed order or Richardson
extrapolated field is substituted for actual fields. Two successful adjacent
comparisons do not prove whole-history or spatial convergence.

Every checkpoint has accepted=false and restart_allowed=false. Acceptance time
advanced is always zero. A physically completed trajectory is distinguished
from passing all temporal screens. Neither outcome authorizes a process restart
without a separately reviewed publication procedure and independent evidence.

Before each stage, write an immutable protocol binding sources, exact parent,
the independent failed-window readback, all saved reference artifacts, step
count and cost formula. Measured per-candidate cost is the maximum of all
complete candidate times and whole-stage wall/count from the original window,
the previous parent window and (for level eight) the new four-step stage,
including unpublished candidates. Add one second per candidate. Require
30 + 1.5 * count * per_candidate_cost <= 450 s. Each stage is single-threaded,
internally capped at 540 s and externally at 600 s. If the eight-step forecast
does not fit, record BLOCKED and do not run it; do not shorten the study after
the fact. Outputs must be new. No old source, state, failed artifact or contract
is overwritten. Independent saved-only readback is required before using the
result to choose any further process work.

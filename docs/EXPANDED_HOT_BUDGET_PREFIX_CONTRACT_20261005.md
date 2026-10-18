# Seven-pair budget-selected fine prefix of an eight-pair coarse window

This contract reduces the planned work before any fine step is run. It changes
neither physics nor tolerances. The original continuation coarse stage completed
all eight 7.8125 us candidates from the real 15.625 us / global-step-2 parent to
78.125 us. Every coarse candidate remains unaccepted and nonrestartable. The
measured largest complete candidate cost raised the original sixteen-step fine
forecast to about 483.7 s, above the registered 450 s scheduling cap. Therefore
the original fine directory remains absent and its status is NOT_RUN. Its full
eight-pair objective has not been completed.

The new controller writes only
`outputs/expanded_hot_budget_prefix_20261005/window_001/fine` and selects exactly
the first seven already qualified coarse references. It computes fourteen
3.90625 us fine increments from the same real parent, targeting 70.3125 us /
global step 16. The eighth original coarse candidate is explicitly recorded as
`UNUSED_NOT_ACCEPTED`; it is neither deleted nor included in the accepted path.
No coarse candidate is recomputed, and no larger original stage is relabelled
as successful by shortening its history after execution.

Parent checkpoint, parent readback and branch-origin hashes are the same as in
`EXPANDED_HOT_CONTINUATION_CONTRACT_20261005.md`. The full two-record initial
energy prefix, all physical fields, identities, free axes, observer eligibility
and cumulative local-screen seeds are retained. Before selecting the seven
references, verify all original eight candidates' source maps, paths, hashes,
immediate-parent chain, exact times/steps, complete energy prefixes and their
successful externally supervised coarse manifest. Verify the old coarse map
against its actual original controller/contract; it cannot contain this later
controller. The new fine source map is that complete original map plus this
new controller and contract. Bind all eight original candidate artifacts as
immutable reference inputs even though only seven are compared with fine.

Reuse the frozen `expanded_hot_continuation_candidate` and acceptance helper.
Every solve still uses real neighborhoods, full material/contact/thermal maps,
fresh force and unchanged raw/all-prefix energy gates, scalar/batched identity,
finite-shell and retained-support path certificates. Every fine pair retains
the original same-parent global/local comparisons and strict crossing-ID check.
No interval is reset; no backward-Euler work defect is subtracted or compensated.
Only a passing pair permits publication, with its true direct-parent and window
hashes. A failed/incomplete pair remains nonrestartable and its evidence stays.

Use the original formula with all eight measured coarse costs:
`per_step = max(startup_complete_candidate_costs, all_eight_coarse_costs) + 1 s`;
`forecast = 30 s + 1.5 * fine_step_count * per_step`.
Register both the rejected sixteen-step forecast and selected fourteen-step
forecast before this new worker runs. Selection is permitted only if the former
exceeds 450 s and the latter does not. Do not lower the measured cost or omit
the expensive last reference to pass the cap. Internal/external limits remain
540/600 s, single-threaded. Neither a stage overrun nor an energy failure permits
automatic shortening, threshold relaxation, hidden time reset or another run.

Successful publication status is
`PASS_EXPANDED_BUDGET_SELECTED_FINE_PREFIX_NOT_FULL_WINDOW`. It means only the
registered seven-pair budget prefix passed. Record `original_full_window_completed=false`
and the unused eighth reference in every terminal outcome. A new independent
reader must audit all eight old coarse references and only the fourteen actual
fine candidates/accepted states, preserving the two inherited energy records.
The raw eight-pair fine stage remains NOT_RUN. Future continuation from a newly
audited true endpoint requires a separately registered time-step plan. This is
still an engineering numerical screen, not completed deposition or validation.

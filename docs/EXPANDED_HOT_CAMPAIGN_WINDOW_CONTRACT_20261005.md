# Explicit, bounded single windows for the independent 918-particle model

This controller replaces repeated controller copies with one parameterized
single-window protocol. It does not create an unattended campaign, change the
physical model, or authorize production. Each invocation runs one coarse or one
fine stage. A later window requires a new real accepted parent, independent
readback and explicit command. A repeated failure of the same local comparison
requires method review, not an automatic sequence of smaller steps.

Coarse registration requires explicit parent checkpoint, parent readback,
coarse step and requested pair-count upper bound from 1 through 8, plus a new
output directory under project outputs. Allowed coarse steps are 7.8125,
15.625, 31.25, 62.5 and 125 us. Fine is exactly half the registered coarse step.
A coarse step may not exceed four times the last actually accepted fine step
(equivalently twice the previous window's coarse step). This growth bound is
permission to test, not evidence that the larger step is accurate or affordable.
No physical parameter, material/contact/thermal law, mass/time/modulus scale,
energy tolerance, clipping rule, free axis, observer or passive label changes.

The parent must be the exact last published checkpoint in a successful saved-
evidence readback from the budget-prefix reader or this campaign-window reader.
Its window may have failed after a real accepted prefix. Require true accepted
and restart flags, the independent 918-particle class/scope and original initial
SHA, exact endpoint path/hash/time/global step, unchanged initial qualification,
all source hashes, the full global energy ledger, and the actual passing pair
screen bound by that checkpoint. Later failed screens must never supply its
seed set or be treated as accepted. Preserve their saved evidence and outcomes.
Neither a status string alone nor a coarse/nonrestartable candidate is a parent.

Before any coarse worker, write a new immutable `window_protocol.json`. It binds
all inputs, parent readback evidence, the parent outcome (including failures),
the last accepted fine step, inherited passing-screen seeds, proposed step,
requested and budget-selected pair counts, source map and scheduling formula.
The source map is the parent's implementation map plus this controller/contract,
the frozen continuation candidate/helper, true parent/readback, and the parent
ledger/screen/history/manifest/supervision files used to establish ancestry and
cost. Parent failure evidence hashes are additionally bound in the protocol.
The protocol's own hash is bound separately, avoiding a circular source hash.

Use complete candidate costs from both stages of the latest audited window,
including complete unaccepted candidates. For each stage also use its whole
wall time divided by its complete candidate count, so failed attempts and I/O
are not discarded. Set measured cost to the maximum of these values and add
one second of allowance per candidate. Forecast a stage as
`30 + 1.5 * planned_candidate_count * per_candidate_cost` seconds. Before coarse,
choose the largest pair count no greater than the explicit upper bound for which
both coarse and doubled fine forecasts are <=450 s. If none fits, record a
budget-blocked protocol and run no worker. A --plan-only coarse registration may
write this exact protocol without physics; a later command consumes it unchanged.

After coarse terminates, verify every fully qualified reference, its complete
source/field/artifact chain and full energy prefix; keep all failed/partial
evidence. Before any fine worker write immutable `fine_selection.json`. Recompute
cost using all new coarse complete candidates and whole-stage amortized cost,
never dropping an expensive unused candidate. Choose a budget-fitting contiguous
qualified prefix, preserving the original target/count and explicit unused
reference hashes. If none fits or none qualified, write a blocked selection and
run no fine worker. The same physical parent and time origin are mandatory for
both branches; fine cannot restart from a later coarse endpoint. No shortening
after execution converts a failed run into a completed original window.

The frozen continuation candidate/helper supplies the unchanged fresh complete
force, raw energy and all-prefix energy, scalar/batched, penetration, thermal,
finite-shell and retained-support path gates. Each fine pair additionally uses
the unchanged global/local comparison with the fine endpoint as reference,
inherited cumulative event/blocked/crossing seeds and one-hop raw/retained graph
neighborhoods. Only passing pairs publish both fine halves. All initial energy
records remain present; no backward-Euler residual is cleared, subtracted, or
converted to compensating heat. Tool pose, real fields, original failures and
unused coarse trajectories retain their actual identities and status.

`QUALIFIED_EXPANDED_CAMPAIGN_COARSE_NOT_ACCEPTED` is diagnostic only.
`PASS_EXPANDED_CAMPAIGN_SELECTED_FINE_SCOPE_NOT_PHYSICAL_VALIDATION` means the
preselected fine scope passed. Requested-target, registered-window and selected-
scope completion are separate booleans. Any accepted prefix of a failed stage
remains such until an independent saved-evidence reader audits it. Published
states bind the window protocol and fine selection hashes in addition to their
true direct-parent, branch-origin, window-origin and pair-screen hashes.

Each stage is single-threaded, forecast <=450 s, internally limited to 540 s
and externally to 600 s. No retry/next-window loop is provided. An optional
planning-only command never calls NNPS or a physical map. A new independently
implemented reader must examine both branches, all discarded references,
protocol/selection arithmetic, all published fields and the complete inherited
ledger before any further window is scheduled. These windows remain engineering
numerical screens, not sustained deposition, spatial convergence, three-case
production or quantitative reproduction of the paper.

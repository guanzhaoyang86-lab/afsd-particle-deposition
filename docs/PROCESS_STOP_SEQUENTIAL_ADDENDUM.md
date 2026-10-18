# Sequential measured events — numerical repair 1

Registered 2026-10-03 BEFORE the changed oracle, unit tests or trajectory.
The original source, rejected root states and report are retained under
`logs/process_stop_20261003/`. The original numerical bounds and physical
inputs remain unchanged. The initial test-assertion correction was not a
numerical model repair; this is the first oracle algorithm repair.

The coarse history has distinct individual roots: its zero mean contains
opposite 5.24774668e-13 m/s residuals exceeding the frozen 512eps full-vector
bound. Therefore simultaneous event grouping is not justified. The repair
measures both individual full-vector roots using private incoming-branch
candidates and selects the earliest measured time, not a particle-ID order.
Each root retains the same <=5e-19 s time bracket (or original floating-time
floor) and <=512eps*max(1,speed) full-vector check. Root brackets must establish
order. At most the two physically discovered mirror contacts may switch.

At the first event switch ONLY that contact after its own static disk check;
state x/v/S/rho/T and ledger are continuous. Locate the still-incoming contact
again from this unchanged event state, using actual changed forces, not the
obsolete second candidate. A positive incoming signed slip may already be
inside the root roundoff bound: that does NOT silently synchronize events;
it still receives its own time root and full-vector certification. No
velocity snap, event-time prescription, conditional simultaneous switch or
threshold increase is permitted.

A previously located contact can have a genuinely positive outgoing speed
smaller than the zero-localization roundoff bound between these two events.
When its actual vector has positive dot product with its measured outgoing
direction and nonzero norm, evaluate the actual normalized Coulomb vector,
not a zero-speed exception or a fixed force. The right-limit direction is
still permitted only inside the original event-time/norm roundoff window.
Later full-zero events remain guarded by dense full-vector checks; this is
not a general sticking/contact-event solver. In private root probes only
the unprocessed incoming branches are formally continued beyond their zero.

Record both individual events, each physical state/ledger before and after,
static feasibility, left/right forces, rejected private probes and measured
order. Compare corresponding identity event times across coarse/fine with
the original 1e-12 s bound. All shared-state uncertainty limits are unchanged.
The source/contract/parent receipts bind both the original design and this
addendum. Add independent close-separated-event/order and no-false-zero
regressions before a single rerun to8us or the next unsupported branch.
Prediction remains <20s, combined internal55s/external60s. Another unsupported
branch will be retained rather than repeatedly refined to obtain PASS.

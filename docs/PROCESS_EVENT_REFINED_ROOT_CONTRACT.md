# Bounded derivative refinement root strategy — before new cloud

2026-10-03. Original sources, tests, contracts and response_probe failure stay
unchanged. This is one numerical repair, not a change of physical time step,
contact law, solver acceptance tolerance or a continuation of a rejected state.
See the separately preregistered PROCESS_EVENT_REFINEMENT_CONTRACT for component
rules. The extra wrapper can halve force-difference h at most twice, with all
original derivative/branch/raw-heat/replay requirements. It must reject an
unresolved cusp; it never substitutes a free-particle contact response.

The fixed six-level read-only PID1 diagnostic distinguishes O(h) one-sided
truncation from a visible nonsmooth response: Cyx side gaps halve from2.27405e-4
to7.10641e-6 s/kg while central errors decrease quadratically. All13 kicks
have122 intervals, no stop events and minimum slip>=.0038149 m/s. Its13 local
kicks took.511 s; this is a cost sample, not a universal performance bound.
The anchor has96 contacted particles (88 one-face,8 two-face) plus246 free.

All139 old/new tests, including the retained continuous-reference FAIL, run
again; new refinement and strategy-wiring tests must pass. Retained FAIL remains
FAIL, not all-green. New named driver installs only response_from_trace and
records that strategy and additional source hashes in the existing driver.
The physical ProcessThermalTrial, initial predictor, event kick, raw tolerances,
32-evaluation root budget,500ns step and original accepted parent are unchanged.

One fresh refined cloud only: at most39 local kicks/contacted particle. PID1's
sample suggests roughly147 s plus evidence IO if every contacted particle
needed all3 attempts. Actual contact heterogeneity makes this uncertain; keep
the original180-s internal/210-s external hard cap. Incomplete cloud is rejected,
cost/fullroot NOT_RUN. No further h-halving/threshold changes after this attempt.

If complete: same-source cost probe<=80/90s, then only if measured forecast
1.5*(build+32*trial_proposal_IO)+30<=450s, one fullroot<=540/570s. No new time step
or extended iteration budget. All source/parent/anchor/state rules from
PROCESS_EVENT_ROOT_CONTRACT remain mandatory. Accepted data still needs the
unchanged full acceptance and independent readback; this is not process accuracy
or transport completion even if a numerical root succeeds.

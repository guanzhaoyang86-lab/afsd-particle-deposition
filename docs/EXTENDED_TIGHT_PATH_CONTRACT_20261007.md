# Extended process-time continuation, registered 2026-10-07

The user requested a longer physical simulation to observe more particle passage
and deposition. After the original 51.1875 ms path-budget stop, four separately
registered 0.25 ms fine windows all passed, including independent saved readback.
The true parent is now 823 accepted steps at 52.1875 ms, with six eligible outlet
crossings (2.025 mg). It is not a completed-deposition result.

Only this full, audited checkpoint may start the new series:
`outputs/tight_path_refinement_20261007/window_004/fine/step_032/accepted.npz`
SHA256 `861bf0120c6321d638d47bc80f0cab57ead09eb16380f9978ff64ee27667f93e`.
Its `window_004/readback/summary.json` SHA256 is `7edea3183aa9a985de655ad1a311a20244d372bfdab2be8657f978c3472ab149`.
The parent screen is `fine/pair_016_temporal_screen.json`. All 823 energy records,
particle IDs/fields, material/contact history, free constraints, crossing events
and original branch identity are preserved, with no resetting or inserted mass.

The controller is a static adaptation of `run_tight_path_process_v2.py`.
The previous qualified 62.5/31.25 us coarse/fine schedule is attempted again only
through the unchanged complete candidate, local/global temporal screens, and
independent readback. This schedule is not assumed stable at the new state.
Each new window is 16 coarse and 32 fine steps, at most 1 ms. At most 48 windows
are registered, targeting 100.1875 ms. Every next window requires the immediately
previous complete audited window. A failure stops the series; there is no blind
retry, timestep change, tolerance relaxation, failed-state promotion or skipped
window. A failed coarse run can supply only its already qualified contiguous
prefix for independent fine comparison; the original failure remains unchanged.

All original force/material/thermal inputs, geometry, initial finite inventory,
256-query/32-depth geometric certificate, fresh/scalar equivalence, continuous
shell, full-prefix raw energy and temporal thresholds remain unchanged. Local
seeds, one-hop comparison, identical crossing IDs and the original first-window
TOP observation requirement remain. The original 342-particle frozen PySPH/CPU
backend is used without scaling. Wider-domain or continuous-inflow experiments
are not substituted or appended to this physical history.

Forecasts use measured parent full-stage/root cost and must remain <=450 s;
per-stage internal/external limits stay 540/600 s. Numerical threads remain 1.
The separate orchestration supervisor has a maximum five-hour wall budget,
630 s wrapper timeout for the original 600 s physical supervisor, 600 s for saved
readback, cleanup of its own process tree, and a STOP_REQUESTED boundary stop.
It refuses to overwrite or silently resume existing stage outputs. It records
commands, PIDs, elapsed time, exit codes and the last independently audited
endpoint. A running target is not a completed result or deposition validation.

A separate static adaptation of the saved-evidence reader checks the new exact
parent, sources, time grid, all particle fields, geometric certificates, raw
energy and temporal comparisons. Original physics and prior component source
hashes must still match. The new reader's parent wiring is verified before
physical execution. No claim of rerunning the full historical test suite is made.

The 4 mm substrate remains smaller than the 5.4 mm shoulder and the bar remains
a finite inventory. Existing engineering/material/spatial/experimental limitations
and every prior failed result remain. These limits prevent promising that merely
reaching the time target will establish sustained or validated deposition.

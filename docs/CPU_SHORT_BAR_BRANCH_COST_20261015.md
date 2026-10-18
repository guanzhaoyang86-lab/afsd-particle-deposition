# Explicit CPU continuation using measured full-branch costs

The user requested continuation after the previous campaign stopped on a cost
forecast. This new administrative campaign starts at the complete, screened
75.1875 ms / step 400 short-bar state (continuation_004, SHA
5f96360bcb151dd91281d630136c3859d8732afb07299a43e88ca92c1413a3b6).
It repeats the preceding five-hour duration: 2026-10-15 20:10:09 through
2026-10-16 01:10:09 Asia/Shanghai, with preparation included. It can end earlier
at the unchanged 220 ms stationary-phase target, user pause or a genuine failure.
No paid compute, GPU switch, new particles or initial-state reset is involved.

The old run and its single-solve-maximum times twelve forecast remain immutable.
This named variant changes only administrative cost prediction. The new forecast
uses actual supervised full-branch wall time, including geometry proof, output,
worker startup and supervision. Each complete branch with n physical solves is
projected to the existing twelve-solve cap as wall_seconds * 12 / n. All verified
event and no-event branches are included; the historical maximum only increases.
The estimate also has a floor of twice the historical maximum single-solve cost.
Forecast = 30 seconds + 1.5 * max(projected full-branch cost, single-solve floor).
This is an empirical scheduling estimate, not a deterministic runtime bound.
Only forecasts <=450 seconds can start. Internal watchdog540 / external600 second
limits are unchanged, with >=650 seconds remaining required before every child.
The original physical driver's Newton retry reserve still receives the actual
maximum single-solve cost, not an averaged value. The maximum is updated after
each verified branch, before scheduling the next one, and never decreases.

Startup fully reads the eight saved pairs (replay, three tolerance siblings and
continuations001-004) with the unchanged original event receipt reader, verifies
their source/artifact hashes and complete original global/local/crossing-ID
screens, and checks the selected continuation chain. Only existing complete
evidence calibrates the estimate; no physical work is repeated for calibration.
All187 inherited sources are retained unchanged. New administrative sources,
this contract and tests form an explicit source-map extension. Frozen CPU
material/thermal/contact/Newton/event integration, tolerance .25, one event and
at most12 solves per branch /25 per pair remain exactly the same.

The three independently bounded branches still form one complete coarse/fine
screened pair. Only its fully read-back final fine state is a continuation parent.
Partial pairs and probe states never become restart parents. Checkpoints persist
across later failure or normal budget exit. Cost, pause, deadline and physical
failure receive separate statuses. There is no unconditional restart of failed
physics, no removal of checks, and no promise of uninterrupted progress through
unknown future events. The old event repair and its sensitivity evidence remain
the qualification of the unchanged physical method.

All output remains ENGINEERING_APPROXIMATE diagnostic evidence with accepted=false
and restart_allowed=false; only explicitly screened continuation is authorized.
Continuation and deposition results remain NOT_RUN until actual receipts exist.

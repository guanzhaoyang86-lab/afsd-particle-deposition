# Continuation using the qualified cancellation-aware visibility bound

Registered before the first physical run. This series starts only from the
independently audited 49.9375 ms / global step 655 accepted checkpoint at
`outputs/visible_local_refinement_20261004/window_004/fine/step_014/accepted.npz`
(SHA256 `396080e1fdd6d9f211fb521d2c44e47ecb19c17c0d327d781b3bb5f61156d4f0`).
All 655 interval records, particle fields, constraints, contact history,
material history, spatial transition and branch origin remain unchanged.

Only the geometric certification algorithm changes. The separate qualification
`outputs/visible_path_tight_20261004/qualification_v1/summary.json`, SHA256
`3d461593a2fab30900a83bcde5b2835b353bec3d66783393b3500ce493b750bc`, binds
26 passing component cases and the saved failed path replay. Its entire source,
input, output and execution-log evidence must be verified before a physical run.
The new bound follows `VISIBLE_SUPPORT_TIGHT_BOUND_CONTRACT_20261004.md` and
keeps the original 256-query / 32-depth caps, strict intersections, numerical
domain rejection and all retained graph edges. It changes no mechanics,
material, heat, contact, geometry, loading or acceptance tolerance. Its rounding
envelope is an explicitly approximate binary64 engineering bound, not a formal
interval-arithmetic proof. The old 50 ms failure and unaccepted candidate remain
immutable; a saved replay alone never promotes them to accepted status.

Each numbered window is 1 ms: 16 x 62.5 us nonrestartable coarse references and
32 x 31.25 us fine increments from the same true parent. All original fresh-map,
full-history energy, finite-shell, scalar/batch and time-screen gates remain.
Local seeds include inherited events, blocked and crossed IDs plus 270, 339 and
341; both branches' raw and retained one-hop neighbors enter comparisons.
The thresholds remain 5 um position, 0.5 K temperature, 10% stress (1 MPa local
floor), 20% plastic strain (1e-5 local floor), and identical strict crossing IDs.
No acceptance is permitted solely because the new path certificate succeeds.

The cost forecast is 30 s + 1.5 * planned steps * max(previous full-stage cost
per qualified candidate, slowest candidate root), plus 1 s per planned step in
the first window for the new geometry certificate. It must be <=450 s. Every
stage retains 540 s internal / 600 s external limits, one CPU thread, the fixed
PySPH commit, and no scaling. A failed coarse stage can supply only its already
qualified contiguous prefix for fine comparison; its original FAIL is retained.
Only an independently audited full window may parent the next numbered window.
At most 100 such separately budgeted windows are registered; a physical or
numerical domain failure stops the series for diagnosis.

New accepted states explicitly identify the path qualification, contract and
component-receipt hashes. A separately implemented readback verifies all saved
fields, energy history, source maps, certificates, publication metadata and
same-parent comparisons. This limited 342-particle continuation is not full
deposition validation: the wider independently initialized domain, replenishment,
spatial/time sensitivity and comparison to the paper remain separate work.

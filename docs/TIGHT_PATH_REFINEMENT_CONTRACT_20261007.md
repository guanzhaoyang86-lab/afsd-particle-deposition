# Bounded continuation after the 51.1875 ms path-budget stop

Registered before new physical evaluations on 2026-10-07, following the user's
request to extend the simulation and observe further particle passage.

The unchanged 342-particle process was extended from 50.9375 to 51.1875 ms,
global step 695, with eight newly accepted fine increments. Independent saved
readback passed, but window 002 remains incomplete. Its ninth fine candidate
failed the 256-query retained-support path check. No failed candidate may be
promoted. The original failure, sources and all prior evidence remain unchanged.

The only starting checkpoint is
`outputs/tight_path_process_v2_20261004/window_002/fine/step_008/accepted.npz`,
SHA256 `857d1fa88b2b74939d9c308306506918611c823847fde13f04684c24324ae408`.
Its independent `window_002/readback/summary.json` has SHA256
`69984500c9337f4f2b4eb4e089755a38290c9b7fa48217e8dafd7f197b932212`.
The parent's qualified last screen is `fine/pair_004_temporal_screen.json`.
All 695 energy records, IDs, fields, constraints, contact/material history,
crossing observations and original branch identity are inherited without reset.

The new controller statically adapts `run_tight_path_process_v2.py`; it calls
the unchanged `visible_tight_path_candidate.candidate`. Only the explicit
parent, provenance links, output series and time schedule change. Each window
uses 16 x 15.625 us coarse references and 32 x 7.8125 us fine increments from
the same parent, adding at most 0.25 ms. At most four windows are registered
(at most 52.1875 ms); each subsequent window requires the preceding complete
window and independent readback. A failure stops the series, with no automatic
retry at another timestep. A fully qualified coarse prefix may supply equal-time
references for fine qualification, without changing the original coarse FAIL.

All original mechanics, constitutive, thermal, geometry and acceptance settings
remain: fixed PySPH/CPU, no scaling, 256 path queries, 32 subdivision levels,
fresh-map and scalar/batch checks, complete-prefix energy, continuous shell,
5 um local position, 0.5 K temperature, 10% stress (1 MPa local floor), 20%
plastic strain (1e-5 floor), and identical crossing IDs. Inherited events,
blocked/crossed IDs and seeds 270/339/341 retain the original one-hop comparison.
Previous TOP contact qualification is retained, with the existing first-window
TOP observation requirement also unchanged.

Before each stage the original measured-cost forecast must be <=450 s; internal
and external limits remain 540/600 s, one numerical thread. A new independent
reader replays the same saved arithmetic and fields with the new explicit
parent and time grid. Controller and reader adaptations are inspected before
execution; prior qualified component source hashes must still match. No claim
of new physics-component or full-regression testing follows from this adaptation.

The limited 4 mm substrate and finite feedstock are unchanged. Additional
accepted time or outlet crossings are not complete deposition validation.
The longer 100 ms objective remains unachieved until actually run and audited.

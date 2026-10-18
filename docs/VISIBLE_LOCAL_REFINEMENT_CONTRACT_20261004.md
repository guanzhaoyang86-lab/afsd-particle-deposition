# Local refinement after the 46.625 ms time comparison failure

Registered before new physical runs. The original TOP window 002 remains
FAILED_WITH_ACCEPTED_PREFIX: its fine candidates 027 and 028 were not published.
At 46.625 ms PID270 has 10.07210375058861% local stress difference against the
same-parent 125 us coarse branch, above the unchanged 10% limit. Global fields
and strict crossing IDs passed. This failure is not rounded down or relabeled.

Use the last true accepted state at 46.500 ms / global step 545:
`outputs/visible_top_event_20261004/window_002/fine/step_026/accepted.npz`, SHA256
`b395158b765378f4dcd64f2fa7a98d236efb8a260be79a385d3acdf71d4be3ce`.
The independent readback of the failed window must confirm this exact prefix
before this controller may begin. All 545 energy records, passive fields,
contact and material history, old branch origin and spatial transition remain.

The first local refinement window is 1 ms, with 16 x 62.5 us nonrestartable
coarse references and 32 x 31.25 us fine increments from this same true parent.
PID270 is included in cumulative local seeds from the start, in addition to all
inherited contact, blocked and crossed IDs. Raw and filtered one-hop neighbors
from both branches enter every equal-time comparison. New fine states publish
only after each pair completes all physical/path gates and global/local/time
screens. Tolerances remain unchanged: 5 um, 0.5 K, 10% stress with 1 MPa local
floor, 20% plastic strain with 1e-5 local floor, and identical strict crossing IDs.

This is a physical time-step refinement with unchanged equations and speeds,
not time scaling or a material/model restart. It tests local continuation from
the accepted state; it does not erase the earlier accumulated temporal error
or establish convergence of the entire 0--46.5 ms prefix. The 250 us TOP-creation
failure remains separate, and larger-step studies without complete references
are NOT_RUN. No OUTER_CYLINDER/BASE event extension is authorized by this file.

The previously qualified TOP candidate implementation is reused unchanged,
including fresh force/thermal/material acceptance, all-history energy checks,
full scalar/batch comparison, continuous finite-shell enclosure, and retained
support edges certified over the whole increment. Parent and all sources are
frozen by hashes. Independent readback must verify actual candidate/publication
arrays, fine/coarse comparisons, complete history and all saved certificates.

Measured full prior-stage cost (including failed pair work) and slowest root
give forecast 30 s + 1.5 * planned steps * max(full cost/qualified candidates,
slowest root), required <=450 s. Each physical stage still has 540 s internal
and 600 s external limits, one CPU thread, unchanged fixed PySPH and parameters.
No GPU or mass/modulus scaling is introduced.

At most 100 numbered 1 ms refined windows can continue, each separately costed
and audited; only a fully completed window with independent PASS can parent
the next number. Failure stops that series with its true accepted prefix saved.
Stop for assessment when newly crossed material develops sustained shoulder
spreading. The 342-particle finite domain cannot certify a full stable production
track; wider substrate, replenishment, spatial sensitivity and paper comparison
remain distinct required work.

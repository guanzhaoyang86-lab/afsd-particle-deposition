# Local subdivision without a cumulative refinement quota

The user explicitly requested removal of the cumulative four-halving stop, then
clarified: fix and test only, do not run. This revision must not launch physical
simulation during repair. The old v3 evidence and199 frozen sources remain intact.

The prepared continuation parent is the complete v3 continuation0023_h2 at
79.6962890625ms /454 steps, SHA
03fefe0f56537c3793be81a1d11b500375b16eec78d5482175295f7b248c2c60.
Its last fully verified macro dt,118.65234375us, is the nominal interval. This
revision does not silently increase it to the earlier474.609375us.

Each nominal interval is represented by ordered pending absolute end times. On
the unchanged, narrowly qualified root multiple-opening failure, subdivide only
the first pending interval at its midpoint. Preserve the original end boundary
and all later pending boundaries. Recompute the entire coarse/fine pair from the
last complete parent. A pair must pass all original receipt/physical screens and
reach its requested boundary before it can be promoted and that boundary removed.
An incomplete pair, failed fine1/fine2 or diagnostic chord cannot become a parent.

There is no campaign-wide or per-interval refinement-count limit. Diagnostic
counters never authorize stopping. The existing midpoint and finer-step times
must remain representable as strictly increasing finite floating-point values;
otherwise record NUMERICAL_TIME_RESOLUTION_EXHAUSTED. This is a numerical domain
limit, not a replacement arbitrary attempt count. After all boundaries of an
interval are completed, the next interval starts with the nominal118.65234375us
again. This local recovery never skips pending material evolution or changes the
original physics or coarse/fine thresholds.

For the prepared starting point, existing fully supervised v3 failures0024_h2,
h3 andh4 can be read and reclassified to seed the three known local subdivisions.
Their SHA-bound ordinary/union/strict/shell/contact records and actual parent must
still qualify. No physical work is replayed merely to rediscover those failures.
They select smaller trial intervals only, never a state or a physical event time.

The original CPU worker, event locator, common graph, material/thermal/contact
providers and complete-pair validation stay frozen. Unknown failures, after-event
failures outside the original classifier, source corruption, pair-screen failure
and the event integrator's own domain bounds remain stopping conditions. Physical
12-solve/1-event branch and25-solve pair limits remain, as do450s cost prediction,
540/600s hard child limits, user pause, and the stationary phase end at220ms.
Removing the administrative counter cannot guarantee progress through all future
physical failures, and this repair does not claim deposition has been achieved.

A future launch requires an explicit wall duration. Only actual launch computes
start+duration once, saves a new immutable configuration and gives that same
absolute deadline to every child. It must not retain the old01:10 deadline or
grant a fresh duration to each worker. This repair creates neither a launch
configuration nor a running process. The inspect command remains read-only.

Validation must include more than four consecutive mock failures followed by
successful completion, local remainder consumption and nominal-step restoration,
same-parent whole-pair rollback, invalid-time/float exhaustion rejection, unchanged
physical worker and pair validation, and launcher tests with process creation
mocked. Actual new physical intervals remain NOT_RUN until separately requested.

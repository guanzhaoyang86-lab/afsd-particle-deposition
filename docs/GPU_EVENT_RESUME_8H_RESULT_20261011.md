# 2026-10-11 eight-hour continuation result

## Stopped session, checked at 18:55 +08:00

The session stopped at 17:22:12 with STOPPED_CHILD_FAILURE; supervisor 49000 no
longer exists. It completed 35 passing stages out of 40 attempts. The latest
complete checkpoint is attempt_0039 at 0.08958959960937454 s, step 599, SHA256
`aedf8bfed9c67f38ec12bf38e1e9ec77e37cbba39914f63064dfedbfb537c474`.
Direct hash read matches both progress and manifest; original pair screen and
crossing IDs pass, supervision exit is zero, source matches, no interruption.
Four supplied particles have crossed (1.35 mg), with zero under the shoulder.

Attempt 40 failed in its first coarse proposal: original Broyden exhausted its
Armijo refresh budget and original Newton fallback raised
`RefinableCandidateError('complete residual line search/evaluation budget failed')`.
Supervision took 66.4487011 s, exited 2, and reports no interruption and unchanged
source inputs. This is an unconverged nonlinear solve, not a wall-time expiry or
geometry-query exhaustion. No candidate from that attempt is a valid parent.
The controller routed it to NON_TEMPORAL_FAILURE_REQUIRES_DIAGNOSIS.

The original CPU adaptive_pair contains same-parent halving for this error class,
but a narrow GPU administrative recovery and real smaller-step result remain
unimplemented/unverified in this session. This status turn did not change frozen
files or restart computation. Any later restoration must reconstruct the final
STOPPED adaptive policy, and investigate why final audit unique/locked files are
30284 versus 30283 expected files before reusing an old equality assertion.
The original authorization cutoff remains 23:41:47; the eight-hour continuous
run was not achieved.

## Live observation at 17:01 +08:00

The session is RUNNING_GPU_CHILD, with supervisor 49000 and child 52572 alive.
There are 24 complete passing stages; attempt 29 is running. The latest complete
stage is attempt_0028 at 0.08436889648437465 s, step 577, checkpoint SHA256
`4fa75fbd2b057ef5200a5e8b973ecb860f8a08e25a4fc94f17e76cbb44135d3c`.
This advances 3.025634766 ms and 50 saved fine steps from the session parent,
past the formerly rejected interval. Original global/local and crossing-ID
checks report PASS. Four supplied particles have crossed (1.35 mg), with zero
under the shoulder; crossing has not yet established deposition.
Current macro timestep has recovered to 474.609375 microseconds. Initial history
preflight reports 765.2767125 seconds. Deadline remains 23:41:47. This status
check did not restart the simulation or change any frozen implementation.

## Authorization and actual start

User request: continue the existing run for eight hours without avoidable interruption.
Authorization window: 2026-10-11 15:41:47 to 23:41:47 +08:00; preparation and verification count.
The detached supervisor started at 16:02:31, PID 49000, output
`outputs/gpu_acceleration_20261007/event_resume_8h_20261011_v1`.
Bootstrap independently reports PASS_DETACHED_JOB_POLICY and kill_on_job_close=false.
The singleton and process-local idle-sleep prevention remain enabled.
At 16:03 the process is alive in INITIAL_FULL_HISTORY_VERIFICATION. There is no
new GPU child or new completed physical stage yet. Full initial history verification
runs once, followed by the continuation using its live immutable history pool.

## Complete starting state and previous stop

The prior 90-minute run ended at 00:55:35 with STOPPED_CHILD_FAILURE, 13 complete
passing stages and a rejected attempt 14. The complete parent is
`event_resume_90m_20261010_v1/attempt_0013/fine_2/NOT_RESTARTABLE_candidate.npz`,
SHA256 `0079328d9ea46263d62cf975299497ea6774136d22e5bdf0ffba5d6cdcb91c97`,
physical time 0.08134326171874973 s, step 527, 1381 particles.
Four supplied particles have crossed (1.35 mg), shoulder occupancy is zero.
Nearest eligible uncrossed PID 1307 is 34.94029918145545 micrometres above the outlet.

Failed attempt 14 had two opening support edges, 1059–1087 and 1111–1243, in one
coarse proposal. Original retained and union routing found both without unresolved
or closing events. The continuous shell certificate passed; the failure used only
2 of 256 path queries. The frozen v5 event integrator only accepts one opening at
a time. The failed proposal is not a continuation parent. Its one physical solve
took 30.7044 s, worker 59.978 s, supervisor 63.324 s. All 30 failed-stage artifacts,
the supervision hashes and the valid parent were independently checked.

## Administrative changes and validation

The exact 636-entry physical implementation remains unchanged. A new strict
administrative classifier recognizes multiple localized opening events using five
SHA-bound reports, including union-edge routing. It permits another real solve
from the same complete parent with half the macro timestep; it does not certify
the rejected path or accept a predicted event state. The initial restored policy
is READY at 237.3046875 microseconds, with the original eight-halving limit,
minimum timestep and cooldown. Unknown, closing, unresolved and hidden-contact
failures do not qualify.

Complete cost history remains bound: 24 completed pairs, maximum supervised pair
590.5055829000012 s, forecast 915.7583743500018 s, maximum actual single solve
113.85019830000238 s. Forecasts no longer end the campaign merely by exceeding a
fixed 900/1200-second ceiling. Per-stage soft budget is
min(3600,max(1260,forecast+60)); hard adds 30 s and outer adds 120 s. All are capped
by the actual session deadline. A limited final best-effort stage is permitted
while the measured single-solve reserve remains. Explicit administrative hard
timeouts with no recorded physical failure can trigger bounded same-parent retry.
Recorded physical/gate errors retain the original conservative handling.
Shoulder occupancy is a milestone rather than a five-millisecond stop condition.
The registered stationary phase and 0.22 s physical target remain scope limits.

Unified component check: PASS_EVENT_RESUME_8H_COMPONENTS, 50 tests (30 controller /
worker / restoration, 20 multi-event policy). Independent review passed both
groups. The actual old saved-state fixture reconstructs 48 frames, 24 completed
cost rows and the exact pending failed attempt, then applies the new recovery.
The root checker verified 1,932 historical artifacts and inherited source closure.
Component receipt: `outputs/gpu_acceleration_20261007/event_resume_8h_components_20261011_v1.json`.
All seven new administrative files and their receipt are frozen after this check.

Launch command:

```powershell
.envs/pysph/python.exe -B scripts/run_gpu_event_resume_8h_20261011.py --output outputs/gpu_acceleration_20261007/event_resume_8h_20261011_v1 --wall-seconds 28800 --max-stages 100000 --detach
```

At this entry, real crossing of the rejected interval remains NOT_RUN in the new
session. Unit tests and a successful background launch are not proof of a new
physical PASS or of shoulder deposition. Existing diagnostic-only qualifications
remain in force.

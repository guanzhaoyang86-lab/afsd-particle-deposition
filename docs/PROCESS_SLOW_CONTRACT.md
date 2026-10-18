# Slow loading / full-startup HHT feasibility scout — preregistered

2026-10-03, before any new trajectory. Classification:
LOCAL_SLOW_LOAD_FULL_STARTUP_FEASIBILITY_NOT_AFSD. This is one bounded
work package under PROCESS_HHT_NEXT_WORK_PACKAGE, not relaxation of the
previous velocity accuracy FAIL and not production authorization.

## Exact physical inputs and provenance

Use the literal27node1cm block, tensor-trapezoid masses totaling.0027kg,
dx=.005m,h=.0065m,rho0=2700kg/m3,E=70GPa,nu=0,G=35GPa,K=70GPa/3,
nonzero S=diag(20000,-10000,-10000)Pa,v=(0,0,.02)m/s,T=400+2*x/.01K,
ep=0 and actual original FKS assumptions from PROCESS_HHT_FIXTURE_CONTRACT.
This is manufactured elastic material provenance, NOT a calibrated AA6061
process. All81 material DOFs remain free; no counterforce or prescribed free
velocity cancels startup. No old checkpoint is altered or reinitialized.

Opposite x planes retain2um initial overlap, radiusdx/2, kn=7e6N/m per face
node, mu=.2, same9nodes each, cp=900J/kg/K,k=167W/m/K,H=.025W/K per face
node, baths410/420K, material friction partition.5 and plastic heat fraction1.
Keep original mass/moduli, room/melt limits, label/ID roles and thermal law.

The NEW named rigid-wall programme is:

    d(t)=1e-7*(1-cos(2*pi*5*t)) m
    xwall_left=-dx/2+2e-6+d(t)
    xwall_right=.01+dx/2-2e-6-d(t)
    w_left=(d'(t),0,+.2) m/s; w_right=(-d'(t),0,+.2) m/s

Tangential wall motion is constant +.2m/s; normal motion is5Hz. No1us jump
or reversal is imposed. The5Hz scale is motivated by300rpm, but these planes
are not a rotating AFSD tool. Choosing a different physical manufactured
boundary is explicit; physical time/mass/moduli are NOT scaled, and old fast
fixture output is not stretched or reused as a solved state.

## Horizons, schedules and costs fixed before execution

Desired next feasibility horizon1ms is initially NOT_RUN: old independently
measured RK4 cost1.333s per2us for both grids predicts666.5s to1ms, over the
60s reference cap. Do not start that reference or declare its accuracy known.
First inspect full startup to20us (10x previous window), or its first genuine
unsupported branch. This remains a startup scout, not a material-transport run.

Candidate schedules for a successful ordinary20us path are fixed:

|Schedule|0–2us|2–4us|4–20us|Total steps|
|---|---:|---:|---:|---:|
|S0|250ns|500ns|1000ns|28|
|S1|125ns|250ns|500ns|56|
|S2|62.5ns|125ns|250ns|112|

Separate invocation for each schedule, internal55s/external60s. One bounded
250ns cost pilot may precede the scout (internal55/external60). Launch a
schedule only if2*planned_steps*max(measured_whole_step,.01)+10<=50s.
Original measured .139s pilot suggests41.1s for S2, not a guarantee. Actual
cost includes all currentNNPS/material/root/contact/heat/guards/output.
Stop a schedule on first numerical/physical/domain failure; do not halve its
time step or expand24nonlinear evaluations to get past that failure. After
S0 fails an unimplemented physical branch, later schedules are NOT_RUN unless
they are needed for a separately registered event-location check, not a hunt
for an accepted microstep. Only all3completed20us permits sensitivity claims.

## Independent reference and event semantics

Reuse independent all-pair formulas/literal inputs only; no tested HHT/path,
mechanics/material/thermal functions may generate the oracle. Independently
code actual wall programme, full2D tangential Coulomb force, Jaumann ODE,
continuity and heat. Fixed RK4 grids7.8125ns and3.90625ns (unscaled) predict
~13.33s for20us on ordinary smooth paths. Common snapshots every62.5ns;
save all accepted states and source/config hashes. Internal55/external60s
for both references combined; stop on the deadline.

The old signed-z slip>.05m/s rule was a fast-fixture domain guard, NOT a
physical stop criterion; it is NOT inherited. A nonzero vector slip whose
z component changes sign is valid curved/reverse sliding, not automatic FAIL.
Record signed z and complete tangential vector separately. Physical branch
boundaries are normal opening, actual full-vector zero slip, yielding, or
invalid geometry/temperature/density. Check intermediate stages and dense
RK4 cubic paths as well as endpoints. A candidate crossing a boundary is
REJECTED, never the accepted physical endpoint. A scalar signed-z root can
certify full zero only with independently checked transverse symmetry at
512eps*max(1, relevant speed) roundoff, not by prescribing transverse motion.

For event localization, at most12bisections/halvings from the last accepted
state are allowed. Keep every rejected candidate's reason/margins/time/state.
Never fabricate a post-event branch or clip negative Fn/Q/speed. If a true
root cannot be certified, report EVENT_APPROACH_BLOCKED with last accepted
vector slip, gap, yield margin and rejected interval; do not call it STOP.
Two reference grids must agree on class/contact and event bracket within
max(one coarse reference step, their retained bracket widths); otherwise
event identity/time stays BLOCKED. Shared accepted-prefix state comparison
must use identical sample times, not different rejected endpoints.

Reference shared-time differences must satisfy the existing1% oracle bounds:
RMSx1e-11m,v2e-6m/s,fullstress70Pa,maxT1e-9K,faceforce.007N,
impulse1.4e-9Ns,facewallwork/rawQ2e-9J,reservoirheat2e-10J.
Complete20us HHT finest accuracy retains existing absolute bounds1e-9m,
2e-4m/s,7000Pa,1e-7K,.7N,1.4e-7Ns,2e-7J,2e-7J,2e-8J
and all errors decrease by factor<=.7 above10x reference uncertainty.
Old fast-fixture FAIL is separate; no tolerances are fitted to new results.

## Driver scope and transactions

Use unmodified hht_coupled current-geometry elastic FKS branch and raw residual
gates; any plastic response is an unimplemented scope boundary, not zeroed.
Use exact raw physical forces, algorithmic work separate from heat, actual
saved-temperature BE balance, ID/mass integrity and atomic publication.
New slow inputs and physical boundary must be independently unit-checked;
all previous4067frozen entries remain intact except the3document indices.

The startup scout may diagnose a candidate interval containing a reference
domain boundary, but must not publish a state beyond the reference's certified
prefix as validated. Save such a candidate explicitly REJECTED/UNREFERENCED;
do not count it as a successful process step. Numerical failure before a
physical event is classified separately. No first-stage interpolation is
evidence of the event cause without the independent reference.

Retain original density<=.005relative, overlap<=5um, solid temperature/finite
checks and node momentum/work/BE roundoff gates. This distinct longer-horizon
fixture permits total translation<=250um (wall travel at1ms plus50um margin),
not old AFSD's50um guard; maximum relative pair-length change<=.005 guards
the manufactured small-strain scope. No value is clipped to enforce a guard.

At most two evidence-based implementation repairs, retain every FAIL. If a
contact branch or accuracy/cost overlap blocks the scout, finish the diagnostic
package with that evidence and identify the exact needed next implementation.
Do not respond by adding equivalent microsteps, deleting inertia, clearing
startup velocity or moving to a different backend/material silently.

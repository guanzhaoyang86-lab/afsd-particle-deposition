# Curved-contact scout: independent causal and cost review

2026-10-03. Read-only state/source inspection and independent algebra only.
No ODE integration, solver step, time refinement, restart or physical state
publication was performed by this review. Local27node manufactured benchmark;
not AFSD transport and not a completed process feasibility result.

## Evidence inspected

- `outputs/process_curved_20261003/saved_red_initial/report.json`:
  `PASS_PRIVATE_SAVED_RED_REGRESSION_NOT_TRAJECTORY`; 12 evaluations, velocity
  residual2.78861e-10m/s. This is the same saved [2,2.5]us interval that the
  old endpoint-chord method rejected; old failure evidence remains intact.
- `outputs/process_curved_20261003/S0_initial/report.json`, accepted history,
  `stage_013/accepted.npz` and `REJECTED_UNREFERENCED.npz`.
- Independent `logs/process_slow_20261003/reference_repair1/` report and both
  shared-grid/fine histories. Original reference source and histories unchanged.

Auditing script `logs/process_curved_20261003/audit_causal_state.py` reads NPZ
directly for diagnostics, including the rejected envelope, without invoking a
solver or treating it as restartable. It recomputes instantaneous geometry,
elastic force, FKS margin and contact slip with the independent all-pair
reference algebra, reconstructs the rejected interval's slip polynomial
independently, and records all input hashes. Output is
`causal_state_audit.json` and stdout/elapsed/exit status in
`causal_state_audit.log`. Audit elapsed .044430s; exit0. Direct inspection of a
REJECTED NPZ does not override its `restart_allowed=false` envelope.

## The actual new stop condition

The new S0 run published **13 local steps through5us**, including one1us step
[4,5]us. The next [5,6]us joint candidate returns from the nonlinear solver.
It is rejected by the independent-reference-prefix publication guard, because
its endpoint6us exceeds the currently certified common prefix5.9375us.

The stored rejected candidate satisfies its raw mechanical root: independently
recomputed `h*max|R/M|=3.43364e-10m/s`. This is not a convergence-budget failure
and not a wall-time timeout. It does not make the candidate a valid physical
endpoint. The parent remains the saved5us accepted state; the6us envelope is
explicitly REJECTED and non-restartable.

The numerical interval contains **no zero slip in its own HHT path**. Its
smallest full-vector slip is .0226522082m/s at PID4, fraction1; PID22 agrees.
For PID4 the independently reconstructed tangent-z polynomial is

    uz(s)=-.107901655751284 + .085431597468668*s
          -.000182149945850*s^2, 0<=s<=1.

It remains negative throughout [5,6]us. Transverse motion is measured at
roundoff, not imposed. At the independent reference's last pre-event time,
5.974655151367073us, this rejected interpolant still has
uz=-.0248083430m/s. Therefore the result cannot be described as resolving,
splitting, sticking at, reversing through, or validating the true reference
event near5.97466us. It is a delayed coarse numerical trajectory, not evidence
that the reference's event can be ignored.

## Same-time physical discrepancy before the event

At exactly5us, compare the accepted state to the independent reference at5us,
using the unchanged trapezoid reference-volume weights. No different failure
endpoints are compared:

|Field|Measured S0 error at5us|
|---|---:|
|weighted RMS position|2.11583e-8m|
|weighted RMS velocity|.0132493766m/s|
|weighted RMS full stress|408995.29Pa|
|maximum temperature|4.96245e-7K|
|maximum density|.0248724kg/m3|

These are **coarse S0 diagnostics**, not a completed three-grid convergence
study or a claim that the preregistered finest-grid accuracy gate has been
run. They are nevertheless much larger than root noise or independent-oracle
uncertainty. A root/ledger PASS must not be substituted for dynamic accuracy.

At5us, PID4 slip is -.107901656m/s in S0 versus -.0995602163m/s in the reference.
Its independent instantaneous internal +z force is4.19956N in S0 versus
5.41261N in the reference. Both are still sliding. These differences support
a temporal-response/phase limitation in the current coarse schedule; this
review does not isolate a single guilty substep, integrator coefficient or
material-rate discretization. It does not authorize tuning a coefficient
until the event time matches.

## True reference event and static-feasibility evidence

The unchanged independent two-grid event approach intervals overlap:

- coarse [5.974655151367281,5.974670410156344]us;
- fine [5.974655151367073,5.974662780761604]us.

At the last accepted independent pre-event state, IDs4/22 have full slip
1.2623e-7m/s and transverse velocity around1e-18m/s. PID4 internal +z force is
4.67484334N, Coulomb capacity2.73923931N, giving instantaneous required holding
force/capacity1.706621. This supports infeasible sticking and forward passage
into reverse sliding once the actual zero is located. It is still pre-event
evidence, not an exact event-state solve or a simulated post-event branch.

No opening or yield precedes it: minimum overlap1.61882um, minimum yield
margin313.363MPa. On the rejected6us coarse state, the instantaneous holding
ratio is also >1 (1.76232), but its velocity is not zero. That ratio alone
cannot justify snapping it to wall speed or flipping its traction early.

## Evidence-ranked conclusions

**CONFIRMED**

- The original nonzero-turning chord failure was overcome on the same saved
  interval using a mechanical impulse correction, not only a heat substitution.
- Current accepted local extent is5us, not6us and not an AFSD process.
- Immediate publication rejection is missing certified post-event reference.
- The rejected coarse path has no zero even where the fine reference reaches
  one; same-time dynamic discrepancies already exist at5us.
- The prescribed guard retained the parent and rejected candidate separately.

**SUPPORTED, not independently isolated to one term**

- Coarse temporal response and accumulated numerical history delay the event.
  Tightening the nonlinear root alone cannot explain/remove the measured gap.
- The first true reference event should restart sliding with opposite direction
  after an actual feasibility check, rather than stick indefinitely.

**RULED_OUT as descriptions of this result**

- A resolved zero-slip/reversal event, a completed20us run, a time-convergence
  demonstration, or proof of accurate process-time transport.
- Claiming that removing the reference-prefix guard would solve the physical
  problem; it would only publish an unverified candidate across a known event.
- Claiming that the present stop is the old1% chord veto or a55s timeout.

**NOT_RUN / unresolved**

- Independent post-event continuation; a located joint event/restart variant;
  S1/S2/common20us accuracy, spatial resolution and1ms economics.
- Positive plastic coupling, full finite rotating square bore/actuator/inlet,
  actual AFSD shoulder crossing or three-case production.

## Measured cost, without an unjustified production extrapolation

The complete S0 command reports8.710573s; trajectory body8.671028s. Sum of
accepted nonlinear-stage times6.269790s, maximum accepted root15 evaluations.
Peak process working set81,080,320bytes (77.324MiB). The largest accepted
step is1us, costing .769345s for its root, before per-stage publication cost.
The full S0 elapsed also includes the final rejected candidate and output.

Two rough **local27node**1ms extrapolations illustrate the current cost risk:
holding the last1us root cost gives769s; using the whole observed scout rate
gives1742s. Both exceed600s, so neither licenses a1ms job. They are not bounds,
not costs for full AFSD geometry, and do not prove later larger accurate steps
are impossible. The current measurements establish neither an accuracy/cost
overlap nor stable post-event step growth. GPU/parallel optimization is not
the next missing scientific task.

## Recommended next bounded package

First extend the independent oracle through its **located** zero with correct
static feasibility and reverse sliding, preserving x/v/S/rho/T and raw work.
Use a new source/result, do not overwrite this pre-event oracle. Validate the
event machinery with the independent1D reversal values and nonzero-stress
material state, then attempt one preregistered short post-event window within
the existing60s reference cap. Stop if another unsupported branch appears.
This supplies the missing evidence needed to judge any solver continuation.

In parallel, use the saved current5us/6us data to define an accuracy/cost
decision for the actual transient, not a campaign of equivalent smaller steps.
The coupled event must be located on the solver's own evolving trajectory;
do not impose the reference event time, snap velocity or move material to
make it agree. Require independently measured event timing/impulse/heat and
same-time response, with tolerances fixed before new tests. The S0 discrepancy
must remain visible even if a later event implementation works. If a bounded
set of registered schedules cannot deliver useful accuracy/cost overlap,
report that result and design the next named numerical/model alternative;
do not keep refining until a green number appears or erase physical inertia.

No further trajectory was run by this review, and no old threshold, source,
accepted state or failure log was changed.
